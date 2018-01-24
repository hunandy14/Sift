#include <iostream>
#include <vector>
using namespace std;

#include <opencv2/opencv.hpp>
using namespace cv;

#include "xform.hpp"
#include "utils.hpp"

/********************************** Structures *******************************/
struct ransac_data {
	void* orig_feat_data; //保存此特徵點的feature_data域的以前的值.  
	int sampled;          //標識位，值爲1標識此特徵點是否被選爲樣本.
};
/******************************* Defs and macros *****************************/
/* RANSAC error tolerance in pixels
RANSAC算法的容錯度 對於匹配點對<pt,mpt> 以及變換矩陣H
如果pt經H變換後的點和mpt之間的距離的平方小於RANSAC_ERR_TOL，則可將其加入當前一致集 */
#define RANSAC_ERR_TOL 3
/** pessimistic estimate of fraction of inlers for RANSAC */
#define RANSAC_INLIER_FRAC_EST 0.25 //內點數目佔樣本總數目的百分比的最小值.
/** estimate of the probability that a correspondence supports a bad model */
#define RANSAC_PROB_BAD_SUPP 0.10 //一個匹配點對支持錯誤模型的概率（不知道是幹什麼用的）.
/* extracts a feature's RANSAC data */
//定義了一個帶參數的函數宏，用來提取參數feat中的feature_data成員並轉換爲ransac_data格式的指針.
#define feat_ransac_data(feat) ((ransac_data*)(feat)->feature_data)
/************************* Local Function Prototypes *************************/
static int get_matched_features(Feature*, int, Feature***);
static int calc_min_inliers(int, int, float, float);
static Feature** draw_ransac_sample(Feature**, int, int);
static void extract_corresp_pts(Feature**, int, fpoint**, fpoint**);
static vector<float> lsq_homog(fpoint *, fpoint *, int n);
static int find_consensus(Feature**, int, vector<float>&, float, Feature***);
static void release_mem(fpoint *, fpoint *, Feature**);
/********************** Functions prototyped in model.h **********************/
vector<float> ransac_xform(
	Feature *features, int n,
	int m, float p_badxform, float err_tol,
	Feature*** inliers, int &n_in)
{
	srand((unsigned int)time(NULL));
	Feature **matched = nullptr, **sample = nullptr;
	Feature **consensus = nullptr, **consensus_max = nullptr;
	ransac_data* rdata = nullptr;
	fpoint *pts = nullptr, *mpts = nullptr;
	vector<float> M;
	float p, in_frac = (float)RANSAC_INLIER_FRAC_EST;
	int i, nm, in_min, k = 0, in = 0, in_max = 0;

	// nm 幾個點是有匹配的.
	nm = get_matched_features(features, n, &matched);
	if(nm < m) {
		cout << "Warning: not enough matches to compute xform" << endl;
		goto end;
	}

	//計算保證RANSAC最終計算出的轉換矩陣錯誤的概率小於p_badxform所需的最小內點數目.
	in_min = calc_min_inliers(nm, m, (float)RANSAC_PROB_BAD_SUPP, p_badxform);
	//當前計算出的模型的錯誤概率,內點所佔比例in_frac越大，錯誤概率越小；迭代次數k越大，錯誤概率越小.
	p = pow(1.f - pow(in_frac, m), k);

	int testC = 0;
	while(p > p_badxform) {
		// 從樣本集matched中隨機抽選一個RANSAC樣本(即一個4個特徵點的數組)，放到樣本變量sample中.
		sample = draw_ransac_sample(matched, nm, m);
		// 從樣本中獲取特徵點和其對應匹配點的二維座標，分別放到輸出參數pts和mpts中.
		extract_corresp_pts(sample, m, &pts, &mpts);
		// 計算將m個點的數組pts變換爲mpts的矩陣，返回變換矩陣給M.
		M = lsq_homog(pts, mpts, m);
		if(M.empty())
			goto iteration_end;
		//給定特徵點集，變換矩陣，誤差函數，計算出當前一致集consensus，返回一致集中元素個數給in  
		in = find_consensus(matched, nm, M, err_tol, &consensus);
		//cout << "in=" << in << endl;
		//若當前一致集大於歷史最優一致集，即當前一致集爲最優，則更新最優一致集consensus_max  
		if(in > in_max) {
			if(consensus_max) { //若之前有最優值，釋放其空間.
				delete[] consensus_max;
			}
			consensus_max = consensus; //更新最優一致集.
			in_max = in;
			in_frac = (float)in_max / (float)nm; //最優一致集中元素個數佔樣本總個數的百分比.

			testC++;
		}
		else { //若當前一致集小於歷史最優一致集，釋放當前一致集.
			delete[] consensus;
		}
		M.clear();

	iteration_end:
		release_mem(pts, mpts, sample);
		p = pow(1.f - pow(in_frac, m), ++k);
	}


	// 根據最優一致集計算最終的變換矩陣.
	/* calculate final transform based on best consensus set */
	// 若最優一致集中元素個數大於最低標準，即符合要求.
	if(in_max >= in_min) {
		//從最優一致集中獲取特徵點和其對應匹配點的二維座標，分別放到輸出參數pts和mpts中.
		extract_corresp_pts(consensus_max, in_max, &pts, &mpts);
		//調用參數中傳入的函數xform_fn，計算將in_max個點的數組pts變換爲mpts的矩陣，返回變換矩陣給M.
		M = lsq_homog(pts, mpts, in_max);
		/***********下面會再進行一次迭代**********/
		//根據變換矩陣M從樣本集matched中計算出一致集consensus，返回一致集中元素個數給in.
		in = find_consensus(matched, nm, M, err_tol, &consensus);
		M.clear();
		release_mem(pts, mpts, consensus_max);
		//從一致集中獲取特徵點和其對應匹配點的二維座標，分別放到輸出參數pts和mpts中.
		extract_corresp_pts(consensus, in, &pts, &mpts);
		fpoint r;
		for(int a = 0; a < in - 1; a++) {
			for(int b = a + 1; b < in; b++) {
				if(pts[a].x > pts[b].x) {
					r = pts[a];
					pts[a] = pts[b];
					pts[b] = r;
					r = mpts[a];
					mpts[a] = mpts[b];
					mpts[b] = r;
				}
			}
		}
		for(int a = 0; a < in - 1; a++) {
			for(int b = a + 1; b < in; b++) {
				if(pts[a].y > pts[b].y) {
					r = pts[a];
					pts[a] = pts[b];
					pts[b] = r;
					r = mpts[a];
					mpts[a] = mpts[b];
					mpts[b] = r;
				}
			}
		}
		//調用參數中傳入的函數xform_fn，計算將in個點的數組pts變換爲mpts的矩陣，返回變換矩陣給M.
		M = lsq_homog(pts, mpts, in);
		if(inliers) {
			*inliers = consensus; //將最優一致集賦值給輸出參數：inliers，即內點集合.
			consensus = NULL;
		}
		if(n_in == 0)
			n_in = in; //將最優一致集中元素個數賦值給輸出參數：n_in，即內點個數.
		release_mem(pts, mpts, consensus);
	} //in_max >= in_min

	  //沒有計算出符合要求的一致集 .
	else if(consensus_max) {
		if(inliers)
			*inliers = NULL;
		if(n_in == 0)
			n_in = 0;
		delete[] consensus_max;
	}
end:
	//RANSAC算法結束：恢復特徵點中被更改的數據域feature_data，並返回變換矩陣M.
	for(i = 0; i < nm; i++) {
		//利用宏feat_ransac_data來提取matched[i]中的feature_data成員並轉換爲ransac_data格式的指針.
		rdata = feat_ransac_data(matched[i]);
		//恢復feature_data成員的以前的值
		matched[i]->feature_data = rdata->orig_feat_data;
		delete rdata;
	}
	delete[] matched;
	return M;
}

/************************ Local funciton definitions *************************/
// 獲取匹配點.
Feature* get_match(Feature* feat) {
	return feat->fwd_match;
}

// 變換矩陣
vector<float> lsq_homog(fpoint * pts, fpoint * mpts, int n) {
	vector<float> H(9, 0.f);
	CvMat *A, *B, X;
	float x[9];

	int i;

	if(n < 4) {
		cout << "Warning: too few points in lsq_homog()" << endl;
		return{};
	}

	A = cvCreateMat(2 * n, 8, CV_32FC1);
	B = cvCreateMat(2 * n, 1, CV_32FC1);
	X = cvMat(8, 1, CV_32FC1, x);
	cvZero(A);
	for(i = 0; i < n; i++) {
		cvmSet(A, i, 0, pts[i].x);
		cvmSet(A, i + n, 3, pts[i].x);
		cvmSet(A, i, 1, pts[i].y);
		cvmSet(A, i + n, 4, pts[i].y);
		cvmSet(A, i, 2, 1.0);
		cvmSet(A, i + n, 5, 1.0);
		cvmSet(A, i, 6, (-pts[i].x * mpts[i].x));
		cvmSet(A, i, 7, (-pts[i].y * mpts[i].x));
		cvmSet(A, i + n, 6, (-pts[i].x * mpts[i].y));
		cvmSet(A, i + n, 7, (-pts[i].y * mpts[i].y));
		cvmSet(B, i, 0, mpts[i].x);
		cvmSet(B, i + n, 0, mpts[i].y);
	}

	cvSolve(A, B, &X, CV_SVD);

	H[8] = 1.f;
	for(i = 0; i < 8; i++) {
		H[i] = X.data.fl[i];
	}

	cvReleaseMat(&A);
	cvReleaseMat(&B);

	return H;
}

// 取出幾個點是有匹配的.
int get_matched_features(Feature *features, int n, Feature ***matched) {
	Feature **_matched;
	ransac_data* rdata;
	int i, m = 0;

	_matched = new Feature*[n];
	if(!_matched) {
		cout << "memory error" << endl;
	}
	else {
		for(i = 0; i < n; i++) {
			if(get_match(features + i)) {
				rdata = new ransac_data{};
				rdata->orig_feat_data = features[i].feature_data;
				_matched[m] = features + i;
				_matched[m]->feature_data = rdata;
				m++;
			}
		}
	}
	*matched = _matched;
	return m;
}

float _log_factorial(int n, int m) {
	float f = 0.f;
	for(int i = n; i <= m; i++){
		f += log((float)i);
	}
	return f;
}
float log_factorial(int n) {
	float f = 0.f;
	for(int i = 1; i <= n; i++){
		f += log((float)i);
	}
	return f;
}
//計算保證RANSAC最終計算出的轉換矩陣錯誤的概率小於p_badxform所需的最小內點數目.
int calc_min_inliers(int n, int m, float p_badsupp, float p_badxform) {
	float pi = 0.f, sum = 0.f;
	int i = 0, j = 0;
	float log_n = log_factorial(n);
	float log_nm = log_n - _log_factorial((n - m), n);
	vector<float> temp(n + 1, 0.f);
	float total_s = 0.f;
	for(i = m + 1; i <= n; i++) {
		pi = (float)(i - m) * log(p_badsupp) + (float)(n - i + m) * log(1.f - p_badsupp) + log_nm - log_factorial(i - m) - (log_n - _log_factorial(n - i, n));
		temp[i] = exp(pi);
		total_s += temp[i];
	}
	sum = total_s;
	for(j = m + 1; j <= n; j++) {
		sum -= temp[j];
		if(sum < p_badxform)
			break;
	}
	temp.clear();
	return j;
}

// 從樣本集matched中隨機抽選一個RANSAC樣本(即一個4個特徵點的數組)，放到樣本變量sample中.
Feature** draw_ransac_sample(Feature** features, int n, int m) {
	Feature** sample, *feat;
	struct ransac_data* rdata;
	int i, x;

	for(i = 0; i < n; i++) {
		rdata = feat_ransac_data(features[i]);
		rdata->sampled = 0;
	}

	sample = new Feature*[m];
	for(i = 0; i < m; i++) {
		do {
			x = rand() % n;
			feat = features[x];
			rdata = feat_ransac_data(feat);
		} while(rdata->sampled);
		sample[i] = feat;
		rdata->sampled = 1;
	}

	return sample;
}
//從一致集中獲取特徵點和其對應匹配點的二維座標，分別放到輸出參數pts和mpts中.
void extract_corresp_pts(Feature** features, int n, fpoint** pts, fpoint** mpts) {
	// 提取相應的資訊出來.
	Feature *match = nullptr;
	fpoint *_pts, *_mpts;

	_pts = new fpoint[n];
	_mpts = new fpoint[n];

	for(int i = 0; i < n; i++) {
		// 獲取這一點，對應的匹配點.
		match = get_match(features[i]);
		if(!match)
			cout << "feature does not have match" << endl;
		else if(match) {
			// 原本寫法.
			//_pts[i] = features[i]->img_pt;
			//_mpts[i] = match->img_pt;
			// 無偏移的方式.
			_pts[i] = fpoint(features[i]->rX(), features[i]->rY());
			_mpts[i] = fpoint(match->rX(), match->rY());
		}
	}
	*pts = _pts;
	*mpts = _mpts;
}

/* 找出一致集. */
fpoint persp_xform_pt(fpoint pt, vector<float> &T)
{
	CvMat* _T;
	_T = cvCreateMat(3, 3, CV_32FC1);
	for(int i = 0; i < 9; i++) {
		_T->data.fl[i] = T[i];
	}
	CvMat XY, UV;
	float xy[3] = { pt.x, pt.y, 1.f }, uv[3] = { 0.f };

	cvInitMatHeader(&XY, 3, 1, CV_32FC1, xy, CV_AUTOSTEP);
	cvInitMatHeader(&UV, 3, 1, CV_32FC1, uv, CV_AUTOSTEP);
	cvMatMul(_T, &XY, &UV);

	fpoint rslt = fpoint(uv[0] / uv[2], uv[1] / uv[2]);

	cvReleaseMat(&_T);
	return rslt;
}
float homog_xfer_err(fpoint pt, fpoint mpt, vector<float> &H) {
	fpoint xpt = persp_xform_pt(pt, H);
	return sqrt(dist_sq_2D(xpt, mpt));
}
//對於給定的模型和錯誤度量函數，從特徵點集和中找出一致集.
int find_consensus(Feature** features, int n, vector<float> &M, float err_tol, Feature*** consensus) {
	// 這裡帶入的 features 已經過濾過了一定有匹配點
	Feature** _consensus;
	Feature* match;
	fpoint pt, mpt;
	float err;
	int i, in = 0;
	_consensus = new Feature*[n];
	for(i = 0; i < n; i++) {
		match = get_match(features[i]);
		if(!match) {
			cout << "feature does not have match" << endl;
		}
		else if(match) {
			// 原本寫法.
			//pt = features[i]->img_pt;
			//mpt = match->img_pt;
			// 無偏移的方式.
			pt = fpoint(features[i]->rX(), features[i]->rY());
			mpt = fpoint(match->rX(), match->rY());
			//與變形後的點算距離.
			err = homog_xfer_err(pt, mpt, M);

			if(err <= err_tol) {
				_consensus[in++] = features[i];
			}
		}
	}
	*consensus = _consensus;
	return in;
}

// 釋放內存.
void release_mem(fpoint *pts1, fpoint *pts2, Feature** features) {
	delete[] pts1;
	delete[] pts2;
	if(features)
		delete[] features;
}