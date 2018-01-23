#include <iostream>
#include <time.h>
#include <vector>
#include <opencv2/opencv.hpp>

#include "xform.hpp"
#include "utils.hpp"

using namespace std;
using namespace cv;
/************************* Local Function Prototypes *************************/
static vector<float> lsq_homog(fpoint *, fpoint *, int n);
static float homog_xfer_err(fpoint pt, fpoint mpt, vector<float> &H);
static fpoint persp_xform_pt(fpoint pt, vector<float> &T);
static Feature* get_match(Feature*);
static int get_matched_features(Feature*, int, Feature***);
static int calc_min_inliers(int, int, float, float);
static float _log_factorial(int, int);
static float log_factorial(int);
static Feature** draw_ransac_sample(Feature**, int, int);
static void extract_corresp_pts(Feature**, int, fpoint**, fpoint**);
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

	// 幾個點是有匹配的
	nm = get_matched_features(features, n, &matched);
	cout << "nm = " << nm << endl;
	if (nm < m) {
		cout << "Warning: not enough matches to compute xform" << endl;
		goto end;
	}

	//计算保证RANSAC最终计算出的转换矩阵错误的概率小于p_badxform所需的最小内点数目  
	in_min = calc_min_inliers(nm, m, (float)RANSAC_PROB_BAD_SUPP, p_badxform);
	//当前计算出的模型的错误概率,内点所占比例in_frac越大，错误概率越小；迭代次数k越大，错误概率越小 
	p = pow(1.f - pow(in_frac, m), k);

	int testC = 0;
	while (p > p_badxform)
	{
		//从样本集matched中随机抽选一个RANSAC样本(即一个4个特征点的数组)，放到样本变量sample中  
		sample = draw_ransac_sample(matched, nm, m);
		//从样本中获取特征点和其对应匹配点的二维坐标，分别放到输出参数pts和mpts中  
		extract_corresp_pts(sample, m, &pts, &mpts);
		//调用参数中传入的函数xform_fn，计算将m个点的数组pts变换为mpts的矩阵，返回变换矩阵给M.

		cout << "## x=" << pts->x << ", y=" << pts->y << endl;

		M = lsq_homog(pts, mpts, m);

		if (M.empty())
			goto iteration_end;

		//给定特征点集，变换矩阵，误差函数，计算出当前一致集consensus，返回一致集中元素个数给in  
		in = find_consensus(matched, nm, M, err_tol, &consensus);
		//cout << "in=" << in << endl;
		//若当前一致集大于历史最优一致集，即当前一致集为最优，则更新最优一致集consensus_max  
		if (in > in_max)
		{
			if(consensus_max) { //若之前有最优值，释放其空间  
				delete[] consensus_max;
			}
			consensus_max = consensus; //更新最优一致集
			in_max = in;
			in_frac = (float)in_max / (float)nm; //最优一致集中元素个数占样本总个数的百分比  

			testC++;
		}
		else { //若当前一致集小于历史最优一致集，释放当前一致集  
			delete[] consensus;
		}
		M.clear();

	iteration_end:
		release_mem(pts, mpts, sample);
		p = pow(1.f - pow(in_frac, m), ++k);
	}
	cout << "testC=" << testC << ", in_max=" << in_max;
	cout << endl;

	// 根据最优一致集计算最终的变换矩阵  
	/* calculate final transform based on best consensus set */  
	// 若最优一致集中元素个数大于最低标准，即符合要求  .
	if (in_max >= in_min) {
		//从最优一致集中获取特征点和其对应匹配点的二维坐标，分别放到输出参数pts和mpts中  .
		extract_corresp_pts(consensus_max, in_max, &pts, &mpts);
		//调用参数中传入的函数xform_fn，计算将in_max个点的数组pts变换为mpts的矩阵，返回变换矩阵给M  .
		M = lsq_homog(pts, mpts, in_max);

		/***********下面会再进行一次迭代**********/  
		//根据变换矩阵M从样本集matched中计算出一致集consensus，返回一致集中元素个数给in.
		cout << "==========最後===============" << in << endl;
		in = find_consensus(matched, nm, M, err_tol, &consensus);
		cout << "	lsq_homog in=" << in << endl;



		/*測試*/
		int eff_num = 0;
		for(size_t i = 0; i < in; i++) {
			if(consensus[i]->fwd_match) {
				eff_num++;
			}
		}
		cout<< "eff_num=" << eff_num << endl;


		M.clear();
		release_mem(pts, mpts, consensus_max);
		//从一致集中获取特征点和其对应匹配点的二维坐标，分别放到输出参数pts和mpts中  .
		extract_corresp_pts(consensus, in, &pts, &mpts);
		fpoint r;
		for (int a = 0; a < in - 1; a++)
		{
			for (int b = a + 1; b < in; b++)
			{
				if (pts[a].x > pts[b].x)
				{
					r = pts[a];
					pts[a] = pts[b];
					pts[b] = r;
					r = mpts[a];
					mpts[a] = mpts[b];
					mpts[b] = r;
				}
			}
		}
		for (int a = 0; a < in - 1; a++)
		{
			for (int b = a + 1; b < in; b++)
			{
				if (pts[a].y > pts[b].y)
				{
					r = pts[a];
					pts[a] = pts[b];
					pts[b] = r;
					r = mpts[a];
					mpts[a] = mpts[b];
					mpts[b] = r;
				}
			}
		}
		//调用参数中传入的函数xform_fn，计算将in个点的数组pts变换为mpts的矩阵，返回变换矩阵给M  .
		M = lsq_homog(pts, mpts, in);
		if (inliers)
		{
			*inliers = consensus; //将最优一致集赋值给输出参数：inliers，即内点集合 .
			consensus = NULL;
		}
		if (n_in == 0)
			n_in = in; //将最优一致集中元素个数赋值给输出参数：n_in，即内点个数 .
		release_mem(pts, mpts, consensus);
	} //in_max >= in_min


	//没有计算出符合要求的一致集 .
	else if (consensus_max){
		if (inliers)
			*inliers = NULL;
		if (n_in == 0)
			n_in = 0;
		delete[] consensus_max;
	}
end:
	//RANSAC算法结束：恢复特征点中被更改的数据域feature_data，并返回变换矩阵M .
	for (i = 0; i < nm; i++)
	{
		//利用宏feat_ransac_data来提取matched[i]中的feature_data成员并转换为ransac_data格式的指针 .
		rdata = feat_ransac_data(matched[i]);
		//恢复feature_data成员的以前的值
		matched[i]->feature_data = rdata->orig_feat_data;
		delete rdata;
	}
	delete[] matched;
	return M;
}

/************************ Local funciton definitions *************************/
// 輸入兩點匯出變換矩陣; n 是最少點數輸入是 4.
vector<float> lsq_homog(fpoint * pts, fpoint * mpts, int n)
{
	vector<float> H(9, 0.f);
	CvMat *A, *B, X;
	float x[9];

	int i;

	if (n < 4)
	{
		cout << "Warning: too few points in lsq_homog()" << endl;
		return{};
	}

	A = cvCreateMat(2 * n, 8, CV_32FC1);
	B = cvCreateMat(2 * n, 1, CV_32FC1);
	X = cvMat(8, 1, CV_32FC1, x);
	cvZero(A);
	for (i = 0; i < n; i++)
	{
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
	for (i = 0; i < 8; i++)
	{
		H[i] = X.data.fl[i];
	}

	cvReleaseMat(&A);
	cvReleaseMat(&B);

	vector<float> Hmat(9, 0.f);
	for(size_t j=0, idx=0; j < 3; j++) {
		for(size_t i=0; i < 3; i++, idx++) {
			//cout << H[j*3+i] << ", ";
		} //cout << endl;
	} 

	return H;

	/*
	CvMat* H, * A, * B, X;
	double x[9];
	int i;

	if(n < 4)
	{
		fprintf(stderr, "Warning: too few points in lsq_homog(), %s line %d\n",
			__FILE__, __LINE__);
		return vector<float>();
	}

	// set up matrices so we can unstack homography into X; AX = B
	A = cvCreateMat(2*n, 8, CV_64FC1);
	B = cvCreateMat(2*n, 1, CV_64FC1);
	X = cvMat(8, 1, CV_64FC1, x);
	H = cvCreateMat(3, 3, CV_64FC1);
	cvZero(A);
	for(i = 0; i < n; i++)
	{
		cvmSet(A, i, 0, pts[i].x);
		cvmSet(A, i+n, 3, pts[i].x);
		cvmSet(A, i, 1, pts[i].y);
		cvmSet(A, i+n, 4, pts[i].y);
		cvmSet(A, i, 2, 1.0);
		cvmSet(A, i+n, 5, 1.0);
		cvmSet(A, i, 6, -pts[i].x * mpts[i].x);
		cvmSet(A, i, 7, -pts[i].y * mpts[i].x);
		cvmSet(A, i+n, 6, -pts[i].x * mpts[i].y);
		cvmSet(A, i+n, 7, -pts[i].y * mpts[i].y);
		cvmSet(B, i, 0, mpts[i].x);
		cvmSet(B, i+n, 0, mpts[i].y);
	}
	cvSolve(A, B, &X, CV_SVD);
	x[8] = 1.0;
	X = cvMat(3, 3, CV_64FC1, x);
	cvConvert(&X, H);

	
	vector<float> Hmat(9, 0.f);
	for(size_t j=0, idx=0; j < 3; j++) {
		for(size_t i=0; i < 3; i++, idx++) {
			Hmat[idx] = (H->data.fl[j*3+i]);
			//cout << x[j*3+i] << ", ";
		} //cout << endl;
	} 
	//cout << endl;

	cvReleaseMat(&A);
	cvReleaseMat(&B);
	cvReleaseMat(&H);
	return Hmat;*/
}
float homog_xfer_err(fpoint pt, fpoint mpt, vector<float> &H)
{
	fpoint xpt = persp_xform_pt(pt, H);
	return sqrt(dist_sq_2D(xpt, mpt));
}
fpoint persp_xform_pt(fpoint pt, vector<float> &T)
{
	CvMat* _T;
	_T = cvCreateMat(3, 3, CV_32FC1);
	for (int i = 0; i < 9; i++)
	{
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
// 獲得這一點的相對匹配點(也可以檢測是否有效).
Feature* get_match(Feature* feat)
{
	return feat->fwd_match;
}
// 
int get_matched_features(Feature *features, int n, Feature ***matched)
{
	Feature **_matched;
	ransac_data* rdata;
	int i, m = 0;

	_matched = new Feature*[n];
	if (!_matched){
		cout << "memory error" << endl;
	} else {
		for (i = 0; i < n; i++)
		{
			if (get_match(features + i))
			{
				//cout << (features[i].feature_data==nullptr) << endl;
				//system("pause");
				rdata = new ransac_data{};
				rdata->orig_feat_data = features[i].feature_data;
				_matched[m] = features + i;
				_matched[m]->feature_data = rdata;
				m++;
			}
		}
	}
	cout << "m=" << m << endl;
	//system("pause");

	*matched = _matched;
	return m;
}
int calc_min_inliers(int n, int m, float p_badsupp, float p_badxform)
{
	float pi = 0.f, sum = 0.f;
	int i = 0, j = 0;
	float log_n = log_factorial(n);
	float log_nm = log_n - _log_factorial((n - m), n);
	vector<float> temp(n + 1, 0.f);
	float total_s = 0.f;
	for (i = m + 1; i <= n; i++)
	{
		pi = (float)(i - m) * log(p_badsupp) + (float)(n - i + m) * log(1.f - p_badsupp) + log_nm - log_factorial(i - m) - (log_n - _log_factorial(n - i, n));
		temp[i] = exp(pi);
		total_s += temp[i];
	}
	sum = total_s;
	for (j = m + 1; j <= n; j++)
	{
		sum -= temp[j];
		if (sum < p_badxform)
			break;
	}
	temp.clear();
	return j;
}
float _log_factorial(int n, int m)
{
	float f = 0.f;
	int i;

	for (i = n; i <= m; i++)
		f += log((float)i);

	return f;
}
float log_factorial(int n)
{
	float f = 0.f;
	int i;

	for (i = 1; i <= n; i++)
		f += log((float)i);

	return f;
}
Feature** draw_ransac_sample(Feature** features, int n, int m)
{
	Feature** sample, *feat;
	struct ransac_data* rdata;
	int i, x;

	for (i = 0; i < n; i++)
	{
		rdata = feat_ransac_data(features[i]);
		rdata->sampled = 0;
	}

	sample = new Feature*[m];
	for (i = 0; i < m; i++)
	{
		do
		{
			x = rand() % n;
			feat = features[x];
			rdata = feat_ransac_data(feat);
		} while (rdata->sampled);
		sample[i] = feat;
		rdata->sampled = 1;
	}

	return sample;
}

//从一致集中获取特征点和其对应匹配点的二维坐标，分别放到输出参数pts和mpts中  
void extract_corresp_pts(Feature** features, int n, fpoint** pts, fpoint** mpts)
{
	// 提取相應的資訊出來.
	Feature *match = nullptr;
	fpoint *_pts, *_mpts;

	_pts = new fpoint[n];
	_mpts = new fpoint[n];


	int ii = 0;
	for (int i = 0; i < n; i++)
	{
		// 獲取這一點，對應的匹配點.
		match = get_match(features[i]);
		if (!match)
			cout << "feature does not have match" << endl;
		else if(match) {
			// 原本寫法.
			//_pts[i] = features[i]->img_pt;
			//_mpts[i] = match->img_pt;
			// 無偏移的方式.
			_pts[i] = fpoint(features[i]->rX(), features[i]->rY());
			_mpts[i] = fpoint(match->rX(), match->rY());
			ii++;
		}
	}
	cout << "ii = " << ii;
	cout << endl;

	*pts = _pts;
	*mpts = _mpts;
}
int find_consensus(Feature** features, int n, vector<float> &M, float err_tol, Feature*** consensus)
{
	// 這裡帶入的 features 已經過濾過了一定有匹配點
	Feature** _consensus;
	Feature* match;
	fpoint pt, mpt;
	float err;
	int i, in = 0;

	int ii = 0;
	_consensus = new Feature*[n];

	for (i = 0; i < n; i++)
	{
		match = get_match(features[i]);
		if(!match) {
			cout << "feature does not have match" << endl;
		}
		else if (match) {
			// 原本寫法.
			//pt = features[i]->img_pt;
			//mpt = match->img_pt;
			// 無偏移的方式.
			pt = fpoint(features[i]->rX(), features[i]->rY());
			mpt = fpoint(match->rX(), match->rY());

			err = homog_xfer_err(pt, mpt, M); //與變形後的點算距離
			//cout << "err=" << err << ", err_tol=" << err_tol << endl;


			if(err <= err_tol) {
				//cout << "err <= err_tol" << endl;
				++ii;
				_consensus[in++] = features[i];
			}
		}
	}

	cout << "最後加入有幾點" << ii << endl;
	*consensus = _consensus;
	return in;
}
void release_mem(fpoint *pts1, fpoint *pts2, Feature** features)
{
	delete[] pts1;
	delete[] pts2;
	if (features)
		delete[] features;
}