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
	Feature **matched = nullptr, **sample = nullptr;
	Feature **consensus = nullptr, **consensus_max = nullptr;
	ransac_data* rdata = nullptr;
	fpoint *pts = nullptr, *mpts = nullptr;
	vector<float> M;
	float p, in_frac = (float)RANSAC_INLIER_FRAC_EST;
	int i, nm, in_min, k = 0, in = 0, in_max = 0;

	nm = get_matched_features(features, n, &matched);
	if (nm < m)
	{
		cout << "Warning: not enough matches to compute xform" << endl;
		goto end;
	}

	srand((unsigned int)time(NULL));

	in_min = calc_min_inliers(nm, m, (float)RANSAC_PROB_BAD_SUPP, p_badxform);

	p = pow(1.f - pow(in_frac, m), k);

	while (p > p_badxform)
	{
		sample = draw_ransac_sample(matched, nm, m);
		extract_corresp_pts(sample, m, &pts, &mpts);
		//cout << pts[1].x << ", " << pts[1].y << ", " << mpts[1].x << ", " << mpts[1].y <<endl;
		M = lsq_homog(pts, mpts, m);
		if (M.empty())
			goto iteration_end;
		in = find_consensus(matched, nm, M, err_tol, &consensus);
		if (in > in_max)
		{
			if (consensus_max)
				delete[] consensus_max;
			consensus_max = consensus;
			in_max = in;
			in_frac = (float)in_max / (float)nm;
		}
		else
			delete[] consensus;
		M.clear();

	iteration_end:
		release_mem(pts, mpts, sample);
		p = pow(1.f - pow(in_frac, m), ++k);
	}

	if (in_max >= in_min)
	{
		extract_corresp_pts(consensus_max, in_max, &pts, &mpts);
		M = lsq_homog(pts, mpts, in_max);
		in = find_consensus(matched, nm, M, err_tol, &consensus);
		M.clear();
		release_mem(pts, mpts, consensus_max);
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
		M = lsq_homog(pts, mpts, in);
		if (inliers)
		{
			*inliers = consensus;
			consensus = NULL;
		}
		if (n_in == 0)
			n_in = in;
		release_mem(pts, mpts, consensus);
	}
	else if (consensus_max)
	{
		if (inliers)
			*inliers = NULL;
		if (n_in == 0)
			n_in = 0;
		delete[] consensus_max;
	}
end:
	for (i = 0; i < nm; i++)
	{
		rdata = feat_ransac_data(matched[i]);
		matched[i]->feature_data = rdata->orig_feat_data;
		delete rdata;
	}
	delete[] matched;
	return M;
}
/*
// 沒用到
vector<float> _ransac_xform(struct feature *features, int n, int m, float p_badxform, float err_tol)
{
	struct feature **matched, **sample, **consensus, **consensus_max = NULL;
	struct ransac_data* rdata;
	fpoint *pts, *mpts;
	vector<float> M;
	float p, in_frac = (float)RANSAC_INLIER_FRAC_EST;
	int i, nm, in_min, k = 0, in = 0, in_max = 0;

	nm = get_matched_features(features, n, &matched);
	if (nm < m)
	{
		cout << "Warning: not enough matches to compute xform" << endl;
		goto end;
	}

	srand((unsigned int)time(NULL));

	in_min = calc_min_inliers(nm, m, (float)RANSAC_PROB_BAD_SUPP, p_badxform);

	p = pow(1.f - pow(in_frac, m), k);

	while (p > p_badxform)
	{
		sample = draw_ransac_sample(matched, nm, m);
		extract_corresp_pts(sample, m, &pts, &mpts);
		M = lsq_homog(pts, mpts, m);
		if (M.empty())
			goto iteration_end;
		in = find_consensus(matched, nm, M, err_tol, &consensus);
		if (in > in_max)
		{
			if (consensus_max)
				delete[] consensus_max;
			consensus_max = consensus;
			in_max = in;
			in_frac = (float)in_max / (float)nm;
		}
		else
			delete[] consensus;
		M.clear();

	iteration_end:
		release_mem(pts, mpts, sample);
		p = pow(1.f - pow(in_frac, m), ++k);
	}

	if (in_max >= in_min)
	{
		extract_corresp_pts(consensus_max, in_max, &pts, &mpts);
		M = lsq_homog(pts, mpts, in_max);
		in = find_consensus(matched, nm, M, err_tol, &consensus);
		M.clear();
		release_mem(pts, mpts, consensus_max);
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
		M = lsq_homog(pts, mpts, in);
		release_mem(pts, mpts, consensus);
	}
	else if (consensus_max)
	{
		delete[] consensus_max;
	}
end:
	for (i = 0; i < nm; i++)
	{
		rdata = feat_ransac_data(matched[i]);
		matched[i]->feature_data = rdata->orig_feat_data;
		delete rdata;
	}
	delete[] matched;
	return M;
}
*/

/************************ Local funciton definitions *************************/
// 顏自己寫的，輸入兩點匯出變換矩陣 n 是最少點數輸入是 4
vector<float> lsq_homog(fpoint * pts, fpoint * mpts, int n)
{
	/*
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

	return H;
	*/

	CvMat* H, * A, * B, X;
	double x[9];
	int i;

	if(n < 4)
	{
		fprintf(stderr, "Warning: too few points in lsq_homog(), %s line %d\n",
			__FILE__, __LINE__);
		return vector<float>();
	}

	/* set up matrices so we can unstack homography into X; AX = B */
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
	return Hmat;
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
Feature* get_match(Feature* feat)
{
	return feat->fwd_match;
}
int get_matched_features(Feature *features, int n, Feature ***matched)
{
	Feature **_matched;
	struct ransac_data* rdata;
	int i, m = 0;

	_matched = new Feature*[n];
	if (!_matched){
		cout << "memory error" << endl;
	} else {
		for (i = 0; i < n; i++)
		{
			if (get_match(features + i))
			{
				rdata = new ransac_data;
				memset(rdata, 0, sizeof(ransac_data));
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
void extract_corresp_pts(Feature** features, int n, fpoint** pts, fpoint** mpts)
{
	// 提取相應的資訊出來
	Feature *match = nullptr;
	fpoint *_pts, *_mpts;

	_pts = new fpoint[n];
	_mpts = new fpoint[n];

	for (int i = 0; i < n; i++)
	{
		// 獲取這一點，對應的匹配點
		match = get_match(features[i]);
		if (!match)
			cout << "feature does not have match" << endl;
		// 我的寫法沒有寫這個東西
		//_pts[i] = features[i]->img_pt;
		//_mpts[i] = match->img_pt;
		// 補救方式 
		_pts[i] = fpoint(features[i]->rX(), features[i]->rY());
		_mpts[i] = fpoint(match->rX(), match->rY());

	}

	*pts = _pts;
	*mpts = _mpts;
}
int find_consensus(Feature** features, int n, vector<float> &M, float err_tol, Feature*** consensus)
{
	Feature** _consensus;
	Feature* match;
	fpoint pt, mpt;
	float err;
	int i, in = 0;

	_consensus = new Feature*[n];

	for (i = 0; i < n; i++)
	{
		match = get_match(features[i]);
		if (!match)
			cout << "feature does not have match" << endl;
		// 我的寫法沒有寫這個東西
		//pt = features[i]->img_pt;
		//mpt = match->img_pt;
		// 補救方式 
		pt = fpoint(features[i]->rX(), features[i]->rY());
		mpt = fpoint(match->rX(), match->rY());

		err = homog_xfer_err(pt, mpt, M);
		if (err <= err_tol)
			_consensus[in++] = features[i];
	}

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