#pragma once

#include "imagedata.hpp"


using std::vector;
//using pyData = vector<vector<float>>;
using pyData = vector<ImgRaw>;


//计算三维偏导数
vector<float> deriv_3D(pyData &dog_pyr, 
	int Col, int intvl, int r, int c)
{
    float dx = (
		dog_pyr[(intvl)][(r)* Col + (c + 1)] -
        dog_pyr[(intvl)][(r)* Col + (c - 1)]
	)/ 2.f;
    float dy = (
		dog_pyr[(intvl)][(r + 1) * Col + (c)] -
        dog_pyr[(intvl)][(r - 1) * Col + (c)]
	) / 2.f;
    float ds = (
		dog_pyr[(intvl + 1)][(r)* Col + (c)] -
        dog_pyr[(intvl - 1)][(r)* Col + (c)]
	) / 2.f;

    vector<float> dI(3, 0.f);
    dI[0] = dx;
    dI[1] = dy;
    dI[2] = ds;
    return dI;
}
vector<float> hessian_3D(pyData &dog_pyr, 
    int Col, int intvl, int r, int c)
{
    float v   =  dog_pyr[intvl][r * Col + c];
    float dxx = (dog_pyr[(intvl)][(r)* Col + (c + 1)] +
        dog_pyr[(intvl)][(r)* Col + (c - 1)] - 2.f * v);
    float dyy = (dog_pyr[(intvl)][(r + 1) * Col + (c)] +
        dog_pyr[(intvl)][(r - 1) * Col + (c)] - 2.f * v);
    float dss = (dog_pyr[(intvl + 1)][(r)* Col + (c)] +
        dog_pyr[(intvl - 1)][(r)* Col + (c)] - 2.f * v);
    float dxy = (dog_pyr[(intvl)][(r + 1) * Col + (c + 1)] -
        dog_pyr[(intvl)][(r + 1) * Col + (c - 1)] -
        dog_pyr[(intvl)][(r - 1) * Col + (c + 1)] +
        dog_pyr[(intvl)][(r - 1) * Col + (c - 1)]) / 4.f;
    float dxs = (dog_pyr[(intvl + 1)][(r)* Col + (c + 1)] -
        dog_pyr[(intvl + 1)][(r)* Col + (c - 1)] -
        dog_pyr[(intvl - 1)][(r)* Col + (c + 1)] +
        dog_pyr[(intvl - 1)][(r)* Col + (c - 1)]) / 4.f;
    float dys = (dog_pyr[(intvl + 1)][(r + 1) * Col + (c)] -
        dog_pyr[(intvl + 1)][(r - 1) * Col + (c)] -
        dog_pyr[(intvl - 1)][(r + 1) * Col + (c)] +
        dog_pyr[(intvl - 1)][(r - 1) * Col + (c)]) / 4.f;
    vector<float> H(3 * 3, 0.f);
    H[0] = dxx;
    H[1] = dxy;
    H[2] = dxs;
    H[3] = dxy;
    H[4] = dyy;
    H[5] = dys;
    H[6] = dxs;
    H[7] = dys;
    H[8] = dss;
    return H;
}
vector<float> Invert3X3(vector<float> &v)
{
    vector<float> m(9, 0.f);
    float det_v  = v[0] * (v[4] * v[8] - v[5] * v[7]);
    det_v -= v[3] * (v[1] * v[8] - v[2] * v[7]);
    det_v += v[6] * (v[1] * v[5] - v[2] * v[4]);

    m[0] =  1 * (v[4] * v[8] - v[5] * v[7]) / det_v;
    m[1] = -1 * (v[1] * v[8] - v[2] * v[7]) / det_v;
    m[2] =  1 * (v[1] * v[5] - v[2] * v[4]) / det_v;
    m[3] = -1 * (v[3] * v[8] - v[5] * v[6]) / det_v;
    m[4] =  1 * (v[0] * v[8] - v[2] * v[6]) / det_v;
    m[5] = -1 * (v[0] * v[5] - v[2] * v[3]) / det_v;
    m[6] =  1 * (v[3] * v[7] - v[4] * v[6]) / det_v;
    m[7] = -1 * (v[0] * v[7] - v[1] * v[6]) / det_v;
    m[8] =  1 * (v[0] * v[4] - v[1] * v[3]) / det_v;

    return m;
}
vector<float> GEMM_3X1(vector<float> &src1, vector<float> &src2, float alpha)
{
    vector<float> data(3, 0.f);

    data[0] = (src1[0] * src2[0] + src1[1] * src2[1] + src1[2] * src2[2]) * alpha;
    data[1] = (src1[3] * src2[0] + src1[4] * src2[1] + src1[5] * src2[2]) * alpha;
    data[2] = (src1[6] * src2[0] + src1[7] * src2[1] + src1[8] * src2[2]) * alpha;

    return data;
}
float GEMM_1X1(vector<float> &src1, vector<float> &src2)
{
    float data = src1[0] * src2[0] + src1[1] * src2[1] + src1[2] * src2[2];
    return data;
}
void interp_step(pyData &dog_pyr, int Col, int intvl, int r, int c, float &xi, float &xr, float &xc)
{
	// intvl=第幾張模糊的圖, col=圖寬度
	vector<float> dD = deriv_3D(dog_pyr, Col, intvl, r, c);
	vector<float> H = hessian_3D(dog_pyr, Col, intvl, r, c);
	vector<float> H_inv = Invert3X3(H);
	vector<float> x = GEMM_3X1(H_inv, dD, -1.0);
	xi = x[2];
	xr = x[1];
	xc = x[0];
}
float interp_contr(pyData &dog_pyr, int Col, int intvl, int r, int c, float xi, float xr, float xc)
{
	float t = 0.0;
	vector<float> x(3, 0.0);
	x[0] = xc;
	x[1] = xr;
	x[2] = xi;
	vector<float> dD = deriv_3D(dog_pyr, Col, intvl, r, c);
	t = GEMM_1X1(dD, x);

	return dog_pyr[intvl][r * Col + c] + t * 0.5f;
}


/*特徵點定位 DoG函數*/
static Feature* new_feature(void)
{
	Feature* feat;
	detection_data* ddata;
	feat = new Feature{};
	ddata = new detection_data{};
	feat->feature_data = ddata;
	return feat;
}
/* 回傳特徵點的相關檢測數據 */
#define feat_detection_data(f) ((detection_data*)(f->feature_data))

Feature* interp_extremum(pyData &dog_pyr, 
	int COL, int ROW, 
	int octv, int intvl, 
	int r, int c, 
	int intvl_Len, double contr_thr)
{
	/*
	COL		= 所在圖寬度
	ROW		= 所在圖像高
	octv	= 金字塔y 第幾個大小 的idx
	intvl	= 金字塔x 第幾個模糊 的idx
	*/
	Feature* feat = nullptr;
	detection_data* ddata = nullptr;
	float xi = 0.0, xr = 0.0, xc = 0.0, contr = 0.0;

	int i = 0;
	while (i < SIFT_MAX_INTERP_STEPS)
	{
		interp_step(dog_pyr, COL, intvl, r, c, xi, xr, xc);
		if (fabs(xi) < 0.5  &&  fabs(xr) < 0.5  &&  fabs(xc) < 0.5)
			break;

		c	  += (int)round(xc);
		r	  += (int)round(xr);
		intvl += (int)round(xi);

		if (intvl < 1				|| intvl >  intvl_Len						||
			c	  < SIFT_IMG_BORDER || c	 >= COL - SIFT_IMG_BORDER ||
			r	  < SIFT_IMG_BORDER || r	 >= ROW - SIFT_IMG_BORDER )
		{
			return nullptr;
		}

		i++;
	}

	/* ensure convergence of interpolation */
	if (i >= SIFT_MAX_INTERP_STEPS)
		return nullptr;

	contr = interp_contr(dog_pyr, COL, intvl, r, c, xi, xr, xc);
	//if (fabs(contr) < (float)(contr_thr / intvls))
	if (fabs(contr) < (float)(contr_thr))
		return nullptr;

	feat = new_feature();
	ddata = feat_detection_data(feat);
	feat->img_pt.x = ((float)c + xc) * (float)pow(2.0, octv) *.5; // 有放大所以除2
	feat->img_pt.y = ((float)r + xr) * (float)pow(2.0, octv) *.5; // 有放大所以除2
	ddata->r = r;
	ddata->c = c;
	ddata->octv = octv;
	ddata->intvl = intvl;
	ddata->subintvl = xi;
	return feat;
}