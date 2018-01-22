/*****************************************************************
Name :
Date : 2017/07/05
By   : CharlotteHonG
Final: 2017/07/05


                       _oo0oo_
                      o8888888o
                      88" . "88
                      (| -_- |)
                      0\  =  /0
                    ___/`---'\___
                  .' \\|     |// '.
                 / \\|||  :  |||// \
                / _||||| -:- |||||- \
               |   | \\\  -  /// |   |
               | \_|  ''\---/''  |_/ |
               \  .-\__  '-'  ___/-. /
             ___'. .'  /--.--\  `. .'___
          ."" '<  `.___\_<|>_/___.' >' "".
         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
         \  \ `_.   \_ __\ /__ _/   .-` /  /
     =====`-.____`.___ \_____/___.-`___.-'=====
                       `=---='


     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

               �򯪫O��         �õLBUG

*****************************************************************/
#pragma warning(disable : 4819)
#pragma once

#include "imglib\imglib.hpp"
#include "Raw2Img\Raw2Img.hpp"
#include "Imgraw.hpp"
#include "imagedata.hpp"

// �פ夤��s
#define SIFT_Sacle 3
// �䷥�ȱ˥h�e���i + �t���Ϥ�1�i
#define SIFT_SacleDiff 3
// �����ҽk����l�Y��
#define SIFT_GauSigma 1.6f
// �h����í�w�S�x�I
#define SIFT_Dx 0.03f
// ��פ嵹�� 3*1.5*sigma
#define SIFT_DESCR_SCL_FCTR 3.f
// �S�x�y�z�l8�Ө���
#define SIFT_DESCR_HIST_BINS 8
// �S�x�I�V�q�������֭�
#define SIFT_DESCR_MAG_THR 0.2
// �N�B�I���ഫ�� uchar
#define SIFT_INT_DESCR_FCTR 512.0
// �D���v�S�x�I��һ֭�(���I�����ɪ��֭�)
#define SIFT_CURV_THR 10

/* maximum steps of keypoint interpolation before failure */
#define SIFT_MAX_INTERP_STEPS 5
/* ��������t�e�� */
#define SIFT_IMG_BORDER 5

/** �S�x�I����q�{�֭� |D(x)| */
#define SIFT_CONTR_THR 0.05


// �ǰt�ۦP���S�x�I
/* the maximum number of keypoint NN candidates to check during BBF search */
#define KDTREE_BBF_MAX_NN_CHKS 200
/* threshold on squared ratio of distances between NN and 2nd NN */
#define NN_SQ_DIST_RATIO_THR 0.4


// �S�x�I���c

//----------------------------------------------------------------
// �ؤo�j�p���X
class Size_error : public std::runtime_error {
public:
    Size_error(const std::string& str): std::runtime_error(str) {}
};
// �W�X�Ϥ����
class Out_of_imgRange : public std::runtime_error {
public:
	Out_of_imgRange(const std::string& str): std::runtime_error(str) {}
};
//----------------------------------------------------------------
class Sift {
private:
    using types = float;
	using Desc = vector<vector<vector<float>>>;
public:
    Sift(ImgRaw img, size_t intvls=6);

public:
	void pyramid();
    void comp(vector<ImgRaw>& pyrs, string name="");
	void drawArrow(string name="feaArrow.bmp", float ratio = 10000.f);
private:
	bool findMaxMin(vector<ImgRaw>& gauDog_imgs, size_t scale_idx, size_t curr_Width, size_t y, size_t x);
	// �y�z�S�x�l
	static bool calc_grad_mag_ori(const vector<float> &img, int &COL, int &ROW, int r, int c, float &mag, float &ori);
	static Desc descr_hist(vector<float> &img, int &COL, int &ROW, int r, int c, float ori, float scl, int d, int n);
	static void interp_hist_entry(Desc &hist, float rbin, float cbin, float obin, float mag, int d, int n);
	static void hist_to_descr(const Desc &hist, int d, int n, Feature* feat);
	static void normalize_descr(Feature* feat);
	static void FeatureDescrip(vector<ImgRaw>& kaidaImag, Feature* FeatureNow);

	void getHistogramMS(Feature* NweFeat, const ImgRaw& doImage, float Insize, size_t scale, float sigma, 
		size_t Iny, size_t Inx, size_t InWidth, size_t Inr);
	void FeatAppend(Feature* NweFeat, int Inx, int Iny, float Insize, int kai, int sigmaOCT, float Inm, int Insita);
	void getHistogramMS(const ImgRaw& doImage, float Insize, size_t scale, float sigma, 
		size_t Iny, size_t Inx, size_t InWidth, size_t Inr);
	void FeatAppend(int Inx, int Iny, float Insize, int kai, int sigmaOCT, float Inm, int Insita);
private:
	bool harris(const vector<float>& p, size_t w, size_t y, size_t x, float r=SIFT_CURV_THR);
	static void DescripNomal(Desc& descripgroup);
public:
    ImgRaw raw_img;    // ���
	size_t pyWidth=6;  // ��(��j�Y�p)
    vector<vector<ImgRaw>> pyrs;
	vector<vector<ImgRaw>> pyrs_dog;
	Feature* FeatStart; // �S�x�I
	Feature* FeatEnd;   // �аO�~���αq�Y�A��
	size_t feaNum=0;    // �S�x�I�ƶq
};
//-----------------------------------------------------------------
class Stitching{
private:
	using Desc = std::vector<std::vector<std::vector<float>>>;
public:
	Stitching(const Sift& desc1, const Sift& desc2);
	static float EuclDist(const Desc& point1, const Desc& point2); // �y�z�l�ڦ��Z��
	float Stitching::EuclDist2(float point1[128], float point2[128]); // �y�z�l�ڦ��Z��
	void Check(float matchTh = NN_SQ_DIST_RATIO_THR); // �ˬd�O�_���ۦP���S�x�y�z�l
private:
	Feature *feat1, *feat2;
	size_t Width, Height;
	ImgRaw matchImg;
	size_t feat1_Count, feat2_Count;
};
