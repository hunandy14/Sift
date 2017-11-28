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

               佛祖保佑         永無BUG

*****************************************************************/
#pragma warning(disable : 4819)
#pragma once

#include "imglib\imglib.hpp"
#include "Raw2Img\Raw2Img.hpp"
#include "Imgraw.hpp"

// 論文中的s
#define SIFT_Sacle 3
// 找極值捨去前後兩張 + 差分圖少1張
#define SIFT_SacleDiff 3
// 高斯模糊的初始係數
#define SIFT_GauSigma 1.6f
// 去除不穩定特徵點
#define SIFT_Dx 0.03f
// 角點偵測 r = 10
#define SIFT_HarrisR 10

// 特徵點結構
struct Feature {
	float size;//階
	int kai;//層
	float sigmaOCT;//高斯模糊係數
	int x, y;//各所在階層的座標
	float mm;//強度
	int sita;//包含主方向與負方向的角度
	vector<vector<vector<float>>> descrip;// 描述子(舊方法的)
	Feature* nextptr = nullptr; // 下一點

	float descr[128] = {};// 統計完成後的描述子(robbs的方法)
	int d; // 特徵點長度
};
// 尺寸大小不合
class Size_error : public std::runtime_error {
public:
    Size_error(const std::string& str): std::runtime_error(str) {}
};
// 超出圖片邊界
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
	void drawArrow(string name="feaArrow.bmp");
private:
	bool findMaxMin(vector<ImgRaw>& gauDog_imgs, size_t scale_idx, size_t curr_Width, size_t y, size_t x);
	void FeatureDescrip(vector<ImgRaw>& kaidaImag, Feature* FeatureNow);
	void FeatureDescrip2(vector<ImgRaw>& kaidaImag, Feature* FeatureNow);

	
	static bool calc_grad_mag_ori(const vector<float> &img, int &COL, int &ROW, int r, int c, float &mag, float &ori);
	static Desc descr_hist(vector<float> &img, int &COL, int &ROW, int r, int c, float ori, float scl, int d, int n);
	static void interp_hist_entry(Desc &hist, float rbin, float cbin, float obin, float mag, int d, int n);
	static void hist_to_descr(Desc &hist, int d, int n, Feature* feat);
	static void normalize_descr(Feature* feat);
	static void FeatureDescrip3(vector<ImgRaw>& kaidaImag, Feature* FeatureNow);


	void FeatureDescrip_ori(vector<ImgRaw>& kaidaImag, Feature* FeatureNow);
	
	void getHistogramMS(const ImgRaw& doImage, float Insize, size_t scale, float sigma, 
		size_t Iny, size_t Inx, size_t InWidth, size_t Inr);
	void AddnewFeaturestruct(int Inx, int Iny, float Insize, int kai, int sigmaOCT, float Inm, int Insita);
private:
	static void DescripNomal(Desc& descripgroup);
public:
    ImgRaw raw_img;  //原圖
	size_t pyWidth=6;  //塔高(放大縮小)
    vector<vector<ImgRaw>> pyrs;
	vector<vector<ImgRaw>> pyrs_dog;
	// 特徵點
	Feature* FeatStart;
	Feature* FeatEnd;
};
//----------------------------------------------------------------
struct Fea_point {
	Fea_point(size_t o, size_t s, size_t y, size_t x, 
		float gau_r ,float sigma, float m=0.f, float sida=0.f):
		o(o), s(s), y(y), x(x), gau_r(gau_r), sigma(sigma), m(m), sida(sida){}
	size_t o;
	size_t s;
	size_t y;
	size_t x;
	size_t gau_r;
	float sigma;
	float m;
	float sida;
};

//-----------------------------------------------------------------
class Stitching{
private:
	using Desc = std::vector<std::vector<std::vector<float>>>;
public:
	Stitching(const Sift& desc1, const Sift& desc2);
	static float EuclDist(const Desc& point1, const Desc& point2); // 描述子歐式距離
	float Stitching::EuclDist2(float point1[128], float point2[128]); // 描述子歐式距離
	void Check(float matchTh=0.6); // 檢查是否有相同的特徵描述子
private:
	int Width, Height;
	Feature *FeatureStart1, *FeatureStart2;
	ImgRaw matchImg;
};
