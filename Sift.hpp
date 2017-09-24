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

#if defined(_MSC_VER)
    #define or ||
    #define and &&
    #define OR ||
    #define AND &&
#endif

#include "imglib\imglib.hpp"
#include "Raw2Img\Raw2Img.hpp"

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
struct Feature
{
	int x, y;//各所在階層的座標
	float mm;//強度
	int sita;//包含主方向與負方向的角度
	float size;//階
	int kai;//層
	float sigmaOCT;
	vector <vector <vector <float>>> descrip;
	Feature *nextptr = nullptr;
};
typedef struct Feature* Featureptr;

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
class ImgRaw {
private:
    using types = float;
public:
	// 初始化
    ImgRaw(vector<types> img, size_t width, size_t height) :
        raw_img(img), width(width), height(height) {}
	ImgRaw(size_t width, size_t height) :
		raw_img(width*height), width(width), height(height){}
	ImgRaw(size_t size = 0) :raw_img(size){}
	ImgRaw(string bmpname, bool gray_tran);
    // 隱式轉換
    operator vector<types>&() { return raw_img; }
	operator const vector<types>&() const { return raw_img; }
    operator vector<unsigned char>() {
        vector<unsigned char> img(raw_img.size());
        for(unsigned i = 0; i < raw_img.size(); ++i) {
            img[i] = (unsigned char)(raw_img[i]*255);
        }
        return img;
    }
    // 重載下標符號
    types & operator[](size_t idx) {
        return const_cast<types&>(static_cast<const ImgRaw&>(*this)[idx]);
    }
    const types & operator[](size_t idx) const {
        return raw_img[idx];
    }
    // 二維讀取
    types & at2d(size_t y, size_t x) {
        return const_cast<types&>(static_cast<const ImgRaw&>(*this).at2d(y, x));
    }
    const types & at2d(size_t y, size_t x) const {
        return raw_img[y*width + x];
    }
	// 重設大小
	void resize(size_t width, size_t height) {
		raw_img.resize(width*height);
		this->width=width;
		this->height=height;
	}
public:
    // 大小是否相等
    friend bool operator!=(const ImgRaw& lhs, const ImgRaw& rhs);
    friend bool operator==(const ImgRaw& lhs, const ImgRaw& rhs);
    // 差分圖
    ImgRaw& operator-=(const ImgRaw& rhs) {
        // 尺寸不同
        if ((*this) != rhs) {
            throw Size_error("Error size is diff.");
        }
        for (unsigned i = 0; i < width*height; ++i) {
            raw_img[i] -= rhs.raw_img[i];
        }
        return (*this);
    }
    friend ImgRaw operator-(ImgRaw lhs, const ImgRaw& rhs) {
        return lhs -= rhs;
    }
    // 寫 BMP 檔
    ImgRaw& bmp(string name, size_t bits=0) {
		if (bits == 0) {
			bits = this->bitCount;
		}
		// 有重載轉換函式(感覺不是很好，有空修掉)
        vector<unsigned char> img = (*this);
        Raw::raw2bmp(name, img, width, height, bits);
        return (*this);
    }

public:
    static void zero(ImgRaw& tar, ImgRaw& sou, float z) {
        Scaling::zero(tar, sou, sou.width, sou.height, z);
        tar.width = (size_t)(sou.width*z);
        tar.height = (size_t)(sou.height*z);
    }
    static void first(ImgRaw& tar, ImgRaw& sou, float z) {
        Scaling::first(tar, sou, sou.width, sou.height, z);
        tar.width = (size_t)(sou.width*z);
        tar.height = (size_t)(sou.height*z);
    }
    static void cubic(ImgRaw& tar, ImgRaw& sou, float z) {
        Scaling::cubic(tar, sou, sou.width, sou.height, z);
        tar.width = (size_t)(sou.width*z);
        tar.height = (size_t)(sou.height*z);
    }
    static void gauBlur(ImgRaw& tar, ImgRaw& sou, float p) {
        Gaus::GauBlur(tar, sou, sou.width, sou.height, p);
        tar.width = (size_t)(sou.width);
        tar.height = (size_t)(sou.height);
    }
public:
    vector<types> raw_img;
    uint32_t width;
	uint32_t height;
	uint16_t bitCount = 0;
	float sigma = 0;
};
// 大小是否相等
inline bool operator!=(const ImgRaw& lhs, const ImgRaw& rhs) {
    return !(lhs == rhs);
}
inline bool operator==(const ImgRaw& lhs, const ImgRaw& rhs) {
    if (lhs.width == rhs.width &&  lhs.height == rhs.height) {
        return 1;
    } return 0;
}
//----------------------------------------------------------------
class Sift {
private:
    using types = float;
public:
    Sift(ImgRaw img);
    Sift(vector<types> raw_img, size_t width, size_t height):
        raw_img(raw_img, width, height) {}
public:
    void pyramid(size_t s = 3); // 3 為論文中所給的
	void pyramid2(); // 3 為論文中所給的
    void comp(vector<ImgRaw>& pyrs, string name="");
    vector<ImgRaw> dog_gau(ImgRaw& img, size_t s, size_t o=1);
	float* getFea(ImgRaw& img, size_t y, size_t x, float sigma, size_t r);
	void display();
private:	
	float fea_m(ImgRaw& img, size_t y, size_t x);
	float fea_sida(ImgRaw& img, size_t y, size_t x);
	bool findMaxMin(vector<ImgRaw>& gauDog_imgs, size_t scale_idx, size_t curr_Width, size_t y, size_t x);
	void ZoomInOut(ImgRaw& doImage, int InWidth, int InHeight);
	void getHistogramMS(const ImgRaw& doImage, size_t Iny, size_t Inx,size_t Inr,float Insize, size_t InWidth, float sigma, int scale);
	void AddnewFeaturestruct(int Inx, int Iny, float Insize, int kai, int sigmaOCT, float Inm, int Insita);
private:
    ImgRaw raw_img;  //原圖
	size_t pyWidth=5;  //塔高(放大縮小)
	size_t pyheight=5; //塔寬(模糊幾次)
    vector<vector<ImgRaw>> pyrs;
	vector<vector<ImgRaw>> pyrs_dog;
	Featureptr FeatureStart, FeatureEnd, FeatureNow; // 特徵點
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
// 畫線
class Draw {
public:
	static void draw_line(ImgRaw& img, size_t y, size_t x, float line_len, float sg);
};