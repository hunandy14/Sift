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

#if defined(_MSC_VER)
    #define or ||
    #define and &&
    #define OR ||
    #define AND &&
#endif

#include <memory>

#include "imglib\imglib.hpp"
#include "Raw2Img\Raw2Img.hpp"

// �פ夤��s
#define SIFT_Sacle 3
// �䷥�ȱ˥h�e���i + �t���Ϥ�1�i
#define SIFT_SacleDiff 3
// �����ҽk����l�Y��
#define SIFT_GauSigma 1.6f
// �h����í�w�S�x�I
#define SIFT_Dx 0.03f
// ���I���� r = 10
#define SIFT_HarrisR 10

// �S�x�I���c
struct Feature
{
	float size;//��
	int kai;//�h
	float sigmaOCT;//�����ҽk�Y��
	int x, y;//�U�Ҧb���h���y��
	float mm;//�j��
	int sita;//�]�t�D��V�P�t��V������
	vector<vector<vector<float>>> descrip;// �y�z�l
	Feature* nextptr = nullptr;
};
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
class ImgRaw {
private:
    using types = float;
public:
	// ��l��
	ImgRaw() = default;
    ImgRaw(vector<types> img, size_t width, size_t height) :
        raw_img(img), width(width), height(height) {}
	ImgRaw(size_t width, size_t height, size_t bits) :
		raw_img(width*height * (bits/8)), 
		width(width), height(height), bitCount(bits){}
	ImgRaw(size_t width, size_t height) :
		raw_img(width*height), width(width), height(height){}
	ImgRaw(string bmpname);
    // �����ഫ
    operator vector<types>&() { return raw_img; }
	operator const vector<types>&() const { return raw_img; }
    operator vector<unsigned char>() {
		const vector<unsigned char> img = static_cast<const ImgRaw&>(*this);
		return const_cast<vector<unsigned char>&>(img);
    }
	operator const vector<unsigned char>() const {
		vector<unsigned char> img(raw_img.size());
		for(unsigned i = 0; i < raw_img.size(); ++i)
			img[i] = (unsigned char)(raw_img[i]*255);
		return img;
	}
    // �����U�вŸ�
    types & operator[](size_t idx) {
        return const_cast<types&>(static_cast<const ImgRaw&>(*this)[idx]);
    }
    const types & operator[](size_t idx) const {
        return raw_img[idx];
    }
    // �G��Ū�� (ı�o���طQ����)
    /*types & at2d(size_t y, size_t x) {
        return const_cast<types&>(static_cast<const ImgRaw&>(*this).at2d(y, x));
    }
    const types & at2d(size_t y, size_t x) const {
        return raw_img[y*width + x];
    }*/
    // �j�p�O�_�۵�
    friend bool operator!=(const ImgRaw& lhs, const ImgRaw& rhs);
    friend bool operator==(const ImgRaw& lhs, const ImgRaw& rhs);
	// ��o�j�p
	size_t size() {
		return this->raw_img.size();
	}
	// ���]�j�p
	void resize(size_t width, size_t height, size_t bits) {
		raw_img.resize(width*height * bits/8);
		this->width=width;
		this->height=height;
		this->bitCount=bits;
	}
public:
	// �ର�Ƕ�
	ImgRaw ConverGray() const {
		if (bitCount == 24) {
			ImgRaw gray(this->width, this->height, 8);
			for (size_t i = 0; i < gray.size(); i++) {
				gray[i] =
					raw_img[i*3+0]*0.299 +
					raw_img[i*3+1]*0.587 +
					raw_img[i*3+2]*0.114;
			} return gray;
		} return (*this);
	}
    // �g BMP ��
    void bmp(string name, size_t bits=0) {
		if (bits == 0) { bits = this->bitCount; }
        vector<unsigned char> img = (*this);// �������ഫ�禡
        Raw::raw2bmp(name, img, width, height, bits);
    }
public: // ��j�Y�p (ı�o���طQ����)
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
};
// �j�p�O�_�۵�
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
    Sift(ImgRaw img, size_t intvls=6);
public:
	void pyramid2();
    void comp(vector<ImgRaw>& pyrs, string name="");
	void drawArrow(string name="feaArrow.bmp");
private:
	bool findMaxMin(vector<ImgRaw>& gauDog_imgs, size_t scale_idx, size_t curr_Width, size_t y, size_t x);
	void FeatureDescrip(vector<ImgRaw>& kaidaImag, Feature* FeatureNow);
	void getHistogramMS(const ImgRaw& doImage, float Insize, size_t scale, float sigma, 
		size_t Iny, size_t Inx, size_t InWidth, size_t Inr);
	void AddnewFeaturestruct(int Inx, int Iny, float Insize, int kai, int sigmaOCT, float Inm, int Insita);
public:
    ImgRaw raw_img;  //���
	size_t pyWidth=6;  //��(��j�Y�p)
    vector<vector<ImgRaw>> pyrs;
	vector<vector<ImgRaw>> pyrs_dog;
	// �S�x�I
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
// �e�u
class Draw {
public:
	static void drawLine_p(ImgRaw& img, int y, int x, int y2, int x2);
	static void drawLine_s(ImgRaw& img, int y, int x, float line_len, float sg);
	static void draw_arrow(ImgRaw& img, int y, int x, float line_len, float sg);

	static void drawLineRGB_p(ImgRaw& img, int y, int x, int y2, int x2);
	static void drawLineRGB_s(ImgRaw& img, int y, int x, float line_len, float sg);
	static void draw_arrowRGB(ImgRaw& img, int y, int x, float line_len, float sg);
};
//-----------------------------------------------------------------
class Stitching{
private:
	using Desc = std::vector<std::vector<std::vector<float>>>;
public:
	Stitching(const Sift& desc1, const Sift& desc2);
	static float EuclDist(const Desc& point1, const Desc& point2); // �y�z�l�ڦ��Z��
	void Check(float matchTh=0.6); // �ˬd�O�_���ۦP���S�x�y�z�l
	void Link(int x1, int y1, int x2, int y2);// �N�a�J�����I�۳s
private:
	int Width, Height;
	Feature *FeatureStart1, *FeatureStart2;
	ImgRaw matchImg;
};
