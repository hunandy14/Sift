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
    ImgRaw(vector<types> img, size_t width, size_t height) :
        raw_img(img), width(width), height(height) {}
	ImgRaw(size_t width, size_t height) :
		raw_img(width*height), width(width), height(height){}
	ImgRaw(size_t size = 0) :raw_img(size){}
	ImgRaw(string bmpname, bool gray_tran);
    // �����ഫ
    operator vector<types>&() { return raw_img; }
	operator const vector<types>&() const { return raw_img; }
    operator vector<unsigned char>() {
        vector<unsigned char> img(raw_img.size());
        for(unsigned i = 0; i < raw_img.size(); ++i) {
            img[i] = (unsigned char)(raw_img[i]*255);
        }
        return img;
    }
    // �����U�вŸ�
    types & operator[](size_t idx) {
        return const_cast<types&>(static_cast<const ImgRaw&>(*this)[idx]);
    }
    const types & operator[](size_t idx) const {
        return raw_img[idx];
    }
    // �G��Ū��
    types & at2d(size_t y, size_t x) {
        return const_cast<types&>(static_cast<const ImgRaw&>(*this).at2d(y, x));
    }
    const types & at2d(size_t y, size_t x) const {
        return raw_img[y*width + x];
    }
	// ���]�j�p
	void resize(size_t width, size_t height) {
		raw_img.resize(width*height);
		this->width=width;
		this->height=height;
	}
public:
    // �j�p�O�_�۵�
    friend bool operator!=(const ImgRaw& lhs, const ImgRaw& rhs);
    friend bool operator==(const ImgRaw& lhs, const ImgRaw& rhs);
    // �g BMP ��
    void bmp(string name, size_t bits=0) {
		if (bits == 0) { bits = this->bitCount; }
        vector<unsigned char> img = (*this);// // �������ഫ�禡
        Raw::raw2bmp(name, img, width, height, bits);
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
    Sift(ImgRaw img);
    Sift(vector<types> raw_img, size_t width, size_t height):
        raw_img(raw_img, width, height) {}
public:
	void pyramid2();
    void comp(vector<ImgRaw>& pyrs, string name="");
	void addArrow(string name="feaArrow.bmp");
private:
	bool findMaxMin(vector<ImgRaw>& gauDog_imgs, size_t scale_idx, size_t curr_Width, size_t y, size_t x);
	void FeatureDescrip(vector<ImgRaw>& kaidaImag);
	void getHistogramMS(const ImgRaw& doImage, float Insize, size_t scale, float sigma, 
		size_t Iny, size_t Inx, size_t InWidth, size_t Inr);
	void AddnewFeaturestruct(int Inx, int Iny, float Insize, int kai, int sigmaOCT, float Inm, int Insita);
public:
    ImgRaw raw_img;  //���
	size_t pyWidth=5;  //��(��j�Y�p)
	size_t pyheight=5; //��e(�ҽk�X��)
    vector<vector<ImgRaw>> pyrs;
	vector<vector<ImgRaw>> pyrs_dog;
	// �S�x�I
	Feature* FeatureStart;
	Feature* FeatureEnd;
	Feature* FeatureNow;
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
	static void draw_line(ImgRaw& img, size_t y, size_t x, float line_len, float sg);
};
//-----------------------------------------------------------------
class Stitching{
private:
	int Width, Height;
	Feature *FeatureStart1, *FeatureStart2;
	ImgRaw matchImg;
	//RGBTRIPLE** color;
	//BITMAPFILEHEADER FileHeader;
	//BITMAPINFOHEADER InfoHeader;
public:
	//*** �e�ⶵ�ѼƬ����Y�ɡA�ĤT�B�|�����Ϥ���T�A���B�������S�x�I��T ***//
	Stitching(
		//BITMAPFILEHEADER Filedata, 
		//BITMAPINFOHEADER Infodata, 
		//RGBTRIPLE** image1, 
		//RGBTRIPLE** image2, 
		Feature* inFeatureptr1, 
		Feature* inFeatureptr2,
		string name1,
		string name2
	);
	~Stitching();
	void OutBMP(std::string outname);
	void OutRAW(std::string outname);
	//void WriteImage(RGBTRIPLE** image1, RGBTRIPLE** image2);
	//*** �ˬd�O�_���ۦP���S�x�y�z�l ***//
	void Check(void);
	//*** �N�a�J�����I�۳s ***//
	void Link(int x1, int y1, int x2, int y2);
	//*** �p���Ӵy�z�l�������ڦ��Z�� ***//
	float EuclideanDistance(std::vector <std::vector <std::vector <float>>> point1, std::vector <std::vector <std::vector <float>>> point2);
};
