/*****************************************************************
Name : imgraw
Date : 2017/11/16
By   : CharlotteHonG
Final: 2017/11/16
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

//-----------------------------------------------------------------
// �ֳtatan�B��
float fastAtan2f(float dy, float dx);



//-----------------------------------------------------------------
class ImgRaw {
private:
	using types = float;
public:
	// ��l��
	ImgRaw() = default;
	ImgRaw(vector<types> img, size_t width, size_t height, size_t bits) :
		raw_img(img), width(width), height(height), bitCount(bits) {}
	ImgRaw(size_t width, size_t height, size_t bits) :raw_img(width*height * (bits/8)), 
		width(width), height(height), bitCount(bits){}
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
	const types & operator[](size_t idx) const { return raw_img[idx]; }
	// �G��Ū��
	types & at2d(size_t y, size_t x) {
		return const_cast<types&>(static_cast<const ImgRaw&>(*this).at2d(y, x));
	}
	const types & at2d(size_t y, size_t x) const {
		return raw_img[y*width + x];
	}
	// �G�����u�ʹB��Ū��
	const types atBilinear(float y, float x) const;
	// �j�p�O�_�۵�
	friend bool operator!=(const ImgRaw& lhs, const ImgRaw& rhs);
	friend bool operator==(const ImgRaw& lhs, const ImgRaw& rhs);
	// ��o�j�p
	size_t size() { return this->raw_img.size(); }
	// ���]�j�p
	void resize(size_t width, size_t height, size_t bits) {
		raw_img.resize(width*height * bits/8);
		this->width=width;
		this->height=height;
		this->bitCount=bits;
	}
public:
	// �ର�Ƕ�
	ImgRaw ConverGray() const;
	// �g BMP ��
	void bmp(string name, size_t bits=0);
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
// �g BMP ��
inline void ImgRaw::bmp(string name, size_t bits) {
	if (bits == 0) { bits = this->bitCount; }
	vector<unsigned char> img = (*this);// �������ഫ�禡
	Raw::raw2bmp(name, img, width, height, bits);
}
// �ର�Ƕ�
inline ImgRaw ImgRaw::ConverGray() const {
	if (bitCount == 24) {
		ImgRaw gray(this->width, this->height, 8);
		for (size_t i = 0; i < gray.size(); i++) {
			const types& R = raw_img[i*3+0];
			const types& G = raw_img[i*3+1];
			const types& B = raw_img[i*3+2];
			gray[i] = R*0.299 + G*0.587 + B*0.114;
		} return gray;
	} else if (bitCount == 8) {
		return (*this);
	}
}



//-----------------------------------------------------------------
// �e�u
class Draw {
public:
	static void drawLine_p(ImgRaw& img, int y, int x, int y2, int x2);
	static void drawLine_s(ImgRaw& img, int y, int x, float line_len, float sg);
	static void draw_arrow(ImgRaw& img, int y, int x, float line_len, float sg);

	static void drawLineRGB_p(ImgRaw& img, int y, int x, int y2, int x2, 
		float r=1, float g=0, float b=0);
	static void drawLineRGB_s(ImgRaw& img, int y, int x, float line_len, float sg);
	static void draw_arrowRGB(ImgRaw& img, int y, int x, float line_len, float sg);
};



//----------------------------------------------------------------