/*****************************************************************
Name : imgraw
Date : 2017/11/16
By   : CharlotteHonG
Final: 2017/11/16
*****************************************************************/
#pragma warning(disable : 4819)
#pragma once

#include "imglib\imglib.hpp"
//-----------------------------------------------------------------
// �ֳtatan�B��
float fastAtan2f(float dy, float dx);
float fastAtan2f_rad(float dy, float dx);
float fastAtanf(float dy);
float fastAtanf_rad(float dy);


//-----------------------------------------------------------------
#define ImgRawPixMax 255.0

class ImgRaw {
private:
	using types = float;
public:
	// ��l��
	ImgRaw() = default;
	ImgRaw(vector<types> img, uint32_t width, uint32_t height, uint16_t bits) :
		raw_img(img), width(width), height(height), bitCount(bits) {}
	ImgRaw(vector<unsigned char> img, uint32_t width, uint32_t height, uint16_t bits) :
		width(width), height(height), bitCount(bits)
	{
		raw_img.resize(img.size());
		for(unsigned i = 0; i < img.size(); ++i) {
			if(nomal) {
				raw_img[i] = img[i]/ImgRawPixMax;
			} else {
				raw_img[i] = img[i];
			}
		}
	}
	ImgRaw(uint32_t width, uint32_t height, uint16_t bits) :raw_img(width*height * (bits/8)), 
		width(width), height(height), bitCount(bits){}
	ImgRaw(string bmpname, string path="");
	ImgRaw(string bmpname, string path, bool nomal);
	// �ƻs�禡
	/*ImgRaw& operator=(const ImgRaw& other) {
		if (this != &other) {
			raw_img = other.raw_img;
			width = other.width;
			height = other.height;
			bitCount = other.bitCount;
		}
		return *this;
	}*/
	// �����ഫ
	operator vector<types>&() { return raw_img; }
	operator const vector<types>&() const { return raw_img; }
	operator vector<unsigned char>() {
		vector<unsigned char> img(raw_img.size());
		for(unsigned i = 0; i < raw_img.size(); ++i) {
			if(nomal) {
				img[i] = (unsigned char)(raw_img[i]*ImgRawPixMax);
			} else {
				img[i] = (unsigned char)(raw_img[i]);
			}
		}
		return img;
		//const vector<unsigned char> img = static_cast<const ImgRaw&>(*this);
		//return const_cast<vector<unsigned char>&>(img);
	}
	operator const vector<unsigned char>() const {
		vector<unsigned char> img(raw_img.size());
		for(unsigned i = 0; i < raw_img.size(); ++i) {
			if(nomal) {
				img[i] = (unsigned char)(raw_img[i]*ImgRawPixMax);
			} else {
				img[i] = (unsigned char)(raw_img[i]);
			}
		}
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
	const size_t size() const { return this->raw_img.size(); }
	// ���]�j�p
	void resize(uint32_t width, uint32_t height, uint16_t bits) {
		raw_img.resize(width*height * bits/8);
		this->width=width;
		this->height=height;
		this->bitCount=bits;
	}
public:
	// �ର�Ƕ�
	ImgRaw ConverGray() const;
	// �g BMP ��
	void bmp(string name, uint32_t bits=0);
	void bmp(string name, uint32_t bits=0) const;
	// ���X����᪺�Ϥ�
	ImgRaw rotateImg(size_t x, size_t y, float radius, float sita);
public: // ��j�Y�p (ı�o���طQ����)
	static void zero(ImgRaw& tar, ImgRaw& sou, float z);
	static void first(ImgRaw& tar, ImgRaw& sou, float z);
	static void cubic(ImgRaw& tar, ImgRaw& sou, float z);
	static void gauBlur(ImgRaw& tar, ImgRaw& sou, float p);
public:
	vector<types> raw_img;
	uint32_t width;
	uint32_t height;
	uint16_t bitCount = 0;
	bool nomal = 1;
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



//-----------------------------------------------------------------
// �e�u
#define DrawPixNomal 255.0

class Draw {
public:
	static void drawLine_p(ImgRaw& img, int y, int x, int y2, int x2, float val=200);
	static void drawLine_s(ImgRaw& img, int y, int x, float line_len, float sg, float val=200);
	static void draw_arrow(ImgRaw& img, int y, int x, float line_len, float sg, float val=200);

	static void drawLineRGB_p(ImgRaw& img, int y, int x, int y2, int x2, 
		float r, float, float);
	static void drawLineRGB_p(ImgRaw& img, int y, int x, int y2, int x2);
	static void drawLineRGB_s(ImgRaw& img, int y, int x, float line_len, float sg);
	static void draw_arrowRGB(ImgRaw& img, int y, int x, float line_len, float sg);
};



//----------------------------------------------------------------