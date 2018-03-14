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
// 快速atan運算
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
	// 初始化
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
	// 複製函式
	/*ImgRaw& operator=(const ImgRaw& other) {
		if (this != &other) {
			raw_img = other.raw_img;
			width = other.width;
			height = other.height;
			bitCount = other.bitCount;
		}
		return *this;
	}*/
	// 隱式轉換
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
	// 重載下標符號
	types & operator[](size_t idx) {
		return const_cast<types&>(static_cast<const ImgRaw&>(*this)[idx]);
	}
	const types & operator[](size_t idx) const { return raw_img[idx]; }
	// 二維讀取
	types & at2d(size_t y, size_t x) {
		return const_cast<types&>(static_cast<const ImgRaw&>(*this).at2d(y, x));
	}
	const types & at2d(size_t y, size_t x) const {
		return raw_img[y*width + x];
	}
	// 二維雙線性運算讀取
	const types atBilinear(float y, float x) const;
	// 大小是否相等
	friend bool operator!=(const ImgRaw& lhs, const ImgRaw& rhs);
	friend bool operator==(const ImgRaw& lhs, const ImgRaw& rhs);
	// 獲得大小
	const size_t size() const { return this->raw_img.size(); }
	// 重設大小
	void resize(uint32_t width, uint32_t height, uint16_t bits) {
		raw_img.resize(width*height * bits/8);
		this->width=width;
		this->height=height;
		this->bitCount=bits;
	}
public:
	// 轉為灰階
	ImgRaw ConverGray() const;
	// 寫 BMP 檔
	void bmp(string name, uint32_t bits=0);
	void bmp(string name, uint32_t bits=0) const;
	// 取出旋轉後的圖片
	ImgRaw rotateImg(size_t x, size_t y, float radius, float sita);
public: // 放大縮小 (覺得累贅想拿掉)
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
// 大小是否相等
inline bool operator!=(const ImgRaw& lhs, const ImgRaw& rhs) {
	return !(lhs == rhs);
}
inline bool operator==(const ImgRaw& lhs, const ImgRaw& rhs) {
	if (lhs.width == rhs.width &&  lhs.height == rhs.height) {
		return 1;
	} return 0;
}



//-----------------------------------------------------------------
// 畫線
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