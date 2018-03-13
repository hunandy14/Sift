/*****************************************************************
Name : imgraw
Date : 2017/11/16
By   : CharlotteHonG
Final: 2017/11/16
*****************************************************************/
#pragma warning(disable : 4819)
#pragma once

//-----------------------------------------------------------------
// �ֳtatan�B��
float fastAtan2f(float dy, float dx);
float fastAtan2f_rad(float dy, float dx);
float fastAtanf(float dy);
float fastAtanf_rad(float dy);


//-----------------------------------------------------------------
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
		for(unsigned i = 0; i < img.size(); ++i)
			raw_img[i] = img[i]/255.0;
	}
	ImgRaw(uint32_t width, uint32_t height, uint16_t bits) :raw_img(width*height * (bits/8)), 
		width(width), height(height), bitCount(bits){}
	ImgRaw(string bmpname, string path="");
	// �ƻs�禡
	ImgRaw& operator=(const ImgRaw& other) {
		if (this != &other) {
			raw_img = other.raw_img;
			width = other.width;
			height = other.height;
			bitCount = other.bitCount;
		}
		return *this;
	}
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
inline void ImgRaw::bmp(string name, uint32_t bits) {
	if (bits == 0) { bits = this->bitCount; }
	vector<unsigned char> img = (*this);// �������ഫ�禡
	Raw2Img::raw2bmp(name, img, width, height, bits);
}
inline void ImgRaw::bmp(string name, uint32_t bits) const {
	if (bits == 0) { bits = this->bitCount; }
	vector<unsigned char> img = (*this);// �������ഫ�禡
	Raw2Img::raw2bmp(name, img, width, height, bits);
}
// �ର�Ƕ�
inline ImgRaw ImgRaw::ConverGray() const {
	if (bitCount == 24) {
		ImgRaw gray(this->width, this->height, 8);
		for (size_t i = 0; i < gray.size(); i++) {
			const types& R = raw_img[i*3+0];
			const types& G = raw_img[i*3+1];
			const types& B = raw_img[i*3+2];
			gray[i] = (float)(R*0.299 + G*0.587 + B*0.114);
		} return gray;
	} else if (bitCount == 8) {
		return (*this);
	}
}



//-----------------------------------------------------------------
// �e�u
class Draw {
public:
	static void drawLine_p(ImgRaw& img, int y, int x, int y2, int x2, float val=200/255.0);
	static void drawLine_s(ImgRaw& img, int y, int x, float line_len, float sg, float val=200/255.0);
	static void draw_arrow(ImgRaw& img, int y, int x, float line_len, float sg, float val=200/255.0);

	static void drawLineRGB_p(ImgRaw& img, int y, int x, int y2, int x2, 
		float r, float, float);
	static void drawLineRGB_p(ImgRaw& img, int y, int x, int y2, int x2);
	static void drawLineRGB_s(ImgRaw& img, int y, int x, float line_len, float sg);
	static void draw_arrowRGB(ImgRaw& img, int y, int x, float line_len, float sg);
};



//----------------------------------------------------------------