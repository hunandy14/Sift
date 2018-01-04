/*****************************************************************
Name : imgraw
Date : 2017/11/16
By   : CharlotteHonG
Final: 2017/11/16
*****************************************************************/
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "imglib\imglib.hpp"
#include "Raw2Img\Raw2Img.hpp"
#include "Imgraw.hpp"

#if defined(_MSC_VER)
	#define or ||
	#define and &&
	#define OR ||
	#define AND &&
#endif

#define M_PI 3.14159265358979323846



// �ֳtatan�B��
float fastAtan2f(float dy, float dx){
	static const float atan2_p1 = 0.9997878412794807f*(float)(180/M_PI);
	static const float atan2_p3 = -0.3258083974640975f*(float)(180/M_PI);
	static const float atan2_p5 = 0.1555786518463281f*(float)(180/M_PI);
	static const float atan2_p7 = -0.04432655554792128f*(float)(180/M_PI);
	static const float atan2_DBL_EPSILON = 2.2204460492503131e-016;

	float ax = std::abs(dx), ay = std::abs(dy);
	float a, c, c2;
	if (ax >= ay) {
		c = ay/(ax + static_cast<float>(atan2_DBL_EPSILON));
		c2 = c*c;
		a = (((atan2_p7*c2 + atan2_p5)*c2 + atan2_p3)*c2 + atan2_p1)*c;
	} else {
		c = ax/(ay + static_cast<float>(atan2_DBL_EPSILON));
		c2 = c*c;
		a = 90.f - (((atan2_p7*c2 + atan2_p5)*c2 + atan2_p3)*c2 + atan2_p1)*c;
	}
	if (dx < 0)
		a = 180.f - a;
	if (dy < 0)
		a = 360.f - a;
	return a;
}



// ImgRaw �غc�l
ImgRaw::ImgRaw(string bmpname){
	vector<unsigned char> img;
	uint32_t width, height;
	uint16_t bits;
	// Ū���Ϥ�
	Raw2Img::read_bmp(img, bmpname, &width, &height, &bits);
	this->width    = width;
	this->height   = height;
	this->bitCount = bits;
	// ��l��(�t���W��)
	raw_img.resize(img.size());
	for (size_t i = 0; i < img.size(); i++) {
		raw_img[i] = (float)img[i] / 255.0;
	}
}
// �G�����u�ʹB��Ū��
const ImgRaw::types ImgRaw::atBilinear(float y, float x) const {
	if (x < 0 || y < 0 || x>=width-1 || y >=height-1){
		std::cerr << "atBilinear(x, y) : [x, y] our of range. \n";
		std::cerr << "x=" << x << "y=" << y << "\n";
		return 0;
	}
	return Scaling::bilinear(raw_img, width, y, x);
}
// ���X����᪺�Ϥ�
ImgRaw ImgRaw::rotateImg(size_t x, size_t y, float radius, float sita) {
	ImgRaw& Img = *this;
	// �e�b�Y�аO
	//Draw::draw_arrow(Img, y, x, radius, sita, 1);
	//Img.bmp("rotate_test.bmp");
	// �s�Ϫ��e�b�|
	float maxRadius = radius;
	int rotH = floor(maxRadius);
	int rotW = floor(maxRadius);
	ImgRaw rotate(rotW*2, rotH*2, 8);
	// �w��
	sita *= -1; // ��s����^0��
	float cos_t = cosf(sita*(float)(M_PI/180));  
	float sin_t = sinf(sita*(float)(M_PI/180));
	// �]�s��
	for (int j = -rotH; j < rotH; j++) {
		for (int i = -rotW; i < rotW; i++) {
			// ��J�s�Ϯy�Ъ�^�¹Ϯy��(�w0, 0����߱���)
			float r_rot = (j)*sin_t + (i)*cos_t; // ���X�y��
			float c_rot = (j)*cos_t - (i)*sin_t; // ���Y�y��
												 // �B���^���w��m
			float rbin = r_rot + x; 
			float cbin = c_rot + y;
			// �h����ϥ~���I
			if (cbin < Img.height - 1 and cbin > 0) {
				if (rbin < Img.width - 1 and rbin > 0) {
					// ���u�m����
					rotate.at2d(j+rotH, i+rotW) = Img.atBilinear(cbin, rbin); 
				}
			}
		}
	}
	return rotate;
}




// �e�u
void Draw::drawLine_p(ImgRaw& img, int y, int x, int y2, int x2, float val) {
	// ���I�������Z���t
	float dx = x2-x;
	float dy = y2-y;
	// �HY�b���D
	float sita=fastAtan2f(dy, dx);
	if (sita>45 and sita<135 or sita>225 and sita<315) {
		float slopeY = dx/dy; // �ײv
		for (int i = 0; i < abs(dy); i++) {
			int iFix = dy>0? i:-i;
			int currPos = iFix*slopeY + x;

			int distX = currPos;
			int distY = y+iFix;

			if (distX<0 or distX>=img.width or distY<0 or distY>=img.height) {
				return;
			}
			img.raw_img[distY*img.width + distX] = val;
		}
	} 
	// �HX�b���D
	else {
		float slopeX = dy/dx; // �ײv
		for (int i = 0; i < abs(dx); i++) {
			int iFix = dx>0? i:-i;
			int currPos = iFix*slopeX + y;

			int distX = x+iFix;
			int distY = currPos;

			if (distX<0 or distX>=img.width or distY<0 or distY>=img.height) {
				return;
			}
			img.raw_img[distY*img.width + distX] = val;
		}
	}
}
void Draw::drawLine_s(ImgRaw& img, int y, int x, float line_len, float sg, float val) {
	// ���b
	if (line_len < 0) {
		return;
	}
	if (line_len==1) {
		img[x*img.width + y] = val;
		return;
	}
	// ���Y��
	int x2 = x + line_len*cos(sg * M_PI/180.0);
	int y2 = y + line_len*sin(sg * M_PI/180.0);
	// �e�u
	drawLine_p(img, y, x, y2, x2, val);
}
void Draw::draw_arrow(ImgRaw& img, int y, int x, float line_len, float sg, float val) {
	// ���b
	if (line_len < 0) {
		return;
	}
	if (line_len==1) {
		img[x*img.width + y] = val;
		return;
	}
	// ���Y��
	int x2 = x + line_len*cos(sg * M_PI/180.0);
	int y2 = y + line_len*sin(sg * M_PI/180.0);
	// �e�u
	drawLine_p(img, y, x, y2, x2, val);
	// �e�Y
	drawLine_s(img, y2, x2, 10, sg-150, val);
	drawLine_s(img, y2, x2, 10, sg+150, val);
}
// �e�uRGB
void Draw::drawLineRGB_p(ImgRaw& img, int y, int x, int y2, int x2, 
	float r, float g, float b) {
	// ���I�������Z���t
	float dx = x2-x;
	float dy = y2-y;
	// �HY�b���D
	float sita=fastAtan2f(dy, dx);
	if (sita>45 and sita<135 or sita>225 and sita<315) {
		float slopeY = dx/dy; // �ײv
		for (int i = 0; i < abs(dy); i++) {
			int iFix = dy>0? i:-i;
			int currPos = iFix*slopeY + x;

			int distX = currPos;
			int distY = y+iFix;

			if (distX<0 or distX>=img.width or distY<0 or distY>=img.height) {
				return;
			}
			size_t posi = distY*img.width + distX;
			img.raw_img[posi*3 + 0] = r;
			img.raw_img[posi*3 + 1] = g;
			img.raw_img[posi*3 + 2] = b;
		}
	} 
	// �HX�b���D
	else {
		float slopeX = dy/dx; // �ײv
		for (int i = 0; i < abs(dx); i++) {
			int iFix = dx>0? i:-i;
			int currPos = iFix*slopeX + y;

			int distX = x+iFix;
			int distY = currPos;

			if (distX<0 or distX>=img.width or distY<0 or distY>=img.height) {
				return;
			}
			size_t posi = distY*img.width + distX;
			img.raw_img[posi*3 + 0] = r;
			img.raw_img[posi*3 + 1] = g;
			img.raw_img[posi*3 + 2] = b;
		}
	}
}
void Draw::drawLineRGB_p(ImgRaw& img, int y, int x, int y2, int x2) {
	// �H���C��
	auto random_num = [] {
		return ((rand() / (RAND_MAX+1.0)) * (1 - 0) + 0);
	};
	float rVal = random_num();
	float gVal = random_num();
	float bVal = random_num();
	drawLineRGB_p(img, y, x, y2, x2, rVal, gVal, bVal);
}
void Draw::drawLineRGB_s(ImgRaw& img, int y, int x, float line_len, float sg) {
	float value = 200 /255.0;
	float endvalue = 255 /255.0;
	// ���b
	if (line_len < 0) {
		return;
	}
	if (line_len==1) {
		img[x*img.width + y] = value;
		return;
	}
	// ���Y��
	int x2 = x + line_len*cos(sg * M_PI/180.0);
	int y2 = y + line_len*sin(sg * M_PI/180.0);
	// �e�u
	drawLineRGB_p(img, y, x, y2, x2, 242/255.f, 66/255.f, 54/255.f);
}
void Draw::draw_arrowRGB(ImgRaw& img, int y, int x, float line_len, float sg) {
	float value = 200 /255.0;
	float endvalue = 255 /255.0;
	// ���b
	if (line_len < 0) {
		return;
	}
	if (line_len==1) {
		img[x*img.width + y] = value;
		return;
	}
	// ���Y��
	int x2 = x + line_len*cos(sg * M_PI/180.0);
	int y2 = y + line_len*sin(sg * M_PI/180.0);
	// �e�u
	drawLineRGB_p(img, y, x, y2, x2, 242/255.f, 66/255.f, 54/255.f);
	// �e�Y
	size_t head_len = 6;
	drawLineRGB_s(img, y2, x2, head_len, sg-150);
	drawLineRGB_s(img, y2, x2, head_len, sg+150);
}
