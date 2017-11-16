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

#include "Imgraw.hpp"
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
	Raw::read_bmp(img, bmpname, &width, &height, &bits);
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
	if (x < 0 || y < 0 || 
		x>=width-1 || y >=height-1)
	{
		std::cerr << "atBilinear(x, y) : [x, y] our of range. \n";
		std::cerr << "x=" << x << "y=" << y << "\n";
		return 0;
	}
	// ����F�I(����� 1+)
	size_t x0 = floor(x);
	size_t x1 = ceil(x);
	size_t y0 = floor(y);
	size_t y1 = ceil(y);
	// ������(�u��� 1-)
	float dx1 = x - x0;
	float dx2 = 1 - dx1;
	float dy1 = y - y0;
	float dy2 = 1 - dy1;
	// ����I
	const float& A = raw_img[y0*width + x0];
	const float& B = raw_img[y0*width + x1];
	const float& C = raw_img[y1*width + x0];
	const float& D = raw_img[y1*width + x1];
	// ���X���(�n��e)
	float AB = A*dx2 + B*dx1;
	float CD = C*dx2 + D*dx1;
	float X = AB*dy2 + CD*dy1;
	return X;
}



// �e�u
void Draw::drawLine_p(ImgRaw& img, int y, int x, int y2, int x2) {
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
			img.raw_img[distY*img.width + distX] = 0.5;
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
			img.raw_img[distY*img.width + distX] = 0.5;
		}
	}
}
void Draw::drawLine_s(ImgRaw& img, int y, int x, float line_len, float sg) {
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
	drawLine_p(img, y, x, y2, x2);
}
void Draw::draw_arrow(ImgRaw& img, int y, int x, float line_len, float sg) {
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
	drawLine_p(img, y, x, y2, x2);
	// �e�Y
	drawLine_s(img, y2, x2, 10, sg-150);
	drawLine_s(img, y2, x2, 10, sg+150);
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
	drawLineRGB_p(img, y, x, y2, x2);
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
	drawLineRGB_p(img, y, x, y2, x2);
	// �e�Y
	size_t head_len = 6;
	drawLineRGB_s(img, y2, x2, head_len, sg-150);
	drawLineRGB_s(img, y2, x2, head_len, sg+150);
}
