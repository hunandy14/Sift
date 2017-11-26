/*****************************************************************
Name : 
Date : 2017/07/04
By   : CharlotteHonG
Final: 2017/07/04
*****************************************************************/
#pragma warning(disable : 4819)
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <memory>
//#include <windows.h>

using namespace std;
#include "Sift.hpp"
#include <opencv2/opencv.hpp>
using namespace cv;

#if defined(_MSC_VER)
	#define or ||
	#define and &&
	#define OR ||
	#define AND &&
#endif
#define M_PI 3.14159265358979323846
#define SQUARE2 1.4142135623730951f


ImgRaw rotateImg(const ImgRaw& sou, size_t x, size_t y, float radius, float sita) {
	ImgRaw test = sou;
	//test.at2d(x, y) = 1;
	Draw::draw_arrow(test, y, x, radius, sita, 1);
	test.bmp("rotate_test.bmp");
	// 新圖長寬半徑
	float maxRadius = radius;
	int rotH = floor(maxRadius);
	int rotW = floor(maxRadius);
	ImgRaw rotate(rotW*2, rotH*2, 8);
	// 預算
	sita *= -1; // 把新圖轉回0度
	float cos_t = cosf(sita*(float)(M_PI/180));  
	float sin_t = sinf(sita*(float)(M_PI/180));
	// 跑新圖
	for (int j = -rotH; j < rotH; j++) {
		for (int i = -rotW; i < rotW; i++) {
			// 輸入新圖座標返回舊圖座標(已0, 0為圓心旋轉)
			float r_rot = (j)*sin_t + (i)*cos_t; // 原圖X座標
			float c_rot = (j)*cos_t - (i)*sin_t; // 原圖Y座標
			// 矯正回指定位置
			float rbin = r_rot + x; 
			float cbin = c_rot + y;
			// 去除原圖外的點
			if (cbin < sou.height - 1 and cbin > 0) {
				if (rbin < sou.width - 1 and rbin > 0) {
					// 雙線姓插補
					rotate.at2d(j+rotH, i+rotW) = test.atBilinear(cbin, rbin); 
				}
			}
		}
	}
	return rotate;
}

ImgRaw rotateImg2(const ImgRaw& sou, size_t x, size_t y, float radius, float sita) {
	ImgRaw test = sou;
	//test.at2d(x, y) = 1;
	Draw::draw_arrow(test, y, x, radius, sita, 1);
	test.bmp("rotate_test.bmp");
	// 新圖長寬半徑
	float maxRadius = radius;
	int rotH = floor(maxRadius);
	int rotW = floor(maxRadius);
	ImgRaw rotate(rotW*2, rotH*2, 8);
	// 預算
	sita *= -1; // 把新圖轉回0度
	float cos_t = cosf(sita*(float)(M_PI/180));  
	float sin_t = sinf(sita*(float)(M_PI/180));
	// 跑新圖
	for (int j = -rotH; j < rotH; j++) {
		for (int i = -rotW; i < rotW; i++) {
			// 輸入新圖座標返回舊圖座標(已0, 0為圓心旋轉)
			float r_rot = (j)*sin_t + (i)*cos_t; // 原圖X座標
			float c_rot = (j)*cos_t - (i)*sin_t; // 原圖Y座標
												 // 矯正回指定位置
			float rbin = r_rot + x; 
			float cbin = c_rot + y;
			// 去除原圖外的點
			if (cbin < sou.height - 1 and cbin > 0) {
				if (rbin < sou.width - 1 and rbin > 0) {
					// 雙線姓插補
					rotate.at2d(j+rotH, i+rotW) = test.atBilinear(cbin, rbin); 
				}
			}
		}
	}
	return rotate;
}
//================================================================
int main(int argc, char const *argv[]){
	//srand((unsigned)time(NULL)); rand();
	// 讀取圖片
	string name1="iori\\iori01.bmp";
	string name2="iori\\iori02.bmp";
	ImgRaw img1(name1);
	ImgRaw img2(name2);
	// 旋轉圖片
	ImgRaw sou = img1.ConverGray();
	ImgRaw rotate = rotateImg(sou, 80, 95, 100/2.0, 45);
	//rotate.bmp("rotate.bmp");
	// 放大圖片
	ImgRaw first = sou.ConverGray();
	//ImgRaw::first(first, sou, 0.3595);
	//first.bmp("first.bmp", 8);


#define testpoint1
#ifdef testpoint1
clock_t start,end;
start = clock();
    // 金字塔1
    Sift fea1(img1);
	fea1.pyramid();
	fea1.drawArrow("feaArrow1.bmp");
	// 金字塔2
	Sift fea2(img2);
	fea2.pyramid();
	fea2.drawArrow("feaArrow2.bmp");
	// 匹配特徵點(兩張大小要一樣)
	Stitching match(fea1, fea2);
	match.Check(0.4);
end = clock();
cout << "time is:" << (end - start)/1000.0 << "s" << endl;
#endif // testpoint1
	return 0;
}

//================================================================
