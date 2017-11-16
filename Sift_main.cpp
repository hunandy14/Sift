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

#define M_PI 3.14159265358979323846
#define SQUARE2 1.4142135623730951f

ImgRaw rotateImg(const ImgRaw& sou, float radius, float sita) {
	int ry = radius * SQUARE2;
	int rx = radius * SQUARE2;
	ImgRaw rotate(rx*2, ry*2, 8);

	float cos_t = cosf(sita*(float)(M_PI/180));  
	float sin_t = sinf(sita*(float)(M_PI/180));

	for (int j = -ry; j < ry; j++) {
		for (int i = -rx; i < rx; i++) {
			float r_rot = j*sin_t + i*cos_t; // 原圖X座標(未轉換)
			float c_rot = j*cos_t - i*sin_t; // 原圖Y座標(未轉換)

			float rbin = r_rot + radius; 
			float cbin = c_rot + radius;

			if (cbin < radius*2 - 1 and cbin > 0) {
				if (rbin < radius*2 - 1 and rbin > 0) {
					// 四捨五入
					//rotate.at2d(j+ry, i+rx) = sou.at2d(round(cbin), round(rbin));
					// 雙線姓插補
					rotate.at2d(j+ry, i+rx) = sou.atBilinear(cbin, rbin);
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
	string name1="Lena.bmp";
	string name2="Lena5.bmp";
	ImgRaw img1(name1);
	ImgRaw img2(name2);
	// 旋轉圖片

	/*ImgRaw sou = img1.ConverGray();
	ImgRaw rotate = rotateImg(sou, sou.width/2, 45);
	rotate.bmp("rotate.bmp");*/

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
	match.Check(0.7);
	


end = clock();
cout << "time is:" << (end - start)/1000.0 << "s" << endl;
#endif // testpoint1
	

	// 畫線測試
	/*ImgRaw img_line(1280, 720, 8);
	for (size_t i = 0; i < 36; i++) {
		Draw::draw_arrow(img_line, 200, 200, 100.f*sqrt(2), i*10);
	} img_line.bmp("line.bmp");*/
	
/*
	ImgRaw img_lineRGB(1280, 720, 24);
	for (size_t i = 0; i < 36; i++) {
		Draw::draw_arrowRGB(img_lineRGB, 200, 200, 100.f*sqrt(2), i*10);
	} img_lineRGB.bmp("lineRGB.bmp");
	*/
	return 0;
}

//================================================================
