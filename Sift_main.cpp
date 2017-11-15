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

using namespace std;
#include "Sift.hpp"
#include <opencv2/opencv.hpp>
using namespace cv;

#define M_PI 3.14159265358979323846


static inline void ImgRotate(float& deltaX, float& deltaY, float sita) {
	// 座標轉換
	float X = cos(sita)*deltaX - sin(sita)*deltaY;
	deltaY = sin(sita)*deltaX + cos(sita)*deltaY;
	deltaX = X;
}
static inline void ImgRotate2(float& x, float& y, float sita) {
	// 座標轉換
	float x2;
	x2 = x*cos(sita) + y*sin(sita);
	y  = y*cos(sita) - x*sin(sita);
	x  = x2;
}
//================================================================
int main(int argc, char const *argv[]){
	//srand((unsigned)time(NULL)); rand();
	// 讀取圖片
	string name1="Lena.bmp";
	string name2="02.bmp";
	ImgRaw img1(name1);
	ImgRaw img2(name2);

	ImgRaw sou = img1.ConverGray();
	ImgRaw rotate(sou.width, sou.height, sou.bitCount);
	
	float ori_X = 120;
	float ori_Y = 120;
	float si= 180;

	// 轉正
	cout << "X=" << ori_X << ", Y=" << ori_Y << endl;
	if (ori_X > sou.width) {
		ori_X = sou.width/2.0 - ori_X;
	} else {
		ori_X = ori_X - sou.width/2.0;
	}
	if (ori_Y > sou.height) {
		ori_Y = sou.height/2.0 - ori_Y;
	} else {
		ori_Y = ori_Y - sou.height/2.0;
	}
	ImgRotate(ori_X, ori_Y, si * M_PI/180.0);
	ori_X += sou.width/2.0 - 1, ori_Y += sou.height/2.0 - 1;
	if (ori_Y < 0) {ori_Y=0;}
	if (ori_X < 0) {ori_X=0;}
	cout << "rat=" << si << ", X=" << ori_X << ", Y=" << ori_Y << endl;
	

	// for 要跑新圖才能插值運算
	for (float j = 0; j < rotate.height; j++) {
		for (float i = 0; i < rotate.width; i++) {
			float X=i, Y=j;
			
			// 轉正
			//cout << "X=" << (int)X << ", Y=" << (int)Y << endl;
			if (X > sou.width) {
				X = sou.width/2.0 - X;
			} else {
				X = X - sou.width/2.0;
			}
			if (Y > sou.height) {
				Y = sou.height/2.0 - Y;
			} else {
				Y = Y - sou.height/2.0;
			}

			ImgRotate2(X, Y, 45 * M_PI/180.0);

			X += sou.width/2.0 - 1, Y += sou.height/2.0 - 1;
			if (Y < 0) {Y=0;}
			if (X < 0) {X=0;}
			//cout << "rat=" << si << ", X=" << (int)X << ", Y=" << (int)Y << endl;

			if ((int)Y < rotate.height and (int)Y > 0) {
				if ((int)X < rotate.width and (int)X > 0) {
					rotate.at2d(j, i) = sou.at2d(Y, X);
				}
			}

			
			//rotate.at2d(j, i) = 0.5;
		}
	}
	rotate.bmp("rotate.bmp");
	



// #define testpoint1
#ifdef testpoint1
clock_t start,end;
start = clock();


    // 金字塔1
    Sift fea1(img1);
	fea1.pyramid();
	fea1.drawArrow2("feaArrow1.bmp");
	// 金字塔2
	Sift fea2(img2);
	fea2.pyramid();
	fea2.drawArrow2("feaArrow2.bmp");
	// 匹配特徵點(兩張大小要一樣)
	Stitching match(fea1, fea2);
	match.Check(0.75);
	


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
