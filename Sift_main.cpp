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



//================================================================
int main(int argc, char const *argv[]){
	//srand((unsigned)time(NULL)); rand();
	// 讀取圖片
	string name1="Lena.bmp";
	string name2="Lena3.bmp";
	ImgRaw img1(name1);
	ImgRaw img2(name2);


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
	fea2.drawArrow2("feaArrow2.bmp");
	// 匹配特徵點(兩張大小要一樣)
	Stitching match(fea1, fea2);
	match.Check(0.6);
	


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
