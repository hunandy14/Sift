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
#include <Timer.hpp>

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



//================================================================
int main(int argc, char const *argv[]){
	//srand((unsigned)time(NULL)); rand();
	// 讀取圖片

	ImgRaw f("kanna.bmp");
	ImgRaw f2;
	f.first(f2, f, 4);

#define testpoint1
#ifdef testpoint1	
	string name1="sd01.bmp", name2="sd02.bmp";
	//string name1="Lena.bmp", name2="Lena4.bmp";
	
	ImgRaw img1(name1);
	ImgRaw img2(name2);


    // 金字塔1
	Timer t1;
	Sift fea1(img1);
	fea1.pyramid();
	fea1.drawArrow("feaArrow1.bmp");
	t1.print("feat1");

	// 金字塔2
	t1.start();
	Sift fea2(img2);
	fea2.pyramid();
	fea2.drawArrow("feaArrow2.bmp");
	t1.print("feat2");

	// KDtree


	// 匹配特徵點(兩張大小要一樣)
	Stitching match(fea1, fea2);
	match.Check(0.4);

#endif // testpoint1
	return 0;
}

//================================================================
