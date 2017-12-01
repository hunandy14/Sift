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



//================================================================
int main(int argc, char const *argv[]){
	//srand((unsigned)time(NULL)); rand();
	// 讀取圖片
	string name1="Lena6.bmp";
	string name2="Lena135.bmp";
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
	fea2.drawArrow("feaArrow2.bmp");
	// 匹配特徵點(兩張大小要一樣)
	Stitching match(fea1, fea2);
	match.Check(0.5);
end = clock();
cout << "time is:" << (end - start)/1000.0 << "s" << endl;
#endif // testpoint1
	return 0;
}

//================================================================
