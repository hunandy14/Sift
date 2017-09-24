/*****************************************************************
Name : 
Date : 2017/07/04
By   : CharlotteHonG
Final: 2017/07/04
*****************************************************************/
#define M_PI 3.14159265358979323846

#pragma warning(disable : 4819)
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;
#include "Sift.hpp"
#include <opencv2/opencv.hpp>
using namespace cv;



//================================================================
int main(int argc, char const *argv[]){
#define testpoint1
#ifdef testpoint1
	// 讀取圖片
	// ImgRaw img("en.bmp");
	// ImgRaw img("ball_01.bmp");
	ImgRaw img("kanna.bmp", 1);
    //ImgRaw input_img(0, 0);
    //ImgRaw::first(input_img, img, 1);
    // 金字塔
    Sift fea(img);
	fea.pyramid2();
	//fea.pyramid(2);
	fea.display();

#endif // testpoint1
	
	/*ImgRaw img_line("kanna.bmp", 1);
	for (size_t i = 0; i < 36; i++) {
		Draw::draw_line(img_line, 400, 300, 100.f*sqrt(2), i*10);
	}

	ImgRaw img_t(720, 480);
	Draw::draw_line(img_t, 40, 30, 100, 10);
	img_t.bmp("line2.bmp", 8);

	img_line.bmp("line.bmp", 8);
	*/
	return 0;

}

//================================================================
