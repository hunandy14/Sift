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
	// 轉浮點數+正規化
    auto raw2f = [](vector<unsigned char>& img) {
        vector<float> temp(size(img));
        for(unsigned i = 0; i < temp.size(); ++i)
            temp[i] = (float)img[i] / 255.0;
        return temp;
    };

//#define testpoint1
#ifdef testpoint1
	// 讀取圖片
	// ImgRaw temp("en.bmp");
	// ImgRaw temp("ball_01.bmp");
	ImgRaw temp("kanna.bmp");
    ImgRaw input_img(0, 0);
    ImgRaw::first(input_img, temp, 1);
    // 金字塔
    Sift img(input_img);
	img.pyramid();

#endif // testpoint1
	
	ImgRaw img_line("kanna.bmp", 1);
	for (size_t i = 0; i < 36; i++) {
		Draw::draw_line(img_line, 400, 300, 100.f*sqrt(2), i*10);
	}

	img_line.bmp("line.bmp", 8);
	return 0;

}

//================================================================
