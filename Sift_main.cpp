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

#define testpoint1
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


//#define harri_corner
#ifdef harri_corner
	ImgRaw fea_pix_ori("harri_ori.bmp");
	ImgRaw fea_pix(420, 200);
	/*
	// 黑白相反
	for (size_t j = 1; j < fea_pix_ori.height - 1; j++) {
		for (size_t i = 1; i < fea_pix_ori.width - 1; i++) {
			if (fea_pix_ori.at2d(j, i) == 1){
				fea_pix_ori.at2d(j, i) = 0;
			}
			else{
				fea_pix_ori.at2d(j, i) = 1;
			}
		}
	}
	fea_pix_ori.bmp("harri.bmp", 8);
	// 偵測
	for (size_t j = 1; j < fea_pix.height -1; j++) {
		for (size_t i = 1; i < fea_pix.width-1; i++) {
			if (fea_pix_ori.at2d(j, i)==1) { // 偵測白點
				fea_pix.at2d(j, i) = 0.5;
				if (Corner::harris(fea_pix_ori, fea_pix_ori.width, j, i) == 1) {
					fea_pix.at2d(j, i) = 0;
				}
			}
		}
	}
	fea_pix.bmp("harri.bmp", 8);
	*/

	// 黑白相反
	for (size_t j = 0; j < fea_pix_ori.height; j++) {
		for (size_t i = 0; i < fea_pix_ori.width; i++) {
			if (fea_pix_ori.at2d(j, i) == 1){
				fea_pix_ori.at2d(j, i) = 1 /255.0; // 背景
			}
			else{
				fea_pix_ori.at2d(j, i) = 200 /255.0; // 圖形
				//cout << fea_pix_ori.at2d(j, i) << ", ";
			}
		}
	}
	fea_pix_ori.bmp("harri.bmp", 8);
	// 偵測
	for (size_t j = 1; j < fea_pix.height -1; j++) {
		for (size_t i = 1; i < fea_pix.width-1; i++) {
			if (fea_pix_ori.at2d(j, i)*255 == 200) { // 偵測的圖形
				fea_pix.at2d(j, i) = 128 /255.0; // 圖形顏色
				if (Corner::harris(fea_pix_ori, fea_pix_ori.width, j, i) == 1) {
					fea_pix.at2d(j, i) = 20 /255.0; // 角點
				}
			} else {
				fea_pix.at2d(j, i) = 255 /255.0; // 背景
			}
		}
	}
	fea_pix.bmp("harri2.bmp", 8);


	/*
	ImgRaw cor(fea_pix_ori.width, fea_pix_ori.height, 255);
	for (size_t j = 1; j < fea_pix_ori.height -1; j++) {
		for (size_t i = 1; i < fea_pix_ori.width-1; i++) {
			if (fea_pix_ori.at2d(j, i)==0) { // 偵測黑點
				cor.at2d(j, i) = 0.5;
				if (Corner::harris(fea_pix_ori, fea_pix_ori.width, j, i) == 1) {
					cor.at2d(j, i) = 0;
				}
			}
		}
	}
	cor.bmp("harri.bmp", 8);
	*/

#endif // harri_corner

	return 0;

}
//================================================================
