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

#define testpoint1
#ifdef testpoint1
	// 讀取圖片
	ImgRaw temp("kanna.bmp");
    ImgRaw input_img(0, 0);
    ImgRaw::first(input_img, temp, 1);
    // 金字塔
    Sift img(input_img);
	img.pyramid(1);

#endif // testpoint1

//#define harri_corner
#ifdef harri_corner
	vector<unsigned char> fea_img;
	Raw::read_bmp(fea_img, "harri_ori.bmp");
	ImgRaw fea_pix_ori(raw2f(fea_img), 420, 200);
	ImgRaw fea_pix(raw2f(fea_img), 420, 200);

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
	for (size_t j = 1; j < fea_pix.height -1; j++) {
		for (size_t i = 1; i < fea_pix.width-1; i++) {
			if (fea_pix_ori.at2d(j, i)==1) {
				fea_pix.at2d(j, i) = 0.5;
				if (Corner::harris(fea_pix_ori, fea_pix_ori.width, j, i) == 1) {
					fea_pix.at2d(j, i) = 0;
				}
			}
			else {
				//fea_pix.at2d(j, i) = 0.5;
			}
		}
	}
	fea_pix.bmp("harri.bmp", 8);
#endif // harri_corner

	return 0;

}
//================================================================
