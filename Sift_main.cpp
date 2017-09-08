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

// #define testpoint1
#ifndef testpoint1

	//IplImage*  -> 轉換 Mat
	//IplImage* img01;
	//img01 = cvLoadImage("kanna.png", 1);

	//Mat kanna = imread("kanna.png");
	//kanna.at<Vec3b>(30, 20)[0] = 255;


    // 讀取圖片
    vector<unsigned char> raw_img;
    vector<unsigned char> raw_img2x;
    // 圖片1
    Raw::read_bmp(raw_img, "kanna.bmp");
    Raw::raw2gray(raw_img); // 轉灰階
	// 圖片2
	//Raw::read_bmp(raw_img, "en.bmp");

	//Mat kama_test(1000, 1334, 0);
	//namedWindow("Image", WINDOW_AUTOSIZE);
	//resizeWindow("Image", 1000, 1334);
	//kama_test.data = raw_img.data();
	//imshow("Image", kama_test);
	//cvWaitKey(0);


    // 圖片2
    // Raw::read_raw(raw_img, "en_420x420.raw");


    // 創建結構(寬, 長)
    //ImgRaw temp(raw2f(raw_img), 420, 420);
	ImgRaw temp(raw2f(raw_img), 850, 602);
    ImgRaw input_img(0, 0);
    ImgRaw::first(input_img, temp, 1);
    // 金字塔
    Sift img(input_img);
	img.pyramid(1);

#endif // testpoint1


	
	vector<unsigned char> fea_img;
	/*
	Raw::read_bmp(fea_img, "fea/fea0-1.bmp");
	ImgRaw fea_pix(raw2f(fea_img), 850, 602);
	*/
	
	Raw::read_bmp(fea_img, "tt.bmp");
	ImgRaw fea_pix_ori(raw2f(fea_img), 420, 200);
	ImgRaw fea_pix(raw2f(fea_img), 420, 200);
	
	auto harri = [&](ImgRaw& p, size_t y, size_t x) {
		constexpr float r = 10;
		float Dxx = 2*p.at2d(y, x) - p.at2d(y, x-2) - p.at2d(y, x+1);
		float Dyy = 2*p.at2d(y, x) - p.at2d(y-1, x) - p.at2d(y+1, x);
		float Dxy = p.at2d(y+1, x+1) + p.at2d(y-1, x-1)
			- p.at2d(y-1, x+1) - p.at2d(y+1, x-1);
		Dxy /= 4;
		float Tr = Dxx + Dyy;
		float Det = Dxx * Dyy - Dxy*Dxy;

		if ((Tr*Tr/Det) < ((r+1)*(r+1))/r){
			return 1;
		} return 0;
	};

	for (size_t j = 1; j < fea_pix_ori.height - 1; j++) {
		for (size_t i = 1; i < fea_pix_ori.width - 1; i++) {
			if (fea_pix_ori.at2d(j, i) == 1)
			{
				fea_pix_ori.at2d(j, i) = 0;
			}
			else
			{
				fea_pix_ori.at2d(j, i) = 1;
			}
		}
	}
	for (size_t j = 1; j < fea_pix.height -1; j++){
		for (size_t i = 1; i < fea_pix.width-1; i++){
			if (fea_pix_ori.at2d(j, i)==1) {
				fea_pix.at2d(j, i) = 0.5;
				if (harri(fea_pix_ori, j, i) == 1)
				{
					fea_pix.at2d(j, i) = 0;
				}
			}
			else {
				//fea_pix.at2d(j, i) = 0.5;
			}
		}
	}
	fea_pix.bmp("harri.bmp", 8);
	
	return 0;

}
//================================================================
