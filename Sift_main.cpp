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
	

	//Mat kama_test(1000, 1334, 0);
	//namedWindow("Image", WINDOW_AUTOSIZE);
	//resizeWindow("Image", 1000, 1334);
	//kama_test.data = raw_img.data();
	//imshow("Image", kama_test);
	//cvWaitKey(0);


    // 圖片2
    // Raw::read_raw(raw_img, "en_420x420.raw");

    // 轉浮點數
	
    auto raw2f = [](vector<unsigned char>& img) {
        vector<float> temp(size(img));
        for(unsigned i = 0; i < temp.size(); ++i)
            temp[i] = (float)img[i] / 255.0;
        return temp;
    };
    // 創建結構(寬, 長)
    ImgRaw temp(raw2f(raw_img), 850 , 602);
    ImgRaw input_img(0, 0);
    ImgRaw::first(input_img, temp, 1);
    // 金字塔
    Sift img(input_img);
    img.pyramid(1);
	
	return 0;
}
//================================================================
