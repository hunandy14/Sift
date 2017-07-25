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

//================================================================
int main(int argc, char const *argv[]){
    // 讀取圖片
    vector<unsigned char> raw_img;
    vector<unsigned char> raw_img2x;


    // 圖片1
    // Raw::read_bmp(raw_img, "Lena.bmp");
    // Raw::raw2gray(raw_img); // 轉灰階
    // 圖片2
    Raw::read_raw(raw_img, "en_420x420.raw");
    // 轉浮點數
    auto raw2f = [](vector<unsigned char>& img) {
        vector<float> temp(size(img));
        for(unsigned i = 0; i < temp.size(); ++i)
            temp[i] = (float)img[i] / 255.0;
        return temp;
    };
    // 轉浮char
    // auto f2raw = [](vector<float>& img) {
    //     vector<unsigned char> temp(size(img));
    //     for(unsigned i = 0; i < temp.size(); ++i)
    //         temp[i] = img[i];
    //     return temp;
    // };

    // 創建結構
    ImgRaw temp(raw2f(raw_img), 420 , 420);
    ImgRaw input_img(0, 0);
    ImgRaw::first(input_img, temp, 0.5);
    // 金字塔
    Sift img(input_img);
    img.pyramid(1);
	return 0;
}
//================================================================
