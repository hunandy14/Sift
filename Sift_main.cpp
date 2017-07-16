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

string rawName = "Rena_256x256_8bit.raw";
string rawName2 = "Seymour_Park_960x540_24bit.raw";
string bmpName = "Rena.bmp";
string bmpName2 = "Rena2.bmp";
//================================================================
int main(int argc, char const *argv[]){
    // 讀取圖片
    vector<unsigned char> raw_img;
    vector<unsigned char> raw_img2x;
    Raw::read_bmp(raw_img, "conna.bmp");
    Raw::raw2gray(raw_img); // 轉灰階
    // 轉浮點數
    auto tof = [](vector<unsigned char>& img) {
        vector<float> f_img(size(img));
        for(unsigned i = 0; i < f_img.size(); ++i)
            f_img[i] = img[i];
        return f_img;
    };
    // 轉浮char
    auto toch = [](vector<float>& img) {
        vector<unsigned char> f_img(size(img));
        for(unsigned i = 0; i < f_img.size(); ++i)
            f_img[i] = img[i];
        return f_img;
    };
    // 創建結構
    vector<float> f_img = tof(raw_img);
    ImgRaw s_img(f_img, 850, 602);
    ImgRaw s_img2(850*2, 602*2);

// 測試成功zero可以轉
    // Scaling::zero(raw_img2x, raw_img, 850, 602, 2);
    // Scaling::zero(s_img2, s_img, 850, 602, 2);


    ImgRaw::first(s_img2, s_img, 2);
    // vector<unsigned char> b=toch(s_img2);
    // Raw::raw2bmp("2X.bmp", b, 1700, 1204, 8);
    s_img2.bmp("2x.bmp", 8);









    vector<float> raw_img2;
    // raw_img.resize(raw_img2.size());
    // for(unsigned i = 0; i < raw_img2.size(); ++i) {
    //     raw_img[i] = raw_img2[i];
    // }
    // Raw::raw2bmp("2x.bmp", raw_img, 850*2, 602*2, 24);


    // s_img.bmp("2x.bmp", 24);
    // raw_img=s_img;
    // ImgRaw::zero(s_img2, s_img, 2);
    // s_img2.bmp("2x.bmp", 24);
    // vector<unsigned char> i = s_img;
    // Sift img(raw_img, 420, 420);
    // img.pyramid();
	return 0;
}
//================================================================
