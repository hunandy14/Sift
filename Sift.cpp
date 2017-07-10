/*****************************************************************
Name : 
Date : 2017/07/05
By   : CharlotteHonG
Final: 2017/07/05
*****************************************************************/
#include <iostream>
#include <vector>
#include <cmath>
#include "Sift.hpp"
using namespace std;

void Sift::pyramid(size_t s){
    s+=3;
    // size_t octvs = log(min(width, height)) / log(2.0) - 2;
    
    // 測試層數先設3
    size_t octvs = 3;

    vector<v_uch> pyr_img;
    pyr_img.resize(octvs*s);
    // 放大
    v_uch temp;
    Scaling::cubic(temp, raw_img, width, height, 2);
    height*=2, width*=2;
    raw_img=temp;
    // 模糊
    float p=1.6;
    GauBlur::raw2GauBlur(pyr_img[0], raw_img, width, height, p);
    for(unsigned i = 1; i < s; ++i) {
        GauBlur::raw2GauBlur(pyr_img[i], pyr_img[i-1], 
            width, height, p*= pow(2.0, 1.0/s));
    }
    // 合成圖片
    v_uch big_img(width*height*s);
    size_t idx=0;
    for(unsigned j = 0; j < height; ++j) {
        for(unsigned n = 0; n < s; ++n) {
            for(unsigned i = 0; i < width; ++i) {
                big_img[idx++] = pyr_img[n][j*width + i];
            }
        }
    }
    Raw::raw2bmp("big.bmp", 
            big_img, width*s, height, 8);
    system("big.bmp");
}