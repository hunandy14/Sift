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
    // size_t octvs = log(std::min(height, width)) / log(2.0) - 2;
    size_t octvs = 3;
    vector<v_uch> pyr_img;
    pyr_img.resize(octvs*s);
    // 放大
    v_uch temp;
    Scaling::cubic(temp, raw_img, height, width, 2);
    height*=2, width*=2;
    raw_img=temp;
    // 模糊
    float p=1.6;
    GauBlur::raw2GauBlur(pyr_img[0], raw_img, height, width, p);
    for(unsigned i = 1; i < s; ++i) {
        p *= pow(2.0, 1.0/s);
        GauBlur::raw2GauBlur(pyr_img[i], pyr_img[i-1], height, width, p);
    }
    

    // Raw::raw2bmp("2.bmp", temp, width*2, height*2, 8);

    // vector<unsigned char> img1;
    float r=2;
    // Scaling::cubic(img1, raw_img, 256, 256, r);
    Raw::raw2bmp("0.bmp", pyr_img[0], 256*r, 256*r, 8);
    Raw::raw2bmp("1.bmp", pyr_img[1], 256*r, 256*r, 8);
    Raw::raw2bmp("2.bmp", pyr_img[2], 256*r, 256*r, 8);
    Raw::raw2bmp("3.bmp", pyr_img[3], 256*r, 256*r, 8);
    Raw::raw2bmp("4.bmp", pyr_img[4], 256*r, 256*r, 8);
    Raw::raw2bmp("5.bmp", pyr_img[5], 256*r, 256*r, 8);
}