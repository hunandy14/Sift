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

    // 高斯金字塔
    vector<v_uch> pyr_img;
    pyr_img.resize(octvs*s);

    
    auto fog_gau = [&](v_uch raw_img) {
        // 模糊
        float p=1.6;
        GauBlur::raw2GauBlur(pyr_img[0], raw_img, width, height, p);
        for(unsigned i = 1; i < s; ++i) {
            GauBlur::raw2GauBlur(pyr_img[i], pyr_img[i-1], 
                width, height, p*= pow(2.0, 1.0/s));
        }
        // 高斯差分
        for(unsigned j = 0; j < s-1; ++j) {
            for (unsigned i = 0; i < width*height; ++i){
                pyr_img[j][i] = pyr_img[j+1][i]-pyr_img[j][i] + uch(128);
            }
        }
        // 合成圖片
        size_t s2 = s-1;
        v_uch big_img(width*height * s2);
        for(unsigned j=0, idx=0; j < height; ++j) {
            for(unsigned n = 0; n < s2; ++n) {
                for(unsigned i = 0; i < width; ++i) {
                    big_img[idx++] = pyr_img[n][j*width + i];
                }
            }
        }
        // 輸出圖片並開啟
        Raw::raw2bmp("big.bmp", 
                big_img, width*s2, height, 8);
        system("big.bmp");
    };


    // 放大
    v_uch temp;

    Scaling::cubic(temp, raw_img, width, height, 2);
    height*=2, width*=2;
    // fog_gau(temp);

    // Scaling::first(temp, pyr_img[s/2], width, height, 0.5);
    // height*=0.5, width*=0.5;
    // fog_gau(temp);

    // Scaling::first(temp, pyr_img[s/2], width, height, 0.5);
    // height*=0.5, width*=0.5;
    // fog_gau(temp);
    
    auto dog_gau2 = [&](ImgRaw img) {
        // 初始化
        vector<ImgRaw> pyrs;
        for(unsigned i = 0; i < s; ++i){
            pyrs.emplace_back(v_uch(), img.width, img.height);
        }
        // 高斯模糊
        float p=1.6;
        ImgRaw::gauBlur(pyrs[0], img, p);
        for(unsigned i = 1; i < s; ++i) {
            ImgRaw::gauBlur(pyrs[i], pyrs[i-1], p*= pow(2.0, 1.0/s));
        }
        // 高斯差分
        for(unsigned j = 0; j < s-1; ++j) {
            pyrs[j] = pyrs[j+1] - pyrs[j];
        }
        // 合成圖片
        size_t s2 = s-1;
        v_uch big_img(width*height * s2);
        for(unsigned j=0, idx=0; j < img.height; ++j) {
            for(unsigned n = 0; n < s2; ++n) {
                for(unsigned i = 0; i < img.width; ++i) {
                    big_img[idx++] = pyrs[n][j*width + i];
                }
            }
        }
        // 輸出圖片並開啟
        Raw::raw2bmp("big.bmp", 
                big_img, width*s2, height, 8);
        system("big.bmp");
        return pyrs;
    };
    
    // 原圖
    ImgRaw sou(temp, width, height);
    auto pyrs = dog_gau2(sou);

}