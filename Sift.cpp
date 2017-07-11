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
    s += 3;
    size_t w = raw_img.width, h = raw_img.height;
    // size_t octvs = log(min(w, h)) / log(2.0) - 2;
    
    // 測試層數先設3
    size_t octvs = 3;

    auto dog_gau = [s](ImgRaw img, string name="") {
        size_t w = img.width, h = img.height;
        // 初始化
        vector<ImgRaw> pyrs;
        for(unsigned i = 0; i < s; ++i){
            pyrs.emplace_back(v_uch(), w, h);
        }
        // 高斯模糊
        float p=1.6;
        ImgRaw::gauBlur(pyrs[0], img, p);
        for(unsigned i = 1; i < s; ++i) {
            ImgRaw::gauBlur(pyrs[i], pyrs[i-1], p*= pow(2.0, 1.0/s));
        }
        // 高斯差分
        // for(unsigned j = 0; j < s-1; ++j) {
        //     pyrs[j] = pyrs[j+1] - pyrs[j];
        // }
        // 合成圖片
        size_t s2 = s-1;
        ImgRaw big_img(w*s2, h);
        for(unsigned j=0, idx=0; j < h; ++j) {
            for(unsigned n = 0; n < s2; ++n) {
                for(unsigned i = 0; i < w; ++i) {
                    big_img[idx++] = pyrs[n][j*w + i];
                }
            }
        }
        // 輸出圖片並開啟
        if(name!="") {
            name += ".bmp";
            big_img.bmp(name);
            system(name.c_str());
        }
        return pyrs;
    };
    // 輸入圖
    ImgRaw temp(raw_img.width, raw_img.height);
    ImgRaw::cubic(temp, raw_img, 2);
    // 高斯金字塔
    vector<ImgRaw> pyrs;
    pyrs = dog_gau(temp, "0");
    for(unsigned i = 0; i < octvs-1; ++i) {
        ImgRaw::first(temp, pyrs[s/2], 0.5);
        pyrs = dog_gau(temp, to_string(i+1));
    }

}