/*****************************************************************
Name : 
Date : 2017/07/05
By   : CharlotteHonG
Final: 2017/07/05
*****************************************************************/
#pragma warning(disable : 4819)
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include "Sift.hpp"
using namespace std;

void Sift::pyramid(size_t s){
    s += 3;
    size_t octvs = 3;
    octvs = (size_t)(log(min(raw_img.width, raw_img.height)) / log(2.0)-2);
    // 合成圖片
    auto comp = [](vector<ImgRaw>& pyrs, string name=""){
        // 合成圖片
        size_t w = pyrs[0].width, h = pyrs[0].height;
        size_t s2 = pyrs.size();
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
            // system(name.c_str());
        }
    };
    // 高斯模糊
    auto dog_gau = [s](ImgRaw img, string name="") {
        size_t w = img.width, h = img.height;
        // 初始化
        vector<ImgRaw> pyrs;
        for(unsigned i = 0; i < s; ++i){
            pyrs.emplace_back(v_uch(), w, h);
        }
        // 高斯模糊
        float p = float(1.6);
        ImgRaw::gauBlur(pyrs[0], img, p);
        for(unsigned i = 1; i < s; ++i) {
            p *= (float)pow(2.0, 1.0/s);
            ImgRaw::gauBlur(pyrs[i], pyrs[i-1], p);
        }
        return pyrs;
    };
    // 差分圖
    auto dog_diff = [](vector<ImgRaw>& pyrs) {
        // 高斯差分
        for(unsigned j = 0; j < pyrs.size()-1; ++j) {
            pyrs[j] = pyrs[j+1] - pyrs[j];
        } 
		pyrs.erase(--(pyrs.end()));
    };
    // 輸入圖
    ImgRaw temp(raw_img.width, raw_img.height);
    ImgRaw::first(temp, raw_img, 2);
    // 高斯金字塔
    vector<vector<ImgRaw>> pyrs(octvs);
    pyrs[0] = dog_gau(temp, "0");
    comp(pyrs[0], "Sift-gau_"+to_string(0));
    for(unsigned i = 1; i < octvs; ++i) {
        ImgRaw::first(temp, pyrs[i-1][s/2], 0.5);
		pyrs[i] = dog_gau(temp, to_string(i));
		comp(pyrs[i], "Sift-gau_" + to_string(i));
    }
    // 高斯金字塔差分圖
    for(unsigned i = 0; i < octvs; ++i) {
        dog_diff(pyrs[i]);
		comp(pyrs[i], "Sift-diff_" + to_string(i));
    }
}