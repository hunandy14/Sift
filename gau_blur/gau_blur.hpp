/*****************************************************************
Name : 
Date : 2017/07/04
By   : CharlotteHonG
Final: 2017/07/04
*****************************************************************/
#pragma once

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
// 來源相同的例外
class file_same : public std::runtime_error {
public:
    file_same(const std::string& str): std::runtime_error(str) {}
};
// 高斯模糊
class GauBlur{
public:
    static void raw2GauBlur(vector<unsigned char>& img_gau, 
        vector<unsigned char>& img_ori, 
        size_t width, size_t height, float p=1.69864);
private:
    static vector<float> gau_matrix(float p=1.69864);
    static float g_onedim(size_t r, float p=1.69864);
};

