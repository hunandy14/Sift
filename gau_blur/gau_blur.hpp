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

class GauBlur{
public:
    static void raw2GauBlur(vector<unsigned char>& img_gau, 
        vector<unsigned char>& img_ori, 
        size_t width, size_t height, float p=1.69864);
private:
    static vector<float> gau_matrix(float p=1.69864);
    static float g_onedim(size_t r, float p=1.69864);
};

