/*****************************************************************
Name : 
Date : 2017/07/04
By   : CharlotteHonG
Final: 2017/07/04
*****************************************************************/
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#include "gau_blur.hpp"
constexpr auto M_PI = 3.14159265358979323846;

// 高斯模糊
void GauBlur::raw2GauBlur(vector<unsigned char>& img_gau, 
    vector<unsigned char>& img_ori, 
    size_t width, size_t height, float p)
{
    // 來源相同例外
    if(&img_gau == &img_ori) {
        throw file_same("## Erroe! in and out is same.");
    }
    // 設定正確的大小
    vector<unsigned char> img_gauX(img_ori.size());
    img_gau.resize(img_ori.size());
    // 高斯矩陣與半徑
    vector<float> gau_mat = gau_matrix(p);
    const int r = gau_mat.size()/2;
    // 高斯模糊 X 軸
    for(unsigned j = 0; j < height; ++j) {
        for(unsigned i = 0; i < width; ++i) {
            size_t sum = 0;
            for(unsigned k = 0; k < gau_mat.size(); ++k) {
                int idx = (i-r+k);
                if(idx < 0) { idx=0; }
                sum += img_ori[j*width+idx]*gau_mat[k];
            }
            img_gauX[j*width+i] = sum;
        }
    }
    // 高斯模糊 Y 軸
    for(unsigned j = 0; j < height; ++j) {
        for(unsigned i = 0; i < width; ++i) {
            size_t sum = 0;
            for(unsigned k = 0; k < gau_mat.size(); ++k) {
                int idx = (j-r+k);
                if(idx < 0) { idx=0; }
                sum += img_gauX[i*height+idx]*gau_mat[k];
            }
            img_gau[i*height+j] = sum;
        }
    }
}

//----------------------------------------------------------------
// 高斯公式
float GauBlur::gau_meth(size_t r, float p) {
    float two = 2;
    float num = exp(-pow(r, two) / (two*pow(p, two)));
    num /= sqrt(two*M_PI)*p;
    return num;
}
// 高斯矩陣
vector<float> GauBlur::gau_matrix(float p){
    vector<float> gau_mat;
    // 計算矩陣長度
    int mat_len = ((p-0.8) / 0.3+1.0) * 2.0;
    if (mat_len % 2 == 0) {++mat_len;}
    // 一維高斯矩陣
    gau_mat.resize(mat_len);
    float sum=0;
    for(int i=0, j=mat_len/2; j < mat_len; ++i, ++j) {
        float temp;
        if(i) {
            temp = gau_meth(i);
            gau_mat[j] = temp;
            gau_mat[mat_len-j-1] = temp;
            sum += temp += temp;
        } else {
            gau_mat[j]=gau_meth(i);
            sum += gau_mat[j];
        }
    }
    // 歸一化
    for(auto&& i : gau_mat) { i /= sum; }
    return gau_mat;
}
//----------------------------------------------------------------
