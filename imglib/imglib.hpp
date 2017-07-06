/*****************************************************************
Name :
Date : 2017/07/04
By   : CharlotteHonG
Final: 2017/07/04
*****************************************************************/
#pragma once

using std::vector;
using std::string;
// 來源相同的例外
class file_same : public std::runtime_error {
public:
    file_same(const string& str): std::runtime_error(str) {}
};
// 高斯模糊
class GauBlur{
public:
    static void raw2GauBlur(vector<unsigned char>& img_gau,
        vector<unsigned char>& img_ori,
        size_t width, size_t height, float p=1.69864);
private:
    static vector<float> gau_matrix(float p=1.69864);
    static float gau_meth(size_t r, float p=1.69864);
};
// 圖像縮放
class Scaling{
public:
    // ZroOrder調整大小
    static void zero(vector<unsigned char>& img,
        vector<unsigned char>& img_ori, 
        size_t width, size_t height, float Ratio);
    // FisrtOrder調整大小
    static void first(vector<unsigned char>& img,
        vector<unsigned char>& img_ori, 
        size_t width, size_t height, float Ratio);
    // Bicubic調整大小
    static void cubic(vector<unsigned char>& img,
        vector<unsigned char>& img_ori, size_t width, 
        size_t height, float Ratio);
private:
    // Bicubic 插值核心運算
    static unsigned char cubicInterpolate (
        unsigned char* p, double x)
    {
        double temp = (double)(p[1] + 0.5 * 
            x*(p[2] - p[0] +x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - 
                p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0]))));
        if (temp > 255) { temp = 255; }
        else if (temp < 0) { temp = 0; }
        return static_cast<unsigned char>(temp);
    }
    // Bicubic 輸入16點與插入位置，取得目標值
    static unsigned char bicubicInterpolate (
        unsigned char** p, double y, double x)
    {
        unsigned char* arr = new unsigned char[4];
        for (int i = 0; i < 4; ++i)
            arr[i] = cubicInterpolate(p[i], x);
        return cubicInterpolate(arr, y);
    }
};

