/*****************************************************************
Name :
Date : 2017/07/04
By   : CharlotteHonG
Final: 2017/07/04
*****************************************************************/
#pragma warning(disable : 4819)
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
private:
    using types = float;
public:
    static void raw2GauBlur(vector<types>& img_gau,
        vector<types>& img_ori,
        size_t width, size_t height, float p);
private:
    static vector<types> gau_matrix(float p);
    static types gau_meth(size_t r, float p);
};
// 圖像縮放
class Scaling{
private:
    using types = float;
public:
    // ZroOrder調整大小
    static void zero(vector<types>& img,
        vector<types>& img_ori, 
        size_t width, size_t height, float Ratio);
    // FisrtOrder調整大小
    static void first(vector<types>& img,
        vector<types>& img_ori, 
        size_t width, size_t height, float Ratio);
    // Bicubic調整大小
    static void cubic(vector<types>& img,
        vector<types>& img_ori, size_t width, 
        size_t height, float Ratio);
private:
    // Bicubic 插值核心運算
    static float cubicInterpolate (
        float* p, float x)
    {
        float temp = (float)(p[1] + 0.5 * 
            x*(p[2] - p[0] +x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - 
                p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0]))));
        if (temp > 255) { temp = 255; }
        else if (temp < 0) { temp = 0; }
        return temp;
    }
    // Bicubic 輸入16點與插入位置，取得目標值
    static float bicubicInterpolate (
        float* p, float y, float x)
    {
        float arr[4];
        for (int i = 0; i < 4; ++i){
            arr[i] = cubicInterpolate((i*4 + p), x);
        } return cubicInterpolate(arr, y);
    }
};

