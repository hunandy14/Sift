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
class Gaus{
private:
    using types = float;
public:
    static void GauBlur(vector<types>& img_gau,
		const vector<types>& img_ori,
        size_t width, size_t height, float p);
	static void GauDog(vector<types>& img_dog,
		vector<types>& img_gau, size_t width, size_t height);
    static vector<types> gau_matrix(float p, size_t mat_len = 3);
	static vector<types> gau_matrix2d(vector<types>& gau_mat2d, types p, size_t mat_len=0);
public:
	static void regularization(vector<types>& img, vector<unsigned char>& img_ori);
	static void unregularization(vector<unsigned char>& img, vector<types>& img_ori);
private:
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
        float* p, float x);
    // Bicubic 輸入16點與插入位置，取得目標值
    static float bicubicInterpolate (
        float* p, float y, float x);
};
// 角點偵測
class Corner {
public:
	static bool harris(const vector<float>& p,
		size_t w, size_t y, size_t x, float r);
};
