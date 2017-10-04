/*****************************************************************
Name :
Date : 2017/07/04
By   : CharlotteHonG
Final: 2017/07/04
*****************************************************************/
#pragma warning(disable : 4819)
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

#include "imglib.hpp"
constexpr auto M_PI = 3.14159265358979323846;

// 高斯模糊(3x3)
void Gaus::GauBlur(vector<types>& img_gau,
    const vector<types>& img_ori,
    size_t width, size_t height, float p)
{
    // 來源不可相同
    if (&img_gau == &img_ori) {
        throw file_same("## Erroe! in and out is same.");
    }

	/*
    constexpr size_t GauMat_R = 3;
    // 設定正確的大小
    img_gau.resize(img_ori.size());
    // 高斯矩陣
    vector<types> gau_mat = gau_matrix(p, GauMat_R);
    // 高斯橫向
	vector<types> gau(width*height);
    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width; i++) {
            float mat[GauMat_R];
            if (i == 0) {
                mat[0] = img_ori[j*width + i+0] * gau_mat[0];
                mat[1] = img_ori[j*width + i+0] * gau_mat[1];
                mat[2] = img_ori[j*width + i+1] * gau_mat[2];
            } else if (i == width-1) {
                mat[0] = img_ori[j*width + i-1] * gau_mat[0];
                mat[1] = img_ori[j*width + i+0] * gau_mat[1];
                mat[2] = img_ori[j*width + i+0] * gau_mat[2];
            } else {
                mat[0] = img_ori[j*width + i-1] * gau_mat[0];
                mat[1] = img_ori[j*width + i+0] * gau_mat[1];
                mat[2] = img_ori[j*width + i+1] * gau_mat[2];
            }
			gau[j*width + i] = (mat[0] + mat[1] + mat[2]);
        }
    }
    // 高斯縱向
    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width; i++) {
            float mat[GauMat_R];
            if (j == 0) {
                mat[0] = gau[(j+0)*width + i] * gau_mat[0];
				mat[1] = gau[(j+0)*width + i] * gau_mat[1];
				mat[2] = gau[(j+1)*width + i] * gau_mat[2];
            } else if (j == height-1) {
                mat[0] = gau[(j-1)*width + i] * gau_mat[0];
				mat[1] = gau[(j+0)*width + i] * gau_mat[1];
				mat[2] = gau[(j+0)*width + i] * gau_mat[2];
            } else {
                mat[0] = gau[(j-1)*width + i] * gau_mat[0];
				mat[1] = gau[(j+0)*width + i] * gau_mat[1];
				mat[2] = gau[(j+1)*width + i] * gau_mat[2];
            }
            img_gau[j*width+i] = (mat[0] + mat[1] + mat[2]);
        }
    }*/
	
	
	vector<types> gau_mat = gau_matrix(p, 0);
	img_gau.resize(height*width);
    // 緩存
    vector<types> img_gauX(img_ori.size());
    // 高斯模糊 X 軸
    const size_t r = gau_mat.size() / 2;
    for (unsigned j = 0; j < height; ++j) {
        for (unsigned i = 0; i < width; ++i) {
            double sum = 0;
            for (unsigned k = 0; k < gau_mat.size(); ++k) {
                int idx = i-r + k;
                // idx超出邊緣處理
                if (idx < 0) {
                  idx = 0;
                } else if (idx >(int)(width-1)) {
                  idx = (width-1);
                }
                sum += img_ori[j*width + idx] * gau_mat[k];
            }
            img_gauX[j*width + i] = sum;
        }
    }
    // 高斯模糊 Y 軸
    for (unsigned j = 0; j < height; ++j) {
        for (unsigned i = 0; i < width; ++i) {
            double sum = 0;
            for (unsigned k = 0; k < gau_mat.size(); ++k) {
                int idx = j-r + k;
                // idx超出邊緣處理
                if (idx < 0) {
                  idx = 0;
                } else if (idx > (int)(height-1)) {
                  idx = (height-1);
                }
                sum += img_gauX[idx*width + i] * gau_mat[k];
            }
            img_gau[j*width + i] = sum;
        }
    }
}
// 高斯公式
float Gaus::gau_meth(size_t r, float p) {
    float two = 2.0;
    float num = exp(-pow(r, two) / (two*pow(p, two)));
    num /= sqrt(two*M_PI)*p;
    return num;
}
// 高斯差分
void Gaus::GauDog(vector<types>& img_dog,
	vector<types>& img_gau, size_t width, size_t height)
{
	
}

// 高斯矩陣 (mat_len defa=3)
vector<float> Gaus::gau_matrix(float p, size_t mat_len) {
    vector<float> gau_mat;
    // 計算矩陣長度
    if (mat_len == 0) {
        mat_len = (int)(((p - 0.8) / 0.3 + 1.0) * 2.0);// (顏瑞穎給的公式)
    }
    // 奇數修正
    if (mat_len % 2 == 0) { ++mat_len; }
    // 一維高斯矩陣
    gau_mat.resize(mat_len);
    float sum = 0;
    for (int i = 0, j = mat_len / 2; j < mat_len; ++i, ++j) {
        float temp;
        if (i) {
            temp = gau_meth(i, p);
            gau_mat[j] = temp;
            gau_mat[mat_len - j - 1] = temp;
            sum += temp += temp;
        }
        else {
            gau_mat[j] = gau_meth(i, p);
            sum += gau_mat[j];
        }
    }
    // 歸一化
    for (auto&& i : gau_mat) { i /= sum; }
    return gau_mat;
}
vector<Gaus::types> Gaus::gau_matrix2d(vector<types>& gau_mat2d, types p, size_t mat_len) {
    // 高斯2d矩陣
    gau_mat2d.resize(mat_len*mat_len);
    // 二維讀取
    auto at2d = [&](int y, int x)->float& {
        return gau_mat2d[y*mat_len + x];
    };

    vector<float> gau_mat1d = Gaus::gau_matrix(p, mat_len);
    // 做 X
    for (size_t j = 0; j < mat_len; j++) {
        for (size_t i = 0; i < mat_len; i++) {
            at2d(j, i) = gau_mat1d[i];
        }
    }
    // 做 Y
    for (size_t j = 0; j < mat_len; j++) {
        for (size_t i = 0; i < mat_len; i++) {
            at2d(i, j) *= gau_mat1d[i];
        }
    }
    // 檢查
/*
float sum=0;
for (size_t j = 0; j < mat_len; j++)
    for (size_t i = 0; i < mat_len; i++)
        sum+=at2d(i, j);
if (sum == 1) cout << "gau_matrix2d -> check ok" << endl;
*/
    return gau_mat2d;
}

// 正規化
void Gaus::regularization(vector<types>& img, vector<unsigned char>& img_ori) {
	img.resize(img_ori.size());
	for (size_t i = 0; i < img_ori.size(); i++) {
		img[i] = img_ori[i]/255.0;
	}
}
void Gaus::unregularization(vector<unsigned char>& img, vector<types>& img_ori) {
	img.resize(img_ori.size());
	for (size_t i = 0; i < img_ori.size(); i++) {
		img[i] = img_ori[i]*255.0;
	}
}
//----------------------------------------------------------------
// ZroOrder調整大小
void Scaling::zero(vector<types>& img,
    vector<types>& img_ori, size_t width,
    size_t height, float Ratio)
{
    int w = (int)floor(width * Ratio);
    int h = (int)floor(height * Ratio);
    img.resize(w*h);
    for (int j = 0; j < h; ++j) {
        for (int i = 0; i < w; ++i) {
            img[j*w + i] =
                img_ori[(int)(j/Ratio)*width + (int)(i/Ratio)];
        }
    }
}
// FisrtOrder調整大小
void Scaling::first(vector<types>& img,
    vector<types>& img_ori, size_t width,
    size_t height, float Ratio)
{
    int w = (int)floor(width * Ratio);
    int h = (int)floor(height * Ratio);
    img.resize(h*w);
    for (int j = 0; j < h; ++j) {
        for (int i = 0; i < w; ++i) {
            // 對應到原圖的座標
            int oy = (int)floor(j / Ratio);
            int ox = (int)floor(i / Ratio);
            // 附近的四個點
            size_t xp = (ox+1) > (int)(width-1)? width-1: (ox+1);
            size_t yp = (oy+1) > (int)(height-1)? height-1: (oy+1);
            types A = img_ori[oy*width + ox];
            types B = img_ori[oy*width + xp];
            types C = img_ori[yp*width + ox];
            types D = img_ori[yp*width + xp];
            // 公式的 a 與 b
            float a = (i - ox*Ratio) / (Ratio);
            float b = (j - oy*Ratio) / (Ratio);
            types AB = (A*(1.0 - a)) + (B*a);
            types CD = (C*(1.0 - a)) + (D*a);
            types X = (AB*(1.0 - b)) + (CD*b);
            img[j*w + i] = X;
        }
    }
}
// Bicubic調整大小
void Scaling::cubic(vector<types>& img,
    vector<types>& img_ori, size_t width,
    size_t height, float Ratio)
{
    /*
    using uch = types;
    // Bicubic 取得周圍16點
    auto getMask = [&](uch* mask, size_t oy, size_t ox) {
        // 取得周圍16點
        int foy, fox; // 修復後的原始座標
        for (int j = 0, idx = 0; j < 4; ++j) {
            for (int i = 0; i < 4; ++i, ++idx) {
                foy = (int)(oy+(j-1)), fox = (int)(ox+(i-1));
                // 超過左邊界修復
                if (foy<0) { foy = 1; }
                // 超過上邊界修復
                if (fox<0) { fox = 1; }
                // 超過下邊界修復
                if (foy == (int)height) { foy -= 2; }
                if (foy == (int)height - 1) { foy -= 1; }
                // 超過右邊界修復
                if (fox == (int)width) { fox -= 2; }
                if (fox == (int)width - 1) { fox -= 1; }
                // 紀錄對應的指標
                mask[idx] = img_ori[foy*width + fox];
            }
        }
    };

    int w = (int)floor(width * Ratio);
    int h = (int)floor(height * Ratio);
    img.resize(h*w);
    for (int j = 0; j < h; ++j) {
        for (int i = 0; i < w; ++i) {
            // 對應到原圖的座標
            int oy = (int)floor(j / Ratio);
            int ox = (int)floor(i / Ratio);
            double a = (i - ox*Ratio) / (Ratio);
            double b = (j - oy*Ratio) / (Ratio);
            // 取得周圍16點
            uch mask[16];
            // getMask(mask, oy, ox);
            // 導入周圍16點與插入的比例位置
            // uch X = bicubicInterpolate(mask, b, a);
            // cout << "X=" << X << endl;
            // 寫入暫存內
            img[j*w + i] = 100;
        }
    }
    */
}
// Bicubic 插值核心運算
float Scaling::cubicInterpolate (
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
float Scaling::bicubicInterpolate (
    float* p, float y, float x)
{
    float arr[4];
    for (int i = 0; i < 4; ++i){
        arr[i] = cubicInterpolate((i*4 + p), x);
    } return cubicInterpolate(arr, y);
}
//----------------------------------------------------------------
bool Corner::harris(const vector<float>& p,
    size_t w, size_t y, size_t x, float r)
{
    // 閥值
    float thre = ((r + 1)*(r + 1)) / r;
    // 二維讀取
    auto at2d = [&](int y, int x)->float {
        return p[y*w + x];
    };
    // 公式
    float Dxx = 2.f*at2d(y, x) - at2d(y, x - 1) - at2d(y, x + 1);
    float Dyy = 2.f*at2d(y, x) - at2d(y - 1, x) - at2d(y + 1, x);
    float Dxy = at2d(y + 1, x + 1) + at2d(y - 1, x - 1)
        - at2d(y - 1, x + 1) - at2d(y + 1, x - 1);
    Dxy /= 4.f;
    float Tr = Dxx + Dyy;
    float Det = Dxx*Dyy - Dxy*Dxy;
    // 判斷閥值
    float val = (Tr*Tr / Det);
    if (val < thre) {
        return 1;
    }
    // 不成立則刪除
    return 0;
}




