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

// 高斯模糊
void GauBlur::raw2GauBlur(vector<types>& img_gau,
    vector<types>& img_ori,
    size_t width, size_t height, float p)
{
    // 來源相同例外
    if (&img_gau == &img_ori) {
        throw file_same("## Erroe! in and out is same.");
    }
    // 設定正確的大小
    img_gau.resize(img_ori.size());
    // 緩存
    vector<types> img_gauX(img_ori.size());
    // 高斯矩陣與半徑
    vector<types> gau_mat = gau_matrix(p);
    const size_t r = gau_mat.size() / 2;
    // 高斯模糊 X 軸
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
                sum += img_gauX[i*height + idx] * gau_mat[k];
            }
            img_gau[i*height + j] = sum;
        }
    }
}
// 高斯公式
float GauBlur::gau_meth(size_t r, float p) {
    float two = 2.0;
    float num = exp(-pow(r, two) / (two*pow(p, two)));
    num /= sqrt(two*M_PI)*p;
    return num;
}
// 高斯矩陣
vector<float> GauBlur::gau_matrix(float p) {
    vector<float> gau_mat;
    // 計算矩陣長度
    int mat_len = (int)(((p - 0.8) / 0.3 + 1.0) * 2.0);
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
//----------------------------------------------------------------
bool Corner::harris(const vector<float>& p,
    size_t w, size_t y, size_t x)
{
    // 閥值
    constexpr float r = 10;
    constexpr float thre = ((r + 1)*(r + 1)) / r;
    // 二維讀取
    auto at2d = [&](int y, int x) {return p[y*w + x]; };
    // 公式
    float Dxx = 2 * at2d(y, x) - at2d(y, x - 1) - at2d(y, x + 1);
    float Dyy = 2 * at2d(y, x) - at2d(y - 1, x) - at2d(y + 1, x);
    float Dxy = at2d(y + 1, x + 1) + at2d(y - 1, x - 1)
        - at2d(y - 1, x + 1) - at2d(y + 1, x - 1);
    Dxy /= 4;
    float Tr = Dxx + Dyy;
    float Det = Dxx * Dyy - Dxy*Dxy;
    // 判斷閥值
    if ((Tr*Tr / Det) < thre) {
        return 1;
    } return 0;

}