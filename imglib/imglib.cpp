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

#include "imglib.hpp"
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

// ZroOrder調整大小
void Scaling::zero(vector<unsigned char>& img,
    vector<unsigned char>& img_ori, size_t width,
    size_t height, float Ratio)
{
    int w=floor(width * Ratio);
    int h=floor(height * Ratio);
    img.resize(w*h);
    for(int j = 0; j < h; ++j) {
        for(int i = 0; i < w; ++i) {
            img[j*w+i] = img_ori[(int)(j/Ratio)*width + (int)(i/Ratio)];
        }
    }
}
// FisrtOrder調整大小
void Scaling::first(vector<unsigned char>& img,
    vector<unsigned char>& img_ori, size_t width,
    size_t height, float Ratio)
{
    int w=floor(width * Ratio);
    int h=floor(height * Ratio);
    img.resize(h*w);

    unsigned char A, B, C, D;// 附近的四個點
    unsigned char AB, CD, X;
    int oy, ox;// 對應到原圖的座標
    int a, b;// 公式的 a 與 b

    for(int j=0; j < h; ++j) {
        for(int i=0; i < w; ++i) {
            oy=j/Ratio; ox=i/Ratio;

            A = img_ori[oy*width + ox];
            B = img_ori[oy*width + ox+1];
            C = img_ori[oy+1*width + ox];
            D = img_ori[oy+1*width + ox+1];

            a = (i-ox*Ratio)/(Ratio);
            b = (j-oy*Ratio)/(Ratio);
            AB = (A*(1.0-a)) + (B*a);
            CD = (C*(1.0-a)) + (D*a);
            X = ((AB*(1.0-b)) + (CD*b));

            img[j*w + i] = X;
        }
    }
}
// Bicubic調整大小
void Scaling::cubic(vector<unsigned char>& img,
    vector<unsigned char>& img_ori, size_t width,
    size_t height, float Ratio)
{
    // Bicubic 取得周圍16點
    auto getMask = [&](size_t oy, size_t ox){
        unsigned char** mask = new unsigned char*[4];
        for (unsigned i = 0; i < 4; ++i)
            mask[i] = new unsigned char[4];
        // 取得周圍16點
        int foy,fox; // 修復後的原始座標
        for (int j = 0; j < 4; ++j){
            for (int i = 0; i < 4; ++i){
                foy=oy+(j-1); fox=ox+(i-1);
                // 超過左邊界修復
                if (foy<0){foy=1;}
                // 超過上邊界修復
                if (fox<0){fox=1;}
                // 超過下邊界修復
                if(foy==(int)height){foy-=2;}
                if(foy==(int)height-1){foy-=1;}
                // 超過右邊界修復
                if (fox==(int)width){fox-=2;}
                if (fox==(int)width-1){fox-=1;}
                // 紀錄對應的指標
                mask[j][i] = img_ori[foy*width + fox];
            }
        }
        // 釋放記憶體
        // for (int i = 0; i < 4; ++i)
        //     delete [] mask[i];
        // delete [] mask;
        return mask;
    };

    int w=floor(width * Ratio);
    int h=floor(height * Ratio);
    img.resize(h*w);
    int oy, ox; // 對應到原圖的座標
    double a, b;// 插入的比例位置
    unsigned char** mask;// 遮罩(周圍的16點)
    unsigned char X;     // 暫存
    for(int j = 0; j < h; ++j) {
        for(int i = 0; i < w; ++i) {
            oy=(int)j/Ratio; ox=(int)i/Ratio;
            a = (i-ox*Ratio)/(Ratio);
            b = (j-oy*Ratio)/(Ratio);
            // 取得周圍16點
            mask = getMask(oy, ox);
            // 導入周圍16點與插入的比例位置
            X = bicubicInterpolate(mask, b, a);
            // 寫入暫存內
            img[j*w + i] = X;
            // 釋放記憶體
            for (int i = 0; i < 4; ++i)
                delete [] mask[i];
            delete [] mask;
        }
    }
}
//----------------------------------------------------------------
