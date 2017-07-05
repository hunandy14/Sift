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
    static float gau_meth(size_t r, float p=1.69864);
};
// 圖像縮放
class Scaling{
public:
    // ZroOrder調整大小
    static void zero(vector<unsigned char>& img,
        vector<unsigned char>& img_ori, size_t width, size_t height,
        float Ratio)
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
    static void first(vector<unsigned char>& img,
        vector<unsigned char>& img_ori, size_t width, size_t height,
        float Ratio)
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
    // void imgraw::resize_bicubic(float Ratio) {
    //     if(Ratio <= 0) {
    //         cout << "Ratio more than the zero." << endl;
    //         return;
    //     }
    //     int w=floor(this->width * Ratio);
    //     int h=floor(this->high * Ratio);
    //     imgraw img2(h, w);
    //     int oy, ox; // 對應到原圖的座標
    //     double a, b;// 插入的比例位置
    //     unsigned char** mask;// 遮罩(周圍的16點)
    //     unsigned char X;     // 暫存
    //     for(int j = 0; j < h; ++j) {
    //         for(int i = 0; i < w; ++i) {
    //             oy=(int)j/Ratio; ox=(int)i/Ratio;
    //             a = (i-ox*Ratio)/(Ratio);
    //             b = (j-oy*Ratio)/(Ratio);
    //             // 取得周圍16點
    //             mask = this->getMask(oy, ox);
    //             // 導入周圍16點與插入的比例位置
    //             X = bicubicInterpolate(mask, b, a);
    //             // 寫入暫存內
    //             img2.point_write(j, i, X);

    //             // 釋放記憶體
    //             for (int i = 0; i < 4; ++i)
    //                 delete [] mask[i];
    //             delete [] mask;
    //         }
    //     }
    //     // 輸出暫存
    //     *this = img2;
    // }
private:
    /*
    // Bicubic 取得周圍16點
    unsigned char** getMask(int oy, int ox){
        unsigned char** mask = new unsigned char*[4];
        for (int i = 0; i < 4; ++i)
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
                if(foy==this->high){foy-=2;}
                if(foy==this->high-1){foy-=1;}
                // 超過右邊界修復
                if (fox==this->width){fox-=2;}
                if (fox==this->width-1){fox-=1;}
                // 紀錄對應的指標
                mask[j][i] = this->point_read(foy, fox);
            }
        }
        // 釋放記憶體
        // for (int i = 0; i < 4; ++i)
        //     delete [] mask[i];
        // delete [] mask;
        return mask;
    }
    // Bicubic 插值核心運算
    unsigned char cubicInterpolate (unsigned char* p, double x) {
        double temp = (double)(p[1] + 0.5 * 
            x*(p[2] - p[0] +x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - 
                p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0]))));
        if (temp > 255) { temp = 255; }
        else if (temp < 0) { temp = 0; }
        return (unsigned char)temp;
    }
    // Bicubic 輸入16點與插入位置，取得目標值
    unsigned char bicubicInterpolate (unsigned char** p, double y, double x) {
        unsigned char* arr = new unsigned char[4];
        for (int i = 0; i < 4; ++i)
            arr[i] = cubicInterpolate(p[i], x);
        return cubicInterpolate(arr, y);
    }
    */
};

