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
#include "Sift.hpp"
#include <opencv2/opencv.hpp>
using namespace cv;

//================================================================
int main(int argc, char const *argv[]) {
    // 轉浮點數+正規化
    auto raw2f = [](vector<unsigned char>& img) {
        vector<float> temp(size(img));
        for (unsigned i = 0; i < temp.size(); ++i)
            temp[i] = (float)img[i] / 255.0;
        return temp;
    };

#define testpoint1
#ifdef testpoint1
    // 讀取圖片
    vector<unsigned char> raw_img;
    vector<unsigned char> raw_img2;
    // 圖片1
    Raw::read_bmp(raw_img, "kanna.bmp");
    Raw::raw2gray(raw_img); // 轉灰階
                            // 圖片2
                            //Raw::read_bmp(raw_img, "en.bmp");

                            // 創建結構(寬, 長)
                            //ImgRaw temp(raw2f(raw_img), 420, 420);
    ImgRaw temp(raw2f(raw_img), 850, 602);
    ImgRaw input_img(0, 0);
    ImgRaw::first(input_img, temp, 1);
    // 金字塔
    Sift img(input_img);
    img.pyramid(1);

#endif // testpoint1
    /*
    vector<unsigned char> fea_img;
    Raw::read_bmp(fea_img, "fea\fea0-1.bmp");

    ImgRaw fea_pix(raw2f(fea_img), 850, 602);
    ImgRaw fea_pix_ori(raw2f(fea_img), 850, 602);

    for (size_t j = 1; j < fea_pix.height - 1; j++) {
    for (size_t i = 1; i < fea_pix.width - 1; i++) {
    if (temp.at2d(j, i) == 1) {
    fea_pix.at2d(j, i) = 0.5;
    if (Corner::harri(temp, temp.width, j, i) == 1) {
    fea_pix.at2d(j, i) = 0;
    }
    }
    else {
    //fea_pix.at2d(j, i) = 0.5;
    }
    }
    }
    */
    //fea_pix.bmp("harri.bmp", 8);
















    // #define harri_corner
#ifdef harri_corner
    vector<unsigned char> fea_img;
    Raw::read_bmp(fea_img, "tt.bmp");
    ImgRaw fea_pix_ori(raw2f(fea_img), 420, 200);
    ImgRaw fea_pix(raw2f(fea_img), 420, 200);

    for (size_t j = 1; j < fea_pix_ori.height - 1; j++) {
        for (size_t i = 1; i < fea_pix_ori.width - 1; i++) {
            if (fea_pix_ori.at2d(j, i) == 1) {
                fea_pix_ori.at2d(j, i) = 0;
            }
            else {
                fea_pix_ori.at2d(j, i) = 1;
            }
        }
    }
    for (size_t j = 1; j < fea_pix.height - 1; j++) {
        for (size_t i = 1; i < fea_pix.width - 1; i++) {
            if (fea_pix_ori.at2d(j, i) == 1) {
                fea_pix.at2d(j, i) = 0.5;
                if (Corner::harri(fea_pix_ori, fea_pix_ori.width, j, i) == 1) {
                    fea_pix.at2d(j, i) = 0;
                }
            }
            else {
                //fea_pix.at2d(j, i) = 0.5;
            }
        }
    }
    fea_pix.bmp("harri.bmp", 8);
#endif // harri_corner

    return 0;

}
//================================================================
