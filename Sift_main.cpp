/*****************************************************************
Name : 
Date : 2017/07/04
By   : CharlotteHonG
Final: 2017/07/04
*****************************************************************/
#include <iostream>
using namespace std;
#include "Sift.hpp"
// #include "OpenRAW_fun\OpenRAW.hpp"

string rawName = "Rena_256x256_8bit.raw";
string rawName2 = "Seymour_Park_960x540_24bit.raw";
string bmpName = "Rena.bmp";
string bmpName2 = "Rena2.bmp";
//================================================================
int main(int argc, char const *argv[]){
    // 讀取圖片
    vector<unsigned char> raw_img;
    Raw::read_raw(raw_img, rawName2);
    Raw::raw2gray(raw_img);
    // 高斯模糊
    vector<unsigned char> img_gau;
    vector<unsigned char> img1;
    float r=2;
    GauBlur::raw2GauBlur(img_gau, raw_img, 960, 540);
    // 縮放大小
    Scaling::cubic(img1, img_gau, 960, 540, r);
    // 輸出BMP
    Raw::raw2bmp(bmpName, img1, 960*r, 540*r, 8);
    // Raw::raw2bmp(bmpName2, raw_img, 960, 540, 8);
    Raw::write_raw("a.raw", img1);
    system(bmpName.c_str());
    // system(bmpName2.c_str());
    return 0;
}
//================================================================
