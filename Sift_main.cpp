/*****************************************************************
Name : 
Date : 2017/07/04
By   : CharlotteHonG
Final: 2017/07/04
*****************************************************************/
#include <iostream>
using namespace std;
#include "gau_blur\gau_blur.hpp"
#include "Raw2Img\Raw2Img.hpp"

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
    // 高斯暫存
    vector<unsigned char> gau_img;
    GauBlur::raw2GauBlur(gau_img, raw_img, 960, 540);
    // 輸出BMP
    Raw::raw2bmp(bmpName, gau_img, 960, 540, 8);
    Raw::raw2bmp(bmpName2, raw_img, 960, 540, 8);
    system(bmpName.c_str());
    system(bmpName2.c_str());
    return 0;
}
//================================================================
