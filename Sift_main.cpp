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
string bmpName = "Rena.bmp";
//================================================================
int main(int argc, char const *argv[]){
    // 讀取圖片
    vector<unsigned char> raw_img;
    Raw::read_raw(raw_img, rawName);
    // 高斯暫存
    vector<unsigned char> gau_img;
    GauBlur::raw2GauBlur(gau_img, raw_img, 256, 256);
    // 輸出BMP
    Raw::raw2bmp(bmpName, gau_img, 256, 256, 8);
    system(bmpName.c_str());
    return 0;
}
//================================================================
