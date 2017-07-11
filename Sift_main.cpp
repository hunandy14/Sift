/*****************************************************************
Name : 
Date : 2017/07/04
By   : CharlotteHonG
Final: 2017/07/04
*****************************************************************/
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
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
    Raw::read_raw(raw_img, rawName);

    // 特徵點
    Sift img(raw_img, 256, 256);
    img.pyramid();
    return 0;
}
//================================================================
