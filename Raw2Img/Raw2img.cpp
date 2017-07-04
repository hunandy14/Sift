/*****************************************************************
Name :
Date : 2017/06/14
By   : CharlotteHonG
Final: 2017/06/14
*****************************************************************/
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Raw2img.hpp"

using namespace std;
using uch = unsigned char;
//----------------------------------------------------------------
// 寫檔
void Raw::raw2bmp(
    string name, vector<uch>& raw,
    uint32_t width, uint32_t height, uint16_t bits)
{
    // 檔案資訊
    FileHeader file_h = makeFH(width, height, bits);
    // 圖片資訊
    InfoHeader info_h = makeIH(width, height, bits);
    // 寫檔
    fstream img(name, ios::out | ios::binary);
    img << file_h << info_h;
    // 寫調色盤
    if(bits == 8) {
        for(unsigned i = 0; i < 256; ++i) {
            img << uch(i) << uch(i) << uch(i) << uch(0);
        }
    }
    // 寫 Raw 檔
    size_t alig = (4 - width%4)%4;
    for(int j = height-1; j >= 0; --j) {
        for(unsigned i = 0; i < width; ++i) {
            if(bits==24) {
                img << raw[(j*width+i)*3 + B];
                img << raw[(j*width+i)*3 + G];
                img << raw[(j*width+i)*3 + R];
            } else if(bits==8) {
                img << raw[(j*width+i)];
            }
        }
        // 對齊4byte
        for(unsigned i = 0; i < alig; ++i)
            img << uch(0);
    }
    img.close();
}
//----------------------------------------------------------------
// 讀 Bmp 檔案
void Raw::read_bmp(vector<uch>& raw, string name) {
    fstream file(name.c_str(), ios::in | ios::binary);
    file.seekg(0, ios::beg);
    // 讀檔頭
    FileHeader file_h;
    file >> file_h;
    InfoHeader info_h;
    file >> info_h;
    // 讀 Raw
    file.seekg(file_h.headSize, ios::beg);
    raw.resize(info_h.imagesize);
    size_t alig = (4 - info_h.width%4)%4;
    char* p = reinterpret_cast<char*>(raw.data());
    for(int j = info_h.height-1; j >= 0; --j) {
        for(unsigned i = 0; i < info_h.width; ++i) {
            // 來源是 rgb
            if(info_h.bits == 24) {
                file.read(p + j*info_h.width*3+i*3 + G, 1);
                file.read(p + j*info_h.width*3+i*3 + B, 1);
                file.read(p + j*info_h.width*3+i*3 + R, 1);
            }
            // 來源是 gray
            else if(info_h.bits == 8) {
                file.read(p + j*info_h.width+i, 1);
            }
        }
        file.seekg(alig, ios::cur); // 跳開 4bite 對齊的空格
    }
    file.close();
}
