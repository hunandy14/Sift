/*****************************************************************
Name :
Date : 2017/06/12
By   : CharlotteHonG
Final: 2017/06/14
*****************************************************************/
#pragma once
//----------------------------------------------------------------
#include <fstream>
#include "Raw2img_type.hpp"
/*
     ######
     ##   ##
     ##   ##   ######  ##   ##
     ######   ##   ##  ## # ##
     ## ##    ##   ##  ## # ##
     ##  ##   ##  ###  ## # ##
     ##   ##   ### ##   ## ##

*/
class Raw {
private:
    using uch = unsigned char;
    // RGB 轉灰階公式
    static uch rgb2gray(uch* p) {
        return ((
            19595 * (*(p+R))+
            38469 * (*(p+G))+
            7472  * (*(p+B))) >> 16);
    }
    // 創建檔頭
    static BmpFileHeader makeFH(
        uint32_t width, uint32_t height, uint16_t bits)
    {
        BmpFileHeader file_h;
        file_h.bfSize = file_h.bfOffBits + width*height * bits/8;
        if(bits==8) {file_h.bfSize += 1024, file_h.bfOffBits += 1024;}
        return file_h;
    }
    static BmpInfoHeader makeIH(
        uint32_t width, uint32_t height, uint16_t bits)
    {
        BmpInfoHeader info_h;
        info_h.biWidth = width;
        info_h.biHeight = height;
        info_h.biBitCount = bits;
        info_h.biSizeImage = width*height * bits/8;
        if(bits==8) {info_h.biClrUsed=256;}
        return info_h;
    }
public:
    // 轉灰階
    static void raw2gray(std::vector<uch>& raw) {
        raw2gray(raw, raw);
    }
    static void raw2gray(std::vector<uch>& gray, std::vector<uch>& raw) {
        // 判定相等
        if(&gray!=&raw)
            gray.resize(raw.size()/3);
        // 轉換
        for(unsigned i = 0; i < raw.size()/3; ++i)
            gray[i] = rgb2gray(&raw[i*3]);
        gray.resize(raw.size()/3);
    }
public:
    // 讀 Bmp 檔案
    static void read_bmp(std::vector<uch>& raw, std::string name,
        uint32_t* width=nullptr, uint32_t* height=nullptr, 
        uint16_t* bits=nullptr);
    // 讀 Raw 檔
    static void read_raw(std::vector<uch>& raw, std::string name);
    // 寫 Bmp 檔
    static void raw2bmp(std::string name, std::vector<uch>& raw,
        uint32_t width, uint32_t height, uint16_t bits=24);
    // 寫 Raw 檔
    static void write_raw(std::string name, std::vector<uch>& raw);
};
