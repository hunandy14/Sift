/*****************************************************************
Name :
Date : 2017/06/12
By   : CharlotteHonG
Final: 2017/06/14
*****************************************************************/
#pragma warning(disable : 4819)
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
    static FileHeader makeFH(
        uint32_t width, uint32_t height, uint16_t bits)
    {
        FileHeader file_h;
        file_h.size = file_h.headSize + width*height * bits/8;
        if(bits==8) {
            file_h.size += 1024;
            file_h.headSize += 1024;
        }
        return file_h;
    }
    static InfoHeader makeIH(
        uint32_t width, uint32_t height, uint16_t bits)
    {
        InfoHeader info_h;
        info_h.width = width;
        info_h.height = height;
        info_h.bits = bits;
        info_h.imagesize = width*height * bits/8;
        if(bits==8) {
            info_h.ncolours=256;
        }
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
    static void read_bmp(std::vector<uch>& raw, std::string name);
    // 讀檔
    static void read_raw(std::vector<uch>& raw, std::string name) {
        std::fstream file(name.c_str(),
                          std::ios::in | std::ios::binary | std::ios::ate);
        raw.resize(static_cast<size_t>(file.tellg()));
        file.seekg(0, std::ios::beg);
        file.read(reinterpret_cast<char*>(raw.data()), raw.size());
        file.close();
    }
    // 寫檔
    static void raw2bmp(std::string name, std::vector<uch>& raw,
                        uint32_t width, uint32_t height, uint16_t bits=24);
    // 寫檔
    static void write_raw(std::string name, std::vector<uch>& raw) {
        std::fstream img(name.c_str(), std::ios::out | std::ios::binary);
        img.write(reinterpret_cast<char*>(raw.data()), raw.size());
    }
};
