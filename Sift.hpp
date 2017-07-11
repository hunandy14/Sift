/*****************************************************************
Name : 
Date : 2017/07/05
By   : CharlotteHonG
Final: 2017/07/05
*****************************************************************/
#pragma once

#include "imglib\imglib.hpp"
#include "Raw2Img\Raw2Img.hpp"


class ImgRaw{
private:
    using uch = unsigned char;
public:
    ImgRaw(vector<uch> img, size_t width, size_t height):
        raw_img(img), width(width), height(height){}
    ImgRaw(size_t width, size_t height):
        raw_img(width*height), width(width), height(height){}
    operator vector<uch>&(){ return raw_img; }
    // 重載下標符號
    uch & operator[](size_t idx){
        return const_cast<uch&>(static_cast<const ImgRaw&>(*this)[idx]);
    }
    const uch & operator[](size_t idx) const{
        return raw_img[idx];
    }
public:
    // 大小是否相等
    friend bool operator!=(const ImgRaw& lhs, const ImgRaw& rhs);
    friend bool operator==(const ImgRaw& lhs, const ImgRaw& rhs);
    // 差分圖
    ImgRaw& operator-=(const ImgRaw& rhs){
        if((*this) != rhs) { // 尺寸不同
            cout << "Error size is diff." << endl;
            return (*this);
        }
        for(unsigned i = 0; i < width*height; ++i) {
            raw_img[i] -= rhs.raw_img[i]+128;
        }
        return (*this);
    }
    friend ImgRaw operator-(ImgRaw lhs, const ImgRaw& rhs){
        return lhs -= rhs;
    }
    // 寫 BMP 檔
    ImgRaw& bmp(string name){
        Raw::raw2bmp(name, raw_img, width, height, 8);
        return (*this);
    }

public:
    static void zero(ImgRaw& tar, ImgRaw& sou, float z){
        Scaling::zero(tar, sou, sou.width, sou.height, z);
        tar.width = sou.width*z, tar.height = sou.width*z;
    }
    static void first(ImgRaw& tar, ImgRaw& sou, float z){
        Scaling::first(tar, sou, sou.width, sou.height, z);
        tar.width = sou.width*z, tar.height = sou.width*z;
    }
    static void cubic(ImgRaw& tar, ImgRaw& sou, float z){
        Scaling::cubic(tar, sou, sou.width, sou.height, z);
        tar.width = sou.width*z, tar.height = sou.width*z;
    }
    static void gauBlur(ImgRaw& tar, ImgRaw& sou, float p){
        GauBlur::raw2GauBlur(tar, sou, sou.width, sou.height, p);
    }
public:
    vector<uch> raw_img;
    size_t width;
    size_t height;
};

// 大小是否相等
inline bool operator!=(const ImgRaw& lhs, const ImgRaw& rhs){
    return !(lhs == rhs);
}
inline bool operator==(const ImgRaw& lhs, const ImgRaw& rhs){
    if(lhs.width==rhs.width &&  lhs.height == rhs.height) {
        return 1;
    } return 0;
}
//----------------------------------------------------------------
class Sift{
public:
    Sift(vector<unsigned char> raw_img, size_t width, size_t height): 
        // raw_img(raw_img), width(width), height(height)
        raw_img(raw_img, width, height)
    {

    }
public:
    void pyramid(size_t s=3);
private:
    using uch = unsigned char;
    using v_uch = vector<unsigned char>;
    // v_uch raw_img;
    // size_t width;
    // size_t height;
    ImgRaw raw_img;
};