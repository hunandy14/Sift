/*****************************************************************
Name :
Date : 2017/07/05
By   : CharlotteHonG
Final: 2017/07/05
*****************************************************************/
#pragma warning(disable : 4819)
#pragma once

#if defined(_MSC_VER)
    #define or ||
    #define and &&
    #define OR ||
    #define AND &&
#endif

#include "imglib\imglib.hpp"
#include "Raw2Img\Raw2Img.hpp"

// 尺寸大小不合
class Size_error : public std::runtime_error {
public:
    Size_error(const std::string& str): std::runtime_error(str) {}
};
//----------------------------------------------------------------
class ImgRaw {
private:
    using types = float;
public:
    ImgRaw(vector<types> img, size_t width, size_t height) :
        raw_img(img), width(width), height(height) {}
    ImgRaw(size_t width, size_t height) :
        raw_img(width*height), width(width), height(height) {}
    // 隱式轉換
    operator vector<types>&() {
        return raw_img;
    }
    operator vector<unsigned char>() {
        vector<unsigned char> img(raw_img.size());
        for(unsigned i = 0; i < raw_img.size(); ++i) {
            img[i] = (unsigned char)(raw_img[i]*255);
        }
        return img;
    }
    // 重載下標符號
    types & operator[](size_t idx) {
        return const_cast<types&>(static_cast<const ImgRaw&>(*this)[idx]);
    }
    const types & operator[](size_t idx) const {
        return raw_img[idx];
    }
    // 二維讀取
    types & at2d(size_t y, size_t x) {
        return const_cast<types&>(static_cast<const ImgRaw&>(*this).at2d(y, x));
    }
    const types & at2d(size_t y, size_t x) const {
        return raw_img[y*width + x];
    }
public:
    // 大小是否相等
    friend bool operator!=(const ImgRaw& lhs, const ImgRaw& rhs);
    friend bool operator==(const ImgRaw& lhs, const ImgRaw& rhs);
    // 差分圖
    ImgRaw& operator-=(const ImgRaw& rhs) {
        // 尺寸不同
        if ((*this) != rhs) {
            throw Size_error("Error size is diff.");
        }
        for (unsigned i = 0; i < width*height; ++i) {
            raw_img[i] -= rhs.raw_img[i];
        }
        return (*this);
    }
    friend ImgRaw operator-(ImgRaw lhs, const ImgRaw& rhs) {
        return lhs -= rhs;
    }
    // 寫 BMP 檔
    ImgRaw& bmp(string name, size_t bits) {
        vector<unsigned char> img = (*this);
        Raw::raw2bmp(name, img, width, height, bits);
        return (*this);
    }

public:
    static void zero(ImgRaw& tar, ImgRaw& sou, float z) {
        Scaling::zero(tar, sou, sou.width, sou.height, z);
        tar.width = (size_t)(sou.width*z);
        tar.height = (size_t)(sou.height*z);
    }
    static void first(ImgRaw& tar, ImgRaw& sou, float z) {
        Scaling::first(tar, sou, sou.width, sou.height, z);
        tar.width = (size_t)(sou.width*z);
        tar.height = (size_t)(sou.height*z);
    }
    static void cubic(ImgRaw& tar, ImgRaw& sou, float z) {
        Scaling::cubic(tar, sou, sou.width, sou.height, z);
        tar.width = (size_t)(sou.width*z);
        tar.height = (size_t)(sou.height*z);
    }
    static void gauBlur(ImgRaw& tar, ImgRaw& sou, float p) {
        GauBlur::raw2GauBlur(tar, sou, sou.width, sou.height, p);
        tar.width = (size_t)(sou.width);
        tar.height = (size_t)(sou.height);
    }
public:
    vector<types> raw_img;
    size_t width;
    size_t height;
};
// 大小是否相等
inline bool operator!=(const ImgRaw& lhs, const ImgRaw& rhs) {
    return !(lhs == rhs);
}
inline bool operator==(const ImgRaw& lhs, const ImgRaw& rhs) {
    if (lhs.width == rhs.width &&  lhs.height == rhs.height) {
        return 1;
    } return 0;
}
//----------------------------------------------------------------
class Sift {
private:
    using types = float;
public:
    Sift(ImgRaw img): raw_img(img) {}
    Sift(vector<float> raw_img, size_t width, size_t height):
        raw_img(raw_img, width, height) {}
public:
    void pyramid(size_t s = 3);
    void comp(vector<ImgRaw>& pyrs, string name="");
    vector<ImgRaw> dog_gau(ImgRaw& img, size_t s);
	
private:
    ImgRaw raw_img;
    vector<vector<ImgRaw>> pyrs;
};