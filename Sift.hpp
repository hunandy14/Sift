/*****************************************************************
Name : 
Date : 2017/07/05
By   : CharlotteHonG
Final: 2017/07/05
*****************************************************************/
#pragma once

#include "imglib\imglib.hpp"
#include "Raw2Img\Raw2Img.hpp"

class Sift{
public:
    Sift(vector<unsigned char> raw_img, size_t width, size_t height): 
        raw_img(raw_img), width(width), height(height)
    {

    }
public:
    void pyramid(size_t s=3);
private:
    using uch = unsigned char;
    using v_uch = vector<unsigned char>;
    v_uch raw_img;
    size_t width;
    size_t height;
};