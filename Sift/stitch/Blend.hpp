//***********************************************************//
//Company       :   KUAS                                     //
//***********************************************************//
//File          :   Blend.hpp                                //
//Author        :   YAN,RUI-YING                             //
//Language      :   C++                                      //
//IDE           :   Microsoft Visual Studio 2017             //
//OS            :   Windows 10                               //
//Date          :   2017/05/16                               //
//***********************************************************//
#pragma once
#include "imagedata.hpp"
#define RIGHT 1
#define LEFT 0

extern std::vector<unsigned char> MultiBandBlending(std::vector<unsigned char> left, std::vector<unsigned char> right, int width, int height);
extern void multiBandBlend(Raw &limg, Raw &rimg, int dx, int dy);

extern void multiBandBlend_G(unsigned char *limg, unsigned char *rimg, int Xsize, int Ysize, int dx, int dy);