#pragma once

/*************************** Function Prototypes *****************************/
std::vector<unsigned char> StitchingImg(Raw &img, Raw &xformed, int &W, int &H);
void WarpPerspective(const Raw &src, Raw &dst, const std::vector<float> &H);
void _WarpPerspective(const Raw &src, Raw &dst, const std::vector<float> &H);
Raw DealWithImgData(Raw &srcdata, int width, int height, float R);
void warping(const std::vector<Raw> &inputArrays, float FL2, 
	std::vector<Raw> &Output, std::vector<fpoint> &upedge, std::vector<fpoint> &downedge);
