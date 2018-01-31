#pragma once

/*************************** Function Prototypes *****************************/
extern std::vector<unsigned char> StitchingImg(Raw &img, Raw &xformed, int &W, int &H);
extern void WarpPerspective(const Raw &src, Raw &dst, const std::vector<float> &H);
extern void _WarpPerspective(const Raw &src, Raw &dst, const std::vector<float> &H);
extern Raw DealWithImgData(Raw &srcdata, int width, int height, float R);
extern void warping(const std::vector<Raw> &inputArrays, float FL2, 
	std::vector<Raw> &Output, std::vector<fpoint> &upedge, std::vector<fpoint> &downedge);
