#pragma once

// Blend
struct Blend_Image
{
	int width;
	int height;
	float *RGB;
};

void multiBandBlend(Raw &limg, Raw &rimg, int dx, int dy);

float getFocal(const vector<float>& HomogMat, size_t img1Size, size_t img2Size);

void blen2img(const ImgRaw & img1, const ImgRaw & img2, ImgRaw & dst, const vector<float>& HomogMat, const Feature * const * RANSAC_feat, int RANSAC_num);
