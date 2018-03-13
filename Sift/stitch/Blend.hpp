#pragma once

// Blend
struct Blend_Image
{
	int width;
	int height;
	float *RGB;
};

void multiBandBlend(Raw &limg, Raw &rimg, int dx, int dy);