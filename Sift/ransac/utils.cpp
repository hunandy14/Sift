#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>

#include "imagedata.hpp"
#include "utils.hpp"
/*************************** Function Definitions ****************************/
int array_double(void** array, int n, int size)
{
	void* tmp;

	tmp = realloc(*array, 2 * n * size);
	if (!tmp)
	{
		std::cout << "Warning: unable to allocate memory in array_double()" << std::endl;
		if (*array)
			free(*array);
		*array = NULL;
		return 0;
	}
	*array = tmp;
	return n * 2;
}
float dist_sq_2D(fpoint p1, fpoint p2)
{
	float x_diff = p1.x - p2.x;
	float y_diff = p1.y - p2.y;

	return x_diff * x_diff + y_diff * y_diff;
}

/*
void KeyPointImg(char *name, Image_Data &src, struct feature **feat, int n)
{
	Image_Data KeyPointRGB(src.getCol(), src.getRow());
	memcpy(&KeyPointRGB, &src, sizeof(src));
	int col = 0, row = 0;
	int i = 0, j = 0;
	for(i = 0; i < n; ++i)
	{
		col = (int)(*feat)[i].img_pt.x;
		row = (int)(*feat)[i].img_pt.y;
		for (i = 0; i < 3; i++)
		{
			KeyPointRGB.RGB[(row * KeyPointRGB.getCol() + (col + i)) * 3] = (unsigned char)255;
			KeyPointRGB.RGB[(row * KeyPointRGB.getCol() + (col + i)) * 3 + 1] = (unsigned char)255;
			KeyPointRGB.RGB[(row * KeyPointRGB.getCol() + (col + i)) * 3 + 2] = (unsigned char)0;

			KeyPointRGB.RGB[(row * KeyPointRGB.getCol() + (col - i)) * 3] = (unsigned char)255;
			KeyPointRGB.RGB[(row * KeyPointRGB.getCol() + (col - i)) * 3 + 1] = (unsigned char)255;
			KeyPointRGB.RGB[(row * KeyPointRGB.getCol() + (col - i)) * 3 + 2] = (unsigned char)0;

			KeyPointRGB.RGB[((row + i) * KeyPointRGB.getCol() + col) * 3] = (unsigned char)255;
			KeyPointRGB.RGB[((row + i) * KeyPointRGB.getCol() + col) * 3 + 1] = (unsigned char)255;
			KeyPointRGB.RGB[((row + i) * KeyPointRGB.getCol() + col) * 3 + 2] = (unsigned char)0;

			KeyPointRGB.RGB[((row - i) * KeyPointRGB.getCol() + col) * 3] = (unsigned char)255;
			KeyPointRGB.RGB[((row - i) * KeyPointRGB.getCol() + col) * 3 + 1] = (unsigned char)255;
			KeyPointRGB.RGB[((row - i) * KeyPointRGB.getCol() + col) * 3 + 2] = (unsigned char)0;
		}
	}
	std::cout << "KeyPointImg->Size : " << KeyPointRGB.getCol() << " X " << KeyPointRGB.getRow() << std::endl;
	std::fstream img_kp;
	img_kp.open(name, std::ios::binary | std::ios::out);
	img_kp.write((char*)&KeyPointRGB.RGB[0], KeyPointRGB.RGB.size());
}
*/

Raw stack_Img(Raw &img1, Raw &img2)
{
	Raw stacked((img1.getCol() + img2.getCol()), std::max(img1.getRow(), img2.getRow()));

	int Row = stacked.getRow();
	int Col = stacked.getCol();
	for (int i = 0; i < Row; i++)
	{
		int j;
		for (j = 0; j < img1.getCol() * 3; j++)
		{
			if (i < img1.getRow())
				stacked.RGB[i * Col * 3 + j] = img1.RGB[i * img1.getCol() * 3 + j];
			else
				stacked.RGB[i * Col * 3 + j] = (unsigned char)0;
		}
		for (; j < Col * 3; j++)
		{
			if (i < img2.getRow())
				stacked.RGB[i * Col * 3 + j] = img2.RGB[i * img2.getCol() * 3 + (j - (img1.getCol() * 3))];
			else
				stacked.RGB[i * Col * 3 + j] = (unsigned char)0;
		}
	}
	return stacked;
}