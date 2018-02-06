#include <iostream>
#include <vector>
#include <algorithm>
#include <timer.hpp>
using namespace std;

#include "imagedata.hpp"
#include "stitch.hpp"

#define M_PI 3.14159265358979323846

/************************* Local Function Prototypes *************************/
// bilinear補點
struct Color {
	unsigned char R;
	unsigned char G;
	unsigned char B;
};
static Color bilinear(const Raw &src, float _x, float _y)
{
	Color color;
	int x, y;
	x = (int)_x;
	y = (int)_y;

	float l_x = 0.f, r_x = 0.f;
	float t_y = 0.f, b_y = 0.f;

	l_x = _x - (float)x;
	r_x = 1.f - l_x;

	t_y = _y - (float)y;
	b_y = 1.f - t_y;

	float R = 0.f, G = 0.f, B = 0.f;
	int x1 = (x + 1) > (src.getCol() - 1) ? (x + 1) : (src.getCol() - 1);
	int y1 = (y + 1) > (src.getRow() - 1) ? (y + 1) : (src.getRow() - 1);

	if (x >= 0 && x < src.getCol() - 1 && y >= 0 && y < src.getRow() - 1)
	{
		R = (float)src.RGB[((y)* src.getCol() + (x)) * 3 + 0] * (r_x * b_y);
		G = (float)src.RGB[((y)* src.getCol() + (x)) * 3 + 1] * (r_x * b_y);
		B = (float)src.RGB[((y)* src.getCol() + (x)) * 3 + 2] * (r_x * b_y);

		R += (float)src.RGB[((y)* src.getCol() + (x + 1)) * 3 + 0] * (l_x * b_y);
		G += (float)src.RGB[((y)* src.getCol() + (x + 1)) * 3 + 1] * (l_x * b_y);
		B += (float)src.RGB[((y)* src.getCol() + (x + 1)) * 3 + 2] * (l_x * b_y);

		R += (float)src.RGB[((y + 1) * src.getCol() + (x)) * 3 + 0] * (r_x * t_y);
		G += (float)src.RGB[((y + 1) * src.getCol() + (x)) * 3 + 1] * (r_x * t_y);
		B += (float)src.RGB[((y + 1) * src.getCol() + (x)) * 3 + 2] * (r_x * t_y);

		R += (float)src.RGB[((y + 1) * src.getCol() + (x + 1)) * 3 + 0] * (l_x * t_y);
		G += (float)src.RGB[((y + 1) * src.getCol() + (x + 1)) * 3 + 1] * (l_x * t_y);
		B += (float)src.RGB[((y + 1) * src.getCol() + (x + 1)) * 3 + 2] * (l_x * t_y);
	}

	color.R = (unsigned char)(R > 255.0 ? 255 : R < 0.0 ? 0 : R);
	color.G = (unsigned char)(G > 255.0 ? 255 : G < 0.0 ? 0 : G);
	color.B = (unsigned char)(B > 255.0 ? 255 : B < 0.0 ? 0 : B);

	return color;
}
/************************ Functions prototyped here **************************/
// 仿射投影
void WarpPerspective(const Raw &src, Raw &dst, const vector<float> &H) {
	//-------------------------------------
	Color color;
	int iiii = 0, j = 0;
	int d_Row = dst.getRow();
	int d_Col = dst.getCol();
	int s_Row = src.getRow();
	int s_Col = src.getCol();
	float x = 0.0, y = 0.0;
	unsigned char *dst_RGB;
	//-------------------------------------
	// 由投影後的點去回推原圖點，並利用 bilinear 補齊像素
	// http://www.bubufx.com/detail-1797795.html
	//for (iiii = 0; iiii < d_Row; ++iiii)
	for (iiii = -d_Row / 4; iiii < d_Row / 4 * 3; ++iiii)
	{
		dst_RGB = &dst.RGB[(iiii + d_Row / 4) * d_Col * 3];
		//dst_RGB = &dst.RGB[iiii * d_Col * 3];
		//dst_gray = &dst[iiii * d_Col];
		//for (idx = -d_Col / 2; idx < d_Col / 2; ++idx)
		for (j = 0; j < d_Col; ++j)
		{
			x  = (H[2] - H[8] * (float)j) * (H[4] - H[7] * (float)iiii) - (H[1] - H[7] * (float)j) * (H[5] - H[8] * (float)iiii);
			x /= (H[1] - H[7] * (float)j) * (H[3] - H[6] * (float)iiii) - (H[0] - H[6] * (float)j) * (H[4] - H[7] * (float)iiii);
			y  = (H[0] - H[6] * (float)j) * (H[5] - H[8] * (float)iiii) - (H[2] - H[8] * (float)j) * (H[3] - H[6] * (float)iiii);
			y /= (H[1] - H[7] * (float)j) * (H[3] - H[6] * (float)iiii) - (H[0] - H[6] * (float)j) * (H[4] - H[7] * (float)iiii);

			if (x < (float)s_Col && x >= 0.0 && y < (float)s_Row && y >= 0.0)
			{
				color = bilinear(src, x, y);

				dst_RGB[j * 3 + 0] = color.R;
				dst_RGB[j * 3 + 1] = color.G;
				dst_RGB[j * 3 + 2] = color.B;

				//cout << "test=" << (int)dst_RGB[j * 3 + 1] << endl;
			}
		}
	}
}
void _WarpPerspective(const Raw &src, Raw &dst, const vector<float> &H)
{
	//-------------------------------------
	Color color;
	int i = 0, j = 0;
	int d_Row = dst.getRow();
	int d_Col = dst.getCol();
	int s_Row = src.getRow();
	int s_Col = src.getCol();
	float x = 0.0, y = 0.0;
	unsigned char *dst_RGB;
	//-------------------------------------
	// 由投影後的點去回推原圖點，並利用 bilinear 補齊像素
	// http://www.bubufx.com/detail-1797795.html
	//for (iiii = 0; iiii < d_Row; ++iiii)
	for (i = 0; i < d_Row; ++i)
	{
		dst_RGB = &dst.RGB[i * d_Col * 3];
		for (j = 0; j < d_Col; ++j)
		{
			x = (H[2] - H[8] * (float)j) * (H[4] - H[7] * (float)i) - (H[1] - H[7] * (float)j) * (H[5] - H[8] * (float)i);
			x /= (H[1] - H[7] * (float)j) * (H[3] - H[6] * (float)i) - (H[0] - H[6] * (float)j) * (H[4] - H[7] * (float)i);
			y = (H[0] - H[6] * (float)j) * (H[5] - H[8] * (float)i) - (H[2] - H[8] * (float)j) * (H[3] - H[6] * (float)i);
			y /= (H[1] - H[7] * (float)j) * (H[3] - H[6] * (float)i) - (H[0] - H[6] * (float)j) * (H[4] - H[7] * (float)i);

			if (x < (float)s_Col && x >= 0.0 && y < (float)s_Row && y >= 0.0)
			{
				color = bilinear(src, x, y);

				dst_RGB[j * 3 + 0] = color.R;
				dst_RGB[j * 3 + 1] = color.G;
				dst_RGB[j * 3 + 2] = color.B;
			}
		}
	}
}
//-------------------------------------
// 圓柱投影
Raw DealWithImgData(Raw &srcdata, int width, int height, float R)
{
	//-------------------------------------
	int hnum = 0, wnum = 0;
	float x = 0.f, y = 0.f;
	Raw drcdata(width, height * 2);
	Color color;
	unsigned char *drcdata_RGB;
	//-------------------------------------
	// 由投影後的點去回推原圖點，並利用 bilinear 補齊像素
	// http://blog.csdn.net/weixinhum/article/details/50611750
	float  k = 0.f;
	for (hnum = -(height / 2); hnum < (height * 3 / 2); hnum++)
	{
		drcdata_RGB = &drcdata.RGB[(hnum + (height / 2)) * width * 3];
		for (wnum = 0; wnum < width; wnum++)
		{
			float k = R / sqrt(R * R + ((float)wnum - (float)width / 2.f) * ((float)wnum - (float)width / 2.f));
			x = ((float)wnum - (float)width / 2.f) / k + (float)width / 2.f;
			y = ((float)hnum - (float)height / 2.f) / k + (float)height / 2.f;

			if (x >= 0 && y >= 0 && x < width - 1 && y < height - 1)
			{
				color = bilinear(srcdata, x, y);

				drcdata_RGB[wnum * 3 + 0] = color.R;
				drcdata_RGB[wnum * 3 + 1] = color.G;
				drcdata_RGB[wnum * 3 + 2] = color.B;
			}
		}
	}
	return drcdata;
}

void warping(const vector<Raw> &inputArrays, float FL2, 
	vector<Raw> &Output, vector<fpoint> &upedge, vector<fpoint> &downedge)
{
	float FL = FL2;
	Output.resize(inputArrays.size());

	Timer t;
	for (int idx = 0; idx < inputArrays.size(); idx++)
	{
		const Raw& image = inputArrays[idx];
		Output[idx].resize(image.getCol(), image.getRow());
		float mid_x = (float)image.getCol() / 2.f;
		float mid_y = (float)image.getRow() / 2.f;
		size_t curr_h = image.getRow();
		size_t curr_w = image.getCol();

		vector<float> k_core(curr_w);
		Timer t;
		for(int j = 0; j < curr_h; j++){
			for(int i = 0; i < curr_w; i++){
				if(j == 0)
					k_core[i] = sqrtf(FL * FL + (i - mid_x) * (i - mid_x)) / FL;
				float k = k_core[i];
				float x = (i-mid_x)*k + mid_x;
				float y = (j-mid_y)*k + mid_y;

				if(x >= 0 && y >= 0 &&
					x < curr_w && y < curr_h)
				{
					Color color = bilinear(image, x, y);
					Output[idx].RGB[(j*curr_w + i) * 3 + 0] = color.R;
					Output[idx].RGB[(j*curr_w + i) * 3 + 1] = color.G;
					Output[idx].RGB[(j*curr_w + i) * 3 + 2] = color.B;
				}
				if(j == 0) { // 圖上邊邊界.
					upedge.emplace_back(fpoint((float)x, (float)y));
				} else if(j == curr_h - 1) { // 圖下邊邊界.
					downedge.emplace_back(fpoint((float)x, (float)y));
				}
			}
		}
	}
	//t.print("warping time = ");
}




