
// #include "gaussian.hpp"
// #include "Bmp.hpp"
// #include "bicubic.hpp"
// #include "First_order_GPU.hpp"
// #include "Convolution.hpp"

#include <iostream>
#include <algorithm>
#include <fstream>
#include <time.h>
#include <vector>
#include <string>
using namespace std;

#include "imglib\imglib.hpp"
#include "Raw2Img\Raw2Img.hpp"
#include "imglib\imglib.hpp"
#include "Imgraw.hpp"

#include "imagedata.hpp"
#include "Blend.hpp"

#define RIGHT 1
#define LEFT 0
//----------------------------------------
// �h��
#define PYR_OCTAVE 5
//----------------------------------------
// �����b�|�A�P�̫᪺�����v�����b�|���@��
#define PYR_R 2
#define M_PI 3.14159265358979323846

struct X_S
{
	int m;
	int n;
	int middle;
	int l;
};

using namespace std;
/************************* Local Function Prototypes *************************/
/*static vector<float> getGaussianKernel(int);
static Blend_Image sub(const Blend_Image&, const Blend_Image&);
static Blend_Image add(const Blend_Image&, const Blend_Image&);
static vector<struct X_S> getXsize(Blend_Image);

vector<bool> buildLaplacianMap(const Raw &, vector<Blend_Image> &, int, int, int);
vector<float> getGaussianKernel_rr(int, int, int, int);
void blendImg(Raw &, const Blend_Image &, int, int, int, const vector<bool> &, const vector<bool> &);*/


/****************************** Local Function *******************************/
// �إ߼e�פj�p�Ӱ���kernel
static vector<float> getGaussianKernel(int x)
{
	vector<float> kernel(x, 0.f);
	float sigma; // sigma �ѥb�|�h�p��

	float half = ((float)x - 1.f) / 2.f;
	sigma = sqrt((-1.f) * pow((float)x - 1.f - half, 2.f) / (2.f * log(0.5f)));

	for(int i = 0; i < x; i++){
		float g = 0.f;
		if(i <= (x - half)){
			g = exp((-1.f) * (float)(i * i) / (2.f * sigma * sigma));
		} else{
			g = 1.f - exp((-1.f) * pow((float)(x - i) - 1.f, 2.f) / (2.f * sigma * sigma));
		}
		kernel[i] = g;
	}
	return kernel;
}
//----------------------------------------
// ��ϼh�۴�
static Blend_Image sub(const Blend_Image& m, const Blend_Image& n)
{
	Blend_Image out;
	out.width = m.width;
	out.height = m.height;
	out.RGB = new float[out.width * out.height * 3];
	int i = 0, j = 0;
	for (i = 0; i < out.height; i++)
	{
		for (j = 0; j < out.width; j++)
		{
			if (i < n.height && j < n.width)
			{
				out.RGB[(i * out.width + j) * 3 + 0] = m.RGB[(i * m.width + j) * 3 + 0] - n.RGB[(i * n.width + j) * 3 + 0];
				out.RGB[(i * out.width + j) * 3 + 1] = m.RGB[(i * m.width + j) * 3 + 1] - n.RGB[(i * n.width + j) * 3 + 1];
				out.RGB[(i * out.width + j) * 3 + 2] = m.RGB[(i * m.width + j) * 3 + 2] - n.RGB[(i * n.width + j) * 3 + 2];
			}
		}
	}
	return out;
}
//----------------------------------------
// ��ϼh�ۥ[
static Blend_Image add(const Blend_Image& m, const Blend_Image& n)
{
	Blend_Image out;
	out.width = m.width;
	out.height = m.height;
	out.RGB = new float[out.width * out.height * 3];
	int i = 0, j = 0;
	//float *out_RGB, *m_RGB, *n_RGB;
	for (i = 0; i < out.height; i++)
	{
		const auto& out_RGB = &out.RGB[i * out.width * 3];
		const auto& m_RGB = &m.RGB[i * m.width * 3];
		const auto& n_RGB = &n.RGB[i * n.width * 3];

		for (j = 0; j < out.width; j++)
		{
			out_RGB[j * 3 + 0] = m_RGB[j * 3 + 0] + n_RGB[j * 3 + 0];
			out_RGB[j * 3 + 1] = m_RGB[j * 3 + 1] + n_RGB[j * 3 + 1];
			out_RGB[j * 3 + 2] = m_RGB[j * 3 + 2] + n_RGB[j * 3 + 2];
			out_RGB[j * 3 + 0] = out_RGB[j * 3 + 0] < 0.f ? 0.f : out_RGB[j * 3 + 0]>255.f ? 255.f : out_RGB[j * 3 + 0];
			out_RGB[j * 3 + 1] = out_RGB[j * 3 + 1] < 0.f ? 0.f : out_RGB[j * 3 + 1]>255.f ? 255.f : out_RGB[j * 3 + 1];
			out_RGB[j * 3 + 2] = out_RGB[j * 3 + 2] < 0.f ? 0.f : out_RGB[j * 3 + 2]>255.f ? 255.f : out_RGB[j * 3 + 2];
		}
	}
	return out;
}
//----------------------------------------
// �p��C�Ӱ��׸̪�X�e��
static vector<struct X_S> getXsize(Blend_Image input)
{
	vector<struct X_S> temp;
	temp.resize(input.height);
	int m = 0, n = 0;
	int i = 0, j = 0;
	for (i = 0; i < input.height; i++)
	{
		m = 0;
		n = 0;
		for (j = 0; j < input.width; j++)
		{
			if ((input.RGB[(i * input.width + j) * 3 + 0] != 0.f) ||
				(input.RGB[(i * input.width + j) * 3 + 1] != 0.f) ||
				(input.RGB[(i * input.width + j) * 3 + 2] != 0.f) )
			{
				m = j;
				break;
			}
		}
		for (j = input.width - 1; j >= 0; j--)
		{
			if ((input.RGB[(i * input.width + j) * 3 + 0] != 0.f) ||
				(input.RGB[(i * input.width + j) * 3 + 1] != 0.f) ||
				(input.RGB[(i * input.width + j) * 3 + 2] != 0.f))
			{
				n = j;
				break;
			}
		}
		temp[i].m = m;
		temp[i].n = n;
		if (m == 0 && n == 0)
		{
			temp[i].l = 0;
			temp[i].middle = 0;
		}
		else
		{
			temp[i].l = (n - m) + 1;
			temp[i].middle = (n + m) / 2;
		}
	}
	return temp;
}
//----------------------------------------
// �����x�}
static vector<float> getGaussianKernel_rr(int x, int y, int dx, int dy = 0)
{
	vector<float> kernel(x * y, 0.f);
	float half = (dx - 1) / 2.f;
	float sigma = sqrt((-1.f) * pow((float)x - 1.f - half, 2.f) / (2.f * std::log(0.5f)));
	for(int i = (x - dx); i < x; i++){
		float g;
		if(i <= (x - half)){
			g = exp((-1) * i * i / (2 * sigma * sigma));
		} else{
			g = 1.f - exp((-1.f) * pow(x - i - 1.f, 2.f) / (2.f * sigma * sigma));
		}

		for(int j = 0; j < y; j++){
			kernel[j * x + i] = g;
		}
	}
	return kernel;
}
//----------------------------------------


static float BiCubicPoly(float x) {
	float abs_x = abs(x);
	float a = -0.5f;
	if(abs_x <= 1.f){
		return (a + 2.f) * pow(abs_x, 3) - (a + 3.f) * pow(abs_x, 2) + 1.f;
	} else if(abs_x < 2.0){
		return a * pow(abs_x, 3) - 5.f * a * pow(abs_x, 2) + 8.f * a * abs_x - 4.f * a;
	}
	return 0.f;
}
// ��g�ഫ+�u�ʸɭ�.
static Blend_Image Bicubic(const Blend_Image& src, float TransMat[3][3])
{
	Blend_Image dst;
	// calculate margin point of dst image  
	float left = 0;
	float right = 0;
	float top = 0;
	float down = 0;

	float x = src.width * 1.0f;
	float y = 0.0f;
	float u1 = x * TransMat[0][0] + y * TransMat[0][1];
	float v1 = x * TransMat[1][0] + y * TransMat[1][1];
	x = src.width * 1.0f;
	y = src.height * 1.0f;
	float u2 = x * TransMat[0][0] + y * TransMat[0][1];
	float v2 = x * TransMat[1][0] + y * TransMat[1][1];
	x = 0.0f;
	y = src.width * 1.0f;
	float u3 = x * TransMat[0][0] + y * TransMat[0][1];
	float v3 = x * TransMat[1][0] + y * TransMat[1][1];

	left = min(min(min(0.0f, u1), u2), u3);
	right = max(max(max(0.0f, u1), u2), u3);
	top = min(min(min(0.0f, v1), v2), v3);
	down = max(max(max(0.0f, v1), v2), v3);

	// create dst image  
	int t_size = int(abs(right - left)) * int(abs(down - top));
	dst.width = int(abs(right - left));
	dst.height = int(abs(down - top));
	dst.RGB = new float[t_size * 3];

	int i, j;
	float* p;
	float* q0;
	float* q1;
	float* q2;
	float* q3;

	int x0 = 0, x1 = 0, x2 = 0, x3 = 0;
	int y0 = 0, y1 = 0, y2 = 0, y3 = 0;
	float dist_x0 = 0.f, dist_x1 = 0.f, dist_x2 = 0.f, dist_x3 = 0.f;
	float dist_y0 = 0.f, dist_y1 = 0.f, dist_y2 = 0.f, dist_y3 = 0.f;
	float dist_x0y0 = 0.f, dist_x0y1 = 0.f, dist_x0y2 = 0.f, dist_x0y3 = 0.f;
	float dist_x1y0 = 0.f, dist_x1y1 = 0.f, dist_x1y2 = 0.f, dist_x1y3 = 0.f;
	float dist_x2y0 = 0.f, dist_x2y1 = 0.f, dist_x2y2 = 0.f, dist_x2y3 = 0.f;
	float dist_x3y0 = 0.f, dist_x3y1 = 0.f, dist_x3y2 = 0.f, dist_x3y3 = 0.f;
	float thre = 0.f;

	for (i = 0; i < dst.height; ++i)
	{
		p = &dst.RGB[i * dst.width * 3];
		for (j = 0; j < dst.width; ++j)
		{  
			x = (j + left) / TransMat[0][0];
			y = (i + top) / TransMat[1][1];

			x0 = int(x) - 1;
			y0 = int(y) - 1;
			x1 = int(x);
			y1 = int(y);
			x2 = int(x) + 1;
			y2 = int(y) + 1;
			x3 = int(x) + 2;
			y3 = int(y) + 2;

			if ((x0 >= 0) && (x3 < src.width) && (y0 >= 0) && (y3 < src.height))
			{      
				q0 = &src.RGB[y0 * src.width * 3];
				q1 = &src.RGB[y1 * src.width * 3];
				q2 = &src.RGB[y2 * src.width * 3];
				q3 = &src.RGB[y3 * src.width * 3];
				dist_x0 = BiCubicPoly(x - x0);
				dist_x1 = BiCubicPoly(x - x1);
				dist_x2 = BiCubicPoly(x - x2);
				dist_x3 = BiCubicPoly(x - x3);
				dist_y0 = BiCubicPoly(y - y0);
				dist_y1 = BiCubicPoly(y - y1);
				dist_y2 = BiCubicPoly(y - y2);
				dist_y3 = BiCubicPoly(y - y3);

				dist_x0y0 = dist_x0 * dist_y0;
				dist_x0y1 = dist_x0 * dist_y1;
				dist_x0y2 = dist_x0 * dist_y2;
				dist_x0y3 = dist_x0 * dist_y3;
				dist_x1y0 = dist_x1 * dist_y0;
				dist_x1y1 = dist_x1 * dist_y1;
				dist_x1y2 = dist_x1 * dist_y2;
				dist_x1y3 = dist_x1 * dist_y3;
				dist_x2y0 = dist_x2 * dist_y0;
				dist_x2y1 = dist_x2 * dist_y1;
				dist_x2y2 = dist_x2 * dist_y2;
				dist_x2y3 = dist_x2 * dist_y3;
				dist_x3y0 = dist_x3 * dist_y0;
				dist_x3y1 = dist_x3 * dist_y1;
				dist_x3y2 = dist_x3 * dist_y2;
				dist_x3y3 = dist_x3 * dist_y3;

				p[3 * j + 0] = (q0[3 * x0] * dist_x0y0 +
					q1[3 * x0] * dist_x0y1 +
					q2[3 * x0] * dist_x0y2 +
					q3[3 * x0] * dist_x0y3 +
					q0[3 * x1] * dist_x1y0 +
					q1[3 * x1] * dist_x1y1 +
					q2[3 * x1] * dist_x1y2 +
					q3[3 * x1] * dist_x1y3 +
					q0[3 * x2] * dist_x2y0 +
					q1[3 * x2] * dist_x2y1 +
					q2[3 * x2] * dist_x2y2 +
					q3[3 * x2] * dist_x2y3 +
					q0[3 * x3] * dist_x3y0 +
					q1[3 * x3] * dist_x3y1 +
					q2[3 * x3] * dist_x3y2 +
					q3[3 * x3] * dist_x3y3);

				p[3 * j + 1] = (q0[3 * x0 + 1] * dist_x0y0 +
					q1[3 * x0 + 1] * dist_x0y1 +
					q2[3 * x0 + 1] * dist_x0y2 +
					q3[3 * x0 + 1] * dist_x0y3 +
					q0[3 * x1 + 1] * dist_x1y0 +
					q1[3 * x1 + 1] * dist_x1y1 +
					q2[3 * x1 + 1] * dist_x1y2 +
					q3[3 * x1 + 1] * dist_x1y3 +
					q0[3 * x2 + 1] * dist_x2y0 +
					q1[3 * x2 + 1] * dist_x2y1 +
					q2[3 * x2 + 1] * dist_x2y2 +
					q3[3 * x2 + 1] * dist_x2y3 +
					q0[3 * x3 + 1] * dist_x3y0 +
					q1[3 * x3 + 1] * dist_x3y1 +
					q2[3 * x3 + 1] * dist_x3y2 +
					q3[3 * x3 + 1] * dist_x3y3);

				p[3 * j + 2] = (q0[3 * x0 + 2] * dist_x0y0 +
					q1[3 * x0 + 2] * dist_x0y1 +
					q2[3 * x0 + 2] * dist_x0y2 +
					q3[3 * x0 + 2] * dist_x0y3 +
					q0[3 * x1 + 2] * dist_x1y0 +
					q1[3 * x1 + 2] * dist_x1y1 +
					q2[3 * x1 + 2] * dist_x1y2 +
					q3[3 * x1 + 2] * dist_x1y3 +
					q0[3 * x2 + 2] * dist_x2y0 +
					q1[3 * x2 + 2] * dist_x2y1 +
					q2[3 * x2 + 2] * dist_x2y2 +
					q3[3 * x2 + 2] * dist_x2y3 +
					q0[3 * x3 + 2] * dist_x3y0 +
					q1[3 * x3 + 2] * dist_x3y1 +
					q2[3 * x3 + 2] * dist_x3y2 +
					q3[3 * x3 + 2] * dist_x3y3);

				thre = 198.f;
				if ((abs(p[3 * j] - q1[3 * x1]) > thre) || (abs(p[3 * j + 1] - q1[3 * x1 + 1]) > thre) ||
					(abs(p[3 * j + 2] - q1[3 * x1 + 2]) > thre))
				{
					p[3 * j] = q1[3 * x1];
					p[3 * j + 1] = q1[3 * x1 + 1];
					p[3 * j + 2] = q1[3 * x1 + 2];
				}
			}
			else
			{
				p[j * 3 + 0] = 0.f;
				p[j * 3 + 1] = 0.f;
				p[j * 3 + 2] = 0.f;
			}
		}
	}
	return dst;
}

//�@��������.
#define KERNEL_RADIUS 8
#define KERNEL_LENGTH (2 * KERNEL_RADIUS + 1)
float *GaussianKernel1D(int length, float sigma)
{
	cout << sigma << endl;
	float *h_Kernel;
	h_Kernel = new float[KERNEL_LENGTH];

	double s2 = sigma * sigma;
	double m = 1.0 / (sqrt(2.0 * M_PI) * sigma);
	double v = 0.0;
	int i = 0;

	for (i = 0; i < KERNEL_LENGTH; ++i)
	{
		h_Kernel[i] = 0.0;
	}

	int radius = length / 2;

	for (i = 0; i <= radius; ++i)
	{
		v = m * exp(-(1.0 * i * i) / (2.0 * s2));
		h_Kernel[radius + i] = (float)v;
		h_Kernel[radius - i] = (float)v;
	}
	float total = 0.0;
	for (i = 0; i < length; ++i)
	{
		total += h_Kernel[i];
	}
	for (i = 0; i < length; ++i)
	{
		h_Kernel[i] /= total;
		//	cout << h_Kernel[i] << endl;
	}
	//system("pause");
	return h_Kernel;
}

// �����ҽk
static void GauBlur3d(vector<float>& img_gau, const vector<float>& img_ori,
	size_t width, size_t height, float p, size_t mat_len)
{
	// �ӷ����i�ۦP
	if (&img_gau == &img_ori) {
		throw file_same("## Erroe! in and out is same.");
	}
	vector<float> gau_mat = Gaus::gau_matrix(p, mat_len);
	img_gau.resize(height*width*3);
	// �w�s
	vector<float> img_gauX(img_ori.size());
	// �����ҽk X �b
	const size_t r = gau_mat.size() / 2;
	for (unsigned j = 0; j < height; ++j) {
		for (unsigned i = 0; i < width; ++i) {
			double sumR = 0;
			double sumG = 0;
			double sumB = 0;
			for (unsigned k = 0; k < gau_mat.size(); ++k) {
				int idx = i-r + k;
				// idx�W�X��t�B�z
				if (idx < 0) {
					idx = 0;
				} else if (idx >(int)(width-1)) {
					idx = (width-1);
				}
				sumR += img_ori[(j*width + idx)*3 + 0] * gau_mat[k];
				sumG += img_ori[(j*width + idx)*3 + 1] * gau_mat[k];
				sumB += img_ori[(j*width + idx)*3 + 2] * gau_mat[k];
			}
			img_gauX[(j*width + i)*3 + 0] = sumR;
			img_gauX[(j*width + i)*3 + 1] = sumG;
			img_gauX[(j*width + i)*3 + 2] = sumB;
		}
	}
	// �����ҽk Y �b
	for (unsigned j = 0; j < height; ++j) {
		for (unsigned i = 0; i < width; ++i) {
			double sumR = 0;
			double sumG = 0;
			double sumB = 0;
			for (unsigned k = 0; k < gau_mat.size(); ++k) {
				int idx = j-r + k;
				// idx�W�X��t�B�z
				if (idx < 0) {
					idx = 0;
				} else if (idx > (int)(height-1)) {
					idx = (height-1);
				}
				sumR += img_gauX[(idx*width + i)*3 + 0] * gau_mat[k];
				sumG += img_gauX[(idx*width + i)*3 + 2] * gau_mat[k];
				sumB += img_gauX[(idx*width + i)*3 + 3] * gau_mat[k];

			}
			img_gau[(j*width + i)*3 + 0] = sumR;
			img_gau[(j*width + i)*3 + 1] = sumG;
			img_gau[(j*width + i)*3 + 2] = sumB;
		}
	}
}
static Blend_Image BlurImage(Blend_Image &img, float sigma, int r) {
	size_t w = img.width;
	size_t h = img.height;

	vector<float> img_ori(w*h*3);
	for(size_t i = 0; i < img_ori.size(); i++) {
		img_ori[i] = img.RGB[i];
	}

	vector<float> img_gau;
	GauBlur3d(img_gau, img_ori, w, h, sigma, (r*2 + 1));

	Blend_Image dst;
	dst.width = img.width;
	dst.height = img.height;
	dst.RGB = new float[dst.width * dst.height * 3];
	for(size_t i = 0; i < img_gau.size(); i++) {
		dst.RGB[i] = img_gau[i] ;
	}
	return dst;
}

/*********************** Functions prototyped in Blend.h **********************/
// ��B�z�n�����|�ϼg�J���.
static void blendImg(Raw &inputArray, const Blend_Image &overlap_area, int dx, int dy, int lr, const vector<bool> &bol, const vector<bool> &bor)
{
	int disx = (lr == RIGHT) ? 0 : (inputArray.getCol() - dx);
	int disy = dy >= 0 ? ((lr == RIGHT) ? 0 : (inputArray.getRow() - dy)) : ((lr == RIGHT) ? (inputArray.getRow() + dy) : 0);

	if (disy < 0) { disy = 0; }
	if (disx < 0) { disx = 0; }

	for (int j = 0; j < overlap_area.height; j++){
		auto&& currW = overlap_area.width;
		for (int i = 0; i < currW - 1; i++)
		{
			// ��X���|�I.
			if (bol[j*currW + i]     == true && bor[j*currW + i]     == true 
				&& bol[j*currW + (i+1)] == true && bor[j*currW + (i+1)] == true)
			{
				//todo 0203��z��o��
				if (j + disy < inputArray.getRow() && i + disx < inputArray.getCol())
				{
					/*if (flag == false)
					{
					inputArray.RGB[((j + disy) * inputArray.getCol() + (i + disx)) * 3 + 0] = 0xff;
					inputArray.RGB[((j + disy) * inputArray.getCol() + (i + disx)) * 3 + 0] = 0x00;
					inputArray.RGB[((j + disy) * inputArray.getCol() + (i + disx)) * 3 + 0] = 0x00;
					}*/
					for (int c = 0; c < 3; c++)
					{
						inputArray.RGB[((j + disy)*inputArray.getCol() + (i + disx)) * 3 + c] = (unsigned char)
							// �ץ��W�X�W�U���� (�O���ܦP�@��ѨM�A���O���ܬݤU�@��).
							((overlap_area.RGB[(j*currW + i)*3 + c] < 0.f)? 0:
							(overlap_area.RGB[(j*currW + i)*3 + c] > 255.f)? 255:
								overlap_area.RGB[(j*currW + i)*3 + c]);
					}
				}
			}
		}
	}
}
// �إ� Laplacian ���r��.
static vector<bool> buildLaplacianMap(const Raw &inputArray, vector<Blend_Image> &outputArrays, int dx, int dy, int chkLR)
{
	vector<bool> OverlapBool;
	Blend_Image tmp;
	tmp.height = abs(dy);
	tmp.width = dx;
	tmp.RGB = new float[tmp.width * tmp.height * 3];

	auto&& imgW = inputArray.getCol();
	auto&& imgH = inputArray.getRow();

	// ���|�ϰ_�l��m (�ҥH�k�Ϥ~�O�q0, 0�}�l)
	int disx = (chkLR == RIGHT) ? 0 : (inputArray.getCol() - dx);
	int disy = (dy >= 0) ?
		((chkLR == RIGHT) ? 0 : (inputArray.getRow() - dy)) : 
		((chkLR == RIGHT) ? (inputArray.getRow() + dy) : 0);
	//cout << "dis (x, y) = " << disx << ", " << disy << endl;;

	auto&& Blend_to_imgraw = [](const Blend_Image& src){
		vector<unsigned char> temp(src.width*src.height*3);
		for(size_t i = 0; i < temp.size(); i++){
			temp[i]=src.RGB[i];
		}
		ImgRaw dst(temp, src.width, src.height, 24);
		return dst;
	};

	if (disx < 0) { disx = 0; }
	// copy�X���|��
	for(int j = 0; j < tmp.height; j++) {
		for(int i = 0; i < tmp.width; i++) {
			bool flag = false;
			if( (j + disy) < imgH && 
				(i + disx) < imgW)
			{
				for(int c = 0; c < 3; c++) {
					tmp.RGB[(j*tmp.width + i) * 3 + c] = (float)inputArray.RGB[((j+disy)*imgW + (i+disx)) * 3 + c];
					if(inputArray.RGB[((j + disy) * imgW + (i + disx)) * 3 + c] != 0) {
						flag = true;
					}
				}
			}
			OverlapBool.push_back(flag);
		}
	}
	//Blend_to_imgraw(tmp).bmp("overlap0 .bmp");
	//system("pause");
	

	// ��¦�ϰ�ɦ�����(���C�W�U����W���������۵M)
	int set_r = 1;
	int ur = 0, dr = 0;
	for(int j = 0; j < tmp.width; j++) {
		for(int i = 0; i < tmp.height; i++) {
			if(tmp.RGB[(i * tmp.width + j) * 3 + 0] != 0.f || tmp.RGB[(i * tmp.width + j) * 3 + 1] != 0.f || tmp.RGB[(i * tmp.width + j) * 3 + 2] != 0.f) {
				ur = i;
				break;
			}
		}
		for(int i = ur; i >= 0; i--) {
			tmp.RGB[(i * tmp.width + j) * 3 + 0] = tmp.RGB[((ur + set_r) * tmp.width + j) * 3 + 0];
			tmp.RGB[(i * tmp.width + j) * 3 + 1] = tmp.RGB[((ur + set_r) * tmp.width + j) * 3 + 1];
			tmp.RGB[(i * tmp.width + j) * 3 + 2] = tmp.RGB[((ur + set_r) * tmp.width + j) * 3 + 2];
		}
		for(int i = tmp.height - 1; i >= 0; i--) {
			if(tmp.RGB[(i * tmp.width + j) * 3 + 0] != 0.f || tmp.RGB[(i * tmp.width + j) * 3 + 1] != 0.f || tmp.RGB[(i * tmp.width + j) * 3 + 2] != 0.f) {
				dr = i;
				break;
			}
		}
		dr = dr < 1 ? 1 : dr;
		for(int i = dr; i < tmp.height; i++) {
			tmp.RGB[(i * tmp.width + j) * 3 + 0] = tmp.RGB[((dr - set_r) * tmp.width + j) * 3 + 0];
			tmp.RGB[(i * tmp.width + j) * 3 + 1] = tmp.RGB[((dr - set_r) * tmp.width + j) * 3 + 1];
			tmp.RGB[(i * tmp.width + j) * 3 + 2] = tmp.RGB[((dr - set_r) * tmp.width + j) * 3 + 2];
		}
		//-------------------------------------------------------------------------------------------------------------------------------------------------
	}
	for(int i = 0; i < tmp.height; i++) {
		for(int j = 0; j < tmp.width; j++) {
			if(tmp.RGB[(i * tmp.width + j) * 3 + 0] != 0.f || tmp.RGB[(i * tmp.width + j) * 3 + 1] != 0.f || tmp.RGB[(i * tmp.width + j) * 3 + 2] != 0.f) {
				ur = j;
				break;
			}
		}
		for(int j = ur; j >= 0; j--) {
			tmp.RGB[(i * tmp.width + j) * 3 + 0] = tmp.RGB[(i * tmp.width + (ur + set_r)) * 3 + 0];
			tmp.RGB[(i * tmp.width + j) * 3 + 1] = tmp.RGB[(i * tmp.width + (ur + set_r)) * 3 + 1];
			tmp.RGB[(i * tmp.width + j) * 3 + 2] = tmp.RGB[(i * tmp.width + (ur + set_r)) * 3 + 2];
		}
		for(int j = tmp.width - 1; j >= 0; j--) {
			if(tmp.RGB[(i * tmp.width + j) * 3 + 0] != 0.f || tmp.RGB[(i * tmp.width + j) * 3 + 1] != 0.f || tmp.RGB[(i * tmp.width + j) * 3 + 2] != 0.f) {
				dr = j;
				break;
			}
		}
		dr = dr < 1 ? 1 : dr;
		for(int j = dr; j < tmp.width; j++) {
			tmp.RGB[(i * tmp.width + j) * 3 + 0] = tmp.RGB[(i * tmp.width + (dr - set_r)) * 3 + 0];
			tmp.RGB[(i * tmp.width + j) * 3 + 1] = tmp.RGB[(i * tmp.width + (dr - set_r)) * 3 + 1];
			tmp.RGB[(i * tmp.width + j) * 3 + 2] = tmp.RGB[(i * tmp.width + (dr - set_r)) * 3 + 2];
		}
		//-------------------------------------------------------------------------------------------------------------------------------------------------
	}
	

	// TODO 0207 �o�̴N�X�ƤF2 �ݤ@�U temp �A�ĤT����
	outputArrays.clear();
	outputArrays.resize(PYR_OCTAVE);
	outputArrays[0] = tmp;

	float sigma = sqrt((-1) * pow(2.f, 2.f) / (2.f * log(0.5f)));
	// �Y�p���ܴ��x�}.
	float transMat_0_5[3][3] = {
		{ 0.5f,  0.f, 0.f },
		{  0.f, 0.5f, 0.f },
		{  0.f,  0.f, 1.f }
	};
	//Blend_to_imgraw(outputArrays[0]).bmp("overlap1 .bmp");

	//system("pause");


	for(int i = 1; i < PYR_OCTAVE; i++) {
		// �ҽk�Ϥ�.
		const Blend_Image&& blurImg = BlurImage(outputArrays[i - 1], sigma, PYR_R);

		// TODO �o�̴N�X�ƤF
		// Blend_to_imgraw(blurImg).bmp("blurImg .bmp");
		// �Y�p�Ϥ�.
		const Blend_Image&& scalImg = Bicubic(blurImg, transMat_0_5);
		outputArrays[i] = scalImg;
	}
	// ��j���ܴ��x�}.
	float transMat_2[3][3] = {
		{ 2.0f,  0.f, 0.f },
		{  0.f, 2.0f, 0.f },
		{  0.f,  0.f, 1.f }
	};
	for(int a = 0; a < PYR_OCTAVE - 1; a++) {
		// �Y�p�Ϥ�.
		const Blend_Image&& scalImg = Bicubic(outputArrays[a+1], transMat_2);
		outputArrays[a] = sub(outputArrays[a], scalImg);
	}

	return OverlapBool;
}
// �h�W�q�V�X�A��J��L�P��R�A�٦���R�������q.
void multiBandBlend(Raw &limg, Raw &rimg, int dx, int dy)
{
	// �ڧ�L�g�X�h�F�A�@�}�l����N�ץ��n
	/*if(dx % 2 == 0) {
		if(dx + 1 <= limg.getCol() && dx + 1 <= rimg.getCol()) {
			dx += 1;
		} else {
			dx -= 1;
		}
	}
	if(dy % 2 == 0) {
		if(dy + 1 <= limg.getRow() && dy + 1 <= rimg.getRow()) {
			dy += 1;
		} else {
			dy -= 1;
		}
	} else if(dy % 2 == 1){
		dy += 1;
	}*/


	vector<Blend_Image> llpyr, rlpyr;
	vector<bool> bol, bor; // ��Ϲ������|�ϰ�A�άO�X��
	bol = buildLaplacianMap(limg, llpyr, dx, dy, LEFT);
	bor = buildLaplacianMap(rimg, rlpyr, dx, dy, RIGHT);

	int center = 0;
	int i, c;
	vector<Blend_Image> LS(PYR_OCTAVE);
	vector<float> k = getGaussianKernel_rr(llpyr[llpyr.size() - 1].width, llpyr[llpyr.size() - 1].height, llpyr[llpyr.size() - 1].width, 0);
	for (int a = 0; a < llpyr.size(); a++)
	{
		LS[a].width = llpyr[a].width;
		LS[a].height = llpyr[a].height;
		LS[a].RGB = new float[LS[a].width * LS[a].height * 3];
		center = (int)(llpyr[a].width / 2);

		for(int j = 0; j < LS[a].height; j++) {
			for(i = 0; i < LS[a].width; i++) {
				for(c = 0; c < 3; c++) {
					if(a == llpyr.size() - 1) {
						LS[a].RGB[(j * LS[a].width + i) * 3 + c] = llpyr[a].RGB[(j * llpyr[a].width + i) * 3 + c] * k[j * llpyr[a].width + i];
						LS[a].RGB[(j * LS[a].width + i) * 3 + c] += rlpyr[a].RGB[(j * rlpyr[a].width + i) * 3 + c] * (1.f - k[j * llpyr[a].width + i]);
					} else {
						if(i == center) {
							LS[a].RGB[(j * LS[a].width + i) * 3 + c] = (llpyr[a].RGB[(j * llpyr[a].width + i) * 3 + c] + rlpyr[a].RGB[(j * rlpyr[a].width + i) * 3 + c]) / 2.f;
						} else if(i > center) {
							LS[a].RGB[(j * LS[a].width + i) * 3 + c] = rlpyr[a].RGB[(j * rlpyr[a].width + i) * 3 + c];
						} else {
							LS[a].RGB[(j * LS[a].width + i) * 3 + c] = llpyr[a].RGB[(j * llpyr[a].width + i) * 3 + c];
						}
					}
				}
			}
		}
	}
	float transMat_2[3][3] = { { 2.0f, 0.f, 0.f },{ 0.f, 2.0f, 0.f },{ 0.f, 0.f, 1.f } };
	Blend_Image result;
	for(int a = PYR_OCTAVE - 1; a > 0; a--) {
		result = Bicubic(LS[a], transMat_2);

		for(int j = 0; j < LS[a - 1].height; j++) {
			for(i = 0; i < LS[a - 1].width; i++) {
				for(c = 0; c < 3; c++) {
					if(j < result.height && i < result.width) {
						LS[a - 1].RGB[(j * LS[a - 1].width + i) * 3 + c] = 
							LS[a - 1].RGB[(j * LS[a - 1].width + i) * 3 + c] + result.RGB[(j * result.width + i) * 3 + c];
						LS[a - 1].RGB[(j * LS[a - 1].width + i) * 3 + c] = 
							LS[a - 1].RGB[(j * LS[a - 1].width + i) * 3 + c] < 0.f ? 
								0.f : LS[a - 1].RGB[(j * LS[a - 1].width + i) * 3 + c] > 255.f ? 
								255.f : LS[a - 1].RGB[(j * LS[a - 1].width + i) * 3 + c];
					}
				}
			}
		}
	}

	result = LS[0];

	// �ഫ�Ψ禡.
	auto&& raw_to_imgraw = [](const Raw& src){
		ImgRaw dst(src.RGB, src.getCol(), src.getRow(), 24);
		return dst;
	};
	auto&& imgraw_to_raw = [](const ImgRaw& src){
		Raw dst(src.width, src.height);
		dst.RGB = src; // �o�̷|�I�s�����禡��uch
		return dst;
	};
	auto&& Blend_to_imgraw = [](const Blend_Image& src){
		vector<unsigned char> temp(src.width*src.height*3);
		for(size_t i = 0; i < temp.size(); i++){
			temp[i]=src.RGB[i];
		}
		ImgRaw dst(temp, src.width, src.height, 24);
		return dst;
	};

	//Blend_to_imgraw(LS[0]).bmp("overlap.bmp");

	blendImg(limg, result, dx, dy, LEFT, bol, bor);
	blendImg(rimg, result, dx, dy, RIGHT, bol, bor);
}

/* �o�өγ\�i�H�ݩԴ������r���򰵪�
vector<unsigned char> MultiBandBlending(vector<unsigned char> left, vector<unsigned char> right, int width, int height)
{
	int Col = width;
	int Row = height;
	vector<unsigned char> I(Col * Row * 3, 0);
	vector<bool> R(Col * Row, false);
	int o = 0, i = 0, j = 0;
	//----------------------------------------
	// bicubic ���v�x�} X.Y.Z
	float transMat_2[3][3] = { { 2.0f, 0.f, 0.f },{ 0.f, 2.0f, 0.f },{ 0.f, 0.f, 1.f } };
	float transMat_0_5[3][3] = { { 0.5f, 0.f, 0.f },{ 0.f, 0.5f, 0.f },{ 0.f, 0.f, 1.f } };
	//----------------------------------------
	// ���X���|�϶�
	vector<struct X_S> M(Row);
	vector<struct X_S> U(Col);
	int L_x = INT_MAX, R_x = INT_MIN;
	int U_y = INT_MAX, D_y = INT_MIN;
	for (i = 0; i < Row; ++i)
	{
		for (j = 0; j < Col; j++)
		{
			if ((left[(i * Col + j) * 3 + 0] != 0 || left[(i * Col + j) * 3 + 1] != 0 || left[(i * Col + j) * 3 + 2] != 0) &&
				(right[(i * Col + j) * 3 + 0] != 0 || right[(i * Col + j) * 3 + 1] != 0 || right[(i * Col + j) * 3 + 2] != 0))
			{
				right[(i * Col + j) * 3 + 0] = 0;
				right[(i * Col + j) * 3 + 1] = 0;
				right[(i * Col + j) * 3 + 2] = 0;
				M[i].m = j + 1;
				break;
			}
		}
		for (j = Col - 1; j >= 0; j--)
		{
			if ((left[(i * Col + j) * 3 + 0] != 0 || left[(i * Col + j) * 3 + 1] != 0 || left[(i * Col + j) * 3 + 2] != 0) &&
				(right[(i * Col + j) * 3 + 0] != 0 || right[(i * Col + j) * 3 + 1] != 0 || right[(i * Col + j) * 3 + 2] != 0))
			{
				left[(i * Col + j) * 3 + 0] = 0;
				left[(i * Col + j) * 3 + 1] = 0;
				left[(i * Col + j) * 3 + 2] = 0;
				M[i].n = j - 1;
				break;
			}
		}
		if (M[i].m == 0 && M[i].n == 0)
		{
			M[i].l = 0;
			M[i].middle = 0;
		}
		else
		{
			M[i].middle = (M[i].m + M[i].n) / 2;
			M[i].l = M[i].n - M[i].m + 1;
			if (M[i].m < L_x)
			{
				L_x = M[i].m;
			}
			if (M[i].n > R_x)
			{
				R_x = M[i].n;
			}
			if (i < U_y)
			{
				U_y = i;
			}
			D_y = i;
		}
	}
	for (j = 0; j < Col; j++)
	{
		for (i = 0; i < Row; ++i)
		{
			if (left[(i * Col + j) * 3 + 0] != 0 || left[(i * Col + j) * 3 + 1] != 0 || left[(i * Col + j) * 3 + 2] != 0)
			{
				left[(i * Col + j) * 3 + 0] = 0;
				left[(i * Col + j) * 3 + 1] = 0;
				left[(i * Col + j) * 3 + 2] = 0;
				break;
			}
		}
		for (i = 0; i < Row; ++i)
		{
			if (right[(i * Col + j) * 3 + 0] != 0 || right[(i * Col + j) * 3 + 1] != 0 || right[(i * Col + j) * 3 + 2] != 0)
			{
				right[(i * Col + j) * 3 + 0] = 0;
				right[(i * Col + j) * 3 + 1] = 0;
				right[(i * Col + j) * 3 + 2] = 0;
				break;
			}
		}
		for (i = Row - 1; i >= 0; --i)
		{
			if (left[(i * Col + j) * 3 + 0] != 0 || left[(i * Col + j) * 3 + 1] != 0 || left[(i * Col + j) * 3 + 2] != 0)
			{
				left[(i * Col + j) * 3 + 0] = 0;
				left[(i * Col + j) * 3 + 1] = 0;
				left[(i * Col + j) * 3 + 2] = 0;
				break;
			}
		}
		for (i = Row - 1; i >= 0; --i)
		{
			if (right[(i * Col + j) * 3 + 0] != 0 || right[(i * Col + j) * 3 + 1] != 0 || right[(i * Col + j) * 3 + 2] != 0)
			{
				right[(i * Col + j) * 3 + 0] = 0;
				right[(i * Col + j) * 3 + 1] = 0;
				right[(i * Col + j) * 3 + 2] = 0;
				break;
			}
		}
	}
	//----------------------------------------
	// �̷Ӷ��h�����Y�p�A�]��size�i�ର�_�Ʒ|�y���~�t�A�ҥH�N���|�ϰ����j�� octave^2
	// �ťզa��ɳ̫�@�I����
	int size_a = (int)pow(2.f, PYR_OCTAVE);
	int re_Col = (R_x - L_x + 1) % size_a == 0 ? (R_x - L_x + 1) : ((R_x - L_x + 1) / size_a + 1) * size_a;
	//int re_Row = Row % size_a == 0 ? Row : (Row / size_a + 1) * size_a;
	int re_Row = (D_y - U_y + 1) % size_a == 0 ? (D_y - U_y + 1) : ((D_y - U_y + 1) / size_a + 1) * size_a;
	//----------------------------------------
	Blend_Image left_t, right_t;
	//----------------------------------------
	//left
	left_t.width = re_Col;
	left_t.height = re_Row;
	left_t.RGB = new float[left_t.width * left_t.height * 3];
	//----------------------------------------
	//right
	right_t.width = re_Col;
	right_t.height = re_Row;
	right_t.RGB = new float[right_t.width * right_t.height * 3];
	//----------------------------------------
	//�k0
	for (i = 0; i < re_Col * re_Row; ++i)
	{
		left_t.RGB[i * 3 + 0] = 0.f;
		left_t.RGB[i * 3 + 1] = 0.f;
		left_t.RGB[i * 3 + 2] = 0.f;
		right_t.RGB[i * 3 + 0] = 0.f;
		right_t.RGB[i * 3 + 1] = 0.f;
		right_t.RGB[i * 3 + 2] = 0.f;
	}
	//----------------------------------------
	// ���X���|���a��
	for (i = U_y; i < U_y + re_Row; i++)
	{
		for (j = L_x; j < L_x + re_Col; j++)
		{
			left_t.RGB[((i - U_y) * re_Col + (j - L_x)) * 3 + 0] = (float)left[(i * Col + j) * 3 + 0];
			left_t.RGB[((i - U_y) * re_Col + (j - L_x)) * 3 + 1] = (float)left[(i * Col + j) * 3 + 1];
			left_t.RGB[((i - U_y) * re_Col + (j - L_x)) * 3 + 2] = (float)left[(i * Col + j) * 3 + 2];

			right_t.RGB[((i - U_y) * re_Col + (j - L_x)) * 3 + 0] = (float)right[(i * Col + j) * 3 + 0];
			right_t.RGB[((i - U_y) * re_Col + (j - L_x)) * 3 + 1] = (float)right[(i * Col + j) * 3 + 1];
			right_t.RGB[((i - U_y) * re_Col + (j - L_x)) * 3 + 2] = (float)right[(i * Col + j) * 3 + 2];
		}
	}
	int ur = 0, dr = 0;
	int set_r = 1;

	for (j = 0; j < re_Col; j++)
	{
		for (i = 0; i < re_Row; i++)
		{
			if (left_t.RGB[(i * re_Col + j) * 3 + 0] != 0.f || left_t.RGB[(i * re_Col + j) * 3 + 1] != 0.f || left_t.RGB[(i * re_Col + j) * 3 + 2] != 0.f)
			{
				ur = i;
				break;
			}
		}
		for (i = ur; i >= 0; i--)
		{
			left_t.RGB[(i * re_Col + j) * 3 + 0] = left_t.RGB[((ur + set_r) * re_Col + j) * 3 + 0];
			left_t.RGB[(i * re_Col + j) * 3 + 1] = left_t.RGB[((ur + set_r) * re_Col + j) * 3 + 1];
			left_t.RGB[(i * re_Col + j) * 3 + 2] = left_t.RGB[((ur + set_r) * re_Col + j) * 3 + 2];
		}
		for (i = re_Row - 1; i >= 0; i--)
		{
			if (left_t.RGB[(i * re_Col + j) * 3 + 0] != 0.f || left_t.RGB[(i * re_Col + j) * 3 + 1] != 0.f || left_t.RGB[(i * re_Col + j) * 3 + 2] != 0.f)
			{
				dr = i;
				break;
			}
		}
		for (i = dr; i < re_Row; i++)
		{
			left_t.RGB[(i * re_Col + j) * 3 + 0] = left_t.RGB[((dr - set_r) * re_Col + j) * 3 + 0];
			left_t.RGB[(i * re_Col + j) * 3 + 1] = left_t.RGB[((dr - set_r) * re_Col + j) * 3 + 1];
			left_t.RGB[(i * re_Col + j) * 3 + 2] = left_t.RGB[((dr - set_r) * re_Col + j) * 3 + 2];
		}
		//-------------------------------------------------------------------------------------------------------------------------------------------------
		for (i = 0; i < re_Row; i++)
		{
			if (right_t.RGB[(i * re_Col + j) * 3 + 0] != 0.f || right_t.RGB[(i * re_Col + j) * 3 + 1] != 0.f || right_t.RGB[(i * re_Col + j) * 3 + 2] != 0.f)
			{
				ur = i;
				break;
			}
		}
		for (i = ur; i >= 0; i--)
		{
			right_t.RGB[(i * re_Col + j) * 3 + 0] = right_t.RGB[((ur + set_r) * re_Col + j) * 3 + 0];
			right_t.RGB[(i * re_Col + j) * 3 + 1] = right_t.RGB[((ur + set_r) * re_Col + j) * 3 + 1];
			right_t.RGB[(i * re_Col + j) * 3 + 2] = right_t.RGB[((ur + set_r) * re_Col + j) * 3 + 2];
		}
		for (i = re_Row - 1; i >= 0; i--)
		{
			if (right_t.RGB[(i * re_Col + j) * 3 + 0] != 0.f || right_t.RGB[(i * re_Col + j) * 3 + 1] != 0.f || right_t.RGB[(i * re_Col + j) * 3 + 2] != 0.f)
			{
				dr = i;
				break;
			}
		}
		for (i = dr; i < re_Row; i++)
		{
			right_t.RGB[(i * re_Col + j) * 3 + 0] = right_t.RGB[((dr - set_r) * re_Col + j) * 3 + 0];
			right_t.RGB[(i * re_Col + j) * 3 + 1] = right_t.RGB[((dr - set_r) * re_Col + j) * 3 + 1];
			right_t.RGB[(i * re_Col + j) * 3 + 2] = right_t.RGB[((dr - set_r) * re_Col + j) * 3 + 2];
		}
	}
	for (i = 0; i < re_Row; i++)
	{
		for (j = 0; j < re_Col; j++)
		{
			if (left_t.RGB[(i * re_Col + j) * 3 + 0] != 0.f || left_t.RGB[(i * re_Col + j) * 3 + 1] != 0.f || left_t.RGB[(i * re_Col + j) * 3 + 2] != 0.f)
			{
				ur = j;
				break;
			}
		}
		for (j = ur; j >= 0; j--)
		{
			left_t.RGB[(i * re_Col + j) * 3 + 0] = left_t.RGB[(i * re_Col + (ur + set_r)) * 3 + 0];
			left_t.RGB[(i * re_Col + j) * 3 + 1] = left_t.RGB[(i * re_Col + (ur + set_r)) * 3 + 1];
			left_t.RGB[(i * re_Col + j) * 3 + 2] = left_t.RGB[(i * re_Col + (ur + set_r)) * 3 + 2];
		}
		for (j = re_Col - 1; j >= 0; j--)
		{
			if (left_t.RGB[(i * re_Col + j) * 3 + 0] != 0.f || left_t.RGB[(i * re_Col + j) * 3 + 1] != 0.f || left_t.RGB[(i * re_Col + j) * 3 + 2] != 0.f)
			{
				dr = j;
				break;
			}
		}
		for (j = dr; j < re_Col; j++)
		{
			left_t.RGB[(i * re_Col + j) * 3 + 0] = left_t.RGB[(i * re_Col + (dr - set_r)) * 3 + 0];
			left_t.RGB[(i * re_Col + j) * 3 + 1] = left_t.RGB[(i * re_Col + (dr - set_r)) * 3 + 1];
			left_t.RGB[(i * re_Col + j) * 3 + 2] = left_t.RGB[(i * re_Col + (dr - set_r)) * 3 + 2];
		}
		//-------------------------------------------------------------------------------------------------------------------------------------------------
		for (j = 0; j < re_Col; j++)
		{
			if (right_t.RGB[(i * re_Col + j) * 3 + 0] != 0.f || right_t.RGB[(i * re_Col + j) * 3 + 1] != 0.f || right_t.RGB[(i * re_Col + j) * 3 + 2] != 0.f)
			{
				ur = j;
				break;
			}
		}
		for (j = ur; j >= 0; j--)
		{
			right_t.RGB[(i * re_Col + j) * 3 + 0] = right_t.RGB[(i * re_Col + (ur + set_r)) * 3 + 0];
			right_t.RGB[(i * re_Col + j) * 3 + 1] = right_t.RGB[(i * re_Col + (ur + set_r)) * 3 + 1];
			right_t.RGB[(i * re_Col + j) * 3 + 2] = right_t.RGB[(i * re_Col + (ur + set_r)) * 3 + 2];
		}
		for (j = re_Col - 1; j >= 0; j--)
		{
			if (right_t.RGB[(i * re_Col + j) * 3 + 0] != 0.f || right_t.RGB[(i * re_Col + j) * 3 + 1] != 0.f || right_t.RGB[(i * re_Col + j) * 3 + 2] != 0.f)
			{
				dr = j;
				break;
			}
		}
		for (j = dr; j < re_Col; j++)
		{
			right_t.RGB[(i * re_Col + j) * 3 + 0] = right_t.RGB[(i * re_Col + (dr - set_r)) * 3 + 0];
			right_t.RGB[(i * re_Col + j) * 3 + 1] = right_t.RGB[(i * re_Col + (dr - set_r)) * 3 + 1];
			right_t.RGB[(i * re_Col + j) * 3 + 2] = right_t.RGB[(i * re_Col + (dr - set_r)) * 3 + 2];
		}
	}

	// �إ߰����h
	vector<Blend_Image> L_G(PYR_OCTAVE), R_G(PYR_OCTAVE);
	//----------------------------------------
	// ��0�h�����
	L_G[0] = left_t;
	R_G[0] = right_t;
	//----------------------------------------
	// ��l�ϼh���A��0�h�������Y�p�⭿
	float sigma = sqrt((-1) * pow(2.f, 2.f) / (2.f * log(0.5f)));
	for (i = 1; i < PYR_OCTAVE; i++)
	{
		L_G[i] = Bicubic(BlurImage(L_G[i - 1], sigma, PYR_R), transMat_0_5);
		R_G[i] = Bicubic(BlurImage(R_G[i - 1], sigma, PYR_R), transMat_0_5);
	}

	//----------------------------------------
	// �إ� L �h
	vector<Blend_Image> L_L(PYR_OCTAVE), R_L(PYR_OCTAVE);
	//----------------------------------------
	// �̰��h���m���T�A��l���h����ɸ�T
	for (i = 0; i < PYR_OCTAVE - 1; i++)
	{
		L_L[i] = sub(L_G[i], Bicubic(L_G[i + 1], transMat_2));
		R_L[i] = sub(R_G[i], Bicubic(R_G[i + 1], transMat_2));
	}
	L_L[PYR_OCTAVE - 1] = L_G[PYR_OCTAVE - 1];
	R_L[PYR_OCTAVE - 1] = R_G[PYR_OCTAVE - 1];
	//----------------------------------------
	// �N��ӭ��|�϶��� L �h�@���X�A��ɸ�T�q�����I���}�U�e�@�b
	// ��m��T�̷Ӱ����v�������t�A���� kernel ���e�פj�p
	vector<Blend_Image> L_T(PYR_OCTAVE);
	vector<float> kernel = getGaussianKernel(L_L[PYR_OCTAVE - 1].width);
	int middle = 0;

	for (o = 0; o < PYR_OCTAVE; o++)
	{
		L_T[o].width = L_L[o].width;
		L_T[o].height = L_L[o].height;
		L_T[o].RGB = new float[L_T[o].width * L_T[o].height * 3];
		middle = L_T[o].width % 2 == 0 ? L_T[o].width / 2 + 1 : L_T[o].width / 2;
		for (i = 0; i < L_T[o].height; ++i)
		{
			for (j = 0; j < L_T[o].width; ++j)
			{
				if (o == PYR_OCTAVE - 1)
				{
					L_T[o].RGB[(i * L_T[o].width + j) * 3 + 0] = kernel[j] * L_L[o].RGB[(i * L_L[o].width + j) * 3 + 0] + (1.f - kernel[j]) * R_L[o].RGB[(i * R_L[o].width + j) * 3 + 0];
					L_T[o].RGB[(i * L_T[o].width + j) * 3 + 1] = kernel[j] * L_L[o].RGB[(i * L_L[o].width + j) * 3 + 1] + (1.f - kernel[j]) * R_L[o].RGB[(i * R_L[o].width + j) * 3 + 1];
					L_T[o].RGB[(i * L_T[o].width + j) * 3 + 2] = kernel[j] * L_L[o].RGB[(i * L_L[o].width + j) * 3 + 2] + (1.f - kernel[j]) * R_L[o].RGB[(i * R_L[o].width + j) * 3 + 2];
				}
				else
				{
					if (j < middle)
					{
						L_T[o].RGB[(i * L_T[o].width + j) * 3 + 0] = L_L[o].RGB[(i * L_L[o].width + j) * 3 + 0];
						L_T[o].RGB[(i * L_T[o].width + j) * 3 + 1] = L_L[o].RGB[(i * L_L[o].width + j) * 3 + 1];
						L_T[o].RGB[(i * L_T[o].width + j) * 3 + 2] = L_L[o].RGB[(i * L_L[o].width + j) * 3 + 2];
					}
					else if (j == middle)
					{
						L_T[o].RGB[(i * L_T[o].width + j) * 3 + 0] = (L_L[o].RGB[(i * L_L[o].width + j) * 3 + 0] + R_L[o].RGB[(i * R_L[o].width + j) * 3 + 0]) / 2;
						L_T[o].RGB[(i * L_T[o].width + j) * 3 + 1] = (L_L[o].RGB[(i * L_L[o].width + j) * 3 + 1] + R_L[o].RGB[(i * R_L[o].width + j) * 3 + 1]) / 2;
						L_T[o].RGB[(i * L_T[o].width + j) * 3 + 2] = (L_L[o].RGB[(i * L_L[o].width + j) * 3 + 2] + R_L[o].RGB[(i * R_L[o].width + j) * 3 + 2]) / 2;
					}
					else
					{
						L_T[o].RGB[(i * L_T[o].width + j) * 3 + 0] = R_L[o].RGB[(i * R_L[o].width + j) * 3 + 0];
						L_T[o].RGB[(i * L_T[o].width + j) * 3 + 1] = R_L[o].RGB[(i * R_L[o].width + j) * 3 + 1];
						L_T[o].RGB[(i * L_T[o].width + j) * 3 + 2] = R_L[o].RGB[(i * R_L[o].width + j) * 3 + 2];
					}
				}
			}
		}
	}

	//----------------------------------------
	// �٭���r��A�q�̰��h����m��T�̧ǩ�j�֥[��ɸ�T�A�̩��h���V�X�L�᪺���G
	for (o = PYR_OCTAVE - 2; o >= 0; o--)
	{
		L_T[o] = add(L_T[o], Bicubic(L_T[o + 1], transMat_2));
	}
	//----------------------------------------
	// ��X�̫ᵲ�G�A�N�����|���϶��ɦ^
	for (i = 0; i < Row; i++)
	{
		for (j = 0; j < Col; j++)
		{
			if ((left[(i * Col + j) * 3 + 0] == 0 && left[(i * Col + j) * 3 + 1] == 0 && left[(i * Col + j) * 3 + 2] == 0) ||
				(right[(i * Col + j) * 3 + 0] == 0 && right[(i * Col + j) * 3 + 1] == 0 && right[(i * Col + j) * 3 + 2] == 0))
			{
				I[(i * Col + j) * 3 + 0] = left[(i * Col + j) * 3 + 0] + right[(i * Col + j) * 3 + 0];
				I[(i * Col + j) * 3 + 1] = left[(i * Col + j) * 3 + 1] + right[(i * Col + j) * 3 + 1];
				I[(i * Col + j) * 3 + 2] = left[(i * Col + j) * 3 + 2] + right[(i * Col + j) * 3 + 2];
			}
			else
			{
				I[(i * Col + j) * 3 + 0] = 255;
				I[(i * Col + j) * 3 + 1] = 255;
				I[(i * Col + j) * 3 + 2] = 255;
			}
		}
	}

	for (i = U_y; i < U_y + re_Row; i++)
	{
		for (j = L_x; j < L_x + re_Col; j++)
		{
			if (I[(i * Col + j) * 3 + 0] == 255 && I[(i * Col + j) * 3 + 1] == 255 && I[(i * Col + j) * 3 + 2] == 255)
			{
				I[(i * Col + j) * 3 + 0] = (unsigned char)L_T[0].RGB[((i - U_y) * re_Col + (j - L_x)) * 3 + 0];
				I[(i * Col + j) * 3 + 1] = (unsigned char)L_T[0].RGB[((i - U_y) * re_Col + (j - L_x)) * 3 + 1];
				I[(i * Col + j) * 3 + 2] = (unsigned char)L_T[0].RGB[((i - U_y) * re_Col + (j - L_x)) * 3 + 2];
			}
		}
	}
	//RawToBmp("test", I, Col, Row);
	//system("test.bmp");
	//----------------------------------------
	// ����O����
	for (o = 0; o < PYR_OCTAVE; o++)
	{
		delete[] L_G[o].RGB;
		delete[] R_G[o].RGB;
		if (o != PYR_OCTAVE - 1)
		{
			delete[] L_L[o].RGB;
			delete[] R_L[o].RGB;
		}
		delete[] L_T[o].RGB;
	}
	L_G.clear();
	R_G.clear();
	L_L.clear();
	R_L.clear();
	L_T.clear();
	return I;
}
*/