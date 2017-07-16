/*****************************************************************
Name :
Date : 2017/07/04
By   : CharlotteHonG
Final: 2017/07/04
*****************************************************************/
#pragma warning(disable : 4819)
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

#include "imglib.hpp"
constexpr auto M_PI = 3.14159265358979323846;

// �����ҽk
void GauBlur::raw2GauBlur(vector<unsigned char>& img_gau,
	vector<unsigned char>& img_ori,
	size_t width, size_t height, float p)
{
	// �ӷ��ۦP�ҥ~
	if (&img_gau == &img_ori) {
		throw file_same("## Erroe! in and out is same.");
	}
	// �]�w���T���j�p
	img_gau.resize(img_ori.size());
	// �w�s
	vector<double> img_gauX(img_ori.size());
	// �����x�}�P�b�|
	vector<float> gau_mat = gau_matrix(p);
	const size_t r = gau_mat.size() / 2;
	// �����ҽk X �b
	for (unsigned j = 0; j < height; ++j) {
		for (unsigned i = 0; i < width; ++i) {
			double sum = 0;
			for (unsigned k = 0; k < gau_mat.size(); ++k) {
				int idx = (int)(i - r + k);
				if (idx < 0) { idx = 0; }
				else if (idx >(int)(width - 1)) { idx = (int)(width-1); }
				sum += (float)img_ori[j*width + idx] * gau_mat[k];
			}
			img_gauX[j*width + i] = sum;
		}
	}
	// �����ҽk Y �b
	for (unsigned j = 0; j < height; ++j) {
		for (unsigned i = 0; i < width; ++i) {
			double sum = 0;
			for (unsigned k = 0; k < gau_mat.size(); ++k) {
				int idx = (int)(j-r+k);
				if (idx < 0) { idx = 0; }
				else if (idx > (int)(width-1)) { idx = (int)(width-1); }
				sum += (float)img_gauX[i*height + idx] * gau_mat[k];
			}
			img_gau[i*height + j] = (unsigned char)sum;
		}
	}
}
// ��������
float GauBlur::gau_meth(size_t r, float p) {
	double two = 2.0;
	double num = exp(-pow(r, two) / (two*pow(p, two)));
	num /= sqrt(two*M_PI)*p;
	return (float)num;
}
// �����x�}
vector<float> GauBlur::gau_matrix(float p) {
	vector<float> gau_mat;
	// �p��x�}����
	int mat_len = (int)(((p - 0.8) / 0.3 + 1.0) * 2.0);
	if (mat_len % 2 == 0) { ++mat_len; }
	// �@�������x�}
	gau_mat.resize(mat_len);
	float sum = 0;
	for (int i = 0, j = mat_len / 2; j < mat_len; ++i, ++j) {
		float temp;
		if (i) {
			temp = gau_meth(i, p);
			gau_mat[j] = temp;
			gau_mat[mat_len - j - 1] = temp;
			sum += temp += temp;
		}
		else {
			gau_mat[j] = gau_meth(i, p);
			sum += gau_mat[j];
		}
	}
	// �k�@��
	for (auto&& i : gau_mat) { i /= sum; }
	return gau_mat;
}
//----------------------------------------------------------------
// ZroOrder�վ�j�p
void Scaling::zero(vector<types>& img,
	vector<types>& img_ori, size_t width,
	size_t height, float Ratio)
{
	int w = (int)floor(width * Ratio);
	int h = (int)floor(height * Ratio);
	img.resize(w*h);
	for (int j = 0; j < h; ++j) {
		for (int i = 0; i < w; ++i) {
			img[j*w + i] =
				img_ori[(int)(j/Ratio)*width + (int)(i/Ratio)];
		}
	}
}
// FisrtOrder�վ�j�p
void Scaling::first(vector<types>& img,
	vector<types>& img_ori, size_t width,
	size_t height, float Ratio)
{
	using uch = types;
	int w = (int)floor(width * Ratio);
	int h = (int)floor(height * Ratio);
	img.resize(h*w);
	for (int j = 0; j < h; ++j) {
		for (int i = 0; i < w; ++i) {
			// �������Ϫ��y��
			int oy = (int)floor(j / Ratio);
			int ox = (int)floor(i / Ratio);
			// ���񪺥|���I
			size_t xp = (ox+1) > (int)(width-1)? width-1: (ox+1);
			size_t yp = (oy+1) > (int)(height-1)? height-1: (oy+1);
			uch A = img_ori[oy*width + ox];
			uch B = img_ori[oy*width + xp];
			uch C = img_ori[yp*width + ox];
			uch D = img_ori[yp*width + xp];
			// ������ a �P b
			float a = (i - ox*Ratio) / (Ratio);
			float b = (j - oy*Ratio) / (Ratio);
			uch AB = (uch)(A*(1.0 - a)) + (B*a);
			uch CD = (uch)(C*(1.0 - a)) + (D*a);
			uch X = (uch)((AB*(1.0 - b)) + (CD*b));
			img[j*w + i] = X;
		}
	}
}
// Bicubic�վ�j�p
void Scaling::cubic(vector<unsigned char>& img,
	vector<unsigned char>& img_ori, size_t width,
	size_t height, float Ratio)
{
	using uch = unsigned char;
	// Bicubic ���o�P��16�I
	auto getMask = [&](uch* mask, size_t oy, size_t ox) {
		// ���o�P��16�I
		int foy, fox; // �״_�᪺��l�y��
		for (int j = 0, idx = 0; j < 4; ++j) {
			for (int i = 0; i < 4; ++i, ++idx) {
				foy = (int)(oy+(j-1)), fox = (int)(ox+(i-1));
				// �W�L����ɭ״_
				if (foy<0) { foy = 1; }
				// �W�L�W��ɭ״_
				if (fox<0) { fox = 1; }
				// �W�L�U��ɭ״_
				if (foy == (int)height) { foy -= 2; }
				if (foy == (int)height - 1) { foy -= 1; }
				// �W�L�k��ɭ״_
				if (fox == (int)width) { fox -= 2; }
				if (fox == (int)width - 1) { fox -= 1; }
				// ��������������
				mask[idx] = img_ori[foy*width + fox];
			}
		}
	};

	int w = (int)floor(width * Ratio);
	int h = (int)floor(height * Ratio);
	img.resize(h*w);
	for (int j = 0; j < h; ++j) {
		for (int i = 0; i < w; ++i) {
			// �������Ϫ��y��
			int oy = (int)floor(j / Ratio);
			int ox = (int)floor(i / Ratio);
			double a = (i - ox*Ratio) / (Ratio);
			double b = (j - oy*Ratio) / (Ratio);
			// ���o�P��16�I
			uch mask[16];
			getMask(mask, oy, ox);
			// �ɤJ�P��16�I�P���J����Ҧ�m
			uch X = bicubicInterpolate(mask, b, a);
			// �g�J�Ȧs��
			img[j*w + i] = X;
		}
	}
}
//----------------------------------------------------------------
