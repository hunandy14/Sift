/*****************************************************************
Name :
Date : 2017/07/04
By   : CharlotteHonG
Final: 2017/07/04
*****************************************************************/
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
    if(&img_gau == &img_ori) {
        throw file_same("## Erroe! in and out is same.");
    }
    // �]�w���T���j�p
    vector<unsigned char> img_gauX(img_ori.size());
    img_gau.resize(img_ori.size());
    // �����x�}�P�b�|
    vector<float> gau_mat = gau_matrix(p);
    const int r = gau_mat.size()/2;
    // �����ҽk X �b
    for(unsigned j = 0; j < height; ++j) {
        for(unsigned i = 0; i < width; ++i) {
            size_t sum = 0;
            for(unsigned k = 0; k < gau_mat.size(); ++k) {
                int idx = (i-r+k);
                if(idx < 0) { idx=0; }
                sum += img_ori[j*width+idx]*gau_mat[k];
            }
            img_gauX[j*width+i] = sum;
        }
    }
    // �����ҽk Y �b
    for(unsigned j = 0; j < height; ++j) {
        for(unsigned i = 0; i < width; ++i) {
            size_t sum = 0;
            for(unsigned k = 0; k < gau_mat.size(); ++k) {
                int idx = (j-r+k);
                if(idx < 0) { idx=0; }
                sum += img_gauX[i*height+idx]*gau_mat[k];
            }
            img_gau[i*height+j] = sum;
        }
    }
}
// ��������
float GauBlur::gau_meth(size_t r, float p) {
    float two = 2;
    float num = exp(-pow(r, two) / (two*pow(p, two)));
    num /= sqrt(two*M_PI)*p;
    return num;
}
// �����x�}
vector<float> GauBlur::gau_matrix(float p){
    vector<float> gau_mat;
    // �p��x�}����
    int mat_len = ((p-0.8) / 0.3+1.0) * 2.0;
    if (mat_len % 2 == 0) {++mat_len;}
    // �@�������x�}
    gau_mat.resize(mat_len);
    float sum=0;
    for(int i=0, j=mat_len/2; j < mat_len; ++i, ++j) {
        float temp;
        if(i) {
            temp = gau_meth(i);
            gau_mat[j] = temp;
            gau_mat[mat_len-j-1] = temp;
            sum += temp += temp;
        } else {
            gau_mat[j]=gau_meth(i);
            sum += gau_mat[j];
        }
    }
    // �k�@��
    for(auto&& i : gau_mat) { i /= sum; }
    return gau_mat;
}
//----------------------------------------------------------------
// ZroOrder�վ�j�p
void Scaling::zero(vector<unsigned char>& img,
    vector<unsigned char>& img_ori, size_t width,
    size_t height, float Ratio)
{
    int w=floor(width * Ratio);
    int h=floor(height * Ratio);
    img.resize(w*h);
    for(int j = 0; j < h; ++j) {
        for(int i = 0; i < w; ++i) {
            img[j*w+i] = img_ori[(int)(j/Ratio)*width + (int)(i/Ratio)];
        }
    }
}
// FisrtOrder�վ�j�p
void Scaling::first(vector<unsigned char>& img,
    vector<unsigned char>& img_ori, size_t width,
    size_t height, float Ratio)
{
    int w=floor(width * Ratio);
    int h=floor(height * Ratio);
    img.resize(h*w);

    unsigned char A, B, C, D;// ���񪺥|���I
    unsigned char AB, CD, X;
    int oy, ox;// �������Ϫ��y��
    int a, b;// ������ a �P b

    for(int j=0; j < h; ++j) {
        for(int i=0; i < w; ++i) {
            oy = floor(j/Ratio);
            ox = floor(i/Ratio);

            A = img_ori[oy*width + ox];
            B = img_ori[oy*width + ox+1];
            C = img_ori[oy+1*width + ox];
            D = img_ori[oy+1*width + ox+1];

            a = (i-ox*Ratio)/(Ratio);
            b = (j-oy*Ratio)/(Ratio);
            AB = (A*(1.0-a)) + (B*a);
            CD = (C*(1.0-a)) + (D*a);
            X = ((AB*(1.0-b)) + (CD*b));

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
    auto getMask = [&](size_t oy, size_t ox){
        uch** mask = new uch*[4];
        for (unsigned i = 0; i < 4; ++i)
            mask[i] = new uch[4];
        // ���o�P��16�I
        int foy,fox; // �״_�᪺��l�y��
        for (int j = 0; j < 4; ++j){
            for (int i = 0; i < 4; ++i){
                foy=oy+(j-1); fox=ox+(i-1);
                // �W�L����ɭ״_
                if (foy<0){foy=1;}
                // �W�L�W��ɭ״_
                if (fox<0){fox=1;}
                // �W�L�U��ɭ״_
                if(foy==(int)height){foy-=2;}
                if(foy==(int)height-1){foy-=1;}
                // �W�L�k��ɭ״_
                if (fox==(int)width){fox-=2;}
                if (fox==(int)width-1){fox-=1;}
                // ��������������
                mask[j][i] = img_ori[foy*width + fox];
            }
        }
        // ����O����
        // for (int i = 0; i < 4; ++i)
        //     delete [] mask[i];
        // delete [] mask;
        return mask;
    };

    int w=floor(width * Ratio);
    int h=floor(height * Ratio);
    img.resize(h*w);
    int oy, ox; // �������Ϫ��y��
    double a, b;// ���J����Ҧ�m
    uch** mask;// �B�n(�P��16�I)
    uch X;     // �Ȧs
    for(int j = 0; j < h; ++j) {
        for(int i = 0; i < w; ++i) {
            oy = floor(j/Ratio);
            ox = floor(i/Ratio);
            a = (i-ox*Ratio)/(Ratio);
            b = (j-oy*Ratio)/(Ratio);
            // ���o�P��16�I
            mask = getMask(oy, ox);
            // �ɤJ�P��16�I�P���J����Ҧ�m
            X = bicubicInterpolate(mask, b, a);
            // �g�J�Ȧs��
            img[j*w + i] = X;
            // ����O����
            // for (int i = 0; i < 4; ++i)
            //     delete [] mask[i];
            // delete [] mask;
        }
    }
}
//----------------------------------------------------------------
