/*****************************************************************
Name :
Date : 2017/07/05
By   : CharlotteHonG
Final: 2017/07/05
*****************************************************************/
#pragma warning(disable : 4819)
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <cmath>

#include "Sift.hpp"
using namespace std;

#define M_PI 3.14159265358979323846

// (多一個 gray_tran 參數判定感覺不是很好有機會修掉)
ImgRaw::ImgRaw(string bmpname, bool gray_tran){
	vector<unsigned char> img;
	uint32_t width, height;
	uint16_t bits;
	// 讀取圖片
	Raw::read_bmp(img, bmpname, &width, &height, &bits);
	this->width    = width;
	this->height   = height;
	this->bitCount = bits;
	// 轉換灰階
	if (gray_tran == 1 and bits == 24) {
		Raw::raw2gray(img);
		this->bitCount = 8;
	}
	// 初始化(含正規化)
	raw_img.resize(img.size());
	for (size_t i = 0; i < img.size(); i++) {
		raw_img[i] = (float)img[i] / 255.0;
	}
}
// Sift建構子
Sift::Sift(ImgRaw img): raw_img(img) {
	FeatureStart = new Feature;
	FeatureEnd = FeatureStart;
}
// 合成圖片
void Sift::comp(vector<ImgRaw>& pyrs, string name) {
    // 合成圖片
    size_t w = pyrs[0].width, h = pyrs[0].height;
    size_t s2 = pyrs.size();
    ImgRaw big_img(w*s2, h);
    for (unsigned j = 0, idx = 0; j < h; ++j) {
        for (unsigned n = 0; n < s2; ++n) {
            for (unsigned i = 0; i < w; ++i) {
                big_img[idx++] = pyrs[n][j*w + i];
            }
        }
    }
    // 輸出圖片並開啟
    if (name != "") {
        name += ".bmp";
        big_img.bmp("pic/" + name, 8);
        // system(name.c_str());
    }
};
// 獲得特徵點直方圖統計
static inline float getMag(float dx, float dy){
	// 獲得強度
	return sqrtf(dx*dx + dy*dy);
}
static inline float getSita(float dx, float dy){
	// 獲得角度
	float s = atan2f(dy, dx);
	return (s > 0? s: (s+M_PI*2));
}
void Sift::AddnewFeaturestruct(int Inx, int Iny, float Insize, int kai, int sigmaOCT, float Inm, int Insita)
{
	Featureptr newnode = new Feature;
	newnode->x = Inx;
	newnode->y = Iny;
	newnode->mm = Inm;
	newnode->sita = Insita;
	newnode->size = Insize;
	newnode->kai = kai;
	newnode->sigmaOCT = sigmaOCT;
	newnode->nextptr = NULL;
	FeatureEnd->nextptr = newnode;
	FeatureEnd = newnode;
}
void Sift::getHistogramMS(const ImgRaw& doImage, float Insize, size_t scale, float sigma, 
	size_t Iny, size_t Inx, size_t InWidth, size_t Inr)
{
	const size_t matLen = Inr*2 + 1; // 遮罩長度
	vector<float> mag(matLen * matLen);  // 儲存強度
	vector<float> sita(matLen * matLen); // 儲存角度
	//做m、sita計算
	for (int j = (Iny - Inr), idx = 0; j <= (Iny + Inr); j++) {
		for (int i = (Inx - Inr); i <= (Inx + Inr); i++, idx++) {
			const float& dx = doImage[(j+0)*InWidth + (i+1)] - doImage[(j+0)*InWidth + (i-1)];
			const float& dy = doImage[(j+1)*InWidth + (i+0)] - doImage[(j-1)*InWidth + (i+0)];
			mag[idx] = getMag(dx, dy);
			sita[idx] = getSita(dx, dy) * 180.0/M_PI;
		}
	}
	// 高斯矩陣相乘
	vector<types> gauMat = Gaus::gau_matrix(sigma, matLen); // 一維高斯矩陣
	for (size_t j = 0; j < matLen; j++) {
		for (size_t i = 0; i < matLen; i++) {
			mag[j*matLen + i] *= gauMat [i];
		}
	}
	for (size_t j = 0; j < matLen; j++) {
		for (size_t i = 0; i < matLen; i++) {
			mag[j*matLen + i] *= gauMat [j];
		}
	}
	// 遮罩內的強度依照角度分類累加
	vector<float> magSum(36); // 強度累加陣列
	for (size_t i = 0; i < mag.size(); i++) {
		magSum[(sita[i]/10.0)] += mag[i];
	}
	// 新增主方向與副方向
	float maxMag = *max_element(magSum.begin(), magSum.end());
	for (size_t i = 0; i < 36; i++) {
		if (magSum[i] >= maxMag*0.8) {
			AddnewFeaturestruct(Inx, Iny, Insize, scale, sigma, magSum[i], i*10);
		}
	}
}
// 高斯金字塔
static inline size_t getPyramidOctv(size_t width, size_t height) {
	// 獲得金字塔八度的長
	int n = 1;
	for (float Limit = 1.0; Limit < width || Limit < height; ++n){
		Limit = powf(2.0, n);
	}
	if (n < 7) {n = 1;}
	else {n -= 6;}
	return n;
}
void Sift::pyramid2() {
	ImgRaw first_img(raw_img.width, raw_img.height);
	size_t octv = getPyramidOctv(raw_img.width, raw_img.height);
	float curr_size=2;
	for (size_t ph = 0; ph < octv; ph++, curr_size /= 2) { // 金字塔的sigma還沒修復第二層沒*2
		const size_t curr_Width = raw_img.width*curr_size;
		const size_t curr_Height = raw_img.height*curr_size;
		ImgRaw currSize_img(curr_Width ,curr_Height);
		//ImgRaw::first(currSize_img, raw_img, curr_size);
		ZoomInOut(currSize_img, curr_Width, curr_Height);

		// 高斯模糊
		vector<ImgRaw> gau_imgs(pyWidth);
		ImgRaw::gauBlur(gau_imgs[0], currSize_img, SIFT_GauSigma * (ph+1));
		for (size_t i = 1; i < pyWidth; i++) {
			const float curr_Sigma = SIFT_GauSigma * powf(sqrt(2.0),i) * (ph+1); // 這裡的2是 S+3 的 S
			gau_imgs[i].raw_img.resize(curr_Width*curr_Height);
			Gaus::GauBlur(gau_imgs[i], gau_imgs[i-1], curr_Width, curr_Height, curr_Sigma);
			//gau_imgs[i].bmp("_gau.bmp", 8);
		}
		// 高斯差分
		vector<ImgRaw> gauDog_imgs(pyWidth-1);
		for (size_t i = 0; i < pyWidth-1; i++) {
			gauDog_imgs[i].resize(curr_Width, curr_Height);
			for (size_t idx = 0; idx < curr_Width*curr_Height; idx++) {
				gauDog_imgs[i][idx] = gau_imgs[i][idx] - gau_imgs[i+1][idx];
			}
		}
		// 尋找特徵點並累加直方圖
		FeatureNow = FeatureEnd;
		for (size_t scale_idx = 1; scale_idx < pyWidth-1-1; scale_idx++) {
			ImgRaw MaxMin_img(curr_Width, curr_Height);
			ImgRaw Herris_img(curr_Width, curr_Height);
			// 特徵直方圖累加半徑
			const size_t r = curr_size * 3. * 1.5 * (scale_idx+1);
			for (size_t j = r+1; j < curr_Height-r-1; j++) {
				for (size_t i = r+1; i < curr_Width-r-1; i++) {
					// 尋找極值點
					if (findMaxMin(gauDog_imgs, scale_idx, curr_Width, j, i)) {
						const ImgRaw& currImg = gauDog_imgs[scale_idx];
						const float currSigma = SIFT_GauSigma * powf(sqrt(2.0),scale_idx) * (ph+1);
						if (Corner::harris(currImg, curr_Width, j, i, r)) {
							getHistogramMS(currImg, curr_size, scale_idx, currSigma, j, i, curr_Width, r);
							//Herris_img[j*curr_Width + i] = 255 /255.0;
						}
						//MaxMin_img[j*curr_Width + i] = 255 /255.0;
					}
				}
			}
			//MaxMin_img.bmp("MaxMin/MaxMin_"+to_string(scale_idx)+".bmp", 8);
			//Herris_img.bmp("Herris/Herris_"+to_string(scale_idx)+".bmp", 8);
		}
	}
}
// 畫線
void Draw::draw_line(ImgRaw& img, size_t y, size_t x, float line_len, float sg) {
	float value = 200 /255.0;
	float endvalue = 255 /255.0;

	if (line_len==1) {
		img.at2d(x, y) = value;
		return;
	}


	if (sg > 180) {
		sg -= 360.f;
	}
	sg*=-1; // 轉正圖片上下顛倒
	size_t x2 = x + line_len*cos(sg * M_PI/180);
	size_t y2 = y + line_len*sin(sg * M_PI/180);
	float m = tan(sg * M_PI/180);


	auto draw_line = [&](size_t yValue, size_t xValue) {
		if (yValue) {

		}
		//img.at2d(yValue, xValue) = value;
		if (yValue > 0 and xValue > 0 and 
			xValue < img.width-1 and yValue < img.height-1)
		{
			img.raw_img[yValue*img.width + xValue] = value;
		}
	};

	for (int i = 0; i < abs((int)(x-x2)); i++) {
		if (sg > 0 && abs(sg>90.f)) {
			int yValue = i*m + std::max(y, y2);
			int xValue = i+ std::min(x, x2);
			//img.at2d(yValue, xValue) = value;
			draw_line(yValue, xValue);
		}
		else if (sg < 0 && abs(sg<-90.f)) {
			int yValue = i*m + std::min(y, y2);
			int xValue = i+ std::min(x, x2);
			draw_line(yValue, xValue);
		}
		else if (sg > 0 && abs(sg<90.f)) {
			int yValue = i*m + std::min(y, y2);
			int xValue = i+ std::min(x, x2);
			draw_line(yValue, xValue);
		}
		else if (sg < 0 && abs(sg>-90.f)) {
			int yValue = i*m + std::max(y, y2);
			int xValue = i+ std::min(x, x2);
			draw_line(yValue, xValue);
		}
		else if (sg==0) {
			int yValue = y;
			int xValue = i+x;
			draw_line(yValue, xValue);
		}
	}
	// 垂直處理
	if (sg == 90 || sg == -90) {
		for (size_t i = 0; i < abs((int)(y-y2)); i++) {
			size_t yValue = i+std::min(y, y2);
			size_t xValue = x;
			draw_line(yValue, xValue);
		}
	}
	// 頭尾
	draw_line(y, x);
	draw_line(y2, x2);
}
// 找極值
bool Sift::findMaxMin(vector<ImgRaw>& gauDog_imgs, size_t scale_idx, size_t curr_Width, size_t y, size_t x) {
	const float& Val = gauDog_imgs[scale_idx][y*curr_Width + x];

	bool maxc, minc;
	maxc = minc = true;
	for (int k = (scale_idx - 1); k <= (scale_idx + 1); k++)
	{
		for (int j = (y - 1); j <= (y + 1); j++)
		{
			for (int i = (x - 1); i < (x + 1); i++)
			{
				if (gauDog_imgs[k][j*curr_Width + i] > Val)
					maxc = false;
				if (gauDog_imgs[k][j*curr_Width + i] < Val)
					minc = false;
				if (!(maxc | minc))
					return false;
			}
		}
	}

	if (abs(Val) < 0.015)
		return false;
	else
		return true;


	/*for (int i = -1; i < 2; i++) {
	if (fabs(Val) > SIFT_Dx/2.f) {
	if (
	(Val>0 and
	Val>=gauDog_imgs[scale_idx+i][(y-1)*curr_Width+x-1] and
	Val>=gauDog_imgs[scale_idx+i][(y-1)*curr_Width+x+0] and
	Val>=gauDog_imgs[scale_idx+i][(y-1)*curr_Width+x+1] and
	Val>=gauDog_imgs[scale_idx+i][(y+0)*curr_Width+x-1] and
	Val>=gauDog_imgs[scale_idx+i][(y+0)*curr_Width+x+0] and
	Val>=gauDog_imgs[scale_idx+i][(y+0)*curr_Width+x+1] and
	Val>=gauDog_imgs[scale_idx+i][(y+1)*curr_Width+x-1] and
	Val>=gauDog_imgs[scale_idx+i][(y+1)*curr_Width+x+0] and
	Val>=gauDog_imgs[scale_idx+i][(y+1)*curr_Width+x+1])
	){
	return true;
	}
	}
	}
	return false;*/
};
// 放大縮小
void Sift::ZoomInOut(ImgRaw& doImage, int InWidth, int InHeight)
{
	// 原圖資訊(取自資料成員)
	const size_t Height=raw_img.height;
	const size_t Width=raw_img.width;
	const vector<float>& gray=raw_img;

	float XX1, XX2, XXYY;
	float dx, dy;
	float Io, Jo;
	for (int j = 0; j < InHeight - 1; ++j)
	{
		Jo = j * (Height - 1) / (InHeight - 1);
		dy = 1 - ((j * (Height - 1) / (InHeight - 1)) - Jo);
		for (int i = 0; i < InWidth - 1; ++i)
		{
			Io = i * (Width - 1) / (InWidth - 1);
			dx = 1 - ((i * (Width - 1) / (InWidth - 1)) - Io);
			XX1 = gray[(int)((Io + 0) + (Jo + 0)*Width)] * dx + gray[(int)((Io + 1) + (Jo + 0)*Width)] * (1 - dx);
			XX2 = gray[(int)((Io + 0) + (Jo + 1)*Width)] * dx + gray[(int)((Io + 1) + (Jo + 1)*Width)] * (1 - dx);
			XXYY = XX1*dy + XX2*(1 - dy);
			doImage[j*InWidth + i] = XXYY;
		}
		XX1 = gray[(int)((Jo + 1)*Width - 1)];
		XX2 = gray[(int)((Jo + 2)*Width - 1)];
		XXYY = XX1*dy + XX2*(1 - dy);
		doImage[(j + 1)*InWidth - 1] = XXYY;
	}
	for (int i = 0; i < InWidth - 1; i++)
	{
		Io = i * (Width - 1) / (InWidth - 1);
		dx = 1 - ((i * (Width - 1) / (InWidth - 1)) - Io);
		XX1 = gray[(int)((Height - 1)*Width + (Io + 0))];
		XX2 = gray[(int)((Height - 1)*Width + (Io + 1))];
		XXYY = XX1*dx + XX2*(1 - dx);
		doImage[(InHeight - 1)*InWidth + i] = XXYY;
	}
	doImage[InHeight*InWidth - 1] = gray[Height*Width - 1];
}
// 印出箭頭
void Sift::addArrow()
{
	
	//箭頭倍率
	int mag = 100000.0;

	Featureptr pictureS = FeatureStart;
	size_t Height=raw_img.height;
	size_t Width=raw_img.width;

	
	
	// 臨時圖庫
	/*ImgRaw img(Width*3 * Height);
	img.width=Width;
	img.height=Height;*/

	ImgRaw img("kanna.bmp", 0);

	ImgRaw img_gray(Width, Height);

	//for (size_t j = 0; j < Height; j++) {
	//	for (size_t i = 0; i < Width; i++) {
	//		//img[(j*Width+i) + 0]=0.5;
	//		img[(j*Width+i)*3 + 0] = 0.5;
	//		img[(j*Width+i)*3 + 1] = 0;
	//		img[(j*Width+i)*3 + 2] = 0;
	//	}
	//}
	//img.bmp("feaArrow.bmp", 24);


	while (pictureS->nextptr != NULL)
	{
		//cout << pictureS->x  << ", "<< pictureS->y  << ", "<< pictureS->mag << ", " << pictureS->sita << endl;


		pictureS = pictureS->nextptr;
		int x, y;
		float roominout = pictureS->size;
		x = 1;

		/*
		//cout<< roominout << endl;
		x=pictureS->x/roominout;
		y=pictureS->y/roominout;
		//cout << x  << ", "<< y  << ", "<< roominout << ", " << pictureS->sita << endl;
		Draw::draw_line(img, y, x, 10, -pictureS->sita);
		*/

		
		//劃出箭頭的直線
		while (true)
		{
			
			if (pictureS->sita == 90 || pictureS->sita == 270) y = 0;
			else y = x * tan((pictureS->sita)*M_PI / 180.0);

			if ((x*x + y*y) <= (pictureS->mm)*mag)
			{
				int dx, dy;
				if (pictureS->sita > 90 && pictureS->sita <= 270)
				{
					dx = (int)(pictureS->x / roominout) - x;
					dy = (int)(pictureS->y / roominout) - y;
				}
				else
				{
					dx = (int)(pictureS->x / roominout) + x;
					dy = (int)(pictureS->y / roominout) + y;
				}

				if (dx >= 0 && dx < Width && dy >= 0 && dy < Height)
				{
					//cout << "123" << endl;
					img[(dy*Width+dx)*3 + 0] = 255 /255.0;
					img[(dy*Width+dx)*3 + 1] = 0 /255.0;
					img[(dy*Width+dx)*3 + 2] = 0 /255.0;
					//color[dy][dx].rgbtBlue = 255.0;
					//color[dy][dx].rgbtGreen = 0.0;
					//color[dy][dx].rgbtRed = 0.0;
				}
			}
			else
			{
				break;
			}
			++x;
		}
		if (pictureS->sita > 90 && pictureS->sita <= 270)
		{
			x = -x;
			y = -y;
		}
		//畫出箭頭的兩條斜線
		int xR, yR;
		int xL, yL;
		xR = xL = 1;
		int sitaR = pictureS->sita - 150;
		int sitaL = pictureS->sita + 150;
		if (sitaR < 0) sitaR = 360 + sitaR;
		if (sitaL >= 360) sitaL = sitaL - 360;
		while (true)//右下線
		{
			if (sitaR == 90 || sitaR == 270) yR = 0;
			else yR = xR * tan((sitaR)*M_PI / 180.0);

			if ((xR*xR + yR*yR) <= (pictureS->mm)*mag / 4.0)
			{
				int dx, dy;
				if (sitaR > 90 && sitaR <= 270)
				{
					dx = pictureS->x / roominout + x - xR;
					dy = pictureS->y / roominout + y - yR;
				}
				else
				{
					dx = pictureS->x / roominout + x + xR;
					dy = pictureS->y / roominout + y + yR;
				}

				if (dx >= 0 && dx < Width && dy >= 0 && dy < Height)
				{
					img[(dy*Width+dx)*3 + 0] = 255 /255.0;
					img[(dy*Width+dx)*3 + 1] = 0 /255.0;
					img[(dy*Width+dx)*3 + 2] = 0 /255.0;
					
					//color[dy][dx].rgbtBlue = 255.0;
					//color[dy][dx].rgbtGreen = 0.0;
					//color[dy][dx].rgbtRed = 0.0;
					
				}
			}
			else
			{
				break;
			}
			++xR;
		}

		while (true)//左下線
		{
			if (sitaL == 90 || sitaL == 270) yL = 0;
			else yL = xL * tan((sitaL)*M_PI / 180.0);

			if ((xL*xL + yL*yL) <= (pictureS->mm)*mag / 4.0)
			{
				int dx, dy;
				if (sitaL > 90 && sitaL <= 270)
				{
					dx = pictureS->x / roominout + x - xL;
					dy = pictureS->y / roominout + y - yL;
				}
				else
				{
					dx = pictureS->x / roominout + x + xL;
					dy = pictureS->y / roominout + y + yL;
				}

				if (dx >= 0 && dx < Width && dy >= 0 && dy < Height)
				{
					img[(dy*Width+dx)*3 + 0] = 255 /255.0;
					img[(dy*Width+dx)*3 + 1] = 0 /255.0;
					img[(dy*Width+dx)*3 + 2] = 0 /255.0;
					
					//color[dy][dx].rgbtBlue = 255.0;
					//color[dy][dx].rgbtGreen = 0.0;
					//color[dy][dx].rgbtRed = 0.0;
					
				}
			}
			else
			{
				break;
			}
			++xL;
		}
	}

	/*輸出圖片*/
	img.bmp("feaArrow.bmp", 24);
	//img_gray.bmp("feaArrow.bmp", 24);
}