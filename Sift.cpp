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

// 初始化
ImgRaw::ImgRaw(size_t width, size_t height, float val) :
	raw_img(width*height), width(width), height(height)
{
	if (val) {
		std::fill_n (raw_img.begin(), raw_img.size(),val/255.0);
	}
}
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
Sift::Sift(ImgRaw img): raw_img(img) {

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
// 高斯模糊(第一圖, 模糊幾次) / 回傳整個度(octvs)回去
vector<ImgRaw> Sift::dog_gau(ImgRaw& img, size_t s, size_t o) {
	// 初始化
	vector<ImgRaw> pyrs;
	for (unsigned i = 0; i < s; ++i) {
		pyrs.emplace_back(img.width, img.height);
	}
	// 高斯模糊
	float p = SIFT_GauSigma;
	float k = pow(2.f, 1.f / SIFT_Sacle);
	pyrs[0].sigma = p*o;
	//cout << pyrs[0].sigma << ", ";
	ImgRaw::gauBlur(pyrs[0], img, p);
	for (unsigned i = 1; i < s; ++i) {
		p *= k;
		pyrs[i].sigma = p * o; // 儲存kp多*o
		//cout << pyrs[i].sigma << ", ";
		ImgRaw::gauBlur(pyrs[i], pyrs[i - 1], p); // 圖片是取上一層的，不用再*o
	} //cout << endl;
	return pyrs;
}

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

void Sift::Gusto(ImgRaw& doImage, ImgRaw& doImage_ori, int InWidth, int InHeight, float sigma)
{
	vector<float> musk = Gaus::gau_matrix(sigma, 3);
	float block[3];
	float total;

	doImage.raw_img.reserve(InWidth*InHeight);
	float* gau_temp = new float[InWidth*InHeight];
	//做橫向的高斯矩陣
	for (int j = 0; j < InHeight; j++)
	{
		for (int i = 0; i < InWidth; i++)
		{
			if (i == 0)
			{
				block[0] = doImage_ori[j*InWidth + i + 0];
				block[1] = doImage_ori[j*InWidth + i + 0];
				block[2] = doImage_ori[j*InWidth + i + 1];
			}
			else if (i == (InWidth - 1))
			{
				block[0] = doImage_ori[j*InWidth + i - 1];
				block[1] = doImage_ori[j*InWidth + i + 0];
				block[2] = doImage_ori[j*InWidth + i + 0];
			}
			else
			{
				block[0] = doImage_ori[j*InWidth + i - 1];
				block[1] = doImage_ori[j*InWidth + i + 0];
				block[2] = doImage_ori[j*InWidth + i + 1];
			}
			total = 0;
			for (int v = 0; v < 3; v++)
			{
				total += block[v] * musk[v];
			}
			gau_temp[j*InWidth + i] = total;
		}
	}
	//做縱向的高斯矩陣
	for (int j = 0; j < InHeight; j++)
	{
		for (int i = 0; i < InWidth; i++)
		{
			if (j == 0)
			{
				block[0] = gau_temp[(j + 0)*InWidth + i];
				block[1] = gau_temp[(j + 0)*InWidth + i];
				block[2] = gau_temp[(j + 1)*InWidth + i];
			}
			else if (j == (InHeight - 1))
			{
				block[0] = gau_temp[(j - 1)*InWidth + i];
				block[1] = gau_temp[(j + 0)*InWidth + i];
				block[2] = gau_temp[(j + 0)*InWidth + i];
			}
			else
			{
				block[0] = gau_temp[(j - 1)*InWidth + i];
				block[1] = gau_temp[(j + 0)*InWidth + i];
				block[2] = gau_temp[(j + 1)*InWidth + i];
			}
			total = 0;
			for (int v = 0; v < 3; v++)
			{
				total += block[v] * musk[v];
			}
			doImage[j*InWidth + i] = total;
		}
	}
	delete[] gau_temp;
}
void Sift::GusB(vector<ImgRaw>& doImage,int inz, int InWidth, int InHeight, float sigma)
{
	vector<float> musk = Gaus::gau_matrix(sigma, 3);
	float block[3];
	float total;
	float* gau_temp = new float[InWidth*InHeight];
	//做橫向的高斯模糊
	for (int j = 0; j < InHeight; j++)
	{
		for (int i = 0; i < InWidth; i++)
		{
			if (i == 0)
			{
				block[0] = doImage[inz - 1][j*InWidth + i + 0];
				block[1] = doImage[inz - 1][j*InWidth + i + 0];
				block[2] = doImage[inz - 1][j*InWidth + i + 1];
			}
			else if (i == (InWidth - 1))
			{
				block[0] = doImage[inz - 1][j*InWidth + i - 1];
				block[1] = doImage[inz - 1][j*InWidth + i + 0];
				block[2] = doImage[inz - 1][j*InWidth + i + 0];
			}
			else
			{
				block[0] = doImage[inz - 1][j*InWidth + i - 1];
				block[1] = doImage[inz - 1][j*InWidth + i + 0];
				block[2] = doImage[inz - 1][j*InWidth + i + 1];
			}
			total = 0;
			for (int v = 0; v < 3; v++)
			{
				total += block[v]*musk[v];
			}
			gau_temp[j*InWidth + i] = total;
		}
	}

	//做縱向的高斯模糊
	for (int j = 0; j < InHeight; j++)
	{
		for (int i = 0; i < InWidth; i++)
		{
			if (j == 0)
			{
				block[0] = gau_temp[(j + 0)*InWidth + i];
				block[1] = gau_temp[(j + 0)*InWidth + i];
				block[2] = gau_temp[(j + 1)*InWidth + i];
			}
			else if (j == (InHeight - 1))
			{
				block[0] = gau_temp[(j - 1)*InWidth + i];
				block[1] = gau_temp[(j + 0)*InWidth + i];
				block[2] = gau_temp[(j + 0)*InWidth + i];
			}
			else
			{
				block[0] = gau_temp[(j - 1)*InWidth + i];
				block[1] = gau_temp[(j + 0)*InWidth + i];
				block[2] = gau_temp[(j + 1)*InWidth + i];
			}
			total = 0;
			for (int v = 0; v < 3; v++)
			{
				total += block[v] * musk[v];
			}
			doImage[inz][j*InWidth + i] = total;
		}
	}

}
// 高斯金字塔
void Sift::pyramid2() {
	ImgRaw first_img;
	float Limit;
	/*金字塔高*/
	int n = 1;
	for (Limit = 1.0;Limit < first_img.width || Limit < first_img.height; ++n){
		Limit = powf(2.0, n);
	}
	if (n < 7) {
		n = 1;
	} else {
		n -= 6;
	}
	Limit = pow(0.5, n);
	/*做金字塔*/
	float curr_size=2;
	for (size_t ph = 0; ph < 1; ph++) {
		//ImgRaw::first(currSize_img, raw_img, curr_size);
		const size_t curr_Width = raw_img.width*curr_size;
		const size_t curr_Height = raw_img.height*curr_size;
		ImgRaw currSize_img(curr_Width ,curr_Height);
		ZoomInOut(currSize_img, curr_Width, curr_Height);

		// 高斯模糊
		vector<ImgRaw> gau_imgs(pyWidth);
		ImgRaw::gauBlur(gau_imgs[0], currSize_img, SIFT_GauSigma);
		//Gusto(gau_imgs[0], currSize_img, curr_Width, curr_Height, SIFT_GauSigma);
		for (size_t i = 1; i < pyWidth; i++) {
			const float curr_Sigma = SIFT_GauSigma*pow(sqrt(2.0),i);
			gau_imgs[i].raw_img.resize(curr_Width*curr_Height);
			Gaus::raw2GauBlur(gau_imgs[i], gau_imgs[i-1], curr_Width, curr_Height, curr_Sigma);
			//Gusto(gau_imgs[i], gau_imgs[i-1], curr_Width, curr_Height, curr_Sigma);
			//gau_imgs[i].bmp("_gau.bmp", 8);
		}
		// 高斯差分
		vector<ImgRaw> gauDog_imgs(pyWidth-1);
		for (size_t i = 0; i < pyWidth-1; i++) {
			ImgRaw img(curr_Width, curr_Height);
			for (size_t idx = 0; idx < curr_Width*curr_Height; idx++) {
				img[idx] = gau_imgs[i][idx] - gau_imgs[i+1][idx];
			}
			gauDog_imgs[i]=img;
			gauDog_imgs[i].bmp("gauc/gauDog"+to_string(i)+".bmp", 8);
		}
		// 尋找特徵點並累加直方圖
		for (size_t scale_idx = 1; scale_idx < pyWidth-1-1; scale_idx++) {
			ImgRaw teest_img(curr_Width, curr_Height);
			// 特徵直方圖累加半徑
			const size_t r = curr_size * 3. * 1.5 * (scale_idx+1);
			for (size_t j = r+1; j < curr_Height-r-1; j++) {
				for (size_t i = r+1; i < curr_Width-r-1; i++) {
					// 尋找極值點
					const float& Val = gauDog_imgs[scale_idx][j*curr_Width + i];
					if (findMaxMin(gauDog_imgs, scale_idx, curr_Width, j, i)) {
						teest_img[j*curr_Width + i] = 255 /255.0;
					}
				}
			}
			teest_img.bmp("MaxMin/MaxMin"+to_string(scale_idx)+".bmp", 8);
		}




		curr_size /= 2; // 每次遞減
		cout << endl;
	}
	return;
}

void Sift::pyramid(size_t s) {
	size_t sacle = s+SIFT_SacleDiff; // 找極值捨去前後兩張 + 差分圖少1張
    size_t octvs = 2; // 不同大小圖片(自訂)
    //octvs = (size_t)(log(min(raw_img.width, raw_img.height)) / log(2.0)-2);
    pyrs.resize(octvs);

    // 輸入圖
    ImgRaw temp(raw_img.width, raw_img.height);
	ImgRaw::first(temp, raw_img, 2);
	temp.sigma = SIFT_GauSigma;
    // 高斯金字塔
    pyrs[0] = dog_gau(temp, sacle);
    comp(pyrs[0], "Sift-gau_" + to_string(0));
    for (unsigned i = 1; i < octvs; ++i) {
        // 縮小一張圖(中間那張的p倍率正好會是下一排的第一張)
        ImgRaw::first(temp, pyrs[i - 1][SIFT_SacleDiff], 0.5);
        // 回傳 s 張模糊圖
        pyrs[i] = dog_gau(temp, sacle, i*2);
        // 輸出高斯模糊圖
        comp(pyrs[i], "Sift-gau_" + to_string(i));
    }

// 測試 sigma 數值
//for (size_t o = 0; o < pyrs.size(); o++) {
//	for (size_t s = 0; s < pyrs[0].size(); s++) {
//		cout << pyrs[o][s].sigma << ", ";
//	} cout << endl;
//} cout << "==================================" << endl;

	// 差分圖金字塔
	pyrs_dog = pyrs;
    // 高斯金字塔差分圖
    for (unsigned j = 0; j < octvs; ++j) {
        // 高斯差分
        for (unsigned i = 0; i < pyrs[j].size() - 1; ++i) {
			pyrs_dog[j][i] = pyrs[j][i] - pyrs[j][i+1];
        }
		// 移除最後一張
		pyrs_dog[j].erase(--(pyrs_dog[j].end()));
		// 輸出差分圖
        comp(pyrs_dog[j], "Sift-diff_" + to_string(j));
    }

// 測試 sigma 數值
//for (size_t o = 0; o < pyrs.size(); o++) {
//	for (size_t s = 0; s < pyrs[0].size(); s++) {
//		cout << pyrs[o][s].sigma << ", ";
//	} cout << endl;
//}

    // 取得遮罩(沒有邊緣防呆)
    auto getMask = [](vector<types>& mask,
        ImgRaw& img, size_t y, size_t x)
    {
        // 複製遮罩
        mask.resize(9);
        mask[0] = img.at2d(y - 1, x - 1);
        mask[1] = img.at2d(y - 1, x + 0);
        mask[2] = img.at2d(y - 1, x - 1);

        mask[3] = img.at2d(y + 0, x - 1);
        mask[4] = img.at2d(y + 0, x + 0);
        mask[5] = img.at2d(y + 0, x - 1);

        mask[6] = img.at2d(y + 1, x - 1);
        mask[7] = img.at2d(y + 1, x + 0);
        mask[8] = img.at2d(y + 1, x + 1);
    };
    // 尋找 cubic 極值
    vector<types> mask; // 臨時方塊27點找極值
	vector<Fea_point> feap;
    for (unsigned py = 0; py < pyrs_dog.size(); ++py) {
        for (unsigned px = 1; px < pyrs_dog[py].size() - 1; ++px) { // 捨棄前後兩張
			size_t gau_r;// 當前這張圖的高斯半徑
			gau_r = (size_t)(2.f*1.f/powf(2.f, py) * (px+1) * 1.5f * 3.f);
//cout << 2.f*1.f/powf(2.f, py) << ", " << (px+1) << "----";
//cout << gau_r << ", ";
            ImgRaw fea(pyrs_dog[py][px].width, pyrs_dog[py][px].height); // 暫存畫布
			// 選定層級找極值(#防呆去除邊緣)(移除半徑會超出邊緣的點)
            for (unsigned j = 1+gau_r; j < pyrs_dog[py][px].height - 2 - gau_r; ++j) {
                for (unsigned i = 1+gau_r; i < pyrs_dog[py][px].width - 2 - gau_r; ++i) {
                    // 取得上下層的九宮格
                    getMask(mask, pyrs_dog[py][px - 1], j, i);
                    getMask(mask, pyrs_dog[py][px], j, i);
                    getMask(mask, pyrs_dog[py][px + 1], j, i);
					// 找兩極值
                    float max = *std::max_element(mask.begin(), mask.end());
                    float min = *std::min_element(mask.begin(), mask.end());
                    size_t mid = (mask.size() - 1) / 2;

					// 保留極值
					if (mask[mid] == max or mask[mid] == min) {
                        constexpr float thre = SIFT_Dx / 2.f;
						// 保留變動大於 0.03
						if (abs(mask[mid]) > thre) { 
							// 角點偵測
							if (Corner::harris(pyrs_dog[py][px],
								pyrs_dog[py][px].width, j, i, SIFT_HarrisR) == 1)
							{
								// 計算特徵點的統計方向與強度
								float fea_m = 0.f;
								float fea_sida = 0.f;
								float* fea_p = getFea(pyrs[py][px], j, i, pyrs_dog[py][px].sigma, gau_r);
								// 建立特徵點的統計方向與強度
								// feap.emplace_back(py, px, j, i, gau_r, pyrs_dog[py][px].sigma, fea_p[0], fea_p[1]);
								float rate=1.0;
								//cout << fea_p[0]*rate << ", " << fea_p[1] << endl;
								delete fea_p;
								//Draw::draw_line(fea, j, i, fea_p[0]*rate, fea_p[1]);
								//Draw::draw_line(fea, j, i, 10, fea_p[1]);

								fea.at2d(j, i) = 1; // 畫點
							}
                        }
                    }
                }
            }
            // 特徵點圖片預覽
            fea.bmp("fea/fea" + to_string(py) + "-" + to_string(px) + ".bmp", 8);
			
        }
//cout << endl;
    }



	/* pyrs到這裡已經是特徵點圖 */

	/*
	vector<types> getGau;
	auto getGau = [](vector<types>& mask,
		ImgRaw& img, size_t y, size_t x)
	{
		// 複製遮罩
		mask.resize();
	};*/

	// feap存了各階層的特徵點




} // 高斯金字塔

// 計算特徵描述長度
float Sift::fea_m(ImgRaw& img, size_t y, size_t x) {
	// 檢查邊界
	if (x<=0 or y<=0 or x>= img.width-2 or y>= img.height-2) {
		Out_of_imgRange("超出圖片邊界");
	}
	// 公式
	float m = sqrtf(
		powf(img.at2d(y, x+1)-img.at2d(y, x-1), 2)+
		powf(img.at2d(y+1, x)-img.at2d(y-1, x), 2)
	);
	return m;
}
// 計算特徵描述角度
float Sift::fea_sida(ImgRaw& img, size_t y, size_t x) {
	// 檢查邊界
	if (x<=0 or y<=0 or x>= img.width-2 or y>= img.height-2) {
		Out_of_imgRange("超出圖片邊界");
	}
	// 公式
	float valY = img.at2d(y+1, x)-img.at2d(y-1, x);
	float valX = img.at2d(y, x+1)-img.at2d(y, x-1);
	float sida = atan2(valY, valX) * (180/M_PI);
	// 修正負號
	if (sida < 0) {
		sida += 360.f;
	}
	return sida;
}
// 輸入原圖與點，計算強度與角度
float* Sift::getFea(ImgRaw& img, size_t y, size_t x, float sigma, size_t r) {
	size_t diam = r*2;
	// 高斯矩陣
	vector<float> gau_2dmat;
	Gaus::gau_matrix2d(gau_2dmat, sigma, diam+1);
	// 獲取特徵點計算半徑內的所有 s, m
	map<int, float> fea_hist;
	for (size_t j = 0; j < diam+1; j++) {
		for (size_t i = 0; i < diam+1; i++) {
			// 當前點
			float m = fea_m(img, y-r+j, x-r+i) * gau_2dmat[j*diam + i];
			float s = fea_sida(img, y-r+j, x-r+i);
			fea_hist[floor(s/10)] += m;
//cout << m << ",	" << floor(s/10) << endl;
		}
	}

// 找出主方向
//cout << "*********" << endl;
	float value = fea_hist[0];
	float mapidx = 0;
	for (size_t i = 1; i < 36; i++) {
//cout << fea_hist[i] << endl;
		if (value < fea_hist[i]) {
			value = fea_hist[i];
			mapidx = i;
		}
	}
	float* fea_value = new float[2]{value, mapidx*10};
	return fea_value;
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
