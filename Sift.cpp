/*****************************************************************
Name :
Date : 2017/07/05
By   : CharlotteHonG
Final: 2017/07/05
*****************************************************************/
#pragma warning(disable : 4819)
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <random>

#if defined(_MSC_VER)
	#define or ||
	#define and &&
	#define OR ||
	#define AND &&
#endif
#define M_PI 3.14159265358979323846

#include "Sift.hpp"
using namespace std;

// Sift建構子
Sift::Sift(ImgRaw img, size_t intvls): raw_img(img), pyWidth(intvls) {
	FeatStart = new Feature;
	FeatEnd = FeatStart;
}
// 合成圖片
void Sift::comp(vector<ImgRaw>& pyrs, string name) {
    // 合成圖片
    size_t w = pyrs[0].width, h = pyrs[0].height;
    size_t s2 = pyrs.size();
    ImgRaw big_img(w*s2, h, 8);
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
void Sift::AddnewFeaturestruct(int Inx, int Iny, float Insize, int scale_idx, int sigmaOCT, float Inm, int Insita)
{
	Feature* newnode = new Feature;
	newnode->x = Inx;
	newnode->y = Iny;
	newnode->mm = Inm;
	newnode->sita = Insita;
	newnode->size = Insize;
	newnode->kai = scale_idx;
	newnode->sigmaOCT = sigmaOCT;
	newnode->nextptr = nullptr;

	FeatEnd->nextptr = newnode;// 加入該點
	FeatEnd = newnode; // 更新結尾點的標記
}
// 過濾特徵點並找極值
bool Sift::findMaxMin(vector<ImgRaw>& gauDog_imgs, size_t scale_idx, size_t curr_Width, size_t y, size_t x) {
	const float& Val = gauDog_imgs[scale_idx][y*curr_Width + x];
	// 過濾變化率過小的特徵點
	if (abs(Val) > SIFT_Dx/2.0) {
		// 找極值
		for (int k = (scale_idx - 1); k <= (scale_idx + 1); k++) {
			const ImgRaw& currImg = gauDog_imgs[k];
			for (int j = (y - 1); j <= (y + 1); j++) {
				size_t currIdx = j*curr_Width+x;
				for (int i = - 1; i <= 1; i++) {
					if ((Val>0 and Val<currImg[currIdx +i]) or 
						(Val<0 and Val>currImg[currIdx +i]))
					{
						return false;
					}
				}
				//if ((Val>0 and Val<currImg[currIdx -1]) or (Val<0 and Val>currImg[currIdx -1])) {return false;}
				//if ((Val>0 and Val<currImg[currIdx +0]) or (Val<0 and Val>currImg[currIdx +0])) {return false;}
				//if ((Val>0 and Val<currImg[currIdx +1]) or (Val<0 and Val>currImg[currIdx +1])) {return false;}
			}
		} return true;
	} return false;
};
// 獲得特徵點直方圖統計

static inline float getMag(float dx, float dy){
	// 獲得強度
	return sqrtf(dx*dx + dy*dy);
}
static inline float getSita(float dx, float dy){
	// 獲得角度
	return fastAtan2f(dy, dx);
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
			sita[idx] = getSita(dx, dy);
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
// 產生特徵描述子
static inline void CoordinateChange(int& deltaX, int& deltaY, float sita) {
	// 座標轉換
	float deltaX2 = deltaY*sin(sita) + deltaX*cos(sita);
	deltaY = deltaY*cos(sita) - deltaX*sin(sita);
	deltaX = deltaX2;
}
void Sift::FeatureDescrip(vector<ImgRaw>& kaidaImag, Feature* FeatureNow) {
	//新增一個4*4*8的特徵描述空間
	vector <vector <vector <float>>> descripgroup(4);//儲存特徵描述子
	for (int v = 0; v < 4; v++) {
		descripgroup[v].resize(4);
		for (int j = 0; j < 4; j++) {
			descripgroup[v][j].resize(8);
		}
	}
	// 從當前極值產出的點開始做
	for (Feature* FeatureS = FeatureNow; FeatureS->nextptr != nullptr;) {
		FeatureS = FeatureS->nextptr;
		// 轉制
		int radius = (FeatureS->sigmaOCT * 3 * sqrt(2)*(5) + 1) / 2;
		int inHeight = raw_img.height*FeatureS->size;
		int inWidth = raw_img.width*FeatureS->size;
		// 直徑x直徑
		for (int j = FeatureS->y - radius; j < FeatureS->y + radius; j++) {
			for (int i = FeatureS->x - radius; i < FeatureS->x + radius; i++) {
				if (j < 1 || j >= inHeight - 1 || i < 1 || i >= inWidth - 1) {
					continue; // 邊緣不算
				}

				int newx, newy;
				newx = i - FeatureS->x;
				newy = j - FeatureS->y;
				CoordinateChange(newx, newy, (FeatureS->sita * M_PI / 180.0));
				float blockx = newx / ((radius*sqrt(2) / 2.0)) * 2.0 + 3.0;// 取範圍落在1~4內的值
				float blocky = newy / ((radius*sqrt(2) / 2.0)) * 2.0 + 3.0;// 取範圍落在1~4內的值

				
				float sita, mm;
				const float* curryImg = kaidaImag[FeatureS->kai].raw_img.data();

				float dx = curryImg[j*inWidth + i+1] - curryImg[j*inWidth + i-1];
				float dy = curryImg[(j+1)*inWidth + i]- curryImg[(j-1)*inWidth + i];
				mm = getMag(dx, dy);
				sita = getSita(dx, dy);
				
				const float unit = 360.0/8.0;
				float weix0, weix1;
				float weiy0, weiy1;
				weix1 = blockx - (int)blockx; weix0 = 1 - weix1;
				weiy1 = blocky - (int)blocky; weiy0 = 1 - weiy1;

				//****** 以firstorder將值分別存入descripgroup容器裡面
				int bx0, bx1;
				int by0, by1;
				bx0 = blockx;
				bx1 = bx0 + 1;
				by0 = blocky;
				by1 = by0 + 1;

				float RU, RD, LU, LD;
				float SitaGroup;

				SitaGroup = sita / unit;
				
				if (bx0 >= 1 && bx0 <= 4 && by0 >= 1 && by0 <= 4)//RU
				{
					RU = weix0 * weiy0 * mm;
					descripgroup[by0 - 1][bx0 - 1][SitaGroup] += RU;
				}
				if (bx0 >= 1 && bx0 <= 4 && by1 >= 1 && by1 <= 4)//RD
				{
					RD = weix0 * weiy1 * mm;
					descripgroup[by1 - 1][bx0 - 1][SitaGroup] += RD;
				}
				if (bx1 >= 1 && bx1 <= 4 && by0 >= 1 && by0 <= 4)//LU
				{
					LU = weix1 * weiy0 * mm;
					descripgroup[by0 - 1][bx1 - 1][SitaGroup] += LU;
				}
				if (bx1 >= 1 && bx1 <= 4 && by1 >= 1 && by1 <= 4)//LD
				{
					LD = weix1 * weiy1 * mm;
					descripgroup[by1 - 1][bx1 - 1][SitaGroup] += LD;
				}
				
			}
		}
		//****** 限定門檻值
		float add = 0.0;
		for (int j = 0; j < 4; j++) {
			for (int i = 0; i < 4; i++) {
				for (int v = 0; v < 8; v++) {
					if (descripgroup[j][i][v] > 0.2) 
						descripgroup[j][i][v] = 0.2;
					add += descripgroup[j][i][v];
				}
			}
		}
		//****** 向量正規化
		add = sqrt(add);
		for (int j = 0; j < 4; j++) {
			for (int i = 0; i < 4; i++) {
				for (int v = 0; v < 8; v++) {
					descripgroup[j][i][v] /= add;
				}
			}
		}
		FeatureS->descrip = descripgroup;
	}
}
// 高斯金字塔
static inline void ZoomInOut(ImgRaw& doImage, const ImgRaw& gray, int InWidth, int InHeight) {
// 放大縮小
	// 原圖資訊(取自資料成員)
	const size_t Height = gray.height;
	const size_t Width = gray.width;

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
void Sift::pyramid() {
	ImgRaw gray_img = raw_img.ConverGray();
	size_t octv = getPyramidOctv(gray_img.width, gray_img.height);
	float curr_size=2.f, curr_sigmaCoef=1.f;
	for (size_t ph = 0; ph < octv; ph++, curr_size /= 2.f, curr_sigmaCoef*=2.f) { // 金字塔的sigma還沒修復第二層沒*2
		const size_t curr_Width = gray_img.width*curr_size;
		const size_t curr_Height = gray_img.height*curr_size;
		ImgRaw currSize_img(curr_Width ,curr_Height, 8);
		ImgRaw::first(currSize_img, gray_img, curr_size);
		//ZoomInOut(currSize_img, gray_img, curr_Width, curr_Height);

		// 高斯模糊
		vector<ImgRaw> gau_imgs(pyWidth);
		ImgRaw::gauBlur(gau_imgs[0], currSize_img, SIFT_GauSigma * (curr_sigmaCoef));
		//gau_imgs[0].bmp("gau/gau"+to_string(0)+".bmp", 8);
		for (size_t i = 1; i < pyWidth; i++) {
			const float curr_Sigma = SIFT_GauSigma * powf(sqrt(2.0),i) * (curr_sigmaCoef); // 這裡的2是 S+3 的 S
			gau_imgs[i].resize(curr_Width, curr_Height, 8);
			Gaus::GauBlur(gau_imgs[i], gau_imgs[i-1], curr_Width, curr_Height, curr_Sigma);
			//gau_imgs[i].bmp("gau/gau"+to_string(i)+".bmp", 8);
		}
		// 高斯差分
		vector<ImgRaw> gauDog_imgs(pyWidth-1);
		for (size_t i = 0; i < pyWidth-1; i++) {
			gauDog_imgs[i].resize(curr_Width, curr_Height, 8);
			for (size_t idx = 0; idx < curr_Width*curr_Height; idx++) {
				gauDog_imgs[i][idx] = gau_imgs[i][idx] - gau_imgs[i+1][idx];
			}
			//gauDog_imgs[i].bmp("gauDog/gauDog"+to_string(i)+".bmp", 8);
		}
		// 紀錄當前指針位置
		Feature* FeatureNow = FeatEnd;
		// 尋找特徵點並累加直方圖
		for (size_t scale_idx = 1; scale_idx < pyWidth-1-1; scale_idx++) {
			ImgRaw MaxMin_img(curr_Width, curr_Height, 8);
			ImgRaw Herris_img(curr_Width, curr_Height, 8);
			// 特徵直方圖累加半徑
			const size_t r = curr_size * 3. * 1.5 * (scale_idx+1);
			for (size_t j = r+1; j < curr_Height-r-1; j++) {
				for (size_t i = r+1; i < curr_Width-r-1; i++) {
					// 尋找極值點
					if (findMaxMin(gauDog_imgs, scale_idx, curr_Width, j, i)) {
						const ImgRaw& currImg = gauDog_imgs[scale_idx];
						const float currSigma = SIFT_GauSigma * powf(sqrt(2.0),scale_idx) * (curr_sigmaCoef);
						if (Corner::harris(currImg, curr_Width, j, i, r)) {
							getHistogramMS(currImg, curr_size, scale_idx, currSigma, j, i, curr_Width, r);
							Herris_img[j*curr_Width + i] = 255 /255.0;
						}
						MaxMin_img[j*curr_Width + i] = 255 /255.0;
					}
				}
			}
			//MaxMin_img.bmp("MaxMin/MaxMin_"+to_string(scale_idx)+".bmp", 8);
			//Herris_img.bmp("Herris/Herris_"+to_string(scale_idx)+".bmp", 8);

		}
		// 描述剛剛才找到的新特徵點(尋找前先記錄位置)
		FeatureDescrip(gau_imgs, FeatureNow);
	}
}
// 印出箭頭
void Sift::drawArrow(string name){
	//箭頭倍率
	float mag = 10000.f;
	// 臨時圖庫
	ImgRaw img = raw_img;
	for (Feature* feaPoint = FeatStart;
		feaPoint->nextptr != NULL;
		feaPoint = feaPoint->nextptr)
	{
		size_t x = feaPoint->x/feaPoint->size;
		size_t y = feaPoint->y/feaPoint->size;
		Draw::draw_arrowRGB(img, y, x, (feaPoint->mm)*mag, feaPoint->sita);
	}
	/*輸出圖片*/
	img.bmp(name);
}


/************/
/*** 匹配 ****/
/************/
Stitching::Stitching(const Sift& desc1, const Sift& desc2) {
	FeatureStart1 = desc1.FeatStart;
	FeatureStart2 = desc2.FeatStart;;
	const ImgRaw& img1=desc1.raw_img;
	const ImgRaw& img2=desc2.raw_img;

	Width=img1.width+img2.width;
	Height=img1.height;
	// 合併兩張圖
	matchImg.resize(img1.width*2, img2.height, 24);
	for (size_t j = 0; j < img1.height; j++) {
		for (size_t i = 0; i < img1.width; i++) {
			matchImg[(j*matchImg.width+i)*3 + 0] = img1[(j*img1.width+i)*3 + 0];
			matchImg[(j*matchImg.width+i)*3 + 1] = img1[(j*img1.width+i)*3 + 1];
			matchImg[(j*matchImg.width+i)*3 + 2] = img1[(j*img1.width+i)*3 + 2];
		}
		for (size_t i = img1.width; i < img2.width+img1.width; i++) {
			matchImg[(j*matchImg.width+i)*3 + 0] = img2[(j*img1.width+i-img1.width)*3 + 0];
			matchImg[(j*matchImg.width+i)*3 + 1] = img2[(j*img1.width+i-img1.width)*3 + 1];
			matchImg[(j*matchImg.width+i)*3 + 2] = img2[(j*img1.width+i-img1.width)*3 + 2];
		}
	}
	// matchImg.bmp("matchImg.bmp", 24); // 測試
	cout << endl;
}
// 檢查是否有相同的特徵描述子(歐式距離計算)
float Stitching::EuclDist(const Desc& point1, const Desc& point2) {
	float sum = 0.f;
	for (size_t k = 0; k < 4; k++) {
		for (size_t j = 0; j < 4; j++) {
			for (size_t i = 0; i < 8; i++) {
				const float reduce = point1[k][j][i] - point2[k][j][i];
				sum += reduce * reduce;
			}
		}
	}
	return sqrtf(sum);
}
void Stitching::Check(float matchTh) {
	Feature *startlink1, *startlink2;
	startlink1 = FeatureStart1;
	startlink2 = FeatureStart2;
	int pc = 0;
	while (startlink1->nextptr != nullptr) {
		startlink1 = startlink1->nextptr;

		float distance1, distance2;
		distance1 = distance2 = 0xffffffff;
		float useX1, useY1;
		useX1 = useY1 = 0xffffffff;

		startlink2 = FeatureStart2;
		while (startlink2->nextptr != nullptr) {
			startlink2 = startlink2->nextptr;
			//******* 找處距離最近的兩個點
			float value = EuclDist(startlink1->descrip, startlink2->descrip);
			if (value < distance1) {
				distance2 = distance1;
				distance1 = value;

				useX1 = startlink2->x / startlink2->size;
				useY1 = startlink2->y / startlink2->size;
			} else if (value < distance2) {
				distance2 = value;
			}
		}
		//****** 確認是否匹配成功並畫線
		if ((distance1 / distance2) < matchTh) {
			//閥值越小找到的點越少但越可靠，論文建議 0.6~0.8
			int x1, y1;
			x1 = startlink1->x / startlink1->size;
			y1 = startlink1->y / startlink1->size;
			
			// 隨機顏色
			auto random_num = [] {
				return ((rand() / (RAND_MAX+1.0)) * (1 - 0) + 0);
			};
			float rVal = random_num();
			float gVal = random_num();
			float bVal = random_num();
			// 畫線
			Draw::drawLineRGB_p(matchImg, y1, x1, useY1,useX1 + (Width / 2), rVal, gVal, bVal);

			++pc;
		}
	}
	matchImg.bmp("_matchImg.bmp", 24);
}

