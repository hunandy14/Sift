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
#include <limits>

#if defined(_MSC_VER)
	#define or ||
	#define and &&
	#define OR ||
	#define AND &&
#endif
#define M_PI 3.14159265358979323846
#define SQUARE2 1.4142135623730951f

#include "Sift.hpp"
using namespace std;

/*
     ######  #### ######## ########
    ##    ##  ##  ##          ##
    ##        ##  ##          ##
     ######   ##  ######      ##
          ##  ##  ##          ##
    ##    ##  ##  ##          ##
     ######  #### ##          ##
*/
// Sift建構子
Sift::Sift(ImgRaw img, size_t intvls): raw_img(img), pyWidth(intvls) {
	FeatStart = new Feature; // 第一點為空
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
void Sift::FeatAppend(int Inx, int Iny, float Insize, int scale_idx, int sigmaOCT, float Inm, int Insita)
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
void Sift::pyramid() {
	ImgRaw gray_img = raw_img.ConverGray();
	size_t octv = getPyramidOctv(gray_img.width, gray_img.height);
	float curr_size=2.f; // 尺度倍率
	float curr_sigmaCoef=1.f; // 模糊度倍率
	for (size_t ph = 0; ph < octv; ph++, curr_size /= 2.f, curr_sigmaCoef*=2.f) {
		const size_t curr_Width = gray_img.width*curr_size;
		const size_t curr_Height = gray_img.height*curr_size;
		ImgRaw currSize_img(curr_Width ,curr_Height, 8);
		ImgRaw::first(currSize_img, gray_img, curr_size);
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
		// 紀錄當前指針位置開始的位置
		Feature* CurrFeature = FeatEnd;
		// 尋找特徵點並累加直方圖
		for (size_t scale_idx = 1; scale_idx < pyWidth-1-1; scale_idx++) {
			//ImgRaw MaxMin_img(curr_Width, curr_Height, 8);
			//ImgRaw Herris_img(curr_Width, curr_Height, 8);
			// 特徵直方圖累加半徑
			const size_t r = curr_size * 3. * 1.5 * (scale_idx+1);
			for (size_t j = r+1; j < curr_Height-r-1; j++) {
				for (size_t i = r+1; i < curr_Width-r-1; i++) {
					// 尋找極值點
					if (findMaxMin(gauDog_imgs, scale_idx, curr_Width, j, i)) {
						const ImgRaw& currImg = gauDog_imgs[scale_idx];
						const float currSigma = SIFT_GauSigma * powf(sqrt(2.0),scale_idx) * (curr_sigmaCoef);
						if (harris(currImg, curr_Width, j, i, r)) { // r=SIFT_CURV_THR=10
							// 這裡不知道為什麼設置為r效果更好論文所提設置是10
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
		// 描述剛剛才找到的那堆新特徵點(當前的點是舊的)
		FeatureDescrip(gau_imgs, CurrFeature->nextptr);
	} // 不同尺度大小的for()
}
bool Sift::harris(const vector<float>& p,
	size_t w, size_t y, size_t x, float r)
{
	// 閥值 (論文r=10)
	float thre = ((r + 1)*(r + 1)) / r;
	// 二維讀取
	auto at2d = [&](int y, int x)->float {
		return p[y*w + x];
	};
	// 公式
	float Dxx = 2.f*at2d(y, x) - at2d(y, x - 1) - at2d(y, x + 1);
	float Dyy = 2.f*at2d(y, x) - at2d(y - 1, x) - at2d(y + 1, x);
	float Dxy = at2d(y + 1, x + 1) + at2d(y - 1, x - 1)
		- at2d(y - 1, x + 1) - at2d(y + 1, x - 1);
	Dxy /= 4.f;
	float Tr = Dxx + Dyy;
	float Det = Dxx*Dyy - Dxy*Dxy;
	// 判斷閥值
	float val = (Tr*Tr / Det);
	if (val < thre) {
		return 1;
	}
	// 不成立則刪除
	return 0;
}
void Sift::getHistogramMS(const ImgRaw& doImage, float Insize, size_t scale, float sigma, 
	size_t Iny, size_t Inx, size_t InWidth, size_t Inr)
{
	const size_t matLen = Inr*2 + 1;     // 遮罩長度
	vector<float> mag(matLen * matLen);  // 儲存強度
	vector<float> sita(matLen * matLen); // 儲存角度
	//做m、sita計算
	for (int j = (Iny - Inr), idx = 0; j <= (Iny + Inr); j++) {
		for (int i = (Inx - Inr); i <= (Inx + Inr); i++, idx++) {
			const float& dx = doImage[(j+0)*InWidth + (i+1)] - doImage[(j+0)*InWidth + (i-1)];
			const float& dy = doImage[(j+1)*InWidth + (i+0)] - doImage[(j-1)*InWidth + (i+0)];
			mag[idx] = sqrtf(dx*dx + dy*dy); // 獲得強度
			sita[idx] = fastAtan2f(dy, dx);  // 獲得角度
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
		if (magSum[i] >= maxMag*0.8) {// 大於80%的都存下來
			FeatAppend(Inx, Iny, Insize, scale, sigma, magSum[i], i*10);
		}
	}
}


// 描述特徵子
bool Sift::calc_grad_mag_ori(const vector<float> &img, int &COL, int &ROW, int r, int c, float &mag, float &ori) {
	float dx = 0.0, dy = 0.0;
	int Col = COL;
	if (r > 0 && r < ROW - 1 && c > 0 && c < COL - 1) {
		dx = img[(r)* Col + (c + 1)] - img[(r)* Col + (c - 1)];
		dy = img[(r - 1) * Col + (c)] - img[(r + 1) * Col + (c)];
		mag = sqrt(dx * dx + dy * dy);
		//mag = dx > dy ? max(0.875f*dx+0.5f*dy, dx) : max(0.875f*dy + 0.5f*dx, dy);
		ori = atan2(dy, dx);
		return true;
	} else {
		return false;
	}
}
Sift::Desc Sift::descr_hist(vector<float> &img, int &COL, int &ROW, int r, int c, float ori, float scl, int d, int n) {
	vector<vector<vector<float>>> hist;
	float cos_t, sin_t, hist_width, exp_denom, r_rot, c_rot, grad_mag,
		grad_ori, w, rbin, cbin, obin, bins_per_rad, PI2 = 2.f * M_PI;
	int radius, i, j;

	hist.resize(d);
	for (i = 0; i < d; i++) {
		hist[i].resize(d);
		for (j = 0; j < d; j++) {
			hist[i][j].resize(n);
		}
	}

	cos_t = cos(ori);
	sin_t = sin(ori);
	bins_per_rad = n / PI2;
	exp_denom = d * d * 0.5f;
	hist_width = (float)SIFT_DESCR_SCL_FCTR * scl;
	radius = (int)(hist_width * sqrt(2.f) * (d + 1.f) * 0.5f + 0.5f);

	for (i = -radius; i <= radius; i++) {
		for (j = -radius; j <= radius; j++) {
			// 把位置縮放到 [0~4)
			c_rot = (j * cos_t - i * sin_t) / hist_width;
			r_rot = (j * sin_t + i * cos_t) / hist_width;
			rbin = r_rot + d / 2.f - 0.5f;
			cbin = c_rot + d / 2.f - 0.5f;
			/*cout << hist_width << endl;
			cout << radius << endl;
			cout << c_rot << " " << r_rot << endl;
			cout << cbin << " " << rbin << endl;
			system("pause");*/

			// 確保位置 [0~4) 刪除超出原圖指定半徑外的點
			if (rbin > -1.f  &&  rbin < d  &&  cbin > -1.f  &&  cbin < d) {
				// 計算成功回傳1 並修改數值
				if (calc_grad_mag_ori(img, COL, ROW, r + i, c + j, grad_mag, grad_ori)) {
					grad_ori -= ori;
					// 修正0~360
					while (grad_ori < 0.f) {
						grad_ori += PI2;
					}
					while (grad_ori >= PI2) {
						grad_ori -= PI2;
					}

					obin = grad_ori * bins_per_rad;// 高斯加權分布
					w = exp(-(c_rot * c_rot + r_rot * r_rot) / exp_denom); // 計算權值
					interp_hist_entry(hist, rbin, cbin, obin, grad_mag * w, d, n); //插補
				}
			}
		}
	}
	return hist;
}
void Sift::interp_hist_entry(Desc &hist, float rbin, float cbin, float obin, float mag, int d, int n) {
	float d_r, d_c, d_o;
	int r0, c0, o0;

	r0 = (int)floor(rbin);
	c0 = (int)floor(cbin);
	o0 = (int)floor(obin);
	d_r = rbin - r0;
	d_c = cbin - c0;
	d_o = obin - o0;

	for (int r = 0; r <= 1; r++) {
		int rb = r0 + r;
		if (rb >= 0 && rb < d) {
			float v_r = mag * ((r == 0) ? 1.f - d_r : d_r);
			vector<vector<float>>& row = hist[rb];
			for (int c = 0; c <= 1; c++) {
				int cb = c0 + c;
				if (cb >= 0 && cb < d) {
					float v_c = v_r * ((c == 0) ? 1.f - d_c : d_c);
					vector<float>& h = row[cb];
					for (int o = 0; o <= 1; o++) {
						int ob = (o0 + o) % n;
						float v_o = v_c * ((o == 0) ? 1.f - d_o : d_o);
						h[ob] += v_o;
					}
				}
			}
		}
	}
}
void Sift::hist_to_descr(const Desc &hist, int d, int n, Feature* feat) {
	int int_val, i, r, c, o, k = 0;

	for (r = 0; r < d; r++) {
		for (c = 0; c < d; c++) {
			for (o = 0; o < n; o++) {
				feat->descr[k++] = hist[r][c][o];
			}
		}
	}
	feat->d = k;
	normalize_descr(feat);
	for (i = 0; i < k; i++) {
		if (feat->descr[i] > SIFT_DESCR_MAG_THR) {
			feat->descr[i] = (float)SIFT_DESCR_MAG_THR;
		}
	}
	normalize_descr(feat);

	// convert floating-point descriptor to integer valued descriptor
	for (i = 0; i < k; i++) {
		int_val = (int)((SIFT_INT_DESCR_FCTR) * feat->descr[i]);
		feat->descr[i] = (float)min(255, int_val);
	}
}
void Sift::normalize_descr(Feature* feat) {
	float cur, len_inv, len_sq = 0.f;
	int i, d = feat->d;

	for (i = 0; i < d; i++) {
		cur = feat->descr[i];
		len_sq += cur*cur;
	}
	len_inv = 1.f / sqrt(len_sq);
	for (i = 0; i < d; i++)
		feat->descr[i] *= len_inv;
}
void Sift::FeatureDescrip(vector<ImgRaw>& kaidaImag, Feature* FeatureNow) {
	//新增一個4*4*8的特徵描述空間
	vector <vector <vector <float>>> hist(4);//儲存特徵描述子
	for (int v = 0; v < 4; v++) {
		hist[v].resize(4);
		for (int j = 0; j < 4; j++) {
			hist[v][j].resize(8); // 8個角度
		}
	}
	// 從當前極值產出的點開始做
	for (Feature* FeatureS = FeatureNow;
		 FeatureS != nullptr;
		 FeatureS = FeatureS->nextptr)
	{
		//const int radius = (FeatureS->sigmaOCT*3.f * SQUARE2*5.f) * 0.5f;
		ImgRaw& curryImg = kaidaImag[FeatureS->kai];
		int col = curryImg.width;
		int row = curryImg.height;
		int ori = FeatureS->sita;
		int d=4;
		int r=FeatureS->y;
		int c=FeatureS->x;
		int n=SIFT_DESCR_HIST_BINS;
		int scl_octv = FeatureS->sigmaOCT;


		// 描述特徵子
		hist = descr_hist(curryImg, col, row, r, c, ori, scl_octv, d, n);
			// rob規一化方法
			hist_to_descr(hist, d, n, FeatureS);
			// 舊規一化方法(不知道為什麼效果比較好)
			Sift::DescripNomal(hist);
			FeatureS->descrip = hist;
	}
	// cout << endl;
}
// 描述特徵子 (舊方法有誤)
static inline void CoordinateChange(int& deltaX, int& deltaY, float sita) {
	// 座標轉換
	float deltaX2 = deltaY*sin(sita) + deltaX*cos(sita);
	deltaY = deltaY*cos(sita) - deltaX*sin(sita);
	deltaX = deltaX2;

	//float r_rot = (j)*sin_t + (i)*cos_t; // 原圖X座標
	//float c_rot = (j)*cos_t - (i)*sin_t; // 原圖Y座標
}
void Sift::FeatureDescrip_ori(vector<ImgRaw>& kaidaImag, Feature* FeatureNow) {
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
		int radius = (FeatureS->sigmaOCT * 3 * sqrt(2)*(5)) / 2;
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
				mm = sqrtf(dx*dx + dy*dy);
				sita = fastAtan2f(dy, dx);

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
		Sift::DescripNomal(descripgroup);
		FeatureS->descrip = descripgroup;
	}
}
void Sift::DescripNomal(Desc& descripgroup) {
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
}
// 印出箭頭
void Sift::drawArrow(string name, float ratio){
	// 臨時圖庫
	ImgRaw img = raw_img;
	for (Feature* feaPoint = FeatStart;
		feaPoint->nextptr != NULL;
		feaPoint = feaPoint->nextptr)
	{
		size_t x = feaPoint->x/feaPoint->size;
		size_t y = feaPoint->y/feaPoint->size;
		Draw::draw_arrowRGB(img, y, x, (feaPoint->mm)*ratio, feaPoint->sita);
	}
	/*輸出圖片*/
	img.bmp(name);
}

/*
    ##     ##    ###    ########  ######  ##     ##
    ###   ###   ## ##      ##    ##    ## ##     ##
    #### ####  ##   ##     ##    ##       ##     ##
    ## ### ## ##     ##    ##    ##       #########
    ##     ## #########    ##    ##       ##     ##
    ##     ## ##     ##    ##    ##    ## ##     ##
    ##     ## ##     ##    ##     ######  ##     ##
*/
Stitching::Stitching(const Sift& desc1, const Sift& desc2):
	feat1(desc1.FeatStart), feat2(desc2.FeatStart)
{
	const ImgRaw& img1 = desc1.raw_img;
	const ImgRaw& img2 = desc2.raw_img;
	this->Width  = img1.width+img2.width;
	this->Height = img1.height;
	// 合併兩張圖
	this->matchImg.resize(img1.width*2, img2.height, 24);
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
float Stitching::EuclDist2(float point1[128], float point2[128]) {
	float sum = 0.f;
	for (size_t i = 0; i < 128; i++) {
		const float reduce = point1[i] - point2[i];
		sum += reduce * reduce;
	}
	return sqrtf(sum);
}
// 匹配相同的特徵點
void Stitching::Check(float matchTh) {
	for (Feature* startlink1 = feat1->nextptr;
		startlink1 != nullptr;
		startlink1 = startlink1->nextptr)
	{
		// 初始化最大值
		float dist1, dist2;
		dist1 = dist2 = numeric_limits<float>::max();
		float useX1, useY1;
		useX1 = useY1 = numeric_limits<float>::max();
		// 找處距離最近的兩個點
		for (Feature* startlink2 = feat2->nextptr;
			startlink2 != nullptr;
			startlink2 = startlink2->nextptr)
		{
			//float value = EuclDist2(startlink1->descr, startlink2->descr); // rob方法
			float value = EuclDist(startlink1->descrip, startlink2->descrip); // 舊函式
			if (value < dist1) {
				dist2 = dist1;
				dist1 = value;

				useX1 = startlink2->rX();
				useY1 = startlink2->rY();
			} else if (value < dist2) {
				dist2 = value;
			}
		}
		// 確認是否匹配成功並畫線
		if ((dist1 / dist2) < matchTh) {
			// 閥值越小找到的點越少但越可靠，論文建議 0.4~0.6
			int x1, y1;
			x1 = startlink1->rX();
			y1 = startlink1->rY();
			// 畫線
			Draw::drawLineRGB_p(matchImg, y1, x1, useY1,useX1 + (Width / 2));
		}
	}
	matchImg.bmp("_matchImg.bmp", 24);
}

