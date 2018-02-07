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
#include "timer.hpp"
using namespace std;

#if defined(_MSC_VER)
	#define or ||
	#define and &&
	#define OR ||
	#define AND &&
#endif
#define M_PI 3.14159265358979323846
#define SQUARE2 1.4142135623730951f

#include "opencv2/opencv.hpp"
using namespace cv;

#include "imagedata.hpp"
#include "Sift.hpp"

#include "kdtree.hpp"
#include "xform.hpp"
#include "stitch.hpp"
#include "Blend.hpp"



/*
     ######  #### ######## ########
    ##    ##  ##  ##          ##
    ##        ##  ##          ##
     ######   ##  ######      ##
          ##  ##  ##          ##
    ##    ##  ##  ##          ##
     ######  #### ##          ##
*/
// Sift建構子.
Sift::Sift(ImgRaw img, size_t intvls): raw_img(img), pyWidth(intvls+SIFT_INTVLS_DIFF) {
	FeatStart = new Feature; // 第一點為空
	FeatEnd = FeatStart;
}
// 合成圖片.
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
// 印出箭頭.
void Sift::drawArrow(string name, float ratio){
	ImgRaw img = raw_img;
	for(size_t i = 0; i < feaNum; i++) {
		const Feature& fea = FeatStart[i];
		Draw::draw_arrowRGB(img, fea.rY(), fea.rX(), fea.mm*ratio, fea.sita);
	} img.bmp(name);
}


/* 高斯金字塔. */
static bool findMaxMin(vector<ImgRaw>& gauDog_imgs, size_t intvl_idx, size_t curr_Width, size_t y, size_t x) {
	const float& Val = gauDog_imgs[intvl_idx][y*curr_Width + x];
	// 過濾變化率過小的特徵點
	if (fabs(Val) > SIFT_Dx/2.0) {
		// 找極值
		for (int k = (intvl_idx - 1); k <= (intvl_idx + 1); k++) {
			const ImgRaw& currImg = gauDog_imgs[k];
			for (int j = (y - 1); j <= (y + 1); j++) {
				size_t currIdx = j*curr_Width+x;
				for (int i = - 1; i <= 1; i++) {
					if (Val>0 and Val<currImg[currIdx +i])
						return false;
					if(Val<0 and Val>currImg[currIdx +i])
						return false;
				}
			}
		} return true;
	} return false;
};
static bool harris_corner(const vector<float>& p,
	size_t w, size_t y, size_t x, float r=SIFT_CURV_THR)
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
static size_t getPyramidOctv(size_t width, size_t height) {
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
	size_t octvs = getPyramidOctv(gray_img.width, gray_img.height);
	float curr_size=2.f; // 尺度倍率.
	float curr_sigmaCoef=1.f; // 模糊度倍率.
	for (size_t ph = 0; ph < octvs; ph++, curr_size /= 2.f, curr_sigmaCoef*=2.f) {
		const size_t curr_Width = gray_img.width*curr_size;
		const size_t curr_Height = gray_img.height*curr_size;
		ImgRaw currSize_img(curr_Width ,curr_Height, 8);
		ImgRaw::first(currSize_img, gray_img, curr_size);
		// 高斯模糊.
		vector<ImgRaw> gau_imgs(pyWidth);
		ImgRaw::gauBlur(gau_imgs[0], currSize_img, SIFT_GauSigma * (curr_sigmaCoef));
		//gau_imgs[0].bmp("gau/gau"+to_string(0)+".bmp", 8);
		for (size_t i = 1; i < pyWidth; i++) {
			const float curr_Sigma = SIFT_GauSigma * powf(sqrt(2.0),i) * (curr_sigmaCoef); // 這裡的2是 S+3 的 S
			gau_imgs[i].resize(curr_Width, curr_Height, 8);
			Gaus::GauBlur(gau_imgs[i], gau_imgs[i-1], curr_Width, curr_Height, curr_Sigma);
			//gau_imgs[j].bmp("gau/gau"+to_string(j)+".bmp", 8);
		}
		// 高斯差分
		vector<ImgRaw> gauDog_imgs(pyWidth-1);
		for (size_t i = 0; i < pyWidth-1; i++) {
			gauDog_imgs[i].resize(curr_Width, curr_Height, 8);
			for (size_t idx = 0; idx < curr_Width*curr_Height; idx++) {
				gauDog_imgs[i][idx] = gau_imgs[i][idx] - gau_imgs[i+1][idx];
			}
			//gauDog_imgs[j].bmp("gauDog/gauDog"+to_string(j)+".bmp", 8);
		}
	/* 尋找特徵點. */
		// 紀錄當前指針位置開始的位置.
		Feature* CurrFeature = FeatEnd;
		// 尋找特徵點並累加直方圖.
		for (size_t intvl_idx = 1; intvl_idx < pyWidth-1-1; intvl_idx++) {
			//ImgRaw MaxMin_img(curr_Width, curr_Height, 8);
			//ImgRaw Herris_img(curr_Width, curr_Height, 8);
			// 特徵直方圖累加半徑.
			const size_t r = curr_size * 3. * 1.5 * (intvl_idx+1);
			for (size_t j = r+1; j < curr_Height-r-1; j++) {
				for (size_t i = r+1; i < curr_Width-r-1; i++) {
					// 尋找極值點.
					if (findMaxMin(gauDog_imgs, intvl_idx, curr_Width, j, i)) {
						const ImgRaw& currImg = gauDog_imgs[intvl_idx];
						const float currSigma = SIFT_GauSigma * powf(sqrt(2.0),intvl_idx) * (curr_sigmaCoef);
						if (harris_corner(currImg, curr_Width, j, i, r)) { // r=SIFT_CURV_THR=10
							// 這裡不知道為什麼設置為r效果更好論文所提設置是10.
							Feature* NweFeat = nullptr;
							//NweFeat = interp_extremum(gauDog_imgs, curr_Height, curr_Width, ph, intvl_idx, j, j, pyWidth-3, SIFT_CONTR_THR);
							if(NweFeat) {
								//getHistogramMS(NweFeat, currImg, curr_size, intvl_idx, currSigma, j, j, r);
							}
							getHistogramMS(nullptr, currImg, curr_size, intvl_idx, currSigma, j, i, r);
							//Herris_img[j*curr_Width + j] = 255 /255.0;
						}
						//MaxMin_img[j*curr_Width + j] = 255 /255.0;
					}
				}
			}
			//MaxMin_img.bmp("MaxMin/MaxMin_"+to_string(intvl_idx)+".bmp", 8);
			//Herris_img.bmp("Herris/Herris_"+to_string(intvl_idx)+".bmp", 8);
		}
		// 描述剛剛才找到的那堆新特徵點(當前的點是舊的).
		FeatureDescrip(gau_imgs, CurrFeature->nextptr);
	} // 不同尺度大小的for()



	// 複製到連續的位置上.
	Feature* temp = new Feature[feaNum];
	int i = 0;
	for(Feature *preFea = nullptr, 
		*fea = FeatStart->nextptr; // 跳過第一點空點.
		fea != nullptr;
		fea = fea->nextptr
	){
		memcpy(&temp[i++], fea, sizeof(*fea));
		if(preFea) {
			delete preFea; // 刪除上一點.
			// temp[j-1].nextptr = &temp[j]; // next 改成連續記憶體上的.
		}
		preFea = fea; // 迴圈記住上一點.
	} delete FeatStart;
	//　置換過去.
	FeatStart = temp;
	temp = nullptr;
}


/* 獲取特徵點. */
static Feature* FeatAppend(Feature* NweFeat, int Inx, int Iny, float Insize, int intvl_idx, int sigmaOCT, float Inm, int Insita) {
	Feature* newnode;
	if(NweFeat) {
		newnode = NweFeat;
	}
	else {
		newnode = new Feature;
		// 新增原點
		//newnode->img_pt.x = Inx/Insize;
		//newnode->img_pt.y = Iny/Insize;
	}
	newnode->x = Inx;
	newnode->y = Iny;
	newnode->mm = Inm;
	newnode->sita = Insita;
	newnode->size = Insize;
	newnode->kai = intvl_idx;
	newnode->sigmaOCT = sigmaOCT;
	newnode->nextptr = nullptr;
	return newnode;
}
void Sift::getHistogramMS(Feature* NweFeat, const ImgRaw& img, float Insize, size_t intvl_idx, float sigma, 
	size_t Iny, size_t Inx, size_t Inr)
{
	const size_t matLen = Inr*2 + 1;     // 遮罩長度
	vector<float> mag(matLen * matLen);  // 儲存強度
	vector<float> sita(matLen * matLen); // 儲存角度
	const size_t InWidth = img.width;// 圖片寬度
	//做curr_m、sita計算
	for (int j = (Iny - Inr), idx = 0; j <= (Iny + Inr); j++) {
		for (int i = (Inx - Inr); i <= (Inx + Inr); i++, idx++) {
			const float& dx = img[(j+0)*InWidth + (i+1)] - img[(j+0)*InWidth + (i-1)];
			const float& dy = img[(j+1)*InWidth + (i+0)] - img[(j-1)*InWidth + (i+0)];
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
			Feature* newnode = FeatAppend(NweFeat, Inx, Iny, Insize, intvl_idx, sigma, magSum[i], i*10);
			if(newnode) {
				FeatEnd->nextptr = newnode;// 加入該點
				FeatEnd = newnode; // 更新結尾點的標記
				++feaNum; // 計數幾點
			}
		}
	}
}


/* 描述特徵子(rob方法). */
static bool calc_grad_mag_ori(const ImgRaw& img, int r, int c, float &mag, float &ori) {
	int col = img.width;
	int row = img.height;

	if (r > 0 && r < row - 1 && c > 0 && c < col - 1) {
		float dx = img[(r  )*col + (c+1)] - img[(r  )*col + (c-1)];
		float dy = img[(r-1)*col + (c  )] - img[(r+1)*col + (c  )];
		mag = sqrt(dx*dx + dy*dy);
		//mag = Align_dx > Align_dy ? max(0.875f*Align_dx+0.5f*Align_dy, Align_dx) : max(0.875f*Align_dy + 0.5f*Align_dx, Align_dy);
		ori = atan2(dy, dx);
		return true;
	} else {
		return false;
	}
}
static void interp_hist_entry(Desc &hist, float rbin, float cbin, float obin, float mag, int d, int n) {
	int r0 = (int)floor(rbin);
	int c0 = (int)floor(cbin);
	int o0 = (int)floor(obin);
	float d_r = rbin - r0;
	float d_c = cbin - c0;
	float d_o = obin - o0;

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
Sift::Desc Sift::descr_hist(const ImgRaw &img, int r, int c, float ori, float scl, int d, int n) {
	vector<vector<vector<float>>> hist(d);
	for (int i = 0; i < d; i++) {
		hist[i].resize(d);
		for (int j = 0; j < d; j++) {
			hist[i][j].resize(n);
		}
	}

	float cos_t = cos(ori);
	float sin_t = sin(ori);
	float PI2 = 2.f * M_PI;
	float bins_per_rad = n / PI2;
	float exp_denom = d * d * 0.5f;
	float hist_width = (float)SIFT_DESCR_SCL_FCTR * scl;
	int radius = (int)(hist_width * sqrt(2.f) * (d + 1.f) * 0.5f + 0.5f);

	for (int i = -radius; i <= radius; i++) {
		for (int j = -radius; j <= radius; j++) {
			// 把位置縮放到 [0~4)
			float c_rot = (j * cos_t - i * sin_t) / hist_width;
			float r_rot = (j * sin_t + i * cos_t) / hist_width;
			float rbin = r_rot + d / 2.f - 0.5f;
			float cbin = c_rot + d / 2.f - 0.5f;
			// 確保位置 [0~4) 刪除超出原圖指定半徑外的點
			if (rbin > -1.f  &&  rbin < d  &&  cbin > -1.f  &&  cbin < d) {
				// 計算成功回傳1 並修改數值
				float grad_mag, grad_ori;
				if (calc_grad_mag_ori(img, r + i, c + j, grad_mag, grad_ori)) {
					grad_ori -= ori;
					// 修正0~360
					while (grad_ori < 0.f) {
						grad_ori += PI2;
					}
					while (grad_ori >= PI2) {
						grad_ori -= PI2;
					}
					float obin = grad_ori * bins_per_rad;// 高斯加權分布
					float w = exp(-(c_rot * c_rot + r_rot * r_rot) / exp_denom); // 計算權值
					interp_hist_entry(hist, rbin, cbin, obin, grad_mag * w, d, n); //插補
				}
			}
		}
	}
	return hist;
}
// rob hess normalize.
static void normalize_descr(Feature* feat) {
	float len_sq = 0.f;
	for (int i = 0; i < feat->d; i++) {
		float cur = feat->descr[i];
		len_sq += cur*cur;
	}
	float len_inv = 1.f / sqrt(len_sq);
	for (int i = 0; i < feat->d; i++){
		feat->descr[i] *= len_inv;
	}
}
void Sift::hist_to_descr(const Desc &hist, int d, int n, Feature* feat) {
	int featSize = 0;
	for (int r = 0; r < d; r++) {
		for (int c = 0; c < d; c++) {
			for (int o = 0; o < n; o++) {
				feat->descr[featSize++] = hist[r][c][o];
			}
		}
	}
	normalize_descr(feat);
	for (int i = 0; i < featSize; i++) {
		if (feat->descr[i] > SIFT_DESCR_MAG_THR) {
			feat->descr[i] = (float)SIFT_DESCR_MAG_THR;
		}
	}
	normalize_descr(feat);
	// convert floating-point descriptor to integer valued descriptor
	for (int i = 0; i < featSize; i++) {
		int int_val = (int)((SIFT_INT_DESCR_FCTR) * feat->descr[i]);
		feat->descr[i] = (float)min(255, int_val);
	} feat->d = featSize;
}
// 彥誠 normalize.
static void DescripNomal(Desc& descripgroup, Feature* feat) {
	// 限定門檻值
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
	// 向量正規化
	add = sqrt(add);
	int idx = 0;
	for (int j = 0; j < 4; j++) {
		for (int i = 0; i < 4; i++) {
			for (int v = 0; v < 8; v++) {
				descripgroup[j][i][v] /= add;
				feat->descr[idx++] = descripgroup[j][i][v];
			}
		}
	}
	// 儲存特徵描述子長度
	feat->d = idx;
}
// 描述特徵子.
void Sift::FeatureDescrip(vector<ImgRaw>& kaidaImag, Feature* FeatureNow) {
	//新增一個4*4*8的特徵描述空間.
	vector <vector <vector <float>>> hist(4);//儲存特徵描述子.
	for (int v = 0; v < 4; v++) {
		hist[v].resize(4);
		for (int j = 0; j < 4; j++) {
			hist[v][j].resize(8); // 8個角度.
		}
	}
	// 從當前極值產出的點開始做.
	for (Feature* FeatureS = FeatureNow;
		 FeatureS != nullptr;
		 FeatureS = FeatureS->nextptr)
	{
		// 轉換為 rob hess 函式可用的數據
		ImgRaw& curryImg = kaidaImag[FeatureS->kai];
		int r = FeatureS->y;
		int c = FeatureS->x;
		int ori = FeatureS->sita;
		int scl = FeatureS->sigmaOCT;
		int d = SIFT_DESCR_HIST_RADIUS;
		int n = SIFT_DESCR_HIST_BINS; 

		// 描述這當前一點特徵子.
		hist = descr_hist(curryImg, r, c, ori, scl, d, n);
		// rob hess normalize.
		hist_to_descr(hist, d, n, FeatureS);
		// normalize2 (不知道為什麼，match比較好).
		//DescripNomal(hist, FeatureS);
	}
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
	feat1(desc1.FeatStart), feat2(desc2.FeatStart), 
	feat1_Count(desc1.feaNum), feat2_Count(desc2.feaNum)
{
	img1 = desc1.raw_img;
	img2 = desc2.raw_img;
	this->Width  = img1.width+img2.width;
	this->Height = img1.height;
	// 合併兩張圖
	this->stackImg.resize(img1.width*2, img2.height, 24);
	for (size_t j = 0; j < img1.height; j++) {
		for (size_t i = 0; i < img1.width; i++) {
			stackImg[(j*stackImg.width+i)*3 + 0] = img1[(j*img1.width+i)*3 + 0];
			stackImg[(j*stackImg.width+i)*3 + 1] = img1[(j*img1.width+i)*3 + 1];
			stackImg[(j*stackImg.width+i)*3 + 2] = img1[(j*img1.width+i)*3 + 2];
		}
		for (size_t i = img1.width; i < img2.width+img1.width; i++) {
			stackImg[(j*stackImg.width+i)*3 + 0] = img2[(j*img1.width+i-img1.width)*3 + 0];
			stackImg[(j*stackImg.width+i)*3 + 1] = img2[(j*img1.width+i-img1.width)*3 + 1];
			stackImg[(j*stackImg.width+i)*3 + 2] = img2[(j*img1.width+i-img1.width)*3 + 2];
		}
	}
}


// 匹配相同的特徵點.
static float descr_dist_sq(const Feature* f1, const Feature* f2) {
	int f1d = f1->d;
	if(f2->d != f1d) 
		return FLT_MAX;

	const float* descr1 = f1->descr;
	const float* descr2 = f2->descr;

	float dsq = 0.0;
	for (int i = 0; i < f1d; i++){
		float diff = descr1[i] - descr2[i];
		dsq += diff * diff;
	}
	return dsq;
}
static void forceMatchFeat(Feature* feat1, size_t feat1_Count, Feature* feat2, size_t feat2_Count, float matchTh) {
	for(size_t j = 0; j < feat1_Count; j++) {
		Feature* currFeat1 = &feat1[j];
		Feature* featMatch = nullptr;
		// 初始化最大值.
		float dist1, dist2;
		dist1 = dist2 = numeric_limits<float>::max();
		// 找處距離最近的兩個點.
		for(size_t j = 0; j < feat2_Count; j++) {
			Feature* currFeat2 = &feat2[j];
			// 歐式距離
			float value = 0;
			for (size_t i = 0; i < 128; i++) {
				const float reduce = currFeat1->descr[i] - currFeat2->descr[i];
				value += reduce * reduce;
			}
			// 判定距離
			if (value < dist1) {
				dist2 = dist1;
				dist1 = value;
				featMatch = currFeat2;
			} else if (value < dist2) {
				dist2 = value;
			}
		}
		// 儲存匹配點.
		if ((dist1 / dist2) < matchTh) {
			currFeat1->fwd_match = featMatch;
		}
	}
}
static void kdtreeMatchFeat(Feature* feat1, size_t feat1_Count, Feature* feat2, size_t feat2_Count, float matchTh) {
	kd_node *kd_root = kdtree_build(feat2, feat2_Count);
	for(int i = 0; i < feat1_Count; i++) {
		Feature **nbrs, *feat_st = feat1 + i;
		int k = kdtree_bbf_knn(kd_root, feat_st, 2, &nbrs, KDTREE_BBF_MAX_NN_CHKS);
		if(k == 2) {
			float d0 = descr_dist_sq(feat_st, nbrs[0]);
			float d1 = descr_dist_sq(feat_st, nbrs[1]);
			if(d0 < d1 * matchTh) {
				feat1[i].fwd_match = nbrs[0]; // 把匹配點存入
			}
		} else {
			feat1[i].fwd_match = nullptr;
		} delete[] nbrs;
	}
	kdtree_release(kd_root);
}
static void featDrawLine(string name, const ImgRaw& stackImg, const Feature* feat , size_t featNum) {
	ImgRaw outImg = stackImg;
	for(int i = 0; i < featNum; i++) {
		if(feat[i].fwd_match) {
			/*
			fpoint pt11 = fpoint(round(feat[j].rX()), round(feat[j].rY()));
			fpoint pt22 = fpoint(round(feat[j].fwd_match->rX()), round(feat[j].fwd_match->rY()));
			*/
			fpoint pt11 = fpoint(round(feat[i].fwd_match->rX()), round(feat[i].fwd_match->rY()));
			fpoint pt22 = fpoint(round(feat[i].rX()), round(feat[i].rY()));

			const int& x1 = pt11.x;
			const int& y1 = pt11.y;
			const int& x2 = pt22.x + (outImg.width *.5);
			const int& y2 = pt22.y;
			Draw::drawLineRGB_p(outImg, y1, x1, y2, x2);
		} else {
			// cerr << "feat[j].fwd_match is nullptr" << endl;
		}
	}
	outImg.bmp(name, 24);
}
static void featDrawLine(string name, const ImgRaw& stackImg, Feature const* const* RANfeat , size_t RANfeatNum) {
	ImgRaw outImg = stackImg;
	for(int i = 0; i < RANfeatNum; i++) {
		if(RANfeat[i]->fwd_match) {
			/*
			fpoint pt11 = fpoint(round(RANfeat[j]->rX()), round(RANfeat[j]->rY()));
			fpoint pt22 = fpoint(round(RANfeat[j]->fwd_match->rX()), round(RANfeat[j]->fwd_match->rY()));
			*/
			fpoint pt11 = fpoint(round(RANfeat[i]->fwd_match->rX()), round(RANfeat[i]->fwd_match->rY()));
			fpoint pt22 = fpoint(round(RANfeat[i]->rX()), round(RANfeat[i]->rY()));

			const int& x1 = pt11.x;
			const int& y1 = pt11.y;
			const int& x2 = pt22.x + (outImg.width *.5);
			const int& y2 = pt22.y;
			Draw::drawLineRGB_p(outImg, y1, x1, y2, x2);
		} else {
			cerr << "RANSAC feat[i].fwd_match is nullptr" << endl;
		}
	}
	outImg.bmp(name, 24);
}






// 取出旋轉後的圖片
ImgRaw imgWarpAffine(const ImgRaw& Img, size_t x, size_t y, float radius, float sita) {
	// 新圖長寬半徑
	float maxRadius = radius;
	int rotH = floor(maxRadius);
	int rotW = floor(maxRadius);
	ImgRaw rotate(rotW*2, rotH*2, 8);
	// 預算
	sita *= -1; // 把新圖轉回0度
	float cos_t = cosf(sita*(float)(M_PI/180));  
	float sin_t = sinf(sita*(float)(M_PI/180));
	// 跑新圖
	for (int j = -rotH; j < rotH; j++) {
		for (int i = -rotW; i < rotW; i++) {
			// 輸入新圖座標返回舊圖座標(已0, 0為圓心旋轉)
			float r_rot = (j)*sin_t + (i)*cos_t; // 原圖X座標
			float c_rot = (j)*cos_t - (i)*sin_t; // 原圖Y座標
												 // 矯正回指定位置
			float rbin = r_rot + x; 
			float cbin = c_rot + y;
			// 去除原圖外的點
			if (cbin < Img.height - 1 and cbin > 0) {
				if (rbin < Img.width - 1 and rbin > 0) {
					// 雙線姓插補
					rotate.at2d(j+rotH, i+rotW) = Img.atBilinear(cbin, rbin); 
				}
			}
		}
	}
	return rotate;
}

static void focalsFromHomography(const vector<float> &Hmg, float &f0, float &f1, bool &f0_ok, bool &f1_ok)
{
	float d1, d2; // Denominators. 分母. 
	float v1, v2; // Focal squares value candidates. Focal平方值候選.

	f1_ok = true;
	d1 = Hmg[6] * Hmg[7];
	d2 = (Hmg[7] - Hmg[6]) * (Hmg[7] + Hmg[6]);
	v1 = -(Hmg[0] * Hmg[1] + Hmg[3] * Hmg[4]) / d1;
	v2 = (Hmg[0] * Hmg[0] + Hmg[3] * Hmg[3] - Hmg[1] * Hmg[1] - Hmg[4] * Hmg[4]) / d2;

	if (v1 < v2) {
		swap(v1, v2);
	}
	if (v1 > 0 && v2 > 0) {
		f1 = sqrt(abs(d1) > abs(d2) ? v1 : v2);
	} else if (v1 > 0) {
		f1 = sqrt(v1);
	} else {
		f1_ok = false;
	}

	f0_ok = true;
	d1 = Hmg[0] * Hmg[3] + Hmg[1] * Hmg[4];
	d2 = Hmg[0] * Hmg[0] + Hmg[1] * Hmg[1] - Hmg[3] * Hmg[3] - Hmg[4] * Hmg[4];
	v1 = -Hmg[2] * Hmg[5] / d1;
	v2 = (Hmg[5] * Hmg[5] - Hmg[2] * Hmg[2]) / d2;

	if (v1 < v2) {
		swap(v1, v2);
	}
	if (v1 > 0 && v2 > 0) {
		f0 = sqrt(abs(d1) > abs(d2) ? v1 : v2);
	} else if (v1 > 0) {
		f0 = sqrt(v1);
	} else {
		f0_ok = false;
	}
}
static void alignMatch(
	const Raw &img1, const Raw &img2,
	Feature const* const* good_match, int gm_num,
	vector<int> &x, vector<int> &y, 
	float FL)
{
	int cal_dx = 0;
	int cal_dy = 0;
	// 中間值.
	const float&& mid_x1 = (float)img1.getCol() / 2.f;
	const float&& mid_x2 = (float)img2.getCol() / 2.f;
	const float&& mid_y1 = (float)img1.getRow() / 2.f;
	const float&& mid_y2 = (float)img2.getRow() / 2.f;
	// 先算平方.
	const float& fL1 = FL;
	const float& fL2 = FL;
	const float&& fL1_pow = pow(fL1, 2);
	const float&& fL2_pow = pow(fL2, 2);
	
	Timer t;
	for (int i = 0; i < gm_num-1; i++) {
		Feature const* const& curr_m = good_match[i];
		const float&& imgX1 = curr_m->rX();
		const float&& imgY1 = curr_m->rY();
		const float&& imgX2 = curr_m->fwd_match->rX();
		const float&& imgY2 = curr_m->fwd_match->rY();
		// 圖1
		float theta1 = fastAtanf_rad((imgX2 - mid_x1) / fL1);
		float h1 = imgY2 - mid_y1;
		h1 /= sqrt(pow((imgX2 - mid_x1), 2) + fL1_pow);
		int x1 = (int)round(fL1*theta1 + mid_x1);
		int y1 = (int)round(fL1*h1 + mid_y1);
		// 圖2
		float theta2 = fastAtanf_rad((imgX1 - mid_x2) / fL2);
		float h2 = imgY1 - mid_y2;
		h2 /= sqrt(pow((imgX1 - mid_x2), 2) + fL2_pow);
		int x2 = (int)round(fL2*theta2 + mid_x2 + img1.getCol());
		int y2 = (int)round(fL2*h2 + mid_y2);
		// 累加座標.
		int distX = x2 - x1;
		int distY = img1.getRow() - y1 + y2;
		cal_dx += distX;
		cal_dy += distY;
	}
	//t.print("match time.");

	// 平均座標.
	int avg_dx = (float)cal_dx / (float)(gm_num-1);
	int avg_dy = (float)cal_dy / (float)(gm_num-1);


	// 修正座標(猜測是4捨5入哪裡怎樣沒寫好才變成這樣).
	if(avg_dx % 2 == 0){
		if(avg_dx + 1 <= img1.getCol() && avg_dx + 1 <= img2.getCol()){
			avg_dx += 1;
		} else{
			avg_dx -= 1;
		}
	}
	if(avg_dy % 2 == 0){
		if(avg_dy + 1 <= img1.getRow() && avg_dy + 1 <= img2.getRow()){
			avg_dy += 1;
		} else{
			avg_dy -= 1;
		}
	} else if(avg_dy % 2 == 1){
		avg_dy+=0;
	}

	// 輸出座標.
	x.emplace_back(avg_dx);
	y.emplace_back(avg_dy);



}


void Stitching::Check(float matchTh) {

	Timer t1;
/* kdtree */
	t1.start();
	kdtreeMatchFeat(feat2, feat2_Count, feat1, feat1_Count, matchTh);
	//forceMatchFeat(feat1, feat1_Count, feat2, feat2_Count, matchTh);
	t1.print("KD-tree match time");
	// 畫出連線.
	//featDrawLine("_matchImg_kdLinkImg.bmp", stackImg, feat2, feat2_Count);
/* ransac */
	// 獲得矩陣.
	Feature** RANSAC_feat = nullptr;
	int RANSAC_num = 0;
	t1.start();
	const vector<float> HomogMat = ransac_xform(feat2, feat2_Count, 4, 0.005f, 3.f, &RANSAC_feat, RANSAC_num);
	t1.print("ransac time");
	// 查看矩陣.
	cout << "RANSAC_num=" << RANSAC_num << endl;
	for(size_t j = 0; j < 3; j++) {
		for(size_t i = 0; i < 3; i++) {
			cout << HomogMat[j*3 + i] << ", ";
		} cout << endl;
	} cout << endl;
	// 畫出連線.
	//featDrawLine("_matchImg_RANSACImg.bmp", stackImg, RANSAC_feat, RANSAC_num);
// 縫合圖片.
	//------------------------------------------------------------------------
	// 轉換用函式.
	auto&& raw_to_imgraw = [](const Raw& src){
		ImgRaw dst(src.RGB, src.getCol(), src.getRow(), 24);
		return dst;
	};
	auto&& imgraw_to_raw = [](const ImgRaw& src){
		Raw dst(src.width, src.height);
		dst.RGB = src; // 這裡會呼叫重載函式轉uch
		return dst;
	};
	// 宣告所需資料項目.
	vector<Raw> InputImage = {
		imgraw_to_raw(img1),
		imgraw_to_raw(img2)
	};
	//------------------------------------------------------------------------
	// WarpPerspective
	// todo 這裡的長寬沒有算出來，有空讓他自動算出完整的大小
	Raw warpImg(InputImage[1].getCol()*1.6, InputImage[1].getRow()*1.1);
	//Raw warpImg(InputImage[1].getCol(), InputImage[1].getRow());
	_WarpPerspective(InputImage[1], warpImg, HomogMat);

	ImgRaw warpImg2 = raw_to_imgraw(warpImg);
	//warpImg2.bmp("_warpImg2.bmp");

	ImgRaw matchImg=warpImg2;
	for(size_t j = 0; j < img1.height; j++) {
		int start = img1.width;
		int end = img1.width;
		for(size_t i = 0; i <= (img1.width-1); i++) {
			if(warpImg2[j*warpImg2.width*3 + i*3+0] == 0 and 
				warpImg2[j*warpImg2.width*3 + i*3+1] == 0 and
				warpImg2[j*warpImg2.width*3 + i*3+2] == 0)
			{
				// 這裡要補原圖a的.
				matchImg[j*matchImg.width*3 + i*3+0] = img1[j*img1.width*3 + i*3+0];
				matchImg[j*matchImg.width*3 + i*3+1] = img1[j*img1.width*3 + i*3+1];
				matchImg[j*matchImg.width*3 + i*3+2] = img1[j*img1.width*3 + i*3+2];
			} else {
				// 這裡是重疊處.
				// 比例.
				if(start==end) {
					start=i;
					matchImg[j*matchImg.width*3 + i*3+0] = 100;
					matchImg[j*matchImg.width*3 + i*3+1] = 100;
					matchImg[j*matchImg.width*3 + i*3+2] = 100;
				}
				if(start<end) {
					float len = end-start;
					float ratioR = (i-start)/len;
					float ratioL = 1.0 - ratioR;
					matchImg[j*matchImg.width*3 + i*3+0] = 
						img1[j*img1.width*3 + i*3+0]*ratioL + warpImg2[j*warpImg2.width*3 + i*3+0]*ratioR;
					matchImg[j*matchImg.width*3 + i*3+1] = 
						img1[j*img1.width*3 + i*3+1]*ratioL + warpImg2[j*warpImg2.width*3 + i*3+1]*ratioR;
					matchImg[j*matchImg.width*3 + i*3+2] = 
						img1[j*img1.width*3 + i*3+2]*ratioL + warpImg2[j*warpImg2.width*3 + i*3+2]*ratioR;
				}
				
				
			}
		}
	}
	//matchImg.bmp("_matchImg_Warp.bmp");
	//------------------------------------------------------------------------
	// 獲得焦距(所有圖共用一個).
	int img_total = 2;
	float f0 = 0.f, f1 = 0.f, ft = 0.f;
	bool f0ok = false, f1ok = false;
	vector<float> all_focals;
	if(!HomogMat.empty()) {
		focalsFromHomography(HomogMat, f0, f1, f0ok, f1ok);
		if(f0ok && f1ok) {
			all_focals.push_back(sqrtf(f0 * f1));
		}
	}
	if(all_focals.size() >= img_total - 1) {
		sort(all_focals.begin(), all_focals.end());
		if(all_focals.size() % 2 == 1) {
			ft = all_focals[all_focals.size() / 2];
		} else {
			ft = (all_focals[all_focals.size() / 2 - 1] + all_focals[all_focals.size() / 2]) * 0.5f;
		}
	} else {
		float focals_sum = 0.f;
		focals_sum += img1.size()+img2.size();
		ft = focals_sum / (float)img_total;
	}
	cout << "ft = " << ft << endl;

	// Warping 圓柱投影.
	vector<fpoint> upedge;
	vector<fpoint> downedge;
	vector<Raw> warpingImg;
	warping(InputImage, ft, warpingImg, upedge, downedge);
	cout << "warp=" << warpingImg.size() << endl;

	raw_to_imgraw(warpingImg[0]).bmp("_Warp1.bmp");
	raw_to_imgraw(warpingImg[1]).bmp("_Warp2.bmp");


	//-------------------------------------------------------------------------
	//Align 對齊.
	vector<int> Align_dx, Align_dy; // 第二張圖整張的偏移量
	alignMatch(InputImage[0], InputImage[1], RANSAC_feat, RANSAC_num, Align_dx, Align_dy, ft);
	cout << "dxSize = " << Align_dx.size() << ", dx = " << Align_dx[0] << endl;
	cout << "dySize = " << Align_dy.size() << ", dy = " << Align_dy[0] << endl;
	

	//-------------------------------------------------------------------------
	// Blend 多頻段混合.
	for(int num = 0; num < warpingImg.size() - 1; num++) {
		if(Align_dy[num] > warpingImg[num].getRow()) { // 假如 y 的偏移量大於圖片高
			int dyy = -(warpingImg[num].getRow() - abs(warpingImg[num].getRow() - Align_dy[num]));
			cout << "dy--->dyy = " << Align_dy[num] << ", " << dyy << endl;
			multiBandBlend(warpingImg[num], warpingImg[num + 1], Align_dx[num], dyy);
		} else { // 通常情況
			multiBandBlend(warpingImg[num], warpingImg[num + 1], Align_dx[num], Align_dy[num]);
		}
	}
	// 輸出混合後的個別兩張圖.
	raw_to_imgraw(warpingImg[0]).bmp("_WarpBlend1.bmp");
	raw_to_imgraw(warpingImg[1]).bmp("_WarpBlend2.bmp");


	//-------------------------------------------------------------------------
	// Writing 輸出圖片.
	int cols = 0;
	int rows = 0;
	for(int i = 0; i < warpingImg.size(); i++) {
		if(i != warpingImg.size() - 1) {
			if(Align_dy[i] > warpingImg[i].getRow()) {
				rows += warpingImg[i].getRow() - warpingImg[i].getRow();
			} else
				rows += warpingImg[i].getRow() - Align_dy[i];
		} else {
			rows += warpingImg[i].getRow();
		}

		if(i != warpingImg.size() - 1)
			cols += warpingImg[i].getCol() - Align_dx[i];
		else
			cols += warpingImg[i].getCol();
	}

	//Raw result(cols, rows * 2);
	int rows2 = rows+150;
	vector<unsigned char> result(cols*rows2*3);
	int distancey = warpingImg[0].getRow() / 10;


	double minb = warpingImg[0].getRow() / 2;
	double mina = warpingImg[0].getCol() / 2;
	double maxb = warpingImg[0].getRow() / 2;
	double maxa = warpingImg[0].getCol() / 2;

	// 可能是4個邊界點.
	double upa, upb, downa, downb;
	upa = upb = downa = downb = 0;

	int distance = 0;
	for(int j = distancey; j < distancey + warpingImg[0].getRow(); j++)
	{
		for(int i = distance; i < warpingImg[0].getCol()+distance; i++)
		{
			if(j - distancey >= 0 && j - distancey < warpingImg[0].getRow() && 
				i - distance >= 0 && i - distance < warpingImg[0].getCol())
			{
				if(j < rows2 && j >= 0 && i >= 0 && i <= cols)
				{
					const auto&  currImg = warpingImg[0].RGB;
					const auto&& imgW = warpingImg[0].getCol();
					const auto&& imgY = (j-distancey);
					const auto&& imgX = (i-distance);
					// 寫入新圖.
					result[(j*cols + i)*3 + 0] = currImg[(imgY*imgW + imgX)*3 + 0];
					result[(j*cols + i)*3 + 1] = currImg[(imgY*imgW + imgX)*3 + 1];
					result[(j*cols + i)*3 + 2] = currImg[(imgY*imgW + imgX)*3 + 2];
				}
			}
		}
	}
	//system("pause");
	for(int i = 1; i <= Align_dx.size(); i++) {
		distance = distance + warpingImg[i - 1].getCol() - Align_dx[i - 1];

		if(Align_dy[i - 1] > warpingImg[i - 1].getRow()) {
			distancey = distancey - abs(warpingImg[i - 1].getRow() - Align_dy[i - 1]);
		} else {
			distancey = distancey + abs(warpingImg[i - 1].getRow() - Align_dy[i - 1]);
		}

		const auto& currImg = warpingImg[i];
		const auto& currRgb = warpingImg[i].RGB;
		const auto&& currW = warpingImg[i].getCol();
		const auto&& currH = warpingImg[i].getRow();
		for(int b = distancey; b < distancey+currH; b++) {
			for(int a = distance; a < currW + distance; a++) {
				if(b >= 0) {
					//if (result.RGB[(b * result.getCol() + a) * 3 + 0] == 0 && result.RGB[(b * result.getCol() + a) * 3 + 1] == 0 && result.RGB[(b * result.getCol() + a) * 3 + 2] == 0)
					//{
					//if (b - distancey >= 0 && b - distancey<warpingImg[i].getRow() && a - distance >= 0 && a - distance<warpingImg[i].getCol())
					//{
					//result.RGB[(b * result.getCol() + a) * 3 + 0] = warpingImg[i].RGB[((b - distancey) * warpingImg[i].getCol() + (a - distance)) * 3 + 0];
					//result.RGB[(b * result.getCol() + a) * 3 + 1] = warpingImg[i].RGB[((b - distancey) * warpingImg[i].getCol() + (a - distance)) * 3 + 1];
					//result.RGB[(b * result.getCol() + a) * 3 + 2] = warpingImg[i].RGB[((b - distancey) * warpingImg[i].getCol() + (a - distance)) * 3 + 2];
					//}
					//}


					const auto&& imgY = b-distancey;
					const auto&& imgX = a-distance;
					const auto&& imgIdx = (imgY*currW + imgX)*3;

					if(imgY >= 0 && imgY < currH && 
						a-distance >= 0 && a-distance < currW)
					{
						if( currRgb[imgIdx + 0] != 0 ||
							currRgb[imgIdx + 1] != 0 ||
							currRgb[imgIdx + 2] != 0) {
							result[(b*cols + a)*3 + 0] = currRgb[imgIdx + 0];
							result[(b*cols + a)*3 + 1] = currRgb[imgIdx + 1];
							result[(b*cols + a)*3 + 2] = currRgb[imgIdx + 2];
						}
					}

					if(b < minb) {
						minb = upb = b;
						mina = upa = a;
					}
					if(b>maxb) {
						maxb = b;
						maxa = a;
						downb = b - currW;
						downa = a;
					}
				}
			}
		}
	}

	ImgRaw maru_match(result, cols, rows2, 24);
	maru_match.bmp("__maru_match.bmp");
	//------------------------------------------------------------------------
	// get UpDw limit
	struct posiMaxComp {
		bool operator() (fpoint i,fpoint j) { return i.y>j.y; }
	} posiMaxComp;
	struct posiMinComp {
		bool operator() (fpoint i,fpoint j) { return i.y<j.y; }
	} posiMinComp;

	float posiUp = max_element(upedge.begin(), upedge.end(), posiMaxComp)->y;
	//posiUp = (min_element(upedge.begin(), upedge.end(), posiMinComp)->y);

	float posiDw = max_element(downedge.begin(), downedge.end(), posiMaxComp)->y;

	//posiUp += minb;
	//posiDw += distancey;

	//posiUp=89;
	//posiDw=673;
	cout << "minb=" << minb << endl;
	cout << "distancey=" << distancey << endl;
	cout << "posiUp=" << posiUp << endl;
	cout << "posiDw=" << posiDw << endl;

	
	//------------------------------------------------------------------------
	// draw cut img
	int cutImgsizeW = cols;
	int cutImgsizeH = posiDw - posiUp;
	cutImgsizeH = rows2;

	ImgRaw cutImage(cols, rows2, 24);
	/*for (int j = posiUp; j < posiDw; j++)
	{
		for (int i = 0; i < cutImgsizeW; i++)
		{
			cutImage[((j-posiUp)*cutImgsizeW + i) * 3 + 0] = result[(j*cols + i) * 3 + 0];
			cutImage[((j-posiUp)*cutImgsizeW + i) * 3 + 1] = result[(j*cols + i) * 3 + 1];
			cutImage[((j-posiUp)*cutImgsizeW + i) * 3 + 2] = result[(j*cols + i) * 3 + 2];
		}
	}*/
	for(size_t j = 0; j < cutImgsizeH; j++){
		for(size_t i = 0; i < cutImgsizeW; i++){
			cutImage[(j*cutImgsizeW+i)*3 +0] = maru_match[(j*cols+i)*3 +0];
			cutImage[(j*cutImgsizeW+i)*3 +1] = maru_match[(j*cols+i)*3 +1];
			cutImage[(j*cutImgsizeW+i)*3 +2] = maru_match[(j*cols+i)*3 +2];
		}
	}

	cutImage.bmp("_cutImg.bmp");
}
