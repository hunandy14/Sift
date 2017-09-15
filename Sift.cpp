/*****************************************************************
Name :
Date : 2017/07/05
By   : CharlotteHonG
Final: 2017/07/05
*****************************************************************/
#pragma warning(disable : 4819)
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include "Sift.hpp"
using namespace std;

// 初始化
ImgRaw::ImgRaw(size_t width, size_t height, float val) :
	raw_img(width*height), width(width), height(height)
{
	if (val) {
		std::fill_n (raw_img.begin(), raw_img.size(),val/255.0);
	}
}
ImgRaw::ImgRaw(string bmpname, bool gray_tran){
	vector<unsigned char> img;
	uint32_t width, height;
	uint16_t bits;
	// 讀取圖片
	Raw::read_bmp(img, bmpname, &width, &height, &bits);
	if (gray_tran == 1 and bits == 24)
		Raw::raw2gray(img);
	// 初始化
	raw_img.resize(img.size());
	for (size_t i = 0; i < img.size(); i++)
		raw_img[i] = (float)img[i] / 255.0;
	this->width = width;
	this->height = height;
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
// 高斯金字塔
void Sift::pyramid(size_t s) {
	size_t sacle = s+SIFT_SacleDiff; // 找極值捨去前後兩張 + 差分圖少1張
    size_t octvs = 2; // 不同大小圖片(自訂)
    //octvs = (size_t)(log(min(raw_img.width, raw_img.height)) / log(2.0)-2);
    pyrs.resize(octvs);

    // 輸入圖
    ImgRaw temp(raw_img.width, raw_img.height);
	ImgRaw::first(temp, raw_img, 1);
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
			pyrs_dog[j][i] = pyrs[j][i + 1] - pyrs[j][i];
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
			float gau_r;// 當前這張圖的高斯半徑
			gau_r =  2.f*1.f/powf(2.f, py) * (px+1) * 1.5f * 3.f;
//cout << 2.f*1.f/powf(2.f, py) << ", " << (px+1) << "----";
//cout << gau_r << ", ";
            ImgRaw fea(pyrs_dog[py][px].width, pyrs_dog[py][px].height); // 暫存畫布
			// 選定層級找極值(#防呆去除邊緣)
            for (unsigned j = 1; j < pyrs_dog[py][px].height - 2; ++j) {
                for (unsigned i = 1; i < pyrs_dog[py][px].width - 2; ++i) {
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
								// 移除半徑會超出邊緣的點
								if (j>=gau_r and j<=pyrs_dog[py][px].height-gau_r-1
									and i>=gau_r and i<=pyrs_dog[py][px].width-gau_r-1)
								{
									// 計算特徵點的統計方向與強度
									float fea_m = 0.f;
									float fea_sida = 0.f;
									
									// 建立特徵點的統計方向與強度
									feap.emplace_back(py, px, j, i, gau_r, pyrs_dog[py][px].sigma);
									fea.at2d(j, i) = 1;
								}
								
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


// 測試矩陣
//for (size_t j = 0; j < gr; j++) {
//	for (size_t i = 0; i < gr; i++) {
//		cout << gau_mat.at2d(j, i);
//	} cout << endl;
//}

	ImgRaw gau_2dmat(0, 0); // 高斯矩陣
	ImgRaw fea_Mmat(0, 0);
	ImgRaw fea_Smat(0, 0);
	//for (size_t i = 0; i < feap.size(); i++) {
		int i=0;
		size_t gr = feap[i].gau_r;
		float gs = feap[i].sigma;
		GauBlur::gau_matrix2d(gau_2dmat, gs, gr*2+1);
		auto& cur_img = pyrs_dog[feap[i].o][feap[i].s];

		for (size_t j = 0; j < gr*2+1; j++) {
			for (size_t i = 0; i < gr*2+1; i++) {
				float temp = cur_img.at2d(j, i);
				float m = fea_m(cur_img, );
				float s = fea_m(temp);
			}
		}
	//}

	
} // 高斯金字塔

// 計算特徵描述長度
float Sift::fea_m(ImgRaw& img, size_t y, size_t x) {
	// 檢查邊界
	if (x<=0 or y<=0 or x>= img.width-2 or y>= img.height-2) {
		Out_of_imgRange("超出圖片邊界");
	}
	// 公式
	float m = sqrt(
		pow(img.at2d(x+1, y)-img.at2d(x-1, y), 2)+
		pow(img.at2d(x, y+1)-img.at2d(x, y-1), 2)
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
	float sida = 
		(float)(img.at2d(x, y+1)-img.at2d(x, y-1))/
		(float)(img.at2d(x+1, y)-img.at2d(x-1, y));
	return sida;
}



