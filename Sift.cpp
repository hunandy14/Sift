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

// 合成圖片
void Sift::comp(vector<ImgRaw>& pyrs, string name){
    // 合成圖片
    size_t w = pyrs[0].width, h = pyrs[0].height;
    size_t s2 = pyrs.size();
    ImgRaw big_img(w*s2, h);
    for(unsigned j=0, idx=0; j < h; ++j) {
        for(unsigned n = 0; n < s2; ++n) {
            for(unsigned i = 0; i < w; ++i) {
                big_img[idx++] = pyrs[n][j*w + i];
            }
        }
    }
    // 輸出圖片並開啟
    if(name!="") {
        name += ".bmp";
        big_img.bmp("pic/"+name, 8);
        // system(name.c_str());
    }
};

// 高斯模糊
vector<ImgRaw> Sift::dog_gau(ImgRaw& img, size_t s){
    size_t w = img.width, h = img.height;
    // 初始化
    vector<ImgRaw> pyrs;
    for(unsigned i = 0; i < s; ++i){
        pyrs.emplace_back(w, h);
    }
    // 高斯模糊
    float p = float(1.6);
    ImgRaw::gauBlur(pyrs[0], img, p);
    for(unsigned i = 1; i < s; ++i) {
        p *= (float)pow(2.0, 1.0/s);
        ImgRaw::gauBlur(pyrs[i], pyrs[i-1], p);
    }
    return pyrs;
}

// 高斯金字塔
void Sift::pyramid(size_t s){
    // 長寬不同的圖會出問題
    s += 3;
    size_t octvs = 3;
    //octvs = (size_t)(log(min(raw_img.width, raw_img.height)) / log(2.0)-2);
    pyrs.resize(octvs);
    // 輸入圖
    ImgRaw temp(raw_img, raw_img.width, raw_img.height);
    // 高斯金字塔
    pyrs[0] = dog_gau(temp, s);
    comp(pyrs[0], "Sift-gau_"+to_string(0));
    for(unsigned i = 1; i < octvs; ++i) {
        // 縮小一張圖
        ImgRaw::first(temp, pyrs[i-1][s/2], 0.5);
        // 回傳 s 張模糊圖
        pyrs[i] = dog_gau(temp, s);
        // 合成圖片
        comp(pyrs[i], "Sift-gau_" + to_string(i));
    }
    // 高斯金字塔差分圖
    for(unsigned j = 0; j < octvs; ++j) {
        // 高斯差分
        for(unsigned i = 0; i < pyrs[j].size()-1; ++i) {
            pyrs[j][i] = pyrs[j][i+1] - pyrs[j][i];
        } pyrs[j].erase(--(pyrs[j].end()));
        comp(pyrs[j], "Sift-diff_" + to_string(j));
    }
    // 取得遮罩(沒有邊緣防呆)
    auto getMask = [](vector<types>& mask,
            ImgRaw& img, size_t y, size_t x)
    {
        // 複製遮罩
		mask.resize(9);
		mask[0]= img.at2d(y-1, x-1);
		mask[1]= img.at2d(y-1, x+0);
		mask[2]= img.at2d(y-1, x-1);

		mask[3]= img.at2d(y+0, x-1);
		mask[4]= img.at2d(y+0, x+0);
		mask[5]= img.at2d(y+0, x-1);

		mask[6]= img.at2d(y+1, x-1);
		mask[7]= img.at2d(y+1, x+0);
		mask[8]= img.at2d(y+1, x+1);
    };

    // 尋找 cubic 極值
    vector<types> mask; // 臨時方塊27點找極值
	for (unsigned py = 0; py < pyrs.size() - 1; ++py) { // 
		for (unsigned px = 1; px < pyrs[py].size() - 1; ++px) { // 
			ImgRaw fea(pyrs[py][px].width, pyrs[py][px].height); // 暫存畫布
			// 選定層級找極值(#沒邊緣防呆)
			for (unsigned j = 1; j < pyrs[py][px].height-2; ++j) { // j = 圖片 Y
				for (unsigned i = 1; i < pyrs[py][px].width-2; ++i) { // i = 圖片 X
					// 取得上下層的九宮格
					getMask(mask, pyrs[py][px - 1], j, i);
					getMask(mask, pyrs[py][px], j, i);
					getMask(mask, pyrs[py][px + 1], j, i);

					float max = *std::max_element(mask.begin(), mask.end());
					float min = *std::min_element(mask.begin(), mask.end());

					size_t mid = (mask.size() - 1) / 2;
					if (mask[mid] == max or mask[mid] == min) { // 保留極值
						constexpr float thre = 0.03/2;
						if (mask[mid] > thre or mask[mid] < -thre) { // 保留保留變動大於 0.03
							// 保留角點偵測(#沒邊緣防呆)
							ImgRaw& pic_temp = pyrs[py][px];
							
							fea.at2d(j, i) = 1;
						}
					}
				}
			}
			// 特徵點圖片預覽
			fea.bmp("fea/fea" + to_string(py) + "-" + to_string(px) + ".bmp", 8);
		}
	}



	cout << "end" << endl;
}



