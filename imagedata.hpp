/*****************************************************************
Name : 
Date : 2018/01/04
By   : CharlotteHonG
Final: 2018/01/04
*****************************************************************/
#pragma once
using std::vector;

struct fpoint{
public:
	fpoint() :x(0.0), y(0.0) {}
	fpoint(float p_x, float p_y) : x(p_x), y(p_y) {}
	float x;
	float y;
};
class Raw
{
public:
	Raw() {}
	Raw(size_t w, size_t h) : RGB(w * h * 3, 0), col(w), row(h) {}
	~Raw()
	{
		RGB.clear();
	}
	int getCol()
	{
		return (int)col;
	}
	int getRow()
	{
		return (int)row;
	}
	std::vector<unsigned char> RGB;
protected:
	size_t col;
	size_t row;
};
class Image_Data {
public:
	Image_Data() {}
	Image_Data(size_t w, size_t h) : col(w), row(h), array(w * h, 0.0) {}
	~Image_Data(){
		array.clear();
	}
	float& operator[](size_t idx){
		return array[idx];
	}
	const float& operator[](size_t idx) const{
		return array[idx];
	}
	size_t getArraysize(){
		return array.size();
	}
	int getCol(){
		return (int)col;
	}
	int getRow(){
		return (int)row;
	}
private:
	vector<float> array;
	size_t col;
	size_t row;
};

struct Feature {
	float size;//階
	int kai;//層
	float sigmaOCT;//高斯模糊係數
	int x, y;//各所在階層的座標
	float mm;//強度
	int sita;//包含主方向與負方向的角度
	vector<vector<vector<float>>> descrip;// 描述子(舊方法的)
	Feature* nextptr = nullptr; // 下一點
	// rob
	float descr[128] = {};// 統計完成後的描述子(robbs的方法)
	int d = 0; // 特徵點長度
	// 顏
	float d_l;// 這個不知道幹嘛的，可能是測試用的
	Feature* fwd_match = nullptr;  /**< matching feature from forward image */
	fpoint img_pt;				   /**< location in image 這個只是xy位置，因為有存了就棄用了*/
	fpoint mdl_pt;				   /**< location in model 在xform裡面我把它用 rX rY 替換掉，*/
										/*到時候平時運算可能需要"預"算，先留著*/
	void* feature_data = nullptr;  /**< user-definable data */


public:
	float rX() {return x/size;}
	float rY() {return y/size;}
};

/** 保存檢測的相關數據 */
struct detection_data
{
	int r;
	int c;
	int octv;
	int intvl;
	float subintvl;
	float scl_octv;
};