/*****************************************************************
Name : 
Date : 2018/01/04
By   : CharlotteHonG
Final: 2018/01/04
*****************************************************************/
#pragma once

struct fpoint{
public:
	fpoint() :x(0.0), y(0.0) {}
	fpoint(float p_x, float p_y) : x(p_x), y(p_y) {}
	float x;
	float y;
};

class Raw {
public:
	Raw() = default;
	Raw(size_t w, size_t h) :
		RGB(w*h * 3), col(w), row(h) {}
	void resize(size_t w, size_t h){
		RGB.resize(w*h*3);
		col=w, row=h;
	}
	int getCol() const { return (int)col; }
	int getRow() const { return (int)row; }
public:
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
	std::vector<float> array;
	size_t col;
	size_t row;
};

struct Feature {
private:
	using Desc = std::vector<std::vector<std::vector<float>>>;
public:
	float size;                     // ��
	int kai;                        // �h
	float sigmaOCT;                 // �����ҽk�Y��
	int x, y;                       // �U�Ҧb���h���y��
	float mm;                       // �j��
	int sita;                       // �]�t�D��V�P�t��V������
	Feature* nextptr = nullptr;     // �U�@���쵲�I
									// �y�z�l����
	float descr[128] = {};          // �έp�����᪺�y�z�l
	int d = 0;                      // �y�z�l����                    
	Feature* fwd_match = nullptr;   // �ǰt�I
	fpoint img_pt;                  // �Y��^��Ϫ��I
	fpoint mdl_pt;                  // 
	void* feature_data = nullptr;   // rob�禡�B��ɪ��Ȧs
									// �C�o
	float d_l;                      // �o�Ӥ����D�F�����A�i��O���եΪ�
public:
	float rX() const {return x/size;}
	float rY() const {return y/size;}
};

// �O�s�˴��������ƾ�
struct detection_data
{
	int r;
	int c;
	int octv;
	int intvl;
	float subintvl;
	float scl_octv;
};


// Blend
struct Blend_Image
{
	int width;
	int height;
	float *RGB;
};