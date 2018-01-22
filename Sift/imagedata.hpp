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
	float size;//��
	int kai;//�h
	float sigmaOCT;//�����ҽk�Y��
	int x, y;//�U�Ҧb���h���y��
	float mm;//�j��
	int sita;//�]�t�D��V�P�t��V������
	vector<vector<vector<float>>> descrip;// �y�z�l(�¤�k��)
	Feature* nextptr = nullptr; // �U�@�I
	// rob
	float descr[128] = {};// �έp�����᪺�y�z�l(robbs����k)
	int d = 0; // �S�x�I����
	// �C
	float d_l;// �o�Ӥ����D�F�����A�i��O���եΪ�
	Feature* fwd_match = nullptr;  /**< matching feature from forward image */
	fpoint img_pt;				   /**< location in image �o�ӥu�Oxy��m�A�]�����s�F�N��ΤF*/
	fpoint mdl_pt;				   /**< location in model �bxform�̭��ڧ⥦�� rX rY �������A*/
										/*��ɭԥ��ɹB��i��ݭn"�w"��A���d��*/
	void* feature_data = nullptr;  /**< user-definable data */


public:
	float rX() {return x/size;}
	float rY() {return y/size;}
};

/** �O�s�˴��������ƾ� */
struct detection_data
{
	int r;
	int c;
	int octv;
	int intvl;
	float subintvl;
	float scl_octv;
};