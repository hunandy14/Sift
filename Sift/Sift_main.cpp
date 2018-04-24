/*****************************************************************
Name : 
Date : 2017/07/04
By   : CharlotteHonG
Final: 2017/07/04
*****************************************************************/
#pragma warning(disable : 4819)
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <memory>
#include <timer.hpp>

using namespace std;
#include "Sift.hpp"
#include <opencv2/opencv.hpp>
using namespace cv;

#if defined(_MSC_VER)
	#define or ||
	#define and &&
	#define OR ||
	#define AND &&
#endif
#define M_PI 3.14159265358979323846
#define SQUARE2 1.4142135623730951f


void StitchImg(string  name1, string name2, string path="testImg") {
	ImgRaw img1(name1, path);
	ImgRaw img2(name2, path);


    // 金字塔1
	Timer t1;
	Sift fea1(img1);
	fea1.pyramid();
	//fea1.drawArrow("feaArrow1.bmp");
	t1.print("feat1");

	// 金字塔2
	t1.start();
	Sift fea2(img2);
	fea2.pyramid();
	//fea2.drawArrow("feaArrow2.bmp");
	t1.print("feat2");


	// 匹配特徵點(兩張大小要一樣)
	Stitching match(fea1, fea2);
	match.Check(0.5);
	return;
}
//================================================================
int main(int argc, char const *argv[]){

	//srand((unsigned)time(NULL)); rand();

	//string name1="kanna.bmp", name2="kanna2.bmp";
	//string name1="ball_01.bmp", name2="ball_02.bmp";
	//string name1="ball_03.bmp", name2="ball_04.bmp";
	//string name1="sd01.bmp", name2="sd02.bmp";
	//string name1="d01.bmp", name2="d02.bmp";
	//string name1="kanna.bmp", name2="kanna2.bmp";
	//string name1="kanna_L.bmp", name2="kanna_R.bmp";
	//string name1="sm01.bmp", name2="sm02.bmp";
	//string name1="mm01.bmp", name2="mm02.bmp";
	//string name1="mm03.bmp", name2="mm04.bmp";
	//string name1="mm05.bmp", name2="mm06.bmp";
	//string name1="mm07.bmp", name2="mm08.bmp";
	//string name1="mori1.bmp", name2="mori2.bmp";
	//string name1="morii1.bmp", name2="morii2.bmp";
	//string name1="mori2.bmp", name2="mori3.bmp";

	//string name1="l01.bmp", name2="l02.bmp";
	//string name1="dk05.bmp", name2="dk06.bmp";
	//string name1="dk13.bmp", name2="dk14.bmp";
	//string name1="sc01.bmp", name2="sc02.bmp";
	string name1="sc02.bmp", name2="sc03.bmp";
	
	//string name1="lib01.bmp", name2="lib02.bmp";
	//string name1="lib03.bmp", name2="lib04.bmp";

	//string name1="b01.bmp", name2="b02.bmp";
	//string name1 = "tf01.bmp", name2= "tf02.bmp";
	//string name1 = "tf11.bmp", name2= "tf12.bmp";
	//string name1 = "tfv01.bmp", name2= "tfv02.bmp";
	//string name1 = "tfv03.bmp", name2= "tfv04.bmp";
	
	//StitchImg(name1, name2);
	StitchImg("DSC_2936.bmp", "DSC_2937.bmp", "data");
	StitchImg("DSC_2938.bmp", "DSC_2939.bmp", "data");
	StitchImg("DSC_2940.bmp", "DSC_2941.bmp", "data");
	StitchImg("DSC_2942.bmp", "DSC_2943.bmp", "data");
	StitchImg("DSC_2944.bmp", "DSC_2945.bmp", "data");
	StitchImg("DSC_2946.bmp", "DSC_2947.bmp", "data");
	StitchImg("DSC_2950.bmp", "DSC_2951.bmp", "data");

	return 0;
}

//================================================================
