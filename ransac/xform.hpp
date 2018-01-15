#ifndef XFORM_H
#define XFORM_H

#include "imagedata.hpp"
/********************************** Structures *******************************/
/** holds feature data relevant to ransac */
struct ransac_data  
{  
	void* orig_feat_data; //保存此特征点的feature_data域的以前的值  
	int sampled;		  //标识位，值为1标识此特征点是否被选为样本  
};  
/******************************* Defs and macros *****************************/
/* RANSAC error tolerance in pixels */
/*RANSAC算法的容错度 
对于匹配点对<pt,mpt>，以及变换矩阵H， 
如果pt经H变换后的点和mpt之间的距离的平方小于RANSAC_ERR_TOL，则可将其加入当前一致集 */  
#define RANSAC_ERR_TOL 3

//内点数目占样本总数目的百分比的最小值
/** pessimistic estimate of fraction of inlers for RANSAC */
#define RANSAC_INLIER_FRAC_EST 0.25

/** estimate of the probability that a correspondence supports a bad model */
//一个匹配点对支持错误模型的概率（不知道是干什么用的）
#define RANSAC_PROB_BAD_SUPP 0.10

/* extracts a feature's RANSAC data */
//定义了一个带参数的函数宏，用来提取参数feat中的feature_data成员并转换为ransac_data格式的指针 
#define feat_ransac_data(feat) ((ransac_data*)(feat)->feature_data)


/***************************** Function Prototypes ***************************/
std::vector<float> ransac_xform(Feature *features, int n, int m, float p_badxform, float err_tol, Feature*** inliers, int &n_in);
//extern std::vector<float> _ransac_xform(struct feature *features, int n, int m, float p_badxform, float err_tol);
#endif
