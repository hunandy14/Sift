#ifndef XFORM_H
#define XFORM_H

#include "imagedata.hpp"
/********************************** Structures *******************************/
/** holds feature data relevant to ransac */
struct ransac_data
{
	void* orig_feat_data;
	int sampled;
};
/******************************* Defs and macros *****************************/
/* RANSAC error tolerance in pixels */
#define RANSAC_ERR_TOL 3
/** pessimistic estimate of fraction of inlers for RANSAC */
#define RANSAC_INLIER_FRAC_EST 0.25
/** estimate of the probability that a correspondence supports a bad model */
#define RANSAC_PROB_BAD_SUPP 0.10
/* extracts a feature's RANSAC data */
#define feat_ransac_data(feat) ((struct ransac_data*)(feat)->feature_data)
/***************************** Function Prototypes ***************************/
extern std::vector<float> ransac_xform(struct feature *features, int n, int m, float p_badxform, float err_tol, struct feature*** inliers, int &n_in);
extern std::vector<float> _ransac_xform(struct feature *features, int n, int m, float p_badxform, float err_tol);
#endif
