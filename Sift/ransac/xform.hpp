#pragma once
#include "imagedata.hpp"

using std::vector;
vector<float> ransac_xform(Feature *features, int n, 
	int m, float p_badxform, float err_tol, 
	Feature*** inliers, int &n_in);
