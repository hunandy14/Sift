#pragma once

/********************************* Structures ********************************/
struct Feature;
/** a node in a k-d tree */
struct kd_node
{
	int ki;                     /**< partition key index */ //關鍵點直方圖方差最大向量系列位置 
	float kv;                   /**< partition key value */ //直方圖方差最大向量系列中最中間模值
	int leaf;                   /**< 1 if node is a leaf, 0 otherwise */
	Feature* features;			/**< features at this node */
	int n;                      /**< number of features */
	kd_node* kd_left;			/**< left child */
	kd_node* kd_right;			/**< right child */
};
/*************************** Function Prototypes *****************************/
extern kd_node* kdtree_build(Feature* features, int n);
extern int kdtree_bbf_knn(kd_node* kd_root, Feature* feat, int k, Feature*** nbrs, int max_nn_chks);
extern void kdtree_release(kd_node* kd_root);
