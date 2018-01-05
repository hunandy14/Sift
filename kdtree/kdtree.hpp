#pragma once

/********************************* Structures ********************************/
struct Feature;
/** a node in a k-d tree */
struct kd_node
{
	int ki;                     /**< partition key index */ //�����I����Ϥ�t�̤j�V�q�t�C��m 
	float kv;                   /**< partition key value */ //����Ϥ�t�̤j�V�q�t�C���̤����ҭ�
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
