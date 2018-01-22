
#include <iostream>
#include <vector>
#include "imagedata.hpp"
#include "minpq.hpp"
#include "kdtree.hpp"

using namespace std;

struct bbf_data
{
	float d;  //此特征点到目标点的欧式距离值=
	void* old_data; //保存此特征点的feature_data域的以前的值=
};
/************************* Local Function Prototypes *************************/
static kd_node* kd_node_init(Feature*, int);
static void expand_kd_node_subtree(kd_node*);
static void assign_part_key(kd_node*);
static float median_select(float *, int);
static float rank_select(float *, int, int);
static void insertion_sort(float *, int);
static int partition_array(float *, int, float);
static void partition_features(kd_node*);
static kd_node* explore_to_leaf(kd_node*, Feature*, min_pq*);
static int insert_into_nbr_array(Feature*, Feature**, int, int);
/******************** Functions prototyped in kdtree.hpp **********************/
kd_node* kdtree_build(Feature* features, int n)
{
	kd_node* kd_root;

	if (!features || n <= 0)
	{
		cout << "Warning: kdtree_build(): no features" << endl;
		return NULL;
	}
	// 初始化填預設數據
	kd_root = kd_node_init(features, n);

	expand_kd_node_subtree(kd_root);
	return kd_root;
}
// 計算歐式距離
float descr_dist_sq(Feature* f1, Feature* f2)
{
	float diff, dsq = 0.0;
	float* descr1, *descr2;
	int i = 0, d = 0;

	d = f1->d;
	if (f2->d != d)
		return FLT_MAX;
	descr1 = f1->descr;
	descr2 = f2->descr;

	for (i = 0; i < d; i++)
	{
		diff = descr1[i] - descr2[i];
		dsq += diff * diff;
	}
	return dsq;
}

int kdtree_bbf_knn(kd_node* kd_root, Feature* feat, int k, Feature*** nbrs, int max_nn_chks)
{
	kd_node *expl;
	min_pq *min_pq;
	Feature *tree_feat, **_nbrs = nullptr;
	bbf_data *bbf_data;
	int i, t = 0, n = 0;

	if (!nbrs || !feat || !kd_root)
	{
		cout << "Warning: NULL pointer error" << endl;
		return -1;
	}

	_nbrs = new Feature*[k];
	min_pq = minpq_init();
	minpq_insert(min_pq, kd_root, 0);
	while (min_pq->n > 0 && t < max_nn_chks)
	{
		expl = (kd_node*)minpq_extract_min(min_pq);
		if (!expl)
		{
			cout << "Warning: PQ unexpectedly empty" << endl;
			goto fail;
		}
		expl = explore_to_leaf(expl, feat, min_pq);
		if (!expl)
		{
			cout << "Warning: PQ unexpectedly empty" << endl;
			goto fail;
		}
		for (i = 0; i < expl->n; i++)
		{
			tree_feat = &expl->features[i];
			bbf_data = new struct bbf_data;
			if (!bbf_data)
			{
				cout << "Warning: unable to allocate memory" << endl;
				goto fail;
			}
			bbf_data->old_data = tree_feat->feature_data;
			bbf_data->d = descr_dist_sq(feat, tree_feat);
			tree_feat->feature_data = bbf_data;
			n += insert_into_nbr_array(tree_feat, _nbrs, n, k);
		}
		t++;
	}

	minpq_release(&min_pq);
	for (i = 0; i < n; i++)
	{
		bbf_data = (struct bbf_data*)_nbrs[i]->feature_data; // 這裡已經是空值
		_nbrs[i]->feature_data = bbf_data->old_data;

		if(_nbrs[i]->feature_data==nullptr) {
			//cout << "_nbrs[i]->feature_data = nullptr" << endl;
		}
		else {
			//cout << "_nbrs[i]->feature_data != nullptr" << endl;
		}

		delete(bbf_data);
	}
	*nbrs = _nbrs;
	return n;

fail:
	minpq_release(&min_pq);
	for (i = 0; i < n; i++)
	{
		bbf_data = (struct bbf_data*)_nbrs[i]->feature_data;
		_nbrs[i]->feature_data = bbf_data->old_data;
		delete(bbf_data);
	}
	delete(_nbrs);
	*nbrs = NULL;
	return -1;
}
void kdtree_release(kd_node* kd_root)
{
	if (!kd_root)
		return;
	kdtree_release(kd_root->kd_left);
	kdtree_release(kd_root->kd_right);
	delete kd_root;
}
/************************ Functions prototyped here **************************/
// 初始化填入預設數據
static kd_node* kd_node_init(Feature* features, int n)
{
	kd_node* kd_node = nullptr;

	kd_node = new struct kd_node;
	memset(kd_node, 0, sizeof(kd_node));
	kd_node->ki = -1;
	kd_node->features = features;
	kd_node->n = n;

	// e04 rob為什麼在這裡沒有歸零，尼瑪出了個bug找好久
	// 我把它改到 partition_features 內減少入次數
	kd_node->kd_left = nullptr;
	kd_node->kd_right = nullptr;

	//cout << "n=" << n << endl;
	//cout << endl;

	return kd_node;
}
// 擴展子樹
static void expand_kd_node_subtree(kd_node* kd_node)
{
	/* base case: leaf node */
	// cout << "n=" << kd_node->n << endl;

	if (kd_node->n == 1 || kd_node->n == 0)
	{
		kd_node->leaf = 1;
		return;
	}

	assign_part_key(kd_node);
	
	partition_features(kd_node);

	
	if (kd_node->kd_left != nullptr){
		expand_kd_node_subtree(kd_node->kd_left);
	}
	if(kd_node->kd_right != nullptr) {
		expand_kd_node_subtree(kd_node->kd_right);
	}
	
}
// 獲得 ki(位置) 與 kv(中間值)
static void assign_part_key(kd_node* kd_node)
{
	Feature* features;
	float kv, x, mean, var, var_max = 0;
	float *tmp;
	int d, n, i, j, ki = 0;

	features = kd_node->features;

	n = kd_node->n;     // 總共找到幾個特徵點
	d = features[0].d;  // 特徵點內有128個數據

	/* partition key index is that along which descriptors have most variance */
	for (j = 0; j < d; j++)
	{
		mean = 0;
		var = 0;
		for(i = 0; i < n; i++) {
			mean += features[i].descr[j];
			//cout << features[i].descr[j] << endl;
		}

		mean /= n;
		for (i = 0; i < n; i++)
		{
			x = features[i].descr[j] - mean;
			var += x * x;
		}
		var /= n;

		if (var > var_max)
		{
			ki = j;
			var_max = var;
		}
	}
	/* partition key value is median of descriptor values at ki */
	tmp = new float[n];
	for(i = 0; i < n; i++) {
		tmp[i] = features[i].descr[ki];
	}
	kv = median_select(tmp, n); 
	delete tmp;

	kd_node->ki = ki;
	kd_node->kv = kv;
}
static float median_select(float *array, int n)
{
	return rank_select(array, n, (n - 1) / 2);
}
static float rank_select(float *array, int n, int r)
{
	float med;
	float *tmp;
	int gr_5, gr_tot, rem_elts, i, j;
	/* base case */
	if (n == 1)
		return array[0];
	/* divide array into groups of 5 and sort them */
	gr_5 = n / 5;
	gr_tot = (int)ceil(n / 5.0);
	rem_elts = n % 5;
	tmp = array;
	for (i = 0; i < gr_5; i++)
	{
		insertion_sort(tmp, 5);
		tmp += 5;
	}
	insertion_sort(tmp, rem_elts);
	tmp = NULL;
	/* recursively find the median of the medians of the groups of 5 */
	tmp = new float[gr_tot];
	for (i = 0, j = 2; i < gr_5; i++, j += 5)
		tmp[i] = array[j];
	if (rem_elts)
		tmp[i++] = array[n - 1 - rem_elts / 2];
	med = rank_select(tmp, i, (i - 1) / 2);
	delete tmp;
	/* partition around median of medians and recursively select if necessary */
	j = partition_array(array, n, med);
	if (r == j)
		return med;
	else if (r < j)
		return rank_select(array, j, r);
	else
	{
		array += j + 1;
		return rank_select(array, (n - j - 1), (r - j - 1));
	}
}
static void insertion_sort(float *array, int n)
{
	float k;
	int i, j;

	for (i = 1; i < n; i++)
	{
		k = array[i];
		j = i - 1;
		while (j >= 0 && array[j] > k)
		{
			array[j + 1] = array[j];
			j -= 1;
		}
		array[j + 1] = k;
	}
}
static int partition_array(float *array, int n, float pivot)
{
	float tmp;
	int p, i, j;

	i = -1;
	for (j = 0; j < n; j++)
		if (array[j] <= pivot)
		{
			tmp = array[++i];
			array[i] = array[j];
			array[j] = tmp;
			if (array[i] == pivot)
				p = i;
		}
	array[p] = array[i];
	array[i] = pivot;

	return i;
}
static void partition_features(kd_node* kd_node)
{
	/*在指定的k-d树节点上划分特征点集 
	使得特征点集的前半部分是第ki维小于枢轴的点，后半部分是第ki维大于枢轴的点 
	*/
	Feature *features, tmp;
	float kv;
	int n, ki, p, i, j = -1;

	features = kd_node->features;
	n = kd_node->n;
	ki = kd_node->ki;
	kv = kd_node->kv;
	for (i = 0; i < n; i++)
	{
		if (features[i].descr[ki] <= kv)
		{
			tmp = features[++j];
			features[j] = features[i];
			features[i] = tmp;
			if (features[j].descr[ki] == kv)
				p = j;
		}
	}
	tmp = features[p];
	features[p] = features[j];
	features[j] = tmp;

	/* if all records fall on same side of partition, make node a leaf */
	if (j == n - 1)
	{
		kd_node->leaf = 1;
		// 修復bug歸零指針 (rob沒歸零)
		kd_node->kd_left = nullptr;
		kd_node->kd_right = nullptr;
		return;
	}

	int nl = j + 1;
	int nr = (n - j - 1);
	// cout << "                  nl=" << nl << ", nr=" << nr << endl;

	kd_node->kd_left = kd_node_init(features, nl);
	kd_node->kd_right = kd_node_init(features+(j + 1), nr);
}
// 搜索路径和优先级队列的生成函数.
static kd_node* explore_to_leaf(kd_node* kd_node, Feature* feat, min_pq* min_pq)
{
	struct kd_node *unexpl; //unexpl中存放着优先级队列的候选特征点.
	struct kd_node *expl = kd_node; // expl为开始搜索节点.
	float kv; //kv是分割维度的数据.
	int ki; //ki是分割维度序号.

	while (expl && !expl->leaf)
	{
		ki = expl->ki; //获得分割节点的ki，kv数据.
		kv = expl->kv;

		if (ki >= feat->d)
		{
			cout << "Warning: comparing imcompatible descriptors" << endl;
			return NULL;
		}
		if (feat->descr[ki] <= kv) //目标特征和分割节点分割维上的数据比较.
		{
			unexpl = expl->kd_right; //小于右子树根节点成为候选节点.
			expl = expl->kd_left; //并进入左子树搜索.
		}
		else
		{
			unexpl = expl->kd_left;  //大于左子树根节点成为候选节点.
			expl = expl->kd_right;   //并进入右子树搜索.
		}
		//将候选节点unexpl根据目标与分割超平面的距离插入到优先级队列中.
		if (minpq_insert(min_pq, unexpl, (int)fabs(kv - feat->descr[ki])))
		{
			cout << "Warning: unable to insert into PQ" << endl;
			return NULL;
		}
	}

	return expl; //返回搜索路径中最后的叶子节点.
}


/*插入一个特征点到最近邻数组，使数组中的点按到目标点的距离升序排列 
参数： 
feat：要插入的特征点，其feature_data域应是指向bbf_data结构的指针，其中的d值时feat和目标点的距离的平方 
nbrs：最近邻数组 
n：已在最近邻数组中的元素个数 
k：最近邻数组元素个数的最大值 
返回值：若feat成功插入，返回1，否则返回0 
*/ 
static int insert_into_nbr_array(Feature* feat, Feature** nbrs, int n, int k )  
{  
    struct bbf_data* fdata, * ndata;//fdata是要插入的点的bbf结构，ndata是最近邻数组中的点的bbf结构  
    double dn, df; //dn是最近邻数组中特征点的bbf结构中的距离值，df是要插入的特征点的bbf结构中的距离值  
    int i, ret = 0;  
//原最近邻数组为空  
    if( n == 0 )  
    {  
        nbrs[0] = feat;  
        return 1;//插入成功，返回1  
    }  
/* check at end of array */  
    fdata = (struct bbf_data*)feat->feature_data;//要插入点的bbf结构  
    df = fdata->d;//要插入的特征点的bbf结构中的距离值  
	//cout << "df=" << df << endl;
    ndata = (struct bbf_data*)nbrs[n-1]->feature_data;//最近邻数组中的点的bbf结构  
    dn = ndata->d;//最近邻数组中最后一个特征点的bbf结构中的距离值  
	//cout << "dn=" << dn << endl;
	//df>=dn，说明要插入在最近邻数组的末尾  
    if( df >= dn )  
    {  
        //最近邻数组中元素个数已达到最大值，无法插入  
        if( n == k )  
        {  
            feat->feature_data = fdata->old_data;//不明白这里是干什么？  
            free( fdata );  
            return 0;//插入失败，返回0  
        }  
        nbrs[n] = feat;//插入到末尾  
        return 1;//插入成功，返回1  
    }  

/* find the right place in the array */  
//运行到此处说明插入位置不在数组末尾
    if( n < k )//最近邻数组中元素个数小于最大值，可插入  
    {  
        nbrs[n] = nbrs[n-1];//原数组最后一个点后移  
        ret = 1;//插入结果设为1  
    }  
    else//最近邻数组中元素个数大于或等于最大值，无法插入，插入结果ret还是0  
    {//其实也不是无法插入，而是最近邻数组中元素个数不变，但值会更新  
        nbrs[n-1]->feature_data = ndata->old_data;  
        free( ndata );  
    }  
    i = n-2;  
    //在最近邻数组中查找要插入的位置  
    while( i >= 0 )  
    {  
        ndata = (struct bbf_data*)nbrs[i]->feature_data;  
        dn = ndata->d;  
        if( dn <= df )//找到插入点  
            break;  
        nbrs[i+1] = nbrs[i];//一次后移  
        i--;  
    }  
    i++;  
    nbrs[i] = feat;//插入  
    return ret;//返回结果  
}  
