
#include <iostream>
#include <vector>
#include "imagedata.hpp"
#include "minpq.hpp"
#include "kdtree.hpp"

using namespace std;

struct bbf_data
{
	float d;
	void* old_data;
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
	Feature *tree_feat, **_nbrs;
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
		bbf_data = (struct bbf_data*)_nbrs[i]->feature_data;
		_nbrs[i]->feature_data = bbf_data->old_data;
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
	kd_node* kd_node;

	kd_node = new struct kd_node;
	memset(kd_node, 0, sizeof(kd_node));
	kd_node->ki = -1;
	kd_node->features = features;
	kd_node->n = n;

	return kd_node;
}
// 擴展子樹
static void expand_kd_node_subtree(kd_node* kd_node)
{
	/* base case: leaf node */
	if (kd_node->n == 1 || kd_node->n == 0)
	{
		kd_node->leaf = 1;
		return;
	}

	assign_part_key(kd_node);
	
	partition_features(kd_node);
	if (kd_node->kd_left)
		expand_kd_node_subtree(kd_node->kd_left);
	if (kd_node->kd_right)
		expand_kd_node_subtree(kd_node->kd_right);
	
}
// 獲得 ki(位置) 與 kv(中間值)
static void assign_part_key(kd_node* kd_node)
{
	Feature* features;
	float kv, x, mean, var, var_max = 0;
	float *tmp;
	int d, n, i, j, ki = 0;

	features = kd_node->features;// 這裡怎麼是這樣做!

	n = kd_node->n;     // 總共找到幾個特徵點
	d = features[0].d;  // 特徵點內有128個數據
	//d = 128; // 不知道為什麼數值怪怪的，直接寫死

	cout << "n=" << n << ", d=" << d << endl;
	cout << endl;
	/* partition key index is that along which descriptors have most variance */
	for (j = 0; j < d; j++)
	{
		mean = 0;
		var = 0;
		for(i = 0; i < n; i++) {

			// 這裡怎麼是這樣做!特徵點應該是new出來的怎麼用[]取
			// 顏瑞穎在算完之後有重新複製到，新的陣列
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
		//tmp[i];
		//features[i].descr[ki];
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
		return;
	}

	kd_node->kd_left = kd_node_init(features, j + 1);
	kd_node->kd_right = kd_node_init(features + (j + 1), (n - j - 1));
}
static kd_node* explore_to_leaf(kd_node* kd_node, Feature* feat, min_pq* min_pq)
{
	struct kd_node *unexpl, *expl = kd_node;
	float kv;
	int ki;

	while (expl && !expl->leaf)
	{
		ki = expl->ki;
		kv = expl->kv;

		if (ki >= feat->d)
		{
			cout << "Warning: comparing imcompatible descriptors" << endl;
			return NULL;
		}
		if (feat->descr[ki] <= kv)
		{
			unexpl = expl->kd_right;
			expl = expl->kd_left;
		}
		else
		{
			unexpl = expl->kd_left;
			expl = expl->kd_right;
		}

		if (minpq_insert(min_pq, unexpl, (int)fabs(kv - feat->descr[ki])))
		{
			cout << "Warning: unable to insert into PQ" << endl;
			return NULL;
		}
	}

	return expl;
}
static int insert_into_nbr_array(Feature *feat, Feature **nbrs, int n, int k)
{
	bbf_data *fdata, *ndata;
	float dn, df;
	int i, ret = 0;

	if (n == 0)
	{
		nbrs[0] = feat;
		return 1;
	}

	/* check at end of array */
	fdata = (bbf_data*)feat->feature_data;
	df = fdata->d;
	ndata = (bbf_data*)nbrs[n - 1]->feature_data;
	dn = ndata->d;
	if (df >= dn)
	{
		if (n == k)
		{
			feat->feature_data = fdata->old_data;
			delete(fdata);
			return 0;
		}
		nbrs[n] = feat;
		return 1;
	}

	/* find the right place in the array */
	if (n < k)
	{
		nbrs[n] = nbrs[n - 1];
		ret = 1;
	}
	else
	{
		nbrs[n - 1]->feature_data = ndata->old_data;
		delete(ndata);
	}
	i = n - 2;
	while (i >= 0)
	{
		ndata = (bbf_data*)nbrs[i]->feature_data;
		dn = ndata->d;
		if (dn <= df)
			break;
		nbrs[i + 1] = nbrs[i];
		i--;
	}
	i++;
	nbrs[i] = feat;

	return ret;
}
