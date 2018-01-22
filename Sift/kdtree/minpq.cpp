#include <iostream>
#include <vector>

#include "minpq.hpp"
#include "imagedata.hpp"
#include "utils.hpp"

/************************* Local Function Prototypes *************************/
static void restore_minpq_order(struct pq_node*, int, int);
static void decrease_pq_node_key(struct pq_node*, int, int);
/************************** Local Inline Functions ***************************/
/* returns the array index of element i's parent */
static inline int parent(int i)
{
	return (i - 1) / 2;
}
/* returns the array index of element i's right child */
static inline int right(int i)
{
	return 2 * i + 2;
}
/* returns the array index of element i's left child */
static inline int left(int i)
{
	return 2 * i + 1;
}
/********************** Functions prototyped in minpq.h **********************/
struct min_pq* minpq_init()
{
	struct min_pq *min_pq;
	min_pq = new struct min_pq;
	min_pq->pq_array = new struct pq_node[MINPQ_INIT_NALLOCD];
	min_pq->nallocd = MINPQ_INIT_NALLOCD;
	min_pq->n = 0;

	return min_pq;
}
int minpq_insert(struct min_pq *min_pq, void *data, int key)
{
	int n = min_pq->n;

	/* double array allocation if necessary */
	if (min_pq->nallocd == n)
	{
		min_pq->nallocd = array_double((void**)&min_pq->pq_array, min_pq->nallocd, sizeof(struct pq_node));
		if (!min_pq->nallocd)
		{
			std::cout << "Warning: unable to allocate memory" << std::endl;
			return 1;
		}
	}

	min_pq->pq_array[n].data = data;
	min_pq->pq_array[n].key = INT_MAX;
	decrease_pq_node_key(min_pq->pq_array, min_pq->n, key);
	min_pq->n++;

	return 0;
}
void* minpq_get_min(struct min_pq *min_pq)
{
	if (min_pq->n < 1)
	{
		std::cout << "Warning: PQ empty" << std::endl;
		return NULL;
	}
	return min_pq->pq_array[0].data;
}
void* minpq_extract_min(struct min_pq* min_pq)
{
	void* data;

	if (min_pq->n < 1)
	{
		std::cout << "Warning: PQ empty" << std::endl;
		return NULL;
	}
	data = min_pq->pq_array[0].data;
	min_pq->n--;
	min_pq->pq_array[0] = min_pq->pq_array[min_pq->n];
	restore_minpq_order(min_pq->pq_array, 0, min_pq->n);

	return data;
}
void minpq_release(struct min_pq** min_pq)
{
	if (!min_pq)
	{
		std::cout << "Warning: NULL pointer error" << std::endl;
		return;
	}
	if (*min_pq && (*min_pq)->pq_array)
	{
		delete((*min_pq)->pq_array);
		delete(*min_pq);
		*min_pq = NULL;
	}
}
static void decrease_pq_node_key(struct pq_node* pq_array, int i, int key)
{
	struct pq_node tmp;

	if (key > pq_array[i].key)
		return;

	pq_array[i].key = key;
	while (i > 0 && pq_array[i].key < pq_array[parent(i)].key)
	{
		tmp = pq_array[parent(i)];
		pq_array[parent(i)] = pq_array[i];
		pq_array[i] = tmp;
		i = parent(i);
	}
}
static void restore_minpq_order(struct pq_node* pq_array, int i, int n)
{
	struct pq_node tmp;
	int l, r, min = i;

	l = left(i);
	r = right(i);
	if (l < n)
		if (pq_array[l].key < pq_array[i].key)
			min = l;
	if (r < n)
		if (pq_array[r].key < pq_array[min].key)
			min = r;

	if (min != i)
	{
		tmp = pq_array[min];
		pq_array[min] = pq_array[i];
		pq_array[i] = tmp;
		restore_minpq_order(pq_array, min, n);
	}
}
