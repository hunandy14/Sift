#pragma once

/******************************* Defs and macros *****************************/
/* initial # of priority queue elements for which to allocate space */
#define MINPQ_INIT_NALLOCD 512
/********************************** Structures *******************************/
/** an element in a minimizing priority queue */
struct pq_node
{
	void* data;
	int key;
};
/** a minimizing priority queue */
struct min_pq
{
	struct pq_node* pq_array;    /* array containing priority queue */
	int nallocd;                 /* number of elements allocated */
	int n;                       /**< number of elements in pq */
};
/*************************** Function Prototypes *****************************/
extern struct min_pq* minpq_init();
extern int minpq_insert(struct min_pq* min_pq, void* data, int key);
extern void* minpq_get_min(struct min_pq* min_pq);
extern void* minpq_extract_min(struct min_pq* min_pq);
extern void minpq_release(struct min_pq** min_pq);
