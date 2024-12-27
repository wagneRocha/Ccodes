#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ADT_vectors.h"

/* *********************************************************************************** */
/* ---------------------------------- TAD int vectors -------------------------------- */
/* *********************************************************************************** */
int *alloc_vec_d(int n){
	
	int *v;
	v = (int*)malloc(n*sizeof(int));
	if (v == NULL){
		printf("Error: Insufficient Memory.\n");
		exit(1);
	} 
	return v;
}
/* *********************************************************************************** */
int *zeros_vec_d(int n){

	int k;
	int *v = alloc_vec_d(n);
	for(k = 0; k < n; k++)
		v[k] = 0;
	return v;
}
/* *********************************************************************************** */
int *ones_vec_d(int n){

	int k;
	int *v = alloc_vec_d(n);
	for(k = 0; k < n; k++)
		v[k] = 1;
	return v;
}
/* *********************************************************************************** */
void vec_e_vec_d(int *v, int *u, int n){

	int k;
	for(k = 0; k < n; k++)
		v[k] = u[k];
}
/* *********************************************************************************** */
void vec_p_vec_d(int *w, int *u, int *v, int n){

	int k;
	for(k = 0; k < n; k++)
		w[k] = u[k] + v[k];
}
/* *********************************************************************************** */
void print_vec_d(int *v, int n){

	int k;
	for(k = 0; k < n; k++)
		printf("[%d]\n", v[k]);
	printf("\n");
}
/* *********************************************************************************** */
void print_vec_d_t(int *v, int n){

	int k;
	printf("[");
	for(k = 0; k < n; k++)
		printf("%d ", v[k]);
	printf("]\n");
}
/* *********************************************************************************** */
int *find_val_vec_d(int *v, int n, int val){

	int *ans = zeros_vec_d(n);
	int k;
	for(k = 0; k < n; k++)
		if(v[k] == val)
			ans[k] = 1;
	return ans;
}
/* *********************************************************************************** */
int max_vec_d(int *v, int n){

	int max = 0;
	
	if((v != NULL) && (n > 0)){
		int k;
		max = v[0];
		for(k = 1; k < n; k++)
			if(max < v[k])
				max = v[k];
	}
	
	return max;
}
/* *********************************************************************************** */
int *find_ones_position_vec_d(int *v, int n, int num1){

	int i, j = 0;
	int *ans = alloc_vec_d(num1);
	
	for(i = 0; i < n; i++)
		if(v[i] == 1){
			ans[j++] = i;
			}
	
	return ans;
}
/* *********************************************************************************** */
int sum_vec_d(int *v, int n){

	int ans = 0;
	int k;
	for(k = 0; k < n; k++)
		ans = ans + v[k];
		
	return ans;
}

/* *********************************************************************************** */
int minAB_d(int a, int b){

	if(a <= b)
		return a;
	else
		return b;
}
/* *********************************************************************************** */
int maxAB_d(int a, int b){

	if(a >= b)
		return a;
	else
		return b;
}
/* *********************************************************************************** */
/* ---------------------------------- TAD double vectors ----------------------------- */
/* *********************************************************************************** */
double *alloc_vec_lf(int n){
	
	double *v;
	v = (double*)malloc(n*sizeof(double));
 	if (v == NULL){
		printf("Error: Insufficient Memory.\n");
		exit(1);
	} 
	return v;
}
/* *********************************************************************************** */
double *zeros_vec_lf(int n){
 
	int k;
	double *v = alloc_vec_lf(n);
	for(k = 0; k < n; k++)
		v[k] = 0.0;
	return v;
}
/* *********************************************************************************** */
double *ones_vec_lf(int n){
 
	int k;
	double *v = alloc_vec_lf(n);
	for(k = 0; k < n; k++)
		v[k] = 1.0;
	return v;
}
/* *********************************************************************************** */
void vec_e_vec_lf(double *v, double *u, int n){

	int k;
	for(k = 0; k < n; k++)
		v[k] = u[k];
}
/* *********************************************************************************** */
void vec_p_vec_lf(double *w, double *u, double *v, int n){

	int k;
	for(k = 0; k < n; k++)
		w[k] = u[k] + v[k];
}
/* *********************************************************************************** */
void vec_m_vec_lf(double *w, double *u, double *v, int n){

	int k;
	for(k = 0; k < n; k++)
		w[k] = u[k] - v[k];
}
/* *********************************************************************************** */
void l_t_vec_lf(double *v, double l, double *u, int n){

	int k;
	for(k = 0; k < n; k++)
		v[k] = l*u[k];
}
/* *********************************************************************************** */
void print_vec_lf(double *v, int n){

	int k;
	for(k = 0; k < n; k++)
		printf("[%.6f]\n", v[k]);
	printf("\n");
}
/* *********************************************************************************** */
double dot_product_lf(double *u, double *v, int n){

	int k;
	double p = 0;
	for(k = 0; k < n; k++)
		p = p + u[k]*v[k];

	return p;
}
/* *********************************************************************************** */
double norm2_lf(double *v, int n){

	return sqrt(dot_product_lf(v, v, n));
}
/* *********************************************************************************** */
void vec_lf_zeros(double *v, int n){

	int i;
	for(i = 1; i < n; i++)
		v[i] = 0.0;
}
/* *********************************************************************************** */
int *find_val_vec_lf(double *v, int n, double val){

	int *ans = zeros_vec_d(n);
	int k;
	for(k = 0; k < n; k++)
		if(fabs(v[k] - val) < 0.000001)
			ans[k] = 1;
	
	return ans;
}
/* *********************************************************************************** */
double min_val_vec_lf(double *v, int n){

	double ans = 0.0;

	if((v != NULL) && (n > 0)){
		int k;
		ans = v[0];
		for(k = 1; k < n; k++)
			if(v[k] <= ans)
				ans = v[k];
	}
	
	return ans;
}
/* *********************************************************************************** */
double max_val_vec_lf(double *v, int n){

	double ans = 0.0;

	if((v != NULL) && (n > 0)){
		int k;
		ans = v[0];
		for(k = 1; k < n; k++)
			if(v[k] >= ans)
				ans = v[k];
	}
	
	return ans;
}
/* *********************************************************************************** */
double sum_vec_lf(double *v, int n){

	double ans = 0;
	int k;
	for(k = 0; k < n; k++)
		ans += v[k];
		
	return ans;
}
/* *********************************************************************************** */
double minAB_lf(double a, double b){

	if((a < b) || (fabs(a - b) < 0.000001))
		return a;
	else
		return b;
}
/* *********************************************************************************** */
double maxAB_lf(double a, double b){

	if((a > b) || (fabs(a - b) < 0.000001))
		return a;
	else
		return b;
}
/* *********************************************************************************** */
int *vec_lf_2_vec_d(double *vec_lf, int n){

	int i;

	int *vec_d = alloc_vec_d(n);
	
	for(i = 0; i < n; i++)
		vec_d[i] = (int) vec_lf[i];
		
	return vec_d;
}
/* *********************************************************************************** */
void swap_lf(double *a, double *b){

	double t = (*a);
	(*a) = (*b);
	(*b) = t;
}
/* *********************************************************************************** */
void swap_d(int *a, int *b){

	int t = (*a);
	(*a) = (*b);
	(*b) = t;
}
/* *********************************************************************************** */
int partition(double array[], int array_index[], int low, int high){
	
	// select the rightmost element as pivot
	double pivot = array[high];
	// pointer for greater element
	int i = (low - 1);
	int j;
	// traverse each element of the array
	// compare them with the pivot
	for(j = low; j < high; j++){
		if((array[j] < pivot) || (fabs(array[j] - pivot) < 0.000001)){
			// if element smaller than pivot is found
			// swap it with the greater element pointed by i
			i++;
			// swap element at i with element at j
			swap_lf(&array[i], &array[j]);
			swap_d(&array_index[i], &array_index[j]);
		}
	}

	// swap the pivot element with the greater element at i
	swap_lf(&array[i+1], &array[high]);
	swap_d(&array_index[i+1], &array_index[high]);
	// return the partition point
	return (i + 1);
}
/* *********************************************************************************** */
void quicksort(double array[], int array_index[], int low, int high){
	if(low < high){
		// find the pivot element such that
		// elements smaller than pivot are on left of pivot
		// elements greater than pivot are on right of pivot
		int p_i = partition(array, array_index, low, high);

		// recursive call on the left of pivot
		quicksort(array, array_index, low, p_i - 1);

		// recursive call on the right of pivot
		quicksort(array, array_index, p_i + 1, high);
	}
}
/* *********************************************************************************** */
/* -------------------------------- TAD 3D double vectors ---------------------------- */
/* *********************************************************************************** */
double *alloc_vec_lf3(){
	
	return alloc_vec_lf(3);
}
/* *********************************************************************************** */
double *zeros_vec_lf3(){
 
	return zeros_vec_lf(3);
}
/* *********************************************************************************** */
void vec_e_vec_lf3(double *v, double *u){

	vec_e_vec_lf(v, u, 3);
}
/* *********************************************************************************** */
void vec_p_vec_lf3(double *w, double *u, double *v){
	
	vec_p_vec_lf(w, u, v, 3);
}
/* *********************************************************************************** */
void vec_m_vec_lf3(double *w, double *u, double *v){
	
	vec_m_vec_lf(w, u, v, 3);
}
/* *********************************************************************************** */
void l_t_vec_lf3(double *v, double l, double *u){
	
	l_t_vec_lf(v, l, u, 3);
}
/* *********************************************************************************** */
void print_vec_lf3(double *v){

	print_vec_lf(v, 3);
}
/* *********************************************************************************** */
double dot_product_lf3(double *u, double *v){
	
	return dot_product_lf(u, v, 3);
}
/* *********************************************************************************** */
double norm2_lf3(double *v){

	return norm2_lf(v, 3);	
}
/* *********************************************************************************** */
void cross_product0_lf3(double *w, double *u, double *v){

	w[0] = u[1]*v[2] - u[2]*v[1];
	w[1] = u[2]*v[0] - u[0]*v[2];
	w[2] = u[0]*v[1] - u[1]*v[0];
}
/* *********************************************************************************** */
double *cross_product_lf3(double *u, double *v){

	double *w = alloc_vec_lf3();
	w[0] = u[1]*v[2] - u[2]*v[1];
	w[1] = u[2]*v[0] - u[0]*v[2];
	w[2] = u[0]*v[1] - u[1]*v[0];

	return w;
}
/* *********************************************************************************** */
void vec_lf3_zeros(double *v){

	vec_lf_zeros(v, 3);
}

/* *********************************************************************************** */
/* ---------------------------------- TAD string vectors ----------------------------- */
/* *********************************************************************************** */
char *alloc_vec_s(int n){
	
	char *str;
	str = (char*)malloc(n*sizeof(char));
 	if (str == NULL){
		printf("Error: Insufficient Memory.\n");
		exit(1);
	} 
	return str;
}
/* *********************************************************************************** */
void blank_string(char *str, int L){

	memset(str, ' ', L);
	str[L] = '\0';
}
/* *********************************************************************************** */
void remove_spaces(char *input, char *output){
	
	int i, j = 0;
	for (i = 0; input[i] != '\0'; i++){
		if (input[i] != ' '){
			output[j] = input[i];
			j++;
		}
	}
	output[j] = '\0';
}
/* *********************************************************************************** */
/* -------------------------------- TAD dcfile vectors --------------------------- */
/* *********************************************************************************** */
dcinputfile *alloc_vec_dcInputFile(int n){
	
	dcinputfile *dc_vec = (dcinputfile*)malloc(n*sizeof(dcinputfile));
 	if (dc_vec == NULL){
		printf("Error: Insufficient Memory.\n");
		exit(1);
	} 
	return dc_vec;
}
/* *********************************************************************************** */
/* ---------------------------------- TAD protein strucute --------------------------- */
/* *********************************************************************************** */
proteinstructure *alloc_vec_proteinstructure(int n){
	
	proteinstructure *protein = (proteinstructure*)malloc(n*sizeof(proteinstructure));
	if (protein == NULL){
		printf("Error: Insufficient Memory.\n");
		exit(1);
	} 
	return protein;
}
/* *********************************************************************************** */
/* ---------------------------------- TAD edges set ---------------------------------- */
/* *********************************************************************************** */
prune_edges_set *alloc_vec_pruneedgesset(int n){
	
	prune_edges_set *pruneeges = (prune_edges_set*)malloc(n*sizeof(prune_edges_set));
	if (pruneeges == NULL){
		printf("Error: Insufficient Memory.\n");
		exit(1);
	} 
	return pruneeges;
}
/* *********************************************************************************** */
ddgp_parameters *alloc_vec_ddgparameters(int n){

	ddgp_parameters *ddgpparam = (ddgp_parameters*)malloc(n*sizeof(ddgp_parameters));
 	if (ddgpparam == NULL){
		printf("Error: Insufficient Memory.\n");
		exit(1);
	} 
	return ddgpparam;
}
