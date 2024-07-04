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

	if(n > 0){
		int k;
		int max = v[0];
		
		for(k = 1; k < n; k++)
			if(max < v[k])
				max = v[k];
		
		return max;
	}
	else
		return 0;
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
