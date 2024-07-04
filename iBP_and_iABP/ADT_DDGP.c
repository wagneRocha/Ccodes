#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ADT_matrices.h"
#include "ADT_vectors.h"
#include "ADT_strings.h"
#include "ADT_stack.h"
#include "ADT_DDGP.h"
#define PI 3.14159265359
#define TWOPI 6.28318530718
#define HALFPI 1.57079632680

/* *********************************************************************************** */
int *adjacent_predecessors_cardinality(int *vertices, int m, int n){
	
	int k, *vec_aux;
	int *ans = zeros_vec_d(n);
	
	for(k = 0; k < n; k++){
		vec_aux = find_val_vec_d(vertices, m, k+1);
		ans[k] = sum_vec_d(vec_aux, m);
		free(vec_aux);
	}
		
	return ans;
}
/* *********************************************************************************** */
void dcinputfile_2_instance(dcinputfile *dc_vec, double ****instance_adjlist, proteinstructure **protein, int **kv, int m, int *n){

	int k, i, *adjpred;
	
	int *vertices = alloc_vec_d(m);
	for(k = 0; k < m; k++)
		vertices[k] = dc_vec[k].i;	
	
	(*n) = max_vec_d(vertices, m);
	
	(*kv) = adjacent_predecessors_cardinality(vertices, m, (*n));
	
	(*instance_adjlist) = alloc_3Darray_lf((*n), (*kv));
	
	(*protein) = alloc_vec_proteinstructure((*n));
	int cardUk, *pos_adjpred_vk;

	(*protein)[0].v_k = dc_vec[0].j;
	(*protein)[0].res_seq = dc_vec[0].rj;
	strcpy((*protein)[0].atom_name, dc_vec[0].aj);
	strcpy((*protein)[0].res_name, dc_vec[0].rtj);
	
	for(k = 1; k < (*n); k++){
		adjpred = find_val_vec_d(vertices, m, k+1);
		cardUk = sum_vec_d(adjpred, m);
		pos_adjpred_vk = find_ones_position_vec_d(adjpred, m, cardUk);

		(*protein)[k].v_k = dc_vec[pos_adjpred_vk[0]].i;
		(*protein)[k].res_seq = dc_vec[pos_adjpred_vk[0]].ri;
		strcpy((*protein)[k].atom_name, dc_vec[pos_adjpred_vk[0]].ai);
		strcpy((*protein)[k].res_name, dc_vec[pos_adjpred_vk[0]].rti);
		
		for(i = 0; i < cardUk; i++){
			(*instance_adjlist)[k][i][0] = dc_vec[pos_adjpred_vk[i]].j;
			(*instance_adjlist)[k][i][1] = dc_vec[pos_adjpred_vk[i]].dl;
			(*instance_adjlist)[k][i][2] = dc_vec[pos_adjpred_vk[i]].du;
		}
		free(adjpred);
		free(pos_adjpred_vk);
	}
	
	free(vertices);
}
/* *********************************************************************************** */
double dij_from_adjlist(double ***adjlist, int *cardUi, int i, int j, char l_or_u[]){

	int pos;
	if(!strcmp(l_or_u, "l"))
		pos = 1;
	else
		pos = 2;

	double dij, *Ui, *Uj;
	int *vec1, num1, *pos_dij;
	
	Ui = get_column_mat_lf(adjlist[i-1], cardUi[i-1], 1);
	vec1 = find_val_vec_lf(Ui, cardUi[i-1], (double) j);
	num1 = sum_vec_d(vec1, cardUi[i-1]);
	if(num1 > 0){
		pos_dij = find_ones_position_vec_d(vec1, cardUi[i-1], num1);
		dij = adjlist[i-1][pos_dij[0]][pos];
		free(pos_dij);
		free(vec1);
		free(Ui);
	}
	else{
		free(Ui);
		free(vec1);
		
		Uj = get_column_mat_lf(adjlist[j-1], cardUi[j-1], 1);
		vec1 = find_val_vec_lf(Uj, cardUi[j-1], (double) i);
		num1 = sum_vec_d(vec1, cardUi[j-1]);
		int *pos_dji = find_ones_position_vec_d(vec1, cardUi[j-1], num1);
		
		dij = adjlist[j-1][pos_dji[0]][pos];
		
		free(Uj);
		free(vec1);
		free(pos_dji);
	}
	
	return dij;
}
/* *********************************************************************************** */
double dij_from_discretizable_edges_set(double ***discretizableEdges, int i, int j, char l_or_u[]){

	int pos;
	if(!strcmp(l_or_u, "l"))
		pos = 1;
	else
		pos = 2;

	double dij, *Ui, *Uj;
	int *vec1, num1, *pos_dij;
	
	Ui = get_column_mat_lf(discretizableEdges[i-1], 3, 1);
	vec1 = find_val_vec_lf(Ui, 3, (double) j);
	num1 = sum_vec_d(vec1, 3);
	if(num1 > 0){
		pos_dij = find_ones_position_vec_d(vec1, 3, num1);
		dij = discretizableEdges[i-1][pos_dij[0]][pos];
		free(pos_dij);
		free(vec1);
		free(Ui);
	}
	else{
		free(vec1);
		free(Ui);
		
		Uj = get_column_mat_lf(discretizableEdges[j-1], 3, 1);
		vec1 = find_val_vec_lf(Uj, 3, (double) i);
		num1 = sum_vec_d(vec1, 3);
		int *pos_dji = find_ones_position_vec_d(vec1, 3, num1);
		
		dij = discretizableEdges[j-1][pos_dji[0]][pos];
		free(pos_dji);
		free(vec1);
		free(Uj);
	}
	
	return dij;
}
/* *********************************************************************************** */
void referential_x1_x2_x3(double ***X, double ***discretationEdges){
	
	double d12 = discretationEdges[1][0][1];//[2][1][2]-1
	double d13 = discretationEdges[2][1][1];//[3][2][2]-1
	double d23 = discretationEdges[2][0][1];//[3][1][2]-1
	
	double d12d12   = d12*d12;
	double d12cosTh = 0.5*(d12d12 + d23*d23 - d13*d13)/d23;
	double d12sinTh = sqrt(d12d12 - d12cosTh*d12cosTh);
	
	(*X)[0][0] = -d12sinTh;
	(*X)[0][1] = d12cosTh - d23;
	(*X)[1][1] = -d23;
}
/* *********************************************************************************** */
void adjlist_2_discretation_edges(double ***instance_adjlist, int *kv, int n, int **cliques, double ****discretationEdges){

	int k, i;
	(*discretationEdges) = alloc_fix3Darray_lf(n, 3, 3);
	
	(*discretationEdges)[1][0][0] = cliques[1][1];
	(*discretationEdges)[1][0][1] = dij_from_adjlist(instance_adjlist, kv, 2, 1, "l");
	(*discretationEdges)[1][0][2] = dij_from_adjlist(instance_adjlist, kv, 2, 1, "u");
	
	(*discretationEdges)[2][0][0] = cliques[2][1];
	(*discretationEdges)[2][0][1] = dij_from_adjlist(instance_adjlist, kv, 3, 2, "l");
	(*discretationEdges)[2][0][2] = dij_from_adjlist(instance_adjlist, kv, 3, 2, "u");
	
	(*discretationEdges)[2][1][0] = cliques[2][2];
	(*discretationEdges)[2][1][1] = dij_from_adjlist(instance_adjlist, kv, 3, 1, "l");
	(*discretationEdges)[2][1][2] = dij_from_adjlist(instance_adjlist, kv, 3, 1, "u");
	
	for(k = 3; k < n; k++)
		for(i = 0; i < 3; i++){
			(*discretationEdges)[k][i][0] = cliques[k][i+1];
			(*discretationEdges)[k][i][1] = dij_from_adjlist(instance_adjlist, kv, k+1, cliques[k][i+1], "l");
			(*discretationEdges)[k][i][2] = dij_from_adjlist(instance_adjlist, kv, k+1, cliques[k][i+1], "u");
		}
}
/* *********************************************************************************** */
void adjlist_2_prune_edges(double ***instance_adjlist, int *kv, int n, int **cliques, prune_edges_set **pruneEdges){

	int i, j;
	double *Ui, dl, du;
	int *vec1i1, *vec1i2, *vec1i3, *vec1i1i2i3, numPd, numId, numPd_i, numId_i;
	
	(*pruneEdges) = alloc_vec_pruneedgesset(n);
	
	(*pruneEdges)[0].cardUkP = 0;
	(*pruneEdges)[0].cardUkI = 0;
	(*pruneEdges)[1].cardUkP = 0;
	(*pruneEdges)[1].cardUkI = 0;
	(*pruneEdges)[2].cardUkP = 0;
	(*pruneEdges)[2].cardUkI = 0;
	
	for(i = 3; i < n; i++)	
		if(kv[i] > 3){
			numPd = 0;
			numId = 0;
			Ui = get_column_mat_lf(instance_adjlist[i], kv[i], 1);
			vec1i1 = find_val_vec_lf(Ui, kv[i], cliques[i][1]);
			vec1i2 = find_val_vec_lf(Ui, kv[i], cliques[i][2]);
			vec1i3 = find_val_vec_lf(Ui, kv[i], cliques[i][3]);
			vec1i1i2i3 = alloc_vec_d(kv[i]);
			vec_p_vec_d(vec1i1i2i3, vec1i1, vec1i2, kv[i]);
			vec_p_vec_d(vec1i1i2i3, vec1i1i2i3, vec1i3, kv[i]);
			free(vec1i1);
			free(vec1i2);
			free(vec1i3);
			
			for(j = 0; j < kv[i]; j++)
				if(!vec1i1i2i3[j]){
					if(fabs(instance_adjlist[i][j][1] - instance_adjlist[i][j][2]) < 0.000001)
						numPd++;
					else
						numId++;
				}
			
			(*pruneEdges)[i].cardUkP = numPd;
			(*pruneEdges)[i].cardUkI = numId;
			if(numPd > 0)
				(*pruneEdges)[i].precise = alloc_mat_lf(numPd, 3);
			if(numId > 0)
				(*pruneEdges)[i].interval = alloc_mat_lf(numId, 3);
			
			numPd_i = 0;
			numId_i = 0;
			
			for(j = 0; j < kv[i]; j++)
				if(!vec1i1i2i3[j]){
					dl = dij_from_adjlist(instance_adjlist, kv, i+1, Ui[j], "l");
					du = dij_from_adjlist(instance_adjlist, kv, i+1, Ui[j], "u");
					if(fabs(dl - du) < 0.000001){
						(*pruneEdges)[i].precise[numPd_i][0] = instance_adjlist[i][j][0];
						(*pruneEdges)[i].precise[numPd_i][1] = instance_adjlist[i][j][1];
						(*pruneEdges)[i].precise[numPd_i][2] = instance_adjlist[i][j][2];
						numPd_i++;
					}
					else{
						(*pruneEdges)[i].interval[numId_i][0] = instance_adjlist[i][j][0];
						(*pruneEdges)[i].interval[numId_i][1] = instance_adjlist[i][j][1];
						(*pruneEdges)[i].interval[numId_i][2] = instance_adjlist[i][j][2];
						numId_i++;
					}
				}
			free(Ui);
			free(vec1i1i2i3);
		}
		else{
			(*pruneEdges)[i].cardUkP = 0;
			(*pruneEdges)[i].cardUkI = 0;
		}
}
/* *********************************************************************************** */
void adjlist_2_graph_parameters(double ***instance_adjlist, int *kv, int n, int **cliques, double ****discretationEdges, prune_edges_set **pruneEdges){

	adjlist_2_discretation_edges(instance_adjlist, kv, n, cliques, &(*discretationEdges));
	adjlist_2_prune_edges(instance_adjlist, kv, n, cliques, &(*pruneEdges));
}
/* *********************************************************************************** */
double kappai_d2(double dii2_2, double di1i2, double dii1_2){

	return (dii2_2 + di1i2*di1i2 - dii1_2)/(2.0*di1i2);
}
/* *********************************************************************************** */
double rhoi2_d2(double dii2_2, double ki){

	return (dii2_2 - ki*ki);
}
/* *********************************************************************************** */
double round_double(double v, int n){

	long int num_digits = pow(10,n);
	return (double) round(num_digits*v)/num_digits;
}
/* *********************************************************************************** */
double abs_torsion_angle_with_constants_d2(double p, double q, double dii3_2){
	
	double tau = 0.0;
	double twoq = 2.0*q;
	double dmin2 = p - twoq;
	double dmax2 = p + twoq;
	
	if((dmin2 < dii3_2) && (dii3_2 < dmax2))
		tau = acos(round_double((p - dii3_2)/twoq, 6));
	else if((dmax2 < dii3_2) || (fabs(dmax2 - dii3_2) < 0.000001))
		tau = PI;

	return tau;
}
/* *********************************************************************************** */
double abs_torsion_angle_with_distances_d2(double d12_2, double d13_2, double d14_2, double d23, double d24_2, double d34_2){
	
	double kbar 	= kappai_d2(d12_2, d23, d13_2);
	double rhobar2 	=  rhoi2_d2(d12_2, kbar);
	double k 	= kappai_d2(d24_2, d23, d34_2);
	double rho2 	=  rhoi2_d2(d24_2, k);
	
	double kmkbar = k - kbar;
	double p = kmkbar*kmkbar + rhobar2 + rho2;
	double q = sqrt(rhobar2*rho2);

	return abs_torsion_angle_with_constants_d2(p, q, d14_2);
}
/* *********************************************************************************** */
double torsion_angle_with_points_d2(double *x1, double *x2, double *x3, double *x4){
	
	double d12_2, d13_2, d14_2, d23, d24_2, d34_2;
	
	d12_2 = d2xixj_lf3(x1, x2);
	d13_2 = d2xixj_lf3(x1, x3);
	d14_2 = d2xixj_lf3(x1, x4);	
	d24_2 = d2xixj_lf3(x2, x4);
	d34_2 = d2xixj_lf3(x3, x4);
	
	d23 = dxixj_lf3(x2, x3);
	
	return sign_torsion_angle(x1, x2, x3, x4)*abs_torsion_angle_with_distances_d2(d12_2, d13_2, d14_2, d23, d24_2, d34_2);
}
/* *********************************************************************************** */
int sign_torsion_angle(double *x1, double *x2, double *x3, double *x4){

	double *v = alloc_vec_lf3();
	double *r = alloc_vec_lf3();
	double *s = alloc_vec_lf3();
	vec_m_vec_lf3(v, x1, x2);
	vec_m_vec_lf3(r, x3, x2);
	vec_m_vec_lf3(s, x4, x2);

	double *sperp = cross_product_lf3(r,v);
	int sign0 = sign_double(dot_product_lf3(sperp, s));
	
	free(v);
	free(r);
	free(s);
	free(sperp);

	return sign0;
}
/* *********************************************************************************** */
double d2xixj_lf3(double *x1, double *x2){ 
	
	double *x;
	double norm2;
	x = alloc_vec_lf3();
	vec_m_vec_lf3(x, x1, x2);
	
	norm2 = dot_product_lf3(x, x);
	
	free(x);
	
	return norm2;
}
/* *********************************************************************************** */
double dxixj_lf3(double *x1, double *x2){ 
	
	double *x;
	double norm;
	x = alloc_vec_lf3();
	vec_m_vec_lf3(x, x1, x2);
	
	norm = norm2_lf3(x);
	
	free(x);
	
	return norm;
}
/* *********************************************************************************** */
int sign_double(double v){
	
	if(v >= 0)
		return 1;
	else
		return -1;
}
/* *********************************************************************************** */
double changeReferential(double phase, double tau_ref1){
	
	double tau_ref0 = phase + tau_ref1;
	
	if(tau_ref0 < -PI)
		tau_ref0 = TWOPI + tau_ref0;
	else if(PI < tau_ref0)
		tau_ref0 = -TWOPI + tau_ref0;
	
	return tau_ref0;
}
/* *********************************************************************************** */
double round_lf(double val, int n){
	
	int i;
	double factor = 1.0;
	
	for(i = 0; i < n; i++)
		factor *= 10.0;
		
	return (double)((long long)(val*factor + 0.5))/factor;
}
/* *********************************************************************************** */
void cap2intervals(double *interA, double *interB, int *pos_interC, double ***interC, double myZero){
	
	double aA, bA, aB, bB;
	
	// aA <= aB always!
	if(interA[0] < interB[0]){
		aA = interA[0];
		bA = interA[1];
		aB = interB[0];
		bB = interB[1];
	}
	else{
		aA = interB[0];
		bA = interB[1];
		aB = interA[0];
		bB = interA[1];
	}
	
	int aA_e_bA = (fabs(bA - aA) < myZero) ? 1 : 0;
	int aB_e_bB = (fabs(bB - aB) < myZero) ? 1 : 0;
	int bA_e_aB = (fabs(bA - aB) < myZero) ? 1 : 0;
	
	if((aB < bA) || bA_e_aB){
		if(!aA_e_bA){
			(*interC)[(*pos_interC)][0] = aB;
			if(aB_e_bB || bA_e_aB)
				(*interC)[(*pos_interC)][1] = aB;
			else
				(*interC)[(*pos_interC)][1] = minAB_lf(bA, bB);
		}
		else{
			(*interC)[(*pos_interC)][0] = aA;
			(*interC)[(*pos_interC)][1] = aA;
		}
		(*pos_interC)++;
	}
}
/* *********************************************************************************** */
void interIntervals(ddgp_interval *tau, ddgp_interval tau_k, double myZero){
	
	int k, l;
	double **cap_res;
	
	int num_inter, num_max_inter;
	
	// positive interval intersection
	if((*tau).num_pos_inter != 0){
		num_max_inter = (*tau).num_pos_inter + tau_k.num_pos_inter;
		cap_res = alloc_mat_lf(num_max_inter, 2);
		num_inter = 0;
		for(k = 0; k < (*tau).num_pos_inter; k++)
			for(l = 0; l < tau_k.num_pos_inter; l++)
				if(fabs((tau_k.pos_inter[l][1] - tau_k.pos_inter[l][0]) - PI) < myZero){
					cap_res[num_inter][0] = (*tau).pos_inter[k][0];
					cap_res[num_inter][1] = (*tau).pos_inter[k][1];
					num_inter++;
				}
				else
					cap2intervals(&(*tau).pos_inter[k][0], &tau_k.pos_inter[l][0], &num_inter, &cap_res, myZero);
		
		dealloc_mat_lf((*tau).pos_inter, (*tau).num_pos_inter);
		
		if(num_inter > 0){
			(*tau).num_pos_inter = num_inter;
			(*tau).pos_inter = alloc_mat_lf(num_inter, 2);
			l = 0;
			for(k = 0; k < num_inter; k++){
				(*tau).pos_inter[l][0] = cap_res[l][0];
				(*tau).pos_inter[l][1] = cap_res[l][1];
				l++;
			}
		}
		else
			(*tau).num_pos_inter = 0;
		
		dealloc_mat_lf(cap_res, num_max_inter);
	}
	
	// negative interval intersection
	if((*tau).num_neg_inter != 0){
		num_max_inter = (*tau).num_neg_inter + tau_k.num_neg_inter;
		cap_res = alloc_mat_lf(num_max_inter, 2);
		num_inter = 0;
		for(k = 0; k < (*tau).num_neg_inter; k++)
			for(l = 0; l < tau_k.num_neg_inter; l++)
				if(fabs((tau_k.neg_inter[l][1] - tau_k.neg_inter[l][0]) - PI) < myZero){
					cap_res[num_inter][0] = (*tau).neg_inter[k][0];
					cap_res[num_inter][1] = (*tau).neg_inter[k][1];
					num_inter++;
				}
				else
					cap2intervals(&(*tau).neg_inter[k][0], &tau_k.neg_inter[l][0], &num_inter, &cap_res, myZero);
		
		dealloc_mat_lf((*tau).neg_inter, (*tau).num_neg_inter);
		
		if(num_inter > 0){
			(*tau).num_neg_inter = num_inter;
			(*tau).neg_inter = alloc_mat_lf(num_inter, 2);
			l = 0;
			for(k = 0; k < num_inter; k++){
				(*tau).neg_inter[l][0] = cap_res[l][0];
				(*tau).neg_inter[l][1] = cap_res[l][1];
				l++;
			}
		}
		else
			(*tau).num_neg_inter = 0;
		
		dealloc_mat_lf(cap_res, num_max_inter);
	}
}
/* *********************************************************************************** */
void sample_interval_0(double **J, int num_inter, int *vec_pos, double **sample, int sample_size, double tolerance){
	
	int j, k, l, num_points, n_points;
	int *vec_ni, *vec_ni_sample;
	double *vec_wi;
	double maxW, h, step;
	
	// analise
	//tolerance = 1000.0*tolerance;
	
	(*vec_pos) = sample_size;
	if(num_inter > 0){
		vec_wi = zeros_vec_lf(num_inter);
		maxW = 0.0;
		for(k = 0; k < num_inter; k++){
			vec_wi[k] = J[k][1] - J[k][0];
			if(maxW < vec_wi[k])
				maxW = vec_wi[k];
		}
		vec_ni = ones_vec_d(num_inter);
		for(k = 0; k < num_inter; k++)
			if(vec_wi[k] > tolerance)
				vec_ni[k] = maxAB_d((int) round(((double) sample_size)*vec_wi[k]/maxW), 1);
			
		num_points = sum_vec_d(vec_ni, num_inter);
		vec_ni_sample = zeros_vec_d(num_inter);
		
		if(num_points <= sample_size)
			vec_e_vec_d(vec_ni_sample, vec_ni, num_inter);
		else{
			int pos;
			int pos_j = vec_ni[0]-1;
			l = 0;
			j = 0;
			step = ((double) (num_points - 1))/((double) (sample_size - 1));
			for(k = 0; k < sample_size; k++){
				pos = (int) round(k*step);
				if(pos <= pos_j)
					vec_ni_sample[j]++;
				else{
					l++;
					pos_j += vec_ni[l];
					j++;
					while(pos > pos_j){
						//printf("num_inter = %d, l = %d\n", num_inter, l);
						l++;
						pos_j += vec_ni[l];
						j++;
					}
					vec_ni_sample[j]++;
				}
			}
		}
		
		(*sample) = zeros_vec_lf(sample_size);
		n_points = num_points - sample_size;
		if((n_points < 0) || (num_inter == sample_size)){
			(*vec_pos) = -n_points;
			k = (*vec_pos);
			for(l = 0; l < num_inter; l++)
				(*sample)[k++] = J[l][0];
		}
		else{
			(*vec_pos) = 0;
			j = 0;
			for(k = 0; k < num_inter; k++){
				if(vec_ni_sample[k] == 1)
					(*sample)[j++] = (J[k][0] + J[k][1])/2.0;
					
				else if(vec_ni_sample[k] > 1){
					h = vec_wi[k]/((double) (vec_ni_sample[k] - 1));
					for(l = 0; l < vec_ni_sample[k]; l++)
						(*sample)[j++] = J[k][0] + h*l;
				}
			}
		}
		free(vec_wi);
		free(vec_ni);
		free(vec_ni_sample);
	}
}
/* *********************************************************************************** */
void sampleInterval(int i, double ***branches, int **branchNum, ddgp_interval tau, int sampleSize, double tolerance){
	
	int k, l;
	int branchNum_pos, branchNum_neg;
	double *sample, *sample_pos, *sample_neg;
	int twoSampleSize, num_points_neg, num_points_pos;
	
	// positive interval
	sample_interval_0(tau.pos_inter, tau.num_pos_inter, &branchNum_pos, &sample_pos, sampleSize, tolerance);
	// negative interval
	sample_interval_0(tau.neg_inter, tau.num_neg_inter, &branchNum_neg, &sample_neg, sampleSize, tolerance);
	
	num_points_neg = sampleSize - branchNum_neg;
	num_points_pos = sampleSize - branchNum_pos;
	
	sample = zeros_vec_lf(num_points_neg + num_points_pos);
	l = 0;
	
	for(k = branchNum_neg; k < sampleSize; k++)
		sample[l++] = sample_neg[k];
		
	for(k = branchNum_pos; k < sampleSize; k++)
		sample[l++] = sample_pos[k];
	
	if(num_points_pos > 0)
		free(sample_pos);
	if(num_points_neg > 0)
		free(sample_neg);
		
	(*branchNum)[i] = branchNum_pos + branchNum_neg;
	twoSampleSize = 2*sampleSize;
	l = 0;
	for(k = (*branchNum)[i]; k < twoSampleSize; k++)
		(*branches)[i][k] = sample[l++];
	
	free(sample);
}
/* *********************************************************************************** */
void sampleSimpleInterval(int i, double ***branches, int **branchNum, ddgp_interval tau, int sampleSize, double myZero){
	
	double w = tau.pos_inter[0][1] - tau.pos_inter[0][0];
	if(w > myZero){
		(*branchNum)[i] = 0;
		double h = w/((double) (sampleSize - 1));
		int k, l = 0;
		for(k = 0; k < sampleSize; k++)
			(*branches)[i][l++] = tau.neg_inter[0][0] + h*k;
		for(k = 0; k < sampleSize; k++)
			(*branches)[i][l++] = tau.pos_inter[0][0] + h*k;
	}
	else{
		int twoSampleSize = 2*sampleSize;
		if(tau.num_neg_inter > 0){
			(*branchNum)[i] = twoSampleSize - 2;
			(*branches)[i][(*branchNum)[i]] = tau.neg_inter[0][0];
			(*branches)[i][(*branchNum)[i] + 1] = tau.pos_inter[0][0];
		}
		else{
			(*branchNum)[i] = twoSampleSize - 1;
			(*branches)[i][(*branchNum)[i]] = tau.pos_inter[0][0];
		}
	}
}
/* *********************************************************************************** */
void circunference_parameters(int i, double ***A, double ***B, double ***C, double *xi3, double *xi2, double *xi1, double kappa_210, double kappa_213, double rho2_210, double q_30, double di1i2){
	
	double *vec_aux1 = alloc_vec_lf3();
	double *vec_aux2 = alloc_vec_lf3();
	double *vec_aux3 = alloc_vec_lf3();
	double *vec_aux4 = alloc_vec_lf3();
	double *vec_aux5 = alloc_vec_lf3();
	double *vec_aux6 = alloc_vec_lf3();
	
	vec_m_vec_lf3(vec_aux1, xi1, xi2);			// vec_aux1 = r = xi1 - xi2
	l_t_vec_lf3(vec_aux2, 1/di1i2, vec_aux1);		// vec_aux2 = rhat = r/||r||
	
	vec_m_vec_lf3(vec_aux3, xi3, xi2);			// vec_aux3 = v = xi3 - xi2
	
	l_t_vec_lf3(vec_aux4, kappa_210, vec_aux2); 		// vec_aux4 = kappa_210*rhat
	l_t_vec_lf3(vec_aux5, kappa_213, vec_aux2); 		// vec_aux5 = kappa_213*rhat
	
	vec_m_vec_lf3(vec_aux6, vec_aux3, vec_aux5); 		// vec_aux6 = v - kappa_213*rhat
	
	vec_p_vec_lf3(&(*A)[i][0], xi2, vec_aux4); 		// A = xi2 + kappa_210*rhat
	l_t_vec_lf3(&(*B)[i][0], rho2_210/q_30, vec_aux6);	// B = (rho2_210/q_30)(v - kappa_213*rhat)
	cross_product0_lf3(&(*C)[i][0], vec_aux2, &(*B)[i][0]);	// C = rhat x B
	
	free(vec_aux1);
	free(vec_aux2);
	free(vec_aux3);
	free(vec_aux4);
	free(vec_aux5);
	free(vec_aux6);
}
/* *********************************************************************************** */
void vertex_embedding(int i, double ***X, double *A, double *B, double *C, double tauAngle){

	double *vec_aux1 = alloc_vec_lf3();
	double *vec_aux2 = alloc_vec_lf3();
	double *vec_aux3 = alloc_vec_lf3();
	
	l_t_vec_lf3(vec_aux1, cos(tauAngle), B);	// vec_aux1 = B * cos(tau)
	l_t_vec_lf3(vec_aux2, sin(tauAngle), C);	// vec_aux2 = C * sin(tau)
	
	vec_p_vec_lf3(vec_aux3, vec_aux1, vec_aux2);	// vec_aux3 = B * cos(tau) + C * sin(tau)
	
	vec_p_vec_lf3(&(*X)[i][0], A, vec_aux3);	// X = A + B * cos(tau) + C * sin(tau)
	
	free(vec_aux1);
	free(vec_aux2);
	free(vec_aux3);
}
/* *********************************************************************************** */
void change_referential_interval_case(ddgp_interval *tau_k, double tauk_l, double tauk_u, double phase, double myZero){

	double tauk_pos_u, tauk_pos_l, tauk_neg_u, tauk_neg_l;
	int sign_tauk_pos_u, sign_tauk_pos_l, sign_tauk_neg_u, sign_tauk_neg_l;
	
	if(fabs(tauk_u - PI) < myZero){
	// tauk_u = PI
		if(tauk_l < myZero){
		// tauk_u = PI & tauk_l = 0
			(*tau_k).num_pos_inter = 1;
			(*tau_k).num_neg_inter = 1;
							
			(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
			(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);

			(*tau_k).pos_inter[0][0] = 0.0;
			(*tau_k).pos_inter[0][1] = PI;
							
			(*tau_k).neg_inter[0][0] = -PI;
			(*tau_k).neg_inter[0][1] = 0.0;
		}
		else{
		// tauk_u = PI & tauk_l != 0
			tauk_pos_l = round_lf(changeReferential(phase,  tauk_l), 6);
			tauk_neg_l = round_lf(changeReferential(phase, -tauk_l), 6);
							
			sign_tauk_pos_l = sign_double(tauk_pos_l);
			sign_tauk_neg_l = sign_double(tauk_neg_l);
			if((sign_tauk_pos_l > 0) && (sign_tauk_neg_l > 0)){
				if((tauk_u - tauk_l) > HALFPI){
					(*tau_k).num_pos_inter = 2;
					(*tau_k).num_neg_inter = 1;								
					(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
					(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
					
					(*tau_k).pos_inter[0][0] = 0.0;
					(*tau_k).pos_inter[0][1] = tauk_neg_l;
					(*tau_k).pos_inter[1][0] = tauk_pos_l;
					(*tau_k).pos_inter[1][1] = PI;
					
					(*tau_k).neg_inter[0][0] = -PI;
					(*tau_k).neg_inter[0][1] = 0.0;
				}
				else{
					(*tau_k).num_pos_inter = 1;
					(*tau_k).num_neg_inter = 0;
					(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);

					(*tau_k).pos_inter[0][0] = tauk_pos_l; //minAB_lf(tauk_neg_l, tauk_pos_l);
					(*tau_k).pos_inter[0][1] = tauk_neg_l; //maxAB_lf(tauk_neg_l, tauk_pos_l);
				}
			}
			else if((sign_tauk_pos_l < 0) && (sign_tauk_neg_l < 0)){
				if((tauk_u - tauk_l) > HALFPI){
					(*tau_k).num_pos_inter = 1;
					(*tau_k).num_neg_inter = 2;
					(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
					(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);

					(*tau_k).pos_inter[0][0] = 0.0;
					(*tau_k).pos_inter[0][1] = PI;

					(*tau_k).neg_inter[0][0] = -PI;
					(*tau_k).neg_inter[0][1] = tauk_neg_l;
					(*tau_k).neg_inter[1][0] = tauk_pos_l;
					(*tau_k).neg_inter[1][1] = 0.0;
				}
				else{
					(*tau_k).num_pos_inter = 0;
					(*tau_k).num_neg_inter = 1;
					(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
					
					(*tau_k).neg_inter[0][0] = tauk_pos_l; //minAB_lf(tauk_neg_l, tauk_pos_l);
					(*tau_k).neg_inter[0][1] = tauk_neg_l; //maxAB_lf(tauk_neg_l, tauk_pos_l);
				}
			}
			else if((sign_tauk_pos_l > 0) && (sign_tauk_neg_l < 0)){
				(*tau_k).num_pos_inter = 1;
				(*tau_k).num_neg_inter = 1;
				(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
				(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);

				(*tau_k).pos_inter[0][0] = tauk_pos_l;
				(*tau_k).pos_inter[0][1] = PI;
									
				(*tau_k).neg_inter[0][0] = -PI;
				(*tau_k).neg_inter[0][1] = tauk_neg_l;
			}
			else if((sign_tauk_pos_l < 0) && (sign_tauk_neg_l > 0)){
				(*tau_k).num_pos_inter = 1;
				(*tau_k).num_neg_inter = 1;
				(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
				(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
				
				(*tau_k).pos_inter[0][0] = 0.0;
				(*tau_k).pos_inter[0][1] = tauk_neg_l;
				
				(*tau_k).neg_inter[0][0] = tauk_pos_l;
				(*tau_k).neg_inter[0][1] = 0.0;
			}
		}
	}
	else{
	// tauk_u != PI
		if(tauk_l < myZero){
		// tauk_u != PI & tauk_l = 0
			tauk_pos_u = round_lf(changeReferential(phase,  tauk_u), 6);
			tauk_neg_u = round_lf(changeReferential(phase, -tauk_u), 6);
							
			sign_tauk_pos_u = sign_double(tauk_pos_u);
			sign_tauk_neg_u = sign_double(tauk_neg_u);
			if((sign_tauk_pos_u > 0) && (sign_tauk_neg_u > 0)){
				if((tauk_u - tauk_l) > HALFPI){
					(*tau_k).num_pos_inter = 2;
					(*tau_k).num_neg_inter = 1;
					(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
					(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);

					(*tau_k).pos_inter[0][0] = 0.0;
					(*tau_k).pos_inter[0][1] = tauk_pos_u;
					(*tau_k).pos_inter[1][0] = tauk_neg_u;
					(*tau_k).pos_inter[1][1] = PI;
									
					(*tau_k).neg_inter[0][0] = -PI;
					(*tau_k).neg_inter[0][1] = 0.0;
				}
				else{
					(*tau_k).num_pos_inter = 1;
					(*tau_k).num_neg_inter = 0;
					(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);

					(*tau_k).pos_inter[0][0] = tauk_neg_u; //minAB_lf(tauk_neg_u, tauk_pos_u);
					(*tau_k).pos_inter[0][1] = tauk_pos_u; //maxAB_lf(tauk_neg_u, tauk_pos_u);
				}
			}
			else if((sign_tauk_pos_u < 0) && (sign_tauk_neg_u < 0)){
				if((tauk_u - tauk_l) > HALFPI){
					(*tau_k).num_pos_inter = 1;
					(*tau_k).num_neg_inter = 2;
					(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
					(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);

					(*tau_k).pos_inter[0][0] = 0.0;
					(*tau_k).pos_inter[0][1] = PI;

					(*tau_k).neg_inter[0][0] = -PI;
					(*tau_k).neg_inter[0][1] = tauk_pos_u;
					(*tau_k).neg_inter[1][0] = tauk_neg_u;
					(*tau_k).neg_inter[1][1] = 0.0;
				}
				else{
					(*tau_k).num_pos_inter = 1;
					(*tau_k).num_neg_inter = 0;
					(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);

					(*tau_k).pos_inter[0][0] = tauk_neg_u; //minAB_lf(tauk_neg_u, tauk_pos_u);
					(*tau_k).pos_inter[0][1] = tauk_pos_u; //maxAB_lf(tauk_neg_u, tauk_pos_u);
				}
			}
			else if((sign_tauk_pos_u > 0) && (sign_tauk_neg_u < 0)){
				(*tau_k).num_pos_inter = 1;
				(*tau_k).num_neg_inter = 1;
				(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
				(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);

				(*tau_k).pos_inter[0][0] = 0.0;
				(*tau_k).pos_inter[0][1] = tauk_pos_u;
									
				(*tau_k).neg_inter[0][0] = tauk_neg_u;
				(*tau_k).neg_inter[0][1] = 0.0;
			}
			else if((sign_tauk_pos_u < 0) && (sign_tauk_neg_u > 0)){
				(*tau_k).num_pos_inter = 1;
				(*tau_k).num_neg_inter = 1;
				(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
				(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);

				(*tau_k).pos_inter[0][0] = tauk_neg_u;
				(*tau_k).pos_inter[0][1] = PI;
									
				(*tau_k).neg_inter[0][0] = -PI;
				(*tau_k).neg_inter[0][1] = tauk_pos_u;
			}
		}
		else{
		// tauk_u != PI & tauk_l != 0
			tauk_pos_u = round_lf(changeReferential(phase,  tauk_u), 6);
			tauk_pos_l = round_lf(changeReferential(phase,  tauk_l), 6);
			tauk_neg_u = round_lf(changeReferential(phase, -tauk_u), 6);
			tauk_neg_l = round_lf(changeReferential(phase, -tauk_l), 6);
			
			sign_tauk_pos_l = sign_double(tauk_pos_l);
			sign_tauk_pos_u = sign_double(tauk_pos_u);
			
			sign_tauk_neg_l = sign_double(tauk_neg_l);
			sign_tauk_neg_u = sign_double(tauk_neg_u);
			/* ******************************************************** *
			label:
			sign_tauk_pos_u sign_tauk_pos_l sign_tauk_neg_u sign_tauk_neg_l
			0 means > 0
			1 means < 0
			
			0 0 0 0
			0 0 0 1
			0 0 1 0
			0 0 1 1
			
			0 1 0 0
			0 1 0 1
			0 1 1 0 XXX
			0 1 1 1
			
			1 0 0 0
			1 0 0 1 XXX
			1 0 1 0
			1 0 1 1
			
			1 1 0 0
			1 1 0 1
			1 1 1 0
			1 1 1 1
			* ******************************************************** */
			if(sign_tauk_pos_u > 0){
				if(sign_tauk_pos_l > 0){
					if(sign_tauk_neg_u > 0){
						if(sign_tauk_neg_l > 0){
							// 0 0 0 0
							(*tau_k).num_pos_inter = 2;
							(*tau_k).num_neg_inter = 0;
							(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
											
							if(tauk_pos_l < tauk_neg_u){
								(*tau_k).pos_inter[0][0] = tauk_pos_l;
								(*tau_k).pos_inter[0][1] = tauk_pos_u;	
								(*tau_k).pos_inter[1][0] = tauk_neg_u;
								(*tau_k).pos_inter[1][1] = tauk_neg_l;
							}
							else{
								(*tau_k).pos_inter[0][0] = tauk_neg_u;
								(*tau_k).pos_inter[0][1] = tauk_neg_l;
								(*tau_k).pos_inter[1][0] = tauk_pos_l;
								(*tau_k).pos_inter[1][1] = tauk_pos_u;
							}
						}
						else{
							// 0 0 0 1
							(*tau_k).num_pos_inter = 2;
							(*tau_k).num_neg_inter = 1;
							(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
							(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
							
							(*tau_k).pos_inter[0][0] = tauk_pos_l;
							(*tau_k).pos_inter[0][1] = tauk_pos_u;
							(*tau_k).pos_inter[1][0] = tauk_neg_u;
							(*tau_k).pos_inter[1][1] = PI;
							
							(*tau_k).neg_inter[0][0] = -PI;
							(*tau_k).neg_inter[0][1] = tauk_neg_l;
						}
					}
					else{
						if(sign_tauk_neg_l > 0){
							// 0 0 1 0
							(*tau_k).num_pos_inter = 2;
							(*tau_k).num_neg_inter = 1;
							(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
							(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
							
							(*tau_k).pos_inter[0][0] = 0.0;
							(*tau_k).pos_inter[0][1] = tauk_neg_l;
							(*tau_k).pos_inter[1][0] = tauk_pos_l;
							(*tau_k).pos_inter[1][1] = tauk_pos_u;
							
							(*tau_k).neg_inter[0][0] = tauk_neg_u;
							(*tau_k).neg_inter[0][1] = 0.0;
						}
						else{
							// 0 0 1 1
							(*tau_k).num_pos_inter = 1;
							(*tau_k).num_neg_inter = 1;
							(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
							(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
							
							(*tau_k).pos_inter[0][0] = tauk_pos_l;
							(*tau_k).pos_inter[0][1] = tauk_pos_u;
							
							(*tau_k).neg_inter[0][0] = tauk_neg_u;
							(*tau_k).neg_inter[0][1] = tauk_neg_l;
						}
					}
				}
				else{
					if(sign_tauk_neg_u > 0){
						if(sign_tauk_neg_l > 0){
							// 0 1 0 0
							(*tau_k).num_pos_inter = 2;
							(*tau_k).num_neg_inter = 1;
							(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
							(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
							
							(*tau_k).pos_inter[0][0] = 0.0;
							(*tau_k).pos_inter[0][1] = tauk_pos_u;
							(*tau_k).pos_inter[1][0] = tauk_neg_u;
							(*tau_k).pos_inter[1][1] = tauk_neg_l;
							
							(*tau_k).neg_inter[0][0] = tauk_pos_l;
							(*tau_k).neg_inter[0][1] = 0.0;
						}
						else{
							// 0 1 0 1
							(*tau_k).num_pos_inter = 2;
							(*tau_k).num_neg_inter = 2;
							(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
							(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
										
							(*tau_k).pos_inter[0][0] = 0.0;
							(*tau_k).pos_inter[0][1] = tauk_pos_u;
							(*tau_k).pos_inter[1][0] = tauk_neg_u;
							(*tau_k).pos_inter[1][1] = PI;
							
							(*tau_k).neg_inter[0][0] = -PI;
							(*tau_k).neg_inter[0][1] = tauk_neg_l;
							(*tau_k).neg_inter[1][0] = tauk_pos_l;
							(*tau_k).neg_inter[1][1] = 0.0;
						}
					}
					else{
						if(sign_tauk_neg_l > 0){
							// 0 1 1 0
							printf("Impossible case\n");
						}
						else{
							// 0 1 1 1
							(*tau_k).num_pos_inter = 1;
							(*tau_k).num_neg_inter = 2;
							(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
							(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
							
							(*tau_k).pos_inter[0][0] = 0.0;
							(*tau_k).pos_inter[0][1] = tauk_pos_u;
							
							(*tau_k).neg_inter[0][0] = tauk_neg_u;
							(*tau_k).neg_inter[0][1] = tauk_neg_l;
							(*tau_k).neg_inter[1][0] = tauk_pos_l;
							(*tau_k).neg_inter[1][1] = 0.0;
						}
					}
				}
			}
			else{
				if(sign_tauk_pos_l > 0){
					if(sign_tauk_neg_u > 0){
						if(sign_tauk_neg_l > 0){
							// 1 0 0 0
							(*tau_k).num_pos_inter = 2;
							(*tau_k).num_neg_inter = 1;
							(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
							(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
							
							(*tau_k).pos_inter[0][0] = tauk_neg_u;
							(*tau_k).pos_inter[0][1] = tauk_neg_l;
							(*tau_k).pos_inter[1][0] = tauk_pos_l;
							(*tau_k).pos_inter[1][1] = PI;
							
							(*tau_k).neg_inter[0][0] = -PI;
							(*tau_k).neg_inter[0][1] = tauk_pos_u;
						}
						else{
							// 1 0 0 1
							printf("Impossible case\n");
						}
					}
					else{
						if(sign_tauk_neg_l > 0){
							// 1 0 1 0
							(*tau_k).num_pos_inter = 2;
							(*tau_k).num_neg_inter = 2;
							(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
							(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
							
							(*tau_k).pos_inter[0][0] = 0.0;
							(*tau_k).pos_inter[0][1] = tauk_neg_l;
							(*tau_k).pos_inter[1][0] = tauk_pos_l;
							(*tau_k).pos_inter[1][1] = PI;
							
							(*tau_k).neg_inter[0][0] = -PI;
							(*tau_k).neg_inter[0][1] = tauk_pos_u;
							(*tau_k).neg_inter[1][0] = tauk_neg_u;
							(*tau_k).neg_inter[1][1] = 0.0;
						}
						else{
							// 1 0 1 1
							(*tau_k).num_pos_inter = 1;
							(*tau_k).num_neg_inter = 2;
							(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
							(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
							
							(*tau_k).pos_inter[0][0] = tauk_pos_l;
							(*tau_k).pos_inter[0][1] = PI;
							
							(*tau_k).neg_inter[0][0] = -PI;
							(*tau_k).neg_inter[0][1] = tauk_pos_u;
							(*tau_k).neg_inter[1][0] = tauk_neg_u;
							(*tau_k).neg_inter[1][1] = tauk_neg_l;
						}
					}
				}
				else{
					if(sign_tauk_neg_u > 0){
						if(sign_tauk_neg_l > 0){
							// 1 1 0 0
							(*tau_k).num_pos_inter = 1;
							(*tau_k).num_neg_inter = 1;
							(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
							(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
							
							(*tau_k).pos_inter[0][0] = tauk_neg_u;
							(*tau_k).pos_inter[0][1] = tauk_neg_l;
							
							(*tau_k).neg_inter[0][0] = tauk_pos_l;
							(*tau_k).neg_inter[0][1] = tauk_pos_u;
						}
						else{
							// 1 1 0 1
							(*tau_k).num_pos_inter = 1;
							(*tau_k).num_neg_inter = 2;
							(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
							(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
											
							(*tau_k).pos_inter[0][0] = tauk_neg_u;
							(*tau_k).pos_inter[0][1] = PI;
											
							(*tau_k).neg_inter[0][0] = -PI;
							(*tau_k).neg_inter[0][1] = tauk_neg_l;
							(*tau_k).neg_inter[1][0] = tauk_pos_l;
							(*tau_k).neg_inter[1][1] = tauk_pos_u;
						}
					}
					else{
						if(sign_tauk_neg_l > 0){
							// 1 1 1 0
							(*tau_k).num_pos_inter = 1;
							(*tau_k).num_neg_inter = 2;
							(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
							(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
											
							(*tau_k).pos_inter[0][0] = 0.0;
							(*tau_k).pos_inter[0][1] = tauk_neg_l;
											
							(*tau_k).neg_inter[0][0] = tauk_pos_l;
							(*tau_k).neg_inter[0][1] = tauk_pos_u;
							(*tau_k).neg_inter[1][0] = tauk_neg_u;
							(*tau_k).neg_inter[1][1] = 0.0;
						}
						else{
							// 1 1 1 1
							(*tau_k).num_pos_inter = 0;
							(*tau_k).num_neg_inter = 2;
							(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
											
							if(tauk_pos_l < tauk_neg_u){
								(*tau_k).neg_inter[0][0] = tauk_pos_l;
								(*tau_k).neg_inter[0][1] = tauk_pos_u;	
								(*tau_k).neg_inter[1][0] = tauk_neg_u;
								(*tau_k).neg_inter[1][1] = tauk_neg_l;
							}
							else{
								(*tau_k).neg_inter[0][0] = tauk_neg_u;
								(*tau_k).neg_inter[0][1] = tauk_neg_l;
								(*tau_k).neg_inter[1][0] = tauk_pos_l;
								(*tau_k).neg_inter[1][1] = tauk_pos_u;
							}
						}
					}
				}
			}
		}
	}
}
/* *********************************************************************************** */
void change_referential_precise_case(ddgp_interval *tau_k, double tauk_m, double phase){
	
	double tau_k1, tau_k2;
	double tauk_pos = round_lf(changeReferential(phase,  tauk_m), 6);
	double tauk_neg = round_lf(changeReferential(phase, -tauk_m), 6);
					
	int sign_tauk_pos = sign_double(tauk_pos);
	int sign_tauk_neg = sign_double(tauk_neg);
					
	if((sign_tauk_pos > 0) && (sign_tauk_neg > 0)){
		(*tau_k).num_pos_inter = 2;
		(*tau_k).num_neg_inter = 0;
		(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
		
		tau_k1 = minAB_lf(tauk_pos, tauk_neg);
		tau_k2 = maxAB_lf(tauk_pos, tauk_neg);
		
		(*tau_k).pos_inter[0][0] = tau_k1;
		(*tau_k).pos_inter[0][1] = tau_k1;
		(*tau_k).pos_inter[1][0] = tau_k2;
		(*tau_k).pos_inter[1][1] = tau_k2;
	}
	else if((sign_tauk_pos < 0) && (sign_tauk_neg < 0)){
		(*tau_k).num_pos_inter = 0;
		(*tau_k).num_neg_inter = 2;
		(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
		
		tau_k1 = minAB_lf(tauk_pos, tauk_neg);
		tau_k2 = maxAB_lf(tauk_pos, tauk_neg);
		
		(*tau_k).neg_inter[0][0] = tau_k1;
		(*tau_k).neg_inter[0][1] = tau_k1;
		(*tau_k).neg_inter[1][0] = tau_k2;
		(*tau_k).neg_inter[1][1] = tau_k2;
	}
	else if((sign_tauk_pos > 0) && (sign_tauk_neg < 0)){
		(*tau_k).num_pos_inter = 1;
		(*tau_k).num_neg_inter = 1;
		(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
		(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
		
		(*tau_k).pos_inter[0][0] = tauk_pos;
		(*tau_k).pos_inter[0][1] = tauk_pos;
		
		(*tau_k).neg_inter[0][0] = tauk_neg;
		(*tau_k).neg_inter[0][1] = tauk_neg;
	}
	else if((sign_tauk_pos < 0) && (sign_tauk_neg > 0)){
		(*tau_k).num_pos_inter = 1;
		(*tau_k).num_neg_inter = 1;
		(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
		(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
		
		(*tau_k).pos_inter[0][0] = tauk_neg;
		(*tau_k).pos_inter[0][1] = tauk_neg;
		
		(*tau_k).neg_inter[0][0] = tauk_pos;
		(*tau_k).neg_inter[0][1] = tauk_pos;
	}
}
/* *********************************************************************************** */
void tree_backtracking(int *i, int **exploredVertex, double ***branches, int **branchNum, int sampleSize){
	
	(*exploredVertex)[(*i)--] = 0;
	while(1){
		(*branches)[(*i)][(*branchNum)[(*i)]] = 404.0;
		if((*branchNum)[(*i)] < sampleSize)
			break;
		else
			(*exploredVertex)[(*i)--] = 0;
	}
}
/* *********************************************************************************** */
void tauii3_with_parameters_d2(ddgp_interval *tau, double p30, double q30, double dii3_2_l, double dii3_2_u, double myZero){
	
	double tau_l, tau_u;
	
	// case dl = du
	if(fabs(dii3_2_l - dii3_2_u) < myZero){
		tau_l = round_lf(abs_torsion_angle_with_constants_d2(p30, q30, dii3_2_l), 6);
		
		if((tau_l < myZero) || (fabs(tau_l - PI) < myZero)){
			(*tau).num_pos_inter = 1;
			(*tau).num_neg_inter = 0;
			
			(*tau).pos_inter = alloc_mat_lf((*tau).num_pos_inter, 2);
			
			(*tau).pos_inter[0][0] = tau_l;
			(*tau).pos_inter[0][1] = tau_l;
		}
		else{
			(*tau).num_pos_inter = 1;
			(*tau).num_neg_inter = 1;
			
			(*tau).pos_inter = alloc_mat_lf((*tau).num_pos_inter, 2);
			(*tau).neg_inter = alloc_mat_lf((*tau).num_neg_inter, 2);
			
			(*tau).pos_inter[0][0] = tau_l;
			(*tau).pos_inter[0][1] = tau_l;
			
			(*tau).neg_inter[0][0] = -tau_l;
			(*tau).neg_inter[0][1] = -tau_l;
		}
	}
	// case dl < du
	else{
		tau_l = round_lf(abs_torsion_angle_with_constants_d2(p30, q30, dii3_2_l), 6);
		tau_u = round_lf(abs_torsion_angle_with_constants_d2(p30, q30, dii3_2_u), 6);
		
		(*tau).num_pos_inter = 1;
		(*tau).num_neg_inter = 1;
			
		(*tau).pos_inter = alloc_mat_lf((*tau).num_pos_inter, 2);
		(*tau).neg_inter = alloc_mat_lf((*tau).num_neg_inter, 2);
		
		(*tau).pos_inter[0][0] = tau_l;
		(*tau).pos_inter[0][1] = tau_u;
			
		(*tau).neg_inter[0][0] = -tau_u;
		(*tau).neg_inter[0][1] = -tau_l;
	}
}
/* *********************************************************************************** */
void compute_tauii3(int i, double **X, double ***discretationEdges, double myZero, double **xi1, double **xi2, double **xi3, double *di1i2, double *kappa_210, double *rho2_210, double *kappa_213, double *rho2_213, double *q30, ddgp_interval *tau){
	
	int i1, i2, i3;
	double di1i3_2, di2i3_2, dii1_2, dii2_2, dii3_2_l, dii3_2_u, kmkb, p30;
	
	i1 = (int) (discretationEdges[i][0][0] - 1.0);
	i2 = (int) (discretationEdges[i][1][0] - 1.0);
	i3 = (int) (discretationEdges[i][2][0] - 1.0);
	
	vec_e_vec_lf3((*xi1), &X[i1][0]);
	vec_e_vec_lf3((*xi2), &X[i2][0]);
	vec_e_vec_lf3((*xi3), &X[i3][0]);
	
	di1i3_2  = d2xixj_lf3((*xi3), (*xi1));
	(*di1i2) =  dxixj_lf3((*xi2), (*xi1));
	di2i3_2  = d2xixj_lf3((*xi3), (*xi2));
	
	dii1_2   = discretationEdges[i][0][1]*discretationEdges[i][0][1];
	dii2_2   = discretationEdges[i][1][1]*discretationEdges[i][1][1];
	dii3_2_l = discretationEdges[i][2][1]*discretationEdges[i][2][1];
	dii3_2_u = discretationEdges[i][2][2]*discretationEdges[i][2][2];
	
	(*kappa_210) = kappai_d2(dii2_2, (*di1i2), dii1_2);
	(*rho2_210)  = rhoi2_d2(dii2_2, (*kappa_210));
	
	(*kappa_213) = kappai_d2(di2i3_2, (*di1i2), di1i3_2);
	(*rho2_213)  = rhoi2_d2(di2i3_2, (*kappa_213));
	
	kmkb = (*kappa_210) - (*kappa_213);
	p30  = kmkb*kmkb + (*rho2_210) + (*rho2_213);
	(*q30)  = sqrt((*rho2_210)*(*rho2_213));
	
	tauii3_with_parameters_d2(&(*tau), p30, (*q30), dii3_2_l, dii3_2_u, myZero);
}
/* *********************************************************************************** */
void compute_tauiik_precise(int i, int k, double **X, prune_edges_set *pruneEdges, double kappa_210, double rho2_210, double kappa_213, double rho2_213, double di1i2, double *xi3, double *xi2, double *xi1, ddgp_interval *tau_k){
	
	int j, ik;
	double *xik = zeros_vec_lf3();
	double di1ik_2, di2ik_2, di3ik_2, diik_2_l;
	double kappa_21k, rho2_21k, kmkb, pk3, qk3, phase, pk0, qk0;
	double tauk_m;
	
	ik = (int) (pruneEdges[i].precise[k][0] - 1.0);
	
	vec_e_vec_lf3(xik, &X[ik][0]);
	
	di1ik_2 = d2xixj_lf3(xik, xi1);
	di2ik_2 = d2xixj_lf3(xik, xi2);
	di3ik_2 = d2xixj_lf3(xik, xi3);
		
	kappa_21k = kappai_d2(di2ik_2, di1i2, di1ik_2);
	rho2_21k  = rhoi2_d2(di2ik_2, kappa_21k);
	
	kmkb = kappa_213 - kappa_21k;
	pk3  = kmkb*kmkb + rho2_21k + rho2_213;
	qk3  = sqrt(rho2_21k*rho2_213);
	
	phase = sign_torsion_angle(xi3, xi2, xi1, xik)*abs_torsion_angle_with_constants_d2(pk3, qk3, di3ik_2);
	
	kmkb = kappa_210 - kappa_21k;
	pk0  = kmkb*kmkb + rho2_210 + rho2_21k;
	qk0  = sqrt(rho2_210*rho2_21k);
						
	diik_2_l = pruneEdges[i].precise[k][1]*pruneEdges[i].precise[k][1];
	
	tauk_m = abs_torsion_angle_with_constants_d2(pk0, qk0, diik_2_l);
	
	change_referential_precise_case(&(*tau_k), tauk_m, phase);
	
	for(j = 0; j < (*tau_k).num_pos_inter; j++){
		(*tau_k).pos_inter[j][0] = maxAB_lf((*tau_k).pos_inter[j][0] - 0.0088, 0.0);
		(*tau_k).pos_inter[j][1] = minAB_lf((*tau_k).pos_inter[j][1] + 0.0088, PI);
	}
	for(j = 0; j < (*tau_k).num_neg_inter; j++){
		(*tau_k).neg_inter[j][0] = maxAB_lf((*tau_k).neg_inter[j][0] - 0.0088, -PI);
		(*tau_k).neg_inter[j][1] = minAB_lf((*tau_k).neg_inter[j][1] + 0.0088, 0.0);
	}
		
	free(xik);
}
/* *********************************************************************************** */
void compute_tauiik_interval(int i, int k, double **X, prune_edges_set *pruneEdges, double kappa_210, double rho2_210, double kappa_213, double rho2_213, double di1i2, double *xi3, double *xi2, double *xi1, double myZero, ddgp_interval *tau_k, int *isTauIntervalEmpty){
	
	int ik;
	double *xik = zeros_vec_lf3();
	double di1ik_2, di2ik_2, di3ik_2, diik_2_l, diik_2_u;
	double kappa_21k, rho2_21k, kmkb, pk3, qk3, phase, pk0, qk0;
	double tauk_l, tauk_u;
		
	ik = (int) (pruneEdges[i].interval[k][0] - 1.0);
	
	vec_e_vec_lf3(xik, &X[ik][0]);
	
	di1ik_2 = d2xixj_lf3(xik, xi1);
	di2ik_2 = d2xixj_lf3(xik, xi2);
	di3ik_2 = d2xixj_lf3(xik, xi3);
		
	kappa_21k = kappai_d2(di2ik_2, di1i2, di1ik_2);
	rho2_21k  = rhoi2_d2(di2ik_2, kappa_21k);
	
	kmkb = kappa_213 - kappa_21k;
	pk3  = kmkb*kmkb + rho2_21k + rho2_213;
	qk3  = sqrt(rho2_21k*rho2_213);
	
	phase = sign_torsion_angle(xi3, xi2, xi1, xik)*abs_torsion_angle_with_constants_d2(pk3, qk3, di3ik_2);
	
	kmkb = kappa_210 - kappa_21k;
	pk0  = kmkb*kmkb + rho2_210 + rho2_21k;
	qk0  = sqrt(rho2_210*rho2_21k);
			
	diik_2_l = pruneEdges[i].interval[k][1]*pruneEdges[i].interval[k][1];
	diik_2_u = pruneEdges[i].interval[k][2]*pruneEdges[i].interval[k][2];
	
	tauk_l = abs_torsion_angle_with_constants_d2(pk0, qk0, diik_2_l);
	tauk_u = abs_torsion_angle_with_constants_d2(pk0, qk0, diik_2_u);
	
	if(fabs(tauk_u - tauk_l) < myZero)
		(*isTauIntervalEmpty) = 1;
	else				
		change_referential_interval_case(&(*tau_k), tauk_l, tauk_u, phase, myZero);
		
	free(xik);
}
/* *********************************************************************************** */
void dealloc_tau(ddgp_interval tau){
	
	if(tau.num_pos_inter > 0)
		dealloc_mat_lf(tau.pos_inter, tau.num_pos_inter);
	tau.num_pos_inter = 0;
	
	if(tau.num_neg_inter > 0)
		dealloc_mat_lf(tau.neg_inter, tau.num_neg_inter);
	tau.num_neg_inter = 0;
}
/* *********************************************************************************** */
void iABP(double ***X, int n, int *lev, double *cpu_time_used, double *numIt, double ***discretationEdges, prune_edges_set *pruneEdges, int sampleSize, double myZero, double timeMAX, double tolerance){
	
	// ---------- Variables initialization ----------
	
	(*lev) = 2;
	(*numIt) = 0.000003;
	
	clock_t start, end;
	
	int k, isTauIntervalEmpty;
	int i = 3;
	int twoSampleSize = 2*sampleSize;
	int *exploredVertex = zeros_vec_d(n);
	int *branchNum = zeros_vec_d(n);
	
	double kappa_210, rho2_210, kappa_213, rho2_213, q30, di1i2;
    	double *xi1 = zeros_vec_lf3();
	double *xi2 = zeros_vec_lf3();
	double *xi3 = zeros_vec_lf3();
	double *xik = zeros_vec_lf3();
	double **A  = zeros_mat_lf(n, 3);
	double **B  = zeros_mat_lf(n, 3);
	double **C  = zeros_mat_lf(n, 3);
	double **Xr = zeros_mat_lf(n, 3);
	double **branches = ones_mat_lf(n, twoSampleSize);
	l_t_mat_lf(404.0, branches, n, twoSampleSize);
    	
	ddgp_interval tau, tau_k;
    	
	// ----------------------------------------------
    	
    	mat_e_mat_lf((*X), Xr, n, 3);
    	
    	start = clock();
    	clock_t startTime = time(NULL);
    	
    	referential_x1_x2_x3(&Xr, discretationEdges);
    	
    	while((i < n) && (difftime(time(NULL), startTime) < timeMAX)){
	    	//printf("***********************************************\n");
    		//printf("i = %d\n", i);
    		//printf("***********************************************\n");
    		if(exploredVertex[i] == 0){
    			// ************************************************************
    			// Compute tau_{i,i3}
    			// ************************************************************
    			compute_tauii3(i, Xr, discretationEdges, myZero, &xi1, &xi2, &xi3, &di1i2, &kappa_210, &rho2_210, &kappa_213, &rho2_213, &q30, &tau);
			
			isTauIntervalEmpty = 0;
			
			//printf("\There are %d precise prune edges\n\n", pruneEdges[i].cardUkP);
			for(k = 0; k < pruneEdges[i].cardUkP; k++){
				// ************************************************************
		    		// Compute tau_{i,ik}
		    		// ************************************************************
		    		compute_tauiik_precise(i, k, Xr, pruneEdges, kappa_210, rho2_210, kappa_213, rho2_213, di1i2, xi3, xi2, xi1, &tau_k);
				
				interIntervals(&tau, tau_k, myZero);
				
				dealloc_tau(tau_k);
				
				if((tau.num_pos_inter + tau.num_neg_inter) == 0){
					isTauIntervalEmpty = 1;
					break;
				}
			}
			
			if(isTauIntervalEmpty == 0){
				//printf("\There are %d interval prune edges\n\n", pruneEdges[i].cardUkI);
				for(k = 0; k < pruneEdges[i].cardUkI; k++){
					// ************************************************************
		    			// Compute tau_{i,ik}
		    			// ************************************************************
		    			compute_tauiik_interval(i, k, Xr, pruneEdges, kappa_210, rho2_210, kappa_213, rho2_213, di1i2, xi3, xi2, xi1, myZero, &tau_k, &isTauIntervalEmpty);
		    			if(isTauIntervalEmpty == 1){
						dealloc_tau(tau);
						break;
					}
					
					interIntervals(&tau, tau_k, myZero);
					
					dealloc_tau(tau_k);
					
					if((tau.num_pos_inter + tau.num_neg_inter) == 0){
						isTauIntervalEmpty = 1;
						break;
					}
				}
			}
			
			if(isTauIntervalEmpty == 0){
				sampleInterval(i, &branches, &branchNum, tau, sampleSize, tolerance);
				
				dealloc_tau(tau);
				
				exploredVertex[i] = 1;
				circunference_parameters(i, &A, &B, &C, xi3, xi2, xi1, kappa_210, kappa_213, rho2_210, q30, di1i2);
			}
			else{
				// backtracking in the tree
				tree_backtracking(&i, &exploredVertex, &branches, &branchNum, twoSampleSize);
				if(i == 2){
					printf("iABP: ERROR! The instance solution was not found\n");
					printf("iABP: The whole search space was explored\n");
					break;
				}
				continue;
			}
    		}
    		else
			branches[i][branchNum[i]++] = 404.0;
		if(branchNum[i] < twoSampleSize){
			vertex_embedding(i, &Xr, &A[i][0], &B[i][0], &C[i][0], branches[i][branchNum[i]]);
			
			(*numIt)+= 0.000001;
			i++;
			if((*lev) < i)
				(*lev) = i;
		}
		else{
			tree_backtracking(&i, &exploredVertex, &branches, &branchNum, twoSampleSize);
			if(i == 2){
				printf("iABP: ERROR! The instance solution was not found\n");
				printf("iABP: The whole search space was explored\n");
				break;
			}
			continue;
		}
    	}
    	
    	if(i == n)
    		mat_e_mat_lf((*X), Xr, n, 3);
    	
	free(xi1);
	free(xi2);
    	free(xi3);
    	free(xik);
    	
    	dealloc_mat_lf(A, n);
    	dealloc_mat_lf(B, n);
    	dealloc_mat_lf(C, n);
    	dealloc_mat_lf(Xr, n);
    	
    	dealloc_mat_lf(branches, n);
    	free(exploredVertex);
	free(branchNum);
    	    	
	end = clock();
	(*cpu_time_used) = ((double) (end - start))/CLOCKS_PER_SEC;
	
	if((*cpu_time_used) > timeMAX){
		printf("iABP: Time Limit Reached.\n");
		printf("iABP: ERROR! The instance solution was not found\n");
	}
}
/* *********************************************************************************** */
void iBP(double ***X, int n, int *lev, double *cpu_time_used, double *numIt, double ***discretationEdges, prune_edges_set *pruneEdges, int sampleSize, double myZero, double timeMAX, double tolerance){
	
	// ---------- Variables initialization ----------
	
	(*lev) = 2;
	(*numIt) = 0.000003;
	
	clock_t start, end;
	
	int k, ik;
	int i = 3;
	int isTauIntervalEmpty;
	int twoSampleSize = 2*sampleSize;
	int *exploredVertex = zeros_vec_d(n);
	int *branchNum = zeros_vec_d(n);
	
	double kappa_210, rho2_210, kappa_213, rho2_213, q30;
	double di1i2, diik_2_l, diik_2_0, diik_2_u;
	double **A  = zeros_mat_lf(n,3);
	double **B  = zeros_mat_lf(n,3);
	double **C  = zeros_mat_lf(n,3);
	double **Xr = zeros_mat_lf(n,3);
	double *xi1 = zeros_vec_lf3();
	double *xi2 = zeros_vec_lf3();
	double *xi3 = zeros_vec_lf3();
	double **branches = ones_mat_lf(n, twoSampleSize);
	l_t_mat_lf(404.0, branches, n, twoSampleSize);
	
	ddgp_interval tau;
	
    	// ----------------------------------------------
    	
	mat_e_mat_lf((*X), Xr, n, 3);
	
    	start = clock();
    	clock_t startTime = time(NULL);
    	
    	referential_x1_x2_x3(&Xr, discretationEdges);
    	
	while((i < n) && (difftime(time(NULL), startTime) < timeMAX)){
		//printf("***********************************************\n");
		//printf("i = %d\n", i);
		//printf("***********************************************\n");
		if(exploredVertex[i] == 0){
    			// ************************************************************
    			// Compute tau_{i,i3}
    			// ************************************************************
    			compute_tauii3(i, Xr, discretationEdges, myZero, &xi1, &xi2, &xi3, &di1i2, &kappa_210, &rho2_210, &kappa_213, &rho2_213, &q30, &tau);
						
			sampleSimpleInterval(i, &branches, &branchNum, tau, sampleSize, myZero);	
			
			dealloc_tau(tau);
			
			exploredVertex[i] = 1;
			circunference_parameters(i, &A, &B, &C, xi3, xi2, xi1, kappa_210, kappa_213, rho2_210, q30, di1i2);
		}
    		else
			branches[i][branchNum[i]++] = 404.0;
		
		if(branchNum[i] < twoSampleSize){
			vertex_embedding(i, &Xr, &A[i][0], &B[i][0], &C[i][0], branches[i][branchNum[i]]);
			(*numIt)+= 0.000001;
			
		}
		else{
			tree_backtracking(&i, &exploredVertex, &branches, &branchNum, twoSampleSize);
			if(i == 2){
				printf("iBP: The whole search space was explored\n");
				printf("iBP: ERROR! The instance solution was not found\n");
				break;
			}
			continue;
		}
				
		isTauIntervalEmpty = 0;
		
		if(pruneEdges[i].cardUkP > 0){
			//printf("\There are %d precise prune edges\n\n", pruneEdges[i].cardUkP);
			for(k = 0; k < pruneEdges[i].cardUkP; k++){
				ik = (int) (pruneEdges[i].precise[k][0] - 1.0);
				
		    		diik_2_0 = d2xixj_lf3(&Xr[i][0], &Xr[ik][0]);
				
				diik_2_l = pruneEdges[i].precise[k][1]*pruneEdges[i].precise[k][1];
			    	
				if(fabs(diik_2_l - diik_2_0) > tolerance){
					isTauIntervalEmpty = 1;
					break;
				}
			}
		}
		
		if((isTauIntervalEmpty == 0) && (pruneEdges[i].cardUkI > 0)){
			//printf("\There are %d interval prune edges\n\n", pruneEdges[i].cardUkI);
			for(k = 0; k < pruneEdges[i].cardUkI; k++){
				ik = (int) (pruneEdges[i].interval[k][0] - 1.0);
					
				diik_2_0 = d2xixj_lf3(&Xr[i][0], &Xr[ik][0]);
				
				diik_2_l = pruneEdges[i].interval[k][1]*pruneEdges[i].interval[k][1];
		    		diik_2_u = pruneEdges[i].interval[k][2]*pruneEdges[i].interval[k][2];
			    						
				if((diik_2_0 < fabs(diik_2_l - tolerance)) || (fabs(diik_2_u + tolerance) < diik_2_0)){
					isTauIntervalEmpty = 1;
					break;
				}
			}
		}
		
		if(isTauIntervalEmpty == 0){
			i++;
			if((*lev) < i)
				(*lev) = i;
		}
    	}
    	
    	if(i == n)
    		mat_e_mat_lf((*X), Xr, n, 3);

	free(xi1);
	free(xi2);
    	free(xi3);
    	
    	dealloc_mat_lf(A, n);
    	dealloc_mat_lf(B, n);
    	dealloc_mat_lf(C, n);
    	dealloc_mat_lf(Xr, n);
    	
    	dealloc_mat_lf(branches, n);
    	free(exploredVertex);
	free(branchNum);
    	    	
	end = clock();
	(*cpu_time_used) = ((double) (end - start))/CLOCKS_PER_SEC;
	
	if((*cpu_time_used) > timeMAX){
		printf("iBP: Time Limit Reached.\n");
		printf("iBP: ERROR! The instance solution was not found\n");
	}
}
/* *********************************************************************************** */
void satisfiesInstanceQM(double **X, dcinputfile *dc_vec, int m, double resolution){

	int k, i, j, sum, ndd;
	int *flags = ones_vec_d(m);

	double dL_2, dU_2, dij_2, dL, dU, dij;
	double *xi = zeros_vec_lf3();
	double *xj = zeros_vec_lf3();
	
	ndd = 6;
	for(k = 0; k < m; k++){
		i = dc_vec[k].i - 1;
		j = dc_vec[k].j - 1;

		vec_e_vec_lf3(xi, &X[i][0]);
		vec_e_vec_lf3(xj, &X[j][0]);

		dij_2 = round_lf(d2xixj_lf3(xi, xj), ndd);

		dL_2 = round_lf(dc_vec[k].dl*dc_vec[k].dl, ndd);
		dU_2 = round_lf(dc_vec[k].du*dc_vec[k].du, ndd);
		
		if(((dL_2 < dij_2) || (fabs(dL_2 - dij_2) < resolution)) &&
		   ((dij_2 < dU_2) || (fabs(dij_2 - dU_2) < resolution)))
			flags[k] = 0;
	}
	
	sum = sum_vec_d(flags, m);
	
	if(sum == 0)
		printf("X Satisfies the Instance\n");
	else{
	//if(sum > 0){
		printf("X does not Satisfies the Instance: %d wrong edges were detected\n", sum);
		for(k = 0; k < m; k++)
			if(flags[k] == 1){
				i = dc_vec[k].i - 1;
				j = dc_vec[k].j - 1;

				vec_e_vec_lf3(xi, &X[i][0]); //xi = X(i, :);
				vec_e_vec_lf3(xj, &X[j][0]); //xj = X(j, :);

				dij = dxixj_lf3(xi, xj);

				dL = dc_vec[k].dl;
				dU = dc_vec[k].du;
				printf("line %.3d: d_%.2d,%.2d_L = %.6lf <= %.6lf <= %.6lf = d_%.2d,%.2d_U\n", k+1, i+1, j+1, dL, dij, dU, i+1, j+1);
			}
		
	}
	
	free(flags);
	free(xi);
	free(xj);
}
