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
#define MYZERO 0.000001
#define CF_RAD2DEG 57.29577951308
#define CF_DEG2RAD 0.01745329252

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

	int k, i;
	int cardUk, *adjpred, *pos_adjpred_vk;
	int *vertices = alloc_vec_d(m);
	
	for(k = 0; k < m; k++)
		vertices[k] = dc_vec[k].i;	
	
	(*n) = max_vec_d(vertices, m);
	
	(*kv) = adjacent_predecessors_cardinality(vertices, m, (*n));
	
	(*instance_adjlist) = alloc_3Darray_lf((*n), (*kv));
	
	(*protein) = alloc_vec_proteinstructure((*n));
	
	(*protein)[0].v_k = dc_vec[0].j;
	(*protein)[0].res_seq = dc_vec[0].rj;
	strcpy((*protein)[0].atom_name, dc_vec[0].aj);
	strcpy((*protein)[0].res_name, dc_vec[0].rtj);
	(*protein)[0].x = 0.0;
	(*protein)[0].y = 0.0;
	(*protein)[0].z = 0.0;
	for(k = 1; k < (*n); k++){
		adjpred = find_val_vec_d(vertices, m, k+1);
		cardUk = sum_vec_d(adjpred, m);
		pos_adjpred_vk = find_ones_position_vec_d(adjpred, m, cardUk);

		(*protein)[k].v_k = dc_vec[pos_adjpred_vk[0]].i;
		(*protein)[k].res_seq = dc_vec[pos_adjpred_vk[0]].ri;
		strcpy((*protein)[k].atom_name, dc_vec[pos_adjpred_vk[0]].ai);
		strcpy((*protein)[k].res_name, dc_vec[pos_adjpred_vk[0]].rti);
		
		(*protein)[k].x = 0.0;
		(*protein)[k].y = 0.0;
		(*protein)[k].z = 0.0;
		
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
void referential_x1_x2_x3(double ***X, double ***discretizationEdges_2){
	
	double d12_2 = discretizationEdges_2[1][0][1];//[2][1][2]-1
	double d13_2 = discretizationEdges_2[2][1][1];//[3][2][2]-1
	double d23_2 = discretizationEdges_2[2][0][1];//[3][1][2]-1
	double d23 = sqrt(d23_2);
	
	double d12cosTh = 0.5*(d12_2 + d23_2 - d13_2)/d23;
	double d12sinTh = sqrt(d12_2 - d12cosTh*d12cosTh);
	
	(*X)[0][0] = -d12sinTh;
	(*X)[0][1] = d12cosTh - d23;
	(*X)[1][1] = -d23;
}
/* *********************************************************************************** */
void adjlist_2_discretization_edges_2(double ***instance_adjlist, int *kv, int n, int **cliques, double ****discretizationEdges_2){

	int k, i;
	double dl, du;
	(*discretizationEdges_2) = alloc_fix3Darray_lf(n, 3, 3);
	
	(*discretizationEdges_2)[1][0][0] = cliques[1][1];
	dl = dij_from_adjlist(instance_adjlist, kv, 2, 1, "l");
	du = dij_from_adjlist(instance_adjlist, kv, 2, 1, "u");
	(*discretizationEdges_2)[1][0][1] = dl*dl;
	(*discretizationEdges_2)[1][0][2] = du*du;
	
	(*discretizationEdges_2)[2][0][0] = cliques[2][1];
	dl = dij_from_adjlist(instance_adjlist, kv, 3, 2, "l");
	du = dij_from_adjlist(instance_adjlist, kv, 3, 2, "u");
	(*discretizationEdges_2)[2][0][1] = dl*dl;
	(*discretizationEdges_2)[2][0][2] = du*du;
	
	(*discretizationEdges_2)[2][1][0] = cliques[2][2];
	dl = dij_from_adjlist(instance_adjlist, kv, 3, 1, "l");
	du = dij_from_adjlist(instance_adjlist, kv, 3, 1, "u");
	(*discretizationEdges_2)[2][1][1] = dl*dl;
	(*discretizationEdges_2)[2][1][2] = du*du;
	
	for(k = 3; k < n; k++)
		for(i = 0; i < 3; i++){
			(*discretizationEdges_2)[k][i][0] = cliques[k][i+1];
			dl = dij_from_adjlist(instance_adjlist, kv, k+1, cliques[k][i+1], "l");
			du = dij_from_adjlist(instance_adjlist, kv, k+1, cliques[k][i+1], "u");
			(*discretizationEdges_2)[k][i][1] = dl*dl;
			(*discretizationEdges_2)[k][i][2] = du*du;
		}
}
/* *********************************************************************************** */
void adjlist_2_prune_edges_2(double ***instance_adjlist, int *kv, int n, int **cliques, prune_edges_set **pruneEdges_2){

	int i, j;
	double *Ui, dl, du;
	int *vec1i1, *vec1i2, *vec1i3, *vec1i1i2i3, numPd, numId, numPd_i, numId_i;
	
	(*pruneEdges_2) = alloc_vec_pruneedgesset(n);
	
	(*pruneEdges_2)[0].cardUkP = 0;
	(*pruneEdges_2)[0].cardUkI = 0;
	(*pruneEdges_2)[1].cardUkP = 0;
	(*pruneEdges_2)[1].cardUkI = 0;
	(*pruneEdges_2)[2].cardUkP = 0;
	(*pruneEdges_2)[2].cardUkI = 0;
	
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
					if(fabs(instance_adjlist[i][j][1] - instance_adjlist[i][j][2]) < MYZERO)
						numPd++;
					else
						numId++;
				}
			
			(*pruneEdges_2)[i].cardUkP = numPd;
			(*pruneEdges_2)[i].cardUkI = numId;
			if(numPd > 0)
				(*pruneEdges_2)[i].precise = alloc_mat_lf(numPd, 3);
			if(numId > 0)
				(*pruneEdges_2)[i].interval = alloc_mat_lf(numId, 3);
			
			numPd_i = 0;
			numId_i = 0;
			
			for(j = 0; j < kv[i]; j++)
				if(!vec1i1i2i3[j]){
					dl = dij_from_adjlist(instance_adjlist, kv, i+1, Ui[j], "l");
					du = dij_from_adjlist(instance_adjlist, kv, i+1, Ui[j], "u");
					if(fabs(dl - du) < MYZERO){
						(*pruneEdges_2)[i].precise[numPd_i][0] = instance_adjlist[i][j][0];
						(*pruneEdges_2)[i].precise[numPd_i][1] = instance_adjlist[i][j][1]*instance_adjlist[i][j][1];
						(*pruneEdges_2)[i].precise[numPd_i][2] = instance_adjlist[i][j][2]*instance_adjlist[i][j][2];
						numPd_i++;
					}
					else{
						(*pruneEdges_2)[i].interval[numId_i][0] = instance_adjlist[i][j][0];
						(*pruneEdges_2)[i].interval[numId_i][1] = instance_adjlist[i][j][1]*instance_adjlist[i][j][1];
						(*pruneEdges_2)[i].interval[numId_i][2] = instance_adjlist[i][j][2]*instance_adjlist[i][j][2];
						numId_i++;
					}
				}
			free(Ui);
			free(vec1i1i2i3);
			
			double **auxMat = alloc_mat_lf(numId_i, 3);
			mat_e_mat_lf(auxMat, (*pruneEdges_2)[i].interval, numId_i, 3);
			
			double *vecW = alloc_vec_lf(numId_i);
			for(j = 0; j < numId_i; j++)
				vecW[j] = (*pruneEdges_2)[i].interval[j][2] - (*pruneEdges_2)[i].interval[j][1];
				
			int *vecIndices = alloc_vec_d(numId_i);
			for(j = 0; j < numId_i; j++)
				vecIndices[j] = j;
				
			quicksort(vecW, vecIndices, 0, numId_i - 1);
			
				for(j = 0; j < numId_i; j++){
					(*pruneEdges_2)[i].interval[j][0] = auxMat[vecIndices[j]][0];
					(*pruneEdges_2)[i].interval[j][1] = auxMat[vecIndices[j]][1];
					(*pruneEdges_2)[i].interval[j][2] = auxMat[vecIndices[j]][2];
				}
			
			dealloc_mat_lf(auxMat, numId_i);
			free(vecW);
		}
		else{
			(*pruneEdges_2)[i].cardUkP = 0;
			(*pruneEdges_2)[i].cardUkI = 0;
		}
}
/* *********************************************************************************** */
void adjlist_2_graph_parameters(double ***instance_adjlist, int *kv, int n, int **cliques, double ****discretizationEdges_2, prune_edges_set **pruneEdges_2){

	adjlist_2_discretization_edges_2(instance_adjlist, kv, n, cliques, &(*discretizationEdges_2));
	adjlist_2_prune_edges_2(instance_adjlist, kv, n, cliques, &(*pruneEdges_2));
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
double abs_torsion_angle_with_constants_d2(double p, double q, double dii3_2){
	
	double absTau = 0.0;
	double twoq = 2.0*q;
	double dmin2 = p - twoq;
	double dmax2 = p + twoq;
	
	if((dmin2 < dii3_2) && (dii3_2 < dmax2))
		absTau = acos(round_lf((p - dii3_2)/twoq, 6));
	else if((dmax2 < dii3_2) || (fabs(dmax2 - dii3_2) < MYZERO))
		absTau = PI;

	return absTau;
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
double change_referential_1_to_0(double phase, double tau_ref1){
	
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
void interA_cap_interB(double *interA, double *interB, double **interC, int *isCapEmpty){
	
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
	
	int aA_e_bA = (fabs(bA - aA) < MYZERO) ? 1 : 0;
	int aB_e_bB = (fabs(bB - aB) < MYZERO) ? 1 : 0;
	int bA_e_aB = (fabs(bA - aB) < MYZERO) ? 1 : 0;
	
	(*isCapEmpty) = 1;
	
	if((aB < bA) || bA_e_aB){
		(*isCapEmpty) = 0;
		if(!aA_e_bA){
			(*interC)[0] = aB;
			if(aB_e_bB || bA_e_aB)
				(*interC)[1] = aB;
			else
				(*interC)[1] = minAB_lf(bA, bB);
		}
		else{
			(*interC)[0] = aA;
			(*interC)[1] = aA;
		}
	}
}
/* *********************************************************************************** */
void interA_cap_interB_matrix(double *interA, double *interB, double ***interC_matrix, int *matrix_row){
	
	int isEmpty;
	double *interC = alloc_vec_lf(2);
	
	interA_cap_interB(interA, interB, &interC, &isEmpty);
	
	if(isEmpty == 0){
		(*interC_matrix)[(*matrix_row)][0] = interC[0];
		(*interC_matrix)[(*matrix_row)][1] = interC[1];
		(*matrix_row)++;
	}
	
	free(interC);
}
/* *********************************************************************************** */
void cap_interAk(ddgp_interval *tau, ddgp_interval tau_k){
	
	int k, l;
	double **Ak_cap_Al_matrix;
	
	int num_inter, num_max_inter;
	
	// positive interval intersection
	if((*tau).num_pos_inter > 0){
		num_max_inter = (*tau).num_pos_inter + tau_k.num_pos_inter;
		Ak_cap_Al_matrix = alloc_mat_lf(num_max_inter, 2);
		num_inter = 0;
		for(k = 0; k < (*tau).num_pos_inter; k++)
			for(l = 0; l < tau_k.num_pos_inter; l++)
				if(fabs((tau_k.pos_inter[l][1] - tau_k.pos_inter[l][0]) - PI) < MYZERO){
					Ak_cap_Al_matrix[num_inter][0] = (*tau).pos_inter[k][0];
					Ak_cap_Al_matrix[num_inter][1] = (*tau).pos_inter[k][1];
					num_inter++;
				}
				else
					interA_cap_interB_matrix(&(*tau).pos_inter[k][0], &tau_k.pos_inter[l][0], &Ak_cap_Al_matrix, &num_inter);
		
		dealloc_mat_lf((*tau).pos_inter, (*tau).num_pos_inter);
		
		if(num_inter > 0){
			(*tau).num_pos_inter = num_inter;
			(*tau).pos_inter = alloc_mat_lf(num_inter, 2);
			l = 0;
			for(k = 0; k < num_inter; k++){
				(*tau).pos_inter[l][0] = Ak_cap_Al_matrix[l][0];
				(*tau).pos_inter[l][1] = Ak_cap_Al_matrix[l][1];
				l++;
			}
		}
		else
			(*tau).num_pos_inter = 0;
			
		dealloc_mat_lf(Ak_cap_Al_matrix, num_max_inter);
	}
	
	// negative interval intersection
	if((*tau).num_neg_inter > 0){
		num_max_inter = (*tau).num_neg_inter + tau_k.num_neg_inter;
		Ak_cap_Al_matrix = alloc_mat_lf(num_max_inter, 2);
		num_inter = 0;
		for(k = 0; k < (*tau).num_neg_inter; k++)
			for(l = 0; l < tau_k.num_neg_inter; l++)
				if(fabs((tau_k.neg_inter[l][1] - tau_k.neg_inter[l][0]) - PI) < MYZERO){
					Ak_cap_Al_matrix[num_inter][0] = (*tau).neg_inter[k][0];
					Ak_cap_Al_matrix[num_inter][1] = (*tau).neg_inter[k][1];
					num_inter++;
				}
				else
					interA_cap_interB_matrix(&(*tau).neg_inter[k][0], &tau_k.neg_inter[l][0], &Ak_cap_Al_matrix, &num_inter);
				
		dealloc_mat_lf((*tau).neg_inter, (*tau).num_neg_inter);
		
		if(num_inter > 0){
			(*tau).num_neg_inter = num_inter;
			(*tau).neg_inter = alloc_mat_lf(num_inter, 2);
			l = 0;
			for(k = 0; k < num_inter; k++){
				(*tau).neg_inter[l][0] = Ak_cap_Al_matrix[l][0];
				(*tau).neg_inter[l][1] = Ak_cap_Al_matrix[l][1];
				l++;
			}
		}
		else
			(*tau).num_neg_inter = 0;
			
		dealloc_mat_lf(Ak_cap_Al_matrix, num_max_inter);
	}
}
/* *********************************************************************************** */
void sample_interval_union(double **J, int numInter, double **sample, int *sampleSize, int maxSampleSize, double resolution){
	
	int ij, k, l, fullSampleSize;
	int *vec_nk, numPoints_k;
	double *vec_wk, *fullSample, *sample_k, *sampleAux;
	double maxW;
	
	(*sampleSize) = 0;
	if(numInter > 0){
		vec_wk = alloc_vec_lf(numInter);
		maxW = 0.0;
		for(k = 0; k < numInter; k++){
			vec_wk[k] = J[k][1] - J[k][0];
			if(maxW < vec_wk[k])
				maxW = vec_wk[k];
		}
		vec_nk = alloc_vec_d(numInter);
		for(k = 0; k < numInter; k++)
			if(vec_wk[k] > resolution)
				vec_nk[k] = minAB_d(maxAB_d((int) round(((double) maxSampleSize)*vec_wk[k]/maxW), 1), (int) floor(vec_wk[k]/resolution));
			else
				vec_nk[k] = 1;
		
		fullSampleSize = sum_vec_d(vec_nk, numInter);
		fullSample = alloc_vec_lf(fullSampleSize);
		
		l = 0;
		for(k = 0; k < numInter; k++){
			numPoints_k = 0;
			sample_simple_interval(&sample_k, &J[k][0], vec_nk[k], resolution, &numPoints_k);
			for(ij = 0; ij < numPoints_k; ij++)
				fullSample[l++] = sample_k[ij];
			
			free(sample_k);
		}
		
		if(fullSampleSize <= maxSampleSize){
			(*sampleSize) = fullSampleSize;
			(*sample) = alloc_vec_lf((*sampleSize));
			sort_sample_indices(&(*sample), fullSample, (*sampleSize));
		}
		else{
			(*sampleSize) = maxSampleSize;
			sampleAux = alloc_vec_lf((*sampleSize));
			int pos;
			double a, b;
			double step = ((double) (fullSampleSize - 1))/((double) maxSampleSize);
			a = 0.0;
			for(k = 0; k < (*sampleSize); k++){
				b = a + step;
				pos = (int) round((a + b)/2.0);
				a = b;
				sampleAux[k] = fullSample[pos];
			}
			
			(*sample) = alloc_vec_lf((*sampleSize));
			sort_sample_indices(&(*sample), sampleAux, (*sampleSize));
			free(sampleAux);
		}
		
		free(vec_wk);
		free(vec_nk);
		free(fullSample);
		
	}
}
/* *********************************************************************************** */
void sample_interval_union_DDGP(int i, double ***branches, int **branchNum, ddgp_interval tau, int sampleSize, double resolution){
	
	int k, l, fullSampleSize;
	double *sample, *samplePos, *sampleNeg;
	int sampleSizePos = 0;
	int sampleSizeNeg = 0;
	int twoSampleSize = 2*sampleSize;
	
	// positive interval
	sample_interval_union(tau.pos_inter, tau.num_pos_inter, &samplePos, &sampleSizePos, sampleSize, resolution);
	// negative interval
	sample_interval_union(tau.neg_inter, tau.num_neg_inter, &sampleNeg, &sampleSizeNeg, sampleSize, resolution);
	
	fullSampleSize = sampleSizeNeg + sampleSizePos;
	
	sample = alloc_vec_lf(fullSampleSize);
	
	l = 0;
	for(k = 0; k < sampleSizeNeg; k++)
		sample[l++] = sampleNeg[k];
	for(k = 0; k < sampleSizePos; k++)
		sample[l++] = samplePos[k];
	
	if(sampleSizeNeg > 0)
		free(sampleNeg);
	if(sampleSizePos > 0)
		free(samplePos);
	
	(*branchNum)[i] = twoSampleSize - fullSampleSize;
	l = 0;
	for(k = (*branchNum)[i]; k < twoSampleSize; k++)
		(*branches)[i][k] = sample[l++];
	
	free(sample);
}
/* *********************************************************************************** */
void sample_simple_DDGPinterval(int i, double ***branches, int **branchNum, ddgp_interval tau, int sampleSize, double resolution){
	
	int k, l;
	int numNegPoints = 0;
	int numPosPoints = 0;
	double *sampleNeg, *samplePos;
	
	if(tau.num_pos_inter > 0)
		sample_simple_interval(&samplePos, &tau.pos_inter[0][0], sampleSize, resolution, &numPosPoints);
	
	if(tau.num_neg_inter > 0)
		sample_simple_interval(&sampleNeg, &tau.neg_inter[0][0], sampleSize, resolution, &numNegPoints);	
	
	(*branchNum)[i] = 2*sampleSize - (numPosPoints + numNegPoints);
	
	l = (*branchNum)[i];
	for(k = 0; k < numNegPoints; k++)
		(*branches)[i][l++] = sampleNeg[k];
	for(k = 0; k < numPosPoints; k++)
		(*branches)[i][l++] = samplePos[k];
		
	if(numNegPoints > 0)
		free(sampleNeg);
	if(numPosPoints > 0)
		free(samplePos);
}
/* *********************************************************************************** */
void sample_simple_interval(double **sample, double *interval, int sampleSize, double tolerance, int *numPoints){
	
	double w = interval[1] - interval[0];
	
	if((w > tolerance) && (sampleSize > 1)){
		int k;
		double h;
		double *sampleAux;
		
		(*numPoints) = minAB_d(sampleSize, (int) floor(w/tolerance));
		sampleAux = zeros_vec_lf((*numPoints));
		
		h = w/((double) ((*numPoints) - 1));
		sampleAux[0] = interval[0];
		for(k = 1; k < (*numPoints); k++)
			sampleAux[k] = sampleAux[k-1] + h;
		
		(*sample) = alloc_vec_lf((*numPoints));
		sort_sample_indices(&(*sample), sampleAux, (*numPoints));
		
		free(sampleAux);
	}
	else{
		(*numPoints) = 1;
		(*sample) = alloc_vec_lf((*numPoints));
		(*sample)[0] = (interval[0] + interval[1])/2.0;
	}	
}
/* *********************************************************************************** */
void sort_sample_indices(double **sortedSample, double *sample, int sampleSize){
	
	int k, centralIndex, lastIndexLeft, lastIndexRight;
	
	if(sampleSize == 1)
		(*sortedSample)[0] = sample[0];
	else if(sampleSize % 2 == 0){
		lastIndexRight = sampleSize/2;
		lastIndexLeft  = lastIndexRight-1;
		(*sortedSample)[0] = sample[lastIndexLeft];
		(*sortedSample)[1] = sample[lastIndexRight];
		
		for(k = 2; k < sampleSize; k+=2){
			(*sortedSample)[k]   = sample[--lastIndexLeft];
			(*sortedSample)[k+1] = sample[++lastIndexRight];		
		}
	}
	else{
		centralIndex = (int) floor(((double) sampleSize)/2.0);
		(*sortedSample)[0] = sample[centralIndex];
		lastIndexLeft  = centralIndex - 1;
		lastIndexRight = centralIndex + 1;
		(*sortedSample)[1] = sample[lastIndexLeft];
		(*sortedSample)[2] = sample[lastIndexRight];
		
		for(k = 3; k < sampleSize; k+=2){
			(*sortedSample)[k]   = sample[--lastIndexLeft];
			(*sortedSample)[k+1] = sample[++lastIndexRight];
		
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
void change_referential_interval_case(ddgp_interval *tau_k, double absTauk_l, double absTauk_u, double phase){

	double tauk_pos_u, tauk_pos_l, tauk_neg_u, tauk_neg_l;
	int sign_tauk_pos_u, sign_tauk_pos_l, sign_tauk_neg_u, sign_tauk_neg_l;
	
	if(fabs(absTauk_u - PI) < MYZERO){
	// absTauk_u = PI
		if(absTauk_l < MYZERO){
		// absTauk_u = PI & absTauk_l = 0
			(*tau_k).num_neg_inter = 1;
			(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
			(*tau_k).neg_inter[0][0] = -PI;
			(*tau_k).neg_inter[0][1] = 0.0;
			
			(*tau_k).num_pos_inter = 1;				
			(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
			(*tau_k).pos_inter[0][0] = 0.0;
			(*tau_k).pos_inter[0][1] = PI;
		}
		else{
		// absTauk_u = PI & absTauk_l != 0
			tauk_pos_l = round_lf(change_referential_1_to_0(phase,  absTauk_l), 6);
			tauk_neg_l = round_lf(change_referential_1_to_0(phase, -absTauk_l), 6);
							
			sign_tauk_pos_l = sign_double(tauk_pos_l);
			sign_tauk_neg_l = sign_double(tauk_neg_l);
			if((sign_tauk_pos_l > 0) && (sign_tauk_neg_l > 0)){
				if((absTauk_u - absTauk_l) > HALFPI){
					(*tau_k).num_neg_inter = 1;								
					(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
					(*tau_k).neg_inter[0][0] = -PI;
					(*tau_k).neg_inter[0][1] = 0.0;
					
					(*tau_k).num_pos_inter = 2;
					(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
					(*tau_k).pos_inter[0][0] = 0.0;
					(*tau_k).pos_inter[0][1] = tauk_neg_l;
					(*tau_k).pos_inter[1][0] = tauk_pos_l;
					(*tau_k).pos_inter[1][1] = PI;
				}
				else{
					(*tau_k).num_neg_inter = 0;
					
					(*tau_k).num_pos_inter = 1;
					(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
					(*tau_k).pos_inter[0][0] = tauk_pos_l; //minAB_lf(tauk_neg_l, tauk_pos_l);
					(*tau_k).pos_inter[0][1] = tauk_neg_l; //maxAB_lf(tauk_neg_l, tauk_pos_l);
				}
			}
			else if((sign_tauk_pos_l < 0) && (sign_tauk_neg_l < 0)){
				if((absTauk_u - absTauk_l) > HALFPI){
					(*tau_k).num_neg_inter = 2;
					(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
					(*tau_k).neg_inter[0][0] = -PI;
					(*tau_k).neg_inter[0][1] = tauk_neg_l;
					(*tau_k).neg_inter[1][0] = tauk_pos_l;
					(*tau_k).neg_inter[1][1] = 0.0;
					
					(*tau_k).num_pos_inter = 1;
					(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
					(*tau_k).pos_inter[0][0] = 0.0;
					(*tau_k).pos_inter[0][1] = PI;
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
				(*tau_k).num_neg_inter = 1;
				(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
				(*tau_k).neg_inter[0][0] = -PI;
				(*tau_k).neg_inter[0][1] = tauk_neg_l;
				
				(*tau_k).num_pos_inter = 1;
				(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
				(*tau_k).pos_inter[0][0] = tauk_pos_l;
				(*tau_k).pos_inter[0][1] = PI;
			}
			else if((sign_tauk_pos_l < 0) && (sign_tauk_neg_l > 0)){
				(*tau_k).num_neg_inter = 1;
				(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
				(*tau_k).neg_inter[0][0] = tauk_pos_l;
				(*tau_k).neg_inter[0][1] = 0.0;
				
				(*tau_k).num_pos_inter = 1;
				(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
				(*tau_k).pos_inter[0][0] = 0.0;
				(*tau_k).pos_inter[0][1] = tauk_neg_l;
			}
		}
	}
	else{
	// absTauk_u != PI
		if(absTauk_l < MYZERO){
		// absTauk_u != PI & absTauk_l = 0
			tauk_pos_u = round_lf(change_referential_1_to_0(phase,  absTauk_u), 6);
			tauk_neg_u = round_lf(change_referential_1_to_0(phase, -absTauk_u), 6);
							
			sign_tauk_pos_u = sign_double(tauk_pos_u);
			sign_tauk_neg_u = sign_double(tauk_neg_u);
			if((sign_tauk_pos_u > 0) && (sign_tauk_neg_u > 0)){
				if((absTauk_u - absTauk_l) > HALFPI){
					(*tau_k).num_neg_inter = 1;
					(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);				
					(*tau_k).neg_inter[0][0] = -PI;
					(*tau_k).neg_inter[0][1] = 0.0;
					
					(*tau_k).num_pos_inter = 2;
					(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
					(*tau_k).pos_inter[0][0] = 0.0;
					(*tau_k).pos_inter[0][1] = tauk_pos_u;
					(*tau_k).pos_inter[1][0] = tauk_neg_u;
					(*tau_k).pos_inter[1][1] = PI;
				}
				else{
					(*tau_k).num_neg_inter = 0;
	
					(*tau_k).num_pos_inter = 1;
					(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
					(*tau_k).pos_inter[0][0] = tauk_neg_u; //minAB_lf(tauk_neg_u, tauk_pos_u);
					(*tau_k).pos_inter[0][1] = tauk_pos_u; //maxAB_lf(tauk_neg_u, tauk_pos_u);
				}
			}
			else if((sign_tauk_pos_u < 0) && (sign_tauk_neg_u < 0)){
				if((absTauk_u - absTauk_l) > HALFPI){
					(*tau_k).num_neg_inter = 2;
					(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
					(*tau_k).neg_inter[0][0] = -PI;
					(*tau_k).neg_inter[0][1] = tauk_pos_u;
					(*tau_k).neg_inter[1][0] = tauk_neg_u;
					(*tau_k).neg_inter[1][1] = 0.0;
					
					(*tau_k).num_pos_inter = 1;
					(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
					(*tau_k).pos_inter[0][0] = 0.0;
					(*tau_k).pos_inter[0][1] = PI;
				}
				else{
					(*tau_k).num_pos_inter = 0;
					
					(*tau_k).num_neg_inter = 1;
					(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
					(*tau_k).neg_inter[0][0] = tauk_neg_u; //minAB_lf(tauk_neg_u, tauk_pos_u);
					(*tau_k).neg_inter[0][1] = tauk_pos_u; //maxAB_lf(tauk_neg_u, tauk_pos_u);
				}
			}
			else if((sign_tauk_pos_u > 0) && (sign_tauk_neg_u < 0)){
				(*tau_k).num_neg_inter = 1;
				(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);					
				(*tau_k).neg_inter[0][0] = tauk_neg_u;
				(*tau_k).neg_inter[0][1] = 0.0;
				
				(*tau_k).num_pos_inter = 1;
				(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
				(*tau_k).pos_inter[0][0] = 0.0;
				(*tau_k).pos_inter[0][1] = tauk_pos_u;
			}
			else if((sign_tauk_pos_u < 0) && (sign_tauk_neg_u > 0)){
				(*tau_k).num_neg_inter = 1;
				(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
				(*tau_k).neg_inter[0][0] = -PI;
				(*tau_k).neg_inter[0][1] = tauk_pos_u;
				
				(*tau_k).num_pos_inter = 1;
				(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
				(*tau_k).pos_inter[0][0] = tauk_neg_u;
				(*tau_k).pos_inter[0][1] = PI;
			}
		}
		else{
		// absTauk_u != PI & absTauk_l != 0
			tauk_pos_u = round_lf(change_referential_1_to_0(phase,  absTauk_u), 6);
			tauk_pos_l = round_lf(change_referential_1_to_0(phase,  absTauk_l), 6);
			tauk_neg_u = round_lf(change_referential_1_to_0(phase, -absTauk_u), 6);
			tauk_neg_l = round_lf(change_referential_1_to_0(phase, -absTauk_l), 6);
			
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
							(*tau_k).num_neg_inter = 0;
							
							(*tau_k).num_pos_inter = 2;
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
							(*tau_k).num_neg_inter = 1;
							(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
							(*tau_k).neg_inter[0][0] = -PI;
							(*tau_k).neg_inter[0][1] = tauk_neg_l;
							
							(*tau_k).num_pos_inter = 2;
							(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
							(*tau_k).pos_inter[0][0] = tauk_pos_l;
							(*tau_k).pos_inter[0][1] = tauk_pos_u;
							(*tau_k).pos_inter[1][0] = tauk_neg_u;
							(*tau_k).pos_inter[1][1] = PI;
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
void change_referential_precise_case(ddgp_interval *tau_k, double absTauk_m, double phase, double resolution){
	
	double tau_k1, tau_k2;
	double tauk_pos = round_lf(change_referential_1_to_0(phase,  absTauk_m), 6);
	double tauk_neg = round_lf(change_referential_1_to_0(phase, -absTauk_m), 6);
					
	int sign_tauk_pos = sign_double(tauk_pos);
	int sign_tauk_neg = sign_double(tauk_neg);
					
	if((sign_tauk_pos > 0) && (sign_tauk_neg > 0)){
		(*tau_k).num_neg_inter = 0;
		
		tau_k1 = minAB_lf(tauk_pos, tauk_neg);
		tau_k2 = maxAB_lf(tauk_pos, tauk_neg);
		
		if(fabs(tau_k2 - tau_k1) < 2*resolution){
			(*tau_k).num_pos_inter = 1;
			(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
			
			(*tau_k).pos_inter[0][0] = maxAB_lf(minAB_lf(tau_k1 - resolution, tau_k2 - resolution), 0.0);
			(*tau_k).pos_inter[0][1] = minAB_lf(maxAB_lf(tau_k1 + resolution, tau_k2 + resolution), PI);
		}
		else{
			(*tau_k).num_pos_inter = 2;
			(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
			
			(*tau_k).pos_inter[0][0] = maxAB_lf(tau_k1 - resolution, 0.0);
			(*tau_k).pos_inter[0][1] = minAB_lf(tau_k1 + resolution, PI);
			(*tau_k).pos_inter[1][0] = maxAB_lf(tau_k2 - resolution, 0.0);
			(*tau_k).pos_inter[1][1] = minAB_lf(tau_k2 + resolution, PI);
		}
	}
	else if((sign_tauk_pos < 0) && (sign_tauk_neg < 0)){
		(*tau_k).num_pos_inter = 0;
		
		tau_k1 = minAB_lf(tauk_pos, tauk_neg);
		tau_k2 = maxAB_lf(tauk_pos, tauk_neg);
		
		if(fabs(tau_k2 - tau_k1) < 2*resolution){
			(*tau_k).num_neg_inter = 1;
			(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
			
			(*tau_k).neg_inter[0][0] = maxAB_lf(minAB_lf(tau_k1 - resolution, tau_k2 - resolution), -PI);
			(*tau_k).neg_inter[0][1] = minAB_lf(maxAB_lf(tau_k1 + resolution, tau_k2 + resolution), 0.0);
		}
		else{
			(*tau_k).num_neg_inter = 2;
			(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
			
			(*tau_k).neg_inter[0][0] = maxAB_lf(tau_k1 - resolution, -PI);
			(*tau_k).neg_inter[0][1] = minAB_lf(tau_k1 + resolution, 0.0);
			(*tau_k).neg_inter[1][0] = maxAB_lf(tau_k2 - resolution, -PI);
			(*tau_k).neg_inter[1][1] = minAB_lf(tau_k2 + resolution, 0.0);
		}
	}
	else if((sign_tauk_pos > 0) && (sign_tauk_neg < 0)){
		(*tau_k).num_pos_inter = 1;
		(*tau_k).num_neg_inter = 1;
		(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
		(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
		
		(*tau_k).pos_inter[0][0] = maxAB_lf(tauk_pos - resolution, 0.0);
		(*tau_k).pos_inter[0][1] = minAB_lf(tauk_pos + resolution, PI);
		
		(*tau_k).neg_inter[0][0] = maxAB_lf(tauk_neg - resolution, -PI);
		(*tau_k).neg_inter[0][1] = minAB_lf(tauk_neg + resolution, 0.0);
	}
	else if((sign_tauk_pos < 0) && (sign_tauk_neg > 0)){
		(*tau_k).num_pos_inter = 1;
		(*tau_k).num_neg_inter = 1;
		(*tau_k).pos_inter = alloc_mat_lf((*tau_k).num_pos_inter, 2);
		(*tau_k).neg_inter = alloc_mat_lf((*tau_k).num_neg_inter, 2);
		
		(*tau_k).pos_inter[0][0] = maxAB_lf(tauk_neg - resolution, 0.0);
		(*tau_k).pos_inter[0][1] = minAB_lf(tauk_neg + resolution, PI);
		
		(*tau_k).neg_inter[0][0] = maxAB_lf(tauk_pos - resolution, -PI);
		(*tau_k).neg_inter[0][1] = minAB_lf(tauk_pos + resolution, 0.0);
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
void tauii3_with_parameters_d2(ddgp_interval *tau, double p30, double q30, double dii3_2_l, double dii3_2_u){
	
	double absTau_l, absTau_u;
	
	// case dl = du
	if(fabs(dii3_2_l - dii3_2_u) < MYZERO){
		absTau_l = round_lf(abs_torsion_angle_with_constants_d2(p30, q30, dii3_2_l), 6);
		
		if((absTau_l < MYZERO) || (fabs(absTau_l - PI) < MYZERO)){
			(*tau).num_pos_inter = 1;
			(*tau).pos_inter = alloc_mat_lf((*tau).num_pos_inter, 2);
			(*tau).pos_inter[0][0] = absTau_l;
			(*tau).pos_inter[0][1] = absTau_l;
			
			(*tau).num_neg_inter = 0;
		}
		else{
			(*tau).num_neg_inter = 1;
			(*tau).neg_inter = alloc_mat_lf((*tau).num_neg_inter, 2);
			(*tau).neg_inter[0][0] = -absTau_l;
			(*tau).neg_inter[0][1] = -absTau_l;
				
			(*tau).num_pos_inter = 1;
			(*tau).pos_inter = alloc_mat_lf((*tau).num_pos_inter, 2);
			(*tau).pos_inter[0][0] = absTau_l;
			(*tau).pos_inter[0][1] = absTau_l;
		}
	}
	// case dl < du
	else{
		absTau_l = round_lf(abs_torsion_angle_with_constants_d2(p30, q30, dii3_2_l), 6);
		absTau_u = round_lf(abs_torsion_angle_with_constants_d2(p30, q30, dii3_2_u), 6);
	        
		(*tau).num_neg_inter = 1;
		(*tau).neg_inter = alloc_mat_lf((*tau).num_neg_inter, 2);
		(*tau).neg_inter[0][0] = -absTau_u;
		(*tau).neg_inter[0][1] = -absTau_l;	
			
		(*tau).num_pos_inter = 1;
		(*tau).pos_inter = alloc_mat_lf((*tau).num_pos_inter, 2);
		(*tau).pos_inter[0][0] = absTau_l;
		(*tau).pos_inter[0][1] = absTau_u;
	}
}
/* *********************************************************************************** */
void compute_tauii3_parameters(int i, double **X, double ***discretizationEdges_2, int *i1, int *i2, int *i3, double *di1i2, double *kappa_210, double *rho2_210, double *kappa_213, double *rho2_213, double *q30){
	
	double di1i3_2, di2i3_2, dii1_2, dii2_2;
	
	(*i1) = (int) (discretizationEdges_2[i][0][0] - 1.0);
	(*i2) = (int) (discretizationEdges_2[i][1][0] - 1.0);
	(*i3) = (int) (discretizationEdges_2[i][2][0] - 1.0);
	
	di1i3_2  = d2xixj_lf3(&X[(*i3)][0], &X[(*i1)][0]);
	(*di1i2) =  dxixj_lf3(&X[(*i2)][0], &X[(*i1)][0]);
	di2i3_2  = d2xixj_lf3(&X[(*i3)][0], &X[(*i2)][0]);
	
	dii1_2   = discretizationEdges_2[i][0][1];
	dii2_2   = discretizationEdges_2[i][1][1];
	
	(*kappa_210) = kappai_d2(dii2_2, (*di1i2), dii1_2);
	(*rho2_210)  = rhoi2_d2(dii2_2, (*kappa_210));
	
	(*kappa_213) = kappai_d2(di2i3_2, (*di1i2), di1i3_2);
	(*rho2_213)  = rhoi2_d2(di2i3_2, (*kappa_213));
	
	(*q30)  = sqrt((*rho2_210)*(*rho2_213));
}
/* *********************************************************************************** */
void compute_tauii3(int i, double **X, double ***discretizationEdges_2, int *i1, int *i2, int *i3, double *di1i2, double *kappa_210, double *rho2_210, double *kappa_213, double *rho2_213, double *q30, ddgp_interval *tau){
	
	double di1i3_2, di2i3_2, dii1_2, dii2_2, dii3_2_l, dii3_2_u, kmkb, p30;
	
	(*i1) = (int) (discretizationEdges_2[i][0][0] - 1.0);
	(*i2) = (int) (discretizationEdges_2[i][1][0] - 1.0);
	(*i3) = (int) (discretizationEdges_2[i][2][0] - 1.0);
	
	di1i3_2  = d2xixj_lf3(&X[(*i3)][0], &X[(*i1)][0]);
	(*di1i2) =  dxixj_lf3(&X[(*i2)][0], &X[(*i1)][0]);
	di2i3_2  = d2xixj_lf3(&X[(*i3)][0], &X[(*i2)][0]);
	
	dii1_2   = discretizationEdges_2[i][0][1];
	dii2_2   = discretizationEdges_2[i][1][1];
	dii3_2_l = discretizationEdges_2[i][2][1];
	dii3_2_u = discretizationEdges_2[i][2][2];
	
	(*kappa_210) = kappai_d2(dii2_2, (*di1i2), dii1_2);
	(*rho2_210)  = rhoi2_d2(dii2_2, (*kappa_210));
	
	(*kappa_213) = kappai_d2(di2i3_2, (*di1i2), di1i3_2);
	(*rho2_213)  = rhoi2_d2(di2i3_2, (*kappa_213));
	
	kmkb = (*kappa_210) - (*kappa_213);
	p30  = kmkb*kmkb + (*rho2_210) + (*rho2_213);
	(*q30)  = sqrt((*rho2_210)*(*rho2_213));
	
	tauii3_with_parameters_d2(&(*tau), p30, (*q30), dii3_2_l, dii3_2_u);
}
/* *********************************************************************************** */
void compute_tauiik_precise(int i, int k, double **X, prune_edges_set *pruneEdges_2, double kappa_210, double rho2_210, double kappa_213, double rho2_213, double di1i2, int i3, int i2, int i1, ddgp_interval *tau_k, double resolution){
	
	int ik;
	double di1ik_2, di2ik_2, di3ik_2, diik_2_l;
	double kappa_21k, rho2_21k, kmkb, pk3, qk3, phase, pk0, qk0;
	double absTauk_m;
	
	ik = (int) (pruneEdges_2[i].precise[k][0] - 1.0);
	
	di1ik_2 = d2xixj_lf3(&X[ik][0], &X[i1][0]);
	di2ik_2 = d2xixj_lf3(&X[ik][0], &X[i2][0]);
	di3ik_2 = d2xixj_lf3(&X[ik][0], &X[i3][0]);
	
	kappa_21k = kappai_d2(di2ik_2, di1i2, di1ik_2);
	rho2_21k  = rhoi2_d2(di2ik_2, kappa_21k);
	
	kmkb = kappa_213 - kappa_21k;
	pk3  = kmkb*kmkb + rho2_21k + rho2_213;
	qk3  = sqrt(rho2_21k*rho2_213);
	
	phase = sign_torsion_angle(&X[i3][0], &X[i2][0], &X[i1][0], &X[ik][0])*abs_torsion_angle_with_constants_d2(pk3, qk3, di3ik_2);
	
	kmkb = kappa_210 - kappa_21k;
	pk0  = kmkb*kmkb + rho2_210 + rho2_21k;
	qk0  = sqrt(rho2_210*rho2_21k);
						
	diik_2_l = pruneEdges_2[i].precise[k][1];
	
	absTauk_m = abs_torsion_angle_with_constants_d2(pk0, qk0, diik_2_l);
	
	change_referential_precise_case(&(*tau_k), absTauk_m, phase, resolution);
}
/* *********************************************************************************** */
void compute_tauiik_interval(int i, int k, double **X, prune_edges_set *pruneEdges_2, double kappa_210, double rho2_210, double kappa_213, double rho2_213, double di1i2, int i3, int i2, int i1, ddgp_interval *tau_k, int *isTauIntervalEmpty){
	
	int ik;
	double di1ik_2, di2ik_2, di3ik_2, diik_2_l, diik_2_u;
	double kappa_21k, rho2_21k, kmkb, pk3, qk3, phase, pk0, qk0;
	double absTauk_l, absTauk_u;
		
	ik = (int) (pruneEdges_2[i].interval[k][0] - 1.0);
	
	di1ik_2 = d2xixj_lf3(&X[ik][0], &X[i1][0]);
	di2ik_2 = d2xixj_lf3(&X[ik][0], &X[i2][0]);
	di3ik_2 = d2xixj_lf3(&X[ik][0], &X[i3][0]);
		
	kappa_21k = kappai_d2(di2ik_2, di1i2, di1ik_2);
	rho2_21k  = rhoi2_d2(di2ik_2, kappa_21k);
	
	kmkb = kappa_213 - kappa_21k;
	pk3  = kmkb*kmkb + rho2_21k + rho2_213;
	qk3  = sqrt(rho2_21k*rho2_213);
	
	phase = sign_torsion_angle(&X[i3][0], &X[i2][0], &X[i1][0], &X[ik][0])*abs_torsion_angle_with_constants_d2(pk3, qk3, di3ik_2);
	
	kmkb = kappa_210 - kappa_21k;
	pk0  = kmkb*kmkb + rho2_210 + rho2_21k;
	qk0  = sqrt(rho2_210*rho2_21k);
			
	diik_2_l = pruneEdges_2[i].interval[k][1];
	diik_2_u = pruneEdges_2[i].interval[k][2];
	
	absTauk_l = abs_torsion_angle_with_constants_d2(pk0, qk0, diik_2_l);
	absTauk_u = abs_torsion_angle_with_constants_d2(pk0, qk0, diik_2_u);
	
	if(fabs(absTauk_u - absTauk_l) < MYZERO)
		(*isTauIntervalEmpty) = 1;
	else
		change_referential_interval_case(&(*tau_k), absTauk_l, absTauk_u, phase);
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
void print_interval(ddgp_interval tau, int opt1, int opt2, int k){
	
	int l;
	if(opt1){ // tau
		if(opt2){ // radians
			print("");
			if(tau.num_neg_inter > 0)
				for(l = 0; l < tau.num_neg_inter; l++)
					printf("tau^- = [%9.6lf %9.6lf]\n", tau.neg_inter[l][0], tau.neg_inter[l][1]);
			else
				printf("tau^- = {}\n");

			if(tau.num_pos_inter > 0)		
				for(l = 0; l < tau.num_pos_inter; l++)
					printf("tau^+ = [%9.6lf %9.6lf]\n", tau.pos_inter[l][0], tau.pos_inter[l][1]);
			else
				printf("tau^+ = {}\n");
			print("");
		}
		else{ // degree
			print("");
			if(tau.num_neg_inter > 0)
				for(l = 0; l < tau.num_neg_inter; l++)
					printf("tau^- = [%9.6lf %9.6lf]\n", tau.neg_inter[l][0]*180/PI, tau.neg_inter[l][1]*180/PI);
			else
				printf("tau^- = {}\n");

			if(tau.num_pos_inter > 0)		
				for(l = 0; l < tau.num_pos_inter; l++)
					printf("tau^+ = [%9.6lf %9.6lf]\n", tau.pos_inter[l][0]*180/PI, tau.pos_inter[l][1]*180/PI);
			else
				printf("tau^+ = {}\n");
			print("");
		}
	}
	else{ // tau_k
		if(opt2){ // radians
			print("");
			if(tau.num_neg_inter > 0)
				for(l = 0; l < tau.num_neg_inter; l++)
					printf("tau_%d^- = [%9.6lf %9.6lf]\n", k, tau.neg_inter[l][0], tau.neg_inter[l][1]);
			else
				printf("tau_%d^- = {}\n", k);

			if(tau.num_pos_inter > 0)		
				for(l = 0; l < tau.num_pos_inter; l++)
					printf("tau_%d^+ = [%9.6lf %9.6lf]\n", k, tau.pos_inter[l][0], tau.pos_inter[l][1]);
			else
				printf("tau_%d^+ = {}\n", k);
			print("");
		}
		else{ // degree
			print("");
			if(tau.num_neg_inter > 0)
				for(l = 0; l < tau.num_neg_inter; l++)
					printf("tau_%d^- = [%9.6lf %9.6lf]\n", k, tau.neg_inter[l][0]*180/PI, tau.neg_inter[l][1]*180/PI);
			else
				printf("tau_%d^- = {}\n", k);

			if(tau.num_pos_inter > 0)		
				for(l = 0; l < tau.num_pos_inter; l++)
					printf("tau_%d^+ = [%9.6lf %9.6lf]\n", k, tau.pos_inter[l][0]*180/PI, tau.pos_inter[l][1]*180/PI);
			else
				printf("tau_%d^+ = {}\n", k);
			print("");
		}
	}
}
/* *********************************************************************************** */
void get_torsion_angle(ddgp_interval *tau, double tauAst, double deltaTau){

	double tau_neg = tauAst - deltaTau;
	double tau_pos = tauAst + deltaTau;
	
	int sign_tau_neg = sign_double(tau_neg);
	int sign_tau_pos = sign_double(tau_pos);
	
	if((sign_tau_pos > 0) && (sign_tau_neg > 0)){
		if(tau_pos > PI){
			(*tau).num_neg_inter = 1;
			(*tau).neg_inter = alloc_mat_lf((*tau).num_neg_inter, 2);
			(*tau).neg_inter[0][0] = -PI;
			(*tau).neg_inter[0][1] = tau_pos - TWOPI;
			
			(*tau).num_pos_inter = 1;
			(*tau).pos_inter = alloc_mat_lf((*tau).num_pos_inter, 2);
			(*tau).pos_inter[0][0] = tau_neg;
			(*tau).pos_inter[0][1] = PI;
		}
		else{
			(*tau).num_neg_inter = 0;
			
			(*tau).num_pos_inter = 1;
			(*tau).pos_inter = alloc_mat_lf((*tau).num_pos_inter, 2);
			(*tau).pos_inter[0][0] = tau_neg;
			(*tau).pos_inter[0][1] = tau_pos;
		}
	}
	else if((sign_tau_pos < 0) && (sign_tau_neg < 0)){
		if(tau_neg < -PI){
			(*tau).num_neg_inter = 1;
			(*tau).neg_inter = alloc_mat_lf((*tau).num_neg_inter, 2);
			(*tau).neg_inter[0][0] = -PI;
			(*tau).neg_inter[0][1] = tau_pos;
			
			(*tau).num_pos_inter = 1;
			(*tau).pos_inter = alloc_mat_lf((*tau).num_pos_inter, 2);
			(*tau).pos_inter[0][0] = tau_neg + TWOPI;
			(*tau).pos_inter[0][1] = PI;
		}
		else{
			(*tau).num_neg_inter = 1;
			(*tau).neg_inter = alloc_mat_lf((*tau).num_neg_inter, 2);
			(*tau).neg_inter[0][0] = tau_neg;
			(*tau).neg_inter[0][1] = tau_pos;
			
			(*tau).num_pos_inter = 0;
		}
	}
	else if((sign_tau_pos > 0) && (sign_tau_neg < 0)){
		(*tau).num_neg_inter = 1;
		(*tau).neg_inter = alloc_mat_lf((*tau).num_neg_inter, 2);
		(*tau).neg_inter[0][0] = tau_neg;
		(*tau).neg_inter[0][1] = 0.0;
		
		(*tau).num_pos_inter = 1;
		(*tau).pos_inter = alloc_mat_lf((*tau).num_pos_inter, 2);
		(*tau).pos_inter[0][0] = 0.0;
		(*tau).pos_inter[0][1] = tau_pos;
	}
	else if((sign_tau_pos < 0) && (sign_tau_neg > 0))
		printf("Erro\n");
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
	
	//if(sum == 0)
	//	printf("X Satisfies the Instance\n");
	//else{
	if(sum > 0){
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
/* *********************************************************************************** */
void compute_mde_and_lde(double **X, dcinputfile *dc_vec, int m, double *mde, double *lde){

	int i, j, k, m0, k0;
	double dij, dijL, dijU;
	int *flagsVec = zeros_vec_d(m);
	m0 = 0;
	for(k = 0; k < m; k++)
		if((dc_vec[k].dl < dc_vec[k].du) && (dc_vec[k].du < 900.00)){
			flagsVec[k] = 1;
			m0++;
		}
	
	(*lde) = 0.0;
	(*mde) = 0.0;
	
	if(m0 > 0){
		double *deVec = zeros_vec_lf(m0);
		k0 = 0;
		for(k = 0; k < m; k++)
			if(flagsVec[k]){
				dijU = dc_vec[k].du;
				dijL = dc_vec[k].dl;
				
				i = dc_vec[k].i - 1;
				j = dc_vec[k].j - 1;
				dij = dxixj_lf3(&X[i][0], &X[j][0]);
				
				deVec[k0++] = maxAB_lf(fabs(dijL - dij), fabs(dijU - dij));
			}
		
		(*lde) = max_val_vec_lf(deVec, m0);
		(*mde) = sum_vec_lf(deVec, m0)/m0;
		
		free(deVec);
	}
	
	free(flagsVec);
}
/* *********************************************************************************** */
double compute_rmsd(double **X, double **Y, int n){ // X = X0, Y = Xr;

	double **X0 = alloc_mat_lf(n, 3);
	mat_e_mat_lf(X0, X, n, 3);
	
	double **Xr = alloc_mat_lf(n, 3);
	mat_e_mat_lf(Xr, Y, n, 3);

	double cm_x_Xr, cm_y_Xr, cm_z_Xr, cm_x_X0, cm_y_X0, cm_z_X0;
	double *vec_aux;
	
	// compute the CM of Xr
	vec_aux = get_column_mat_lf(Xr, n, 1);
	cm_x_Xr = sum_vec_lf(vec_aux, n)/n;
	free(vec_aux);
	
	vec_aux = get_column_mat_lf(Xr, n, 2);
	cm_y_Xr = sum_vec_lf(vec_aux, n)/n;
	free(vec_aux);
	
	vec_aux = get_column_mat_lf(Xr, n, 3);
	cm_z_Xr = sum_vec_lf(vec_aux, n)/n;
	free(vec_aux);
	
	// compute the CM of X0
	vec_aux = get_column_mat_lf(X0, n, 1);
	cm_x_X0 = sum_vec_lf(vec_aux, n)/n;
	free(vec_aux);
	
	vec_aux = get_column_mat_lf(X0, n, 2);
	cm_y_X0 = sum_vec_lf(vec_aux, n)/n;
	free(vec_aux);
	
	vec_aux = get_column_mat_lf(X0, n, 3);
	cm_z_X0 = sum_vec_lf(vec_aux, n)/n;
	free(vec_aux);
	
	// subtract the CM of the two structures
	int i;
	for(i = 0; i < n; i++){
		Xr[i][0] -= cm_x_Xr;
		Xr[i][1] -= cm_y_Xr;
		Xr[i][2] -= cm_z_Xr;
		
		X0[i][0] -= cm_x_X0;
		X0[i][1] -= cm_y_X0;
		X0[i][2] -= cm_z_X0;
	}
	
	// compute Xr^tX0
	double **Xrt = trans_mat_lf(Xr, n, 3);
	double **Xrt_x_X0 = mat_times_mat_lf(Xrt, X0, 3, n, 3);
	dealloc_mat_lf(Xrt, 3);
		
	// compute USV^t = svd(Xr^tX0)
	double **U = alloc_mat_lf(3, 3);
	double **S = alloc_mat_lf(3, 3);
	double **Vt = alloc_mat_lf(3, 3);
	
	svd_row_major_matrix(Xrt_x_X0, 3, 3, U, S, Vt);
	dealloc_mat_lf(Xrt_x_X0, 3);
	dealloc_mat_lf(S, 3);
	// compute Q = V*Ut
	double **Ut = trans_mat_lf(U, 3, 3);
	dealloc_mat_lf(U, 3);
	double **V  = trans_mat_lf(Vt, 3, 3);
	dealloc_mat_lf(Vt, 3);
	double **Q = mat_times_mat_lf(V, Ut, 3, 3, 3);
	dealloc_mat_lf(Ut, 3);
	dealloc_mat_lf(V, 3);
	
	// compute XrQ^t
	double **Qt = trans_mat_lf(Q, 3, 3);
	dealloc_mat_lf(Q, 3);
	double **Xr_x_Qt = mat_times_mat_lf(Xr, Qt, n, 3, 3);
	dealloc_mat_lf(Qt, 3);

	// compute X0 - XrQ^t
	double **X0_m_XrQt = zeros_mat_lf(n, 3);
	mat_m_mat_lf(X0_m_XrQt, X0, Xr_x_Qt, n, 3);
	dealloc_mat_lf(Xr_x_Qt, n);
	
	double rmsd = frobenius_norm(X0_m_XrQt, n, 3)/sqrt(n);
	
	dealloc_mat_lf(X0_m_XrQt, n);
	dealloc_mat_lf(X0, n);
	dealloc_mat_lf(Xr, n);
	
	return rmsd;
}
/* *********************************************************************************** */
void iABP(int n, double ***discretizationEdges_2, int *tauSign, double *givenTau, double *givenTauDeviation, prune_edges_set *pruneEdges_2, int sampleSize, double resolution, double tolerance, double timeLimit, proteinstructure **protein, int GivenNumOfSols, solution_metrics *solutionMetrics, dcinputfile *dc_vec, int num_dc, int is_X0_known, double **X0){
	
	// ---------- Variables initialization ----------
	clock_t startTime, start, end;
	
	int j, k, isTauIntervalEmpty;
	int i = 3;
	int lev = 2;
	int time_limit_not_reached = 1;
	int twoSampleSize = 2*sampleSize;
	int *exploredVertex = zeros_vec_d(n);
	int *branchNum = zeros_vec_d(n);
	int i1, i2, i3;
	
	long nos = 0;
	
	double kappa_210, rho2_210, kappa_213, rho2_213, q30, di1i2;
	double mde, lde, rmsd;
	double npp = 0.000003;
	double rmsdMin = 1000.0;
	double **A  = zeros_mat_lf(n, 3);
	double **B  = zeros_mat_lf(n, 3);
	double **C  = zeros_mat_lf(n, 3);
	double **Xr = zeros_mat_lf(n, 3);
	double **branches = ones_mat_lf(n, twoSampleSize);
	l_t_mat_lf(404.0, branches, n, twoSampleSize);
	
	ddgp_interval tau, tau_k;

	(*solutionMetrics).RMSD = stack_create();
	(*solutionMetrics).MDE = stack_create();
	(*solutionMetrics).LDE = stack_create();

	referential_x1_x2_x3(&Xr, discretizationEdges_2);
	
	start = clock();
	startTime = time(NULL);
	while(time_limit_not_reached){
		if(difftime(time(NULL), startTime) > timeLimit)
			time_limit_not_reached = 0;
		
		if(i == n){
			satisfiesInstanceQM(Xr, dc_vec, num_dc, tolerance);
			
			if(is_X0_known){
				rmsd = compute_rmsd(X0, Xr, n);
				stack_push((*solutionMetrics).RMSD, rmsd);
				
				if(rmsdMin > rmsd){
					rmsdMin = rmsd;
					for(j = 0; j < n; j++){
						(*protein)[j].x = Xr[j][0];
						(*protein)[j].y = Xr[j][1];
						(*protein)[j].z = Xr[j][2];
					}

				}
			}
			compute_mde_and_lde(Xr, dc_vec, num_dc , &mde, &lde);
			stack_push((*solutionMetrics).MDE, mde);
			stack_push((*solutionMetrics).LDE, lde);
			nos++;
			
			if((nos == GivenNumOfSols) && (GivenNumOfSols != 0))
				break;
			
			i--;
			tree_backtracking(&i, &exploredVertex, &branches, &branchNum, twoSampleSize);
			if(i == 2){
				printf("iABP: ERROR! The instance solution was not found\n");
				printf("iABP: The whole search space was explored\n");
				break;
			}
		}
		
		if(exploredVertex[i] == 0){
			isTauIntervalEmpty = 0;
			// Compute tau_{i,i3}
			if(tauSign[i] == 0){
				compute_tauii3(i, Xr, discretizationEdges_2, &i1, &i2, &i3, &di1i2, &kappa_210, &rho2_210, &kappa_213, &rho2_213, &q30, &tau);
				for(k = 0; k < pruneEdges_2[i].cardUkP; k++){
					// Compute tau_{i,ik}
					compute_tauiik_precise(i, k, Xr, pruneEdges_2, kappa_210, rho2_210, kappa_213, rho2_213, di1i2, i3, i2, i1, &tau_k, resolution);
					
					cap_interAk(&tau, tau_k);
					
					dealloc_tau(tau_k);
					
					if((tau.num_pos_inter + tau.num_neg_inter) == 0){
						isTauIntervalEmpty = 1;
						break;
					}
				}
			}
			else{
				compute_tauii3_parameters(i, Xr, discretizationEdges_2, &i1, &i2, &i3, &di1i2, &kappa_210, &rho2_210, &kappa_213, &rho2_213, &q30);
				get_torsion_angle(&tau, tauSign[i]*givenTau[i], givenTauDeviation[i]);
			}
			
			if(isTauIntervalEmpty == 0){
				for(k = 0; k < pruneEdges_2[i].cardUkI; k++){
					// Compute tau_{i,ik}
					compute_tauiik_interval(i, k, Xr, pruneEdges_2, kappa_210, rho2_210, kappa_213, rho2_213, di1i2, i3, i2, i1, &tau_k, &isTauIntervalEmpty);
					
					if(isTauIntervalEmpty == 1){
						dealloc_tau(tau);
						break;
					}
					
					cap_interAk(&tau, tau_k);
					
					dealloc_tau(tau_k);
					
					if((tau.num_pos_inter + tau.num_neg_inter) == 0){
						isTauIntervalEmpty = 1;
						break;
					}
				}
			}
			
			if(isTauIntervalEmpty == 0){
				sample_interval_union_DDGP(i, &branches, &branchNum, tau, sampleSize, resolution);
				
				dealloc_tau(tau);
				
				exploredVertex[i] = 1;
				
				circunference_parameters(i, &A, &B, &C, &Xr[i3][0], &Xr[i2][0], &Xr[i1][0], kappa_210, kappa_213, rho2_210, q30, di1i2);
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
			
			npp+= 0.000001;
			i++;
			if(lev < i)
				lev = i;
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
	
	dealloc_mat_lf(A, n);
	dealloc_mat_lf(B, n);
	dealloc_mat_lf(C, n);
	dealloc_mat_lf(Xr, n);
	
	dealloc_mat_lf(branches, n);
	free(exploredVertex);
	free(branchNum);
	
	end = clock();
	
	(*solutionMetrics).cpu_time_used = ((double) (end - start))/CLOCKS_PER_SEC;
	(*solutionMetrics).lev = lev;
	(*solutionMetrics).npp = 1000000*npp;
	(*solutionMetrics).nos = nos;
	
	if((*solutionMetrics).cpu_time_used > timeLimit - 0.5){
		printf("iABP: Time Limit Reached.\n");
		if(lev < n)
			printf("iABP: ERROR! The instance solution was not found\n");
		
	}
}
/* *********************************************************************************** */
void iBP(int n, double ***discretizationEdges_2, prune_edges_set *pruneEdges_2, int sampleSize, double tolerance, double timeLimit, proteinstructure **protein, int GivenNumOfSols, solution_metrics *solutionMetrics, dcinputfile *dc_vec, int num_dc, int is_X0_known, double **X0){
	// ---------- Variables initialization ----------
	clock_t startTime, start, end;
	
	int j, k, satisfiedPruneEdges;
	int i = 3;
	int lev = 2;
	int time_limit_not_reached = 1;
	int twoSampleSize = 2*sampleSize;
	int *exploredVertex = zeros_vec_d(n);
	int *branchNum = zeros_vec_d(n);
	int i1, i2, i3, ik;
	
	long nos = 0;
	
	double kappa_210, rho2_210, kappa_213, rho2_213, q30, di1i2, diik_2_0, diik_2_l, diik_2_u;
	double mde, lde, rmsd;
	double npp = 0.000003;
	double rmsdMin = 1000.0;
	double **A  = zeros_mat_lf(n, 3);
	double **B  = zeros_mat_lf(n, 3);
	double **C  = zeros_mat_lf(n, 3);
	double **Xr = zeros_mat_lf(n, 3);
	double **branches = ones_mat_lf(n, twoSampleSize);
	l_t_mat_lf(404.0, branches, n, twoSampleSize);
	
	ddgp_interval tau;
	
	(*solutionMetrics).RMSD = stack_create();
	(*solutionMetrics).MDE = stack_create();
	(*solutionMetrics).LDE = stack_create();
	
	referential_x1_x2_x3(&Xr, discretizationEdges_2);
	
	start = clock();
	startTime = time(NULL);
	while(time_limit_not_reached){
		if(difftime(time(NULL), startTime) > timeLimit)
		time_limit_not_reached = 0;
	
		if(i == n){
			satisfiesInstanceQM(Xr, dc_vec, num_dc, tolerance);
			
			if(is_X0_known){
				rmsd = compute_rmsd(X0, Xr, n);
				stack_push((*solutionMetrics).RMSD, rmsd);
				
				if(rmsdMin > rmsd){
					rmsdMin = rmsd;
					for(j = 0; j < n; j++){
						(*protein)[j].x = Xr[j][0];
						(*protein)[j].y = Xr[j][1];
						(*protein)[j].z = Xr[j][2];
					}
				}
			}
			compute_mde_and_lde(Xr, dc_vec, num_dc , &mde, &lde);
			stack_push((*solutionMetrics).MDE, mde);
			stack_push((*solutionMetrics).LDE, lde);
			nos++;
			
			if((nos == GivenNumOfSols) && (GivenNumOfSols != 0))
				break;
			
			i--;
			tree_backtracking(&i, &exploredVertex, &branches, &branchNum, twoSampleSize);
			if(i == 2){
				printf("iBP: ERROR! The instance solution was not found\n");
				printf("iBP: The whole search space was explored\n");
				break;
			}
		}
		
		if(exploredVertex[i] == 0){
			// Compute tau_{i,i3}
			compute_tauii3(i, Xr, discretizationEdges_2, &i1, &i2, &i3, &di1i2, &kappa_210, &rho2_210, &kappa_213, &rho2_213, &q30, &tau);
		
			sample_simple_DDGPinterval(i, &branches, &branchNum, tau, sampleSize, tolerance);	
			
			dealloc_tau(tau);
			
			exploredVertex[i] = 1;
			circunference_parameters(i, &A, &B, &C, &Xr[i3][0], &Xr[i2][0], &Xr[i1][0], kappa_210, kappa_213, rho2_210, q30, di1i2);
		}
		else
			branches[i][branchNum[i]++] = 404.0;
		
		if(branchNum[i] < twoSampleSize){
			vertex_embedding(i, &Xr, &A[i][0], &B[i][0], &C[i][0], branches[i][branchNum[i]]);
			npp += 0.000001;
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
				
		satisfiedPruneEdges = 1;
		
		for(k = 0; k < pruneEdges_2[i].cardUkP; k++){
			ik = (int) (pruneEdges_2[i].precise[k][0] - 1.0);
				
			diik_2_0 = d2xixj_lf3(&Xr[i][0], &Xr[ik][0]);
				
			diik_2_l = pruneEdges_2[i].precise[k][1];
			    	
			if(fabs(diik_2_l - diik_2_0) > tolerance){
				satisfiedPruneEdges = 0;
				break;
			}
		}
		
		if(satisfiedPruneEdges){
			for(k = 0; k < pruneEdges_2[i].cardUkI; k++){
				ik = (int) (pruneEdges_2[i].interval[k][0] - 1.0);
				
				diik_2_0 = d2xixj_lf3(&Xr[i][0], &Xr[ik][0]);
				
				diik_2_l = pruneEdges_2[i].interval[k][1];
				diik_2_u = pruneEdges_2[i].interval[k][2];
				
				if((diik_2_0 < fabs(diik_2_l - tolerance)) || (fabs(diik_2_u + tolerance) < diik_2_0)){
					satisfiedPruneEdges = 0;
					break;
				}
			}
		}
		
		if(satisfiedPruneEdges){
			i++;
			if(lev < i)
				lev = i;
		}
	}
	
	dealloc_mat_lf(A, n);
	dealloc_mat_lf(B, n);
	dealloc_mat_lf(C, n);
	dealloc_mat_lf(Xr, n);
	
	dealloc_mat_lf(branches, n);
	free(exploredVertex);
	free(branchNum);

	end = clock();
	
	(*solutionMetrics).cpu_time_used = ((double) (end - start))/CLOCKS_PER_SEC;
	(*solutionMetrics).lev = lev;
	(*solutionMetrics).npp = 1000000*npp;
	(*solutionMetrics).nos = nos;
	
	if((*solutionMetrics).cpu_time_used > timeLimit - 0.5){
		printf("iBP: Time Limit Reached.\n");
		if(lev < n)
			printf("iBP: ERROR! The instance solution was not found\n");
		
	}
}
/* *********************************************************************************** */
void copy_protein_structure(proteinstructure **p1, proteinstructure *p0, int n){

	int k;
	(*p1) = alloc_vec_proteinstructure(n);
	for(k = 0; k < n; k++){
		(*p1)[k].v_k = p0[k].v_k;
		(*p1)[k].res_seq = p0[k].res_seq;
		strcpy((*p1)[k].atom_name, p0[k].atom_name);
		strcpy((*p1)[k].res_name, p0[k].res_name);
		(*p1)[k].x = p0[k].x;
		(*p1)[k].y = p0[k].y;
		(*p1)[k].z = p0[k].z;
	}
}