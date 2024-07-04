#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <lapacke.h>
#include <cblas.h>
#include "ADT_vectors.h"
#include "ADT_matrices.h"
#include "ADT_files.h"
#include "ADT_strings.h"
#include "ADT_stack.h"
#include "ADT_DDGP.h"
#define PI 3.14159265359
#define myZero 0.000001
#define MAX_LINE_LENGTH 128

//int main(int argc, char *argv[]){
int main(){

	char *hcOrder[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
	char hco[MAX_LINE_LENGTH];
	strcpy(hco, hcOrder[8]);
	
	int sampleSize = 5;
	double timeMAX = 1.0*1.0*10.0;
	double tolerance = 0.0001;
	double resolution = 0.01;
	
	char dirname[] = "dataset/";
	
	/*
	char *idCodes[] = {	"1TOS", "1UAO", "1KUW",  // n.n.a = 10
				"1ID6", "1DNG", "1O53",  // n.n.a = 15
				"1DU1", "1DPK", "1HO7",  // n.n.a = 20
				"1CKZ", "1LFC", "1A11",  // n.n.a = 25
				"1HO0", "1MMC", "1D0R",  // n.n.a = 30
				"1ZWD", "1D1H", "1SPF",  // n.n.a = 35
				"1AML", "1BA4", "1C56"}; // n.n.a = 40
	*/
	char *idCodes[] = {"1TOS"};
	
	int num_structures = sizeof(idCodes)/sizeof(idCodes[0]);
	
	char *dirname_i, *pstr, *strI, *strX, *strT;
	char aux_str[MAX_LINE_LENGTH];
	char instanceFile[MAX_LINE_LENGTH];
	char initialStructureFile[MAX_LINE_LENGTH];
	char cliquesFile[MAX_LINE_LENGTH];
		
	int lev;
	double numIt;
	double cpu_time_used;
	int vec_cliques[] = {1, 2, 3, 4};
	
	int m, n, i, k;
	int *kv;
	int **cliques;
	
	double **T0, **cliques_lf, **X0, **Xiabp, **Xibp;
	double ***instance_adjlist, ***discretationEdges;
	
	dcinputfile *dc_vec;
	
	prune_edges_set *pruneEdges;
	
	proteinstructure *protein;//, *protein_iabp, *protein_ibp;
	
	for(i = 0; i < num_structures; i++){
		strcpy(aux_str, dirname);
		dirname_i = strcat(strcat(aux_str, idCodes[i]), "/dIw_1_aIw_40/");
		
		strcpy(instanceFile, dirname_i);
		strcpy(initialStructureFile, dirname_i);
		strcpy(cliquesFile, dirname_i);
		
		strcpy(aux_str, idCodes[i]);
		pstr = strcat(strcat(strcat(aux_str, "_model1_chainA_ddgpHCorder"), hco), ".dat");
		
		strI = strcat(strcat(instanceFile, "I_"), pstr);
		strX = strcat(strcat(initialStructureFile, "X_"), pstr);
		strT = strcat(strcat(cliquesFile, "T_"), pstr);
		
		printf("%s\n", strX);
		printf("%s\n", strT);
		printf("%s\n", strI);
		
		X0 = file_2_mat_lf(initialStructureFile, ' ');
		T0 = file_2_mat_lf(cliquesFile, ' ');
		
		m = num_lines_file(instanceFile);
		
		dc_vec = dcInputFile_2_dcDataVector(instanceFile);
		
		dcinputfile_2_instance(dc_vec, &instance_adjlist, &protein, &kv, m, &n);
		
		cliques_lf = get_columns_mat_lf(T0, n, vec_cliques, 4);
		dealloc_mat_lf(T0, n);
		
		cliques = mat_lf_2_mat_d(cliques_lf, n, 4);
		dealloc_mat_lf(cliques_lf, n);
		
		adjlist_2_graph_parameters(instance_adjlist, kv, n, cliques, &discretationEdges, &pruneEdges);
		
		dealloc_3Darray_lf(instance_adjlist, kv, n);
		free(kv);
		
		printf("\nidcode = %s\n", idCodes[i]);
		
		print("****************** iABP ******************");
		Xiabp = alloc_mat_lf(n, 3);
		iABP(&Xiabp, n, &lev, &cpu_time_used, &numIt, discretationEdges, pruneEdges, sampleSize, myZero, timeMAX, tolerance);
		printf("cpu_time_used = %.8lf\n", cpu_time_used);
		printf("l.e.v = %d/%d\n", lev, n);
		numIt *= 1000000;
		printf("numIt = %.0lf\n", numIt);
		if(Xiabp[0][0] < 0){
			satisfiesInstanceQM(Xiabp, dc_vec, m, resolution);
			print_mat_lf(Xiabp, n, 3);
		}
		dealloc_mat_lf(Xiabp, n);
		
		print("****************** iBP ******************");
		Xibp = alloc_mat_lf(n, 3);
		iBP(&Xibp, n, &lev, &cpu_time_used, &numIt, discretationEdges, pruneEdges, sampleSize, myZero, timeMAX, tolerance);
		printf("cpu_time_used = %.8lf\n", cpu_time_used);
		printf("l.e.v = %d/%d\n", lev, n);
		numIt *= 1000000;
		printf("numIt = %.0lf\n", numIt);
		if(Xibp[0][0] < 0){
			satisfiesInstanceQM(Xibp, dc_vec, m, resolution);
			print_mat_lf(Xibp, n, 3);
		}
		dealloc_mat_lf(Xibp, n);
		
		print("-----------------------------------------");
		free(dc_vec);
		
		free(protein);
		
		for(k = 0; k < n; k++){
			if(pruneEdges[k].cardUkP > 0)
				dealloc_mat_lf(pruneEdges[k].precise, pruneEdges[k].cardUkP);
			if(pruneEdges[k].cardUkI > 0)
				dealloc_mat_lf(pruneEdges[k].interval, pruneEdges[k].cardUkI);
		}
		
		free(pruneEdges);
		
		dealloc_fix3Darray_lf(discretationEdges, n, 3);
	}
	
	dealloc_mat_lf(X0, n);
	
	dealloc_mat_d(cliques, n);
	
	return 0;
}
