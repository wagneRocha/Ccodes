#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
#define MYZERO 0.000001
#define CF_RAD2DEG 57.29577951308
#define CF_DEG2RAD 0.01745329252
#define MAX_LINE_LENGTH 512

int main(int argc, char *argv[]){

	if(argc < 3){
		printf("Usage: %s "
		"<%%s: input_filename_path> "
		"<%%s: output_folder_path> \n", argv[0]);
		return 0;
	}

	char *fname = argv[1];
	char *outputFolder = argv[2];
	char *structure_id = NULL;
	char *structure_chain = NULL;
	char *method = NULL;
	char *instanceFile = NULL;
	char *cliquesFile = NULL;
	char *initialStructureFile = NULL;
	double timeLimit = 0.0;
	double tolerance = 0.0;
	double angularResolution = 0.0;
	int numOfSols = 0; 
	int sampleSize = 0;

	read_input_file(fname, &structure_id, &structure_chain, &method, &instanceFile, &cliquesFile, &initialStructureFile, &timeLimit, &tolerance, &angularResolution, &numOfSols, &sampleSize);
	
	char aux_str[MAX_LINE_LENGTH];
	
	int m, n, k;
	int is_X0_knownQM = 0;
	int *kv, *tauSign;
	int **cliques;
	int vec_cliques[] = {1, 2, 3, 4};
	
	double *mde_vec, *lde_vec, *rmsd_vec; 
	double *tauSign_lf, *givenTau, *givenTauDeviation;
	double **T0, **cliques_lf, **X0;
	double ***instance_adjlist, ***discretizationEdges_2;
	
	dcinputfile *dc_vec;
	
	prune_edges_set *pruneEdges_2;
	
	proteinstructure *protein;
	
	solution_metrics solutionMetrics;
	
	m = num_lines_file(instanceFile);
	
	dc_vec = dcInputFile_2_dcDataVector(instanceFile);
	
	dcinputfile_2_instance(dc_vec, &instance_adjlist, &protein, &kv, m, &n);

	if(initialStructureFile != NULL){
		is_X0_knownQM = 1;

		X0 = file_2_mat_lf(initialStructureFile, ' ');

		proteinstructure *protein_pdb;

		copy_protein_structure(&protein_pdb, protein, n);
		for(k = 0; k < n; k++){
			protein_pdb[k].x = X0[k][0];
			protein_pdb[k].y = X0[k][1];
			protein_pdb[k].z = X0[k][2];
		}

		strcpy(aux_str, outputFolder);
		save_protein_PDBformat("PDB", protein_pdb, n, structure_id, "1", structure_chain, strcat(strcat(strcat(aux_str, "/"), structure_id), "_0.pdb"));

		free(protein_pdb);
	}
	else
		X0 = zeros_mat_lf(n,3);
	
	if(cliquesFile != NULL)
		T0 = file_2_mat_lf(cliquesFile, ' ');
	else{
		printf("Error: the path to the cliques file could not be found\n");
		return 0;
	}
	
	cliques_lf = get_columns_mat_lf(T0, n, vec_cliques, 4);
	cliques = mat_lf_2_mat_d(cliques_lf, n, 4);
	dealloc_mat_lf(cliques_lf, n);
	
	tauSign_lf = get_column_mat_lf(T0, n, 5);
	tauSign = vec_lf_2_vec_d(tauSign_lf, n);
	//tauSign = zeros_vec_d(n);
	free(tauSign_lf);
	
	givenTau = get_column_mat_lf(T0, n, 6);
	l_t_vec_lf(givenTau, CF_DEG2RAD, givenTau, n);
	
	givenTauDeviation = get_column_mat_lf(T0, n, 7);
	l_t_vec_lf(givenTauDeviation, CF_DEG2RAD, givenTauDeviation, n);
	
	dealloc_mat_lf(T0, n);
	
	adjlist_2_graph_parameters(instance_adjlist, kv, n, cliques, &discretizationEdges_2, &pruneEdges_2);
	dealloc_3Darray_lf(instance_adjlist, kv, n);
	free(kv);
	
	free(instanceFile);
	free(initialStructureFile);
	free(cliquesFile);
	
	if(strcmp(method, "iabp") == 0){
		print("****************** iABP ******************");
		iABP(n, discretizationEdges_2, tauSign, givenTau, givenTauDeviation, pruneEdges_2, sampleSize, angularResolution, tolerance, timeLimit, &protein, numOfSols, &solutionMetrics, dc_vec, m, is_X0_knownQM, X0);
	}
	else if(strcmp(method, "ibp") == 0){
		print("****************** iBP ******************");
		iBP(n, discretizationEdges_2, pruneEdges_2, sampleSize, tolerance, timeLimit, &protein, numOfSols, &solutionMetrics, dc_vec, m, is_X0_knownQM, X0);
	}
	else{
		printf("Error: the chosen method is neither iABP nor iBP\n");
		return 0;
	}
	
	strcpy(aux_str, outputFolder);
	FILE *fp = fopen(strcat(aux_str, "/results.txt"), "w");
	if(!fp){
		printf("Error: file reading error.\n");
		exit(1);
	}
	
	fprintf(fp, "CPU time = %.8lf\n", solutionMetrics.cpu_time_used);
	fprintf(fp, "Last embedded vertex = %d/%d\n", solutionMetrics.lev, n);
	fprintf(fp, "Number of placed vertices = %.0lf\n", solutionMetrics.npp);
	fprintf(fp, "Number of found solutions = %ld\n", solutionMetrics.nos);
	
	fclose(fp);
	
	if(solutionMetrics.nos > 0){
		strcpy(aux_str, outputFolder);
		save_protein_PDBformat(method, protein, n, structure_id, "1", structure_chain, strcat(strcat(strcat(aux_str, "/"), structure_id), ".pdb"));
		
		mde_vec = stack_2_vec_lf(solutionMetrics.MDE);
		strcpy(aux_str, outputFolder);
		vec_lf_2_file(mde_vec, solutionMetrics.nos, strcat(aux_str, "/MDE.dat"));
		free(mde_vec);
		
		lde_vec = stack_2_vec_lf(solutionMetrics.LDE);
		strcpy(aux_str, outputFolder);
		vec_lf_2_file(lde_vec, solutionMetrics.nos, strcat(aux_str, "/LDE.dat"));
		free(lde_vec);
		
		rmsd_vec = stack_2_vec_lf(solutionMetrics.RMSD);
		strcpy(aux_str, outputFolder);
		vec_lf_2_file(rmsd_vec, solutionMetrics.nos, strcat(aux_str, "/RMSD.dat"));
		free(rmsd_vec);
	}
	else{
		stack_delete(solutionMetrics.MDE);
		stack_delete(solutionMetrics.LDE);
		stack_delete(solutionMetrics.RMSD);
	}
	
	free(method);
	free(dc_vec);
	free(protein);
	dealloc_mat_lf(X0, n);
	dealloc_mat_d(cliques, n);
	free(tauSign);
	free(givenTau);
	free(givenTauDeviation);
	dealloc_fix3Darray_lf(discretizationEdges_2, n, 3);
	for(k = 0; k < n; k++){
		if(pruneEdges_2[k].cardUkP > 0)
			dealloc_mat_lf(pruneEdges_2[k].precise, pruneEdges_2[k].cardUkP);
		if(pruneEdges_2[k].cardUkI > 0)
			dealloc_mat_lf(pruneEdges_2[k].interval, pruneEdges_2[k].cardUkI);
	}
	free(pruneEdges_2);
	
	return 0;
}