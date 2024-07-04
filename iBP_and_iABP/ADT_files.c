#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ADT_vectors.h"
#include "ADT_matrices.h"
#include "ADT_files.h"

#define MAX_LINE_LENGTH 128
/* *********************************************************************************** */
/* ---------------------------------- TAD files -------------------------------------- */
/* *********************************************************************************** */
int num_lines_file(char filename[]){
	
	int m = 0;
	char c;

	FILE *fp = fopen(filename, "r");
	if(fp == NULL){
		printf("File reading error.\n");
		exit(1);
	}
	
	for (c = getc(fp); c != EOF; c = getc(fp))
			if (c == '\n')
				m++;
	fclose(fp);
	
	return m;
}
/* *********************************************************************************** */
int num_columns_file_lines(char filename[], char sep_char){
	
	int n = 0;
	char c;
	
	FILE *fp = fopen(filename, "r");
	if(fp == NULL){
		printf("File reading error.\n");
		exit(1);
	}
	
	for (c = getc(fp); c != '\n'; c = getc(fp))
			if (c == sep_char)
				n++;
	fclose(fp);
	
	return ++n;
}
/* *********************************************************************************** */
void mat_lf_2_file(double **M, int m, int n, char filename[], char sep_char){
	
	int i, j;
	
	FILE *fp = fopen(filename, "w");
	if(fp == NULL){
		printf("File reading error.\n");
		exit(1);
	}

	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++)
			if(j <= n-2)
				fprintf(fp, "%e%c", M[i][j], sep_char);
			else
				fprintf(fp, "%e\n", M[i][j]);
	}
	fclose(fp);
}
/* *********************************************************************************** */
double **file_2_mat_lf(char filename[], char sep_char){
	
	int status;

	FILE *fp = fopen(filename, "r");
	if(fp == NULL){
		printf("File reading error.\n");
		exit(1);
	}
	
	int m = num_lines_file(filename);
	int n = num_columns_file_lines(filename, sep_char);
	
	int i, j;
	double **T = alloc_mat_lf(m, n);
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++){
			status = fscanf(fp,"%lf", &T[i][j]);
			if(status != 1){
				printf("Error reading file at line %d, column %d.\n", i+1, j+1);
				fclose(fp);
				dealloc_mat_lf(T, m);
				exit(1);
			}
		}
	fclose(fp);
	
	return T;
}
/* *********************************************************************************** */
int **file_2_mat_d(char filename[], char sep_char){

	int status;
	
	FILE *fp = fopen(filename, "r");
	if(fp == NULL){
		printf("File reading error.\n");
		exit(1);
	}
	int m = num_lines_file(filename);
	int n = num_columns_file_lines(filename, sep_char);

	int i, j;
	int **T = alloc_mat_d(m, n);
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++){
			status = fscanf(fp,"%d", &T[i][j]);
			if(status != 1){
				printf("Error reading file at line %d, column %d.\n", i+1, j+1);
				fclose(fp);
				dealloc_mat_d(T, m);
				exit(1);
			}
		}
	fclose(fp);
	
	return T;
}
/* *********************************************************************************** */
dcinputfile *dcInputFile_2_dcDataVector(char *fname){
	
	char line_i[MAX_LINE_LENGTH], maxsubstring[22];
	char sepChar = ' ';
	int m = num_lines_file(fname);
	
	int *vecspaces;
	
	dcinputfile *dcvector = alloc_vec_dcInputFile(m);
	
	int i, len;
		
	FILE *fp = fopen(fname, "r");
	if(fp == NULL){
		printf("File reading error.\n");
		exit(1);
	}
	
	for(i = 0; i < m; i++){
		if (fgets(line_i, MAX_LINE_LENGTH, fp) == NULL) {
			// If the file ends before reaching the desired line
			fprintf(stderr, "Error: File has fewer lines than expected.\n");
			break;
		}
		replaceMultipleSpaces(line_i);
		vecspaces = findAllPositions(line_i, sepChar);
		//printf("Line %d: %s", i, line_i);
		// vertex i ************************************************
		strcpy(maxsubstring, "                     ");
		len = vecspaces[0];
		strncpy(maxsubstring, line_i, len);
		maxsubstring[len] = '\0';
		dcvector[i].i  = atoi(maxsubstring);
		//printf("%d\n", dcvector[i].i);
		// vertex j ************************************************
		strcpy(maxsubstring, "                     ");
		len = vecspaces[1] - vecspaces[0]-1;
		strncpy(maxsubstring, line_i + vecspaces[0] + 1, len);
		maxsubstring[len] = '\0';
		dcvector[i].j  = atoi(maxsubstring);
		//printf("%d\n", dcvector[i].j);
		// residue vertex i ************************************************
		strcpy(maxsubstring, "                     ");
		len = vecspaces[2] - vecspaces[1]-1;
		strncpy(maxsubstring, line_i + vecspaces[1] + 1, len);
		maxsubstring[len] = '\0';
		dcvector[i].ri = atoi(maxsubstring);
		//printf("%d\n", dcvector[i].ri);
		// residue vertex j ************************************************
		strcpy(maxsubstring, "                     ");
		len = vecspaces[3] - vecspaces[2]-1;
		strncpy(maxsubstring, line_i + vecspaces[2] + 1, len);
		maxsubstring[len] = '\0';
		dcvector[i].rj = atoi(maxsubstring);
		//printf("%d\n", dcvector[i].rj);
		// lower bound ************************************************
		strcpy(maxsubstring, "                     ");
		len = vecspaces[4] - vecspaces[3]-1;
		strncpy(maxsubstring, line_i + vecspaces[3] + 1, len);
		maxsubstring[len] = '\0';
		dcvector[i].dl = strtod(maxsubstring, NULL);
		//printf("%.16lf\n", dcvector[i].dl);
		// upper bound ************************************************
		strcpy(maxsubstring, "                     ");
		len = vecspaces[5] - vecspaces[4]-1;
		strncpy(maxsubstring, line_i + vecspaces[4] + 1, len);
		maxsubstring[len] = '\0';
		dcvector[i].du = strtod(maxsubstring, NULL);
		//printf("%.16lf\n", dcvector[i].du);
		// atom i ************************************************
		strcpy(maxsubstring, "                     ");
		len = vecspaces[6] - vecspaces[5]-1;
		strncpy(maxsubstring, line_i + vecspaces[5] + 1, len);
		maxsubstring[len] = '\0';
		strcpy(dcvector[i].ai, maxsubstring);
		//printf("%s\n", dcvector[i].ai);
		// atom j ************************************************
		strcpy(maxsubstring, "                     ");
		len = vecspaces[7] - vecspaces[6]-1;
		strncpy(maxsubstring, line_i + vecspaces[6] + 1, len);
		maxsubstring[len] = '\0';
		strcpy(dcvector[i].aj, maxsubstring);
		//printf("%s\n", dcvector[i].aj);
		// res type atom i ************************************************
		strcpy(maxsubstring, "                     ");
		len = 3;
		strncpy(maxsubstring, line_i + vecspaces[7] + 1, len);
		maxsubstring[len] = '\0';
		strcpy(dcvector[i].rti, maxsubstring);
		//printf("%s\n", dcvector[i].rti);
		// res type atom j ************************************************
		strcpy(maxsubstring, "                     ");
		len = 3;
		strncpy(maxsubstring, line_i + vecspaces[8] + 1, len);
		maxsubstring[len] = '\0';
		strcpy(dcvector[i].rtj, maxsubstring);
		//printf("%s\n", dcvector[i].rtj);
		
		free(vecspaces);
	}
	
	fclose(fp);
			
	return dcvector;
}
