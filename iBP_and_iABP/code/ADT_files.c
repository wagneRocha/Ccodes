#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ADT_vectors.h"
#include "ADT_matrices.h"
#include "ADT_files.h"
#include "ADT_strings.h"
#define PI 3.14159265359

#define MAX_LINE_LENGTH 512
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
void vec_lf_2_file(double *v, int m, char filename[]){
	
	int i;
	
	FILE *fp = fopen(filename, "w");
	if(fp == NULL){
		fprintf(stderr, "File opening error for writing.\n");
		exit(1);
	}

	for(i = 0; i < m; i++)
		fprintf(fp, "%.8lf\n", v[i]);
		
	fclose(fp);
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
		replace_multiple_spaces(line_i);
		vecspaces = find_all_positions(line_i, sepChar);
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
/* *********************************************************************************** */
void read_input_file(char *fname, char **structure_id, char **structure_chain, char **method, char **fname_dc, char **fname_T0, char **fname_X0, double *timeLimit, double *tolerance, double *angularResolution, int *numOfSols, int *sampleSize){
	
	FILE *fp = fopen(fname, "r");
	if(fp == NULL){
		printf("Error: file reading error.\n");
		exit(1);
	}

	int numLines = num_lines_file(fname);
	
	if(numLines == 11){
		char line_i[MAX_LINE_LENGTH];
		char *content[numLines];
		int i = 0;
		
		for(i = 0; i < numLines; i++){
			if (fgets(line_i, MAX_LINE_LENGTH, fp) == NULL){
				fprintf(stderr, "Error: file has fewer lines than expected.\n");
				break;
			}
			char *colon_position = strchr(line_i, ':');
			if(colon_position){
				char *value = colon_position + 1;
				while(*value == ' ')
					value++;

				content[i] = custom_strdup(value);
				if(!content[i]){
					perror("Error: memory allocation failed");
					fclose(fp);
					exit(1);
				}
				
				content[i][strcspn(content[i], "\n")] = '\0';
			}
		}
		
		fclose(fp);
		
		*structure_id = malloc(strlen(content[0]) + 1);
		strcpy(*structure_id, content[0]);
		free(content[0]);
		
		*structure_chain = malloc(strlen(content[1]) + 1);
		strcpy(*structure_chain, content[1]);
		free(content[1]);
		
		*method = malloc(strlen(content[2]) + 1);
		strcpy(*method, content[2]);
		free(content[2]);
		
		*fname_dc = malloc(strlen(content[3]) + 1);
		strcpy(*fname_dc, content[3]);
		free(content[3]);
		
		*fname_T0 = malloc(strlen(content[4]) + 1);
		strcpy(*fname_T0, content[4]);
		free(content[4]);
		
		*fname_X0 = malloc(strlen(content[5]) + 1);
		strcpy(*fname_X0, content[5]);
		free(content[5]);
		
		double days, hours, minutes, seconds;
		int status = sscanf(content[6], "%lf-%lf:%lf:%lf", &days, &hours, &minutes, &seconds);
		if(status != 4)
			printf("Error: the time limit format is incorrect, the correct form is 'days-hours:minutes:seconds'\n");
		else
			*timeLimit = 86400.0*days + 3600.0*hours + 60.0*minutes + seconds;
		free(content[6]);
		
		*tolerance = strtod(content[7], NULL);
		free(content[7]);
		
		*angularResolution = strtod(content[8], NULL)*PI/180;
		free(content[8]);
		
		*numOfSols = (int) strtod(content[9], NULL);
		free(content[9]);
		
		*sampleSize = (int) strtod(content[10], NULL);
		free(content[10]);
	}
	else
		printf("Error: the number of lines in the input file is incorrect; it must contain 11 lines\n");
}
/* *********************************************************************************** */
void save_protein_PDBformat(char method[], proteinstructure *protein, int n, char structure_id[], char model[], char structure_chain[], char fname[]){
	
	FILE *fp = fopen(fname, "w");
	if(!fp){
		printf("Error: file reading error.\n");
		exit(1);
	}

	time_t t = time(NULL);
	struct tm tm = *localtime(&t);
	char date[11];
	snprintf(date, sizeof(date), "%04d-%02d-%02d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday);

	// Write HEADER
	fprintf(fp, "HEADER                                                %s   %s         \n", date, structure_id);

	// Write TITLE
	if(strcmp(method, "PDB") != 0)
		fprintf(fp, "TITLE     %-4s solution                                                          \n", method);
	else
		fprintf(fp, "TITLE     PDB structure                                                         \n");
	
	// Write REMARKS
	fprintf(fp, "REMARK   2                                                                      \n");
	fprintf(fp, "REMARK   2 RESOLUTION.    2.00 ANGSTROMS.                                       \n");
	fprintf(fp, "MODEL   %6s                                                                  \n", model);

	char atomType[7] = "";
	int serial;
	char name[5] = "";
	char name_aux[5] = "";
	const char altLoc = ' ';
	char resName[4] = "";
	char chainID[2] = "";
	strcpy(chainID, structure_chain);
	int resSeq;
	const char iCode = ' ';
	double x;
	double y;
	double z;
	double occupancy = 1.0;
	double tempFactor = 1.0;
	char element[2];
	element[1] = '\0';
	const char *charge = "  ";
	
	// Write ATOM entries
	strcpy(atomType, "ATOM  ");
	int k;
	for(k = 0; k < n; k++){
		serial = protein[k].v_k;
		strcpy(name_aux, protein[k].atom_name);
		strcpy(resName, protein[k].res_name);
		resSeq = protein[k].res_seq;
		x = protein[k].x;
		y = protein[k].y;
		z = protein[k].z;
		element[0] = protein[k].atom_name[0];
		
		if(strlen(name_aux) == 4)
			strcpy(name, name_aux);
		else if(strlen(name_aux) == 3){
			strcpy(name, " ");
			strcat(name, name_aux);
		}
		else if(strlen(name_aux) == 2){
			strcpy(name, " ");
			strcat(strcat(name, name_aux), " ");
		}
		else if(strlen(name_aux) == 1){
			strcpy(name, " ");
			strcat(strcat(name, name_aux), "  ");
		}
		
		fprintf(fp, "%6s%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
			atomType, serial, name, altLoc, resName, *chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, element, charge);
	}

	// Write TER entry
	strcpy(atomType, "TER   ");
	serial = n + 1;
	strcpy(name, "    ");
	strcpy(resName, protein[n-1].res_name);
	resSeq = protein[n-1].res_seq;

	fprintf(fp, "%6s%5d %4s%c%3s %c%4d%c\n",
		atomType, serial, name, altLoc, resName, *chainID, resSeq, iCode);
	
	fclose(fp);
}
