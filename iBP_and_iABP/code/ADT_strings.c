#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ADT_vectors.h"
#include "ADT_matrices.h"
#include "ADT_files.h"
#include "ADT_strings.h"

#define MAX_LINE_LENGTH 128
/* *********************************************************************************** */
/* ---------------------------------- TAD strings ------------------------------------ */
/* *********************************************************************************** */
void replace_multiple_spaces(char *str){

	char result[1000];  // Adjust the size based on your needs
	char *token;
	char delimiter[] = " ";

	// Initialize the result string
	result[0] = '\0';

	// Tokenize the original string using space as a delimiter
	token = strtok(str, delimiter);

	// Construct the result string with a single space between each word
	while (token != NULL){
		strcat(result, token);
		strcat(result, " ");  // Add a single space after each word
		token = strtok(NULL, delimiter);
	}

	// Remove the trailing space, if any
	if (result[strlen(result) - 1] == ' '){
		result[strlen(result) - 1] = '\0';
	}

	// Copy the result back to the original string
	strcpy(str, result);
}
/* *********************************************************************************** */
int *find_all_positions(char *str, char sepChar){
	
	int i, j = 0;
	int *vecSep = alloc_vec_d(9);
		
	for(i = 0; str[i] != '\0'; i++)
		if (str[i] == sepChar)
			vecSep[j++] = i;
	
	return vecSep;
}
/* *********************************************************************************** */
char *custom_strdup(char *str){
	
	size_t len = strlen(str) + 1;
	char *dup = malloc(len);
	
	if(dup)
		memcpy(dup, str, len);
	
	return dup;
}
/* *********************************************************************************** */
void print(char *str){
	
	printf("%s\n", str);
}
