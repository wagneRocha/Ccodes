#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ADT_vectors.h"
#include "ADT_stack.h"

/* *********************************************************************************** */
/* ---------------------------------- TAD stack double ------------------------------- */
/* *********************************************************************************** */
node *node_create(double value){

	node *no = (node*)malloc(sizeof(node));
	no -> value = value;
	no -> next  = NULL;

	return no;
}
/* *********************************************************************************** */
stack *stack_create(){
	
	stack *s = (stack*)malloc(sizeof(stack));
	s -> head = node_create(404.0);
	s -> size = 0;

	return s;
}
/* *********************************************************************************** */
void stack_push(stack *s, double val){

	node *new_no = node_create(val);
	
	new_no -> next = (s -> head) -> next;
	(s -> head) -> next = new_no;			
	s -> size = s -> size + 1;	
}
/* *********************************************************************************** */
void stack_pop(stack *s){

	node *aux;
	
	if(s -> size == 0)
		printf("Empty Stack!\n");
	else{
		aux = (s -> head) -> next;
		(s -> head) -> next = aux -> next;
		free(aux);
		s -> size = s -> size - 1;
	}
}
/* *********************************************************************************** */
double stack_top(stack *s){

	if(s -> size == 0){
		printf("Empty Stack!\n");
		return 0.0; 
	}
	return (s -> head) -> next -> value;
}
/* *********************************************************************************** */
int stack_size(stack *s){

	return (s -> size);
}
/* *********************************************************************************** */
void stack_delete(stack *s){

	node *aux, *trash;

	if(stack_size(s) > 0){
		aux = (s -> head) -> next;	
		while(aux){
			trash = aux;
			aux = aux -> next;
			free(trash);
		}
		free(aux);
		s -> size = 0;
	}
	else{
		free(s -> head);
		free(s);
	}
}
/* *********************************************************************************** */
double *stack_2_vec_lf(stack *s){

	int i = 0;
	int n = stack_size(s);
	double *v = zeros_vec_lf(n);
	
	while(s -> size != 0){
		v[i++] = stack_top(s);
		stack_pop(s);
	}
	
	stack_delete(s);
	
	return v;
}