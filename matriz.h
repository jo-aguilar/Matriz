#ifndef MATRIZ_H
#define MATRIZ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct Matriz Matriz;

struct Matriz {
	int colunas;
	int linhas;
	double* matriz;
	double (*ret)(Matriz*, int, int);
	double (*transposta)(Matriz*);
	void (*print)(Matriz*);
};

Matriz matriz(double* mat, int linhas, int colunas);
double ret(Matriz* self, int l, int c);
Matriz transposta(Matriz *self);
void print(Matriz *self);
Matriz pm(Matriz*, Matriz*, char*);
Matriz sm (Matriz*, Matriz*);
Matriz sub(Matriz*, Matriz*);
Matriz ident(int);
Matriz ret_lin(Matriz*, int);
Matriz mult(Matriz*, double);
void subst(Matriz*, Matriz* , int);
Matriz extend(Matriz*, Matriz*);
Matriz mat_mult(Matriz*, Matriz*);
double int_prod(Matriz*, Matriz*);
Matriz inv(Matriz *);
Matriz zeros(int, int); 

#endif
