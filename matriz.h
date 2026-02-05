#ifndef MATRIZ_H
#define MATRIZ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define BUFF 100
#define MATRIZ_NULA() matriz((double[BUFF]){0}, 0, 0)


typedef struct Matriz Matriz;

struct Matriz {
	int colunas;
	int linhas;
	double* matriz;
	double (*ret)(Matriz*, int, int);
	double (*transposta)(Matriz*);
	void (*print)(Matriz*);
};
void cop_mat(Matriz* m1, Matriz* m2);
Matriz matriz(double* mat, int linhas, int colunas);
Matriz matriz_nula();
void mat_pot(Matriz* mat, Matriz* res, int pot);
double ret(Matriz* self, int l, int c);
void transposta(Matriz* self, Matriz* res);
void print(Matriz *self);
void pm(Matriz* m1, Matriz* m2, Matriz* res, char* op);
void sm (Matriz* m1, Matriz* m2, Matriz* res);
void sub(Matriz* m1, Matriz* m2, Matriz* res);
void ident(Matriz* m, int t);
void ret_lin(Matriz* m, Matriz* dest, int l);
void mult(Matriz* m, Matriz* destino, double v);
void subst(Matriz* m_linha, Matriz* dest , int indice);
void extend(Matriz* m, Matriz* n, Matriz* dest);
void mat_mult(Matriz* a, Matriz* b, Matriz* res);
double int_prod(Matriz* l1, Matriz* l2);
void inv(Matriz *origin, Matriz* dest);
Matriz zeros(int l, int c); 

#endif
