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
Matriz sm (Matriz* m1, Matriz* m2){ return pm(m1, m2, "mais");}
Matriz sub(Matriz* m1, Matriz* m2){ return pm(m1, m2, "menos");}

Matriz matriz(double* mat, int linhas, int colunas) {
//Construtor do pseudo-objeto Matriz
	Matriz m;
	m.matriz = mat;
	m.colunas = colunas;
	m.linhas = linhas;
	m.ret = ret;
	m.print = print;
	return m;
}


int main(){
	
	double m1[2][2] = {{1, 2}, {3,4}};
	double m2[2][2] = {{5, 6}, {7,8}};
	Matriz m11 = matriz(&m1[0][0], 2, 2);
	Matriz m22 = matriz(&m2[0][0], 2, 2);
	Matriz res1 = sm(&m11, &m22);
	Matriz res2 = sub(&m11, &m22);
	res1.print(&res1);
	res2.print(&res2);
}

double ret(Matriz* self, int l, int c){
	return self->matriz[l*self->colunas + c];
}

void print(Matriz *self){
	const int l = self->linhas;
	const int c = self->colunas;
	double* mat = self->matriz;
	printf("[");
	
	for(int i = 0; i < l; i++){
		for(int j = 0; j < c; j++){
			if(j == c-1) printf("%f", self->ret(self, i, j));
			else printf("%f \t", self->ret(self, i, j));
		}
		if(i == l-1) printf("]\n");
		else printf("\n");
	}	
	printf("\n");
}

Matriz transposta(Matriz* self){
	const int l = self->linhas;
	const int c = self->colunas;
	double* nova_mat = (double *)malloc(l*c*sizeof(double));
	for(int i = 0; i < c; i++){
		for(int j = 0; j < l; j++){
			nova_mat[i*l + j] = self->ret(self, j, i);
		}
	}
	return matriz(nova_mat, c, l);
}

Matriz pm (Matriz* m1, Matriz* m2, char* op){
	if(m1->linhas == m2->linhas && m1->colunas==m2->colunas){
		const int l = m1->linhas;
		const int c = m1->colunas;
		double* r = (double*) malloc(l*c*sizeof(double));
		for(int i = 0; i < l; i++){
			for(int j = 0; j < c; j++){
				if(strcmp(op, "mais")==0){
					r[i*c + j] = m1->ret(m1, i, j) + m2->ret(m2, i, j);}
				else if (strcmp(op, "menos")==0){
					r[i*c + j] = m1->ret(m1, i, j) - m2->ret(m2, i, j);}
				}	
			}
		return matriz(r, m1->linhas, m1->colunas);
	}
	else {
		//Erro: Matrizes de tamanhos diferentes
		//Retorna uma matriz |0
		double m[1][1] = {0};
		return matriz(&m[0][0], 1, 1);
	}				       
}
