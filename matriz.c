#include <stdio.h>

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
	double matx[3][3] = {
		{1,2,3},{4,5,6},{7,8,9}
	};
	Matriz matj = matriz(&matx[0][0], 3, 3);
	printf("Tamanho da linha: %i\n", matj.linhas);
	printf("Tamanho da coluna: %i\n", matj.colunas);
	printf("[][]: %f\n", matj.ret(&matj, 1,1));
	matj.print(&matj);
	Matriz nova = transposta(&matj);
	nova.print(&nova);
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

