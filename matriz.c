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
Matriz ident(int);

Matriz matriz(double* mat, int linhas, int  colunas) {
//Construtor do pseudo-objeto Matriz
	Matriz m;
	m.matriz = mat;
	m.colunas = colunas;
	m.linhas = linhas;
	m.ret = ret;
	m.print = print;
	return m;
}

Matriz ret_lin(Matriz*, int);
Matriz mult(Matriz*, double);
void subst(Matriz*, Matriz* , int);
Matriz extend(Matriz* m, Matriz*n)
{
	int colunas = m->colunas;
	int linhas  = m->linhas;
	double *original = m->matriz;
	double *identidade = n->matriz;
	double *retorno = malloc(2*linhas*colunas*sizeof(double));
	printf("[DEBUG] Quant. colunas: %d\n",colunas);
	printf("[DEBUG] Quant. linhas:  %d\n", linhas);	
	for(int i = 0; i < linhas; i++){	
		double* nova_linha = malloc(2*colunas*sizeof(double));
		memcpy(nova_linha,           original   + i*colunas, colunas*sizeof(double));
		memcpy(nova_linha + colunas, identidade + i*colunas, colunas*sizeof(double));
		memcpy(retorno + i*2*colunas, nova_linha, 2*colunas*sizeof(double));
		free(nova_linha);
	}
	return matriz(retorno, linhas, 2*colunas);

}


int main(){	
	//double m1[3][3] = {{1, 2, 3}, {4, 5, 6},{7, 8, 9}};
	//Matriz m = matriz(&m1[0][0], 3,3);
	Matriz m1 = matriz(((double[]){1, 2, 3, 4, 5, 6, 7, 8, 9}), 3, 3);
	Matriz m2 = matriz(((double[]){10, 11, 12, 13, 14, 15, 16, 17, 18}), 3, 3);	

	//m.print(&m);
	//i.print(&i);
	Matriz m3 = extend(&m1, &m2);
	m3.print(&m3);
	
}

Matriz ret_lin(Matriz* m, int l){
//Retorna uma linha qualquer completa de dentro da matriz fornecida
	int colunas = m->colunas;
	double* linha = malloc(colunas*sizeof(double));
	linha = memcpy(linha, &(m->matriz[l*colunas]) , colunas*sizeof(double));
	return matriz(linha, 1, colunas); 
}

Matriz mult(Matriz* m, double v){
//Devolve uma linha qualquer completa de uma matriz multiplicada 
//por um fator v fornecido pelo usuário
	double* mat = m->matriz;
	for(int i = 0; i < m->colunas; i++)
		mat[i] = v*mat[i]; 
	return matriz(mat, m->linhas, m->colunas);
}

void subst(Matriz* m, Matriz* l, int indice){
//Substitui uma linha arbitrária de uma matriz por uma nova linha
	int colunas = l->colunas;
	memcpy(&(m->matriz[indice*colunas]), l->matriz, indice*colunas*sizeof(double));
}

Matriz ident(int t){
//Retorna uma matriz identidade de dimensão txt
	double* arr = malloc(t*t*sizeof(double));
	memset(arr, 0, sizeof(arr));
	for(int i = 0; i < t*t; i++) { arr[(t+1)*i] = 1;}
	return matriz(arr, t, t);
}

double ret(Matriz* self, int l, int c){
//Retorna um elemento definido pelo índice (a,b) de uma matriz
//de qualquer dimensão
	return self->matriz[l*self->colunas + c];
}

void print(Matriz *self){
//Apresenta uma matriz formatada para a visualização do usuário
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
//Retorna uma matriz transposta de uma matriz retangular
//de qualquer dimensão
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
//Função superior para soma e subtração de matrizes
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
