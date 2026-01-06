#include "matriz.h"

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

Matriz mat_pot(Matriz* mat, double pot) {
	double* m1 = mat->matriz;
	double* m = malloc(mat->linhas*mat->colunas*sizeof(double));
	memcpy(m, m1, mat->linhas*mat->colunas*sizeof(double));
	for(int i = 0; i < mat->colunas*mat->linhas; i++) {
		m[i] = pow(m[i], pot);
	}
	return matriz(m, mat->linhas, mat->colunas);
}


Matriz zeros(int l, int c){
//Retorna uma matriz de zeros do tamanho especificado pelo
//usuário.
	double* mat = malloc(l*c*sizeof(double));
	memset(mat, 0, l*c*sizeof(double));
	return matriz(mat, l, c);
}


Matriz inv(Matriz *origin){
//Utiliza Gauss-Jordan para a obtenção da matriz inversa de uma 
//matriz quadrada fornecida pelo usuário
	int linhas  = origin->linhas;
	int colunas = origin->colunas; 
	Matriz id = ident(linhas);
	Matriz est = extend(origin, &id);
	
	for(int i = 0; i < linhas; i++){
		int n = i;
		double fator = 1/ret(&est, n, n);
		Matriz l1 = ret_lin(&est, n);
		l1 = mult(&l1, fator);
		subst(&est, &l1, i);
		
		for(int j = 0; j < linhas; j++){
			if(i!=j){
				Matriz l2 = ret_lin(&est, j);
				double novo_fator = ret(&est, j, i);
				l1 = ret_lin(&est, i);
				Matriz novo_l1 = mult(&l1, novo_fator);
				Matriz novo_l2 = sub(&l2, &novo_l1);
				subst(&est, &novo_l2, j);
			}
		}
	}
	double* retorno = malloc(linhas*colunas*sizeof(double));
	for(int i = 0; i < linhas; i++){
		Matriz linha = ret_lin(&est, i);
		double* l = linha.matriz;
		memcpy(&retorno[i*colunas], &l[colunas] , colunas*sizeof(double));
	}
	return matriz(retorno, linhas, colunas);
}

double int_prod(Matriz* l1, Matriz *l2){
//A partir de duas matrizes-linha, retorna o produto
//interno de ambas
	double prod = 0;
	for(int i=0; i < l1->colunas; i++)
		prod = prod + (l1->ret(l1, 0, i))*(l2->ret(l2, 0,i));
	return prod;
}

Matriz mat_mult(Matriz* a, Matriz* b){
//Multiplica duas matrizes fornecidas pelo usuário
	int linhas = a->linhas;
	int colunas = b->colunas;
	double* res = malloc(linhas*colunas*sizeof(double));
	Matriz b_transp = transposta(b);
	
	int indice = 0;
	for(int i = 0; i < linhas; i++){
		for(int j = 0; j < colunas; j++){
			Matriz l1 = ret_lin(a, i);
			Matriz l2 = ret_lin(&b_transp, j);
			res[indice++] = int_prod(&l1, &l2);
		}
	}
	return matriz(res, linhas, colunas);
}

Matriz extend(Matriz* m, Matriz*n){
//Entrega uma matriz estendida [A|B] a partir de duas matrizes 
//A e B fornecidas pelo usuário (com a mesma quantidade de linhas)
	int colunas = m->colunas;
	int linhas  = m->linhas;
	int n_col   = n->colunas;
	int n_lin   = n->linhas;
	double *original = m->matriz;
	double *identidade = n->matriz;
	double *retorno = malloc(linhas*(colunas+n_col)*sizeof(double));
	for(int i = 0; i < linhas; i++){	
		double* nova_linha = malloc((colunas + n_col)*sizeof(double));
		memcpy(nova_linha,           original   + i*colunas, colunas*sizeof(double));
		memcpy(nova_linha + colunas, identidade + i*n_col, n_col*sizeof(double));
		memcpy(retorno + i*(n_col + colunas), nova_linha, (n_col+colunas)*sizeof(double));
		free(nova_linha);
	}
	return matriz(retorno, linhas, n_col + colunas);
}

Matriz mult(Matriz* m, double v){
//Devolve uma linha qualquer completa de uma matriz multiplicada 
//por um fator v fornecido pelo usuário
	double* mat = m->matriz;
	for(int i = 0; i < m->colunas*m->linhas; i++)
		mat[i] = v*mat[i]; 
	return matriz(mat, m->linhas, m->colunas);
}

void subst(Matriz* m, Matriz* l, int indice){
//Substitui uma linha arbitrária de uma matriz por uma nova linha
	int colunas = l->colunas;
	memcpy(&(m->matriz[indice*colunas]), l->matriz, colunas*sizeof(double));
}

Matriz ret_lin(Matriz* m, int l){
//Retorna uma linha qualquer completa de dentro da matriz fornecida
	int colunas = m->colunas;
	double* linha = malloc(colunas*sizeof(double));
	linha = memcpy(linha, &(m->matriz[l*colunas]) , colunas*sizeof(double));
	return matriz(linha, 1, colunas); 
}

Matriz ident(int t){
//Retorna uma matriz identidade de dimensão txt
	double* arr = malloc(t*t*sizeof(double));
	memset(arr, 0, t*t*sizeof(arr));
	for(int i = 0; i < t; i++) { arr[(t+1)*i] = 1;}
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
			if(j == c-1) printf("%lg", self->ret(self, i, j));
			else printf("%8.5lg \t", self->ret(self, i, j));
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

Matriz sm (Matriz* m1, Matriz* m2){
//Retorna a soma entre duas matrizes
	return pm(m1, m2, "mais");
}

Matriz sub(Matriz* m1, Matriz* m2){
//Retorna a subtração entre duas matrizes
	return pm(m1, m2, "menos");
}
