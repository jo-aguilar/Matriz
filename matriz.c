#include "matriz.h"

void cop_mat(Matriz* m1, Matriz* m2){
	m2->linhas = m1->linhas;
	m2->colunas = m1->colunas;
	int contador = 0;
	for(int linhas = 0; linhas < m2->linhas; linhas++){
		for(int colunas = 0; colunas < m2->colunas; colunas++){
			m2->matriz[contador++] = m1->ret(m1, linhas, colunas);
		}
	}
}

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

void mat_pot(Matriz* mat, Matriz* res, int pot) {
//Entrega uma potência qualquer de uma matriz fornecida pelo usuário
	if(pot == 0){
		Matriz I = MATRIZ_NULA(); ident(&I, mat->colunas);
		res->matriz = I.matriz;
		res->colunas = I.colunas;
		res->linhas = I.linhas;
		return;
	}

	if(pot == 1){
		res->linhas = mat->linhas;
		res->colunas = mat->colunas;
		res->matriz = mat->matriz;
		return;
	}

	Matriz r = MATRIZ_NULA();
	for(int i = 0; i < pot-1; i++){
		if(i == 0) mat_mult(mat, mat, &r);
		else mat_mult(mat, &r, &r); 
	}
	res->linhas  = r.linhas;
	res->colunas = r.colunas;
	res->matriz  = r.matriz;
}

Matriz zeros(int l, int c){
//Retorna uma matriz de zeros do tamanho especificado pelo
//usuário.
	double* mat = malloc(l*c*sizeof(double));
	memset(mat, 0, l*c*sizeof(double));
	return matriz(mat, l, c);
} 

void inv(Matriz *origin, Matriz* dest){
//Utiliza Gauss-Jordan para a obtenção da matriz inversa de uma 
//matriz quadrada fornecida pelo usuário
	int linhas  = origin->linhas;
	int colunas = origin->colunas; 
	Matriz id  = MATRIZ_NULA(); ident(&id, linhas);
	Matriz est = MATRIZ_NULA(); extend(origin, &id, &est); 
	Matriz l1 = MATRIZ_NULA();
	Matriz l2 = MATRIZ_NULA();
	Matriz novo_l1 = MATRIZ_NULA();
	Matriz novo_l2 = MATRIZ_NULA();
	Matriz temp = MATRIZ_NULA();

	for(int i = 0; i < linhas; i++){
	//Opera Gauss-Jordan sobre uma matriz estendida [A|B], transformando
	//A em uma matriz identidade e B na inversa
		int n = i;
		double fator = 1/ret(&est, n, n);
		ret_lin(&est, &l1, n);
		mult(&l1, &temp, fator);
		subst(&temp, &est, i+1);
		
		for(int j = 0; j < linhas; j++){
			if(i!=j){
				ret_lin(&est, &l2, j);
				double novo_fator = ret(&est, j, i);
				ret_lin(&est, &l1, i);
				mult(&l1, &novo_l1, novo_fator);
				sub(&l2, &novo_l1, &novo_l2);
				subst(&novo_l2, &est, j+1);
			}
		}
	}
	int contador = 0;
	Matriz Lin_temp = MATRIZ_NULA();
	for(int i = 0; i < linhas; i++){
	//Extrai a matriz B como inversa e passa seu resultado para a
	//matriz de destino
		ret_lin(&est, &Lin_temp, i);
		for(int j = colunas; j < 2*colunas; j++)
			dest->matriz[contador++] = Lin_temp.ret(&Lin_temp, 0, j);
	}
	dest->linhas = origin->linhas;
	dest->colunas = origin->colunas;
}

double int_prod(Matriz* l1, Matriz *l2){
//A partir de duas matrizes-linha, retorna o produto
//interno de ambas
	double prod = 0;
	for(int i=0; i < l1->colunas; i++)
		prod = prod + (l1->ret(l1, 0, i))*(l2->ret(l2, 0,i));
	return prod;
}


void mat_mult(Matriz* a, Matriz* b, Matriz* res){
//Multiplica duas matrizes fornecidas pelo usuário, contanto que 
//ambas possuam as mesmas dimensões linhas/colunas, colunas/linhas
	int linhas = a->linhas;
	int colunas = b->colunas;
	Matriz b_transp = MATRIZ_NULA();
	transposta(b, &b_transp);
	int indice = 0;
	res->linhas = linhas;
	res->colunas = colunas;

	for(int i = 0; i < linhas; i++){
		Matriz l1 = MATRIZ_NULA();
		ret_lin(a, &l1, i);
		for(int j = 0; j < colunas; j++){
			Matriz l2 = MATRIZ_NULA();
			ret_lin(&b_transp, &l2, j);
			res->matriz[indice++] = int_prod(&l1, &l2);
		}
	}
}


void extend(Matriz* m, Matriz*n, Matriz* dest){
//Entrega uma matriz estendida [A|B] a partir de duas matrizes 
//A e B fornecidas pelo usuário (com a mesma quantidade de linhas)
	if(m->linhas!=n->linhas){
		printf("ERRO!!! Dimensões de matrizes não conformantes.\nRetornando...\n");
		return;
	}

	int colunas = m->colunas;
	int linhas  = m->linhas;
	int contador = 0;
	for(int i = 0; i < linhas; i++){
		for(int j = 0; j < m->colunas; j++)
			dest->matriz[contador++] = m->ret(m, i, j);
		for(int j = 0; j < n->colunas; j++)
			dest->matriz[contador++] = n->ret(n, i, j);
	}
	dest->linhas = linhas;
	dest->colunas = m->colunas + n->colunas;
}




void mult(Matriz* m, Matriz* destino, double v){
//Devolve uma linha qualquer completa de uma matriz multiplicada 
//por um fator v fornecido pelo usuário
	destino->linhas  = m->linhas;
	destino->colunas = m->colunas;
	for(int i = 0; i < m->colunas*m->linhas; i++)
		destino->matriz[i] = v*m->matriz[i];
}


void subst(Matriz* m_linha, Matriz* dest, int indice){
//Substitui uma linha arbitrária de uma matriz por uma nova linha
	if(m_linha->colunas != dest->colunas){
		printf("ERRO!!! Colunas de linha e de destino de tamanhos diferentes"
			"\nRetornando em subst(...)...\n");
		return;
	}
	int colunas = m_linha->colunas;
	//printf("%f\n", retorno);
	for(int i = 0; i < colunas; i++) 
		dest->matriz[(indice-1)*colunas + i] = 
				m_linha->ret(m_linha, 0, i);
}

void ret_lin(Matriz* m, Matriz* dest, int l){
//Retorna uma linha qualquer completa de dentro da matriz fornecida
	int colunas = m->colunas;
	dest->linhas = 1;
	dest->colunas = colunas;
	for(int i = 0; i < colunas; i++)
		dest->matriz[i] = m->ret(m, l, i);
}


void ident(Matriz* m,  int t){
//Retorna uma matriz identidade de dimensão txt
	m->linhas = t;
	m->colunas = t;
	memset(m->matriz, 0, sizeof(m->matriz));
	for(int i = 0; i < t; i++) { m->matriz[(t+1)*i] = 1.0;}
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
			if(j == c-1) printf("%3lg", self->ret(self, i, j));
			else printf("%8.3lg \t", self->ret(self, i, j));
		}
		if(i == l-1) printf("]\n");
		else printf("\n");
	}	
	printf("\n");
}


void transposta(Matriz* self, Matriz* res){
//Retorna uma matriz transposta de uma matriz retangular
//de qualquer dimensão
	const int l = self->linhas;
	const int c = self->colunas;
	res->linhas  = c;
	res->colunas = l;
	for(int i = 0; i < c; i++){
		for(int j = 0; j < l; j++){
			res->matriz[i*l + j] = self->ret(self, j, i);
		}
	}
}

void pm (Matriz* m1, Matriz* m2, Matriz* res, char* op){
//Função superior para soma e subtração de matrizes
	res->linhas = m2->linhas;
	res->colunas = m2->colunas;
	if(m1->linhas == m2->linhas && m1->colunas==m2->colunas){
		const int l = m1->linhas;
		const int c = m1->colunas;
		for(int i = 0; i < l; i++){
			for(int j = 0; j < c; j++){
				if(strcmp(op, "mais")==0){
					res->matriz[i*c + j] =
						m1->ret(m1, i, j) + m2->ret(m2, i, j);}
				else if (strcmp(op, "menos")==0){
					res->matriz[i*c + j] =
						m1->ret(m1, i, j) - m2->ret(m2, i, j);}
				}	
			}
	}
	else {
		//Erro: Matrizes de tamanhos diferentes
		//Retorna uma matriz |0
		return;
	}				       
}

void sm (Matriz* m1, Matriz* m2, Matriz* res){
//Retorna a soma entre duas matrizes
	pm(m1, m2, res, "mais");
}

void sub(Matriz* m1, Matriz* m2, Matriz* res){
//Retorna a subtração entre duas matrizes
	pm(m1, m2, res, "menos");
}
