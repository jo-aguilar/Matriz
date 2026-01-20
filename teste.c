#include <stdio.h>
#include <math.h>
#include "dmpc.h"


void malha_fechada_teste(Matriz*, Matriz*, Matriz*, int, int, 
		  Matriz*, Matriz*, double, Matriz*, double);

int main(){	
	
	double a = 0.8; 
	double b = 0.1;
	Matriz A = matriz((double[]){a, 0, a, 1}, 2, 2);
	Matriz B = matriz((double[]){b, b}, 2, 1);
	Matriz C = matriz((double[]){0, 1}, 1, 2);
	Matriz A2 = matriz((double[]){a, a, b, b}, 2, 2);
	int Np = 10;
	int Nc = 4;
	Matriz Phi = MATRIZ_NULA(); func_Phi(&A, &B, &C, Np, Nc, &Phi);
	Matriz F = MATRIZ_NULA(); func_F(&A, &C, &F, Np);

	double r_w = 0.033620009; //encontrado empiricamente
	double r_ki = 1;

	Matriz x_ki = matriz((double[]){0.1, 0.2}, 2, 1);

	double* rs_barra = malloc(Np*sizeof(double));
	for(int i = 0; i < Np; i++) rs_barra[i] = 1.0;
	Matriz Rs_barra = matriz(rs_barra, Np, 1);
	Matriz M0 = MATRIZ_NULA(); ident(&M0, Nc);
	Matriz R_barra = MATRIZ_NULA(); mult(&M0, &R_barra, r_w);
	
	Matriz R = MATRIZ_NULA();
	//func_deltau(&Phi, &R_barra, &Rs_barra, &F, r_ki, &x_ki, &R);
	//R.print(&R);
	
	
	for(int i = 0; i < 10; i++){
		double du = malha_fechada(Nc, Np, &F, &Phi, &R_barra, 
			&Rs_barra, r_ki, r_w, &x_ki);
		Matriz Res = MATRIZ_NULA(); mat_mult(&C, &x_ki, &Res);
		printf("Valor de D_U em main(): %f\n", du);
	}
	
}


/*
void malha_fechada_teste(Matriz* A, Matriz* B, Matriz* C, int Nc,
		int Np, Matriz* F, Matriz* Phi, double r_ki, 
		Matriz* x_ki, double rw){
	double* rs_barra = malloc(Np*sizeof(double));
	for(int i = 0; i < Np; i++) rs_barra[i] = 1.0;
	Matriz Rs_barra = matriz(rs_barra, Np, 1);
	Matriz R_barra = ident(Nc);
	R_barra = mult(&R_barra, rw);

	Matriz Delta_U = func_deltau(Phi, &R_barra, &Rs_barra,
				     F, r_ki, x_ki);

	Matriz Um = matriz((double[]){1}, 1, 1);
	Matriz Zeros = zeros(1, Nc - 1);
	Matriz Primeiro_fator = extend(&Um, &Zeros);
	Matriz Phi_T = transposta(Phi);
	Matriz K = mat_mult(&Phi_T, Phi);
	K = sm(&K, &R_barra);
	K = inv(&K);
	K = mat_mult(&K, &Phi_T);
	K = mat_mult(&Primeiro_fator, &K);

	Matriz Ky   = mat_mult(&K, &Rs_barra);
	Matriz Kmpc = mat_mult(&K, F);
	Matriz M0 = mat_mult(B, &Kmpc);
	M0 = sub(A, &M0);
	M0 = mat_mult(&M0, x_ki);
	Matriz M1 = mat_mult(B, &Ky);
	M1 = mult(&M1, r_ki);
	
	//retorno substitutivo
	Matriz Ret = sm(&M0, &M1);
	*x_ki = Ret;
}
*/
