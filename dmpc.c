#include "dmpc.h"

Matriz func_Phi(Matriz* A, Matriz* B, Matriz* C, int Np, int Nc){
//Matriz deslizada de Toeplitz formada a partir de um vetor de fatores de
//potência a serem deslizados
	Matriz Phi = matriz((double[]){}, 0, 0);
	Matriz Phi_primario = matriz((double[]){}, 0, 0);
	
	//Crição do vetor (1xNp) a ser deslizado para a criação de Toeplitz
	for(int i = 1; i <= Np; i++){
		if((Np-i) != 0) {
			int pot = Np - i;
			Matriz MA = mat_pot(A, pot);
			MA = mat_mult(C, &MA);
			MA = mat_mult(&MA, B);
			Phi_primario = extend(&MA, &Phi_primario);
		} else {
			Matriz MA = mat_mult(C, B);
			Phi_primario = extend(&MA, &Phi_primario); 
		}
	}

	double* flip(double* array, int t){
	//Função interna para a fazer a inversão do array 
	//da matriz primária
		double* array_invertido = malloc(t*sizeof(double));
		for(int i = 0; i < t; i++)
			array_invertido[t - 1 - i] = array[i];
		return array_invertido;
	}

	int tamanho = Phi_primario.linhas*Phi_primario.colunas;
	Phi_primario = matriz(flip(Phi_primario.matriz, tamanho),
				   Phi_primario.linhas,
				   Phi_primario.colunas);
	
	Phi_primario.print(&Phi_primario);
	double *phi_mat = Phi_primario.matriz;
	
	//Primeira linha da matriz de Toeplitz é copiada como a primeira linha
	//da resposta para que a resposta final não seja uma matriz nula	
	Matriz res = zeros(1, Nc);
	memcpy(res.matriz, phi_mat+Np-1, 1*sizeof(double));

	for(int i = 1; i < Np; i++){
	//A quantidade de colunas deve ser garantidamente do tamanho de Nc, sendo
	//garantido tanto a partir de ind_max
		int ind_max = i+1;
		if(ind_max > Nc) ind_max = Nc; //quant. cols. não maior que Nc
		Matriz Toep = zeros(1, Nc);
		memcpy(Toep.matriz, phi_mat+Np-i-1, (ind_max)*sizeof(double));
		res = extend(&res, &Toep);
		free(Toep.matriz);
	}
	free(Phi_primario.matriz);
	return res;
}

Matriz func_F(int Np, Matriz* C, Matriz* A){
//Entrega a matriz compacta F (representativa das operações inerentes
//ao sistema
	Matriz F = matriz((double[]){}, 0, 0);	
	for(int i = 0; i < Np; i++){
		Matriz M0 = mat_pot(A, i+1);
		Matriz M1 = mat_mult(C, &M0);
		extend(&F, &M1);
		free(M0.matriz);
		free(M1.matriz);
	}
	return F;
}

void model_estend(Matriz* Am, Matriz* Bm, Matriz* Cm, Matriz *A, Matriz* B, Matriz*C){
//Entrega versões estenidas das matrizes A, B e C a partir das matrizes do modelo
//Am, Bm e Cm
	Matriz M1 = zeros(Am->linhas, 1);          
	Matriz M2 = extend(Am, &M1);               
	Matriz CmAm = mat_mult(Cm, Am);            
	Matriz Intg = matriz((double[]){1}, 1, 1);
	Matriz M3 = extend(&CmAm, &Intg);          
	
	double* a_mat = malloc(M2.colunas*(M3.linhas + M2.linhas)*sizeof(double));
	memcpy(a_mat, M2.matriz, M2.linhas*M2.colunas*sizeof(double));
	memcpy(a_mat + M2.linhas*M2.colunas, M3.matriz,
	       M3.linhas*M3.colunas*sizeof(double));
	*A = matriz(a_mat, M3.linhas + M2.linhas, M2.colunas);
	free(M1.matriz);
	free(M2.matriz);
	free(CmAm.matriz);
	free(M3.matriz);

	Matriz CmBm = mat_mult(Cm, Bm);
	Matriz Bm_t = transposta(Bm);
	CmBm = extend(&Bm_t, &CmBm);
	*B = transposta(&CmBm);
	free(CmBm.matriz);
	free(Bm_t.matriz);

	Matriz M4 = zeros(1, Cm->colunas);
	Matriz M5 = matriz((double[]){1}, 1, 1);
	*C = extend(&M4, &M5);
	free(M4.matriz);
}

Matriz ret_A(Matriz* Am, Matriz* Cm){
//Entrega da versão estendida de A para CPDM
	Matriz M0 = zeros(Am->linhas, 1);
	Matriz CmAm = mat_mult(Cm, Am);
	Matriz Intg = matriz((double[]){1}, 1, 1);

	Matriz M1 = extend(Am, &M0);
	Matriz M2 = extend(&CmAm, &Intg);
	int tm = M1.linhas*M1.colunas + M2.linhas*M2.colunas;
	double *mat_res = malloc(tm*sizeof(double));
	memcpy(mat_res, M1.matriz, M1.linhas*M1.colunas*sizeof(double));
	memcpy(mat_res + M1.linhas*M1.colunas, M2.matriz,
	       M2.linhas*M2.colunas*sizeof(double));
	free(M1.matriz);
	free(M2.matriz);
	free(M0.matriz);
	free(CmAm.matriz);
	return matriz(mat_res, M1.linhas + M2.linhas, M2.colunas);
}

Matriz ret_B(Matriz* Bm, Matriz* Cm){
//Entrega da versão estendida de B para CPDM
	Matriz CmBm = mat_mult(Cm, Bm);
	double *res = malloc((CmBm.linhas + Bm->linhas)*sizeof(double));
	memcpy(res, Bm->matriz, Bm->linhas*sizeof(double));
	memcpy(res + Bm->linhas, CmBm.matriz, CmBm.linhas*sizeof(double));
	return matriz(res, Bm->linhas + CmBm.linhas ,1);
}

Matriz ret_C(Matriz* Cm){
//Entrega da versão estendida de C para CPDM
	Matriz M0 = zeros(1, Cm->colunas + 1);
	M0.matriz[Cm->colunas] = 1;
	return M0;
}

void discret(double T, Matriz* Am, Matriz* Bm, Matriz* Ad, Matriz* Bd){
//Modifica os valores internos das matrizes Bd, e Ad a partir da discretização
//de Tustin feita sobre as Matrizes Am, Bm e sobre o valor de tempo amostrado T
	Matriz I = ident(Am->colunas);  
	Matriz A_modif = mult(Am, T/2); //(a*A)
	Matriz M0 = sub(&I, &A_modif);  //(I - a*A);
	Matriz M1 = inv(&M0);           //inv(I - a*A);
	Matriz M2 = sm(&I, &A_modif);   //(I + a*A);
	*Ad = mat_mult(&M1, &M2);

	Matriz M3 = mult(Bm, T);        //B*T
	*Bd = mat_mult(&M1, &M3);       //inv(I - a*A)*B*T
	
	free(I.matriz);
	free(M0.matriz);
	free(M1.matriz);
	free(M2.matriz);
}

