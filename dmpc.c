#include "dmpc.h"

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




