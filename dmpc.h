#ifndef DMPC_H
#define DMPC_H
#include "matriz.h"

Matriz func_deltau(Matriz*, Matriz*, Matriz*, Matriz*, double, Matriz*);
Matriz func_Phi(Matriz*, Matriz*, Matriz*, int, int);
Matriz func_F(int, Matriz*, Matriz*);
Matriz ret_A(Matriz* Am, Matriz* Cm);
Matriz ret_B(Matriz* Bm, Matriz* Cm);
Matriz ret_C(Matriz* Cm);
void discret(double, Matriz*, Matriz*, Matriz*, Matriz*);
void model_estend(Matriz*, Matriz*, Matriz*, Matriz*, Matriz*, Matriz*);

#endif
