#include "bitwise.h"

#define QUBIT_TYPE uint64_t
#define QUBIT_LENGHT sizeof(QUBIT_TYPE)

void bit_wise_test(){
    const uint64_t qubit_quant = 100;
    const QUBIT_TYPE* qubits = (QUBIT_TYPE*) malloc(ceil((2*qubit_quant+1)/QUBIT_LENGHT));



    free(qubits);
    return;
}