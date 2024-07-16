#include "circuits.hpp"

const std::string str_error = "Circuito não contem numero minimo de qubits";

void test(CliffordSimulator& c){
    if(c.getSize() < 1){
        std::cout << str_error;
        return;
    }

    c.H(0);
    for(int i = 1; i < c.getSize(); i++){
        c.CNOT(0, i);
    }
}

void ghz(CliffordSimulator& c){
    if(c.getSize() < 3){
        std::cout << str_error;
        return;
    }

    c.H(0);
    c.H(1);
    c.CNOT(0, 2);
    c.CNOT(1, 2);
    c.S(0);
    c.S(1);
    c.S(2);

    c.H(0);
    c.H(1);
}