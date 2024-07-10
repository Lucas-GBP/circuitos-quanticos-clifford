#include "main.hpp"

using namespace std;

int main() {
    srand(time(0));  //Seed para numeros aleatórios
    CliffordSimulator circuit = CliffordSimulator(100);

    test(circuit);

    cout << "Circuito concluído, carregando representações...";

    //cout << circuit;
    //cout << circuit.strPaulli();
    cout << circuit.strKet();
    cout << endl << endl;    

    return 0;
}