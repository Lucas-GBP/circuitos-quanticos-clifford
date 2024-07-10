#include "qliff.hpp"

#define TOTAL_SIZE (2*n)    //Tamanho da matriz quadrada G e do vetor
#define BUFFER_INDEX (2*n)  //Index do buffer da matrix e vetor
#define X_INDEX(Q) Q        //Index da coluna X do qubit q
#define Z_INDEX(Q) Q+n      //Index da coluna Z do qubit q
#define D_INDEX(R) R        //Index da linha R dos desestabilizadores
#define S_INDEX(R) R+n      //Index da linha R dos estabilizadores
#define IS_X(I, Q) (G[I][X_INDEX(Q)] && !G[I][Z_INDEX(Q)])
#define IS_Y(I, Q) (G[I][X_INDEX(Q)] && G[I][Z_INDEX(Q)])
#define IS_Z(I, Q) (!G[I][X_INDEX(Q)] && G[I][Z_INDEX(Q)])

CliffordSimulator::CliffordSimulator(QBIT_TYPE numQubits): 
    n(numQubits), 
    G(2 * numQubits + 1, std::vector<MATRIX_TYPE>(2 * numQubits)), 
    F(2 * numQubits + 1) 
    {
    // Inicializar G e F para o estado |0>^n
    for (int i = 0; i < n; ++i) {
        G[D_INDEX(i)][X_INDEX(i)] = 1;  // X
        G[S_INDEX(i)][Z_INDEX(i)] = 1;  // Z
    }
}
void CliffordSimulator::H(QBIT_TYPE qubit) {
    // H troca X e Z
    MATRIX_TYPE tmp;
    for(QBIT_TYPE i = 0; i < TOTAL_SIZE; ++i){
        // Troca as colunas X e Z do qubit
        tmp = G[i][X_INDEX(qubit)];
        G[i][X_INDEX(qubit)] = G[i][Z_INDEX(qubit)];
        G[i][Z_INDEX(qubit)] = tmp;

        // Adiciona a fase se X e Z estiverem presentes (caso Y)
        if (IS_Y(i, qubit)) addPhase(i);
    }
}
void CliffordSimulator::S(QBIT_TYPE qubit) {
    // S mapeia Z -> Z e X -> Y
    for (QBIT_TYPE i = 0; i < TOTAL_SIZE; ++i) {
        // Adiciona a fase se X e Z estiverem presentes (caso Y)
        if (IS_Y(i, qubit)) addPhase(i);

        G[i][Z_INDEX(qubit)] ^= G[i][X_INDEX(qubit)]; //z[i][qubit] ^= x[i][qubit]
    }
}
void CliffordSimulator::CNOT(QBIT_TYPE control, QBIT_TYPE target) {
    // CNOT mapeia X -> XX e Z -> IZ
    for (QBIT_TYPE i = 0; i < TOTAL_SIZE; ++i) {
        G[i][X_INDEX(target)] ^= G[i][X_INDEX(control)];
        G[i][Z_INDEX(control)] ^= G[i][Z_INDEX(target)];

        if ((G[i][X_INDEX(control)] && G[i][Z_INDEX(target)] &&
             G[i][X_INDEX(target)]  && G[i][Z_INDEX(control)])||
            (G[i][X_INDEX(control)] && G[i][Z_INDEX(target)] &&
            !G[i][X_INDEX(target)] && !G[i][Z_INDEX(control)])
        ) {
            addPhase(i); // Adiciona a fase se X e Z estiverem presentes (caso Y)
        }
    }
}
MeasureReturns CliffordSimulator::Measure(QBIT_TYPE qubit, bool suppress) {
    //TODO: - otimizar e criar opção para medir passivamente
    bool indeterminated = false;//O qubit é indeterminado?
    QBIT_TYPE s_pivot;          //linha pivo dos estabilizadores
    QBIT_TYPE d_pivot;          //linha pivo dos desestabilizadores

    for(s_pivot = 0; s_pivot < n; s_pivot++){
        if(G[S_INDEX(s_pivot)][X_INDEX(qubit)]){
            indeterminated = true;
            break;
        }
    }

    if(indeterminated){
        copyRows(D_INDEX(s_pivot), S_INDEX(s_pivot));
        setRow(S_INDEX(s_pivot), S_INDEX(qubit));
        F[Z_INDEX(s_pivot)] = 2*(rand()%2);

        for(QBIT_TYPE i = 0; i < 2*n; i++){
            if(i!=s_pivot && G[i][X_INDEX(qubit)]){
                multRow(i, s_pivot);
            }
        }

        if(F[Z_INDEX(s_pivot)]) return RandomOne;
        return RandomZero;
    }

    if(!indeterminated && !suppress){
        for(d_pivot = 0; d_pivot < n; d_pivot++){
            if(G[d_pivot][X_INDEX(qubit)]) break;
        }
        
        copyRows(BUFFER_INDEX, d_pivot+n);
        for(QBIT_TYPE i = d_pivot+1; i < n; i++){
            if(G[i][qubit]){
                multRow(BUFFER_INDEX, i+n);
            }
        }

        if(F[BUFFER_INDEX]) return AlwaysOne;
        return AlwaysZero;
    }

    return AlwaysZero;
}
void CliffordSimulator::X(QBIT_TYPE qubit){
    H(qubit);
    Z(qubit);
    H(qubit);
}
void CliffordSimulator::Y(QBIT_TYPE qubit){
    S(qubit);
    X(qubit);
    St(qubit);
}
void CliffordSimulator::Z(QBIT_TYPE qubit){
    S(qubit);
    S(qubit);
}
void CliffordSimulator::St(QBIT_TYPE qubit){
    Z(qubit);
    S(qubit);
}
QBIT_TYPE CliffordSimulator::getSize() {
    return n;
}
void CliffordSimulator::minusX(QBIT_TYPE qubit){
    Z(qubit);
    X(qubit);
    Z(qubit);
}
void CliffordSimulator::minusY(QBIT_TYPE qubit){
    Z(qubit);
    Y(qubit);
    Z(qubit);
}
std::string CliffordSimulator::strMatrix(){
    std::string str = "";
    for (int i = 0; i < 2*n; ++i) {
        for (int j = 0; j < 2*n; ++j) {
            str += std::to_string((int) G[i][j]) + " ";
        }
        str += "| " + std::to_string((int) F[i]) + "\n";
    }
    str += "\n";

    return str;
}
std::string CliffordSimulator::strPaulli(){
    std::string str = "";
    int i = 0;
    str += "Destabilizer:{\n";
    for (; i < n; ++i) {
        str += "\t" + paulliRepre(i) + ",\n";
    }

    str += "},\nStabilizer:{\n";
    for(; i < 2*n; i++) {
        str += "\t" + paulliRepre(i) + ",\n";
    }

    str += "}.\n";
    return str;
}
std::string CliffordSimulator::strMesurement(MeasureReturns m){
    return measurementStates[m];
}
std::string CliffordSimulator::strKet(){
    QBIT_TYPE gauss = gaussian();
    QBIT_TYPE states_quant = 1 << gauss;
    std::string str = "\n";
    str += std::to_string(states_quant);
    str += " possiveis estados\n";

    seed(gauss);
    str += strBaseState();

    for(QBIT_TYPE i = 0 ; i < states_quant-1; i++){
        QBIT_TYPE i2 = i ^ (i+1);
        for(QBIT_TYPE j = 0; j < gauss; j++){
            if(i2 & (1<<j)){
                multRow(BUFFER_INDEX, n+j);
            }
        }
        str += strBaseState();
    }

    return str;
}

inline void CliffordSimulator::addPhase(QBIT_TYPE qubit){
    F[qubit] = (F[qubit] + PHASE_QUANT/2)%PHASE_QUANT;
}
inline void CliffordSimulator::swapRows(QBIT_TYPE row1, QBIT_TYPE row2){
    copyRows(BUFFER_INDEX, row2);
    copyRows(row1, row2);
    copyRows(row1, BUFFER_INDEX);
    return;
}
void CliffordSimulator::copyRows(QBIT_TYPE target, QBIT_TYPE control){
    for (int qubit = 0; qubit < n; qubit++){
        G[target][X_INDEX(qubit)] = G[control][X_INDEX(qubit)];
        G[target][Z_INDEX(qubit)] = G[control][Z_INDEX(qubit)];
    }
    F[target] = F[control];
    return;
}
void CliffordSimulator::multRow(QBIT_TYPE target_row, QBIT_TYPE control_row){
    // Ajusta a Fase
    QBIT_TYPE e = 0;// expoente que i esta elevado
    for(QBIT_TYPE i = 0; i < n; i++){
        if(IS_X(control_row, i)){       // X
            if(IS_Y(target_row, i)) e++; // XY=iZ
            else 
            if(IS_Z(target_row, i)) e--;// XZ=-iY
        } else 
        if(IS_Y(control_row, i)){ // Y
            if(IS_Z(target_row, i)) e++;// YZ=iX
            else 
            if(IS_X(target_row, i)) e--;// YX=-iZ
        } else 
        if(IS_Z(control_row, i)){ // Z
            if(IS_X(target_row, i)) e++;// ZX=iY
            else
            if(IS_Y(target_row, i)) e--; // ZY=-iX
        }
    }

    e = (e+F[target_row]+F[control_row])%PHASE_QUANT;
    if(!(e >= 0)){
        e+=PHASE_QUANT;
    }
    F[target_row] = e;

    // Realiza a multiplicação
    for(QBIT_TYPE qubit = 0; qubit < n; qubit++){
        G[target_row][X_INDEX(qubit)] ^= G[control_row][X_INDEX(qubit)];
        G[target_row][Z_INDEX(qubit)] ^= G[control_row][Z_INDEX(qubit)];
    }
}
void CliffordSimulator::setRow(QBIT_TYPE row, QBIT_TYPE obs){    
    // Reseta a linha
    for(QBIT_TYPE i = 0; i < n; i++){
        G[row][X_INDEX(i)] = 0; //X
        G[row][Z_INDEX(i)] = 0; //Z
    }
    F[row] = 0;                 //Y

    G[row][obs] = 1;
    return;
}
void CliffordSimulator::seed(QBIT_TYPE gauss){
    //TODO: otimizar
    QBIT_TYPE min;

    // Limpa o buffer
    F[BUFFER_INDEX] = 0;
    for(QBIT_TYPE i = 0; i < n; i++){
        G[BUFFER_INDEX][X_INDEX(i)] = 0;
        G[BUFFER_INDEX][Z_INDEX(i)] = 0;
    }

    
    for(QBIT_TYPE i = BUFFER_INDEX-1; i >= n+gauss; i--){
        PHASE_TYPE f = F[i];
        for(QBIT_TYPE j = n-1; j >= 0; j--){
            if(G[i][Z_INDEX(j)]){
                min = j;
                if(G[BUFFER_INDEX][X_INDEX(j)]){
                    f = (f+2)%4;
                }
            }
        }
        if(f == 2){
            G[BUFFER_INDEX][X_INDEX(min)] = 1;
        }
    }

    return;
}
void CliffordSimulator::printBuffer(){
    for(QBIT_TYPE i = 0; i < TOTAL_SIZE; i++){
        std::cout << (int) G[BUFFER_INDEX][i] << " "; 
    }
    std::cout << "| " << (int) F[BUFFER_INDEX] << "\n"; 
    return;
}
QBIT_TYPE CliffordSimulator::gaussian(){
    QBIT_TYPE i = n;
    QBIT_TYPE result;

    for(QBIT_TYPE j = 0; j < n; j++){
        for(QBIT_TYPE k = i; k < TOTAL_SIZE; k++){
            if (G[k][X_INDEX(j)]){
                swapRows(i, k);
                swapRows(i-n, k-n);

                for(QBIT_TYPE z = i+1; z < TOTAL_SIZE; z++){
                    if(G[z][X_INDEX(j)]){
                        multRow(z, i);
                        multRow(i-n, z-n);
                    }
                }

                i++;
                break;
            };
        }
    }
    result = i-n;

    for(QBIT_TYPE j = 0; j < n; j++){
        for(QBIT_TYPE k = i; k < TOTAL_SIZE; k++){
            if(G[k][Z_INDEX(j)]){
                swapRows(i, k);
                swapRows(i-n, k-n);

                for(QBIT_TYPE z = i+1; z < TOTAL_SIZE; z++){
                    if(G[z][Z_INDEX(j)]){
                        multRow(z, i);
                        multRow(i-n, z-n);
                    }
                }

                i++;
                break;
            }
        }
    }
    return result;
}
std::string CliffordSimulator::paulliRepre(QBIT_TYPE row){
    std::string str = "";
    str += phaseStates[F[row]];

    for (int qubit = 0; qubit < n; qubit++){
        str += paulliGates[G[row][X_INDEX(qubit)]][G[row][Z_INDEX(qubit)]];
    }

    return str;
}
std::string CliffordSimulator::strBaseState(){
    int e = F[BUFFER_INDEX];
    std::string str = "";

    for(QBIT_TYPE i = 0; i < n; i++){
        if(IS_Y(BUFFER_INDEX, i)){
            e = (e+1)%4;
        }
    }
    str += phaseStates[e];
    str += "|";

    for(QBIT_TYPE i = 0; i < n; i++){
        str += std::to_string(G[BUFFER_INDEX][X_INDEX(i)]);
    }
    str += ">";

    return str;
}
std::ostream& operator << (std::ostream& os, CliffordSimulator& c) {
    std::cout << c.strMatrix();
    return os;
}