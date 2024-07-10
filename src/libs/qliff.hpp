#ifndef QLIFF_HPP
#define QLIFF_HPP

#include <iostream>
#include <stdio.h>
#include <inttypes.h>
#include <vector>
#include <string>
#include <bitset>

#define MATRIX_TYPE int_fast8_t //Tipo dos elementos da matriz estabilizadaora
#define PHASE_TYPE int_fast8_t  //TIpo do vetor de fase
#define QBIT_TYPE int64_t       //Tipo da quantidade de qubits
#define PHASE_QUANT 4           //Quantidade de possiveis fases

const std::string phaseStates[PHASE_QUANT] = {  //string das fases possiveis
    " +", "+i", " -", "-i"
};
const char paulliGates[2][2] = { //String das portas de Paulli
    {'I', 'Z'},
    {'X', 'Y'}
};
const std::string measurementStates[4] = {  //representação das possiveis medições
    "0", "1", "0 (random)", "1 (random)"
};

enum MeasureReturns {
    AlwaysZero, //Qubit sempre será 0;
    AlwaysOne,  //Qubit sempre será 1;
    RandomZero, //Qubit esta indeterminado 0;
    RandomOne   //Qubit esta indeterminado 1.
};

/**
 * @class CliffordSimulator
 * @brief Classe que retorna um circuito quântico
 */
class CliffordSimulator {
public:
    /** 
     *  @brief Construtor que inicializa o simulador com o número de qubits especificado
     *  @param numQubits quantidade de qubits totais no circuito
    */
    explicit CliffordSimulator(QBIT_TYPE numQubits);
    /** 
     *  @brief Aplica a porta Hadamard em um qubit específico
     *  @param qubit index do qubit aplicado
    */
    void H(QBIT_TYPE qubit);
    /** 
     *  @brief Aplica a porta de fase S em um qubit específico
     *  @param qubit index do qubit aplicado
    */
    void S(QBIT_TYPE qubit);
    /** 
     *  @brief Aplica a porta CNOT com qubit de controle e alvo especificados
     *  @param control index do qubit de controle
     *  @param target index do qubit alvo
    */
    void CNOT(QBIT_TYPE control, QBIT_TYPE target);
    /**
     *  @brief Faz a medição em um qubit específico
     *  @param qubit index do qubit aplicado
     *  @return 
     *  0 se qubit sempre for 0;
     *  1 se qubit sempre será 1;
     *  2 se qubit esta indeterminado 0;
     *  3 se qubit esta indeterminado 1.
     */
    MeasureReturns Measure(QBIT_TYPE qubit, bool suppress = false);
    /** 
     *  @brief Aplica a porta X (NOT) em um qubit específico
     *  @param qubit index do qubit aplicado
    */
    void X(QBIT_TYPE qubit);
    /** 
     *  @brief Aplicar a porta Y em um qubit específico
     *  @param qubit index do qubit aplicado
    */
    void Y(QBIT_TYPE qubit);
    /** 
     *  @brief Aplica a porta Z em um qubit específico
     *  @param qubit index do qubit aplicado
    */
    void Z(QBIT_TYPE qubit);
    /** 
     *  @brief Aplicar a porta S^\dagger em um qubit especifico
     *  @param qubit index do qubit aplicado
    */
    void St(QBIT_TYPE qubit);
    /** 
     *  @brief Aplicar a porta -X em um qubit especifico
     *  @param qubit index do qubit aplicado
    */
    void minusX(QBIT_TYPE qubit);
    /** 
     *  @brief Aplicar a porta -Y em um qubit especifico
     *  @param qubit index do qubit aplicado
    */
    void minusY(QBIT_TYPE qubit);
    void seed(QBIT_TYPE gauss);
    /** 
     *  @brief Retorna o numero de qbits no circuito 
    */
    QBIT_TYPE getSize();
    /** 
     *  @brief matrix G e vetor F estabilizadores
     *  @return retorna uma string contendo a matrix e o vetor
    */
    std::string strMatrix();
    /** 
     *  @brief estabilizadores e desestabilizadores em matrizes de Palli 
     *  @return retorna uma string contendo os conjuntos estabilizadores e desestabilizadores em matrizes de Paulli
    */
    std::string strPaulli();
    /**
     * @brief resultado da medida
     * @param m mensuração
     * @return retorna uma string com o resultado da medida, e se ela for determinada ou aleatória
     */
    std::string strMesurement(MeasureReturns m);
    //TODO: documetação strKet()
    /**
     * 
     */
    std::string strKet();
private:
    QBIT_TYPE n;                            // Número de qubits
    std::vector<std::vector<MATRIX_TYPE>> G;// Matriz G para a representação estabilizadora
    std::vector<PHASE_TYPE> F;              // Vetor F para a representação estabilizadora
    /**
     * @brief adiciona a fase em um qubit
     * @param qubit index do qubit
     */
    inline void addPhase(QBIT_TYPE qubit);
    /**
     * @brief Troca duas linhas da matriz
     * @param row1 index de uma das linhas
     * @param row2 index de uma das linhas
    */
    inline void swapRows(QBIT_TYPE row1, QBIT_TYPE row2);
    /**
     * @brief copia um linha da matriz em outra
     * @param control index da linha a ser copiada
     * @param target index da linha que abrigará a copia
    */
    void copyRows(QBIT_TYPE control, QBIT_TYPE target);
    /**
     *  @brief multiplica a linha target pela linha control
     *  @param target_row index da linha target
     *  @param control_row index da linha de controle
    */
    void multRow(QBIT_TYPE target_row, QBIT_TYPE control_row);
    /**
     *  @brief defina a linha row com o observavel da posição obs fazendo com que ela seja composta por zeros em todas as posições exceto no index obs onde será 1
     *  @param row index da linha
     *  @param obs index do observavel
    */
    void setRow(QBIT_TYPE row, QBIT_TYPE obs);
    //TODO: documetação printBuffer()
    /**
     * 
     */
    void printBuffer();
    //TODO: documentação gaussian()
    /**
     * @brief gausian( )
     * @return a
    */
    QBIT_TYPE gaussian();
    /** 
     *  @brief Função auxiliar: retorna As matrizes de pauli correspondestes ao row da matriz 
     *  @param row index da linha
    */
    std::string paulliRepre(QBIT_TYPE row);
    /**
     * @brief estado em notação Ket armazenado no buffer
     * 
     * @return retorna uma string com o estado armazenado no buffer em notação ket com a fase
    */
    std::string strBaseState();
    /** 
     *  @brief Função para std::cout 
    */
    friend std::ostream& operator << (std::ostream&, CliffordSimulator&);
};

#endif