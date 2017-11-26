#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
#include "mex.h"
/*
 * LDPC Soft Decoder written for Advanced Digital Communications
 * Authors: Biniyam Zewede and Shem Kikamaze
 *
 */

using namespace std;

class LDPC {
private:
    int **parityMatrix; // parity matrix
    double *inputData; // input vector
    int k;
    int n;

    struct CheckNode{
        int *variableNodes;
        double *values;
        int size;
    };

    struct VariableNode{
        int *checkNodes;
        double *values;
        int size;
    };

    CheckNode *checkNodes;
    VariableNode *variableNodes;

    /* data */
public:
    LDPC (int **parityMatrix, double *inputData, int n, int k);
    ~LDPC();
    void display(int ** parityMatrix,int n,int k);
    void setVariableNodes(double *inputData);
    void setCheckNodes(double * inputData);
    double getValueFromCheckNode(int j, int i);
    double setNewInputValue(int i,double *&input);
    double getValueFromVariableNode(int i, int j);
};

class LDPC_H_Matrix{
private:
    int **parityMatrix;
    int column_weight;
    int row_weight;
    int k;
    int n;

public:
    LDPC_H_Matrix(int c_w, int n, int k);
    void createMatrix();
    void display();
};
