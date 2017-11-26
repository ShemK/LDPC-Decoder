#include "soft_decoder.h"
/*
 * LDPC Soft Decoder written for Advanced Digital Communications
 * Authors: Biniyam Zewede and Shem Kikamaze
 *
 */

LDPC::LDPC(int **parityMatrix, double *inputData, int n, int k){
    //this->parityMatrix = parityMatrix;
    this->inputData = inputData;
    this->n = n;
    this->k = k;
    checkNodes = new CheckNode[n-k];
    variableNodes = new VariableNode[n];
    // Setup the graph edges
    // get the number of variable nodes connected to each check node
    for(int j = 0; j < n-k;j++) {
        int sumf = 0;
        for(int i = 0; i < n; i++) {
            sumf = sumf + parityMatrix[j][i];
        }
        checkNodes[j].size = sumf;
        checkNodes[j].variableNodes = new int[sumf];
        checkNodes[j].values = new double[sumf];


        int count = 0;
        for(int i = 0; i < n; i++) {
            if(parityMatrix[j][i] == 1) {
                checkNodes[j].variableNodes[count] = i;
                count++;
            }
        }
    }

    for(int i = 0; i < n; i++) {
        variableNodes[i].size = 0;
        for(int j = 0; j < n-k; j++) {
            variableNodes[i].size = parityMatrix[j][i] + variableNodes[i].size;
        }
        int count = 0;
        variableNodes[i].checkNodes = new int[variableNodes[i].size];
        variableNodes[i].values = new double[variableNodes[i].size];
        for(int j = 0; j < n-k;j++) {
            if(parityMatrix[j][i] == 1) {
                variableNodes[i].checkNodes[count] = j;
                count++;
            }
        }
    }
    // initialize check node values qij;
    for(int j = 0; j < n-k; j++) {
        for(int i = 0; i < checkNodes[j].size; i++){
            checkNodes[j].values[i] = 0; //qij
        }
    }

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < variableNodes[i].size; j++){
            variableNodes[i].values[j] = 0;
        }
    }

}

LDPC::~LDPC(){
    delete[]checkNodes;
    delete[]variableNodes;
}
void LDPC::setVariableNodes(double *inputData) {
    for(int j = 0; j < n-k; j++) {
        cout << "Check Node " << j << endl;
        for(int i = 0; i < checkNodes[j].size; i++){
            checkNodes[j].values[i] = inputData[checkNodes[j].variableNodes[i]] +
                        getValueFromVariableNode(checkNodes[j].variableNodes[i],j); //qij
            cout << checkNodes[j].values[i] << endl;
        }
    }
}

void LDPC::setCheckNodes(double * inputData) {
    for(int i = 0; i < n; i++) {
        cout << "Variable Node " << i << " with size = " << variableNodes[i].size<< endl;
        for(int j = 0; j < variableNodes[i].size; j++){
            variableNodes[i].values[j] = getValueFromCheckNode(variableNodes[i].checkNodes[j],i); //rji
            cout << variableNodes[i].values[j] << endl;
        }
        inputData[i] = setNewInputValue(i,inputData);
    }
}

double LDPC::getValueFromCheckNode(int j, int i) {
    double newValue = 0;
    double multi_val = 1;

    for(int m = 0; m < checkNodes[j].size; m++){
        if(checkNodes[j].variableNodes[m] != i) {
            multi_val = multi_val*tanh(checkNodes[j].values[m]/2);
        }
    }
    
    if(multi_val==1){
        newValue = 2*19.07;
        return newValue;
    }else if(multi_val==-1) {
        newValue = -2*19.07;
        return newValue;
    } else{
        newValue = 2*atanh(multi_val);;
        return newValue;
    }
}

double LDPC::setNewInputValue(int i,double *&input) {
    double sum = this->inputData[i];
    for(int j = 0; j < variableNodes[i].size; j++) {
        sum = sum + variableNodes[i].values[j];
    }
    return sum;
}



void LDPC::display(int ** parityMatrix,int n,int k) {
    std::cout << "-----------------------" << std::endl;
    for(int i = 0; i<n-k;i++){
        for(int j = 0; j < n; j++) {
            std::cout << parityMatrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

double LDPC::getValueFromVariableNode(int i, int j) {
    double sum = 0;
    for(int m = 0; m < variableNodes[i].size; m++){
        if(variableNodes[i].checkNodes[m] != j) {
            sum = sum + variableNodes[i].values[m];
        }
    }
    return sum;
}

/* The gateway function */

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{

    int n = (int) mxGetScalar(prhs[2]);
    int k = (int) mxGetScalar(prhs[3]);
    int iter = (int) mxGetScalar(prhs[4]);
    double *outMatrix;
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)n,mxREAL);
    outMatrix = mxGetPr(plhs[0]);
    double *in_parity = mxGetPr(prhs[0]);
    // outMatrix = mxGetPr(plhs[0]);
    int **parityMatrix; // parity matrix
    double *input = mxGetPr(prhs[1]); // input vector
    parityMatrix = (int **) malloc(sizeof(double)*(n-k));
    for(int r = 0; r<n-k;r++){
        *(parityMatrix+r) = (int*) malloc (n*sizeof(int));
        for(int c = 0; c < n; c++) {
            parityMatrix[r][c] =  (int) in_parity[c*(n-k)+r];
        }
    }
    LDPC myLDPC(parityMatrix,input,n,k);
    
    double *output = new double[n];
    for(int i = 0; i < iter; i++){
        myLDPC.setVariableNodes(input);
        myLDPC.setCheckNodes(outMatrix);
    }
}
