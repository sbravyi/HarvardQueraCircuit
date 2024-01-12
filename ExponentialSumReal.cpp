#include <cstdlib>
#include <iostream>
#include <assert.h>     /* assert */

using std::endl;
using std::cout;
using std::vector;
using std::pair;

// we only consider Clifford circuits of the form H-CZ-Z-H
// where H stands for the bitwise Hadamard
// CZ is any layer of CZ gates
// Z is any layer of Z gates
struct clifford_circuit
{ 
    vector<pair<unsigned int,unsigned int> > CZ;// pairs of control,target qubit for each CZ gate
    vector<unsigned int> Z;// qubits acted by Z-gates
};


void print_circuit(clifford_circuit circ)
{   
    cout<<"H-CZ-Z-H circuit:"<<endl;
    cout<<"Bitwise Hadamard"<<endl;
    for (vector<pair<unsigned int, unsigned int> >::iterator it=circ.CZ.begin(); it!=circ.CZ.end(); it++)
        cout<<"CZ "<<(*it).first<<", "<<(*it).second<<endl;
    for (vector<unsigned int>::iterator it=circ.Z.begin(); it!=circ.Z.end(); it++)
        cout<<"Z "<<(*it)<<endl;
    cout<<"Bitwise Hadamard"<<endl;
}



// this function computes zero-zero amplitude <0^n|circ|0^n> where circ is n-qubit H-CZ-Z-H circuit 
// To compute an amplitude <v|circ|0^n> for some n-bit string v
// toggle Z gates in the support of v
// Example for n=2:
// <01|H0 H1 CZ[0,1] H0 H1 |00> = <00|H0 H1 CZ[0,1] Z[1] H0 H1 |00>
//
// we parameterize the amplitude by a triple of variables (bool outZero, unsigned int outP, int outSign) such that
// if (!outZero) amplitude = (1.0/(1<<(n-outP)))*outSign;
// if (outZero) amplitude = 0.0;
// the integer outP is always between 0 and n
// see example in the main() function
void ExponentialSumReal(unsigned int n, clifford_circuit circ, bool &outZero, unsigned int &outP, int &outSign)
{
    // this is an optimized version of the algorithm ExponentialSum, see page 12 of
    // https://arxiv.org/pdf/1601.07601.pdf
    // the algorithm is adapted to the real case (no phase gates)



    long unsigned one=1ul;
    
    
    if (n>64) {
       cout<<"ExponentialSumReal:error, expect n<=64"<<endl;
       exit(1);
    }
    
    size_t pow2=0;
    bool sigma=0;
    bool isZero=0;
    
    // vector L parameterizes Z-gates
    long unsigned L=0ul;
    for (vector<unsigned>::iterator it = circ.Z.begin(); it!=circ.Z.end(); it++)
        L^= (one << (*it));
    
    
    // matrix M parameterizes CZ-gates
    // convert each column to long integer
    long unsigned M[64];
    for (size_t j=0; j<n; j++) M[j]=0ul;
    
    for (vector<pair<unsigned,unsigned> >::iterator it = circ.CZ.begin(); it!=circ.CZ.end(); it++)
    {
        unsigned int con = (*it).first;
        unsigned int tar = (*it).second;
        assert((con<tar));
        M[con]^= (one<<tar);
    }
    
    
    bool active[64];
    for (size_t j=0; j<n; j++)
        active[j]=true;
    
    int nActive=n;
    
    while (nActive>=1)
    {
      // find the first active variable
      size_t i1;
      for (i1=0; i1<n; i1++)
          if (active[i1])
              break;
      
      // find i2 such that M(i1,i2)!=M(i2,i1)
      size_t i2;
      bool isFound=false;
      for (i2=0; i2<n; i2++)
      {
          isFound = ( ((M[i1]>>i2) & one) != ((M[i2]>>i1) & one) );
          if (isFound) break;
      }
      
      bool L1 = ((L>>i1) & one) ^ ((M[i1]>>i1) & one);
      
      // take care of the trivial cases
      if (!isFound)
      {
         // the form is linear in the variable i1
         if (L1)
         {
             outP=0;
             outSign=1;
             outZero=true;
             return;
         }
         else
         {
           pow2+=1;
           nActive-=1;
           // set column i1 to zero
           M[i1]=0ul;
           // set row i1 to zero
           for (size_t j=0; j<n; j++)
               M[j]&=~(one<<i1);
           L&=~(one<<i1);
           active[i1]=0;
           continue;
         }
      }
      
      // Do the recursion
      bool L2 = ( ((L>>i2) & one) ^ ((M[i2]>>i2) & one) );
      L&=~(one<<i1);
      L&=~(one<<i2);
      
      // Extract rows i1 and i2 of M
      long unsigned m1=0ul;
      long unsigned m2=0ul;
      for (size_t j=0; j<n; j++)
      {
         m1^=((M[j]>>i1) & one)<<j;
         m2^=((M[j]>>i2) & one)<<j;
      }
      m1^=M[i1];
      m2^=M[i2];
      
      m1&=~(one<<i1);
      m1&=~(one<<i2);
      m2&=~(one<<i1);
      m2&=~(one<<i2);
      
      // set columns i1,i2 to zero
      M[i1]=0ul;
      M[i2]=0ul;
      // set rows i1,i2 to zero
      for (size_t j=0; j<n; j++)
      {
         M[j]&=~(one<<i1);
         M[j]&=~(one<<i2);
      }
      
      
      // compute new M and L
     if (L1)
        L^=m2;
     
     if (L2)
        L^=m1;
    
     for (size_t j=0; j<n; j++)
         if ((m2 >> j) & one)
             M[j]^=m1;
      
      pow2+=1;
      sigma^=L1 & L2;
      active[i1]=0;
      active[i2]=0;      
      nActive-=2;  
    }// while
    
   outP=pow2;
   outSign=1-2*sigma;
   outZero=false;
   return;
 

}



int main()
{
// number of qubits
unsigned int n = 48;

assert(n<=64);

srand(17);

// print the circuit if true
bool verbosity = true;

// number of circuits to simulate
long unsigned num_circuits = 1;

// generated circuits
vector<clifford_circuit>  C;

cout<<"Number of Clifford circuits="<<num_circuits<<endl;
cout<<"Generating random circuits..."<<endl;
for (size_t i=0; i<num_circuits; i++)
{

    // define H-CZ-Z-H circuit on n qubits
    // bitwise Hadamard on the left and on the right
    // arbitrary layer of CZ and Z gates in the middle

    // pick random CZ gates
    vector<pair<unsigned int,unsigned int> > CZ;
    vector<unsigned int> Z;
    for (size_t con=0; con<n; con++)
    for (size_t tar=(con+1); tar<n; tar++)
    {
        if (rand() % 2) CZ.push_back(pair<unsigned int,unsigned int>(con,tar));

    }
    // pick random Z gates
    for (size_t qubit=0; qubit<n; qubit++)
         if (rand() % 2) Z.push_back(qubit);
     clifford_circuit circ;
     circ.CZ = CZ;
     circ.Z = Z;
     C.push_back(circ);
}
cout<<"Done."<<endl;

bool outZero;
unsigned int outP;
int outSign;

cout<<"Begin simulation"<<endl;
for (vector<clifford_circuit>::iterator it=C.begin(); it!=C.end(); it++)
{
// compute the amplitude
ExponentialSumReal(n, *it, outZero, outP, outSign);
assert(outP<=n);
double amplitude = 0.0;
if (!outZero) amplitude = (1.0/(1<<(n-outP)))*outSign;

if (verbosity)
{
print_circuit(*it);
cout<<"amplitude="<<amplitude<<endl;
}
}
cout<<"Done"<<endl;

}


