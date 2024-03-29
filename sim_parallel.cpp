// computes output amplitude <s|U|00...0> where U is QuEra-Harward circuit
// acting on 3*2^k qubits
//
// define k at line 22; define s at line 358


#include <cstdlib>
#include <iostream>
#include <cmath>
#include <set>
#include <vector>
#include <chrono>
#include <cassert>
#include <thread>
#include <future>
#include <algorithm>

using std::endl;
using std::cout;
using std::array;
using std::set;
using std::vector;
using std::pair;

// dimension of the hypercube
#define k 5

const long unsigned one = 1ul;


// number of nodes in the hypercube = 2^k
#define num_nodes (1<<k)

// total number of qubits in the QuEra circuit = 3*2^k
#define num_qubits 3*num_nodes

// number of qubits in the Clifford simulator = 2^{k+1}
#define num_qubits_clif 2*num_nodes

// data structure to describe degree-3 polynomials with binary variables and binary coefficients
// each element of the outer set defines a monomial, say, x_0*x_2*x_5
// each element of the inner set defines a variable in the monomial
typedef set<set<unsigned> > phase_poly;

// toggle the coefficient of a given monomial in the phase polynomial P
void toggle(phase_poly &P, set<unsigned> monomial)
{
    if (P.find(monomial)==P.end()) P.insert(monomial); else P.erase(monomial);
}


// apply CCZ to a triple of qubits q1,q2,q3
void apply_ccz(phase_poly &P, unsigned q1, unsigned q2, unsigned q3)
{   
    set<unsigned> S;
    S.insert(q1);
    S.insert(q2);
    S.insert(q3);
    toggle(P,S);
}

// apply CZ to a pair of qubits q1,q2
void apply_cz(phase_poly &P, unsigned q1, unsigned q2)
{   
    set<unsigned> S;
    S.insert(q1);
    S.insert(q2);
    toggle(P,S);
}

// apply Z to a qubits q1
void apply_z(phase_poly &P, unsigned q1)
{   
    set<unsigned> S;
    S.insert(q1);
    toggle(P,S);
}

// apply CNOT with control=con and target=tar
void apply_cnot(phase_poly &P, unsigned con, unsigned tar)
{   
    for (set<set<unsigned> >::iterator it=P.begin(); it!=P.end(); ++it)
    {
        // check if the monomial *it contains the target qubit
        if ((*it).find(tar)!=(*it).end())
        {
            set<unsigned> T(*it);
            T.erase(tar);
            // assert(T.find(con)==T.end());// QuEra circuit never has control and target in the same monomial
            T.insert(con);
            toggle(P,T);
        }
    }
}


// data structure for amplitudes of -H-CZ-Z-H- circuits  
// such circuits always have real amplitudes of the form sign*2^{pow2} 
// where sign is 0 or 1 or -1; pow2 is an integer such that -n<=pow2<=0 
// if sign=0 then pow2 can be ignored
struct clifford_amplitude
{
    int sign;
    int pow2;
};

// Clifford circuit of the form H-CZ-Z-H with n<=64 qubits
// Here H denotes biwise Hadamard
// L parameterizes -Z- layer; Apply Z to i-th qubit if ((L>>i) & 1)==1
// M parameterizes -CZ- layer; Apply CZ to i-th and j-th qubit if ((M[i]>>j) & 1)==1
struct clifford_circuit
{   
    
    long unsigned M[num_qubits_clif];// CZ layer
    long unsigned L;// Z layer

    clifford_circuit()
    {   
        L=0ul;
        for (unsigned i=0; i<num_qubits_clif; i++) M[i]=0ul;
    }

};

// Compute amplitude <0^n|C|0^n> where C is n-qubit H-CZ-Z-H circuit 
// To compute an amplitude <v|C|0^n> for some n-bit string v
// toggle Z gates in the support of v
// Example for n=2:
// <01|H0 H1 CZ[0,1] H0 H1 |00> = <00|H0 H1 CZ[0,1] Z[1] H0 H1 |00>
//
// We implement the algorithm described on pages 25,26 of
// https://arxiv.org/pdf/1808.00128.pdf
// M and L data of clifford_circuit encode matrix M and vector L defined in the above paper
clifford_amplitude ExponentialSumReal(clifford_circuit C)
{    
    unsigned n = num_qubits_clif;

    clifford_amplitude a_out;
    a_out.sign = 0;
    a_out.pow2 = 0;

    long unsigned one=1ul;
    
    size_t pow2=0;
    bool sigma=0;
    bool isZero=0;
     
    bool active[64];
    for (size_t j=0; j<n; j++)
        active[j]=true;
    
    unsigned nActive=n;
    
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
          isFound = ( ((C.M[i1]>>i2) & one) != ((C.M[i2]>>i1) & one) );
          if (isFound) break;
      }
      
      bool L1 = ((C.L>>i1) & one) ^ ((C.M[i1]>>i1) & one);
      
      // take care of the trivial cases
      if (!isFound)
      {
         // the form is linear in the variable i1
         if (L1)
         {  
             a_out.sign = 0;
             a_out.pow2 = 0;
             return a_out;
         }
         else
         {
           pow2+=1;
           nActive-=1;
           // set column i1 to zero
           C.M[i1]=0ul;
           // set row i1 to zero
           for (size_t j=0; j<n; j++)
               C.M[j]&=~(one<<i1);
           C.L&=~(one<<i1);
           active[i1]=0;
           continue;
         }
      }
      
      // Do the recursion
      bool L2 = ( ((C.L>>i2) & one) ^ ((C.M[i2]>>i2) & one) );
      C.L&=~(one<<i1);
      C.L&=~(one<<i2);
      
      // Extract rows i1 and i2 of M
      long unsigned m1=0ul;
      long unsigned m2=0ul;
      for (size_t j=0; j<n; j++)
      {
         m1^=((C.M[j]>>i1) & one)<<j;
         m2^=((C.M[j]>>i2) & one)<<j;
      }
      m1^=C.M[i1];
      m2^=C.M[i2];
      
      m1&=~(one<<i1);
      m1&=~(one<<i2);
      m2&=~(one<<i1);
      m2&=~(one<<i2);
      
      // set columns i1,i2 to zero
      C.M[i1]=0ul;
      C.M[i2]=0ul;
      // set rows i1,i2 to zero
      for (size_t j=0; j<n; j++)
      {
         C.M[j]&=~(one<<i1);
         C.M[j]&=~(one<<i2);
      }
      
      
      // compute new M and L
     if (L1)
        C.L^=m2;
     
     if (L2)
        C.L^=m1;
    
     for (size_t j=0; j<n; j++)
         if ((m2 >> j) & one)
             C.M[j]^=m1;
      
      pow2+=1;
      sigma^=L1 & L2;
      active[i1]=0;
      active[i2]=0;      
      nActive-=2;  
    }// while

    a_out.sign = 1-2*sigma;
    a_out.pow2 = pow2 - n;
    // assert(a_out.pow2<=0);
   return a_out;

}
void print_phase_poly(phase_poly &P)
{

    for (set<set<unsigned> >::iterator it=P.begin(); it!=P.end(); ++it)
    {
        cout<<"Monomial=(";
        for (set<unsigned>::iterator it1=(*it).begin(); it1!=(*it).end(); ++it1)
            cout<<(*it1)<<",";
        cout<<")"<<endl;
    }
}

// convert qubit index ranging between 0 and 3*2^{2^k}
// to index ranging between 0 and 2^{2^k} for red qubits
// and index ranging between 0 and 2*2^{2^k} for blue and green qubits
// Note: we use red qubits for slicing QuEra circuit
// each slice defines a -H-CZ-Z-H- circuit acting on blue and green qubits
unsigned qubit_index(unsigned qubit)
{   
    if ((qubit % 3)==0) return unsigned((qubit/3));// red qubit
    if ((qubit % 3)==1) return unsigned(((qubit-1)/3));// blue qubit
    // green qubit
    unsigned out = unsigned(((qubit-2)/3));
    out+= num_nodes;
    return out;
}


double exponential_task(std::tuple<unsigned long, unsigned long> boundaries, clifford_circuit C, unsigned long sR, const unsigned long* P1, const unsigned long (*P2)[num_qubits_clif])
{ 
    unsigned long start = std::get<0>(boundaries);
    unsigned long end = std::get<1>(boundaries);
    double amplitude = 0.0;
    // initial circuit population
    if (start != 1) {
        long unsigned before_start = start - 1;
        long unsigned starting_gray_code = before_start ^ (before_start>>1);
        for (unsigned bit_idx = 0; bit_idx < num_nodes; bit_idx++) {
            long unsigned starting_code_bitval = starting_gray_code & (one << bit_idx);
            if (starting_code_bitval) {
                for (unsigned q=0; q<num_qubits_clif; q++) {
                    unsigned long pval = P2[bit_idx][q];
                    C.M[q]^= pval;
                }
                C.L^= P1[bit_idx];
            }
        }
    }
    for (long unsigned x = start; x<end; x++) {
        // y = gray code encoding of x
        long unsigned y = x ^ (x>>1);
        long unsigned xprev = x - one;
        long unsigned yprev = xprev ^ (xprev>>1);
        // u = bit where gray_code(x) and gray_code(x-1) differ
        unsigned u = __builtin_ffs (y ^ yprev);  
        // assert(u>=1);
        u-=1;
        // assert(u>=0);
        // assert(u<num_nodes);
        for (unsigned q=0; q<num_qubits_clif; q++) {
            unsigned long pval = P2[u][q];
            C.M[q]^= pval;
        }
        C.L^= P1[u];

        // quick test that can detect -H-CZ-Z-H- circuit with zero amplitude
        bool test1 = ((__builtin_popcountl(y & C.L) % 2)==0);
        bool test2 = ((__builtin_popcountl(y & (C.L>>num_nodes)) % 2)==0);
        if (test1 && test2)
        {
            clifford_amplitude a = ExponentialSumReal(C);// this is likely to be the most expensive step
            // assert((num_nodes-a.pow2)>=0);
            if (a.sign!=0) {
                int overlap = (__builtin_popcountl(sR & y) % 2);
                double amp_inc = ((a.sign)*(1-2*overlap)*(1.0/double(one<<(num_nodes-a.pow2))));
                // cout << "amplitude change on " << x << " from:(" << amplitude << ")" << endl;
                // cout << "amp inc(" << amp_inc << "," << x << ")" << endl;
                amplitude += amp_inc;
            }
        }
    }
    return amplitude;
}

int main()
{

    // std::vector<bool> s = {true, true, true};
    // cout<<"Qubits="<<num_qubits<<endl;
    // cout<<"output string s=0b";
    // for (std::vector<bool>::reverse_iterator it = s.rbegin(); it != s.rend(); ++it) {
    //     std::cout << (*it ? "1" : "0") << " ";
    // }
    // cout<<endl;
    auto begin = std::chrono::high_resolution_clock::now();


// partition 3*2^k qubits into red, blue, and green. There are 2^k qubits of each color.
vector<unsigned> Red;
vector<unsigned> Blue;
vector<unsigned> Green;
for (unsigned i=0; i<num_nodes; i++)
{
    Red.push_back(3*i);
    Blue.push_back(3*i+1);
    Green.push_back(3*i+2);
}



// phase polynomial describing the output state of QuEra circuit immediately before
// the final H-layer
phase_poly P;


// apply the initial layer of "A-rectangles", see page 29 in 
// https://arxiv.org/pdf/2312.03982.pdf
for (unsigned i=0; i<num_nodes; i++)
{
    apply_ccz(P,Red[i],Blue[i],Green[i]);
    apply_cz(P,Red[i],Blue[i]);
    apply_cz(P,Blue[i],Green[i]);
    apply_cz(P,Red[i],Green[i]);
    // we ignore pauli Z gates since they can be absorbed into a Pauli frame
}


for (unsigned direction=0; direction<k; direction++)
{
    // apply CNOTs oriented along this direction on the cube
    // cube nodes with even pariry = control qubits
    // cube nodes with odd parity = target qubits
    for (unsigned x=0; x<num_nodes; x++)
    {
        if ((__builtin_popcount(x) % 2)==0)
        {
            unsigned y = x ^ (1<<direction);
            apply_cnot(P,Red[x],Red[y]);
            apply_cnot(P,Blue[x],Blue[y]);
            apply_cnot(P,Green[x],Green[y]);
        }
    }

    // alternate between layers of A or B rectangles, see page 29 in 
    // https://arxiv.org/pdf/2312.03982.pdf
    // some A/B rectangles acting on nodes with even parity cancel each other
    for (unsigned i=0; i<num_nodes; i++)
    {
        apply_ccz(P,Red[i],Blue[i],Green[i]);
        apply_cz(P,Red[i],Blue[i]);
        apply_cz(P,Blue[i],Green[i]);
        if (direction % 2) apply_cz(P,Red[i],Green[i]);
        // we ignore pauli Z gates since they can be absorbed into a Pauli frame
    }
    
}




// // project s onto red, blue, and green qubits
// long unsigned sR = 0ul;
// long unsigned sB = 0ul;
// long unsigned sG = 0ul;
// unsigned s_index = 0;
// for (unsigned i = 0; s_index < s.size(); i++)
// {
//     sR^= (s[s_index] & one)<<i;
//     sB^= (s[s_index] & one)<<i;
//     sG^= (s[s_index] & one)<<i;
//     s_index += 3;
// }

// define output basis vector |s> of the QuEra circuit
long unsigned s = 123;
cout<<"Qubits="<<num_qubits<<endl;
cout<<"output string s="<<s<<endl;

// project s onto red, blue, and green qubits
long unsigned sR = 0ul;
long unsigned sB = 0ul;
long unsigned sG = 0ul;
for (unsigned i=0; i<num_nodes; i++)
{
    sR^= ((s>>(3*i)) & one)<<i;
    sB^= ((s>>(3*i+1)) & one)<<i;
    sG^= ((s>>(3*i+2)) & one)<<i;
}

// initial -H-CZ-Z-H- circuit on blue+green qubits. All red qubits are set to zero.
clifford_circuit C;
C.L = sB ^ (sG<<num_nodes);


// repackage the phase polynomial 
// group monomials that contain a given red variable
long unsigned P2[num_nodes][num_qubits_clif] = { {0ul} }; 
long unsigned  P1[num_nodes] = {0ul};

for (set<set<unsigned> >::iterator it=P.begin(); it!=P.end(); ++it)
{   
    // we should not get linear terms
    //assert((*it).size()>=1);
    unsigned red=0;
    unsigned blue=0;
    unsigned green=0;
    bool has_red = false;
    bool has_blue = false;
    bool has_green = false;
    for (set<unsigned>::iterator it1=(*it).begin(); it1!=(*it).end(); it1++)
    {
        if (((*it1) % 3)==0) {red=qubit_index(*it1); has_red=true;}
        if (((*it1) % 3)==1) {blue=qubit_index(*it1); has_blue=true;}
        if (((*it1) % 3)==2) {green=qubit_index(*it1); has_green=true;}
    }
    
    //assert(has_red || has_blue || has_green);

    if ( (has_red) && (has_blue) && (has_green) ) P2[red][blue]^= (1<<green);
    if ( (!has_red) && (has_blue) && (has_green) ) C.M[blue]^= (1<<green);
    if ( (has_red) && (!has_blue) && (has_green) ) P1[red]^= (1<<green);
    if ( (has_red) && (has_blue) && (!has_green) ) P1[red]^= (1<<blue);
}




// output amplitude is a sum over 2^{2^k} -H-CZ-Z-H- circuits on blue+green qubits
// iterate over basis vectors on red qubits
long unsigned N = one<<num_nodes;

clifford_amplitude a = ExponentialSumReal(C);
double amplitude = 0.0;
if (a.sign!=0) amplitude = 1.0*(a.sign)/(one<<(num_nodes-a.pow2));
// iterate over gray code index of bit strings of length num_nodes
// has to be a power of two to evenly divide the set
const unsigned long N_TASKS = std::min(N/4, (1UL << 7));
std::future<double> futures[N_TASKS];
for (unsigned long i = 0; i < N_TASKS; ++i) {
    unsigned long n_multiple = N / N_TASKS;
    unsigned long start;
    unsigned long end = n_multiple * (i + 1);
    if (i == 0) {
        start = 1;
    } else {
        start = n_multiple * i;
    }
    std::tuple<long unsigned, long unsigned> bitstring_boundaries = std::make_tuple(
        start,
        end
    );
    futures[i] = std::async(std::launch::async, [=] {
        return exponential_task(bitstring_boundaries, C, sR, std::ref(P1), std::ref(P2));
    });
}
// Wait for all the tasks to complete
for (unsigned i = 0; i < N_TASKS; ++i) {
    futures[i].wait();
    amplitude += futures[i].get();
}
auto end = std::chrono::high_resolution_clock::now();
auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
printf("Time measured: %.5f seconds.\n", elapsed.count() * 1e-9);
cout<<"output amplitude="<<amplitude<<endl;

}


