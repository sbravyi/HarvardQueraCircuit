import numpy as np

# we consider Harvard/QuEra circuit defined on the Boolean cube {0,1}^k
# each node of the cube contains a triple of logical qubits (Red,Blue,Green)
# see page 29 of https://arxiv.org/pdf/2312.03982.pdf
# for the k=4 example

# Boolean cube dimension
k = 2

# number of cube nodes (denoted m in my notes)
nodes = 1<<k

# total number of qubits
n = 3*nodes

print('Qubits n=',n)

# partition of qubits into Red, Blue, Green
colors = ['Red','Blue','Green']
def my_color(i):
	return colors[i % 3]

Red = [i for i in range(n) if my_color(i)=='Red']
Blue = [i for i in range(n) if my_color(i)=='Blue']
Green = [i for i in range(n) if my_color(i)=='Green']

assert(len(Red)==nodes)
assert(len(Blue)==nodes)
assert(len(Green)==nodes)

# extract indices of red/blue/green variables contained in h
def separate_colors(h):
	red_vars = [i//3 for i in h if my_color(i)=='Red']
	blue_vars = [i//3 for i in h if my_color(i)=='Blue']
	green_vars = [i//3 for i in h if my_color(i)=='Green']
	# when we only care about variables of a fixed color, variable index ranges between 0 and nodes-1
	# hence the division by three
	return red_vars,blue_vars,green_vars


### Functions of manipulations with bit strings
def bit_count(x):
    return bin(x).count("1")

def to_binary(num_bits,x):
	xb = [(x>>q) & 1 for q in range(num_bits)]
	return np.array(xb,dtype=int)

def make_gray_code(number_of_bits):
	assert(number_of_bits>=1)
	if number_of_bits==1:
		return [0]
	C = make_gray_code(number_of_bits-1)
	return C+[number_of_bits-1]+C

# gray_code[i] = which bit to flip to go from i-th to (i+1)-th bit string in the
# Gray code order
# WARNING: exponentially big list !
gray_code = make_gray_code(nodes) + [nodes-1]

### Linear system solver
# solve(A,b) solves a linear system Ax=b (mod 2)
# returns
# status,x,Aperp
# where status='OK' if the system is feasible and 'NULL' otherwise
# x is some solution of the system (if any)
# Aperp is a matrix whose columns span the nullspace of A
# if A is full-rank then Aperp=[[],[]] is array with shape (0,0)
def solve(A,b):
	num_eqs,num_vars = A.shape
	if np.count_nonzero(b)>0:
		# add b to A as the last column
		B = np.hstack((A,np.reshape(b,(num_eqs,1))))
		num_vars+=1
	else:
		B = A.copy()
	# compute the nullspace of B
	# add one linear equation per time
	# initially there are no equations
	# any basis vector is in the nullspace of B
	X = np.identity(num_vars,dtype=int)
	for i in range(num_eqs):
		syndrome = (B[i,:] @ X) % 2
		# columns of X that have inner product=0 with the i-th row of A
		good = X[:,syndrome==0]
		# columns of X that have inner product=1 with the i-th row of A
		bad = X[:, syndrome==1]
		if bad.shape[1]>0:
			ones = np.ones(bad.shape[1],dtype=int)
			# add the first bad column to all other bad columns
			bad^= np.outer(bad[:,0],ones)
			# remove the first bad column
			bad = np.delete(bad, 0, axis=1)
		
		X = np.hstack((good, bad))
		if X.shape[1]==0:
			if np.count_nonzero(b)>0:
				return 'NULL',[],[]
			else:
				# if b is all-zeros vectors then we can get here only if 
				# A is full-rank matrix
				# in this case all-zero vector x is the only solution of the system
				x = np.zeros(num_vars,dtype=int)
				Aperp = np.array([[],[]])
				return 'OK',x,Aperp
	
	if np.count_nonzero(b)==0:
		return 'OK',X[:,0],X

	# any column of X that ends with '1' is a solution x of the system Ax=b (mod 2)
	good = np.nonzero(X[-1,:])[0]
	if len(good)==0:
		# the system is not feasible
		return 'NULL',[],[]
	col = good[0]
	# a solution of the system
	x = X[:,col]
	x = np.delete(x,-1)
	X^= np.outer(X[:,col],X[-1,:])
	
	Aperp = np.delete(X,col,axis=1)
	Aperp = np.delete(Aperp,-1,axis=0)

	return 'OK',x,Aperp


### Functions to simulate CCZ/CZ/CNOT gates via the phase polynomial formalism
# We only consider phase polynomials with degree-2 and degree-3 terms
# A phase polynomial is described by a list M of tuples of variable indices
# Example:
# f(x_0,x_1,x_2,x_3) = x_0*x_1*x_3 + x_1*x_2 (mod 2) 
# is described by M = [(0,1,3),(1,2)]
# A phase polynomial M can be converted to n-qubit state |psi> as follows:
#   N = 1<<n
#	psi = np.ones(N)
#	for i in range(N):
#		x = to_binary(n,i)
#		for h in M:
#			if all([x[i] for i in h]):
#				psi[i]*=-1
#	psi = psi/np.sqrt(N)
def ccz(M,q1,q2,q3):
	assert(not(my_color(q1)==my_color(q2)))
	assert(not(my_color(q1)==my_color(q3)))
	assert(not(my_color(q2)==my_color(q3)))
	supp = [q1,q2,q3]
	supp.sort()
	supp = tuple(supp)
	if supp in M:
		M.remove(supp)
	else:
		M.add(supp)

def cz(M,q1,q2):
	assert(not(my_color(q1)==my_color(q2)))
	supp = [q1,q2]
	supp.sort()
	supp = tuple(supp)
	if supp in M:
		M.remove(supp)
	else:
		M.add(supp)

def cnot(M,con,tar):
	assert(not(con==tar))
	assert(my_color(con)==my_color(tar))
	M1 = M.copy()
	for h in M1:
		if tar in h:
			h1 = [i for i in h if not(i==tar)]+[con]
			h1.sort()
			h1 = tuple(h1)
			if h1 in M:
				M.remove(h1)
			else:
				M.add(h1)

# this is only used for debugging to compare our fast algorithm and brute force simulation
def apply_hadamard(psi,q):
	N = psi.shape[0] 
	Nq = 1<<q
	psiZ = [psi[x]*((-1)**((x>>q) & 1)) for x in range(N)] # Z[q]*psi
	psiX = [psi[x ^ Nq] for x in range(N)] # X[q]*psi
	return (np.array(psiZ) + np.array(psiX))/np.sqrt(2)


# compute the phase polynomial generated by the initial Hadamard layer and -CNOT-CCZ-CZ- layers
# note: for simplicity we assume that all n qubits are initialized in |0>
# Actual Harvard/QuEra circuit initializes red and green qubits in |1>
# while blue qubits are initialized in |0>
# This is equivalent to initializing all n qubits in |0> and flipping red and green qubits
# after the last Hadamard layer
print('Computing the phase polynomial...')

# initial phase polynomial
M = set([])

# apply the initial layer of "A-rectangles", see page 29 in 
# https://arxiv.org/pdf/2312.03982.pdf
for i in range(nodes):
	ccz(M,Red[i],Blue[i],Green[i])
	cz(M,Red[i],Blue[i])
	cz(M,Blue[i],Green[i])
	cz(M,Red[i],Green[i])
	# we ignore pauli Z gates since they can be absorbed into a Pauli frame

for direction in range(k):
	# apply CNOTs oriented along this direction on the cube
	# cube nodes with even pariry = control qubits
	# cube nodes with odd parity = target qubits
	for x in range(nodes):
		if bit_count(x) % 2:
			continue
		y = x ^ (1<<direction)
		cnot(M,Red[x],Red[y])
		cnot(M,Blue[x],Blue[y])
		cnot(M,Green[x],Green[y])

	# alternate between layers of A or B rectangles, see page 29 in 
	# https://arxiv.org/pdf/2312.03982.pdf
	# some A/B rectangles acting on nodes with even parity cancel each other
	for i in range(nodes):
		ccz(M,Red[i],Blue[i],Green[i])
		cz(M,Red[i],Blue[i])
		cz(M,Blue[i],Green[i])
		if direction % 2:
			cz(M,Red[i],Green[i])

# at this point we have the final phase polynomial M

# collect monomials of each type contained in M
M_RBG = {} # degree-3 monomials with red,blue,green variables
M_RB = {} # degree-2 monomials with red,blue variables
M_RG = {} # degree-2 monomials with red,green variables
M_BG = [] # degree-2 monomials with blue,green variables
# monomials that contain a red variable a grouped into bins according to the label
# of the red variable. 
for h in M:
	assert(len(h)>=2) # we should not get degree-1 monomials
	red_vars,blue_vars,green_vars = separate_colors(h)
	assert(len(red_vars)<=1)
	assert(len(blue_vars)<=1)
	assert(len(green_vars)<=1)
	has_red = False
	has_blue = False
	has_green = False
	if len(red_vars)==1:
		has_red = True
		hRed = red_vars[0]
	if len(blue_vars)==1:
		has_blue = True
		hBlue = blue_vars[0]
	if len(green_vars)==1:
		has_green = True
		hGreen = green_vars[0]
	
	if has_red and has_blue and has_green:
		if hRed in M_RBG:
			M_RBG[hRed].append((hBlue,hGreen))
		else:
			M_RBG[hRed]= [(hBlue,hGreen)]
		continue

	if has_red and has_blue and not(has_green):
		if hRed in M_RB:
			M_RB[hRed].append(hBlue)
		else:
			M_RB[hRed] = [hBlue]
		continue

	if has_red and has_green and not(has_blue):
		if hRed in M_RG:
			M_RG[hRed].append(hGreen)
		else:
			M_RG[hRed] = [hGreen]
		continue

	if has_blue and has_green and not(has_red):
		M_BG.append((hBlue,hGreen))
		continue

	# we should not get here
	assert(False)

print('Done')

# compute the output amplitude <s|H^{otimes n}|psi>
# where |psi> is the n-qubit state described by the phase polynomial M



# pick random output basis vector |s>
s = np.array([1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0], dtype=int)
sR = s[Red]
sB = s[Blue]
sG = s[Green]

print('Computing output amplitude <s|U|00...0>')
print('Output bit string s=',s)
print('s(Red)=',sR)
print('s(Blue)=',sB)
print('s(Green)=',sG)
print('This can take a while...')

# Below we use matrices Gamma, delta^B, and delta^G defined on page 2 of my note
amplitude = 0.0
Gamma = np.zeros((nodes,nodes),dtype=int)
for h in M_BG:
	Gamma[h[0],h[1]]^=1
deltaB = np.zeros(nodes,dtype=int)
deltaG = np.zeros(nodes,dtype=int)
xR = np.zeros(nodes,dtype=int)

# WARNING: exponentially big loop ! 
# We iterate over (n/3)-bit strings xR 
# This is where most of the time is spent
for i, flip_bit in enumerate(gray_code):
	# quick test that does not require solving the linear system
	status = (np.sum(xR*(sB^deltaB)) % 2)==0 and (np.sum(xR*(sG^deltaG)) % 2)==0
	if status:
		# I suspect that 99% of the time is spent in this line:
		status,xG,Gperp = solve(Gamma,sB^deltaB)

	# if status=='OK' then 
	# xG is a solution of the linear system Gamma @ xG = sB (mod 2)
	# and columns of Gperp span the nullspace of Gamma
	if status=='OK':
		# check if the vector sG belongs to the nullspace of Gamma
		if Gperp.shape[1]>0:
			syndrome = ((sG^deltaG) @ Gperp) % 2
			# rank of Gamma
			rk = nodes - Gperp.shape[1]
		else:
			# Gamma is a full rank matrix
			syndrome = 0
			rk = nodes

		if np.count_nonzero(syndrome)==0:
			# we got a nonzero contribution to the amplitude from this boolean pattern xR
			phase = (-1)**(np.sum((sG^deltaG)*xG) + np.sum(sR*xR) % 2)
			amplitude+= phase/(1<<rk)
	# update xR and G
	xR[flip_bit]^=1
	for h in M_RBG[flip_bit]:
		Gamma[h[0],h[1]]^=1
	# update deltaB and deltaG
	if flip_bit in M_RB:
		for h in M_RB[flip_bit]:
			deltaB[h]^=1
	if flip_bit in M_RG:
		for h in M_RG[flip_bit]:
			deltaG[h]^=1

# the last flip in the Gray code should bring us back to the all-zero string
assert(np.count_nonzero(xR)==0)
	
amplitude/=(1<<nodes)

print('Amplitude <s|U|00...0>=',amplitude)

# test: compute the same amplitude by the brute force algorithm
if n<=24:
	N = 1<<n
	psi = np.ones(N)
	for i in range(N):
		x = to_binary(n,i)
		for h in M:
			if all([x[i] for i in h]):
				psi[i]*=-1
	psi = psi/np.sqrt(N)
	for q in range(n):
		psi = apply_hadamard(psi,q)
	s_integer = 0
	for q in range(n):
		if s[q]:
			s_integer^= 1<<q
	amplitude_exact = psi[s_integer]
	print('Brute force simulation: amplitude=',amplitude_exact)
	assert(np.abs(amplitude-amplitude_exact)<1e-10)