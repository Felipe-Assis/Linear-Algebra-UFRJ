from math import *

def f1(V,i):
    x,y,z = V[0][0], V[1][0], V[2][0]
    if (i==1):
        resultado = 16.0*(x**4) + 16.0*(y**4) + 1.0*(z**4) - 16.0
        return resultado
    
    if (i==2):
        resultado = 1.0*(x**2) + 1.0*(y**2) + 1.0*(z**2) - 3.0
        return resultado
    
    if (i==3):
        resultado = 1.0*(x**3) - (1.0*y) + 1.0*z - 1.0
        return resultado


def df1(V,i):
    x,y,z = V[0][0], V[1][0], V[2][0]
    if (i==1):
        res1 = 64.0*(x**3)
        res2 = 64.0*(y**3)
        res3 = 4.0*(z**3)
        return [res1, res2, res3]

    if (i==2):
        res1 = 2.0*(x)
        res2 = 2.0*(y)
        res3 = 2.0*z
        return [res1, res2, res3]

    if (i==3):
        res1 = 3.0*(x**2)
        res2 = -1.0
        res3 = 1.0
        return [res1, res2, res3]


def f2(V,i):
	#initialize variables
	c2, c3, c4 = V[0][0], V[1][0], V[2][0]
	teta1, teta2 = 0, 11.667

	if (i==1):
		ret = 2*(c3**2) + c2**2 + 6*(c4**2) - 1
		return ret
	if (i==2):
		ret = 8*(c3**3) + 6*c3*(c2**2) + 36*c3*c2*c4 + 108*c3*(c4**2) - teta1
		return ret
	if (i==3):
		ret = 60*(c4**4) + 60*(c3**2)*(c2**2) + 576*(c3**2)*c2*c4 + 2232*(c3**2)*(c4**2) + 252*(c4**2)*(c2**2) + 1296*(c4**3)*c2 + 3348*(c4**4) + 24*(c2**3)*c4 + 3*c2 - teta2
		return ret


def df2(V,i):
    c2, c3, c4 = V[0][0], V[1][0], V[2][0]
    teta1, teta2 = 0, 11667.0

    #derivadas em relação a c2
    if (i==1):
        res1 = 2*c2
        res2 = 4*c3
        res3 = 12*c4
        return [res1, res2, res3]
    
    #derivadas em relação a c3
    if (i==2):
        res1 = 12*c3*c2 + 36*c3*c4
        res2 = 24*(c3**2) + 6*(c2**2) + 108*(c4**2)
        res3 = 36*c3*c2 + 216*c3*c4
        return [res1, res2, res3]
    
    #derivadas em relação a c4
    if (i==3):
        res1 = 20*(c3**2)*c2 + 576*(c3**2)*c4 + 504*(c4**2)*c2 + 1296*(c4**3)
        + 72*(c2**2)*c4 + 3
        res2 = 240*(c3**3) + 120*c3*(c2**2) + 1152*c3*c2*c4 + 4464*c3*(c4**2)
        res3 = 576*(c3**2)*c2 + 4464*(c3**2)*c4 + 504*c4*(c2**2) + 3888*(c4**2)*c2 + 13392*(c4**3) + 24*(c2**3)
        return [res1, res2, res3]


#Função a ser ajustada
def f3(V,i):
    b0,b1,b2 = V[0][0], V[1][0], V[2][0]
    x = [1.0,2.0,3.0]
    y = [1.0,2.0,9.0]
    ret = b0 + b1*(x[i]**b2) - y[i]
    print ret
    return ret

def df3(V,i):
    b0,b1,b2 = V[0][0], V[1][0], V[2][0]
    x = [1.0,2.0,3.0]
    ret1 = 1.0
    ret2 = x[i]**b2
    ret3 = (b1*(x[i]**b2)*log(x[i]))
    return [ret1, ret2, ret3]

def euclidean(X):
    total = 0;
    for x in X:
        total += 1.0*x[0]**2
    return sqrt(total)


def det2x2(M):
        ret = (M[0][0]*M[1][1] - M[0][1]*M[1][0])
        return ret
   

def det3x3(M):
    a = M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
    b = M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
    c = M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0])
    return (a-b+c)


def inverse3x3(M):
    a, b, c = M[0][0], M[0][1], M[0][2]
    d, e, f = M[1][0], M[1][1], M[1][2]
    g, h, i = M[2][0], M[2][1], M[2][2]
    A = [[e*i-f*h, c*h-b*i, b*f-c*e],
        [f*g-d*i, a*i-c*g, c*d-a*f],
        [d*h-e*g, b*g-a*h, a*e-b*d]]

    det = det3x3(M)

    for i in range(len(M)):
        for j in range(len(M)):
            A[i][j] = A[i][j]*(1.0/det)
    return A


def cholesky(M):
        n = len(M)
        L = [[0.0]*n for i in xrange(n)]
        for i in xrange(n):
                for k in range(i+1):
                    temp = sum(L[i][j]*L[k][j] for j in xrange(k))

                if (i==k):
                        L[i][k] = sqrt(M[i][i] - temp)
                
                else:
                        L[i][k] = (1.0/L[k][k] * (M[i][k] - temp))
        print L

                  
def transposeMatrix(M):
    T = []
    for i in range(len(M)):
        tRow = []
        for j in range(len(M[i])):
            if i == j:
                tRow.append(M[i][j])
            else:
                tRow.append(M[j][i])
        T.append(tRow)
    return T

def transposeVector(V):
    T = []
    for i in range(len(V)):
        T += [V[i][0]]
    T = [T]
    return T

def getMatrixMinor(M,i,j):
    ret = [row[:j] + row[j+1:] for row in (M[:i]+M[i+1:])]
    return ret

def getDeternminant(M):
    #base case for 2x2 matrix
    if len(M) == 2:
        return M[0][0]*M[1][1]-M[0][1]*M[1][0]

    determinant = 0
    for i in range(len(M)):
        determinant += ((-1)**i)*M[0][i]*getDeternminant(getMatrixMinor(M,0,i))
    return determinant

def getInverse(M):
    determinant = getDeternminant(M)
    #special case for 2x2 matrix:
    if len(M) == 2:
        return [[M[1][1]/determinant, -1*M[0][1]/determinant],
                [-1*M[1][0]/determinant, M[0][0]/determinant]]

    #find matrix of cofactors
    cofactors = []
    for i in range(len(M)):
        cofactorRow = []
        for j in range(len(M)):
            minor = getMatrixMinor(M,i,j)
            cofactorRow.append(((-1)**(i+j)) * getDeternminant(minor))
        cofactors.append(cofactorRow)
    cofactors = transposeMatrix(cofactors)
    
    for i in range(len(cofactors)):
        for j in range(len(cofactors)):
            cofactors[i][j] = cofactors[i][j]/determinant
    return cofactors

def mult(A,B):
	#height: qtd de linhas
	#width: qtd de colunas
	heightA, widthA = len(A), len(A[0])
	heightB, widthB = len(B), len(B[0])

	#verificar dimenses das matrizes
	if (widthA != heightB):
		print "Dimensões das matrizes são incompatíveis"
		return 0

	#inicializar matriz Resultado
	R = [[0.0 for x in range(widthB)] for y in range(heightA)] 

	for i in range(heightA):
		for j in range(widthB):
			for k in range(widthA):
				R[i][j] += A[i][k] * B[k][j]
	return R

def multVector(A,V):
    R = [0.0 for i in range(len(A))]
    for i in range(len(A)):
        for j in range(len(A)):
            R[i] += A[i][j]*V[j]
    return R


def add(A,B):
    R = []
    for i in range(len(A)):
        row = []
        for j in range(len(A[i])):
            row.append(A[i][j] + B[i][j])
        R.append(row)
    return R


def addVector(X,Y):
    R = [0.0 for i in range(len(X))]
    for i in range(len(X)):
        R[i] = X[i] + Y[i]
    

def scalar(a,A):
	ret = [[(a)* x for x in A[i]] for i in range(len(A))]
	return ret

def cout(M):
    for i in range(len(M)):
        print M[i]

def NewtonMethod(n, niter):
    X0 = []
    tol = 10.0**(-5)
    tolk = 0

    X0 = [[1.0] for i in range(n)]

    for k in range (niter):
        J, F, deltaX = [], [], []
            
        for i in range(n):
            F += [[f2(X0,i+1)]]
            J += [df2(X0,i+1)]

        Ji = getInverse(J)

        deltaX = mult(Ji,F)


        X0 = add(X0,scalar(-1,deltaX))

        tolk = 1.0*euclidean(deltaX)/euclidean(X0)
        print ""
        print X0
        if (tolk < tol):
            print "Convergence REACHED"
            return

    print "Convergence not reached"
    

def BroydenMethod(n, niter):
    X0 = []
    x0 = 1.0
    tol = 10.0**(-5)
    tolk = 0.0
    dx = 10.0**(-4)

    #Inicializar matriz Jacobianos e vetor de entrada
    B = [[0.0]*n for i in range(n)]
    X0 = [[x0] for i in range(n)]
    for i in range(n):
        for j in range(n):
            Xt = X0[:]
            Xt[j] = [x0+dx]
            B[i][j] = (f2(Xt,i+1)-f2(X0,i+1))/dx
    #print B

    for k in range (niter):
        Y, F, deltaX = [], [], []
        J = B[:]
        
        for i in range(n):
            F += [[f2(X0,i+1)]]

        #inverter a matriz para resolver o sistema
        Ji = getInverse(J)

        #multiplicar matriz inversa Ji pelo vetor V0
        deltaX = scalar(-1,mult(Ji,F))

        #salvar valor antigo de F
        F0 = F    

        #atualizar X0
        X0 = add(X0,deltaX)

        #calcular F com X0 atualizado
        F = []
        for i in range(n):
            F+= [[f2(X0,i+1)]]

        #calcular Y
        Y = add(F,scalar(-1,F0))

        #calcular o erro e comparar com tolerÃ¢ncia
        tolk = 1.0*euclidean(deltaX)/euclidean(X0)
        #print X0
        #print ""
        
        if (tolk < tol):
            print X0
            print "Convergence REACHED"
            return
        
        else:
            aux = mult(B,deltaX)
            aux = add(Y,scalar(-1.0,aux))
            aux = mult(aux,transposeVector(deltaX))
            aux2 = mult(transposeVector(deltaX), deltaX)
            B = add(B,scalar((1.0/aux2[0][0]),aux))
        
    print "Convergence NOT REACHED"



def NL_MinimumSquare(n, niter):
    B0 = []
    tol = 10.0**(-4)
    tolk = 0
    B0 = [[1.0] for i in range(n)]
    #B0 = [[1.0],[2.0],[3.0]]
    
    for k in range (niter):
        J, F, deltaB = [], [], []

        for i in range(n):
            F += [[f3(B0,i)]]
            J += [df3(B0,i)]

        #calcular deltaB: -inv(J'*J)*J'*F
        a = mult(transposeMatrix(J),J)
        a = getInverse(a)
        b = mult(transposeMatrix(J),F)

        deltaB = scalar(-1,mult(a,b))

        #Calcular B0 como no método de Newton
        B0 = add(B0,deltaB)
    
        tolk = 1.0*euclidean(deltaB)/euclidean(B0)
        #print ""
        print B0

        #return 0
        if (tolk < tol):
            print "Convergence REACHED"
            return

    print "Convergence not reached"


#main()
#NewtonMethod(3,100)
BroydenMethod(3,100)
#NL_MinimumSquare(3,500)
