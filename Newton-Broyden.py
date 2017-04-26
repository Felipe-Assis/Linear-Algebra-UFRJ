from math import *

def f1(V,i):
    if (i==1):
        resultado = 16.0*(V[0]**4) + 16.0*(V[1]**4) + 1.0*V[2]**4 - 16.0
        return resultado
    
    if (i==2):
        resultado = (1.0*V[0]**2) + (1.0*V[1]**2) + 1.0*V[2]**2 - 3.0
        return resultado
    
    if (i==3):
        resultado = (1.0*V[0]**3) - (1.0*V[1]) + 1.0*V[2] - 1.0
        return resultado


def df1(V,i):
    if (i==1):
        res1 = 64.0*(V[0]**3)
        res2 = 64.0*(V[1]**3)
        res3 = 4.0*V[2]**3
        return [res1, res2, res3]

    if (i==2):
        res1 = 2.0*(V[0])
        res2 = 2.0*(V[1])
        res3 = 2.0*V[2]
        return [res1, res2, res3]

    if (i==3):
        res1 = 3.0*(V[0]**2)
        res2 = -1.0
        res3 = 1.0
        return [res1, res2, res3]


def f2(V,i):
	#initialize variables
	c2, c3, c4 = V[0], V[1], V[2]
	teta1, teta2 = 0, 0

	if (i==1):
		ret = 2*(c3**2) + c2**2 + 6*(c4**2) - 1
		return ret
	if (i==2):
		ret = 8*(c3**3) + 6*c3*(c2**2) + 36*c3*c2*c4 + 108*c3*(c4**2) - teta1
		return ret
	if (i==3):
		ret = 60*(c3**4) + 60*(c3**2)*(c2**2) + 576*(c3**2)*c2*c4 + 2232*(c3**2)*(c4**2)
		+ 252*(c4**2)*(c2**2) + 1296*(c4**3)*c2 + 3348*(c4**4) + 24*(c2**3) + 3*c2 - teta2
		return ret

def euclidean(X):
    total = 0;
    for x in X:
        total += 1.0*x**2
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

                  
def transposeMatrix(m):
    t = []
    for r in range(len(m)):
        tRow = []
        for c in range(len(m[r])):
            if c == r:
                tRow.append(m[r][c])
            else:
                tRow.append(m[c][r])
        t.append(tRow)
    return t

def getMatrixMinor(m,i,j):
    ret = [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]
return ret

def getMatrixDeternminant(m):
    #base case for 2x2 matrix
    if len(m) == 2:
        return m[0][0]*m[1][1]-m[0][1]*m[1][0]

    determinant = 0
    for c in range(len(m)):
        determinant += ((-1)**c)*m[0][c]*getMatrixDeternminant(getMatrixMinor(m,0,c))
    return determinant

def getMatrixInverse(m):
    determinant = getMatrixDeternminant(m)
    #special case for 2x2 matrix:
    if len(m) == 2:
        return [[m[1][1]/determinant, -1*m[0][1]/determinant],
                [-1*m[1][0]/determinant, m[0][0]/determinant]]

    #find matrix of cofactors
    cofactors = []
    for r in range(len(m)):
        cofactorRow = []
        for c in range(len(m)):
            minor = getMatrixMinor(m,r,c)
            cofactorRow.append(((-1)**(r+c)) * getMatrixDeternminant(minor))
        cofactors.append(cofactorRow)
    cofactors = transposeMatrix(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c] = cofactors[r][c]/determinant
    return cofactors

def mult(A,B):
	#height: qtd de linhas
	#width: qtd de colunas
	heightA, widthA = len(A), len(A[0])
	heightB, widthB = len(B), len(B[0])

	#verificar dimensões das matrizes
	if (widthA != heightB):
		print "Dimensões das matrizes não são incompatíveis"
		return 0

	#inicializar matriz Resultado
	R = [[0.0 for x in range(heightA)] for y in range(widthB)] 

	for i in range(heightA):
		for j in range(widthB):
			for k in range(heightB):
				R[i][j] += A[i][k] * B[k][j]
	return R

def add(A,B):
	#height: qtd de linhas
	#width: qtd de colunas
	heightA, widthA = len(A), len(A[0])
	heightB, widthB = len(B), len(B[0])

	#verificar dimensões das matrizes
	if ((widthA != widthB) or (heightA != heightB)):
		print "Dimensões das matrizes não são incompatíveis"
		return 0

	R = [[0.0 for x in range(heightA)] for y in range(widthA)]
	for i in range(heightA):
		for j in range(widthA):
			R [i][j] = A[i][j] + B[i][j]
	return ret

def scalar(a,A):
	ret = [[(a)* x for x in A[i]] for i in range(len(A))]
	return ret

def NewtonMethod(n, niter):
    X0 = []
    tol = 10.0**(-5)
    tolk = 0


    X0 = [1 for i in range(n)]

    for k in range (niter):
        J, F, deltaX = [], [], []
        
        for i in range(n):
            F += [f(X0,i+1)]
            J += [df(X0,i+1)]

        Ji = getMatrixInverse(J)

        for i in range(n):
            temp = 0
            for j in range(n):
                temp += Ji[i][j]*F[j]
            deltaX += [round(temp,5)]


        for i in range(n):
            X0[i] = X0[i]-deltaX[i]

        tolk = 1.0*euclidean(deltaX)/euclidean(X0)
        print X0
        print "============================"
        if (tolk < tol):
            print "Convergence REACHED"
            print X0
            return

    print "Convergence not reached"

def BroydenMethod(n, niter):
    X0 = []
    tol = 10.0**(-5)
    tolk = 0

    #Inicializar matriz Jacobianos e vetor de entrada
    B = [[1.0]*n for i in range(n)]
    X0 = [1.0 for i in range(n)]

    for k in range (niter):
        F, deltaX = [], []
        J = B
        
        for i in range(n):
            F += [f2(X0,i+1)]

        #inverter a matriz para resolver o sistema
        Ji = getMatrixInverse(J)

        #multiplicar matriz inversa Ji pelo vetor V0
        for i in range(n):
            temp = 0.0
            for j in range(n):
                temp += Ji[i][j]*F[j]
            deltaX += [round(temp,5)]

        #salvar valor antigo de F
        F0 = F    

        #atualizar X0
        for i in range(n):
            X0[i] = X0[i]-deltaX[i]

        #calcular F com V0 atualizado
        F = []
        for i in range(n):
        	F+= [f2(X0,i+1)]

        #calcular Y
        for i in range(n):
        	Y += [F[i]-F0[i]]

        #calcular o erro e comparar com tolerância
        tolk = 1.0*euclidean(deltaX)/euclidean(X0)
        print X0
        print "============================"
        if (tolk < tol):
            print "Convergence REACHED"
            print X0
            return
        else:
            aux = mult(B,transposeMatrix(deltaX))
            aux = add(Y,scalar((-1),B))
            aux2 = mult(deltaX, transposeMatrix(deltaX))
            B = add(B,(aux/aux2)
                    
    print "Convergence not reached"

#main()
BroydenMethod(3,50)
