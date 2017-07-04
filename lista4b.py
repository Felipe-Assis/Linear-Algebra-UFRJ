from math import *

def f1(V,i): #retorna valores para uma funcao i do sistema do exercicio 3
    x,y,z = V[0][0], V[1][0], V[2][0] #le as variaveis do vetor de entrada
    if (i==1): #calcula valor da funcao 1 do sistema
        resultado = 16.0*(x**4) + 16.0*(y**4) + 1.0*(z**4) - 16.0
        return resultado
    if (i==2): #calcula valor da funcao 2 do sistema
        resultado = 1.0*(x**2) + 1.0*(y**2) + 1.0*(z**2) - 3.0
        return resultado
    if (i==3): #calcula valor da funcao 3 do sistema
        resultado = 1.0*(x**3) - (1.0*y) + 1.0*z - 1.0
        return resultado

def df1(V,i): #retorna o vetor gradiente de uma funcao i do sistema do exercicio 3
    x,y,z = V[0][0], V[1][0], V[2][0] #le as variaveis do vetor de entrada
    if (i==1): #calcula o gradiente da funcao 1
        res1 = 64.0*(x**3)
        res2 = 64.0*(y**3)
        res3 = 4.0*(z**3)
        return [res1, res2, res3]
    if (i==2): #calcula o gradiente da funcao 2
        res1 = 2.0*(x)
        res2 = 2.0*(y)
        res3 = 2.0*z
        return [res1, res2, res3]
    if (i==3): #calcula o gradiente da funcao 3
        res1 = 3.0*(x**2)
        res2 = -1.0
        res3 = 1.0
        return [res1, res2, res3]

def f2(V,i,teta1=0,teta2=11.667): #retorna valores para uma funcao i do sistema do exercicio 4
    c2, c3, c4 = V[0][0], V[1][0], V[2][0] #le as variaveis do vetor de entrada
    if (i==1): #calcula valor da funcao 1 do sistema
	    ret = 2*(c3**2) + c2**2 + 6*(c4**2) - 1
	    return ret
    if (i==2): #calcula valor da funcao 2 do sistema
            ret = 8*(c3**3) + 6*c3*(c2**2) + 36*c3*c2*c4 + 108*c3*(c4**2) - teta1
            return ret
    if (i==3): #calcula valor da funcao 3 do sistema
            ret = 60*(c3**4) + 60*(c3**2)*(c2**2) + 576*(c3**2)*c2*c4 + 2232*(c3**2)*(c4**2) + 252*(c4**2)*(c2**2) + 1296*(c4**3)*c2 + 3348*(c4**4) + 24*(c2**3)*c4 + 3*c2 - teta2
            return ret

def df2(V,i): #retorna o vetor gradiente de uma funcao i do sistema do exercicio 4
    c2, c3, c4 = V[0][0], V[1][0], V[2][0] #le as variaveis do vetor de entrada
    if (i==1): #calcula o gradiente da funcao 1
        res1 = 2*c2
        res2 = 4*c3
        res3 = 12*c4
        return [res1, res2, res3]
    if (i==2): #calcula o gradiente da funcao 2
        res1 = 12*c3*c2 + 36*c3*c4
        res2 = 24*(c3**2) + 6*(c2**2) + 36*c2*c4 + 108*(c4**2)
        res3 = 36*c3*c2 + 216*c3*c4
        return [res1, res2, res3]
    if (i==3): #calcula o gradiente da funcao 3
        res1 = 120*(c3**2)*c2 + 576*(c3**2)*c4 + 504*(c4**2)*c2 + 1296*(c4**3)
        + 72*(c2**2)*c4 + 3
        res2 = 240*(c3**3) + 120*c3*(c2**2) + 1152*c3*c2*c4 + 4464*c3*(c4**2)
        res3 = 576*(c3**2)*c2 + 4464*(c3**2)*c4 + 504*c4*(c2**2) + 3888*(c4**2)*c2 + 13392*(c4**3) + 24*(c2**3)
        return [res1, res2, res3]

def f3(V,i): #retorna valores para a funcao do sistema do exercicio 5 com coeficientes V=[b0,b1,b2]
    b0,b1,b2 = V[0][0], V[1][0], V[2][0] #le as variaveis do vetor de entrada
    x = [1.0,2.0,3.0]
    y = [1.0,2.0,9.0]
    ret = b0 + b1*(x[i-1]**b2) - y[i-1]
    return ret

def df3(V,i): #retorna vetor com derivadas parciais da funcao do exercicio 5 em relacao a b0, b1 e b2 
    b0,b1,b2 = V[0][0], V[1][0], V[2][0] #le as variaveis do vetor de entrada
    x = [1.0,2.0,3.0]
    ret1 = 1.0
    ret2 = x[i-1]**b2
    ret3 = (b1*(x[i-1]**b2)*log(x[i-1]))
    return [ret1, ret2, ret3]

def funcao(nfuncao, V, i): #decodifica quando funcao deve ser chamada
    if nfuncao<1 or nfuncao>3: #verifica se nfuncao identifica uma funcao valida
        return "Numero da funcao invalido"
    if (nfuncao==1): #efetua chamada para o exercicio 3
        return f1(V, i)
    if (nfuncao==2): #efetua chamada para o exercicio 4
        return f2(V, i)
    if (nfuncao==3): #efetua chamada para o exercicio 5
        return f3(V, i)

def dfuncao(nfuncao, V, i): #decodifica quando funcao deve ser chamada para calcular derivadas
    if nfuncao<1 or nfuncao>3: #verifica se nfuncao identifica uma funcao valida
        return "Numero da funcao invalido"
    if (nfuncao==1): #efetua chamada para o exercicio 3
        return df1(V, i)
    if (nfuncao==2): #efetua chamada para o exercicio 4
        return df2(V, i)
    if (nfuncao==3): #efetua chamada para o exercicio 5
        return df3(V, i)

def euclidean(X): #retorna a norma euclidiana de um vetor X
    total = 0;
    for x in X:
        total += 1.0*x[0]**2
    return sqrt(total)
                  
def transpose(M): #retorna a transposta de uma matriz (ou de um vetor) M
    return [list(i) for i in zip(*M)] #unzipping para uma lista

def getMatrixMinor(M,i,j): #retorna submatriz de M retirando a linha i e a coluna j
    ret = [row[:j] + row[j+1:] for row in (M[:i]+M[i+1:])]
    return ret

def getDeterminant(M): #retorna o determinante da matriz M
    if len(M)==1: #se matriz 1x1, retorna o valor de sua unica posicao
        return M[0][0]
    #caso especial para matrizes 2x2
    if len(M) == 2:
        return M[0][0]*M[1][1]-M[0][1]*M[1][0]
    #
    determinant = 0
    for i in range(len(M)): #para cada linha da matriz
        determinant += ((-1)**i)*M[0][i]*getDeterminant(getMatrixMinor(M,0,i))
    return determinant

def getInverse(M): #retorna a inversa de uma matriz M
    determinant = getDeterminant(M)
    if determinant==0: #se determinante for 0, matriz nao invertivel
        print "Matrix cannot be inversed. Its determinant is 0."
        return "Matrix cannot be inversed. Its determinant is 0."
    if len(M)==1: #se matriz 1x1, inversa e o inverso de seu unico elemento
        return [[1/M[0][0]]]
    #caso especial para matrizes 2x2
    if len(M) == 2:
        return [[M[1][1]/determinant, -1*M[0][1]/determinant],
                [-1*M[1][0]/determinant, M[0][0]/determinant]]
    #encontrar matriz de cofatores
    cofactors = []
    for i in range(len(M)): #para cada linha i de M
        cofactorRow = []
        for j in range(len(M)): #para cada coluna j de M
            minor = getMatrixMinor(M,i,j) #obtem submatriz excluindo linha i e coluna j
            cofactorRow+=[((-1)**(i+j)) * getDeterminant(minor)]
        cofactors+=[cofactorRow] #adiciona a matriz de cofatores
    cofactors = transpose(cofactors) #obtem transposta da matriz de cofatores
    for i in range(len(cofactors)): #para cada linha i da matriz de cofatores
        for j in range(len(cofactors)): #para cada linha i da matriz de cofatores
            cofactors[i][j] = cofactors[i][j]/determinant
    return cofactors

def mult(A,B): #retorna matriz resultante da multiplicacao da matriz A pela matriz B
    #height: qtd de linhas
    #width: qtd de colunas
    heightA, widthA = len(A), len(A[0])
    heightB, widthB = len(B), len(B[0])
    #verificar dimensoes das matrizes
    if (widthA != heightB):
        return "Dimensoes das matrizes sao incompativeis"
    #inicializar matriz Resultado
    R = [[0.0 for x in range(widthB)] for y in range(heightA)] 
    #
    for i in range(heightA): #para cada linha de A
        for j in range(widthB): # para cada coluna de B
            for k in range(widthA): # para cada elemento da linha i de A e da coluna j de B
                R[i][j] += A[i][k] * B[k][j] #produto elemento a elemento
    return R

def add(A,B): #retorna a matriz resultante da soma entre matirzes A e B
    R = []
    for i in range(len(A)): #para cada linha
        row = []
        for j in range(len(A[i])): #para cada coluna
            row+=[A[i][j] + B[i][j]] #soma elemento a elemento
        R+=[row]
    return R

def scalar(a,A): #retorna a matriz resultante do produto entre um escalar a e uma matriz A
    ret = [[(a)* x for x in A[i]] for i in range(len(A))] #produto para cada posicao de cada linha de A
    return ret

def NewtonMethod(nfuncao, X0=[[1.0],[1.0],[1.0]], n=3, niter=1000, tol=10.0**(-4)): #X0=[[-1.0],[0.0],[1.0]] #executa o metodo de Newton para solucao do sistema
    tolk = 0.0
    for k in range (niter): #por no maximo niter iteracoes
        J, F, deltaX = [], [], []
        for i in range(n):
            F += [[funcao(nfuncao,X0,i+1)]] #prepara vetor F aplicando X0
            J += [dfuncao(nfuncao,X0,i+1)] #prepara Jacobiana aplicando X0
        Ji = getInverse(J)
        deltaX = scalar(-1,mult(Ji,F))
        X0 = add(X0,deltaX)
        tolk = 1.0*euclidean(deltaX)/euclidean(X0) #obtem variacao relativa entre passos
        if (tolk < tol): #convergencia alcancada
            return X0
    return "Convergencia NAO ALCANCADA"

def BroydenMethod(nfuncao, X0=[[1.0],[1.0],[1.0]], calculaJacobiana=True, B=[[0,-3,-1],[3,0,5],[1,5,0]], niter=1000, tol=10.0**(-4)): #aplica metodo de Broyden para solucao do sistema 
    tolk = 0.0
    dx = 10.0**(-4)
    n=len(B)
    if calculaJacobiana: #se for escolhido calcular a Jacobiana automaticamente
        #Inicializar matriz Jacobiana
        B = [[0.0]*n for i in range(n)]
        for i in range(n): #para cada funcao
            for j in range(n): #para cada variavel
                Xt = X0[:]
                Xt[j] = [X0[j][0]+dx]
                B[i][j] = (funcao(nfuncao,Xt,i+1)-funcao(nfuncao,X0,i+1))/dx #obtem derivada numericamente
    for k in range (niter): #por no maximo niter iteracoes
        Y, F, deltaX = [], [], []
        J = B
        for i in range(n):
            F += [[funcao(nfuncao,X0,i+1)]] #calcula vetor F aplicando X0
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
            F += [[funcao(nfuncao,X0,i+1)]]
        #calcular Y
        Y = add(F,scalar(-1,F0))
        #calcular o erro e comparar com tolerancia
        tolk = 1.0*euclidean(deltaX)/euclidean(X0)
        if (tolk < tol): #convergencia alcancada
            return X0
        else: #calcula novo B
            aux = mult(B,deltaX)
            aux = add(Y,scalar(-1.0,aux))
            aux = mult(aux,transpose(deltaX))
            aux2 = mult(transpose(deltaX), deltaX)
            B = add(B,scalar(getInverse(aux2)[0][0], aux))
    return "Convergencia NAO ALCANCADA"

def NL_MinimumSquares(nfuncao, B0=[[3.0],[2.0],[1.0]], n=3, niter=1000, tol=10.0**(-4)): #B0=[[[10.0],[10.0],[10.0]] #calcula coeficientes B que melhor ajustam uma funcao a um conjunto de pontos
    tolk = 0.0
    #B0 = [[1.0],[2.0],[3.0]]
    for k in range (niter): #por no maximo niter iteracoes
        J, F, deltaB = [], [], []
        for i in range(n):
            F += [[funcao(nfuncao,B0,i+1)]] #calcula F aplicando B0
            J += [dfuncao(nfuncao,B0,i+1)] #calcula Jacobiana aplicando B0
        #calcular deltaB: -inv(J'*J)*J'*F
        a = mult(transpose(J),J)
        a = getInverse(a)
        b = mult(transpose(J),F)
        deltaB = scalar(-1,mult(a,b))
        #Calcular B0 como no metodo de Newton
        B0 = add(B0,deltaB)
        tolk = 1.0*euclidean(deltaB)/euclidean(B0) #obtem variacao relativa entre passos
        if (tolk < tol): #convergencia alcancada
            return B0
    return "Convergencia NAO ALCANCADA"
