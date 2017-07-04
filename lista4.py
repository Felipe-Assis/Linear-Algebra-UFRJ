from math import *

def bissecao(nFuncao, a, b, tol=0.0001): #aplica metodo da bissecao para resolver equacao nao-linear
    """Lista 4: (1,250,300) (1,-300,-250) (2,-10,-5) (2,0,2)"""
    if a==b: #caso a=b 
        if funcao(nFuncao, a-tol)<0<funcao(nFuncao, a+tol):
            return a #raiz e o proprio 'a'; verificado que a funcao corta o eixo x 
        else:
            return "Parametro nao e raiz. Insira dois parametros diferentes" #funcao nao corta o eixo x
    #verificar se f(a) > f(b) e trocar os valores de a e b caso verdadeiro
    if funcao(nFuncao, a)>funcao(nFuncao, b): 
        a,b=b,a 
    while abs(b-a)>=tol: #verificar se |a-b| e maior ou igual a tolerancia
        x=(a+b)/2.0 #ponto x em que a funcao sera calculada na proxima etapa
        f=funcao(nFuncao, x) #calcula a funcao no ponto x
        if f>0: #se f(x)>0, b recebe valor de x; isso se deve a f(b) ser maior que 0
            b=x
        else:
            a=x #se f(x)<0, a recebe valor de x
    return x

def newtonRaiz(nFuncao, x0, secante=False, ITERMAX=1000000, tol=0.0001): #aplica metodo de newton para resolver equacao NL
    try:
        if not secante: #verifica se flag de secante = False
            for i in range(ITERMAX):
                x=x0-funcao(nFuncao, x0)/dfuncao(nFuncao, x0) #x = x0 -  f(x)/f'(x)
                if abs(x-x0)<tol:
                    return x
                x0=x
        else:
            #derivada calculada numericamente por diferencas finitas
            delta=0.001 #define um acrescimo delta
            x=x0+delta #x* = x + h
            fa=funcao(nFuncao, x0)
            for i in range(ITERMAX):
                fi=funcao(nFuncao, x) #calcular f em x
                xprox=x-fi*(x-x0)/(fi-fa) #xk+1 = xk - f(xk)/df(xk);
                if abs(xprox-x)<tol: #se diferenca entre xk+1 - xk menor que tolerancia
                    return xprox #retornar xprox; solucao
                fa=fi
                x0=x
                x=xprox #atualizar x com o xprox obtido anteriormente e repetir o processo
        return "Convergencia nao foi alcancada"
    except:
        return "Ponto inadequado"

def interpolacaoInversa(nFuncao, X, ITERMAX=1000000, tol=0.0001): #constituir um polinomio quadratico
    if len(X)<3: #numero de pontos
        return "Insira mais pontos"
    if len(set(X))<len(X): #set(X) retornar lista ordenada sem elementos repetidos; #len(set(X)) retornar tamanho
        return "Todos os pontos devem ser distintos entre si" #se len(set(X)) < len(X) significa que os pontos nao sao todos diferentes 
    X.sort()
    x0=None
    for i in range(ITERMAX):
        Y=[] #inicializa o vetor Y
        for i in X:
            Y+=[funcao(nFuncao, i)] #preenche Y com a funcao calculada em cada elemento i de X
        x=0 #inicializa variavel que recebera o produtorio
        for i in range(len(X)):
            numerador=1
            denominador=1
            for k in range(len(Y)):
                if k!=i:
                    numerador*=Y[k] #calcula o produto entre todos Y[k] exceto Y[i]
                    denominador*=Y[i]-Y[k]  #calcula o produto entre todos (Y[i]-Y[k])
            x+=X[i]*numerador/denominador  #calcula numerador/denominador e multiplica pelo respectivo X[i] somando ao x que será a solução
        if x0!=None and abs(x-x0)<tol: #verificar diferenca entre valor absoluto de x-x0
            return x
        x0=x
        if abs(max(Y))>abs(min(Y)): #verifica se valor absoluto do maior elemento de Y é maior do que o menor
            i=Y.index(max(Y)) #retorna indice do maior elemento de Y
        else:
            i=Y.index(min(Y)) #retorna indice do menor elemento de Y
        X[i]=x #atualiza o x na posicao do xi equivalente ao fi de Y
        X.sort()
    return "Convergencia nao foi alcancada"

def funcao(n, x): #retornar valores para uma funcao do exercicio 1 ou 2
    if n<1 or n>2:
        return "Funcao desconhecida"
    if n==1:
        g=9.806
        k=0.00341
        return log(cosh(x*sqrt(g*k)))-50
    if n==2:
        return 4*cos(x)-exp(2*x)

def dfuncao(n, x): #retornar valores da derivada para uma funcao do exercicio 1 ou 2
    if n<1 or n>2:
        return "Funcao desconhecida"
    if n==1:
        g=9.806
        k=0.00341
        return sqrt(g*k)*tanh(x*sqrt(g*k))
    if n==2:
        return -2*(exp(2*x)+2*sin(x))

def funcao1(x, g=9.806, k=0.00341):
    return log(cosh(x*sqrt(g*k)))-50

def funcao2(x):
    return 4*cos(x)-exp(2*x)

def dfuncao1(x, g=9.806, k=0.00341):
    return sqrt(g*k)*tanh(x*sqrt(g*k))

def dfuncao2(x):
    return -2*(exp(2*x)+2*sin(x))
