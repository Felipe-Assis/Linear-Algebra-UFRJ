from math import *

#questao 1; lista 5
def integracao(nfuncao, a, b, N=3, deGauss=True): #funcao que realiza a integracao numerica de uma funcao f(x) no intervalo [a,b]
    if N<1 or N>10: #numero de pontos de integracao deve estar entre 1 e 10
        return "Numero de pontos invalido"
    L=b-a #define L, o comprimento do intervalo de integracao
    if not deGauss: #integracao polinomial
        if N>3: #se N>3, e necessario definir um delta para acrescer aos pontos de integracao
            delta=1.0*L/(N-1)
        P=[]
        if N==1: #caso apenas 1 ponto de integracao, o valor medio da funcao e calculado; regra do Ponto Medio
            P+=[(a+b)/2.0]
            W=[L]
        elif N==2: #regra do Trapezio
            P+=[a,b]
            W=[L/2.0,L/2.0]
        elif N==3: #regra de Simpson
            P+=[a,(a+b)/2.0,b]
            W=[L/6.0,2.0*L/3,L/6.0]
        #calculo das coodernadas e pesos para demais numeros de pontos
        elif N==4:
            for i in range(N):
                P+=[a+i*delta]
            W=[L/8.0,3.0*L/8,3.0*L/8,L/8.0]
        elif N==5:
            for i in range(N):
                P+=[a+i*delta]
            W=[7.0*L/90,16.0*L/45,2.0*L/15,16.0*L/45,7.0*L/90]
        elif N==6:
            for i in range(N):
                P+=[a+i*delta]
            W=[19.0*L/288,25.0*L/96,25.0*L/144,25.0*L/144,25.0*L/96,19.0*L/288]
        elif N==7:
            for i in range(N):
                P+=[a+i*delta]
            W=[41.0*L/840,9.0*L/35,9.0*L/280,34.0*L/105,9.0*L/280,9.0*L/35,41.0*L/840]
        elif N==8:
            for i in range(N):
                P+=[a+i*delta]
            W=[751.0*L/17280,3577.0*L/17280,49.0*L/640,2989.0*L/17280,2989.0*L/17280,49.0*L/640,3577.0*L/17280,751.0*L/17280]
        elif N==9:
            for i in range(N):
                P+=[a+i*delta]
            W=[989.0*L/28350,2944.0*L/14175,-464.0*L/14175,5248.0*L/14175,-454.0*L/2835,5248.0*L/14175,-464.0*L/14175,2944.0*L/14175,989.0*L/28350]
        elif N==10:
            for i in range(N):
                P+=[a+i*delta]
            W=[2857.0*L/89600,15741.0*L/89600,27.0*L/2240,1209.0*L/5600,2889.0*L/44800,2889.0*L/44800,1209.0*L/5600,27.0*L/2240,15741.0*L/89600,2857.0*L/89600]
    #Metodo da Quadratura de Gauss
    #Pesos da Quadratura de Gauss sao tabelados e ja definidos
    else: #http://keisan.casio.com/exec/system/1329114617
        if N==1:
            Z=[0]
            W=[2]
        elif N==2:
            Z=[-0.5773502691896258,0.5773502691896258]
            W=[1,1]
        elif N==3:
            Z=[-0.7745966692414834,0,0.7745966692414834]
            W=[0.555555555555556,0.8888888888888889,0.555555555555556]
        elif N==4:
            Z=[-0.8611363115940526,-0.3399810435848563,0.3399810435848563,0.8611363115940526]
            W=[0.3478548451374539,0.652145154862546,0.652145154862546,0.3478548451374539]
        elif N==5:
            Z=[-0.9061798459386640,-0.5384693101056831,0,0.5384693101056831,0.9061798459386640]
            W=[0.2369268850561891,0.4786286704993665,0.568888888888889,0.4786286704993665,0.2369268850561891]
        elif N==6:
            Z=[-0.9324695142031520,-0.6612093864662645,-0.2386191860831969,0.2386191860831969,0.6612093864662645,0.9324695142031520]
            W=[0.1713244923791703,0.3607615730481386,0.4679139345726910,0.4679139345726910,0.3607615730481386,0.1713244923791703]
        elif N==7:
            Z=[-0.949107912342759,-0.741531185599394,-0.4058451513773972,0,0.4058451513773972,0.741531185599394,0.949107912342759]
            W=[0.1294849661688697,0.279705391489277,0.3818300505051189,0.4179591836734694,0.3818300505051189,0.279705391489277,0.1294849661688697]
        elif N==8:
            Z=[-0.960289856497536,-0.796666477413627,-0.5255324099163290,-0.1834346424956498,0.1834346424956498,0.5255324099163290,0.796666477413627,0.960289856497536]
            W=[0.1012285362903763,0.222381034453374,0.3137066458778873,0.3626837833783620,0.3626837833783620,0.3137066458778873,0.222381034453374,0.1012285362903763]
        elif N==9:
            Z=[-0.968160239507626,-0.836031107326636,-0.6133714327005904,-0.3242534234038089,0,0.3242534234038089,0.6133714327005904,0.836031107326636,0.968160239507626]
            W=[0.08127438836157441,0.1806481606948574,0.2606106964029355,0.312347077040003,0.3302393550012598,0.312347077040003,0.2606106964029355,0.1806481606948574,0.08127438836157441]
        elif N==10:
            Z=[-0.973906528517172,-0.8650633666889845,-0.6794095682990244,-0.4333953941292472,-0.1488743389816312,0.1488743389816312,0.4333953941292472,0.6794095682990244,0.8650633666889845,0.973906528517172]
            W=[0.0666713443086881,0.1494513491505806,0.2190863625159820,0.2692667193099964,0.295524224714753,0.295524224714753,0.2692667193099964,0.2190863625159820,0.1494513491505806,0.0666713443086881]
        P=[] #inicializa lista de pontos
        for i in Z:
            P+=[(a+b+i*L)/2.0] #calcula os pontos em que a funcao sera avaliada para posterior integracao
    integral=0
    for i in range(N):
        integral+=W[i]*funcao(nfuncao,P[i]) #integra a funcao multiplicando a funcao f calculada no ponto P[i] por W[i] (respectivo peso)
    if deGauss:
        integral*=L/2.0 #na quadratura de gauss multiplica-se ainda o termo obtido por L/2.0 para concluir a integracao  
    return integral

def funcao(n, x): #retorna valores das funcoes da lista 5
    if n<1 or n>5:
        return "Funcao desconhecida"
    if n==1: #retorna valor da funcao do exericio 2
        return exp(-x**2/2)/sqrt(2*pi)
    if n==2: #retorna valor da funcao do exericio 3
        omegaN=1.0
        sN=2.0
        csi=0.05
        return sN/((1-(x/omegaN)**2 )**2+(2*csi*x/omegaN)**2)
    if n==2.5: #retorna valor da segunda funcao do exercicio 3
        omegaN=1.0
        sN=2.0
        csi=0.05
        return (x**2)*sN/((1-(x/omegaN)**2)**2+(2*csi*x/omegaN)**2)
    if n==3: #retorna valor da funcao do exercicio 4
        hs=3.0
        tz=5.0
        return (funcao(2,x)/2)*exp(-16*pi**3/((x*tz)**4))*4*(pi**3)*(hs**2)/((x**5)*(tz**4))
    if n==3.5: #retorna valor da segunda funcao do exericio 4
        hs=3.0
        tz=5.0
        return (x**2)*(funcao(2,x)/2)*exp(-16*pi**3/((x*tz)**4))*4*(pi**3)*(hs**2)/((x**5)*(tz**4))
    if n==4: #retorna valor da funcao do exericio 5
        return 2+2*x-x**2+3*x**3
    if n==5: #retorna valor da funcao do exericio 6
        return 1.0/(1+x**2)

def funcao1(x):
    return exp(-x**2/2)/sqrt(2*pi)

def funcao2(x,omegaN=1.0,sN=2.0,csi=0.05):
    return sN/((1-(x/omegaN)**2)**2+(2*csi*x/omegaN)**2)

def funcao2b(x,omegaN=1.0,sN=2.0,csi=0.05):
    return (x**2)*sN/((1-(x/omegaN)**2)**2+(2*csi*x/omegaN)**2)

def funcao3(x,hs=3.0,tz=5.0):
    return (funcao2(x)/2)*exp(-16*pi**3/((x*tz)**4))*4*(pi**3)*(hs**2)/((x**5)*(tz**4))

def funcao3b(x,hs=3.0,tz=5.0):
    return (x**2)*(funcao2(x)/2)*exp(-16*pi**3/((x*tz)**4))*4*(pi**3)*(hs**2)/((x**5)*(tz**4))

def funcao4(x):
    return 2+2*x-x**2+3*x**3

def funcao5(x):
    return 1/(1+x**2)
