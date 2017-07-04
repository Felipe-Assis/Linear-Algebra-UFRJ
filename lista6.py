# "conventional comment"
## code regarding .csv writing
### code regarding printing and writing only some points from all points iterated

from math import *
##import csv

def f1(x, y): #retorna valores para a funcao descrita no exercicio 1
    return -2*x*(y**2)

def f2(x, y, y_1, m=1.0, c=0.2,k=1, w=0.5): #retorna valores para a funcao descrita no exercicio 2
    return (2*sin(w*x)+sin(2*w*x)+cos(3*w*x)-k*y-c*y_1)/m

def f3(x, y, y_1, kd=1, g=9.81): #retorna valores para a funcao descrita no exercicio 3
    return -g-kd*y_1*abs(y_1)

def EulerEDO(h, x=0, y=1, lim_sup=2): #retorna pontos da solucao para uma EDO de primeira ordem usando o metodo de Euler
    ##fd=open('EulerEDO.csv', 'wb')
    ##wr=csv.writer(fd)
    ##wr.writerow([x,y])
    #pontos=[[x,y]]
    print [x,y]
    while x+h-0.001<=lim_sup:
        y+=h*f1(x,y)
        x+=h
        #pontos+=[[x,y]]
        ##wr.writerow([x,y])
        print [x,y]
    ##fd.close()
    #return pontos

def RK2EDO(h, x=0, y=1, lim_sup=2): #retorna pontos da solucao para uma EDO de primeira ordem usando o metodo de Runge-Kutta de segunda ordem
    ##fd=open('RK2EDO.csv', 'wb')
    ##wr=csv.writer(fd)
    ##wr.writerow([x,y])
    #pontos=[[x,y]]
    print [x,y]
    while x+h-0.001<=lim_sup:
        K1=f1(x,y)
        K2=f1(x+h,y+h*K1)
        y+=(h/2.0)*(K1+K2)
        x+=h
        #pontos+=[[x,y]]
        ##wr.writerow([x,y])
        print [x,y]
    ##fd.close()
    #return pontos

def RK4EDO(h, x=0, y=1, lim_sup=2): #retorna pontos da solucao para uma EDO de primeira ordem usando o metodo de Runge-Kutta de quarta ordem
    ##fd=open('RK4EDO.csv', 'wb')
    ##wr=csv.writer(fd)
    ##wr.writerow([x,y])
    #pontos=[[x,y]]
    print [x,y]
    while x+h-0.001<=lim_sup:
        K1=f1(x,y)
        K2=f1(x+(h/2),y+(h/2)*K1)
        K3=f1(x+(h/2),y+(h/2)*K2)
        K4=f1(x+h,y+h*K3)
        y+=(h/6.0)*(K1+2*K2+2*K3+K4)
        x+=h
        #pontos+=[[x,y]]
        ##wr.writerow([x,y])
        print [x,y]
    ##fd.close()
    #return pontos

def TaylorEDO2(h, x=0, y=0, y_1=0, lim_sup=100): #retorna pontos da solucao para uma EDO de segunda ordem usando o metodo de expansao em series de Taylor
    ##fd=open('TaylorEDO2.csv', 'wb')
    ##wr=csv.writer(fd)
    ##wr.writerow([x,y])
    #pontos=[[x,y]]
    print [x,y]
    while x+h-0.001<=lim_sup:
        y_2=f2(x,y,y_1) #f3
        y+=y_1*h+(y_2/2.0)*(h**2)
        y_1+=y_2*h
        x+=h
        #pontos+=[[x,y]]
        
        ###if round(x,1).is_integer(): ###prints integers
        ###    wr.writerow([x,y])
        ###    print [x,y]
        
        ###if round(x,1).is_integer() or (round(x,1)+0.5).is_integer() or (round(x,1)-0.5).is_integer(): ###prints integers and ending with .5
        ###    wr.writerow([x,y])
        ###    print [x,y]
        
        ##wr.writerow([x,y])
        print [x,y]
    ##fd.close()
    #return pontos

def RKNEDO2(h, x=0, y=0, y_1=0, lim_sup=100): #retorna pontos da solucao para uma EDO de primeira ordem usando o metodo de Runge-Kutta-Nystrom
    ##fd=open('RKNEDO2.csv', 'wb')
    ##wr=csv.writer(fd)
    ##wr.writerow([x,y])
    #pontos=[[x,y]]
    print [x,y]
    while x+h-0.001<=lim_sup:
        K1=(h/2.0)*f2(x,y,y_1) #f3
        Q=(h/2.0)*(y_1+K1/2.0)
        K2=(h/2.0)*f2(x+h/2.0,y+Q,y_1+K1) #f3
        K3=(h/2.0)*f2(x+h/2.0,y+Q,y_1+K2) #f3
        L=h*(y_1+K3)
        K4=(h/2.0)*f2(x+h,y+L,y_1+2*K3) #f3
        y+=h*(y_1+(K1+K2+K3)/3.0)
        y_1+=(K1+2*K2+2*K3+K4)/3.0
        x+=h
        #pontos+=[[x,y]]
        
        ###if round(x,1).is_integer(): ###prints integers
        ###    wr.writerow([x,y])
        ###    print [x,y]
        
        ###if round(x,1).is_integer() or (round(x,1)+0.5).is_integer() or (round(x,1)-0.5).is_integer():  ###prints integers and ending with .5
        ###    wr.writerow([x,y])
        ###    print [x,y]
        
        ##wr.writerow([x,y])
        print [x,y]
    ##fd.close()
    #return pontos
