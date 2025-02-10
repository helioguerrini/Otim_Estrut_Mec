# -*- coding: utf-8 -*-
"""
Spyder Editor
Criado por:
    Hélio Guerrini Filho
"""
########### PROGRAMA PARA MINIMIZAR MASSA DE TREM DE ENGRENAGENS ###########
from numpy import array
from scipy import optimize
#from scipy.optimize import minimize
from math import pi
import time
###################### DADOS DE ENTRADA ###########################
# Potencia
Pot = 800000 # W (Watts)
# Rotacao
n = 1600. # rpm
# Coeficiente elastico
ff = 2.86e11 # Pa
# Vida em horas
h = 15000.
# Dureza Brinell
HB = 6270. # MPa
# Relacao de transmissao total
iT = 8
# Parametro de Penalidade
plam = 20.
# Tensão admissivel do aco
sigmaadm = 200. # MPa
# Fator de serviço
fie = 1.0

###################### TORQUE DE ENTRADA ##########################
MT1 = Pot/(2*(pi)*(n/60))

'''
####################### PROJETO INICIAL ###########################
z1, z3, i1 = 19.0, 19.0, 6.0
i2 = iT/i1
m1, m3, b1, b3 = 0.00275 , 0.004, 0.012, 0.075
'''
####################### PROJETO INICIAL ###########################
z1, z3, i1 = 46.0, 44.0, 1.6
i2 = iT/i1
m1, m3, b1, b3 = 0.006, 0.007, 0.082, 0.075

# Volume inicial dos pares engrenados
print ('Volume inicial:')
Vi = (pi/4)*(b1*(m1*z1)**2 + b1*(m1*z1*i1)**2 + b3*(m3*z3)**2 + b3*(m3*z3*i2)**2)
print (Vi)
#################### VETOR DE PARAMETROS ##########################
param = (MT1, n, ff, h, HB, iT, plam, Vi, sigmaadm,fie)

################### VETOR DE VARIÁVEIS DE PROJETO #################
dStart = array([z1, z3, i1, i2, m1, m3, b1, b3])

######################### RESTRICOES ##############################
# Restrição da relacao de transmissao total
def h1(i1,i2,iT):
    restrh1 = (((i1*i2)/iT) - 1)*10
    return restrh1

########## Restricao da pressao entre os flancos dos dentes das engrenagens ###########
def g1(ne, he, HBe, i1e, MT1e, ffe, b1e, m1e, z1e, fie):
    # Fator de durabilidade
    w = (60*ne*he)/(10**6)
    # Pressao admissivel
    Padm = (0.487*HBe)/(w**(1/6))
    # Pressao de contato
    Pc = (((fie*2*ffe*MT1e*(i1e+1))/((i1e+0.14)*b1e*(m1e*z1e)**2))**0.5)*1.05
    restrg1 = ((Pc*10**(-6))/Padm) - 1
    return restrg1

############### Restricao de tensao no pe do dente das engrenagens #################
def g2(MT1, z1, m1, b1, sigmaadm):
    ########## Forca tangencial ##########
    WT = 2*MT1/(m1*z1)
    ########## Fator de forma ##########
    if z1 < 17:
        q = 1000
    if z1 >= 17:
        qq1 = 3.6 + (-1.00000000e-01)*(z1-17) + (8.33333333e-03)*(z1-17)*(z1-18)
        qq2 = (-3.96825397e-04)*(z1-17)*(z1-18)*(z1-21)
        q = qq1 +qq2
        if z1 >= 24:
            q = 0
            qq3 = 3.2 + (-2.50000000e-02)*(z1-24) + (8.33333333e-04)*(z1-24)*(z1-28)
            qq4 = (-5.20833333e-05)*(z1-24)*(z1-28)*(z1-34)
            q = qq3 + qq4
            if z1 >= 40:
                q = 0
                q = 2.9 + (-1.00000000e-02)*(z1-40) + (1.33333333e-04)*(z1-40)*(z1-50)
                if z1 >= 65:
                    q = 0
                    q = 2.7 + (-6.66666667e-03)*(z1-65) + (1.90476190e-04)*(z1-65)*(z1-80)
                    if z1 > 100:
                        q = 2.45
    # Tensão maxima                   
    sigmamax = (q*WT)/(b1*m1)
    restrg2 = (sigmamax*10**(-6))/(sigmaadm) - 1
    return restrg2

########### Restrições de limites de variáveis ############
def g3(z1,m1,b1,z3,m3,b3,i1):
    restrgA = (z1/25)-1
    restrgB = (z3/25)-1
    restrgC = (m1/0.0003) - 1
    restrgD = (m3/0.0003) - 1
    restrgE = (m1/0.075) - 1
    restrgF = (m3/0.075) - 1
    restrgG = (b1/0.005) - 1
    restrgH = (b3/0.005) - 1
    restrgI = (b1/(1.2*z1*m1)) - 1
    restrgJ = (b3/(1.2*z3*m3)) - 1
    
    restrgK = (m3/m1) - 1
   
    restrgM = (b3/b1) - 1
    
    listarestr = [restrgA,restrgB,restrgC,restrgD,restrgE,restrgF,restrgG,restrgH,restrgI,restrgJ,restrgK,restrgM]
    return listarestr
###################### FUNÇÃO OBJETIVO A SER MINIMIZADA ####################
def F(x, *param):
    z1, z3, i1, i2 = x[0], x[1], x[2], x[3] 
    m1, m3, b1, b3 = x[4], x[5], x[6], x[7]
    MT1, n, ff, h, HB, iT, plam, Vi, sigmaadm, fie = param
    #print (n)

    # Função de penalidade
    Pn1 = (h1(i1,i2,iT))**2
    Pn2 = (max(0,g1(n, h, HB, i1, MT1, ff, b1, m1, z1, fie)))**2
    Pn3 = (max(0,g2(MT1, z1, m1, b1, sigmaadm)))**2
    n2, MT2 = n/i1, i1*MT1
    Pn4 = (max(0,g1(n2, h, HB, i2, MT2, ff, b3, m3, z3, fie)))**2
    Pn5 = (max(0,g2(MT2, z3, m3, b3, sigmaadm)))**2
    
    listarestricao = g3(z1,m1,b1,z3,m3,b3,i1)
    Pn6 = (min(0,listarestricao[0]))**2
    Pn7 = (min(0,listarestricao[1]))**2
    Pn8 = (min(0,listarestricao[2]))**2
    Pn9 = (min(0,listarestricao[3]))**2
    Pn10 = (max(0,listarestricao[4]))**2
    Pn11 = (max(0,listarestricao[5]))**2
    Pn12 = (min(0,listarestricao[6]))**2
    Pn13 = (min(0,listarestricao[7]))**2
    Pn14 = (max(0,listarestricao[8]))**2
    Pn15 = (max(0,listarestricao[9]))**2
    Pn16 = (max(0,listarestricao[10]))**2
    Pn17 = (max(0,listarestricao[11]))**2
    
    
    Pn = Pn1 + Pn2 + Pn3 + Pn4 + Pn5 + Pn6 + Pn7 + Pn8 + Pn9 + Pn10 + Pn11 + Pn12 + Pn13 + Pn14 + Pn15 + Pn16 + Pn17
    fobj = (1/Vi)*(pi/4)*(b1*(m1*z1)**2 + b1*(m1*z1*i1)**2 + b3*(m3*z3)**2 + b3*(m3*z3*i2)**2)
    Pfobj = fobj + plam*Pn
    #print (fobj)


    return Pfobj

####################### PROJETO FINAL #############################
res1 = optimize.fmin(F, dStart, args=param, maxiter= 3000)
#res = minimize(F, dStart, args=param)
#print (res)

print ('PROJETO FINAL:')
print ('z1 = ' + str(res1[0]))
print ('z3 = ' + str(res1[1]))
print ('i1 = ' + str(res1[2]))
print ('i2 = ' + str(res1[3]))
print ('m1 = ' + str(res1[4]))
print ('m3 = ' + str(res1[5]))
print ('b1 = ' + str(res1[6]))
print ('b3 = ' + str(res1[7]))

Vfinall = (pi/4)*(res1[6]*(res1[4]*res1[0])**2 + res1[6]*(res1[4]*res1[0]*res1[2])**2 + res1[7]*(res1[5]*res1[1])**2 + res1[7]*(res1[5]*res1[1]*res1[3])**2)
print ('Volume final:')
print (Vfinall)

#time.sleep(60)  # Pausa o código por 5 segundos

