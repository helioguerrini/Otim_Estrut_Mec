# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 21:45:42 2022

@author: Me. Hélio Guerrini Filho
"""

#Importar os módulos
from scipy import optimize

# Parâmetros iniciais de seção
wi, ti = 200.0, 35.0 #mm

# Carga P
P = 20000.0 #N

# Comprimento da viga
L = 2000.0 #mm

# Momento fletor
M = P*L

# Limites das dimensões de seção
winf, wsup = 60.0, 300.0 #mm
tinf, tsup = 10.0, 40.0 #mm

# Módulo de elasticida
E = 210000.0 #MPa

# deflexão máxima
qsup = 10.0 #mm

# Tensão admissível
sigmaAdm = 125 #MPa

# Volume inicial
Vi = L*(wi**2 - (wi - 2*ti)**2)
print(Vi)

#Parâmetro de penalidade
lam = 10

#xstart = array([[wi, ti]],float)
xstart = [wi, ti]

def Pfobj(x):
    w = x[0]
    t = x[1]
    V = L*(w**2 - (w - 2*t)**2)
    f = V/Vi
    
    #Restrições
    g1 = w/winf - 1
    g2 = w/wsup - 1
    g3 = t/tinf - 1
    g4 = t/tsup - 1
    P1 = (min(0.0,g1))**2
    P2 = (max(0.0,g2))**2
    P3 = (min(0.0,g3))**2
    P4 = (max(0.0,g4))**2
    
    # Momento de inércia
    I = (w**4/12) - ((w - 2*t)**4/12)
    # Deflexão da extremidade livre
    q = (P*L**3)/(3*E*I)
    # Restrição de deflexão
    g5 = q/qsup - 1
    P5 = (max(0.0,g5))**2
    
    # Tensão normal
    sigma = (M*w)/(2*I)
    # Restrição de tensão
    g6 = sigma/sigmaAdm - 1
    P6 = (max(0.0,g6))**2
    
    #Pseudo funcão objetivo
    Flinha = f + lam*(P1 + P2 + P3 + P4 + P5 + P6)
    
    return Flinha

res = optimize.fmin(Pfobj, xstart)
print(res)
print('w = ' + str(res[0]))
print('t = ' + str(res[1]))




