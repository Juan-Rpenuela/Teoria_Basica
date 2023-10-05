import numpy as np
import Libcplxmat as lm
import math
def modulo(c):
    return c.real**2+c.imag**2


def normal(v):
    b = 0
    for i in range(len(v)):
        b += modulo(v[i])
    b = math.sqrt(b.real)
    return b
#El sistema debe calcular la probabilidad de encontrarlo en una posición en particular.
def probSisLi(estado,pos):
    a = normal(estado)
    return ((modulo(estado[pos])**2/normal(estado)**2))
#El sistema si se le da otro vector Ket debe buscar la probabilidad de transitar del primer vector al segundo.
def transicion (ket1,ket2):
    a = lm.producintern(ket2,ket1)
    a = a/(normal(ket1)*normal(ket2))
    return modulo(a)

def prob_transi(v1,v2):
    p = lm.producintern(v2,v1)/ (normal(v1)*normal(v2))
    return modulo(p)
def media(mat, ket):
    if normal(ket) != 1:
        ket = normal(ket)
    resp = lm.hermitiana(mat)
    if resp == "No es hermitiana" or resp == "Tamaño incorrecto":
        return resp
    else:
        final = lm.prod_interno(lm.sobrevector(mat, ket), ket)
        return final.real
    
def varianza(mat, ket):
    if normal(ket) != 1:
        ket = normal(ket)
    resp = lm.hermitiana(mat)
    if resp == "No es hermitiana" or resp == "Tamaño incorrecto":
        return resp
    else:
        identidad = [[0 for i in range(len(mat))] for j in range(len(mat))]
        for i in range(len(identidad)):
            for j in range(len(identidad)):
                if i == j:
                    identidad[i][j] = 1
        operador = lm.esc_mat(media(mat, ket), lm.invAd_mat(identidad))
        delta = lm.suma_mat(mat, operador)
        cuadrado = lm.productoma(delta,delta)
        return media(cuadrado,ket)

def val_prop(mat):
    mat = np.array(mat)
    valores = np.linalg.eigvals(mat)
    return valores

def vect_prop(mat):
    mat = np.array(mat)
    valores, vectores = np.linalg.eig(mat)
    return vectores

def transitar_vect_prop(mat, ket):
    if normal(ket) != 1:
        ket = normal(ket)
    vectores = vect_prop(mat)
    prob = []
    for i in range(len(vectores)):
        proba = prob_transi(ket,vectores[i])
        prob.append(proba)
    return prob

def final(seq, ket):
    estado = ket
    chack = True
    for i in range(len(seq)):
        unit = lm.unitario(seq[i])
        if unit == "No es unit" or unit == "Tamaño incorrecto":
            chack = False
            break
    if chack:
        for i in range(len(seq)):
            estado = lm.accionmsobrev(seq[i],estado)
        return estado
    else:
        return "No unitarias algunas matrices."

print(probSisLi([2+1j,-1+2j,1j,1,3-1j,2,-2j,-2+1j,1-3j,-1j],7))
print(modulo([2+1j,-1+2j,1j,1,3-1j,2,-2j,-2+1j,1-3j,-1j]))

