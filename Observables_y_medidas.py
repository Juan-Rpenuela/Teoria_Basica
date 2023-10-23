import numpy as np
import Libcplxmat as lm
import math
import saltoalocuant as sac
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
def normalizar(v):
    norma = normal(v)
    for i in range(len(v)):
        v[i] = v[i]*(1/norma)
    return v
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
        ket = normalizar(ket)
    resp = lm.hermitiana(mat)
    if resp == "No es hermitiana" or resp == "Tamaño incorrecto":
        return resp
    else:
        final = lm.producintern(lm.accionmsobrev(mat, ket), ket)
        return final.real
    
def varianza(mat, ket):
    if normal(ket) != 1:
        ket = normalizar(ket)
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
        ket = normalizar(ket)
    vectores = vect_prop(mat)
    prob = []
    for i in range(len(vectores)):
        proba = prob_transi(ket,vectores[i])
        prob.append(proba)
    return prob

def final(seq, ket):
    estado = ket
    check = True
    for i in range(len(seq)):
        unit = lm.unitario(seq[i])
        if unit == "No es unit" or unit == "Tamaño incorrecto":
            check = False
            break
    if check:
        for i in range(len(seq)):
            estado = lm.accionmsobrev(seq[i],estado)
        return estado
    else:
        return "No unitarias algunas matrices."
# 4.3.1
print(vect_prop([[0,1],[1,0]]))

# 4.3.2
print(transitar_vect_prop([[0, 1], [1, 0]], [0, 1]))
print(val_prop([[0, 1], [1, 0]]))
print(media([[0, 1], [1, 0]], [0, 1]))

# 4.4.1
A =[[0,1],[1,0]]
c = (2**(1/2))/2
B = [[c,c],[c,-c]]
multi = lm.productoma(A,B)
print(lm.unitario(A))
print(lm.unitario(B))
print(lm.unitario(multi))

# 4.4.2
c = 1/math.sqrt(2)
billar = [[0,c,c,0,],[1j*c,0,0,c],[c,0,0,1j*c],[0,c,-c,0]]
billar1 = [[0,c,c,0,],[c,0,0,-c],[c,0,0,c],[0,-c,c,0]]
state = sac.cplx(billar,[1,0,0,0],3)
state1 = sac.cplx(billar1,[1,0,0,0],3)
print(state)
print(state1)


