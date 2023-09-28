import numpy as np 
import math
def adicc_vec(v1,v2):
    res = []
    if len(v1)== len (v2):
        for j in range (len(v1)):
            res.append(v1[j]+v2[j])
        return res
    else :
        return "Vectores con tamaños diferentes"

def inv_adit(v1):
    res = []
    for j in range (len(v1)):
        res.append(v1[j]*-1)
    return res

def esc_vec(e,v1):
    res = []
    for j in range (len(v1)):
        res.append(v1[j]*e)
    return res

def suma_mat(m1,m2):
    matR= []
    if len(m1) == len (m2):
        for j in range (len(m1)):
            fila = []
            for k in range (len(m1[0])):
                fila.append(m1[j][k]+m2[j][k])
            matR.append(fila)
        return matR
    else :
        return "matrices con tamaños diferentes"

def invAd_mat(m1):
    matR = [[0 for y in range(len(m1[0]))] for x in range(len(m1))]
    for j in range (len(m1)):
        for k in range(len(m1[0])):
            matR[j][k] = -1*(m1[j][k])
    return matR

def esc_mat(e,m1):
    matR = [[0 for y in range(len(m1[0]))] for x in range(len(m1))]
    for j in range (len(m1)):
        for k in range(len(m1[0])):
            matR[j][k] = e*(m1[j][k])
    return matR

def transpuesta(m1):
    matR = [[0 for y in range(len(m1))] for x in range(len(m1[0]))]
    for j in range (len(m1)):
        for k in range(len(m1[0])):
            matR[k][j] = m1[j][k]
    return matR

def conjugada(m1):
    fil = len(m1)
    col = len(m1[0])
    for i in range(fil):
        for j in range(col):
            m1[i][j] = m1[i][j].conjugate()
    return m1

def adjunta(mat):
    mat = transpuesta(mat)
    mat = conjugada(mat)
    return mat

def productoma(mat1,mat2):
    fil_mat1 = len(mat1) 
    col_mat1 = len(mat1[0])
    fil_mat2= len(mat2) 
    col_mat2= len(mat2[0]) 
    if col_mat1 == fil_mat2:
        matriz = []
        for k in range(col_mat2):
            fila = []
            for i in range(fil_mat1):
                suma = 0
                for j in range(col_mat1):
                    suma += mat1[k][j]*mat2[j][i]
                fila.append(suma)
            matriz.append(fila)
    else: 
        matriz = "Tamaños erroneos"
    return matriz

def accionmsobrev(mat,v):
    fil_mat = len(mat) 
    col_mat = len(mat[0]) 
    fil_v = len(v) 
    if col_mat == fil_v:
        vector = []
        for i in range(fil_mat):
            suma = 0
            for j in range(col_mat):
                suma += mat[i][j]*v[j]
            vector.append(suma)
    else: vector = "verificar entrada"
    return vector

def producintern(v1,v2):
    if len(v1) == len(v2):
        suma = 0
        for i in range(len(v1)):
            suma += v1[i].conjugate()*v2[i]
    else: suma = "Error, tamaños erroneos"
    return suma

def normalvec(v):
    b = producintern(v,v)
    b = round(math.sqrt(b.real),2)
    return b

def distvec(u,v):
    h = adicc_vec(u,inv_adit(v))
    dis = normalvec(h)
    return dis

def val_vec_propios(a):
    a = np.asarray(a)
    tam = a.shape
    fil = tam[0]
    col = tam[1]
    if fil == col:
        valores , vectores = np.linalg.eig(a)
        vectores = np.transpose(vectores)
        valor = np.round(valores,2)
        vector = np.round(vectores,2)
        return  f"valores: {valor}. vectores: {vectores[0]}"
    else:
        return "Error, matriz no cuadrada"
    
def unitario(mat):
    if len(mat) == len(mat[0]):
            adjunt = adjunta(mat)
            comprobador = productoma(mat,adjunt)
            verificador = True
            for i in range(len(comprobador)):
                for j in range(len(comprobador)):
                    if (i == j and comprobador[i][j] != 1) or (i != j and comprobador[i][j] != 0):
                        verificador = False
                        break
            if verificador:
                resp = "Es unitaria"
            else:
                resp = "NO es unitaria"
    else:
        resp = "Error, tamaño de matriz incorrecto"
    return resp

def hermitiana(mat):
    if len(mat) == len(mat[0]):
            adjunt = adjunta(mat)
            verificador = True
            for i in range(len(adjunt)):
                for j in range(len(adjunt)):
                    if mat[i][j] != adjunt[i][j]:
                        verificador = False
                        break
            if verificador:
                resp = "Es hermitiana"
            else:
                resp = "NO es hermitiana"
    else:
        resp = "Error, tamaño de matriz incorrecto"
    return resp

def producto_tensor(a,b):
    m = len(a)
    n = len(b)
    m_1 = len(a[0])
    n_1 = len(b[0])
    fil = m*n
    col = m_1*n_1
    result = [[0 for j in range(col)] for k in range(fil)]
    for j in range(fil):
        for k in range(col):
            result[j][k] = (a[j//n][k//n_1])*(b[j%n][k%n_1])
    return result
    

