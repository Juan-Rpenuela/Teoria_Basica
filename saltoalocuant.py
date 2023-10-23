#Juan Andres Rodriguez
import numpy as np
import math
import Libcplxmat as lm
import matplotlib.pyplot as plot


def ModuloCuadrado(c):
    return round(c.real**2+c.imag**2,2)

def Calculo(mtrx,vecto,clicks):
    if clicks == 0:
        return vecto
    elif clicks == 1:
        return lm.accionmsobrev(mtrx,vecto)
    else:
        return lm.accionmsobrev(mtrx,Calculo(mtrx,vecto,clicks-1))

def Canicas(mtrx,vecto,clics):
    if len(mtrx) != len(mtrx[0]):
        result = "Tamaño incorrecto de matriz"
    for i in range(len(mtrx)):
        for j in range(len(mtrx[0])):
            if mtrx[i][j] != 0 or mtrx[i][j] != 1:
                result = "Matriz no booleana"
                break
    else:
        result = Calculo(mtrx,vecto,clics)
    return result

def Fracciones(mtrx,vecto,clics):
    if len(mtrx) != len(mtrx[0]):
        result = "Tamaño incorrecto de matriz"
    for i in range(len(mtrx)):
        suma = 0
        for j in range(len(mtrx[0])):
            suma += mtrx[i][j]
        if suma != 1:
            result = "Matriz no doblemente estocástica"
            break
    for i in range(len(lm.transpuesta(mtrx))):
        suma = 0
        for j in range(len(lm.transpuesta(mtrx)[0])):
            suma += mtrx[i][j]
        if suma != 1:
            result = "Matriz no doblemente estocástica"
            break
    else:
        result = Calculo(mtrx,vecto,clics)
        result = [round(x,2) for x in result]
    return result

def ExperimentoRendijas(blancos,rendijas):

    matriz = [[0 for x in range(blancos+rendijas+1)] for x in range(blancos+rendijas+1)]
    cada = blancos//rendijas if blancos%2==0 else (blancos+1)//(rendijas)

    j=rendijas+1

    for k in range(1,rendijas+1):
        matriz[k][0]=round(1/rendijas,2)
        for i in range(j,j+cada):
            matriz[i][k]= round(1/cada,2)
        j+=cada-1

        for l in range(len(matriz)-blancos,len(matriz)):
            matriz[l][l]=1
    vector = [0 for x in range(len(matriz))]
    vector[0] = 1
    matriz1 = lm.productoma(matriz,matriz)
    vector = lm.accionmsobrev(matriz1,vector)

    return matriz,vector


def cplx(mtrx,vecto,clics):
    if len(mtrx) != len(mtrx[0]):
        result = "Tamaño incorrecto de matriz"
    for i in range(len(mtrx)):
        suma = 0
        for j in range(len(mtrx[0])):
            suma += ModuloCuadrado(mtrx[i][j])
        if suma != 1:
            result = "Matriz no doblemente estocástica"
            break
    for i in range(len(lm.transpuesta(mtrx))):
        suma = 0
        for j in range(len(lm.transpuesta(mtrx)[0])):
            suma += ModuloCuadrado(mtrx[i][j])
        if suma != 1:
            result = "Matriz no doblemente estocástica"
            break
    else:
        result = Calculo(mtrx,vecto,clics)
    return result

def ExperimentoRendijasCua(mtrx):
    
    vector = [0 for x in range(len(mtrx))]
    vector[0] = 1
    matriz = lm.productoma(mtrx,mtrx)
    vector = lm.accionmsobrev(mtrx,vector)

    return vector

print(ExperimentoRendijasCua([[0,0,0,0,0,0,0,0],[1/math.sqrt(2),0,0,0,0,0,0,0],[1/math.sqrt(2),0,0,0,0,0,0,0],[0,(-1+1j)/math.sqrt(6),0,1,0,0,0,0],[0,(-1-1j)/math.sqrt(6),0,0,1,0,0,0],[0,(1-1j)/math.sqrt(6),(-1+1j)/math.sqrt(6),0,0,1,0,0],[0,0,(-1-1j)/math.sqrt(6),0,0,0,1,0],[0,0,(1-1j)/math.sqrt(6),0,0,0,0,1]]))

def graficar(vector):
    data =len(vector)
    x = np.array([x for x in range(data)])
    y = np.array([round(vector[x] * 100, 2) for x in range(data)])

    plot.bar(x, y, color='b', align='center')
    plot.title('Probabilidad del vector en unidades de %')
    plot.savefig("Caso.jpg")
    plot.show()


