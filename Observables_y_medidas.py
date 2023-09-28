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



