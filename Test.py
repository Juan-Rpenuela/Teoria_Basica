import unittest
import numpy as np
import Observables_y_medidas as om
import math

#Dos pruebas por funcion
class TestSalto_Observables_Medidas(unittest.TestCase):
        def test_prob_1(self):
                v1 = [-3-1j,-2j,1j,2]
                pos = 2
                resultm = 0.05263157894736841
                resultc = om.probSisLi(v1,pos)
                self.assertEqual(resultm, resultc)

        def test_prob_2(self):
                v1 = [2-1j,2j,1-1j,1,-2j,2]
                pos = 3
                resultm = 0.04999999999999999
                resultc = om.probSisLi(v1, pos)
                self.assertEqual(resultm, resultc)

        def test_transicion_1(self):
                v1 = [(2/math.sqrt(2))*1j,-(2/math.sqrt(2))]
                v2 = [(2/math.sqrt(2)),(2/math.sqrt(2))*-1j]
                resultm = 0
                resultc = om.transicion(v1,v2)
                self.assertEqual(resultm, resultc)

        def test_transicion_2(self):
                v1 = [1,-1j]
                v2 = [1j,1]
                resultm = 0.9999999999999996
                resultc = om.transicion(v1, v2)
                self.assertEqual(resultm, resultc)
        
        def test_media_1(self):
                matrix = [[1,-1j],[1j,2]]
                ket = [math.sqrt(2)/2,(math.sqrt(2)/2)*1j]
                resultm = 2.5000000000000004
                resultc = om.media(matrix, ket)
                self.assertEqual(resultm, resultc)

        def test_media_2(self):
                matrix = [[3,1+2j],[1-2j,-1]]
                ket = [math.sqrt(2)/2,(math.sqrt(2)/2)*-1]
                resultm = 0
                resultc = om.media(matrix, ket)
                self.assertEqual(resultm, resultc)

        def test_varianza_1(self):
                matrix = [[1,-1j],[1j,2]]
                ket = [math.sqrt(2)/2,(math.sqrt(2)/2)*1j]
                resultm = 0.25
                resultc = om.varianza(matrix, ket)
                self.assertEqual(resultm, resultc)

        def test_varianza_2(self):
                matrix = [[3,1+2j],[1-2j,-1]]
                ket = [math.sqrt(2)/2,(math.sqrt(2)/2)*-1]
                resultm = 8
                resultc = om.varianza(matrix, ket)
                self.assertEqual(resultm, resultc)

        def test_val_prop_1(self):
                matrix = [[1, 0], [0, -1]]
                resultm = "[ 1. -1.]"
                resultc = str(om.val_prop(matrix))
                self.assertEqual(resultm, resultc)

        def test_val_prop_2(self):
                matrix = [[3, 1 + 2j], [1 - 2j, -1]]
                resultm = "[ 4.+1.85037171e-17j -2.+9.25185854e-17j]"
                resultc = str(om.val_prop(matrix))
                self.assertEqual(resultm, resultc)

        def test_transitar_a_vectores_propios_1(self):
                matrix = [[0, 1], [1, 0]]
                ket = [0, 1]
                resultm =[0.5000000000000001, 0.5000000000000001]
                resultc = om.transitar_vect_prop(matrix,ket)
                self.assertEqual(resultm, resultc)

        def test_transitar_a_vectores_propios_2(self):
                matrix = [[-1,-1j],[1j,1]]
                ket = [1/2,1/2]
                resultm =[0.5, 0.5]
                resultc = om.transitar_vect_prop(matrix,ket)
                self.assertEqual(resultm, resultc)

        def test_finalstate_1(self):
                seq = [[[(1 + 1j)/2, (1 - 1j)/2], [(1 - 1j)/2, (1 + 1j)/2]],[[1, 0], [0, 1j]]]
                ket = [1j, 2j]
                resultm = [(0.5+1.5j), (-1.5-0.5j)]
                resultc = om.final(seq,ket)
                self.assertEqual(resultm, resultc)

        def test_finalstate_2(self):
                seq = [[[0,1],[1,0]],[[(math.sqrt(2))/2,(math.sqrt(2))/2],[(math.sqrt(2))/2,-1*(math.sqrt(2))/2]]]
                ket = [1+6j,3-2j]
                resultm = [(2.8284271247461903+2.8284271247461907j), (1.4142135623730954-5.656854249492381j)]
                resultc = om.final(seq,ket)
                self.assertEqual(resultm,resultc)
    
if __name__ == "__main__":
    unittest.main()