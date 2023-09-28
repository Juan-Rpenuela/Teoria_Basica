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
    
if __name__ == "__main__":
    unittest.main()