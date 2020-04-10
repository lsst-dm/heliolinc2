import numpy as np

import ..vector 

def test_vector_module():
    """Unit tests for vector module"""
    a=np.array([-4, 2.2, 0.5])
    check=np.testing.assert_array_almost_equal_nulp
    
    # TEST norm
    tst=v.norm(a)
    correct_solution=np.array([4.592385001282014])
    check(tst,correct_solution)
    
    # TEST unit_vector
    tst=v.unitVector(a)
    correct_solution=np.array([-0.8710071126186844, 0.4790539119402764, 0.10887588907733554])
    check(tst,correct_solution)
    
    # TEST rotate_vector
    axis=np.array([-5.6, -12.5, 10.2])
    angle=np.pi/2
    tst=v.rotateVector(angle, axis, a, deg=False)
    correct_solution=np.array([-1.6799623851846146, -2.225115741966377, -3.6491898168248587])
    check(tst,correct_solution,nulp=2)

    return
    
    
    


