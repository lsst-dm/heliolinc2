import numpy as np

from .. import vector as v

def test_vector_module():
    """Unit tests for vector module"""
    a=np.array([-4, 2.2, 0.5])
    check=np.testing.assert_array_almost_equal_nulp
    check2=np.testing.assert_array_almost_equal
    
    # TEST norm
    tst=v.norm(a)
    correct_solution=np.array([4.592385001282014])
    check(tst,correct_solution)
    
    # TEST unit_vector
    tst=v.unitVector(a)
    correct_solution=np.array([-0.8710071126186844, 0.4790539119402764, 0.10887588907733554])
    check(tst,correct_solution)
    
    # TEST rotate_vector
    axs=np.array([-5.6, -12.5, 10.2])
    angle=np.pi/2
    tst=v.rotateVector(angle, axs, a, deg=False)
    correct_solution=np.array([-1.6799623851846146, -2.225115741966377, -3.6491898168248587])
    check(tst,correct_solution,nulp=2)

    phi=np.array([45,90,180])
    axs=np.array([[1,2,1]/np.linalg.norm([1,2,1]),[0,0,1],[0,0,1]])
    vec=np.array([[1,0,0],[1,1,0],[1,0,0]])
    tst=v.rotateVector(phi, axs, vec, deg=True, axis=1)
    correct_solution=np.array([[7.559223176554566e-01,  3.863062075326305e-01,-5.285347327207170e-01],
                               [-1., 1., 0.],
                               [-1., 1.2246468e-16,  0.]])
    check2(tst,correct_solution)
    
#     print(correct_solution)
#     for i in range(3):
#         print(tst[i,:],correct_solution[i,:])
#         check2(tst,correct_solution)
#     return
    
    
    


