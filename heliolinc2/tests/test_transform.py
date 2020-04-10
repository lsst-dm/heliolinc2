import numpy as np

# time scale transforms
from astropy.time import Time

import ..transform as tr

def test_orbital_element_conversion():
    
    check2=np.testing.assert_allclose
    
    epoch=Time(2456165.500000000, format='jd', scale='tdb').value
    
    duende_ecliptic_state = np.array([6.832993744130379E-01, -7.578357776615274E-01,  4.887250343702505E-02,
                                      1.115044267163202E-02,  1.227682492646761E-02, -2.98353508862975E-03])
    # q,e,i,node,w,tp
    duende_cometary_orbital_elements = np.array([8.934926246676773E-01,1.081293521931094E-01,1.033729958423475E+01,
                                                 1.472624792437615E+02, 2.710862614793172E+02, 2456261.553285983857]) 
    # a,e,i,w,node,M
    duende_keplerian_orbital_elements = np.array([1.001818623434660E+00, 1.081293521931094E-01,  1.033729958423475E+01,
                                                  2.710862614793172E+02,1.472624792437615E+02, 2.655868148489221E+02])
    
    # q,e,i,node,w,tp
    ceres_cometary_orbital_elements = np.array([2.556737030655413E+00,7.705127985019911E-02, 1.059239607081042E+01,
                                                8.033159124626791E+01, 7.187631305619260E+01,2456549.349030518439])
    # a,e,i,w,node,M
    ceres_keplerian_orbital_elements = np.array([2.770183190936585E+00,7.705127985019911E-02, 1.059239607081042E+01,
                                                 7.187631305619260E+01,8.033159124626791E+01,2.779456163962251E+02])
    
    
     # q,e,i,node,w,tp
    siegfried_cometary_orbital_elements = np.array([3.208390751707551E+00,1.209336623887380E-02,8.941381299848775E+00,
                                                    5.801891370459320E+01,2.884583753517642E+02, 2456831.134322461672])
    # a,e,i,w,node,M
    siegfried_keplerian_orbital_elements = np.array([3.247665965651702E+00,1.209336623887380E-02,8.941381299848775E+00,
                                                     2.884583753517642E+02,5.801891370459320E+01,2.480719216306982E+02])
    

    cart2com = tr.cartesian2cometary(epoch, duende_ecliptic_state, frame='ecliptic')
    
    cart2kep = tr.cartesian2keplerian(epoch, duende_ecliptic_state, frame='ecliptic')
    
    kep2cart = tr.keplerian2cartesian(epoch, duende_keplerian_orbital_elements, frame='ecliptic')
    
    print(cart2com)
    print(cart2com[0])
        
    com2kep = tr.cometary2keplerian(epoch, cart2com[0])
        
    com2cart = tr.cometary2cartesian(epoch, cart2com[0])
        
    print('Ecliptic states')
    print(duende_ecliptic_state)
    print(com2cart)
    print(kep2cart)  
    correct_solution = duende_ecliptic_state
    
    tst = kep2cart
    check2(tst,correct_solution,rtol=1e-10, atol=1e-11)
    
    print('Cometary elements')
    print(duende_cometary_orbital_elements)
    print(cart2com)
    correct_solution = duende_cometary_orbital_elements
    tst = cart2com[0]
    tst[5]=tst[5]+cart2com[2]
    print(tst)
    check2(tst[0:4],correct_solution[0:4],rtol=1e-9, atol=1e-9)

    
    print('Keplerian elements')
    print(duende_keplerian_orbital_elements)
    print(com2kep)
    print(cart2kep)
    
    correct_solution = duende_keplerian_orbital_elements
    tst = cart2kep[0]
    check2(tst,correct_solution,rtol=1e-7, atol=1e-7)
    
    return