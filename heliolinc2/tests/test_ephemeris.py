import numpy as np

# time scale transforms
from astropy.time import Time

import ..transform as tr

import ..ephemeris as ephem


def test_transforms():
    """Unit tests for transforms module"""
    
    # Checks
    check=np.testing.assert_array_almost_equal_nulp
    check2=np.testing.assert_allclose
    
    au2km=149597870.700
    # Dynamical states and orbital elements in TDB timescale
    epoch=Time(2456165.500000000, format='jd', scale='tdb')
    tdb_minus_ut=epoch.tdb.mjd-epoch.utc.mjd
        
    # 367943 Duende
    duende_ecliptic_state=np.array([6.832993744130379E-01, -7.578357776615274E-01,  4.887250343702505E-02,
                                    1.115044267163202E-02,  1.227682492646761E-02, -2.98353508862975E-03])
    
    duende_icrf_state=np.array([6.832993744130379E-01, -7.147410974191527E-01, -2.566101150697676E-01,
                                1.115044267163202E-02,  1.245054875137605E-02, 2.146100577752525E-03])
    # q,e,i,node,w,tp
    duende_cometary_orbital_elements=np.array([8.934926246676773E-01,1.081293521931094E-01,1.033729958423475E+01,
                                                1.472624792437615E+02, 2.710862614793172E+02, 2456261.553285983857]) 
    # a,e,i,w,node,M
    duende_keplerian_orbital_elements=np.array([1.001818623434660E+00, 1.081293521931094E-01,  1.033729958423475E+01,
                                                2.710862614793172E+02,1.472624792437615E+02, 2.655868148489221E+02])                                                             
    # Right Ascension and Declination in UTC
    duende_radec_ut=np.array([313.706696703, -14.548583617])
    
    # WARNING: JPL dRA/dt and dDEC/dt are not in ICRF! They can only serve as order of magnitude checks.
    duende_dradec_dt_ut=np.array([140.2116,  -14.4983])/3600*24
    
    duende_r_rdot=np.array([1.02157736197894, -3.1025058/au2km*86400])
    
    # 1 Ceres
    ceres_ecliptic_state=np.array([1.307798047465823E+00,  2.421468434364515E+00 ,-1.650429115858807E-01,
                                   -9.298551562868396E-03,  4.238614777929977E-03,  1.847323238873149E-03])
    
    ceres_icrf_state=np.array([1.307798047465823E+00,  2.287304152373543E+00,  8.117809161487053E-01,
                               -9.298551562868396E-03,  3.154030042726210E-03,  3.380910065965624E-03])
    # q,e,i,node,w,tp
    ceres_cometary_orbital_elements=np.array([2.556737030655413E+00,7.705127985019911E-02, 1.059239607081042E+01,
                                              8.033159124626791E+01, 7.187631305619260E+01,2456549.349030518439])
    # a,e,i,w,node,M
    ceres_keplerian_orbital_elements=np.array([2.770183190936585E+00,7.705127985019911E-02, 1.059239607081042E+01,
                                                7.187631305619260E+01,8.033159124626791E+01,2.779456163962251E+02])
    ceres_radec_ut=np.array([60.237463150,  17.122942759])
    # WARNING: JPL dRA/dt and dDEC/dt are not in ICRF! They can only serve as order of magnitude checks.
    ceres_dradec_dt_ut=np.array([32.23178,  1.685305 ])
    
    ceres_r_rdot=np.array([2.75701848771017,  -1.3828025/au2km*86400])
                                                
    # 15147 Siegfried
    siegfried_ecliptic_state=np.array([-1.955900170259298E+00, -2.611267177109618E+00,  4.342623471545343E-02,
                                       7.563883190764837E-03, -5.557598139717739E-03, -1.472566354935759E-03])
    
    siegfried_icrf_state=np.array([-1.955900170259298E+00, -2.413064758406049E+00, -9.988596397143862E-01,
                                   7.563883190764837E-03, -4.513243344792708E-03, -3.561738797689575E-03])
                                   
    siegfried_radec_ut=np.array([230.970732638, -17.825060926]) 
    # WARNING: JPL dRA/dt and dDEC/dt are not in ICRF! They can only serve as order of magnitude checks.
    siegfried_dradec_dt_ut=np.array([231.131603719, -17.891296611])/3600*24  
    
    siegfried_r_rdot=np.array([3.26284541961338,  -0.1834848/au2km*86400])
    # q,e,i,node,w,tp
    siegfried_cometary_orbital_elements=np.array([3.208390751707551E+00,1.209336623887380E-02,8.941381299848775E+00,
                                            5.801891370459320E+01,2.884583753517642E+02, 2456831.134322461672])
    # a,e,i,w,node,M
    siegfried_keplerian_orbital_elements=np.array([3.247665965651702E+00,1.209336623887380E-02,8.941381299848775E+00,
                                                   2.884583753517642E+02,5.801891370459320E+01,2.480719216306982E+02])
                                   
    # TESTS FOR SINGLE STATES
    
    # TEST ecliptic2icrf
    tst=tr.ecliptic2icrf(duende_ecliptic_state)
    correct_solution=duende_icrf_state
    check2(tst,correct_solution,rtol=1e-15, atol=1e-14)
    
    # TEST icrf2eclip
    tst=tr.icrf2ecliptic(duende_icrf_state)
    correct_solution=duende_ecliptic_state
    check2(tst,correct_solution,rtol=1e-15, atol=1e-14)
    
    # TEST icrf2ephemeris
#     tst=tr.icrf2ephemeris(epoch,duende_icrf_state,timescale_epoch='utc',timescale_state='tdb',deg=True,lttc=True)
#     correct_solution=np.hstack([duende_radec_ut,duende_dradec_dt_ut,duende_r_rdot])
#     tst=tr.icrf2ephemeris(epoch,ceres_icrf_state,timescale_epoch='utc',timescale_state='tdb',deg=True,lttc=True)
#     correct_solution=np.hstack([ceres_radec_ut,ceres_dradec_dt_ut,ceres_r_rdot])
    
    tst=ephem.icrf2ephemeris(epoch,siegfried_icrf_state,timescale_epoch='utc',timescale_state='tdb',deg=True,lttc=True)
    correct_solution=np.hstack([siegfried_radec_ut,siegfried_dradec_dt_ut,siegfried_r_rdot])
    
    print(tst)
    print(correct_solution)
    print('O-C: RA & DEC [mas], dRA/dt*cos(DEC) & dDEC/dt [arcsec/hr], r [km], rdot [km/s]')
    print((tst-correct_solution)[0:2]*3600*1000,(tst-correct_solution)[2:4]*3600/24., 
          (tst-correct_solution)[4]*au2km,(tst-correct_solution)[5]*au2km/86400)
    # 5mas tolerance for RA & DEC
    check2(tst[0:1],correct_solution[0:1],rtol=1e-15, atol=1/3600/1000*10)
    # 1% tolerance for sky plane motion
    #check2(np.abs(tst[2:3]),np.abs(correct_solution[2:3]),rtol=1e-15, atol=0.01)
    
    check2(np.abs(tst[4:5]),np.abs(correct_solution[4:5]),rtol=1e-15, atol=10/au2km)
   
    # TESTS FOR MULTIPLE STATES
    
    # TEST ecliptic2icrf
    ecliptic_states=np.array([duende_ecliptic_state,ceres_ecliptic_state,siegfried_ecliptic_state])                                            
    icrf_states=np.array([duende_icrf_state,ceres_icrf_state,siegfried_icrf_state])     
    
    tst=tr.ecliptic2icrf(ecliptic_states[:,0:3])
    correct_solution=icrf_states[:,0:3]  
    #print(tst)
    #print(correct_solution)
    #print((tst-correct_solution))
    check2(tst,correct_solution,rtol=1e-15, atol=1e-14)
    
    # TEST icrf2eclip
    tst=tr.ecliptic2icrf(ecliptic_states)
    correct_solution=icrf_states  
    #print(tst)
    #print(correct_solution)
    #print((tst-correct_solution))
    check2(tst,correct_solution,rtol=1e-15, atol=1e-14)
    
    # TEST icrf2ephemeris
    tst=ephem.icrf2ephemeris(epoch,icrf_states,timescale_epoch='utc',timescale_state='tdb',deg=True,lttc=True)
    correct_solution=np.vstack([np.hstack([duende_radec_ut,duende_dradec_dt_ut, duende_r_rdot]),
                               np.hstack([ceres_radec_ut,ceres_dradec_dt_ut,ceres_r_rdot]),
                               np.hstack([siegfried_radec_ut,siegfried_dradec_dt_ut,siegfried_r_rdot])])
                                          
    print('icrf2ephemeris')
    print(tst)
    print(correct_solution)
    print('O-C: RA & DEC [mas], dRA/dt*cos(DEC) & dDEC/dt [arcsec/hr], r [km], rdot [km/s]')
    print((tst-correct_solution)[:,0:2]*3600*1000)
    print((tst-correct_solution)[:,2:4]*3600/24.)
    print((tst-correct_solution)[:,4]*au2km,(tst-correct_solution)[:,5]*au2km/86400)
          

    # 5mas tolerance for RA & DEC
    check2(tst[:,0:1],correct_solution[:,0:1],rtol=1e-15, atol=1/3600/1000*10)
    # 1% tolerance for sky plane motion
    #check2(np.abs(tst[2:3]),np.abs(correct_solution[2:3]),rtol=1e-15, atol=0.01)
    check2(np.abs(tst[:,4:5]),np.abs(correct_solution[:,4:5]),rtol=1e-15, atol=10/au2km)
    
#     tst=tr.icrf2radec(epoch,duende_icrf_state[0:6],deg=True,lttc=True)+duende_dradec_dt_ut*tdb_minus_ut
#     print(tst)
#     print(correct_solution)
#     print((tst-correct_solution)*3600)
#     check(tst,correct_solution,nulp=1)
    
    return icrf_states