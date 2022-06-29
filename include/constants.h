// Constants used in makeTracklets
#ifndef HELA_CONSTANTS_H
#define HELA_CONSTANTS_H

#define DEGPRAD (180.0L/M_PI) /*Degrees per radian*/
#define LSQUARE(x) (long(x)*long(x))
#define DSQUARE(x) (double(x)*double(x))
#define LDSQUARE(x) (((long double)x)*((long double)x))

#define SOLARDAY 86400.0L
#define NEPRA 270.0L //Right ascension of the North Ecliptic Pole.
#define NEPDEC 66.560708333333L // Declination of the North Ecliptic Pole
#define NGPRA 192.728L //Right ascension of the North Galactic Pole (Karim and Mamajek 2017)
#define NGPDEC 26.863L // Declination of the North Galactic Pole (Karim and Mamajek 2017)
#define NCPGAL_LON 122.928L // Galactic longitude of the North Celestial Pole (Karim and Mamajek 2017)
#define MJDOFF 2400000.5L // Offset from Julian Days to Modified Julian days.
#define EARTHEQUATRAD 6378.140L // Earth's equatorial radius.
#define AU_KM 1.495978700e8L /*One AU in km*/
#define AU 1.49597870700e11L /*One AU in meters*/
#define EFOLDS_PER_MAG 0.921034037197618 // E-foldings per astronomical magnitude
#define ZET0 2.5976176L /*precession formula constants in units of arcsec*/
#define ZET1 2306.0809506L
#define ZET2 0.3019015L
#define ZET3 0.0179663L
#define ZET4 -32.7e-6L
#define ZET5 -0.2e-6L
#define Z0 -2.5976176L
#define Z1 2306.0803226L
#define Z2 1.0947790L
#define Z3 0.0182273L
#define Z4 47.0e-6L
#define Z5 -0.3e-6L
#define THET1 2004.1917476L
#define THET2 -0.4269353L
#define THET3 -0.0418251L
#define THET4 -60.1e-6L
#define THET5 -0.1e-6L
#define SMALLANG 1.0e-7L // angle in radians within which declinations are
                         // collapsed to the pole in precess01 routine.
#define VSMALLANG 1.0e-9L // angle difference below which angles are collapsed
                          // to zero in distradec routines.*/

#define TTDELTAT 70.0L // Time in seconds by which TT differs from UT. Thus,
                      // a light-travel-time corrected dynamical position at
                      // TT =  T + TTDELTAT corresponds to the celestial
                      // coordinates observed at UT = T.*/
#define CLIGHT 2.99792458e8L // Speed of light in m/sec
#define CLIGHT_AUDAY CLIGHT*SOLARDAY/AU // Speed of light in AU/day
#define KEPTRANSITMAX 50 // Maximum number of iterations to use in Newton's
                         // method solution of the trancendental Kepler Equation.
#define KEPTRANSTOL 1e-15L // Maximum error for an acceptable solution of the
                           // trancendental Kepler Equation.
#define LARGERR 1e30L // Large number supposed to be a safe initialization
                      // for most minimum-finding problems.
#define GMSUN_KM3_SEC2 132712440041.279419L // GM for the Sun: that is, the Universal
                                            // Gravitational Constant times the solar mass,
                                            // in units of km^3/sec^2.
#define KCONST 0.0172020989484485L // Mean daily motion in radians for an orbit around
                                   // the sun with a semimajor axis of 1AU.
#define MINSTRINGLEN 5 // Minimum size of character array we use: e.g., for filter bandpass or obscode.
#define SHORTSTRINGLEN 20 // Standard size for a short-ish string, used, e.g. for detection idstring
#define MEDSTRINGLEN 80 // Medium string length, should hold most file paths
#define LONGSTRINGLEN 256 // Should hold any reasonable file path, use if not pressed for memory.
#define RAND_MAX_64 18446744073709551616.0L

#define PHASECONST_A1 3.33L // Constants from H-G system phase equation
#define PHASECONST_A2 1.87L // for calculating asteroid apparent magnitudes.
#define PHASECONST_B1 0.63L
#define PHASECONST_B2 1.22L

#endif // HELA_CONSTANTS_H