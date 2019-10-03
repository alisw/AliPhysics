/***************************************************************************
 *
 * $Id$
 *
 * Author: blasiuk adapted from CLHEP
 ***************************************************************************
 *
 * Description:  This file is based on the SystemOfUnits provided
 *               in the CLHEP library v1.2:  The units remain the same.
 *               It is just the naming conventions that are different:
 *
 * 1) No single letter unit:
 *    : m --> meter
 *    : s --> second
 *    : g --> gram
 *
 * 2) All prefixes are spelled out explicitly (except electron Volt):
 *    : ns --> nanosecond
 *    : mm --> millimeter
 *
 * 3) All units with proper names follow the international standard
 *    of being lower case:
 *    : farad --> farad
 *    : volt  --> volt
 *
 * The basic units are :
 *              centimeter              (centimeter)
 *              second                  (second)
 *              Giga electron Volt      (GeV)
 *              positron charge         (eplus)
 *              degree Kelvin           (kelvin)
 *              the amount of substance (mole)
 *              radian                  (radian)
 *              steradian               (steradian)
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.2.1  2007/10/05 09:38:17  akisiel
 * Fix stray colons
 *
 * Revision 1.1  2007/05/16 10:22:12  akisiel
 * Making the directory structure of AliFemto flat. All files go into one common directory
 *
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.5  2003/09/02 17:59:35  perev
 * gcc 3.2 updates + WarnOff
 *
 * Revision 1.4  1999/03/22 16:21:38  fisyak
 * Add anti CINT flags
 *
 * Revision 1.3  1999/03/11 14:53:07  ullrich
 * Added definition of inch.
 *
 * Revision 1.2  1999/03/02 20:15:08  ullrich
 * Added millivolt.
 *
 * Revision 1.1  1999/01/30 03:59:06  fisyak
 * Root Version of StarClassLibrary
 *
 * Revision 1.1  1999/01/23 00:28:08  ullrich
 * Initial Revision
 *
 **************************************************************************/
#ifndef HEP_SYSTEM_OF_UNITS_H
#define HEP_SYSTEM_OF_UNITS_H


#ifndef M_PI
#define M_PI 3.14159265358979312
#endif


namespace units {
    // new macro for CLHEP SystemOfUnits: at end of file
    //  ST_ADD_OLD_CLHEP_SYSTEM_OF_UNITS
    // 
    // Length [L]
    //
    static const double millimeter  = 0.1;
    static const double millimeter2 = millimeter*millimeter;
    static const double millimeter3 = millimeter*millimeter*millimeter;

    static const double centimeter  = 10*millimeter;
    static const double centimeter2 = centimeter*centimeter;
    static const double centimeter3 = centimeter*centimeter*centimeter;

    static const double meter       = 100.*centimeter;
    static const double meter2      = meter*meter;
    static const double meter3      = meter*meter*meter;

    static const double kilometer   = 1000.*meter;
    static const double kilometer2  = kilometer*kilometer;
    static const double kilometer3  = kilometer*kilometer*kilometer;

    static const double micrometer  = 1.e-6*meter;
    static const double nanometer   = 1.e-9*meter;
    static const double femtometer  = 1.e-15*meter;
    static const double fermi       = 1*femtometer;
    
    static const double      barn   = 1.e-28*meter2;
    static const double millibarn   = 1.e-3*barn;
    static const double microbarn   = 1.e-6*barn;
    static const double  nanobarn   = 1.e-9*barn;
    static const double      inch   = 2.54*centimeter;
    
    //
    // Angle
    //
    static const double      radian = 1.;
    static const double milliradian = 1.e-3*radian;
#ifndef __CINT__
    static const double      degree = (M_PI/180.0)*radian;
#endif    
    static const double   steradian = 1.;

    //
    // Time [T]
    //
    static const double      second = 1;
    static const double millisecond = 1.e-3*second;
    static const double microsecond = 1.e-3*millisecond;
    static const double  nanosecond = 1.e-3*microsecond;
    
    static const double     hertz   = 1./second;
    static const double kilohertz   = 1.e+3*hertz;
    static const double Megahertz   = 1.e+6*hertz;
    
    // but these are also unambiguous and unlikely to be used as variable!
    static const double  Hz         = 1*hertz;
    static const double kHz         = 1*kilohertz;
    static const double MHz         = 1*Megahertz;

    //
    // Electric charge [Q]
    //
    static const double eplus   = 1. ;		        // positron charge
    static const double e_SI    = 1.60217733e-19;	// positron charge in coulomb
    static const double coulomb = eplus/e_SI;
    
    //
    // Energy [E]
    //
    static const double Gigaelectronvolt = 1.;
    static const double Megaelectronvolt = 1.e-3*Gigaelectronvolt;
    static const double     electronvolt = 1.e-6*Megaelectronvolt;
    static const double kiloelectronvolt = 1.e+3*electronvolt;
    static const double Teraelectronvolt = 1.e+3*Gigaelectronvolt;
    
    // but these are also unambiguous and unlikely to be used as variables
    static const double MeV     = Megaelectronvolt;
    static const double  eV     =     electronvolt;
    static const double keV     = kiloelectronvolt;
    static const double GeV     = Gigaelectronvolt;
    static const double TeV     = Teraelectronvolt;
    
    static const double joule   = electronvolt/e_SI;
    
    //
    // Mass [E][T^2][L^-2]
    //
    static const double  kilogram = joule*second*second/(meter*meter);
    static const double      gram = 1.e-3*kilogram;
    static const double milligram = 1.e-3*gram;

    //
    // Power [E][T^-1]
    //
    static const double watt    = joule/second;
    
    //
    // Force [E][L^-1]
    //
    static const double newton  = joule/meter;

    //
    // Pressure [E][L^-3]
    //
#ifndef __CINT__    
#define pascal hep_pascal       // a trick to avoid warnings 
    static const double hep_pascal = newton/meter2;
#else
    static const double pascal     = newton/meter2;
#endif
    static const double bar        = 100000*pascal;
    static const double atmosphere = 101325*pascal;

    //
    // Electric current [Q][T^-1]
    //
    static const double ampere   = coulomb/second;
    
    //
    // Electric potential [E][Q^-1]
    //
    static const double Megavolt = MeV/eplus;
    static const double kilovolt = 1.e-3*Megavolt;
    static const double     volt = 1.e-6*Megavolt;
    static const double millivolt = 1.e-3*volt;
    
    //
    // Electric resistance [E][T][Q^-2]
    //
    static const double ohm = volt/ampere;
    
    //
    // Electric capacitance [Q^2][E^-1]
    //
    static const double farad = coulomb/volt;
    static const double millifarad = 1.e-3*farad;
    static const double microfarad = 1.e-6*farad;
    static const double  nanofarad = 1.e-9*farad;
    static const double  picofarad = 1.e-12*farad;
    
    //
    // Magnetic Flux [T][E][Q^-1]
    //
    static const double weber = volt*second;
    
    //
    // Magnetic Field [T][E][Q^-1][L^-2]
    //
    static const double tesla     = volt*second/meter2;
    
    static const double gauss     = 1.e-4*tesla;
    static const double kilogauss = 1.e-1*tesla;

    //
    // Inductance [T^2][E][Q^-2]
    //
    static const double henry = weber/ampere;

    //
    // Temperature
    //
    static const double kelvin = 1.;

    //
    // Amount of substance
    //
    static const double mole = 1.;
    
    //
    // Activity [T^-1]
    //
    static const double becquerel = 1./second;
    static const double curie = 3.7e+10 * becquerel;
    
    //
    // Absorbed dose [L^2][T^-2]
    //
    static const double gray = joule/kilogram ;

    //
    // Miscellaneous
    //
    static const double perCent     = 0.01 ;
    static const double perThousand = 0.001;
    static const double perMillion  = 0.000001;

#ifdef ST_ADD_OLD_CLHEP_SYSTEM_OF_UNITS

    static const double mm  = 0.1;         // millimeter
    static const double mm2 = mm*mm;
    static const double mm3 = mm*mm*mm;
    
    static const double cm  = 10.*mm;      // centimeter
    static const double cm2 = cm*cm;
    static const double cm3 = cm*cm*cm;
    
    static const double m  = 1000.*mm;     // meter
    static const double m2 = m*m;
    static const double m3 = m*m*m;
    
    static const double km = 1000.*m;      // kilometer
    static const double km2 = km*km;
    static const double km3 = km*km*km;

    static const double microm = 1.e-6*m;  // micro meter
    static const double  nanom = 1.e-9*m;
    //static const double  fermi = 1.e-15*m;

    //
    // Angle
    //
    static const double  rad = 1.;        // radian 
    static const double mrad = 1.e-3*rad; // milliradian
    static const double  deg = (M_PI/180.0)*rad;

    static const double   st = 1.;	    // steradian

    //
    // Time [T]
    //
    static const double  s = 1;           // second
    static const double ns = 1.e-9*s;     // nano second
    static const double ms = 1.e-3*s;     // milli second

    // Mass [E][T^2][L^-2]
    //
    static const double kg = joule*second*second/(meter*meter);	// kg = 6.24150 e+24 * MeV*ns*ns/(mm*mm)   
    static const double  g = 1.e-3*kg;
    static const double mg = 1.e-3*g;

#endif

}
using namespace units;
#endif /* HEP_SYSTEM_OF_UNITS_H */
