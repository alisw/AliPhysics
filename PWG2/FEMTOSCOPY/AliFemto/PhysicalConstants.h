/***************************************************************************
 *
 * $Id$
 *
 * Author: CLHEP (see below)
 ***************************************************************************
 *
 * Description:  Taken as-is from CLHEP.
 *               Modified original CVS-Id to retain version info. 
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.2  1999/02/22 16:52:47  didenko
 * updates from Gene
 *
 * Revision 1.1  1999/01/30 03:58:59  fisyak
 * Root Version of StarClassLibrary
 *
 * Revision 1.1  1999/01/23 00:27:34  ullrich
 * Initial Revision
 *
 **************************************************************************/

#ifndef HEP_PHYSICAL_CONSTANTS_H
#define HEP_PHYSICAL_CONSTANTS_H

#include "SystemOfUnits.h"

#ifndef ST_NO_NAMESPACES
using namespace std;
using namespace units;
#endif

#include <math.h>
#include <cmath>

/*
//
#ifndef M_PI
#define M_PI 3.14159265358979312
#endif
//
//#ifdef MACOSX
//extern const double     pi;    // from <math.h>
//#else
static const double     pi  = M_PI;    // from <math.h>
//#endif
static const double  twopi  = 2*pi;
static const double halfpi  = pi/2;
static const double    pi2  = pi*pi;

//
// 
//
static const double Avogadro = 6.0221367e+23/mole;
*/
//
// c   = 299.792458 mm/ns
// c^2 = 898.7404 (mm/ns)^2 
//
static const double c_light   = 2.99792458e+8 * meter/second;
static const double c_squared = c_light * c_light;
/*
//
// h     = 4.13566e-12 MeV*ns
// hbar  = 6.58212e-13 MeV*ns
// hbarc = 197.32705e-12 MeV*mm
//
static const double h_Planck      = 6.6260755e-34 * joule*second;
static const double hbar_Planck   = h_Planck/twopi;
static const double hbarc         = hbar_Planck * c_light;
static const double hbarc_squared = hbarc * hbarc;

//
//
//
static const double electron_charge = - eplus; // see SystemOfUnits.h
static const double e_squared = eplus * eplus;

*/
//
// amu_c2 - atomic equivalent mass unit
// amu    - atomic mass unit
//
static const double electron_mass_c2 = 0.51099906 * MeV;
static const double   proton_mass_c2 = 938.27231 * MeV;
static const double  neutron_mass_c2 = 939.56563 * MeV;
static const double           amu_c2 = 931.49432 * MeV;
//static const double              amu = amu_c2/c_squared;

static const double kaon_0_short_mass_c2 = 497.672  * MeV;
static const double    pion_plus_mass_c2 = 139.5700 * MeV;
static const double   pion_minus_mass_c2 = 139.5700 * MeV;
static const double       lambda_mass_c2 = 1115.684 * MeV;
static const double   antilambda_mass_c2 = 1115.684 * MeV;
static const double     xi_minus_mass_c2 = 1321.32  * MeV;

/*
//
// permeability of free space mu0    = 2.01334e-16 Mev*(ns*eplus)^2/mm
// permittivity of free space epsil0 = 5.52636e+10 eplus^2/(MeV*mm)
//
static const double mu0      = 4*pi*1.e-7 * henry/meter;
static const double epsilon0 = 1./(c_squared*mu0);

//
// electromagnetic coupling = 1.43996e-12 MeV*mm/(eplus^2)
//
static const double elm_coupling           = e_squared/(4*pi*epsilon0);
static const double fine_structure_const   = elm_coupling/hbarc;
static const double classic_electr_radius  = elm_coupling/electron_mass_c2;
static const double electron_Compton_length = hbarc/electron_mass_c2;
static const double Bohr_radius = electron_Compton_length/fine_structure_const;

static const double alpha_rcl2 = fine_structure_const
                                   *classic_electr_radius
                                   *classic_electr_radius;

static const double twopi_mc2_rcl2 = twopi*electron_mass_c2
                                             *classic_electr_radius
                                             *classic_electr_radius;
//
//
//
static const double k_Boltzmann = 8.617385e-11 * MeV/kelvin;

//
//
//
static const double STP_Temperature = 273.15*kelvin;
static const double STP_Pressure    = 1.*atmosphere;
static const double kGasThreshold   = 1.e-2*gram/centimeter3;
*/

#endif /* HEP_PHYSICAL_CONSTANTS_H */





