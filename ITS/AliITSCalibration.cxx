/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//////////////////////////////////////////////////////
//  Calibration class for set:ITS                   //
//  Specific subdetector implementation is done in  //
//  AliITSCalibrationSPD                            //
//  AliITSCalibrationSDD                            //
//  AliITSCalibrationSSD                            //
//////////////////////////////////////////////////////

#include <Riostream.h>
#include <TMath.h>

#include "AliITSCalibration.h"
#include "AliLog.h"

ClassImp(AliITSCalibration)

//______________________________________________________________________
AliITSCalibration::AliITSCalibration():
TObject(),
fDataType(),
fdv(0.000375),
fN(0.),
fT(300.),
fGeVcharge(0.),
fResponse(){
    // Default Constructor (300 microns and 80 volts)

    SetGeVToCharge();
    fResponse = 0;
}
//______________________________________________________________________
AliITSCalibration::AliITSCalibration(Double_t thickness):
TObject(),
fDataType(),
fdv(0.),
fN(0.),
fT(300.),
fGeVcharge(0.),
fResponse(){
    // Default Constructor

    fdv = thickness/80.0;   // 80 volts.
    SetGeVToCharge();
    fResponse = 0;
}

//______________________________________________________________________
AliITSCalibration::AliITSCalibration(const AliITSCalibration &ob):
TObject(ob),
fDataType(ob.fDataType),
fdv(ob.fdv),
fN(ob.fN),
fT(ob.fT),
fGeVcharge(ob.fGeVcharge),
fResponse(ob.fResponse)
{
  // Copy constructor

}
/*
//______________________________________________________________________________
AliITSCalibration& AliITSCalibration::operator= (const AliITSCalibration& source)
{
  // Asignment operator

  this->~AliITSCalibration();
  new(this) AliITSCalibration(source);
  return *this;


}
*/
//______________________________________________________________________
Double_t AliITSCalibration::MobilityElectronSiEmp() const {
    // Computes the electron mobility in cm^2/volt-sec. Taken from SILVACO
    // International ATLAS II, 2D Device Simulation Framework, User Manual 
    // Chapter 5 Equation 5-6. An empirical function for low-field mobiliity 
    // in silicon at different tempeatures.
    // Inputs:
    //    none.
    // Output:
    //    none.
    // Return:
    //    The Mobility of electrons in Si at a give temprature and impurity
    //    concentration. [cm^2/Volt-sec]
    const Double_t km0  = 55.24; // cm^2/Volt-sec
    const Double_t km1  = 7.12E+08; // cm^2 (degree K)^2.3 / Volt-sec
    const Double_t kN0  = 1.072E17; // #/cm^3
    const Double_t kT0  = 300.; // degree K.
    const Double_t keT0 = -2.3; // Power of Temp.
    const Double_t keT1 = -3.8; // Power of Temp.
    const Double_t keN  = 0.73; // Power of Dopent Consentrations
    Double_t m;
    Double_t tT = fT,nN = fN;

    if(nN<=0.0){ // Simple case.
	if(tT==300.) return 1350.0; // From Table 5-1 at consentration 1.0E14.
	m = km1*TMath::Power(tT,keT0);
	return m;
    } // if nN<=0.0
    m = km1*TMath::Power(tT,keT0) - km0;
    m /= 1.0 + TMath::Power(tT/kT0,keT1)*TMath::Power(nN/kN0,keN);
    m += km0;
    return m;
}
//______________________________________________________________________
Double_t AliITSCalibration::MobilityHoleSiEmp() const {
    // Computes the Hole mobility in cm^2/volt-sec. Taken from SILVACO
    // International ATLAS II, 2D Device Simulation Framework, User Manual 
    // Chapter 5 Equation 5-7 An empirical function for low-field mobiliity 
    // in silicon at different tempeatures.
    // Inputs:
    //    none.
    // Output:
    //    none.
    // Return:
    //    The Mobility of Hole in Si at a give temprature and impurity
    //    concentration. [cm^2/Volt-sec]
    const Double_t km0a = 49.74; // cm^2/Volt-sec
    const Double_t km0b = 49.70; // cm^2/Volt-sec
    const Double_t km1  = 1.35E+08; // cm^2 (degree K)^2.3 / Volt-sec
    const Double_t kN0  = 1.606E17; // #/cm^3
    const Double_t kT0  = 300.; // degree K.
    const Double_t keT0 = -2.2; // Power of Temp.
    const Double_t keT1 = -3.7; // Power of Temp.
    const Double_t keN  = 0.70; // Power of Dopent Consentrations
    Double_t m;
    Double_t tT = fT,nN = fN;

    if(nN<=0.0){ // Simple case.
	if(tT==300.) return 495.0; // From Table 5-1 at consentration 1.0E14.
	m = km1*TMath::Power(tT,keT0) + km0a-km0b;
	return m;
    } // if nN<=0.0
    m = km1*TMath::Power(tT,keT0) - km0b;
    m /= 1.0 + TMath::Power(tT/kT0,keT1)*TMath::Power(nN/kN0,keN);
    m += km0a;
    return m;
}
//______________________________________________________________________
Double_t AliITSCalibration::DiffusionCoefficientElectron() const {
    // Computes the Diffusion coefficient for electrons in cm^2/sec. Taken 
    // from SILVACO International ATLAS II, 2D Device Simulation Framework, 
    // User Manual Chapter 5 Equation 5-53. Einstein relations for diffusion 
    // coefficient. Note: 1 cm^2/sec = 10 microns^2/nanosec.
    // Inputs:
    //    none.
    // Output:
    //    none.
    // Return:
    //    The Diffusion Coefficient of electrons in Si at a give temprature
    //    and impurity concentration. [cm^2/sec]
    // const Double_t kb = 1.3806503E-23; // Joules/degree K
    // const Double_t qe = 1.60217646E-19; // Coulumbs.
    const Double_t kbqe = 8.617342312E-5; // Volt/degree K
    Double_t m = MobilityElectronSiEmp();
    Double_t tT = fT;

    return m*kbqe*tT;  // [cm^2/sec]
}
//______________________________________________________________________
Double_t AliITSCalibration::DiffusionCoefficientHole() const {
    // Computes the Diffusion coefficient for Holes in cm^2/sec. Taken 
    // from SILVACO International ATLAS II, 2D Device Simulation Framework, 
    // User Manual Chapter 5 Equation 5-53. Einstein relations for diffusion 
    // coefficient. Note: 1 cm^2/sec = 10 microns^2/nanosec.
    // Inputs:
    //    none.
    // Output:
    //    none.
    // Return:
    //    The Defusion Coefficient of Hole in Si at a give temprature and 
    //    impurity concentration. [cm^2/sec]
    //    and impurity concentration. [cm^2/sec]
    // const Double_t kb = 1.3806503E-23; // Joules/degree K
    // const Double_t qe = 1.60217646E-19; // Coulumbs.
    const Double_t kbqe = 8.617342312E-5; // Volt/degree K
    Double_t m = MobilityHoleSiEmp();
    Double_t tT = fT;

    return m*kbqe*tT;  // [cm^2/sec]
}
//______________________________________________________________________
Double_t AliITSCalibration::SpeedElectron() const {
    // Computes the average speed for electrons in Si under the low-field 
    // approximation. [cm/sec].
    // Inputs:
    //    none.
    // Output:
    //    none.
    // Return:
    //    The speed the holes are traveling at due to the low field applied. 
    //    [cm/sec]
    Double_t m = MobilityElectronSiEmp();

    return m/fdv;  // [cm/sec]
}
//______________________________________________________________________
Double_t AliITSCalibration::SpeedHole() const {
    // Computes the average speed for Holes in Si under the low-field 
    // approximation.[cm/sec].
    // Inputs:
    //    none.
    // Output:
    //    none.
    // Return:
    //    The speed the holes are traveling at due to the low field applied. 
    //    [cm/sec]
    Double_t m = MobilityHoleSiEmp();

    return m/fdv;  // [cm/sec]
}
//______________________________________________________________________
Double_t AliITSCalibration::SigmaDiffusion3D(Double_t l) const {
    // Returns the Gaussian sigma^2 == <x^2+y^2+z^2> [cm^2] due to the
    // defusion of electrons or holes through a distance l [cm] caused 
    // by an applied voltage v [volt] through a distance d [cm] in any
    //  material at a temperature T [degree K]. The sigma diffusion when
    //  expressed in terms of the distance over which the diffusion 
    // occures, l=time/speed, is independent of the mobility and therefore
    //  the properties of the material. The charge distributions is given by 
    // n = exp(-r^2/4Dt)/(4piDt)^1.5. From this <r^2> = 6Dt where D=mkT/e
    // (m==mobility, k==Boltzman's constant, T==temparature, e==electric 
    // charge. and vel=m*v/d. consiquently sigma^2=6kTdl/ev.
    // Inputs:
    //    Double_t l   Distance the charge has to travel.
    // Output:
    //    none.
    // Return:
    //    The Sigma due to the diffution of electrons. [cm]
    const Double_t kcon = 5.17040258E-04; // == 6k/e [J/col or volts]

    return TMath::Sqrt(kcon*fT*fdv*l);  // [cm]
}
//______________________________________________________________________
Double_t AliITSCalibration::SigmaDiffusion2D(Double_t l) const {
    // Returns the Gaussian sigma^2 == <x^2+z^2> [cm^2] due to the defusion 
    // of electrons or holes through a distance l [cm] caused by an applied
    // voltage v [volt] through a distance d [cm] in any material at a
    // temperature T [degree K]. The sigma diffusion when expressed in terms
    // of the distance over which the diffusion occures, l=time/speed, is 
    // independent of the mobility and therefore the properties of the
    // material. The charge distributions is given by 
    // n = exp(-r^2/4Dt)/(4piDt)^1.5. From this <x^2+z^2> = 4Dt where D=mkT/e
    // (m==mobility, k==Boltzman's constant, T==temparature, e==electric 
    // charge. and vel=m*v/d. consiquently sigma^2=4kTdl/ev.
    // Inputs:
    //    Double_t l   Distance the charge has to travel.
    // Output:
    //    none.
    // Return:
    //    The Sigma due to the diffution of electrons. [cm]
    const Double_t kcon = 3.446935053E-04; // == 4k/e [J/col or volts]

    return TMath::Sqrt(kcon*fT*fdv*l);  // [cm]
}
//______________________________________________________________________
Double_t AliITSCalibration::SigmaDiffusion1D(Double_t l) const {
    // Returns the Gaussian sigma^2 == <x^2> [cm^2] due to the defusion 
    // of electrons or holes through a distance l [cm] caused by an applied
    // voltage v [volt] through a distance d [cm] in any material at a
    // temperature T [degree K]. The sigma diffusion when expressed in terms
    // of the distance over which the diffusion occures, l=time/speed, is 
    // independent of the mobility and therefore the properties of the
    // material. The charge distributions is given by 
    // n = exp(-r^2/4Dt)/(4piDt)^1.5. From this <r^2> = 2Dt where D=mkT/e
    // (m==mobility, k==Boltzman's constant, T==temparature, e==electric 
    // charge. and vel=m*v/d. consiquently sigma^2=2kTdl/ev.
    // Inputs:
    //    Double_t l   Distance the charge has to travel.
    // Output:
    //    none.
    // Return:
    //    The Sigma due to the diffution of electrons. [cm]
    const Double_t kcon = 1.723467527E-04; // == 2k/e [J/col or volts]

    return TMath::Sqrt(kcon*fT*fdv*l);  // [cm]
}
//----------------------------------------------------------------------
Double_t AliITSCalibration::DepletedRegionThicknessA(Double_t dopCons,
                                                 Double_t voltage,
                                                 Double_t elecCharge,
                                                 Double_t voltBuiltIn)const{
    // Computes the thickness of the depleted region in Si due to the 
    // application of an external bias voltage. From the Particle Data
    // Book, 28.8 Silicon semiconductor detectors equation 28.19 (2004)
    // Physics Letters B "Review of Particle Physics" Volume 592, Issue 1-4
    // July 15 2004, ISSN 0370-2693 page 263. First equation.
    // Inputs:
    //    Double_t dopCons           "N" doping concentration
    //    Double_t voltage           "V" external bias voltage
    //    Double_t elecCharge        "e" electronic charge
    //    Double_t voltBuiltIn=0.5   "V_bi" "built-in" Voltage (~0.5V for
    //                               resistivities typically used in detectors)
    // Output:
    //    none.
    // Return:
    //    The thickness of the depleted region

    return TMath::Sqrt(2.0*(voltage+voltBuiltIn)/(dopCons*elecCharge));
}
//----------------------------------------------------------------------
Double_t AliITSCalibration::DepletedRegionThicknessB(Double_t resist,
                                                 Double_t voltage,
                                                 Double_t mobility,
                                                 Double_t voltBuiltIn,
                                                 Double_t dielConst)const{
    // Computes the thickness of the depleted region in Si due to the 
    // application of an external bias voltage. From the Particle Data
    // Book, 28.8 Silicon semiconductor detectors equation 28.19 (2004)
    // Physics Letters B "Review of Particle Physics" Volume 592, Issue 1-4
    // July 15 2004, ISSN 0370-2693 page 263. Second Equation.
    // Inputs:
    //    Double_t resist            "rho" resistivity (typically 1-10 kOhm cm)
    //    Double_t voltage           "V" external bias voltage
    //    Double_t mobility          "mu" charge carrier mobility
    //                                  (electons 1350, holes 450 cm^2/V/s)
    //    Double_t voltBuiltIn=0.5   "V_bi" "built-in" Voltage (~0.5V for
    //                               resistivities typically used in detectors)
    //    Double_t dielConst=1.E-12  "epsilon" dielectric constant = 11.9 *
    //                                (permittivity of free space) or ~ 1 pF/cm
    // Output:
    //    none.
    // Return:
    //    The thickness of the depleted region

    return TMath::Sqrt(2.8*resist*mobility*dielConst*(voltage+voltBuiltIn));
}
//----------------------------------------------------------------------
Double_t AliITSCalibration::ReverseBiasCurrent(Double_t temp,
                                            Double_t revBiasCurT1,
                                            Double_t tempT1,
                                            Double_t energy)const{
    // Computes the temperature dependance of the reverse bias current
    // of Si detectors. From the Particle Data
    // Book, 28.8 Silicon semiconductor detectors equation 28.21 (2004)
    // Physics Letters B "Review of Particle Physics" Volume 592, Issue 1-4
    // July 15 2004, ISSN 0370-2693 page 263.
    // Inputs:
    //    Double_t temp         The temperature at which the current is wanted
    //    Double_t revBiasCurT1 The reference bias current at temp T1
    //    Double_t tempT1       The temperature correstponding to revBiasCurT1 
    //    Double_t energy=1.2   Some energy [eV]
    // Output:
    //    none.
    // Return:
    //    The reverse bias current at the tempeature temp.
    const Double_t kBoltz = 8.617343E-5; //[eV/K]

    return revBiasCurT1*(temp*temp/(tempT1*tempT1))*
        TMath::Exp(-0.5*energy*(tempT1-temp)/(kBoltz*tempT1*temp));
}
//----------------------------------------------------------------------
void AliITSCalibration::Print(ostream *os) const {
  // Standard output format for this class.
  // Inputs:
    *os << fdv << " " << fN << " " << fT << " ";
    *os << fGeVcharge;    
  //    printf("%-10.6e  %-10.6e %-10.6e %-10.6e \n",fdv,fN,fT,fGeVcharge);
    return;
}
//----------------------------------------------------------------------
void AliITSCalibration::Read(istream *is) {
  // Standard input format for this class.
  // Inputs:
  //    ostream *is  Pointer to the output stream
  // Outputs:
  //    none:
  // Return:
  //    none.

    *is >> fdv >> fN >> fT >> fGeVcharge;
    return;
}
//----------------------------------------------------------------------

ostream &operator<<(ostream &os,AliITSCalibration &p){
  // Standard output streaming function.
  // Inputs:
  //    ostream *os  Pointer to the output stream
  // Outputs:
  //    none:
  // Return:
  //    none.

    p.Print(&os);
    return os;
}

//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSCalibration &r){
  // Standard input streaming function.
  // Inputs:
  //    ostream *os  Pointer to the output stream
  // Outputs:
  //    none:
  // Return:
  //    none.

    r.Read(&is);
    return is;
}
//----------------------------------------------------------------------
