/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class to store the parameters used in   //
// the simulation of SPD, SDD and SSD detectors                  //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSSimuParam.h"
#include <TMath.h>

const Float_t  AliITSSimuParam::fgkSPDBiasVoltageDefault = 18.182;
const Double_t AliITSSimuParam::fgkSPDThreshDefault = 3000.;
const Double_t AliITSSimuParam::fgkSPDSigmaDefault = 250.;
const TString  AliITSSimuParam::fgkSPDCouplingOptDefault = "old";
const Double_t AliITSSimuParam::fgkSPDCouplColDefault = 0.;
const Double_t AliITSSimuParam::fgkSPDCouplRowDefault = 0.055;
const Float_t  AliITSSimuParam::fgkSPDEccDiffDefault = 0.85;
const Float_t  AliITSSimuParam::fgkSPDLorentzHoleWeightDefault = 1.0;
const Float_t  AliITSSimuParam::fgkSDDDiffCoeffDefault = 3.23;
const Float_t  AliITSSimuParam::fgkSDDDiffCoeff1Default = 30.;
const Float_t  AliITSSimuParam::fgkSDDJitterErrorDefault = 20.; // 20 um from beam test 2001
const Float_t  AliITSSimuParam::fgkSDDDynamicRangeDefault = 1400./2.5; // mV/MOhm = nA
const Int_t    AliITSSimuParam::fgkSDDMaxAdcDefault = 1024;
const Float_t  AliITSSimuParam::fgkSDDChargeLossDefault = 0.;
const Float_t  AliITSSimuParam::fgkSDDTrigDelayDefault = 54.3;
const Float_t  AliITSSimuParam::fgkSDDMapPrecDefault = 20.; // 20 um from laser tests
const Float_t  AliITSSimuParam::fgkSDDkeVtoADCDefault = 3.42;
const Double_t AliITSSimuParam::fgkSSDCouplingPRDefault = 0.01;
const Double_t AliITSSimuParam::fgkSSDCouplingPLDefault = 0.01;
const Double_t AliITSSimuParam::fgkSSDCouplingNRDefault = 0.01;
const Double_t AliITSSimuParam::fgkSSDCouplingNLDefault = 0.01;
const Int_t    AliITSSimuParam::fgkSSDZSThresholdDefault = 3;

const Float_t AliITSSimuParam::fgkNsigmasDefault = 3.;
const Int_t AliITSSimuParam::fgkNcompsDefault = 121;

ClassImp(AliITSSimuParam)

//______________________________________________________________________
AliITSSimuParam::AliITSSimuParam():
  TObject(),
  fGeVcharge(0.),
  fDOverV(0.),
//fSPDBiasVoltage(fgkSPDBiasVoltageDefault),
//fSPDThresh(fgkSPDThreshDefault),
//fSPDSigma(fgkSPDSigmaDefault),
  fSPDCouplOpt(0),
  fSPDCouplCol(fgkSPDCouplColDefault),
  fSPDCouplRow(fgkSPDCouplRowDefault),
  fSPDEccDiff(0.),
  fSPDLorentzDrift(kTRUE),
  fSPDLorentzHoleWeight(fgkSPDLorentzHoleWeightDefault),
  fSPDAddNoisyFlag(kFALSE),
  fSPDRemoveDeadFlag(kFALSE),
  fSDDElectronics(0),
  fSDDDiffCoeff(0.),
  fSDDDiffCoeff1(0.),
  fSDDJitterError(fgkSDDJitterErrorDefault),
  fSDDDynamicRange(fgkSDDDynamicRangeDefault),
  fSDDMaxAdc(0.),
  fSDDChargeLoss(fgkSDDChargeLossDefault),
  fSDDTrigDelay(fgkSDDTrigDelayDefault),
  fSDDMapPrec(fgkSDDMapPrecDefault),
  fSDDkeVtoADC(fgkSDDkeVtoADCDefault),
  fSDDRawFormat(7),
  fSSDLorentzDrift(kTRUE),
  fSSDCouplingPR(0),
  fSSDCouplingPL(0),
  fSSDCouplingNR(0),
  fSSDCouplingNL(0),
  fSSDZSThreshold(fgkSSDZSThresholdDefault),
  fNsigmas(fgkNsigmasDefault),
  fNcomps(fgkNcompsDefault),
  fGaus(),
  fN(0.),
  fT(300.)
{  
  // default constructor
  SetSPDBiasVoltageAll(fgkSPDBiasVoltageDefault);
  SetSPDThresholdsAll(fgkSPDThreshDefault,fgkSPDSigmaDefault);
  SetSPDNoiseAll(0,0);
  SetGeVToCharge();
  SetDistanceOverVoltage();
  SetSPDCouplingOption(fgkSPDCouplingOptDefault);
  SetSPDSigmaDiffusionAsymmetry(fgkSPDEccDiffDefault);
  SetSDDElectronics();
  SetSDDDiffCoeff(fgkSDDDiffCoeffDefault,fgkSDDDiffCoeff1Default);
  SetSDDMaxAdc((Double_t)fgkSDDMaxAdcDefault);
  SetSSDCouplings(fgkSSDCouplingPRDefault,fgkSSDCouplingPLDefault,fgkSSDCouplingNRDefault,fgkSSDCouplingNLDefault);
  fSPDHitStrobe[0]=-1; // 100 ns before the collision time ( 300 ns readout strobe in total)
  fSPDHitStrobe[1]=2; // 200 ns after the collision time  (300 ns readout strobe in total)
  fSPDFoStrobe[0]=0;  // coincidence with collision time
  fSPDFoStrobe[1]=1;  // 100 ns after the collision time
}
//______________________________________________________________________
AliITSSimuParam::AliITSSimuParam(const AliITSSimuParam &simpar):
  TObject(),
  fGeVcharge(simpar.fGeVcharge),
  fDOverV(simpar.fDOverV),
  //fSPDBiasVoltage(simpar.fSPDBiasVoltage),
  //fSPDThresh(simpar.fSPDThresh),
  //fSPDSigma(simpar.fSPDSigma),
  fSPDCouplOpt(simpar.fSPDCouplOpt),
  fSPDCouplCol(simpar.fSPDCouplCol),
  fSPDCouplRow(simpar.fSPDCouplRow),
  fSPDEccDiff(simpar.fSPDEccDiff),
  fSPDLorentzDrift(simpar.fSPDLorentzDrift),
  fSPDLorentzHoleWeight(simpar.fSPDLorentzHoleWeight),
  fSPDAddNoisyFlag(simpar.fSPDAddNoisyFlag),
  fSPDRemoveDeadFlag(simpar.fSPDRemoveDeadFlag),
  fSDDElectronics(simpar.fSDDElectronics),
  fSDDDiffCoeff(simpar.fSDDDiffCoeff),
  fSDDDiffCoeff1(simpar.fSDDDiffCoeff1),
  fSDDJitterError(simpar.fSDDJitterError),
  fSDDDynamicRange(simpar.fSDDDynamicRange),
  fSDDMaxAdc(simpar.fSDDMaxAdc),
  fSDDChargeLoss(simpar.fSDDChargeLoss),
  fSDDTrigDelay(simpar.fSDDTrigDelay),
  fSDDMapPrec(simpar.fSDDMapPrec),
  fSDDkeVtoADC(simpar.fSDDkeVtoADC),
  fSDDRawFormat(simpar.fSDDRawFormat),
  fSSDLorentzDrift(simpar.fSSDLorentzDrift),
  fSSDCouplingPR(simpar.fSSDCouplingPR),
  fSSDCouplingPL(simpar.fSSDCouplingPL),
  fSSDCouplingNR(simpar.fSSDCouplingNR),
  fSSDCouplingNL(simpar.fSSDCouplingNL),
  fSSDZSThreshold(simpar.fSSDZSThreshold),
  fNsigmas(simpar.fNsigmas),
  fNcomps(simpar.fNcomps),
  fGaus(),
  fN(simpar.fN),
  fT(simpar.fT){
  // copy constructor
  for (Int_t i=0;i<240;i++) {
    fSPDBiasVoltage[i]=simpar.fSPDBiasVoltage[i];
    fSPDThresh[i]=simpar.fSPDThresh[i];
    fSPDSigma[i]=simpar.fSPDSigma[i];
    fSPDNoise[i]=simpar.fSPDNoise[i];
    fSPDBaseline[i]=simpar.fSPDBaseline[i];
  }
 for (Int_t j=0; j<2; j++){
  fSPDHitStrobe[j]=simpar.fSPDHitStrobe[j]; 
  fSPDFoStrobe[j]=simpar.fSPDFoStrobe[j]; 
 }
}

//______________________________________________________________________
AliITSSimuParam& AliITSSimuParam::operator=(const AliITSSimuParam& source){
  // Assignment operator. 
  this->~AliITSSimuParam();
  new(this) AliITSSimuParam(source);
  return *this;
  
}


//______________________________________________________________________
AliITSSimuParam::~AliITSSimuParam() {
  // destructor
  if(fGaus) delete fGaus;
}
//________________________________________________________________________
void AliITSSimuParam::SetNLookUp(Int_t p1){
  // Set number of sigmas over which cluster disintegration is performed
  fNcomps=p1;
  if (fGaus) delete fGaus;
  fGaus = new TArrayF(fNcomps+1);
  for(Int_t i=0; i<=fNcomps; i++) {
    Float_t x = -fNsigmas + (2.*i*fNsigmas)/(fNcomps-1);
    (*fGaus)[i] = exp(-((x*x)/2));
  }
}
//________________________________________________________________________
void AliITSSimuParam::PrintParameters() const{
  // Dump all parameters
  printf("GeVToCharge               = %G\n",fGeVcharge);
  printf("DistanveOverVoltage       = %f \n",fDOverV);
  printf("\n");
  printf("=====  SPD parameters  =====\n");
  printf("Bias Voltage              = %f \n",fSPDBiasVoltage[0]);
  printf("Threshold and sigma       = %f %f\n",fSPDThresh[0],fSPDSigma[0]);
  printf("Coupling Option           = %s\n",fSPDCouplOpt.Data());
  printf("Coupling value (column)   = %f\n",fSPDCouplCol);
  printf("Coupling value (row)      = %f\n",fSPDCouplRow);
  printf("Eccentricity in diffusion = %f\n",fSPDEccDiff);
  printf("Flag to add Lorentz Drift = %d\n",fSPDLorentzDrift);
  printf("Weight of Holes in Lor.Drift = %f\n",fSPDLorentzHoleWeight);
  printf("Flag to add noisy         = %d\n",fSPDAddNoisyFlag);
  printf("Flag to remove dead       = %d\n",fSPDRemoveDeadFlag);
  printf("Hit Strobe params         = %d %d\n",fSPDHitStrobe[0],fSPDHitStrobe[1]);
  printf("FO  Strobe params         = %d %d\n",fSPDFoStrobe[0],fSPDFoStrobe[1]);
  printf("\n");
  printf("=====  SDD parameters  =====\n");
  printf("Electronic chips          = %d\n",fSDDElectronics);
  printf("Diffusion Coefficients    = %f %f\n",fSDDDiffCoeff,fSDDDiffCoeff1);
  printf("Jitter Error              = %f um\n",fSDDJitterError);
  printf("Dynamic Range             = %f\n",fSDDDynamicRange);
  printf("Max. ADC                  = %f\n",fSDDMaxAdc);
  printf("Charge Loss               = %f\n",fSDDChargeLoss);  
  printf("Trigger Delay (ns)        = %f\n",fSDDTrigDelay);  
  printf("Smear from map (um)       = %f\n",fSDDMapPrec);
  printf("keV->ADC conv. fact.        = %f\n",fSDDkeVtoADC);
  printf("Raw Data Format           = %d\n",fSDDRawFormat);  
  printf("\n");
  printf("=====  SSD parameters  =====\n");
  printf("Flag to add Lorentz Drift = %d\n",fSSDLorentzDrift);
  printf("Coupling PR               = %f\n",fSSDCouplingPR);
  printf("Coupling PL               = %f\n",fSSDCouplingPL);
  printf("Coupling NR               = %f\n",fSSDCouplingNR);
  printf("Coupling NL               = %f\n",fSSDCouplingNL);
  printf("Zero Supp threshold       = %d\n",fSSDZSThreshold);
}
//______________________________________________________________________
Double_t AliITSSimuParam::MobilityElectronSiEmp() const {
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
Double_t AliITSSimuParam::MobilityHoleSiEmp() const {
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
Double_t AliITSSimuParam::DiffusionCoefficientElectron() const {
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
Double_t AliITSSimuParam::DiffusionCoefficientHole() const {
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
Double_t AliITSSimuParam::LorentzAngleHole(Double_t B) const {
   // Computes the Lorentz angle for electrons in Si
   // Input: magnetic Field in KGauss
   // Output: Lorentz angle in radians (positive if Bz is positive)
   // Main Reference: NIM A 497 (2003) 389–396.
   // "An algorithm for calculating the Lorentz angle in silicon detectors", V. Bartsch et al.
   //
   const Double_t krH=0.70; // Hall scattering factor for Hole
   const Double_t kT0  = 300.;       // reference Temperature (degree K).
   const Double_t kmulow0 = 470.5;   // cm^2/Volt-sec
   const Double_t keT0 = -2.5;       // Power of Temp.
   const Double_t beta0 = 1.213;     // beta coeff. at T0=300K
   const Double_t keT1 = 0.17;       // Power of Temp. for beta
   const Double_t kvsat0 = 8.37E+06; // saturated velocity at T0=300K (cm/sec)
   const Double_t keT2 = 0.52;       // Power of Temp. for vsat
   Double_t tT = fT;
   Double_t eE= 1./fDOverV;
   Double_t muLow=kmulow0*TMath::Power(tT/kT0,keT0);
   Double_t beta=beta0*TMath::Power(tT/kT0,keT1);
   Double_t vsat=kvsat0*TMath::Power(tT/kT0,keT2);
   Double_t mu=muLow/TMath::Power(1+TMath::Power(muLow*eE/vsat,beta),1/beta);
   Double_t angle=TMath::ATan(krH*mu*B*1.E-05); // Conversion Factor
   return angle;
}
//______________________________________________________________________
Double_t AliITSSimuParam::LorentzAngleElectron(Double_t B) const {
   // Computes the Lorentz angle for electrons in Si
   // Input: magnetic Field in KGauss
   // Output: Lorentz angle in radians (positive if Bz is positive)
   // Main Reference: NIM A 497 (2003) 389–396.
   // "An algorithm for calculating the Lorentz angle in silicon detectors", V. Bartsch et al.
   //
   const Double_t krH=1.15; // Hall scattering factor for Electron
   const Double_t kT0  = 300.;       // reference Temperature (degree K).
   const Double_t kmulow0 = 1417.0;  // cm^2/Volt-sec
   const Double_t keT0 = -2.2;       // Power of Temp.
   const Double_t beta0 = 1.109;     // beta coeff. at T0=300K
   const Double_t keT1 = 0.66;       // Power of Temp. for beta
   const Double_t kvsat0 = 1.07E+07; // saturated velocity at T0=300K (cm/sec)
   const Double_t keT2 = 0.87;       // Power of Temp. for vsat
   Double_t tT = fT;
   Double_t eE= 1./fDOverV;
   Double_t muLow=kmulow0*TMath::Power(tT/kT0,keT0);
   Double_t beta=beta0*TMath::Power(tT/kT0,keT1);
   Double_t vsat=kvsat0*TMath::Power(tT/kT0,keT2);
   Double_t mu=muLow/TMath::Power(1+TMath::Power(muLow*eE/vsat,beta),1/beta);
   Double_t angle=TMath::ATan(krH*mu*B*1.E-05);
   return angle;
}
//______________________________________________________________________
Double_t AliITSSimuParam::SpeedElectron() const {
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

    return m/fDOverV;  // [cm/sec]
}
//______________________________________________________________________
Double_t AliITSSimuParam::SpeedHole() const {
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

    return m/fDOverV;  // [cm/sec]
}
//______________________________________________________________________
Double_t AliITSSimuParam::SigmaDiffusion3D(Double_t l) const {
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

    return TMath::Sqrt(kcon*fT*fDOverV*l);  // [cm]
}
//______________________________________________________________________
Double_t AliITSSimuParam::SigmaDiffusion2D(Double_t l) const {
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

    return TMath::Sqrt(kcon*fT*fDOverV*l);  // [cm]
}
//______________________________________________________________________
Double_t AliITSSimuParam::SigmaDiffusion1D(Double_t l) const {
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

    return TMath::Sqrt(kcon*fT*fDOverV*l);  // [cm]
}
//----------------------------------------------------------------------
Double_t AliITSSimuParam::DepletedRegionThicknessA(Double_t dopCons,
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
Double_t AliITSSimuParam::DepletedRegionThicknessB(Double_t resist,
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
Double_t AliITSSimuParam::ReverseBiasCurrent(Double_t temp,
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
//______________________________________________________________________
 void   AliITSSimuParam::SPDThresholds(const Int_t mod, Double_t& thresh, Double_t& sigma) const {
   // Get SPD threshold values
    if(mod<0 || mod>239) {
       thresh=0;
       sigma=0; 
       return;
     } 
     thresh=fSPDThresh[mod];
     sigma=fSPDSigma[mod];
     return;
}
//_______________________________________________________________________
 void   AliITSSimuParam::SPDNoise(const Int_t mod,Double_t &noise, Double_t &baseline) const {
   //Get SPD noise and baseline values
     if(mod<0 || mod>239) {
       noise=0;
       baseline=0; 
       return;
     } 
     noise=fSPDNoise[mod];
     baseline=fSPDBaseline[mod];
     return;
}

