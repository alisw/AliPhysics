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

/* $Id: AliITSUSimuParam.cxx 48165 2011-03-07 17:48:57Z masera $ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class to store the parameters used in   //
// the simulation of SPD, SDD and SSD detectors                  //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////
#include "AliITSUSimuParam.h"
#include "AliLog.h"


const Double_t  AliITSUSimuParam::fgkPixUpgBiasVoltageDefault = 18.182;
const Double_t  AliITSUSimuParam::fgkPixUpgThreshDefault = 3000.;
const Double_t  AliITSUSimuParam::fgkPixUpgThrSigmaDefault = 250.;
const UInt_t    AliITSUSimuParam::fgkPixUpgCouplingOptDefault = AliITSUSimuParam::kOldCouplingPixUpg;
const Double_t  AliITSUSimuParam::fgkPixUpgCouplColDefault = 0.;
const Double_t  AliITSUSimuParam::fgkPixUpgCouplRowDefault = 0.055;
const Double_t  AliITSUSimuParam::fgkPixUpgEccDiffDefault = 0.85;
const Double_t  AliITSUSimuParam::fgkPixUpgLorentzHoleWeightDefault = 1.0;
const Double_t  AliITSUSimuParam::fgkGeVtoChargeDefault = 3.6e-9;
const Double_t  AliITSUSimuParam::fgkDOverVDefault = 0.000375;
const Double_t  AliITSUSimuParam::fgkTDefault = 300;

const Double_t  AliITSUSimuParam::fgkNsigmasDefault = 3.;
const Int_t    AliITSUSimuParam::fgkNcompsDefault = 121;

ClassImp(AliITSUSimuParam)

//______________________________________________________________________
AliITSUSimuParam::AliITSUSimuParam()
:  fGeVcharge(fgkGeVtoChargeDefault)
  ,fDOverV(fgkDOverVDefault)
  ,fT(fgkTDefault)
  //
  ,fNPixUpg(0)
  ,fPixUpgCouplOpt(kOldCouplingPixUpg)
  ,fPixUpgCouplCol(fgkPixUpgCouplColDefault)
  ,fPixUpgCouplRow(fgkPixUpgCouplRowDefault)
  ,fPixUpgEccDiff(fgkPixUpgEccDiffDefault)
  ,fPixUpgLorentzDrift(kTRUE)
  ,fPixUpgLorentzHoleWeight(fgkPixUpgLorentzHoleWeightDefault)
  ,fPixUpgAddNoisyFlag(kFALSE)
  ,fPixUpgRemoveDeadFlag(kFALSE)
  //
  ,fPixUpgThreshDef(fgkPixUpgThreshDefault)
  ,fPixUpgThrSigmaDef(fgkPixUpgThrSigmaDefault)
  ,fPixUpgBiasVoltageDef(fgkPixUpgBiasVoltageDefault)
  ,fPixUpgNoiseDef(0)
  ,fPixUpgBaselineDef(0)
  //
  ,fPixUpgThresh(0)
  ,fPixUpgThrSigma(0)
  ,fPixUpgBiasVoltage(0)
  ,fPixUpgSigma(0)
  ,fPixUpgNoise(0)
  ,fPixUpgBaseline(0)
{  
  // default constructor
}

//______________________________________________________________________
AliITSUSimuParam::AliITSUSimuParam(UInt_t nPixUpg)
  :fGeVcharge(fgkGeVtoChargeDefault)
  ,fDOverV(fgkDOverVDefault)
  ,fT(fgkTDefault)
    //
  ,fNPixUpg(nPixUpg)
  ,fPixUpgCouplOpt(kOldCouplingPixUpg)
  ,fPixUpgCouplCol(fgkPixUpgCouplColDefault)
  ,fPixUpgCouplRow(fgkPixUpgCouplRowDefault)
  ,fPixUpgEccDiff(fgkPixUpgEccDiffDefault)
  ,fPixUpgLorentzDrift(kTRUE)
  ,fPixUpgLorentzHoleWeight(fgkPixUpgLorentzHoleWeightDefault)
  ,fPixUpgAddNoisyFlag(kFALSE)
  ,fPixUpgRemoveDeadFlag(kFALSE)
  //
  ,fPixUpgThreshDef(fgkPixUpgThreshDefault)
  ,fPixUpgThrSigmaDef(fgkPixUpgThrSigmaDefault)
  ,fPixUpgBiasVoltageDef(fgkPixUpgBiasVoltageDefault)
  ,fPixUpgNoiseDef(0)
  ,fPixUpgBaselineDef(0)
  //
  ,fPixUpgThresh(0)
  ,fPixUpgThrSigma(0)
  ,fPixUpgBiasVoltage(0)
  ,fPixUpgSigma(0)
  ,fPixUpgNoise(0)
  ,fPixUpgBaseline(0)
{  
  // regular constructor
  if (fNPixUpg>0) {
    fPixUpgBiasVoltage = new Double_t[fNPixUpg];
    fPixUpgThresh      = new Double_t[fNPixUpg];
    fPixUpgThrSigma    = new Double_t[fNPixUpg];
    fPixUpgNoise       = new Double_t[fNPixUpg];
    fPixUpgBaseline    = new Double_t[fNPixUpg];
  }
  SetPixUpgThreshold(fgkPixUpgThreshDefault,fgkPixUpgThrSigmaDefault);
  SetPixUpgNoise(0.,0.);
  SetPixUpgBiasVoltage(fgkPixUpgBiasVoltageDefault);
  //
}

//______________________________________________________________________
AliITSUSimuParam::AliITSUSimuParam(const AliITSUSimuParam &simpar)
  :TObject(simpar)
  ,fGeVcharge(simpar.fGeVcharge)
  ,fDOverV(simpar.fDOverV)
  ,fT(simpar.fT)
   //
  ,fNPixUpg(simpar.fNPixUpg)
  ,fPixUpgCouplOpt(simpar.fPixUpgCouplOpt)
  ,fPixUpgCouplCol(simpar.fPixUpgCouplCol)
  ,fPixUpgCouplRow(simpar.fPixUpgCouplRow)
  ,fPixUpgEccDiff(simpar.fPixUpgEccDiff)
  ,fPixUpgLorentzDrift(simpar.fPixUpgLorentzDrift)
  ,fPixUpgLorentzHoleWeight(simpar.fPixUpgLorentzHoleWeight)
  ,fPixUpgAddNoisyFlag(simpar.fPixUpgAddNoisyFlag)
  ,fPixUpgRemoveDeadFlag(simpar.fPixUpgRemoveDeadFlag)
   //
  ,fPixUpgThreshDef(simpar.fPixUpgThreshDef)
  ,fPixUpgThrSigmaDef(simpar.fPixUpgThrSigmaDef)
  ,fPixUpgBiasVoltageDef(simpar.fPixUpgBiasVoltageDef)
  ,fPixUpgNoiseDef(simpar.fPixUpgNoiseDef)
  ,fPixUpgBaselineDef(simpar.fPixUpgBaselineDef)
  //
  ,fPixUpgThresh(0)
  ,fPixUpgThrSigma(0)
  ,fPixUpgBiasVoltage(0)
  ,fPixUpgSigma(0)
  ,fPixUpgNoise(0)
  ,fPixUpgBaseline(0)   
   //
{
  // copy constructor
  if (fNPixUpg) {
    fPixUpgBiasVoltage = new Double_t[fNPixUpg];
    fPixUpgThresh      = new Double_t[fNPixUpg];
    fPixUpgThrSigma    = new Double_t[fNPixUpg];
    fPixUpgNoise       = new Double_t[fNPixUpg];
    fPixUpgBaseline    = new Double_t[fNPixUpg];
  }
  for (Int_t i=fNPixUpg;i--;) {
    fPixUpgBiasVoltage[i] = simpar.fPixUpgBiasVoltage[i];
    fPixUpgThresh[i]      = simpar.fPixUpgThresh[i];
    fPixUpgThrSigma[i]    = simpar.fPixUpgThrSigma[i];
    fPixUpgNoise[i]       = simpar.fPixUpgNoise[i];
    fPixUpgBaseline[i]    = simpar.fPixUpgBaseline[i];
  }
  //
}

//______________________________________________________________________
AliITSUSimuParam& AliITSUSimuParam::operator=(const AliITSUSimuParam& source)
{
  // Assignment operator. 
  if (this==&source) return *this;
  this->~AliITSUSimuParam();
  new(this) AliITSUSimuParam(source);
  return *this;
  //
}

//______________________________________________________________________
AliITSUSimuParam::~AliITSUSimuParam() 
{
  // destructor
  delete[] fPixUpgBiasVoltage;
  delete[] fPixUpgThresh;
  delete[] fPixUpgThrSigma;
  delete[] fPixUpgNoise;
  delete[] fPixUpgBaseline;
}

//________________________________________________________________________
void AliITSUSimuParam::Print(Option_t *) const
{
  // Dump all parameters
  Dump();
}

//_______________________________________________________________________
Double_t AliITSUSimuParam::ApplyPixUpgBaselineAndNoise(UInt_t mod) const 
{
  // generate random noise 
  double base,noise;
  if (mod>=fNPixUpg) {
    if (fNPixUpg>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPixUpg));}
    base = fPixUpgBaselineDef;
    noise = fPixUpgNoiseDef;
  }
  else {
    base  = fPixUpgBaseline[mod];
    noise = fPixUpgNoise[mod];    
  }
  return base+noise*gRandom->Gaus();
}

//_______________________________________________________________________
Double_t AliITSUSimuParam::CalcProbNoiseOverThreshold(UInt_t mod) const 
{
  // calculate probability of noise exceeding the threshold
  double base,noise,thresh;
  if (mod>=fNPixUpg) {
    if (fNPixUpg>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPixUpg));}
    base   = fPixUpgBaselineDef;
    noise  = fPixUpgNoiseDef;
    thresh = fPixUpgThreshDef;
  }
  else {
    base   = fPixUpgBaseline[mod];
    noise  = fPixUpgNoise[mod];
    thresh = fPixUpgThresh[mod];
  }
  if (noise<1e-12) {
    if (base>thresh) return 1;
    else             return 0;
  }
  return CalcProbNoiseOverThreshold(base, noise, thresh);
}

//_______________________________________________________________________
void AliITSUSimuParam::SetPixUpgThreshold(Double_t thresh, Double_t sigma, int mod)
{
  // set threshold params
  if (mod<0) {
    fPixUpgThreshDef = thresh;
    fPixUpgThrSigmaDef  = sigma;
    for (int i=fNPixUpg;i--;) {
      fPixUpgThresh[i] = thresh;
      fPixUpgThrSigma[i]  = sigma;
    }
  }
  else if (mod>=(int)fNPixUpg) {
    if (fNPixUpg>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPixUpg));}
    fPixUpgThreshDef = thresh;
    fPixUpgThrSigmaDef  = sigma;
  }
  else {
    fPixUpgThresh[mod] = thresh;
    fPixUpgThrSigma[mod]  = sigma;
  }
  //
}

//_______________________________________________________________________
Double_t AliITSUSimuParam::GetPixUpgThreshold(UInt_t mod) const
{
  // obtain threshold
  if (mod>=fNPixUpg) {
    if (fNPixUpg>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPixUpg));}
    return fPixUpgThreshDef;
  }
  else return fPixUpgThresh[mod];
}

//_______________________________________________________________________
void AliITSUSimuParam::GetPixUpgThreshold(UInt_t mod, Double_t &thresh, Double_t &sigma) const
{
  // obtain thresholds
  if (mod>=fNPixUpg) {
    if (fNPixUpg>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPixUpg));}
    thresh = fPixUpgThreshDef;
    sigma  = fPixUpgThrSigmaDef;
  }
  else {
    thresh = fPixUpgThresh[mod];
    sigma  = fPixUpgThrSigma[mod];
  }
}

//_______________________________________________________________________
void AliITSUSimuParam::SetPixUpgBiasVoltage(Double_t val, int mod)
{
  // set threshold params
  if (mod<0) {
    fPixUpgBiasVoltageDef = val;
    for (int i=fNPixUpg;i--;) fPixUpgBiasVoltage[i] = val;
  }
  else if (mod>=(int)fNPixUpg) {
    if (fNPixUpg>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPixUpg));}
    fPixUpgBiasVoltageDef = val;
  }
  else fPixUpgBiasVoltage[mod] = val;
  //
}

//_______________________________________________________________________
Double_t AliITSUSimuParam::GetPixUpgBiasVoltage(UInt_t mod) const
{
  // obtain threshold
  if (mod>=fNPixUpg) {
    if (fNPixUpg>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPixUpg));}
    return fPixUpgBiasVoltageDef;
  }
  else return fPixUpgBiasVoltage[mod];
}

//_______________________________________________________________________
void AliITSUSimuParam::SetPixUpgNoise(Double_t noise, Double_t baseline, int mod)
{
  // set noise params
  if (mod<0) {
    fPixUpgNoiseDef = noise;
    fPixUpgBaselineDef  = baseline;
    for (int i=fNPixUpg;i--;) {
      fPixUpgNoise[i] = noise;
      fPixUpgBaseline[i]  = baseline;
    }
  }
  else if (mod>=(int)fNPixUpg) {
    if (fNPixUpg>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPixUpg));}
    fPixUpgNoiseDef = noise;
    fPixUpgBaselineDef  = baseline;
  }
  else {
    fPixUpgNoise[mod] = noise;
    fPixUpgBaseline[mod]  = baseline;
  }
  //
}

//_______________________________________________________________________
void AliITSUSimuParam::GetPixUpgNoise(UInt_t mod, Double_t &noise, Double_t &baseline) const
{
  // obtain noise
  if (mod>=fNPixUpg) {
    if (fNPixUpg>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPixUpg));}
    noise     = fPixUpgNoiseDef;
    baseline  = fPixUpgBaselineDef;
  }
  else {
    noise     = fPixUpgNoise[mod];
    baseline  = fPixUpgBaseline[mod];
  }
}

//_______________________________________________________________________
void AliITSUSimuParam::SetPixUpgCouplingOption(UInt_t opt)
{
  // set coupling option
  if (opt>=kMaxCouplingOptPixUpg) AliFatal(Form("Coupling option %d should be less than %d",opt,kMaxCouplingOptPixUpg));
  fPixUpgCouplOpt = opt;
}


//______________________________________________________________________
Double_t AliITSUSimuParam::LorentzAngleHole(Double_t B) const 
{
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
Double_t AliITSUSimuParam::LorentzAngleElectron(Double_t B) const 
{
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
Double_t AliITSSimuParam::SigmaDiffusion3D(Double_t l) const 
{
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
Double_t AliITSUSimuParam::SigmaDiffusion2D(Double_t l) const 
{
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
Double_t AliITSSimuParam::SigmaDiffusion1D(Double_t l) const 
{
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
