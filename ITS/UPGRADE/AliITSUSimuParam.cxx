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
#include "AliParamList.h"

using namespace TMath;


const Double_t  AliITSUSimuParam::fgkPixBiasVoltageDefault = 18.182;
const Double_t  AliITSUSimuParam::fgkPixThreshDefault = 100.;
const Double_t  AliITSUSimuParam::fgkPixThrSigmaDefault = 10.;
const Double_t  AliITSUSimuParam::fgkPixMinElToAddDefault = 2.;
const UInt_t    AliITSUSimuParam::fgkPixCouplingOptDefault = AliITSUSimuParam::kOldCouplingPix;
const Double_t  AliITSUSimuParam::fgkPixCouplColDefault = 0.;
const Double_t  AliITSUSimuParam::fgkPixCouplRowDefault = 0.055;
const Double_t  AliITSUSimuParam::fgkPixEccDiffDefault = 0.85;
const Double_t  AliITSUSimuParam::fgkPixLorentzHoleWeightDefault = 1.0;
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
  ,fNPix(0)
  ,fPixCouplOpt(kOldCouplingPix)
  ,fPixCouplCol(fgkPixCouplColDefault)
  ,fPixCouplRow(fgkPixCouplRowDefault)
  ,fPixLorentzDrift(kTRUE)
  ,fPixLorentzHoleWeight(fgkPixLorentzHoleWeightDefault)
  ,fPixAddNoisyFlag(kFALSE)
  ,fPixRemoveDeadFlag(kFALSE)
  //
  ,fPixThreshDef(fgkPixThreshDefault)
  ,fPixThrSigmaDef(fgkPixThrSigmaDefault)
  ,fPixBiasVoltageDef(fgkPixBiasVoltageDefault)
  ,fPixNoiseDef(0)
  ,fPixBaselineDef(0)
  ,fPixMinElToAddDef(fgkPixMinElToAddDefault)
  //
  ,fPixThresh(0)
  ,fPixThrSigma(0)
  ,fPixBiasVoltage(0)
  ,fPixSigma(0)
  ,fPixNoise(0)
  ,fPixBaseline(0)
  ,fRespFunParam(0)
{  
  // default constructor
  fRespFunParam.SetOwner(kTRUE);
}

//______________________________________________________________________
AliITSUSimuParam::AliITSUSimuParam(UInt_t nPix)
  :fGeVcharge(fgkGeVtoChargeDefault)
  ,fDOverV(fgkDOverVDefault)
  ,fT(fgkTDefault)
    //
  ,fNPix(nPix)
  ,fPixCouplOpt(kOldCouplingPix)
  ,fPixCouplCol(fgkPixCouplColDefault)
  ,fPixCouplRow(fgkPixCouplRowDefault)
  ,fPixLorentzDrift(kTRUE)
  ,fPixLorentzHoleWeight(fgkPixLorentzHoleWeightDefault)
  ,fPixAddNoisyFlag(kFALSE)
  ,fPixRemoveDeadFlag(kFALSE)
  //
  ,fPixThreshDef(fgkPixThreshDefault)
  ,fPixThrSigmaDef(fgkPixThrSigmaDefault)
  ,fPixBiasVoltageDef(fgkPixBiasVoltageDefault)
  ,fPixNoiseDef(0)
  ,fPixBaselineDef(0)
  ,fPixMinElToAddDef(fgkPixMinElToAddDefault)
  //
  ,fPixThresh(0)
  ,fPixThrSigma(0)
  ,fPixBiasVoltage(0)
  ,fPixSigma(0)
  ,fPixNoise(0)
  ,fPixBaseline(0)
  ,fRespFunParam(0)
{  
  // regular constructor
  if (fNPix>0) {
    fPixBiasVoltage = new Double_t[fNPix];
    fPixThresh      = new Double_t[fNPix];
    fPixThrSigma    = new Double_t[fNPix];
    fPixNoise       = new Double_t[fNPix];
    fPixBaseline    = new Double_t[fNPix];
  }
  SetPixThreshold(fgkPixThreshDefault,fgkPixThrSigmaDefault);
  SetPixNoise(0.,0.);
  SetPixBiasVoltage(fgkPixBiasVoltageDefault);
  fRespFunParam.SetOwner(kTRUE);
  //
}

//______________________________________________________________________
AliITSUSimuParam::AliITSUSimuParam(const AliITSUSimuParam &simpar)
  :TObject(simpar)
  ,fGeVcharge(simpar.fGeVcharge)
  ,fDOverV(simpar.fDOverV)
  ,fT(simpar.fT)
   //
  ,fNPix(simpar.fNPix)
  ,fPixCouplOpt(simpar.fPixCouplOpt)
  ,fPixCouplCol(simpar.fPixCouplCol)
  ,fPixCouplRow(simpar.fPixCouplRow)
  ,fPixLorentzDrift(simpar.fPixLorentzDrift)
  ,fPixLorentzHoleWeight(simpar.fPixLorentzHoleWeight)
  ,fPixAddNoisyFlag(simpar.fPixAddNoisyFlag)
  ,fPixRemoveDeadFlag(simpar.fPixRemoveDeadFlag)
   //
  ,fPixThreshDef(simpar.fPixThreshDef)
  ,fPixThrSigmaDef(simpar.fPixThrSigmaDef)
  ,fPixBiasVoltageDef(simpar.fPixBiasVoltageDef)
  ,fPixNoiseDef(simpar.fPixNoiseDef)
  ,fPixBaselineDef(simpar.fPixBaselineDef)
  ,fPixMinElToAddDef(simpar.fPixMinElToAddDef)
  //
  ,fPixThresh(0)
  ,fPixThrSigma(0)
  ,fPixBiasVoltage(0)
  ,fPixSigma(0)
  ,fPixNoise(0)
  ,fPixBaseline(0)   
  ,fRespFunParam(0)
   //
{
  // copy constructor
  if (fNPix) {
    fPixBiasVoltage = new Double_t[fNPix];
    fPixThresh      = new Double_t[fNPix];
    fPixThrSigma    = new Double_t[fNPix];
    fPixNoise       = new Double_t[fNPix];
    fPixBaseline    = new Double_t[fNPix];
  }
  for (Int_t i=fNPix;i--;) {
    fPixBiasVoltage[i] = simpar.fPixBiasVoltage[i];
    fPixThresh[i]      = simpar.fPixThresh[i];
    fPixThrSigma[i]    = simpar.fPixThrSigma[i];
    fPixNoise[i]       = simpar.fPixNoise[i];
    fPixBaseline[i]    = simpar.fPixBaseline[i];
  }
  //
  for (int i=0;i<simpar.fRespFunParam.GetEntriesFast();i++) {
    AliParamList* pr = (AliParamList*)simpar.fRespFunParam[0];
    if (pr) fRespFunParam.AddLast(new AliParamList(*pr));
  }
  fRespFunParam.SetOwner(kTRUE);
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
  delete[] fPixBiasVoltage;
  delete[] fPixThresh;
  delete[] fPixThrSigma;
  delete[] fPixNoise;
  delete[] fPixBaseline;
}

//________________________________________________________________________
void AliITSUSimuParam::Print(Option_t *) const
{
  // Dump all parameters
  Dump();
}

//_______________________________________________________________________
Double_t AliITSUSimuParam::ApplyPixBaselineAndNoise(UInt_t mod) const 
{
  // generate random noise 
  double base,noise;
  if (mod>=fNPix) {
    if (fNPix>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPix));}
    base = fPixBaselineDef;
    noise = fPixNoiseDef;
  }
  else {
    base  = fPixBaseline[mod];
    noise = fPixNoise[mod];    
  }
  return base+noise*gRandom->Gaus();
}

//_______________________________________________________________________
Double_t AliITSUSimuParam::CalcProbNoiseOverThreshold(UInt_t mod) const 
{
  // calculate probability of noise exceeding the threshold
  double base,noise,thresh;
  if (mod>=fNPix) {
    if (fNPix>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPix));}
    base   = fPixBaselineDef;
    noise  = fPixNoiseDef;
    thresh = fPixThreshDef;
  }
  else {
    base   = fPixBaseline[mod];
    noise  = fPixNoise[mod];
    thresh = fPixThresh[mod];
  }
  if (noise<1e-12) {
    if (base>thresh) return 1;
    else             return 0;
  }
  return CalcProbNoiseOverThreshold(base, noise, thresh);
}

//_______________________________________________________________________
void AliITSUSimuParam::SetPixThreshold(Double_t thresh, Double_t sigma, int mod)
{
  // set threshold params
  if (mod<0) {
    fPixThreshDef = thresh;
    fPixThrSigmaDef  = sigma;
    for (int i=fNPix;i--;) {
      fPixThresh[i] = thresh;
      fPixThrSigma[i]  = sigma;
    }
  }
  else if (mod>=(int)fNPix) {
    if (fNPix>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPix));}
    fPixThreshDef = thresh;
    fPixThrSigmaDef  = sigma;
  }
  else {
    fPixThresh[mod] = thresh;
    fPixThrSigma[mod]  = sigma;
  }
  //
}

//_______________________________________________________________________
Double_t AliITSUSimuParam::GetPixThreshold(UInt_t mod) const
{
  // obtain threshold
  if (mod>=fNPix) {
    if (fNPix>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPix));}
    return fPixThreshDef;
  }
  else return fPixThresh[mod];
}

//_______________________________________________________________________
void AliITSUSimuParam::GetPixThreshold(UInt_t mod, Double_t &thresh, Double_t &sigma) const
{
  // obtain thresholds
  if (mod>=fNPix) {
    if (fNPix>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPix));}
    thresh = fPixThreshDef;
    sigma  = fPixThrSigmaDef;
  }
  else {
    thresh = fPixThresh[mod];
    sigma  = fPixThrSigma[mod];
  }
}

//_______________________________________________________________________
void AliITSUSimuParam::SetPixBiasVoltage(Double_t val, int mod)
{
  // set threshold params
  if (mod<0) {
    fPixBiasVoltageDef = val;
    for (int i=fNPix;i--;) fPixBiasVoltage[i] = val;
  }
  else if (mod>=(int)fNPix) {
    if (fNPix>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPix));}
    fPixBiasVoltageDef = val;
  }
  else fPixBiasVoltage[mod] = val;
  //
}

//_______________________________________________________________________
Double_t AliITSUSimuParam::GetPixBiasVoltage(UInt_t mod) const
{
  // obtain threshold
  if (mod>=fNPix) {
    if (fNPix>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPix));}
    return fPixBiasVoltageDef;
  }
  else return fPixBiasVoltage[mod];
}

//_______________________________________________________________________
void AliITSUSimuParam::SetPixNoise(Double_t noise, Double_t baseline, int mod)
{
  // set noise params
  if (mod<0) {
    fPixNoiseDef = noise;
    fPixBaselineDef  = baseline;
    for (int i=fNPix;i--;) {
      fPixNoise[i] = noise;
      fPixBaseline[i]  = baseline;
    }
  }
  else if (mod>=(int)fNPix) {
    if (fNPix>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPix));}
    fPixNoiseDef = noise;
    fPixBaselineDef  = baseline;
  }
  else {
    fPixNoise[mod] = noise;
    fPixBaseline[mod]  = baseline;
  }
  //
}

//_______________________________________________________________________
void AliITSUSimuParam::GetPixNoise(UInt_t mod, Double_t &noise, Double_t &baseline) const
{
  // obtain noise
  if (mod>=fNPix) {
    if (fNPix>0) {AliFatal(Form("Wrong module %d, NPidUpg=%d",mod,fNPix));}
    noise     = fPixNoiseDef;
    baseline  = fPixBaselineDef;
  }
  else {
    noise     = fPixNoise[mod];
    baseline  = fPixBaseline[mod];
  }
}

//_______________________________________________________________________
void AliITSUSimuParam::SetPixCouplingOption(UInt_t opt)
{
  // set coupling option
  if (opt>=kMaxCouplingOptPix) AliFatal(Form("Coupling option %d should be less than %d",opt,kMaxCouplingOptPix));
  fPixCouplOpt = opt;
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
  Double_t muLow=kmulow0*Power(tT/kT0,keT0);
  Double_t beta=beta0*Power(tT/kT0,keT1);
  Double_t vsat=kvsat0*Power(tT/kT0,keT2);
  Double_t mu=muLow/Power(1+Power(muLow*eE/vsat,beta),1/beta);
  Double_t angle=ATan(krH*mu*B*1.E-05); // Conversion Factor
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
  Double_t muLow=kmulow0*Power(tT/kT0,keT0);
  Double_t beta=beta0*Power(tT/kT0,keT1);
  Double_t vsat=kvsat0*Power(tT/kT0,keT2);
  Double_t mu=muLow/Power(1+Power(muLow*eE/vsat,beta),1/beta);
  Double_t angle=ATan(krH*mu*B*1.E-05);
  return angle;
}

//___________________________________________________________
const AliParamList* AliITSUSimuParam::FindRespFunParams(Int_t detId) const
{
  // find parameters list for detID
  for (int i=fRespFunParam.GetEntriesFast();i--;) {
    const AliParamList* pr = GetRespFunParams(i);
    if (int(pr->GetUniqueID())==detId) return pr;
  }
  return 0;
}

//___________________________________________________________
void AliITSUSimuParam::AddRespFunParam(AliParamList* pr)
{
  // add spread parameterization data
  fRespFunParam.AddLast(pr);
}
