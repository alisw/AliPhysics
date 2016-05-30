#ifndef ALIITSSIMUPARAM_H
#define ALIITSSIMUPARAM_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to store the parameters used in the simulation of       //
// SPD, SDD and SSD detectors                                    //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////
#include <TRandom.h>
#include<TObject.h>
#include <TString.h>
#include <TArrayF.h>

class AliITSSimuParam : public TObject {

 public:
  AliITSSimuParam();
  AliITSSimuParam(const AliITSSimuParam& simpar);
  // assignment operator 
  AliITSSimuParam& operator=(const AliITSSimuParam& source);
  ~AliITSSimuParam();

 
  void SetGeVToCharge(Double_t gc=3.6e-9){fGeVcharge = gc;}
  Double_t GetGeVToCharge() const {return fGeVcharge;}
  Double_t GeVToCharge(Double_t gev) const {return gev/fGeVcharge;}
  
  void SetDistanceOverVoltage(Double_t d,Double_t v){fDOverV = d/v;}
  void SetDistanceOverVoltage(Double_t dv=0.000375){fDOverV = dv;}
  Double_t GetDistanceOverVoltage() const {return fDOverV;}



  void    SetSPDBiasVoltageAll(Double_t bias=18.182) {for(Int_t i=0;i<240;i++) fSPDBiasVoltage[i]=bias;}
  void    SetSPDBiasVoltage(Int_t mod, Double_t bias=18.182) {if(mod<0 || mod>239) return; fSPDBiasVoltage[mod]=bias;}
  Double_t  GetSPDBiasVoltage(Int_t mod=0) const {if(mod<0 || mod>239) return 0;  return fSPDBiasVoltage[mod];}

  void   SetSPDThresholdsAll(Double_t thresh, Double_t sigma)
        {for(Int_t i=0;i<240;i++) {fSPDThresh[i]=thresh; fSPDSigma[i]=sigma;}}
  void   SetSPDThresholds(Int_t mod,Double_t thresh, Double_t sigma)
        {if(mod<0 || mod>239) return; fSPDThresh[mod]=thresh; fSPDSigma[mod]=sigma; }
  void   SPDThresholds(const Int_t mod, Double_t& thresh, Double_t& sigma) const;
  void   SetSPDNoiseAll(Double_t noise, Double_t baseline)
        {for(Int_t i=0;i<240;i++) {fSPDNoise[i]=noise; fSPDBaseline[i]=baseline;}}
  void   SetSPDNoise(Int_t mod,Double_t noise, Double_t baseline)
        {if(mod<0 || mod>239) return; fSPDNoise[mod]=noise; fSPDBaseline[mod]=baseline; }
  void   SPDNoise(const Int_t mod,Double_t &noise, Double_t &baseline) const;
  // Applies a random noise and addes the baseline
  Double_t ApplySPDBaselineAndNoise(Int_t mod=0) const 
    {if (mod<0 || mod>239) mod=0; return fSPDBaseline[mod]+fSPDNoise[mod]*gRandom->Gaus();}


  void SetSPDCouplingOption(const char *opt) {fSPDCouplOpt=opt;}
  void GetSPDCouplingOption(char *opt) const {strncpy(opt,fSPDCouplOpt.Data(),fSPDCouplOpt.Sizeof());}

  void SetSPDCouplingParam(Double_t col, Double_t row)
        {fSPDCouplCol = col; fSPDCouplRow = row;}
  void GetSPDCouplingParam(Double_t &col, Double_t &row) const
        {col = fSPDCouplCol; row = fSPDCouplRow;}

  void   SetSPDSigmaDiffusionAsymmetry(Double_t ecc) {fSPDEccDiff=ecc;}   
  void   GetSPDSigmaDiffusionAsymmetry(Double_t &ecc) const {ecc=fSPDEccDiff;}

  void    SetSPDLorentzDrift(Bool_t ison) {fSPDLorentzDrift=ison;}
  Bool_t  GetSPDLorentzDrift() const {return fSPDLorentzDrift;}
  void    SetSPDLorentzHoleWeight(Float_t weight) {fSPDLorentzHoleWeight=weight;}
  Float_t GetSPDLorentzHoleWeight() const {return fSPDLorentzHoleWeight;}
  
  void   SetSPDAddNoisyFlag(Bool_t value) {fSPDAddNoisyFlag = value;}
  Bool_t GetSPDAddNoisyFlag() const {return fSPDAddNoisyFlag;}
  void   SetSPDRemoveDeadFlag(Bool_t value) {fSPDRemoveDeadFlag = value;}
  Bool_t GetSPDRemoveDeadFlag() const {return fSPDRemoveDeadFlag;}
  
  void SetSDDElectronics(Int_t p1=1) {fSDDElectronics=p1;   }
  Int_t GetSDDElectronics()  const {return fSDDElectronics;}

  void  SetSDDDiffCoeff(Float_t p1, Float_t p2) {
      fSDDDiffCoeff=p1; fSDDDiffCoeff1=p2;}
  void  GetSDDDiffCoeff(Float_t &diff,Float_t &diff1) const {
      diff=fSDDDiffCoeff; diff1=fSDDDiffCoeff1;}

  void  SetSDDJitterError(Float_t jitter) {fSDDJitterError=jitter;}
  Float_t  GetSDDJitterError() const {return fSDDJitterError;}

  void    SetSDDDynamicRange(Double_t p1) {fSDDDynamicRange = p1;}
  Float_t GetSDDDynamicRange() const {return fSDDDynamicRange;}

  void    SetSDDMaxAdc(Double_t p1) {fSDDMaxAdc=p1;}
  Float_t GetSDDMaxAdc() const  {return fSDDMaxAdc;}

  void    SetSDDChargeLoss(Double_t p1) {fSDDChargeLoss=p1;}
  Float_t GetSDDChargeLoss() const {return fSDDChargeLoss;}

  void    SetSDDTrigDelay(Double_t p1) {fSDDTrigDelay=p1;}
  Float_t GetSDDTrigDelay() const {return fSDDTrigDelay;}

  void    SetSDDCorrMapPrecision(Double_t p1) {fSDDMapPrec=p1;}
  Float_t GetSDDCorrMapPrecision() const {return fSDDMapPrec;}

  void    SetSDDkeVtoADC(Double_t p1) {fSDDkeVtoADC=p1;}
  Float_t GetSDDkeVtoADC() const {return fSDDkeVtoADC;}

  void    SetSDDRawDataFormatCarlos() {fSDDRawFormat=7;}
  void    SetSDDRawDataFormatFixLen8bitEncoded() {fSDDRawFormat=0;}
  Char_t  GetSDDRawDataFormat() const {return fSDDRawFormat;}

  // Use Lorentz's angle
  void    SetSSDLorentzDrift(Bool_t ison) {fSSDLorentzDrift=ison;}
  Bool_t  GetSSDLorentzDrift() const {return fSSDLorentzDrift;}

  Int_t GetSSDZSThreshold() const { // ZS threshold
    return fSSDZSThreshold; }
  virtual void SetSSDZSThreshold(Int_t zsth) { fSSDZSThreshold = zsth; }
  
  void SetSSDCouplings(Double_t pr, Double_t pl, Double_t nr, Double_t nl) {
      fSSDCouplingPR=pr; fSSDCouplingPL=pl; fSSDCouplingNR=nr; fSSDCouplingNL=nl; }
  Double_t  GetSSDCouplingPR() const {return fSSDCouplingPR;}
  Double_t  GetSSDCouplingPL() const {return fSSDCouplingPL;}
  Double_t  GetSSDCouplingNR() const {return fSSDCouplingNR;}
  Double_t  GetSSDCouplingNL() const {return fSSDCouplingNL;}

  void    SetNSigmaIntegration(Double_t p1) {fNsigmas=p1;}
  Float_t GetNSigmaIntegration() const {return fNsigmas;}
  void    SetNLookUp(Int_t p1);
  Int_t   GetGausNLookUp() const {return fNcomps;}
  Float_t GetGausLookUp(Int_t i)  {
    if (!fGaus) SetNLookUp(fgkNcompsDefault);
    if(i<0 || i>=fNcomps) return 0.;return fGaus->At(i);
  }

  // Set the impurity concentrations in [#/cm^3]
  void SetImpurity(Double_t n=0.0){fN = n;}
  // Returns the impurity consentration in [#/cm^3]
  Double_t Impurity() const {return fN;}

  // Electron mobility in Si. [cm^2/(Volt Sec)]. T in degree K, N in #/cm^3
  Double_t MobilityElectronSiEmp() const ;
  // Hole mobility in Si. [cm^2/(Volt Sec)]  T in degree K, N in #/cm^3
  Double_t MobilityHoleSiEmp() const ;
  // Einstein relation for Diffusion Coefficient of Electrons. [cm^2/sec]
  //  T in degree K, N in #/cm^3
  Double_t DiffusionCoefficientElectron() const ;
  // Einstein relation for Diffusion Coefficient of Holes. [cm^2/sec]
  //  T in [degree K], N in [#/cm^3]
  Double_t DiffusionCoefficientHole() const ;
  // Electron <speed> under an applied electric field E=Volts/cm. [cm/sec]
  // d distance-thickness in [cm], v in [volts], T in [degree K],
  // N in [#/cm^3]
  Double_t SpeedElectron() const ;
  // Holes <speed> under an applied electric field E=Volts/cm. [cm/sec]
  // d distance-thickness in [cm], v in [volts], T in [degree K],
  // N in [#/cm^3]
  Double_t SpeedHole() const ;
  // Returns the Gaussian sigma == <x^2+z^2> [cm^2] due to the defusion of
  // electrons or holes through a distance l [cm] caused by an applied
  // voltage v [volt] through a distance d [cm] in any material at a
  // temperature T [degree K].
  Double_t SigmaDiffusion3D(Double_t  l) const;
  // Returns the Gaussian sigma == <x^2 +y^2+z^2> [cm^2] due to the
  // defusion of electrons or holes through a distance l [cm] caused by an
  // applied voltage v [volt] through a distance d [cm] in any material at a
  // temperature T [degree K].
  Double_t SigmaDiffusion2D(Double_t l) const;
  // Returns the Gaussian sigma == <x^2+z^2> [cm^2] due to the defusion of
  // electrons or holes through a distance l [cm] caused by an applied
  // voltage v [volt] through a distance d [cm] in any material at a
  // temperature T [degree K].
  Double_t SigmaDiffusion1D(Double_t l) const;
  // Computes the Lorentz angle for Electron and Hole, under the Magnetic field bz (in kGauss)
  Double_t LorentzAngleElectron(Double_t bz) const;
  Double_t LorentzAngleHole(Double_t bz) const;
  // Compute the thickness of the depleted region in a Si detector, version A
  Double_t DepletedRegionThicknessA(Double_t dopCons,
                                    Double_t voltage,
                                    Double_t elecCharge,
                                    Double_t voltBuiltIn=0.5)const;
  // Compute the thickness of the depleted region in a Si detector, version B
  Double_t DepletedRegionThicknessB(Double_t resist,Double_t voltage,
                                    Double_t mobility,
                                    Double_t voltBuiltIn=0.5,
                                    Double_t dielConst=1.E-12)const;
  // Computes the temperature dependance of the reverse bias current
  Double_t ReverseBiasCurrent(Double_t temp,Double_t revBiasCurT1,
  	                      Double_t tempT1,Double_t energy=1.2)const;


  void PrintParameters() const; 
  void SetSPDReadoutStrobe(Int_t win[2]){fSPDHitStrobe[0]=win[0];fSPDHitStrobe[1]=win[1];}
  void SetSPDFastOrStrobe(Int_t win[2]){fSPDFoStrobe[0]=win[0];fSPDFoStrobe[1]=win[1];}
  const Int_t* GetSPDHitStrobe() const {return fSPDHitStrobe;}
  const Int_t* GetSPDFoStrobe() const {return fSPDFoStrobe;}

 protected:

  static const Float_t fgkSPDBiasVoltageDefault;//default for fSPDBiasVoltage
  static const Double_t fgkSPDThreshDefault; //default for fThresh
  static const Double_t fgkSPDSigmaDefault; //default for fSigma
  static const TString fgkSPDCouplingOptDefault;  // type of pixel Coupling (old or new)
  static const Double_t fgkSPDCouplColDefault; //default for fSPDCouplCol
  static const Double_t fgkSPDCouplRowDefault; //default for fSPDCouplRow
  static const Float_t fgkSPDEccDiffDefault;//default for fSPDEccDiff
  static const Float_t fgkSPDLorentzHoleWeightDefault;//default for fSPDLorentzHoleWeight
  static const Float_t fgkSDDDiffCoeffDefault; // default for fSDDDiffCoeff
  static const Float_t fgkSDDDiffCoeff1Default; // default for fSDDDiffCoeff1 
  static const Float_t fgkSDDJitterErrorDefault; // default for fSDDJitterError
  static const Float_t fgkSDDDynamicRangeDefault; // default for fSDDDynamicRange
  static const Int_t fgkSDDMaxAdcDefault; // default for fSDDMaxAdc
  static const Float_t fgkSDDChargeLossDefault; // default for fSDDChargeLoss
  static const Float_t fgkSDDTrigDelayDefault; // default for fSDDTrigDelay
  static const Float_t fgkSDDMapPrecDefault; // default for fSDDTrigDelay
  static const Float_t fgkSDDkeVtoADCDefault; // default for keV->ADC conv.

  static const Double_t fgkSSDCouplingPRDefault;  // default values
  static const Double_t fgkSSDCouplingPLDefault;  // for the
  static const Double_t fgkSSDCouplingNRDefault;  // various SSD
  static const Double_t fgkSSDCouplingNLDefault;  // couplings
  static const Int_t fgkSSDZSThresholdDefault;  // default for fSSDZSThreshold

  static const Float_t fgkNsigmasDefault; //default for fNsigmas
  static const Int_t fgkNcompsDefault; //default for fNcomps

 private:
  Double_t fGeVcharge;      // Energy to ionize (free an electron) in GeV
  Double_t fDOverV;  // The parameter d/v where d is the disance over which the
                     // the potential v is applied d/v [cm/volts]

  
  Double_t fSPDBiasVoltage[240]; // Bias Voltage for the SPD
  Double_t fSPDThresh[240];      // SPD Threshold value
  Double_t fSPDSigma[240];       // SPD threshold fluctuations spread
  Double_t fSPDNoise[240];       // SPD electronic noise: sigma
  Double_t fSPDBaseline[240];    // SPD electronic noise: baseline
  TString  fSPDCouplOpt;    // SPD Coupling Option
  Double_t fSPDCouplCol;    // SPD Coupling parameter along the cols
  Double_t fSPDCouplRow;    // SPD Coupling parameter along the rows
  Int_t    fSPDHitStrobe[2]; // Hit read out strobe (left and right readout windows wrt the collision)
  Int_t    fSPDFoStrobe[2];  // FastOr read out strobe (left and right readout windows wrt the collision)
  Float_t  fSPDEccDiff;     // Eccentricity (i.e. asymmetry parameter) in the 
                            // Gaussian diffusion for SPD  
  Bool_t   fSPDLorentzDrift;     // Flag to decide whether to simulate the Lorentz Drift or not in SPD
  Float_t  fSPDLorentzHoleWeight;// Lorentz Angle is computed for SPD as average of Hole and Electron
                                 // this parameter gives the relative weights between the two
  Bool_t   fSPDAddNoisyFlag;     // Flag saying whether noisy pixels should be added to digits
  Bool_t   fSPDRemoveDeadFlag;   // Flag saying whether dead pixels should be removed from digits

  Int_t    fSDDElectronics;  // SDD Electronics Pascal (1) or OLA (2)
  Float_t  fSDDDiffCoeff;    // SDD Diffusion Coefficient (scaling the time)
  Float_t  fSDDDiffCoeff1;   // SDD Diffusion Coefficient (constant term)
  Float_t  fSDDJitterError;  // SDD jitter error
  Float_t  fSDDDynamicRange; // SDD Dynamic Range 
  Float_t  fSDDMaxAdc;       // SDD ADC saturation value
  Float_t  fSDDChargeLoss;   // Set Linear Coefficient for Charge Loss 
  Float_t  fSDDTrigDelay;    // SDD time-zero
  Float_t  fSDDMapPrec;      // SDD maps precision
  Float_t  fSDDkeVtoADC;     // SDD keV->ADC conv. factor
  Char_t   fSDDRawFormat;    // Index for SDD RawFormat
  
  Bool_t   fSSDLorentzDrift;     // Flag to decide whether to simulate the Lorentz Drift or not in SSD

  Double_t fSSDCouplingPR;  // SSD couplings
  Double_t fSSDCouplingPL;  // SSD couplings
  Double_t fSSDCouplingNR;  // SSD couplings
  Double_t fSSDCouplingNL;  // SSD couplings   
  Int_t    fSSDZSThreshold; // SSD threshold for the zero suppresion

  Float_t  fNsigmas;   // Number of sigmas over which charge disintegration
                       // is performed
  Int_t      fNcomps;  // Number of samplings along the gaussian
  TArrayF   *fGaus;    // Gaussian lookup table for signal generation

  Double_t fN;  // the impurity concentration of the material in #/cm^3  (NOT USED!)
  Float_t fT;   // The temperature of the Si in Degree K.

  ClassDef(AliITSSimuParam,7);
};
#endif
