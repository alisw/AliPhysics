#ifndef ALIITSSIMUPARAM_H
#define ALIITSSIMUPARAM_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id:$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to store the parameters used in the simulation of       //
// SPD, SDD and SSD detectors                                    //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

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



  void    SetSPDBiasVoltage(Double_t bias=18.182) {fSPDBiasVoltage=bias;}
  Double_t  GetSPDBiasVoltage() const {return fSPDBiasVoltage;}

  void   SetSPDThresholds(Double_t thresh, Double_t sigma)
	{fSPDThresh=thresh; fSPDSigma=sigma;}
  void   Thresholds(Double_t &thresh, Double_t &sigma) const
	{thresh=fSPDThresh; sigma=fSPDSigma;}

  void SetSPDCouplingOption(const char *opt) {fSPDCouplOpt=opt;}
  void GetSPDCouplingOption(char *opt) const {strcpy(opt,fSPDCouplOpt.Data());}

  void SetSPDCouplingParam(Double_t col, Double_t row)
        {fSPDCouplCol = col; fSPDCouplRow = row;}
  void GetSPDCouplingParam(Double_t &col, Double_t &row) const
        {col = fSPDCouplCol; row = fSPDCouplRow;}

  void   SetSPDSigmaDiffusionAsymmetry(Double_t ecc) {fSPDEccDiff=ecc;}   
  void   GetSPDSigmaDiffusionAsymmetry(Double_t &ecc) const {ecc=fSPDEccDiff;}

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


  void SetSSDADCpereV(Double_t a=120./24888.9){fSSDADCpereV = a;}
  Double_t GetSSDDEvToADC(Double_t eV) const {return eV*fSSDADCpereV;}
  Int_t GetSSDIEvToADC(Double_t eV) const { 
      return ((Int_t) GetSSDDEvToADC(eV)); }

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

  void PrintParameters() const; 

 protected:

  static const Float_t fgkSPDBiasVoltageDefault;//default for fSPDBiasVoltage
  static const Double_t fgkSPDThreshDefault; //default for fThresh
  static const Double_t fgkSPDSigmaDefault; //default for fSigma
  static const TString fgkSPDCouplingOptDefault;  // type of pixel Coupling (old or new)
  static const Double_t fgkSPDCouplColDefault; //default for fSPDCouplCol
  static const Double_t fgkSPDCouplRowDefault; //default for fSPDCouplRow
  static const Float_t fgkSPDEccDiffDefault;//default for fSPDEccDiff
  static const Float_t fgkSDDDiffCoeffDefault; // default for fSDDDiffCoeff
  static const Float_t fgkSDDDiffCoeff1Default; // default for fSDDDiffCoeff1 
  static const Float_t fgkSDDJitterErrorDefault; // default for fSDDJitterError
  static const Float_t fgkSDDDynamicRangeDefault; // default for fSDDDynamicRange
  static const Int_t fgkSDDMaxAdcDefault; // default for fSDDMaxAdc
  static const Float_t fgkSDDChargeLossDefault; // default for fSDDChargeLoss

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

  
  Double_t fSPDBiasVoltage; // Bias Voltage for the SPD
  Double_t fSPDThresh;      // SPD Threshold value
  Double_t fSPDSigma;       // SPD Noise + threshold fluctuations value  
  TString  fSPDCouplOpt;    // SPD Coupling Option
  Double_t fSPDCouplCol;    // SPD Coupling parameter along the cols
  Double_t fSPDCouplRow;    // SPD Coupling parameter along the rows
  Float_t  fSPDEccDiff;     // Eccentricity (i.e. asymmetry parameter) in the 
                            // Gaussian diffusion for SPD

  Int_t    fSDDElectronics;  // SDD Electronics Pascal (1) or OLA (2)
  Float_t  fSDDDiffCoeff;    // SDD Diffusion Coefficient (scaling the time)
  Float_t  fSDDDiffCoeff1;   // SDD Diffusion Coefficient (constant term)
  Float_t  fSDDJitterError;  // SDD jitter error
  Float_t  fSDDDynamicRange; // SDD Dynamic Range 
  Float_t  fSDDMaxAdc;       // SDD ADC saturation value
  Float_t  fSDDChargeLoss;   // Set Linear Coefficient for Charge Loss 
  
  Double_t fSSDADCpereV;    // Constant to convert eV to ADC for SSD.
  Double_t fSSDCouplingPR;  // SSD couplings
  Double_t fSSDCouplingPL;  // SSD couplings
  Double_t fSSDCouplingNR;  // SSD couplings
  Double_t fSSDCouplingNL;  // SSD couplings   
  Int_t    fSSDZSThreshold; // SSD threshold for the zero suppresion

  Float_t  fNsigmas;   // Number of sigmas over which charge disintegration
                       // is performed
  Int_t      fNcomps;  // Number of samplings along the gaussian
  TArrayF   *fGaus;    // Gaussian lookup table for signal generation

  ClassDef(AliITSSimuParam,1);
};
#endif
