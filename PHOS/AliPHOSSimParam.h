#ifndef ALIPHOSSIMPARAM_H
#define ALIPHOSSIMPARAM_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id: AliPHOSSimParam.h 23530 2008-01-25 06:46:13Z prsnko $ */
                                              
// Base class for the PHOS simulation parameters.
// Do not use in the reconstruction; use derivative classes instead.

#include "TNamed.h"
#include "TMath.h"

class AliPHOSSimParam : public TNamed {

public:

  AliPHOSSimParam();
  AliPHOSSimParam(const AliPHOSSimParam& recoParam);
  AliPHOSSimParam& operator = (const AliPHOSSimParam& recoParam);
  virtual ~AliPHOSSimParam() {}

  static AliPHOSSimParam * GetInstance() ;

  //Parameters used in conversion of deposited energy to APD response [see AliPHOSv1 for details]
  void SetMeanLightYield(Float_t ly=47000.){fLightYieldMean=ly; //recalculate now dependencies
                                            fLightFactor = fLightYieldMean * fIntrinsicAPDEfficiency ; 
                                            fAPDFactor   = (13.418/fLightYieldMean/100.) * 300. ; }
  void SetAPDEfficiency(Float_t eff=0.02655){fIntrinsicAPDEfficiency=eff;
                                             fLightFactor = fLightYieldMean * fIntrinsicAPDEfficiency ;}
  Float_t GetLightFactor(void) const { return fLightFactor ;}
  Float_t GetAPDFactor(void) const { return fAPDFactor ;}

  //Parameters use in EMC noise calculation [see AliPHOSDigitizer for details]
  Float_t GetAPDNoise() const { return fAPDNoise;  }          //RMS of APD noise
  void SetAPDNoise(Float_t noise=0.012){fAPDNoise = noise;  }

  //Parameters to apply non-lineary on cell level
  Bool_t IsCellNonlinearityOn() const {return fCellNonLineaityOn;}
  void SetCellNonLinearity(Bool_t on=kTRUE){fCellNonLineaityOn=on;} //default: on=kFALSE
  Double_t GetCellNonLineairyA(void) const {return fCellNonLineaityA; }
  Double_t GetCellNonLineairyB(void) const {return fCellNonLineaityB; }
  void SetCellNonLineairyA(Double_t a=0.30) {fCellNonLineaityA = a; }
  void SetCellNonLineairyB(Double_t b=0.109){fCellNonLineaityB = b; }



  Float_t GetEmcDigitsThreshold() const { return fEMCDigitThreshold ; }  //Minimal energy to keep digit
  void SetEMCDigitsThreshold(Float_t thresh=0.01){fEMCDigitThreshold=thresh;} 

  //Parameters for energy digitization [see AliPHOSDigitizer for details]
  void SetADCchannelW(Float_t width=0.005){fEMCADCchannel=width ;} //EMC channel width
  Float_t GetADCchannelW(void) const {return fEMCADCchannel; }

  void SetEDigitizationOn(Bool_t on=kTRUE){fDigitizeE=on ;}   //Use digitization in simulation or left it
  Bool_t IsEDigitizationOn(void)const {return fDigitizeE ;}   //for Digits2Raw->Digits procedure

  //Parameters for EMC TOF resolution [see AliPHOSDigitizer::TimeResolution()]
  Float_t GetTOFa()const{return fTOFa ;}  //constant term
  Float_t GetTOFb()const{return fTOFb ;}  //stohastic term
  void SetTOFparameters(Float_t a=0.5e-9, Float_t b=1.5e-9){fTOFa=a; fTOFb=b; }

  //Parameters for CPV noise and digitization  [see AliPHOSDigitizer for details] 
  Float_t GetCPVNoise() const {return fCPVNoise ;}           //RMS of CPV noise in
  void SetCPVNoise(Float_t noise=0.01){ fCPVNoise = noise ;} //CPV popugais

  Float_t GetCpvDigitsThreshold() const {return fCPVDigitThreshold ;}           //Minimal energy to keep digit in
  void SetCpvDigitsThreshold(Float_t thresh=0.09){fCPVDigitThreshold = thresh;} //CPV popugais

  Float_t GetADCpedestalCpv() const {return fADCpedestalCpv ;}      //CPV pedestal value 
  void SetADCpedestalCpv(Float_t ped=0.012){ fADCpedestalCpv=ped ;} //in CPV popugais

  Float_t GetADCchanelCpv() const {return fADCchanelCpv;}     //Price of one ADC channel 
  void SetADCchanelCpv(Float_t w=0.0012) {fADCchanelCpv=w; }  //for CPV

  Int_t GetNADCcpv() const {return fNADCcpv ;}                         //Max number of channels
  void SettNADCcpv(Int_t n=12) { fNADCcpv=(Int_t)TMath::Power(2,n) ; } //in CPV  ADC

  //Mark streams for mixing as streams contaning Digits (true) or SDigits (false)
  //Streamt numbering the same as in StreamManager
  void    SetStreamDigits(Int_t i){if(i<10)fDStream[i]=kTRUE;}
  Bool_t  IsStreamDigits(Int_t i){return fDStream[i]; }

  //Parameters for RAW embedding
  void SetEMCSubtractPedestals(Bool_t subtract) { fEMCSubtractPedestals = subtract;}
  Bool_t  EMCSubtractPedestals()      const { return fEMCSubtractPedestals;    }

  void SetGlobalAltroOffset(Int_t offset)  { fGlobalAltroOffset = offset ; }
  Int_t   GetGlobalAltroOffset()      const { return fGlobalAltroOffset ;  }

  void SetGlobalAltroThreshold(Int_t ZSth) { fGlobalAltroThreshold = ZSth; }
  Int_t   GetGlobalAltroThreshold()   const { return fGlobalAltroThreshold;}

  void SetSampleQualityCut(Float_t qcut) { fEMCSampleQualityCut=qcut; }
  Float_t GetEMCSampleQualityCut()    const { return fEMCSampleQualityCut; }

private:

  AliPHOSSimParam(Int_t i); //True constructor which should be called by GetInstance()

private:

  //Parameters used in conversion of deposited energy to APD response (AliPHOSv1)
  Float_t  fLightYieldMean ;        //Average number of photoelectrons per GeV
  Float_t  fIntrinsicAPDEfficiency; //APD efficiency including geometric coverage
  Float_t  fLightFactor ;           //Average number of photons collected by APD per GeV deposited energy
  Float_t  fAPDFactor ;             //factor relating light yield and APD response 
 
  //Parameters used in electronic noise calculation and thresholds for EMC (AliPHOSDigitizer)
  Float_t fAPDNoise;            //RMS of APD noise
  Float_t fEMCDigitThreshold ;  //minimal energy to keep digit 
  Float_t fEMCADCchannel ;      //width of ADC channel in GeV
  Float_t fTOFa  ;              //constant term of TOF resolution 
  Float_t fTOFb  ;              //stohastic term of TOF resolution 
  Float_t fCellNonLineaityA ;   //Amp of cel non-linearity
  Float_t fCellNonLineaityB ;   //Energy scale of cel non-linearity

  //Parameters used for RAW embedding
  Bool_t  fEMCSubtractPedestals;   // true if pedestal should be subtracted (in non-ZS)
  Int_t   fGlobalAltroOffset ;     // Offset used in ALTRO chips in SZ runs
  Int_t   fGlobalAltroThreshold ;  // Threshold used in ALTRO chips in SZ runs
  Float_t fEMCSampleQualityCut;    // Cut on pulse shape fit quality

  //CPV parameters
  Float_t fADCpedestalCpv ;    //Pedestal value
  Float_t fADCchanelCpv ;      //ADC channel width
  Float_t fCPVNoise ;          //RMS of CPV noise 
  Float_t fCPVDigitThreshold ; //Minimal energy to keep digit 
  Int_t   fNADCcpv ;           //Max number of channels in CPV ADC
 
  Bool_t fDStream[10] ;   //Mark mixing stream contains digits or SDigits
  Bool_t fDigitizeE ;     //Use energy digitization in simulation or left to Digits2Raw()
  Bool_t fCellNonLineaityOn ;  //Model scintillator non-linearity in AliPHOSDigitizer
  
  static AliPHOSSimParam * fgSimParam ; // pointer to the unique instance of the class

  ClassDef(AliPHOSSimParam,2)
};

#endif
