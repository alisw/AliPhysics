#ifndef ALIJETAODFILLUNITARRAYEMCALDIGITS_H
#define ALIJETAODFILLUNITARRAYEMCALDIGITS_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Jet Fill Unit Array 
// Called by ESD Reader for jet analysis
// Author: Magali Estienne (magali.estienne@subatech.in2p3.fr)
//---------------------------------------------------------------------

#include "AliJetFillUnitArray.h"

class AliESDCaloCluster;
class AliJetReader;
class AliJetAODReader;

class AliJetAODFillUnitArrayEMCalDigits : public AliJetFillUnitArray
{
 public: 
  AliJetAODFillUnitArrayEMCalDigits();
  AliJetAODFillUnitArrayEMCalDigits(AliAODEvent *fAOD);
  virtual ~AliJetAODFillUnitArrayEMCalDigits();
  
  // Setter
  void SetApplyElectronCorrection(Int_t flag = 1)     {fApplyElectronCorrection = flag;}
  void SetApplyFractionHadronicCorrection(Bool_t val) {fApplyFractionHadronicCorrection = val;}
  void SetFractionHadronicCorrection(Double_t val)    {fFractionHadronicCorrection = val;}
  void SetAOD(AliAODEvent* const aod)                 {fAOD = aod;}
  void SetInitMult(Int_t mult)                        {fNDigitEmcal = mult;}
  void SetInitMultCut(Int_t multcut)                  {fNDigitEmcalCut = multcut;}

  // Getter
  Int_t         GetMult()      const {return fNDigitEmcal;}
  Int_t         GetMultCut()   const {return fNDigitEmcalCut;}

  // For calibration 
  //  virtual Float_t Calibrate(Int_t amp, Int_t cellId) ;  // Tranforms Amp to energy

  // Other
  void Exec(Option_t* const option);

 protected:
  AliAODEvent  *fAOD;                 // ESD
  Int_t        fNIn;                  // Number of Array filled in UnitArray
  Int_t        fCluster;              // Use all cells or cells in clusters for jet finding 
  Int_t        fNCEMCAL;              // Number of clusters in EMCAL
  Int_t        fNCPHOS;               // Number of clusters in PHOS
  Int_t        fNCCalo;               // Number of cluster in EMCAL + PHOS calorimeters

  Bool_t       fApplyElectronCorrection;          // Electron correction flag
  Bool_t       fApplyFractionHadronicCorrection;  // Fraction hadronic correction flag
  Double_t     fFractionHadronicCorrection;       // Fraction hadronic correction 

  //Track-matching (mth)
// Temporarily commented - will be updated removing the AliAODpid dependence
//  void GetTracksPointingToCell(TArrayS *arr, Double_t eta, Double_t phi, Double_t res);

  // geometry info
  AliAODCaloCluster *fClus;           //! 
  Int_t fNDigitEmcal;                 //!
  Int_t fNDigitEmcalCut;              //!
  //Calibration parameters... to be replaced by database
//  AliEMCALCalibData *fCalibData;     //! Calibration database if aval
//  Float_t            fADCchannelECA; // width of one ADC channel for EC section (GeV)
//  Float_t            fADCpedestalECA;// pedestal of ADC for EC section (GeV)


 private:
  AliJetAODFillUnitArrayEMCalDigits(const AliJetAODFillUnitArrayEMCalDigits &det);
  AliJetAODFillUnitArrayEMCalDigits &operator=(const AliJetAODFillUnitArrayEMCalDigits &det);

//  void    GetCalibrationParameters(void) ;
  
  ClassDef(AliJetAODFillUnitArrayEMCalDigits,1) // Fill Unit Array with tpc and/or emcal information
};

#endif
