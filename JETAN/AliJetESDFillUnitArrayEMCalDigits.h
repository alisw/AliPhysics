#ifndef ALIJETESDFILLUNITARRAYEMCALDIGITS_H
#define ALIJETESDFILLUNITARRAYEMCALDIGITS_H
 
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
class AliJetESDReader;
class AliJetESDReaderHeader;
//class AliEMCALCalibData ;

class AliJetESDFillUnitArrayEMCalDigits : public AliJetFillUnitArray
{
 public: 
  AliJetESDFillUnitArrayEMCalDigits();
  AliJetESDFillUnitArrayEMCalDigits(AliESDEvent *fESD);
  virtual ~AliJetESDFillUnitArrayEMCalDigits();
  
  // Setter
  void SetApplyElectronCorrection(Int_t flag = 1)     {fApplyElectronCorrection = flag;}
  void SetApplyFractionHadronicCorrection(Bool_t val) {fApplyFractionHadronicCorrection = val;}
  void SetFractionHadronicCorrection(Double_t val)    {fFractionHadronicCorrection = val;}
  void SetESD(AliESDEvent* const esd)                 {fESD = esd;}
  void SetInitMult(Int_t mult)                        {fNDigitEmcal = mult;}
  void SetInitMultCut(Int_t multcut)                  {fNDigitEmcalCut = multcut;}

  // Getter
  Int_t         GetMult()    const {return fNDigitEmcal;}
  Int_t         GetMultCut() const {return fNDigitEmcalCut;}

  // Other
  void Exec(Option_t* const option);
  // For calibration 
  //  virtual Float_t Calibrate(Int_t amp, Int_t cellId) ;  // Tranforms Amp to energy

 protected:
  AliESDEvent  *fESD;                 // ESD
  Int_t        fNIn;                  // Number of Array filled in UnitArray
  Int_t        fOpt;                  // Detector to be used for jet reconstruction
  Int_t        fCluster;              // Use all cells or cells in clusters for jet finding 
  Int_t        fDebug;                // Debug option
  Int_t        fNCEMCAL;              // Number of clusters in EMCAL
  Int_t        fNCPHOS;               // Number of clusters in PHOS
  Int_t        fNCCalo;               // Number of cluster in EMCAL + PHOS calorimeters

  Bool_t       fApplyElectronCorrection;          // Electron correction flag
  Bool_t       fApplyFractionHadronicCorrection;  // Fraction hadronic correction flag
  Bool_t       fFractionHadronicCorrection;       // Fraction hadronic correction 

  AliESDCaloCluster *fClus;           //! 
  Int_t              fNDigitEmcal;    //!
  Int_t              fNDigitEmcalCut; //!
  // Calibration parameters... to be replaced by database
/*   AliEMCALCalibData *fCalibData;      //! Calibration database if aval */
/*   Float_t            fADCchannelECA;  // width of one ADC channel for EC section (GeV) */
/*   Float_t            fADCpedestalECA; // pedestal of ADC for EC section (GeV) */


 private:
  AliJetESDFillUnitArrayEMCalDigits(const AliJetESDFillUnitArrayEMCalDigits &det);
  AliJetESDFillUnitArrayEMCalDigits &operator=(const AliJetESDFillUnitArrayEMCalDigits &det);

//  void    GetCalibrationParameters(void) ;
  
  ClassDef(AliJetESDFillUnitArrayEMCalDigits,1) // Fill Unit Array with tpc and/or emcal information
};

#endif
