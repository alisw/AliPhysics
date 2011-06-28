#ifndef ALITRDCALONLINEGAINTABLEMCM_H
#define ALITRDCALONLINEGAINTABLEMCM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

//////////////////////////////////////////////////////////////////////////////////////////////
//
// Data structure to store gaintables of the online calibration in the OCDB
// consisting of three classes:
// AliTRDCalOnlineGainTable 
// AliTRDCalOnlineGainTableROC 
// AliTRDCalOnlineGainTableMCM
//
// AliTRDCalOnlineGainTable is the main class from which all stored data can be accessed.
// The two sub-classes AliTRDCalOnlineGainTableROC and AliTRDCalOnlineGainTableMCM
// contain the gaintables on ROC level and on the MCM level respectively.
//
// The online calibration is used to compensate gain deviations on the pad level.
// For the offline reconstruction the online calibration has to be undone. 
// The corresponding gain correction factor that was used by the online gain filter can be accessed 
// via the functions AliTRDCalOnlineGainTable::GetGainCorrectionFactor(Int_t det, Int_t row, Int_t col) 
// and AliTRDCalOnlineGainTable::GetGainCorrectionFactor(Int_t sector, Int_t stack, Int_t layer, Int_t row, Int_t col).
//
// With the class AliTRDCalOnlineGainTablesMCM all values used for the 
// online calibration can be set and accessed on the MCM/channel level
//
//////////////////////////////////////////////////////////////////////////////////////////////

class AliTRDCalOnlineGainTableMCM: public TObject
{
public:

  AliTRDCalOnlineGainTableMCM();
  ~AliTRDCalOnlineGainTableMCM();

  Float_t GetGainCorrectionFactor(Int_t channel);
  Float_t GetMCMGain();
  Short_t GetAdcdac();
  Short_t GetFGAN(Int_t channel);
  Short_t GetFGFN(Int_t channel);
  
  void SetAdcdac(Short_t x) {fAdcdac = x;} // Sets fAdcdac to the given value
  void SetMCMGain(Float_t x) {fMCMGain = x;} // Sets fMCMGain to the given value

  void SetFGFN(Short_t ch, Short_t x) {fFGFN[ch] = x;} // Sets fFGFN to the given value
  void SetFGAN(Short_t ch, Short_t x) {fFGAN[ch] = x;} // Sets fFGAN to the given value

protected:

  Short_t fAdcdac;   // Reference voltage of the ADCs  U_Ref =  (1.05V + (fAdcdac/31)*0.4V

  Short_t fFGFN[21]; // Gain Correction Filter Factor
  Short_t fFGAN[21]; // Gain Correction Filter Additive

  Float_t fMCMGain;  // Gain Factor which would lead to a Correction Factor of 1.0 within the MCM


  ClassDef(AliTRDCalOnlineGainTableMCM,1); // TRD online gain table of a MCM

};

#endif
