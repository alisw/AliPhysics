#ifndef ALITRDCALONLINEGAINTABLEROC_H
#define ALITRDCALONLINEGAINTABLEROC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include "AliTRDCalOnlineGainTableMCM.h"

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
// AliTRDCalOnlineGainTableROC is a class to allocate MCM Gain Tables 
// and to access all stored calibration values from the ROC level by indicating row and col
//
//////////////////////////////////////////////////////////////////////////////////////////////

class AliTRDCalOnlineGainTableROC: public TObject
{
public:

  AliTRDCalOnlineGainTableROC(); 
   AliTRDCalOnlineGainTableROC(const AliTRDCalOnlineGainTableROC& other);
   AliTRDCalOnlineGainTableROC& operator=(const AliTRDCalOnlineGainTableROC& other);
  ~AliTRDCalOnlineGainTableROC();

  Float_t GetGainCorrectionFactor(Int_t row, Int_t col); 
  Short_t GetAdcdac(Int_t row, Int_t col);
  Float_t GetMCMGain(Int_t row, Int_t col); 
  Short_t GetFGAN(Int_t row, Int_t col);
  Short_t GetFGFN(Int_t row, Int_t col);

  void AllocateGainTableMCM(Int_t rob, Int_t mcm);

  AliTRDCalOnlineGainTableMCM* GetGainTableMCM(Int_t index) const
  { 
    //returns the Gain Table of the given MCM
    return fMCMGainTables[index]; 
  }

  AliTRDCalOnlineGainTableMCM* GetGainTableMCM(Int_t rob, Int_t mcm) const
  { 
    //returns the Gain Table of the given MCM
    return GetGainTableMCM(16*rob+mcm); 
  }

protected:
  
  AliTRDCalOnlineGainTableMCM* fMCMGainTables[128]; // Array of gain tables for MCMs
  
  ClassDef(AliTRDCalOnlineGainTableROC,1);          // TRD online gain table of a ROC

};

#endif
