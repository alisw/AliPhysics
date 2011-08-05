#ifndef ALITRDCALONLINEGAINTABLE_H
#define ALITRDCALONLINEGAINTABLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include "AliTRDCalOnlineGainTableROC.h"

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
//////////////////////////////////////////////////////////////////////////////////////////////

class AliTRDCalOnlineGainTable: public TObject
{
public:

  AliTRDCalOnlineGainTable();
  AliTRDCalOnlineGainTable(const AliTRDCalOnlineGainTable& other);
  AliTRDCalOnlineGainTable& operator=(const AliTRDCalOnlineGainTable& other);
  ~AliTRDCalOnlineGainTable();

  Float_t GetGainCorrectionFactor(Int_t det, Int_t row, Int_t col) const;
  Float_t GetGainCorrectionFactor(Int_t sector, Int_t stack, Int_t layer, Int_t row, Int_t col) const;

  Short_t GetAdcdac(Int_t det, Int_t row, Int_t col);
  Short_t GetAdcdac(Int_t sector, Int_t stack, Int_t layer, Int_t row, Int_t col);

  Float_t GetMCMGain(Int_t det, Int_t row, Int_t col);
  Float_t GetMCMGain(Int_t sector, Int_t stack, Int_t layer, Int_t row, Int_t col);

  Short_t GetFGAN(Int_t det, Int_t row, Int_t col);
  Short_t GetFGAN(Int_t sector, Int_t stack, Int_t layer, Int_t row, Int_t col);

  Short_t GetFGFN(Int_t det, Int_t row, Int_t col);
  Short_t GetFGFN(Int_t sector, Int_t stack, Int_t layer, Int_t row, Int_t col);

  void AllocateGainTableROC(Int_t det);

  AliTRDCalOnlineGainTableROC *GetGainTableROC(Int_t det) const
  {
    // returns the Gain Table for the given detector
    return fROCGainTables[det];
  }

  AliTRDCalOnlineGainTableROC *GetGainTableROC(Int_t sector, Int_t stack, Int_t layer) const
  {
    // returns the Gain Table for the given detector
    return GetGainTableROC(sector*30 + stack*6 + layer);
  } 

  static const Float_t UnDef;

protected:

  AliTRDCalOnlineGainTableROC* fROCGainTables[540]; // Array of gain tables for all ROCs
  
  ClassDef(AliTRDCalOnlineGainTable,2);             // TRD online gain tables

};

#endif
