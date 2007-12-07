#ifndef ALIPMDCALIBGAIN_H
#define ALIPMDCALIBGAIN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TObject.h"

class TH1F;  
class AliRawReader;

class AliPMDCalibGain : public TObject
{
 public:
  AliPMDCalibGain() ;              // ctor
  AliPMDCalibGain(const AliPMDCalibGain &pmdcalibgain);  // copy constructor
  AliPMDCalibGain &operator=(const AliPMDCalibGain &pmdcalibgain); // assignment op

  virtual ~AliPMDCalibGain() ;     // dtor

  Bool_t ProcessEvent(AliRawReader *rawReader);  //Looks for iso cells

  void Analyse(TTree *gaintree);
  
 private:

  enum
      {
	  kDet    = 2,   // Number of Planes
	  kMaxSMN = 24,  // Number of Modules
	  kMaxRow = 48,  // Number of Rows
	  kMaxCol = 96   // Number of Columns
      };

  Float_t fSMIso[kDet][kMaxSMN];
  Float_t fSMCount[kDet][kMaxSMN]; // counter
  Float_t fCellIso[kDet][kMaxSMN][kMaxRow][kMaxCol]; // adc of iso cells
  Float_t fCellCount[kDet][kMaxSMN][kMaxRow][kMaxCol]; // counter

ClassDef(AliPMDCalibGain,2)        // description 
};
#endif // ALIPMDCALIBGAIN_H
