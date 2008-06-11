#ifndef ALIPMDCALIBPEDESTAL_H
#define ALIPMDCALIBPEDESTAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TObject.h"

class TH1F;
class AliRawReader;
class TTree;

class AliPMDCalibPedestal : public TObject {

public:
  AliPMDCalibPedestal();
  AliPMDCalibPedestal(const AliPMDCalibPedestal &ped);
  AliPMDCalibPedestal& operator = (const  AliPMDCalibPedestal &source);
  virtual ~AliPMDCalibPedestal();

  Bool_t ProcessEvent(AliRawReader  *rawReader, TObjArray *pmdddlcont);
  void   Analyse(TTree *pedtree);

private:

  enum
      {
	  kDet    = 2,   // Number of Planes
	  kMaxSMN = 24,  // Number of Modules
	  kMaxRow = 48,  // Number of Rows
	  kMaxCol = 96   // Number of Columns
      };

  Float_t fPedVal[kDet][kMaxSMN][kMaxRow][kMaxCol];
  Float_t fPedValSq[kDet][kMaxSMN][kMaxRow][kMaxCol];
  Float_t fPedCount[kDet][kMaxSMN][kMaxRow][kMaxCol];
  Int_t   fPedChain[6][51][25][64];
  Int_t   fRunNumber;
  Int_t   fEventNumber;

  ClassDef(AliPMDCalibPedestal,5)
};



#endif

