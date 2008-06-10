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

  Bool_t ProcessEvent(AliRawReader  *rawReader);
  void   Analyse(TTree *pedtree);
  void   ConvertDDL(Int_t det, Int_t smn, Int_t &ddlno);

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
  UInt_t  fPedChain[kDet][kMaxSMN][kMaxRow][kMaxCol];
  

  ClassDef(AliPMDCalibPedestal,4)
};



#endif

