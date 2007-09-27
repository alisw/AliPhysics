#ifndef ALIPMDCALIBPEDESTAL_H
#define ALIPMDCALIBPEDESTAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TObject.h"

class TH1F;
class AliRawReader;

class AliPMDCalibPedestal : public TObject {

public:
  AliPMDCalibPedestal();
  AliPMDCalibPedestal(const AliPMDCalibPedestal &ped);
  AliPMDCalibPedestal& operator = (const  AliPMDCalibPedestal &source);
  virtual ~AliPMDCalibPedestal();

  Bool_t ProcessEvent(AliRawReader  *rawReader);
  void   Analyse(TTree *pedtree);

private:

  TH1F *fPedHisto[2][24][48][96];


  ClassDef(AliPMDCalibPedestal,2)
};



#endif

