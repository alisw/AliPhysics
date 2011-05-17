#ifndef ALIQATHRESHOLDS_H
#define ALIQATHRESHOLDS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Container for the parameters that might be needed by the QA.              //
// Each detector has 1 of these objects with the list of parameters.         //
// The parameters are created online and are passed to the offline through   //
// the shuttle to be added to the OCDB in the GRP.                           //
// Created by Barthelemy.von.Haller@cern.ch  30/11/2010                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TObjArray.h"
#include "TParameter.h"

class AliQAThresholds: public TObject {

 public:

  AliQAThresholds(Int_t detId);
  virtual ~AliQAThresholds();

  Int_t GetDetectorId();
  void SetDetectorId(Int_t i);
  void AddThreshold(TParameter<long>* item);
  void AddThreshold(TParameter<int>* item);
  void AddThreshold(TParameter<double>* item);
  void AddThreshold(TParameter<float>* item);
  void AddThresholdAt(TParameter<long>* item, Int_t index);
  void AddThresholdAt(TParameter<int>* item, Int_t index);
  void AddThresholdAt(TParameter<double>* item, Int_t index);
  void AddThresholdAt(TParameter<float>* item, Int_t index);
  TObject* GetThreshold(Int_t i);
  Int_t GetSize();

 private:

  TObjArray fThresholds;
  Int_t fDetectorId;       // in the sense of the class AliDAQ

  ClassDef(AliQAThresholds, 2)
};

#endif
