#ifndef ALITRDCLUSTERIZERV1_H
#define ALITRDCLUSTERIZERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliTRDclusterizer.h"

///////////////////////////////////////////////////////
//  Finds and handles cluster (slow simulation)      //
///////////////////////////////////////////////////////

class AliTRDdigitsManager;

class AliTRDclusterizerV1 : public AliTRDclusterizer {

 public:

  AliTRDclusterizerV1();
  AliTRDclusterizerV1(const Text_t* name, const Text_t* title);
  AliTRDclusterizerV1(const AliTRDclusterizerV1 &c);
  virtual ~AliTRDclusterizerV1();
  AliTRDclusterizerV1 &operator=(const AliTRDclusterizerV1 &c);

  virtual void    Copy(TObject &c);
  virtual void    Init();
  virtual Bool_t  MakeClusters();
  virtual Bool_t  ReadDigits();

  virtual void    SetClusMaxThresh(Int_t thresh)          { fClusMaxThresh = thresh; };
  virtual void    SetClusSigThresh(Int_t thresh)          { fClusSigThresh = thresh; };

  virtual Int_t GetClusMaxThresh() const                  { return fClusMaxThresh; };
  virtual Int_t GetClusSigThresh() const                  { return fClusSigThresh; };

 protected:

  AliTRDdigitsManager *fDigitsManager; //! TRD digits manager

  Int_t              fClusMaxThresh; // Threshold value for cluster maximum
  Int_t              fClusSigThresh; // Threshold value for cluster signal

 private:

  virtual Float_t  Unfold(Float_t eps, Float_t *padSignal);
  virtual Float_t  PadResponse(Float_t x);

  ClassDef(AliTRDclusterizerV1,2)      // TRD-Cluster finder, slow simulator

};

#endif
