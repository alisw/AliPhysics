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
class AliTRDparameter;

class AliTRDclusterizerV1 : public AliTRDclusterizer {

 public:

  AliTRDclusterizerV1();
  AliTRDclusterizerV1(const Text_t* name, const Text_t* title);
  AliTRDclusterizerV1(const AliTRDclusterizerV1 &c);
  virtual ~AliTRDclusterizerV1();
  AliTRDclusterizerV1 &operator=(const AliTRDclusterizerV1 &c);

  virtual void     Copy(TObject &c);
  virtual Bool_t   MakeClusters();
  virtual Bool_t   ReadDigits();
  virtual void     SetParameter(AliTRDparameter *par)      { fPar           = par; };

  AliTRDparameter *GetParameter()                    const { return fPar;          };

 protected:

  AliTRDdigitsManager *fDigitsManager;      //! TRD digits manager
  AliTRDparameter     *fPar;                //  TRD digitization parameter object

 private:

  virtual Float_t  Unfold(Float_t eps, Int_t plane, Float_t *padSignal);

  ClassDef(AliTRDclusterizerV1,4)           // TRD-Cluster finder, slow simulator

};

#endif
