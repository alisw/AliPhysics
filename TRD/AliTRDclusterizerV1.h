#ifndef ALITRDCLUSTERIZERV1_H
#define ALITRDCLUSTERIZERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliTRDclusterizer.h"

///////////////////////////////////////////////////////
//  Finds and handles cluster (slow simulation)      //
///////////////////////////////////////////////////////

class AliTRDdataArrayI;
class AliTRDdigitsManager;
class AliTRDparameter;
class AliRawReader;

class AliTRDclusterizerV1 : public AliTRDclusterizer {

 public:

  AliTRDclusterizerV1();
  AliTRDclusterizerV1(const Text_t* name, const Text_t* title);
  AliTRDclusterizerV1(const AliTRDclusterizerV1 &c);
  virtual ~AliTRDclusterizerV1();
  AliTRDclusterizerV1 &operator=(const AliTRDclusterizerV1 &c);

  virtual void     Copy(TObject &c) const;
  virtual Bool_t   MakeClusters();
  virtual Bool_t   ReadDigits();
  virtual Bool_t   ReadDigits(AliRawReader* rawReader);

 protected:

  AliTRDdigitsManager *fDigitsManager;      //! TRD digits manager

 private:
  void DeConvExp(Double_t *source, Double_t *target, Int_t nTimeTotal, Int_t nexp);
  void Transform(AliTRDdataArrayI* digitsIn, AliTRDdataArrayI* digitsOut,
		 Int_t idet, Int_t nRowMax, Int_t nColMax, Int_t nTimeTotal, Float_t ADCthreshold);
  virtual Double_t Unfold(Double_t eps, Int_t plane, Double_t *padSignal);
  Double_t GetCOG(Double_t signal[5]);      // get COG position

  ClassDef(AliTRDclusterizerV1,5)           // TRD-Cluster finder, slow simulator

};

#endif
