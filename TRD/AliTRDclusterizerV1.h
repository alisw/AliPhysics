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
          void    UseLUT()                                { fUseLUT        = kTRUE;  };
          void    UseCOG()                                { fUseLUT        = kFALSE; };

          void    SetClusMaxThresh(Int_t thresh)          { fClusMaxThresh = thresh; };
          void    SetClusSigThresh(Int_t thresh)          { fClusSigThresh = thresh; };

          Int_t   GetClusMaxThresh() const                { return fClusMaxThresh; };
          Int_t   GetClusSigThresh() const                { return fClusSigThresh; };

 protected:
 
  enum { 
    kNlut = 128                        //  Dimension of the lookup table
  };                    

  AliTRDdigitsManager *fDigitsManager; //! TRD digits manager

  Int_t                fClusMaxThresh; //  Threshold value for cluster maximum
  Int_t                fClusSigThresh; //  Threshold value for cluster signal
  Bool_t               fUseLUT;        //  Switch for the lookup table method
  Float_t              fLUT[kNlut];    //  The lookup table

 private:

  virtual Float_t  Unfold(Float_t eps, Float_t *padSignal);
  virtual Float_t  PadResponse(Float_t x);

  ClassDef(AliTRDclusterizerV1,2)      // TRD-Cluster finder, slow simulator

};

#endif
