#ifndef MUONPROTO_H
#define MUONPROTO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMUON.h"

class AliMUONproto : public AliMUON {
  public:
    AliMUONproto();
    AliMUONproto(const char *name, const char *title);
    virtual ~AliMUONproto() {}
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual void   Init();
    virtual Int_t  IsVersion() const {return 2;}
    virtual void   StepManager();
//    void GetRawDigits(Int_t iev, Int_t* ptr, Int_t len);
    void SetPadSize(Int_t id, Int_t isec, Float_t p1, Float_t p2);
//    void SetThreshold();
    void SetNsigma(Int_t nsig) {fNsigma = nsig;};
    void BuildGeometry();
    void SetChargeSlope(Int_t id, Float_t p1);
    void SetChargeSpread(Int_t id, Float_t p1, Float_t p2);
    void SetMaxAdc(Int_t id, Float_t p1);
    void FindClusters(Int_t, Int_t);

  private:
    ClassDef(AliMUONproto,1) // MUON Detector class for test beam prototypes
  protected:
    Int_t fNsigma;
    Float_t fThreshold[100];
};
#endif
