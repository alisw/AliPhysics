#ifndef ALIPXCONEJETHEADER_H
#define ALIPXCONEJETHEADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
 
//---------------------------------------------------------------------
// Jet header class for Pxcone algorithm
// Stores the parameters of the Pxcone jet algorithm
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------
 
#include "AliJetHeader.h"
 
class AliPxconeJetHeader : public AliJetHeader
{
 public:
 
  AliPxconeJetHeader();
  virtual ~AliPxconeJetHeader() { }

  // Getters

  Int_t    GetMode()    const {return fMode;}
  Double_t GetRadius()  const {return fRadius;}
  Double_t GetMinPt()   const {return fMinPt;}
  Double_t GetOverlap() const {return fOverlap;}
   
  // Setters
  void SetMode(Int_t m=2) {fMode=m;}
  void SetRadius(Double_t r=0.3) {fRadius=r;}
  void SetMinPt(Double_t p=10) {fMinPt=p;}
  void SetOverlap(Double_t o=0.75) {fOverlap=o;}

  // others
  void PrintParameters() const;
   
protected:
  Int_t fMode;           // ee or pp mode
  Double_t fRadius;      // jet radius
  Double_t fMinPt;       // min pt of jets  
  Double_t fOverlap;     // fraction of overlap energy

  ClassDef(AliPxconeJetHeader,1)
};
 
#endif
