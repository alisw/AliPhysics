/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef AliFlowLYZHist2_H
#define AliFlowLYZHist2_H

#include "TComplex.h"
#include "TProfile.h"   //no forward declaration possible because of inline functions

class TProfile2D;


// Description: Class to organise histograms for Flow
//              by the LeeYangZeros method in the second run.
//              Also contains methods to get values from the histograms
//              which are called in AliFlowLeeYandZerosMaker::Finish().


class AliFlowLYZHist2 {

 public:

  AliFlowLYZHist2(Int_t i1, Int_t i2);          //constructor
  AliFlowLYZHist2(const AliFlowLYZHist2&);      //copy constructor (dummy)
  virtual  ~AliFlowLYZHist2();                  //destructor

  AliFlowLYZHist2& operator=(const AliFlowLYZHist2&);  //assignment operator (dummy)

  void Fill(Float_t f1,Float_t f2, TComplex C); //fill the histograms
  Int_t GetNbinsX()                {Int_t fMaxEtaBins = fHistProReNumer->GetNbinsX();  return fMaxEtaBins;}     
  Int_t GetNbinsXPt()              {Int_t fMaxPtBins = fHistProReNumerPt->GetNbinsX(); return fMaxPtBins;}
  Float_t GetBinCenter(Int_t i)    {Float_t fEta = fHistProReNumer->GetXaxis()->GetBinCenter(i);  return fEta;}
  Float_t GetBinCenterPt(Int_t i)  {Float_t fPt = fHistProReNumerPt->GetXaxis()->GetBinCenter(i); return fPt;}
  TComplex GetfNumer(Int_t i);     //get numerator for diff. flow (eta)
  TComplex GetfNumerPt(Int_t i);   //get numerator for diff. flow (pt)
  Int_t GetNprime(Int_t i);        //get number of entries in bin (eta)
  Int_t GetNprimePt(Int_t i);      //get number of entries in bin (pt) 


 private:
  
  TProfile* fHistProReNumer;                     //!
  TProfile* fHistProImNumer;                     //!
  TProfile* fHistProReNumerPt;                   //!
  TProfile* fHistProImNumerPt;                   //!
  TProfile2D* fHistProReNumer2D;                 //!
  TProfile2D* fHistProImNumer2D;                 //!
  




  ClassDef(AliFlowLYZHist2,0)                    // macro for rootcint
    };
 
     
#endif

