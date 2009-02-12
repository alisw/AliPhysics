/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliFlowLYZHist2.h 29627 2008-10-31 10:30:41Z snelling $ */

#ifndef AliFlowLYZHist2_H
#define AliFlowLYZHist2_H

#include "TProfile.h"   //no forward declaration possible because of inline functions

class TComplex;
class TProfile2D;
class TCollection;
class TList;
class TBrowser;

// Description: Class to organise histograms for Flow
//              by the LeeYangZeros method in the second run.
//              Also contains methods to get values from the histograms
//              which are called in AliFlowLeeYandZerosMaker::Finish().


class AliFlowLYZHist2: public TNamed  {

 public:

  AliFlowLYZHist2();                                //default constructor
  AliFlowLYZHist2(Int_t theta, const char *name = "AliFlowLYZHist2" ,const char *title = "AliFlowLYZHist2");                     //constructor
  virtual  ~AliFlowLYZHist2();                      //destructor
  
  Bool_t  IsFolder() const {return kTRUE;};
  void Browse(TBrowser *b); 

  void     Fill(Double_t d1,Double_t d2, TComplex c);   //fill the histograms
  Int_t    GetNbinsX()                
    {Int_t iMaxEtaBins = fHistProReNumer->GetNbinsX();  return iMaxEtaBins;}     
  Int_t    GetNbinsXPt()              
    {Int_t iMaxPtBins = fHistProReNumerPt->GetNbinsX(); return iMaxPtBins;}
  Double_t GetBinCenter(Int_t i)   
    {Double_t dEta = fHistProReNumer->GetXaxis()->GetBinCenter(i);  return dEta;}
  Double_t GetBinCenterPt(Int_t i) 
    {Double_t dPt = fHistProReNumerPt->GetXaxis()->GetBinCenter(i); return dPt;}
  TComplex GetNumerEta(Int_t i);                   //get numerator for diff. flow (eta)
  TComplex GetNumerPt(Int_t i);                    //get numerator for diff. flow (pt)

  TProfile*   GetHistProReNumer()   {return this->fHistProReNumer;}         
  TProfile*   GetHistProImNumer()   {return this->fHistProImNumer;}        
  TProfile*   GetHistProReNumerPt() {return this->fHistProReNumerPt;} 
  TProfile*   GetHistProImNumerPt() {return this->fHistProImNumerPt;} 
  TProfile2D* GetHistProReNumer2D() {return this->fHistProReNumer2D;}       
  TProfile2D* GetHistProImNumer2D() {return this->fHistProImNumer2D;}
  TList*      GetHistList()         {return this->fHistList;}   
  
  virtual Double_t Merge(TCollection *aList);  //merge function
  void Print(Option_t* option = "") const;     //method to print stats

 private:
 
  AliFlowLYZHist2(const AliFlowLYZHist2& aAnalysis);             //copy constructor
  AliFlowLYZHist2& operator=(const AliFlowLYZHist2& aAnalysis);  //assignment operator
  
  TProfile*   fHistProReNumer;         //holds Re of Numerator(eta)
  TProfile*   fHistProImNumer;         //holds Im of Numerator(eta)
  TProfile*   fHistProReNumerPt;       //holds Re of Numerator(pt)
  TProfile*   fHistProImNumerPt;       //holds Im of Numerator(pt)
  TProfile2D* fHistProReNumer2D;       //holds Re of Numerator
  TProfile2D* fHistProImNumer2D;       //holds Im of Numerator
  TList*      fHistList;               //list to hold all histograms  

  ClassDef(AliFlowLYZHist2,1)          // macro for rootcint
    };
 
     
#endif
