/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#ifndef AliFlowCommonHistResults_H
#define AliFlowCommonHistResults_H

class TH1F;
class TH1D;

// AliFlowCommonHistResults:
// Class to organize the common histograms for Flow Analysis
// authors: N. van der Kolk (kolk@nikhef.nl) and A. Bilandzic (anteb@nikhef.nl)

class AliFlowCommonHistResults : public TObject {

 public:

  AliFlowCommonHistResults(TString input);  //constructor
  virtual ~AliFlowCommonHistResults();      //destructor

  //make fill methods here
  Bool_t FillIntegratedFlow(Double_t aV, Double_t anError);                //fill fHistIntFlow
  Bool_t FillDifferentialFlow(Int_t aBin, Double_t av, Double_t anError);  //fill fHistDiffFlow
  Bool_t FillChi(Double_t aChi);                                           //fill fHistChi

  //make get methods here
  TH1D*    GetfHistDiffFlow()               {return fHistDiffFlow; } ; 
  TH1D*    GetfHistChi()                    {return fHistChi; } ;
  TH1D*    GetfHistIntFlow()                {return fHistIntFlow; } ;
 
 private:

  AliFlowCommonHistResults(const AliFlowCommonHistResults& aSetOfResultHists);            //copy constructor
  AliFlowCommonHistResults& operator=(const AliFlowCommonHistResults& aSetOfResultHists); //assignment operator

  //integrated flow
  TH1D*     fHistIntFlow;      
  //differential flow
  TH1D*     fHistDiffFlow; 
  //resolution
  TH1D*     fHistChi;
  
  ClassDef(AliFlowCommonHistResults,0)                 // macro for rootcint
};
 
#endif

