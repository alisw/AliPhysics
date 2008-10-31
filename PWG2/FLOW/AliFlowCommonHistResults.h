/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#ifndef AliFlowCommonHistResults_H
#define AliFlowCommonHistResults_H

class TH1F;
class TH1D;
class TCollection;
class TList;
class TBrowser;

// AliFlowCommonHistResults:
// Class to organize the common histograms for Flow Analysis
// authors: N. van der Kolk (kolk@nikhef.nl) and A. Bilandzic (anteb@nikhef.nl)

class AliFlowCommonHistResults : public TNamed {

 public:
  AliFlowCommonHistResults();               //default constructor
  AliFlowCommonHistResults(const char *name,const char *title = "AliFlowCommonHist");  //constructor
  virtual ~AliFlowCommonHistResults();      //destructor

  Bool_t  IsFolder() const {return kTRUE;};
  void Browse(TBrowser *b); 

  //make fill methods here
  Bool_t FillIntegratedFlow(Double_t aV, Double_t anError);                //fill fHistIntFlow
  Bool_t FillDifferentialFlow(Int_t aBin, Double_t av, Double_t anError);  //fill fHistDiffFlow
  Bool_t FillChi(Double_t aChi);                                           //fill fHistChi

  //make get methods here
  TH1D*    GetHistDiffFlow()               {return fHistDiffFlow; } ; 
  TH1D*    GetHistChi()                    {return fHistChi; } ;
  TH1D*    GetHistIntFlow()                {return fHistIntFlow; } ;
  TList*   GetHistList()                   {return fHistList;} ;  

  virtual Double_t  Merge(TCollection *aList);  //merge function
  void Print(Option_t* option = "") const;      //method to print stats

 private:

  AliFlowCommonHistResults(const AliFlowCommonHistResults& aSetOfResultHists);            //copy constructor
  AliFlowCommonHistResults& operator=(const AliFlowCommonHistResults& aSetOfResultHists); //assignment operator

  TH1D*     fHistIntFlow;      //integrated flow
  TH1D*     fHistDiffFlow;     //differential flow
  TH1D*     fHistChi;          //resolution
  TList*    fHistList;         //list to hold all histograms

  ClassDef(AliFlowCommonHistResults,1)                 // macro for rootcint
};
 
#endif

