/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliFlowLYZHist1.h 29627 2008-10-31 10:30:41Z snelling $ */

#ifndef AliFlowLYZHist1_H
#define AliFlowLYZHist1_H


#include "TComplex.h"  //explicitly called in void Fill(Float_t f, TComplex C);

class TH1D;
class TCollection;
class TList;
class TBrowser;

// Description: Class to organise histograms for Flow 
//              by the LeeYangZeros method in the first run.
//              Also includes methods to find the first minimum R0
//              in the generating function.


class AliFlowLYZHist1: public TNamed  {

 public:

  AliFlowLYZHist1();                          //default constructor
  AliFlowLYZHist1(Int_t theta, const char *name ,const char *title = "AliFlowLYZHist1");               //constructor
  virtual  ~AliFlowLYZHist1();                //destructor
  
  Bool_t  IsFolder() const {return kTRUE;};
  void Browse(TBrowser *b); 

  void Fill(Double_t f, TComplex c);          //fill the histograms
  TH1D* FillGtheta();                         //fills fHistGtheta
  Double_t GetR0();                           //get R0
  Double_t GetBinCenter(Int_t i);             //Get a bincentre of fHistGtheta
  Int_t GetNBins();                           //Gets Nbins

  TH1D*     GetHistGtheta()      {return this->fHistGtheta;} ;
  TProfile* GetHistProReGtheta() {return this->fHistProReGtheta;} ;
  TProfile* GetHistProImGtheta() {return this->fHistProImGtheta;} ;
  TList*    GetHistList()        {return this->fHistList;} ;  

  virtual Double_t  Merge(TCollection *aList);  //merge function
  void Print(Option_t* option = "") const;      //method to print stats

 
private:

  AliFlowLYZHist1(const AliFlowLYZHist1& aAnalysis);             //copy constructor
  AliFlowLYZHist1& operator=(const AliFlowLYZHist1& aAnalysis);  //assignment operator

  TH1D*     fHistGtheta;             //holds |Gtheta|^2(r)
  TProfile* fHistProReGtheta;        //holds Re of Gtheta(r)
  TProfile* fHistProImGtheta;        //holds Im of Gtheta(r)
  TList*    fHistList;               //list to hold all histograms  


  ClassDef(AliFlowLYZHist1,1)        // macro for rootcint
    };
 
     
#endif
