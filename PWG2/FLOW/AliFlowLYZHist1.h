/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef AliFlowLYZHist1_H
#define AliFlowLYZHist1_H


#include "TComplex.h"  //explicitly called in void Fill(Float_t f, TComplex C);

class TProfile;
class TH1D;

// Description: Class to organise histograms for Flow 
//              by the LeeYangZeros method in the first run.
//              Also includes methods to find the first minimum R0
//              in the generating function.


class AliFlowLYZHist1 {

 public:

  AliFlowLYZHist1(Int_t i1, Int_t i2);        //constructor
  AliFlowLYZHist1(const AliFlowLYZHist1&);    //copy constructor (dummy)
  virtual  ~AliFlowLYZHist1();                //destructor

  AliFlowLYZHist1& operator=(const AliFlowLYZHist1&); //assignment operator (dummy)

  void Fill(Float_t f, TComplex C);           //fill the histograms
  TH1D* FillGtheta();                         //fills fHistGtheta
  Float_t GetR0();                            //get R0
  Float_t GetBinCenter(Int_t i);              //Get a bincentre of fHistGtheta
  Int_t GetNBins();                           //Gets fNbins

  //void Save();                              //save the histograms 

 
private:

  TH1D* fHistGtheta;                          //!
  TProfile* fHistProReGtheta;                 //!
  TProfile* fHistProImGtheta;                 //!
  

  ClassDef(AliFlowLYZHist1,0)                 // macro for rootcint
    };
 
     
#endif

