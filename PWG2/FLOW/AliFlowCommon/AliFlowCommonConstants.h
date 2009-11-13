/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef ALIFLOWCOMMONCONSTANTS_H
#define ALIFLOWCOMMONCONSTANTS_H

// AliFlowCommonConstants:
// Description: constants for the Common Histograms in the Flow Analysis
// Author: Naomi van der Kolk (kolk@nikhef.nl)
// modified: Mikolaj Krzewicki, Nikhef (mikolaj.krzewicki@cern.ch)

//class TNamed;

#include <TNamed.h>

class AliFlowCommonConstants: public TNamed {

 public:
  AliFlowCommonConstants();
  virtual ~AliFlowCommonConstants();
  static AliFlowCommonConstants* GetMaster();

  Int_t GetNbinsMult() { return fNbinsMult; }
  Int_t GetNbinsPt()   { return fNbinsPt; }
  Int_t GetNbinsPhi()  { return fNbinsPhi; }
  Int_t GetNbinsEta()  { return fNbinsEta; }
  Int_t GetNbinsQ()    { return fNbinsQ; }
   
  Double_t GetMultMin() { return fMultMin; }
  Double_t GetMultMax() { return fMultMax; }
  Double_t GetPtMin()   { return fPtMin; }
  Double_t GetPtMax()   { return fPtMax; }
  Double_t GetPhiMin()  { return fPhiMin; }
  Double_t GetPhiMax()  { return fPhiMax; }
  Double_t GetEtaMin()  { return fEtaMin; }
  Double_t GetEtaMax()  { return fEtaMax; }
  Double_t GetQMin()    { return fQMin; }
  Double_t GetQMax()    { return fQMax; }
  
  void SetNbinsMult( Int_t i ) { fNbinsMult = i; }
  void SetNbinsPt( Int_t i )   { fNbinsPt = i; }
  void SetNbinsPhi( Int_t i )  { fNbinsPhi = i; }
  void SetNbinsEta( Int_t i )  { fNbinsEta = i; }
  void SetNbinsQ( Int_t i )    { fNbinsQ = i; }
   
  void SetMultMin( Double_t i ) { fMultMin = i; }
  void SetMultMax( Double_t i ) { fMultMax = i; }
  void SetPtMin( Double_t i )   { fPtMin = i; }
  void SetPtMax( Double_t i )   { fPtMax = i; }
  void SetPhiMin( Double_t i )  { fPhiMin = i; }
  void SetPhiMax( Double_t i )  { fPhiMax = i; }
  void SetEtaMin( Double_t i )  { fEtaMin = i; }
  void SetEtaMax( Double_t i )  { fEtaMax = i; }
  void SetQMin( Double_t i )    { fQMin = i; }
  void SetQMax( Double_t i )    { fQMax = i; }
  
 private:
  AliFlowCommonConstants& operator= (const AliFlowCommonConstants& c);
  AliFlowCommonConstants(const AliFlowCommonConstants& a);
  
  //histogram sizes
  Int_t  fNbinsMult; //
  Int_t  fNbinsPt;   //
  Int_t  fNbinsPhi;  //
  Int_t  fNbinsEta;  //
  Int_t  fNbinsQ;    //
 
  // Histograms limits
  Double_t  fMultMin;  //          
  Double_t  fMultMax;  //
  Double_t  fPtMin;    //
  Double_t  fPtMax;    //
  Double_t  fPhiMin;	 //  
  Double_t  fPhiMax;   //
  Double_t  fEtaMin;	 //  
  Double_t  fEtaMax;	 //    
  Double_t  fQMin;	   //  
  Double_t  fQMax;     //
 
  static AliFlowCommonConstants* fgPMasterConfig;
  
  ClassDef(AliFlowCommonConstants,1)  // macro for rootcint
};

#endif


