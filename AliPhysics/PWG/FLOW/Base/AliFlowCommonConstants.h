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

  enum ERefMultSource {kExternal, kRP, kPOI};

  Int_t GetNbinsMult() const { return fNbinsMult; }
  Int_t GetNbinsPt()   const { return fNbinsPt; }
  Int_t GetNbinsPhi()  const { return fNbinsPhi; }
  Int_t GetNbinsEta()  const { return fNbinsEta; }
  Int_t GetNbinsQ()    const { return fNbinsQ; }
  Int_t GetNbinsMass() const { return fNbinsMass; }
   
  Double_t GetMultMin() const { return fMultMin; }
  Double_t GetMultMax() const { return fMultMax; }
  Double_t GetPtMin()   const { return fPtMin; }
  Double_t GetPtMax()   const { return fPtMax; }
  Double_t GetPhiMin()  const { return fPhiMin; }
  Double_t GetPhiMax()  const { return fPhiMax; }
  Double_t GetEtaMin()  const { return fEtaMin; }
  Double_t GetEtaMax()  const { return fEtaMax; }
  Double_t GetQMin()    const { return fQMin; }
  Double_t GetQMax()    const { return fQMax; }
  Double_t GetMassMin()    const { return fMassMin; }
  Double_t GetMassMax()    const { return fMassMax; }
  Double_t GetHistWeightvsPhiMax() const {return fHistWeightvsPhiMax;}
  Double_t GetHistWeightvsPhiMin() const {return fHistWeightvsPhiMin;}
  
  void SetNbinsMult( Int_t i ) { fNbinsMult = i; }
  void SetNbinsPt( Int_t i )   { fNbinsPt = i; }
  void SetNbinsPhi( Int_t i )  { fNbinsPhi = i; }
  void SetNbinsEta( Int_t i )  { fNbinsEta = i; }
  void SetNbinsQ( Int_t i )    { fNbinsQ = i; }
  void SetNbinsMass( Int_t i ) { fNbinsMass = i; }
   
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
  void SetMassMin( Double_t i )    { fMassMin = i; }
  void SetMassMax( Double_t i )    { fMassMax = i; }
  void SetHistWeightvsPhiMax( Double_t d ) {fHistWeightvsPhiMax=d;}
  void SetHistWeightvsPhiMin( Double_t d ) {fHistWeightvsPhiMin=d;}
  
 private:
  AliFlowCommonConstants& operator= (const AliFlowCommonConstants& c);
  AliFlowCommonConstants(const AliFlowCommonConstants& a);
  
  //histogram sizes
  Int_t  fNbinsMult; // histogram size
  Int_t  fNbinsPt;   // histogram size
  Int_t  fNbinsPhi;  // histogram size
  Int_t  fNbinsEta;  // histogram size
  Int_t  fNbinsQ;    // histogram size
  Int_t  fNbinsMass; // histogram size
 
  // Histograms limits
  Double_t  fMultMin;  // histogram limit 
  Double_t  fMultMax;  // histogram limit
  Double_t  fPtMin;    // histogram limit
  Double_t  fPtMax;    // histogram limit
  Double_t  fPhiMin;	 // histogram limit
  Double_t  fPhiMax;   // histogram limit
  Double_t  fEtaMin;	 // histogram limit
  Double_t  fEtaMax;	 // histogram limit
  Double_t  fQMin;	   // histogram limit
  Double_t  fQMax;     // histogram limit
  Double_t  fMassMin;  // histogram limit 
  Double_t  fMassMax;  // histogram limit
  Double_t  fHistWeightvsPhiMin; // histogram limit
  Double_t  fHistWeightvsPhiMax; // histogram limit
 
  static AliFlowCommonConstants* fgPMasterConfig; //master object
  
  ClassDef(AliFlowCommonConstants,2) //ClassDef
};

#endif


