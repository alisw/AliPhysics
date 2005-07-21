#ifndef ALIJETKINEREADERHEADER_H
#define ALIJETKINEREADERHEADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Jet Kinematics Reader Header
// Header for the Kine reader in the jet analysis
// Author: Andreas Morsch (andreas.morsch@cern.ch)

#include "AliJetReaderHeader.h"
 
class AliJetKineReaderHeader : public AliJetReaderHeader
{

 public:
  AliJetKineReaderHeader();
  virtual ~AliJetKineReaderHeader();
  
  // Setters
  void SetFastSimTPC(Bool_t flag = kTRUE) {fFastSimTPC = flag;}
  void SetFastSimEMCAL(Bool_t flag = kTRUE) {fFastSimEMCAL = flag;}
  virtual void SetPtCut(Float_t par = 2.0) {fPtCut = par;}
  // Getter
  Bool_t  FastSimTPC() const  {return fFastSimTPC;}
  Bool_t  FastSimEMCAL() const  {return fFastSimEMCAL;}
  Float_t GetPtCut()   const  {return fPtCut;}
	  
 protected:
  //parameters set by user
  Bool_t fFastSimTPC;
  Bool_t fFastSimEMCAL;
  Float_t fPtCut;
  ClassDef(AliJetKineReaderHeader,1);
};
 
#endif
