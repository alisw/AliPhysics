#ifndef ALIEMCALJETFINDERALGO_H
#define ALIEMCALJETFINDERALGO_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *  * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Base Class for JetFinder Algorithms     
//                  
//*-- Author: Mark Horner (LBL/UCT)


#include "TTask.h"
#include "AliEMCALJetFinderInput.h"
#include "AliEMCALJetFinderOutput.h"


class AliEMCALJetFinderAlgo : public TTask
{
public:
  AliEMCALJetFinderAlgo();
  virtual ~AliEMCALJetFinderAlgo();
  virtual void FindJets() = 0 ;
  void InitInput(AliEMCALJetFinderInput* input);
  void SetDebug(Int_t debug = 0){fDebug = debug; fOutputObject.SetDebug(debug);}
  AliEMCALJetFinderOutput* GetOutput() { return &fOutputObject;}
  Float_t PropagatePhi(Float_t pt, Float_t charge, Bool_t& curls);  	
protected:
   
  AliEMCALJetFinderInput*        fInputPointer;  // pointer to the input object 
  AliEMCALJetFinderOutput        fOutputObject;  // output object for results
  Int_t				 fDebug; 	 // debug level

  ClassDef(AliEMCALJetFinderAlgo,2)

};
#endif

