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
  void SetDebug(Int_t debug = 0){fDebug = debug; if (fOutputPointer) fOutputPointer->SetDebug(debug);}
  AliEMCALJetFinderOutput* GetOutput() { return fOutputPointer;}
  void SetOutput(AliEMCALJetFinderOutput *output);
  Float_t PropagatePhi(Float_t pt, Float_t charge, Bool_t& curls);  	
  Bool_t GetPythiaComparison(){return fPythiaComparison;}
  void SetPythiaComparison(Bool_t value){fPythiaComparison=value;}

  AliEMCALJetFinderAlgo (const AliEMCALJetFinderAlgo&);
  AliEMCALJetFinderAlgo & operator = (const AliEMCALJetFinderAlgo & ) {
    Fatal("operator =", "not implemented") ;
    return *this ;
  }

protected:
   
  AliEMCALJetFinderInput*        fInputPointer;  // pointer to the input object 
  AliEMCALJetFinderOutput*       fOutputPointer;  // output object for results
  Bool_t                         fOutputAllocated; // flags if ouputobject is woned by JetFinderAlgo
  Int_t				 fDebug; 	 // debug level
  Bool_t 			 fPythiaComparison; // for pyclus comparison

  ClassDef(AliEMCALJetFinderAlgo,4)

};
#endif

