#ifndef ALIANACALOTRIGGER_H
#define ALIANACALOTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the trigger data in ESD
// Creates an ntuple for 2x2 and NxN triggers
// Each ntuple connects the maximum trigger amplitudes 
// and its positions with reconstructed clusters
// and if MC stack available, with pt of parent.
//
//*-- Yves Schutz (CERN) & Gustavo Conesa Balbastre (INFN-LNF)
//////////////////////////////////////////////////////////////////////////////


#include "AliAnalysisTaskSE.h"  
class TNtuple ;

class AliESDEvent ; 

class AliAnaCaloTrigger : public AliAnalysisTaskSE {

public:
  AliAnaCaloTrigger() ;
  AliAnaCaloTrigger(const char *name) ;
  AliAnaCaloTrigger(const AliAnaCaloTrigger & trig) ;
  AliAnaCaloTrigger & operator=(const AliAnaCaloTrigger& source);
  virtual ~AliAnaCaloTrigger() ;
   
  virtual void UserExec(Option_t * opt = "") ;
  virtual void UserCreateOutputObjects();
//  virtual void Terminate(Option_t * opt = "") const ;

  TString GetCalorimeter()     const   {return fCalorimeter ; }
  void    SetCalorimeter(TString calo) {fCalorimeter = calo ; }

private:
 
  TList * fOutputContainer ; //! output data container
  TString fCalorimeter ; // "PHOS" or "EMCAL"

  // Histograms
  TNtuple * fNtTrigger22 ; //Ntuple with 2x2 max dig amplitude and cluster energy, and positions.
  TNtuple * fNtTriggerNN ; //Ntuple with NxN max dig amplitude and cluster energy, and positions.

  ClassDef(AliAnaCaloTrigger, 2); // a trigger analysis task 
};
#endif // ALIANACALOTRIGGER_H
