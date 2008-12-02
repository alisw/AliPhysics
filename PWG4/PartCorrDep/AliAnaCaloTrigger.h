#ifndef ALIANACALOTRIGGER_H
#define ALIANACALOTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the trigger data in ESD
// Creates an ntuple for 2x2 and NxN triggers
// Each ntuple connects the maximum trigger amplitudes 
// and its positions with reconstructed clusters
//
//*-- Yves Schutz (CERN) & Gustavo Conesa Balbastre (INFN-LNF)
//////////////////////////////////////////////////////////////////////////////


#include "AliAnalysisTask.h"  
class TFile ;
class TNtuple ;
class TH1D ; 
class TH1I ; 
class TChain;

class AliAnalysisManager ;
class AliESDEvent ; 

class AliAnaCaloTrigger : public AliAnalysisTask {

public:
  AliAnaCaloTrigger() ;
  AliAnaCaloTrigger(const char *name) ;
  AliAnaCaloTrigger(const AliAnaCaloTrigger & trig) ;
  AliAnaCaloTrigger & operator=(const AliAnaCaloTrigger& source);
  virtual ~AliAnaCaloTrigger() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
//  virtual void Terminate(Option_t * opt = "") const ;

  TString GetCalorimeter()     const   {return fCalorimeter ; }
  void    SetCalorimeter(TString calo) {fCalorimeter = calo ; }

private:
  TChain       * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESDEvent  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container

  TString fCalorimeter ; // "PHOS" or "EMCAL"

  // Histograms
  TNtuple * fNtTrigger22 ; //Ntuple with 2x2 max dig amplitude and cluster energy, and positions.
  TNtuple * fNtTriggerNN ; //Ntuple with NxN max dig amplitude and cluster energy, and positions.

  ClassDef(AliAnaCaloTrigger, 1); // a photon analysis task 
};
#endif // ALIANACALOTRIGGER_H
