#ifndef ALIANACALOTRIGGERMC_H
#define ALIANACALOTRIGGERMC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the trigger data in ESD with MC data
// Creates an ntuple for 2x2 and NxN triggers
// Each ntuple connects the maximum trigger amplitudes 
// and its positions with reconstructed clusters and MC
//
//*-- Yves Schutz (CERN) & Gustavo Conesa Balbastre (INFN-LNF)
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"  

class TFile ;
class TNtuple ;
class TH1D ; 
class TH1I ; 
class TChain;

class AliMCEvent ;
class AliESDEvent ; 

class AliAnaCaloTriggerMC : public AliAnalysisTask {

public:
  AliAnaCaloTriggerMC();
  AliAnaCaloTriggerMC(const char *name) ;
  AliAnaCaloTriggerMC(const AliAnaCaloTriggerMC & trig) ;
  AliAnaCaloTriggerMC & operator=(const AliAnaCaloTriggerMC& source);
  virtual ~AliAnaCaloTriggerMC() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
//  virtual void Terminate(Option_t * opt = "") const ;

  TString GetCalorimeter()     const   {return fCalorimeter ; }
  void    SetCalorimeter(TString calo) {fCalorimeter = calo ; }

private:
  TChain   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESDEvent  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container

  TString fCalorimeter ; // "PHOS" or "EMCAL"

  // Histograms
  TNtuple * fNtTrigger22 ; //Ntuple with 2x2 max dig amplitude and MC particle and cluster energy, and positions.
  TNtuple * fNtTriggerNN ; //Ntuple with NxN max dig amplitude  and MC particle and cluster energy, and positions.

  ClassDef(AliAnaCaloTriggerMC, 1); // a photon analysis task 
};
#endif // ALIANACALOTRIGGERMC_H
