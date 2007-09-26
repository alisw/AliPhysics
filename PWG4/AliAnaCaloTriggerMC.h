#ifndef ALIANACALOTRIGGERMC_H
#define ALIANACALOTRIGGERMC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the PHOS photon data in simulated data
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"  

class AliESDEvent ; 
class TNtuple ;
class TH1D ; 
class TH1I ; 
class TChain;

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
  virtual void Terminate(Option_t * opt = "") ;

  TString GetCalorimeter()     const   {return fCalorimeter ; }
  void    SetCalorimeter(TString calo) {fCalorimeter = calo ; }

private:
  TChain   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESDEvent  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container

  TString fCalorimeter ; // "PHOS" or "EMCAL"

  // Histograms
  TNtuple * fNtTrigger22 ;
  TNtuple * fNtTriggerNN ;

  ClassDef(AliAnaCaloTriggerMC, 0); // a PHOS photon analysis task 
};
#endif // ALIANACALOTRIGGERMC_H
