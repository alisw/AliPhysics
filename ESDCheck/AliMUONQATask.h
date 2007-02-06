#ifndef ALIMUONQATASK_H
#define ALIMUONQATASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the MUON data in simulated data
//
//*-- Ivana Hrivnacova
//////////////////////////////////////////////////////////////////////////////

#include <TTree.h> 
#include "AliAnalysisTask.h"  

class AliESD ; 
class TH1F ;
 
class AliMUONQATask : public AliAnalysisTask {

public:
  AliMUONQATask(const char *name) ;
  virtual ~AliMUONQATask() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "") ;

private:
  TTree   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESD  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container

  Int_t fnTrackTrig ; //!
  Int_t ftracktot   ; //!
  Int_t fnevents    ; //!
  Int_t fSPLowpt    ; //!
  Int_t fSPHighpt   ; //!
  Int_t fSPAllpt    ; //!
  Int_t fSMLowpt    ; //!
  Int_t fSMHighpt   ; //!
  Int_t fSMAllpt    ; //!
  Int_t fSULowpt    ; //!
  Int_t fSUHighpt   ; //!
  Int_t fSUAllpt    ; //!
  Int_t fUSLowpt    ; //!
  Int_t fUSHighpt   ; //!
  Int_t fUSAllpt    ; //! 
  Int_t fLSLowpt    ; //!
  Int_t fLSHighpt   ; //! 
  Int_t fLSAllpt    ; //!

  // Histograms
  TH1F * fhMUONVertex ; //! 
  TH1F * fhMUONMult   ; //!
   
  ClassDef(AliMUONQATask, 0); // a MUON photon analysis task 
};
#endif // ALIMUONQATASK_H
