#ifndef ALITOFPROB_H
#define ALITOFPROB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//_________________________________________________________________________
//  Task Class for Probabilities in TOF      
//                  
//-- Author: F. Pierella


#include "TTask.h"
#include "TCutG.h"
#include "TString.h"
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>

class AliTOFProb: public TTask {

public:
  AliTOFProb() ;          // ctor
  AliTOFProb(char* headerFile) ; 
  virtual ~AliTOFProb() ; // dtor   
  virtual void  Exec(const Option_t *dummyOpt); // do the main work
  void          Init();

  TNtuple*      GetInputNtuple() const {return fNtuple;}

  Bool_t   operator == (const AliTOFProb & tofprob) const ;

private:

  TNtuple *fNtuple; // pointer to ntuple 
  TFile *fhfile;
  TFile *fgen;
  char* foutfileName; // destination file name for histos
  Int_t fTask;

 protected:

  ClassDef(AliTOFProb,0)  // Task class for PID probabilities from TOF

};

#endif // AliTOFPROB_H
