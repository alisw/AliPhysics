#ifndef ALITOFPID_H
#define ALITOFPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//_________________________________________________________________________
//  Task Class for PID in TOF      
//                  
//-- Authors: B. Zagreev, F. Pierella


#include "TTask.h"
#include "TCutG.h"
#include "TString.h"
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>

class AliTOFPID: public TTask {

public:
  AliTOFPID() ;          // ctor
  AliTOFPID(char* headerFile, char *cutsFile, const Option_t* opt="pp") ; 
  virtual ~AliTOFPID() ; // dtor   
  virtual void  Exec(const Option_t *eventType, const Option_t *outputmode, const Option_t *outputsavemode); // do the main work
  void          Init(const Option_t* opt);
  void          SetNEvents(Int_t Nevents) {fNevents = Nevents;}
  Int_t         GetNEvents() const {return fNevents;}
  void  SetDbg(Int_t dbgflag)                        {fdbg=dbgflag;}
  // getter for pointers to cuts
  TCutG*        GetElectronCut() const {return felectron;}
  TCutG*        GetPionCut()     const {return fpion;}
  TCutG*        GetKaonCut()     const {return fkaon;}
  TCutG*        GetProtonCut()   const {return fproton;}

  TFile *       GetFileCuts()    const {return fcut;}
  TNtuple*      GetInputNtuple() const {return fNtuple;}

  Bool_t   operator == (const AliTOFPID & tofpid) const ;

private:
  TCutG *felectron; // pointer to cut for electron
  TCutG *fpion;     // pointer to cut for pion
  TCutG *fkaon;     // pointer to cut for kaon
  TCutG *fproton;   // pointer to cut for proton
  TFile *fcut;      // pointer to file containing cuts
  TNtuple *fNtuple; // pointer to ntuple 
  TFile *fhfile;
  TFile *fgen;
  char* foutfileName; // destination file name for histos
  Int_t fNevents; // number of events
  Int_t fdbg;
  Int_t fTask;

 protected:

  ClassDef(AliTOFPID,0)  // Task class for TOF pid

};

#endif // AliTOFPID_H
