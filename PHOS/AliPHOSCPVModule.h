#ifndef ALIPHOSCPVMODULE_H
#define ALIPHOSCPVMODULE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for one CPV module          //
//                                            //
//  Author: Yuri Kharlov, IHEP, Protvino      //
//  e-mail: Yuri.Kharlov@cern.ch              //
//  Last modified: 7 November 2000            //
////////////////////////////////////////////////
 
#include <assert.h> 

// --- ROOT system ---
#include <TClonesArray.h> 
#include <TLorentzVector.h>

// --- galice header files ---

class AliPHOSCPVModule : public TObject {

public:

  virtual ~AliPHOSCPVModule(void);
           AliPHOSCPVModule(void);
           AliPHOSCPVModule(const AliPHOSCPVModule & cpv) {
	     // cpy ctor requested by Coding Convention but not yet needed
	     assert(0==1) ;
	   } 
  
  AliPHOSCPVModule & operator = (const AliPHOSCPVModule  & rvalue) {
    // assignement operator requested by coding convention but not needed
    assert(0==1) ;
    return *this ; 
  }
  void     Clear(Option_t *opt="");
  void     Print(Option_t *opt="");
  void     AddHit(TLorentzVector p, Float_t *xy, Int_t ipart);
  void     MakeBranch(Int_t i);
  void     SetTreeAddress(Int_t i);
  
  TClonesArray *Hits(void) {return fHits;}

private:
  
  TClonesArray *fHits;              // List of hits in the Module

  ClassDef(AliPHOSCPVModule,1)      // CPV Module
};
 
#endif
