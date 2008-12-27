#ifndef AliTRDQADataMakerSim_H
#define AliTRDQADataMakerSim_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Produces the data needed to calculate the quality assurance.          //
//  All data must be mergeable objects.                                   //
//                                                                        //
//  Author:                                                               //
//    Sylwester Radomski (radomski@physi.uni-heidelberg.de)               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
class TH1F ; 
class TH1I ; 

// --- AliRoot header files ---
class AliExternalTrackParam;

#include "AliQADataMakerSim.h"

class AliTRDQADataMakerSim: public AliQADataMakerSim {

 public:

  AliTRDQADataMakerSim() ;          // ctor
  AliTRDQADataMakerSim(const AliTRDQADataMakerSim& qadm) ;   
  AliTRDQADataMakerSim& operator = (const AliTRDQADataMakerSim& qadm) ;
  virtual ~AliTRDQADataMakerSim() {;} // dtor

 private:

  virtual void EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray ** list) ;
  virtual void InitHits() ; 
  virtual void InitDigits() ; 
  virtual void InitSDigits() ;

  virtual void MakeHits(TTree * hitTree);
  virtual void MakeHits(TClonesArray * hits);

  virtual void MakeSDigits(TTree *sdigitTree);
  virtual void MakeSDigits(TClonesArray * sigits); 

  virtual void MakeDigits(TTree *digitTree);
  virtual void MakeDigits(TClonesArray * digits); 

  virtual void StartOfDetectorCycle() ; 
  Int_t    CheckPointer(TObject *obj, const char *name);

  ClassDef(AliTRDQADataMakerSim,1)   // Creates the TRD QA data

};
#endif // AliTRDQADataMakerSim_H
