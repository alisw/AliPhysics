#ifndef AliHMPIDQADataMakerSim_H
#define AliHMPIDQADataMakerSim_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.
//  A. Mastroserio



// --- ROOT system ---


class TH1F ; 
class TH2F ;
class TProfile;
class TObjArray;
// --- Standard library ---
#include <TString.h>
// --- AliRoot header files ---

#include "AliQADataMakerSim.h"

class AliHMPIDQADataMakerSim: public AliQADataMakerSim {

public:
  AliHMPIDQADataMakerSim() ;          // ctor
  AliHMPIDQADataMakerSim(const AliHMPIDQADataMakerSim& qadm) ;   
  AliHMPIDQADataMakerSim& operator = (const AliHMPIDQADataMakerSim& qadm) ;
  virtual ~AliHMPIDQADataMakerSim() {;} // dtor

private:

  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX_t, TObjArray ** obj) ;
  virtual void   InitHits() ;      //book hit QA histo 
  virtual void   InitDigits() ;    //book Digit QA histo
  virtual void   InitSDigits() ;   //book SDigits QA histo
  virtual void   MakeHits(TClonesArray * hits) ;       //Fill hit QA histo
  virtual void   MakeHits(TTree * hits) ;       // access to hit tree
  virtual void   MakeDigits(TClonesArray * digits) ;   //Fill Digit QA histo
  virtual void   MakeDigits(TTree * digits) ;   //access to digit tree
  virtual void   MakeSDigits(TClonesArray * sdigits) ; //Fill SDigit QA histo
  virtual void   MakeSDigits(TTree * sdigits) ; //access to SDigits tree
  virtual void   StartOfDetectorCycle() ;

  ClassDef(AliHMPIDQADataMakerSim,1)  // description 

};

#endif // AliHMPIDQADataMakerSim_H
