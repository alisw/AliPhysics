#ifndef ALIHMPIDQADATAMAKER_H
#define ALIHMPIDQADATAMAKER_H
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

#include "AliQADataMaker.h"

class AliHMPIDQADataMaker: public AliQADataMaker {

public:
  AliHMPIDQADataMaker() ;          // ctor
  AliHMPIDQADataMaker(const AliHMPIDQADataMaker& qadm) ;   
  AliHMPIDQADataMaker& operator = (const AliHMPIDQADataMaker& qadm) ;
  virtual ~AliHMPIDQADataMaker() {;} // dtor

private:

  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX, TObjArray * obj) ;
  virtual void   InitHits() ;      //book hit QA histo 
  virtual void   InitDigits() ;    //book Digit QA histo
  virtual void   InitSDigits() ;   //book SDigits QA histo
  virtual void   InitRecPoints();  //book cluster QA histo
  virtual void   InitRaws();     //book raw QA histo
  virtual void   InitESDs() ;      //book ESD QA histo 
  virtual void   MakeHits(TClonesArray * hits) ;       //Fill hit QA histo
  virtual void   MakeHits(TTree * hits) ;       // access to hit tree
  virtual void   MakeDigits(TClonesArray * digits) ;   //Fill Digit QA histo
  virtual void   MakeDigits(TTree * digits) ;   //access to digit tree
  virtual void   MakeSDigits(TClonesArray * sdigits) ; //Fill SDigit QA histo
  virtual void   MakeSDigits(TTree * sdigits) ; //access to SDigits tree
  virtual void   MakeRecPoints(TTree * clusters)    ;  //Fill cluster QA histo
  virtual void   MakeRaws(AliRawReader* rawReader);
  virtual void   MakeESDs(AliESDEvent * esd) ;         //Fill hit QA histo
  virtual void   StartOfDetectorCycle() ;

  ClassDef(AliHMPIDQADataMaker,1)  // description 

};

#endif // AliHMPIDQADataMaker_H
