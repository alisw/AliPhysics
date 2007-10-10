#ifndef ALIPHOSQADATAMAKER_H
#define ALIPHOSQADATAMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  Y. Schutz CERN July 2007
*/


// --- ROOT system ---
class TH1F ; 
class TH1I ; 
class TList ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQADataMaker.h"

class AliPHOSQADataMaker: public AliQADataMaker {

public:
  AliPHOSQADataMaker() ;          // ctor
  AliPHOSQADataMaker(const AliPHOSQADataMaker& qadm) ;   
  AliPHOSQADataMaker& operator = (const AliPHOSQADataMaker& qadm) ;
  virtual ~AliPHOSQADataMaker() {;} // dtor
  
private:
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX, TList * list) ;
  virtual void   InitHits() ; 
  virtual void   InitESDs() ; 
  virtual void   InitDigits() ; 
  //virtual void   InitRecParticles() ; 
  virtual void   InitRecPoints() ; 
  virtual void   InitRaws() ; 
  //virtual void   InitTrackSegments() ; 
  virtual void   InitSDigits() ; 
  virtual void   MakeESDs(AliESDEvent * esd) ;
  virtual void   MakeHits(TClonesArray * hits) ;
  virtual void   MakeDigits(TClonesArray * digits) ; 
  // virtual void   MakeRecParticles(TTree * recpar) ; 
  virtual void   MakeRecPoints(TTree * recpo) ; 
  virtual void   MakeRaws(AliRawReader* rawReader) ; 
  virtual void   MakeSDigits(TClonesArray * sigits) ; 
  //virtual void   MakeTrackSegments(TTree *ts ) ; 
  virtual void   StartOfDetectorCycle() ; 

  ClassDef(AliPHOSQADataMaker,1)  // description 

};

#endif // AliPHOSQADataMaker_H
