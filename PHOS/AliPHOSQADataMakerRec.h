#ifndef ALIPHOSQADataMakerRec_H
#define ALIPHOSQADataMakerRec_H
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
class TObjArray ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQADataMakerRec.h"

class AliPHOSQADataMakerRec: public AliQADataMakerRec {

public:
  AliPHOSQADataMakerRec() ;          // ctor
  AliPHOSQADataMakerRec(const AliPHOSQADataMakerRec& qadm) ;   
  AliPHOSQADataMakerRec& operator = (const AliPHOSQADataMakerRec& qadm) ;
  virtual ~AliPHOSQADataMakerRec() {;} // dtor
  
private:
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX, TObjArray * list) ;
  virtual void   InitESDs() ; 
  //virtual void   InitRecParticles() ; 
  virtual void   InitRecPoints() ; 
  virtual void   InitRaws() ; 
  //virtual void   InitTrackSegments() ; 
  virtual void   MakeESDs(AliESDEvent * esd) ;
  // virtual void   MakeRecParticles(TTree * recpar) ; 
  virtual void   MakeRecPoints(TTree * recpoTree) ; 
  virtual void   MakeRaws(AliRawReader* rawReader) ; 
  //virtual void   MakeTrackSegments(TTree *ts ) ; 
  virtual void   StartOfDetectorCycle() ; 

  ClassDef(AliPHOSQADataMakerRec,1)  // description 

};

#endif // AliPHOSQADataMakerRec_H
