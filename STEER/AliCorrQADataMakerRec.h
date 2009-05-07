#ifndef ALICORRQADataMakerRec_H
#define ALICORRQADataMakerRec_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id: AliCORRQADataMakerRec.h 27570 2008-07-24 21:49:27Z cvetan $ */

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

class AliCorrQADataMakerRec: public AliQADataMakerRec {

public:
  AliCorrQADataMakerRec(AliQADataMaker **) ;          // ctor
  AliCorrQADataMakerRec(const AliCorrQADataMakerRec& qadm) ;   
  AliCorrQADataMakerRec& operator = (const AliCorrQADataMakerRec& qadm) ;
  virtual ~AliCorrQADataMakerRec() ; // dtor
  
private:

  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** list) ;
  virtual void   InitESDs() ; 
  virtual void   InitRecPoints() ; 
  virtual void   InitRaws() ; 
  virtual void   MakeESDs(AliESDEvent * esd) ;
  virtual void   MakeRecPoints(TTree * recpoTree) ; 
  virtual void   MakeRaws(AliRawReader *) ; 
  virtual void   StartOfDetectorCycle() ; 

  Int_t fMaxRawVar ;              //! number of raw parameters in the ntuple
  AliQADataMaker **    fqadm ;    //! array of detectors QA data makers pointers
  ClassDef(AliCorrQADataMakerRec,1)  // description 

};

#endif // AliCORRQADataMakerRec_H
