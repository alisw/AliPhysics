#ifndef AliHMPIDQADataMakerRec_H
#define AliHMPIDQADataMakerRec_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.
//  A. Mastroserio



// --- ROOT system ---
class TTree ;
class AliRawReader;
class AliESDEvent ;
// --- AliRoot header files ---

#include "AliQADataMakerRec.h"

class AliHMPIDQADataMakerRec: public AliQADataMakerRec {

public:
  AliHMPIDQADataMakerRec() ;          // ctor
  AliHMPIDQADataMakerRec(const AliHMPIDQADataMakerRec& qadm) ;   
  AliHMPIDQADataMakerRec& operator = (const AliHMPIDQADataMakerRec& qadm) ;
  virtual ~AliHMPIDQADataMakerRec() {;} // dtor

private:
  virtual void   InitRecPoints();  //book cluster QA histo
  virtual void   InitRaws();     //book raw QA histo
  virtual void   InitESDs() ;      //book ESD QA histo 
  virtual void   MakeRecPoints(TTree * clusters)    ;  //Fill cluster QA histo
  virtual void   MakeRaws(AliRawReader* rawReader);
  virtual void   MakeESDs(AliESDEvent * esd) ;         //Fill hit QA histo
  virtual void   StartOfDetectorCycle() ;
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** obj) ;
  Int_t   fEvtRaw;    //internal event counter for raw 

  ClassDef(AliHMPIDQADataMakerRec,1)  // description 

};

#endif // AliHMPIDQADataMakerRec_H
