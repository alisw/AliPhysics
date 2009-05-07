#ifndef AliTRDQADATAMAKER_H
#define AliTRDQADATAMAKER_H
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

#include "AliQADataMaker.h"

class AliTRDQADataMaker: public AliQADataMaker {

 public:

  AliTRDQADataMaker() ;          // ctor
  AliTRDQADataMaker(const AliTRDQADataMaker& qadm) ;   
  AliTRDQADataMaker& operator = (const AliTRDQADataMaker& qadm) ;
  virtual ~AliTRDQADataMaker() {;} // dtor

 private:

  virtual void EndOfDetectorCycle(AliQAv1::TASKINDEX task, TObjArray * list) ;
  virtual void InitHits() ; 
  virtual void InitESDs() ; 
  virtual void InitDigits() ; 
  virtual void InitRecPoints() ; 
  virtual void InitRaws() ; 
  virtual void InitSDigits() ;

  virtual void MakeHits(TTree * hitTree);
  virtual void MakeHits(TClonesArray * hits);

  virtual void MakeSDigits(TTree *sdigitTree);
  virtual void MakeSDigits(TClonesArray * sigits); 

  virtual void MakeDigits(TTree *digitTree);
  virtual void MakeDigits(TClonesArray * digits); 

  virtual void MakeRaws(AliRawReader* rawReader); 
  virtual void MakeRecPoints(TTree * recpo); 
  virtual void MakeESDs(AliESDEvent * esd);

  virtual void StartOfDetectorCycle() ; 
  Int_t    CheckPointer(TObject *obj, const char *name);

  // internal methods
  Int_t    GetSector(const Double_t alpha) const;
  Double_t GetExtZ(const AliExternalTrackParam *paramIn) const;

  ClassDef(AliTRDQADataMaker,1)   // Creates the TRD QA data

};
#endif // AliTRDQADATAMAKER_H
