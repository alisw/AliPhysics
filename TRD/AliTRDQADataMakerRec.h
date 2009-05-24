#ifndef AliTRDQADataMakerRec_H
#define AliTRDQADataMakerRec_H
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

#include "AliQADataMakerRec.h"

class AliTRDQADataMakerRec: public AliQADataMakerRec {

 public:

  AliTRDQADataMakerRec() ;          // ctor
  AliTRDQADataMakerRec(const AliTRDQADataMakerRec& qadm) ;   
  AliTRDQADataMakerRec& operator = (const AliTRDQADataMakerRec& qadm) ;
  virtual ~AliTRDQADataMakerRec() {;} // dtor

 private:

  virtual void EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list) ;
  virtual void InitESDs() ; 
  virtual void InitDigits() ; 
  virtual void InitRecPoints() ; 
  virtual void InitRaws() ; 

  virtual void MakeRaws(AliRawReader* rawReader); 
  virtual void MakeDigits(TTree * dig); 
  virtual void MakeDigits(TClonesArray * digits); 
  virtual void MakeRecPoints(TTree * recpo); 
  virtual void MakeESDs(AliESDEvent * esd);

  virtual void StartOfDetectorCycle() ; 
  Int_t    CheckPointer(TObject *obj, const char *name);

  // internal methods
  Int_t    GetSector(Double_t alpha) const;
  Int_t    GetStack(const AliExternalTrackParam *paramOut) const;
  Double_t GetExtZ(const AliExternalTrackParam *paramIn) const;
  Int_t    FillBits(TH1D *hist, Int_t code, Int_t offset);

  void BuildRatio(TH1D *ratio, TH1D *histN, TH1D *histD);

  ClassDef(AliTRDQADataMakerRec,1)   // Creates the TRD QA data

};
#endif // AliTRDQADataMakerRec_H
