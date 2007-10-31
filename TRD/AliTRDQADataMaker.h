#ifndef AliTRDQADatamaker_H
#define AliTRDQADatamaker_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
Produces the data needed to calculate the quality assurance. 
All data must be mergeable objects.
S.Radomski Uni-Heidelberg October 2007
*/

// --- ROOT system ---
class TH1F ; 
class TH1I ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQADataMaker.h"

class AliTRDQADataMaker: public AliQADataMaker {

public:
  AliTRDQADataMaker() ;          // ctor
  AliTRDQADataMaker(const AliTRDQADataMaker& qadm) ;   
  AliTRDQADataMaker& operator = (const AliTRDQADataMaker& qadm) ;
  virtual ~AliTRDQADataMaker() {;} // dtor
  
private:
  virtual void EndOfDetectorCycle() ;
  virtual void InitHits() ; 
  virtual void InitESDs() ; 
  virtual void InitDigits() ; 
  virtual void InitRecPoints() ; 
  virtual void InitRaws() ; 
  virtual void InitSDigits() ;
 
  virtual void MakeHits(TTree * hitTree);
  virtual void MakeHits(TClonesArray * hits);

  //virtual void MakeSDigits(TTree *sdigitTree);
  virtual void MakeSDigits(TClonesArray * sigits); 

  //virtual void MakeDigits(TTree *digitTree);
  virtual void MakeDigits(TClonesArray * digits); 

  virtual void MakeRaws(AliRawReader* rawReader); 
  virtual void MakeRecPoints(TTree * recpo); 
  virtual void MakeESDs(AliESDEvent * esd);
  
  virtual void StartOfDetectorCycle() ; 
  Int_t CheckPointer(TObject *obj, const char *name);

  ClassDef(AliTRDQADataMaker,1)  // description 

};

#endif // AliTRDQADatamaker_H
