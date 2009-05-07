#ifndef ALITOFQADATAMAKER_H
#define ALITOFQADATAMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////
//                                                                // 
//  Produces the data needed to calculate the quality assurance.  //
//    All data must be mergeable objects.                         //
//    S. Arcelli                                                  //
//                                                                // 
////////////////////////////////////////////////////////////////////


#include "AliQADataMaker.h"
class AliTOFQADataMaker: public AliQADataMaker {

public:
  AliTOFQADataMaker() ;          // ctor
  AliTOFQADataMaker(const AliTOFQADataMaker& qadm) ;   
  AliTOFQADataMaker& operator = (const AliTOFQADataMaker& qadm) ;
  virtual ~AliTOFQADataMaker() {;} // dtor
  
private:
  virtual void   InitHits() ; 
  virtual void   InitESDs() ; 
  virtual void   InitDigits() ; 
  virtual void   InitRecPoints() ; 
  virtual void   InitRaws() ; 
  virtual void   InitSDigits() ; 
  virtual void   MakeESDs(AliESDEvent * esd) ;
  virtual void   MakeHits(TClonesArray * hits) ;
  virtual void   MakeHits(TTree * hitTree);
  virtual void   MakeDigits(TClonesArray * digits) ; 
  virtual void   MakeDigits(TTree * digTree);
  virtual void   MakeSDigits(TClonesArray * sdigits) ; 
  virtual void   MakeSDigits(TTree * sdigTree);
  virtual void   MakeRecPoints(TTree * recTree) ; 
  virtual void   MakeRaws(AliRawReader* rawReader) ; 
  virtual void   StartOfDetectorCycle() ; 
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX, TObjArray * list) ;
  virtual void   GetMapIndeces(Int_t *in, Int_t *out) ; 

  ClassDef(AliTOFQADataMaker,1)  // description 

};

#endif // AliTOFQADataMaker_H
