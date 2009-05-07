#ifndef ALIZDCQADATAMAKER_H
#define ALIZDCQADATAMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////
//  						    //
//  Produces the data needed for quality assurance. //
//  C. Oppedisano Chiara.Oppedisano@to.infn.it      //
//  						    //
//////////////////////////////////////////////////////


#include "AliQADataMaker.h"

class AliZDCQADataMaker: public AliQADataMaker {

public:
  AliZDCQADataMaker();          // ctor
  AliZDCQADataMaker(const AliZDCQADataMaker& qadm);   
  AliZDCQADataMaker& operator = (const AliZDCQADataMaker& qadm);
  virtual ~AliZDCQADataMaker() {;} // dtor

private:

  virtual void  InitHits();	 
  virtual void  InitDigits();	 
  virtual void  InitSDigits() ;   
  virtual void  InitRecPoints();  
  virtual void  InitRaws();	  
  virtual void  InitESDs();	 
  virtual void  MakeHits(TClonesArray * hits = 0);       
  virtual void  MakeHits(TTree * hits);       
  virtual void  MakeDigits(TClonesArray * digits  = 0);   
  virtual void  MakeDigits(TTree * digits);   
  virtual void  MakeSDigits(TClonesArray * /*sdigits*/) {;} 
  virtual void  MakeSDigits(TTree * /*sdigTree*/) {;}
  virtual void  MakeRecPoints(const TTree * const clusters) {;}  
  virtual void  MakeRaws(AliRawReader* rawReader);
  virtual void  MakeESDs(AliESDEvent * esd);	     
  virtual void  StartOfDetectorCycle();
  virtual void  EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray * list);
  //
  TClonesArray * fHits; 	//! Array containing ZDC hits
  TClonesArray * fDigits; 	//! Array containing ZDC digits

  ClassDef(AliZDCQADataMaker,2)  // description 

};

#endif // AliZDCQADataMaker_H
