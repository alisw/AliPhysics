#ifndef AliZDCQADataMakerRec_H
#define AliZDCQADataMakerRec_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////
//  						    //
//  Produces the data needed for quality assurance. //
//  C. Oppedisano Chiara.Oppedisano@to.infn.it      //
//  						    //
//////////////////////////////////////////////////////


#include "AliQADataMakerRec.h"
class AliZDCQADataMakerRec: public AliQADataMakerRec {

public:
  AliZDCQADataMakerRec() ;          // ctor
  AliZDCQADataMakerRec(const AliZDCQADataMakerRec& qadm) ;   
  AliZDCQADataMakerRec& operator = (const AliZDCQADataMakerRec& qadm) ;
  virtual ~AliZDCQADataMakerRec() {;} // dtor
  
private:
  virtual void   InitESDs(); 
  virtual void   InitDigits() ; 
  virtual void   InitRecPoints() {;} 
  virtual void   InitRaws(); 
  virtual void   MakeDigits(TTree * /*recTree*/) ;
  virtual void   MakeDigits(TClonesArray* /*digits*/)  {return;}
  virtual void   MakeRecPoints(TTree * /*recTree*/) {;} 
  virtual void   MakeRaws(AliRawReader* rawReader) ; 
  virtual void   MakeESDs(AliESDEvent * esd) ;
  virtual void   StartOfDetectorCycle() ; 
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list) ;

  ClassDef(AliZDCQADataMakerRec,1)  // description 

};

#endif // AliZDCQADataMakerRec_H
