#ifndef AliZDCQADataMakerSim_H
#define AliZDCQADataMakerSim_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////
//  						    //
//  Produces the data needed for quality assurance. //
//  C. Oppedisano Chiara.Oppedisano@to.infn.it      //
//  						    //
//////////////////////////////////////////////////////

class AliZDCDigit;

#include "AliQADataMakerSim.h"
class AliZDCQADataMakerSim: public AliQADataMakerSim {

public:
  AliZDCQADataMakerSim();          // ctor
  AliZDCQADataMakerSim(const AliZDCQADataMakerSim& qadm);   
  AliZDCQADataMakerSim& operator = (const AliZDCQADataMakerSim& qadm) ;
  virtual ~AliZDCQADataMakerSim() {;} // dtor
  
private:
  virtual void   InitHits(); 
  virtual void   InitDigits(); 
  virtual void   InitSDigits() {;} 
  virtual void   MakeHits(TClonesArray * hits = 0);
  virtual void   MakeHits(TTree * hitTree);
  virtual void   MakeDigits(TClonesArray * /*digits*/){;} 
  virtual void   MakeDigits(TTree * digTree);
  virtual void   MakeSDigits(TClonesArray * /*sdigits*/) {;} 
  virtual void   MakeSDigits(TTree * /*sdigTree*/) {;}
  virtual void   StartOfDetectorCycle(); 
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray ** list);
  //
  TClonesArray   *fHits; 	//! Array containing ZDC hits
  AliZDCDigit    *fDigit; 	//! Pointer to digit in tree  
  
  ClassDef(AliZDCQADataMakerSim,2)  // description 

};

#endif // AliZDCQADataMakerSim_H
