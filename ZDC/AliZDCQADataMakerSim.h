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
  virtual void   MakeHits();
  virtual void   MakeHits(TTree * hitTree);
  virtual void   MakeDigits() {;} 
  virtual void   MakeDigits(TTree * digTree);
  virtual void   MakeSDigits() {;} 
  virtual void   MakeSDigits(TTree * /*sdigTree*/) {;}
  virtual void   StartOfDetectorCycle(); 
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list);
  
  ClassDef(AliZDCQADataMakerSim, 3)  // description 

};

#endif // AliZDCQADataMakerSim_H
