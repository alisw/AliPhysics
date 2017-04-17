// -*- C++ -*-
#ifndef ALIADQADATAMAKERSIM_H
#define ALIADQADATAMAKERSIM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved.
 *
 * See cxx source for full Copyright notice
 */

#include "AliQADataMakerSim.h"

class TH1F ;
class TH1I ;
class TList ;

//_____________________________________________________________________
//
// This class implements the AliQADataMakerSim for the AD. Some
// functions are not implemented yet.
// Author : BC
//_____________________________________________________________________

class AliADQADataMakerSim: public AliQADataMakerSim {
public:
  AliADQADataMakerSim() ;          // ctor
  AliADQADataMakerSim(const AliADQADataMakerSim& qadm) ;
  AliADQADataMakerSim& operator = (const AliADQADataMakerSim& qadm) ;
  virtual ~AliADQADataMakerSim() {} // dtor

private:
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** list);
  virtual void   InitHits();
  virtual void   InitSDigits();
  virtual void   InitDigits();
  virtual void   MakeHits() ;
  virtual void   MakeHits(TTree* hitTree) ;
  virtual void   MakeSDigits() ;
  virtual void   MakeSDigits(TTree* digitTree) ;
  virtual void   MakeDigits() ;
  virtual void   MakeDigits(TTree* sdigitTree) ;
  virtual void   StartOfDetectorCycle() ;

  ClassDef(AliADQADataMakerSim,0);
};

#endif // AliADQADataMakerSim_H
//____________________________________________________________________
