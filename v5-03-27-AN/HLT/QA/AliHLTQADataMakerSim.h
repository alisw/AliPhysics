//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTQADATAMAKERSIM_H
#define ALIHLTQADATAMAKERSIM_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTQADataMakerSim.h
    @author Matthias Richter
    @date   2009-05-14
    @brief  Container for the HLT offline QA
*/

#include "AliQADataMakerSim.h"

class AliHLTQADataMakerSim: public AliQADataMakerSim {

public:
  
  AliHLTQADataMakerSim();
  virtual ~AliHLTQADataMakerSim();
  
private:
  /** copy constructor prohibited */
  AliHLTQADataMakerSim(const AliHLTQADataMakerSim&);   
  /** assignment operator prohibited */
  AliHLTQADataMakerSim& operator = (const AliHLTQADataMakerSim&);

  virtual void   StartOfDetectorCycle();
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray** list) ;

  // Digits QA
  virtual void   InitDigits();
  virtual void   MakeDigits(TTree *digitTree);
  virtual void   MakeDigits(TClonesArray *);

  // Hits QA
  virtual void   InitHits();
  virtual void   MakeHits(TTree *hitTree);
  virtual void   MakeHits(TClonesArray *);

  // SDigits QA (empty)
  virtual void   InitSDigits();
  virtual void   MakeSDigits(TTree* );
  virtual void   MakeSDigits(TClonesArray* );

  ClassDef(AliHLTQADataMakerSim,0)  // HLT Quality Assurance Data Maker for simulation
};

#endif // ALIHLTQADATAMAKERSIM_H
