#ifndef AliITSQADataMakerSim_H
#define AliITSQADataMakerSim_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Checks the quality assurance. 
//  By comparing with reference data
//  contained in a DB
//
//
//  W. Ferrarese Nov 2007




#include "AliQADataMakerSim.h"
class TObjArray;
class TH1F;
class TH2D;
class AliRawReader;
class AliITSgeomTGeo;

class AliITSQADataMakerSim: public AliQADataMakerSim {

public:
  AliITSQADataMakerSim();          // ctor
  AliITSQADataMakerSim(Int_t /* ldc */, Bool_t /* kMode =  kFALSE */);
  AliITSQADataMakerSim(const AliITSQADataMakerSim& qadm);
  AliITSQADataMakerSim& operator = (const AliITSQADataMakerSim& qac);
  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list);
   virtual ~AliITSQADataMakerSim() {;} // dtor

private:

  ClassDef(AliITSQADataMakerSim,1)  // description 

};

#endif


