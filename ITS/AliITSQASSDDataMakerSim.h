#ifndef AliITSQASSDDataMakerSim_H
#define AliITSQASSDDataMakerSim_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
//  Checks the quality assurance. 
//  By comparing with reference data
//  contained in a DB
//
//
//  W. Ferrarese + P. Cerello Feb 2008
//  SSD QA part: P. Christakoglou

#include "AliQA.h"
class AliITSQADataMakerSim;
class TObjArray;
class TClonesArray;

class AliITSQASSDDataMakerSim : public TObject {

public:
  AliITSQASSDDataMakerSim(AliITSQADataMakerSim *aliITSQADataMakerSim); //ctor
  AliITSQASSDDataMakerSim(const AliITSQASSDDataMakerSim& qadm);
  AliITSQASSDDataMakerSim& operator = (const AliITSQASSDDataMakerSim& qac);

  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list);
  virtual ~AliITSQASSDDataMakerSim() {;}   // dtor
  virtual void InitDigits();
  virtual void InitSDigits();
  virtual void InitHits();
  virtual void MakeDigits(TClonesArray * /*digits*/){;}
  virtual void MakeSDigits(TClonesArray * /*sdigits*/){;}
  virtual void MakeHits (TClonesArray * /*hits*/){;}
  virtual void MakeDigits(TTree * digits);
  virtual void MakeSDigits(TTree * sdigits);
  virtual void MakeHits(TTree * hits);
  Int_t GetOffset() { return fGenOffset; }
  Int_t GetTaskHisto() { return fSSDhTask; }

private:

  AliITSQADataMakerSim *fAliITSQADataMakerSim; //pointer to the main ctor
  Int_t   fSSDhTask;    //number of booked SSD histograms for each task;
  Int_t   fGenOffset;                         // qachecking offset

  static const Int_t fgkNumberOfPSideStrips = 768; //number of P-side strips

  ClassDef(AliITSQASSDDataMakerSim,1)      // description 
};

#endif


