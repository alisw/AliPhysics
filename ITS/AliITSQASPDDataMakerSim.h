#ifndef AliITSQASPDDataMakerSim_H
#define AliITSQASPDDataMakerSim_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//
//  Checks the quality assurance. 
//  By comparing with reference data
//  contained in a DB
//
//
//  W. Ferrarese + P. Cerello Feb 2008
//


/* $Id$ */

#include "AliQA.h"
class AliITSQADataMakerSim;
class TObjArray;
class TClonesArray;

class AliITSQASPDDataMakerSim : public TObject {

public:
  AliITSQASPDDataMakerSim(AliITSQADataMakerSim *aliITSQADataMakerSim); //ctor
  AliITSQASPDDataMakerSim(const AliITSQASPDDataMakerSim& qadm);
  AliITSQASPDDataMakerSim& operator = (const AliITSQASPDDataMakerSim& qac);

  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list);
  virtual ~AliITSQASPDDataMakerSim() {;}   // dtor
  virtual void InitDigits();
  virtual void InitSDigits();
  virtual void InitHits();
  virtual void MakeDigits(TClonesArray * /*digits*/){;}
  virtual void MakeSDigits(TClonesArray * /*sdigits*/){;}
  virtual void MakeHits(TClonesArray * /*hits*/){;}
  virtual void MakeDigits(TTree * digits);
  virtual void MakeSDigits(TTree * sdigits);
  virtual void MakeHits(TTree * hits);
  Int_t GetOffset() { return fGenOffset; }
  Int_t GetTaskHisto() { return fSPDhTask; }

private:

  AliITSQADataMakerSim *fAliITSQADataMakerSim;//pointer to the main ctor
  Int_t   fSPDhTask;                        //number of booked SPD histograms for each task;
  Int_t   fGenOffset;                         // qachecking offset       
  ClassDef(AliITSQASPDDataMakerSim,2)      // description 

};

#endif

