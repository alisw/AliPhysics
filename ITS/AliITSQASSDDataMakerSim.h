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

#include "AliQAv1.h"
class AliITSQADataMakerSim;
class TObjArray;
class TClonesArray;

class AliITSQASSDDataMakerSim : public TObject {

public:
  AliITSQASSDDataMakerSim(AliITSQADataMakerSim *aliITSQADataMakerSim); //ctor
  AliITSQASSDDataMakerSim(const AliITSQASSDDataMakerSim& qadm);
  AliITSQASSDDataMakerSim& operator = (const AliITSQASSDDataMakerSim& qac);

  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray * list);
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
  Int_t GetOffset(AliQAv1::TASKINDEX_t task);
  Int_t GetTaskHisto(AliQAv1::TASKINDEX_t task);

private:

  AliITSQADataMakerSim *fAliITSQADataMakerSim; //pointer to the main ctor

  Int_t   fSSDhHTask;   //number of booked SSD histograms for the hit task;
  Int_t   fSSDhSTask;   //number of booked SSD histograms for the sdigits task;
  Int_t   fSSDhDTask;   //number of booked SSD histograms for the digit task;
  Int_t   fGenOffsetH;                         // qachecking offset hits
  Int_t   fGenOffsetS;                         // qachecking offset sdigits
  Int_t   fGenOffsetD;                         // qachecking offset digits

  static const Int_t fgkNumberOfPSideStrips = 768; //number of P-side strips

  ClassDef(AliITSQASSDDataMakerSim,2)      // description 
};

#endif


