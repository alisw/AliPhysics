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
  virtual Int_t InitDigits();
  virtual Int_t InitSDigits();
  virtual Int_t InitHits();
  virtual Int_t MakeDigits(){return 0;}
  virtual Int_t MakeSDigits(){return 0;}
  virtual Int_t MakeHits (){return 0;}
  virtual Int_t MakeDigits(TTree * digits);
  virtual Int_t MakeSDigits(TTree * sdigits);
  virtual Int_t MakeHits(TTree * hits);
  Int_t GetOffset(AliQAv1::TASKINDEX_t task);
  void  SetOffset(AliQAv1::TASKINDEX_t task, Int_t offset,Int_t specie = 0);
  Int_t GetTaskHisto(AliQAv1::TASKINDEX_t task);

private:

  AliITSQADataMakerSim *fAliITSQADataMakerSim; //pointer to the main ctor

  Int_t   fSSDhHTask;   //number of booked SSD histograms for the hit task;
  Int_t   fSSDhSTask;   //number of booked SSD histograms for the sdigits task;
  Int_t   fSSDhDTask;   //number of booked SSD histograms for the digit task;
  Int_t   *fGenOffsetH;                         // qachecking offset hits
  Int_t   *fGenOffsetS;                         // qachecking offset sdigits
  Int_t   *fGenOffsetD;                         // qachecking offset digits

  static const Int_t fgkNumberOfPSideStrips = 768; //number of P-side strips

  ClassDef(AliITSQASSDDataMakerSim,3)      // description 
};

#endif


