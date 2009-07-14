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

#include "AliQAv1.h"
class AliITSQADataMakerSim;
class TObjArray;
class TClonesArray;

class AliITSQASPDDataMakerSim : public TObject {

public:
  AliITSQASPDDataMakerSim(AliITSQADataMakerSim *aliITSQADataMakerSim); //ctor
  AliITSQASPDDataMakerSim(const AliITSQASPDDataMakerSim& qadm);
  AliITSQASPDDataMakerSim& operator = (const AliITSQASPDDataMakerSim& qac);

  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray * list);
  virtual ~AliITSQASPDDataMakerSim() {;}   // dtor
  virtual Int_t InitDigits();
  virtual Int_t InitSDigits();
  virtual Int_t InitHits();
  virtual Int_t MakeDigits(){return 0;}
  virtual Int_t MakeSDigits(){return 0;}
  virtual Int_t MakeHits(){return 0;}
  virtual Int_t MakeDigits(TTree * digits);
  virtual Int_t MakeSDigits(TTree * sdigits);
  virtual Int_t MakeHits(TTree * hits);
  Int_t GetOffset(AliQAv1::TASKINDEX_t task);
  void  SetOffset(AliQAv1::TASKINDEX_t task, Int_t offset,Int_t specie = 0);
  Int_t GetTaskHisto(AliQAv1::TASKINDEX_t task);
  //Int_t GetOffsetH() { return fGenOffsetH; }
  //Int_t GetOffsetS() { return fGenOffsetS; }
  //Int_t GetOffsetD() { return fGenOffsetD; }
  //Int_t GetTaskHisto() { return fSPDhTask; }

private:

  AliITSQADataMakerSim *fAliITSQADataMakerSim;//pointer to the main ctor
  //  Int_t   fSPDhTask;                        //number of booked SPD histograms for each task;
  Int_t   fSPDhHTask;                          // number of booked SPD histograms for each task;
  Int_t   fSPDhSTask;                          // number of booked SPD histograms for each task;
  Int_t   fSPDhDTask;                          // number of booked SPD histograms for each task;
  Int_t   *fGenOffsetH;                         // qachecking offset       
  Int_t   *fGenOffsetS;                         // qachecking offset       
  Int_t   *fGenOffsetD;                         // qachecking offset       
  ClassDef(AliITSQASPDDataMakerSim,4)      // description 

};

#endif

