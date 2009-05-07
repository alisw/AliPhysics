#ifndef AliITSQASDDDataMakerSim_H
#define AliITSQASDDDataMakerSim_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//
//  Checks the quality assurance. 
//  By comparing with reference data
//  contained in a DB
//
//
//  W. Ferrarese + P. Cerello Feb 2008

/* $Id$ */

#include "AliQAv1.h"
class AliITSQADataMakerSim;
class AliRunLoader;
class AliRun;
class TObjArray;
class TClonesArray;

class AliITSQASDDDataMakerSim : public TObject {

public:
  AliITSQASDDDataMakerSim(AliITSQADataMakerSim *aliITSQADataMakerSim); //ctor
  AliITSQASDDDataMakerSim(const AliITSQASDDDataMakerSim& qadm);
  AliITSQASDDDataMakerSim& operator = (const AliITSQASDDDataMakerSim& qac);

  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray * list);
  virtual ~AliITSQASDDDataMakerSim() {;}   // dtor
  virtual void InitDigits();
  virtual void InitSDigits();
  virtual void InitHits();
  virtual void MakeDigits(TClonesArray * /*digits*/){;}
  virtual void MakeSDigits(TClonesArray * /*sdigits*/){;}
  virtual void MakeHits(TClonesArray * /*hits*/){;}
  virtual void MakeDigits(TTree * digits);
  virtual void MakeSDigits(TTree * sdigits);
  virtual void MakeHits(TTree * hits);
  Int_t GetOffset(AliQAv1::TASKINDEX_t task);
  Int_t GetTaskHisto(AliQAv1::TASKINDEX_t task);

private:

  AliITSQADataMakerSim *fAliITSQADataMakerSim;//pointer to the main ctor
  Int_t   fSDDhHTask;                        //number of booked SDD histograms for each task;
  Int_t   fSDDhSTask;                        //number of booked SDD histograms for each task;
  Int_t   fSDDhDTask;                        //number of booked SDD histograms for each task;
  Int_t   fGenOffsetH;                         // qachecking offset
  Int_t   fGenOffsetS;                         // qachecking offset
  Int_t   fGenOffsetD;                         // qachecking offset
  ClassDef(AliITSQASDDDataMakerSim,3)      // description 

};

#endif


