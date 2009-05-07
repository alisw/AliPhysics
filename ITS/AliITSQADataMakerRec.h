#ifndef AliITSQADataMakerRec_H
#define AliITSQADataMakerRec_H
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
//
//  ESD QA (Tracking and primary vertex)
//  A. Dainese Jun 2008 

#include "AliQADataMakerRec.h"

class AliITSQASPDDataMakerRec;
class AliITSQASDDDataMakerRec;
class AliITSQASSDDataMakerRec;
class AliRawReader;

class AliITSQADataMakerRec: public AliQADataMakerRec {

friend class AliITSQASPDDataMakerRec;
friend class AliITSQASDDDataMakerRec;
friend class AliITSQASSDDataMakerRec;

public:
  AliITSQADataMakerRec(Bool_t kMode = kFALSE, Short_t subDet = 0, Short_t ldc = 0); // kMode = kFALSE (offline), kTRUE (online); subDet = 0 (ALL), 1 (SPD), 2 (SDD), 3 (SSD)
  AliITSQADataMakerRec(const AliITSQADataMakerRec& qadm);
  AliITSQADataMakerRec& operator = (const AliITSQADataMakerRec& qac);
  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list);
  virtual void EndOfDetectorCycle(const char *fgDataName);
  virtual void InitRaws();
  virtual void InitRecPoints();
  virtual void InitESDs();
  virtual void MakeRaws(AliRawReader *rawReader);
  virtual void MakeRecPoints(TTree *clustersTree);
  virtual void MakeESDs(AliESDEvent *esd);

  void SetHLTMode(Bool_t khltmode=kFALSE){fHLTMode=khltmode;};
  Bool_t GetHLTMode(){return fHLTMode;};
  virtual ~AliITSQADataMakerRec(); // dtor
 Short_t GetSubDet(){return fSubDetector;};
 Int_t GetDetTaskOffset(Int_t subdet,AliQAv1::TASKINDEX_t task);


private:

  Bool_t  fkOnline;                        //online (1) or offline (0) use
  Bool_t  fHLTMode;                        // HLT MODE kTRUE mode C kFALSE mode A
  Short_t fSubDetector;                    // subDetector: 0 (ALL), 1 (SPD), 2 (SDD), 3 (SSD)
  Short_t fLDC;				   // number of LDC: 0 (one LDC for the whole subdetector)

  AliITSQASPDDataMakerRec *fSPDDataMaker;  // SPD Data Maker 
  AliITSQASDDDataMakerRec *fSDDDataMaker;  // SDD Data Maker 
  AliITSQASSDDataMakerRec *fSSDDataMaker;  // SSD Data Maker 

  ClassDef(AliITSQADataMakerRec,4)         // description 

};

#endif

