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
#include "AliDetectorRecoParam.h"
#include "AliReconstructor.h"

class AliITSQASPDDataMakerRec;
class AliITSQASDDDataMakerRec;
class AliITSQASSDDataMakerRec;
class AliITSRecPoint;
class AliRawReader;
class TH2F;

class AliITSQADataMakerRec: public AliQADataMakerRec {

friend class AliITSQASPDDataMakerRec;
friend class AliITSQASDDDataMakerRec;
friend class AliITSQASSDDataMakerRec;

public:
  AliITSQADataMakerRec(Bool_t kMode = kFALSE, Short_t subDet = 0, Short_t ldc = 0); // kMode = kFALSE (offline), kTRUE (online); subDet = 0 (ALL), 1 (SPD), 2 (SDD), 3 (SSD)
  AliITSQADataMakerRec(const AliITSQADataMakerRec& qadm);
  AliITSQADataMakerRec& operator = (const AliITSQADataMakerRec& qac);
  virtual Int_t GetEventSpecie() const { return AliRecoParam::AConvert(fEventSpecie); }
  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list);
  virtual void EndOfDetectorCycle(const char *fgDataName);
  virtual void InitRaws();
  virtual void InitDigits();
  virtual void InitRecPoints();
  virtual void InitESDs();
  virtual void MakeRaws(AliRawReader *rawReader);
  virtual void MakeDigits(){AliWarning("Signature not implemented. A TTree* of digits should be passed as input argument");} 
  virtual void MakeDigits(TTree *digitsTree);
  virtual void MakeRecPoints(TTree *clustersTree);
  virtual void MakeESDs(AliESDEvent *esd);
  virtual void FillRecPoint(AliITSRecPoint rcp);

  virtual ~AliITSQADataMakerRec(); // dtor
 Short_t GetSubDet(){return fSubDetector;};
 Int_t GetDetTaskOffset(Int_t subdet,AliQAv1::TASKINDEX_t task);
 TH2F *GetITSGlobalHisto(Int_t layer);
 static Bool_t IsEqual(Double_t a1, Double_t a2);


private:

  Bool_t  fkOnline;                        //online (1) or offline (0) use
  Short_t fSubDetector;                    // subDetector: 0 (ALL), 1 (SPD), 2 (SDD), 3 (SSD)
  Short_t fLDC;				   // number of LDC: 0 (one LDC for the whole subdetector)

  AliITSQASPDDataMakerRec *fSPDDataMaker;  // SPD Data Maker 
  AliITSQASDDDataMakerRec *fSDDDataMaker;  // SDD Data Maker 
  AliITSQASSDDataMakerRec *fSSDDataMaker;  // SSD Data Maker 

  ClassDef(AliITSQADataMakerRec,5)         // description 

};

#endif

