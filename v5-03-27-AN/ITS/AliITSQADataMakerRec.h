#ifndef ALIITSQADATAMAKERREC_H
#define ALIITSQADATAMAKERREC_H
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

class AliDetectorRecoParam;
class AliReconstructor;
//#include "AliITSDDLModuleMapSDD.h"

class AliQAManager;
class AliITSQASPDDataMakerRec;
class AliITSQASDDDataMakerRec;
class AliITSQASSDDataMakerRec;
class AliITSRecPoint;
class AliRawReader;
class TH2F;
class AliITSDDLModuleMapSDD;

class AliITSQADataMakerRec: public AliQADataMakerRec {

  friend class AliITSQASPDDataMakerRec; //friend class of SPD QA
  friend class AliITSQASDDDataMakerRec; //friend class of SDD QA
  friend class AliITSQASSDDataMakerRec; //friend class of SSD QA

public:
  AliITSQADataMakerRec(Bool_t kMode = kFALSE, Short_t subDet = 0, Short_t ldc = 0); // kMode = kFALSE (offline), kTRUE (online); subDet = 0 (ALL), 1 (SPD), 2 (SDD), 3 (SSD)
  AliITSQADataMakerRec(const AliITSQADataMakerRec& qadm);
  AliITSQADataMakerRec& operator = (const AliITSQADataMakerRec& qac);
  virtual Int_t GetEventSpecie() const { return AliRecoParam::AConvert(fEventSpecie); }
  virtual void StartOfDetectorCycle();
  virtual void StartOfCycle(AliQAv1::TASKINDEX_t task, Int_t run, const Bool_t sameCycle = kFALSE) ;
  virtual void StartOfCycle(Int_t run){AliQADataMakerRec::StartOfCycle(run);}
  virtual void EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list);
  //  virtual void EndOfDetectorCycle(const char *fgDataName);
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
  virtual Bool_t ListExists(AliQAv1::TASKINDEX_t task) const;
  AliQAv1::TASKINDEX_t GetTaskIndexSelected() const {return fSelectedTaskIndex;}

  virtual void ResetDetector(AliQAv1::TASKINDEX_t task);

  virtual ~AliITSQADataMakerRec(); // dtor
  Short_t GetSubDet()const {return fSubDetector;};
  Int_t GetDetTaskOffset(Int_t subdet,AliQAv1::TASKINDEX_t task,Int_t specie=0);
  Int_t GetDetTaskHisto(Int_t subdet,AliQAv1::TASKINDEX_t task);
  TH2F *GetITSGlobalHisto(Int_t layer);
  static Bool_t AreEqual(Double_t a1, Double_t a2);
  
  virtual void SetRunNumber(Int_t runnumber){fRunNumber=runnumber;};
  Int_t GetRunNumber()const {return fRunNumber;};
  
  virtual void SetEventNumber(Int_t eventnumber){fEventNumber=eventnumber;};
  Int_t GetEventNumber() const {return fEventNumber;};
  AliITSDDLModuleMapSDD *GetDDLSDDModuleMap();
  
 private:
  
  Bool_t  fkOnline;                        //online (1) or offline (0) use
  Short_t fSubDetector;                    // subDetector: 0 (ALL), 1 (SPD), 2 (SDD), 3 (SSD)
  Short_t fLDC;				   // number of LDC: 0 (one LDC for the whole subdetector)
  Int_t fRunNumber;                        //run number
  Int_t fEventNumber;                      //Event number (online mode)
  AliQAv1::TASKINDEX_t fSelectedTaskIndex; //Current TaskIndex

  AliITSQASPDDataMakerRec *fSPDDataMaker;  // SPD Data Maker 
  AliITSQASDDDataMakerRec *fSDDDataMaker;  // SDD Data Maker 
  AliITSQASSDDataMakerRec *fSSDDataMaker;  // SSD Data Maker 

  ClassDef(AliITSQADataMakerRec,8)         // description 

};

#endif

