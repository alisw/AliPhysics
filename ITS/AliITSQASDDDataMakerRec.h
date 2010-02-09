#ifndef AliITSQASDDDataMakerRec_H
#define AliITSQASDDDataMakerRec_H
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
#include "AliITSQADataMakerRec.h"

class TObjArray;
class AliITSDDLModuleMapSDD;
class AliITSQASDDDataMakerRec: public TObject {

public:
  AliITSQASDDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode = kFALSE, Short_t ldc = 0);
  AliITSQASDDDataMakerRec(const AliITSQASDDDataMakerRec& qadm);
  AliITSQASDDDataMakerRec& operator = (const AliITSQASDDDataMakerRec& qac);
  virtual Int_t InitRaws();
  virtual Int_t InitDigits();
  virtual Int_t InitRecPoints();
  virtual Int_t MakeRaws(AliRawReader *rawReader);
  virtual Int_t MakeDigits()  {return 0;}
  virtual Int_t MakeDigits(TTree *clustersTree);
  virtual Int_t MakeRecPoints(TTree *clustersTree);
  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray * list);
  virtual void CreateTheMap();

  virtual ~AliITSQASDDDataMakerRec(); // dtor
  Int_t GetOffset(AliQAv1::TASKINDEX_t task,Int_t specie=0);
  void  SetOffset(AliQAv1::TASKINDEX_t task, Int_t offset, Int_t specie = 0);
  Int_t GetTaskHisto(AliQAv1::TASKINDEX_t task);
  virtual void ResetDetector(AliQAv1::TASKINDEX_t task);
  AliITSDDLModuleMapSDD* GetDDLSDDModuleMap(){return fDDLModuleMap; };

private:

  static const Int_t fgknSDDmodules = 260; // number of SDD modules
  static const Int_t fgkmodoffset = 240;   // number of SPD modules
  static const Int_t fgknAnode = 256;      // anode per half-module
  static const Int_t fgknSide =2;          // side per module
  //static const Int_t fgkDDLIDshift = 0;    // necessary option until RawStream Table is complete
  static const Int_t fgkLADDonLAY3 = 14;   // number of ladder on layer 3
  static const Int_t fgkLADDonLAY4 = 22;   // number of ladder on layer 4

  AliITSQADataMakerRec *fAliITSQADataMakerRec;// pointer to the main ctor
  Bool_t  fkOnline;                           // online (1) or offline (0) use
  Int_t   fLDC;                               // LDC number (0 for offline, 1 to 4 for online) 
  Int_t   fSDDhRawsTask;                      // number of histo booked for each the Raws Task SDD
  Int_t   fSDDhDigitsTask;                    // number of histo booked for each the RecPoints Task SDD
  Int_t   fSDDhRecPointsTask;                 // number of histo booked for each the RecPoints Task SDD
  Int_t   *fGenRawsOffset;                     // QAchecking Raws offset       
  Int_t   *fGenDigitsOffset;                   // QAchecking RecPoints offset       
  Int_t   *fGenRecPointsOffset;                // QAchecking RecPoints offset       
  Int_t   fTimeBinSize;			      // time bin width in number of clocks
  AliITSDDLModuleMapSDD  *fDDLModuleMap;      // SDD Detector configuration for the decoding


  ClassDef(AliITSQASDDDataMakerRec,11)         // description 

};

#endif

