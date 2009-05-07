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
//  M.Siciliano Aug 2008 QA RecPoints and HLT mode

/* $Id$ */

#include "AliQAv1.h"
#include "AliITSQADataMakerRec.h"

class TObjArray;
class AliITSDDLModuleMapSDD;
class AliITSHLTforSDD;
class AliITSQASDDDataMakerRec: public TObject {

public:
  AliITSQASDDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode = kFALSE, Short_t ldc = 0);
  AliITSQASDDDataMakerRec(const AliITSQASDDDataMakerRec& qadm);
  AliITSQASDDDataMakerRec& operator = (const AliITSQASDDDataMakerRec& qac);
  virtual void InitRaws();
  virtual void InitRecPoints();
  virtual void MakeRaws(AliRawReader *rawReader);
  virtual void MakeRecPoints(TTree *clustersTree);
  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray * list);


  virtual ~AliITSQASDDDataMakerRec(); // dtor
  Int_t GetOffset(AliQAv1::TASKINDEX_t task);
  Int_t GetTaskHisto(AliQAv1::TASKINDEX_t task);

  void SetHLTMode(Bool_t khltmode=kFALSE){fHLTMode=khltmode;};
  Bool_t GetHLTMode(){return fHLTMode;};
  void SetHLTModeFromEnvironment();

private:

  static const Int_t fgknSDDmodules = 260; // number of SDD modules
  static const Int_t fgkmodoffset = 240;   // number of SPD modules
  static const Int_t fgknAnode = 256;      // anode per half-module
  static const Int_t fgknSide =2;          // side per module
  static const Int_t fgkDDLIDshift = 0;    // necessary option until RawStream Table is complete
  static const Int_t fgkLADDonLAY3 = 14;   // number of ladder on layer 3
  static const Int_t fgkLADDonLAY4 = 22;   // number of ladder on layer 4

  AliITSQADataMakerRec *fAliITSQADataMakerRec;// pointer to the main ctor
  Bool_t  fkOnline;                           // online (1) or offline (0) use
  Int_t   fLDC;                               // LDC number (0 for offline, 1 to 4 for online) 
  Int_t   fSDDhRawsTask;                      // number of histo booked for each the Raws Task SDD
  Int_t   fSDDhRecPointsTask;                 // number of histo booked for each the RecPoints Task SDD
  Int_t   fGenRawsOffset;                     // QAchecking Raws offset       
  Int_t   fGenRecPointsOffset;                // QAchecking RecPoints offset       
  Int_t   fTimeBinSize;			      // time bin width in number of clocks
  AliITSDDLModuleMapSDD  *fDDLModuleMap;      // SDD Detector configuration for the decoding
  Bool_t fHLTMode;                            // kTRUE mode C kFALSE mode A 
                                              // Used in online mode only
  AliITSHLTforSDD *fHLTSDD;                   // used for offline QA as the HLT mode flag
  ClassDef(AliITSQASDDDataMakerRec,7)         // description 

};

#endif

