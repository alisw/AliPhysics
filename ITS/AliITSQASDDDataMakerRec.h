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

#include "AliQA.h"
#include "AliITSQADataMakerRec.h"
class TObjArray;
class TH1D;
class TProfile2D;
class AliRawReader;
class AliITSgeomTGeo;
class AliITSDDLModuleMapSDD;

class AliITSQADataMakerRec;

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
	virtual void EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list);
  virtual ~AliITSQASDDDataMakerRec(); // dtor
  inline Int_t Raws() { return fSDDhRaws; }
  inline Int_t Recs() { return fSDDhRecs; }

private:

  static const Int_t fgknSDDmodules = 260; //number of SDD modules
  static const Int_t fgkmodoffset = 240;   //number of SPD modules
  static const Int_t fgknAnode = 256;      //anode per half-module
  static const Int_t fgknSide =2;          //side per module
  static const Int_t fgkeqOffset = 256;    //DDL offset
  static const Int_t fgkDDLidRange = 24;   //number of DDL:so DDL range is 257-280
  static const Int_t fgkDDLIDshift = 0;    //necessary option until RawStream Table is complete
  static const Int_t fgkLADDonLAY3 = 14;   //number of ladder on layer 3
  static const Int_t fgkLADDonLAY4 = 22;   //number of ladder on layer 4

  AliITSQADataMakerRec *fAliITSQADataMakerRec;//pointer to the main ctor
  Bool_t  fkOnline;                        //online (1) or offline (0) use
  Int_t   fLDC;                            //LDC number (0 for offline, 1 to 4 for online) 
  Int_t   fSDDhRaws;                       // number of histo booked for Raws SDD
  Int_t   fSDDhRecs;                       // number of histo booked for Recs SDD
  Int_t   fRawsOffset;                     // number of histo booked when SDD start
  Int_t   fRecsOffset;                     // number of histo booked when SDD start
  AliITSDDLModuleMapSDD  *fDDLModuleMap;// SDD Detector configuration for the decoding
/*
  TProfile2D *fModuleChargeMap[2*fgknSDDmodules];//module map
  TProfile2D *fModuleChargeMapFSE[2*fgknSDDmodules];//module map for one event 
*/ 
  ClassDef(AliITSQASDDDataMakerRec,3)      // description 

};

#endif


