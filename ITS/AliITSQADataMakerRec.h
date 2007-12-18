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
//  W. Ferrarese Nov 2007




#include "AliQADataMakerRec.h"
class TObjArray;
class TH1F;
class TH2D;
class AliRawReader;
class AliITSgeomTGeo;

class AliITSQADataMakerRec: public AliQADataMakerRec {

public:
  AliITSQADataMakerRec();          // ctor
  AliITSQADataMakerRec(Int_t ldc, Bool_t kMode = kFALSE);
  AliITSQADataMakerRec(const AliITSQADataMakerRec& qadm);
  AliITSQADataMakerRec& operator = (const AliITSQADataMakerRec& qac);
  virtual void StartOfDetectorCycle() const;
  virtual void EndOfDetectorCycle(AliQA::TASKINDEX task, TObjArray * list);
  virtual void EndOfDetectorCycle(const char * fgDataName);
  virtual void InitRaws();
  virtual void InitRecPoints();
  virtual void MakeRaws(AliRawReader *rawReader);
  virtual void MakeRecPoints(TTree *clustersTree);
  virtual ~AliITSQADataMakerRec() {;} // dtor

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
	
  Bool_t  fkOnline;                         //online (1) or offline (0) use
  Int_t   fLDC;                            //LDC number (0 for offline, 1 to 4 for online) 
  Int_t   fnSDDHistos;                     // number of histogrma booked for SDDs
  TH2D *fModuleChargeMap[2*fgknSDDmodules];//module map
  TH1D *fmonoD[2*fgknSDDmodules] ;         //histo used as support
 
  ClassDef(AliITSQADataMakerRec,1)  // description 

};

#endif


