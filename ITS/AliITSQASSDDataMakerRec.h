#ifndef AliITSQASSDDataMakerRec_H
#define AliITSQASSDDataMakerRec_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
//  Checks the quality assurance. 
//  By comparing with reference data
//  contained in a DB
//  -------------------------------------------------------------
//  W. Ferrarese + P. Cerello Feb 2008
//  INFN Torino

#include "AliQA.h"
#include "AliITSQADataMakerRec.h"

class TObjArray;
class TH1D;
class AliRawReader;

class AliITSQADataMakerRec;

class AliITSQASSDDataMakerRec: public TObject {

public:
  AliITSQASSDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode = kFALSE, Int_t ldc=0);  //ctor
  AliITSQASSDDataMakerRec(const AliITSQASSDDataMakerRec& qadm);
  AliITSQASSDDataMakerRec& operator = (const AliITSQASSDDataMakerRec& qac);
  virtual void InitRaws();
  virtual void InitRecPoints();
  virtual void MakeRaws(AliRawReader *rawReader);
  virtual void MakeRecPoints(TTree *clustersTree);
  virtual void StartOfDetectorCycle();
  virtual void EndOfDetectorCycle(AliQA::TASKINDEX task, TObjArray * list);
  virtual ~AliITSQASSDDataMakerRec() {;} // dtor
  inline Int_t Raws() { return fSSDhRaws; }
  inline Int_t Recs() { return fSSDhRecs; }

private:

  Double_t GetSSDOccupancyRaws(TH1 *lHisto); 

  static const Int_t fgkSSDMODULES = 1698;      //total number of SSD modules
  static const Int_t fgkSSDMODULESLAYER5 = 748; //total number of SSD modules - layer5
  static const Int_t fgkSSDMODULESLAYER6 = 950; //total number of SSD modules - layer6

  AliITSQADataMakerRec *fAliITSQADataMakerRec;  //pointer to the main ctor
  Bool_t  fkOnline;                             //online (1) or offline (0) use
  Int_t   fLDC;                                 //LDC number (0 for offline, 1 to 4 for online) 
  Int_t   fSSDhRaws;                            // number of histo booked for Raws SSD
  Int_t   fSSDhRecs;                            // number of histo booked for Recs SSD
  Int_t   fRawsOffset;                          // number of histo booked when SSD start
  Int_t   fRecsOffset;                          // number of histo booked when SSD start
 
  TH1D *fHistSSDRawSignalModule[fgkSSDMODULES]; //raw signal vs strip number - SSD
  ClassDef(AliITSQASSDDataMakerRec,1)           // description 

};

#endif


