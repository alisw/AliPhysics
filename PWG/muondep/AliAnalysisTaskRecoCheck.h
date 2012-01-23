#ifndef AliAnalysisTaskRecoCheck_cxx
#define AliAnalysisTaskRecoCheck_cxx

/* $Id$ */ 

// analysis task for decoding reconstructed tracks and kinematics (AliMUONRecoCheck)
// Authors: Bogdan Vulpescu

class TTree;
class TClonesArray;
class AliESDEvent;

#include "AliAnalysisTask.h"

class AliAnalysisTaskRecoCheck : public AliAnalysisTask {
 public:
  AliAnalysisTaskRecoCheck(const char *name = "AliAnalysisTaskRecoCheck");
  virtual ~AliAnalysisTaskRecoCheck() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetL3Current(Float_t val) { fL3Current = val; }
  
 private:
  AliESDEvent  *fESDEvent;    // ESD object
  TTree        *fTree;        // output tree
  TClonesArray *fArray1Mu;    // per event array of single muons
  TClonesArray *fArray2Mu;    // per event array of muon pairs
  Float_t       fL3Current;   // current in the L3 magnet for field calculation

  AliAnalysisTaskRecoCheck(const AliAnalysisTaskRecoCheck&); // not implemented
  AliAnalysisTaskRecoCheck& operator=(const AliAnalysisTaskRecoCheck&); // not implemented

  ClassDef(AliAnalysisTaskRecoCheck, 1); // ESD and Kine analysis with RecoCheck
};

#endif
