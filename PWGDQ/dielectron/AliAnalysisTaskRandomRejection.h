#ifndef ALIANALYSISTASKRANDOMREJECTION_H
#define ALIANALYSISTASKRANDOMREJECTION_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#########################################################
//#                                                       #
//# based on AliAnalysisTaskMultiDielectron (04.2015)     #
//#                                                       #
//# extension (by Patrick Reichelt)                       #
//#  -  random rejection probability of signal electrons  #
//#     due to the pair prefilter. needs to be included   #
//#     in the efficiency correction.                     #
//#     (see comments in the function implementations)    #
//#                                                       #
//#########################################################

#include "TList.h"
#include "TF1.h"      // added
#include "AliAnalysisTaskMultiDielectron.h"

class AliDielectron;
class TH1D;
class AliAnalysisCuts;
class AliTriggerAnalysis;

class AliAnalysisTaskRandomRejection : public AliAnalysisTaskMultiDielectron {
  
public:
  AliAnalysisTaskRandomRejection();
  AliAnalysisTaskRandomRejection(const char *name);
  virtual ~AliAnalysisTaskRandomRejection();


  virtual void UserExec(Option_t *option);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
  //temporary
//   virtual void NotifyRun(){AliDielectronPID::SetCorrVal((Double_t)fCurrentRunNumber);}
    
  ///
  /// _____ extension compared to AliAnalysisTaskMultiDielectron _____
  ///
  void CalcRandomPairs(AliDielectron* die);
  void FillHistogramsRandomPairs(AliDielectron* die, AliDielectronPair* part, Bool_t wasRejected);
  void FillHistogramsTestpart(AliDielectron* die, TObjArray* arrTracks1, Bool_t* bTracks1);
  void FillHistogramsDataEle(AliDielectron* die, TObjArray* arrTracks1, Bool_t* bTracks1);
  void FillFinalTrackArrays(AliVEvent * const ev, AliDielectron* die);
  void InitTestparticles();
  void SetPtFunc(TF1* func)    { fPtFunc=func; }
  void SetEtaMax(Double_t val) { fRndmEtaMax=val; }
  void SetNPtEtaPhi(UInt_t npt, UInt_t neta, UInt_t nphi) { if(npt>0)fNRndmPt=npt; if(neta>0)fNRndmEta=neta; if(nphi>0)fNRndmPhi=nphi; }

private:
  TF1*        fPtFunc;               //! ("!" to avoid streamer error) function used for random pt-distribution of testparticles [default is "exp(-x/3.)"]
  Double_t    fRndmEtaMax;           // edge of eta sampling (symmetric around eta=0) [default=+-0.9]
  UInt_t      fNRndmPt;              // pt and eta are sampled 'fNRndmPt*fNRndmEta' times
  UInt_t      fNRndmEta;             // 
  UInt_t      fNRndmPhi;             // phi is sampled 'fNRndmPt*fNRndmEta*fNRndmPhi' times
  //                                 // -> also gives the number of testparticles per event [default of 512 good for p-Pb].
  /// _____
  
protected:
  
  ///
  /// _____ extension compared to AliAnalysisTaskMultiDielectron _____
  ///
  TObjArray*  fTestparticles;        //! ("!" not needed!?) testparticles for random pairing
  TObjArray*  fFinalTracks[2];       //! ("!" needed!) track arrays which pass final analysis cuts // (in AliDielectron they are not pointers...)
  /// _____
    
  AliAnalysisTaskRandomRejection(const AliAnalysisTaskRandomRejection &c);
  AliAnalysisTaskRandomRejection& operator= (const AliAnalysisTaskRandomRejection &c);
  
  ClassDef(AliAnalysisTaskRandomRejection, 1); //Analysis Task handling multiple instances of AliDielectron
};
#endif
