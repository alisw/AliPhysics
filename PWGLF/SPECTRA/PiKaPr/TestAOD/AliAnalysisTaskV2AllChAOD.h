#ifndef ALIANALYSISTASKV2ALLCHAOD_H
#define ALIANALYSISTASKV2ALLCHAOD_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliAnalysisTaskV2AllChAOD
//
//
//
//
// Author: Leonardo Milano, CERN
//-------------------------------------------------------------------------

class AliAODEvent;
class AliSpectraAODTrackCuts;
class AliSpectraAODEventCuts;

#include "AliAnalysisTaskSE.h"
#include <TProfile.h>

class AliAnalysisTaskV2AllChAOD : public AliAnalysisTaskSE
{
 public:
  // constructors
 AliAnalysisTaskV2AllChAOD() : AliAnalysisTaskSE(),
    fAOD(0x0),
    fTrackCuts(0x0),
    fEventCuts(0x0),
    fIsMC(0),
    fCharge(0),
    fVZEROside(0),
    fOutput(0x0),
    fOutput_lq(0x0),
    fOutput_sq(0x0),
    fnCentBins(20),
    fnQvecBins(40),
    fIsQvecCalibMode(0),
    fQvecUpperLim(100),
    fCutLargeQperc(9.),
    fCutSmallQperc(10.),
    fEtaGapMin(-0.5),
    fEtaGapMax(0.5)
      {}
  AliAnalysisTaskV2AllChAOD(const char *name);
  virtual ~AliAnalysisTaskV2AllChAOD() {
    Printf("calling detructor of AliAnalysisTaskV2AllChAOD - To be implemented");
  }
  
  void SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; };
  Bool_t GetIsMC()           const           { return fIsMC;};
 
  void SetCharge(Int_t charge = 0)    {fCharge = charge; };
  Int_t GetCharge()           const           { return fCharge;};
  
  void SetVZEROside(Int_t side = 0)    {fVZEROside = side; };
  Int_t GetVZEROside()           const           { return fVZEROside;};
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  AliSpectraAODTrackCuts      * GetTrackCuts()         {  return fTrackCuts; }
  AliSpectraAODEventCuts      * GetEventCuts()         {  return fEventCuts; }
  TList                          * GetOutputList()         { return fOutput; }
  
  void SetTrackCuts(AliSpectraAODTrackCuts * tc)       { fTrackCuts = tc; }
  void SetEventCuts(AliSpectraAODEventCuts * vc)       { fEventCuts = vc; }
  void SetnCentBins(Int_t val)                             { fnCentBins = val; }
  void SetnQvecBins(Int_t val)                             { fnQvecBins = val; }
  void SetQvecCalibMode(Bool_t mode)                  { fIsQvecCalibMode = mode; }
  void SetQvecUpperLimit(Double_t val)                { fQvecUpperLim = val; }
  
  Bool_t GetDCA(const AliAODTrack* trk, Double_t * p);

  void     SetEtaGap(Float_t etamin,Float_t etamax)   { fEtaGapMin = etamin; fEtaGapMax = etamax; }
  
 private:
  
  AliAODEvent                   * fAOD;                         //! AOD object
  AliSpectraAODTrackCuts      * fTrackCuts;                   // Track Cuts
  AliSpectraAODEventCuts      * fEventCuts;                   // Event Cuts
  Bool_t                          fIsMC;                         // true if processing MC
  Int_t                            fCharge;                      // charge to be selected
  Int_t                            fVZEROside;                  // 0: VZERO-A 1: VZERO-C
  TList                          * fOutput;                     // output list
  TList                          * fOutput_lq;                  // output list large Q
  TList                          * fOutput_sq;                  // output list small Q
  Int_t                            fnCentBins;                  // number of bins for the centrality axis
  Int_t                            fnQvecBins;                 // number of bins for the q vector axis
  Bool_t                           fIsQvecCalibMode;          //calib mode for Qvector percentile
  Double_t                         fQvecUpperLim;             //Upper limit for Qvector
  
  Int_t                            fCutLargeQperc; // cut on 10% large Q-vec events
  Int_t                            fCutSmallQperc; // cut on 10% small Q-vec events
  
  Double_t fEtaGapMin;
  Double_t fEtaGapMax;
  
  //output object
  TProfile*     fResSP;             //! resolution
    
  TProfile*     fv2SPGap1A[9];         //! v2{2} eta gap 1 for all events
  TH2F*     fh2v2SPGap1A[9];         //! v2{2} eta gap 1 for all events
  TProfile*     fv2SPGap1B[9];         //! v2{2} eta gap 1 for all events
  TH2F*     fh2v2SPGap1B[9];         //! v2{2} eta gap 1 for all events

  TProfile*     fSinGap1A[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1A[9];      //! <cos> vs pT gap 1
  TProfile*     fSinGap1B[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1B[9];      //! <cos> vs pT gap 1

  //large q
  TProfile*     fResSP_lq;             //! resolution
  TProfile*     fv2SPGap1A_lq[9];         //! v2{2} eta gap 1 for all events
  TProfile*     fv2SPGap1B_lq[9];         //! v2{2} eta gap 1 for all events
  TProfile*     fSinGap1A_lq[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1A_lq[9];      //! <cos> vs pT gap 1
  TProfile*     fSinGap1B_lq[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1B_lq[9];      //! <cos> vs pT gap 1
  
  //small q
  TProfile*     fResSP_sq;             //! resolution
  TProfile*     fv2SPGap1A_sq[9];         //! v2{2} eta gap 1 for all events
  TProfile*     fv2SPGap1B_sq[9];         //! v2{2} eta gap 1 for all events
  TProfile*     fSinGap1A_sq[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1A_sq[9];      //! <cos> vs pT gap 1
  TProfile*     fSinGap1B_sq[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1B_sq[9];      //! <cos> vs pT gap 1
  
  AliAnalysisTaskV2AllChAOD(const AliAnalysisTaskV2AllChAOD&);
  AliAnalysisTaskV2AllChAOD& operator=(const AliAnalysisTaskV2AllChAOD&);
  
  ClassDef(AliAnalysisTaskV2AllChAOD, 1);
};

#endif
