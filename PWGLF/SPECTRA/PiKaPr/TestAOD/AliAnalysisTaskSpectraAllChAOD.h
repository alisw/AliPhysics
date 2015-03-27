#ifndef ALIANALYSISTASKSPECTRAALLCHAOD_H
#define ALIANALYSISTASKSPECTRAALLCHAOD_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliAnalysisTaskSpectraAllChAOD
//
//
//
//
// Author: Leonardo Milano, CERN
//-------------------------------------------------------------------------

class AliAODEvent;
class AliSpectraAODTrackCuts;
class AliSpectraAODEventCuts;
class AliHelperPID;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSpectraAllChAOD : public AliAnalysisTaskSE
{
 public:
  // constructors
 AliAnalysisTaskSpectraAllChAOD() : AliAnalysisTaskSE(),
    fAOD(0x0),
    fTrackCuts(0x0),
    fEventCuts(0x0),
    fHelperPID(0x0),
    fIsMC(0),
    fDoDoubleCounting(0),
    fFillOnlyEvents(0),
    fCharge(0),
    fVZEROside(0),
    fOutput(0x0),
    fnCentBins(20),
    fnQvecBins(40),
    fnNchBins(200),
    fIsQvecCalibMode(0),
    fQvecUpperLim(100),
    fIsAOD160(1),
    fnDCABins(60),
    fDCAmin(-3),
    fDCAmax(3),
    fDCAzCut(0),
    fQst(1),
    fQtrk(0),
    fQgenType(0),
    fDoCentrSystCentrality(0),
    fEtaGap(0.5)
      {}
  AliAnalysisTaskSpectraAllChAOD(const char *name);
  virtual ~AliAnalysisTaskSpectraAllChAOD() {
    Printf("calling detructor of AliAnalysisTaskSpectraAllChAOD - To be implemented");
  }
  
  void SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; };
  Bool_t GetIsMC()           const           { return fIsMC;};
 
  void SetDoDoubleCounting(Bool_t doDoubleCounting = kFALSE)    {fDoDoubleCounting = doDoubleCounting; };
  Bool_t GetDoDoubleCounting()           const           { return fDoDoubleCounting;};
 
  void SetFillOnlyEvents(Bool_t fillOnlyEvents = kFALSE)    {fFillOnlyEvents = fillOnlyEvents; };
  Bool_t GetFillOnlyEvents()           const           { return fFillOnlyEvents;};
 
  void SetCharge(Int_t charge = 0)    {fCharge = charge; };
  Int_t GetCharge()           const           { return fCharge;};
  
  void SetVZEROside(Int_t side = 0)    {fVZEROside = side; };
  Int_t GetVZEROside()           const           { return fVZEROside;};
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  AliSpectraAODTrackCuts      * GetTrackCuts()         {  return fTrackCuts; }
  AliSpectraAODEventCuts      * GetEventCuts()         {  return fEventCuts; }
  AliHelperPID                   * GetHelperPID()          { return fHelperPID; }
  TList                          * GetOutputList()         { return fOutput; }
  
  void SetTrackCuts(AliSpectraAODTrackCuts * tc)       { fTrackCuts = tc; }
  void SetEventCuts(AliSpectraAODEventCuts * vc)       { fEventCuts = vc; }
  void SetHelperPID(AliHelperPID* pid)                     { fHelperPID = pid; }
  void SetnCentBins(Int_t val)                             { fnCentBins = val; }
  void SetnQvecBins(Int_t val)                             { fnQvecBins = val; }
  void SetnNchBins(Int_t val)                             { fnNchBins = val; }
  void SetQvecCalibMode(Bool_t mode)                  { fIsQvecCalibMode = mode; }
  void SetQvecUpperLimit(Double_t val)                { fQvecUpperLim = val; }
  
  void SetIsAOD160(Bool_t aod)                        { fIsAOD160 = aod; }
  void SetnDCABin(Int_t val)                          { fnDCABins = val; } 
  void SetDCAmin(Double_t val)                        { fDCAmin = val; }
  void SetDCAmax(Double_t val)                        { fDCAmax = val; }
  Bool_t GetDCA(const AliAODTrack* trk, Double_t * p);
  void SetDCAzCut(Double_t val)                        { fDCAzCut = val; }
  
  void SetQStack(Int_t val) { fQst = val ; } // type==0 q-reco - type==1 qgen
  void SetQTrack(Int_t val) { fQtrk = val ; } // type==0 q-reco - type==1 qgen
  void SetQgenType(Int_t val) { fQgenType = val ; } // type==0 qgen from tracks - type==1 qgen from vzero
  
  void SetDoCentrSystCentrality(Bool_t val) { fDoCentrSystCentrality = val; } //enable systematic for centrality
  
  void SetEtaGap(Double_t val) { fEtaGap = val; }

 private:
  
  AliAODEvent                   * fAOD;                         //! AOD object
  AliSpectraAODTrackCuts      * fTrackCuts;                   // Track Cuts
  AliSpectraAODEventCuts      * fEventCuts;                   // Event Cuts
  AliHelperPID                   * fHelperPID;                    // points to class for PID
  Bool_t                          fIsMC;                         // true if processing MC
  Bool_t                          fDoDoubleCounting;           // true is double counting for Nsigma accepetd
  Bool_t                          fFillOnlyEvents;               // if true fill only NSparseHistEv
  Int_t                            fCharge;                      // charge to be selected
  Int_t                            fVZEROside;                  // 0: VZERO-A 1: VZERO-C
  TList                          * fOutput;                     // output list
  Int_t                            fnCentBins;                  // number of bins for the centrality axis
  Int_t                            fnQvecBins;                 // number of bins for the q vector axis
  Int_t                            fnNchBins;                 // number of bins for the Nch axis
  Bool_t                           fIsQvecCalibMode;          //calib mode for Qvector percentile
  Double_t                         fQvecUpperLim;             //Upper limit for Qvector
  Bool_t                           fIsAOD160;              // enable DCA for AOD160
  Int_t                            fnDCABins;              // number of bins for DCA axis
  Double_t                         fDCAmin;                // min DCA value
  Double_t                         fDCAmax;                // max DCA value
  Double_t                         fDCAzCut;               //cut on DCA z
  
  Int_t fQst; // type==0 q-reco - type==1 qgen
  Int_t fQtrk; // type==0 q-reco - type==1 qgen
  Int_t fQgenType; // type==0 qgen from tracks - type==1 qgen from vzero
  
  Bool_t fDoCentrSystCentrality; //systematic check on centrality estimation
  
  Double_t fEtaGap;
  
  AliAnalysisTaskSpectraAllChAOD(const AliAnalysisTaskSpectraAllChAOD&);
  AliAnalysisTaskSpectraAllChAOD& operator=(const AliAnalysisTaskSpectraAllChAOD&);
  
  ClassDef(AliAnalysisTaskSpectraAllChAOD, 13);
};

#endif
