#ifndef ALIANALYSISTASKSPECTRAALLCHNANOAOD_H
#define ALIANALYSISTASKSPECTRAALLCHNANOAOD_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliAnalysisTaskSpectraAllChNanoAOD
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

class AliAnalysisTaskSpectraAllChNanoAOD : public AliAnalysisTaskSE
{
 public:
  // constructors
 AliAnalysisTaskSpectraAllChNanoAOD() : AliAnalysisTaskSE(),
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
    fnNchBins(200)
      {}
  AliAnalysisTaskSpectraAllChNanoAOD(const char *name);
  virtual ~AliAnalysisTaskSpectraAllChNanoAOD() {
    Printf("calling detructor of AliAnalysisTaskSpectraAllChNanoAOD - To be implemented");
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


  Int_t GetNanoTrackID(AliVTrack * track) ;
  
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
  AliAnalysisTaskSpectraAllChNanoAOD(const AliAnalysisTaskSpectraAllChNanoAOD&);
  AliAnalysisTaskSpectraAllChNanoAOD& operator=(const AliAnalysisTaskSpectraAllChNanoAOD&);
  
  ClassDef(AliAnalysisTaskSpectraAllChNanoAOD, 6);
};

#endif
