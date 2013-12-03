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
    fOutput(0x0)
      {}
  AliAnalysisTaskSpectraAllChAOD(const char *name);
  virtual ~AliAnalysisTaskSpectraAllChAOD() {}
  
  void SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; };
  Bool_t GetIsMC()           const           { return fIsMC;};
  
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
  
 private:
  
  AliAODEvent                   * fAOD;           //! AOD object
  AliSpectraAODTrackCuts      * fTrackCuts;     // Track Cuts
  AliSpectraAODEventCuts      * fEventCuts;     // Event Cuts
  AliHelperPID                   * fHelperPID;      // points to class for PID
  Bool_t                          fIsMC;           // true if processing MC
  TList                          * fOutput;        // output list
  AliAnalysisTaskSpectraAllChAOD(const AliAnalysisTaskSpectraAllChAOD&);
  AliAnalysisTaskSpectraAllChAOD& operator=(const AliAnalysisTaskSpectraAllChAOD&);
  
  ClassDef(AliAnalysisTaskSpectraAllChAOD, 2);
};

#endif
