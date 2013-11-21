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

class TH1F;
class TH2F;
class AliAODEvent;
class AliSpectraAODTrackCuts;
class AliSpectraAODEventCuts;
#include "AliAnalysisTaskSE.h"


class AliAnalysisTaskSpectraAllChAOD : public AliAnalysisTaskSE
{
 public:
  
  // constructors
 AliAnalysisTaskSpectraAllChAOD() : AliAnalysisTaskSE(), fAOD(0), fTrackCuts(0), fEventCuts(0), fIsMC(0), fOutput(0)
    {}
  AliAnalysisTaskSpectraAllChAOD(const char *name);
  virtual ~AliAnalysisTaskSpectraAllChAOD() {}

  void SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; };
  Bool_t GetIsMC()           const           { return fIsMC;};
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  AliSpectraAODTrackCuts * GetTrackCuts()         {  return fTrackCuts; }
  AliSpectraAODEventCuts * GetEventCuts()         {  return fEventCuts; }
  TList *GetOutputList()       {return fOutput;};
 
  void SetTrackCuts(AliSpectraAODTrackCuts * tc)   {   fTrackCuts = tc;   }
  void SetEventCuts(AliSpectraAODEventCuts * vc)   {   fEventCuts = vc;   }
  
 private:
  
  AliAODEvent           * fAOD;         //! AOD object
  AliSpectraAODTrackCuts      * fTrackCuts;     // Track Cuts
  AliSpectraAODEventCuts      * fEventCuts;     // Event Cuts
  Bool_t          fIsMC;// true if processing MC
  TList          *fOutput;// output list
  AliAnalysisTaskSpectraAllChAOD(const AliAnalysisTaskSpectraAllChAOD&);
  AliAnalysisTaskSpectraAllChAOD& operator=(const AliAnalysisTaskSpectraAllChAOD&);
  
  ClassDef(AliAnalysisTaskSpectraAllChAOD, 1);
};

#endif
