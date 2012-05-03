#ifndef ALIANALYSISTASKSPECTRAAOD_H
#define ALIANALYSISTASKSPECTRAAOD_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliAnalysisTaskSpectraAOD
//
//
//
//
// Author: Michele Floris, CERN
//-------------------------------------------------------------------------

class TH1F;
class TH2F;
class AliAODEvent;
class AliSpectraAODHistoManager;
class AliSpectraAODTrackCuts;
class AliSpectraAODEventCuts;
class AliSpectraAODPID;
#include "AliSpectraAODHistoManager.h"
#include "AliAnalysisTaskSE.h"


class AliAnalysisTaskSpectraAOD : public AliAnalysisTaskSE
{
public:

   // constructors
  AliAnalysisTaskSpectraAOD() : AliAnalysisTaskSE(), fAOD(0), fHistMan(0), fTrackCuts(0), fEventCuts(0), fPID(0), fIsMC(0)
 {}
   AliAnalysisTaskSpectraAOD(const char *name);
   virtual ~AliAnalysisTaskSpectraAOD() {}

   void SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; };
   Bool_t GetIsMC()           const           { return fIsMC;};

   virtual void   UserCreateOutputObjects();
   virtual void   UserExec(Option_t *option);
   virtual void   Terminate(Option_t *);
   void SetTrackCuts(AliSpectraAODTrackCuts * tc)   {   fTrackCuts = tc;   }
   void SetEventCuts(AliSpectraAODEventCuts * vc)   {   fEventCuts = vc;   }
   void SetPID      (AliSpectraAODPID      * pid)   {   fPID       = pid;  }

private:

   AliAODEvent           * fAOD;         //! AOD object
   AliSpectraAODHistoManager      * fHistMan;       // Histogram Manager
   AliSpectraAODTrackCuts      * fTrackCuts;     // Track Cuts
   AliSpectraAODEventCuts      * fEventCuts;     // Event Cuts
   AliSpectraAODPID             * fPID;// PID class
   Bool_t          fIsMC;// true if processing MC
   AliAnalysisTaskSpectraAOD(const AliAnalysisTaskSpectraAOD&);
   AliAnalysisTaskSpectraAOD& operator=(const AliAnalysisTaskSpectraAOD&);

   ClassDef(AliAnalysisTaskSpectraAOD, 1);
};

#endif
