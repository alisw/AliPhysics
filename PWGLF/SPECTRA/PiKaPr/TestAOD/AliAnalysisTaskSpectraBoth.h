#ifndef ALIANALYSISTASKSPECTRABOTH_H
#define ALIANALYSISTASKSPECTRABOTH_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliAnalysisTaskSpectraBoth
//
//
//
//
// Author: Michele Floris, CERN
//-------------------------------------------------------------------------

class TH1F;
class TH2F;
class AliVEvent;
class AliSpectraBothHistoManager;
class AliSpectraBothTrackCuts;
class AliSpectraBothEventCuts;
class AliSpectraBothPID;
class AliESDtrackCuts;
class AliGenEventHeader;

#include "AliSpectraBothHistoManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"

class AliAnalysisTaskSpectraBoth : public AliAnalysisTaskSE
{
public:

   // constructors
  AliAnalysisTaskSpectraBoth() : AliAnalysisTaskSE(), fAOD(0), fHistMan(0), fTrackCuts(0), fEventCuts(0), fPID(0), fIsMC(0), fNRebin(0),fUseMinSigma(0),fCuts(0),fdotheMCLoopAfterEventCuts(0)
 {}
  AliAnalysisTaskSpectraBoth(const char *name);
   virtual ~AliAnalysisTaskSpectraBoth() {}

   void SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; };
   Bool_t GetIsMC()           const           { return fIsMC;};

   virtual void   UserCreateOutputObjects();
   virtual void   UserExec(Option_t *option);
   virtual void   Terminate(Option_t *);

   AliSpectraBothHistoManager * GetHistoManager()         {  return fHistMan; }
   AliSpectraBothTrackCuts * GetTrackCuts()         {  return fTrackCuts; }
   AliSpectraBothEventCuts * GetEventCuts()         {  return fEventCuts; }
   AliSpectraBothPID * GetPID()         {  return fPID; }
   
   void SetTrackCuts(AliSpectraBothTrackCuts * tc)   {   fTrackCuts = tc;   }
   void SetEventCuts(AliSpectraBothEventCuts * vc)   {   fEventCuts = vc;   }
   void SetPID      (AliSpectraBothPID      * pid)   {   fPID       = pid;  }
   void SetNRebin(Int_t nreb){fNRebin=nreb;}
   void SetUseMinSigma (Bool_t flag) {fUseMinSigma=flag;}
   Int_t   GetNRebin() const {return fNRebin;}
   void SetAliESDtrackCuts(AliESDtrackCuts*  cuts ){fCuts=cuts;}
   void SetdotheMCLoopAfterEventCuts (Bool_t flag) {fdotheMCLoopAfterEventCuts=flag;}
   Bool_t GetdotheMCLoopAfterEventCuts () const {return fdotheMCLoopAfterEventCuts;}
private:

   AliVEvent           * fAOD;         //! AOD object
   AliSpectraBothHistoManager      * fHistMan;       // Histogram Manager
   AliSpectraBothTrackCuts      * fTrackCuts;     // Track Cuts
   AliSpectraBothEventCuts      * fEventCuts;     // Event Cuts
   AliSpectraBothPID             * fPID;// PID class
   Bool_t          fIsMC;// true if processing MC
   Int_t      fNRebin; //rebin of histos
   Bool_t fUseMinSigma; // if true use min sigma 
     AliESDtrackCuts *fCuts; // ESD track cuts 
    Bool_t fdotheMCLoopAfterEventCuts; // if true first check the ESD event cuts than loop over MC info , if flase other approach     

   AliAnalysisTaskSpectraBoth(const AliAnalysisTaskSpectraBoth&);
   AliAnalysisTaskSpectraBoth& operator=(const AliAnalysisTaskSpectraBoth&);

   ClassDef(AliAnalysisTaskSpectraBoth, 2);
};

#endif
