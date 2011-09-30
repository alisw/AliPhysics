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
#include "AliSpectraAODHistoManager.h"
#include "AliAnalysisTaskSE.h"


class AliAnalysisTaskSpectraAOD : public AliAnalysisTaskSE
{
public:

   // constructors
   AliAnalysisTaskSpectraAOD() : AliAnalysisTaskSE(), fAOD(0), fHistMan(0), fTrackCuts(0), fEventCuts(0), fIsMC(0), fPIDResponse(0), fNSigmaPID(0), fYCut(0) {}
   AliAnalysisTaskSpectraAOD(const char *name);
   virtual ~AliAnalysisTaskSpectraAOD() {}

   void SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; };
   Bool_t GetIsMC()           const           { return fIsMC;};
   void SetNSigmaForIdentification (Double_t sigma ) { fNSigmaPID = sigma; }
   Double_t GetNSigmaForIdentification () const {return fNSigmaPID; }
   void SetYCut (Double_t y ) { fYCut = y; }
   Double_t GetYCut () const {return fYCut; }

   virtual void   UserCreateOutputObjects();
   Bool_t         CheckYCut(AliSpectraNameSpace::AODParticleSpecies_t species, AliAODTrack* track) const;
   Bool_t         CheckYCut(AliAODMCParticle* particle) const;
   virtual void   UserExec(Option_t *option);
   virtual void   Terminate(Option_t *);
   void SetTrackCuts(AliSpectraAODTrackCuts * tc)   {      fTrackCuts = tc;   }
   void SetEventCuts(AliSpectraAODEventCuts * vc)   {      fEventCuts = vc;   }

private:

   AliAODEvent           * fAOD;         //! AOD object
   AliSpectraAODHistoManager      * fHistMan;       // Histogram Manager
   AliSpectraAODTrackCuts      * fTrackCuts;     // Track Cuts
   AliSpectraAODEventCuts      * fEventCuts;     // Event Cuts
   Bool_t          fIsMC;// true if processing MC
   AliPIDResponse                        *fPIDResponse;     // ! PID response object
   Double_t        fNSigmaPID; // Maximum number of sigmas allowed in particle identification
   Double_t        fYCut; // Maximum rapidity - calculated from identified particle's mass
   AliAnalysisTaskSpectraAOD(const AliAnalysisTaskSpectraAOD&);
   AliAnalysisTaskSpectraAOD& operator=(const AliAnalysisTaskSpectraAOD&);

   ClassDef(AliAnalysisTaskSpectraAOD, 1);
};

#endif
