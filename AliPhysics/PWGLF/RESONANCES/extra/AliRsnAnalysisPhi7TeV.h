//
// Header file for implementation of data analysis aft 900 GeV
//
// Author: A. Pulvirenti
//

#ifndef ALIRSNANALYSISPHI7TEV_H
#define ALIRSNANALYSISPHI7TEV_H

#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"
#include "AliPID.h"

class TH1I;
class TH1F;
class TH3F;

class AliStack;
class AliESDEvent;
class AliESDVertex;
class AliESDpid;
class AliESDtrack;
class AliCFContainer;

class AliRsnAnalysisPhi7TeV : public AliAnalysisTaskSE {
public:

   enum EVertexType {
      kGoodTracksPrimaryVertex = 0,
      kGoodSPDPrimaryVertex    = 1,
      kFarTracksPrimaryVertex  = 2,
      kFarSPDPrimaryVertex     = 3,
      kNoGoodPrimaryVertex     = 4,
      kEmptyEvent              = 5,
      kVertexTypes
   };
   
   enum ETrackType {
      kITSStandAlone,
      kTPConly,
      kTPCwithTOF,
      kTrackTypes
   };

   AliRsnAnalysisPhi7TeV(const char *name = "Phi7TeV", Bool_t isMC = kFALSE);
   AliRsnAnalysisPhi7TeV(const AliRsnAnalysisPhi7TeV& copy);
   AliRsnAnalysisPhi7TeV& operator=(const AliRsnAnalysisPhi7TeV& copy);
   virtual ~AliRsnAnalysisPhi7TeV();

   void             SetMC       (Bool_t yn = kTRUE)     {fIsMC     = yn;}
   void             SetAddITSTPC(Bool_t yn = kTRUE)     {fAddTPC   = yn;}
   void             SetAddITSSA (Bool_t yn = kTRUE)     {fAddITS   = yn;}
   void             SetCheckITS (Bool_t yn = kTRUE)     {fCheckITS = yn;}
   void             SetCheckTPC (Bool_t yn = kTRUE)     {fCheckTPC = yn;}
   void             SetCheckTOF (Bool_t yn = kTRUE)     {fCheckTOF = yn;}
   void             SetStep     (Int_t step)            {fStep     = step;}

   void             SetPhiMass (Double_t value)         {fPhiMass  = value;}
   void             SetKaonMass(Double_t value)         {fKaonMass = value;}

   void             SetMaxVz(Double_t v)                {fMaxVz = v;}

   void             SetITSNSigma           (Double_t v) {fITSNSigma            = v;}
   void             SetTPCNSigmaLow        (Double_t v) {fTPCNSigma[0]         = v;}
   void             SetTPCNSigmaHigh       (Double_t v) {fTPCNSigma[1]         = v;}
   void             SetTPCMomentumThreshold(Double_t v) {fTPCMomentumThreshold = v;}
   void             SetTOFNSigma           (Double_t v) {fTOFNSigma            = v;}
   
   void             SetCutsTPC(AliESDtrackCuts *cuts)   {fESDtrackCutsTPC = cuts;}
   void             SetCutsITS(AliESDtrackCuts *cuts)   {fESDtrackCutsITS = cuts;}
   void             SetESDpid (AliESDpid       *pid )   {fESDpid          = pid ;}
   AliESDtrackCuts* GetCutsTPC()                        {return fESDtrackCutsTPC;}
   AliESDtrackCuts* GetCutsITS()                        {return fESDtrackCutsITS;}
   AliESDpid*       GetESDpid()                         {return fESDpid;         }
   
   virtual void     UserCreateOutputObjects();
   virtual void     UserExec(Option_t *option = "");
   virtual void     Terminate(Option_t *option = "");

   EVertexType      EventEval(AliESDEvent *esd);
   
   Bool_t           IsTPC    (AliESDtrack *track);
   Bool_t           IsITS    (AliESDtrack *track);
   Bool_t           MatchTOF (AliESDtrack *track);
   ETrackType       TrackType(AliESDtrack *track);
   
   Bool_t           OkQualityITS(AliESDtrack *track) {return fESDtrackCutsITS->IsSelected(track);}
   Bool_t           OkQualityTPC(AliESDtrack *track) {return fESDtrackCutsTPC->IsSelected(track);}
   Bool_t           OkQuality   (AliESDtrack *track);
   Bool_t           OkPIDITS    (AliESDtrack *track, AliPID::EParticleType pid = AliPID::kKaon);
   Bool_t           OkPIDTPC    (AliESDtrack *track, AliPID::EParticleType pid = AliPID::kKaon);
   Bool_t           OkPIDTOF    (AliESDtrack *track, AliPID::EParticleType pid = AliPID::kKaon);
   Bool_t           OkPID       (AliESDtrack *track, AliPID::EParticleType pid = AliPID::kKaon);

private:

   void  ProcessESD      (AliESDEvent *esd, AliMCEvent *mc);
   void  ProcessMC       (AliESDEvent *esd, AliMCEvent *mc);
   Int_t MatchedTrack    (AliESDEvent *esd, Int_t label, Int_t &npassed, Int_t start = 0);
   Int_t BestMatchedTrack(AliESDEvent *esd, Int_t label);

   Bool_t           fIsMC;                      // use MC or data?
   Bool_t           fAddTPC;                    // include TPC+ITS tracks?
   Bool_t           fAddITS;                    // add ITS standalone tracks?
   Bool_t           fCheckITS;                  // check ITS PID (only for ITS-SA)?
   Bool_t           fCheckTPC;                  // check TPC PID (only for TPC+ITS)?
   Bool_t           fCheckTOF;                  // check TOF PID (only for TPC+ITS)?
   Int_t            fStep;                      // step for progress message
   
   Double_t         fPhiMass;                   // mass hypothesis for phi
   Double_t         fKaonMass;                  // mass hypothesis for kaon
                   
   Double_t         fMaxVz;                     // range in Z of primary vertex w.r. to origin
                   
   Double_t         fITSNSigma;                 // range for ITS dE/dx band (in sigmas)
   Double_t         fTPCNSigma[2];              // ranges for TPC dE/dx band (in sigmas)
   Double_t         fTPCMomentumThreshold;      // below this, [0] band is used, above, [1] bandd is used
   Double_t         fTOFNSigma;                 // range for TOF time band (in sigmas)
   
   AliESDtrackCuts *fESDtrackCutsTPC;           // ESD standard defined track cuts for TPC tracks
   AliESDtrackCuts *fESDtrackCutsITS;           // ESD standard defined track cuts for ITS-SA tracks
   AliESDpid       *fESDpid;                    // PID manager
                   
   TList           *fOutList;                   // list containing all task outputs
                   
   AliCFContainer  *fCFunlike;                  // CF container for unlike-sign pairs
   AliCFContainer  *fCFlikePP;                  // CF container for like-sign pairs ++
   AliCFContainer  *fCFlikeMM;                  // CF container for like-sign pairs --
   AliCFContainer  *fCFtrues;                   // CF container for true phis
   AliCFContainer  *fCFkaons;                   // CF container for kaons (monitoring)
                   
   TH1I            *fHEvents;                   // histogram of event types
   TH1F            *fVertexX[2];                // histogram of X coordinate of primary vertex ([0] = tracks, [1] = SPD)
   TH1F            *fVertexY[2];                // histogram of Y coordinate of primary vertex ([0] = tracks, [1] = SPD)
   TH1F            *fVertexZ[2];                // histogram of Z coordinate of primary vertex ([0] = tracks, [1] = SPD)

   // ROOT dictionary
   ClassDef(AliRsnAnalysisPhi7TeV, 1)
};

inline Bool_t AliRsnAnalysisPhi7TeV::IsTPC(AliESDtrack *vtrack)
{
//
// Checks if the track has the status flags required for a TPC+ITS track
//

   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   return vtrack->IsOn(AliESDtrack::kTPCin);
}

inline Bool_t AliRsnAnalysisPhi7TeV::IsITS(AliESDtrack *vtrack)
{
//
// Checks if the track has the status flags required for an ITS standalone track
//

   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   Bool_t isTPCin     = vtrack->IsOn(AliESDtrack::kTPCin);
   Bool_t isITSrefit  = vtrack->IsOn(AliESDtrack::kITSrefit);
   Bool_t isITSpureSA = vtrack->IsOn(AliESDtrack::kITSpureSA);
   Bool_t isITSpid    = vtrack->IsOn(AliESDtrack::kITSpid);

   return ((!isTPCin) && isITSrefit && (!isITSpureSA) && isITSpid);
}


inline Bool_t AliRsnAnalysisPhi7TeV::MatchTOF(AliESDtrack *vtrack)
{
//
// Checks if the track has matched the TOF detector
//

   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   // require a minimum length to have meaningful match
   if (vtrack->GetIntegratedLength() < 350.) return kFALSE;

   Bool_t isTOFout = vtrack->IsOn(AliESDtrack::kTOFout);
   Bool_t isTIME   = vtrack->IsOn(AliESDtrack::kTIME);

   return (isTOFout && isTIME);
}

inline AliRsnAnalysisPhi7TeV::ETrackType AliRsnAnalysisPhi7TeV::TrackType(AliESDtrack *track)
{
//
// Assigns the track type according to enum
//

   if (IsITS(track)) 
      return kITSStandAlone;
   else if (IsTPC(track)) {
      if (MatchTOF(track))
         return kTPCwithTOF;
      else
         return kTPConly;
   }
   else
      return kTrackTypes;
}

inline Bool_t AliRsnAnalysisPhi7TeV::OkQuality(AliESDtrack *track)
{
//
// Check track quality cut, which depends on track type
//

   ETrackType type = TrackType(track);
      
   switch (type) {
      case kITSStandAlone:
         return OkQualityITS(track);
      case kTPConly:
      case kTPCwithTOF:
         return OkQualityTPC(track);
      default:
         return kFALSE;
   }
}

inline Bool_t AliRsnAnalysisPhi7TeV::OkPID(AliESDtrack *track, AliPID::EParticleType pid)
{
//
// Check PID cut, which depends on track type
//

   ETrackType type = TrackType(track);
      
   switch (type) {
      case kITSStandAlone:
         if (fCheckITS) return OkPIDITS(track);
         else return kTRUE;
      case kTPConly:
         if (fCheckTPC) return OkPIDTPC(track);
         else return kTRUE;
      case kTPCwithTOF:
         if (fCheckTPC && fCheckTOF) return (OkPIDTPC(track, pid) && OkPIDTOF(track, pid));
         else if (fCheckTOF) return OkPIDTOF(track, pid);
         else if (fCheckTPC) return OkPIDTPC(track, pid);
         else return kTRUE;
      default:
         return kFALSE;
   }
}

#endif
