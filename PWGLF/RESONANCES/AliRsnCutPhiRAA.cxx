//
// This cut implements all the checks done to accept a track as a Kaon
// for the PbPb analysis using 2010 runs.
// It is based on standard cuts on track quality and nsigma cuts
// with respect to the TPC and TOF signals for the PID.
//

#include <Riostream.h>

#include "AliLog.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliESDpid.h"
#include "AliAODpidUtil.h"
#include "AliRsnCutPhiRAA.h"
#include "AliESDtrackCuts.h"

ClassImp(AliRsnCutPhiRAA)

//__________________________________________________________________________________________________
AliRsnCutPhiRAA::AliRsnCutPhiRAA
(const char *name) :
   AliRsnCut(name, AliRsnTarget::kDaughter),
   fMode(k2011_1),
   fCutQuality(Form("%s_quality", name)),
   cut1(AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 1)),
   cut2(AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 0)),
   cut3(AliESDtrackCuts::GetStandardITSTPCTrackCuts2010())
{
//
// Constructor
// Initialize the contained cuts and sets defaults
//

//  cut1
  cut1->SetEtaRange(-0.8, 0.8);
  cut1->SetPtRange(0.30,1e10);
//  cut2
  cut2->SetEtaRange(-0.8, 0.8);
  cut2->SetPtRange(0.30,1e10);
//  cut3
  cut3->SetEtaRange(-0.8, 0.8);
  cut3->SetPtRange(0.30,1e10);

  //V end

}

//__________________________________________________________________________________________________
AliRsnCutPhiRAA::AliRsnCutPhiRAA(const AliRsnCutPhiRAA &copy) :
   AliRsnCut(copy),
   fMode(copy.fMode),
   fCutQuality(copy.fCutQuality),
   cut1(copy.cut1),
   cut2(copy.cut2),
   cut3(copy.cut3)

{
//
// Copy constructor
//
}

//__________________________________________________________________________________________________
AliRsnCutPhiRAA &AliRsnCutPhiRAA::operator=(const AliRsnCutPhiRAA &copy)
{
//
// Assignment operator
//

   AliRsnCut::operator=(copy);
   if (this == &copy)
      return *this;
   fMode = copy.fMode;
   fCutQuality = copy.fCutQuality;

   return *this;
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutPhiRAA::IsSelected(TObject *obj)
{
//
// Global check
//

   // coherence check
   if (!TargetOK(obj)) return kFALSE;

   // check track
   AliVTrack *track = fDaughter->Ref2Vtrack();
   if (!track) {
      if (!fDaughter->GetRef()) AliWarning("NULL ref");
      return kFALSE;
   }

   // initialize check variables
   Bool_t   accept = kFALSE;

   // decide cut result depending on mode
   switch (fMode) {
     case k2010:
       fCutQuality.SetESDtrackCuts(cut3);
       if (fCutQuality.IsSelected(obj)) accept = kTRUE;
       ::Info("AnalysisSetup", "ESD cut 2010!!! ");
     break;
     case k2011_0:
 	   fCutQuality.SetESDtrackCuts(cut2);
       if (fCutQuality.IsSelected(obj)) accept = kTRUE;
     break;
     case k2011_1:
 	   fCutQuality.SetESDtrackCuts(cut1);
       if (fCutQuality.IsSelected(obj)) accept = kTRUE;
       ::Info("AnalysisSetup", "ESD cut 2011 ");
     break;
     case k2011_1_05:
    	  cut1->SetPtRange(0.50,1e10);
 	   fCutQuality.SetESDtrackCuts(cut1);
       if (fCutQuality.IsSelected(obj)) accept = kTRUE;
       ::Info("AnalysisSetup", "ESD cut 2011 pT > 0.5 GeV/c");
     break;
     case k2011_1_075:
    	  cut1->SetPtRange(0.75,1e10);
 	   fCutQuality.SetESDtrackCuts(cut1);
       if (fCutQuality.IsSelected(obj)) accept = kTRUE;
       ::Info("AnalysisSetup", "ESD cut 2011 pT > 0.75 GeV/c ");
     break;
     default:
         AliDebugClass(1, Form("[%s] Wrong mode", GetName()));
         accept = kFALSE;
   }

   AliDebugClass(1, Form("[%s] Track %s", GetName(), (accept ? "accepted" : "rejected")));
   return accept;
}
