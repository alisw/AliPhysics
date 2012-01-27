#include <Riostream.h>

#include <TChain.h>
#include <TH1.h>
#include <TH2.h>

#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"

#include "AliAnalysisTaskEventCount.h"

// analysis task creating basic QA plots for resonance particles
// Author: Ayben Karasu Uysal


ClassImp(AliAnalysisTaskEventCount)

//________________________________________________________________________
AliAnalysisTaskEventCount::AliAnalysisTaskEventCount(const char *name) :
   AliAnalysisTaskSE(name), 
   fAcceptTPC(kFALSE),
   fVz(1E20),
   fOutputList(0),
   fHEvent(0), 
   fTrackCuts(0)
{
//   
// Constructor
//

   // Define input and output slots here
   // Input slot #0 works with a TChain
   DefineInput(0, TChain::Class());
   
   // Output slot #0 id reserved by the base class for AOD
   // Output slot #1 writes into a TH1 container
   DefineOutput(1, TList::Class());
   
   // setup to NULL all remaining histos
   Int_t i;
   for (i = 0; i < kTypes; i++) fCheck[i] = 0;
}

//________________________________________________________________________
void AliAnalysisTaskEventCount::UserCreateOutputObjects()
{
//
// Create histograms (called once)
//

   // standard objects
   fESDpid = new AliESDpid();
   fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
   fTrackCuts->SetPtRange(0.1, 1E20);
   fTrackCuts->SetEtaRange(-0.8, 0.8);

   // create output list
   fOutputList = new TList();
   fOutputList->SetOwner();
   
   // create histograms
   fHEvent = new TH1I("EventStat", "Event statistics", kTypes, 0, kTypes);
   for (Int_t i = 0; i < kTypes; i++) {
      fHEvent->GetXaxis()->SetBinLabel(i + 1, Title((EType)i));
   }
   
   // add histograms to list
   fOutputList->Add(fHEvent);
   
   PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskEventCount::UserExec(Option_t *)
{
//
// Execute computations
//

   // assign first checks
   fCheck[kAll]               = kTRUE;
   fCheck[kPassPhysicsSel]    = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
   fCheck[kGoodPrimaryVertex] = kFALSE;
   fCheck[kHas1TrackAny]      = kFALSE;
   fCheck[kHas1TrackQuality]  = kFALSE;
   
   // try to cast input event to ESD and loop on tracks
   AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
   if (!esd) return;
   
   // check primary vertex
   const AliESDVertex *v0 = esd->GetPrimaryVertexTracks();
   if (v0) {
      if (v0->GetNContributors() > 0) 
         fCheck[kGoodPrimaryVertex] = kTRUE;
      else {
         v0 = fESD->GetPrimaryVertexSPD();
         if (v0) {
            if (v0->GetNContributors() > 0) 
               fCheck[kGoodPrimaryVertex] = kTRUE;
            else {
               if (fAcceptTPC) {
                  v0 = esd->GetPrimaryVertexTPC();
                  if (v0) {
                     if (v0->GetNContributors() > 0) fCheck[kGoodPrimaryVertex] = kTRUE;
                  }
               }
            }
         }
      }
   }
   
   // check track count
   if (esd->GetNumberOfTracks() > 0) fCheck[kHas1TrackAny] = kTRUE;
   if (fTrackCuts->CountAcceptedTracks(esd) > 0) fCheck[kHas1TrackQuality] = kTRUE;
   
   // fill histogram
   for (Int_t i = 0; i < kTypes; i++) {
      if (fCheck[i]) fHEvent->Fill((Double_t)i + 0.01);
   }
   
   PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskEventCount::Terminate(Option_t *)
{
//
// Termination
// Quickly check that the output list is there
//
  
   fOutputList = dynamic_cast<TList*>(GetOutputData(1));
   if (!fOutputList) {
      AliError("Output list not found!!!");
      return;
   }
}
