//
// Analysis task for 'mini' sub-package
// Contains all definitions needed for running an analysis:
// -- global event cut
// -- a list of track cuts (any number)
// -- definitions of output histograms
// -- values to be computed.
// Each one must be defined using the "CREATE" methods, which
// add directly a new element in the task collections, and don't
// need an external object to be passed to the task itself.
//

#include <Riostream.h>

#include <TH1.h>
#include <TList.h>
#include <TTree.h>

#include "AliLog.h"
#include "AliEventplane.h"
#include "AliMultiplicity.h"
#include "AliTriggerAnalysis.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliESDtrackCuts.h"
#include "AliESDUtils.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

#include "AliRsnCutSet.h"
#include "AliRsnMiniPair.h"
#include "AliRsnMiniEvent.h"
#include "AliRsnMiniParticle.h"

#include "AliRsnMiniMonitorTask.h"

ClassImp(AliRsnMiniMonitorTask)

//__________________________________________________________________________________________________
AliRsnMiniMonitorTask::AliRsnMiniMonitorTask() :
   AliAnalysisTaskSE(),
   fUseMC(kFALSE),
   fEvNum(0),
   fUseCentrality(kFALSE),
   fCentralityType("QUALITY"),
   fOutput(0x0),
   fHistograms("AliRsnMiniMonitor", 0),
   fEventCuts(0x0),
   fTrackCuts(0),
   fRsnEvent(),
   fBigOutput(kFALSE)
{
//
// Dummy constructor ALWAYS needed for I/O.
//
}

//__________________________________________________________________________________________________
AliRsnMiniMonitorTask::AliRsnMiniMonitorTask(const char *name, Bool_t useMC) :
   AliAnalysisTaskSE(name),
   fUseMC(useMC),
   fEvNum(0),
   fUseCentrality(kFALSE),
   fCentralityType("QUALITY"),
   fOutput(0x0),
   fHistograms("AliRsnMiniMonitor", 0),
   fEventCuts(0x0),
   fTrackCuts(0),
   fRsnEvent(),
   fBigOutput(kFALSE)
{
//
// Default constructor.
// Define input and output slots here (never in the dummy constructor)
// Input slot #0 works with a TChain - it is connected to the default input container
// Output slot #1 writes into a TH1 container
//

   DefineOutput(1, TList::Class());
}

//__________________________________________________________________________________________________
AliRsnMiniMonitorTask::AliRsnMiniMonitorTask(const AliRsnMiniMonitorTask& copy) :
   AliAnalysisTaskSE(copy),
   fUseMC(copy.fUseMC),
   fEvNum(0),
   fUseCentrality(copy.fUseCentrality),
   fCentralityType(copy.fCentralityType),
   fOutput(0x0),
   fHistograms(copy.fHistograms),
   fEventCuts(copy.fEventCuts),
   fTrackCuts(copy.fTrackCuts),
   fRsnEvent(),
   fBigOutput(copy.fBigOutput)
{
//
// Copy constructor.
// Implemented as requested by C++ standards.
// Can be used in PROOF and by plugins.
//
}

//__________________________________________________________________________________________________
AliRsnMiniMonitorTask& AliRsnMiniMonitorTask::operator=(const AliRsnMiniMonitorTask& copy)
{
//
// Assignment operator.
// Implemented as requested by C++ standards.
// Can be used in PROOF and by plugins.
//

   AliAnalysisTaskSE::operator=(copy);
   
   fUseMC = copy.fUseMC;
   fUseCentrality = copy.fUseCentrality;
   fCentralityType = copy.fCentralityType;
   fHistograms = copy.fHistograms;
   fEventCuts = copy.fEventCuts;
   fTrackCuts = copy.fTrackCuts;
   fBigOutput = copy.fBigOutput;
   
   return (*this);
}

//__________________________________________________________________________________________________
AliRsnMiniMonitorTask::~AliRsnMiniMonitorTask()
{
//
// Destructor. 
// Clean-up the output list, but not the histograms that are put inside
// (the list is owner and will clean-up these histograms). Protect in PROOF case.
//


   if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
   }
}

//__________________________________________________________________________________________________
Int_t AliRsnMiniMonitorTask::AddTrackCuts(AliRsnCutSet *cuts)
{
//
// Add a new cut set for a new criterion for track selection.
// A user can add as many as he wants, and each one corresponds
// to one of the available bits in the AliRsnMiniParticle mask.
// The only check is the following: if a cut set with the same name
// as the argument is there, this is not added.
// Return value is the array position of this set.
//

   TObject *obj = fTrackCuts.FindObject(cuts->GetName());
   
   if (obj) {
      AliInfo(Form("A cut set named '%s' already exists", cuts->GetName()));
      return fTrackCuts.IndexOf(obj);
   } else {
      fTrackCuts.AddLast(cuts);
      return fTrackCuts.IndexOf(cuts);
   }
}

//__________________________________________________________________________________________________
void AliRsnMiniMonitorTask::UserCreateOutputObjects()
{
//
// Initialization of outputs.
// This is called once per worker node.
//

   // reset counter
   fEvNum = 0;
   
   // message
   AliInfo(Form("Selected event characterization: %s (%s)", (fUseCentrality ? "centrality" : "multiplicity"), fCentralityType.Data()));
   
   // create list and set it as owner of its content (MANDATORY)
   if (fBigOutput) OpenFile(1);
   fOutput = new TList();
   fOutput->SetOwner();
   
   // create one histogram per each stored definition (event histograms)
   Int_t i, ndef = fHistograms.GetEntries();
   AliRsnMiniMonitor *def = 0x0;
   for (i = 0; i < ndef; i++) {
      def = (AliRsnMiniMonitor*)fHistograms[i];
      if (!def) continue;
      if (!def->Init(GetName(), fOutput)) {
         AliError(Form("Def '%s': failed initialization", def->GetName()));
         continue;
      }
   }
   
   // post data for ALL output slots >0 here, to get at least an empty histogram
   PostData(1, fOutput);
}

//__________________________________________________________________________________________________
void AliRsnMiniMonitorTask::UserExec(Option_t *)
{
//
// Computation loop.
// In this case, it checks if the event is acceptable, and eventually
// creates the corresponding mini-event and stores it in the buffer.
// The real histogram filling is done at the end, in "FinishTaskOutput".
//

   // event counter
   fEvNum++;
   
   // check current event
   Char_t check = CheckCurrentEvent();
   if (!check) return;
   
   // setup PID response
   AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager(); 
	AliInputEventHandler *inputHandler = (AliInputEventHandler*)man->GetInputEventHandler(); 
	fRsnEvent.SetPIDResponse(inputHandler->GetPIDResponse());
   
   // loop on monitors and fill them
   Int_t it, icut, nTracks = fRsnEvent.GetAbsoluteSum(), nCuts = fTrackCuts.GetEntriesFast();
   AliRsnDaughter cursor;
   AliRsnCutSet *cut = 0x0;
   AliRsnMiniMonitor *mon = 0x0;
   TObjArrayIter next(&fHistograms);
   for (it = 0; it < nTracks; it++) {
      fRsnEvent.SetDaughter(cursor, it, fUseMC);
      next.Reset();
      while ( (mon = (AliRsnMiniMonitor*)next()) ) {
         icut = mon->GetCutID();
         if (icut >= 0 && icut < nCuts) {
            cut = (AliRsnCutSet*)fTrackCuts[icut];
            if (!cut) {
               AliError("Cut not found");
               continue;
            }
            if (!cut->IsSelected(&cursor)) continue;
         }
         mon->Fill(&cursor, &fRsnEvent);
      }
   }
      
   // post data for computed stuff
   PostData(1, fOutput);
}

//__________________________________________________________________________________________________
void AliRsnMiniMonitorTask::Terminate(Option_t *)
{
//
// Draw result to screen, or perform fitting, normalizations
// Called once at the end of the query
//

   fOutput = dynamic_cast<TList*>(GetOutputData(1));
   if (!fOutput) { 
      AliError("Could not retrieve TList fOutput"); 
      return; 
   }
}

//__________________________________________________________________________________________________
Char_t AliRsnMiniMonitorTask::CheckCurrentEvent()
{
//
// This method checks if current event is OK for analysis.
// In case it is, the pointers of the local AliRsnEvent data member
// will point to it, in order to allow cut checking, otherwise the
// function exits with a failure message.
// ---
// ESD events must pass the physics selection, AOD are supposed to do.
// ---
// While checking the event, a histogram is filled to count the number
// of CINT1B, V0AND and CANDLE events, which are needed for normalization
// ---
// Return values can be:
//    -- 'E' if the event is accepted and is ESD
//    -- 'A' if the event is accepted and is AOD
//    --  0  if the event is not accepted
//

   // string to sum messages
   TString msg("");
   
   // check input type
   // exit points are provided in all cases an event is bad
   // if this block is passed, an event can be rejected only
   // if it does not pass the set of event cuts defined in the task
   Char_t output = 0;
   Bool_t isSelected;
   if (fInputEvent->InheritsFrom(AliESDEvent::Class())) {
      // type ESD
      output = 'E';
      // ESD specific check: Physics Selection
      // --> if this is failed, the event is rejected
      isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
      if (!isSelected) {
         AliDebugClass(2, "Event does not pass physics selections");
         fRsnEvent.SetRef(0x0);
         fRsnEvent.SetRefMC(0x0);
         return 0;
      }
      // set reference to input
      fRsnEvent.SetRef(fInputEvent);
      // add MC if requested and available
      if (fUseMC) {
         if (fMCEvent) 
            fRsnEvent.SetRefMC(fMCEvent);
         else {
            AliWarning("MC event requested but not available");
            fRsnEvent.SetRefMC(0x0);
         }
      }
   } else if (fInputEvent->InheritsFrom(AliAODEvent::Class())) {
      // type AOD
      output = 'A';
      // set reference to input
      fRsnEvent.SetRef(fInputEvent);
      // add MC if requested and available (it is in the same object)
      if (fUseMC) {
         fRsnEvent.SetRefMC(fInputEvent);
         if (!fRsnEvent.GetAODList()) {
            AliWarning("MC event requested but not available");
            fRsnEvent.SetRefMC(0x0);
         }
      }
   } else {
      AliError(Form("Bad input event class: %s", fInputEvent->ClassName()));
      // reset pointers in local AliRsnEvent object
      fRsnEvent.SetRef(0x0);
      fRsnEvent.SetRefMC(0x0);
      return 0;
   }
   
   // if event cuts are defined, they are checked here
   // final decision on the event depends on this
   isSelected = kTRUE;
   if (fEventCuts) {
      if (!fEventCuts->IsSelected(&fRsnEvent)) {
         msg += " -- Local cuts = REJECTED";
         isSelected = kFALSE;
      } else {
         msg += " -- Local cuts = ACCEPTED";
         isSelected = kTRUE;
      }
   } else {
      msg += " -- Local cuts = NONE";
      isSelected = kTRUE;
   }
   
   // if the above exit point is not taken, the event is accepted
   AliDebugClass(2, Form("Stats for event %d: %s", fEvNum, msg.Data()));
   if (isSelected) {
      return output;
   } else {
      return 0;
   }
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniMonitorTask::ComputeCentrality(Bool_t isESD)
{
//
// Computes event centrality/multiplicity according to the criterion defined
// by two elements: (1) choice between multiplicity and centrality and
// (2) the string defining what criterion must be used for specific computation.
//

   if (fUseCentrality) {
      AliCentrality *centrality = fInputEvent->GetCentrality();
      if (!centrality) {
         AliError("Cannot compute centrality!");
         return -1.0;
      }
      return centrality->GetCentralityPercentile(fCentralityType.Data());
   } else {
      if (!fCentralityType.CompareTo("TRACKS"))
         return fInputEvent->GetNumberOfTracks();
      else if (!fCentralityType.CompareTo("QUALITY"))
         if (isESD)
            return AliESDtrackCuts::GetReferenceMultiplicity((AliESDEvent*)fInputEvent, kTRUE);
         else {
            Double_t count = 0.;
            Int_t iTrack, ntracksLoop = fInputEvent->GetNumberOfTracks();
            for (iTrack = 0; iTrack < ntracksLoop; iTrack++) {    
               AliVTrack   *track = (AliVTrack*)fInputEvent->GetTrack(iTrack);
               AliAODTrack *aodt  = dynamic_cast<AliAODTrack*>(track);
               if (!aodt) continue;
               if (!aodt->TestFilterBit(5)) continue;
               count++;
            }
            return count;
         }
      else if (!fCentralityType.CompareTo("TRACKLETS")) {
         if (isESD) {
            const AliMultiplicity *mult = ((AliESDEvent*)fInputEvent)->GetMultiplicity();
            Float_t nClusters[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
            for(Int_t ilay = 0; ilay < 6; ilay++) nClusters[ilay] = (Float_t)mult->GetNumberOfITSClusters(ilay);
            return AliESDUtils::GetCorrSPD2(nClusters[1], fInputEvent->GetPrimaryVertex()->GetZ());
         } else {
            AliWarning("Cannot compute multiplicity with SPD tracklets from AOD");
            return 1E20;
         }
      } else {
         AliError(Form("String '%s' does not define a possible multiplicity/centrality computation", fCentralityType.Data()));
         return -1.0;
      }
   }
}
