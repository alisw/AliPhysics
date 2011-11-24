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
#include <TStopwatch.h>

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

#include "AliRsnMiniAnalysisTask.h"


ClassImp(AliRsnMiniAnalysisTask)

//__________________________________________________________________________________________________
AliRsnMiniAnalysisTask::AliRsnMiniAnalysisTask() :
   AliAnalysisTaskSE(),
   fUseMC(kFALSE),
   fEvNum(0),
   fUseCentrality(kFALSE),
   fCentralityType("QUALITY"),
   fContinuousMix(kTRUE),
   fNMix(0),
   fMaxDiffMult(10),
   fMaxDiffVz(1.0),
   fMaxDiffAngle(1E20),
   fOutput(0x0),
   fHistograms("AliRsnMiniOutput", 0),
   fValues("AliRsnMiniValue", 0),
   fHEventStat(0x0),
   fEventCuts(0x0),
   fTrackCuts(0),
   fRsnEvent(),
   fEvBuffer(0x0),
   fTriggerAna(0x0),
   fESDtrackCuts(0x0),
   fMiniEvent(0x0),
   fBigOutput(kFALSE)
{
//
// Dummy constructor ALWAYS needed for I/O.
//
}

//__________________________________________________________________________________________________
AliRsnMiniAnalysisTask::AliRsnMiniAnalysisTask(const char *name, Bool_t useMC) :
   AliAnalysisTaskSE(name),
   fUseMC(useMC),
   fEvNum(0),
   fUseCentrality(kFALSE),
   fCentralityType("QUALITY"),
   fContinuousMix(kTRUE),
   fNMix(0),
   fMaxDiffMult(10),
   fMaxDiffVz(1.0),
   fMaxDiffAngle(1E20),
   fOutput(0x0),
   fHistograms("AliRsnMiniOutput", 0),
   fValues("AliRsnMiniValue", 0),
   fHEventStat(0x0),
   fEventCuts(0x0),
   fTrackCuts(0),
   fRsnEvent(),
   fEvBuffer(0x0),
   fTriggerAna(0x0),
   fESDtrackCuts(0x0),
   fMiniEvent(0x0),
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
AliRsnMiniAnalysisTask::AliRsnMiniAnalysisTask(const AliRsnMiniAnalysisTask &copy) :
   AliAnalysisTaskSE(copy),
   fUseMC(copy.fUseMC),
   fEvNum(0),
   fUseCentrality(copy.fUseCentrality),
   fCentralityType(copy.fCentralityType),
   fContinuousMix(copy.fContinuousMix),
   fNMix(copy.fNMix),
   fMaxDiffMult(copy.fMaxDiffMult),
   fMaxDiffVz(copy.fMaxDiffVz),
   fMaxDiffAngle(copy.fMaxDiffAngle),
   fOutput(0x0),
   fHistograms(copy.fHistograms),
   fValues(copy.fValues),
   fHEventStat(0x0),
   fEventCuts(copy.fEventCuts),
   fTrackCuts(copy.fTrackCuts),
   fRsnEvent(),
   fEvBuffer(0x0),
   fTriggerAna(copy.fTriggerAna),
   fESDtrackCuts(copy.fESDtrackCuts),
   fMiniEvent(0x0),
   fBigOutput(copy.fBigOutput)
{
//
// Copy constructor.
// Implemented as requested by C++ standards.
// Can be used in PROOF and by plugins.
//
}

//__________________________________________________________________________________________________
AliRsnMiniAnalysisTask &AliRsnMiniAnalysisTask::operator=(const AliRsnMiniAnalysisTask &copy)
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
   fContinuousMix = copy.fContinuousMix;
   fNMix = copy.fNMix;
   fMaxDiffMult = copy.fMaxDiffMult;
   fMaxDiffVz = copy.fMaxDiffVz;
   fMaxDiffAngle = copy.fMaxDiffAngle;
   fHistograms = copy.fHistograms;
   fValues = copy.fValues;
   fEventCuts = copy.fEventCuts;
   fTrackCuts = copy.fTrackCuts;
   fTriggerAna = copy.fTriggerAna;
   fESDtrackCuts = copy.fESDtrackCuts;
   fBigOutput = copy.fBigOutput;

   return (*this);
}

//__________________________________________________________________________________________________
AliRsnMiniAnalysisTask::~AliRsnMiniAnalysisTask()
{
//
// Destructor.
// Clean-up the output list, but not the histograms that are put inside
// (the list is owner and will clean-up these histograms). Protect in PROOF case.
//


   if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
      delete fEvBuffer;
   }
}

//__________________________________________________________________________________________________
Int_t AliRsnMiniAnalysisTask::AddTrackCuts(AliRsnCutSet *cuts)
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
void AliRsnMiniAnalysisTask::UserCreateOutputObjects()
{
//
// Initialization of outputs.
// This is called once per worker node.
//

   // reset counter
   fEvNum = -1;

   // message
   AliInfo(Form("Selected event characterization: %s (%s)", (fUseCentrality ? "centrality" : "multiplicity"), fCentralityType.Data()));

   // initialize trigger analysis
   if (fTriggerAna) delete fTriggerAna;
   fTriggerAna = new AliTriggerAnalysis;

   // initialize ESD quality cuts
   if (fESDtrackCuts) delete fESDtrackCuts;
   fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();

   // create list and set it as owner of its content (MANDATORY)
   if (fBigOutput) OpenFile(1);
   fOutput = new TList();
   fOutput->SetOwner();

   // initialize event statistics counter
   fHEventStat = new TH1F("hEventStat", "Event statistics", 4, 0.0, 4.0);
   fHEventStat->GetXaxis()->SetBinLabel(1, "CINT1B");
   fHEventStat->GetXaxis()->SetBinLabel(2, "V0AND");
   fHEventStat->GetXaxis()->SetBinLabel(3, "Candle");
   fHEventStat->GetXaxis()->SetBinLabel(4, "Accepted");
   fOutput->Add(fHEventStat);

   // create temporary tree for filtered events
   if (fMiniEvent) delete fMiniEvent;
   fEvBuffer = new TTree("EventBuffer", "Temporary buffer for mini events");
   fEvBuffer->Branch("events", "AliRsnMiniEvent", &fMiniEvent);

   // create one histogram per each stored definition (event histograms)
   Int_t i, ndef = fHistograms.GetEntries();
   AliRsnMiniOutput *def = 0x0;
   for (i = 0; i < ndef; i++) {
      def = (AliRsnMiniOutput *)fHistograms[i];
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
void AliRsnMiniAnalysisTask::UserExec(Option_t *)
{
//
// Computation loop.
// In this case, it checks if the event is acceptable, and eventually
// creates the corresponding mini-event and stores it in the buffer.
// The real histogram filling is done at the end, in "FinishTaskOutput".
//

   // increment event counter
   fEvNum++;

   // check current event
   Char_t check = CheckCurrentEvent();
   if (!check) return;

   // setup PID response
   AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler *inputHandler = (AliInputEventHandler *)man->GetInputEventHandler();
   fRsnEvent.SetPIDResponse(inputHandler->GetPIDResponse());

   // fill a mini-event from current
   // and skip this event if no tracks were accepted
   FillMiniEvent(check);

   // fill MC based histograms on mothers,
   // which do need the original event
   if (fUseMC) {
      if (fRsnEvent.IsESD() && fMCEvent)
         FillTrueMotherESD(fMiniEvent);
      else if (fRsnEvent.IsAOD() && fRsnEvent.GetAODList())
         FillTrueMotherAOD(fMiniEvent);
   }

   // if the event is not empty, store it
   if (fMiniEvent->IsEmpty()) {
      AliDebugClass(2, Form("Rejecting empty event #%d", fEvNum));
   } else {
      Int_t id = fEvBuffer->GetEntries();
      AliDebugClass(2, Form("Adding event #%d with ID = %d", fEvNum, id));
      fMiniEvent->ID() = id;
      fEvBuffer->Fill();
   }

   // post data for computed stuff
   PostData(1, fOutput);
}

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTask::FinishTaskOutput()
{
//
// This function is called at the end of the loop on available events,
// and then the buffer will be full with all the corresponding mini-events,
// each one containing all tracks selected by each of the available track cuts.
// Here a loop is done on each of these events, and both single-event and mixing are computed
//

   // security code: reassign the buffer to the mini-event cursor
   fEvBuffer->SetBranchAddress("events", &fMiniEvent);
   TStopwatch timer;
   // prepare variables
   Int_t ievt, nEvents = (Int_t)fEvBuffer->GetEntries();
   Int_t idef, nDefs   = fHistograms.GetEntries();
   Int_t imix, iloop, ifill;
   AliRsnMiniOutput *def = 0x0;
   AliRsnMiniOutput::EComputation compType;

   Int_t printNum = 0;
   if (nEvents>1e5) printNum=nEvents/100;
   if (nEvents>1e4) printNum=nEvents/10;
   
   // loop on events, and for each one fill all outputs
   // using the appropriate procedure depending on its type
   // only mother-related histograms are filled in UserExec,
   // since they require direct access to MC event
   timer.Start();
   for (ievt = 0; ievt < nEvents; ievt++) {
      // get next entry
      fEvBuffer->GetEntry(ievt);
      if (printNum&&(ievt%printNum==0)) {
         AliInfo(Form("[%s] Std.Event %d/%d",GetName(), ievt,nEvents));
         timer.Stop(); timer.Print();fflush(stdout); timer.Start(kFALSE);
      }
      // fill
      for (idef = 0; idef < nDefs; idef++) {
         def = (AliRsnMiniOutput *)fHistograms[idef];
         if (!def) continue;
         compType = def->GetComputation();
         // execute computation in the appropriate way
         switch (compType) {
            case AliRsnMiniOutput::kEventOnly:
               //AliDebugClass(1, Form("Event %d, def '%s': event-value histogram filling", ievt, def->GetName()));
               ifill = 1;
               def->FillEvent(fMiniEvent, &fValues);
               break;
            case AliRsnMiniOutput::kTruePair:
               //AliDebugClass(1, Form("Event %d, def '%s': true-pair histogram filling", ievt, def->GetName()));
               ifill = def->FillPair(fMiniEvent, fMiniEvent, &fValues);
               break;
            case AliRsnMiniOutput::kTrackPair:
               //AliDebugClass(1, Form("Event %d, def '%s': pair-value histogram filling", ievt, def->GetName()));
               ifill = def->FillPair(fMiniEvent, fMiniEvent, &fValues);
               break;
            case AliRsnMiniOutput::kTrackPairRotated1:
               //AliDebugClass(1, Form("Event %d, def '%s': rotated (1) background histogram filling", ievt, def->GetName()));
               ifill = def->FillPair(fMiniEvent, fMiniEvent, &fValues);
               break;
            case AliRsnMiniOutput::kTrackPairRotated2:
               //AliDebugClass(1, Form("Event %d, def '%s': rotated (2) background histogram filling", ievt, def->GetName()));
               ifill = def->FillPair(fMiniEvent, fMiniEvent, &fValues);
               break;
            default:
               // other kinds are processed elsewhere
               ifill = 0;
               AliDebugClass(2, Form("Computation = %d", (Int_t)compType));
         }
         // message
         AliDebugClass(1, Form("Event %6d: def = '%15s' -- fills = %5d", ievt, def->GetName(), ifill));
      }
   }

   // if no mixing is required, stop here and post the output
   if (fNMix < 1) {
      AliDebugClass(2, "Stopping here, since no mixing is required");
      PostData(1, fOutput);
      return;
   }

   // initialize mixing counter
   Int_t    nmatched[nEvents];
   TString *smatched = new TString[nEvents];
   for (ievt = 0; ievt < nEvents; ievt++) {
      smatched[ievt] = "|";
      nmatched[ievt] = 0;
   }
   
   
   AliInfo(Form("[%s] Std.Event %d/%d",GetName(), nEvents,nEvents));
   timer.Stop(); timer.Print(); timer.Start();fflush(stdout);

   // search for good matchings
   for (ievt = 0; ievt < nEvents; ievt++) {
      if (printNum&&(ievt%printNum==0)) {
         AliInfo(Form("[%s] EventMixing searching %d/%d",GetName(),ievt,nEvents));
         timer.Stop(); timer.Print(); timer.Start(kFALSE); fflush(stdout);
      }
      if (nmatched[ievt] >= fNMix) continue;
      fEvBuffer->GetEntry(ievt);
      AliRsnMiniEvent evMain(*fMiniEvent);
      for (iloop = 1; iloop < nEvents; iloop++) {
         imix = ievt + iloop;
         if (imix >= nEvents) imix -= nEvents;
         if (imix == ievt) continue;
         // text next entry
         fEvBuffer->GetEntry(imix);
         // skip if events are not matched
         if (!EventsMatch(&evMain, fMiniEvent)) continue;
         // check that the array of good matches for mixed does not already contain main event
         if (smatched[imix].Contains(Form("|%d|", ievt))) continue;
         // check that the found good events has not enough matches already
         if (nmatched[imix] >= fNMix) continue;
         // add new mixing candidate
         smatched[ievt].Append(Form("%d|", imix));
         nmatched[ievt]++;
         nmatched[imix]++;
         if (nmatched[ievt] >= fNMix) break;
      }
      AliDebugClass(1, Form("Matches for event %5d = %d [%s] (missing are declared above)", evMain.ID(), nmatched[ievt], smatched[ievt].Data()));
   }
   
   AliInfo(Form("[%s] EventMixing searching %d/%d",GetName(),nEvents,nEvents));
   timer.Stop(); timer.Print();fflush(stdout); timer.Start();

   // perform mixing
   TObjArray *list = 0x0;
   TObjString *os = 0x0;
   for (ievt = 0; ievt < nEvents; ievt++) {
      if (printNum&&(ievt%printNum==0)) {
         AliInfo(Form("[%s] EventMixing %d/%d",GetName(),ievt,nEvents));
         timer.Stop(); timer.Print(); timer.Start(kFALSE); fflush(stdout);
      }
      ifill = 0;
      fEvBuffer->GetEntry(ievt);
      AliRsnMiniEvent evMain(*fMiniEvent);
      list = smatched[ievt].Tokenize("|");
      TObjArrayIter next(list);
      while ( (os = (TObjString *)next()) ) {
         imix = os->GetString().Atoi();
         fEvBuffer->GetEntry(imix);
         for (idef = 0; idef < nDefs; idef++) {
            def = (AliRsnMiniOutput *)fHistograms[idef];
            if (!def) continue;
            if (!def->IsTrackPairMix()) continue;
            ifill += def->FillPair(&evMain, fMiniEvent, &fValues, kTRUE);
            if (!def->IsSymmetric()) {
               AliDebugClass(2, "Reflecting non symmetric pair");
               ifill += def->FillPair(fMiniEvent, &evMain, &fValues, kFALSE);
            }
         }
      }
      delete list;
   }

   delete [] smatched;

   AliInfo(Form("[%s] EventMixing %d/%d",GetName(),nEvents,nEvents));
   timer.Stop();timer.Print();fflush(stdout);

   /*
   OLD
   ifill = 0;
   for (iloop = 1; iloop < nEvents; iloop++) {
      imix = ievt + iloop;
      // restart from beginning if reached last event
      if (imix >= nEvents) imix -= nEvents;
      // avoid to mix an event with itself
      if (imix == ievt) continue;
      // skip all events already mixed enough times
      if (fNMixed[ievt] >= fNMix) break;
      if (fNMixed[imix] >= fNMix) continue;
      fEvBuffer->GetEntry(imix);
      // skip if events are not matched
      if (TMath::Abs(evMain.Vz()    - fMiniEvent->Vz()   ) > fMaxDiffVz   ) continue;
      if (TMath::Abs(evMain.Mult()  - fMiniEvent->Mult() ) > fMaxDiffMult ) continue;
      if (TMath::Abs(evMain.Angle() - fMiniEvent->Angle()) > fMaxDiffAngle) continue;
      // found a match: increment counter for both events
      AliDebugClass(1, Form("Event %d, def '%s': event mixing (%d with %d)", ievt, def->GetName(), ievt, imix));
      fNMixed[ievt]++;
      fNMixed[imix]++;
      // process mixing
      ifill += def->FillPair(&evMain, fMiniEvent, &fValues);
      // stop if mixed enough times
      if (fNMixed[ievt] >= fNMix) break;
   }
   break;
   // print number of mixings done with each event
   for (ievt = 0; ievt < nEvents; ievt++) {
      AliDebugClass(2, Form("Event %6d: mixed %2d times", ievt, fNMixed[ievt]));
   }
   */

   // post computed data
   PostData(1, fOutput);
}

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTask::Terminate(Option_t *)
{
//
// Draw result to screen, or perform fitting, normalizations
// Called once at the end of the query
//

   fOutput = dynamic_cast<TList *>(GetOutputData(1));
   if (!fOutput) {
      AliError("Could not retrieve TList fOutput");
      return;
   }
}

//__________________________________________________________________________________________________
Char_t AliRsnMiniAnalysisTask::CheckCurrentEvent()
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
      isSelected = (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
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

   // fill counter of accepted events
   fHEventStat->Fill(0.1);

   // check if it is V0AND
   // --> uses a cast to AliESDEvent even if the input is an AliAODEvent
   Bool_t v0A = fTriggerAna->IsOfflineTriggerFired((AliESDEvent *)fInputEvent, AliTriggerAnalysis::kV0A);
   Bool_t v0C = fTriggerAna->IsOfflineTriggerFired((AliESDEvent *)fInputEvent, AliTriggerAnalysis::kV0C);
   if (v0A && v0C) {
      msg += " -- VOAND = YES";
      fHEventStat->Fill(1.1);
   } else {
      msg += " -- VOAND = NO ";
   }

   // check candle
   // --> requires at least one good quality track with Pt > 0.5 and |eta| <= 0.8
   Int_t iTrack, ntracksLoop = fInputEvent->GetNumberOfTracks();
   Bool_t candle = kFALSE;
   for (iTrack = 0; iTrack < ntracksLoop; iTrack++) {
      AliVTrack   *track = (AliVTrack *)fInputEvent->GetTrack(iTrack);
      AliESDtrack *esdt  = dynamic_cast<AliESDtrack *>(track);
      AliAODTrack *aodt  = dynamic_cast<AliAODTrack *>(track);
      if (track->Pt() < 0.5) continue;
      if(TMath::Abs(track->Eta()) > 0.8) continue;
      if (esdt) if (!fESDtrackCuts->AcceptTrack(esdt)) continue;
      if (aodt) if (!aodt->TestFilterBit(5)) continue;
      candle = kTRUE;
      break;
   }
   if (candle) {
      msg += " -- CANDLE = YES";
      fHEventStat->Fill(2.1);
   } else {
      msg += " -- CANDLE = NO ";
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
   AliDebugClass(2, Form("Stats: %s", msg.Data()));
   if (isSelected) {
      fHEventStat->Fill(3.1);
      return output;
   } else {
      return 0;
   }
}

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTask::FillMiniEvent(Char_t evType)
{
//
// Refresh cursor mini-event data member to fill with current event.
// Returns the total number of tracks selected.
//

   // assign event-related values
   if (fMiniEvent) delete fMiniEvent;
   fMiniEvent = new AliRsnMiniEvent;
   fMiniEvent->Vz()    = fInputEvent->GetPrimaryVertex()->GetZ();
   fMiniEvent->Angle() = ComputeAngle();
   fMiniEvent->Mult()  = ComputeCentrality((evType == 'E'));
   AliDebugClass(2, Form("Event %d: type = %c -- vz = %f -- mult = %f -- angle = %f", fEvNum, evType, fMiniEvent->Vz(), fMiniEvent->Mult(), fMiniEvent->Angle()));

   // loop on daughters and assign track-related values
   Int_t ic, ncuts = fTrackCuts.GetEntries();
   Int_t ip, npart = fRsnEvent.GetAbsoluteSum();
   Int_t npos = 0, nneg = 0, nneu = 0;
   AliRsnDaughter cursor;
   AliRsnMiniParticle miniParticle;
   for (ip = 0; ip < npart; ip++) {
      // point cursor to next particle
      fRsnEvent.SetDaughter(cursor, ip);
      // copy momentum and MC info if present
      miniParticle.CopyDaughter(&cursor);
      miniParticle.Index() = ip;
      // switch on the bits corresponding to passed cuts
      for (ic = 0; ic < ncuts; ic++) {
         AliRsnCutSet *cuts = (AliRsnCutSet *)fTrackCuts[ic];
         if (cuts->IsSelected(&cursor)) miniParticle.SetCutBit(ic);
      }
      // if a track passes at least one track cut, it is added to the pool
      if (miniParticle.CutBits()) {
         fMiniEvent->AddParticle(miniParticle);
         if (miniParticle.Charge() == '+') npos++;
         else if (miniParticle.Charge() == '-') nneg++;
         else nneu++;
      }
   }

   // get number of accepted tracks
   AliDebugClass(1, Form("Event %6d: total = %5d, accepted = %4d (pos %4d, neg %4d, neu %4d)", fEvNum, npart, (Int_t)fMiniEvent->Particles().GetEntriesFast(), npos, nneg, nneu));
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniAnalysisTask::ComputeAngle()
{
//
// Get the plane angle
//

   AliEventplane *plane = 0x0;

   if (fInputEvent->InheritsFrom(AliESDEvent::Class()))
      plane = fInputEvent->GetEventplane();
   else if (fInputEvent->InheritsFrom(AliAODEvent::Class())) {
      AliAODEvent *aodEvent = (AliAODEvent *)fInputEvent;
      plane = aodEvent->GetHeader()->GetEventplaneP();
   }

   if (plane)
      return plane->GetEventplane("Q");
   else {
      AliWarning("No event plane defined");
      return 1E20;
   }
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniAnalysisTask::ComputeCentrality(Bool_t isESD)
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
            return AliESDtrackCuts::GetReferenceMultiplicity((AliESDEvent *)fInputEvent, kTRUE);
         else {
            Double_t count = 0.;
            Int_t iTrack, ntracksLoop = fInputEvent->GetNumberOfTracks();
            for (iTrack = 0; iTrack < ntracksLoop; iTrack++) {
               AliVTrack   *track = (AliVTrack *)fInputEvent->GetTrack(iTrack);
               AliAODTrack *aodt  = dynamic_cast<AliAODTrack *>(track);
               if (!aodt) continue;
               if (!aodt->TestFilterBit(5)) continue;
               count++;
            }
            return count;
         }
      else if (!fCentralityType.CompareTo("TRACKLETS")) {
         if (isESD) {
            const AliMultiplicity *mult = ((AliESDEvent *)fInputEvent)->GetMultiplicity();
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

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTask::FillTrueMotherESD(AliRsnMiniEvent *miniEvent)
{
//
// Fills the histograms with true mother (ESD version)
//

   Bool_t okMatch;
   Int_t id, ndef = fHistograms.GetEntries();
   Int_t ip, label1, label2, npart = fMCEvent->GetNumberOfTracks();
   static AliRsnMiniPair miniPair;
   AliMCParticle *daughter1, *daughter2;
   TLorentzVector p1, p2;
   AliRsnMiniOutput *def = 0x0;

   for (id = 0; id < ndef; id++) {
      def = (AliRsnMiniOutput *)fHistograms[id];
      if (!def) continue;
      if (!def->IsMother()) continue;
      for (ip = 0; ip < npart; ip++) {
         AliMCParticle *part = (AliMCParticle *)fMCEvent->GetTrack(ip);
         if (part->Particle()->GetPdgCode() != def->GetMotherPDG()) continue;
         // check that daughters match expected species
         if (part->Particle()->GetNDaughters() < 2) continue;
         label1 = part->Particle()->GetDaughter(0);
         label2 = part->Particle()->GetDaughter(1);
         daughter1 = (AliMCParticle *)fMCEvent->GetTrack(label1);
         daughter2 = (AliMCParticle *)fMCEvent->GetTrack(label2);
         okMatch = kFALSE;
         if (TMath::Abs(daughter1->Particle()->GetPdgCode()) == def->GetPDG(0) && TMath::Abs(daughter2->Particle()->GetPdgCode()) == def->GetPDG(1)) {
            okMatch = kTRUE;
            p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(0));
            p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(1));
         } else if (TMath::Abs(daughter1->Particle()->GetPdgCode()) == def->GetPDG(1) && TMath::Abs(daughter2->Particle()->GetPdgCode()) == def->GetPDG(0)) {
            okMatch = kTRUE;
            p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(1));
            p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(0));
         }
         if (!okMatch) continue;
         // assign momenta to computation object
         miniPair.Sum(0) = miniPair.Sum(1) = (p1 + p2);
         miniPair.FillRef(def->GetMotherMass());
         // do computations
         def->FillMother(&miniPair, miniEvent, &fValues);
      }
   }
}

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTask::FillTrueMotherAOD(AliRsnMiniEvent *miniEvent)
{
//
// Fills the histograms with true mother (AOD version)
//

   Bool_t okMatch;
   TClonesArray *list = fRsnEvent.GetAODList();
   Int_t id, ndef = fHistograms.GetEntries();
   Int_t ip, label1, label2, npart = list->GetEntries();
   static AliRsnMiniPair miniPair;
   AliAODMCParticle *daughter1, *daughter2;
   TLorentzVector p1, p2;
   AliRsnMiniOutput *def = 0x0;

   for (id = 0; id < ndef; id++) {
      def = (AliRsnMiniOutput *)fHistograms[id];
      if (!def) continue;
      if (!def->IsMother()) continue;
      for (ip = 0; ip < npart; ip++) {
         AliAODMCParticle *part = (AliAODMCParticle *)list->At(ip);
         if (part->GetPdgCode() != def->GetMotherPDG()) continue;
         // check that daughters match expected species
         if (part->GetNDaughters() < 2) continue;
         label1 = part->GetDaughter(0);
         label2 = part->GetDaughter(1);
         daughter1 = (AliAODMCParticle *)list->At(label1);
         daughter2 = (AliAODMCParticle *)list->At(label2);
         okMatch = kFALSE;
         if (TMath::Abs(daughter1->GetPdgCode()) == def->GetPDG(0) && TMath::Abs(daughter2->GetPdgCode()) == def->GetPDG(1)) {
            okMatch = kTRUE;
            p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(0));
            p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(1));
         } else if (TMath::Abs(daughter1->GetPdgCode()) == def->GetPDG(1) && TMath::Abs(daughter2->GetPdgCode()) == def->GetPDG(0)) {
            okMatch = kTRUE;
            p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(1));
            p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(0));
         }
         if (!okMatch) continue;
         // assign momenta to computation object
         miniPair.Sum(0) = miniPair.Sum(1) = (p1 + p2);
         miniPair.FillRef(def->GetMotherMass());
         // do computations
         def->FillMother(&miniPair, miniEvent, &fValues);
      }
   }
}

//__________________________________________________________________________________________________
Bool_t AliRsnMiniAnalysisTask::EventsMatch(AliRsnMiniEvent *event1, AliRsnMiniEvent *event2)
{
//
// Check if two events are compatible.
// If the mixing is continuous, this is true if differences in vz, mult and angle are smaller than
// the specified values.
// If the mixing is binned, this is true if the events are in the same bin.
//

   if (!event1 || !event2) return kFALSE;
   Int_t ivz1, ivz2, imult1, imult2, iangle1, iangle2;
   Double_t dv, dm, da;

   if (fContinuousMix) {
      dv = TMath::Abs(event1->Vz()    - event2->Vz()   );
      dm = TMath::Abs(event1->Mult()  - event2->Mult() );
      da = TMath::Abs(event1->Angle() - event2->Angle());
      if (dv > fMaxDiffVz) {
         //AliDebugClass(2, Form("Events #%4d and #%4d don't match due to a too large diff in Vz = %f", event1->ID(), event2->ID(), dv));
         return kFALSE;
      }
      if (dm > fMaxDiffMult ) {
         //AliDebugClass(2, Form("Events #%4d and #%4d don't match due to a too large diff in Mult = %f", event1->ID(), event2->ID(), dm));
         return kFALSE;
      }
      if (da > fMaxDiffAngle) {
         //AliDebugClass(2, Form("Events #%4d and #%4d don't match due to a too large diff in Angle = %f", event1->ID(), event2->ID(), da));
         return kFALSE;
      }
      return kTRUE;
   } else {
      ivz1 = (Int_t)(event1->Vz() / fMaxDiffVz);
      ivz2 = (Int_t)(event2->Vz() / fMaxDiffVz);
      imult1 = (Int_t)(event1->Mult() / fMaxDiffMult);
      imult2 = (Int_t)(event2->Mult() / fMaxDiffMult);
      iangle1 = (Int_t)(event1->Angle() / fMaxDiffAngle);
      iangle2 = (Int_t)(event2->Angle() / fMaxDiffAngle);
      if (ivz1 != ivz2) return kFALSE;
      if (imult1 != imult2) return kFALSE;
      if (iangle1 != iangle2) return kFALSE;
      return kTRUE;
   }
}


