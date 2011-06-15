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

#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliESDtrackCuts.h"
#include "AliESDUtils.h"
#include "AliMultiplicity.h"
#include "AliInputEventHandler.h"
#include "AliTriggerAnalysis.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"

#include "AliRsnCutSet.h"
#include "AliRsnMiniEvent.h"
#include "AliRsnMiniParticle.h"
#include "AliRsnMiniPair.h"

#include "AliRsnMiniAnalysisTask.h"

ClassImp(AliRsnMiniAnalysisTask)

//__________________________________________________________________________________________________
AliRsnMiniAnalysisTask::AliRsnMiniAnalysisTask() :
   AliAnalysisTaskSE(),
   fUseMC(kFALSE),
   fEvNum(0),
   fUseCentrality(kFALSE),
   fCentralityType("QUALITY"),
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
   fNMixed(0),
   fTriggerAna(0x0),
   fESDtrackCuts(0x0)
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
   fNMixed(0),
   fTriggerAna(0x0),
   fESDtrackCuts(0x0)
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
AliRsnMiniAnalysisTask::AliRsnMiniAnalysisTask(const AliRsnMiniAnalysisTask& copy) :
   AliAnalysisTaskSE(copy),
   fUseMC(copy.fUseMC),
   fEvNum(0),
   fUseCentrality(copy.fUseCentrality),
   fCentralityType(copy.fCentralityType),
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
   fNMixed(0),
   fTriggerAna(copy.fTriggerAna),
   fESDtrackCuts(copy.fESDtrackCuts)
{
//
// Copy constructor.
// Implemented as requested by C++ standards.
// Can be used in PROOF and by plugins.
//
}

//__________________________________________________________________________________________________
AliRsnMiniAnalysisTask& AliRsnMiniAnalysisTask::operator=(const AliRsnMiniAnalysisTask& copy)
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
   fEvNum = 0;

   // initialize trigger analysis
   if (fTriggerAna) delete fTriggerAna;
   fTriggerAna = new AliTriggerAnalysis;
   
   // initialize ESD quality cuts
   if (fESDtrackCuts) delete fESDtrackCuts;
   fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();

   // create list and set it as owner of its content (MANDATORY)
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
   AliRsnMiniEvent *mini = 0x0;
   fEvBuffer = new TTree("EventBuffer", "Temporary buffer for mini events");
   fEvBuffer->Branch("events", "AliRsnMiniEvent", &mini);
   
   // create one histogram per each stored definition (event histograms)
   Int_t i, ndef = fHistograms.GetEntries();
   AliRsnMiniOutput *def = 0x0;
   for (i = 0; i < ndef; i++) {
      def = (AliRsnMiniOutput*)fHistograms[i];
      if (!def) continue;
      if (!def->Init(GetName(), fOutput)) {
         AliError(Form("Def '%s': failed initialization", def->GetName()));
         continue;
      }
   }
   
   // setup PID response
   AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager(); 
	AliInputEventHandler *inputHandler = (AliInputEventHandler*)man->GetInputEventHandler(); 
	fRsnEvent.SetPIDResponse(inputHandler->GetPIDResponse());
   
   // post data for ALL output slots >0 here, to get at least an empty histogram
   PostData(1, fOutput);
}

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTask::UserExec(Option_t *)
{
//
// Main computation loop.
// To be precise, this loop only checks if events are OK and eventually
// stores them in the temporary tree used for the analysis.
// The computation of correlations is done only at the end,
// in the automatic function which is called every time
//

   // event counter
   fEvNum++;
   
   // check current event
   Char_t check = CheckCurrentEvent();
   if (!check) return;
      
   // if the check is successful, the mini-event is created and stored
   AliRsnMiniEvent *miniEvent = 0x0;
   fEvBuffer->SetBranchAddress("events", &miniEvent);
   
   // assign event-related values
   miniEvent = new AliRsnMiniEvent;
   miniEvent->Vz() = fInputEvent->GetPrimaryVertex()->GetZ();
   miniEvent->Angle() = 0.0;
   miniEvent->Mult() = ComputeCentrality((check == 'E'));
   AliDebugClass(1, Form("Event %d: vz = %f -- mult = %f -- angle = %f", fEvNum, miniEvent->Vz(), miniEvent->Mult(), miniEvent->Angle()));
   
   // fill all histograms for events only
   Int_t id, ndef = fHistograms.GetEntries();
   AliRsnMiniOutput *def = 0x0;
   for (id = 0; id < ndef; id++) {
      def = (AliRsnMiniOutput*)fHistograms[id];
      if (!def) continue;
      if (!def->IsEventOnly()) continue;
      def->Fill(miniEvent, &fValues);
   }
   
   // fill all histograms for mother only (if MC is present)
   if (fUseMC) {
      if (fRsnEvent.IsESD() && fMCEvent)
         FillTrueMotherESD(miniEvent);
      else if (fRsnEvent.IsAOD() && fRsnEvent.GetAODList())
         FillTrueMotherAOD(miniEvent);
   }
   
   // loop on daughters
   // and store only those that pass at least one cut
   Int_t it, ic, ncuts = fTrackCuts.GetEntries(), ntot = fRsnEvent.GetAbsoluteSum();
   AliRsnDaughter cursor;
   AliRsnMiniParticle miniParticle;
   for (it = 0; it < ntot; it++) {
      fRsnEvent.SetDaughter(cursor, it);
      miniParticle.CopyDaughter(&cursor);
      for (ic = 0; ic < ncuts; ic++) {
         AliRsnCutSet *cuts = (AliRsnCutSet*)fTrackCuts[ic];
         if (cuts->IsSelected(&cursor)) miniParticle.SetCutBit(ic);
      }
      AliDebugClass(3, Form("-- Track %d: REC px = %f -- py = %f -- pz = %f", it, miniParticle.Px(0), miniParticle.Py(0), miniParticle.Pz(0)));
      AliDebugClass(3, Form("-- Track %d: SIM px = %f -- py = %f -- pz = %f", it, miniParticle.Px(1), miniParticle.Py(1), miniParticle.Pz(1)));
      AliDebugClass(3, Form("-- Track %d: Charge = %c -- PDG = %d -- mother = %d -- motherPDG = %d", it, miniParticle.Charge(), miniParticle.PDG(), miniParticle.Mother(), miniParticle.MotherPDG()));
      AliDebugClass(2, Form("-- Track %d: Cutbit = %u", it, miniParticle.CutBits()));
      if (miniParticle.CutBits()) {
         miniEvent->AddParticle(miniParticle);
      }
   }
   AliDebugClass(1, Form("Event %d: selected tracks = %d", fEvNum, miniEvent->Particles().GetEntriesFast()));
   
   // store event
   fEvBuffer->Fill();
   
   // process single event
   ProcessEvents(miniEvent);
   
   // post outputs
   PostData(1, fOutput);
}

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTask::FinishTaskOutput()
{
//
// At the end of execution, when the temporary tree is filled,
// perform mixing with all found events
//

   Int_t i1, i2, imix, nEvents = fEvBuffer->GetEntries();
   AliRsnMiniEvent *evMix = 0x0, evMain;
   fEvBuffer->SetBranchAddress("events", &evMix);
   
   // initialize mixing counter
   fNMixed.Set(nEvents);
   TString msg();
   
   // loop on events
   for (i1 = 0; i1 < nEvents; i1++) {
      if (fNMixed[i1] >= fNMix) continue;
      fEvBuffer->GetEntry(i1);
      evMain = (*evMix);
      for (i2 = 1; i2 < nEvents; i2++) {
         imix = i1 + i2;
         if (imix >= nEvents) imix -= nEvents;
         if (imix == i1) continue;
         if (fNMixed[i1] >= fNMix) break;
         if (fNMixed[imix] >= fNMix) continue;
         fEvBuffer->GetEntry(imix);
         // exit if events are not matched
         if (TMath::Abs(evMain.Vz() - evMix->Vz()) > fMaxDiffVz) continue;
         if (TMath::Abs(evMain.Mult() - evMix->Mult()) > fMaxDiffMult) continue;
         if (TMath::Abs(evMain.Angle() - evMix->Angle()) > fMaxDiffAngle) continue;
         // found a match: increment counter for both events
         fNMixed[i1]++;
         fNMixed[imix]++;
         ProcessEvents(&evMain, evMix);
         // if mixed enough times, stop
      }
   }
   
   // message
   for (Int_t ii = 0; ii < nEvents; ii++) cout << Form("Event #%6d mixed %3d times", ii, fNMixed[ii]) << endl;
   
   // post outputs
   PostData(1, fOutput);
}

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTask::Terminate(Option_t *)
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
      isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
      if (!isSelected) {
         AliDebugClass(1, "Event does not pass physics selections");
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
   Bool_t v0A = fTriggerAna->IsOfflineTriggerFired((AliESDEvent*)fInputEvent, AliTriggerAnalysis::kV0A);
   Bool_t v0C = fTriggerAna->IsOfflineTriggerFired((AliESDEvent*)fInputEvent, AliTriggerAnalysis::kV0C);
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
      AliVTrack   *track = (AliVTrack*)fInputEvent->GetTrack(iTrack);
      AliESDtrack *esdt  = dynamic_cast<AliESDtrack*>(track);
      AliAODTrack *aodt  = dynamic_cast<AliAODTrack*>(track);
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
   AliDebugClass(1, Form("Stats for event %d: %s", fEvNum, msg.Data()));
   if (isSelected) {
      fHEventStat->Fill(3.1);
      return output;
   } else {
      return 0;
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
      def = (AliRsnMiniOutput*)fHistograms[id];
      if (!def) continue;
      if (!def->IsMother()) continue;
      for (ip = 0; ip < npart; ip++) {
         AliMCParticle *part = (AliMCParticle*)fMCEvent->GetTrack(ip);
         if (TMath::Abs(part->Particle()->GetPdgCode()) != def->GetMotherPDG()) continue;
         // check that daughters match expected species
         if (part->Particle()->GetNDaughters() < 2) continue;
         label1 = part->Particle()->GetDaughter(0);
         label2 = part->Particle()->GetDaughter(1);
         daughter1 = (AliMCParticle*)fMCEvent->GetTrack(label1);
         daughter2 = (AliMCParticle*)fMCEvent->GetTrack(label2);
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
         def->Fill(&miniPair, miniEvent, &fValues);
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
      def = (AliRsnMiniOutput*)fHistograms[id];
      if (!def) continue;
      if (!def->IsMother()) continue;
      for (ip = 0; ip < npart; ip++) {
         AliAODMCParticle *part = (AliAODMCParticle*)list->At(ip);
         if (TMath::Abs(part->GetPdgCode()) != def->GetMotherPDG()) continue;
         // check that daughters match expected species
         if (part->GetNDaughters() < 2) continue;
         label1 = part->GetDaughter(0);
         label2 = part->GetDaughter(1);
         daughter1 = (AliAODMCParticle*)list->At(label1);
         daughter2 = (AliAODMCParticle*)list->At(label2);
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
         def->Fill(&miniPair, miniEvent, &fValues);
      }
   }
}

//________________________________________________________________________
void AliRsnMiniAnalysisTask::ProcessEvents(AliRsnMiniEvent *evMain, AliRsnMiniEvent *evMix)
{
//
// This function processes all single tracks.
// Loops on the daughters defined in current event
// and on all of them which pass the single-track cuts defined, computes what needed
//

   // check if it is mixed
   Bool_t isMix;
   if (evMix) {
      isMix = kTRUE;
   } else {
      evMix = evMain;
      isMix = kFALSE;
   }
   AliDebugClass(1, Form("Event %d: processing %s", fEvNum, (isMix ? "mixing" : "single event")));
   
   // 2 nested loops on tracks
   // inner starts from next track or from first, depending if is mixing or not
   Int_t i1, i2, idef, start;
   Int_t n1 = evMain->Particles().GetEntries();
   Int_t n2 = evMix ->Particles().GetEntries();
   Int_t ndef = fHistograms.GetEntries();
   AliRsnMiniParticle *daughter1 = 0x0, *daughter2 = 0x0;
   AliRsnMiniOutput *def = 0x0;
   for (i1 = 0; i1 < n1; i1++) {
      daughter1 = (AliRsnMiniParticle*)evMain->Particles()[i1];
      if (!daughter1) continue;
      if (isMix) start = 0; else start = i1 + 1;
      for (i2 = start; i2 < n2; i2++) {
         daughter2 = (AliRsnMiniParticle*)evMix->Particles()[i2];
         if (!daughter2) continue;
         // loop on definitions
         for (idef = 0; idef < ndef; idef++) {
            def = (AliRsnMiniOutput*)fHistograms[idef];
            if (!def) continue;
            if (def->IsEventOnly()) continue;
            if (def->IsMother()) continue;
            if (isMix && !def->IsTrackPairMix()) continue;
            if (!isMix && def->IsTrackPairMix()) continue;
            // check if definitions here match the current pair
            def->Fill(daughter1, daughter2, evMain, &fValues);
         } // end loop on definitions
      } // end loop on daughter #2
   } // end loop on daughter #1
}
