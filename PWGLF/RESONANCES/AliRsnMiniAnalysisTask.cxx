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
// Author: A. Pulvirenti
// Developers: F. Bellini (fbellini@cern.ch)

#include <Riostream.h>

#include <TH1.h>
#include <TList.h>
#include <TTree.h>
#include <TStopwatch.h>
#include "TRandom.h"

#include "AliLog.h"
#include "AliEventplane.h"
#include "AliMultiplicity.h"
#include "AliTriggerAnalysis.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisUtils.h"

#include "AliESDtrackCuts.h"
#include "AliESDUtils.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

#include "AliMultSelection.h"

#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsQnVector.h"

#include "AliRsnCutSet.h"
#include "AliRsnMiniPair.h"
#include "AliRsnMiniEvent.h"
#include "AliRsnMiniParticle.h"

#include "AliRsnMiniAnalysisTask.h"
#include "AliRsnMiniResonanceFinder.h"
//#include "AliSpherocityUtils.h"

ClassImp(AliRsnMiniAnalysisTask)

//__________________________________________________________________________________________________
/// Default constructor
AliRsnMiniAnalysisTask::AliRsnMiniAnalysisTask() :
   AliAnalysisTaskSE(),
   fEventCut(0x0),
   fUseMC(kFALSE),
   fEvNum(0),
   fTriggerMask(0),
   fSkipTriggerMask(0),
   fUseCentrality(kFALSE),
   fCentralityType("QUALITY"),
   fRefMultiType("GLOBAL"),
   fUseAOD049CentralityPatch(kFALSE),
   fUseCentralityPatchPbPb2011(0),
   fFlowQnVectorMgr(0),
   fFlowQnVectorSubDet("VZEROA"),
   fFlowQnVectorExpStep("latest"),
   fContinuousMix(kTRUE),
   fNMix(0),
   fMaxDiffMult(10),
   fMaxDiffVz(1.0),
   fMaxDiffAngle(1E20),
   fOutput(0x0),
   fHistograms("AliRsnMiniOutput", 0),
   fValues("AliRsnMiniValue", 0),
   fHEventStat(0x0),
   fHAEventsVsMulti(0x0),
   fHAEventsVsTracklets(0x0),
   fHAEventVzCent(0x0),
   fHAEventSpherocityCent(0x0),
   fHAEventMultiCent(0x0),
   fHAEventRefMultiCent(0x0),
   fHAEventPlane(0x0),
   fUseBuiltinEventCuts(kFALSE),
   fEventCuts(0x0),
   fTrackCuts(0),
   fRsnEvent(),
   fEvBuffer(0x0),
   fTriggerAna(0x0),
   fESDtrackCuts(0x0),
   fMiniEvent(0x0),
   fBigOutput(kFALSE),
   fMixPrintRefresh(-1),
   fCheckDecay(kTRUE),
   fMaxNDaughters(-1),
   fCheckP(kFALSE),
   fCheckFeedDown(kFALSE),
   fOriginDselection(kFALSE),
   fKeepDfromB(kFALSE),
   fKeepDfromBOnly(kFALSE),
   fRejectIfNoQuark(kFALSE),
   fMotherAcceptanceCutMinPt(0.0),
   fMotherAcceptanceCutMaxEta(0.9),
   fKeepMotherInAcceptance(kFALSE),
   fRsnTreeInFile(kFALSE),
   fComputeSpherocity(kFALSE),
   fTrackFilter(0x0),
   fSpherocity(-10),
   fResonanceFinders(0)
{
//
// Dummy constructor ALWAYS needed for I/O.
//
}

//__________________________________________________________________________________________________
/// Main constructor
///
/// \param name Name of the task
/// \param useMC Flag for Monte Carlo
/// \param saveRsnTreeInFile Flag for saving tree of rsn daughters in file
/// 
AliRsnMiniAnalysisTask::AliRsnMiniAnalysisTask(const char *name, Bool_t useMC,Bool_t saveRsnTreeInFile) :
   AliAnalysisTaskSE(name),
   fEventCut(0x0),
   fUseMC(useMC),
   fEvNum(0),
   fTriggerMask(AliVEvent::kMB),
   fSkipTriggerMask(0),
   fUseCentrality(kFALSE),
   fCentralityType("QUALITY"),
   fRefMultiType("GLOBAL"),
   fUseAOD049CentralityPatch(kFALSE),
   fUseCentralityPatchPbPb2011(0),
   fFlowQnVectorMgr(0),
   fFlowQnVectorSubDet("VZEROA"),
   fFlowQnVectorExpStep("latest"),
   fContinuousMix(kTRUE),
   fNMix(0),
   fMaxDiffMult(10),
   fMaxDiffVz(1.0),
   fMaxDiffAngle(1E20),
   fOutput(0x0),
   fHistograms("AliRsnMiniOutput", 0),
   fValues("AliRsnMiniValue", 0),
   fHEventStat(0x0),
   fHAEventsVsMulti(0x0),
   fHAEventsVsTracklets(0x0),
   fHAEventVzCent(0x0),
   fHAEventSpherocityCent(0x0),
   fHAEventMultiCent(0x0),
   fHAEventRefMultiCent(0x0),
   fHAEventPlane(0x0),
   fUseBuiltinEventCuts(kFALSE),
   fEventCuts(0x0),
   fTrackCuts(0),
   fRsnEvent(),
   fEvBuffer(0x0),
   fTriggerAna(0x0),
   fESDtrackCuts(0x0),
   fMiniEvent(0x0),
   fBigOutput(kFALSE),
   fMixPrintRefresh(-1),
   fCheckDecay(kTRUE),
   fMaxNDaughters(-1),
   fCheckP(kFALSE),
   fCheckFeedDown(kFALSE),
   fOriginDselection(kFALSE),
   fKeepDfromB(kFALSE),
   fKeepDfromBOnly(kFALSE),
   fRejectIfNoQuark(kFALSE),
   fMotherAcceptanceCutMinPt(0.0),
   fMotherAcceptanceCutMaxEta(0.9),
   fKeepMotherInAcceptance(kFALSE),
   fRsnTreeInFile(saveRsnTreeInFile),
   fComputeSpherocity(kFALSE),
   fTrackFilter(0x0),
   fSpherocity(-10),
   fResonanceFinders(0)
{
//
// Default constructor.
// Define input and output slots here (never in the dummy constructor)
// Input slot #0 works with a TChain - it is connected to the default input container
// Output slot #1 writes into a TH1 container
//

   DefineOutput(1, TList::Class());
   if (fRsnTreeInFile) DefineOutput(2, TTree::Class());
}

//__________________________________________________________________________________________________
/// Copy constructor
AliRsnMiniAnalysisTask::AliRsnMiniAnalysisTask(const AliRsnMiniAnalysisTask &copy) :
   AliAnalysisTaskSE(copy),
   // fEventCut(copy.fEventCut), // <- need to update!  (2020.05.08 blim)
   fUseMC(copy.fUseMC),
   fEvNum(0),
   fTriggerMask(copy.fTriggerMask),
   fSkipTriggerMask(copy.fSkipTriggerMask),
   fUseCentrality(copy.fUseCentrality),
   fCentralityType(copy.fCentralityType),
   fRefMultiType(copy.fRefMultiType),
   fUseAOD049CentralityPatch(copy.fUseAOD049CentralityPatch),
   fUseCentralityPatchPbPb2011(copy.fUseCentralityPatchPbPb2011),
   fFlowQnVectorMgr(copy.fFlowQnVectorMgr),
   fFlowQnVectorSubDet(copy.fFlowQnVectorSubDet),
   fFlowQnVectorExpStep(copy.fFlowQnVectorExpStep),
   fContinuousMix(copy.fContinuousMix),
   fNMix(copy.fNMix),
   fMaxDiffMult(copy.fMaxDiffMult),
   fMaxDiffVz(copy.fMaxDiffVz),
   fMaxDiffAngle(copy.fMaxDiffAngle),
   fOutput(0x0),
   fHistograms(copy.fHistograms),
   fValues(copy.fValues),
   fHEventStat(0x0),
   fHAEventsVsMulti(0x0),
   fHAEventsVsTracklets(0x0),
   fHAEventVzCent(0x0),
   fHAEventSpherocityCent(0x0),
   fHAEventMultiCent(0x0),
   fHAEventRefMultiCent(0x0),
   fHAEventPlane(0x0),
   fEventCuts(copy.fEventCuts),
   fTrackCuts(copy.fTrackCuts),
   fRsnEvent(),
   fEvBuffer(0x0),
   fTriggerAna(copy.fTriggerAna),
   fESDtrackCuts(copy.fESDtrackCuts),
   fMiniEvent(0x0),
   fBigOutput(copy.fBigOutput),
   fMixPrintRefresh(copy.fMixPrintRefresh),
   fCheckDecay(copy.fCheckDecay),
   fMaxNDaughters(copy.fMaxNDaughters),
   fCheckP(copy.fCheckP),
   fCheckFeedDown(copy.fCheckFeedDown),
   fOriginDselection(copy.fOriginDselection),
   fKeepDfromB(copy.fOriginDselection),
   fKeepDfromBOnly(copy.fKeepDfromBOnly),
   fRejectIfNoQuark(copy.fRejectIfNoQuark),
   fMotherAcceptanceCutMinPt(copy.fMotherAcceptanceCutMinPt),
   fMotherAcceptanceCutMaxEta(copy.fMotherAcceptanceCutMaxEta),
   fKeepMotherInAcceptance(copy.fKeepMotherInAcceptance),
   fRsnTreeInFile(copy.fRsnTreeInFile),
   fComputeSpherocity(copy.fComputeSpherocity),
   fTrackFilter(copy.fTrackFilter),
   fSpherocity(copy.fSpherocity),
   fResonanceFinders(copy.fResonanceFinders)
{
//
// Copy constructor.
// Implemented as requested by C++ standards.
// Can be used in PROOF and by plugins.
//
}

//__________________________________________________________________________________________________
/// Assignment Operator 
AliRsnMiniAnalysisTask &AliRsnMiniAnalysisTask::operator=(const AliRsnMiniAnalysisTask &copy)
{
//
// Assignment operator.
// Implemented as requested by C++ standards.
// Can be used in PROOF and by plugins.
//

   AliAnalysisTaskSE::operator=(copy);
   if (this == &copy)
      return *this;
   // fEventCut = copy.fEventCut, // <- need to update! (2020.05.08 blim)
   fUseMC = copy.fUseMC;
   fEvNum = copy.fEvNum;
   fTriggerMask = copy.fTriggerMask;
   fSkipTriggerMask = copy.fSkipTriggerMask;
   fUseCentrality = copy.fUseCentrality;
   fCentralityType = copy.fCentralityType;
   fRefMultiType = copy.fRefMultiType;
   fUseAOD049CentralityPatch = copy.fUseAOD049CentralityPatch;
   fUseCentralityPatchPbPb2011 = copy.fUseCentralityPatchPbPb2011;
   fFlowQnVectorMgr = copy.fFlowQnVectorMgr;
   fFlowQnVectorSubDet = copy.fFlowQnVectorSubDet;
   fFlowQnVectorExpStep = copy.fFlowQnVectorExpStep;
   fContinuousMix = copy.fContinuousMix;
   fNMix = copy.fNMix;
   fMaxDiffMult = copy.fMaxDiffMult;
   fMaxDiffVz = copy.fMaxDiffVz;
   fMaxDiffAngle = copy.fMaxDiffAngle;
   fHistograms = copy.fHistograms;
   fValues = copy.fValues;
   fHEventStat = copy.fHEventStat;
   fHAEventsVsMulti = copy.fHAEventsVsMulti;
   fHAEventsVsTracklets = copy.fHAEventsVsTracklets;
   fHAEventVzCent = copy.fHAEventVzCent;
   fHAEventSpherocityCent = copy.fHAEventSpherocityCent;
   fHAEventMultiCent = copy.fHAEventMultiCent;
   fHAEventRefMultiCent = copy.fHAEventRefMultiCent;
   fHAEventPlane = copy.fHAEventPlane;
   fEventCuts = copy.fEventCuts;
   fTrackCuts = copy.fTrackCuts;
   fTriggerAna = copy.fTriggerAna;
   fESDtrackCuts = copy.fESDtrackCuts;
   fBigOutput = copy.fBigOutput;
   fMixPrintRefresh = copy.fMixPrintRefresh;
   fCheckDecay = copy.fCheckDecay;
   fMaxNDaughters = copy.fMaxNDaughters;
   fCheckP = copy.fCheckP;
   fCheckFeedDown = copy.fCheckFeedDown;
   fOriginDselection = copy.fOriginDselection;
   fKeepDfromB = copy.fOriginDselection;
   fKeepDfromBOnly = copy.fKeepDfromBOnly;
   fRejectIfNoQuark = copy.fRejectIfNoQuark;
   fMotherAcceptanceCutMinPt = copy.fMotherAcceptanceCutMinPt;
   fMotherAcceptanceCutMaxEta = copy.fMotherAcceptanceCutMaxEta;
   fKeepMotherInAcceptance = copy.fKeepMotherInAcceptance;
   fRsnTreeInFile = copy.fRsnTreeInFile;
   fComputeSpherocity = copy.fComputeSpherocity;
   fTrackFilter = copy.fTrackFilter;
   fSpherocity = copy.fSpherocity;
   fResonanceFinders = copy.fResonanceFinders;

   return (*this);
}

//__________________________________________________________________________________________________
/// Destructor
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
/// Add a new cut set for a new criterion for track selection.
/// A user can add as many as (s)he wants, and each one corresponds
/// to one of the available bits in the AliRsnMiniParticle mask.
/// The only check is the following: if a cut set with the same name
/// as the argument is there, this is not added.
///
/// \param cuts Pointer to an AliRsnCutSet object
/// \return A value that is the array position of this set.
///
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
   Int_t v = 0;

   if (obj) {
      AliInfo(Form("A cut set named '%s' already exists", cuts->GetName()));
      v = fTrackCuts.IndexOf(obj);
   } else {
      fTrackCuts.AddLast(cuts);
      v = fTrackCuts.IndexOf(cuts);
   }

   for (Int_t i=0; i<fResonanceFinders.GetEntries(); i++){
      AliRsnMiniResonanceFinder* f = (AliRsnMiniResonanceFinder*) fResonanceFinders[i];
      if(f) f->IncrementResonanceCutID();
   }
    
   return v;
}

//__________________________________________________________________________________________________
/// Create output objects 
/// 
/// It creates a list of custom histograms and event counters. 
/// This is called once per worker node.
/// If required by setting the corresponding flag, 
/// it initialises the flow vector (Qn) from the external task.
///
void AliRsnMiniAnalysisTask::UserCreateOutputObjects()
{
//
// Initialization of outputs.
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

   //initialize quality trackcuts for spherocity
   if(!fTrackFilter){	
     fTrackFilter = new AliAnalysisFilter("trackFilter2015");
     SetTrackCuts(fTrackFilter);
   }

   // create list and set it as owner of its content (MANDATORY)
   if (fBigOutput) OpenFile(1);
   fOutput = new TList();
   fOutput->SetOwner();

   // initialize event statistics counter
   fHEventStat = new TH1F("hEventStat", "Event statistics", 16, 0.0, 16.0);
   fHEventStat->GetXaxis()->SetBinLabel(1, "CINT1B");
   fHEventStat->GetXaxis()->SetBinLabel(2, "V0AND");
   fHEventStat->GetXaxis()->SetBinLabel(3, "Candle");
   fHEventStat->GetXaxis()->SetBinLabel(4, "Accepted");
   fHEventStat->GetXaxis()->SetBinLabel(5, "Not Accepted - Total");
   fHEventStat->GetXaxis()->SetBinLabel(6, "Not Accepted - No Track Vertex");
   fHEventStat->GetXaxis()->SetBinLabel(7, "Not Accepted - Not Enough Contributors");
   fHEventStat->GetXaxis()->SetBinLabel(8, "Not Accepted - No Vertex inside |z| < 10 cm");
   fHEventStat->GetXaxis()->SetBinLabel(9, "Not Accepted - Pile Up Events");
   fHEventStat->GetXaxis()->SetBinLabel(10, "Not Accepted - SPD Vertex z Resolution");
   fHEventStat->GetXaxis()->SetBinLabel(11, "Not Accepted - SPD Vertex Dispersion");
   fHEventStat->GetXaxis()->SetBinLabel(12, "Not Accepted - z Difference vTrack - vSPD");
   fHEventStat->GetXaxis()->SetBinLabel(13, "Not Accepted - Incomplete DAQ");
   fHEventStat->GetXaxis()->SetBinLabel(14, "Not Accepted - Past/Future");
   fHEventStat->GetXaxis()->SetBinLabel(15, "Not Accepted - Cluster-Tracklets Correlation");
   fHEventStat->GetXaxis()->SetBinLabel(16, "Not Accepted - All local cuts");
   
   fOutput->Add(fHEventStat);

   if (!fHAEventsVsMulti) {
      if (fUseCentrality)
         fHAEventsVsMulti = new TH1F("hAEventsVsMulti", "Accepted events vs Centrality", 101, 0, 101.0);
      else
         fHAEventsVsMulti = new TH1F("hAEventsVsMulti", "Accepted events vs Multiplicity",1000, 0, 1000.0);
   }
   fOutput->Add(fHAEventsVsMulti);
   
   if (!fHAEventsVsTracklets)
      fHAEventsVsTracklets = new TH1F("hAEventsVsTracklets", "Accepted events vs Tracklet Number",1000, 0, 1000.0);
   fOutput->Add(fHAEventsVsTracklets);

   if(fHAEventVzCent) fOutput->Add(fHAEventVzCent);
   if(fHAEventMultiCent) fOutput->Add(fHAEventMultiCent);
   if(fHAEventRefMultiCent) fOutput->Add(fHAEventRefMultiCent);
   if(fHAEventSpherocityCent) fOutput->Add(fHAEventSpherocityCent);
   if(fHAEventPlane) fOutput->Add(fHAEventPlane);

   AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections *>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
   if (flowQnVectorTask) {
     AliInfo("Using Flow Qn vector corrections framework task ...");
     fFlowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
   }

   TIter next(&fTrackCuts);
   AliRsnCutSet *cs;
   while ((cs = (AliRsnCutSet *) next())) {
      cs->Init(fOutput);
   }

   // create temporary tree for filtered events
   if (fMiniEvent) SafeDelete(fMiniEvent);
   if (fRsnTreeInFile) OpenFile(2);
   fEvBuffer = new TTree("EventBuffer", "Temporary buffer for mini events");
   fMiniEvent = new AliRsnMiniEvent();
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
   if (fRsnTreeInFile) PostData(2, fEvBuffer);
}

//__________________________________________________________________________________________________
/// Main computation loop, performed for each event.
///
/// It checks if the event is acceptable, and eventually
/// creates the corresponding mini-event and stores it in the buffer.
/// The real histogram filling is done at the end, in "FinishTaskOutput".
///
/// Note: "true" mother-related histograms are filled in UserExec,
/// since they require direct access to MC event
///
void AliRsnMiniAnalysisTask::UserExec(Option_t *)
{
//
// Computation loop.
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

   for (Int_t i=0; i<fResonanceFinders.GetEntries(); i++){
      AliRsnMiniResonanceFinder* f = (AliRsnMiniResonanceFinder*) fResonanceFinders[i];
      if(f) f->RunResonanceFinder(fMiniEvent);
   }

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
   if (fRsnTreeInFile) PostData(2, fEvBuffer);
}

//__________________________________________________________________________________________________
/// Finish task. 
/// This function is called at the end of the loop on available events,
/// and then the buffer will be full with all the corresponding mini-events,
/// each one containing all tracks selected by each of the available track cuts.
/// Here a loop is done on each of these events, and both single-event and mixing are computed
///
void AliRsnMiniAnalysisTask::FinishTaskOutput()
{
//
// FInish task: loop on each of the selected events, and compute both single-event and mixing 
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

   Int_t printNum = fMixPrintRefresh;
   if (printNum < 0) {
      if (nEvents>1e5) printNum=nEvents/100;
      else if (nEvents>1e4) printNum=nEvents/10;
      else printNum = 0;
   }

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
         timer.Stop(); timer.Print(); fflush(stdout); timer.Start(kFALSE);
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
   timer.Stop(); timer.Print(); timer.Start(); fflush(stdout);

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
   timer.Stop(); timer.Print(); fflush(stdout); timer.Start();

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
   timer.Stop(); timer.Print(); fflush(stdout);

   // post computed data
   PostData(1, fOutput);
   if (fRsnTreeInFile) PostData(2, fEvBuffer);
}

//__________________________________________________________________________________________________
/// Terminate function. 
/// Called only once at the end.
/// It only checks that the output list is there.
///
void AliRsnMiniAnalysisTask::Terminate(Option_t *)
{
//
// Called once at the end of the query
//

   fOutput = dynamic_cast<TList *>(GetOutputData(1));
   if (!fOutput) {
      AliError("Could not retrieve TList fOutput");
      return;
   }
}

//__________________________________________________________________________________________________
/// This method checks if current event is OK for analysis.
///
/// In case it is, the pointers of the local AliRsnEvent data member
/// will point to it, in order to allow cut checking, otherwise the
/// function exits with a failure message.
/// ---
/// ESD events must pass the physics selection, AOD are supposed to do.
/// ---
/// While checking the event, a histogram is filled to count the number
/// of CINT1B, V0AND and CANDLE events, which are needed for normalization
/// ---
/// \return A character such as 'E' or 'A' if the event is accepted and is ESD or AOD
/// \return 0 if the event is not accepted   
///
Char_t AliRsnMiniAnalysisTask::CheckCurrentEvent()
{
//
// This method checks if current event is OK for analysis.
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
      isSelected = (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);

      if(isSelected && fSkipTriggerMask){
	      if( (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fSkipTriggerMask) == fSkipTriggerMask) 
            isSelected = kFALSE;
      }

      if (!isSelected && !fUseBuiltinEventCuts) {
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
   if(fUseBuiltinEventCuts){
      isSelected = fEventCut.AcceptEvent(fRsnEvent.GetRef());
      if(isSelected){
         msg += " -- Built-in AliEventCuts = ACCEPTED";
      }
      else{
         msg += " -- Built-in AliEventCuts = REJECTED";
         fHEventStat->Fill(15.1);
      }
   }
   else{
      if (fEventCuts) {
         if (!fEventCuts->IsSelected(&fRsnEvent)) {
            msg += " -- Local cuts = REJECTED";
            isSelected = kFALSE;
      // fill counter
      fHEventStat->Fill(15.1);
         } else {
            msg += " -- Local cuts = ACCEPTED";
            isSelected = kTRUE;
         }
      } else {
         msg += " -- Local cuts = NONE";
         isSelected = kTRUE;
      }
   }

   // if the above exit point is not taken, the event is accepted
   AliDebugClass(2, Form("Stats: %s", msg.Data()));
   if (isSelected) {
     
     fHEventStat->Fill(3.1);
     Double_t multi = ComputeCentrality((output == 'E'));
     Double_t tracklets = ComputeTracklets();
     Double_t refmulti = ComputeReferenceMultiplicity((output == 'E'), fRefMultiType.Data());
     fSpherocity = (fComputeSpherocity) ? ComputeSpherocity() : -10;
     fHAEventsVsMulti->Fill(multi);
     fHAEventsVsTracklets->Fill(tracklets);
     if(fHAEventVzCent) fHAEventVzCent->Fill(multi,fInputEvent->GetPrimaryVertex()->GetZ());
     if(fHAEventSpherocityCent && fComputeSpherocity) fHAEventSpherocityCent->Fill(multi,fSpherocity);
     if(fHAEventMultiCent) fHAEventMultiCent->Fill(multi,ComputeMultiplicity(output == 'E', fHAEventMultiCent->GetYaxis()->GetTitle()));
     if(fHAEventRefMultiCent) fHAEventRefMultiCent->Fill(refmulti, ComputeReferenceMultiplicity(output == 'E', fHAEventRefMultiCent->GetYaxis()->GetTitle()));
     if(fHAEventPlane) fHAEventPlane->Fill(multi,ComputeAngle());
     return output;
     
   } else {

     fHEventStat->Fill(4.1);
     const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
     if (!vertex) {
       fHEventStat->Fill(5.1);
     } else {
       TString title=vertex->GetTitle();
       AliRsnCutEventUtils* evtUtils = new AliRsnCutEventUtils("temporary_cutEventUtils");
       evtUtils->IsSelected(&fRsnEvent);
       AliRsnCutPrimaryVertex* cutPrimaryVertex = new AliRsnCutPrimaryVertex("temporary_cutPrimaryVertex");
       cutPrimaryVertex->IsSelected(&fRsnEvent);
       
       if( (title.Contains("Z")) || (title.Contains("3D")) || vertex->GetNContributors()<1.) {
	 if( (title.Contains("Z")) || (title.Contains("3D")) ) fHEventStat->Fill(5.1);
	 if(vertex->GetNContributors()<1.) fHEventStat->Fill(6.1);
       }
       else if(TMath::Abs(vertex->GetZ())>10.) fHEventStat->Fill(7.1);
       else if(((AliAnalysisUtils *)evtUtils->GetAnalysisUtils())->IsPileUpEvent(fInputEvent)) fHEventStat->Fill(8.1);
       else if(!cutPrimaryVertex->GoodZResolutionSPD()) fHEventStat->Fill(9.1);
       else if(!cutPrimaryVertex->GoodDispersionSPD()) fHEventStat->Fill(10.1);
       else if(!cutPrimaryVertex->GoodZDifferenceSPDTrack()) fHEventStat->Fill(11.1);
       else if(evtUtils->IsIncompleteDAQ()) fHEventStat->Fill(12.1);
       else if(evtUtils->FailsPastFuture()) fHEventStat->Fill(13.1);
       else if(evtUtils->IsSPDClusterVsTrackletBG(fInputEvent)) fHEventStat->Fill(14.1);
       SafeDelete(evtUtils);
       SafeDelete(cutPrimaryVertex);  
     }
     
     return 0;
   }
}

//__________________________________________________________________________________________________
/// Refresh cursor mini-event data member to fill with current event.
/// 
/// It fills the RsnMiniEvent with the event properties
/// and it fills the track pool after applying track cuts.
///
/// \param evType A character for ESD ('E') or AOD ('A') input
/// \return The total number of selected tracks.
/// 
void AliRsnMiniAnalysisTask::FillMiniEvent(Char_t evType)
{
   fMiniEvent->Clear();
   // assign event-related values
   fMiniEvent->SetRef(fRsnEvent.GetRef());
   fMiniEvent->SetRefMC(fRsnEvent.GetRefMC());
   fMiniEvent->Vz()    = fInputEvent->GetPrimaryVertex()->GetZ();
   fMiniEvent->Spherocity()  = fSpherocity;
   fMiniEvent->Angle() = ComputeAngle();
   fMiniEvent->Mult()  = ComputeCentrality((evType == 'E'));
   fMiniEvent->RefMult()  = ComputeReferenceMultiplicity((evType == 'E'), fRefMultiType.Data());
   fMiniEvent->Tracklets() = ComputeTracklets();
   AliDebugClass(2, Form("Event %d: type = %c -- vz = %f -- mult = %f -- angle = %f", fEvNum, evType, fMiniEvent->Vz(), fMiniEvent->Mult(), fMiniEvent->Angle()));

   if (fFlowQnVectorMgr) {
      TList *qnlist = fFlowQnVectorMgr->GetQnVectorList();
      if (qnlist) {
         fMiniEvent->SetQnVector(GetQnVectorFromList(qnlist, fFlowQnVectorSubDet.Data(), fFlowQnVectorExpStep.Data()));
      }
  }
   // loop on daughters and assign track-related values
   Int_t ic, ncuts = fTrackCuts.GetEntries();
   Int_t ip, npart = fRsnEvent.GetAbsoluteSum();
   Int_t npos = 0, nneg = 0, nneu = 0;
   AliRsnDaughter cursor;
  //  AliRsnMiniParticle miniParticle;
   AliRsnMiniParticle *miniParticlePtr;
   for (ip = 0; ip < npart; ip++) {
      // point cursor to next particle
      fRsnEvent.SetDaughter(cursor, ip);
      miniParticlePtr = fMiniEvent->AddParticle();
      miniParticlePtr->CopyDaughter(&cursor);
      miniParticlePtr->Index() = ip;

      AliAODTrack* aodtrack = cursor.Ref2AODtrack();
      if(aodtrack) miniParticlePtr->Index() = aodtrack->GetID();
      
      // copy momentum and MC info if present
      // miniParticle.CopyDaughter(&cursor);
      // miniParticle.Index() = ip;
      // switch on the bits corresponding to passed cuts
      for (ic = 0; ic < ncuts; ic++) {
         AliRsnCutSet *cuts = (AliRsnCutSet *)fTrackCuts[ic];
        //  if (cuts->IsSelected(&cursor)) miniParticle.SetCutBit(ic);
         if (cuts->IsSelected(&cursor)) miniParticlePtr->SetCutBit(ic);
        }
        // continue;
       
        // if a track passes at least one track cut, it is added to the pool
      // if (miniParticle.CutBits()) {
        if (miniParticlePtr->CutBits()) {
          // fMiniEvent->AddParticle(miniParticle);
        // if (miniParticle.Charge() == '+') npos++;
        //  else if (miniParticle.Charge() == '-') nneg++;
        //  else nneu++;
        if (miniParticlePtr->Charge() == '+') npos++;
        else if (miniParticlePtr->Charge() == '-') nneg++;
        else nneu++;
     } else {
        TClonesArray &arr = fMiniEvent->Particles();
        // Printf("B %d",arr.GetEntries());
        arr.RemoveAt(arr.GetEntries()-1);
        // Printf("A %d",arr.GetEntries());
      }
   }

   // get number of accepted tracks
   AliDebugClass(1, Form("Event %6d: total = %5d, accepted = %4d (pos %4d, neg %4d, neu %4d)", fEvNum, npart, (Int_t)fMiniEvent->Particles().GetEntriesFast(), npos, nneg, nneu));
}

//__________________________________________________________________________________________________
/// Compute event plane angle.
///
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
      plane = ((AliVAODHeader*)aodEvent->GetHeader())->GetEventplaneP();
   }

   if (plane)
      return plane->GetEventplane("Q");
   else {
      AliWarning("No event plane defined");
      return 1E20;
   }
}

//__________________________________________________________________________________________________
/// Computes event centrality/multiplicity.
/// 
/// COmputation is carried out according to the criterion defined
/// by two elements: (1) choice between multiplicity and centrality and
/// (2) the string defining what estimator must be used for specific computation.
///
/// \param isESD Flag =1 for ESD
/// \return Centrality/ multiplicity value (or percentile)
///
Double_t AliRsnMiniAnalysisTask::ComputeCentrality(Bool_t isESD)
{
   if (fUseCentrality) {
     if ((!fUseMC) && (fUseCentralityPatchPbPb2011)) {
       return ApplyCentralityPatchPbPb2011();//
    }
     if ((!fUseMC) && (!isESD) && (fUseAOD049CentralityPatch)) {
       return ApplyCentralityPatchAOD049();
     } else {
       AliCentrality *centrality = fInputEvent->GetCentrality();
         if (!centrality) {
	   AliError("Cannot compute centrality!");
	   return -1.0;
         }
         return centrality->GetCentralityPercentile(fCentralityType.Data());
     }
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
      } else if (fCentralityType.Contains("AliMultSelection")) {
	 TObjArray* a=(TObjArray*) fCentralityType.Tokenize("_");
	 if (!a) {
            AliWarning("Problem with AliRsnMiniAnalysisTask::fCentralityType string");
            return 1E20;
	 }

	 TObjString* o=(TObjString*) a->At(1);
	 if (!o) {
            AliWarning("Problem with AliRsnMiniAnalysisTask::fCentralityType string");
            return 1E20;
	 }

	 TString s=o->GetString();

	 AliMultSelection *MultSelection = (AliMultSelection*) fInputEvent -> FindListObject("MultSelection");
	 if (!MultSelection) {
            AliWarning("Could not find MultSelection object");
            return 1E20;
	 }

	 if (s.EqualTo("TEST")) {
	   MultSelection->PrintInfo();
	   return 50.;
	 }

	 return MultSelection->GetMultiplicityPercentile(s.Data());
      } else {
         AliError(Form("String '%s' does not define a possible multiplicity/centrality computation", fCentralityType.Data()));
         return -1.0;
      }
   }
}

//__________________________________________________________________________________________________
/// Computes event reference multiplicity according to the strategy defined by the "type"
///
/// Implementation follows AliPPVsMultUtils::GetStandardReferenceMultiplicity(...)
/// type = "TRACKLETS" uses SPD tracklets only
/// type = "QUALITY" (default) uses global tracks unless in cases when the track vertex is not available. Then tracklets are used anyway.
///
/// \param isESD Flag =1 for ESD
/// \param type Estimator 
/// \return Multiplicity value
/// 
Double_t AliRsnMiniAnalysisTask::ComputeReferenceMultiplicity(Bool_t isESD,TString type)
{

  type.ToUpper();
  Double_t computedRefMulti = -10.0;
  if ( type.CompareTo("TRACKLETS") && type.CompareTo("GLOBAL") ) {
    AliError(Form("String '%s' does not define a supported reference multiplicity computation", type.Data()));
    return computedRefMulti;
  }
  
  if (isESD) {
    AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(fInputEvent);
    if (!esdevent) return kFALSE;
    const AliESDVertex *lPrimaryVtxTracks = esdevent->GetPrimaryVertexTracks();
    if ( (!lPrimaryVtxTracks->GetStatus()) || (!type.CompareTo("TRACKLETS")) ) {
      /*If no track vertex available or if explicitly requested, use kTracklets
	= number of tracklets in the SPD within the specified eta range
      */
      computedRefMulti = AliESDtrackCuts::GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTracklets, 0.8, 0.0);
    } else {
      /* If track vertex available, use combined estimator
	 = number of global tracks/event
	 + number of ITS standalone tracks complementary to TPC for a given event
	 + number of SPD tracklets complementary to global/ITSSA tracks for a given events,
	 all of them within the specified eta range
      */
      computedRefMulti = AliESDtrackCuts::GetReferenceMultiplicity((AliESDEvent *)fInputEvent, AliESDtrackCuts::kTrackletsITSTPC, 0.8, 0.0);
    }
  } else {
    AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(fInputEvent);
    if (!aodevent) return kFALSE;
    AliAODHeader * header = dynamic_cast<AliAODHeader*>(aodevent->GetHeader());
    Long_t lStoredRefMult = header->GetRefMultiplicityComb08();
    // (-1) -> kept, user might discard, but will be killed anyhow by ev. sel.
    // (-2) -> kept, user might discard, but will be killed anyhow by ev. sel.
    // (-3) -> only if old AOD, but then -> recompute
    // (-4) -> kept, user might discard, but will be killed anyhow by ev. sel.
    //Hack: if this is -3, this is an old AOD filtering and we need to fall back on tracklets.
    if( (lStoredRefMult == -3)||(!type.CompareTo("TRACKLETS")) ) {
      //(then -1 and -2 are NOT the case, -4 unchecked but impossible since no track vtx)
      //Get Multiplicity object
      AliAODTracklets *spdmult = aodevent->GetMultiplicity();
      Long_t lNTracklets = 0;
      for (Int_t i=0; i<spdmult->GetNumberOfTracklets(); ++i)
	{
	  if ( TMath::Abs(spdmult->GetEta(i)) < 0.8 ) lNTracklets++;
	}
      computedRefMulti = lNTracklets;
    } else {
      computedRefMulti = lStoredRefMult;
    }
  }
  
  return computedRefMulti;
}

//__________________________________________________________________________________________________
/// Computes event multiplicity according to the specified estimator
/// 
/// Attention: **Deprecated**, better use ComputeCentrality() instead
///
/// \param isESD Flag =1 for ESD
/// \param type Estimator 
/// \return Multiplicity value
/// 
Double_t AliRsnMiniAnalysisTask::ComputeMultiplicity(Bool_t isESD,TString type)
{
   type.ToUpper();

   if (!type.CompareTo("TRACKS"))
      return fInputEvent->GetNumberOfTracks();
   else if (!type.CompareTo("QUALITY"))
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
   else if (!type.CompareTo("TRACKLETS")) {
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
      AliError(Form("String '%s' does not define a possible multiplicity/centrality computation", type.Data()));
      return -1.0;
   }
}

//__________________________________________________________________________________________________
/// Computes number of tracklets
/// 
/// \return Number of tracklets
/// 
Double_t AliRsnMiniAnalysisTask::ComputeTracklets()
{
  Double_t nTr = 100;
  Double_t count = 0.0;

  if (fInputEvent->InheritsFrom(AliESDEvent::Class())){
    AliESDEvent *esdEvent = (AliESDEvent *)fInputEvent;
    const AliMultiplicity *spdmult = esdEvent->GetMultiplicity();
    nTr = 1.0*spdmult->GetNumberOfTracklets();
    for(Int_t iTr=0; iTr<nTr; iTr++){
      Double_t theta=spdmult->GetTheta(iTr);
      Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
      if(eta>-1.0 && eta<1.0) count++;
    }
  }
  else if (fInputEvent->InheritsFrom(AliAODEvent::Class())) {
    AliAODEvent *aodEvent = (AliAODEvent *)fInputEvent;
    AliAODTracklets *spdmult = aodEvent->GetTracklets();
    nTr = 1.0*spdmult->GetNumberOfTracklets();
    for(Int_t iTr=0; iTr<nTr; iTr++){
      Double_t theta=spdmult->GetTheta(iTr);
      Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
      if(eta>-1.0 && eta<1.0) count++;
    }
  }

  return count;
}

//__________________________________________________________________________________________________
/// Computes event spherocity.
/// 
/// For AODs, the filter bit 5 is used to select tracks for the computation. 
/// For tESDs, the track filter set externally is used to select tracks. 
///
/// \return Spherocity value if at least 10 good tracks are used for the computation. If not, returns -10.
/// 
Double_t AliRsnMiniAnalysisTask::ComputeSpherocity()
{
  AliVEvent * evTypeS = InputEvent();
  Int_t ntracksLoop = evTypeS->GetNumberOfTracks();
  Int_t GoodTracks = 0;
  Float_t spherocity = -10.0;
  Float_t pFull = 0;
  Float_t Spherocity = 2;
  Float_t pt[10000],phi[1000];
  
  //computing total pt
  Float_t sumapt = 0;
  for(Int_t i1 = 0; i1 < ntracksLoop; ++i1){
    AliVTrack   *track = (AliVTrack *)evTypeS->GetTrack(i1);
    AliAODTrack *aodt  = dynamic_cast<AliAODTrack *>(track);
    AliESDtrack *esdt  = dynamic_cast<AliESDtrack *>(track);
    if (aodt) if (!aodt->TestFilterBit(5)) continue;
    if (esdt) if (!fTrackFilter->IsSelected(esdt)) continue;
    if (track->Pt() < 0.15) continue;
    if(TMath::Abs(track->Eta()) > 0.8) continue;
    //pt[i1] = track->Pt();
    pt[i1] = 1.0;
    sumapt += pt[i1];
    GoodTracks++;
  }
  if (GoodTracks < 10) return -10.0;
  //Getting thrust
  for(Int_t i = 0; i < 360/0.1; ++i){
	Float_t numerador = 0;
	Float_t phiparam  = 0;
	Float_t nx = 0;
	Float_t ny = 0;
	phiparam=( (TMath::Pi()) * i * 0.1 ) / 180; // parametrization of the angle
	nx = TMath::Cos(phiparam);            // x component of an unitary vector n
	ny = TMath::Sin(phiparam);            // y component of an unitary vector n
	for(Int_t i1 = 0; i1 < ntracksLoop; ++i1){
	  AliVTrack   *track = (AliVTrack *)evTypeS->GetTrack(i1);
	  AliAODTrack *aodt  = dynamic_cast<AliAODTrack *>(track);
	  AliESDtrack *esdt  = dynamic_cast<AliESDtrack *>(track);
	  if (aodt) if (!aodt->TestFilterBit(5)) continue;
	  if (esdt) if (!fTrackFilter->IsSelected(esdt)) continue;
	  if (track->Pt() < 0.15) continue;
	  if(TMath::Abs(track->Eta()) > 0.8) continue;
	  //pt[i1] = track->Pt();
	  pt[i1] = 1.0;
	  phi[i1] = track->Phi();
	  Float_t pxA = pt[i1] * TMath::Cos( phi[i1] );
	  Float_t pyA = pt[i1] * TMath::Sin( phi[i1] );
	  numerador += TMath::Abs( ny * pxA - nx * pyA );//product between p  proyection in XY plane and the unitary vector
	}
	pFull=TMath::Power( (numerador / sumapt),2 );
	if(pFull < Spherocity)//maximization of pFull
	  {
	    Spherocity = pFull;
	  }
  }
  spherocity=((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;
  if (GoodTracks > 9) return spherocity;
  else return -10.0;
}

//__________________________________________________________________________________________________
/// Fills the histograms with mother (ESD version)
///
/// By accessing the Monte Carlo information, it looks for the 
/// desired mother (resonance) and desired decay channel.  
///
/// \param miniEvent Pointer to the miniEvent object to process
/// 
void AliRsnMiniAnalysisTask::FillTrueMotherESD(AliRsnMiniEvent *miniEvent)
{
   Bool_t okMatch;
   Int_t id, ndef = fHistograms.GetEntries();
   Int_t ip, label1, label2, npart = fMCEvent->GetNumberOfTracks();
   static AliRsnMiniPair miniPair;
   AliMCParticle *daughter1, *daughter2;
   TLorentzVector p1, p2;
   AliRsnMiniOutput *def = 0x0;

   for (Int_t i=0; i<fResonanceFinders.GetEntries(); i++){
      AliRsnMiniResonanceFinder* f = (AliRsnMiniResonanceFinder*) fResonanceFinders[i];
      if(f) f->FillMother(fMCEvent, miniEvent);
   }

   for (id = 0; id < ndef; id++) {
      def = (AliRsnMiniOutput *)fHistograms[id];
      if (!def) continue;
      if (!def->IsMother() && !def->IsMotherInAcc() && !def->IsSingle()) continue;
      for (ip = 0; ip < npart; ip++) {
         AliMCParticle *part = (AliMCParticle *)fMCEvent->GetTrack(ip);
	 
         //get mother pdg code
         if (!AliRsnDaughter::IsEquivalentPDGCode(part->Particle()->GetPdgCode() , def->GetMotherPDG())) continue;
         if(def->IsSingle()) {
           def->FillSingle(part, miniEvent, &fValues);
           continue;
         }

         // check that daughters match expected species
         if (part->Particle()->GetNDaughters() < 2) continue;
	      if (fMaxNDaughters > 0 && part->Particle()->GetNDaughters() > fMaxNDaughters) continue;
         //label1 = part->Particle()->GetDaughter(0); // Before Change in accessing MC infor in AliRoot v5-09-46
         //label2 = part->Particle()->GetDaughter(1); // Before Change in accessing MC infor in AliRoot v5-09-46
         label1 = part->GetDaughterLabel(0);
         label2 = part->GetDaughterLabel(1);
          
         daughter1 = (AliMCParticle *)fMCEvent->GetTrack(label1);
         daughter2 = (AliMCParticle *)fMCEvent->GetTrack(label2);
         okMatch = kFALSE;
         if (AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter1->Particle()->GetPdgCode()) , def->GetPDG(0))
	     && AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter2->Particle()->GetPdgCode()) , def->GetPDG(1))) {
            okMatch = kTRUE;
            p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(0));
            p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(1));
         } else if (AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter1->Particle()->GetPdgCode()) , def->GetPDG(1))
		    && AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter2->Particle()->GetPdgCode()) , def->GetPDG(0))) {
            okMatch = kTRUE;
            p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(1));
            p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(0));
         }
         if (fCheckDecay && !okMatch) continue;
	      if (fCheckP && (TMath::Abs(part->Px()-(daughter1->Px()+daughter2->Px()))/(TMath::Abs(part->Px())+1.e-13)) > 0.00001 &&
     				(TMath::Abs(part->Py()-(daughter1->Py()+daughter2->Py()))/(TMath::Abs(part->Py())+1.e-13)) > 0.00001 &&
     				(TMath::Abs(part->Pz()-(daughter1->Pz()+daughter2->Pz()))/(TMath::Abs(part->Pz())+1.e-13)) > 0.00001 ) continue;
         if (fCheckFeedDown){
            Int_t pdgGranma = 0;
            Int_t mother = 0;
            mother = part->GetMother();
            Int_t istep = 0;
            Int_t abspdgGranma =0;
            Bool_t isFromB=kFALSE;
            Bool_t isQuarkFound=kFALSE;
            while (mother >=0 ){
               istep++;
               AliDebug(2,Form("mother at step %d = %d", istep, mother));
               AliMCParticle* mcGranma = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(mother));
               if (mcGranma){
                  pdgGranma = mcGranma->PdgCode();
                  AliDebug(2,Form("Pdg mother at step %d = %d", istep, pdgGranma));
                  abspdgGranma = TMath::Abs(pdgGranma);
                  if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
                  isFromB=kTRUE;
                  }
                  if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
                  mother = mcGranma->GetMother();
               }else{
                  AliError("Failed casting the mother particle!");
                  break;
               }
            }
            if(fRejectIfNoQuark && !isQuarkFound) pdgGranma = -99999;
            if(isFromB){
            if (!fKeepDfromB) pdgGranma = -9999; //skip particle if come from a B meson.
            }
            else{
            if (fKeepDfromBOnly) pdgGranma = -999;
            
            if (pdgGranma == -99999){
               AliDebug(2,"This particle does not have a quark in his genealogy\n");
               continue;
            }
            if (pdgGranma == -9999){
               AliDebug(2,"This particle come from a B decay channel but according to the settings of the task, we keep only the prompt charm particles\n");
               continue;
            }
            
            if (pdgGranma == -999){
               AliDebug(2,"This particle come from a prompt charm particles but according to the settings of the task, we want only the ones coming from B\n");
               continue;
            }

               }
         }
         // assign momenta to computation object
         miniPair.Sum(0) = miniPair.Sum(1) = (p1 + p2);
         miniPair.FillRef(def->GetMotherMass());
         miniPair.P1(1) = p1;
         miniPair.P2(1) = p2;

         // do computations and fill output
         def->FillMother(&miniPair, miniEvent, &fValues);
         if (fKeepMotherInAcceptance){
	         if(daughter1->Pt()<fMotherAcceptanceCutMinPt || daughter2->Pt()<fMotherAcceptanceCutMinPt || TMath::Abs(daughter1->Eta())>fMotherAcceptanceCutMaxEta ||  TMath::Abs(daughter2->Eta())>fMotherAcceptanceCutMaxEta) continue;
	         def->FillMotherInAcceptance(&miniPair, miniEvent, &fValues);
	      }
      }
   }
   return;
}

//__________________________________________________________________________________________________
/// Fills the histograms with mother (AOD version)
///
/// By accessing the Monte Carlo information, it looks for the 
/// desired mother (resonance) and desired decay channel.  
///
/// \param miniEvent Pointer to the miniEvent object to process
/// 
void AliRsnMiniAnalysisTask::FillTrueMotherAOD(AliRsnMiniEvent *miniEvent)
{
   Bool_t okMatch;
   TClonesArray *list = fRsnEvent.GetAODList();
   Int_t id, ndef = fHistograms.GetEntries();
   Int_t ip, label1, label2, npart = list->GetEntries();
   static AliRsnMiniPair miniPair;
   AliAODMCParticle *daughter1, *daughter2;
   TLorentzVector p1, p2;
   AliRsnMiniOutput *def = 0x0;

   for (Int_t i=0; i<fResonanceFinders.GetEntries(); i++){
      AliRsnMiniResonanceFinder* f = (AliRsnMiniResonanceFinder*) fResonanceFinders[i];
      if(f) f->FillMother(list, miniEvent);
   }

   for (id = 0; id < ndef; id++) {
      def = (AliRsnMiniOutput *)fHistograms[id];
      if (!def) continue;
      if (!def->IsMother() && !def->IsMotherInAcc() && !def->IsSingle()) continue;
      for (ip = 0; ip < npart; ip++) {
         AliAODMCParticle *part = (AliAODMCParticle *)list->At(ip);
	 
         if (!AliRsnDaughter::IsEquivalentPDGCode(part->GetPdgCode() , def->GetMotherPDG())) continue;
         if(def->IsSingle()) {
           def->FillSingle(part, miniEvent, &fValues);
           continue;
         }

         // check that daughters match expected species
         if (part->GetNDaughters() < 2) continue;
	 if (fMaxNDaughters > 0 && part->GetNDaughters() > fMaxNDaughters) continue;
         label1 = part->GetDaughterLabel(0);
         label2 = part->GetDaughterLabel(1);
         daughter1 = (AliAODMCParticle *)list->At(label1);
         daughter2 = (AliAODMCParticle *)list->At(label2);
         okMatch = kFALSE;
         if (AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter1->GetPdgCode()) , def->GetPDG(0))
	     && AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter2->GetPdgCode()) , def->GetPDG(1))) {
            okMatch = kTRUE;
            p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(0));
            p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(1));
         } else if (AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter1->GetPdgCode()) , def->GetPDG(1))
		    && AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter2->GetPdgCode()) , def->GetPDG(0))) {
            okMatch = kTRUE;
            p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(1));
            p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(0));
         }
         if (fCheckDecay && !okMatch) continue;
	 if(fCheckP && (TMath::Abs(part->Px()-(daughter1->Px()+daughter2->Px()))/(TMath::Abs(part->Px())+1.e-13)) > 0.00001 &&
     				(TMath::Abs(part->Py()-(daughter1->Py()+daughter2->Py()))/(TMath::Abs(part->Py())+1.e-13)) > 0.00001 &&
     				(TMath::Abs(part->Pz()-(daughter1->Pz()+daughter2->Pz()))/(TMath::Abs(part->Pz())+1.e-13)) > 0.00001 ) continue;
	 if(fCheckFeedDown){
		Int_t pdgGranma = 0;
		Int_t mother = 0;
		mother = part->GetMother();
		Int_t istep = 0;
		Int_t abspdgGranma =0;
		Bool_t isFromB=kFALSE;
		Bool_t isQuarkFound=kFALSE;
		while (mother >=0 ){
			istep++;
			AliDebug(2,Form("mother at step %d = %d", istep, mother));
			AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(list->At(mother));
			if (mcGranma){
				pdgGranma = mcGranma->GetPdgCode();
				AliDebug(2,Form("Pdg mother at step %d = %d", istep, pdgGranma));
				abspdgGranma = TMath::Abs(pdgGranma);
				if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
				  isFromB=kTRUE;
				}
				if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
				mother = mcGranma->GetMother();
			}else{
				AliError("Failed casting the mother particle!");
				break;
			}
		}
		if(fRejectIfNoQuark && !isQuarkFound) pdgGranma = -99999;
		if(isFromB){
		  if (!fKeepDfromB) pdgGranma = -9999; //skip particle if come from a B meson.
		}
		else{
		  if (fKeepDfromBOnly) pdgGranma = -999;
		  }
		  
		if (pdgGranma == -99999){
			AliDebug(2,"This particle does not have a quark in his genealogy\n");
			continue;
		}
		if (pdgGranma == -9999){
			AliDebug(2,"This particle come from a B decay channel but according to the settings of the task, we keep only the prompt charm particles\n");
			continue;
		}
		
		if (pdgGranma == -999){
			AliDebug(2,"This particle come from a prompt charm particles but according to the settings of the task, we want only the ones coming from B\n");
			continue;
		}
	 }
	 // assign momenta to computation object
         miniPair.Sum(0) = miniPair.Sum(1) = (p1 + p2);
         miniPair.FillRef(def->GetMotherMass());
	 miniPair.P1(1) = p1;
	 miniPair.P2(1) = p2;

         // do computations
         def->FillMother(&miniPair, miniEvent, &fValues);
	 if(fKeepMotherInAcceptance){
	      if(daughter1->Pt()<fMotherAcceptanceCutMinPt || daughter2->Pt()<fMotherAcceptanceCutMinPt || TMath::Abs(daughter1->Eta())>fMotherAcceptanceCutMaxEta ||  TMath::Abs(daughter2->Eta())>fMotherAcceptanceCutMaxEta) continue;
	      def->FillMotherInAcceptance(&miniPair, miniEvent, &fValues);
	 }
      }
   }
   return;
}
//___________________________________________________________
/// USed for D0 analysis. Selects the D0 selection criterium.
///
/// \param originDselection Set to 0, 1 or 2 for D0 from c quarks, b quarks or both ,respectively
///
void AliRsnMiniAnalysisTask::SetDselection(UShort_t originDselection)
{
	fOriginDselection = originDselection;
	
	if (fOriginDselection == 0) {
		fKeepDfromB = kFALSE;
		fKeepDfromBOnly = kFALSE;
	}
	
	if (fOriginDselection == 1) {
		fKeepDfromB = kTRUE;
		fKeepDfromBOnly = kTRUE;
	}
	
	if (fOriginDselection == 2) {
		fKeepDfromB = kTRUE;
		fKeepDfromBOnly = kFALSE;
	}
	
	return;
}
//__________________________________________________________________________________________________
/// Check if two events are compatible.
///
/// If the mixing is continuous, this is true if differences in vz, mult and angle are smaller than
/// the specified values.
/// If the mixing is binned, this is true if the events are in the same bin.
///
/// \param event1 First event pointer
/// \param event1 Second event pointer
/// \return Flag = 1 if events are compatible
/// 
Bool_t AliRsnMiniAnalysisTask::EventsMatch(AliRsnMiniEvent *event1, AliRsnMiniEvent *event2)
{
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

//---------------------------------------------------------------------
/// Patch to be used with 2011 Pb-Pb data for flat centrality distribution
///
/// \return Centrality value after flattening (reweighting)
///
Double_t AliRsnMiniAnalysisTask::ApplyCentralityPatchPbPb2011(){
  //This part rejects randomly events such that the centrality gets flat for LHC11h Pb-Pb data
  //for 0-5% and 10-20% centrality bin
  
  if (fCentralityType!="V0M") {
    AliWarning("Wrong value (not centrality from V0).");
    return -999.0;
  }
  
  AliCentrality *centrality = fInputEvent->GetCentrality();
  if (!centrality) {
    AliWarning("Cannot get centrality from AOD event.");
    return -999.0;
  }
  
  Double_t cent = (Float_t)(centrality->GetCentralityPercentile("V0M"));
  Double_t rnd_hc = -1., testf = 0.0, ff = 0, N1 = -1., N2 = -1.;

  if(fUseCentralityPatchPbPb2011==510){
    N1 = 1.9404e+06;
    N2 = 1.56435e+06; //N2 is the reference
    ff = 5.04167e+06 - 1.49885e+06*cent + 2.35998e+05*cent*cent -1.22873e+04*cent*cent*cent;
  } else {
    if(fUseCentralityPatchPbPb2011==1020){
      N2 = 2.0e+05; //N2 is the reference
      N1 = 3.7e+05;
      ff = -1.73979e+06 - 3.05316e+06*cent + 1.05517e+06*cent*cent - 133205*cent*cent*cent + 8187.45*cent*cent*cent*cent - 247.875*cent*cent*cent*cent*cent + 2.9676*cent*cent*cent*cent*cent*cent;
    } else {
      AliError(Form("Patch for the requested centrality (%i) is not available", fUseCentralityPatchPbPb2011));
      return -999.0;
    }
  }
  testf = ( N2 + (N1-ff) ) / N1;
  rnd_hc = gRandom->Rndm();

  //AliDebugClass(1, Form("Flat Centrality %d", fUseCentralityPatchPbPb2011));

  if (rnd_hc < 0 || rnd_hc > 1 )
    {
      AliWarning("Wrong Random number generated");
      return -999.0;
    }
  
  if (rnd_hc < testf)
    return cent;
  else
    return -999.0;
}
//---------------------------------------------------------------------
/// Patch to be used with 2010 Pb-Pb data AOD049 for a correct centrality distribution
///
/// \return Centrality value after re-parameterisation
///
Double_t AliRsnMiniAnalysisTask::ApplyCentralityPatchAOD049()
{
   //
   //Apply centrality patch for AOD049 outliers
   //
   if (fInputEvent->InheritsFrom(AliESDEvent::Class())) {
      AliWarning("Requested patch for AOD049 for ESD. ");
      return -999.0;
   }

   if (fCentralityType!="V0M") {
      AliWarning("Requested patch forAOD049 for wrong value (not centrality from V0).");
      return -999.0;
   }

   AliCentrality *centrality = fInputEvent->GetCentrality();
   if (!centrality) {
      AliWarning("Cannot get centrality from AOD event.");
      return -999.0;
   }

   Float_t cent = (Float_t)(centrality->GetCentralityPercentile("V0M"));

   if(cent>=0.0) {
      Float_t v0 = 0.0;
      AliAODEvent *aodEvent = (AliAODEvent *)fInputEvent;
      AliAODVZERO *aodV0 = (AliAODVZERO *) aodEvent->GetVZEROData();
      v0+=aodV0->GetMTotV0A();
      v0+=aodV0->GetMTotV0C();
      if ( (cent==0) && (v0<19500) ) {
         AliDebug(3, Form("Filtering issue in centrality -> cent = %5.2f",cent));
         return -999.0;
      }
      Float_t tkl = (Float_t)(aodEvent->GetTracklets()->GetNumberOfTracklets());
      Float_t val = 1.30552 +  0.147931 * v0;

      Float_t tklSigma[101] = {176.644, 156.401, 153.789, 153.015, 142.476, 137.951, 136.127, 129.852, 127.436, 124.86,
                               120.788, 115.611, 113.172, 110.496, 109.127, 104.421, 102.479, 99.9766, 97.5152, 94.0654,
                               92.4602, 89.3364, 87.1342, 83.3497, 82.6216, 81.1084, 78.0793, 76.1234, 72.9434, 72.1334,
                               68.0056, 68.2755, 66.0376, 62.9666, 62.4274, 59.65, 58.3776, 56.6361, 54.5184, 53.4224,
                               51.932, 50.8922, 48.2848, 47.912, 46.5717, 43.4114, 43.2083, 41.3065, 40.1863, 38.5255,
                               37.2851, 37.5396, 34.4949, 33.8366, 31.8043, 31.7412, 30.8392, 30.0274, 28.8793, 27.6398,
                               26.6488, 25.0183, 25.1489, 24.4185, 22.9107, 21.2002, 21.6977, 20.1242, 20.4963, 19.0235,
                               19.298, 17.4103, 16.868, 15.2939, 15.2939, 16.0295, 14.186, 14.186, 15.2173, 12.9504, 12.9504,
                               12.9504, 15.264, 12.3674, 12.3674, 12.3674, 12.3674, 12.3674, 18.3811, 13.7544, 13.7544,
                               13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544
                              };

      if ( TMath::Abs(tkl-val) > 6.*tklSigma[(Int_t)cent] )  {
         AliDebug(3, Form("Outlier event in centrality -> cent = %5.2f",cent));
         return -999.0;
      }
   } else {
      //force it to be -999. whatever the negative value was
      cent = -999.;
   }
   return cent;
}

//----------------------------------------------------------------------------------
/// Set event QA histograms
///
/// Used to monitor multiplicity distributions (different estimators),
/// event plane, spherocity, vertex position.
///
/// \param type String defining one specific monitoring histogram among those available
/// \param histo Histogram to be added
///
void AliRsnMiniAnalysisTask::SetEventQAHist(TString type,TH1 *histo)
{
   if(!histo) {
      AliWarning(Form("event QA histogram pointer not defined for slot %s",type.Data()));
      return;
   }

   type.ToLower();
   TString multitype(histo->GetYaxis()->GetTitle());
   multitype.ToUpper();
   
   if(!type.CompareTo("eventsvsmulti")) fHAEventsVsMulti = (TH1F*) histo;
   else if(!type.CompareTo("eventsvstracklets")) fHAEventsVsTracklets = (TH1F*) histo;
   else if(!type.CompareTo("vz")) fHAEventVzCent = (TH2F*) histo;
   else if(!type.CompareTo("spherocitycent")){
      fHAEventSpherocityCent = (TH2F*) histo;
      fComputeSpherocity = kTRUE;
   }
   else if(!type.CompareTo("multicent")) {
      if(multitype.CompareTo("QUALITY") && multitype.CompareTo("TRACKS") && multitype.CompareTo("TRACKLETS")) {
         AliWarning(Form("multiplicity vs. centrality histogram y-axis %s unknown, setting to TRACKS",multitype.Data()));
         histo->GetYaxis()->SetTitle("TRACKS");
      }
      fHAEventMultiCent = (TH2F*) histo;
   }
   else if(!type.CompareTo("refmulti")){
     if ( multitype.CompareTo("GLOBAL") && multitype.CompareTo("TRACKLETS") ) {
       AliWarning(Form("Reference multiplicity vs. centrality histogram y-axis %s unknown, setting to GLOBAL",multitype.Data()));
       histo->GetYaxis()->SetTitle("GLOBAL");
     }
     fHAEventRefMultiCent = (TH2F*) histo;
   }
   else if(!type.CompareTo("eventplane")) fHAEventPlane = (TH2F*) histo;
   else AliWarning(Form("event QA histogram slot %s undefined",type.Data()));

   return;
}

//----------------------------------------------------------------------------------
/// Create a new AliRsnMiniValue in the task
/// and return its ID, which is needed for setting up histograms.
/// If that value was already initialized, returns its ID and does not recreate it.
///
/// \param type Type of value to be created (one of those defined in AliRsnMiniValue)
/// \param useMC Flag = 1 if generated information has to be used (for MC only)
/// \return ID of the value, to be used for setting histograms
///
Int_t AliRsnMiniAnalysisTask::CreateValue(AliRsnMiniValue::EType type, Bool_t useMC)
{
   Int_t valID = ValueID(type, useMC);
   if (valID >= 0 && valID < fValues.GetEntries()) {
      AliInfo(Form("Value '%s' is already created in slot #%d", AliRsnMiniValue::ValueName(type, useMC), valID));
   } else {
      valID = fValues.GetEntries();
      AliInfo(Form("Creating value '%s' in slot #%d", AliRsnMiniValue::ValueName(type, useMC), valID));
      new (fValues[valID]) AliRsnMiniValue(type, useMC);
   }

   return valID;
}

//----------------------------------------------------------------------------------
/// Check if value is initialised
///
/// \param type Type of value to be created (one of those defined in AliRsnMiniValue)
/// \param useMC Flag = 1 if generated information has to be used (for MC only)
/// \return ID of the value, if value is initialised, otherwise returns -1
///
Int_t AliRsnMiniAnalysisTask::ValueID(AliRsnMiniValue::EType type, Bool_t useMC)
{
   const char *name = AliRsnMiniValue::ValueName(type, useMC);
   TObject *obj = fValues.FindObject(name);
   if (obj)
      return fValues.IndexOf(obj);
   else
      return -1;
}

//----------------------------------------------------------------------------------
/// Create a new histogram definition in the task.
/// and returns it to the user for its configuration
///
/// \param name Name of the output histogram
/// \param type Type of the output histogram (out of AliRsnMiniOutput-supported formats)
/// \param src Type of computation (out of AliRsnMiniOutput-supported ones)
/// \return Pointer to the output histogram
///
AliRsnMiniOutput *AliRsnMiniAnalysisTask::CreateOutput(const char *name, AliRsnMiniOutput::EOutputType type, AliRsnMiniOutput::EComputation src)
{
   Int_t n = fHistograms.GetEntries();
   AliRsnMiniOutput *newDef = new (fHistograms[n]) AliRsnMiniOutput(name, type, src);

   return newDef;
}

//----------------------------------------------------------------------------------
/// Create a new histogram definition in the task.
/// and returns it to the user for its configuration
///
/// \param name Name of the output histogram
/// \param type Type of the output histogram (SPARSE, TH1D, TH2D, TH3D)
/// \param src Type of computation (EVENT, PAIR, MIX, TRUE, MOTHER, ROTATE1, ROTATE2)
/// \return Pointer to the output histogram
///
AliRsnMiniOutput *AliRsnMiniAnalysisTask::CreateOutput(const char *name, const char *outType, const char *compType)
{
   Int_t n = fHistograms.GetEntries();
   AliRsnMiniOutput *newDef = new (fHistograms[n]) AliRsnMiniOutput(name, outType, compType);

   return newDef;
}


//----------------------------------------------------------------------------------
/// Gets the flow vector from an external list
///
/// \param list Pointer to the external list (from flow task)
/// \param subdetector String with subdetector to be used for Qn computation
/// \param expectedstep String for desired correction step in Qn framework 
/// \return Pointer to the Qn vector object
///
AliQnCorrectionsQnVector *AliRsnMiniAnalysisTask::GetQnVectorFromList(
  const TList *list, const char *subdetector, const char *expectedstep) const {

  AliQnCorrectionsQnVector *theQnVector = NULL;

  // TList *pQvecList = dynamic_cast<TList *>(list->FindObject(subdetector));
  TList *pQvecList = (TList*)list->FindObject(subdetector);
  if (pQvecList != NULL) {
    /* the detector is present */
    if (TString(expectedstep).EqualTo("latest"))
      theQnVector = (AliQnCorrectionsQnVector *)pQvecList->First();
    else
      theQnVector =
        (AliQnCorrectionsQnVector *)pQvecList->FindObject(expectedstep);
  }
  if (theQnVector != NULL) {
    /* check the Qn vector quality */
    if (!(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0))
      /* not good quality, discarded */
      theQnVector = NULL;
  }
  return theQnVector;
}

//----------------------------------------------------------------------------------
/// Add a new AliRsnMiniResonanceFinder object.
///
/// A user can add as many as he wants, and each one corresponds
/// to one of the available bits in the AliRsnMiniParticle mask.
/// The only check is the following: if a ResonanceFinder set with the same name
/// as the argument is there, this is not added.
///
/// \param f Pointer to the AliRsnMiniResonanceFinder object
/// \return Cut ID for the ResonanceFinder f.
///
Int_t AliRsnMiniAnalysisTask::AddResonanceFinder(AliRsnMiniResonanceFinder* f)
{
   TObject *obj = fResonanceFinders.FindObject(f->GetName());
   Int_t v = 0;

   if (obj) {
      AliInfo(Form("A ResonanceFinder named '%s' already exists", f->GetName()));
      v = fResonanceFinders.IndexOf(obj) + GetNumberOfTrackCuts();
   } else {
      fResonanceFinders.AddLast(f);
      v = fResonanceFinders.IndexOf(f) + GetNumberOfTrackCuts();
      f->SetResonanceCutID(v);
   }

   return v;
}
//----------------------------------------------------------------------------------
/// Defines track cuts for Spherocity computation.
///
/// Uses a track filter to set the quality cuts for the spherocity calculation.
/// Implementation only used for ESDs. AOD filtering is done via the filter bit. 
///
/// \param fTrackFilter Pointer to the AliAnalysisFilter 
///
void AliRsnMiniAnalysisTask::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

	AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts();
	//TPC Only
	esdTrackCuts->SetMinNClustersTPC(50);
	esdTrackCuts->SetMaxChi2PerClusterTPC(4);
	esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
	esdTrackCuts->SetMaxDCAToVertexZ(3.2);
	esdTrackCuts->SetMaxDCAToVertexXY(2.4);
	esdTrackCuts->SetDCAToVertex2D(kTRUE);
	
	esdTrackCuts->SetRequireTPCRefit(kTRUE);// TPC Refit
	esdTrackCuts->SetRequireITSRefit(kTRUE);// ITS Refit
	fTrackFilter->AddCuts(esdTrackCuts);
   return;
}
