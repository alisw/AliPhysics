/**
 * @file   AddTaskCentralTracks.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Fri Jan 28 10:22:26 2011
 * 
 * @brief Class and script to add a multiplicity task for the central
 *        @f$\eta@f$ region
 * 
 * @ingroup pwg2_forward_scripts_tasks
 * 
 */
#include <AliAnalysisTaskSE.h>
#include <AliESDtrackCuts.h>
class TList;
class TH2D;
class TH1D;

/**
 * Task to determine the 
 * @f[
 *   \left.\frac{d^2N_{ch}}{d\eta d\phi}\right|_{central}
 * @f] 
 * from tracks and tracklets. 
 *
 * First, global tracks are investigated. The requirements on global
 * tracks are
 * - @f$ n_{TPC clusters} \ge 70@f$ 
 * - @f$ \chi^2_{TPC cluster} \le 4@f$ 
 * - No daughter kinks 
 * - Re-fit of TPC tracks 
 * - Re-fit of ITS tracks 
 * - No requirement on SPD clusters 
 *
 * Secondly, ITS stand-alone tracks are investigated.  The
 * requirements on ITS standalone tracks are
 * - Re-fit of ITS tracks 
 * - No requirement on SPD clusters 
 * 
 * Tracks that does not meet these quality conditions are flagged as
 * rejected.
 *
 * Both kinds of tracks (global and ITS standalone) have requirements
 * on the distance of closest approach (DCA) to the interaction
 * vertex @f$(v_x,v_y,v_z)@f$. 
 *
 * For tracks with SPD clusters: 
 * - @f$ DCA_{xy} < 0.0182+0.0350/p_t^{1.01}@f$ 
 * - @f$ DCA_{z} < 0.5@f$ 
 *
 * For tracks without SPD clusters 
 * - @f$ DCA_{xy} < 1.5(0.0182+0.0350/p_t^{1.01})@f$ 
 * - @f$ DCA_{z} < 0.5@f$ 
 *
 * Tracks that does not meet these DCA requirements are flagged as
 * secondaries.
 *
 * Thirdly, the number of SPD tracklets are investigated.  If the
 * tracklet is associated with a track, and that track has already
 * been used, then that tracklet is ignored
 * 
 * An @f$(\eta,\phi)@f$ per-event histogram is then filled from these
 * tracks and tracklets, and that histogram is stored in the output
 * AOD event.
 *
 * Furthermore, a number of diagnostics @f$\eta@f$ histograms are filled: 
 * - @f$\eta@f$ of all accepted tracks and tracklets 
 * - @f$\eta@f$ of accepted global tracks
 * - @f$\eta@f$ of accepted ITS tracks
 * - @f$\eta@f$ of accepted SPD tracklets
 *
 * At the end of the job, these histograms are normalized to the
 * number of accepted events and bin width to provide an estimate of
 * the @f$dN_{ch}/d\eta@f$
 *
 * Only minimum bias events with a @f$v_z@f$ within the defined cut
 * are analysed.
 *
 * @ingroup pwg2_forward_aod
 */
class CentralMultTask : public AliAnalysisTaskSE
{
public:
  enum { 
    kValidTrigger = 0x1, 
    kValidVertex  = 0x2
  };
  /** 
   * Constructor 
   * 
   */
  CentralMultTask();
  /** 
   * Constructor
   * 
   * @param name    Name of task 
   * @param maxEta  Set @f$\eta@f$ range
   * @param maxVtx  Set @f$v_z@f$ range
   */
  CentralMultTask(const char* name, 
		  Double_t maxEta=2, 
		  Double_t maxVtx=10);
  /**
   * Destructor
   * 
   */
  virtual ~CentralMultTask();
  void SetUseTracklets(Bool_t use) { fUseTracklets = use; }
  /** 
   * Initialise on master - does nothing
   * 
   */
  virtual void   Init() {}
  /** 
   * Create output objects.  
   *
   * This is called once per slave process 
   */
  virtual void UserCreateOutputObjects();
  /** 
   * Process a single event 
   * 
   * @param option Not used
   */
  virtual void UserExec(Option_t* option);
  /** 
   * Called at end of event processing.
   *
   * This is called once in the master 
   * 
   * @param option Not used 
   */
  virtual void Terminate(Option_t* option);
protected:
  /** 
   * Check if this event is within cuts
   * 
   * @param esd Event
   * @param vz  Vertex @f$ z@f$ coordinate 
   * 
   * @return true if this event is to be considered
   */
  UShort_t CheckEvent(const AliESDEvent& esd, Double_t& vz);

  TH2D*           fHist;      // Per-event d2n/deta/dphi
  TH1D*           fAll;       // Summed d2n/deta/dphi - all track(let)s
  TH1D*           fGlobal;    // Summed d2n/deta/dphi - global tracks
  TH1D*           fITS;       // Summed d2n/deta/dphi - ITS tracks
  TH1D*           fTracklets; // Summed d2n/deta/dphi - SPD tracklets
  TH1D*           fNEventsTr; // Number of triggered events per vertex bin
  TH1D*           fNEventsVtx;// Number of triggered+vertex events per vertex
  TList*          fOutput;    // Output list
  AliESDtrackCuts fQGlo;      // Quality cut on ITS+TPC
  AliESDtrackCuts fQITS;      // Quality cut on ITS standalone (not ITS+TPC)
  AliESDtrackCuts fDCAwSPD;   // DCA for traks with SPD hits
  AliESDtrackCuts fDCAwoSPD;  // DCA for traks without SPD hits
  AliESDtrackCuts fIsPrimary; // Primary particle 
  Bool_t          fUseTracklets; // Whether to count tracklets or not 
  Bool_t          fDebug;     // Debug flag

  ClassDef(CentralMultTask,1); // Determine multiplicity in central area
};

//====================================================================
#include <TMath.h>
#include <TH2D.h>
#include <TH1D.h>
#include <THStack.h>
#include <TList.h>
#include <AliAnalysisManager.h>
#include <AliESDEvent.h>
#include <AliAODHandler.h>
#include <AliAODEvent.h>
#include <AliESDInputHandler.h>
#include <AliESDtrack.h>
#include <AliMultiplicity.h>

//____________________________________________________________________
inline CentralMultTask::CentralMultTask()
  : AliAnalysisTaskSE(), 
    fHist(0),
    fAll(0),
    fGlobal(0),
    fITS(0),
    fTracklets(0),
    fNEventsTr(0),
    fNEventsVtx(0),
    fOutput(0),
    fUseTracklets(false),
    fDebug(false)
{}

//____________________________________________________________________
inline CentralMultTask::CentralMultTask(const char* /* name */,
					Double_t maxEta,
					Double_t maxVtx)
  : AliAnalysisTaskSE("Central"), 
    fHist(0),
    fAll(0),
    fGlobal(0),
    fITS(0),
    fTracklets(0),
    fNEventsTr(0),
    fNEventsVtx(0),
    fOutput(0),
    fUseTracklets(true),
    fDebug(false)
{
  Int_t nBins = 40;
  fHist = new TH2D("Central", "d^{2}N_{ch}/d#etad#phi in central region",
		   nBins, -maxEta, maxEta, 20, 0, 2*TMath::Pi());
  fHist->SetDirectory(0);
  fHist->SetXTitle("#eta");
  fHist->SetYTitle("#phi [radians]");
  fHist->SetZTitle("d^{2}N_{ch}/d#etad#phi");
  fHist->SetStats(0);
  fHist->Sumw2();

  fAll = new TH1D("all", "Central region", nBins, -maxEta, maxEta);
  fAll->SetDirectory(0);
  fAll->SetXTitle("#eta");
  fAll->SetYTitle("dN_{ch}/d#eta");
  fAll->Sumw2();
  fAll->SetFillColor(kGray);
  fAll->SetFillStyle(3001);
  fAll->SetMarkerStyle(28);
  fAll->SetMarkerColor(kGray);
  fAll->SetStats(0);
  
  fGlobal = static_cast<TH1D*>(fAll->Clone("global"));
  fGlobal->SetTitle("Global tracks");
  fGlobal->SetDirectory(0);
  fGlobal->Sumw2();
  fGlobal->SetFillColor(kRed+1);
  fGlobal->SetFillStyle(3001);
  fGlobal->SetMarkerStyle(28);
  fGlobal->SetMarkerColor(kRed+1);
  fGlobal->SetStats(0);

  fITS = static_cast<TH1D*>(fAll->Clone("its"));
  fITS->SetTitle("ITS tracks");
  fITS->SetDirectory(0);
  fITS->Sumw2();
  fITS->SetFillColor(kGreen+1);
  fITS->SetFillStyle(3001);
  fITS->SetMarkerStyle(28);
  fITS->SetMarkerColor(kGreen+1);
  fITS->SetStats(0);

  fTracklets = static_cast<TH1D*>(fAll->Clone("tracklets"));
  fTracklets->SetTitle("SPD tracklets");
  fTracklets->SetDirectory(0);
  fTracklets->Sumw2();
  fTracklets->SetFillColor(kBlue+1);
  fTracklets->SetFillStyle(3001);
  fTracklets->SetMarkerStyle(28);
  fTracklets->SetMarkerColor(kBlue+1);
  fTracklets->SetStats(0);

  fNEventsTr = new TH1D("nEventsTr", 
			"Events per vertex bin", 10, -maxVtx, maxVtx);
  fNEventsTr->SetXTitle("v_{z}");
  fNEventsTr->SetYTitle("Events");
  fNEventsTr->SetDirectory(0);

  fNEventsVtx = new TH1D("nEventsVtx", 
			 "Events per vertex bin", 10, -maxVtx, maxVtx);
  fNEventsVtx->SetXTitle("v_{z}");
  fNEventsVtx->SetYTitle("Events");
  fNEventsVtx->SetDirectory(0);

  // --- Global (ITS+TPC) track quality cuts 
  // TPC  
  fQGlo.SetMinNClustersTPC(70);
  fQGlo.SetMaxChi2PerClusterTPC(4);
  fQGlo.SetAcceptKinkDaughters(kFALSE);
  fQGlo.SetRequireTPCRefit(kTRUE);
  // ITS
  fQGlo.SetRequireITSRefit(kTRUE);
  fQGlo.SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
  fQGlo.SetEtaRange(-maxEta, maxEta);

  // --- ITS standalone track quality cuts 
  fQITS.SetRequireITSRefit(kTRUE);
  fQITS.SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
  fQITS.SetEtaRange(-maxEta, maxEta); 

  // -- Distance-of-Closest-Approach cuts for tracks w/ITS hits 
  fDCAwSPD.SetClusterRequirementITS(AliESDtrackCuts::kSPD,
				    AliESDtrackCuts::kAny);
  fDCAwSPD.SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  fDCAwSPD.SetMaxDCAToVertexZ(0.5);
  fDCAwSPD.SetEtaRange(-maxEta, maxEta);

  // -- Distance-of-Closest-Approach cuts for tracks w/o ITS hits 
  fDCAwoSPD.SetClusterRequirementITS(AliESDtrackCuts::kSPD,
				     AliESDtrackCuts::kNone);
  fDCAwoSPD.SetMaxDCAToVertexXYPtDep("1.5*(0.0182+0.0350/pt^1.01)");
  fDCAwoSPD.SetMaxDCAToVertexZ(0.5);
  fDCAwoSPD.SetEtaRange(-maxEta, maxEta);  

  // -- Primary track cut 
  // https://twiki.cern.ch/twiki/bin/view/ALICE/SelectionOfPrimaryTracksForPpDataAnalysis
  // Quality
  fIsPrimary.SetMinNClustersTPC(70);
  fIsPrimary.SetMaxChi2PerClusterTPC(4);
  fIsPrimary.SetAcceptKinkDaughters(kFALSE);
  fIsPrimary.SetRequireTPCRefit(kTRUE);
  fIsPrimary.SetRequireITSRefit(kTRUE);
  fIsPrimary.SetClusterRequirementITS(AliESDtrackCuts::kSPD,
				      AliESDtrackCuts::kAny);
  //  Dca:
  fIsPrimary.SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  fIsPrimary.SetMaxDCAToVertexZ(2);
  fIsPrimary.SetDCAToVertex2D(kFALSE);
  fIsPrimary.SetRequireSigmaToVertex(kFALSE);  

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class()); 
}

//____________________________________________________________________
inline CentralMultTask::~CentralMultTask()
{
  if (fHist)      delete fHist;
  if (fAll)       delete fAll;
  if (fGlobal)    delete fGlobal;
  if (fITS)       delete fITS;
  if (fTracklets) delete fTracklets;
}

//________________________________________________________________________
inline void 
CentralMultTask::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)

  fOutput = new TList;
  fOutput->SetName(GetName());
  fOutput->SetOwner();

  fOutput->Add(fAll);
  fOutput->Add(fGlobal);
  fOutput->Add(fITS);
  fOutput->Add(fTracklets);
  fOutput->Add(fNEventsTr);
  fOutput->Add(fNEventsVtx);

  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah) AliFatal("No AOD output handler set in analysis manager");

  ah->AddBranch("TH2D", &fHist);
  
  // Post data for ALL output slots >0 here, to get at least an empty histogram
  PostData(1, fOutput); 
}

//____________________________________________________________________
inline UShort_t
CentralMultTask::CheckEvent(const AliESDEvent& esd, Double_t& vz)
{
  // Do some fast cuts first
  vz = 0;

  // Get the analysis manager - should always be there 
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  if (!am) { 
    AliWarning("No analysis manager defined!");
    return 0;
  }

  // Get the input handler - should always be there 
  AliInputEventHandler* ih = 
    static_cast<AliInputEventHandler*>(am->GetInputEventHandler());
  if (!ih) { 
    AliWarning("No input handler");
    return 0;
  }

  // Trigger mask 
  UInt_t    mask     = ih->IsEventSelected();   
  Bool_t   isMinBias = (mask == AliVEvent::kMB) ? 1 : 0;
  UShort_t ret       = 0;
  if (isMinBias) ret |= kValidTrigger;
  
  // check for good reconstructed vertex
  if (!(esd.GetPrimaryVertex()->GetStatus())) {
    if (fDebug)
      AliWarning("Primary vertex has bad status");
    return ret;
  }
  if (!(esd.GetPrimaryVertexSPD()->GetStatus())) { 
    if (fDebug)
      AliWarning("Primary SPD vertex has bad status");
    return ret;
  }

  // if vertex is from spd vertexZ, require more stringent cut
  if (esd.GetPrimaryVertex()->IsFromVertexerZ()) {
    if (esd.GetPrimaryVertex()->GetDispersion() > 0.02 ||  
	esd.GetPrimaryVertex()->GetZRes()       > 0.25) {
      if (fDebug)
	AliWarning(Form("Primary vertex dispersion=%f (0.02) zres=%f (0.05)", 
			esd.GetPrimaryVertex()->GetDispersion(),
			esd.GetPrimaryVertex()->GetZRes()));
      return ret;
    }
  }
  // One can add a cut on the vertex z position here
  vz = esd.GetPrimaryVertex()->GetZ();
  Double_t vl = fNEventsVtx->GetXaxis()->GetXmin();
  Double_t vh = fNEventsVtx->GetXaxis()->GetXmax();
  if (vz <  vl || vz > vh) {
    if (fDebug)
      AliWarning(Form("Primary vertex vz=%f out of range [%f,%f]", vz, vl, vh));
    return ret; 
  }
  ret |= kValidVertex;

  return ret;
}

//____________________________________________________________________
inline void 
CentralMultTask::UserExec(Option_t *) 
{
  // Main loop
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliError("Cannot get the ESD event");
    return;
  }  

  fHist->Reset();

  // Check event 
  Double_t vz         = 0;
  UShort_t eventFlags = CheckEvent(*esd, vz);
  if (!(eventFlags & kValidTrigger)) 
    return; // No trigger
  else {
    fNEventsTr->Fill(vz);
    if (eventFlags & kValidVertex) 
      fNEventsVtx->Fill(vz);
    else 
      return; // No trigger or no vertex
  }
  

  // flags for secondary and rejected tracks
  // set this bit in ESD tracks if it is rejected by a cut
  const int kRejBit = BIT(15); 
  // set this bit in ESD tracks if it is secondary according to a cut
  const int kSecBit = BIT(16); 

  Int_t total      = 0;
  Int_t nESDTracks = esd->GetNumberOfTracks();
  // Loop over ESD tracks 
  for(Int_t i = 0; i < nESDTracks; i++){ // flag the tracks

    AliESDtrack* track = esd->GetTrack(i);
    
    // if track is a secondary from a V0, flag as a secondary
    if (track->IsOn(AliESDtrack::kMultInV0)) {
      track->SetBit(kSecBit); 
      continue;
    } 

    Double_t eta = track->Eta();
    Double_t phi = track->Phi();
    // check tracks with ITS part
    if (track->IsOn(AliESDtrack::kITSin)) {

      // Check ITS pure stand-alone - and if it is, reject it 
      if (track->IsOn(AliESDtrack::kITSpureSA)) {
	track->SetBit(kRejBit);
	continue;
      }

      // Check DCA - if not close enough reject it as secondary 
      if (!fDCAwSPD.AcceptTrack(track) && 
	  !fDCAwoSPD.AcceptTrack(track)) {
	track->SetBit(kSecBit); 
	continue;
      }

      // Check if this is an ITS complementary - no TPC in
      if (!track->IsOn(AliESDtrack::kTPCin)) { 
	if (fQITS.AcceptTrack(track)) { // Check ITS quality 
	  fITS->Fill(eta);
	  fAll->Fill(eta);
	}
	else 
	  track->SetBit(kRejBit);
      }
      else { // Not ITS SA or TPC+ITS
	if (fQGlo.AcceptTrack(track)) { // Check ITS quality 
	  fGlobal->Fill(eta);
	  fAll->Fill(eta);
	}
	else 
	  track->SetBit(kRejBit);
      }
      if (track->IsOn(kSecBit) || track->IsOn(kRejBit)) continue;

      // fill d2n/detadphi
      fHist->Fill(eta, phi);
    }
  }
  
  // Get multiplicity from SPD tracklets 
  const AliMultiplicity* spdmult = esd->GetMultiplicity();
  for (Int_t i=0; i<spdmult->GetNumberOfTracklets(); ++i){
    // if counting tracks+tracklets, 
    // check if clusters were already used in tracks
    Int_t id1,id2;
    spdmult->GetTrackletTrackIDs(i,0,id1,id2);
    AliESDtrack* tr1 = id1>=0 ? esd->GetTrack(id1) : 0;
    AliESDtrack* tr2 = id2>=0 ? esd->GetTrack(id2) : 0;
    //
    if ((tr1 && tr1->TestBit(kSecBit))  || // Flagged as secondary
	(tr2 && tr2->TestBit(kSecBit))  || // Flagged as secondary
	(tr1 && !tr1->TestBit(kRejBit)) || // already accounted for
	(tr2 && !tr2->TestBit(kRejBit)))   // already accounted for 
      continue; 
    ++total;
    Double_t eta = spdmult->GetEta(i);
    Double_t phi = spdmult->GetPhi(i);
    
    // Increment d2n/detadphi from an SPD tracklet 
    if (fUseTracklets) {
      fHist->Fill(eta, phi);
      fAll->Fill(eta);
      total++;
    }
    fTracklets->Fill(eta);
  }
  if (fDebug) AliInfo(Form("A total of %d tracks", total));

  // NEW HISTO should be filled before this point, as PostData puts the
  // information for this iteration of the UserExec in the container
  PostData(1, fOutput);
}

  
//________________________________________________________________________
inline void 
CentralMultTask::Terminate(Option_t *) 
{
  // Draw result to screen, or perform fitting, normalizations Called
  // once at the end of the query
        
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if(!fOutput) {
    AliError("Could not retrieve TList fOutput"); 
    return; 
  }

  TH1D* all       = static_cast<TH1D*>(fOutput->FindObject("all"));
  TH1D* global    = static_cast<TH1D*>(fOutput->FindObject("global"));
  TH1D* its       = static_cast<TH1D*>(fOutput->FindObject("its"));
  TH1D* tracklets = static_cast<TH1D*>(fOutput->FindObject("tracklets"));
  TH1D* eventsTr  = static_cast<TH1D*>(fOutput->FindObject("nEventsTr"));
  TH1D* eventsVtx = static_cast<TH1D*>(fOutput->FindObject("nEventsVtx"));

  Int_t nTriggers  = eventsTr->GetEntries();
  Int_t nVertex    = eventsVtx->GetEntries();
  if (nTriggers <= 0 || nVertex <= 0) {
    AliWarning("No data in the events histogram!");
    nTriggers = 1;
  }

  all      ->Scale(1. / nTriggers, "width");
  global   ->Scale(1. / nTriggers, "width");
  its      ->Scale(1. / nTriggers, "width");
  tracklets->Scale(1. / nTriggers, "width");

  THStack* stack = new THStack("components", "Components");
  if (fUseTracklets) stack->Add(tracklets);
  stack->Add(global);
  stack->Add(its);

  fOutput->Add(stack);
}

//========================================================================
inline AliAnalysisTask*
AddTaskCentralTracks()
{
  // analysis manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  
  // Make our object.  2nd argumenent is absolute max Eta 
  // 3rd argument is absolute max Vz
  CentralMultTask* task = new CentralMultTask("Global", 2, 10);
  // if physics selection performed in UserExec(),
  // this line should be commented 
  // task->SelectCollisionCandidates(AliVEvent::kMB); 
  mgr->AddTask(task);

  // create containers for input/output
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("Central", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  
  // connect input/output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);

  return task;
}

  
//________________________________________________________________________
//
// EOF
// 
