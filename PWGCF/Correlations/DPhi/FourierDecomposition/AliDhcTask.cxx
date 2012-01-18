// Dihadron correlations task - simple task to read ESD or AOD input,
// calculate same- and mixed-event correlations, and fill THnSparse
// output. -A. Adare, Apr 2011

#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TROOT.h"
#include "TTree.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliDhcTask.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "KiddiePoolClasses.h"
#include "AliVParticle.h"

ClassImp(AliDhcTask)

//________________________________________________________________________
AliDhcTask::AliDhcTask(const char *name) 
: AliAnalysisTaskSE(name), fVerbosity(0), fESD(0), fAOD(0), fOutputList(0), 
  fHistPt(0), fHEvt(0), fHTrk(0), fHS(0), fHM(0), fPoolMgr(0), 
  fCentrality(99), fZVertex(99), fZVtxMax(10), fPtMin(0), fPtMax(100), 
  fInputHandler(0), fEsdTrackCutsTPCOnly(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());

  // Initialize cut variables
  fZVtxMax = 10;    // cm
  fPtMin   = 0.50;  // GeV/c
  fPtMax   = 15.;

  /*
  // Not used for anything now...
  fInputHandler = (AliInputEventHandler*)
    ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if (!fInputHandler) {
    Error("AliDhcTask()", "Did not get input handler");
  }
  */

  fEsdTrackCutsTPCOnly = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.,SPDVertex.,TPCVertex.,Tracks "
               "AOD:header,tracks,vertices,";
}

//________________________________________________________________________
void AliDhcTask::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (per slave on PROOF!)

  BookHistos();
  InitEventMixer(); 

  fOutputList = new TList();
  fOutputList->SetOwner(1);

  fOutputList->Add(fHistPt);
  fOutputList->Add(fHS);
  fOutputList->Add(fHM);
  fOutputList->Add(fHEvt);
  fOutputList->Add(fHTrk);

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliDhcTask::BookHistos()
{
  // Setup for THnSparse correlation histos.
  // Important! Order bins[] according to ePairHistAxes.

  const Int_t ndims = 6;
  Int_t nDeta=22, nPtAssc=11, nPtTrig=11, nCent=8, nDphi=36, nZvtx=8;
  Int_t bins[ndims] = {nDeta, nPtAssc, nPtTrig, nCent, nDphi, nZvtx };
  Double_t xmin[ndims] = { -2.2, 0.5, 0.5, 0,  -0.5*TMath::Pi(), -10 };
  Double_t xmax[ndims] = { +2.2, 15., 15., 90, +1.5*TMath::Pi(), +10 };
  fHS = new THnSparseF("fHS", "Same evt Correlations", 
		       ndims, bins, xmin, xmax);

  // Override the uniform binning for some axes. Keep xmin and xmax
  // limits the same, but allow variable bin widths inside.
  Double_t ptt[]  = {0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 10, 15};
  Double_t pta[]  = {0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 10, 15};
  Double_t cent[] = {0, 2, 10, 20, 30, 40, 50, 60, 90};
  Double_t zvtx[] = {-10, -6, -4, -2, 0, 2, 4, 6, 10};
  fHS->SetBinEdges(kPtTrig, ptt);
  fHS->SetBinEdges(kPtAssc, pta);
  fHS->SetBinEdges(kCent, cent);
  fHS->SetBinEdges(kZvtx, zvtx); // Match InitEventMixer - no point in going finer

  fHM = (THnSparse*) fHS->Clone("fHM");
  fHM->SetTitle("Mixed evt Correlations");

  // Event histo
  fHEvt = new TH2F("fHEvt", "Event-level variables", 30, -15, 15, 101, 0, 101);
  // Track histo
  fHTrk = new TH2F("fHTrk", "Track-level variables", 
		   100, 0, TMath::TwoPi(), 100, -2, +2);

  // Left over from the tutorial :)
  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 200, 0., 20.);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);
  
  return;
}

//________________________________________________________________________
void AliDhcTask::InitEventMixer()
{
  // The effective pool size in events is set by trackDepth, so more
  // low-mult events are required to maintain the threshold than
  // high-mult events. Centrality pools are indep. of data histogram
  // binning, no need to match.

  Int_t trackDepth = 1000;   // # tracks to fill pool
  Int_t poolsize   = 200;    // Maximum number of events

  // Centrality pools
  Int_t nCentBins  = 9;
  Double_t centBins[] = {0,1,2,5,10,20,30,40,60,90.1};
 
  //Int_t nCentBins  = 1;
  //Double_t centBins[] = {-1,100.1};
 
  // Z-vertex pools
  Int_t nZvtxBins  = 8;
  Double_t zvtxbin[] = {-10, -6, -4, -2, 0, 2, 4, 6, 10};
  
  fPoolMgr = new KiddiePoolManager(poolsize, trackDepth, nCentBins, 
				   centBins, nZvtxBins, zvtxbin);
  return;
}

//________________________________________________________________________
void AliDhcTask::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  static int iEvent = -1; ++iEvent;

  if (fVerbosity>2) {
    if (iEvent % 10 == 0) 
      cout << iEvent << endl;
  }

  Int_t dType = -1;       // Will be set to kESD or kAOD.
  MiniEvent* sTracks = 0; // Vector of selected MiniTracks.
  Double_t centCL1 = -1;

  LoadBranches();

  // Get event pointers, check for signs of life
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (fESD)
    dType = kESD;
  else if (fAOD)
    dType = kAOD;
  else {
    AliError("Neither fESD nor fAOD available");
    return;
  }

  // Centrality, vertex, other event variables...
  if (dType == kESD) {
    if (!VertexOk(fESD)) {
      if (fVerbosity > 1)
	AliInfo(Form("Event REJECTED. z = %.1f", fZVertex));
      return;
    }
    const AliESDVertex* vertex = fESD->GetPrimaryVertex();
    fZVertex = vertex->GetZ();
    if(fESD->GetCentrality()) {
      fCentrality = 
	fESD->GetCentrality()->GetCentralityPercentile("V0M");
      centCL1 =
	fESD->GetCentrality()->GetCentralityPercentile("CL1");
    }
  }
  if (dType == kAOD) {
    const AliAODVertex* vertex = fAOD->GetPrimaryVertex();
    fZVertex = vertex->GetZ();
    if (!VertexOk(fAOD)) {
      if (fVerbosity > 1)
	Info("Exec()", "Event REJECTED. z = %.1f", fZVertex);
      return;
    }
    const AliCentrality *aodCent = fAOD->GetHeader()->GetCentralityP();
    if (aodCent) {
      fCentrality = aodCent->GetCentralityPercentile("V0M");
      centCL1     = aodCent->GetCentralityPercentile("CL1");
    }
  }
  
  // Fill Event histogram
  fHEvt->Fill(fZVertex, fCentrality);

  if (fCentrality > 90. || fCentrality < 0) {
    AliInfo(Form("Event REJECTED. fCentrality = %.1f", fCentrality));
    return;
  }

  // Get array of selected tracks
  if (dType == kESD) {
    sTracks = GetESDTrax();
  }
  if (dType == kAOD) {
    sTracks = GetAODTrax();
  }

  // Get pool containing tracks from other events like this one
  KiddiePool* pool = fPoolMgr->GetEventPool(fCentrality, fZVertex);
  if (!pool) {
    AliFatal(Form("No pool found. Centrality %f, ZVertex %f", fCentrality, fZVertex));
    return;
  }

  if (!pool->IsReady()) {
    pool->UpdatePool(sTracks);
    return;
  }

  if (pool->IsFirstReady()) {
    // recover events missed before the pool is ready
    Int_t nEvs = pool->GetCurrentNEvents();
    if (nEvs>1) {
      for (Int_t i=0; i<nEvs; ++i) {
	MiniEvent* evI = pool->GetEvent(i);
	for (Int_t j=0; j<nEvs; ++j) {
	  MiniEvent* evJ = pool->GetEvent(j);
	  if (i==j) {
	    Correlate(*evI, *evJ, kSameEvt);
	  } else {
	    Correlate(*evI, *evJ, kDiffEvt, 1.0);
	  }
	}
      }
    }
  } else { /* standard case: same event, then mix*/
    Correlate(*sTracks, *sTracks, kSameEvt);  
    Int_t nMix = pool->GetCurrentNEvents();
    for (Int_t jMix=0; jMix<nMix; ++jMix) {
      MiniEvent* bgTracks = pool->GetEvent(jMix);
      Correlate(*sTracks, *bgTracks, kDiffEvt, 1.0);
    }
  }

  if (fVerbosity>4) {
    cout << "Output of SAME  THnSparse: " << fHS->GetSparseFractionBins() << " " << fHS->GetSparseFractionMem() << endl; 
    cout << "Output of MIXED THnSparse: " << fHM->GetSparseFractionBins() << " " << fHM->GetSparseFractionMem() << endl; 
  }

  pool->UpdatePool(sTracks);
  PostData(1, fOutputList);
  return;
}

//________________________________________________________________________
MiniEvent* AliDhcTask::GetESDTrax() const
{
  // Loop twice: 1. Count sel. tracks. 2. Fill vector.

  Int_t nTrax = fESD->GetNumberOfTracks();
  Int_t nSelTrax = 0;

  if (fVerbosity > 2)
    AliInfo(Form("%d tracks in event",nTrax));

  // Loop 1.
  for (Int_t i = 0; i < nTrax; ++i) {
    AliESDtrack* esdtrack = fESD->GetTrack(i);
    if (!esdtrack) {
      AliError(Form("Couldn't get ESD track %d\n", i));
      continue;
    }
    Bool_t trkOK = fEsdTrackCutsTPCOnly->AcceptTrack(esdtrack);
    if (!trkOK)
      continue;
    Double_t pt = esdtrack->Pt();
    Bool_t ptOK = pt >= fPtMin && pt < fPtMax;
    if (!ptOK)
      continue;
    Double_t eta = esdtrack->Eta();
    if (TMath::Abs(eta) > 1.0)
      continue;
    nSelTrax++;
  }

  MiniEvent* miniEvt = new MiniEvent(0);
  miniEvt->reserve(nSelTrax);

  // Loop 2.  
  for (Int_t i = 0; i < nTrax; ++i) {
    AliESDtrack* esdtrack = fESD->GetTrack(i);
    if (!esdtrack) {
      AliError(Form("Couldn't get ESD track %d\n", i));
      continue;
    }
    Bool_t trkOK = fEsdTrackCutsTPCOnly->AcceptTrack(esdtrack);
    if (!trkOK)
      continue;
    Double_t pt = esdtrack->Pt();
    Bool_t ptOK = pt >= fPtMin && pt < fPtMax;
    if (!ptOK)
      continue;
    Double_t eta  = esdtrack->Eta();
    if (TMath::Abs(eta) > 1.0)
      continue;

    Double_t phi  = esdtrack->Phi();
    Int_t    sign = esdtrack->Charge() > 0 ? 1 : -1;
    miniEvt->push_back(MiniTrack(pt, eta, phi, sign));
  }
  return miniEvt;
}

//________________________________________________________________________
MiniEvent* AliDhcTask::GetAODTrax() const
{
  // Loop twice: 1. Count sel. tracks. 2. Fill vector.

  Int_t nTrax = fAOD->GetNumberOfTracks();
  Int_t nSelTrax = 0;

  if (fVerbosity > 2)
    AliInfo(Form("%d tracks in event",nTrax));

  // Loop 1.
  for (Int_t i = 0; i < nTrax; ++i) {
    AliAODTrack* aodtrack = fAOD->GetTrack(i);
    if (!aodtrack) {
      AliError(Form("Couldn't get AOD track %d\n", i));
      continue;
    }
    // See $ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C
    UInt_t tpcOnly = 1 << 7;
    Bool_t trkOK = aodtrack->TestFilterBit(tpcOnly);
    if (!trkOK)
      continue;
    Double_t pt = aodtrack->Pt();
    Bool_t ptOK = pt >= fPtMin && pt < fPtMax;
    if (!ptOK)
      continue;
    Double_t eta = aodtrack->Eta();
    if (TMath::Abs(eta) > 1.0)
      continue;
    nSelTrax++;
  }

  MiniEvent* miniEvt = new MiniEvent(0);
  miniEvt->reserve(nSelTrax);

  // Loop 2.  
  for (Int_t i = 0; i < nTrax; ++i) {
    AliAODTrack* aodtrack = fAOD->GetTrack(i);
    if (!aodtrack) {
      AliError(Form("Couldn't get AOD track %d\n", i));
      continue;
    }
    
    // See $ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C
    UInt_t tpcOnly = 1 << 7;
    Bool_t trkOK = aodtrack->TestFilterBit(tpcOnly);
    if (!trkOK)
      continue;
    Double_t pt = aodtrack->Pt();
    Bool_t ptOK = pt >= fPtMin && pt < fPtMax;
    if (!ptOK)
      continue;
    Double_t eta  = aodtrack->Eta();
    if (TMath::Abs(eta) > 1.0)
      continue;

    Double_t phi  = aodtrack->Phi();
    Int_t    sign = aodtrack->Charge() > 0 ? 1 : -1;
    miniEvt->push_back(MiniTrack(pt, eta, phi, sign));
  }
  return miniEvt;
}

//________________________________________________________________________
Double_t AliDhcTask::DeltaPhi(Double_t phia, Double_t phib,
			      Double_t rangeMin, Double_t rangeMax) const
{
  Double_t dphi = -999;
  Double_t pi = TMath::Pi();
  
  if (phia < 0)         phia += 2*pi;
  else if (phia > 2*pi) phia -= 2*pi;
  if (phib < 0)         phib += 2*pi;
  else if (phib > 2*pi) phib -= 2*pi;
  dphi = phib - phia;
  if (dphi < rangeMin)      dphi += 2*pi;
  else if (dphi > rangeMax) dphi -= 2*pi;
  
  return dphi;
}

//________________________________________________________________________
Int_t AliDhcTask::Correlate(const MiniEvent &evt1, const MiniEvent &evt2, 
			    Int_t pairing, Double_t weight)
{
  // Triggered angular correlations. If pairing is kSameEvt, particles
  // within evt1 are correlated. If kDiffEvt, correlate triggers from
  // evt1 with partners from evt2.
  
  Int_t iMax = evt1.size();
  Int_t jMax = evt2.size();

  THnSparse *hist = fHM;
  if (pairing == kSameEvt)
    hist = fHS;

  for (Int_t i=0; i<iMax; ++i) {

    // Trigger particles
    const MiniTrack &a(evt1.at(i));

    if (pairing == kSameEvt) {
      fHistPt->Fill(a.Pt());
      fHTrk->Fill(a.Phi(),a.Eta());
    }

    for (int j=0; j<jMax; ++j) {
      // Associated particles
      if (pairing == kSameEvt && i==j)
	continue;

      const MiniTrack &b(evt2.at(j));
      
      Float_t pta  = a.Pt();
      Float_t ptb  = b.Pt();
      if (pta < ptb) 
	continue;

      Float_t dphi = DeltaPhi(a.Phi(), b.Phi());
      Float_t deta = a.Eta() - b.Eta();

      Double_t x[6]; // Match ndims in fHS
      x[kDeta]       = deta;
      x[kPtAssc]     = ptb;
      x[kPtTrig]     = pta;
      x[kCent]       = fCentrality;
      x[kDphi]       = dphi;
      x[kZvtx]       = fZVertex;
      //x[kChargeComb] = a.Sign() * b.Sign(); 
      
      hist->Fill(x, weight);
    }
  }

  return 0;
}

//________________________________________________________________________
void AliDhcTask::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  if (gROOT->IsBatch())
    return;

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    AliError("Output list not available");
    return;
  }
  
  fHistPt = dynamic_cast<TH1F*> (fOutputList->At(0));
  if (!fHistPt) {
    AliError("ERROR: fHistPt not available\n");
    return;
  }
   
  TCanvas *c1 = new TCanvas("AliDhcTask","Pt",10,10,510,510);
  c1->cd(1)->SetLogy();
  fHistPt->DrawCopy("E");
}

//________________________________________________________________________
Bool_t AliDhcTask::VertexOk(TObject* obj) const
{
  // Modified from AliAnalyseLeadingTrackUE::VertexSelection()
  
  Int_t nContributors  = 0;
  Double_t zVertex     = 999;
  TString name("");
  
  if (obj->InheritsFrom("AliESDEvent")) {
    AliESDEvent* esdevt = (AliESDEvent*) obj;
    const AliESDVertex* vtx = esdevt->GetPrimaryVertex();
    if (!vtx)
      return 0;
    nContributors = vtx->GetNContributors();
    zVertex       = vtx->GetZ();
    name          = vtx->GetName();
  }
  else if (obj->InheritsFrom("AliAODEvent")) { 
    AliAODEvent* aodevt = (AliAODEvent*) obj;
    if (aodevt->GetNumberOfVertices() < 1)
      return 0;
    const AliAODVertex* vtx = aodevt->GetPrimaryVertex();
    nContributors = vtx->GetNContributors();
    zVertex       = vtx->GetZ();
    name          = vtx->GetName();
  }
  
  // Reject if TPC-only vertex
  if (name.CompareTo("TPCVertex"))
    return kFALSE;
  
  // Check # contributors and range...
  if( nContributors < 1 || TMath::Abs(zVertex) > fZVtxMax ) {
    return kFALSE;
  }
  
  return kTRUE;
}
