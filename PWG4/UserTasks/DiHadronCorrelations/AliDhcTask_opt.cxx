// Dihadron correlations task - simple task to read ESD or AOD input,
// calculate same- and mixed-event correlations, and fill THnSparse
// output. -A. Adare, Apr 2011

#include "TCanvas.h"
#include "TChain.h"
#include "TFormula.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "THnSparse.h"
#include "TROOT.h"
#include "TTree.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliDhcTask_opt.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "KiddiePoolClasses.h"
#include "AliVParticle.h"

ClassImp(AliDhcTask)

//________________________________________________________________________
AliDhcTask::AliDhcTask(const char *name) 
: AliAnalysisTaskSE(name), fVerbosity(0), fEtaMax(1), fZVtxMax(10), fPtMin(0.25), fPtMax(15), 
  fESD(0), fAOD(0), fOutputList(0), fHistPt(0), fHEvt(0), fHTrk(0), fHSs(0), fHMs(0), 
  fIndex(0), fMeanPtTrg(0), fMeanPtAss(0), fMean2PtTrg(0), fMean2PtAss(0),
  fCentrality(99), fZVertex(99), fEsdTrackCutsTPCOnly(0), fPoolMgr(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());

  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.,SPDVertex.,TPCVertex.,Tracks "
               "AOD:header,tracks,vertices,";
}

//________________________________________________________________________
void AliDhcTask::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (per slave on PROOF!)

  fOutputList = new TList();
  fOutputList->SetOwner(1);

  fEsdTrackCutsTPCOnly = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  //fEsdTrackCutsTPCOnly->SetMinNClustersTPC(70);
  fEsdTrackCutsTPCOnly->SetMinNCrossedRowsTPC(70);
  fEsdTrackCutsTPCOnly->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);

  BookHistos();
  InitEventMixer(); 
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliDhcTask::BookHistos()
{
  // Book histograms.

  Int_t nDeta=20, nPtAssc=12, nPtTrig=12, nCent=12, nDphi=36, nZvtx=8;

  Double_t ptt[]  = {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 15};
  Double_t pta[]  = {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 15};
  Double_t cent[] = {0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 90};
  Double_t zvtx[] = {-10, -6, -4, -2, 0, 2, 4, 6, 10};

  // Event histo
  fHEvt = new TH2F("fHEvt", "Event-level variables; Zvtx; Cent", 30, -15, 15, 101, 0, 101);
  fOutputList->Add(fHEvt);
  // Track histo
  fHTrk = new TH2F("fHTrk", "Track-level variables", 
		   100, 0, TMath::TwoPi(), 100, -fEtaMax, +fEtaMax);
  fOutputList->Add(fHTrk);
  
  // Left over from the tutorial :)
  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 200, 0., fPtMax);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);
  fOutputList->Add(fHistPt);

  fHPtAss = new TH1F("fHPtAss","PtAssoc;P_{T} (GeV/c) [GeV/c]",nPtAssc,pta);
  fOutputList->Add(fHPtAss);
  fHPtTrg = new TH1F("fHPtTrg","PtTrig;P_{T} (GeV/c) [GeV/c]",nPtTrig,ptt);
  fOutputList->Add(fHPtTrg);
  fHCent = new TH1F("fHCent","Cent;bins",nCent,cent);
  fOutputList->Add(fHCent);
  fHZvtx = new TH1F("fHZvtx","Zvertex;bins",nZvtx,zvtx);
  fOutputList->Add(fHZvtx);

  fNbins = nPtTrig*nPtAssc*nCent*nZvtx;
  fHSs = new TH2*[fNbins];
  fHMs = new TH2*[fNbins];

  fIndex = new TFormula("GlobIndex","(t-1)*[0]*[1]*[2]+(z-1)*[0]*[1]+(x-1)*[0]+(y-1)+0*[4]");
  fIndex->SetParameters(nPtTrig,nPtAssc,nZvtx,nCent);
  fIndex->SetParNames("NTrigBins","NAssocBins", "NZvertexBins", "NCentBins");
  //fOutputList->Add(fIndex);

  fMeanPtTrg = new TProfile2D*[nPtTrig*nPtAssc];
  fMeanPtAss = new TProfile2D*[nPtTrig*nPtAssc];
  fMean2PtTrg = new TProfile2D*[nPtTrig*nPtAssc];
  fMean2PtAss = new TProfile2D*[nPtTrig*nPtAssc];
  for (Int_t c=1; c<=nCent; ++c) {
    TString title(Form("cen=%d (%.1f)",c,fHCent->GetBinCenter(c)));
    fMeanPtTrg[c-1]  = new TProfile2D(Form("hMeanPtTrgCen%d",c),title,nPtTrig,ptt,nPtAssc,pta);
    fMeanPtAss[c-1]  = new TProfile2D(Form("hMeanPtAssCen%d",c),title,nPtTrig,ptt,nPtAssc,pta);
    fMean2PtTrg[c-1] = new TProfile2D(Form("hMean2PtTrgCen%d",c),title,nPtTrig,ptt,nPtAssc,pta);
    fMean2PtAss[c-1] = new TProfile2D(Form("hMean2PtAssCen%d",c),title,nPtTrig,ptt,nPtAssc,pta);
    fOutputList->Add(fMeanPtTrg[c-1]);
    fOutputList->Add(fMeanPtAss[c-1]);
    fOutputList->Add(fMean2PtTrg[c-1]);
    fOutputList->Add(fMean2PtAss[c-1]);
  }

  Int_t count = 0;
  for (Int_t c=1; c<=nCent; ++c) {
    for (Int_t z=1; z<=nZvtx; ++z) {
      for (Int_t t=1; t<=nPtTrig; ++t) {
	for (Int_t a=1; a<=nPtAssc; ++a) {
	  fHSs[count] = 0;
	  fHMs[count] = 0;
	  if (a>t) {
	    ++count;
	    continue;
	  }
	  TString title(Form("cen=%d (%.1f), zVtx=%d (%.1f), trig=%d (%.1f), assc=%d (%.1f)",
			     c,fHCent->GetBinCenter(c), z,fHZvtx->GetBinCenter(z),
			     t,fHPtTrg->GetBinCenter(t),a, fHPtAss->GetBinCenter(a)));
	  fHSs[count] = new TH2F(Form("hS%d",count), Form("Signal %s",title.Data()),
				 nDphi,-0.5*TMath::Pi(),1.5*TMath::Pi(),nDeta,-2*fEtaMax,2*fEtaMax);
	  fHMs[count] = new TH2F(Form("hM%d",count), Form("Signal %s",title.Data()),
				 nDphi,-0.5*TMath::Pi(),1.5*TMath::Pi(),nDeta,-2*fEtaMax,2*fEtaMax);
	  fOutputList->Add(fHSs[count]);
	  fOutputList->Add(fHMs[count]);
	  if (fVerbosity>5)
	    cout << count << " " << fIndex->Eval(t,a,z,c) << ": " << title << endl;
	  ++count;
	}
      }
    }
  }

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
  Int_t nCentBins  = 12;
  Double_t centBins[] = {0,1,2,3,4,5,10,20,30,40,50,60,90.1};
 
  //Int_t nCentBins  = 1;
  //Double_t centBins[] = {-1,100.1};
 
  // Z-vertex pools
  Int_t nZvtxBins  = 8;
  Double_t zvtxbin[] = {-10,-6,-4,-2,0,2,4,6,10};

  fPoolMgr = new KiddiePoolManager();
  fPoolMgr->SetTargetTrackDepth(trackDepth);
  if (fVerbosity>4)
    fPoolMgr->SetDebug(1);
  fPoolMgr->InitEventPools(poolsize, nCentBins, centBins, nZvtxBins, zvtxbin);
  
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
    AliWarning(Form("No pool found. Centrality %f, ZVertex %f", fCentrality, fZVertex));
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

  pool->UpdatePool(sTracks);
  PostData(1, fOutputList);
  return;
}

//________________________________________________________________________
MiniEvent* AliDhcTask::GetESDTrax() const
{
  // Loop twice: 1. Count sel. tracks. 2. Fill vector.

  const AliESDVertex *vtxSPD = fESD->GetPrimaryVertexSPD();
  if (!vtxSPD)
    return 0;

  Int_t nTrax = fESD->GetNumberOfTracks();
  if (fVerbosity > 2)
    AliInfo(Form("%d tracks in event",nTrax));

  // Loop 1.
  Int_t nSelTrax = 0;
  TObjArray arr(nTrax);
  arr.SetOwner(1);

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
    if (TMath::Abs(eta) > fEtaMax)
      continue;

    // create a tpc only track
    AliESDtrack *newtrack = AliESDtrackCuts::GetTPCOnlyTrack(fESD,esdtrack->GetID());
    if(!newtrack)
      continue;
    if (newtrack->Pt()<=0) {
      delete newtrack;
      continue;
    }

    AliExternalTrackParam exParam;
    Bool_t relate = newtrack->RelateToVertexTPC(vtxSPD,fESD->GetMagneticField(),kVeryBig,&exParam);
    if (!relate) {
      delete newtrack;
      continue;
    }

    // set the constraint parameters to the track
    newtrack->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),exParam.GetCovariance());

    pt = newtrack->Pt();
    ptOK = pt >= fPtMin && pt < fPtMax;
    if (!ptOK) {
      delete newtrack;
      continue;
    }
    eta  = esdtrack->Eta();
    if (TMath::Abs(eta) > fEtaMax) {
      delete newtrack;
      continue;
    }
    arr.Add(newtrack);
    nSelTrax++;
  }

  MiniEvent* miniEvt = new MiniEvent(0);
  miniEvt->reserve(nSelTrax);

  // Loop 2.
  for (Int_t i = 0; i < nSelTrax; ++i) {
    AliESDtrack* esdtrack = static_cast<AliESDtrack*>(arr.At(i));
    if (!esdtrack) {
      AliError(Form("Couldn't get ESD track %d\n", i));
      continue;
    }
    Double_t pt = esdtrack->Pt();
    Double_t eta  = esdtrack->Eta();
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
    if (TMath::Abs(eta) > fEtaMax)
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
    if (TMath::Abs(eta) > fEtaMax)
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
			    Int_t pairing, Double_t /*weight*/)
{
  // Triggered angular correlations. If pairing is kSameEvt, particles
  // within evt1 are correlated. If kDiffEvt, correlate triggers from
  // evt1 with partners from evt2.

  Int_t cbin = fHCent->FindBin(fCentrality);
  if (fHCent->IsBinOverflow(cbin) ||
      fHCent->IsBinUnderflow(cbin))
    return 0;

  Int_t zbin = fHZvtx->FindBin(fZVertex);
  if (fHZvtx->IsBinOverflow(zbin) ||
      fHZvtx->IsBinUnderflow(zbin))
    return 0;

  Int_t iMax = evt1.size();
  Int_t jMax = evt2.size();

  TH2  **hist = fHMs;
  if (pairing == kSameEvt) {
    hist = fHSs;
    fHCent->AddBinContent(cbin);
    fHZvtx->AddBinContent(zbin);
  }

  Int_t nZvtx = fHZvtx->GetNbinsX();
  Int_t nPtTrig = fHPtTrg->GetNbinsX();
  Int_t nPtAssc = fHPtAss->GetNbinsX();

  Int_t globIndex = (cbin-1)*nZvtx*nPtTrig*nPtAssc+(zbin-1)*nPtTrig*nPtAssc;

  for (Int_t i=0; i<iMax; ++i) {

    // Trigger particles
    const MiniTrack &a(evt1.at(i));

    Float_t pta  = a.Pt();
    Int_t abin = fHPtTrg->FindBin(pta);
    if (fHPtTrg->IsBinOverflow(abin) ||
	fHPtTrg->IsBinUnderflow(abin))
      continue;

    if (pairing == kSameEvt) {
      fHistPt->Fill(pta);
      fHTrk->Fill(a.Phi(),a.Eta());
      fHPtTrg->AddBinContent(abin);
    }

    for (Int_t j=0; j<jMax; ++j) {
      // Associated particles
      if (pairing == kSameEvt && i==j)
	continue;

      const MiniTrack &b(evt2.at(j));
      
      Float_t ptb  = b.Pt();
      if (pta < ptb) 
	continue;

      Int_t bbin = fHPtTrg->FindBin(ptb);
      if (fHPtAss->IsBinOverflow(bbin) ||
	  fHPtAss->IsBinUnderflow(bbin))
	continue;

      Float_t dphi = DeltaPhi(a.Phi(), b.Phi());
      Float_t deta = a.Eta() - b.Eta();

      Int_t index = globIndex+(abin-1)*nPtAssc+(bbin-1);
      hist[index]->Fill(dphi,deta);

      if (pairing == kSameEvt) {
	fHPtAss->AddBinContent(bbin);
	fMeanPtTrg[cbin-1]->Fill(pta,ptb,pta);
	fMeanPtAss[cbin-1]->Fill(pta,ptb,ptb);
	fMean2PtTrg[cbin-1]->Fill(pta,ptb,pta*pta);
	fMean2PtAss[cbin-1]->Fill(pta,ptb,ptb*ptb);
      }
    }
  }

  return 0;
}

//________________________________________________________________________
void AliDhcTask::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  delete fPoolMgr;

  fHCent->SetEntries(fHCent->Integral());
  fHZvtx->SetEntries(fHZvtx->Integral());
  fHPtTrg->SetEntries(fHPtTrg->Integral());
  fHPtAss->SetEntries(fHPtAss->Integral());

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
  if (name.CompareTo("TPCVertex")==0)
    return kFALSE;
  
  // Check # contributors and range...
  if( nContributors < 1 || TMath::Abs(zVertex) > fZVtxMax ) {
    return kFALSE;
  }
  
  return kTRUE;
}
