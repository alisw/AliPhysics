// Dihadron correlations task - simple task to read ESD or AOD input,
// calculate same- and mixed-event correlations, and fill THnSparse
// output. -A. Adare, Apr 2011

#include "TCanvas.h"
#include "TChain.h"
#include "TFormula.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TProfile2D.h"
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
#include "AliPool.h"
#include "AliVParticle.h"

using std::cout;
using std::endl;

ClassImp(AliDhcTask)

//________________________________________________________________________
AliDhcTask::AliDhcTask(const char *name) 
: AliAnalysisTaskSE(name), fVerbosity(0), fEtaMax(1), fZVtxMax(10), fPtMin(0.25), fPtMax(15),
  fTrackDepth(1000), fPoolSize(200), fTracksName(), fDoWeights(0),
  fESD(0), fAOD(0), fOutputList(0), fHEvt(0), fHTrk(0), fHPtAss(0), fHPtTrg(0), fHPtTrg_Evt(0),
  fHCent(0), fHZvtx(0), fNbins(0), fHSs(0), fHMs(0), fHPts(0), fIndex(0), 
  fCentrality(99), fZVertex(99), fEsdTPCOnly(0), fPoolMgr(0),
  fCentMethod("V0M"), fNBdeta(20), fNBdphi(36),
  fBPtT(0x0), fBPtA(0x0), fBCent(0x0), fBZvtx(0x0),
  fMixBCent(0x0), fMixBZvtx(0x0)
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

  Double_t ptt[4] = {0.25, 1.0, 2.0, 15.0};
  fBPtT  = new TAxis(3,ptt); 
  Double_t pta[4] = {0.25, 1.0, 2.0, 15.0};
  fBPtA  = new TAxis(3,pta); 
  Double_t cent[2] = {-100.0, 100.0};
  fBCent = new TAxis(1,cent);
  Double_t zvtx[2] = {-10, 10};
  fBZvtx = new TAxis(1,zvtx);
  Double_t centmix[2] = {-100.0, 100.0};
  fMixBCent = new TAxis(1,centmix);
  Double_t zvtxmix[9] = {-10,-6,-4,-2,0,2,4,6,10};
  fMixBZvtx = new TAxis(8,zvtxmix);
}

//________________________________________________________________________
void AliDhcTask::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (per slave on PROOF!)

  fOutputList = new TList();
  fOutputList->SetOwner(1);

  fEsdTPCOnly = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  //fEsdTPCOnly->SetMinNClustersTPC(70);
  fEsdTPCOnly->SetMinNCrossedRowsTPC(70);
  fEsdTPCOnly->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);

  BookHistos();
  InitEventMixer(); 
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliDhcTask::BookHistos()
{
  // Book histograms.

  if (fVerbosity > 1) {
    AliInfo(Form("Number of pt(t) bins: %d", fBPtT->GetNbins()));
    for (Int_t i=1; i<=fBPtT->GetNbins(); i++) {
      AliInfo(Form("pt(t) bin %d, %f to %f", i, fBPtT->GetBinLowEdge(i), fBPtT->GetBinUpEdge(i)));
    }
    AliInfo(Form("Number of pt(a) bins: %d", fBPtA->GetNbins()));
    for (Int_t i=1; i<=fBPtA->GetNbins(); i++) {
      AliInfo(Form("pt(a) bin %d, %f to %f", i, fBPtA->GetBinLowEdge(i), fBPtA->GetBinUpEdge(i)));
    }
  }

  Int_t nPtAssc=fBPtA->GetNbins();
  Int_t nPtTrig=fBPtT->GetNbins();
  Int_t nCent=fBCent->GetNbins();
  Int_t nZvtx=fBZvtx->GetNbins();
  Double_t ptt[nPtTrig+1];
  ptt[0] = fBPtT->GetBinLowEdge(1);
  for (Int_t i=1; i<=nPtTrig; i++) {
    ptt[i] = fBPtT->GetBinUpEdge(i);
  }
  Double_t pta[nPtAssc+1];
  pta[0] = fBPtA->GetBinLowEdge(1);
  for (Int_t i=1; i<=nPtAssc; i++) {
    pta[i] = fBPtA->GetBinUpEdge(i);
  }
  Double_t cent[nCent+1];
  cent[0] = fBCent->GetBinLowEdge(1);
  for (Int_t i=1; i<=nCent; i++) {
    cent[i] = fBCent->GetBinUpEdge(i);
  }
  Double_t zvtx[nZvtx+1];
  zvtx[0] = fBZvtx->GetBinLowEdge(1);
  for (Int_t i=1; i<=nZvtx; i++) {
    zvtx[i] = fBZvtx->GetBinUpEdge(i);
  }

  fHEvt = new TH2F("fHEvt", "Event-level variables; Zvtx; Cent", nZvtx, zvtx, nCent, cent);
  fOutputList->Add(fHEvt);
  fHTrk = new TH2F("fHTrk", "Track-level variables", 
		   100, 0, TMath::TwoPi(), 100, -fEtaMax, +fEtaMax);
  fOutputList->Add(fHTrk);
  
  fHPtAss = new TH1F("fHPtAss","PtAssoc;P_{T} (GeV/c) [GeV/c]",nPtAssc,pta);
  fOutputList->Add(fHPtAss);
  fHPtTrg = new TH1F("fHPtTrg","PtTrig;P_{T} (GeV/c) [GeV/c]",nPtTrig,ptt);
  fOutputList->Add(fHPtTrg);
  fHPtTrg_Evt = (TH1*) fHPtTrg->Clone("fHPtTrg_Evt");
  fHCent = new TH1F("fHCent","Cent;bins",nCent,cent);
  fOutputList->Add(fHCent);
  fHZvtx = new TH1F("fHZvtx","Zvertex;bins",nZvtx,zvtx);
  fOutputList->Add(fHZvtx);

  fNbins = nPtTrig*nPtAssc*nCent*nZvtx;
  fHSs   = new TH2*[fNbins];
  fHMs   = new TH2*[fNbins];
  fHPts  = new TH1*[fNbins];

  fIndex = new TFormula("GlobIndex","(t-1)*[0]*[1]*[2]+(z-1)*[0]*[1]+(x-1)*[0]+(y-1)+0*[4]");
  fIndex->SetParameters(nPtTrig,nPtAssc,nZvtx,nCent);
  fIndex->SetParNames("NTrigBins","NAssocBins", "NZvertexBins", "NCentBins");
  //fOutputList->Add(fIndex);

  Int_t count = 0;
  for (Int_t c=1; c<=nCent; ++c) {
    for (Int_t z=1; z<=nZvtx; ++z) {
      for (Int_t t=1; t<=nPtTrig; ++t) {
	for (Int_t a=1; a<=nPtAssc; ++a) {
	  fHSs[count]  = 0;
	  fHMs[count]  = 0;
          fHPts[count] = 0;
	  if (a>t) {
	    ++count;
	    continue;
	  }
          if (t==1 && a==1) {
            TString title(Form("cen=%d (%.1f to %.1f), zVtx=%d (%.1f to %.1f)",
                               c, fBCent->GetBinLowEdge(c), fBCent->GetBinUpEdge(c), 
                               z, fBZvtx->GetBinLowEdge(z), fBZvtx->GetBinUpEdge(z)));
            fHPts[count] = new TH1F(Form("hPt%d",count), Form("Ptdist %s;p_{T} (GeV/c)",title.Data()), 200,0,20);
            fOutputList->Add(fHPts[count]);
          }
	  TString title(Form("cen=%d (%.1f to %.1f), zVtx=%d (%.1f to %.1f), trig=%d (%.1f to %.1f), assoc=%d (%.1f to %.1f)",
			     c, fBCent->GetBinLowEdge(c), fBCent->GetBinUpEdge(c), 
                             z, fBZvtx->GetBinLowEdge(z), fBZvtx->GetBinUpEdge(z),
			     t, fBPtT->GetBinLowEdge(t),  fBPtT->GetBinUpEdge(t),
                             a, fBPtA->GetBinLowEdge(a),  fBPtA->GetBinUpEdge(a)));
	  fHSs[count] = new TH2F(Form("hS%d",count), Form("Signal %s;#Delta #varphi;#Delta #eta",title.Data()),
				 fNBdphi,-0.5*TMath::Pi(),1.5*TMath::Pi(),fNBdeta,-2*fEtaMax,2*fEtaMax);
	  fHMs[count] = new TH2F(Form("hM%d",count), Form("Mixed %s;#Delta #varphi;#Delta #eta",title.Data()),
				 fNBdphi,-0.5*TMath::Pi(),1.5*TMath::Pi(),fNBdeta,-2*fEtaMax,2*fEtaMax);
	  fOutputList->Add(fHSs[count]);
	  fOutputList->Add(fHMs[count]);
          Int_t tcount = (Int_t)fIndex->Eval(t,a,z,c);
	  if (fVerbosity>5)
	    cout << count << " " << tcount << ": " << title << endl;
          if (count != tcount)
            AliFatal(Form("Index mismatch: %d %d", count, tcount));
	  ++count;
	}
      }
    }
  }
}

//________________________________________________________________________
void AliDhcTask::InitEventMixer()
{
  // The effective pool size in events is set by trackDepth, so more
  // low-mult events are required to maintain the threshold than
  // high-mult events. Centrality pools are indep. of data histogram
  // binning, no need to match.

  // Centrality pools
  Int_t nCentBins=fMixBCent->GetNbins();
  Double_t centBins[nCentBins+1];
  centBins[0] = fMixBCent->GetBinLowEdge(1);
  for (Int_t i=1; i<=nCentBins; i++) {
    centBins[i] = fMixBCent->GetBinUpEdge(i);
  }
 
  // Z-vertex pools
  Int_t nZvtxBins=fMixBZvtx->GetNbins();
  Double_t zvtxbin[nZvtxBins+1];
  zvtxbin[0] = fMixBZvtx->GetBinLowEdge(1);
  for (Int_t i=1; i<=nZvtxBins; i++) {
    zvtxbin[i] = fMixBZvtx->GetBinUpEdge(i);
  }

  fPoolMgr = new AliEvtPoolManager();
  fPoolMgr->SetTargetTrackDepth(fTrackDepth);
  if (fVerbosity>4)
    fPoolMgr->SetDebug(1);
  fPoolMgr->InitEventPools(fPoolSize, nCentBins, centBins, nZvtxBins, zvtxbin);
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
  MiniEvent* sTracks = new MiniEvent(0); // deleted by pool mgr.

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

  Bool_t mcgen = 0;
  if (fTracksName.Contains("Gen"))
    mcgen = 1;

  // Centrality, vertex, other event variables...
  if (mcgen) {
    fZVertex = 0;
    TList *list = InputEvent()->GetList();
    TClonesArray *tcaTracks = dynamic_cast<TClonesArray*>(list->FindObject(fTracksName));
    if (tcaTracks) 
      fCentrality = tcaTracks->GetEntries();
  } else {
    if (dType == kESD) {
      if (!VertexOk(fESD)) {
        if (fVerbosity > 1)
          AliInfo(Form("Event REJECTED (ESD vertex not OK). z = %.1f", fZVertex));
        return;
      }
      const AliESDVertex* vertex = fESD->GetPrimaryVertex();
      fZVertex = vertex->GetZ();
      if(fESD->GetCentrality()) {
        fCentrality = 
          fESD->GetCentrality()->GetCentralityPercentile(fCentMethod);
      }
    }
    if (dType == kAOD) {
      const AliAODVertex* vertex = fAOD->GetPrimaryVertex();
      fZVertex = vertex->GetZ();
      if (!VertexOk(fAOD)) {
        if (fVerbosity > 1)
          Info("Exec()", "Event REJECTED (AOD vertex not OK). z = %.1f", fZVertex);
        return;
      }
      const AliCentrality *aodCent = fAOD->GetHeader()->GetCentralityP();
      if (aodCent) {
        fCentrality = aodCent->GetCentralityPercentile(fCentMethod);
      }
    }
  }

  // Fill Event histogram
  fHEvt->Fill(fZVertex, fCentrality);
  if (fCentrality > fBCent->GetXmax() || fCentrality < fBCent->GetXmin()) {
    AliInfo(Form("Event REJECTED (centrality out of range). fCentrality = %.1f", fCentrality));
    return;
  }

  // Get array of selected tracks
  if (dType == kESD) {
    GetESDTracks(sTracks);
  }
  if (dType == kAOD) {
    GetAODTracks(sTracks);
  }

  // Get pool containing tracks from other events like this one
  AliEvtPool* pool = fPoolMgr->GetEventPool(fCentrality, fZVertex);
  if (!pool) {
    AliWarning(Form("No pool found. Centrality %f, ZVertex %f", 
		    fCentrality, fZVertex));
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
	    Correlate(*evI, *evJ, kDiffEvt);
	  }
	}
      }
    }
  } else { /* standard case: same event, then mix*/
    Correlate(*sTracks, *sTracks, kSameEvt);  
    Int_t nMix = pool->GetCurrentNEvents();
    for (Int_t jMix=0; jMix<nMix; ++jMix) {
      MiniEvent* bgTracks = pool->GetEvent(jMix);
      Correlate(*sTracks, *bgTracks, kDiffEvt);
    }
  }

  pool->UpdatePool(sTracks);
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliDhcTask::GetESDTracks(MiniEvent* miniEvt)
{
  // Loop twice: 1. Count sel. tracks. 2. Fill vector.

  if (fTracksName.IsNull()) {
    const AliESDVertex *vtxSPD = fESD->GetPrimaryVertexSPD();
    if (!vtxSPD)
      return;
    
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
      Bool_t trkOK = fEsdTPCOnly->AcceptTrack(esdtrack);
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
    return;
  }

  TList *list = InputEvent()->GetList();
  TClonesArray *tcaTracks = dynamic_cast<TClonesArray*>(list->FindObject(fTracksName));
  const Int_t ntracks = tcaTracks->GetEntries();
  if (miniEvt)
    miniEvt->reserve(ntracks);
  else {
    AliError("Ptr to miniEvt zero");
    return;
  }

  for(Int_t itrack = 0; itrack < ntracks; itrack++) {
    AliVTrack *esdtrack = static_cast<AliESDtrack*>(tcaTracks->At(itrack));
    if(!esdtrack) {
      AliError(Form("ERROR: Could not retrieve esdtrack %d",itrack));
      continue;
    }
    Double_t pt = esdtrack->Pt();
    Double_t eta  = esdtrack->Eta();
    Double_t phi  = esdtrack->Phi();
    Int_t    sign = esdtrack->Charge() > 0 ? 1 : -1;
    miniEvt->push_back(AliMiniTrack(pt, eta, phi, sign));
  }
}

//________________________________________________________________________
void AliDhcTask::GetAODTracks(MiniEvent* miniEvt)
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

  if (miniEvt)
    miniEvt->reserve(nSelTrax);
  else {
    AliError("!miniEvt");
    return;
  }
  
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
    miniEvt->push_back(AliMiniTrack(pt, eta, phi, sign));
  }
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
Int_t AliDhcTask::Correlate(const MiniEvent &evt1, const MiniEvent &evt2, Int_t pairing)
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

  Int_t ptindex = (Int_t)fIndex->Eval(1,1,zbin,cbin);

  fHPtTrg_Evt->Reset();
  if (fDoWeights) { // Count trigger particles in this event
    for (Int_t i=0; i<iMax; ++i) {
      const AliMiniTrack &a(evt1.at(i));
      Float_t pta  = a.Pt();
      fHPtTrg_Evt->Fill(pta);
      fHPts[ptindex]->Fill(pta);
    }
  }

  for (Int_t i=0; i<iMax; ++i) {

    // Trigger particles
    const AliMiniTrack &a(evt1.at(i));

    Float_t pta  = a.Pt();
    Int_t abin = fHPtTrg->FindBin(pta);
    if (fHPtTrg->IsBinOverflow(abin) ||
	fHPtTrg->IsBinUnderflow(abin))
      continue;

    if (pairing == kSameEvt) {
      fHTrk->Fill(a.Phi(),a.Eta());
      fHPtTrg->AddBinContent(abin);
    }

    for (Int_t j=0; j<jMax; ++j) {
      // Associated particles
      if (pairing == kSameEvt && i==j)
	continue;

      const AliMiniTrack &b(evt2.at(j));
      
      Float_t ptb  = b.Pt();
      if (pta < ptb) 
	continue;

      Int_t bbin = fHPtAss->FindBin(ptb);
      if (fHPtAss->IsBinOverflow(bbin) ||
	  fHPtAss->IsBinUnderflow(bbin))
	continue;

      Float_t dphi = DeltaPhi(a.Phi(), b.Phi());
      Float_t deta = a.Eta() - b.Eta();

      Int_t index = globIndex+(abin-1)*nPtAssc+(bbin-1);
      Double_t weight = 1.0;
      if (fDoWeights) 
        weight = fHPtTrg_Evt->GetBinContent(abin);

      if (weight==0.0) {
        AliError(Form("Trying to work with weight = %f",weight));
      } else {
        hist[index]->Fill(dphi,deta,1./weight);
      }

      if (pairing == kSameEvt) {
	fHPtAss->AddBinContent(bbin);
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
