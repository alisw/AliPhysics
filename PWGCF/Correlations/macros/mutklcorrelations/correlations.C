#include "TChain.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "AliCFParticle.h"
#include "THn.h"
#include "TH3D.h"
#include "AliEventPoolManager.h"
#include "AliVEvent.h"
#include "Pool.h"
#include "fstream"
#include "TGrid.h"
#include "env.h"

TString sTrackType[5] = {"mc","hyb","its","tkl","mu"};
// Track histograms: deta, pt ass, pt trig, cent, dphi, zvtx 

// muon-tracklet reconstructed bins in multiplicity
//Int_t trg=4; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets; 4 -Muons; 5 - Gen muons
//Int_t ass=3; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets;
//const Float_t etaMinTrg = -4.0;
//const Float_t etaMaxTrg = -2.5;
//const Int_t   ntbins[kTrackNvar] = {  35,  5,    6,       6,      72, 14};
//const Double_t xtmin[kTrackNvar] = {-5.0,  0,    0,       0, -0.5*pi, -7};
//const Double_t xtmax[kTrackNvar] = {-1.5,  5,    4,    1000, +1.5*pi, +7};
//Double_t binsa[] = {0., 1., 2., 3., 4., 5.};
//Double_t binst[] = {0., 0.5, 1.0, 1.5, 2.0, 3.0, 4.0};

// tracklet-tracklet reconstructed
Int_t trg=3; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets; 4 -Muons; 5 - Gen muons
Int_t ass=3; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets;
const Float_t etaMinTrg = -1.0;
const Float_t etaMaxTrg = +1.0;
const Int_t   ntbins[kTrackNvar] = {  40,  5,    5,       6,      72, 14};
const Double_t xtmin[kTrackNvar] = {-2.0,  0,    0,       0, -0.5*pi, -7};
const Double_t xtmax[kTrackNvar] = {+2.0,  5,    5,    1000, +1.5*pi, +7};
Double_t binsa[] = {0., 1., 2., 3., 4., 5.};
Double_t binst[] = {0., 1., 2., 3., 4., 5.};


// central-forward mc
//Int_t trg=0; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets; 4 -Muons; 5 - Gen muons
//Int_t ass=0; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets;
//const Float_t etaMinTrg = -4.0;
//const Float_t etaMaxTrg = -2.5;
//const Int_t   ntbins[kTrackNvar] = {  35,  6,    6,       7,      72,  1};
//const Double_t xtmin[kTrackNvar] = {-5.0,  0,    0,       0, -0.5*pi, -7};
//const Double_t xtmax[kTrackNvar] = {-1.5,  4,    4,   101.1, +1.5*pi, +7};
//Double_t binsa[] = {0., 0.5, 1.0, 1.5, 2.0, 3.0, 4.0};
//Double_t binst[] = {0., 0.5, 1.0, 1.5, 2.0, 3.0, 4.0};

// central-central mc
//Int_t trg=0; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets; 4 -Muons; 5 - Gen muons
//Int_t ass=0; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets;
//const Float_t etaMinTrg = -1.0;
//const Float_t etaMaxTrg = +1.0;
//const Int_t   ntbins[kTrackNvar] = {  40,  6,    6,       7,      72,  1};
//const Double_t xtmin[kTrackNvar] = {-2.0,  0,    0,       0, -0.5*pi, -7};
//const Double_t xtmax[kTrackNvar] = {+2.0,  4,    4,   101.1, +1.5*pi, +7};
//Double_t binsa[] = {0., 0.5, 1.0, 1.5, 2.0, 3.0, 4.0};
//Double_t binst[] = {0., 0.5, 1.0, 1.5, 2.0, 3.0, 4.0};

// muon-tracklet reconstructed
//Int_t trg=4; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets; 4 -Muons; 5 - Gen muons
//Int_t ass=3; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets;
//const Float_t etaMinTrg = -4.0;
//const Float_t etaMaxTrg = -2.5;
//const Int_t   ntbins[kTrackNvar] = {  35,  5,    6,       7,      72, 14};
//const Double_t xtmin[kTrackNvar] = {-5.0,  0,    0,       0, -0.5*pi, -7};
//const Double_t xtmax[kTrackNvar] = {-1.5,  5,    4,   101.1, +1.5*pi, +7};
//Double_t binsa[] = {0., 1., 2., 3., 4., 5.};
//Double_t binst[] = {0., 0.5, 1.0, 1.5, 2.0, 3.0, 4.0};

// tracklet-tracklet reconstructed
//Int_t trg=3; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets; 4 -Muons; 5 - Gen muons
//Int_t ass=3; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets;
//const Float_t etaMinTrg = -1.0;
//const Float_t etaMaxTrg = +1.0;
//const Int_t   ntbins[kTrackNvar] = {  40,  5,    5,       7,      72, 14};
//const Double_t xtmin[kTrackNvar] = {-2.0,  0,    0,       0, -0.5*pi, -7};
//const Double_t xtmax[kTrackNvar] = {+2.0,  5,    5,   101.1, +1.5*pi, +7};
//Double_t binsa[] = {0., 1., 2., 3., 4., 5.};
//Double_t binst[] = {0., 1., 2., 3., 4., 5.};

// muon-tracklet reconstructed eta bins 3.5-4
//Int_t trg=4; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets; 4 -Muons; 5 - Gen muons
//Int_t ass=3; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets;
//const Float_t etaMinTrg = -4.0;
//const Float_t etaMaxTrg = -3.5;
//const Int_t   ntbins[kTrackNvar] = {  25,  5,    6,       7,      72, 14};
//const Double_t xtmin[kTrackNvar] = {-5.0,  0,    0,       0, -0.5*pi, -7};
//const Double_t xtmax[kTrackNvar] = {-2.5,  5,    4,   101.1, +1.5*pi, +7};
//Double_t binsa[] = {0., 1., 2., 3., 4., 5.};
//Double_t binst[] = {0., 0.5, 1.0, 1.5, 2.0, 3.0, 4.0};

// muon-tracklet reconstructed eta bins 3-3.5
//Int_t trg=4; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets; 4 -Muons; 5 - Gen muons
//Int_t ass=3; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets;
//const Float_t etaMinTrg = -3.5;
//const Float_t etaMaxTrg = -3.0;
//const Int_t   ntbins[kTrackNvar] = {  25,  5,    6,       7,      72, 14};
//const Double_t xtmin[kTrackNvar] = {-4.5,  0,    0,       0, -0.5*pi, -7};
//const Double_t xtmax[kTrackNvar] = {-2.0,  5,    4,   101.1, +1.5*pi, +7};
//Double_t binsa[] = {0., 1., 2., 3., 4., 5.};
//Double_t binst[] = {0., 0.5, 1.0, 1.5, 2.0, 3.0, 4.0};

// muon-tracklet reconstructed eta bins 2.5-3
//Int_t trg=4; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets; 4 -Muons; 5 - Gen muons
//Int_t ass=3; // 0 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets;
//const Float_t etaMinTrg = -3.0;
//const Float_t etaMaxTrg = -2.5;
//const Int_t   ntbins[kTrackNvar] = {  25,  5,    6,       7,      72, 14};
//const Double_t xtmin[kTrackNvar] = {-4.0,  0,    0,       0, -0.5*pi, -7};
//const Double_t xtmax[kTrackNvar] = {-1.5,  5,    4,   101.1, +1.5*pi, +7};
//Double_t binsa[] = {0., 1., 2., 3., 4., 5.};
//Double_t binst[] = {0., 0.5, 1.0, 1.5, 2.0, 3.0, 4.0};

// track-track reconstructed
//Int_t trg=1; // 1 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets; 4 -Muons; 5 - Gen muons
//Int_t ass=1; // 1 - MC tracks; 1 - Hybrid tracks; 2 - ITS-SA tracks; 3 - Tracklets;
//const Float_t etaMinTrg = -1.0;
//const Float_t etaMaxTrg = +1.0;
//const Int_t   ntbins[kTrackNvar] = {  40,   4,    4,       7,      72, 14};
//const Double_t xtmin[kTrackNvar] = {-2.0, 0.5,  0.5,       0, -0.5*pi, -7};
//const Double_t xtmax[kTrackNvar] = {+2.0,   4,    4,   101.1, +1.5*pi, +7};
//Double_t binsa[] = {0.5, 1., 2., 3., 4.};
//Double_t binst[] = {0.5, 1., 2., 3., 4.};

const Float_t etaMinAss = -1.0;
const Float_t etaMaxAss = +1.0;
//Double_t binsc[] = {0,5,20,40,60,70,80,101.1}; 
Double_t binsc[] = {0,5,50,100,150,200,1000}; 
Double_t binsz[100];
// Event histograms: pt trig, cent, zvtx
Int_t   nebins[kEventNvar] = {ntbins[kTrackPtTr],ntbins[kTrackCent],ntbins[kTrackZvtx]};
Double_t xemin[kEventNvar] = { xtmin[kTrackPtTr], xtmin[kTrackCent], xtmin[kTrackZvtx]};
Double_t xemax[kEventNvar] = { xtmax[kTrackPtTr], xtmax[kTrackCent], xtmax[kTrackZvtx]};

Bool_t useMCforTracklets=0;

TAxis* axis[6];
TH2D* gFilter;

Int_t SelectTracks(TClonesArray* tracksAll, TClonesArray* tracksSelected, Int_t trackType, Float_t etaMin, Float_t etaMax, Float_t ptMin, Float_t ptMax, Bool_t isHighPtMuonTriggerFired=0, Bool_t apply_pdca=0, Int_t pdg=0);
void FillCorrelations(TClonesArray* vTrg, TClonesArray* vAss, Float_t cent, Float_t zvtx, THnD* hTrackHist, THnD* hEventHist, Bool_t ptOrdering=0, Bool_t mixing=0);
Float_t dphi(Float_t phi1,Float_t phi2);
Float_t phiCorrected(Float_t phi, Float_t dphi);


void correlations(
    TString centMethod = "V0L",
    UInt_t triggerMask = AliVEvent::kMUS7 | AliVEvent::kMUSH7,
    Int_t nMaxEventsInPull=10,
    Bool_t ptOrdering=0,
    Bool_t selectEventsWithMuonsOnly=1,
    Bool_t apply_pdca=1,
    Bool_t efficiency_correction=0,
    Int_t pdg=0,
    Int_t thread=-1,
    Int_t nEventsInThread=1000000
)
{
  if (efficiency_correction) {
    TFile* f = new TFile("Filter.root");
    if (f) gFilter = (TH2D*) f->Get("Filter");
  }

  TChain* fTree = new TChain("CorrelationTree/events");
  for (Int_t i=0;i<12;i++) fTree->AddFile(Form("/data/trees/mu/lhc13b/AnalysisResults.%03i.root",i));
  for (Int_t i=0;i<13;i++) fTree->AddFile(Form("/data/trees/mu/lhc13c/AnalysisResults.%03i.root",i));
  for (Int_t i=0;i<20;i++) fTree->AddFile(Form("/data/trees/mu/lhc13d/AnalysisResults.%03i.root",i));
  for (Int_t i=0;i<26;i++) fTree->AddFile(Form("/data/trees/mu/lhc13e/AnalysisResults.%03i.root",i));
//  for (Int_t i=0;i<64;i++) fTree->AddFile(Form("/data/trees/mu/lhc13f/AnalysisResults.%03i.root",i));

//  if (!TGrid::Connect("alien://")) return;
//  TGridCollection* coll = gGrid->OpenCollection("wn.xml");
//  while (coll->Next()) { fTree->Add(coll->GetTURL());}

  
  Float_t fCentrality[9];
  Float_t fVtxZ;
  Bool_t fVtxTPConly;
  UInt_t fVtxContributors;
  Int_t fRunNumber;
  UInt_t fClassesFired;
  UInt_t fMask;
  Bool_t fIsPileupSPD;
  Int_t fNofITSClusters[6]={0,0,0,0,0,0};
  Float_t fMultV0Cring[4];    //  tree var: 
 
  TClonesArray* fTracks      = new TClonesArray("AliCFParticle",2000);
  TClonesArray* fTracklets   = new TClonesArray("AliCFParticle",2000);
  TClonesArray* fMuons       = new TClonesArray("AliCFParticle",2000);
  TClonesArray* fMcParticles = new TClonesArray("AliCFParticle",2000);

  TClonesArray* fMcMuons = new TClonesArray("AliCFParticle",2000);
  Float_t fMultGenV0S=0;
  TFile* fCent = new TFile("centrality.root");
  TH1D* hCent = (TH1D*) fCent->Get("hCentGenV0S");
  fTree->SetBranchAddress("multGenV0S",&fMultGenV0S);
  fTree->SetBranchAddress("mcmuons",&fMcMuons);

  fTree->SetBranchAddress("cent",&fCentrality);
  fTree->SetBranchAddress("vtxz",&fVtxZ);
  fTree->SetBranchAddress("run",&fRunNumber);
  fTree->SetBranchAddress("classes",&fClassesFired);
  fTree->SetBranchAddress("mask",&fMask);
  fTree->SetBranchAddress("pileupspd",&fIsPileupSPD);
  fTree->SetBranchAddress("tracks",&fTracks);
  fTree->SetBranchAddress("vtxTPConly",&fVtxTPConly);
  fTree->SetBranchAddress("vtxContributors",&fVtxContributors);
  fTree->SetBranchAddress("nofITSClusters",&fNofITSClusters);
  fTree->SetBranchAddress("tracklets",&fTracklets);
  fTree->SetBranchAddress("muons",&fMuons);
  fTree->SetBranchAddress("mcparticles",&fMcParticles);
  fTree->SetBranchAddress("multV0Cring",&fMultV0Cring);

  for (Int_t i=0;i<=ntbins[kTrackZvtx];i++) binsz[i]=xtmin[kTrackZvtx]+(xtmax[kTrackZvtx]-xtmin[kTrackZvtx])*i/ntbins[kTrackZvtx];
  THnD* hTrackHistS = new THnD("hTrackHistS","",kTrackNvar,ntbins,xtmin,xtmax);
  hTrackHistS->GetAxis(kTrackPtAs)->Set(ntbins[kTrackPtAs],binsa);
  hTrackHistS->GetAxis(kTrackPtTr)->Set(ntbins[kTrackPtTr],binst);
  hTrackHistS->GetAxis(kTrackCent)->Set(ntbins[kTrackCent],binsc);
  THnD* hTrackHistM = (THnD*) hTrackHistS->Clone("hTrackHistM");
  
  THnD* hEventHistS = new THnD("hEventHistS","",kEventNvar,nebins,xemin,xemax);
  hEventHistS->GetAxis(kEventPtTr)->Set(nebins[kEventPtTr],binst);
  hEventHistS->GetAxis(kEventCent)->Set(nebins[kEventCent],binsc);
  THnD* hEventHistM = (THnD*) hEventHistS->Clone("hEventHistM");
  
  hTrackHistS->Sumw2();
  hTrackHistM->Sumw2();
  hEventHistS->Sumw2();
  hEventHistM->Sumw2();
  
  for (Int_t i=0;i<kTrackNvar;i++) axis[i] = hTrackHistS->GetAxis(i);
  
  TH3D* hEventCount = new TH3D("hEventCount",";centrality;;",1011,0.0,101.1,1,0,1,1,0,1);
  
  PoolManager* fPoolMgr = new PoolManager(1,nMaxEventsInPull,ntbins[kTrackCent],binsc,ntbins[kTrackZvtx],binsz);
  
  Int_t nEvents = fTree->GetEntries();
  Int_t evStart = 0;
  if (thread>=0) {
    nEvents = TMath::Min(nEvents,(thread+1)*nEventsInThread);
    evStart = thread*nEventsInThread;
  }
  printf("Events=%i\n",nEvents);
  
  Int_t icent = 0;
  if      (centMethod.EqualTo("V0M"))   icent=0;
  else if (centMethod.EqualTo("V0A"))   icent=1;
  else if (centMethod.EqualTo("V0C"))   icent=2;
  else if (centMethod.EqualTo("CL1"))   icent=3;
  else if (centMethod.EqualTo("ZNA"))   icent=4;
  else if (centMethod.EqualTo("ZNC"))   icent=5;
  else if (centMethod.EqualTo("V0A23")) icent=6;
  else if (centMethod.EqualTo("V0C01")) icent=7;
  else if (centMethod.EqualTo("V0S"))   icent=8;
  else if (centMethod.EqualTo("V0L"))   icent=9;
  else { printf("centrality method %s not supported\n",centMethod.Data()); return; }
  TClonesArray* vTrg = new TClonesArray("AliCFParticle",100);
  TClonesArray* vAss = new TClonesArray("AliCFParticle",100);
  
  for (Int_t ev=evStart;ev<nEvents;ev++){
    if (ev%10000==0) printf("Event=%i\n",ev);
    fTree->GetEntry(ev);
    Float_t cent;
    if (icent<9) {
      cent = fCentrality[icent];
      cent = hCent->GetBinContent(hCent->FindFixBin(fMultGenV0S));
      if (cent>100) cent=100;
    } else {
      cent = fMultV0Cring[0]+fMultV0Cring[1];
      if (cent>=1000) cent=999.9;
    }
    
    hEventCount->Fill(cent,"all",Form("%i",fRunNumber),1.);
    if (trg!=0 || ass!=0){
      if (!(fMask & triggerMask)) continue;
      hEventCount->Fill(cent,"after PS check",Form("%i",fRunNumber),1.);
      if (fIsPileupSPD) continue;
      hEventCount->Fill(cent,"after pileup check",Form("%i",fRunNumber),1.);
      if (fVtxTPConly) continue;
      hEventCount->Fill(cent,"after TPConly vertex cut",Form("%i",fRunNumber),1.);
      if (fVtxContributors<1) continue;
      hEventCount->Fill(cent,"after cut on vertex contributors",Form("%i",fRunNumber),1.);
      if (TMath::Abs(fVtxZ)>xtmax[5]) continue;
      hEventCount->Fill(cent,"after zvtx cut",Form("%i",fRunNumber),1.);

      if (selectEventsWithMuonsOnly) {
        Bool_t isTriggerMuon=0;
        for (Int_t iMuons=0;iMuons<fMuons->GetEntriesFast();iMuons++) 
          if (((AliCFParticle*) fMuons->UncheckedAt(iMuons))->Mask()>=1) { isTriggerMuon=1; break; }
        if (!isTriggerMuon) continue;
        hEventCount->Fill(cent,"after cut on muons",Form("%i",fRunNumber),1.);
      }
    }

    Bool_t isHighPtMuonTriggerFired = fMask & AliVEvent::kMUSH7; 
    
    SelectTracks(trg==5 ? fMcMuons : trg==0 ? fMcParticles : (trg==3 ? fTracklets : (trg==4 ? fMuons : fTracks)),vTrg,trg,etaMinTrg,etaMaxTrg,xtmin[kTrackPtTr],xtmax[kTrackPtTr],isHighPtMuonTriggerFired,apply_pdca,pdg);
    SelectTracks(ass==0 ? fMcParticles : (ass==3 ? fTracklets : fTracks)                    ,vAss,ass,etaMinAss,etaMaxAss,xtmin[kTrackPtAs],xtmax[kTrackPtAs],isHighPtMuonTriggerFired,apply_pdca);

    FillCorrelations(vTrg,vAss,cent,fVtxZ,hTrackHistS,hEventHistS,ptOrdering);
    
    Pool* pool = fPoolMgr->GetEventPool(cent, fVtxZ);
    if (!pool) continue;
    if (pool->IsReady()) {
      for (Int_t ev2=0;ev2<pool->GetCurrentNEvents();ev2++){
        FillCorrelations(pool->GetEvent(ev2),vAss,cent,fVtxZ,hTrackHistM,hEventHistM,ptOrdering,1); 
      } // mixed events
    } // pool ready
    //fPoolMgr->PrintInfo();
    pool->UpdatePool(vTrg);
  }
  hEventCount->LabelsDeflate("Z");

  TString outFileName = Form("corr_%s_%s_%s",sTrackType[trg].Data(),sTrackType[ass].Data(),centMethod.Data());
  TFile* fout = new TFile(thread <0 ? Form("%s_bcde.root",outFileName.Data()) : Form("%s_%i.root",outFileName.Data(),thread),"recreate");
  hTrackHistS->Write();
  hTrackHistM->Write();
  hEventHistS->Write();
  hEventHistM->Write();
  hEventCount->Write();
  fout->Close();
}


Float_t dphi(Float_t phi1,Float_t phi2){
  Float_t dphi = phi1-phi2;
  if (dphi > +1.5*pi) dphi -= 2*pi;
  if (dphi < -0.5*pi) dphi += 2*pi;
  return dphi;
}

Float_t phiCorrected(Float_t phi, Float_t dphi){
  phi += 39./34.*dphi;
  if (phi<0  )  phi+=2*pi;
  if (phi>2*pi) phi-=2*pi;
  return phi;
}

void FillCorrelations(TClonesArray* vTrg, TClonesArray* vAss, Float_t cent, Float_t zvtx, THnD* hTrackHist, THnD* hEventHist, Bool_t ptOrdering, Bool_t mixing){
  Int_t ix[kTrackNvar];
  Int_t iv[kEventNvar];
  ix[kTrackCent] = axis[kTrackCent]->FindFixBin(cent);
  ix[kTrackZvtx] = axis[kTrackZvtx]->FindFixBin(zvtx);
  iv[kEventCent] = ix[kTrackCent];
  iv[kEventZvtx] = ix[kTrackZvtx];

  Float_t pt1,eta1,phi1,pt2,eta2,phi2;
  UInt_t id1,id2;
  AliCFParticle* d1;
  AliCFParticle* d2;
  
  for (Int_t i=0;i<vTrg->GetEntriesFast();i++){
    d1 = (AliCFParticle*) vTrg->UncheckedAt(i);
    id1  = d1->GetUniqueID();
    pt1  = d1->Pt();
    eta1 = d1->Eta();
    ix[kTrackPtTr] = axis[kTrackPtTr]->FindFixBin(pt1);
    phi1 = d1->Phi();
    iv[kEventPtTr] = ix[kTrackPtTr];
    Double_t eff = gFilter? gFilter->GetBinContent(gFilter->FindFixBin(pt1,eta1)) : 1.;
    if (eff<1.e-10) continue;
    hEventHist->AddBinContent(iv,1./eff);
    hEventHist->AddBinError2(hEventHist->GetBin(iv),1./eff/eff);

    for (Int_t j=0;j<vAss->GetEntriesFast();j++){
      d2 = (AliCFParticle*) vAss->UncheckedAt(j);
      id2 = d2->GetUniqueID();
      if (!mixing && id1==id2) continue; // check unique id
      pt2 = d2->Pt();
      if (ptOrdering) if(pt2>pt1) continue;
      eta2 = d2->Eta();
      phi2 = d2->Phi();
      ix[kTrackPtAs] = axis[kTrackPtAs]->FindFixBin(pt2);
      ix[kTrackDeta] = axis[kTrackDeta]->FindFixBin(eta1-eta2);
      ix[kTrackDphi] = axis[kTrackDphi]->FindFixBin(dphi(phi1,phi2));
      hTrackHist->AddBinContent(ix,1./eff); // much faster than ->Fill(...), no need to search for bin ids for cent, zvtx and pt trigger
      hTrackHist->AddBinError2(hTrackHist->GetBin(ix),1./eff/eff);
    }
  }
}


Int_t SelectTracks(TClonesArray* tracksAll, TClonesArray* tracksSelected, Int_t trackType, Float_t etaMin, Float_t etaMax, Float_t ptMin, Float_t ptMax, Bool_t isHighPtMuonTriggerFired, Bool_t apply_pdca, Int_t pdg) {
  tracksSelected->Clear();
  Float_t pt,eta,phi,dphi;
  Short_t charge;
  Int_t mask;
  Int_t nSelectedTracks = 0;
  for (Int_t i=0;i<tracksAll->GetEntriesFast();i++){
    AliCFParticle* track = (AliCFParticle*) tracksAll->UncheckedAt(i);
    pt     = track->Pt();
    eta    = track->Eta();
    phi    = track->Phi();
    mask   = track->Mask();
    charge = track->Charge();
    if       (trackType==0) { if (pdg && TMath::Abs(mask)!=pdg) continue; } // select only pi,K,HFM
    else if  (trackType==1) { if (!(mask&(1<<8|1<<9))) continue; }
    else if  (trackType==2) { if (!(mask&(1<<1))     ) continue; }
    else if  (trackType==3) { 
      dphi = pt;
      phi  = phiCorrected(phi,dphi);
      pt   = TMath::Abs(dphi)*1000; // mrad, pseudo pt for tracklets
      if (useMCforTracklets) { 
        if (!mask || !charge) continue; // skip secondaries and fakes
        eta = track->GetAt(1); // etaMC
        phi = track->GetAt(2); // phiMC
      }
    }
    else if  (trackType==4) {
//      if (isHighPtMuonTriggerFired) printf("%i %i %6.2f\n",isHighPtMuonTriggerFired,mask,pt);
      if (isHighPtMuonTriggerFired) { if (mask<2) continue; }
      else                          { if (mask<2) continue; }
      Bool_t pdca = TMath::Nint(track->GetAt(3));
      if (apply_pdca && !pdca) continue;  // skip muons without pdca
    }
    else if  (trackType==5) {
      Int_t absmask = TMath::Abs(mask);
      if (pdg==211) if (absmask!=211) continue; 
      if (pdg==321) if (absmask!=321 && absmask!=310 && absmask!=130) continue; 
    } // select only mu from pi,K
    
    if (eta<etaMin || eta>etaMax) continue;
    if (pt<ptMin || pt>ptMax) continue;
    
    AliCFParticle* part = new ((*tracksSelected)[nSelectedTracks++]) AliCFParticle(pt,eta,phi,charge,mask);
    part->SetUniqueID(10*i+trackType);
  }
  
  // if muons remove unlike-sign events
  return nSelectedTracks;
  //if (trackType!=4) return nSelectedTracks;
  //if (tracksSelected->GetEntriesFast()<2) return nSelectedTracks;
  //charge = ((AliCFParticle*) tracksSelected->UncheckedAt(0))->Charge();
  //for (Int_t i=1;i<tracksSelected->GetEntriesFast();i++){
  //  AliCFParticle* track = (AliCFParticle*) tracksSelected->UncheckedAt(i);
  //  if (charge!=track->Charge()) {charge = 0; break; }
  //}
  //if (charge) return nSelectedTracks;
  //
  //tracksSelected->Clear();
  //return 0;
}
