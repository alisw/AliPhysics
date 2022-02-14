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
#include "utils.h"
#include "env.h"

TAxis* axis[6];
TH2D* gFilter;
TH3D* gFilter3D;

Int_t SelectTracks(TClonesArray* tracksAll, TClonesArray* tracksSelected, Int_t trackType, Float_t etaMin, Float_t etaMax, Float_t ptMin, Float_t ptMax, kFB fb);
void FillCorrelations(TClonesArray* vTrg, TClonesArray* vAssoc, Float_t cent, Float_t zvtx, THnD* hTrackHist, THnD* hEventHist, Bool_t ptOrdering=0, Bool_t mixing=0, Float_t dphimin=-0.5*pi, Float_t dphimax=1.5*pi);
Float_t dphi(Float_t phi1,Float_t phi2,Float_t phimin=-0.5*pi,Float_t phimax=1.5*pi);
Float_t phiCorrected(Float_t phi, Float_t dphi);

void correlations(
    char * period="lhc10de",
    TString centMethod = "V0M",
    Int_t anaType = kTrkTrk,
    kFB fb = fb56,
    UInt_t triggerMask = AliVEvent::kMB,
    UInt_t classfired = 0xffffffff,
    Int_t energy=k7TeV,
    Bool_t efficiency_correction=1,
    char * mcperiod = "lhc10f6a", // needed for efficiency correction
    Bool_t localinputfiles=kTRUE,
    char * pathinputlocalfiles="~/ppcorr/data/", // lhc10de
    Int_t nMaxEventsInPool=100,
    Int_t nEventsInThread=1000000,
    Bool_t ptOrdering=1,
    Int_t thread=-1
)
{

  Float_t etaMinTrg = 0, etaMaxTrg = 0, etaMinAssoc = 0, etaMaxAssoc = 0;
  Int_t assoc=-999,trg=-999 ;
  switch (anaType){
  case kTrkTrk: 
    etaMinTrg = etaMin_trk; etaMaxTrg = etaMax_trk;
    etaMinAssoc = etaMin_trk; etaMaxAssoc = etaMax_trk;
    assoc=kTrk;trg=kTrk;
    break;
  case kTrkTrkGen: 
    etaMinTrg = etaMin_trk; etaMaxTrg = etaMax_trk;
    etaMinAssoc = etaMin_trk; etaMaxAssoc = etaMax_trk;
    assoc=kTrkGen;trg=kTrkGen;
    break;
  case kTrkTrkITS:
    etaMinTrg = etaMin_trkits; etaMaxTrg = etaMax_trkits;
    etaMinAssoc = etaMin_trkits; etaMaxAssoc = etaMax_trkits;
    assoc=kTrkITS;trg=kTrkITS;
    break;
  case kTrkTrkITSGen:
    etaMinTrg = etaMin_trkits; etaMaxTrg = etaMax_trkits;
    etaMinAssoc = etaMin_trkits; etaMaxAssoc = etaMax_trkits;
    assoc=kTrkITSGen;trg=kTrkITSGen;
    break;
  case kTklTkl:
    etaMinTrg = etaMin_tkl; etaMaxTrg = etaMax_tkl;
    etaMinAssoc = etaMin_tkl; etaMaxAssoc = etaMax_tkl;
    assoc=kTkl;trg=kTkl;
    break;
  case kTklTklMC:
    etaMinTrg = etaMin_tkl; etaMaxTrg = etaMax_tkl;
    etaMinAssoc = etaMin_tkl; etaMaxAssoc = etaMax_tkl;
    assoc=kTklMC;trg=kTklMC;
    break;
  case kTklTklGen:
    etaMinTrg = etaMin_tkl; etaMaxTrg = etaMax_tkl;
    etaMinAssoc = etaMin_tkl; etaMaxAssoc = etaMax_tkl;
    assoc=kTklGen;trg=kTklGen;
    break;
  default: cout << "Invalid analysis type choice!" << endl;
  }
  Printf("etaTrig = %.2lf - %.2lf",etaMinTrg,etaMaxTrg);
  Printf("etaAssoc = %.2lf - %.2lf",etaMinAssoc,etaMaxAssoc);

  //todo make it const?
  Int_t   ntbins[kTrackNvar] = {0};
  Double_t xtmin[kTrackNvar] = {0};
  Double_t xtmax[kTrackNvar] = {0};

  for(Int_t ivar=0;ivar<kTrackNvar;ivar++){
    if(anaType==kTklTkl || anaType==kTklTklMC){//tkl-tkl
      ntbins[ivar] = ntbins_tkltkl[ivar];
      xtmin[ivar]  = xtmin_tkltkl[ivar];
      xtmax[ivar]  = xtmax_tkltkl[ivar];
      if(anaType==kTklTklMC){//tkl-tkl - MC info on eta,phi
        useMCforTracklets=kTRUE;
      }
    }
    else if(anaType==kTklTklGen){//tkl-tkl Gen
      ntbins[ivar] = ntbins_tkltklGen[ivar];
      xtmin[ivar]  = xtmin_tkltklGen[ivar];
      xtmax[ivar]  = xtmax_tkltklGen[ivar];
    }
    else if(anaType==kTrkTrk || anaType==kTrkTrkGen){//track-track
      ntbins[ivar] = ntbins_trktrk[ivar];
      xtmin[ivar]  = xtmin_trktrk[ivar];
      xtmax[ivar]  = xtmax_trktrk[ivar];
    }
    else if(anaType==kTrkTrkITS || anaType==kTrkTrkITSGen){//track-track
      ntbins[ivar] = ntbins_trktrkits[ivar];
      xtmin[ivar]  = xtmin_trktrkits[ivar];
      xtmax[ivar]  = xtmax_trktrkits[ivar];
    }
    else {
      cout << "Invalid analysis type choice!" << endl;
      return;
    }
  }

  cout << "number of trigger and associated bins: " << ntbins[kTrackPtAs] << "   " << ntbins[kTrackPtTr] << endl;

  Double_t binsa[ntbins[kTrackPtAs]+1];
  Double_t binst[ntbins[kTrackPtTr]+1];
  if(anaType==kTklTkl || anaType==kTklTklMC){
    for(Int_t ipt=0;ipt<=ntbins[kTrackPtAs];ipt++) binsa[ipt]= bins_tkl[ipt];
    for(Int_t ipt=0;ipt<=ntbins[kTrackPtTr];ipt++) binst[ipt]= bins_tkl[ipt];
    //    Printf("%d",ntbins[kTrackPtTr]);
  }
  else if(anaType==kTklTklGen){
    for(Int_t ipt=0;ipt<=ntbins[kTrackPtAs];ipt++) binsa[ipt]= bins_tklGen[ipt];
    for(Int_t ipt=0;ipt<=ntbins[kTrackPtTr];ipt++) binst[ipt]= bins_tklGen[ipt];
  }
  else if(anaType==kTrkTrk || anaType==kTrkTrkGen){
    for(Int_t ipt=0;ipt<=ntbins[kTrackPtAs];ipt++) binsa[ipt]= bins_trk[ipt];
    for(Int_t ipt=0;ipt<=ntbins[kTrackPtTr];ipt++) binst[ipt]= bins_trk[ipt];
  }
  else if(anaType==kTrkTrkITS || anaType==kTrkTrkITSGen){
    for(Int_t ipt=0;ipt<=ntbins[kTrackPtAs];ipt++) binsa[ipt]= bins_trkits[ipt];
    for(Int_t ipt=0;ipt<=ntbins[kTrackPtTr];ipt++) binst[ipt]= bins_trkits[ipt];
  }
  else {
    cout << "Invalid analysis type choice!" << endl;
    return;
  }
  cout << "trigger pt bins: ";
  for(Int_t ipt=0;ipt<=ntbins[kTrackPtTr];ipt++) cout << binst[ipt] << "   ";
  cout << endl;
  cout << "associated pt bins: ";
  for(Int_t ipt=0;ipt<=ntbins[kTrackPtAs];ipt++) cout << binsa[ipt] << "   ";
  cout << endl;


  Double_t binsz[100];
  // Event histograms:                     pt_trig,              cent,              zvtx
  Int_t   nebins[kEventNvar] = {ntbins[kTrackPtTr],ntbins[kTrackCent],ntbins[kTrackZvtx]};
  Double_t xemin[kEventNvar] = { xtmin[kTrackPtTr], xtmin[kTrackCent], xtmin[kTrackZvtx]};
  Double_t xemax[kEventNvar] = { xtmax[kTrackPtTr], xtmax[kTrackCent], xtmax[kTrackZvtx]};

  if (efficiency_correction) {
    TFile* fEff = 0;
    if(anaType==kTrkTrkITS) fEff = new TFile(Form("../analysis/rootfiles/eff_%s_fbITSsa_final.root",mcperiod),"read");
    else fEff = new TFile(Form("../analysis/rootfiles/eff_%s_%s_final.root",mcperiod,fbname[fb].Data()),"read"); // output of efficiencies.C
    gFilter3D = (TH3D*)fEff->Get("hEffPtEtaVzGen");
  }

  TChain* fTree = new TChain("CorrelationTree/events");
  ChainFiles(fTree,localinputfiles,period,pathinputlocalfiles);
  SetTree(fTree);

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

  TH3D* hEventCount = new TH3D("hEventCount",";centrality;;",100,0.0,10.,1,0,1,1,0,1);
  TH2F* mixedDist = new TH2F("mixedDist", ";centrality;tracks;events", ntbins[kTrackCent], binsc, 200, 0, nMaxEventsInPool*1.5);
  PoolManager* fPoolMgr = new PoolManager(1,nMaxEventsInPool,ntbins[kTrackCent],binsc,ntbins[kTrackZvtx],binsz);

  Int_t nEvents = fTree->GetEntries();
  Int_t evStart = 0;
  if (thread>=0) {
    nEvents = TMath::Min(nEvents,(thread+1)*nEventsInThread);
    evStart = thread*nEventsInThread;
  }
  printf("Events=%i\n",nEvents);

  TClonesArray* vTrg = new TClonesArray("AliCFParticle",100);
  TClonesArray* vAssoc = new TClonesArray("AliCFParticle",100);

  Int_t nRun = 0;

  Bool_t isMC=0;
  if(anaType==kTrkTrkGen || anaType==kTklTklGen || anaType==kTrkTrkITSGen)isMC=1;
  if(fMcParticles)isMC=1; //this is needed to run the same analysis as in the data skipping the pileup selection

  for (Int_t ev=evStart;ev<nEvents;ev++){
    if (ev%100000==0) Printf("Event=%i   %.2lf%% done",ev,Double_t(ev)/Double_t(nEvents)*100.);
    fTree->GetEntry(ev);

    Float_t cent;

    if        (centMethod.EqualTo("V0M"))cent = fMultPercV0Meq;
    else if   (centMethod.EqualTo("CL1"))cent = 1;
    else if   (centMethod.EqualTo("V0A"))cent = fMultPercV0Aeq;
    else if   (centMethod.EqualTo("V0C"))cent = fMultPercV0Ceq;
    else if   (centMethod.EqualTo("TKL"))cent = fMultTKL;
    else if   (centMethod.EqualTo("TPC"))cent = fNchTPC;
    else if   (centMethod.EqualTo("V0Mmc"))cent = 1;
    else if   (centMethod.EqualTo("TPCmc"))cent = 1;
    else if   (centMethod.EqualTo("V0Amc"))cent = 1;
    else if   (centMethod.EqualTo("V0Cmc"))cent = 1;
    else { printf("centrality method %s not supported\n",centMethod.Data()); return; }
    if(cent < xtmin[kTrackCent] || cent > xtmax[kTrackCent]) {/*Printf("rejecting event with cent=%.2lf",cent);*/ continue;}

    //    if(cent > xtmax[kTrackCent]) {Printf("moving event with cent=%.2lf to cent=%.2lf",cent,xtmax[kTrackCent]-0.0001); cent = xtmax[kTrackCent]-0.0001;}
    //    hEventCount->Fill(cent,"all",Form("%i",fCurrentRunNumber),1.);
    //if (trg!=0 || assoc!=0){ // TODO! make this work for MC too!

    if(energy==k7TeV &&  !(IsGoodEvent7TeV(triggerMask ,classfired,xtmin[kTrackZvtx],xtmax[kTrackZvtx],cent,hEventCount)) )continue;
    if(energy==k13TeV && !(IsGoodEvent13TeV(triggerMask,classfired,xtmin[kTrackZvtx],xtmax[kTrackZvtx],cent,hEventCount,isMC)) )continue;

    //}

    //atLeastOnePair = kFALSE;

    switch (anaType){
    case kTrkTrk:
      SelectTracks(fTracks,vTrg,kTrk,etaMinTrg,etaMaxTrg,xtmin[kTrackPtTr],xtmax[kTrackPtTr],fb);
      SelectTracks(fTracks,vAssoc,kTrk,etaMinAssoc,etaMaxAssoc,xtmin[kTrackPtAs],xtmax[kTrackPtAs],fb);
      break;
    case kTrkTrkGen:
      SelectTracks(fMcParticles,vTrg,kTrkGen,etaMinTrg,etaMaxTrg,xtmin[kTrackPtTr],xtmax[kTrackPtTr],fb);
      SelectTracks(fMcParticles,vAssoc,kTrkGen,etaMinAssoc,etaMaxAssoc,xtmin[kTrackPtAs],xtmax[kTrackPtAs],fb);
      break;
    case kTrkTrkITSGen: //no need to create kTrkITSGen tracktype
      SelectTracks(fMcParticles,vTrg,kTrkITSGen,etaMinTrg,etaMaxTrg,xtmin[kTrackPtTr],xtmax[kTrackPtTr],fb);
      SelectTracks(fMcParticles,vAssoc,kTrkITSGen,etaMinAssoc,etaMaxAssoc,xtmin[kTrackPtAs],xtmax[kTrackPtAs],fb);
      break;
    case kTrkTrkITS:
      SelectTracks(fTracks,vTrg,kTrkITS,etaMinTrg,etaMaxTrg,xtmin[kTrackPtTr],xtmax[kTrackPtTr],fb);
      SelectTracks(fTracks,vAssoc,kTrkITS,etaMinAssoc,etaMaxAssoc,xtmin[kTrackPtAs],xtmax[kTrackPtAs],fb);
      break;
    case kTklTkl:
      SelectTracks(fTracklets,vTrg,kTkl,etaMinTrg,etaMaxTrg,xtmin[kTrackPtTr],xtmax[kTrackPtTr],fb);
      SelectTracks(fTracklets,vAssoc,kTkl,etaMinAssoc,etaMaxAssoc,xtmin[kTrackPtAs],xtmax[kTrackPtAs],fb);
      break;
    case kTklTklMC:
      SelectTracks(fTracklets,vTrg,kTkl,etaMinTrg,etaMaxTrg,xtmin[kTrackPtTr],xtmax[kTrackPtTr],fb);
      SelectTracks(fTracklets,vAssoc,kTkl,etaMinAssoc,etaMaxAssoc,xtmin[kTrackPtAs],xtmax[kTrackPtAs],fb);
      break;
    case kTklTklGen:
      SelectTracks(fMcParticles,vTrg,kTklGen,etaMinTrg,etaMaxTrg,xtmin[kTrackPtTr],xtmax[kTrackPtTr],fb);
      SelectTracks(fMcParticles,vAssoc,kTklGen,etaMinAssoc,etaMaxAssoc,xtmin[kTrackPtAs],xtmax[kTrackPtAs],fb);
      break;
    }
    FillCorrelations(vTrg,vAssoc,cent,fVtxZ,hTrackHistS,hEventHistS,ptOrdering,0,-1.*pi,pi);

    Pool* pool = fPoolMgr->GetEventPool(cent, fVtxZ);
    if (!pool) continue;

    if (pool->IsReady()) {
      mixedDist->Fill(cent,pool->GetCurrentNEvents());
      for (Int_t ev2=0;ev2<pool->GetCurrentNEvents();ev2++){
        //if(!atLeastOnePair) continue;
        FillCorrelations(pool->GetEvent(ev2),vAssoc,cent,fVtxZ,hTrackHistM,hEventHistM,ptOrdering,1,-1.*pi,pi);
      } // mixed events
    } // pool ready
    //    fPoolMgr->PrintInfo();
    //not updating the pool if there is only one trigger particle
    //if(!atLeastOnePair) continue;
    pool->UpdatePool(vTrg);
  }

  hEventCount->LabelsDeflate("Z");
  TString outFileName = Form("corr_%s_%s_%s_%s_%s",period,sTrackType[trg].Data(),sTrackType[assoc].Data(),fbname[fb].Data(),centMethod.Data());
  TFile* fout = new TFile(thread <0 ? Form("%s.root",outFileName.Data()) : Form("%s_%i.root",outFileName.Data(),thread),"recreate");
  hTrackHistS->Write();
  hTrackHistM->Write();
  hEventHistS->Write();
  hEventHistM->Write();
  hEventCount->Write();
  mixedDist->Write();
  fout->Close();
}


Float_t dphi(Float_t phi1,Float_t phi2,Float_t phimin,Float_t phimax){
  Float_t dphi = phi1-phi2;
  if (dphi > phimax) dphi -= 2*pi;
  if (dphi < phimin) dphi += 2*pi;
  return dphi;
}

void FillCorrelations(TClonesArray* vTrg, TClonesArray* vAssoc, Float_t cent, Float_t zvtx, THnD* hTrackHist, THnD* hEventHist, Bool_t ptOrdering, Bool_t mixing, Float_t dphimin, Float_t dphimax){
  Int_t ix[kTrackNvar];
  Int_t iv[kEventNvar];
  ix[kTrackCent] = axis[kTrackCent]->FindFixBin(cent);
  ix[kTrackZvtx] = axis[kTrackZvtx]->FindFixBin(zvtx);
  iv[kEventCent] = ix[kTrackCent];
  iv[kEventZvtx] = ix[kTrackZvtx];

  Float_t pt1,eta1,phi1,pt2,eta2,phi2;
  Double_t eff1=1.,eff2=1.;
  UInt_t id1,id2;
  AliCFParticle* d1;
  AliCFParticle* d2;

  for (Int_t i=0;i<vTrg->GetEntriesFast();i++){
    d1 = (AliCFParticle*) vTrg->UncheckedAt(i);
    id1  = d1->GetUniqueID();
    pt1  = d1->Pt();
    eta1 = d1->Eta();
    phi1 = d1->Phi();
    if(gFilter3D){eff1 = gFilter3D->GetBinContent(gFilter3D->FindFixBin(pt1,eta1,zvtx));}
    ix[kTrackPtTr] = axis[kTrackPtTr]->FindFixBin(pt1);
    iv[kEventPtTr] = ix[kTrackPtTr];
    //Double_t eff = gFilter? gFilter->GetBinContent(gFilter->FindFixBin(pt1,eta1)) : 1.;
    if (eff1<1.e-10) continue;
    hEventHist->AddBinContent(iv,1./eff1);
    hEventHist->AddBinError2(hEventHist->GetBin(iv),1./eff1/eff1);

    for (Int_t j=0;j<vAssoc->GetEntriesFast();j++){
      d2 = (AliCFParticle*) vAssoc->UncheckedAt(j);
      id2 = d2->GetUniqueID();
      if (!mixing && id1==id2) continue; // check unique id
      pt2 = d2->Pt();
      if (ptOrdering) if(pt2>pt1) continue;
      eta2 = d2->Eta();
      phi2 = d2->Phi();
      ix[kTrackPtAs] = axis[kTrackPtAs]->FindFixBin(pt2);
      ix[kTrackDeta] = axis[kTrackDeta]->FindFixBin(eta1-eta2);
      ix[kTrackDphi] = axis[kTrackDphi]->FindFixBin(dphi(phi1,phi2,dphimin,dphimax));
      if(gFilter3D){eff2 = gFilter3D->GetBinContent(gFilter3D->FindFixBin(pt2,eta2,zvtx));}
      if (eff2<1.e-10) continue;
      //atLeastOnePair=1;
      hTrackHist->AddBinContent(ix,1./eff1/eff2);
      hTrackHist->AddBinError2(hTrackHist->GetBin(ix),1./(eff1*eff1+eff2*eff2)); // TODO: CHECK THIS!
    }
  }
}


Int_t SelectTracks(TClonesArray* tracksAll, TClonesArray* tracksSelected, Int_t trackType, Float_t etaMin, Float_t etaMax, Float_t ptMin, Float_t ptMax, kFB fb) {
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
    if(charge==0 && trackType!=kTkl) continue;
    if       (trackType==kTrk)
      {
	if (fb==fb56 && !(mask & (1<<5|1<<6))) continue;
	if (fb==fb89 && !(mask & (1<<8|1<<9))) continue;
	if (fb==fb4 && !(mask & (1<<4))) continue;
	if (fb==fb1 && !(mask & (1<<1))) continue;
      }
    else if  (trackType==kTrkITS) { if (!(mask==0)) continue; }
    else if  (trackType==kTkl) { 
      dphi = pt;
      phi  = phiCorrected(phi,dphi);
      pt   = TMath::Abs(dphi)*1000; // mrad, pseudo pt for tracklets
      if (useMCforTracklets) { 
        if (!mask || !charge) continue; // skip secondaries and fakes
        if (mask == 999) continue; // fakes
        eta = track->GetAt(1); // etaMC
        phi = track->GetAt(2); // phiMC
        if(cutPtMCforTracklets){
          Double_t ptMC=track->GetAt(0); // ptMC
          if(ptMC<0.7)continue;
        }
      }
    }
    else if (trackType==kTrkGen || trackType==kTrkITSGen){
      Bool_t isPrimary = (Bool_t)track->GetAt(2);
      if(!isPrimary) continue;
    }

    if (eta<etaMin || eta>etaMax) continue;
    if (pt<ptMin || pt>ptMax) continue;

    AliCFParticle* part = new ((*tracksSelected)[nSelectedTracks++]) AliCFParticle(pt,eta,phi,charge,mask);
    part->SetUniqueID(10*i+trackType);
  }

  return nSelectedTracks;
}
