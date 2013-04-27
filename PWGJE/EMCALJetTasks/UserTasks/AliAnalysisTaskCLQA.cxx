// $Id: AliAnalysisTaskCLQA.cxx 60694 2013-02-04 15:35:56Z morsch $
//
// Constantin's Task
//
// Author: C.Loizides

#include <TChain.h>
#include <TClonesArray.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TNtupleD.h>
#include <TTree.h>

#include "AliESDMuonTrack.h"
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskCLQA.h"
#include "AliAnalysisUtils.h"
#include "AliCentrality.h"
#include "AliEMCALGeoParams.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliEmcalJet.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliTrackerBase.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliVTrack.h"

ClassImp(AliAnalysisTaskCLQA)

//________________________________________________________________________
AliAnalysisTaskCLQA::AliAnalysisTaskCLQA() : 
  AliAnalysisTaskEmcal("AliAnalysisTaskCLQA", kTRUE),
  fDoTracking(1), fDoMuonTracking(0), fDoCumulants(0), 
  fCumPtMin(0.3), fCumPtMax(5.0), fCumEtaMin(-1.0), fCumEtaMax(1.0), fCumMmin(15),
  fCentCL1In(0), fCentV0AIn(0),
  fNtupCum(0), fNtupCumInfo(0), fNtupZdcInfo(0)
{
  // Default constructor.

  for (Int_t i=0;i<1000;++i)
    fHists[i] = 0;
}

//________________________________________________________________________
AliAnalysisTaskCLQA::AliAnalysisTaskCLQA(const char *name) : 
  AliAnalysisTaskEmcal(name, kTRUE),
  fDoTracking(1), fDoMuonTracking(0), fDoCumulants(0), 
  fCumPtMin(0.3), fCumPtMax(5.0), fCumEtaMin(-1.0), fCumEtaMax(1.0), fCumMmin(15),
  fCentCL1In(0), fCentV0AIn(0),
  fNtupCum(0), fNtupCumInfo(0), fNtupZdcInfo(0)
{
  // Standard constructor.

  for (Int_t i=0;i<1000;++i)
    fHists[i] = 0;
}

//________________________________________________________________________
AliAnalysisTaskCLQA::~AliAnalysisTaskCLQA()
{
  // Destructor
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCLQA::FillHistograms()
{
  // Fill histograms.

  AliVEvent *event        = InputEvent();
  AliAnalysisManager *am  = AliAnalysisManager::GetAnalysisManager();

  UInt_t trig = ((AliInputEventHandler*)(am->GetInputEventHandler()))->IsEventSelected();
  for (Int_t i=0;i<31;++i) {
    if (trig & (1<<i))
      fHists[0]->Fill(trig);
  }

  Double_t vz  = event->GetPrimaryVertex()->GetZ();
  fHists[1]->Fill(vz);

  Int_t  run = event->GetRunNumber();
  Int_t  vzn = event->GetPrimaryVertex()->GetNContributors();
  if ((vzn<1)&&(run>0))
    return kTRUE;
  fHists[2]->Fill(vz);

  if (TMath::Abs(vz)>10)
    return kTRUE;

  if ((run>=188356&&run<=188366) || (run>=195344&&run<=197388)) {
    AliAnalysisUtils anau;
    if (anau.IsFirstEventInChunk(event))
      return kFALSE;
    if (!anau.IsVertexSelected2013pA(event))
      return kFALSE;
  }

  // accepted events
  fHists[9]->Fill(1);

  AliCentrality *cent = InputEvent()->GetCentrality();
  Double_t v0acent = cent->GetCentralityPercentile("V0A");
  fHists[10]->Fill(v0acent);
  Double_t znacent = cent->GetCentralityPercentile("ZNA");
  fHists[11]->Fill(znacent);

  if (fDoTracking) {
    const Int_t ntracks = fTracks->GetEntries();
    if (fTracks) {
      for (Int_t i =0; i<ntracks; ++i) {
        AliVParticle *track = dynamic_cast<AliVParticle*>(fTracks->At(i));
        if (!track)
          continue;
        if (track->Charge()==0)
          continue;
        Double_t phi = track->Phi();
        Double_t eta = track->Eta();
        Double_t pt  = track->Pt();
        fHists[20]->Fill(phi,eta);
        fHists[21]->Fill(phi,pt);
        fHists[22]->Fill(eta,pt);
        if (TMath::Abs(eta)<0.8) {
          fHists[23]->Fill(pt,v0acent);
          fHists[24]->Fill(pt,znacent);
        }
      }
    }
  }

  if (fDoMuonTracking) {
    AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
    if (aod) {
      for (Int_t iMu = 0; iMu<aod->GetNumberOfTracks(); ++iMu) {
        AliAODTrack* muonTrack = aod->GetTrack(iMu);
        if (!muonTrack)
          continue;
        if (!muonTrack->IsMuonTrack()) 
          continue;
        Double_t dThetaAbs = TMath::ATan(muonTrack->GetRAtAbsorberEnd()/505.)* TMath::RadToDeg();
        if ((dThetaAbs<2.) || (dThetaAbs>10.)) 
          continue;
        Double_t dEta = muonTrack->Eta();
        if ((dEta<-4.) || (dEta>-2.5)) 
          continue;
        if (0) {
          if (muonTrack->GetMatchTrigger()<0.5) 
            continue;
        }
        Double_t ptMu  = muonTrack->Pt();
        Double_t etaMu = muonTrack->Eta();
        Double_t phiMu = muonTrack->Phi();
        fHists[50]->Fill(phiMu,etaMu);
        fHists[51]->Fill(phiMu,ptMu);
        fHists[52]->Fill(etaMu,ptMu);
        fHists[53]->Fill(ptMu,v0acent);
        fHists[54]->Fill(ptMu,znacent);
      }
    } else {
      AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
      if (esd) {
        for (Int_t iMu = 0; iMu<esd->GetNumberOfMuonTracks(); ++iMu) {
          AliESDMuonTrack* muonTrack = esd->GetMuonTrack(iMu);
          if (!muonTrack)
            continue;
          if (!muonTrack->ContainTrackerData()) 
            continue;
          Double_t thetaTrackAbsEnd = TMath::ATan(muonTrack->GetRAtAbsorberEnd()/505.) * TMath::RadToDeg();
          if ((thetaTrackAbsEnd < 2.) || (thetaTrackAbsEnd > 10.)) 
            continue;
          Double_t eta = muonTrack->Eta();
          if ((eta < -4.) || (eta > -2.5))
            return kFALSE;
          if (0) {
            if (!muonTrack->ContainTriggerData()) 
              continue;
            if (muonTrack->GetMatchTrigger() < 0.5) 
              continue;
          }

        }
      }
    }
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCLQA::RetrieveEventObjects()
{
  // Retrieve event objects.

  if (!AliAnalysisTaskEmcal::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCLQA::Run()
{
  // Run various functions.

  RunCumulants(fCumMmin,fCumPtMin,fCumPtMax,fCumEtaMin,fCumEtaMax);

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskCLQA::RunCumulants(Double_t Mmin, Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax)
{
  // Run cumulant analysis.

  if (!fDoCumulants)
    return;

  if (!fTracks) 
    return;

  Bool_t isMC = 0;
  TString tname(fTracks->GetName());
  if (tname.Contains("mc"))
    isMC = 1;

  const Int_t ntracks = fTracks->GetEntries();
  Int_t Mall=0,M=0,Mall2=0;
  Double_t ptmaxall=0,ptsumall=0,pt2sumall=0,ptsumall2=0;
  Double_t tsa00=0,tsa10=0,tsa11=0;
  Double_t Q2r=0,Q2i=0;
  Double_t Q4r=0,Q4i=0;
  Double_t mpt=0,mpt2=0,ptmaxq=0;
  Double_t ts00=0,ts10=0,ts11=0;
  Double_t v0ach=0, v0cch=0;
  Double_t cl1ch=0;
 
  for (Int_t i =0; i<ntracks; ++i) {
    AliVParticle *track = dynamic_cast<AliVParticle*>(fTracks->At(i));
    if (!track)
      continue;
    if (track->Charge()==0)
      continue;
    Double_t eta = track->Eta();
    if ((eta<5.1)&&(eta>2.8))
      ++v0ach;
    else if ((eta>-3.7)&&(eta<-1.7))
      ++v0cch;
    if (TMath::Abs(eta)<1.4) {
      ++cl1ch;
    }
    if ((eta<etamin) || (eta>etamax))
      continue;
    Double_t pt = track->Pt();
    if (pt>ptmaxall)
      ptmaxall = pt;
    if (pt>1) {
      ptsumall2 += pt;
      ++Mall2;
    }
    ptsumall  +=pt;
    pt2sumall +=pt*pt;
    Double_t px = track->Px();
    Double_t py = track->Py();
    tsa00 += px*px/pt;
    tsa10 += px*py/pt;
    tsa11 += py*py/pt;
    ++Mall;
    if ((pt<ptmin) || (pt>ptmax))
      continue;
    if (pt>ptmaxq)
      ptmaxq = pt;
    Double_t phi = track->Phi();
    ++M;
    mpt  += pt;
    mpt2 += pt*pt;
    ts00 += px*px/pt;
    ts10 += px*py/pt;
    ts11 += py*py/pt;
    Q2r  += TMath::Cos(2*phi);
    Q2i  += TMath::Sin(2*phi);
    Q4r  += TMath::Cos(4*phi);
    Q4i  += TMath::Sin(4*phi);
  }

  if (M<Mmin)
    return;

  Double_t Q2abs = Q2r*Q2r+Q2i*Q2i;
  Double_t Q4abs = Q4r*Q4r+Q4i*Q4i;
  Double_t Q42re = Q4r*Q2r*Q2r-Q4r*Q2i*Q2i+2*Q4i*Q2r*Q2i;
  Double_t tsall = -1;
  Double_t tsax = (tsa00+tsa11)*(tsa00+tsa11)-4*(tsa00*tsa11-tsa10*tsa10);
  if (tsax>=0) {
    Double_t l1 = 0.5*(tsa00+tsa11+TMath::Sqrt(tsax))/ptsumall;
    Double_t l2 = 0.5*(tsa00+tsa11-TMath::Sqrt(tsax))/ptsumall;
    tsall = 2*l2/(l1+l2);
  }

  Double_t ts = -1;
  Double_t tsx = (ts00+ts11)*(ts00+ts11)-4*(ts00*ts11-ts10*ts10);
  if (tsx>=0) {
    Double_t l1 = 0.5*(ts00+ts11+TMath::Sqrt(tsx))/ptsumall;
    Double_t l2 = 0.5*(ts00+ts11-TMath::Sqrt(tsx))/ptsumall;
    ts = 2*l2/(l1+l2);
  }

  if (isMC) {
    fHists[106]->Fill(cl1ch);
    fHists[107]->Fill(v0ach);
    fHists[108]->Fill(v0cch);
    AliCentrality *cent = InputEvent()->GetCentrality();
    if (fCentCL1In) {
      cent->SetCentralityCL1(100*fCentCL1In->GetBinContent(fCentCL1In->FindBin(cl1ch)));
      cent->SetQuality(0);
    }
    if (fCentV0AIn) {
      cent->SetCentralityV0A(100*fCentV0AIn->GetBinContent(fCentV0AIn->FindBin(v0ach)));
      cent->SetQuality(0);
    }
 }

  AliAnalysisUtils anau;
  AliVEvent *event        = InputEvent();
  AliAnalysisManager *am  = AliAnalysisManager::GetAnalysisManager();

  fNtupCumInfo->fTrig     = ((AliInputEventHandler*)(am->GetInputEventHandler()))->IsEventSelected();
  fNtupCumInfo->fRun      = event->GetRunNumber();
  fNtupCumInfo->fVz       = event->GetPrimaryVertex()->GetZ();
  fNtupCumInfo->fIsFEC    = anau.IsFirstEventInChunk(event);
  fNtupCumInfo->fIsVSel   = anau.IsVertexSelected2013pA(event);
  fNtupCumInfo->fIsP      = event->IsPileupFromSPD(3/*minContributors*/,
                                                   0.8/*minZdist*/,
                                                   3./*nSigmaZdist*/,
                                                   2./*nSigmaDiamXY*/,
                                                   5./*nSigmaDiamZ*/);

  fNtupCumInfo->fMall     = Mall;
  fNtupCumInfo->fMall2    = Mall2;
  fNtupCumInfo->fPtMaxall = ptmaxall;
  fNtupCumInfo->fMPtall   = ptsumall/Mall;
  fNtupCumInfo->fMPt2all  = pt2sumall/Mall;
  if (Mall2>0)
    fNtupCumInfo->fMPtall2  = ptsumall2/Mall2;
  else
    fNtupCumInfo->fMPtall2  = -1;
  fNtupCumInfo->fTSall    = tsall;
  fNtupCumInfo->fM        = M;
  fNtupCumInfo->fQ2abs    = Q2abs;
  fNtupCumInfo->fQ4abs    = Q4abs;
  fNtupCumInfo->fQ42re    = Q42re;
  fNtupCumInfo->fCos2phi  = Q2r;
  fNtupCumInfo->fSin2phi  = Q2i;
  fNtupCumInfo->fPtMax    = ptmaxq;
  fNtupCumInfo->fMPt      = mpt/M;
  fNtupCumInfo->fMPt2     = mpt2/M;
  fNtupCumInfo->fTS       = ts;

  if (isMC) {
    fNtupCumInfo->fMV0M     = v0ach + v0cch;
  } else {
    AliVVZERO *vzero = InputEvent()->GetVZEROData();
    fNtupCumInfo->fMV0M     = vzero->GetMTotV0A()+vzero->GetMTotV0C();
  }

  AliCentrality *cent = InputEvent()->GetCentrality();
  fNtupCumInfo->fCl1      = cent->GetCentralityPercentile("CL1");
  fNtupCumInfo->fV0M      = cent->GetCentralityPercentile("V0M");
  fNtupCumInfo->fV0MEq    = cent->GetCentralityPercentile("V0MEq");
  fNtupCumInfo->fV0A      = cent->GetCentralityPercentile("V0A");
  fNtupCumInfo->fV0AEq    = cent->GetCentralityPercentile("V0AEq");
  fNtupCumInfo->fZNA      = cent->GetCentralityPercentile("ZNA");

  AliVZDC *vZDC = InputEvent()->GetZDCData();
  const Double_t *znaTowers = vZDC->GetZNATowerEnergy(); 
  fNtupZdcInfo->fZna0 = znaTowers[0];
  fNtupZdcInfo->fZna1 = znaTowers[1];
  fNtupZdcInfo->fZna2 = znaTowers[2];
  fNtupZdcInfo->fZna3 = znaTowers[3];
  fNtupZdcInfo->fZna4 = znaTowers[4];

  fNtupCum->Fill();

  fHists[109]->Fill(fNtupCumInfo->fCl1);
  fHists[110]->Fill(fNtupCumInfo->fV0A);
  fHists[111]->Fill(fNtupCumInfo->fZNA);

  if ((isMC) || 
      ((TMath::Abs(fNtupCumInfo->fVz)<10) && !fNtupCumInfo->fIsFEC && fNtupCumInfo->fIsVSel)) {
    for (Int_t i =0; i<ntracks; ++i) {
      AliVParticle *track1 = dynamic_cast<AliVParticle*>(fTracks->At(i));
      if (!track1)
        continue;
      Double_t phi1 = track1->Phi();
      Double_t eta1 = track1->Eta();
      Double_t pt1  = track1->Pt();
      ((TH3*)fHists[103])->Fill(pt1,eta1,fNtupCumInfo->fCl1);
      ((TH3*)fHists[104])->Fill(pt1,eta1,fNtupCumInfo->fV0A);
      ((TH3*)fHists[105])->Fill(pt1,eta1,fNtupCumInfo->fZNA);
      if ((eta1<etamin) || (eta1>etamax))
        continue;
      if ((pt1<ptmin) || (pt1>ptmax))
        continue;
      for (Int_t j =0; j<ntracks; ++j) {
        AliVParticle *track2 = dynamic_cast<AliVParticle*>(fTracks->At(j));
        if (!track2)
          continue;
        Double_t eta2 = track2->Eta();
        if ((eta2<etamin) || (eta2>etamax))
          continue;
        Double_t pt2 = track2->Pt();
        if ((pt2<ptmin) || (pt2>ptmax))
          continue;
        Double_t phi2 = track2->Phi();
        Double_t deta = eta1-eta2;
        Double_t dphi = phi1-phi2;
        while (dphi<-TMath::Pi())
          dphi+=TMath::TwoPi();
        while (dphi>3*TMath::Pi()/2)
          dphi-=TMath::TwoPi();
        ((TH3*)fHists[100])->Fill(dphi,deta,fNtupCumInfo->fCl1);
        ((TH3*)fHists[101])->Fill(dphi,deta,fNtupCumInfo->fV0A);
        ((TH3*)fHists[102])->Fill(dphi,deta,fNtupCumInfo->fZNA);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskCLQA::SetCumParams(Double_t Mmin, Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax)
{
  // Set parameters for cumulants.

  fCumMmin   = Mmin;
  fCumPtMin  = ptmin;
  fCumPtMax  = ptmax;
  fCumEtaMin = etamin;
  fCumEtaMax = etamax;
}

//________________________________________________________________________
void AliAnalysisTaskCLQA::UserCreateOutputObjects()
{
  // Create histograms

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  fHists[0] = new TH1D("fTrigBits",";bit",32,-0.5,31.5);
  fOutput->Add(fHists[0]);
  fHists[1] = new TH1D("fVertexZ",";vertex z (cm)",51,-25.5,25.5);
  fOutput->Add(fHists[1]);
  fHists[2] = new TH1D("fVertexZnc",";vertex z (cm)",51,-25.5,25.5);
  fOutput->Add(fHists[2]);
  fHists[9] = new TH1D("fAccepted","",1,0.5,1.5);
  fOutput->Add(fHists[9]);
  fHists[10] = new TH1D("fV0ACent",";percentile",20,0,100);
  fOutput->Add(fHists[10]);
  fHists[11] = new TH1D("fZNACent",";percentile",20,0,100);
  fOutput->Add(fHists[11]);

  if (fDoTracking) {
    fHists[20] = new TH2D("fPhiEtaTracks",";#phi;#eta",60,0,TMath::TwoPi(),20,-2,2);
    fHists[20]->Sumw2();
    fOutput->Add(fHists[20]);
    fHists[21] = new TH2D("fPhiPtTracks",";#phi;p_{T} (GeV/c)",60,0,TMath::TwoPi(),40,0,20);
    fHists[21]->Sumw2();
    fOutput->Add(fHists[21]);
    fHists[22] = new TH2D("fEtaPtTracks",";#eta;p_{T} (GeV/c)",20,-2,2,40,0,20);
    fHists[22]->Sumw2();
    fOutput->Add(fHists[22]);
    fHists[23] = new TH2D("fPtV0ATracks",";#p_{T} (GeV/c);percentile",100,0,20,20,0,100);
    fHists[23]->Sumw2();
    fOutput->Add(fHists[23]);
    fHists[24] = new TH2D("fPtZNATracks",";#p_{T} (GeV/c);percentile",100,0,20,20,0,100);
    fHists[24]->Sumw2();
    fOutput->Add(fHists[24]);
  }

  if (fDoMuonTracking) {
    fHists[50] = new TH2D("fPhiEtaMuonTracks",";#phi;#eta",60,0,TMath::TwoPi(),15,-4,-2.5);
    fHists[50]->Sumw2();
    fOutput->Add(fHists[50]);
    fHists[51] = new TH2D("fPhiPtMuonTracks",";#phi;p_{T} (GeV/c)",60,0,TMath::TwoPi(),200,0,20);
    fHists[51]->Sumw2();
    fOutput->Add(fHists[51]);
    fHists[52] = new TH2D("fEtaPtMuonTracks",";#eta;p_{T} (GeV/c)",15,-4,-2.5,200,0,20);
    fHists[52]->Sumw2();
    fOutput->Add(fHists[52]);
    fHists[53] = new TH2D("fPtV0AMuonTracks",";#p_{T} (GeV/c);percentile",100,0,10,20,0,100);
    fHists[53]->Sumw2();
    fOutput->Add(fHists[53]);
    fHists[54] = new TH2D("fPtZNAMuonTracks",";#p_{T} (GeV/c);percentile",100,0,10,20,0,100);
    fHists[54]->Sumw2();
    fOutput->Add(fHists[54]);
  }

  if (fDoCumulants) {
    fNtupCum = new TTree("NtupCum", "Ntuple for cumulant analysis");
    fHists[100] = new TH3D("fCumPhiEtaCl1",";#Delta#phi;#Delta#eta",32,-TMath::Pi()/2,3*TMath::Pi()/2,60,-3,3,10,0,100);
    fOutput->Add(fHists[100]);
    fHists[101] = new TH3D("fCumPhiEtaV0A",";#Delta#phi;#Delta#eta",32,-TMath::Pi()/2,3*TMath::Pi()/2,60,-3,3,10,0,100);
    fOutput->Add(fHists[101]);
    fHists[102] = new TH3D("fCumPhiEtaZNA",";#Delta#phi;#Delta#eta",32,-TMath::Pi()/2,3*TMath::Pi()/2,60,-3,3,10,0,100);
    fOutput->Add(fHists[102]);
    fHists[103] = new TH3D("fCumPtEtaCl1",";p_{T} (GeV/c);#eta",100,0,25,20,-2,2,10,0,100);
    fOutput->Add(fHists[103]);
    fHists[104] = new TH3D("fCumPtEtaV0A",";p_{T} (GeV/c);#eta",100,0,25,20,-2,2,10,0,100);
    fOutput->Add(fHists[104]);
    fHists[105] = new TH3D("fCumPtEtaZNA",";p_{T} (GeV/c);#eta",100,0,25,20,-2,2,10,0,100);
    fOutput->Add(fHists[105]);
    fHists[106] = new TH1D("fCumCL1MC",";#tracks",2500,0,2500);
    fOutput->Add(fHists[106]);
    fHists[107] = new TH1D("fCumV0AMC",";#tracks",2500,0,2500);
    fOutput->Add(fHists[107]);
    fHists[108] = new TH1D("fCumV0CMC",";#tracks",2500,0,2500);
    fOutput->Add(fHists[108]);
    fHists[109] = new TH1D("fCumCl1Cent",";percentile",10,0,100);
    fOutput->Add(fHists[109]);
    fHists[110] = new TH1D("fCumV0ACent",";percentile",10,0,100);
    fOutput->Add(fHists[110]);
    fHists[111] = new TH1D("fCumZNACent",";percentile",10,0,100);
    fOutput->Add(fHists[111]);

    if (1) {
      fNtupCum->SetDirectory(0);
    } else {
      TFile *f = OpenFile(1); 
      if (f) {
        f->SetCompressionLevel(2);
        fNtupCum->SetDirectory(f);
        fNtupCum->SetAutoFlush(-4*1024*1024);
        fNtupCum->SetAutoSave(0);
      }
    }
    fNtupCumInfo = new AliNtupCumInfo;
    fNtupCum->Branch("cumulants", &fNtupCumInfo, 32*1024, 99);
    fNtupZdcInfo = new AliNtupZdcInfo;
    fNtupCum->Branch("zdc", &fNtupZdcInfo, 32*1024, 99);

    fOutput->Add(fNtupCum);
  }

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}
