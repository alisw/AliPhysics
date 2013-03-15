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
  AliAnalysisTaskEmcalJet("AliAnalysisTaskCLQA", kTRUE),
  fDoCumulants(0), 
  fCumPtMin(0.3), fCumPtMax(5.0), fCumEtaMin(-1.0), fCumEtaMax(1.0), fCumMmin(15),
  fNtupCum(0), fNtupCumInfo(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliAnalysisTaskCLQA::AliAnalysisTaskCLQA(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fDoCumulants(0), 
  fCumPtMin(0.3), fCumPtMax(5.0), fCumEtaMin(-1.0), fCumEtaMax(1.0), fCumMmin(15),
  fNtupCum(0), fNtupCumInfo(0)
{
  // Standard constructor.
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

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCLQA::RetrieveEventObjects()
{
  // Retrieve event objects.

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
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

  const Int_t ntracks = fTracks->GetEntries();
  Int_t Mall=0,M=0,Mall2=0;
  Double_t ptmaxall=0,ptsumall=0,pt2sumall=0,ptsumall2=0;
  Double_t tsa00=0,tsa10=0,tsa11=0;
  Double_t Q2r=0,Q2i=0;
  Double_t Q4r=0,Q4i=0;
  Double_t mpt=0,mpt2=0,ptmaxq=0;
  Double_t ts00=0,ts10=0,ts11=0;
  for (Int_t i =0; i<ntracks; ++i) {
    AliVTrack *track = dynamic_cast<AliVTrack*>(fTracks->At(i));
    Double_t eta = track->Eta();
    if ((eta<etamin) || (eta>etamax))
      continue;
    Double_t pt = track->Pt();
    if (pt>ptmaxall)
      ptmaxall = pt;
    if (pt>2) {
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
  fNtupCumInfo->fMPtall2  = ptsumall2/Mall2;
  fNtupCumInfo->fTSall    = tsall;
  fNtupCumInfo->fM        = M;
  fNtupCumInfo->fQ2abs    = Q2abs;
  fNtupCumInfo->fQ4abs    = Q4abs;
  fNtupCumInfo->fQ42re    = Q42re;
  fNtupCumInfo->fPtMax    = ptmaxq;
  fNtupCumInfo->fMPt      = mpt/M;
  fNtupCumInfo->fMPt2     = mpt2/M;
  fNtupCumInfo->fTS       = ts;
  AliVVZERO *vzero = InputEvent()->GetVZEROData();
  fNtupCumInfo->fMV0M     = vzero->GetMTotV0A()+vzero->GetMTotV0C();
  AliCentrality *cent = InputEvent()->GetCentrality();
  fNtupCumInfo->fCl1      = cent->GetCentralityPercentile("CL1");
  fNtupCumInfo->fV0M      = cent->GetCentralityPercentile("V0M");
  fNtupCumInfo->fV0MEq    = cent->GetCentralityPercentile("V0MEq");
  fNtupCumInfo->fV0A      = cent->GetCentralityPercentile("V0A");
  fNtupCumInfo->fV0AEq    = cent->GetCentralityPercentile("V0AEq");
  fNtupCumInfo->fZNA      = cent->GetCentralityPercentile("ZNA");
  fNtupCum->Fill();
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

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  if (fDoCumulants) {
    fNtupCum = new TTree("NtupCum", "Ntuple for cumulant analysis");
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
    fOutput->Add(fNtupCum);
  }

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}
