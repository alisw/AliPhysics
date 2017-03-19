// $Id: $
//
// Constantin's Task
//
// Author: C.Loizides

#include <complex>

#include "AliAnalysisTaskCLQA.h"
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
#include <TProfile.h>
#include <TTree.h>
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliCentrality.h"
#include "AliEMCALGeoParams.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliEmcalJet.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMultSelection.h"
#include "AliPicoTrack.h"
#include "AliTrackerBase.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"

ClassImp(AliAnalysisTaskCLQA)

//________________________________________________________________________
AliAnalysisTaskCLQA::AliAnalysisTaskCLQA() : 
  AliAnalysisTaskEmcal("AliAnalysisTaskCLQA", kTRUE),
  fDoVertexCut(1),
  fDoTracking(0), fDoMuonTracking(0), fDoCumulants(0), fDoCumNtuple(0), fDoProp(0),
  fCumPtMin(0.3), fCumPtMax(5.0), fCumEtaMin(-1.0), fCumEtaMax(1.0), fCumMmin(15), fCumMbins(250), 
  fDoHet(0), fQC4EG(-1), fHetEtmin(6),
  fCentCL1In(0), fCentV0AIn(0),
  fNtupCum(0), fNtupCumInfo(0), fNtupZdcInfo(0), 
  fNtupHet(0), fNtupHetInfo(0), fCum(0)
{
  // Default constructor.

  for (Int_t i=0;i<1000;++i)
    fHists[i] = 0;
}

//________________________________________________________________________
AliAnalysisTaskCLQA::AliAnalysisTaskCLQA(const char *name) : 
  AliAnalysisTaskEmcal(name, kTRUE),
  fDoVertexCut(1),
  fDoTracking(1), fDoMuonTracking(0), fDoCumulants(0), fDoCumNtuple(0), fDoProp(0),
  fCumPtMin(0.3), fCumPtMax(5.0), fCumEtaMin(-1.0), fCumEtaMax(1.0), fCumMmin(15), fCumMbins(250), 
  fDoHet(0), fQC4EG(-1), fHetEtmin(6),
  fCentCL1In(0), fCentV0AIn(0),
  fNtupCum(0), fNtupCumInfo(0), fNtupZdcInfo(0), 
  fNtupHet(0), fNtupHetInfo(0), fCum(0)
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
Double_t AliAnalysisTaskCLQA::DeltaPhi(Double_t phia, Double_t phib,
                                       Double_t rangeMin, Double_t rangeMax) const
{
  // Calculate Delta Phi.

  Double_t dphi = -999;
  const Double_t tpi = TMath::TwoPi();
  
  if (phia < 0)         phia += tpi;
  else if (phia > tpi) phia -= tpi;
  if (phib < 0)         phib += tpi;
  else if (phib > tpi) phib -= tpi;
  dphi = phib - phia;
  if (dphi < rangeMin)      dphi += tpi;
  else if (dphi > rangeMax) dphi -= tpi;
  
  return dphi;
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
      fHists[0]->Fill(i);
  }

  Double_t vz  = event->GetPrimaryVertex()->GetZ();
  fHists[1]->Fill(vz);

  Int_t  run = event->GetRunNumber();
  Int_t  vzn = event->GetPrimaryVertex()->GetNContributors();
  if ((vzn<1)&&(run>0))
    return kFALSE;
  fHists[2]->Fill(vz);

  if (TMath::Abs(vz)>10)
    return kFALSE;

  if (fDoVertexCut) {
    if ((run>=188356&&run<=188366) || (run>=195344&&run<=197388)) {
      AliAnalysisUtils anau;
      if (anau.IsFirstEventInChunk(event))
	return kFALSE;
      if (!anau.IsVertexSelected2013pA(event))
	return kFALSE;
    } else if (run>=260014) {
      //https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsEventProp#Pileup_removal
      const Double_t SPDZDiffCut=0.8; // (vertices with z separation lower than 8 mm are not considered for tagging)
      Int_t SPDContributorsCut=5; //3 for low multiplicity pp events (higher efficiency, but also higher contamination from false positives), fSPDContributorsCut=5 at high multiplicity 
      Bool_t isPileupFromSPD=event->IsPileupFromSPD(SPDContributorsCut,SPDZDiffCut,3.,2.,5.);
      if (isPileupFromSPD)
	fHists[3]->Fill(1);
      const Double_t minContributors=5; 
      const Double_t minChi2=5.; 
      const Double_t minWeiZDiff=15; 
      const Bool_t checkPlpFromDifferentBC=kFALSE; 
      AliAnalysisUtils anau;
      anau.SetMinPlpContribMV(minContributors);
      anau.SetMaxPlpChi2MV(minChi2);
      anau.SetMinWDistMV(minWeiZDiff);
      anau.SetCheckPlpFromDifferentBCMV(checkPlpFromDifferentBC);
      Bool_t isPileupFromMV = anau.IsPileUpMV(event);
      if (isPileupFromMV)
	fHists[3]->Fill(2);
      if (isPileupFromSPD||isPileupFromMV)
	return kFALSE;
    }
  }

  // accepted events
  fHists[9]->Fill(1,run);

  AliMultSelection *ms = dynamic_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  Double_t v0acent=-1,znacent=-1,v0mcent=-1;
  if (ms) {
    v0acent = ms->GetMultiplicityPercentile("V0A");
    znacent = ms->GetMultiplicityPercentile("ZNA");
    v0mcent = ms->GetMultiplicityPercentile("V0M");
  } else {
    AliCentrality *cent = InputEvent()->GetCentrality();
    v0acent = cent->GetCentralityPercentile("V0A");
    znacent = cent->GetCentralityPercentile("ZNA");
    v0mcent = cent->GetCentralityPercentile("V0M");
  }
  fHists[10]->Fill(v0acent);
  fHists[11]->Fill(znacent);
  fHists[12]->Fill(v0mcent);

  if (fDoTracking) {
    const Int_t ntracks = fTracks->GetEntries();
    if (fTracks) {
      for (Int_t i=0; i<ntracks; ++i) {
        AliVTrack *track = dynamic_cast<AliVTrack*>(fTracks->At(i));
        if (!track)
          continue;
        if (track->Charge()==0)
          continue;

	AliPicoTrack *picot = dynamic_cast<AliPicoTrack*>(track);
	if (picot && picot->GetTrack()) 
	  track = picot->GetTrack();
	
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
	Int_t ttype = AliPicoTrack::GetTrackType(track);
	fHists[25+ttype]->Fill(phi,pt);
	if (fDoProp) {
	  if (track->IsExtrapolatedToEMCAL()) {
	    Double_t dphi = TVector2::Phi_mpi_pi(phi-track->GetTrackPhiOnEMCal());
	    fHists[28]->Fill(dphi,pt);
	  } else {
	    AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(track,440);
	    if (track->IsExtrapolatedToEMCAL()) {
	      Double_t dphi = TVector2::Phi_mpi_pi(phi-track->GetTrackPhiOnEMCal());
	      fHists[29]->Fill(dphi,pt);
	    }
	  }
	  if (track->IsEMCAL() && track->IsExtrapolatedToEMCAL()) {
	    Int_t id = track->GetEMCALcluster();
	    AliVCluster *clus = InputEvent()->GetCaloCluster(id);
	    if (id>=0&&clus) {
	      Float_t pos[3];
	      clus->GetPosition(pos);
	      TVector3 vpos(pos);
	      Double_t dphi = TVector2::Phi_mpi_pi(vpos.Phi()-track->GetTrackPhiOnEMCal());
	      fHists[30]->Fill(dphi,pt);
	    }
	  }
	  if (track->IsExtrapolatedToEMCAL()) {
	    Double_t phi1 = track->GetTrackPhiOnEMCal();
	    AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(track,440);
	    Double_t phi2 = track->GetTrackPhiOnEMCal();
	    Double_t dphi = TVector2::Phi_mpi_pi(phi1-phi2);
	    fHists[31]->Fill(dphi,pt);
	  }
	}
      }
    }
  }

  if (fDoMuonTracking) {
    AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
    if (aod) {
      for (Int_t iMu = 0; iMu<aod->GetNumberOfTracks(); ++iMu) {
        AliAODTrack* muonTrack = static_cast<AliAODTrack*>(aod->GetTrack(iMu));
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

  if (fCum) 
    RunCumulants();
  else 
    RunCumulants(fCumMmin,fCumPtMin,fCumPtMax,fCumEtaMin,fCumEtaMax);

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskCLQA::RunCumulants()
{
  // Run cumulant analysis.

  if (!fDoCumulants)
    return;

  if (!fTracks) 
    return;

  if (!fCum)
    return;

  TObjArray &objs = *fTracks;
  fCum->SetTracks(objs);
  fCum->RunAll();
 
  Int_t M=fCum->GetM();
  if (M<fCumMmin)
    return;
  AliVVZERO *vzero = InputEvent()->GetVZEROData();
  Double_t v0a = vzero->GetMTotV0A();
  Double_t v0c = vzero->GetMTotV0C();
  Double_t v0m = vzero->GetMTotV0A()+vzero->GetMTotV0C();
  fHists[117]->Fill(v0a,M);
  fHists[118]->Fill(v0c,M);
  fHists[119]->Fill(v0m,M);
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
  Double_t Q3r=0,Q3i=0;
  Double_t Q4r=0,Q4i=0;
  Double_t Q6r=0,Q6i=0;
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
    Q3r  += TMath::Cos(3*phi);
    Q3i  += TMath::Sin(3*phi);
    Q4r  += TMath::Cos(4*phi);
    Q4i  += TMath::Sin(4*phi);
    Q6r  += TMath::Cos(6*phi);
    Q6i  += TMath::Sin(6*phi);
  }

  if (M<=1)
    return;

  Int_t pmult=0;
  Double_t v2g=0;
  Double_t v3g=0;
  Int_t pmult14=0;
  Double_t v2g14=0;
  Double_t v3g14=0;
  Int_t pmult18=0;
  Double_t v2g18=0;
  Double_t v3g18=0;
  for (Int_t i=0; i<ntracks; ++i) {
    AliVParticle *track1 = dynamic_cast<AliVParticle*>(fTracks->At(i));
    if (!track1)
      continue;
    if (track1->Charge()==0)
      continue;
    Double_t eta1 = track1->Eta();
    if ((eta1<etamin) || (eta1>etamax))
      continue;
    Double_t pt1 = track1->Pt();
    if ((pt1<ptmin) || (pt1>ptmax))
      continue;
    Double_t phi1 = track1->Phi();
    for (Int_t j = i+1; j<ntracks; ++j) {
      AliVParticle *track2 = dynamic_cast<AliVParticle*>(fTracks->At(j));
      if (!track2)
	continue;
      if (track2->Charge()==0)
	continue;
      Double_t eta2 = track2->Eta();
      if ((eta2<etamin) || (eta2>etamax))
	continue;
      Double_t pt2 = track2->Pt();
      if ((pt2<ptmin) || (pt2>ptmax))
	continue;
      ((TH3*)fHists[128])->Fill(DeltaPhi(phi1,track2->Phi()),eta1-eta2,M);
      fHists[129]->Fill(M);
      Double_t deta=TMath::Abs(eta1-eta2);
      if(deta<1)
	continue;
      Double_t dphi=TVector2::Phi_0_2pi(phi1-track2->Phi());
      pmult++;
      v2g+=TMath::Cos(2*dphi);
      v3g+=TMath::Cos(3*dphi);
      if (deta>1.4) {
        pmult14++;
        v2g14+=TMath::Cos(2*dphi);
        v3g14+=TMath::Cos(3*dphi);
      }
      if (deta>1.8) {
        pmult18++;
        v2g18+=TMath::Cos(2*dphi);
        v3g18+=TMath::Cos(3*dphi);
      }
    }      
  }

  if (pmult>0) {
    v2g/=pmult;
    v3g/=pmult;
  }
  if (pmult14>0) {
    v2g14/=pmult14;
    v3g14/=pmult14;
  }
  if (pmult18>0) {
    v2g18/=pmult18;
    v3g18/=pmult18;
  }

  std::complex<double> q2(Q2r,Q2i);
  std::complex<double> q3(Q3r,Q3i);
  std::complex<double> q4(Q4r,Q4i);
  std::complex<double> q6(Q6r,Q6i);
  Double_t Q22   = std::abs(q2)*std::abs(q2);
  Double_t Q32   = std::abs(q3)*std::abs(q3);
  Double_t Q42   = std::abs(q4)*std::abs(q4);
  Double_t Q62   = std::abs(q6)*std::abs(q6);
  Double_t Q32re = std::real(q6*std::conj(q3)*std::conj(q3));
  Double_t Q42re = std::real(q4*std::conj(q2)*std::conj(q2));
  Double_t Q6are = std::real(q4*q2*std::conj(q2)*std::conj(q2)*std::conj(q2));
  Double_t Q6bre = std::real(q6*std::conj(q2)*std::conj(q2)*std::conj(q2));
  Double_t Q6cre = std::real(q6*std::conj(q4)*std::conj(q2));

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
  fNtupCumInfo->fQ2abs    = Q22;
  fNtupCumInfo->fQ4abs    = Q42;
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

  AliMultSelection *ms = dynamic_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  if (ms) {
    fNtupCumInfo->fCl1      = ms->GetMultiplicityPercentile("CL1");
    fNtupCumInfo->fV0M      = ms->GetMultiplicityPercentile("V0M");
    fNtupCumInfo->fV0MEq    = ms->GetMultiplicityPercentile("V0MEq");
    fNtupCumInfo->fV0A      = ms->GetMultiplicityPercentile("V0A");
    fNtupCumInfo->fV0AEq    = ms->GetMultiplicityPercentile("V0AEq");
    fNtupCumInfo->fZNA      = ms->GetMultiplicityPercentile("ZNA");
  } else {
    AliCentrality *cent = InputEvent()->GetCentrality();
    fNtupCumInfo->fCl1      = cent->GetCentralityPercentile("CL1");
    fNtupCumInfo->fV0M      = cent->GetCentralityPercentile("V0M");
    fNtupCumInfo->fV0MEq    = cent->GetCentralityPercentile("V0MEq");
    fNtupCumInfo->fV0A      = cent->GetCentralityPercentile("V0A");
    fNtupCumInfo->fV0AEq    = cent->GetCentralityPercentile("V0AEq");
    fNtupCumInfo->fZNA      = cent->GetCentralityPercentile("ZNA");
  }

  AliVZDC *vZDC = InputEvent()->GetZDCData();
  const Double_t *znaTowers = vZDC->GetZNATowerEnergy(); 
  fNtupZdcInfo->fZna0 = znaTowers[0];
  fNtupZdcInfo->fZna1 = znaTowers[1];
  fNtupZdcInfo->fZna2 = znaTowers[2];
  fNtupZdcInfo->fZna3 = znaTowers[3];
  fNtupZdcInfo->fZna4 = znaTowers[4];

  if (fDoCumNtuple && (M>=Mmin)) {
    fNtupCum->Fill();
  }

  fHists[109]->Fill(fNtupCumInfo->fCl1);
  fHists[110]->Fill(fNtupCumInfo->fV0A);
  fHists[111]->Fill(fNtupCumInfo->fZNA);

  Bool_t fillCumHist = kTRUE;
  if (fillCumHist) {
    Int_t  run = InputEvent()->GetRunNumber();
    if (fDoVertexCut) {
      if ((run>=188356&&run<=188366) || (run>=195344&&run<=197388)) {
	if (anau.IsFirstEventInChunk(event))
	  fillCumHist = kFALSE;
	if (!anau.IsVertexSelected2013pA(event))
	  fillCumHist = kFALSE;
      }
    }
    Double_t vz = InputEvent()->GetPrimaryVertex()->GetZ();
    if (TMath::Abs(vz)>10)
      fillCumHist = kFALSE;
  }
  if (fillCumHist) {
    AliVVZERO *vzero = InputEvent()->GetVZEROData();
    Double_t v0a = vzero->GetMTotV0A();
    Double_t v0c = vzero->GetMTotV0C();
    Double_t v0m = vzero->GetMTotV0A()+vzero->GetMTotV0C();
    fHists[112]->Fill(Mall);
    fHists[113]->Fill(M);
    fHists[117]->Fill(v0a,M);
    fHists[118]->Fill(v0c,M);
    fHists[119]->Fill(v0m,M);
    if (M>1) {
      fHists[114]->Fill(M,(Q22-M)/M/(M-1));
      fHists[120]->Fill(M,(Q32-M)/M/(M-1));
      fHists[122]->Fill(M,v2g);
      fHists[123]->Fill(M,v3g);
      fHists[124]->Fill(M,v2g14);
      fHists[125]->Fill(M,v3g14);
      fHists[126]->Fill(M,v2g18);
      fHists[127]->Fill(M,v3g18);
      fHists[130]->Fill(M,Q2r/M);
      fHists[131]->Fill(M,Q2i/M);
      fHists[132]->Fill(M,Q3r/M);
      fHists[133]->Fill(M,Q3i/M);
      fHists[134]->Fill(M,Q4r/M);
      fHists[135]->Fill(M,Q4i/M);
      fHists[136]->Fill(M,Q6r/M);
      fHists[137]->Fill(M,Q6i/M);
    }
    if (M>3) {
      Double_t qc4tmp = (Q22*Q22+Q42-2*Q42re-4*(M-2)*Q22+2*M*(M-3));
      fHists[115]->Fill(M,qc4tmp/M/(M-1)/(M-2)/(M-3));
      qc4tmp = (Q32*Q32+Q62-2*Q32re-4*(M-2)*Q32+2*M*(M-3));
      fHists[121]->Fill(M,qc4tmp/M/(M-1)/(M-2)/(M-3));
    }
    if (M>5) {
      Double_t qc6tmp = Q22*Q22*Q22 + 9*Q42*Q22 - 6*Q6are 
                      + 4*Q6bre - 12*Q6cre 
	              + 18*(M-4)*Q42re + 4*Q62
                      - 9*(M-4)*Q22*Q22 - 9*(M-4)*Q42
                      + 18*(M-2)*(M-5)*Q22
                      - 6*M*(M-4)*(M-5);
      fHists[116]->Fill(M,qc6tmp/M/(M-1)/(M-2)/(M-3)/(M-4)/(M-5));
    }
  }

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
void AliAnalysisTaskCLQA::RunHet(Double_t Etmin)
{
  // Run het analysis.

  if (!fDoHet)
    return;
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
void AliAnalysisTaskCLQA::SetHetParams(Double_t etmin)
{
  // Set parameters for het.

  fHetEtmin  = etmin;
}

//________________________________________________________________________
void AliAnalysisTaskCLQA::UserCreateOutputObjects()
{
  // Create histograms

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  AliAnalysisManager *am  = AliAnalysisManager::GetAnalysisManager();
  if (am) {
    AliPhysicsSelectionTask *t=dynamic_cast<AliPhysicsSelectionTask*>(am->GetTopTasks()->At(0));
    if (t) {
      AliPhysicsSelection *ps=t->GetPhysicsSelection();
      fOutput->Add(ps->GetStatistics(""));
    }
  }

  fHists[0] = new TH1D("fTrigBits",";bit",32,-0.5,31.5);
  fOutput->Add(fHists[0]);
  fHists[1] = new TH1D("fVertexZ",";vertex z (cm)",51,-25.5,25.5);
  fOutput->Add(fHists[1]);
  fHists[2] = new TH1D("fVertexZnc",";vertex z (cm)",51,-25.5,25.5);
  fOutput->Add(fHists[2]);
  fHists[3] = new TH1D("fVertexCuts","",1,0.5,2.5);
  fOutput->Add(fHists[3]);
  fHists[9] = new TProfile("fAccepted","",1,0.5,1.5);
  fOutput->Add(fHists[9]);
  fHists[10] = new TH1D("fV0ACent",";percentile",20,0,100);
  fOutput->Add(fHists[10]);
  fHists[11] = new TH1D("fZNACent",";percentile",20,0,100);
  fOutput->Add(fHists[11]);
  fHists[12] = new TH1D("fV0MCent",";percentile",20,0,100);
  fOutput->Add(fHists[12]);

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
    fHists[25] = new TH2D("fPhiPtTracks_type0",";#phi;p_{T} (GeV/c)",60,0,TMath::TwoPi(),40,0,20);
    fHists[25]->Sumw2();
    fOutput->Add(fHists[25]);
    fHists[26] = new TH2D("fPhiPtTracks_type1",";#phi;p_{T} (GeV/c)",60,0,TMath::TwoPi(),40,0,20);
    fHists[26]->Sumw2();
    fOutput->Add(fHists[26]);
    fHists[27] = new TH2D("fPhiPtTracks_type2",";#phi;p_{T} (GeV/c)",60,0,TMath::TwoPi(),40,0,20);
    fHists[27]->Sumw2();
    fOutput->Add(fHists[27]);
    if (fDoProp) {
      fHists[28] = new TH2D("fDPhiPtTracks",";#Delta#phi;p_{T} (GeV/c)",60,-TMath::Pi(),TMath::Pi(),40,0,20);
      fHists[28]->Sumw2();
      fOutput->Add(fHists[28]);
      fHists[29] = new TH2D("fDPhiPtTracks2",";#Delta#phi;p_{T} (GeV/c)",60,-TMath::Pi(),TMath::Pi(),40,0,20);
      fHists[29]->Sumw2();
      fOutput->Add(fHists[29]);
      fHists[30] = new TH2D("fDPhiPtClusTracks",";#Delta#phi;p_{T} (GeV/c)",60,-TMath::Pi(),TMath::Pi(),40,0,20);
      fHists[30]->Sumw2();
      fOutput->Add(fHists[30]);
      fHists[31] = new TH2D("fDPhiPtReTracks",";#Delta#phi;p_{T} (GeV/c)",60,-TMath::Pi(),TMath::Pi(),40,0,20);
      fHists[31]->Sumw2();
      fOutput->Add(fHists[31]);
    }
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
    fCum = new Cumulants("cmhists",fCumMbins,fCumMmin);
    fCum->SetKine(fCumEtaMin,fCumEtaMax,fCumPtMin,fCumPtMax);
    fCum->EnableEG();
    fCum->EnableQC();
    if (fQC4EG<0)
      fCum->EnableQC4withEG();
    else
      fCum->AddQC4withEG(fQC4EG);
    TList *l=fCum->GetList();
    for (Int_t i=0; i<l->GetEntries(); ++i)
      fOutput->Add(l->At(i));
    Int_t v0bins=1000;
    if (fCumMbins>1000)
      v0bins=25000;
    fHists[117] = new TH2D("fCumV0ACentVsM",";v0a;M",v0bins,0,v0bins,fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[117]);
    fHists[118] = new TH2D("fCumV0CCentVsM",";v0c;M",v0bins,0,v0bins,fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[118]);
    fHists[119] = new TH2D("fCumV0MCentVsM",";v0m;M",v0bins,0,v0bins,fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[119]);
  }
  if (!fCum&&fDoCumulants) {
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
    fHists[106] = new TH1D("fCumCL1MC",";#tracks",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[106]);
    fHists[107] = new TH1D("fCumV0AMC",";#tracks",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[107]);
    fHists[108] = new TH1D("fCumV0CMC",";#tracks",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[108]);
    fHists[109] = new TH1D("fCumCl1Cent",";percentile",10,0,100);
    fOutput->Add(fHists[109]);
    fHists[110] = new TH1D("fCumV0ACent",";percentile",10,0,100);
    fOutput->Add(fHists[110]);
    fHists[111] = new TH1D("fCumZNACent",";percentile",10,0,100);
    fOutput->Add(fHists[111]);
    fHists[112] = new TH1D("fCumMall",";mult",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[112]);
    fHists[113] = new TH1D("fCumM",";mult",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[113]);
    fHists[114] = new TProfile("fCumQC2",";qc2",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[114]);
    fHists[115] = new TProfile("fCumQC4",";qc4",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[115]);
    fHists[116] = new TProfile("fCumQC6",";qc6",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[116]);
    Int_t v0bins=1000;
    if (fCumMbins>1000)
      v0bins=25000;
    fHists[117] = new TH2D("fCumV0ACentVsM",";v0a;M",v0bins,0,v0bins,fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[117]);
    fHists[118] = new TH2D("fCumV0CCentVsM",";v0c;M",v0bins,0,v0bins,fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[118]);
    fHists[119] = new TH2D("fCumV0MCentVsM",";v0m;M",v0bins,0,v0bins,fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[119]);
    fHists[120] = new TProfile("fCum3QC2",";qc2",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[120]);
    fHists[121] = new TProfile("fCum3QC4",";qc4",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[121]);
    fHists[122] = new TProfile("fEtaGapV2",";M",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[122]);
    fHists[123] = new TProfile("fEtaGapV3",";M",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[123]);
    fHists[124] = new TProfile("fEtaGapV214",";M",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[124]);
    fHists[125] = new TProfile("fEtaGapV314",";M",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[125]);
    fHists[126] = new TProfile("fEtaGapV218",";M",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[126]);
    fHists[127] = new TProfile("fEtaGapV318",";M",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[127]);
    fHists[128] = new TH3D("fDPhiDEtaTracks",";#Delta#phi;#Delta#eta;M",64,-0.5*TMath::Pi(),1.5*TMath::Pi(),60,-3,3,fCumMbins/10,0,fCumMbins);
    fHists[128]->Sumw2();
    fOutput->Add(fHists[128]);
    fHists[129] = new TH1D("fDPhiDEtaTrigs",";M",fCumMbins/10,0,fCumMbins);
    fHists[129]->Sumw2();
    fOutput->Add(fHists[129]);
    fHists[130] = new TProfile("fCumQ2r",";q2r",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[130]);
    fHists[131] = new TProfile("fCumQ2i",";q2i",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[131]);
    fHists[132] = new TProfile("fCumQ3r",";q3r",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[132]);
    fHists[133] = new TProfile("fCumQ3i",";q3i",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[133]);
    fHists[134] = new TProfile("fCumQ4r",";q4r",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[134]);
    fHists[135] = new TProfile("fCumQ4i",";q4i",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[135]);
    fHists[136] = new TProfile("fCumQ6r",";q6r",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[136]);
    fHists[137] = new TProfile("fCumQ6i",";q6i",fCumMbins,0,fCumMbins);
    fOutput->Add(fHists[137]);

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
    fNtupZdcInfo = new AliNtupZdcInfo;
    fNtupCum->Branch("zdc", &fNtupZdcInfo, 32*1024, 99);
    if (fDoCumNtuple)
      fOutput->Add(fNtupCum);
  }

  if (fDoHet) {
    fNtupHet = new TTree("NtupHet", "Ntuple for het analysis");
    if (1) {
      fNtupHet->SetDirectory(0);
    } else {
      TFile *f = OpenFile(1); 
      if (f) {
	f->SetCompressionLevel(2);
	fNtupHet->SetDirectory(f);
	fNtupHet->SetAutoFlush(-4*1024*1024);
	fNtupHet->SetAutoSave(0);
      }
    }
    fNtupHetInfo = new AliNtupHetInfo;
    fNtupHet->Branch("het", &fNtupHetInfo, 32*1024, 99);
    fOutput->Add(fNtupHet);
  }

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}
