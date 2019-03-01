#include "TString.h"
#include "TMath.h"
#include "AliForwardFlowUtil.h"
#include "TFile.h"

#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TList.h>
#include <THn.h>

#include "AliLog.h"
#include "AliForwardFlowRun2Task.h"
#include "AliForwardQCumulantRun2.h"
#include "AliForwardGenericFramework.h"

#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"

#include "AliAODMCParticle.h"

#include "AliForwardFlowUtil.h"

#include "AliVVZERO.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"

#include "AliESDEvent.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"

#include "AliAnalysisFilter.h"
#include "AliMultSelection.h"
#include "AliMultiplicity.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
//________________________________________________________________________
AliForwardFlowUtil::AliForwardFlowUtil():
fevent(),
fAODevent(),
fMCevent(),
mc(kFALSE),
dNdeta(),
fSettings(),
minpt(0.2),
maxpt(5),
dodNdeta(kTRUE),
fTrackDensity(),
fState(),
fMaxConsequtiveStrips(3),
fLowCutvalue(0)
{
}

void AliForwardFlowUtil::FillData(TH2D*& refDist, TH2D*& centralDist, TH2D*& forwardDist){
    if (!fSettings.mc) {
      AliAODEvent* aodevent = dynamic_cast<AliAODEvent*>(fevent);
      this->fAODevent = aodevent;
      AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(aodevent->FindListObject("Forward"));
      forwardDist = &aodfmult->GetHistogram();
      for (Int_t etaBin = 1; etaBin <= forwardDist->GetNbinsX(); etaBin++) {
        for (Int_t phiBin = 1; phiBin <= forwardDist->GetNbinsX(); phiBin++) {
          if (dodNdeta) dNdeta->Fill(forwardDist->GetXaxis()->GetBinCenter(etaBin),forwardDist->GetBinContent(etaBin, phiBin));
        }
      }

        if (fSettings.ref_mode && fSettings.kSPDref) this->FillFromTracklets(refDist);
        else if (fSettings.ref_mode && fSettings.kITSref) this->FillFromCentralClusters(refDist);
        else if (fSettings.ref_mode && fSettings.kFMDref) {
          AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(aodevent->FindListObject("Forward"));
          refDist = &aodfmult->GetHistogram();
        }
        else {
          this->minpt = 0.2;
          this->maxpt = 5;
          this->FillFromTracks(refDist, fSettings.tracktype);
          this->minpt = fSettings.minpt;
          this->maxpt = fSettings.maxpt;
        }

      if (fSettings.useSPD) this->FillFromTracklets(centralDist);
      else if (fSettings.useITS) this->FillFromCentralClusters(refDist);
      else                  this->FillFromTracks(centralDist, fSettings.tracktype);
    }
    else {

      this->mc = kTRUE;
      Float_t zvertex = this->GetZ();

      if(!fMCevent)
        throw std::runtime_error("Not MC as expected");

      if (fSettings.esd){
        if (fSettings.use_primaries_cen && fSettings.use_primaries_fwd){
          this->FillFromPrimaries(centralDist, forwardDist);
        }
        else if (!fSettings.use_primaries_cen && !fSettings.use_primaries_fwd){
          this->FillFromTrackrefs(centralDist, forwardDist);
        }
        else if (fSettings.use_primaries_cen && !fSettings.use_primaries_fwd){
          this->FillFromPrimaries(centralDist);
          this->FillFromTrackrefs(forwardDist);
        }
        else if (!fSettings.use_primaries_cen && fSettings.use_primaries_fwd){
          this->FillFromTrackrefs(centralDist);
          this->FillFromPrimaries(forwardDist);
        }
      }
      else{ // AOD

        if (fSettings.use_primaries_cen && fSettings.use_primaries_fwd){ //prim central and forward
          if (fSettings.maxpt < 5.0) {
            this->minpt = 0.2;
            this->maxpt = 5.0;
            this->FillFromPrimariesAOD(refDist);
            this->minpt = fSettings.minpt;
            this->maxpt = fSettings.maxpt;
          }

          this->FillFromPrimariesAOD(centralDist, forwardDist);
        }
        else if (fSettings.use_primaries_cen && !fSettings.use_primaries_fwd){ //prim central, AOD forward
          if (fSettings.maxpt < 5.0) {
            this->minpt = 0.2;
            this->maxpt = 5.0;
            this->FillFromPrimariesAOD(refDist);
            this->minpt = fSettings.minpt;
            this->maxpt = fSettings.maxpt;
          }

          this->FillFromPrimariesAOD(centralDist);
          AliAODEvent* aodevent = dynamic_cast<AliAODEvent*>(fevent);

          AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(aodevent->FindListObject("Forward"));
          forwardDist = &aodfmult->GetHistogram();

          for (Int_t etaBin = 1; etaBin <= forwardDist->GetNbinsX(); etaBin++) {
            for (Int_t phiBin = 1; phiBin <= forwardDist->GetNbinsX(); phiBin++) {
              if (dodNdeta) dNdeta->Fill(forwardDist->GetXaxis()->GetBinCenter(etaBin),forwardDist->GetBinContent(etaBin, phiBin));
            }
          }
        }
      }
    }
}


Bool_t AliForwardFlowUtil::ExtraEventCutFMD(TH2D& forwarddNdedp, double cent, Bool_t mc, TH2D* hOutliers){
  Bool_t useEvent = true;
  //if (useEvent) return useEvent;
  Int_t nBadBins = 0;
  Int_t phibins = forwarddNdedp.GetNbinsY();
  Double_t totalFMDpar = 0;

  for (Int_t etaBin = 1; etaBin <= forwarddNdedp.GetNbinsX(); etaBin++) {

    Double_t eta = forwarddNdedp.GetXaxis()->GetBinCenter(etaBin);
    Double_t runAvg = 0;
    Double_t avgSqr = 0;
    Double_t max = 0;
    Int_t nInAvg = 0;

    for (Int_t phiBin = 0; phiBin <= phibins; phiBin++) {
      // if (!mc){
      //   if ( fabs(eta) > 1.7) {
      //     if (phiBin == 0 && forwarddNdedp.GetBinContent(etaBin, 0) == 0) break;
      //   }
      // }
      Double_t weight = forwarddNdedp.GetBinContent(etaBin, phiBin);
      if (!weight){
        weight = 0;
      }
      totalFMDpar += weight;

      // We calculate the average Nch per. bin
      avgSqr += weight*weight;
      runAvg += weight;
      nInAvg++;
      if (weight == 0) continue;
      if (weight > max) {
        max = weight;
      }
    } // End of phi loop

    // Outlier cut calculations
    double fSigmaCut = 4.0;
    if (nInAvg > 0) {
      runAvg /= nInAvg;
      avgSqr /= nInAvg;
      Double_t stdev = (nInAvg > 1 ? TMath::Sqrt(nInAvg/(nInAvg-1))*TMath::Sqrt(avgSqr - runAvg*runAvg) : 0);
      Double_t nSigma = (stdev == 0 ? 0 : (max-runAvg)/stdev);
      hOutliers->Fill(cent,nSigma);
      //std::cout << "sigma = " << nSigma << std::endl;
      if (fSigmaCut > 0. && nSigma >= fSigmaCut && cent < 60) nBadBins++;
      else nBadBins = 0;
      // We still finish the loop, for fOutliers to make sense,
      // but we do no keep the event for analysis
      if (nBadBins > 3) useEvent = false;
     //if (nBadBins > 3) std::cout << "NUMBER OF BAD BINS > 3" << std::endl;
    }
  } // End of eta bin
  if (totalFMDpar < 10) useEvent = false;

  return useEvent;
}


Double_t AliForwardFlowUtil::GetZ(){
  if (this->fSettings.mc) return fMCevent->GetPrimaryVertex()->GetZ();
  else return fevent->GetPrimaryVertex()->GetZ();
}


Double_t AliForwardFlowUtil::GetCentrality(TString centrality_estimator){
  // Get MultSelection
  AliMultSelection *MultSelection;

  MultSelection  = (AliMultSelection*)fevent->FindListObject("MultSelection");

  return MultSelection->GetMultiplicityPercentile(centrality_estimator);
}


void AliForwardFlowUtil::FillFromTrackrefs(TH2D*& cen, TH2D*& fwd) const
{

  Int_t nTracks = fMCevent->Stack()->GetNtrack();

  for (Int_t iTr = 0; iTr < nTracks; iTr++) {
      AliMCParticle* p = static_cast< AliMCParticle* >(fMCevent->GetTrack(iTr));
    if (p->Charge() == 0) continue;

    for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) {
      AliTrackReference* ref = p->GetTrackReference(iTrRef);
      // Check hit on FMD
      if (!ref) continue;
      if (AliTrackReference::kTPC == ref->DetectorId()){
        Double_t x      = ref->X() - fMCevent->GetPrimaryVertex()->GetX();
        Double_t y      = ref->Y() - fMCevent->GetPrimaryVertex()->GetY();
        Double_t z      = ref->Z() - fMCevent->GetPrimaryVertex()->GetZ();
        Double_t rr     = TMath::Sqrt(x * x + y * y);
        Double_t thetaR = TMath::ATan2(rr, z);
        Double_t phiR   = TMath::ATan2(y,x);

        if (phiR < 0) {
          phiR += 2*TMath::Pi();
        }
        if (thetaR < 0) {
          thetaR += 2*TMath::Pi();
        }
        cen->Fill(-TMath::Log(TMath::Tan(thetaR / 2)),phiR);
        if (dodNdeta) dNdeta->Fill(-TMath::Log(TMath::Tan(thetaR / 2)),1);

      }
      else if (AliTrackReference::kFMD == ref->DetectorId()) {
        Double_t x      = ref->X() - fMCevent->GetPrimaryVertex()->GetX();
        Double_t y      = ref->Y() - fMCevent->GetPrimaryVertex()->GetY();
        Double_t z      = ref->Z() - fMCevent->GetPrimaryVertex()->GetZ();
        Double_t rr     = TMath::Sqrt(x * x + y * y);
        Double_t thetaR = TMath::ATan2(rr, z);
        Double_t phiR   = TMath::ATan2(y,x);

        if (phiR < 0) {
          phiR += 2*TMath::Pi();
        }
        if (thetaR < 0) {
          thetaR += 2*TMath::Pi();
        }
        fwd->Fill(-TMath::Log(TMath::Tan(thetaR / 2)),phiR);
        if (dodNdeta) dNdeta->Fill(-TMath::Log(TMath::Tan(thetaR / 2)),1);

      }
    }
  }
}



void AliForwardFlowUtil::FillFromTrackrefs(TH2D*& fwd)
{
  Int_t nTracks   = fMCevent->GetNumberOfTracks();// stack->GetNtrack();

  Int_t nPrim     = fMCevent->GetNumberOfPrimaries();//fAOD->GetNumberOfPrimaries();
  for (Int_t iTr = 0; iTr < nTracks; iTr++) {
    AliMCParticle* particle =
      static_cast<AliMCParticle*>(fMCevent->GetTrack(iTr));

    // Check if this charged and a primary
    if (particle->Charge() == 0) continue;

    Bool_t isPrimary = fMCevent->Stack()->IsPhysicalPrimary(iTr) && iTr < nPrim;

    AliMCParticle* mother = isPrimary ? particle : GetMother(iTr,fMCevent);
    if (!mother) mother = particle;
    // IF the track corresponds to a primary, pass that as both
    // arguments.
    ProcessTrack(particle, mother, fwd);
  } // Loop over tracks
}


AliMCParticle*
AliForwardFlowUtil::GetMother(Int_t iTr, const AliMCEvent* event) const
{
  //
  // Track down primary mother
  //
  Int_t                i         = iTr;
  Bool_t               gammaSeen = false;
  AliMCParticle* candidate = 0;
  do {
    AliMCParticle* p = static_cast<AliMCParticle*>(event->GetTrack(i));
    if (!p) break;
    if (gammaSeen && TMath::Abs(p->PdgCode()) == 111)
      // If we're looking for a mother pi0 of gamma, and we find it
      // here, we return it - irrespective of whether it's flagged as
      // a primary or not.
      return p;

    if (event->IsPhysicalPrimary(i)) {
      candidate = p;
      if (fTrackGammaToPi0 && TMath::Abs(p->PdgCode()) == 22)
	// If we want to track gammas back to a possible pi0, we flag
	// the gamma seen, and store it as a candidate in case we do
	// not find a pi0 in the stack
	gammaSeen = true;
      else
	break;
    }

    // We get here if the current track isn't a primary, or it was a
    // primary gamma and we want to track back to a pi0.
    i = p->GetMother();
  } while (i > 0);

  // Return our candidate (gamma) if we find no mother pi0.  Note, we
  // should never get here with a null pointer, so we issue a warning
  // in that case.
  if (!candidate)
    AliWarningF("No mother found for track # %d", iTr);
  return candidate;
}

Bool_t
AliForwardFlowUtil::ProcessTrack(AliMCParticle* particle,
				    AliMCParticle* mother,TH2D*& fwd)
{
  // Check the returned particle
  //
  // Note: If particle refers to a primary, then particle and mother
  // refers to the same particle (the address are the same)
  //
  if (!particle) return false;

  Int_t              nTrRef = particle->GetNumberOfTrackReferences();
  AliTrackReference* store  = 0;

  BeginTrackRefs();

  // Double_t oTheta= 0;
  for (Int_t iTrRef = 0; iTrRef < nTrRef; iTrRef++) {
    AliTrackReference* ref = particle->GetTrackReference(iTrRef);

    // Check existence
    if (!ref) continue;

    // Check that we hit an Base element
    if (ref->DetectorId() != AliTrackReference::kFMD) continue;

    AliTrackReference* test = ProcessRef(particle, mother, ref,fwd);
    if (test) store = test;

  } // Loop over track references
  if (!store) return true; // Nothing found

  StoreParticle(particle, mother, store,fwd);
  EndTrackRefs();

  return true;
}



AliTrackReference*
AliForwardFlowUtil::ProcessRef(AliMCParticle*       particle,
				  AliMCParticle* mother,
				 AliTrackReference*   ref,TH2D*& fwd)
{
  // Process track references of a track
  //
  // Note: If particle refers to a primary, then particle and mother
  // refers to the same particle (the address are the same)
  //

  // Get the detector coordinates
  UShort_t d, s, t;
  Char_t r;
  AliFMDStripIndex::Unpack(ref->UserId(), d, r, s, t);
  Double_t edep, length, dEdep, dLength;
  AliFMDEncodedEdx::Decode((ref->UserId() >> 19), edep, length, dEdep, dLength);

  Double_t normaldEdx=0.0;
   if(length>0.0)
	normaldEdx=(edep/length)/4.406; // 4.406 mip in Si per 1 cm

  // Calculate distance of previous reference to base of cluster
  UShort_t nT = TMath::Abs(t - fState.startStrip) + 1;

  // Now check if we should flush to output
  Bool_t used = false;

  // If this is a new detector/ring, then reset the other one
  // Check if we have a valid old detectorm ring, and sector
  if (fState.oldDetector >  0 &&
      fState.oldRing     != '\0' &&
      fState.oldSector   != 1024) {
    // New detector, new ring, or new sector
    if (d != fState.oldDetector   ||
	r != fState.oldRing       ||
	s != fState.oldSector) {
      used = true;
    }
    else if (nT > fMaxConsequtiveStrips) {

      used = true;
    }
  }
  if (used) {

    // Int_t nnT   = TMath::Abs(fState.oldStrip - fState.startStrip) + 1;
    StoreParticle(particle, mother, fState.longest,fwd);
    fState.Clear(false);
  }

  if(normaldEdx<fLowCutvalue)
	return 0x0;
  // If base of cluster not set, set it here.
  if (fState.startStrip == 1024) fState.startStrip = t;

  // Calculate distance of previous reference to base of cluster
  fState.nStrips = TMath::Abs(t - fState.startStrip) + 1;

  // Count number of track refs in this sector
  fState.nRefs++;

  fState.oldDetector = d;
  fState.oldRing     = r;
  fState.oldSector   = s;
  fState.oldStrip    = t;



  // The longest passage is determined through the angle
  Double_t ang  = GetTrackRefTheta(ref);
  if (ang > fState.angle) {
    fState.longest = ref;
    fState.angle   = ang;
  }
  return fState.longest;
}


Double_t
AliForwardFlowUtil::GetTrackRefTheta(const AliTrackReference* ref) const
{
  // Get the incidient angle of the track reference.
  const AliVVertex* vertex = this->fMCevent->GetPrimaryVertex();
  // Calculate the vector pointing from the vertex to the track reference on the detector
  Double_t x      = ref->X() - vertex->GetX();
  Double_t y      = ref->Y() - vertex->GetY();
  Double_t z      = ref->Z() - vertex->GetZ();
  Double_t rr   = TMath::Sqrt(x*x+y*y);
  Double_t theta= TMath::ATan2(rr,z);
  Double_t ang  = TMath::Abs(TMath::Pi()-theta);
  return ang;
}

void
AliForwardFlowUtil::StoreParticle(AliMCParticle*       particle,
				     AliMCParticle* mother,
				     AliTrackReference*   ref,TH2D*& fwd)
{
  Double_t eta_mother = mother->Eta();
  Double_t phi_mother = (mother->Phi());//Wrap02pi

  Double_t *etaPhi = new Double_t[2];
  this->GetTrackRefEtaPhi(particle, etaPhi);

  Double_t phi_tr = etaPhi[1]; //Wrap02pi
  Double_t eta_tr = etaPhi[0];

  fwd->Fill(eta_tr,phi_tr,1);
  if (dodNdeta) dNdeta->Fill(eta_tr,1);

  return;
}

void
AliForwardFlowUtil::BeginTrackRefs()
{
  fState.Clear(true);
}


void
AliForwardFlowUtil::EndTrackRefs()
{
  fState.Clear(true);
}


void
AliForwardFlowUtil::State::Clear(Bool_t alsoCount)
{
  angle       = 0;
  oldDetector = 0;
  oldRing     = '\0';
  oldSector   = 1024;
  oldStrip    = 1024;
  startStrip  = 1024;
  nRefs       = 0;
  nStrips     = 0;
  longest     = 0x0;
  if (alsoCount) count = 0;
}

void AliForwardFlowUtil::FillFromPrimariesAOD(TH2D*& cen) const
{
    Int_t nTracksMC = fMCevent->GetNumberOfTracks();
  for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
    AliAODMCParticle* p = static_cast< AliAODMCParticle* >(fMCevent->GetTrack(iTr));
    if (!p->IsPhysicalPrimary()) continue;
    if (p->Charge() == 0) continue;

    Double_t eta = p->Eta();
    if (TMath::Abs(eta) < 1.1) {
      if (p->Pt()>=this->minpt && p->Pt()<=this->maxpt){
        cen->Fill(eta,p->Phi(),1);
        if (dodNdeta) dNdeta->Fill(eta,1);
      }
    }
  }
}



void AliForwardFlowUtil::GetTrackRefEtaPhi(AliMCParticle* p, Double_t* etaPhi) {
  AliTrackReference* ref = 0x0;
  for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) {
    ref = p->GetTrackReference(iTrRef);
    // Check hit on FMD
    if (ref && AliTrackReference::kFMD == ref->DetectorId()) {
      break;
    }
    else {
      ref = 0x0;
    }
  }
  if (!ref) {
    etaPhi = 0x0;
    return;
  }
  const AliVVertex* vertex = fMCevent->GetPrimaryVertex();
  // Calculate the vector pointing from the vertex to the track reference on the detector
  Double_t x      = ref->X() - vertex->GetX();
  Double_t y      = ref->Y() - vertex->GetY();
  Double_t z      = ref->Z() - vertex->GetZ();
  Double_t rr     = TMath::Sqrt(x * x + y * y);
  Double_t thetaR = TMath::ATan2(rr, z);
  Double_t phiR   = TMath::ATan2(y,x);
  // Correct angles
  if (thetaR < 0) {
    thetaR += 2*TMath::Pi();
  }
  if (phiR < 0) {
    phiR += 2*TMath::Pi();
  }
  etaPhi[0] = -TMath::Log(TMath::Tan(thetaR / 2));
  etaPhi[1] = phiR;
  // cout << x << " " << y << " " << z << endl << endl;
  // cout << etaPhi[0] - p->Eta() << " " << etaPhi[1] - p->Phi() << " "
  //      << this->GetDaughters(p).size() << " "
  //      << p->PdgCode() << " "
  //      << endl;
}

void AliForwardFlowUtil::FillFromPrimaries(TH2D*& cen) const
{
  Int_t nTracksMC = fMCevent->GetNumberOfTracks();

  for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
    AliMCParticle* p = static_cast< AliMCParticle* >(fMCevent->GetTrack(iTr));
    if (!p->IsPhysicalPrimary()) continue;
    if (p->Charge() == 0) continue;

    Double_t eta = p->Eta();
    if (TMath::Abs(eta) < 1.1) {
      if (p->Pt()>=this->minpt && p->Pt()<=this->maxpt){
        cen->Fill(eta,p->Phi(),1);
        if (dodNdeta) dNdeta->Fill(eta,1);
      }
    }
  }
}

void AliForwardFlowUtil::FillFromPrimaries(TH2D*& cen, TH2D*& fwd) const
{
  Int_t nTracksMC = fMCevent->GetNumberOfTracks();

  for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
    AliMCParticle* p = static_cast< AliMCParticle* >(fMCevent->GetTrack(iTr));
    if (!p->IsPhysicalPrimary()) continue;
    if (p->Charge() == 0) continue;

    Double_t eta = p->Eta();
    if (TMath::Abs(eta) < 1.1) {
      if (p->Pt()>=this->minpt && p->Pt()<=this->maxpt){
        cen->Fill(eta,p->Phi(),1);
        if (dodNdeta) dNdeta->Fill(eta,1);
      }
    }
    if (eta < 5 /*fwd->GetXaxis()-GetXmax()*/ && eta > -3.5 /*fwd->GetXaxis()-GetXmin()*/) {
      if (TMath::Abs(eta) >= 1.7){
        fwd->Fill(eta,p->Phi(),1);
        if (dodNdeta) dNdeta->Fill(eta,1);
      }
    }
  }
}


void AliForwardFlowUtil::FillFromPrimariesAOD(TH2D*& cen, TH2D*& fwd) const
{
  Int_t nTracksMC = fMCevent->GetNumberOfTracks();

  for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
    AliAODMCParticle* p = static_cast< AliAODMCParticle* >(fMCevent->GetTrack(iTr));
    if (!p->IsPhysicalPrimary()) continue;
    if (p->Charge() == 0) continue;

    Double_t eta = p->Eta();
    if (TMath::Abs(eta) < 1.1) {
      if (p->Pt()>=this->minpt && p->Pt()<=this->maxpt){
        cen->Fill(eta,p->Phi(),1);
        if (dodNdeta) dNdeta->Fill(eta,1);
      }
    }
    if (eta < 5 /*fwd->GetXaxis()-GetXmax()*/ && eta > -3.5 /*fwd->GetXaxis()-GetXmin()*/) {
      if (TMath::Abs(eta) >= 1.7){
        fwd->Fill(eta,p->Phi(),1);
        if (dodNdeta) dNdeta->Fill(eta,1);
      }
    }
  }
}


void AliForwardFlowUtil::FillFromTracklets(TH2D*& cen) const {
  AliAODTracklets* aodTracklets = fAODevent->GetTracklets();

  for (Int_t i = 0; i < aodTracklets->GetNumberOfTracklets(); i++) {
    cen->Fill(aodTracklets->GetEta(i),aodTracklets->GetPhi(i), 1);
    if (dodNdeta) dNdeta->Fill(aodTracklets->GetEta(i),1);
  }
}


void AliForwardFlowUtil::FillFromCentralClusters(TH2D*& cen) const {
  AliAODCentralMult* aodcmult = static_cast<AliAODCentralMult*>(fAODevent->FindListObject("CentralClusters"));
  cen = &aodcmult->GetHistogram();
}

void AliForwardFlowUtil::FillFromTracks(TH2D*& cen, UInt_t tracktype) const {
  Int_t  iTracks(fevent->GetNumberOfTracks());
  std::cout << this->maxpt << std::endl;
  for(Int_t i(0); i < iTracks; i++) {

  // loop  over  all  the  tracks
    AliAODTrack* track = static_cast<AliAODTrack *>(fAODevent->GetTrack(i));

    if (track->TestFilterBit(tracktype) && track->GetTPCNcls() > fSettings.fnoClusters){

      if( fSettings.fCutChargedDCAzMax > 0. || fSettings.fCutChargedDCAxyMax > 0.){
        Double_t dTrackXYZ[3] = {0.,0.,0.};
        Double_t dVertexXYZ[3] = {0.,0.,0.};
        Double_t dDCAXYZ[3] = {0.,0.,0.};
        const AliAODVertex* vertex = fAODevent->GetPrimaryVertex();

        track->GetXYZ(dTrackXYZ);
        vertex->GetXYZ(dVertexXYZ);

        for(Short_t i(0); i < 3; i++) { dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i]; }

        if(fSettings.fCutChargedDCAzMax > 0. && TMath::Abs(dDCAXYZ[2]) > fSettings.fCutChargedDCAzMax) continue;

        if(fSettings.fCutChargedDCAxyMax > 0. && TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]) > fSettings.fCutChargedDCAxyMax) continue;
      }

      if (track->Pt() >= this->minpt && track->Pt() <= this->maxpt){
        cen->Fill(track->Eta(),track->Phi(), 1);
        if (dodNdeta) dNdeta->Fill(track->Eta(),1);
      }
    }
  }
}



void AliForwardFlowUtil::FillFromTrackrefs(TH3D*& cen, TH3D*& fwd, Double_t zvertex) const
{
  Int_t nTracks = fMCevent->Stack()->GetNtrack();

  for (Int_t iTr = 0; iTr < nTracks; iTr++) {
      AliMCParticle* p = static_cast< AliMCParticle* >(fMCevent->GetTrack(iTr));
    if (p->Charge() == 0) continue;

    for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) {
      AliTrackReference* ref = p->GetTrackReference(iTrRef);
      // Check hit on FMD
      if (!ref) continue;
      if (AliTrackReference::kTPC == ref->DetectorId()){
        Double_t x      = ref->X() - fMCevent->GetPrimaryVertex()->GetX();
        Double_t y      = ref->Y() - fMCevent->GetPrimaryVertex()->GetY();
        Double_t z      = ref->Z() - fMCevent->GetPrimaryVertex()->GetZ();
        Double_t rr     = TMath::Sqrt(x * x + y * y);
        Double_t thetaR = TMath::ATan2(rr, z);
        Double_t phiR   = TMath::ATan2(y,x);

        if (phiR < 0) {
          phiR += 2*TMath::Pi();
        }
        if (thetaR < 0) {
          thetaR += 2*TMath::Pi();
        }
        cen->Fill(-TMath::Log(TMath::Tan(thetaR / 2)),phiR, zvertex);
        if (dodNdeta) dNdeta->Fill(-TMath::Log(TMath::Tan(thetaR / 2)),1);

      }
      else if (AliTrackReference::kFMD == ref->DetectorId()) {
        Double_t x      = ref->X() - fMCevent->GetPrimaryVertex()->GetX();
        Double_t y      = ref->Y() - fMCevent->GetPrimaryVertex()->GetY();
        Double_t z      = ref->Z() - fMCevent->GetPrimaryVertex()->GetZ();
        Double_t rr     = TMath::Sqrt(x * x + y * y);
        Double_t thetaR = TMath::ATan2(rr, z);
        Double_t phiR   = TMath::ATan2(y,x);

        if (phiR < 0) {
          phiR += 2*TMath::Pi();
        }
        if (thetaR < 0) {
          thetaR += 2*TMath::Pi();
        }
        fwd->Fill(-TMath::Log(TMath::Tan(thetaR / 2)),phiR, zvertex);
        if (dodNdeta) dNdeta->Fill(-TMath::Log(TMath::Tan(thetaR / 2)),1);
      }
    }
  }
}


void AliForwardFlowUtil::FillFromPrimaries(TH3D*& cen, TH3D*& fwd, Double_t zvertex) const
{
  Int_t nTracksMC = fMCevent->GetNumberOfTracks();

  for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
    AliMCParticle* p = static_cast< AliMCParticle* >(fMCevent->GetTrack(iTr));
    if (!p->IsPhysicalPrimary()) continue;
    if (p->Charge() == 0) continue;

    Double_t eta = p->Eta();
    if (TMath::Abs(eta) < 1.1) {
      if (p->Pt()>=this->minpt && p->Pt()<=this->maxpt){
        cen->Fill(eta,p->Phi(),zvertex,1);
        if (dodNdeta) dNdeta->Fill(eta,1);
      }
    }
    if (eta < 5 /*fwd->GetXaxis()-GetXmax()*/ && eta > -3.5 /*fwd->GetXaxis()-GetXmin()*/) {
      if (TMath::Abs(eta) >= 1.7){
        fwd->Fill(eta,p->Phi(),zvertex,1);
        if (dodNdeta) dNdeta->Fill(eta,1);
      }
    }
  }
}


void AliForwardFlowUtil::FillFromTracklets(TH3D*& cen, Double_t zvertex) const {
  AliAODTracklets* aodTracklets = fAODevent->GetTracklets();

  for (Int_t i = 0; i < aodTracklets->GetNumberOfTracklets(); i++) {
    cen->Fill(aodTracklets->GetEta(i),aodTracklets->GetPhi(i),zvertex, 1);
    if (dodNdeta) dNdeta->Fill(aodTracklets->GetEta(i),1);
  }
}


void AliForwardFlowUtil::FillFromTracks(TH3D*& cen, Int_t tracktype, Double_t zvertex) const {
  Int_t  iTracks(fevent->GetNumberOfTracks());
  for(Int_t i(0); i < iTracks; i++) {

  // loop  over  all  the  tracks
    AliAODTrack* track = static_cast<AliAODTrack *>(fAODevent->GetTrack(i));
    if (track->TestFilterBit(tracktype)){
      if (track->Pt() >= this->minpt && track->Pt() <= this->maxpt){
        cen->Fill(track->Eta(),track->Phi(), zvertex, 1);
        if (dodNdeta) dNdeta->Fill(track->Eta(),1);
      }
    }
  }
}


AliTrackReference* AliForwardFlowUtil::IsHitFMD(AliMCParticle* p) {
  //std::cout << "p->GetNumberOfTrackReferences() = " << p->GetNumberOfTrackReferences() << std::endl;
  for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) {
    AliTrackReference* ref = p->GetTrackReference(iTrRef);
    // Check hit on FMD
    if (!ref || AliTrackReference::kFMD != ref->DetectorId()) {
      continue;
    }
    else {
      return ref;
    }
  }
  return 0x0;
}

AliTrackReference* AliForwardFlowUtil::IsHitTPC(AliMCParticle* p) {
  for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) {
    AliTrackReference* ref = p->GetTrackReference(iTrRef);
    // Check hit on FMD
    if (!ref || AliTrackReference::kTPC != ref->DetectorId()) {
      continue;
    }
    else {
      return ref;
    }
  }
  return 0x0;
}


void AliForwardFlowUtil::MakeFakeHoles(TH2D& forwarddNdedp){
  for (Int_t etaBin = 125; etaBin <= 137; etaBin++){
    forwarddNdedp.SetBinContent(etaBin,17, 0.0);
    forwarddNdedp.SetBinContent(etaBin,18, 0.0);
  }
  for (Int_t etaBin = 168; etaBin <= 185; etaBin++){
    forwarddNdedp.SetBinContent(etaBin,14, 0.0);
  }
}
