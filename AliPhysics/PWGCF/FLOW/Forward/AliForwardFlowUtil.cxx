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
dNdeta(),
fSettings(),
minpt(0.2),
maxpt(5),
dodNdeta(kTRUE),
fTrackDensity(),
fState(),
fLowCutvalue(0),
fTrackGammaToPi0(false),
fStored(0)
{
}

void AliForwardFlowUtil::FillData(TH2D*& refDist, TH2D*& centralDist, TH2D*& forwardDist)
{
  if (!fSettings.mc) {
    AliAODEvent* aodevent = dynamic_cast<AliAODEvent*>(fevent);
    this->fAODevent = aodevent;
      
    // Fill forwardDist
    FillFromForwardClusters(forwardDist);
    // Fill centralDist
    if (fSettings.useSPD) this->FillFromTracklets(centralDist);
    else if (fSettings.useITS) this->FillFromCentralClusters(centralDist);
    else  this->FillFromTracks(centralDist, fSettings.tracktype); //(fSettings.useTPC)

    // Fill refDist
    if (fSettings.ref_mode & fSettings.kSPDref)      this->FillFromTracklets(refDist);
    else if (fSettings.ref_mode & fSettings.kITSref) this->FillFromCentralClusters(refDist);
    else if (fSettings.ref_mode & fSettings.kFMDref) this->FillFromForwardClusters(refDist); 
    else {
      this->minpt = 0.2;
      this->maxpt = 5;
      this->FillFromTracks(refDist, fSettings.tracktype);
      this->minpt = fSettings.minpt;
      this->maxpt = fSettings.maxpt;
    }
  }
  else {
    if(!fMCevent)
      throw std::runtime_error("Not MC as expected");

    if (fSettings.esd){
      // Fill forwardDist
      if (fSettings.use_primaries_fwd) this->FillFromPrimariesFMD(forwardDist);
      else this->FillFromTrackrefsFMD(forwardDist);
      
      // Fill centralDist
      if (fSettings.useITS){
        if (!fSettings.use_primaries_cen) this->FillFromTrackrefsITS(centralDist);
        else this->FillFromPrimariesITS(centralDist);
      }
      else if (fSettings.useSPD && fSettings.use_primaries_cen) this->FillFromPrimariesSPD(centralDist);
      else if (fSettings.useTPC && fSettings.use_primaries_cen) this->FillFromPrimariesTPC(centralDist);
      else std::cout << "No valid central detector chosen." << std::endl;

      // Fill refDist
      if (fSettings.ref_mode & fSettings.kITSref){
        if (!fSettings.use_primaries_cen) this->FillFromTrackrefsITS(refDist);
        else this->FillFromPrimariesITS(refDist);
      }
      else if (fSettings.ref_mode & fSettings.kSPDref && fSettings.use_primaries_cen) this->FillFromPrimariesSPD(centralDist);
      else if (fSettings.ref_mode & fSettings.kTPCref && fSettings.use_primaries_cen) {
        this->minpt = 0.2;
        this->maxpt = 5;
        this->FillFromPrimariesTPC(refDist);
        this->minpt = fSettings.minpt;
        this->maxpt = fSettings.maxpt;
      }
      else std::cout << "No valid reference detector chosen." << std::endl;
    }
    else { // AOD

      // Fill forwardDist
      if (fSettings.use_primaries_fwd) this->FillFromPrimariesAODFMD(forwardDist);
      else this->FillFromForwardClusters(forwardDist);

      // Fill centralDist
      if (fSettings.useITS) {
        if (!fSettings.use_primaries_cen) this->FillFromCentralClusters(centralDist);
        else this->FillFromPrimariesAODITS(centralDist);
      }
      else if (fSettings.useSPD) {
        if (fSettings.use_primaries_cen) this->FillFromPrimariesAODSPD(centralDist); 
        else this->FillFromTracklets(centralDist);
      }
      else if (fSettings.useTPC){
        if (fSettings.use_primaries_cen) this->FillFromPrimariesAODTPC(centralDist); 
        else this->FillFromTracks(centralDist, fSettings.tracktype);
      }
      else std::cout << "No valid central detector chosen." << std::endl;

      // Fill refDist
      if (fSettings.ref_mode & fSettings.kFMDref) {
        if (!fSettings.use_primaries_fwd) this->FillFromForwardClusters(refDist);
        else this->FillFromPrimariesAODFMD(refDist);
      }
      else if (fSettings.ref_mode & fSettings.kITSref) {
        if (!fSettings.use_primaries_cen) this->FillFromCentralClusters(refDist);
        else this->FillFromPrimariesAODITS(refDist);
      }
      else if (fSettings.ref_mode & fSettings.kSPDref) {
        if (fSettings.use_primaries_cen) this->FillFromPrimariesAODSPD(refDist);
        else this->FillFromTracklets(refDist);
      }
      else if (fSettings.ref_mode & fSettings.kTPCref) {
        this->minpt = 0.2;
        this->maxpt = 5;
        if (fSettings.use_primaries_cen) this->FillFromPrimariesAODTPC(refDist);
        else this->FillFromTracks(refDist, fSettings.tracktype);
        this->minpt = fSettings.minpt;
        this->maxpt = fSettings.maxpt;
      }
      else std::cout << "No valid reference detector chosen." << std::endl;
    }
  }
}



void AliForwardFlowUtil::FillDataCentral(TH2D*& centralDist)
{
  if (!fSettings.mc) {
    AliAODEvent* aodevent = dynamic_cast<AliAODEvent*>(fevent);
    this->fAODevent = aodevent;
      
    // Fill centralDist
    if (fSettings.useSPD) this->FillFromTracklets(centralDist);
    else if (fSettings.useITS) this->FillFromCentralClusters(centralDist);
    else  this->FillFromTracks(centralDist, fSettings.tracktype); //(fSettings.useTPC)
  }
  else {
    if(!fMCevent)
      throw std::runtime_error("Not MC as expected");

    if (fSettings.esd){
      
      // Fill centralDist
      if (fSettings.useITS){
        if (!fSettings.use_primaries_cen) this->FillFromTrackrefsITS(centralDist);
        else this->FillFromPrimariesITS(centralDist);
      }
      else if (fSettings.useSPD && fSettings.use_primaries_cen) this->FillFromPrimariesSPD(centralDist);
      else if (fSettings.useTPC && fSettings.use_primaries_cen) this->FillFromPrimariesTPC(centralDist);
      else std::cout << "No valid central detector chosen." << std::endl;
    }
    else { // AOD
      // Fill centralDist
      if (fSettings.useITS) {
        if (!fSettings.use_primaries_cen) this->FillFromCentralClusters(centralDist);
        else this->FillFromPrimariesAODITS(centralDist);
      }
      else if (fSettings.useSPD) {
        if (fSettings.use_primaries_cen) this->FillFromPrimariesAODSPD(centralDist); 
        else this->FillFromTracklets(centralDist);
      }
      else if (fSettings.useTPC){
        if (fSettings.use_primaries_cen) this->FillFromPrimariesAODTPC(centralDist); 
        else this->FillFromTracks(centralDist, fSettings.tracktype);
      }
      else std::cout << "No valid central detector chosen." << std::endl;
    }
  }
}


Double_t AliForwardFlowUtil::GetZ(){
  if (this->fSettings.mc) return fMCevent->GetPrimaryVertex()->GetZ();
  else return fevent->GetPrimaryVertex()->GetZ();
}


Double_t AliForwardFlowUtil::GetCentrality(TString centrality_estimator){
  AliMultSelection *MultSelection;
  MultSelection  = (AliMultSelection*)fevent->FindListObject("MultSelection");
  return MultSelection->GetMultiplicityPercentile(centrality_estimator);
}



void AliForwardFlowUtil::FillFromTrackrefsITS(TH2D*& fwd)
{
  Int_t nTracks   = fMCevent->GetNumberOfTracks();// stack->GetNtrack();

  for (Int_t iTr = 0; iTr < nTracks; iTr++) {
    AliMCParticle* particle =
      static_cast<AliMCParticle*>(fMCevent->GetTrack(iTr));

    // Check if this charged and a primary
    if (particle->Charge() == 0) continue;

    // IF the track corresponds to a primary, pass that as both
    // arguments.
    ProcessTrackITS(particle, fwd);
  } // Loop over tracks
}

void AliForwardFlowUtil::FillFromTrackrefsFMD(TH2D*& fwd) 
{
  Int_t nTracks   = fMCevent->GetNumberOfTracks();// stack->GetNtrack();

  if (this->fSettings.fMaxConsequtiveStrips == 0) {
    for (Int_t iTr = 0; iTr < nTracks; iTr++) {
      AliMCParticle* p = static_cast< AliMCParticle* >(fMCevent->GetTrack(iTr));

      // Ignore things that do not make a signal in the FMD
      AliTrackReference* tr = this->IsHitFMD(p);
      if (tr && p->Charge() != 0){

        Double_t phi_tr = this->GetTrackRefPhi(tr);
        Double_t eta_tr = this->GetTrackRefEta(tr);
        if ((phi_tr > 0) & (eta_tr > -10)) fwd->Fill(eta_tr,phi_tr,1);
      }
    }
  }
  else {
    for (Int_t iTr = 0; iTr < nTracks; iTr++) {
      AliMCParticle* particle =
        static_cast<AliMCParticle*>(fMCevent->GetTrack(iTr));

      // Check if this charged and a primary
      if (particle->Charge() == 0) continue;

      // IF the track corresponds to a primary, pass that as both
      // arguments.
      ProcessTrack(particle, fwd);
    } // Loop over tracks
  }
}



Bool_t
AliForwardFlowUtil::ProcessTrackITS(AliMCParticle* particle,TH2D*& cen)
{

  if (!particle) return false;

  Int_t              nTrRef = particle->GetNumberOfTrackReferences();

  BeginTrackRefs();

  // Double_t oTheta= 0;
  for (Int_t iTrRef = 0; iTrRef < nTrRef; iTrRef++) {
    AliTrackReference* ref = particle->GetTrackReference(iTrRef);

    // Check existence
    if (!ref) continue;

    // Check that we hit an Base element
    if (AliTrackReference::kITS == ref->DetectorId()) {

      // We are interested if it produced a signal, not only a hit in the support structure.
      // This is an envelop around the active area
      if (ref->R() > 3.5 && ref->R() < 4.5 && TMath::Abs(ref->Z()) < 14.1) {
        if (!fStored){
          fStored = ref;

          Double_t phi_tr = this->GetTrackRefPhi(ref);
          Double_t eta_tr = this->GetTrackRefEta(ref);

          
          if ((phi_tr > 0) & (eta_tr > -10)) cen->Fill(eta_tr,phi_tr,1);

          // if (dodNdeta) dNdeta->Fill(eta_tr_spd,1);
        }
      }
    }
  }
  return true;
}



Double_t AliForwardFlowUtil::GetTrackRefPhi(AliTrackReference* ref) {

  const AliVVertex* vertex = fMCevent->GetPrimaryVertex();
  // Calculate the vector pointing from the vertex to the track reference on the detector
  Double_t x      = ref->X() - vertex->GetX();
  Double_t y      = ref->Y() - vertex->GetY();
  Double_t phiR   = TMath::ATan2(y,x);
  if (phiR < 0) {
    phiR += 2*TMath::Pi();
  }

  return phiR;
}


Double_t AliForwardFlowUtil::GetTrackRefEta(AliTrackReference* ref) {

  const AliVVertex* vertex = fMCevent->GetPrimaryVertex();
  // Calculate the vector pointing from the vertex to the track reference on the detector
  Double_t x      = ref->X() - vertex->GetX();
  Double_t y      = ref->Y() - vertex->GetY();
  Double_t z      = ref->Z() - vertex->GetZ();
  Double_t rr     = TMath::Sqrt(x * x + y * y);
  Double_t thetaR = TMath::ATan2(rr, z);
  // Correct angles
  if (thetaR < 0) {
    thetaR += 2*TMath::Pi();
  }
  Double_t eta = -TMath::Log(TMath::Tan(thetaR / 2));

  return eta;
}


Bool_t
AliForwardFlowUtil::ProcessTrack(AliMCParticle* particle,TH2D*& fwd)
{
  // Check the returned particle
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

    AliTrackReference* test = ProcessRef(particle, ref,fwd);
    if (test) store = test;

  } // Loop over track references
  if (!store) return true; // Nothing found

  StoreParticle(particle, store,fwd);
  EndTrackRefs();

  return true;
}



AliTrackReference*
AliForwardFlowUtil::ProcessRef(AliMCParticle* particle, AliTrackReference*ref,TH2D*& fwd)
{
  // Process track references of a track
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
    else if (nT > this->fSettings.fMaxConsequtiveStrips) {

      used = true;
    }
  }
  if (used) {

    // Int_t nnT   = TMath::Abs(fState.oldStrip - fState.startStrip) + 1;
    StoreParticle(particle, fState.longest,fwd);
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
				     AliTrackReference*   ref,TH2D*& fwd)
{
  Double_t phi_tr = this->GetTrackRefPhi(particle); //Wrap02pi
  Double_t eta_tr = this->GetTrackRefEta(particle);

  if ((phi_tr > 0) & (eta_tr > -10)) fwd->Fill(eta_tr,phi_tr,1);
  return;
}



void
AliForwardFlowUtil::BeginTrackRefs()
{
  fState.Clear(true);
  fStored = 0;
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

void AliForwardFlowUtil::FillFromPrimariesAODTPC(TH2D*& cen) const
{
  Int_t nTracksMC = fMCevent->GetNumberOfTracks();
  for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
    AliAODMCParticle* p = static_cast< AliAODMCParticle* >(fMCevent->GetTrack(iTr));
    if (!p->IsPhysicalPrimary()) continue;
    if (p->Charge() == 0) continue;

    Double_t eta = p->Eta();
    if (TMath::Abs(eta) < 1.7) {
      if (p->Pt()>=this->minpt){// && p->Pt()<=this->maxpt){
        cen->Fill(eta,p->Phi(),1);
      }
    }
  }
}



void AliForwardFlowUtil::FillFromPrimariesAODSPD(TH2D*& cen) const
{
  Int_t nTracksMC = fMCevent->GetNumberOfTracks();
  for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
    AliAODMCParticle* p = static_cast< AliAODMCParticle* >(fMCevent->GetTrack(iTr));
    if (!p->IsPhysicalPrimary()) continue;
    if (p->Charge() == 0) continue;

    Double_t eta = p->Eta();
    if (TMath::Abs(eta) < 1.7 && p->Pt() >= 0.05) {
        cen->Fill(eta,p->Phi(),1);
    }
  }
}

void AliForwardFlowUtil::FillFromPrimariesAODITS(TH2D*& cen) const
{
  Int_t nTracksMC = fMCevent->GetNumberOfTracks();
  for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
    AliAODMCParticle* p = static_cast< AliAODMCParticle* >(fMCevent->GetTrack(iTr));
    if (!p->IsPhysicalPrimary()) continue;
    if (p->Charge() == 0) continue;

    Double_t eta = p->Eta();
    if (TMath::Abs(eta) < 1.7) {
        cen->Fill(eta,p->Phi(),1);
        // if (dodNdeta) dNdeta->Fill(eta,1);
    }
  }
}

Double_t AliForwardFlowUtil::GetTrackRefPhi(AliMCParticle* p) {
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
    return -1;
  }
  const AliVVertex* vertex = fMCevent->GetPrimaryVertex();
  // Calculate the vector pointing from the vertex to the track reference on the detector
  Double_t x      = ref->X() - vertex->GetX();
  Double_t y      = ref->Y() - vertex->GetY();
  Double_t phiR   = TMath::ATan2(y,x);
  // Correct angles
  if (phiR < 0) {
    phiR += 2*TMath::Pi();
  }
  return phiR;
}


Double_t AliForwardFlowUtil::GetTrackRefEta(AliMCParticle* p) {
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
    return -10; // dummy value
  }
  const AliVVertex* vertex = fMCevent->GetPrimaryVertex();
  // Calculate the vector pointing from the vertex to the track reference on the detector
  Double_t x      = ref->X() - vertex->GetX();
  Double_t y      = ref->Y() - vertex->GetY();
  Double_t z      = ref->Z() - vertex->GetZ();
  Double_t rr     = TMath::Sqrt(x * x + y * y);
  Double_t thetaR = TMath::ATan2(rr, z);
  // Correct angles
  if (thetaR < 0) {
    thetaR += 2*TMath::Pi();
  }

  Double_t eta = -TMath::Log(TMath::Tan(thetaR / 2));
  return eta;
}

void AliForwardFlowUtil::FillFromPrimariesTPC(TH2D*& cen) const
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
        // if (dodNdeta) dNdeta->Fill(eta,1);
      }
    }
  }
}



void AliForwardFlowUtil::FillFromPrimariesSPD(TH2D*& cen) const
{
  Int_t nTracksMC = fMCevent->GetNumberOfTracks();

  for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
    AliMCParticle* p = static_cast< AliMCParticle* >(fMCevent->GetTrack(iTr));
    if (!p->IsPhysicalPrimary()) continue;
    if (p->Charge() == 0) continue;

    Double_t eta = p->Eta();
    if (TMath::Abs(eta) < 1.7) {
      if (p->Pt()>=0.05){
        cen->Fill(eta,p->Phi(),1);
        // if (dodNdeta) dNdeta->Fill(eta,1);
      }
    }
  }
}



void AliForwardFlowUtil::FillFromPrimariesITS(TH2D*& cen) const
{
  Int_t nTracksMC = fMCevent->GetNumberOfTracks();

  for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
    AliMCParticle* p = static_cast< AliMCParticle* >(fMCevent->GetTrack(iTr));
    if (!p->IsPhysicalPrimary()) continue;
    if (p->Charge() == 0) continue;

    Double_t eta = p->Eta();
    if (TMath::Abs(eta) < 1.7) {
        cen->Fill(eta,p->Phi(),1);
        // if (dodNdeta) dNdeta->Fill(eta,1);
    }
  } 
}

void AliForwardFlowUtil::FillFromPrimariesFMD(TH2D*& fwd) const
{
  Int_t nTracksMC = fMCevent->GetNumberOfTracks();

  for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
    AliMCParticle* p = static_cast< AliMCParticle* >(fMCevent->GetTrack(iTr));
    if (!p->IsPhysicalPrimary()) continue;
    if (p->Charge() == 0) continue;

    Double_t eta = p->Eta();
    if (eta < 5 /*fwd->GetXaxis()-GetXmax()*/ && eta > -3.5 /*fwd->GetXaxis()-GetXmin()*/) {
      if (TMath::Abs(eta) >= 1.7){
        fwd->Fill(eta,p->Phi(),1);
        // if (dodNdeta) dNdeta->Fill(eta,1);
      }
    }
  }
}


void AliForwardFlowUtil::FillFromPrimariesAODFMD(TH2D*& fwd) const
{
  Int_t nTracksMC = fMCevent->GetNumberOfTracks();

  for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
    AliAODMCParticle* p = static_cast< AliAODMCParticle* >(fMCevent->GetTrack(iTr));
    if (!p->IsPhysicalPrimary()) continue;
    if (p->Charge() == 0) continue;

    if (p->Eta() < 5 /*fwd->GetXaxis()-GetXmax()*/ && p->Eta() > -3.5 /*fwd->GetXaxis()-GetXmin()*/) {
      if (TMath::Abs(p->Eta()) >= 1.7){
        fwd->Fill(p->Eta(),p->Phi(),1);
        // if (dodNdeta) dNdeta->Fill(p->Eta(),1);
      }
    }
  }
}


void AliForwardFlowUtil::FillFromTracklets(TH2D*& cen) const {
  AliVMultiplicity* mult = this->fAODevent->GetMultiplicity();
  Int_t nTracklets = mult->GetNumberOfTracklets();

  for (Int_t i = 0; i < nTracklets; i++) {
    // Using a dphi cut in units of mrad; This cut is motivated in
    // https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/lmilano/2017-Aug-11-analysis_note-note.pdf
    auto dphi  = mult->GetDeltaPhi(i);
     if (TMath::Abs(dphi) * 1000 > 5) {
       continue;
     }
    auto eta   = -TMath::Log(TMath::Tan(mult->GetTheta(i)/2));
    // Drop everything outside of -1.7 < eta 1.7 to avoid overlas with the FMD
    if (eta < -1.7 || eta > 1.7) {
      continue;
    }
    auto phi   = mult->GetPhi(i) + 39./34.*dphi;

    cen->Fill(eta,phi, 1);
  }
}


void AliForwardFlowUtil::FillFromCentralClusters(TH2D*& cen) const {
  AliAODCentralMult* aodcmult = static_cast<AliAODCentralMult*>(fAODevent->FindListObject("CentralClusters"));
  cen = &aodcmult->GetHistogram();
}

void AliForwardFlowUtil::FillFromForwardClusters(TH2D*& fwd) {
  AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(fAODevent->FindListObject("Forward"));
  fwd = &aodfmult->GetHistogram();
}

void AliForwardFlowUtil::FillFromTracks(TH2D*& cen, UInt_t tracktype) const {
  Int_t  iTracks(fevent->GetNumberOfTracks());
  for(Int_t i(0); i < iTracks; i++) {

  // loop  over  all  the  tracks
    AliAODTrack* track = static_cast<AliAODTrack *>(fAODevent->GetTrack(i));

    if (track->TestFilterBit(tracktype) && track->GetTPCNcls() > fSettings.fnoClusters){
      if (track->Pt() < this->minpt || track->Pt() > this->maxpt) continue;

      if( fSettings.fCutChargedDCAzMax > 0. || fSettings.fCutChargedDCAxyMax > 0.){
        Double_t dTrackXYZ[3]  = {0.,0.,0.};
        Double_t dVertexXYZ[3] = {0.,0.,0.};
        Double_t dDCAXYZ[3]    = {0.,0.,0.};
        const AliAODVertex* vertex = fAODevent->GetPrimaryVertex();

        track->GetXYZ(dTrackXYZ);
        vertex->GetXYZ(dVertexXYZ);

        for(Short_t i(0); i < 3; i++) { dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i]; }

        if(fSettings.fCutChargedDCAzMax > 0. && TMath::Abs(dDCAXYZ[2]) > fSettings.fCutChargedDCAzMax) continue;
        if(fSettings.fCutChargedDCAxyMax > 0. && TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]) > fSettings.fCutChargedDCAxyMax) continue;
      }
      cen->Fill(track->Eta(),track->Phi(), 1);
    }
  }
}


AliTrackReference* AliForwardFlowUtil::IsHitFMD(AliMCParticle* p) {
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
