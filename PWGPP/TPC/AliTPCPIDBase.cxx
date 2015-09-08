#include "TChain.h"
#include "TF1.h"
#include "TRandom3.h"

#include "AliMCParticle.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliVTrack.h"
#include "AliExternalTrackParam.h"
#include "AliVVertex.h"
#include "AliAnalysisFilter.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliESDv0KineCuts.h"
#include "AliESDv0.h"

#include "AliTPCPIDBase.h"

/*
This class is a base class for all other
analysis tasks that use V0's.
It provides basics for V0 identification.
In addition, some further basic functions are provided.

Class written by Benjamin Hess.
Contact: bhess@cern.ch
*/

ClassImp(AliTPCPIDBase)

Double_t AliTPCPIDBase::fgCutGeo = 1.;   
Double_t AliTPCPIDBase::fgCutNcr = 0.85; 
Double_t AliTPCPIDBase::fgCutNcl = 0.7;  

UShort_t AliTPCPIDBase::fgCutPureNcl = 60;

//________________________________________________________________________
AliTPCPIDBase::AliTPCPIDBase()
  : AliAnalysisTaskSE()
  , fEvent(0x0)
  , fESD(0x0)
  , fMC(0x0)
  , fPIDResponse(0x0)
  , fV0KineCuts(0x0)
  , fAnaUtils(0x0)
  , fIsPbpOrpPb(kFALSE)
  , fUsePhiCut(kFALSE)
  , fTPCcutType(kNoCut)
  , fZvtxCutEvent(10.0)
  , fEtaCut(0.9)
  , fPhiCutLow(0x0)
  , fPhiCutHigh(0x0)
  , fRandom(0x0)
  , fTrackFilter(0x0)
  , fNumTagsStored(0)
  , fV0tags(0x0)
  , fStoreMotherIndex(kFALSE)
  , fV0motherIndex(0x0)
{
  // default Constructor
  
  fRandom = new TRandom3(0); // 0 means random seed


  // Question: Is this the right place to initialize these functions?
  // Will it work on proof? i.e. will they be streamed to the workers?
  // Also one should add getters and setters
  fPhiCutLow  = new TF1("StdPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 100);
  fPhiCutHigh = new TF1("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 100);
}

//________________________________________________________________________
AliTPCPIDBase::AliTPCPIDBase(const char *name)
  : AliAnalysisTaskSE(name)
  , fEvent(0x0)
  , fESD(0x0)
  , fMC(0x0)
  , fPIDResponse(0x0)
  , fV0KineCuts(0x0)
  , fAnaUtils(0x0)
  , fIsPbpOrpPb(kFALSE)
  , fUsePhiCut(kFALSE)
  , fTPCcutType(kNoCut)
  , fZvtxCutEvent(10.0)
  , fEtaCut(0.9)
  , fPhiCutLow(0x0)
  , fPhiCutHigh(0x0)
  , fRandom(0x0)
  , fTrackFilter(0x0)
  , fNumTagsStored(0)
  , fV0tags(0x0)
  , fStoreMotherIndex(kFALSE)
  , fV0motherIndex(0x0)
{
  // Constructor
  
  fRandom = new TRandom3(0); // 0 means random seed

  
  // Question: Is this the right place to initialize these functions?
  // Will it work on proof? i.e. will they be streamed to the workers?
  // Also one should add getters and setters
  fPhiCutLow  = new TF1("StdPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 100);
  fPhiCutHigh = new TF1("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 100);

  DefineInput (0, TChain::Class());
  DefineOutput(0,  TTree::Class());
}


//________________________________________________________________________
AliTPCPIDBase::~AliTPCPIDBase()
{
  // dtor
  
  delete fPhiCutLow;
  fPhiCutLow = 0;
  
  delete fPhiCutHigh;
  fPhiCutHigh = 0;
  
  delete fTrackFilter;
  fTrackFilter = 0;
  
  delete fRandom;
  fRandom = 0;
  
  delete fV0KineCuts;
  fV0KineCuts = 0;
  
  delete [] fV0tags;
  fV0tags = 0;
  fNumTagsStored = 0;
  
  delete [] fV0motherIndex;
  fV0motherIndex = 0;
  
  delete fAnaUtils;
  fAnaUtils = 0;
}


//________________________________________________________________________
void AliTPCPIDBase::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  // Input hander
  AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  
  if (!inputHandler) {
    AliFatal("Input handler needed");
    fPIDResponse = 0x0;
    
    return;
  }

  // PID response object
  fPIDResponse = inputHandler->GetPIDResponse();
  if (!fPIDResponse)
    AliError("PIDResponse object was not created");
  
  // V0 Kine cuts 
  fV0KineCuts = new AliESDv0KineCuts;
  fV0KineCuts->SetGammaCutChi2NDF(5.);
  // Only accept V0el with prod. radius within 45 cm -> PID will by systematically biased for larger values!
  Float_t gammaProdVertexRadiusCuts[2] = { 3.0, 45. }; 
  fV0KineCuts->SetGammaCutVertexR(&gammaProdVertexRadiusCuts[0]);
  
  // Default analysis utils
  fAnaUtils = new AliAnalysisUtils();
  
  // Not used yet, but to be save, forward vertex z cut to analysis utils object
  fAnaUtils->SetMaxVtxZ(fZvtxCutEvent);
}


//________________________________________________________________________
void AliTPCPIDBase::UserExec(Option_t *)
{
  // Main loop
  // Called for each event
}      

//________________________________________________________________________
void AliTPCPIDBase::Terminate(const Option_t *)
{
  // Called once at the end of the query
}


//_____________________________________________________________________________
Double_t AliTPCPIDBase::GetPhiPrime(Double_t phi, Double_t magField, Int_t charge) const
{
  // Get phiPrime which is the cut variable to reject high pT tracks near edges

  if(magField < 0)    // for negatve polarity field
    phi = TMath::TwoPi() - phi;
    
  if(charge < 0) // for negatve charge
    phi = TMath::TwoPi() - phi;
  
  phi += TMath::Pi() / 18.0; // to center gap in the middle
  phi = fmod(phi, TMath::Pi() / 9.0);
  return phi;
}


//______________________________________________________________________________
Bool_t AliTPCPIDBase::PhiPrimeCut(Double_t trackPt, Double_t trackPhi, Short_t trackCharge, Double_t magField) const
{
  // Apply phi' cut to given track parameters
  
  if (trackPt < 2.0)
    return kTRUE;
  
  Double_t phiPrime = GetPhiPrime(trackPhi, magField, trackCharge);
  
  if (phiPrime < fPhiCutHigh->Eval(trackPt) && phiPrime > fPhiCutLow->Eval(trackPt))
    return kFALSE; // reject track
    
    return kTRUE;
}


//______________________________________________________________________________
Bool_t AliTPCPIDBase::PhiPrimeCut(const AliVTrack* track, Double_t magField) const
{
  return PhiPrimeCut(track->Pt(), track->Phi(), track->Charge(), magField);
}


//______________________________________________________________________________
Bool_t AliTPCPIDBase::GetVertexIsOk(AliVEvent* event, Bool_t doVtxZcut) const
{
  // Check whether vertex ful-fills quality requirements.
  // Apply cut on z-position of vertex if doVtxZcut = kTRUE.
  
  AliAODEvent* aod = 0x0;
  AliESDEvent* esd = 0x0;
  
  aod = dynamic_cast<AliAODEvent*>(event);
  if (!aod) {
    esd = dynamic_cast<AliESDEvent*>(event);
    
    if (!esd) {
      AliError("Event seems to be neither AOD nor ESD!");
      return kFALSE;
    }
  }
    
  if (fIsPbpOrpPb) {
    const AliVVertex* trkVtx = (aod ? dynamic_cast<const AliVVertex*>(aod->GetPrimaryVertex()) :
                                      dynamic_cast<const AliVVertex*>(esd->GetPrimaryVertex()));
      
    if (!trkVtx || trkVtx->GetNContributors() <= 0)
      return kFALSE;
      
    TString vtxTtl = trkVtx->GetTitle();
    if (!vtxTtl.Contains("VertexerTracks"))
      return kFALSE;
      
    Float_t zvtx = trkVtx->GetZ();
    const AliVVertex* spdVtx = (aod ? dynamic_cast<const AliVVertex*>(aod->GetPrimaryVertexSPD()) :
                                      dynamic_cast<const AliVVertex*>(esd->GetPrimaryVertexSPD()));
    if (spdVtx->GetNContributors() <= 0)
      return kFALSE;
      
    Double_t cov[6] = {0};
    spdVtx->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (spdVtx->IsFromVertexerZ() && (zRes > 0.25))
      return kFALSE;
      
    if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ()) > 0.5)
      return kFALSE;
    
    if (doVtxZcut) {
      if (TMath::Abs(zvtx) > fZvtxCutEvent) //Default: 10 cm
        return kFALSE;
    }
      
    return kTRUE;
  }
    
  
  // pp and PbPb
  const AliVVertex* primaryVertex = 0x0;
  if (aod) {
    primaryVertex = dynamic_cast<const AliVVertex*>(aod->GetPrimaryVertex());
    if (!primaryVertex || primaryVertex->GetNContributors() <= 0) 
      return kFALSE;
    
    // Reject TPC vertices
    TString primVtxTitle(primaryVertex->GetTitle());
    if (primVtxTitle.Contains("TPCVertex",TString::kIgnoreCase))
      return kFALSE;
  }
  else {
    primaryVertex = dynamic_cast<const AliVVertex*>(esd->GetPrimaryVertexTracks());
    if (!primaryVertex)
      return kFALSE;
    
    if (primaryVertex->GetNContributors() <= 0) {
      // Try SPD vertex
      primaryVertex = dynamic_cast<const AliVVertex*>(esd->GetPrimaryVertexSPD());
      if (!primaryVertex || primaryVertex->GetNContributors() <= 0) 
        return kFALSE;
    }
  }
  
  if (doVtxZcut) {
    if (TMath::Abs(primaryVertex->GetZ()) > fZvtxCutEvent) //Default: 10 cm
      return kFALSE;
  }
  
  return kTRUE;
}


//______________________________________________________________________________
Bool_t AliTPCPIDBase::GetIsPileUp(AliVEvent* event, AliTPCPIDBase::PileUpRejectionType pileUpRejectionType) const
{
  // Check whether event is a pile-up event according to current AnalysisUtils object.
  // If rejection type is kPileUpRejectionOff, kFALSE is returned.
  // In case of errors, the error is displayed and kTRUE is returned.
  
  if (pileUpRejectionType == kPileUpRejectionOff)
    return kFALSE;
  
  if (!event) {
    AliError("No event!");
    return kTRUE;
  }
  
  if (!fAnaUtils) {
    AliError("AnalysisUtils object not available!");
    return kTRUE;
  }
  
  if (pileUpRejectionType == kPileUpRejectionSPD)
    return fAnaUtils->IsPileUpSPD(event);
  else if (pileUpRejectionType == kPileUpRejectionMV)
    return fAnaUtils->IsPileUpMV(event);
  
  return kTRUE;
}


//______________________________________________________________________________
void AliTPCPIDBase::FillV0PIDlist(AliESDEvent* event)
{
  //
  // Fill the PID tag list
  //

  // If no event forwarded as parameter (default), cast current input event.
  // Dynamic cast to ESD events (DO NOTHING for AOD events)
  if (!event)
    event = dynamic_cast<AliESDEvent *>(InputEvent());
  
  // If this fails, do nothing
  if (!event) {
    AliError("Failed to retrieve ESD event. No V0's processed (only works for ESDs by now).");
    return;
  }
  
  if (!fV0KineCuts) {
    AliError("V0KineCuts not available!");
    return;
  }
  
  TString beamType(event->GetBeamType());
  
  if (beamType.CompareTo("Pb-Pb") == 0 || beamType.CompareTo("A-A") == 0) {
    fV0KineCuts->SetMode(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPbPb);
  }
  else {
    fV0KineCuts->SetMode(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPP); 
  }

  // V0 selection
  // set event
  fV0KineCuts->SetEvent(event);

  const Int_t numTracks = event->GetNumberOfTracks();
  fV0tags = new Char_t[numTracks];
  for (Int_t i = 0; i < numTracks; i++)
    fV0tags[i] = 0;
  
  if (fStoreMotherIndex)  {
    fV0motherIndex = new Int_t[numTracks];
    for (Int_t i = 0; i < numTracks; i++)
      fV0motherIndex[i] = -1;
  }
  
  fNumTagsStored = numTracks;
  
  // loop over V0 particles
  for (Int_t iv0 = 0; iv0 < event->GetNumberOfV0s(); iv0++) {
    AliESDv0* v0 = (AliESDv0*)event->GetV0(iv0);
 
    if (!v0)
      continue;
    
    // Reject onFly V0's <-> Only take V0's from offline V0 finder
    if (v0->GetOnFlyStatus())
      continue; 
  
    // Get the particle selection 
    Bool_t foundV0 = kFALSE;
    Int_t pdgV0 = 0, pdgP = 0, pdgN = 0;

    foundV0 = fV0KineCuts->ProcessV0(v0, pdgV0, pdgP, pdgN);
    if (!foundV0)
      continue;
    
    Int_t iTrackP = v0->GetPindex();  // positive track
    Int_t iTrackN = v0->GetNindex();  // negative track

    
    // Fill the Object arrays
    // positive particles
    if (pdgP == -11) {
      fV0tags[iTrackP] = 14;
    }
    else if (pdgP == 211) {
      fV0tags[iTrackP] = 15;
    }
    else if(pdgP == 2212) {
      fV0tags[iTrackP] = 16;
    }
    
    if (fStoreMotherIndex)
      fV0motherIndex[iTrackP] = iv0;

    // negative particles
    if( pdgN == 11){
      fV0tags[iTrackN] = -14;
    }
    else if( pdgN == -211){
      fV0tags[iTrackN] = -15;
    }
    else if( pdgN == -2212){
      fV0tags[iTrackN] = -16;
    }
    
    if (fStoreMotherIndex)
      fV0motherIndex[iTrackN] = iv0;
  }
}


//______________________________________________________________________________
void AliTPCPIDBase::ClearV0PIDlist()
{
  //
  // Clear the PID tag list
  //

  delete [] fV0tags;
  fV0tags = 0;
  
  delete [] fV0motherIndex;
  fV0motherIndex = 0;
  
  fNumTagsStored = 0;
}


//______________________________________________________________________________
Char_t AliTPCPIDBase::GetV0tag(Int_t trackIndex) const
{
  //
  // Get the tag for the corresponding trackIndex. Returns -99 in case of invalid index/tag list.
  //
  
  if (trackIndex < 0 || trackIndex >= fNumTagsStored || !fV0tags)
    return -99;
  
  return fV0tags[trackIndex];
}


//______________________________________________________________________________
Int_t AliTPCPIDBase::GetV0motherIndex(Int_t trackIndex) const
{
  //
  // Get the index of the V0 mother for the corresponding trackIndex. Returns -99 in case of invalid index/mother index list.
  //
  
  if (!fStoreMotherIndex || trackIndex < 0 || trackIndex >= fNumTagsStored || !fV0motherIndex)
    return -99;
  
  return fV0motherIndex[trackIndex];
}


//________________________________________________________________________
Bool_t AliTPCPIDBase::TPCCutMIGeo(const AliVTrack* track, const AliVEvent* evt, TTreeStream* streamer)
{
  //
  // TPC Cut MIGeo
  //

  if (!track || !evt)
    return kFALSE;
  
  const Short_t sign = track->Charge();
  Double_t xyz[50];
  Double_t pxpypz[50];
  Double_t cv[100];

  track->GetXYZ(xyz);
  track->GetPxPyPz(pxpypz);

  AliExternalTrackParam* par = new AliExternalTrackParam(xyz, pxpypz, cv, sign);
  const AliESDtrack dummy;

  const Double_t magField = evt->GetMagneticField();
  Double_t varGeom = dummy.GetLengthInActiveZone(par, 3, 236, magField, 0, 0);
  Double_t varNcr  = track->GetTPCClusterInfo(3, 1);
  Double_t varNcls = track->GetTPCsignalN();

  const Double_t varEval = 130. - 5. * TMath::Abs(1. / track->Pt());
  Bool_t cutGeom   = varGeom > fgCutGeo * varEval;
  Bool_t cutNcr    = varNcr  > fgCutNcr * varEval;
  Bool_t cutNcls   = varNcls > fgCutNcl * varEval;

  Bool_t kout = cutGeom && cutNcr && cutNcls;

  if (streamer) {
    Double_t dedx = track->GetTPCsignal();

    (*streamer)<<"tree"<<
      "param.="<< par<<
      "varGeom="<<varGeom<<
      "varNcr="<<varNcr<<
      "varNcls="<<varNcls<<
      
      "cutGeom="<<cutGeom<<
      "cutNcr="<<cutNcr<<
      "cutNcls="<<cutNcls<<
      
      "kout="<<kout<<
      "dedx="<<dedx<<
      
      "\n";
  }
  
  delete par;
  
  return kout;
}


//________________________________________________________________________
Bool_t AliTPCPIDBase::TPCnclCut(const AliVTrack* track)
{
  //
  // TPC Cut on TPCsignalN() only
  //

  if (!track)
    return kFALSE;
  
  return (track->GetTPCsignalN() >= fgCutPureNcl);
}
