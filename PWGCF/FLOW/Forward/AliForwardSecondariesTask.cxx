
//
// Calculate flow in the forward and central regions using the Q cumulants method.
//
// Inputs:
//  - AliAODEvent
//
// Outputs:
//  - AnalysisResults.root or forward_flow.root
//
#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TList.h>
#include <THn.h>
#include "AliGenEventHeaderTunedPbPb.h"
#include "TGraph.h"
#include "AliFMDStripIndex.h"
#include "AliFMDEncodedEdx.h"
#include <vector>
#include <numeric>
#include <assert.h>

#include "THn.h"
#include "TMath.h"
#include "TCutG.h"
#include "TParticle.h"
#include "AliFMDStripIndex.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"

#include "AliGenEventHeaderTunedPbPb.h"
#include "AliAODMCHeader.h"

#include "AliLog.h"
#include "AliForwardSecondariesTask.h"
#include "AliForwardGenericFramework.h"

#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODEvent.h"

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
#include "AliForwardTaskValidation.h"
using namespace std;
ClassImp(AliForwardSecondariesTask)
#if 0
; // For emacs
#endif

//_____________________________________________________________________
AliForwardSecondariesTask::AliForwardSecondariesTask() : AliAnalysisTaskSE(),
  fOutputList(0),    // output list
  fEventList(0),
  fDeltaList(0),
  fRandom(0),
  fTrackDensity(),
  fnoPrim(),
  fSettings(),
  fUtil(),
  fStored(0),
  fState(),
  fMaxConsequtiveStrips(0),
  fLowCutvalue(0),
  fTrackGammaToPi0(true),
  fStorage(nullptr),
  fdelta_phi_eta(),
  fdelta_eta_phi(),
  fdelta_phi_phi(),
  fdelta_eta_eta() 
  {
    //
    //  Default constructor
    //
  }

//_____________________________________________________________________
AliForwardSecondariesTask::AliForwardSecondariesTask(const char* name) : AliAnalysisTaskSE(name),
  fOutputList(0),    // output list
  fEventList(0),
  fDeltaList(0),
  fRandom(0),
  fTrackDensity(),  fnoPrim(),
  fSettings(),
  fUtil(),
  fStored(0),
  fState(),
  fMaxConsequtiveStrips(0),
  fLowCutvalue(0),
  fTrackGammaToPi0(true),
  fStorage(nullptr),
  fdelta_phi_eta(),
  fdelta_eta_phi(),
  fdelta_phi_phi(),
  fdelta_eta_eta() 
  {
    //
    //  Constructor
    //
    //  Parameters:
    //   name: Name of task
    //
    DefineInput(1, AliForwardTaskValidation::Class());
    //DefineOutput(1, TList::Class());
    DefineOutput(1, AliForwardFlowResultStorage::Class());

  }

//_____________________________________________________________________
void AliForwardSecondariesTask::UserCreateOutputObjects()
  {

  //
  //  Create output objects
  //
  std::cout << "AliForwardSecondariesTask::UserCreateOutputObjects()" << std::endl;

  this->fOutputList = new TList();
  this->fOutputList->SetOwner(kTRUE);

  TRandom fRandom = TRandom();        // random integer to use for creation of samples (used for error bars).
                                      // Needs to be created here, otherwise it will draw the same random number.
  fDeltaList = new TList();
  fDeltaList->SetName("Delta");
  fEventList = new TList();
  fEventList->SetName("EventInfo");

  fOutputList->Add(fEventList);
  fOutputList->Add(fDeltaList);

  Int_t dimensions = 3;

  Int_t bins_phi_eta[3] = {fSettings.fNZvtxBins, fSettings.fNPhiBins+1, fSettings.fNDiffEtaBins} ;
  Double_t xmin_phi_eta[3] = {fSettings.fZVtxAcceptanceLowEdge, -TMath::Pi(), fSettings.fEtaLowEdge};
  Double_t xmax_phi_eta[3] = {fSettings.fZVtxAcceptanceUpEdge, TMath::Pi(), fSettings.fEtaUpEdge}; 

  Int_t bins_eta_phi[3] = {fSettings.fNZvtxBins, 41, fSettings.fNPhiBins} ;
  Double_t xmin_eta_phi[3] = {fSettings.fZVtxAcceptanceLowEdge, -4.92, 0};
  Double_t xmax_eta_phi[3] = {fSettings.fZVtxAcceptanceUpEdge, 4.92, 2*TMath::Pi()}; 

  Int_t bins_eta_eta[3] = {fSettings.fNZvtxBins, 41, fSettings.fNDiffEtaBins} ;
  Double_t xmin_eta_eta[3] = {fSettings.fZVtxAcceptanceLowEdge, -4.92, fSettings.fEtaLowEdge};
  Double_t xmax_eta_eta[3] = {fSettings.fZVtxAcceptanceUpEdge,   4.92, fSettings.fEtaUpEdge}; 


  Int_t bins_phi_phi[3] = {fSettings.fNZvtxBins, fSettings.fNPhiBins + 1, fSettings.fNPhiBins} ;
  Double_t xmin_phi_phi[3] = {fSettings.fZVtxAcceptanceLowEdge, -TMath::Pi(), 0};
  Double_t xmax_phi_phi[3] = {fSettings.fZVtxAcceptanceUpEdge, TMath::Pi(), 2*TMath::Pi()}; 


  fdelta_phi_eta = new THnD("delta_phi_eta", "delta_phi_eta",dimensions,bins_phi_eta, xmin_phi_eta, xmax_phi_eta);
  fdelta_eta_phi = new THnD("delta_eta_phi", "delta_eta_phi",dimensions,bins_eta_phi, xmin_eta_phi, xmax_eta_phi);
  fdelta_eta_eta = new THnD("delta_eta_eta", "delta_eta_eta",dimensions,bins_eta_eta, xmin_eta_eta, xmax_eta_eta);
  fdelta_phi_phi = new THnD("delta_phi_phi", "delta_phi_phi",dimensions,bins_phi_phi, xmin_phi_phi, xmax_phi_phi);


  fdelta_phi_eta->GetAxis(0)->SetName("vertex");
  fdelta_phi_eta->GetAxis(1)->SetName("phi_mother - phi_tr");
  fdelta_phi_eta->GetAxis(2)->SetName("eta");

  fdelta_eta_phi->GetAxis(0)->SetName("vertex");
  fdelta_eta_phi->GetAxis(1)->SetName("eta_mother - eta_tr");
  fdelta_eta_phi->GetAxis(2)->SetName("phi");

  fdelta_eta_eta->GetAxis(0)->SetName("vertex");
  fdelta_eta_eta->GetAxis(1)->SetName("eta_mother - eta_tr");
  fdelta_eta_eta->GetAxis(2)->SetName("eta");

  fdelta_phi_phi->GetAxis(0)->SetName("vertex");
  fdelta_phi_phi->GetAxis(1)->SetName("phi_mother - phi_tr");
  fdelta_phi_phi->GetAxis(2)->SetName("phi");


  TList* list_delta_phi_eta = new TList(); list_delta_phi_eta->SetName("delta_phi_eta"); list_delta_phi_eta->Add(fdelta_phi_eta); fDeltaList->Add(list_delta_phi_eta);
  TList* list_delta_eta_phi = new TList(); list_delta_eta_phi->SetName("delta_eta_phi"); list_delta_eta_phi->Add(fdelta_eta_phi); fDeltaList->Add(list_delta_eta_phi);
  TList* list_delta_eta_eta = new TList(); list_delta_eta_eta->SetName("delta_eta_eta"); list_delta_eta_eta->Add(fdelta_eta_eta); fDeltaList->Add(list_delta_eta_eta);
  TList* list_delta_phi_phi = new TList(); list_delta_phi_phi->SetName("delta_phi_phi"); list_delta_phi_phi->Add(fdelta_phi_phi); fDeltaList->Add(list_delta_phi_phi);


  fEventList->Add(new TH1D("Vertex","Vertex",fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));

  Int_t bins_prim[2] = {fSettings.fNZvtxBins, fSettings.fNDiffEtaBins} ;
  Double_t xmin_prim[2] = {fSettings.fZVtxAcceptanceLowEdge, fSettings.fEtaLowEdge};
  Double_t xmax_prim[2] = {fSettings.fZVtxAcceptanceUpEdge, fSettings.fEtaUpEdge}; //
  Int_t dimensions_prim = 2;

  fnoPrim = new THnD("fnoPrim", "fnoPrim", dimensions_prim, bins_prim, xmin_prim, xmax_prim);

  fnoPrim->GetAxis(0)->SetName("vertex");
  fnoPrim->GetAxis(1)->SetName("eta_mother");

  TList* list_prim = new TList(); list_prim->SetName("prim"); list_prim->Add(fnoPrim); fDeltaList->Add(list_prim);



  fStorage = new AliForwardFlowResultStorage(fSettings.fileName, fOutputList);
  PostData(1, fStorage);
}


//_____________________________________________________________________
void AliForwardSecondariesTask::UserExec(Option_t *)
{
  //
  //  Parameters:
  //   option: Not used
  //

  AliForwardTaskValidation* ev_val = dynamic_cast<AliForwardTaskValidation*>(this->GetInputData(1));
  if (!ev_val->IsValidEvent()){
    PostData(1, this->fOutputList);
    return;
  }
  fUtil.fSettings = fSettings;
  if (fSettings.mc) fUtil.fMCevent = this->MCEvent();
  fUtil.fevent = fInputEvent;
  
  //Double_t cent = fUtil.GetCentrality(fSettings.centrality_estimator);

  AliMCEvent* fAOD = this->MCEvent();




  AliStack* stack = fAOD->Stack();
  if(!fAOD) {
    std::cout << "no mcevent" << std::endl;
    return;
  }
  if (!stack) {
    std::cout << "no stack" << std::endl;
    return;
  }

  // Disregard events without reconstructed vertex
  Double_t event_vtx_z = fUtil.GetZ();
  if (!(TMath::Abs(event_vtx_z) >= 0)) {
    return;
  }


  std::vector<Int_t> listOfMothers;

  Int_t nTracks   = fAOD->GetNumberOfTracks();// stack->GetNtrack();
  
  static_cast<TH1D*>(fEventList->FindObject("Vertex"))->Fill(event_vtx_z);

  if (fMaxConsequtiveStrips == 0){
  for (Int_t iTr = 0; iTr < nTracks; iTr++) {
    AliMCParticle* p = static_cast< AliMCParticle* >(this->MCEvent()->GetTrack(iTr));

    //AliTrackReference* tr = fUtil.IsHitFMD(p);
      for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) {
        //AliTrackReference* tr = p->GetTrackReference(iTrRef);
         // Ignore things that do not make a signal in the FMD
        //if (!tr) continue;    
        AliTrackReference* tr = p->GetTrackReference(iTrRef);

        if (!tr || AliTrackReference::kFMD != tr->DetectorId()) continue;    
        if (tr && p->Charge() == 0) continue;

         AliMCParticle* mother = GetMother(p);
         if (!mother) mother = p;

         Double_t phi_mother = mother->Phi();
         Double_t eta_mother = mother->Eta();

         Double_t *etaPhi = new Double_t[2];
         this->GetTrackRefEtaPhi(tr, etaPhi);

         Double_t phi_tr = etaPhi[1];
         Double_t eta_tr = etaPhi[0];

         Double_t weight = 1;
          // Int_t nuaeta = fSettings.nuaforward->GetXaxis()->FindBin(eta_tr);
          // Int_t nuaphi = fSettings.nuaforward->GetYaxis()->FindBin(phi_tr);
          // Int_t nuavtz = fSettings.nuaforward->GetZaxis()->FindBin(event_vtx_z);
          // weight = weight*fSettings.nuaforward->GetBinContent(nuaeta,nuaphi,nuavtz+10*fSettings.nua_runnumber);
          // if (weight == 0) continue;
         // (samples, vertex,phi_mother - phi_tr ,centrality,eta_mother,eta_tr)

        Double_t dphi_phi[3] = {event_vtx_z, WrapPi(phi_mother - phi_tr), phi_tr};//Wrap02pi
        Double_t dphi_eta[3] = {event_vtx_z, WrapPi(phi_mother - phi_tr), eta_tr};//Wrap02pi
        Double_t deta_phi[3] = {event_vtx_z, eta_tr - eta_mother, phi_tr};//wrappi
        Double_t deta_eta[3] = {event_vtx_z, eta_tr - eta_mother, eta_tr};


         fdelta_phi_eta->Fill(dphi_eta,weight);
         fdelta_eta_phi->Fill(deta_phi,weight);
         fdelta_phi_phi->Fill(dphi_phi,weight);
         fdelta_eta_eta->Fill(deta_eta,weight);


         Double_t x_prim[3] =  {event_vtx_z,eta_tr};
         Bool_t isNewPrimary = AddMotherIfFirstTimeSeen(mother,listOfMothers);
         if (!isNewPrimary){
           listOfMothers.push_back(mother->GetLabel());
           fnoPrim->Fill(x_prim,1);
         }
       }
     }
   }
  else{
    for (Int_t iTr = 0; iTr < nTracks; iTr++) {
      AliMCParticle* particle =
        static_cast<AliMCParticle*>(fAOD->GetTrack(iTr));

      // Check if this charged and a primary
      if (particle->Charge() == 0) continue;

      Bool_t isPrimary = stack->IsPhysicalPrimary(iTr) && iTr < this->MCEvent()->GetNumberOfPrimaries();

      AliMCParticle* mother = isPrimary ? particle : GetMother(iTr,fAOD);
      if (!mother) mother = particle;
      // IF the track corresponds to a primary, pass that as both
      // arguments.

      ProcessTrack(particle, mother,listOfMothers,event_vtx_z);
    }
  }
    PostData(1, fOutputList);
  return;
}


Bool_t
AliForwardSecondariesTask::ProcessTrack(AliMCParticle* particle, AliMCParticle* mother, 
                                        std::vector<Int_t> listOfMothers, Float_t event_vtx_z)
{
  // Check the returned particle
  //
  // Note: If particle refers to a primary, then particle and mother
  // refers to the same particle (the address are the same)
  //
  DGUARD(fDebug,3,"MC track density Process a track");
  if (!particle) return false;

  Int_t              nTrRef = particle->GetNumberOfTrackReferences();
  AliTrackReference* store  = 0;

  BeginTrackRefs();

  // Double_t oTheta= 0;
  for (Int_t iTrRef = 0; iTrRef < nTrRef; iTrRef++) {
    AliTrackReference* ref = particle->GetTrackReference(iTrRef);

    // Check existence
    if (!ref) continue;

    if (ref->DetectorId() != AliTrackReference::kFMD) continue;

    AliTrackReference* test = ProcessRef(particle, mother, ref,listOfMothers, event_vtx_z);
    if (test) store = test;

  } // Loop over track references
  if (!store) return true; // Nothing found

  StoreParticle(particle, mother, store, listOfMothers, event_vtx_z);
  EndTrackRefs();

  return true;
}

//____________________________________________________________________
void
AliForwardSecondariesTask::StoreParticle(AliMCParticle* particle, AliMCParticle* mother, AliTrackReference* ref,
                                         std::vector< Int_t > listOfMothers, Float_t event_vtx_z)
{
  UInt_t packed = ref->UserId();
  UShort_t detector, sector, strip;
  Char_t   ring;
  AliFMDStripIndex::Unpack(packed,detector,ring,sector,strip);
  TString inner = "I";
  //const Char_t* outer = "O";
  // Double_t v0cent = 0.0;
  // if (detector == 1){
  //   v0cent = 10;
  // }
  // if (detector == 2){
  //   if (TString(ring) == inner) v0cent = 30;
  //   else v0cent = 50;
  // }
  // if (detector == 3){
  //   if (ring == inner) v0cent = 70;
  //   else v0cent = 90;
  // }

  Double_t eta_mother = mother->Eta();
  Double_t phi_mother = (mother->Phi());//Wrap02pi

  Double_t *etaPhi = new Double_t[2];
  this->GetTrackRefEtaPhi(ref, etaPhi);//ref, etaPhi

  Double_t phi_tr = (etaPhi[1]); //Wrap02pi
  Double_t eta_tr = etaPhi[0];
  Double_t dphi_phi[3] = {event_vtx_z, WrapPi(phi_mother - phi_tr), phi_tr};//Wrap02pi
  Double_t dphi_eta[3] = {event_vtx_z, WrapPi(phi_mother - phi_tr), eta_tr};//Wrap02pi
  Double_t deta_phi[3] = {event_vtx_z, eta_tr - eta_mother, phi_tr};//wrappi
  Double_t deta_eta[3] = {event_vtx_z, eta_tr - eta_mother, eta_tr};

  fdelta_phi_phi->Fill(dphi_phi,1);
  fdelta_phi_eta->Fill(dphi_eta,1);
  fdelta_eta_phi->Fill(deta_phi,1);
  fdelta_eta_eta->Fill(deta_eta,1);

  Double_t x_prim[4] =  {event_vtx_z,eta_tr};
  Bool_t isNewPrimary = AddMotherIfFirstTimeSeen(mother,listOfMothers);
  if (!isNewPrimary){
    listOfMothers.push_back(mother->GetLabel());
    fnoPrim->Fill(x_prim,1);
  }
  return;
}


void
AliForwardSecondariesTask::State::Clear(Bool_t alsoCount)
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


//____________________________________________________________________
void
AliForwardSecondariesTask::BeginTrackRefs()
{
  fState.Clear(true);
  fStored = 0;
}

void
AliForwardSecondariesTask::EndTrackRefs()
{
  fState.Clear(true);
}

/// Modulo for float numbers
///
/// \param x nominator
/// \param y denominator
///
/// \return Rest of the rounded down division
Double_t AliForwardSecondariesTask::Mod(Double_t x, Double_t y) {
  if (0 == y)
    return x;
  return x - y * floor(x/y);
}


/// Wrap angle around 0 and 2pi
Double_t AliForwardSecondariesTask::Wrap02pi(Double_t angle) {
  const Double_t two_pi = 6.283185307179586;
  Double_t lower_edge = 0;
  Double_t interval = two_pi;
  if (lower_edge <= angle && angle < two_pi) {
    return angle;
  }
  return Mod(angle - lower_edge, interval) + lower_edge;
}


Double_t AliForwardSecondariesTask::WrapPi(Double_t phi){
   if (phi > TMath::Pi()){
     phi = phi - 2*TMath::TwoPi();
   }
   if (phi < -TMath::Pi()){
     phi = phi + 2*TMath::TwoPi();
   }
   return phi;
 }


Bool_t AliForwardSecondariesTask::AddMotherIfFirstTimeSeen(AliMCParticle* p, std::vector<Int_t> v){

  //Checking if v contains elements (is empty):
  if(v.empty()){
    return false;
  }
  Int_t x = p->GetLabel();
  if(std::find(v.begin(), v.end(), x) != v.end()) {
      /* v contains x */
    return true;
  } else {
      /* v does not contain x */
    return false;
  }
}


AliTrackReference*
AliForwardSecondariesTask::ProcessRef(AliMCParticle* particle, AliMCParticle* mother, AliTrackReference* ref,
                                      std::vector< Int_t > listOfMothers, Float_t event_vtx_z)
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
      if (fDebug) Info("Process", "New because new sector");
      used = true;
    }
    else if (nT > fMaxConsequtiveStrips) {
      if (fDebug) Info("Process", "New because too long: %d (%d,%d,%d)",
		       fState.nStrips, t, fState.oldStrip, fState.startStrip);
      used = true;
    }
  }
  if (used) {
    if (fDebug)
      Info("Process", "I=%p L=%p D=%d (was %d), R=%c (was %c), "
	   "S=%2d (was %2d) t=%3d (was %3d) nT=%3d/%4d",
	   ref, fState.longest,
	   d, fState.oldDetector,
	   r, fState.oldRing,
	   s, fState.oldSector,
	   t, fState.oldStrip,
	   fState.nStrips, fMaxConsequtiveStrips);
    StoreParticle(particle, mother, fState.longest,listOfMothers, event_vtx_z);
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

  // Debug output
  if (fDebug) {
    if (t == fState.startStrip)
      Info("Process", "New cluster starting at FMD%d%c[%2d,%3d]",
	   d, r, s, t);
    else
      Info("Process", "Adding to cluster starting at FMD%d%c[%2d,%3d], "
	   "length=%3d (now in %3d, previous %3d)",
	   d, r, s, fState.startStrip, fState.nStrips, t, fState.oldStrip);
  }

  // The longest passage is determined through the angle
  Double_t ang  = GetTrackRefTheta(ref);
  if (ang > fState.angle) {
    fState.longest = ref;
    fState.angle   = ang;
  }
  return fState.longest;
}

Double_t
AliForwardSecondariesTask::GetTrackRefTheta(const AliTrackReference* ref) const
{
  // Get the incidient angle of the track reference.
  const AliVVertex* vertex = this->MCEvent()->GetPrimaryVertex();
  // Calculate the vector pointing from the vertex to the track reference on the detector
  Double_t x      = ref->X() - vertex->GetX();
  Double_t y      = ref->Y() - vertex->GetY();
  Double_t z      = ref->Z() - vertex->GetZ();
  Double_t rr   = TMath::Sqrt(x*x+y*y);
  Double_t theta= TMath::ATan2(rr,z);
  Double_t ang  = TMath::Abs(TMath::Pi()-theta);
  return ang;
}

AliMCParticle*
AliForwardSecondariesTask::GetMother(Int_t iTr, const AliMCEvent* event) const
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


AliMCParticle* AliForwardSecondariesTask::GetMother(AliMCParticle* p) {
  // Recurses until the mother IsPhysicalPrimary
  // Return NULL if no mother was found
  AliMCEvent* event = this->MCEvent();
  // GetLabel() is the index on the Stack!
  // event->Stack()->IsPhysicalPrimary(p->GetLabel());
  Bool_t isPP = this->IsRedefinedPhysicalPrimary(p);
  // Return this particle if it is stable
  if (isPP) {
    return p;
  }
  else {
    // No stable particle found and no mother left !?
    if (p->GetMother() < 0) {
      return 0x0;
    }
    AliMCParticle* ancestor = dynamic_cast< AliMCParticle* >(event->GetTrack(p->GetMother()));
    return GetMother(ancestor);
  }
}


Bool_t AliForwardSecondariesTask::IsRedefinedPhysicalPrimary(AliMCParticle* p) {
  AliMCEvent* event = this->MCEvent();
  // Is this a pi0 which was produced as a primary particle?
  if (TMath::Abs(p->PdgCode()) == 111 /*pi0*/ &&
      p->GetLabel() < event->Stack()->GetNprimary()) {
      std::cout << "found a pi0" << std::endl;

    return true;
  }
  // Is it a Physical Primary by the standard definition?
  Bool_t isPPStandardDef = event->Stack()->IsPhysicalPrimary(p->GetLabel());
  AliMCParticle *pi0Candidate = dynamic_cast< AliMCParticle* >(event->GetTrack(p->GetMother()));
  // Check if this is a primary originating from a pi0
  if (isPPStandardDef && pi0Candidate) {
    if (TMath::Abs(pi0Candidate->PdgCode()) == 111/*pi0*/) {
        std::cout << "found wrong pi0Candidate" << std::endl;

      return false;//false; // Don't allow stable particles stemming from pi0!
    }
  }
  return isPPStandardDef;
}


void AliForwardSecondariesTask::GetTrackRefEtaPhi(AliTrackReference* ref, Double_t* etaPhi) {

  const AliVVertex* vertex = this->MCEvent()->GetPrimaryVertex();
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
}


//_____________________________________________________________________
void AliForwardSecondariesTask::Terminate(Option_t */*option*/)
{
  return;
}


//_____________________________________________________________________
//
//
// EOF
