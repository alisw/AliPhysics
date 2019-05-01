
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
#include "AliForwardQCumulantRun2.h"
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
  fSettings(),
  fUtil(),
  fMultTOFLowCut(),
  fMultTOFHighCut(),
  fMultCentLowCut(),
  fdNdeta(0),
  fPiCheck(0),
  fdNdetaOrigin(0),
  fxray(0),
  fNsecondaries(0),
  fNprimaries(0),
  fITS(0),
  fFMD1(0),
  fFMD2(0),
  fFMD3(0),
  fPipe(0),
  fEarlyDecay(0),
  phihist(),
  fStored(0),
  fState(),
  fMaxConsequtiveStrips(2),
  fLowCutvalue(0),
  fTrackGammaToPi0(true)
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
  fTrackDensity(),
  fSettings(),
  fUtil(),
  fMultTOFLowCut(),
  fMultTOFHighCut(),
  fMultCentLowCut(),
  fdNdeta(0),
  fPiCheck(0),
  fdNdetaOrigin(0),
  fxray(0),
  fNsecondaries(0),
  fNprimaries(0),
  fITS(0),
  fFMD1(0),
  fFMD2(0),
  fFMD3(0),
  fPipe(0),
  fEarlyDecay(0),
  phihist(),
  fStored(0),
  fState(),
  fMaxConsequtiveStrips(2),
  fLowCutvalue(0),
  fTrackGammaToPi0(true)
  {
  //
  //  Constructor
  //
  //  Parameters:
  //   name: Name of task
  //
  DefineInput(1, AliForwardTaskValidation::Class());

    DefineOutput(1, TList::Class());
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

    TRandom fRandom = TRandom();              // random integer to use for creation of samples (used for error bars).
                                        // Needs to be created here, otherwise it will draw the same random number.

    fDeltaList = new TList();
    fDeltaList->SetName("Delta");
    Int_t phibins = 20;
    Int_t etabins = 50;
    fSettings.fnoSamples = 1;
    fSettings.fCentBins = 5;
    Int_t bins_phi_eta[5] = {fSettings.fnoSamples, fSettings.fNZvtxBins, phibins, fSettings.fCentBins, etabins} ;
    Double_t xmin_phi_eta[5] = {0,fSettings.fZVtxAcceptanceLowEdge, -TMath::Pi(), 0, -4};
    Double_t xmax_phi_eta[5] = {10,fSettings.fZVtxAcceptanceUpEdge, TMath::Pi(), 100, 6}; //

    Int_t bins_phi_eta_spd[5] = {fSettings.fnoSamples, fSettings.fNZvtxBins, phibins, 1, etabins} ;
    Double_t xmin_phi_eta_spd[5] = {0,fSettings.fZVtxAcceptanceLowEdge, -TMath::Pi(), 0, -4};
    Double_t xmax_phi_eta_spd[5] = {10,fSettings.fZVtxAcceptanceUpEdge, TMath::Pi(), 100, 6}; //

    Int_t dimensions = 5;

    Int_t bins_eta_phi[5] = {fSettings.fnoSamples, fSettings.fNZvtxBins, 51, 1, 20} ;
    Double_t xmin_eta_phi[5] = {0,fSettings.fZVtxAcceptanceLowEdge, -4, 0, 0.0};
    Double_t xmax_eta_phi[5] = {10,fSettings.fZVtxAcceptanceUpEdge, 6, 100, 2*TMath::Pi()}; //

    Int_t bins_phi[5] = {fSettings.fnoSamples, fSettings.fNZvtxBins, 20, 1, 20} ;
    Double_t xmin_phi[5] = {0,fSettings.fZVtxAcceptanceLowEdge, -TMath::Pi(), 0.0, 0.0};
    Double_t xmax_phi[5] = {10,fSettings.fZVtxAcceptanceUpEdge, TMath::Pi(), 100, 2*TMath::Pi()}; //

    Int_t bins_eta[5] = {fSettings.fnoSamples, fSettings.fNZvtxBins, 51, 1, 50} ;
    Double_t xmin_eta[5] = {0,fSettings.fZVtxAcceptanceLowEdge, -4, 0, -4};
    Double_t xmax_eta[5] = {10,fSettings.fZVtxAcceptanceUpEdge, 6, 100, 6}; //

    fDeltaList->Add(new THnD("delta_phi_eta", "delta_phi_eta",dimensions,bins_phi_eta, xmin_phi_eta, xmax_phi_eta)); // (samples, vertex,phi_mother - phi_tr ,centrality,eta_mother,eta_tr,eta_p)
    fDeltaList->Add(new THnD("delta_phi_eta_spd", "delta_phi_eta_spd",dimensions,bins_phi_eta_spd, xmin_phi_eta_spd, xmax_phi_eta_spd)); // (samples, vertex,phi_mother - phi_tr ,centrality,eta_mother,eta_tr,eta_p)
    fDeltaList->Add(new THnD("delta_eta_phi", "delta_eta_phi",dimensions,bins_eta_phi, xmin_eta_phi, xmax_eta_phi));// (samples, vertex,eta_mother - eta_tr ,centrality,phi_mother,phi_tr,phi_p)
    fDeltaList->Add(new THnD("delta_phi", "delta_phi",dimensions,bins_phi, xmin_phi, xmax_phi)); // (samples, vertex,phi_mother - phi_tr ,centrality,eta_mother,eta_tr,eta_p)
    fDeltaList->Add(new THnD("delta_eta", "delta_eta",dimensions,bins_eta, xmin_eta, xmax_eta));// (samples, vertex,eta_mother - eta_tr ,centrality,phi_mother,phi_tr,phi_p)

    static_cast<THnD*>(fDeltaList->FindObject("delta_phi_eta"))->GetAxis(0)->SetName("samples");
    static_cast<THnD*>(fDeltaList->FindObject("delta_phi_eta"))->GetAxis(1)->SetName("vertex");
    static_cast<THnD*>(fDeltaList->FindObject("delta_phi_eta"))->GetAxis(2)->SetName("phi_mother - phi_tr");
    static_cast<THnD*>(fDeltaList->FindObject("delta_phi_eta"))->GetAxis(3)->SetName("centrality");
    static_cast<THnD*>(fDeltaList->FindObject("delta_phi_eta"))->GetAxis(4)->SetName("eta");

    static_cast<THnD*>(fDeltaList->FindObject("delta_phi_eta_spd"))->GetAxis(0)->SetName("samples");
    static_cast<THnD*>(fDeltaList->FindObject("delta_phi_eta_spd"))->GetAxis(1)->SetName("vertex");
    static_cast<THnD*>(fDeltaList->FindObject("delta_phi_eta_spd"))->GetAxis(2)->SetName("phi_mother - phi_tr");
    static_cast<THnD*>(fDeltaList->FindObject("delta_phi_eta_spd"))->GetAxis(3)->SetName("centrality");
    static_cast<THnD*>(fDeltaList->FindObject("delta_phi_eta_spd"))->GetAxis(4)->SetName("eta");

    static_cast<THnD*>(fDeltaList->FindObject("delta_eta_phi"))->GetAxis(0)->SetName("samples");
    static_cast<THnD*>(fDeltaList->FindObject("delta_eta_phi"))->GetAxis(1)->SetName("vertex");
    static_cast<THnD*>(fDeltaList->FindObject("delta_eta_phi"))->GetAxis(2)->SetName("eta_mother - eta_tr");
    static_cast<THnD*>(fDeltaList->FindObject("delta_eta_phi"))->GetAxis(3)->SetName("centrality");
    static_cast<THnD*>(fDeltaList->FindObject("delta_eta_phi"))->GetAxis(4)->SetName("phi");

    static_cast<THnD*>(fDeltaList->FindObject("delta_phi"))->GetAxis(0)->SetName("samples");
    static_cast<THnD*>(fDeltaList->FindObject("delta_phi"))->GetAxis(1)->SetName("vertex");
    static_cast<THnD*>(fDeltaList->FindObject("delta_phi"))->GetAxis(2)->SetName("phi_mother - phi_tr");
    static_cast<THnD*>(fDeltaList->FindObject("delta_phi"))->GetAxis(3)->SetName("centrality");
    static_cast<THnD*>(fDeltaList->FindObject("delta_phi"))->GetAxis(4)->SetName("phi");

    static_cast<THnD*>(fDeltaList->FindObject("delta_eta"))->GetAxis(0)->SetName("samples");
    static_cast<THnD*>(fDeltaList->FindObject("delta_eta"))->GetAxis(1)->SetName("vertex");
    static_cast<THnD*>(fDeltaList->FindObject("delta_eta"))->GetAxis(2)->SetName("eta_mother - eta_tr");
    static_cast<THnD*>(fDeltaList->FindObject("delta_eta"))->GetAxis(3)->SetName("centrality");
    static_cast<THnD*>(fDeltaList->FindObject("delta_eta"))->GetAxis(4)->SetName("eta");

      //static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(0))->FindObject(Form("cumuRef_v%d", n)))->GetAxis(0)->SetName("samples");

    fEventList = new TList();

    fEventList->Add(new TH1D("Centrality","Centrality",10,0,100));
    fEventList->Add(new TH1D("Vertex","Vertex",fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));
    fEventList->Add(new TH2F("hOutliers","Maximum #sigma from mean N_{ch} pr. bin",
     20, 0., 100., 500, 0., 5.)); //((fFlags & kMC) ? 15. : 5. // Sigma <M> histogram
    fEventList->Add(new TH1D("FMDHits","FMDHits",100,0,10));
    fEventList->SetName("EventInfo");
    phihist = new TH1D("name","name",20,0,2*TMath::Pi());
    phihist->SetDirectory(0);
    Int_t bins_prim[4] = {fSettings.fnoSamples, fSettings.fNZvtxBins, 1, etabins} ;
    Double_t xmin_prim[4] = {0,fSettings.fZVtxAcceptanceLowEdge, 0, -6};
    Double_t xmax_prim[4] = {10,fSettings.fZVtxAcceptanceUpEdge, 100, 6}; //
    Int_t dimensions_prim = 4;
    fDeltaList->Add(new THnD("fnoPrim", "fnoPrim", dimensions_prim, bins_prim, xmin_prim, xmax_prim)); //(samples,vertex, phi, cent, eta)
    static_cast<THnD*>(fDeltaList->FindObject("fnoPrim"))->GetAxis(0)->SetName("samples");
    static_cast<THnD*>(fDeltaList->FindObject("fnoPrim"))->GetAxis(1)->SetName("vertex");
    static_cast<THnD*>(fDeltaList->FindObject("fnoPrim"))->GetAxis(2)->SetName("centrality");
    static_cast<THnD*>(fDeltaList->FindObject("fnoPrim"))->GetAxis(3)->SetName("eta_mother");

    fOutputList->Add(fEventList);
    fOutputList->Add(fDeltaList);

    // do analysis
    PostData(1, fOutputList);
  }







//_____________________________________________________________________
void AliForwardSecondariesTask::UserExec(Option_t *)
{
  //
  //  Parameters:
  //   option: Not used
  //

  std::cout << "getting validation";
  AliForwardTaskValidation* ev_val = dynamic_cast<AliForwardTaskValidation*>(this->GetInputData(1));
   if (!ev_val->IsValidEvent()){
     PostData(1, this->fOutputList);
     return;
   }
   std::cout << "... done" << std::endl;


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
  Float_t event_vtx_z = fAOD->GetPrimaryVertex()->GetZ();
  if (!(TMath::Abs(event_vtx_z) > 0)) {
    return;
  }

  AliMultSelection *MultSelection = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");

  //Int_t nPrim     = stack->GetNprimary();

  Double_t randomInt = static_cast<Double_t>(fRandom.Integer(fSettings.fnoSamples));
  //Float_t event_vtx_z = fAOD->GetPrimaryVertex()->GetZ();
  Double_t v0cent = MultSelection->GetMultiplicityPercentile("V0M");
  bool useEvent = kTRUE;

  if (useEvent){
    static_cast<TH1D*>(fEventList->FindObject("Centrality"))->Fill(v0cent);
    static_cast<TH1D*>(fEventList->FindObject("Vertex"))->Fill(event_vtx_z);
    std::vector< Int_t > listOfMothers;
    Int_t nTracks   = fAOD->GetNumberOfTracks();// stack->GetNtrack();

    Int_t nPrim     = fAOD->GetNumberOfPrimaries();//fAOD->GetNumberOfPrimaries();
    for (Int_t iTr = 0; iTr < nTracks; iTr++) {
      AliMCParticle* particle =
        static_cast<AliMCParticle*>(fAOD->GetTrack(iTr));

      // Check if this charged and a primary
      if (particle->Charge() == 0) continue;

      Bool_t isPrimary = stack->IsPhysicalPrimary(iTr) && iTr < nPrim;

      AliMCParticle* mother = isPrimary ? particle : GetMother(iTr,fAOD);
      if (!mother) mother = particle;
      // IF the track corresponds to a primary, pass that as both
      // arguments.
      ProcessTrack(particle, mother,listOfMothers, randomInt,event_vtx_z,v0cent);

    } // Loop over tracks

    PostData(1, fOutputList);
  } // End of useEvent
  return;
}


Bool_t
AliForwardSecondariesTask::ProcessTrack(AliMCParticle* particle,
				    AliMCParticle* mother, std::vector< Int_t > listOfMothers, Double_t randomInt, Float_t event_vtx_z, Double_t v0cent)
{
  THnD* delta_phi_eta_spd = static_cast<THnD*>(fDeltaList->FindObject("delta_phi_eta_spd")); // (samples, vertex,phi_mother - phi_tr ,centrality,eta_mother,eta_tr,eta_p)

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

    // Check that we hit an Base element
    if (AliTrackReference::kITS == ref->DetectorId()) {

      // We are interested if it produced a signal, not only a hit in the support structure.
      // This is an envelop around the active area
      if (ref->R() > 3.5 && ref->R() < 4.5 && TMath::Abs(ref->Z()) < 14.1) {
        if (!fStored){
          fStored = ref;
        Double_t phi_mother_spd = (mother->Phi());//Wrap02pi

        Double_t *etaPhi_spd = new Double_t[2];
        this->GetTrackRefEtaPhi(ref, etaPhi_spd);

        Double_t phi_tr_spd = etaPhi_spd[1]; //Wrap02pi
        Double_t eta_tr_spd = etaPhi_spd[0];
        Double_t phi[5] = {randomInt,event_vtx_z, WrapPi(phi_mother_spd - phi_tr_spd), v0cent, eta_tr_spd};//Wrap02pi

        delta_phi_eta_spd->Fill(phi,1);
      }
    }
  }




    if (ref->DetectorId() != AliTrackReference::kFMD) continue;

    AliTrackReference* test = ProcessRef(particle, mother, ref,listOfMothers,  randomInt,  event_vtx_z,  v0cent);
    if (test) store = test;

  } // Loop over track references
  if (!store) return true; // Nothing found

  StoreParticle(particle, mother, store, listOfMothers, randomInt,  event_vtx_z,  v0cent);
  EndTrackRefs();

  return true;
}

//____________________________________________________________________
void
AliForwardSecondariesTask::StoreParticle(AliMCParticle*       particle,
				     AliMCParticle* mother,
				     AliTrackReference*   ref,std::vector< Int_t > listOfMothers, Double_t randomInt, Float_t event_vtx_z, Double_t v0cent)
{


  THnD* delta_eta = static_cast<THnD*>(fDeltaList->FindObject("delta_eta")); // (samples, vertex,phi_mother - phi_tr ,centrality,eta_mother,eta_tr,eta_p)
  THnD* delta_phi = static_cast<THnD*>(fDeltaList->FindObject("delta_phi")); // (samples, vertex,phi_mother - phi_tr ,centrality,eta_mother,eta_tr,eta_p)
  THnD* delta_phi_eta = static_cast<THnD*>(fDeltaList->FindObject("delta_phi_eta")); // (samples, vertex,phi_mother - phi_tr ,centrality,eta_mother,eta_tr,eta_p)
  THnD* delta_eta_phi = static_cast<THnD*>(fDeltaList->FindObject("delta_eta_phi")); // (samples, vertex,phi_mother - phi_tr ,centrality,eta_mother,eta_tr,eta_p)
  THnD* fnoPrim = static_cast<THnD*>(fDeltaList->FindObject("fnoPrim"));//->Fill(event_vtx_z,event_vtx_z,event_vtx_z);



    UInt_t packed = ref->UserId();
    UShort_t detector, sector, strip;
    Char_t   ring;
    AliFMDStripIndex::Unpack(packed,detector,ring,sector,strip);
    TString inner = "I";
    //const Char_t* outer = "O";
    if (detector == 1){
      v0cent = 10;
    }

    if (detector == 2){
      if (TString(ring) == inner) v0cent = 30;
      else v0cent = 50;
    }

    if (detector == 3){
      if (ring == inner) v0cent = 70;
      else v0cent = 90;
    }

std::cout << "Detector: " << detector << ", Sector:" << sector << ", Strip: " << strip << ", Ring: " << ring << std::endl;

  Double_t eta_mother = mother->Eta();
  Double_t phi_mother = (mother->Phi());//Wrap02pi

  Double_t *etaPhi = new Double_t[2];
  this->GetTrackRefEtaPhi(ref, etaPhi);


  Double_t phi_tr = etaPhi[1]; //Wrap02pi
  //phi_tr = phihist->GetBinCenter(phihist->FindBin(phi_tr));

  //if (phi_tr < 0) phi_tr += 2*TMath::Pi();
  Double_t eta_tr = etaPhi[0];

  // (samples, vertex,phi_mother - phi_tr ,centrality,eta_mother,eta_tr)
  Double_t phi[5] = {randomInt,event_vtx_z, WrapPi(phi_mother - phi_tr), v0cent, eta_tr};//Wrap02pi
  Double_t eta[5] = {randomInt,event_vtx_z, (eta_mother - eta_tr), v0cent, phi_tr};

  //if (!(fabs(eta_tr - eta_mother) < 0.1)) continue;
  delta_phi_eta->Fill(phi,1);
  delta_eta_phi->Fill(eta,1);

  phi[4] = phi_tr;
  eta[4] = eta_tr;

  delta_phi->Fill(phi,1);
  delta_eta->Fill(eta,1);

  Double_t x_prim[4] =  {randomInt,event_vtx_z,v0cent,eta_tr};
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
AliForwardSecondariesTask::ProcessRef(AliMCParticle*       particle,
				  AliMCParticle* mother,
				 AliTrackReference*   ref,std::vector< Int_t > listOfMothers, Double_t randomInt, Float_t event_vtx_z, Double_t v0cent)
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
    // Int_t nnT   = TMath::Abs(fState.oldStrip - fState.startStrip) + 1;
    StoreParticle(particle, mother, fState.longest,listOfMothers,  randomInt,  event_vtx_z,  v0cent);
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
  Bool_t isPP = event->Stack()->IsPhysicalPrimary(p->GetLabel());
  //Bool_t isPP = this->IsRedefinedPhysicalPrimary(p);
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

Bool_t AliForwardSecondariesTask::hasParticleMaterialInteractionInAncestors(AliMCParticle* p){
  AliMCEvent* event = this->MCEvent();

  Bool_t pIsFromMat = event->Stack()->IsSecondaryFromMaterial(p->GetLabel());
  // If `p` is from material, we don't need to look further
  if (pIsFromMat) {
    return true;
  }
  // `p` has no mother and is not from material interaction
  if ((p->GetMother() < 0)) {
    return false;
  }
  AliMCParticle* ancestor = dynamic_cast< AliMCParticle* >(event->GetTrack(p->GetMother()));
  // Return the ancestor if `p` is from material but `ancestor` is not
  // Recurse if non of the above patterns matched
  return hasParticleMaterialInteractionInAncestors(ancestor);
}


AliMCParticle* AliForwardSecondariesTask::GetChargedMother(AliMCParticle* p) {
  AliMCParticle *mother = this->GetMother(p);
  if (!mother || mother->Charge() == 0) {
    return 0x0;
  }
  return mother;
}

AliMCParticle* AliForwardSecondariesTask::GetNeutralMother(AliMCParticle* p) {
  AliMCParticle *mother = this->GetMother(p);
  if (!mother || mother->Charge() != 0) {
    return 0x0;
  }
  return mother;
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


AliMCParticle* AliForwardSecondariesTask::GetIncidentParticleFromFirstMaterialInteraction(AliMCParticle* p) {
  AliMCEvent* event = this->MCEvent();
  // Is this particle from material, but there is no mother?!
  Bool_t pIsFromMat = event->Stack()->IsSecondaryFromMaterial(p->GetLabel());
  // If `p` is not from material, we don't need to look further
  if (!pIsFromMat) {
    return 0x0;
  }
  // `p` is from material, but has no mother; This should not happen
  if (pIsFromMat && (p->GetMother() < 0)) {
    return 0x0;
  }
  AliMCParticle* ancestor = dynamic_cast< AliMCParticle* >(event->GetTrack(p->GetMother()));
  // Return the ancestor if `p` is from material but `ancestor` is not
  if (pIsFromMat &&
      !event->Stack()->IsSecondaryFromMaterial(ancestor->GetLabel())) {
    return ancestor;
  }
  // Recurse if non of the above patterns matched
  else {
    return GetIncidentParticleFromFirstMaterialInteraction(ancestor);
  }
}

AliMCParticle* AliForwardSecondariesTask::GetFirstNonPrimaryMother(AliMCParticle* p) {
  AliMCEvent* event = this->MCEvent();
  // If `p` is not from material, we don't need to look further
  if (event->Stack()->IsPhysicalPrimary(p->GetLabel())) {
    return 0x0;
  }
  // Are there no more mothers left?
  if (p->GetMother() < 0) {
    return 0x0;
  }
  AliMCParticle* ancestor = dynamic_cast< AliMCParticle* >(event->GetTrack(p->GetMother()));
  // Return p if its ancestor is a primary particle; else recurse
  if (event->Stack()->IsPhysicalPrimary(ancestor->GetLabel())) {
    return p;
  }
  // Recurse if non of the above patterns matched
  else {
    return GetIncidentParticleFromFirstMaterialInteraction(ancestor);
  }
}




AliTrackReference* AliForwardSecondariesTask::IsHitITS(AliMCParticle* p) {
  for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) {
    AliTrackReference* ref = p->GetTrackReference(iTrRef);
    // Check hit on ITS
    if (!ref || AliTrackReference::kITS != ref->DetectorId()) {
      continue;
    }
    // We are interested if it produced a signal, not only a hit in the support structure.
    // This is an envelop around the active area
    if (ref->R() > 3.5 && ref->R() < 4.5 && TMath::Abs(ref->Z()) < 14.1) {
      return ref;
    }
  }
  return 0x0;
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

void AliForwardSecondariesTask::GetTrackRefEtaPhi(AliMCParticle* p, Double_t* etaPhi) {
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

Int_t AliForwardSecondariesTask::GetOriginType(AliMCParticle *p) {
  AliMCEvent* event = this->MCEvent();
  AliStack *stack = event->Stack();
  if (stack->IsPhysicalPrimary(p->GetLabel())) {
    return cOriginType::kPRIMARY;
  }
  Double_t r = TMath::Sqrt(p->Yv() * p->Yv() + p->Xv() * p->Xv());
  if (this->fITS->IsInside(p->Zv(), r)) {
    return cOriginType::kITS;
  }
  if (this->fFMD1->IsInside(p->Zv(), r) ||
      this->fFMD2->IsInside(p->Zv(), r) ||
      this->fFMD3->IsInside(p->Zv(), r)) {
    return cOriginType::kFMD;
  }
  if (this->fPipe->IsInside(p->Zv(), r)) {
    return cOriginType::kPIPE;
  }
  if (this->fEarlyDecay->IsInside(p->Zv(), r)) {
    return cOriginType::kEARLYDECAY;
  }
  return cOriginType::kOTHER;
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
