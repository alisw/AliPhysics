#include "AliBaseMCTrackDensity.h"
#include "AliForwardFlowWeights.h"
#include <AliMCEvent.h>
#include <AliTrackReference.h>
//#include <AliStack.h>
#include <TMath.h>
#include <AliLog.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TList.h>
#include <TROOT.h>
#include <TVector3.h>
#include <iostream>
#include "AliCollisionGeometry.h"
#include "AliGenEventHeader.h"
#include "AliForwardUtil.h"
#include <TF1.h>
#include <TGraph.h>

//____________________________________________________________________
AliBaseMCTrackDensity::AliBaseMCTrackDensity()
  : TNamed(), 
    fUseOnlyPrimary(false), 
    fBinFlow(0), 
    fEtaBinFlow(0),
    fPhiBinFlow(0),
    fNRefs(0),
    fWeights(0),
    fTruthWeights(0),
    fIP(0,0,0), 
    fB(0),
    fPhiR(0),
    fDebug(false),
    fTrackGammaToPi0(false)
{
  // Default constructor 
  DGUARD(fDebug, 3,"Default CTOR of AliBasMCTrackDensity");
}

//____________________________________________________________________
AliBaseMCTrackDensity::AliBaseMCTrackDensity(const char* name)
  : TNamed(name,"mcTrackDensity"), 
    fUseOnlyPrimary(false), 
    fBinFlow(0), 
    fEtaBinFlow(0),
    fPhiBinFlow(0),
    fNRefs(0),
    fWeights(0),
    fTruthWeights(0),
    fIP(0,0,0), 
    fB(0),
    fPhiR(0),
    fDebug(false),
    fTrackGammaToPi0(false)
{
  // Normal constructor constructor 
  DGUARD(fDebug, 3,"Named CTOR of AliBasMCTrackDensity: %s", name);
}

//____________________________________________________________________
AliBaseMCTrackDensity::AliBaseMCTrackDensity(const AliBaseMCTrackDensity& o)
  : TNamed(o),
    fUseOnlyPrimary(o.fUseOnlyPrimary), 
    fBinFlow(o.fBinFlow), 
    fEtaBinFlow(o.fEtaBinFlow),
    fPhiBinFlow(o.fPhiBinFlow),
    fNRefs(o.fNRefs),
    fWeights(o.fWeights),
    fTruthWeights(o.fTruthWeights),
    fIP(o.fIP), 
    fB(o.fB),
    fPhiR(o.fPhiR),
    fDebug(o.fDebug),
    fTrackGammaToPi0(o.fTrackGammaToPi0)
{
  // Normal constructor constructor 
  DGUARD(fDebug, 3,"Copy CTOR of AliBasMCTrackDensity");
}

//____________________________________________________________________
AliBaseMCTrackDensity&
AliBaseMCTrackDensity::operator=(const AliBaseMCTrackDensity& o)
{
  // Assignment operator 
  DGUARD(fDebug,3,"MC track density assignmetn");
  if (&o == this) return *this; 
  TNamed::operator=(o);
  fUseOnlyPrimary       = o.fUseOnlyPrimary;
  fBinFlow              = o.fBinFlow;
  fEtaBinFlow           = o.fEtaBinFlow;
  fPhiBinFlow           = o.fPhiBinFlow;
  fNRefs                = o.fNRefs;
  fDebug                = o.fDebug;
  fWeights              = o.fWeights;
  fTruthWeights         = o.fTruthWeights;
  fIP.SetXYZ(o.fIP.X(), o.fIP.Y(), o.fIP.Z());
  fB                    = o.fB;
  fPhiR                 = o.fPhiR;
  fTrackGammaToPi0      = o.fTrackGammaToPi0;
  return *this;
}

//____________________________________________________________________
void
AliBaseMCTrackDensity::SetWeights(AliBaseMCWeights* w)
{
  if (fWeights) {
    delete fWeights;
    fWeights = 0;
  }
  fWeights = w;
}
//____________________________________________________________________
void
AliBaseMCTrackDensity::SetTruthWeights(AliBaseMCWeights* w)
{
  if (fTruthWeights) {
    delete fTruthWeights;
    fTruthWeights = 0;
  }
  fTruthWeights = w;
}
//____________________________________________________________________
void
AliBaseMCTrackDensity::SetUseFlowWeights(Bool_t use)
{
  if (!use) return;

  SetWeights(new AliForwardFlowWeights());
}
//____________________________________________________________________
void
AliBaseMCTrackDensity::CreateOutputObjects(TList* l)
{
  DGUARD(fDebug,1,"MC track defines output");
  TList* ll = new TList;
  ll->SetName(GetTitle());
  ll->SetOwner();
  l->Add(ll);
  
  fBinFlow = new TH2D("binFlow", "#eta and #varphi bin flow", 
		      200, -5, 5, 40, -180, 180);
  fBinFlow->SetXTitle("#Delta#eta");
  fBinFlow->SetYTitle("#Delta#varphi");
  fBinFlow->SetOption("colz");
  fBinFlow->SetDirectory(0);
  ll->Add(fBinFlow);

  fEtaBinFlow = new TH2D("binFlowEta", "#eta bin flow vs #eta", 
			 200, -4, 6, 200, -5, 5);
  fEtaBinFlow->SetXTitle("#eta");
  fEtaBinFlow->SetYTitle("#Delta#eta");
  fEtaBinFlow->SetOption("colz");
  fEtaBinFlow->SetDirectory(0);
  ll->Add(fEtaBinFlow);

  fPhiBinFlow = new TH2D("binFlowPhi", "#phi bin flow vs #phi", 
			 40, 0, 360, 40, -180, 180);
  fPhiBinFlow->SetXTitle("#varphi");
  fPhiBinFlow->SetYTitle("#Delta#varphi");
  fPhiBinFlow->SetOption("colz");
  fPhiBinFlow->SetDirectory(0);
  ll->Add(fPhiBinFlow);

  fNRefs = new TH1D("nRefs", "# references per track", 21, -.5, 20.5);
  fNRefs->SetXTitle("# references");
  fNRefs->SetFillColor(kMagenta+1);
  fNRefs->SetFillStyle(3001);
  fNRefs->SetDirectory(0);
  ll->Add(fNRefs);

  if(fWeights)  fWeights->Init(ll);
}


//____________________________________________________________________
Double_t
AliBaseMCTrackDensity::StoreParticle(AliMCParticle*       particle, 
				     const AliMCParticle* mother, 
				     AliTrackReference*   ref) const
{
  // Store this particle
  //
  // Note: If particle refers to a primary, then particle and mother
  // refers to the same particle (the address are the same)
  // 
  DGUARD(fDebug,3,"MC track density store particle");
  // Store a particle. 
  if (!ref) return 0;

  Double_t weight = 1;
  if (fWeights) 
    weight = CalculateWeight(mother, (mother == particle));
  

  // Get track-reference stuff 
  Double_t x      = ref->X()-fIP.X();
  Double_t y      = ref->Y()-fIP.Y();
  Double_t z      = ref->Z()-fIP.Z();
  Double_t rr     = TMath::Sqrt(x*x+y*y);
  Double_t thetaR = TMath::ATan2(rr,z);
  Double_t phiR   = TMath::ATan2(y,x);
  Double_t etaR   = -TMath::Log(TMath::Tan(thetaR/2));
  
  // Correct angle and convert to degrees 
  if (thetaR < 0) thetaR += 2*TMath::Pi();
  thetaR *= 180. / TMath::Pi();
  if (phiR < 0) phiR += 2*TMath::Pi();
  phiR *= 180. / TMath::Pi();

  const AliMCParticle* mp = (mother ? mother : particle);
  Double_t dEta = mp->Eta() - etaR;
  Double_t dPhi = mp->Phi() * 180 / TMath::Pi() - phiR;
  if (dPhi >  180) dPhi -= 360;
  if (dPhi < -180) dPhi += 360;
  fBinFlow->Fill(dEta, dPhi);
  fEtaBinFlow->Fill(etaR, dEta);
  fPhiBinFlow->Fill(phiR, dPhi);
  return weight;
}

//____________________________________________________________________
Double_t
AliBaseMCTrackDensity::GetTrackRefTheta(const AliTrackReference* ref) const
{
  // Get the incidient angle of the track reference. 
  Double_t x    = ref->X()-fIP.X();
  Double_t y    = ref->Y()-fIP.Y();
  Double_t z    = ref->Z()-fIP.Z();
  Double_t rr   = TMath::Sqrt(x*x+y*y);
  Double_t theta= TMath::ATan2(rr,z);
  Double_t ang  = TMath::Abs(TMath::Pi()-theta);
  return ang;
}
				    
//____________________________________________________________________
const AliMCParticle*
AliBaseMCTrackDensity::GetMother(Int_t iTr, const AliMCEvent& event) const
{
  // 
  // Track down primary mother 
  // 
  Int_t                i         = iTr;
  Bool_t               gammaSeen = false;
  const AliMCParticle* candidate = 0;
  do { 
    const AliMCParticle* p = static_cast<AliMCParticle*>(event.GetTrack(i));
    if (!p) break;
    if (gammaSeen && TMath::Abs(p->PdgCode()) == 111) 
      // If we're looking for a mother pi0 of gamma, and we find it
      // here, we return it - irrespective of whether it's flagged as
      // a primary or not.
      return p;

    if (const_cast<AliMCEvent&>(event).IsPhysicalPrimary(i)) {
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

//____________________________________________________________________
Bool_t
AliBaseMCTrackDensity::GetCollisionParameters(const AliMCEvent& event)
{ 
  DGUARD(fDebug,3,"MC track density get collision parameters");
  AliCollisionGeometry* hd = 
    dynamic_cast<AliCollisionGeometry*>(event.GenEventHeader());
  fPhiR = (hd ? hd->ReactionPlaneAngle() : 0.);
  fB    = (hd ? hd->ImpactParameter() : -1 );
  return hd != 0;
}

//____________________________________________________________________
Bool_t
AliBaseMCTrackDensity::ProcessTrack(AliMCParticle* particle, 
				    const AliMCParticle* mother)
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
  Int_t nRefs = 0;
  for (Int_t iTrRef = 0; iTrRef < nTrRef; iTrRef++) { 
    AliTrackReference* ref = particle->GetTrackReference(iTrRef);
      
    // Check existence 
    if (!ref) continue;

    // Check that we hit an Base element 
    if (ref->DetectorId() != GetDetectorId()) continue;
    if (!CheckTrackRef(ref)) continue;

    nRefs++;

    AliTrackReference* test = ProcessRef(particle, mother, ref);
    if (test) store = test;

  } // Loop over track references
  if (!store) return true; // Nothing found
  
  StoreParticle(particle, mother, store);

  fNRefs->Fill(nRefs);

  EndTrackRefs(nRefs);

  return true;
}


//____________________________________________________________________
Bool_t
AliBaseMCTrackDensity::ProcessTracks(const AliMCEvent& event,
				     const TVector3&   ip, 
				     TH2D*             primary)
{
  // 
  // Filter the input kinematics and track references, using 
  // some of the ESD information
  // 
  // Parameters:
  //    input   Input ESD event
  //    event   Input MC event
  //    vz      Vertex position 
  //    output  Output ESD-like object
  //    primary Per-event histogram of primaries 
  //
  // Return:
  //    True on succes, false otherwise 
  //
  DGUARD(fDebug,3,"MC track density Process a tracks");
  fIP.SetXYZ(ip.X(), ip.Y(), ip.Z());
  GetCollisionParameters(event);
  
  //  AliStack* stack = const_cast<AliMCEvent&>(event).Stack();
  //  if (!stack) return kFALSE;

  
  Int_t nTracks   = /*stack->GetNtrack();*/event.GetNumberOfTracks();
  // temporary remove constness, since event.GetNumberOfPrimaries() was not const
  Int_t nPrim     = /*stack->GetNprimary();*/((AliMCEvent&)event).GetNumberOfPrimaries();
  for (Int_t iTr = 0; iTr < nTracks; iTr++) { 
    AliMCParticle* particle = 
      static_cast<AliMCParticle*>(event.GetTrack(iTr));
    
    // Check if this charged and a primary 
    if (!particle->Charge() != 0) continue;
    
    Bool_t isPrimary = event.IsPhysicalPrimary(iTr) && iTr < nPrim;
    
    // Fill 'dn/deta' histogram 
    if (isPrimary && primary) {
      Double_t w = 1;
      if (fTruthWeights) w = CalculateTruthWeight(particle);
      primary->Fill(particle->Eta(), particle->Phi(), w);
    }

    // Bail out if we're only processing primaries - perhaps we should
    // track back to the original primary?
    if (fUseOnlyPrimary && !isPrimary) continue;

    const AliMCParticle* mother = isPrimary ? particle : GetMother(iTr, event);
    // IF the track corresponds to a primary, pass that as both
    // arguments.
    ProcessTrack(particle, mother);

  } // Loop over tracks
  return kTRUE;
}

  
//____________________________________________________________________
Double_t
AliBaseMCTrackDensity::CalculateWeight(const AliMCParticle* p,
				       Bool_t isPrimary) const
{
  // Note, it is always the ultimate mother that is passed here.
  if (!p) return 0;
  return fWeights->CalcWeight(p, isPrimary, fPhiR, fB);
}
//____________________________________________________________________
Double_t
AliBaseMCTrackDensity::CalculateTruthWeight(const AliMCParticle* p) const
{
  // Note, it is always the ultimate mother that is passed here.
  return fTruthWeights->CalcWeight(p, true, fPhiR, fB);
}

#define PF(N,V,...)					\
  AliForwardUtil::PrintField(N,V, ## __VA_ARGS__)
#define PFB(N,FLAG)				\
  do {									\
    AliForwardUtil::PrintName(N);					\
    std::cout << std::boolalpha << (FLAG) << std::noboolalpha << std::endl; \
  } while(false)
#define PFV(N,VALUE)					\
  do {							\
    AliForwardUtil::PrintName(N);			\
    std::cout << (VALUE) << std::endl; } while(false)

//____________________________________________________________________
void
AliBaseMCTrackDensity::Print(Option_t* option) const 
{
  AliForwardUtil::PrintTask(*this);
  gROOT->IncreaseDirLevel();  
  PFB("Only primary tracks", fUseOnlyPrimary);
  PFB("Use weights", fWeights);
  if (fWeights) fWeights->Print(option);
  gROOT->DecreaseDirLevel();  
}

//____________________________________________________________________
//
// EOF
//
