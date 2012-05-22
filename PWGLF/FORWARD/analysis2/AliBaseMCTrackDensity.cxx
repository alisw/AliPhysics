#include "AliBaseMCTrackDensity.h"
#include <AliMCEvent.h>
#include <AliTrackReference.h>
#include <AliStack.h>
#include <TMath.h>
#include <AliLog.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TList.h>
#include <TROOT.h>
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
    fUseFlowWeights(false),
    fBinFlow(0), 
    fEtaBinFlow(0),
    fPhiBinFlow(0),
    fNRefs(0),
    fV2Eta(0), 
    fV22Pt(0), 
    fV24Pt(0), 
    fV2B(0),
    fVz(0), 
    fB(0),
    fPhiR(0),
    fDebug(false)
{
  // Default constructor 
  DGUARD(0,0,"MC track density default construction");
}

//____________________________________________________________________
AliBaseMCTrackDensity::AliBaseMCTrackDensity(const char* name)
  : TNamed(name,"mcTrackDensity"), 
    fUseOnlyPrimary(false), 
    fUseFlowWeights(false),
    fBinFlow(0), 
    fEtaBinFlow(0),
    fPhiBinFlow(0),
    fNRefs(0),
    fV2Eta(0), 
    fV22Pt(0), 
    fV24Pt(0), 
    fV2B(0), 
    fVz(0), 
    fB(0),
    fPhiR(0),
    fDebug(false)
{
  // Normal constructor constructor 
  DGUARD(0,0,"MC track density named construction: %s", name);
}

//____________________________________________________________________
AliBaseMCTrackDensity::AliBaseMCTrackDensity(const AliBaseMCTrackDensity& o)
  : TNamed(o),
    fUseOnlyPrimary(o.fUseOnlyPrimary), 
    fUseFlowWeights(o.fUseFlowWeights),
    fBinFlow(o.fBinFlow), 
    fEtaBinFlow(o.fEtaBinFlow),
    fPhiBinFlow(o.fPhiBinFlow),
    fNRefs(o.fNRefs),
    fV2Eta(o.fV2Eta), 
    fV22Pt(o.fV22Pt), 
    fV24Pt(o.fV24Pt), 
    fV2B(o.fV2B), 
    fVz(o.fVz), 
    fB(o.fB),
    fPhiR(o.fPhiR),
    fDebug(o.fDebug)
{
  // Normal constructor constructor 
  DGUARD(0,0,"MC track density copy construction");
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
  fUseFlowWeights       = o.fUseFlowWeights;
  fV2Eta                = o.fV2Eta;
  fV22Pt                = o.fV22Pt;
  fV24Pt                = o.fV24Pt;
  fV2B                  = o.fV2B;
  fVz                   = o.fVz;
  fB                    = o.fB;
  fPhiR                 = o.fPhiR;
  return *this;
}

//____________________________________________________________________
void
AliBaseMCTrackDensity::DefineOutput(TList* l)
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

  SetupWeights(ll);
}

//____________________________________________________________________
void 
AliBaseMCTrackDensity::SetupWeights(TList* l)
{
  DGUARD(fDebug,2,"MC track setup phi weights");
  fV2Eta = new TF1("v2eta", "gaus", -6, 6);
  fV2Eta->SetParameters(20 * 0.1, 0, 9);
  l->Add(fV2Eta);

  Int_t          ptN     = 19;
  const Double_t ptX[]   = {0.00,     0.25,     0.350,    0.45, 
			    0.55,     0.650,    0.75,     0.85, 
			    0.950,    1.10,     1.30,     1.500,
			    1.70,     1.90,     2.250,    2.75, 
			    3.25,     3.750,    4.50};
  { 
    // v2{2} dependence on pt
    const Double_t y[] = {0.00000,  0.043400, 0.059911, 0.073516,
			  0.089756, 0.105486, 0.117391, 0.128199,
			  0.138013, 0.158271, 0.177726, 0.196383,
			  0.208277, 0.216648, 0.242954, 0.249961,
			  0.240131, 0.269006, 0.207796};
    
    fV22Pt = new TGraph(ptN, ptX, y);
    fV22Pt->SetName("v2_2_vs_pt");
    fV22Pt->SetMarkerStyle(20);
    fV22Pt->SetMarkerColor(kRed+1);
    l->Add(fV22Pt);
  }

  {
    const Double_t y[] = {0.000000, 0.038646, 0.049824, 0.066662,
			  0.075856, 0.081583, 0.099778, 0.104674,
			  0.118545, 0.131874, 0.152959, 0.155348,
			  0.169751, 0.179052, 0.178532, 0.198851,
			  0.185737, 0.239901, 0.186098};

    // v2{4} dependence on pt 
    fV24Pt = new TGraph(ptN, ptX, y);
    fV24Pt->SetName("v2_4_vs_pt");
    fV24Pt->SetMarkerStyle(20);
    fV24Pt->SetMarkerColor(kBlue+1);
    l->Add(fV24Pt);
  }
  {
    // V2 dependence on centrality
    Int_t n            = 8;
    const Double_t x[] = {1.75,     4.225,    5.965,    7.765,
			  9.215,    10.46,    11.565,   12.575};
    const Double_t y[] = {0.017855, 0.032440, 0.055818, 0.073137,
			  0.083898, 0.086690, 0.082040, 0.077777};
    fV2B = new TGraph(n, x, y);
    fV2B->SetName("v2_vs_b");
    fV2B->SetMarkerStyle(20);
    fV2B->SetMarkerColor(kGreen+1);
    l->Add(fV2B);
  }
}


//____________________________________________________________________
Double_t
AliBaseMCTrackDensity::StoreParticle(AliMCParticle*       particle, 
				     const AliMCParticle* mother, 
				     AliTrackReference*   ref) const
{
  DGUARD(fDebug,3,"MC track density store particle");
  // Store a particle. 
  if (!ref) return 0;

  Double_t weight = 1;
  if (fUseFlowWeights) {
    Double_t phi = (mother ? mother->Phi() : particle->Phi());
    Double_t eta = (mother ? mother->Eta() : particle->Eta());
    Double_t pt  = (mother ? mother->Pt() : particle->Pt());
    Int_t    id  = (mother ? mother->PdgCode() : 2212);
    weight       = CalculateWeight(eta, pt, phi, id);
  }

  // Get track-reference stuff 
  Double_t x      = ref->X();
  Double_t y      = ref->Y();
  Double_t z      = ref->Z()-fVz;
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
  Double_t x    = ref->X();
  Double_t y    = ref->Y();
  Double_t z    = ref->Z()-fVz;
  Double_t rr   = TMath::Sqrt(x*x+y*y);
  Double_t theta= TMath::ATan2(rr,z);
  Double_t ang  = TMath::Abs(TMath::Pi()-theta);
  return ang;
}
				    
//____________________________________________________________________
const AliMCParticle*
AliBaseMCTrackDensity::GetMother(Int_t     iTr,
				const AliMCEvent& event) const
{
  // 
  // Track down primary mother 
  // 
  Int_t i  = iTr;
  do { 
    const AliMCParticle* p = static_cast<AliMCParticle*>(event.GetTrack(i));
    if (const_cast<AliMCEvent&>(event).Stack()->IsPhysicalPrimary(i)) return p;
    
    i = p->GetMother();
  } while (i > 0);

  return 0;
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
				     Double_t          vz,
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
  fVz = vz;
  GetCollisionParameters(event);
  
  AliStack* stack = const_cast<AliMCEvent&>(event).Stack();
  Int_t nTracks   = stack->GetNtrack();//event.GetNumberOfTracks();
  Int_t nPrim     = stack->GetNprimary();//event.GetNumberOfPrimary();
  for (Int_t iTr = 0; iTr < nTracks; iTr++) { 
    AliMCParticle* particle = 
      static_cast<AliMCParticle*>(event.GetTrack(iTr));
    
    // Check if this charged and a primary 
    if (!particle->Charge() != 0) continue;
    
    Bool_t isPrimary = stack->IsPhysicalPrimary(iTr) && iTr < nPrim;
    
    // Fill 'dn/deta' histogram 
    if (isPrimary && primary) 
      primary->Fill(particle->Eta(), particle->Phi());

    // Bail out if we're only processing primaries - perhaps we should
    // track back to the original primary?
    if (fUseOnlyPrimary && !isPrimary) continue;

    const AliMCParticle* mother = GetMother(iTr, event);
    ProcessTrack(particle, mother);

  } // Loop over tracks
  return kTRUE;
}

  
//____________________________________________________________________
Double_t
AliBaseMCTrackDensity::CalculateWeight(Double_t eta, Double_t pt, 
				       Double_t phi, Int_t id) const
{
  Double_t w = (fV2Eta->Eval(eta) * 0.5 * (fV22Pt->Eval(pt) + fV24Pt->Eval(pt)) 
		* fV2B->Eval(fB) / fV2B->Eval(10.46) 
		* 2 * TMath::Cos(2* (phi - fPhiR)));
  UInt_t aid = TMath::Abs(id);
  switch (aid) { 
  case 211:  w *= 1.3; break; // pions 
  case 2212: w *= 1.0; break; // protons 
  default:   w *= 0.7; break;
  }
  return w;
}
//____________________________________________________________________
void
AliBaseMCTrackDensity::Print(Option_t* /*option*/) const 
{
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << ind << ClassName() << ": " << GetName() << '\n'
	    << std::boolalpha 
	    << ind << " Only primary tracks:    " << fUseOnlyPrimary << '\n'
	    << ind << " Use flow after burner:  " << fUseFlowWeights 
	    << std::noboolalpha << std::endl;
  
}

//____________________________________________________________________
//
// EOF
//
