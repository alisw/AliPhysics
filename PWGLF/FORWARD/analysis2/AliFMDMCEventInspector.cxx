// 
// This class inspects the event 
//
// Input:
//   - AliESDFMD object possibly corrected for sharing
//
// Output:
//   - A histogram of v_z of events with triggers. 
//   - A histogram of v_z of events with vertex and triggers 
//   - A histogram of trigger counters 
// 
// Note, that these are added to the master output list 
//
// Corrections used: 
//   - None
//
#include "AliFMDMCEventInspector.h"
#include "AliLog.h"
#include "AliAODForwardMult.h"
#include "AliForwardUtil.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenHijingEventHeader.h"
// #include "AliGenHydjetEventHeader.h"
#include "AliGenEposEventHeader.h"
#include "AliGenHerwigEventHeader.h"
#include "AliGenGeVSimEventHeader.h"
#include "AliHeader.h"
#include <TList.h>
#include <TH2F.h>
#include <TParticle.h>
#include <TMath.h>

//====================================================================
AliFMDMCEventInspector::AliFMDMCEventInspector()
  : AliFMDEventInspector(), 
    fHVertex(0),
    fHPhiR(0), 
    fHB(0),
    fHBvsPart(0),
    fHBvsBin(0),
    fHBvsCent(0),
    fHVzComp(0),
    fHCentVsPart(0),
    fHCentVsBin(0),
    fProduction("")
{
  // 
  // Constructor 
  //
}

//____________________________________________________________________
AliFMDMCEventInspector::AliFMDMCEventInspector(const char* /* name */)
  : AliFMDEventInspector("fmdEventInspector"), 
    fHVertex(0),
    fHPhiR(0), 
    fHB(0),
    fHBvsPart(0),
    fHBvsBin(0),
    fHBvsCent(0),
    fHVzComp(0),
    fHCentVsPart(0),
    fHCentVsBin(0),
    fProduction("")
{
  // 
  // Constructor 
  // 
  // Parameters:
  //   name Name of object
  //
}

//____________________________________________________________________
AliFMDMCEventInspector::AliFMDMCEventInspector(const AliFMDMCEventInspector& o)
  : AliFMDEventInspector(o), 
    fHVertex(0),
    fHPhiR(0), 
    fHB(0),
    fHBvsPart(0),
    fHBvsBin(0),
    fHBvsCent(0),
    fHVzComp(0),
    fHCentVsPart(0),
    fHCentVsBin(0),
    fProduction("")
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //   o Object to copy from 
  //
}

//____________________________________________________________________
AliFMDMCEventInspector::~AliFMDMCEventInspector()
{
  // 
  // Destructor 
  //
}
//____________________________________________________________________
AliFMDMCEventInspector&
AliFMDMCEventInspector::operator=(const AliFMDMCEventInspector& o)
{
  // 
  // Assignement operator
  // 
  // Parameters:
  //   o Object to assign from 
  // 
  // Return:
  //    Reference to this object
  //
  AliFMDEventInspector::operator=(o);
  return *this;
}

//____________________________________________________________________
void
AliFMDMCEventInspector::Init(const TAxis& vtxAxis)
{
  // 
  // Initialize the object 
  // 
  // Parameters:
  //   vtxAxis Vertex axis in use 
  //
  AliFMDEventInspector::Init(vtxAxis);

  Int_t    maxPart = 450;
  Int_t    maxBin  = 225;
  Int_t    maxB    = 25;
  Int_t    nVtx = vtxAxis.GetNbins();
  Double_t lVtx = vtxAxis.GetXmin();
  Double_t hVtx = vtxAxis.GetXmax();
  fHVertex = new TH1F("vertex", "True vertex distribution", nVtx, lVtx, hVtx);
  fHVertex->SetXTitle("v_{z} [cm]");
  fHVertex->SetYTitle("# of events");
  fHVertex->SetFillColor(kGreen+1);
  fHVertex->SetFillStyle(3001);
  fHVertex->SetDirectory(0);
  // fHVertex->Sumw2();
  fList->Add(fHVertex);

  fHPhiR = new TH1F("phiR", "Event plane", 120, 0, 2*TMath::Pi());
  fHPhiR->SetXTitle("#Phi_{R} [radians]");
  fHPhiR->SetYTitle("# of events");
  fHPhiR->SetFillColor(kGreen+1);
  fHPhiR->SetFillStyle(3001);
  fHPhiR->SetDirectory(0);
  fList->Add(fHPhiR);

  fHB = new TH1F("b", "Impact parameter", 5*maxB, 0, maxB);
  fHB->SetXTitle("b [fm]");
  fHB->SetYTitle("# of events");
  fHB->SetFillColor(kGreen+1);
  fHB->SetFillStyle(3001);
  fHB->SetDirectory(0);
  fList->Add(fHB);

  fHBvsPart = new TH2F("bVsParticipants", "Impact parameter vs Participants",
		       5*maxB, 0, maxB, maxPart, -.5, maxPart-.5);
  fHBvsPart->SetXTitle("b [fm]");
  fHBvsPart->SetYTitle("# of participants");
  fHBvsPart->SetZTitle("Events");
  fHBvsPart->SetDirectory(0);
  fList->Add(fHBvsPart);

  fHBvsBin = new TH2F("bVsBinary", "Impact parameter vs Binary Collisions",
		       5*maxB, 0, maxB, maxBin, -.5, maxBin-.5);
  fHBvsBin->SetXTitle("b [fm]");
  fHBvsBin->SetYTitle("# of binary collisions");
  fHBvsBin->SetZTitle("Events");
  fHBvsBin->SetDirectory(0);
  fList->Add(fHBvsBin);
  
  fHBvsCent = new TH2F("bVsCentrality", "Impact parameter vs Centrality",
		       5*maxB, 0, maxB, fCentAxis->GetNbins(), 
		       fCentAxis->GetXbins()->GetArray());
  fHBvsCent->SetXTitle("b [fm]");
  fHBvsCent->SetYTitle("Centrality [%]");
  fHBvsCent->SetZTitle("Event");
  fHBvsCent->SetDirectory(0);
  fList->Add(fHBvsCent);
  
  
  fHVzComp = new TH2F("vzComparison", "v_{z} truth vs reconstructed",
		      10*nVtx, lVtx, hVtx, 10*nVtx, lVtx, hVtx);
  fHVzComp->SetXTitle("True v_{z} [cm]");
  fHVzComp->SetYTitle("Reconstructed v_{z} [cm]");
  fHVzComp->SetZTitle("Events");
  fHVzComp->SetDirectory(0);
  fList->Add(fHVzComp);

  fHCentVsPart = new TH2F("centralityVsParticipans", 
			  "# of participants vs Centrality",
			  maxPart, -.5, maxPart-.5, fCentAxis->GetNbins(), 
			  fCentAxis->GetXbins()->GetArray());
  fHCentVsPart->SetXTitle("Participants");
  fHCentVsPart->SetYTitle("Centrality [%]");
  fHCentVsPart->SetZTitle("Event");
  fHCentVsPart->SetDirectory(0);
  fList->Add(fHCentVsPart);

  fHCentVsBin = new TH2F("centralityVsBinary", 
			 "# of binary collisions vs Centrality",
			 maxBin, -.5, maxBin-.5, fCentAxis->GetNbins(), 
			 fCentAxis->GetXbins()->GetArray());
  fHCentVsBin->SetXTitle("Binary collisions");
  fHCentVsBin->SetYTitle("Centrality [%]");
  fHCentVsBin->SetZTitle("Event");
  fHCentVsBin->SetDirectory(0);
  fList->Add(fHCentVsBin);
}

//____________________________________________________________________
void
AliFMDMCEventInspector::StoreInformation(Int_t runNo)
{
  // Store information about running conditions in the output list 
  if (!fList) return;
  AliFMDEventInspector::StoreInformation(runNo);
  TNamed* mc = new TNamed("mc", fProduction.Data());
  mc->SetUniqueID(1);
  fList->Add(mc);
}

namespace
{
  TString& AppendField(TString& s, bool test, const char* f)
  {
    if (!test) return s;
    if (!s.IsNull()) s.Append(",");
    s.Append(f);
    return s;
  }
}

//____________________________________________________________________
void
AliFMDMCEventInspector::ReadProductionDetails(AliMCEvent* event)
{
  //Assign MC only triggers : True NSD etc.
  AliHeader*               header          = event->Header();
  AliGenEventHeader*       genHeader       = header->GenEventHeader();
  AliCollisionGeometry*    colGeometry     = 
    dynamic_cast<AliCollisionGeometry*>(genHeader);
  AliGenPythiaEventHeader* pythiaHeader    = 
    dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
  AliGenDPMjetEventHeader* dpmHeader       = 
    dynamic_cast<AliGenDPMjetEventHeader*>(genHeader);
  AliGenGeVSimEventHeader* gevHeader       = 
    dynamic_cast<AliGenGeVSimEventHeader*>(genHeader);
  AliGenHerwigEventHeader* herwigHeader    = 
    dynamic_cast<AliGenHerwigEventHeader*>(genHeader);
  AliGenHijingEventHeader* hijingHeader    = 
    dynamic_cast<AliGenHijingEventHeader*>(genHeader);
  // AliGenHydjetEventHeader* hydjetHeader    = 
  //   dynamic_cast<AliGenHydjetEventHeader*>(genHeader);
  AliGenEposEventHeader*   eposHeader      = 
    dynamic_cast<AliGenEposEventHeader*>(genHeader);

  AppendField(fProduction, colGeometry,  "Geometry");
  AppendField(fProduction, pythiaHeader, "Pythia");
  AppendField(fProduction, dpmHeader,    "DPM");
  AppendField(fProduction, gevHeader,    "GevSim");
  AppendField(fProduction, herwigHeader, "Herwig");
  AppendField(fProduction, hijingHeader, "Hijing");
  // AppendField(fProduction, hydjetHeader, "Hydjet");
  AppendField(fProduction, eposHeader,   "EPOS");
}

//____________________________________________________________________
UInt_t
AliFMDMCEventInspector::ProcessMC(AliMCEvent*       event, 
				  UInt_t&           triggers,
				  UShort_t&         ivz, 
				  Double_t&         vz,
				  Double_t&         b,
				  Int_t&            npart, 
				  Int_t&            nbin,
				  Double_t&         phiR)
{
  // 
  // Process the event 
  // 
  // Parameters:
  //   event     Input event 
  //   triggers  On return, the triggers fired 
  //   ivz       On return, the found vertex bin (1-based).  A zero
  //                  means outside of the defined vertex range
  //   vz        On return, the z position of the interaction
  //   b         On return, the impact parameter - < 0 if not available
  //   phiR      On return, the event plane angle - < 0 if not available 
  // 
  // Return:
  //    0 (or kOk) on success, otherwise a bit mask of error codes 
  //

  // Check that we have an event 
  if (!event) { 
    AliWarning("No MC event found for input event");
    return kNoEvent;
  }

  //Assign MC only triggers : True NSD etc.
  AliHeader*               header          = event->Header();
  AliGenEventHeader*       genHeader       = header->GenEventHeader();
  AliCollisionGeometry*    colGeometry     = 
    dynamic_cast<AliCollisionGeometry*>(genHeader);
  AliGenPythiaEventHeader* pythiaHeader    = 
    dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
  AliGenDPMjetEventHeader* dpmHeader       = 
    dynamic_cast<AliGenDPMjetEventHeader*>(genHeader);
  AliGenGeVSimEventHeader* gevHeader       = 
    dynamic_cast<AliGenGeVSimEventHeader*>(genHeader);
  AliGenHerwigEventHeader* herwigHeader    = 
    dynamic_cast<AliGenHerwigEventHeader*>(genHeader);
  // AliGenHijingEventHeader* hijingHeader    = 
  //   dynamic_cast<AliGenHijingEventHeader*>(genHeader);
  // AliGenHydjetEventHeader* hydjetHeader    = 
  //   dynamic_cast<AliGenHydjetEventHeader*>(genHeader);
  // AliGenEposEventHeader*   eposHeader      = 
  //   dynamic_cast<AliGenEposEventHeader*>(genHeader);
  
  // Check if this is a single diffractive event 
  Bool_t   sd    = kFALSE;
  Double_t phi   = -1111;
  npart          = 0;
  nbin           = 0;
  b              = -1;
  if (colGeometry) { 
    b     = colGeometry->ImpactParameter();
    phi   = colGeometry->ReactionPlaneAngle();
    npart = (colGeometry->ProjectileParticipants() + 
	     colGeometry->TargetParticipants());
    nbin  = colGeometry->NN();
  }
  if(pythiaHeader) {
    Int_t pythiaType = pythiaHeader->ProcessType();
    // 92 and 93 are SD 
    // 94 is DD 
    if (pythiaType==92 || pythiaType==93) sd = kTRUE;
    b     = pythiaHeader->GetImpactParameter();
    npart = 2; // Always 2 protons
    nbin  = 1; // Always 1 binary collision 
  }
  if(dpmHeader) { // Also an AliCollisionGeometry 
    Int_t processType = dpmHeader->ProcessType();
    // 1 & 4 are ND 
    // 5 & 6 are SD 
    // 7 is DD 
    if (processType == 5 || processType == 6)  sd = kTRUE;
  }
  if (gevHeader) { 
    phi  = gevHeader->GetEventPlane();
  }
  if (herwigHeader) {
    Int_t processType = herwigHeader->ProcessType();
    // This is a guess 
    if (processType == 5 || processType == 6)  sd = kTRUE;
    npart = 2; // Always 2 protons
    nbin  = 1; // Always 1 binary collision 
  }
  // if (hijingHeader) { 
  // b    = hijingHeader->ImpactParameter();
  // phi  = hijingHeader->ReactionPlaneAngle();
  // }
  // if (hydjetHeader) {
  //   b    = hydjetHeader->ImpactParameter();
  //   phi  = hydjetHeader->ReactionPlaneAngle();
  // }
  // if (eposHeader) {
  //   b    = eposHeader->ImpactParameter();
  // phi  = eposHeader->ReactionPlaneAngle();
  // }

  // Normalize event plane angle to [0,2pi]
  if (phi <= -1111) phiR = -1;
  else { 
    while (true) {
      if      (phi < 0)             phi += 2*TMath::Pi();
      else if (phi > 2*TMath::Pi()) phi -= 2*TMath::Pi();
      else                          break;
    }
  }

  // Do a check on particles
  sd = IsSingleDiffractive(event->Stack());

  // Set NSD flag
  if(!sd) {
    triggers |= AliAODForwardMult::kMCNSD;
    fHTriggers->Fill(kMCNSD+0.5);
  }
  
  // Get the primary vertex from EG 
  TArrayF vtx;
  genHeader->PrimaryVertex(vtx);
  vz = vtx[2];

  fHVertex->Fill(vz);
  fHPhiR->Fill(phiR);
  fHB->Fill(b);
  fHBvsPart->Fill(b, npart);
  fHBvsBin->Fill(b, nbin);

  // Check for the vertex bin 
  ivz = fVtxAxis.FindBin(vz);
  if (ivz <= 0 || ivz > fHEventsTr->GetXaxis()->GetNbins()) { 
    if (fDebug > 3) {
      AliWarning(Form("Vertex @ %f outside of range [%f,%f]", 
		      vz, fVtxAxis.GetXmin(), fVtxAxis.GetXmax())); }
    ivz = 0;
    return kBadVertex;
  }

  
  return kOk;
}

//____________________________________________________________________
namespace {
  Double_t rapidity(TParticle* p, Double_t mass)
  {
    Double_t pt = p->Pt();
    Double_t pz = p->Pz();
    Double_t ee = TMath::Sqrt(pt*pt+pz*pz+mass*mass);
    if (TMath::Abs(ee - TMath::Abs(pz)) < 1e-9) return TMath::Sign(1e30, pz);
    return .5 * TMath::Log((ee + pz) / (ee - pz)); 
  }
}

//____________________________________________________________________
Bool_t
AliFMDMCEventInspector::IsSingleDiffractive(AliStack* stack,
					    Double_t xiMin, 
					    Double_t xiMax) const
{
  // Re-implementation of AliPWG0Helper::IsHadronLevelSingleDiffrative
  // 
  // This is re-implemented here to be indendent of the PWG0 library. 
  TParticle* p1 = 0;     // Particle with least y 
  TParticle* p2 = 0;     // Particle with largest y 
  Double_t   y1 =  1e10; // y of p1 
  Double_t   y2 = -1e10; // y of p2 
  
  // Loop over primaries 
  for (Int_t i = 0; i < stack->GetNprimary(); i++) { 
    TParticle* p = stack->Particle(i);
    if (!p) continue;
    
    Int_t pdg = TMath::Abs(p->GetPdgCode());
    Int_t c1  = p->GetFirstDaughter();
    Int_t s   = p->GetStatusCode();
    Int_t mfl = Int_t(pdg/TMath::Power(10,Int_t(TMath::Log10(pdg))));

    // Select final state charm and beauty 
    if (c1 > -1 || s != 1) mfl = 0; 

    // Check if this is a primary, pi0, Sigma0, ???, or most
    // significant digit is larger than or equal to 4
    if (!(stack->IsPhysicalPrimary(i) || 
	  pdg == 111  || 
	  pdg == 3212 || 
	  pdg == 3124 || 
	  mfl >= 4)) continue;

    Int_t m1 = p->GetFirstMother();
    if (m1 > 0) { 
      TParticle* pm1  = stack->Particle(m1);
      Int_t      mpdg = TMath::Abs(pm1->GetPdgCode());
      // Check if mother is a p0, Simga0, or ???
      if (mpdg == 111 || mpdg == 3124 || mpdg == 3212) continue;
    }
    
    // Calculate the rapidity of the particle 
    Double_t mm = (pdg != 3124 ? p->GetMass() : 1.5195);
    Double_t yy = rapidity(p, mm);
    
    // Check if the rapidity of this particle is further out than any
    // of the preceding particles
    if (yy < y1) { 
      y1 = yy;
      p1 = p;
    }
    if (yy > y2) { 
      y2 = yy;
      p2 = p;
    }
  }
  if (!p1 || !p2) return false;
  
  // Calculate rapidities assuming protons 
  y1            = TMath::Abs(rapidity(p1, 0.938));
  y2            = TMath::Abs(rapidity(p2, 0.938));

  // Check if both or just one is a proton
  Int_t    pdg1 = p1->GetPdgCode();
  Int_t    pdg2 = p2->GetPdgCode();
  Int_t    arm  = -99999;
  if (pdg1 == 2212 && pdg2 == 2212) 
    arm = (y1 > y2 ? 0 : 1);
  else if (pdg1 == 2212) 
    arm = 0;
  else if (pdg2 == 2212) 
    arm = 1;
  else 
    return false;
  
  // Rapidity shift
  Double_t m02s = 1 - 2 * p1->Energy() / fEnergy; 
  Double_t m12s = 1 - 2 * p2->Energy() / fEnergy;
  
  if (arm == 0 && m02s > xiMin && m02s < xiMax) return true;
  if (arm == 1 && m12s > xiMin && m12s < xiMax) return true;
    
  return false;
}
//____________________________________________________________________
Bool_t
AliFMDMCEventInspector::ReadCentrality(const AliESDEvent* esd, 
				       Double_t& cent, 
				       UShort_t& qual) const
{
  // 
  // Read centrality from event 
  // 
  // Parameters:
  //    esd  Event 
  //    cent On return, the centrality or negative if not found
  // 
  // Return:
  //    False on error, true otherwise 
  //
  cent = -1;
  qual = 0;
  AliCentrality* centObj = const_cast<AliESDEvent*>(esd)->GetCentrality();
  if (!centObj) return true;

  qual = centObj->GetQuality();
  if (qual == 0x8) // Ignore ZDC outliers 
    cent = centObj->GetCentralityPercentileUnchecked("V0M");  
  else
    cent = centObj->GetCentralityPercentile("V0M");        

  return true;
}

//____________________________________________________________________
Bool_t
AliFMDMCEventInspector::CompareResults(Double_t vz,    Double_t trueVz, 
				       Double_t cent,  Double_t b,
				       Int_t    npart, Int_t    nbin)
{
  fHVzComp->Fill(trueVz, vz);
  fHBvsCent->Fill(b, cent);
  fHCentVsPart->Fill(npart, cent);
  fHCentVsBin->Fill(nbin, cent);

  return true;
}  
//
// EOF
//

