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
#include "AliCollisionGeometry.h"
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
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliHeader.h"
#include <TList.h>
#include <TH2F.h>
#include <TParticle.h>
#include <TMath.h>
#include <TParameter.h>

//====================================================================
AliFMDMCEventInspector::AliFMDMCEventInspector()
  : AliFMDEventInspector(), 
    fHVertex(0),
    fHVertexXY(0),
    fHPhiR(0), 
    fHB(0),
    fHMcC(0),
    fHBvsPart(0),
    fHBvsBin(0),
    fHBvsCent(0),
    fHVzComp(0),
    fHCentVsPart(0),
    fHCentVsBin(0),
    fHCentVsMcC(0),
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
    fHVertexXY(0),
    fHPhiR(0), 
    fHB(0),
    fHMcC(0),
    fHBvsPart(0),
    fHBvsBin(0),
    fHBvsCent(0),
    fHVzComp(0),
    fHCentVsPart(0),
    fHCentVsBin(0),
    fHCentVsMcC(0),
    fProduction("")
{
  // 
  // Constructor 
  // 
  // Parameters:
  //   name Name of object
  //
  fMC = true;
}

//____________________________________________________________________
AliFMDMCEventInspector::AliFMDMCEventInspector(const AliFMDMCEventInspector& o)
  : AliFMDEventInspector(o), 
    fHVertex(0),
    fHVertexXY(0),
    fHPhiR(0), 
    fHB(0),
    fHMcC(0),
    fHBvsPart(0),
    fHBvsBin(0),
    fHBvsCent(0),
    fHVzComp(0),
    fHCentVsPart(0),
    fHCentVsBin(0),
    fHCentVsMcC(0),
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
AliFMDMCEventInspector::SetupForData(const TAxis& vtxAxis)
{
  // 
  // Initialize the object 
  // 
  // Parameters:
  //   vtxAxis Vertex axis in use 
  //

  // We temporary disable the displaced vertices so we can initialize
  // the routine ourselves.
  Bool_t saveDisplaced  = AllowDisplaced();
  if (saveDisplaced) SetVertexMethod(kNormal);
  AliFMDEventInspector::SetupForData(vtxAxis);
  if (saveDisplaced) SetVertexMethod(kDisplaced);

  Int_t    maxPart = 450;
  Int_t    maxBin  = 225;
  Int_t    maxB    = 25;
  fHVertex = 0;
  if (vtxAxis.GetXbins() && vtxAxis.GetXbins()->GetArray())
    fHVertex = new TH1F("vertex", "True vertex distribution",
			vtxAxis.GetNbins(), vtxAxis.GetXbins()->GetArray());    
  else 
    fHVertex = new TH1F("vertex", "True vertex distribution",
			vtxAxis.GetNbins(), vtxAxis.GetXmin(),
			vtxAxis.GetXmax());
  fHVertex->SetXTitle("v_{z} [cm]");
  fHVertex->SetYTitle("# of events");
  fHVertex->SetFillColor(kGreen+1);
  fHVertex->SetFillStyle(3001);
  fHVertex->SetDirectory(0);
  // fHVertex->Sumw2();
  fList->Add(fHVertex);

  fHVertexXY = new TH2F("vertexXY", "True vertex distribution",
			100, -1, 1, 100, -1, 1);
  fHVertexXY->SetDirectory(0);
  fHVertexXY->SetXTitle("x [cm]");
  fHVertexXY->SetYTitle("y [cm]");
  fList->Add(fHVertexXY);

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

  fHMcC = static_cast<TH1F*>(fHCent->Clone("mcC"));
  fHMcC->SetFillColor(kCyan+2);
  fHMcC->SetDirectory(0);
  fList->Add(fHMcC);
  
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
  
  
  fHVzComp = 0;
  if (vtxAxis.GetXbins() && vtxAxis.GetXbins()->GetArray())
    fHVzComp = new TH2F("vzComparison", "v_{z} truth vs reconstructed",
			vtxAxis.GetNbins(), vtxAxis.GetXbins()->GetArray(),
			vtxAxis.GetNbins(), vtxAxis.GetXbins()->GetArray());
  else 
    fHVzComp = new TH2F("vzComparison", "v_{z} truth vs reconstructed",
			10*vtxAxis.GetNbins(), vtxAxis.GetXmin(),
			vtxAxis.GetXmax(), 10*vtxAxis.GetNbins(),
			vtxAxis.GetXmin(), vtxAxis.GetXmax());
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

  Int_t    nC = fHCent->GetNbinsX();
  Double_t cL = fHCent->GetXaxis()->GetXmin();
  Double_t cH = fHCent->GetXaxis()->GetXmax();
  fHCentVsMcC = new TH2F("centralityRecoVsMC", 
			 "Centrality from reconstruction vs MC derived", 
			 nC, cL, cH, nC, cL, cH);
  fHCentVsMcC->SetDirectory(0);
  fHCentVsMcC->SetStats(0);
  fHCentVsMcC->SetXTitle("Centralty from Reco [%]");
  fHCentVsMcC->SetYTitle("Centralty derived from Impact Par. [%]");
  fHCentVsMcC->SetZTitle("Events");
  fList->Add(fHCentVsMcC);

  if (AllowDisplaced()) fDisplacedVertex.SetupForData(fList, "", true);
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
				  TVector3&         ip,
				  Double_t&         b,
				  Double_t&         c,
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
  //   c         On return, centrality estimate - < 0 if not available 
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
  // MC events are _always_ B collisions
  triggers |= AliAODForwardMult::kB;
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
  Bool_t   sd    = false;
  Double_t phi   = -1111;
  Bool_t   egSD  = false;
  npart          = 0;
  nbin           = 0;
  b              = -1;
  c              = -1;
  phiR           = -1111;
  if (colGeometry) { 
    b     = colGeometry->ImpactParameter();
    phi   = colGeometry->ReactionPlaneAngle();
    npart = (colGeometry->ProjectileParticipants() + 
	     colGeometry->TargetParticipants());
    nbin  = colGeometry->NN();
  }
  if (fDebug && !colGeometry) { 
    AliWarningF("Collision header of class %s is not a CollisionHeader", 
		genHeader->ClassName());
  }
    
  if(pythiaHeader) {
    Int_t pythiaType = pythiaHeader->ProcessType();
    egSD = true; // We have SD flag in EG
    if (pythiaType <= 100) {
      // Pythia6 
      // 92 and 93 are SD 
      // 94 is DD 
      if (pythiaType==92 || pythiaType==93) sd = true;
    }
    else {
      // Pythia8 
      if (pythiaType==103 || pythiaType==104) sd = true;
    }
    b     = pythiaHeader->GetImpactParameter();
    npart = 2; // Always 2 protons
    nbin  = 1; // Always 1 binary collision 
  }
  if (b >= 0) { 
#if 0
    if      (b <  3.5)  c = 2.5; //0-5%
    else if (b <  4.95) c = 7.5; //5-10%
    else if (b <  6.98) c = 15; //10-20%
    else if (b <  8.55) c = 25; //20-30%
    else if (b <  9.88) c = 35; //30-40%
    else if (b < 11.04) c = 45; //40-50%
    else                c = 55; //50-60%
#else 
    // Updated 4th of November 2014 from 
    // cern.ch/twiki/bin/view/ALICE/CentStudies#Tables_with_centrality_bins_AN1
    Float_t np=0;
    Float_t nc=0;
    if      (0.00 >= b  && b < 1.57)  { c=0.5;  np=403.8; nc=1861; } 
    else if (1.57 >= b  && b < 2.22)  { c=1.5;  np=393.6; nc=1766; } 
    else if (2.22 >= b  && b < 2.71)  { c=2.5;  np=382.9; nc=1678; } 
    else if (2.71 >= b  && b < 3.13)  { c=3.5;  np=372;   nc=1597; }  
    else if (3.13 >= b  && b < 3.50)  { c=4.5;  np=361.1; nc=1520; } 
    else if (3.50 >= b  && b < 4.94)  { c=7.5;  np=329.4; nc=1316; } 
    else if (4.94 >= b  && b < 6.05)  { c=12.5; np=281.2; nc=1032; } 
    else if (6.05 >= b  && b < 6.98)  { c=17.5; np=239;   nc=809.8; }
    else if (6.98 >= b  && b < 7.81)  { c=22.5; np=202.1; nc=629.6; }
    else if (7.81 >= b  && b < 8.55)  { c=27.5; np=169.5; nc=483.7; }
    else if (8.55 >= b  && b < 9.23)  { c=32.5; np=141;   nc=366.7; }
    else if (9.23 >= b  && b < 9.88)  { c=37.5; np=116;   nc=273.4; }
    else if (9.88 >= b  && b < 10.47) { c=42.5; np=94.11; nc=199.4; } 
    else if (10.47 >= b && b < 11.04) { c=47.5; np=75.3;  nc=143.1; } 
    else if (11.04 >= b && b < 11.58) { c=52.5; np=59.24; nc=100.1; }
    else if (11.58 >= b && b < 12.09) { c=57.5; np=45.58; nc=68.46; }
    else if (12.09 >= b && b < 12.58) { c=62.5; np=34.33; nc=45.79; }
    else if (12.58 >= b && b < 13.05) { c=67.5; np=25.21; nc=29.92; }
    else if (13.05 >= b && b < 13.52) { c=72.5; np=17.96; nc=19.08; }
    else if (13.52 >= b && b < 13.97) { c=77.5; np=12.58; nc=12.07; }
    else if (13.97 >= b && b < 14.43) { c=82.5; np=8.812; nc=7.682; }
    else if (14.43 >= b && b < 14.96) { c=87.5; np=6.158; nc=4.904; }
    else if (14.96 >= b && b < 15.67) { c=92.5; np=4.376; nc=3.181; }
    else if (15.67 >= b && b < 20.00) { c=97.5; np=3.064; nc=1.994; }
    // Be careful to round off
    if (npart <= 0) npart = Int_t(np+.5);
    if (nbin  <= 0) nbin  = Int_t(nc+.5)/2;
#endif
  }
  if(dpmHeader) { // Also an AliCollisionGeometry 
    Int_t processType = dpmHeader->ProcessType();
    egSD = true; // We have SD flag in EG
    // 1 & 4 are ND 
    // 5 & 6 are SD 
    // 7 is DD 
    if (processType == 5 || processType == 6)  sd = kTRUE;
#if 0
      // The below - or rather a different implementation with some
      // errors - was proposed by Cvetan - I don't think it's right
      // though.  See also
      //
      //   https://cern.ch/twiki/pub/ALICE/PAPaperCentrality/normalization.pdf
      //   https://cern.ch/twiki/bin/view/ALICE/PAMCProductionStudies
      //
      Int_t nsd1=0, nsd2=0, ndd=0;
      Int_t npP = dpm->ProjectileParticipants();
      Int_t npT = dpm->TargetParticipants();
      // Get the numbeer of single and double diffractive participants
      dpm->GetNDiffractive(nsd1,nsd2,ndd);
      // Check if all partipants are single/double diffractive 
      if      ((ndd == 0) && ((npP == nsd1) || (npT == nsd2)))	sd = true;
      else if (ndd == (npP + npT))                       	dd = true;
      // Printf("Projectile: %3d (%3d) Target: %3d (%3d) DD: %3d Process: %d",
      // 	     npP, nsd1, npT, nsd2, ndd, type);
#endif 
  }
  if (gevHeader) { 
    phi  = gevHeader->GetEventPlane();
  }
  if (herwigHeader) {
    Int_t processType = herwigHeader->ProcessType();
    egSD = true; // We have SD flag in EG
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
    phiR = phi;
  }

  // Do a check on particles - but only if the EG does not give it or
  // the impact parameter is large enough
  if (!egSD && (b < 0 || b > 15)) sd = IsSingleDiffractive(event->Stack());

  // Set NSD flag
  if(!sd) {
    triggers |= AliAODForwardMult::kMCNSD;
    fHTriggers->Fill(kMCNSD+0.5);
  }
  
  // Get the primary vertex from EG 
  TArrayF vtx;
  genHeader->PrimaryVertex(vtx);
  ip.SetXYZ(vtx[0], vtx[1], vtx[2]);

  DMSG(fDebug, 1, "ip=(%f,%f,%f), phiR=%f, b=%f, npart=%d, nbin=%d", 
       vtx[0], vtx[1], vtx[2], phiR, b, npart, nbin);

  fHVertex->Fill(vtx[2]);
  fHVertexXY->Fill(vtx[0], vtx[1]);
  fHPhiR->Fill(phiR);
  fHB->Fill(b);
  fHMcC->Fill(c);
  fHBvsPart->Fill(b, npart);
  fHBvsBin->Fill(b, nbin);

  if(AllowDisplaced()) {
#if 0
    // Put the vertex at fixed locations 
    Double_t zvtx  = vz;
    Double_t ratio = zvtx/37.5;
    if(ratio > 0) ratio = ratio + 0.5;
    if(ratio < 0) ratio = ratio - 0.5;
    Int_t ratioInt = Int_t(ratio);
    zvtx = 37.5*((Double_t)ratioInt);
    if(TMath::Abs(zvtx) > 999) 
      return kBadVertex;
#endif
    if (!fDisplacedVertex.ProcessMC(event)) 
      return kBadVertex;
    if (fDisplacedVertex.IsSatellite())
      ip.SetZ(fDisplacedVertex.GetVertexZ());
  }

  // Check for the vertex bin 
  ivz = fVtxAxis.FindBin(ip.Z());
  if (ivz <= 0 || ivz > fHEventsTr->GetXaxis()->GetNbins()) { 
    if (fDebug > 3) {
      AliWarning(Form("Vertex @ %f outside of range [%f,%f]", 
		      ip.Z(), fVtxAxis.GetXmin(), fVtxAxis.GetXmax())); }
    ivz = 0;
    return kBadVertex;
  }

  
  return kOk;
}
//____________________________________________________________________
Bool_t
AliFMDMCEventInspector::ReadCentrality(const AliESDEvent& esd, 
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
  Bool_t ret = AliFMDEventInspector::ReadCentrality(esd, cent, qual);
  if (qual != 0) {
    AliCentrality* centObj = const_cast<AliESDEvent&>(esd).GetCentrality();
    if (!centObj)  return ret;

    // For MC, we allow `bad' centrality selections 
    cent = centObj->GetCentralityPercentileUnchecked(fCentMethod); 
  }
  return ret;
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
    Double_t mm = (pdg != 3124 && p->GetPDG() ? p->GetMass() : 1.5195);
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
  Double_t m02s = (fEnergy > 0 ? 1 - 2 * p1->Energy() / fEnergy : 0); 
  Double_t m12s = (fEnergy > 0 ? 1 - 2 * p2->Energy() / fEnergy : 0);
  
  if (arm == 0 && m02s > xiMin && m02s < xiMax) return true;
  if (arm == 1 && m12s > xiMin && m12s < xiMax) return true;
    
  return false;
}

//____________________________________________________________________
Bool_t
AliFMDMCEventInspector::CompareResults(Double_t vz,    Double_t trueVz, 
				       Double_t cent,  Double_t mcC, 
				       Double_t b,    
				       Int_t    npart, Int_t    nbin)
{
  fHVzComp->Fill(trueVz, vz);
  fHBvsCent->Fill(b, cent);
  fHCentVsPart->Fill(npart, cent);
  fHCentVsBin->Fill(nbin, cent);
  fHCentVsMcC->Fill(cent, mcC);

  return true;
}  


//
// EOF
//

