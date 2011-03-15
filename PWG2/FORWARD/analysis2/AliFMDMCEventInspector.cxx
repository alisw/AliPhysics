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
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenHijingEventHeader.h"
// #include "AliGenHydjetEventHeader.h"
#include "AliGenEposEventHeader.h"
#include "AliGenHerwigEventHeader.h"
#include "AliGenGeVSimEventHeader.h"
#include "AliHeader.h"
#include <TList.h>
//====================================================================
AliFMDMCEventInspector::AliFMDMCEventInspector()
  : AliFMDEventInspector(), 
    fHVertex(0),
    fHPhiR(0), 
    fHB(0)
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
    fHB(0)
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
    fHB(0)
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

  fHVertex = new TH1F("vertex", "True vertex distribution", 
		      vtxAxis.GetNbins(), 
		      vtxAxis.GetXmin(), 
		      vtxAxis.GetXmax());
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

  fHB = new TH1F("b", "Impact parameter", 125, 0, 25);
  fHB->SetXTitle("b [fm]");
  fHB->SetYTitle("# of events");
  fHB->SetFillColor(kGreen+1);
  fHB->SetFillStyle(3001);
  fHB->SetDirectory(0);
  fList->Add(fHB);
  
}

//____________________________________________________________________
UInt_t
AliFMDMCEventInspector::ProcessMC(AliMCEvent*       event, 
				  UInt_t&           triggers,
				  UShort_t&         ivz, 
				  Double_t&         vz,
				  Double_t&         b,
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
  AliGenPythiaEventHeader* pythiaHeader    = 
    dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
  AliGenDPMjetEventHeader* dpmHeader       = 
    dynamic_cast<AliGenDPMjetEventHeader*>(genHeader);
  AliGenGeVSimEventHeader* gevHeader       = 
    dynamic_cast<AliGenGeVSimEventHeader*>(genHeader);
  AliGenHijingEventHeader* hijingHeader    = 
    dynamic_cast<AliGenHijingEventHeader*>(genHeader);
  AliGenHerwigEventHeader* herwigHeader    = 
    dynamic_cast<AliGenHerwigEventHeader*>(genHeader);
  // AliGenHydjetEventHeader* hydjetHeader    = 
  //   dynamic_cast<AliGenHydjetEventHeader*>(genHeader);
  AliGenEposEventHeader*   eposHeader      = 
    dynamic_cast<AliGenEposEventHeader*>(genHeader);
  
  // Check if this is a single diffractive event 
  Bool_t   sd = kFALSE;
  Double_t phi = -1111;
  b            = -1;
  if(pythiaHeader) {
    Int_t pythiaType = pythiaHeader->ProcessType();
    if (pythiaType==92 || pythiaType==93) sd = kTRUE;
    b = pythiaHeader->GetImpactParameter();
  }
  if(dpmHeader) {
    Int_t processType = dpmHeader->ProcessType();
    if (processType == 5 || processType == 6)  sd = kTRUE;
    b    = dpmHeader->ImpactParameter();
    phi  = dpmHeader->ReactionPlaneAngle();

  }
  if (gevHeader) { 
    phi  = gevHeader->GetEventPlane();
  }
  if (hijingHeader) { 
    b    = hijingHeader->ImpactParameter();
    phi  = hijingHeader->ReactionPlaneAngle();
  }
  if (herwigHeader) {
    Int_t processType = herwigHeader->ProcessType();
    // This is a guess 
    if (processType == 5 || processType == 6)  sd = kTRUE;
  }
  // if (hydjetHeader) {
  //   b    = hydjetHeader->ImpactParameter();
  //   phi  = hydjetHeader->ReactionPlaneAngle();
  // }
  if (eposHeader) {
    b    = eposHeader->ImpactParameter();
    phi  = eposHeader->ReactionPlaneAngle();
  }

  // Normalize event plane angle to [0,2pi]
  if (phi <= -1111) phiR = -1;
  else { 
    while (true) {
      if      (phi < 0)             phi += 2*TMath::Pi();
      else if (phi > 2*TMath::Pi()) phi -= 2*TMath::Pi();
      else                          break;
    }
  }

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

  // Check for the vertex bin 
  ivz = fHEventsTr->GetXaxis()->FindBin(vz);
  if (ivz <= 0 || ivz > fHEventsTr->GetXaxis()->GetNbins()) { 
    if (fDebug > 3) {
      AliWarning(Form("Vertex @ %f outside of range [%f,%f]", 
		      vz, fHEventsTr->GetXaxis()->GetXmin(), 
		      fHEventsTr->GetXaxis()->GetXmax())); }
    ivz = 0;
    return kBadVertex;
  }

  
  return kOk;
}
//____________________________________________________________________
Bool_t
AliFMDMCEventInspector::ReadCentrality(const AliESDEvent* esd, Double_t& cent)
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
  AliCentrality* centObj = const_cast<AliESDEvent*>(esd)->GetCentrality();
  if (centObj) {
    // AliInfo(Form("Got centrality object %p with quality %d", 
    //              centObj, centObj->GetQuality()));
    // centObj->Print();
    if (centObj->GetQuality() == 0x8) 
      cent = centObj->GetCentralityPercentileUnchecked("V0M");  
    else
      cent = centObj->GetCentralityPercentile("V0M");        
  }
  // AliInfo(Form("Centrality is %f", cent));
  fHCent->Fill(cent);

  return true;
}

  
//
// EOF
//

