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
#include "AliFMDEventInspector.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"
#include "AliAODForwardMult.h"
#include "AliForwardUtil.h"
#include "AliCentrality.h"
#include <TH1.h>
#include <TList.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <iostream>
#include <iomanip>

//====================================================================
AliFMDEventInspector::AliFMDEventInspector()
  : TNamed(),
    fHEventsTr(0), 
    fHEventsTrVtx(0),
    fHEventsAccepted(0),
    fHEventsAcceptedXY(0),
    fHTriggers(0),
    fHType(0),
    fHWords(0),
    fHCent(0),
    fHCentVsQual(0),
    fLowFluxCut(1000),
    fMaxVzErr(0.2),
    fList(0),
    fEnergy(0),
    fField(999), 
    fCollisionSystem(kUnknown),
    fDebug(0),
    fCentAxis(0),
    fVtxAxis(10,-10,10),
    fUseFirstPhysicsVertex(true)
{
  // 
  // Constructor 
  //
}

//____________________________________________________________________
AliFMDEventInspector::AliFMDEventInspector(const char* name)
  : TNamed("fmdEventInspector", name),
    fHEventsTr(0), 
    fHEventsTrVtx(0), 
    fHEventsAccepted(0),
    fHEventsAcceptedXY(0),
    fHTriggers(0),
    fHType(0),
    fHWords(0),
    fHCent(0),
    fHCentVsQual(0),
    fLowFluxCut(1000),
    fMaxVzErr(0.2),
    fList(0),
    fEnergy(0),
    fField(999), 
    fCollisionSystem(kUnknown),
    fDebug(0),
    fCentAxis(0),
    fVtxAxis(10,-10,10),
    fUseFirstPhysicsVertex(true)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //   name Name of object
  //
}

//____________________________________________________________________
AliFMDEventInspector::AliFMDEventInspector(const AliFMDEventInspector& o)
  : TNamed(o), 
    fHEventsTr(o.fHEventsTr), 
    fHEventsTrVtx(o.fHEventsTrVtx), 
    fHEventsAccepted(o.fHEventsAccepted),
    fHEventsAcceptedXY(o.fHEventsAcceptedXY),
    fHTriggers(o.fHTriggers),
    fHType(o.fHType),
    fHWords(o.fHWords),
    fHCent(o.fHCent),
    fHCentVsQual(o.fHCentVsQual),
    fLowFluxCut(o.fLowFluxCut),
    fMaxVzErr(o.fMaxVzErr),
    fList(o.fList),
    fEnergy(o.fEnergy),
    fField(o.fField), 
    fCollisionSystem(o.fCollisionSystem),
    fDebug(0),
    fCentAxis(0),
    fVtxAxis(o.fVtxAxis),
    fUseFirstPhysicsVertex(o.fUseFirstPhysicsVertex)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //   o Object to copy from 
  //
}

//____________________________________________________________________
AliFMDEventInspector::~AliFMDEventInspector()
{
  // 
  // Destructor 
  //
  if (fList)         delete fList;
}
//____________________________________________________________________
AliFMDEventInspector&
AliFMDEventInspector::operator=(const AliFMDEventInspector& o)
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
  TNamed::operator=(o);
  fHEventsTr         = o.fHEventsTr;
  fHEventsTrVtx      = o.fHEventsTrVtx;
  fHEventsAccepted   = o.fHEventsAccepted;
  fHEventsAcceptedXY = o.fHEventsAcceptedXY;
  fHTriggers         = o.fHTriggers;
  fHType             = o.fHType;
  fHWords            = o.fHWords;
  fHCent             = o.fHCent;
  fHCentVsQual       = o.fHCentVsQual;
  fLowFluxCut        = o.fLowFluxCut;
  fMaxVzErr          = o.fMaxVzErr;
  fDebug             = o.fDebug;
  fList              = (o.fList ? new TList : 0);
  fEnergy            = o.fEnergy;
  fField             = o.fField;
  fCollisionSystem   = o.fCollisionSystem;
  fVtxAxis.Set(o.fVtxAxis.GetNbins(), o.fVtxAxis.GetXmin(), 
	       o.fVtxAxis.GetXmax());
  
  fUseFirstPhysicsVertex = o.fUseFirstPhysicsVertex;
  
  if (fList) { 
    fList->SetName(GetName());
    if (fHEventsTr)    fList->Add(fHEventsTr);
    if (fHEventsTrVtx) fList->Add(fHEventsTrVtx);
    if (fHTriggers)    fList->Add(fHTriggers);
    if (fHType)        fList->Add(fHType);
    if (fHWords)       fList->Add(fHWords);
    if (fHCent)        fList->Add(fHCent);
    if (fHCentVsQual)  fList->Add(fHCentVsQual);
  }
  return *this;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::FetchHistograms(const TList* d, 
				      TH1I*& hEventsTr, 
				      TH1I*& hEventsTrVtx, 
				      TH1I*& hTriggers) const
{
  // 
  // Fetch our histograms from the passed list 
  // 
  // Parameters:
  //   d             Input
  //   hEventsTr     On return, pointer to histogram, or null
  //   hEventsTrVtx  On return, pointer to histogram, or null
  //   hTriggers     On return, pointer to histogram, or null
  // 
  // Return:
  //    true on success, false otherwise 
  //
  hEventsTr    = 0;
  hEventsTrVtx = 0;
  hTriggers    = 0;
  TList* dd    = dynamic_cast<TList*>(d->FindObject(GetName()));
  if (!dd) return kFALSE;
  
  hEventsTr    = dynamic_cast<TH1I*>(dd->FindObject("nEventsTr"));
  hEventsTrVtx = dynamic_cast<TH1I*>(dd->FindObject("nEventsTrVtx"));
  hTriggers    = dynamic_cast<TH1I*>(dd->FindObject("triggers"));

  if (!hEventsTr || !hEventsTrVtx || !hTriggers) return kFALSE;
  return kTRUE;
}
//____________________________________________________________________
void
AliFMDEventInspector::Init(const TAxis& vtxAxis)
{
  // 
  // Initialize the object 
  // 
  // Parameters:
  //   vtxAxis Vertex axis in use 
  //
  
  // -1.5 -0.5 0.5 1.5 ... 89.5 ... 100.5
  // ----- 92 number --------- ---- 1 ---
  TArrayD limits(93);
  for (Int_t i = 0; i < 92; i++) limits[i] = -1.5 + i;
  limits[92] = 100.5;

  fVtxAxis.Set(vtxAxis.GetNbins(), vtxAxis.GetXmin(), vtxAxis.GetXmax());
  
  fCentAxis  = new TAxis(limits.GetSize()-1, limits.GetArray());
  fHEventsTr = new TH1I("nEventsTr", "Number of events w/trigger", 
			4*vtxAxis.GetNbins(), 
			2*vtxAxis.GetXmin(), 
			2*vtxAxis.GetXmax());
  fHEventsTr->SetXTitle("v_{z} [cm]");
  fHEventsTr->SetYTitle("# of events");
  fHEventsTr->SetFillColor(kRed+1);
  fHEventsTr->SetFillStyle(3001);
  fHEventsTr->SetDirectory(0);
  // fHEventsTr->Sumw2();
  fList->Add(fHEventsTr);

  fHEventsTrVtx = static_cast<TH1I*>(fHEventsTr->Clone("nEventsTrVtx"));
  fHEventsTrVtx->SetTitle("Number of events w/trigger and vertex"); 
  fHEventsTrVtx->SetFillColor(kBlue+1);
  fHEventsTrVtx->SetDirectory(0);
  // fHEventsTrVtx->Sumw2();
  fList->Add(fHEventsTrVtx);

  fHEventsAccepted = new TH1I("nEventsAccepted", 
			      "Number of events  w/trigger and vertex in range",
			      2*vtxAxis.GetNbins(), 
			      2*vtxAxis.GetXmin(), 
			      2*vtxAxis.GetXmax());
  fHEventsAccepted->SetXTitle("v_{z} [cm]");
  fHEventsAccepted->SetYTitle("# of events");
  fHEventsAccepted->SetFillColor(kGreen+1);
  fHEventsAccepted->SetFillStyle(3001);
  fHEventsAccepted->SetDirectory(0);
  // fHEventsAccepted->Sumw2();
  fList->Add(fHEventsAccepted);
			      
  fHEventsAcceptedXY = new TH2D("nEventsAcceptedXY", 
			      "XY vertex w/trigger and Z vertex in range",
				1000,-1,1,1000,-1,1);
  
  fHEventsAcceptedXY->SetXTitle("v_{x} [cm]");
  fHEventsAcceptedXY->SetYTitle("v_{y} [cm]");
  fHEventsAcceptedXY->SetDirectory(0);
  // fHEventsAccepted->Sumw2();
  fList->Add(fHEventsAcceptedXY);

  
  fHTriggers = new TH1I("triggers", "Triggers", kOffline+1, 0, kOffline+1);
  fHTriggers->SetFillColor(kRed+1);
  fHTriggers->SetFillStyle(3001);
  fHTriggers->SetStats(0);
  fHTriggers->SetDirectory(0);
  fHTriggers->GetXaxis()->SetBinLabel(kInel   +1,"INEL");
  fHTriggers->GetXaxis()->SetBinLabel(kInelGt0+1,"INEL>0");
  fHTriggers->GetXaxis()->SetBinLabel(kNSD    +1,"NSD");
  fHTriggers->GetXaxis()->SetBinLabel(kEmpty  +1,"Empty");
  fHTriggers->GetXaxis()->SetBinLabel(kA      +1,"A");
  fHTriggers->GetXaxis()->SetBinLabel(kB      +1,"B");
  fHTriggers->GetXaxis()->SetBinLabel(kC      +1,"C");
  fHTriggers->GetXaxis()->SetBinLabel(kE      +1,"E");
  fHTriggers->GetXaxis()->SetBinLabel(kPileUp +1,"Pileup");
  fHTriggers->GetXaxis()->SetBinLabel(kMCNSD  +1,"NSD_{MC}");
  fHTriggers->GetXaxis()->SetBinLabel(kOffline+1,"Offline");
  fList->Add(fHTriggers);

  fHType = new TH1I("type", Form("Event type (cut: SPD mult>%d)", 
				 fLowFluxCut), 2, -.5, 1.5);
  fHType->SetFillColor(kRed+1);
  fHType->SetFillStyle(3001);
  fHType->SetStats(0);
  fHType->SetDirectory(0);
  fHType->GetXaxis()->SetBinLabel(1,"Low-flux");
  fHType->GetXaxis()->SetBinLabel(2,"High-flux");
  fList->Add(fHType);


  fHWords = new TH1I("words", "Trigger words seen", 1, 0, 0); 
  fHWords->SetFillColor(kBlue+1);
  fHWords->SetFillStyle(3001);
  fHWords->SetStats(0);
  fHWords->SetDirectory(0);
  fHWords->SetBit(TH1::kCanRebin);
  fList->Add(fHWords);

  fHCent = new TH1F("cent", "Centrality", limits.GetSize()-1,limits.GetArray());
  fHCent->SetFillColor(kBlue+1);
  fHCent->SetFillStyle(3001);
  fHCent->SetStats(0);
  fHCent->SetDirectory(0);
  fHCent->SetXTitle("Centrality [%]");
  fHCent->SetYTitle("Events");
  fList->Add(fHCent);

  fHCentVsQual = new TH2F("centVsQuality", "Quality vs Centrality", 
			  5, 0, 5, limits.GetSize()-1, limits.GetArray());
  fHCentVsQual->SetXTitle("Quality");
  fHCentVsQual->SetYTitle("Centrality [%]");
  fHCentVsQual->SetZTitle("Events");
  fHCentVsQual->GetXaxis()->SetBinLabel(1, "OK");
  fHCentVsQual->GetXaxis()->SetBinLabel(2, "Outside v_{z} cut");
  fHCentVsQual->GetXaxis()->SetBinLabel(3, "V0 vs SPD outlier");
  fHCentVsQual->GetXaxis()->SetBinLabel(4, "V0 vs TPC outlier");
  fHCentVsQual->GetXaxis()->SetBinLabel(5, "V0 vs ZDC outlier");
  fList->Add(fHCentVsQual);
}

//____________________________________________________________________
void
AliFMDEventInspector::StoreInformation()
{
  // Write TNamed objects to output list containing information about
  // the running conditions 
  if (!fList) return;

  TNamed* sys = new TNamed("sys", "");
  TNamed* sNN = new TNamed("sNN", "");
  TNamed* fld = new TNamed("field", "");
  sys->SetTitle(AliForwardUtil::CollisionSystemString(fCollisionSystem));
  sNN->SetTitle(AliForwardUtil::CenterOfMassEnergyString(fEnergy));
  fld->SetTitle(AliForwardUtil::MagneticFieldString(fField));
  sys->SetUniqueID(fCollisionSystem);
  sNN->SetUniqueID(fEnergy);
  fld->SetUniqueID(fField);

  fList->Add(sys);
  fList->Add(sNN);
  fList->Add(fld);
}

//____________________________________________________________________
void
AliFMDEventInspector::DefineOutput(TList* dir)
{
  // 
  // Define the output histograms.  These are put in a sub list of the
  // passed list.   The histograms are merged before the parent task calls 
  // AliAnalysisTaskSE::Terminate 
  // 
  //   dir Directory to add to 
  //
  fList = new TList;
  fList->SetName(GetName());
  dir->Add(fList);
}

//____________________________________________________________________
UInt_t
AliFMDEventInspector::Process(const AliESDEvent* event, 
			      UInt_t&            triggers,
			      Bool_t&            lowFlux,
			      UShort_t&          ivz, 
			      Double_t&          vz,
			      Double_t&          cent,
			      UShort_t&          nClusters)
{
  // 
  // Process the event 
  // 
  // Parameters:
  //   event     Input event 
  //   triggers  On return, the triggers fired 
  //   lowFlux   On return, true if the event is considered a low-flux 
  //                  event (according to the setting of fLowFluxCut) 
  //   ivz       On return, the found vertex bin (1-based).  A zero
  //                  means outside of the defined vertex range
  //   vz        On return, the z position of the interaction
  //   cent      On return, the centrality - if not available < 0
  // 
  // Return:
  //    0 (or kOk) on success, otherwise a bit mask of error codes 
  //

  // --- Check that we have an event ---------------------------------
  if (!event) { 
    AliWarning("No ESD event found for input event");
    return kNoEvent;
  }

  // --- Read trigger information from the ESD and store in AOD object
  if (!ReadTriggers(event, triggers, nClusters)) { 
    if (fDebug > 2) {
      AliWarning("Failed to read triggers from ESD"); }
    return kNoTriggers;
  }

  // --- Check if this is a high-flux event --------------------------
  const AliMultiplicity* testmult = event->GetMultiplicity();
  if (!testmult) {
    if (fDebug > 3) {
      AliWarning("No central multiplicity object found"); }
  }
  else 
    lowFlux = testmult->GetNumberOfTracklets() < fLowFluxCut;

  fHType->Fill(lowFlux ? 0 : 1);
  
  // --- Read centrality information 
  cent          = -10;
  UShort_t qual = 0;
  if (!ReadCentrality(event, cent, qual)) {
    if (fDebug > 3) 
      AliWarning("Failed to get centrality");
  }
  fHCent->Fill(cent);
  if (qual == 0) fHCentVsQual->Fill(0., cent);
  else { 
    for (UShort_t i = 0; i < 4; i++) 
      if (qual & (1 << i)) fHCentVsQual->Fill(Double_t(i+1), cent);
  }

  // --- Get the vertex information ----------------------------------
  
  Double_t vx = 0;
  Double_t vy = 0;
  vz          = 0;
  
  Bool_t vzOk = ReadVertex(event, vz,vx,vy);

  fHEventsTr->Fill(vz);
  if (!vzOk) { 
    if (fDebug > 3) {
      AliWarning("Failed to read vertex from ESD"); }
    return kNoVertex;
  }
  fHEventsTrVtx->Fill(vz);
  
  // --- Get the vertex bin ------------------------------------------
  ivz = fVtxAxis.FindBin(vz);
  if (ivz <= 0 || ivz > fVtxAxis.GetNbins()) { 
    if (fDebug > 3) {
      AliWarning(Form("Vertex @ %f outside of range [%f,%f]", 
		      vz, fVtxAxis.GetXmin(), fVtxAxis.GetXmax())); 
    }
    ivz = 0;
    return kBadVertex;
  }
  fHEventsAccepted->Fill(vz);
  fHEventsAcceptedXY->Fill(vx,vy);
  
  // --- Check the FMD ESD data --------------------------------------
  if (!event->GetFMDData()) { 
    if (fDebug > 3) {
      AliWarning("No FMD data found in ESD"); }
    return kNoFMD;
  }

  
  return kOk;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::ReadCentrality(const AliESDEvent* esd, 
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
  if (!centObj)  return true;

  // AliInfo(Form("Got centrality object %p with quality %d", 
  //              centObj, centObj->GetQuality()));
  // centObj->Print();
  cent = centObj->GetCentralityPercentile("V0M");  
  //cent = centObj->GetCentralityPercentile("ZEMvsZDC");  
  qual = centObj->GetQuality();

  return true;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::ReadTriggers(const AliESDEvent* esd, UInt_t& triggers,
				   UShort_t& nClusters)
{
  // 
  // Read the trigger information from the ESD event 
  // 
  // Parameters:
  //   esd        ESD event 
  //   triggers   On return, contains the trigger bits 
  // 
  // Return:
  //    @c true on success, @c false otherwise 
  //
  triggers = 0;

  // Get the analysis manager - should always be there 
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  if (!am) { 
    AliWarning("No analysis manager defined!");
    return kFALSE;
  }

  // Get the input handler - should always be there 
  AliInputEventHandler* ih = 
    static_cast<AliInputEventHandler*>(am->GetInputEventHandler());
  if (!ih) { 
    AliWarning("No input handler");
    return kFALSE;
  }
  
  // Check if this is a collision candidate (MB)
  // Note, that we should use the value cached in the input 
  // handler rather than calling IsCollisionCandiate directly 
  // on the AliPhysicsSelection obejct.  If we called the latter
  // then the AliPhysicsSelection object would overcount by a 
  // factor of 2! :-(
  Bool_t offline  = ih->IsEventSelected() ;
  Bool_t fastonly = (ih->IsEventSelected() & AliVEvent::kFastOnly);
  TString trigStr = esd->GetFiredTriggerClasses();
  
  //If we have the MC input handler,  this must be MC
  Bool_t isMC = am->GetMCtruthEventHandler() != 0;

  // For the 2.76 TeV p+p run, the FMD ran in the slow partition 
  // so it received no triggers from the fast partition. Therefore
  // the fast triggers are removed here but not for MC where all 
  // triggers are fast.
  if(TMath::Abs(fEnergy - 2750.) < 20 && 
     fCollisionSystem == AliForwardUtil::kPP &&
     !isMC)
    if (fastonly) offline = false;
  nClusters = 0;
  
  // MUON triggers are not strictly minimum bias (MB) so they are removed (HHD)
  
  if(offline && trigStr.Contains("CMUS1")) offline = false;
    
  if (offline ) {
    triggers |= AliAODForwardMult::kOffline;
    triggers |= AliAODForwardMult::kInel;
    fHTriggers->Fill(kOffline+0.5);

    // If this is inel, see if we have a tracklet 
    const AliMultiplicity* spdmult = esd->GetMultiplicity();
    if (!spdmult) {
      AliWarning("No SPD multiplicity");
    }
    else { 
      // Check if we have one or more tracklets 
      // in the range -1 < eta < 1 to set the INEL>0 
      // trigger flag. 
      // 
      // Also count tracklets as a single cluster 
      Int_t n = spdmult->GetNumberOfTracklets();
      for (Int_t j = 0; j < n; j++) { 
	if(TMath::Abs(spdmult->GetEta(j)) < 1) { 
	  triggers |= AliAODForwardMult::kInelGt0;
	  nClusters++;
	}
      }
      n = spdmult->GetNumberOfSingleClusters();
      for (Int_t j = 0; j < n; j++) { 
	Double_t eta = -TMath::Log(TMath::Tan(spdmult->GetThetaSingle(j)/2.));
	if (TMath::Abs(eta) < 1) nClusters++;
      }
    }
    if (nClusters > 0) triggers |= AliAODForwardMult::kNClusterGt0;
  }
  
  // Analyse some trigger stuff 
  AliTriggerAnalysis ta;
  if (ta.IsOfflineTriggerFired(esd, AliTriggerAnalysis::kNSD1)) 
    triggers |= AliAODForwardMult::kNSD;
 
  // Check for multiple vertices (pile-up) with at least 3
  // contributors and at least 0.8cm from the primary vertex
  Bool_t pileup = kFALSE;
  if(fCollisionSystem == AliForwardUtil::kPP)
    pileup =  esd->IsPileupFromSPD(3,0.8);
  if (pileup) {
    triggers |= AliAODForwardMult::kPileUp;
    fHTriggers->Fill(kPileUp+.5);
  }
  
  // Get trigger stuff 
  
  //TString trigStr = esd->GetFiredTriggerClasses();
  // AliWarning(Form("Fired trigger classes: %s", trigStr.Data()));
  fHWords->Fill(trigStr.Data(), 1);
#if 0
  if (trigStr.Contains("MB1") || trigStr.Contains("MBBG3"))
    triggers |= AliAOODForwardMult::kB;
  if (trigStr.Contains("COTA")) 
    triggers |= AliAODForwardMult::kA;
  if (trigStr.Contains("COTC")) 
    triggers |= AliAODForwardMult::kC;
#endif
  if (trigStr.Contains("CBEAMB-ABCE-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kEmpty;
    fHTriggers->Fill(kEmpty+.5);
  }
  
  // Check for B triggers
  if (trigStr.Contains("CINT1B-ABCE-NOPF-ALL")   ||   // Early pp
      trigStr.Contains("CINT1-B-NOPF-ALLNOTRD")  ||   // Late pp 
      trigStr.Contains("CINT1-B-NOPF-FASTNOTRD") ||   // Late pp
      //trigStr.Contains("CMUS1-B-NOPF-MUON")      ||   // Late pp -- HHD 160811
      trigStr.Contains("CSMBB-ABCE-NOPF-ALL")    ||   // pp
      trigStr.Contains("CMBACS2-B-NOPF-ALL")     ||   // PbPb
      // trigStr.Contains("C0SMH-B-NOPF-ALL")    ||   // PbPb - high mult
      trigStr.Contains("CMBS2A-B-NOPF-ALL")      ||   // PbPb
      trigStr.Contains("CMBS2C-B-NOPF-ALL")      ||   // PbPb
      trigStr.Contains("CMBAC-B-NOPF-ALL")       ||   // PbPb
      // trigStr.Contains("C0SMH-B-NOPF-ALL")    ||   // PbPb - high mult
      trigStr.Contains("CMBACS2-B-NOPF-ALLNOTRD")     // PbPb
      // trigStr.Contains("C0SMH-B-NOPF-ALLNOTRD")    // PbPb - high mult
      ) {
    Bool_t bTrigger = kTRUE;
    if ( trigStr.Contains("CINT1-B-NOPF-FASTNOTRD") && 
	 !trigStr.Contains("CINT1-B-NOPF-ALLNOTRD") && 
	 TMath::Abs(fEnergy - 2750.) < 20 && 
	 fCollisionSystem == AliForwardUtil::kPP)
      bTrigger = kFALSE;
    if(bTrigger) {
      triggers |= AliAODForwardMult::kB;
      fHTriggers->Fill(kB+.5);
    }
  }
  
  // Check for A triggers
  if (trigStr.Contains("CINT1A-ABCE-NOPF-ALL")    ||   // Early pp
      trigStr.Contains("CINT1-AC-NOPF-ALLNOTRD")  ||   // Late pp
      trigStr.Contains("CINT1-AC-NOPF-FASTNOTRD") ||   // Late pp
      (trigStr.Contains("CSMBA-ABCE-NOPF-ALL") && 
       !(triggers & AliAODForwardMult::kB))       ||   // pp
      trigStr.Contains("CMBACS2-A-NOPF-ALL")      ||   // PbPb
      // trigStr.Contains("C0SMH-A-NOPF-ALL")     ||   // PbPb - high mult
      trigStr.Contains("CMBS2A-A-NOPF-ALL")       ||   // PbPb
      trigStr.Contains("CMBS2C-A-NOPF-ALL")       ||   // PbPb
      trigStr.Contains("CMBAC-A-NOPF-ALL")        ||   // PbPb
      // trigStr.Contains("C0SMH-A-NOPF-ALL")     ||   // PbPb - high mult
      trigStr.Contains("CMBACS2-A-NOPF-ALLNOTRD")     // PbPb
      // trigStr.Contains("C0SMH-A-NOPF-ALLNOTRD")    // PbPb - high mult
      ) {
    triggers |= AliAODForwardMult::kA;
    fHTriggers->Fill(kA+.5);
  }

  // Check for C triggers
  if (trigStr.Contains("CINT1C-ABCE-NOPF-ALL")   ||  // Early pp
      (trigStr.Contains("CSMBC-ABCE-NOPF-ALL") && 
       !(triggers & AliAODForwardMult::kB))      ||   // pp
      trigStr.Contains("CMBACS2-C-NOPF-ALL")     ||   // PbPb
      // trigStr.Contains("C0SMH-B-NOPF-ALL")    ||   // PbPb - high mult
      trigStr.Contains("CMBS2A-C-NOPF-ALL")      ||   // PbPb
      trigStr.Contains("CMBS2C-C-NOPF-ALL")      ||   // PbPb
      trigStr.Contains("CMBAC-C-NOPF-ALL")       ||   // PbPb
      // trigStr.Contains("C0SMH-B-NOPF-ALL")    ||   // PbPb - high mult
      trigStr.Contains("CMBACS2-C-NOPF-ALLNOTRD")     // PbPb
      // trigStr.Contains("C0SMH-B-NOPF-ALLNOTRD")    // PbPb - high mult
      ) {
    triggers |= AliAODForwardMult::kC;
    fHTriggers->Fill(kC+.5);
  }

  // Check for E triggers 
  if (trigStr.Contains("CINT1-E-NOPF-ALL")       ||   // Early pp 
      trigStr.Contains("CINT1-E-NOPF-ALLNOTRD")  ||   // Late pp 
      trigStr.Contains("CINT1-E-NOPF-FASTNOTRD") ||   // Late pp 
      trigStr.Contains("CMBACS2-E-NOPF-ALL")     ||   // PbPb
      // trigStr.Contains("C0SMH-B-NOPF-ALL")    ||   // PbPb - high mult
      trigStr.Contains("CMBS2A-E-NOPF-ALL")      ||   // PbPb
      trigStr.Contains("CMBS2C-E-NOPF-ALL")      ||   // PbPb
      trigStr.Contains("CMBAC-E-NOPF-ALL")       ||   // PbPb
      // trigStr.Contains("C0SMH-B-NOPF-ALL")    ||   // PbPb - high mult
      trigStr.Contains("CMBACS2-E-NOPF-ALLNOTRD")     // PbPb
      // trigStr.Contains("C0SMH-B-NOPF-ALLNOTRD")    // PbPb - high mult
      ) {
    triggers |= AliAODForwardMult::kE;
    fHTriggers->Fill(kE+.5);
  }

  // Now check - if we have a collision - for offline triggers and
  // fill histogram.
  if (triggers & AliAODForwardMult::kB) {
    if (triggers & AliAODForwardMult::kInel) 
      fHTriggers->Fill(kInel);
    
    if (triggers & AliAODForwardMult::kInelGt0)
      fHTriggers->Fill(kInelGt0+.5);
    
    if (triggers & AliAODForwardMult::kNSD)
      fHTriggers->Fill(kNSD+.5);
  }
  
  return kTRUE;
}
//____________________________________________________________________
Bool_t
AliFMDEventInspector::ReadVertex(const AliESDEvent* esd, 
				 Double_t& vz, 
				 Double_t& vx, 
				 Double_t& vy)
{
  // 
  // Read the vertex information from the ESD event 
  // 
  // Parameters:
  //   esd  ESD event 
  //   vz   On return, the vertex Z position 
  // 
  // Return:
  //    @c true on success, @c false otherwise 
  //
  vz = 0;
  vx = 1024;
  vy = 1024;
  if(fUseFirstPhysicsVertex) {
    // This is the code used by the 1st physics people 
    const AliESDVertex* vertex    = esd->GetPrimaryVertex();
    if (!vertex  || !vertex->GetStatus()) {
      if (fDebug > 2) {
	AliWarning(Form("No primary vertex (%p) or bad status %d", 
			vertex, (vertex ? vertex->GetStatus() : -1)));
      }
      return false;
    }
    const AliESDVertex* vertexSPD = esd->GetPrimaryVertexSPD();
    if (!vertexSPD || !vertexSPD->GetStatus()) {
      if (fDebug > 2) {
	AliWarning(Form("No primary SPD vertex (%p) or bad status %d", 
			vertexSPD, (vertexSPD ? vertexSPD->GetStatus() : -1)));
      }
      return false;
    }
    
    // if vertex is from SPD vertexZ, require more stringent cuts 
    if (vertex->IsFromVertexerZ()) { 
      if (vertex->GetDispersion() > fMaxVzErr || 
	  vertex->GetZRes() > 1.25 * fMaxVzErr) {
	if (fDebug > 2) {
	  AliWarning(Form("Dispersion %f > %f or resolution %f > %f",
			  vertex->GetDispersion(), fMaxVzErr,
			  vertex->GetZRes(), 1.25 * fMaxVzErr)); 
	}
	return false;
      }
    }
    vz = vertex->GetZ();
    
    if(!vertex->IsFromVertexerZ()) {
      vx = vertex->GetX();
      vy = vertex->GetY();
    }
    return true;
  }
  else { //Use standard SPD vertex (perhaps preferable for Pb+Pb)
   
    // Get the vertex 
    const AliESDVertex* vertex = esd->GetPrimaryVertexSPD();
    if (!vertex) { 
      if (fDebug > 2) {
	AliWarning("No SPD vertex found in ESD"); }
      return kFALSE;
    }
    
    // Check that enough tracklets contributed 
    if(vertex->GetNContributors() <= 0) {
      if (fDebug > 2) {
	AliWarning(Form("Number of contributors to vertex is %d<=0",
			vertex->GetNContributors())); }
      vz = 0;
      return kFALSE;
    } 
    // Check that the uncertainty isn't too large 
    if (vertex->GetZRes() > fMaxVzErr) { 
      if (fDebug > 2) {
	AliWarning(Form("Uncertaintity in Z of vertex is too large %f > %f", 
		      vertex->GetZRes(), fMaxVzErr)); }
      return kFALSE;
    }
    
    // Get the z coordiante 
    vz = vertex->GetZ();
    const AliESDVertex* vertexXY = esd->GetPrimaryVertex();
    
    if(!vertexXY->IsFromVertexerZ()) {
      vx = vertexXY->GetX();
      vy = vertexXY->GetY();
    }
    return kTRUE;
  } 
}
  
//____________________________________________________________________
Bool_t
AliFMDEventInspector::ReadRunDetails(const AliESDEvent* esd)
{
  // 
  // Read the collision system, collision energy, and L3 field setting
  // from the ESD
  // 
  // Parameters:
  //   esd ESD to get information from 
  // 
  // Return:
  //    true on success, false 
  //
  // AliInfo(Form("Parameters from 1st ESD event: cms=%s, sNN=%f, field=%f",
  // 	       esd->GetBeamType(), 2*esd->GetBeamEnergy(), 
  // 	       esd->GetMagneticField()));
  fCollisionSystem = 
    AliForwardUtil::ParseCollisionSystem(esd->GetBeamType());
  fEnergy          = 
    AliForwardUtil::ParseCenterOfMassEnergy(fCollisionSystem,	
					    2 * esd->GetBeamEnergy());
  fField           = 
    AliForwardUtil::ParseMagneticField(esd->GetMagneticField());

  StoreInformation();
  if (fCollisionSystem   == AliForwardUtil::kUnknown || 
      fEnergy            <= 0                        || 
      TMath::Abs(fField) >  10) 
    return kFALSE;

  return kTRUE;
}

//____________________________________________________________________
void
AliFMDEventInspector::Print(Option_t*) const
{
  // 
  // Print information
  // 
  //   option Not used 
  //
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  TString sNN(AliForwardUtil::CenterOfMassEnergyString(fEnergy));
  sNN.Strip(TString::kBoth, '0');
  sNN.ReplaceAll("GeV", " GeV");
  TString field(AliForwardUtil::MagneticFieldString(fField));
  field.ReplaceAll("p",  "+");
  field.ReplaceAll("m",  "-");
  field.ReplaceAll("kG", " kG");
  
  std::cout << ind << ClassName() << ": " << GetName() << '\n'
	    << ind << " Vertex bins:            " << fVtxAxis.GetNbins() << '\n'
	    << ind << " Vertex range:           [" << fVtxAxis.GetXmin() 
	    << "," << fVtxAxis.GetXmax() << "]\n"
	    << ind << " Low flux cut:           " << fLowFluxCut << '\n'
	    << ind << " Max(delta v_z):         " << fMaxVzErr << " cm\n"
	    << ind << " System:                 " 
	    << AliForwardUtil::CollisionSystemString(fCollisionSystem) << '\n'
	    << ind << " CMS energy per nucleon: " << sNN << '\n'
	    << ind << " Field:                  " <<  field << '\n';
  if (!fCentAxis) { std::cout << std::flush; return; }
  Int_t nBin = fCentAxis->GetNbins();
  std::cout << ind << " Centrality axis:        " << nBin << " bins"
	    << std::flush;
  for (Int_t i = 0; i < nBin; i++) { 
    if ((i % 10) == 0) std::cout << '\n' << ind << "  ";
    std::cout << std::setw(5) << fCentAxis->GetBinLowEdge(i+1) << '-';
  }
  std::cout << std::setw(5) << fCentAxis->GetBinUpEdge(nBin) << std::endl;
}

  
//
// EOF
//

