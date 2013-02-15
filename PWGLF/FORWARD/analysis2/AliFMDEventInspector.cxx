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
#include "AliOADBPhysicsSelection.h"
#include "AliAODForwardMult.h"
#include "AliForwardUtil.h"
#include "AliCentrality.h"
#include <TH1.h>
#include <TList.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TParameter.h>
#include <iostream>
#include <iomanip>
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliVVZERO.h"

//====================================================================
const char* AliFMDEventInspector::fgkFolderName = "fmdEventInspector";

//____________________________________________________________________
AliFMDEventInspector::AliFMDEventInspector()
  : TNamed(),
    fHEventsTr(0), 
    fHEventsTrVtx(0),
    fHEventsAccepted(0),
    fHEventsAcceptedXY(0),
    fHTriggers(0),
    fHTriggerCorr(0),
    fHType(0),
    fHWords(0),
    fHCent(0),
    fHCentVsQual(0),
    fHStatus(0),
    fLowFluxCut(1000),
    fMaxVzErr(0.2),
    fList(0),
    fEnergy(0),
    fField(999), 
    fCollisionSystem(kUnknown),
    fDebug(0),
    fCentAxis(0),
    fVtxAxis(10,-10,10),
    fUseFirstPhysicsVertex(false),
    fUseV0AND(false),
    fMinPileupContrib(3), 
  fMinPileupDistance(0.8),
  fUseDisplacedVertices(false),
  fDisplacedVertex(),
  fCollWords(),
  fBgWords(),
  fCentMethod("V0M"),
  fminCent(-1.0),
  fmaxCent(-1.0)  		
{
  // 
  // Constructor 
  //
  DGUARD(fDebug,1,"Default CTOR of AliFMDEventInspector");
}

//____________________________________________________________________
AliFMDEventInspector::AliFMDEventInspector(const char* name)
  : TNamed(fgkFolderName, name),
    fHEventsTr(0), 
    fHEventsTrVtx(0), 
    fHEventsAccepted(0),
    fHEventsAcceptedXY(0),
    fHTriggers(0),
    fHTriggerCorr(0),
    fHType(0),
    fHWords(0),
    fHCent(0),
    fHCentVsQual(0),
    fHStatus(0),
    fLowFluxCut(1000),
    fMaxVzErr(0.2),
    fList(0),
    fEnergy(0),
    fField(999), 
    fCollisionSystem(kUnknown),
    fDebug(0),
    fCentAxis(0),
    fVtxAxis(10,-10,10),
    fUseFirstPhysicsVertex(false),
    fUseV0AND(false),
    fMinPileupContrib(3), 
    fMinPileupDistance(0.8),
    fUseDisplacedVertices(false),
  fDisplacedVertex(),
  fCollWords(),
  fBgWords(),
  fCentMethod("V0M"),
  fminCent(-1.0),
  fmaxCent(-1.0)  	
{
  // 
  // Constructor 
  // 
  // Parameters:
  //   name Name of object
  //
  DGUARD(fDebug,1,"Named CTOR of AliFMDEventInspector: %s", name);
}

//____________________________________________________________________
AliFMDEventInspector::AliFMDEventInspector(const AliFMDEventInspector& o)
  : TNamed(o), 
    fHEventsTr(o.fHEventsTr), 
    fHEventsTrVtx(o.fHEventsTrVtx), 
    fHEventsAccepted(o.fHEventsAccepted),
    fHEventsAcceptedXY(o.fHEventsAcceptedXY),
    fHTriggers(o.fHTriggers),
    fHTriggerCorr(o.fHTriggerCorr),
    fHType(o.fHType),
    fHWords(o.fHWords),
    fHCent(o.fHCent),
    fHCentVsQual(o.fHCentVsQual),
    fHStatus(o.fHStatus),
    fLowFluxCut(o.fLowFluxCut),
    fMaxVzErr(o.fMaxVzErr),
    fList(o.fList),
    fEnergy(o.fEnergy),
    fField(o.fField), 
    fCollisionSystem(o.fCollisionSystem),
    fDebug(0),
    fCentAxis(0),
    fVtxAxis(o.fVtxAxis),
    fUseFirstPhysicsVertex(o.fUseFirstPhysicsVertex),
    fUseV0AND(o.fUseV0AND),
    fMinPileupContrib(o.fMinPileupContrib), 
    fMinPileupDistance(o.fMinPileupDistance),
    fUseDisplacedVertices(o.fUseDisplacedVertices),
    fDisplacedVertex(o.fDisplacedVertex),
  fCollWords(),
  fBgWords(),
  fCentMethod(o.fCentMethod),
  fminCent(o.fminCent),
  fmaxCent(o.fmaxCent)  	
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //   o Object to copy from 
  //
  DGUARD(fDebug,1,"Copy CTOR of AliFMDEventInspector");
}

//____________________________________________________________________
AliFMDEventInspector::~AliFMDEventInspector()
{
  // 
  // Destructor 
  //
  DGUARD(fDebug,1,"DTOR of AliFMDEventInspector");
  // if (fList)         delete fList;
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
  DGUARD(fDebug,3,"Assignment of AliFMDEventInspector");
  if (&o == this) return *this; 
  TNamed::operator=(o);
  fHEventsTr         = o.fHEventsTr;
  fHEventsTrVtx      = o.fHEventsTrVtx;
  fHEventsAccepted   = o.fHEventsAccepted;
  fHEventsAcceptedXY = o.fHEventsAcceptedXY;
  fHTriggers         = o.fHTriggers;
  fHTriggerCorr      = o.fHTriggerCorr;
  fHType             = o.fHType;
  fHWords            = o.fHWords;
  fHCent             = o.fHCent;
  fHCentVsQual       = o.fHCentVsQual;
  fHStatus           = o.fHStatus;
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
  fUseV0AND              = o.fUseV0AND;
  fMinPileupContrib      = o.fMinPileupContrib;
  fMinPileupDistance     = o.fMinPileupDistance;
  fUseDisplacedVertices  = o.fUseDisplacedVertices;
  fDisplacedVertex       = o.fDisplacedVertex;
  fCentMethod            = o.fCentMethod;
  fminCent		 = o.fminCent;
  fmaxCent		 = o.fmaxCent; 

  if (fList) { 
    fList->SetName(GetName());
    if (fHEventsTr)    fList->Add(fHEventsTr);
    if (fHEventsTrVtx) fList->Add(fHEventsTrVtx);
    if (fHTriggers)    fList->Add(fHTriggers);
    if (fHTriggerCorr) fList->Add(fHTriggerCorr);
    if (fHType)        fList->Add(fHType);
    if (fHWords)       fList->Add(fHWords);
    if (fHCent)        fList->Add(fHCent);
    if (fHCentVsQual)  fList->Add(fHCentVsQual);
    if (fHStatus)      fList->Add(fHStatus);
  }
  return *this;
}

//____________________________________________________________________
void 
AliFMDEventInspector::SetCentralityMethod(ECentMethod m)
{
  switch (m) { 
  case kV0Multiplicity: fCentMethod = "VOM"; break; // VZERO multiplicity 
  case kV0Amplitude:	fCentMethod = "V0A"; break; // VZERO amplitude    
  case kV0Charge: 	fCentMethod = "V0C"; break; // VZERO charge	     
  case kFMDRough: 	fCentMethod = "FMD"; break; // FMD scaled energy l
  case kNTracks: 	fCentMethod = "TRK"; break; // Number of tracks   
  case kLTracks: 	fCentMethod = "TKL"; break; // Number of tracks   
  case kCL0: 		fCentMethod = "CL0"; break; // 		     
  case kCL1: 		fCentMethod = "CL1"; break; // 		     
  case kCND: 		fCentMethod = "CND"; break; // 		     
  case kNParticles:     fCentMethod = "NPA"; break; // Neutral particles  
  case kNeutrons:       fCentMethod = "ZNA"; break; // ZDC neutron amplitu
  case kV0vsFMD: 	fCentMethod = "V0MvsFMD"; break; // VZERO versus FMD   
  case kV0vsNTracks: 	fCentMethod = "TKLvsVOM"; break; // Tracks versus VZERO
  case kZEMvsZDC:	fCentMethod = "ZEMvsZDC"; break; // ZDC		     
  default:              fCentMethod = "V0M"; break;
  }
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::FetchHistograms(const TList* d, 
				      TH1I*& hEventsTr, 
				      TH1I*& hEventsTrVtx, 
				      TH1I*& hEventsAcc,
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
  DGUARD(fDebug,3,"Fetch histograms in AliFMDEventInspector");
  hEventsTr    = 0;
  hEventsTrVtx = 0;
  hEventsAcc   = 0;
  hTriggers    = 0;  
  TList* dd    = dynamic_cast<TList*>(d->FindObject(GetName()));
  if (!dd) return kFALSE;
  
  hEventsTr    = dynamic_cast<TH1I*>(dd->FindObject("nEventsTr"));
  hEventsTrVtx = dynamic_cast<TH1I*>(dd->FindObject("nEventsTrVtx"));
  hEventsAcc   = dynamic_cast<TH1I*>(dd->FindObject("nEventsAccepted"));
  hTriggers    = dynamic_cast<TH1I*>(dd->FindObject("triggers"));

  if (!hEventsTr    || 
      !hEventsTrVtx || 
      !hEventsAcc   ||
      !hTriggers) return kFALSE;
  return kTRUE;
}
//____________________________________________________________________
void
AliFMDEventInspector::CacheConfiguredTriggerClasses(TList& cache, 
						    const TList* classes,
						    AliOADBPhysicsSelection* o)
{
  TIter nextClass(classes);
  TObjString* trigClass = 0;
  // Loop over all trigger classes.  Trigger classes have the format 
  //
  //   class          := positive_words SPACE(s) negative_words 
  //   positive_words := 
  //                  |  '+' words
  //   negative_words := 
  //                  |  '-' words
  //   words          := word 
  //                  |  word ',' words 
  //   
  while ((trigClass = static_cast<TObjString*>(nextClass()))) {
    // Tokenize on space to get positive and negative parts 
    TString     side   = o->GetBeamSide(trigClass->String());
    TObjArray*  parts  = trigClass->String().Tokenize(" ");
    TObjString* part   = 0;
    TIter       nextPart(parts);
    while ((part = static_cast<TObjString*>(nextPart()))) {
      // We only care about the positive ones 
      if (part->GetName()[0] != '+') continue;
      part->String().Remove(0,1);
	
      // Tokenize on a comma to get the words 
      TObjArray*  words = part->String().Tokenize(",");
      TObjString* word  = 0;
      TIter       nextWord(words);
      while ((word = static_cast<TObjString*>(nextWord()))) {
	TNamed* store = new TNamed(word->String(), side);
	cache.Add(store);
	DMSG(fDebug,3,"Caching %s trigger word %s", 
	     store->GetTitle(), store->GetName());
      } // while (word)
      delete words;
    }
    delete parts;
  }
}

//____________________________________________________________________
void
AliFMDEventInspector::SetupForData(const TAxis& vtxAxis)
{
  // 
  // Initialize the object - this is called on the first seen event. 
  // 
  // Parameters:
  //   vtxAxis Vertex axis in use 
  //
  DGUARD(fDebug,1,"Initialize in AliFMDEventInspector");

  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();

  // Get the input handler - should always be there 
  AliInputEventHandler* ih = 
    static_cast<AliInputEventHandler*>(am->GetInputEventHandler());
  if (!ih) { 
    AliWarning("No input handler");
    return;
  }
  // Get the physics selection - should always be there 
  AliPhysicsSelection* ps = 
    static_cast<AliPhysicsSelection*>(ih->GetEventSelection());
  if (!ps) {
    AliWarning("No physics selection");
    return;
  }
  // Get the configured triggers
  AliOADBPhysicsSelection* oadb = 
    const_cast<AliOADBPhysicsSelection*>(ps->GetOADBPhysicsSelection());
  if (!oadb) {
    AliWarning("No OADB physics selection object");
    return;
  }
  // Get the configured trigger words from the physics selection 
  const TList* collTriggClasses = ps->GetCollisionTriggerClasses();
  const TList* bgTriggClasses   = ps->GetBGTriggerClasses();
  if (!collTriggClasses) { 
    AliWarning("No configured collision trigger classes");
    return;
  }
  if (!bgTriggClasses) { 
    AliWarning("No configured background trigger classes");
    return;
  }
  CacheConfiguredTriggerClasses(fCollWords, collTriggClasses, oadb);
  CacheConfiguredTriggerClasses(fBgWords,   bgTriggClasses,   oadb);
  // fCollWords.ls();
  // fBgWords.ls();
  
  
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

  fHTriggerCorr = new TH2I("triggerCorr", "Trigger correlation", 
			   kOffline+1, 0, kOffline+1, 
			   kOffline+1, 0, kOffline+1);
  fHTriggerCorr->SetStats(0);
  fHTriggerCorr->SetDirectory(0);
  fHTriggerCorr->SetXTitle("Requirement");
  fHTriggerCorr->SetYTitle("Companion");

  Int_t binNum[] = { kInel   +1,
		     kInelGt0+1,
		     kNSD    +1,
		     kV0AND  +1,
		     kEmpty  +1,
		     kA      +1,
		     kB      +1,
		     kC      +1,
		     kE      +1,
		     kPileUp +1,
		     kMCNSD  +1,
		     kOffline+1 };
  const char* binLbl[] = { "INEL",	 
			   "INEL>0",
			   "NSD",	 
			   "VOAND",
			   "Empty",	 
			   "A",	 
			   "B",	 
			   "C",	 
			   "E",	 
			   "Pileup",
			   "NSD_{MC}", 
			   "Offline" };
  for (Int_t i = 0; i < kOffline+1; i++) {
    fHTriggers->GetXaxis()->SetBinLabel(binNum[i], binLbl[i]);
    fHTriggerCorr->GetXaxis()->SetBinLabel(binNum[i], binLbl[i]);
    fHTriggerCorr->GetYaxis()->SetBinLabel(binNum[i], binLbl[i]);
  }
  fList->Add(fHTriggers);
  fList->Add(fHTriggerCorr);
  

  fHType = new TH1I("type", Form("Event type (cut: SPD mult>%d)", 
				 fLowFluxCut), 2, -.5, 1.5);
  fHType->SetFillColor(kRed+1);
  fHType->SetFillStyle(3001);
  fHType->SetStats(0);
  fHType->SetDirectory(0);
  fHType->GetXaxis()->SetBinLabel(1,"Low-flux");
  fHType->GetXaxis()->SetBinLabel(2,"High-flux");
  fList->Add(fHType);

#if 0 
  // This histogram disabled as it causes problems in the merge 
  fHWords = new TH1I("words", "Trigger words seen", 1, 0, 0); 
  fHWords->SetFillColor(kBlue+1);
  fHWords->SetFillStyle(3001);
  fHWords->SetStats(0);
  fHWords->SetDirectory(0);
  fHWords->SetBit(TH1::kCanRebin);
  fList->Add(fHWords);
#endif

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
  fHCentVsQual->SetDirectory(0);
  fList->Add(fHCentVsQual);

  fHStatus = new TH1I("status", "Status", 7, 1, 8);
  fHStatus->SetFillColor(kRed+1);
  fHStatus->SetFillStyle(3001);
  fHStatus->SetStats(0);
  fHStatus->SetDirectory(0);
  fHStatus->GetXaxis()->SetBinLabel(1, "OK");
  fHStatus->GetXaxis()->SetBinLabel(2, "No event");
  fHStatus->GetXaxis()->SetBinLabel(3, "No triggers");
  fHStatus->GetXaxis()->SetBinLabel(4, "No SPD");
  fHStatus->GetXaxis()->SetBinLabel(5, "No FMD");
  fHStatus->GetXaxis()->SetBinLabel(6, "No vertex");
  fHStatus->GetXaxis()->SetBinLabel(7, "Bad vertex");
  fList->Add(fHStatus);
}

//____________________________________________________________________
void
AliFMDEventInspector::StoreInformation(Int_t runNo)
{
  // Write TNamed objects to output list containing information about
  // the running conditions 
  DGUARD(fDebug,2,"Store information from AliFMDEventInspector");
  if (!fList) return;

  
  fList->Add(AliForwardUtil::MakeParameter("sys", fCollisionSystem));
  fList->Add(AliForwardUtil::MakeParameter("sNN", fEnergy));
  fList->Add(AliForwardUtil::MakeParameter("field", fField));
  fList->Add(AliForwardUtil::MakeParameter("runNo", runNo));
  fList->Add(AliForwardUtil::MakeParameter("lowFlux", fLowFluxCut));
  fList->Add(AliForwardUtil::MakeParameter("fpVtx",fUseFirstPhysicsVertex));
  fList->Add(AliForwardUtil::MakeParameter("v0and",fUseV0AND));
  fList->Add(AliForwardUtil::MakeParameter("nPileUp", fMinPileupContrib));
  fList->Add(AliForwardUtil::MakeParameter("dPileup", fMinPileupDistance));
  fList->Add(AliForwardUtil::MakeParameter("alirootRev", 
					   AliForwardUtil::AliROOTRevision()));
  fList->Add(AliForwardUtil::MakeParameter("alirootBranch", 
					   AliForwardUtil::AliROOTBranch()));

}

//____________________________________________________________________
void
AliFMDEventInspector::CreateOutputObjects(TList* dir)
{
  // 
  // Define the output histograms.  These are put in a sub list of the
  // passed list.   The histograms are merged before the parent task calls 
  // AliAnalysisTaskSE::Terminate 
  // 
  //   dir Directory to add to 
  //
  DGUARD(fDebug,1,"Define output from AliFMDEventInspector");
  fList = new TList;
  fList->SetName(GetName());
  fList->SetOwner();
  dir->Add(fList);
}

//____________________________________________________________________
UInt_t
AliFMDEventInspector::Process(const AliESDEvent* event, 
			      UInt_t&            triggers,
			      Bool_t&            lowFlux,
			      UShort_t&          ivz, 
			      TVector3&          ip,
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
  DGUARD(fDebug,1,"Process event in AliFMDEventInspector"); 
  // --- Check that we have an event ---------------------------------
  if (!event) { 
    AliWarning("No ESD event found for input event");
    fHStatus->Fill(2);
    return kNoEvent;
  }

  // --- Read trigger information from the ESD and store in AOD object
  if (!ReadTriggers(*event, triggers, nClusters)) { 
    if (fDebug > 2) {
      AliWarning("Failed to read triggers from ESD"); }
    fHStatus->Fill(3);
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

  // --- Process satellite event information is requested ------------
  if (fUseDisplacedVertices) { 
    if (!fDisplacedVertex.Process(event)) 
      AliWarning("Failed to process satellite event");
  }
  
  // --- Read centrality information 
  cent          = -10;
  UShort_t qual = 0;
  if (!ReadCentrality(*event, cent, qual)) {
    if (fDebug > 3) 
      AliWarning("Failed to get centrality");
  }
//check centrality cut
 
  if(fminCent>0.0&&cent<fminCent)
	return  kNoEvent; 
  if(fmaxCent>0.0&&cent>fmaxCent)
	 return  kNoEvent; 
  fHCent->Fill(cent);
  if (qual == 0) fHCentVsQual->Fill(0., cent);
  else { 
    for (UShort_t i = 0; i < 4; i++) 
      if (qual & (1 << i)) fHCentVsQual->Fill(Double_t(i+1), cent);
  }

  // --- Get the interaction point -----------------------------------
  Bool_t vzOk = ReadVertex(*event, ip);
  fHEventsTr->Fill(ip.Z());
  if (!vzOk) { 
    if (fDebug > 3) {
      AliWarning("Failed to read vertex from ESD"); }
    fHStatus->Fill(6);
    return kNoVertex;
  }
  fHEventsTrVtx->Fill(ip.Z());
  
  // --- Get the vertex bin ------------------------------------------
  ivz = fVtxAxis.FindBin(ip.Z());
  if (ivz <= 0 || ivz > fVtxAxis.GetNbins()) { 
    if (fDebug > 3) {
      AliWarning(Form("Vertex @ %f outside of range [%f,%f]", 
		      ip.Z(), fVtxAxis.GetXmin(), fVtxAxis.GetXmax())); 
    }
    ivz = 0;
    fHStatus->Fill(7);
    return kBadVertex;
  }
  fHEventsAccepted->Fill(ip.Z());
  fHEventsAcceptedXY->Fill(ip.X(),ip.Y());
  
  // --- Check the FMD ESD data --------------------------------------
  if (!event->GetFMDData()) { 
    if (fDebug > 3) {
      AliWarning("No FMD data found in ESD"); }
    fHStatus->Fill(5);
    return kNoFMD;
  }

  fHStatus->Fill(1);
  return kOk;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::ReadCentrality(const AliESDEvent& esd, 
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
  DGUARD(fDebug,2,"Read the centrality in AliFMDEventInspector");

  if(fUseDisplacedVertices) {
    Double_t zvtx = fDisplacedVertex.GetVertexZ();
    qual          = 1;
    if(TMath::Abs(zvtx) < 999) {
      cent = fDisplacedVertex.GetCentralityPercentile();
      qual = 0;
    }
    return true;
  }
  
  cent = -1;
  qual = 0;
  AliCentrality* centObj = const_cast<AliESDEvent&>(esd).GetCentrality();
  if (!centObj)  return true;

  cent = centObj->GetCentralityPercentile(fCentMethod);  
  qual = centObj->GetQuality();

  return true;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::CheckpAExtraV0(const AliESDEvent& esd) const
{
  if (fCollisionSystem != AliForwardUtil::kPPb) return true;

   AliVVZERO* esdV0 = esd.GetVZEROData();
   if ((esdV0->GetV0ADecision()!=1) || (esdV0->GetV0CDecision()!=1)) 
     return false;
   return true;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::ReadTriggers(const AliESDEvent& esd, UInt_t& triggers,
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
  DGUARD(fDebug,2,"Read the triggers in AliFMDEventInspector");
  triggers = 0;

  // Get the analysis manager - should always be there 
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  DMSG(fDebug,10,"Got analysis manager %p", am);
  if (!am) { 
    AliWarning("No analysis manager defined!");
    return kFALSE;
  }

  // Get the input handler - should always be there 
  AliInputEventHandler* ih = 
    static_cast<AliInputEventHandler*>(am->GetInputEventHandler());
  DMSG(fDebug,10,"Got input handler %p", ih);
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
  Bool_t  offline  = ih->IsEventSelected();
  Bool_t  fastonly = (ih->IsEventSelected() & AliVEvent::kFastOnly);
  TString trigStr  = esd.GetFiredTriggerClasses();

  if (fHWords) fHWords->Fill(trigStr.Data(), 1);
  
  if(fUseDisplacedVertices) {
    DMSG(fDebug,3,"Using displaced vertex stuff");
    if (TMath::Abs(fDisplacedVertex.GetVertexZ()) >= 999) offline = false;
  }
  
  if (CheckFastPartition(fastonly))     offline = false;
  if (offline && CheckCosmics(trigStr)) offline = false;
  // if (!CheckpAExtraV0(esd))             offline = false;

  DMSG(fDebug,2,"Event is %striggered by off-line", offline ? "" : "NOT ");

  if (offline) {
    triggers |= AliAODForwardMult::kOffline;
    triggers |= AliAODForwardMult::kInel;
    if (!fHTriggers) { 
      AliWarning("Histogram of triggers not defined - has init been called");
      return false;
    }
    // fHTriggers->Fill(kOffline+0.5);
    
    CheckINELGT0(esd, nClusters, triggers);
  }
  
  CheckNSD(esd,triggers);
  CheckPileup(esd, triggers);
  CheckEmpty(trigStr, triggers);
  // if (CheckPileup(esd, triggers)) fHTriggers->Fill(kPileUp+.5);
  // if (CheckEmpty(trigStr, triggers)) fHTriggers->Fill(kEmpty+.5);

  CheckWords(esd, triggers);

#if 0
  // Now check - if we have a collision - for offline triggers and
  // fill histogram.
  if (triggers & AliAODForwardMult::kB) {
    fHTriggers->Fill(kB+.5);
    if (triggers & AliAODForwardMult::kInel) 
      fHTriggers->Fill(kInel+.5);
    
    if (triggers & AliAODForwardMult::kInelGt0)
      fHTriggers->Fill(kInelGt0+.5);
    
    if (triggers & AliAODForwardMult::kNSD)
      fHTriggers->Fill(kNSD+.5);

    if (triggers & AliAODForwardMult::kV0AND)
      fHTriggers->Fill(kV0AND+.5);
  }
  if (triggers & AliAODForwardMult::kA) fHTriggers->Fill(kA+.5);
  if (triggers & AliAODForwardMult::kC) fHTriggers->Fill(kC+.5);
  if (triggers & AliAODForwardMult::kE) fHTriggers->Fill(kE+.5);
#endif
#define TEST_TRIG_BIN(RET,BIN,TRIGGERS) \
  do { switch (BIN) { \
    case kInel:    RET = triggers & AliAODForwardMult::kInel;    break; \
    case kInelGt0: RET = triggers & AliAODForwardMult::kInelGt0; break; \
    case kNSD:     RET = triggers & AliAODForwardMult::kNSD;     break; \
    case kV0AND:   RET = triggers & AliAODForwardMult::kV0AND;   break; \
    case kEmpty:   RET = triggers & AliAODForwardMult::kEmpty;   break; \
    case kA:       RET = triggers & AliAODForwardMult::kA;       break; \
    case kB:       RET = triggers & AliAODForwardMult::kB;       break; \
    case kC:       RET = triggers & AliAODForwardMult::kC;       break; \
    case kE:       RET = triggers & AliAODForwardMult::kE;       break; \
    case kPileUp:  RET = triggers & AliAODForwardMult::kPileUp;  break; \
    case kMCNSD:   RET = triggers & AliAODForwardMult::kMCNSD;   break; \
    case kOffline: RET = triggers & AliAODForwardMult::kOffline; break; \
    default:       RET = false; } } while(false)
      
  
  for (Int_t i = 0; i < kOffline+1; i++) { 
    Bool_t hasX = false;
    TEST_TRIG_BIN(hasX, i, triggers);
    if (!hasX) continue;
    fHTriggers->Fill(i+.5);
    for (Int_t j = 0; j < kOffline+1; j++) { 
      Bool_t hasY = false;
      TEST_TRIG_BIN(hasY, j, triggers);
      if (!hasY) continue;
      
      fHTriggerCorr->Fill(i+.5, j+.5);
    }
  }
  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::CheckFastPartition(bool fastonly) const
{
  // For the 2.76 TeV p+p run, the FMD ran in the slow partition 
  // so it received no triggers from the fast partition. Therefore
  // the fast triggers are removed here but not for MC where all 
  // triggers are fast.
  if (TMath::Abs(fEnergy - 2750.) > 20) return false;
  if (fCollisionSystem != AliForwardUtil::kPP) return false;
  if (fastonly)
    DMSG(fDebug,1,"Fast trigger in pp @ sqrt(s)=2.76TeV removed");

  return fastonly;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::CheckCosmics(const TString& trigStr) const
{
  // MUON triggers are not strictly minimum bias (MB) so they are
  // removed (HHD)
  if(trigStr.Contains("CMUS1")) {
    DMSG(fDebug,1,"Cosmic trigger ins't min-bias, removed");
    return true;
  }
  return false;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::CheckINELGT0(const AliESDEvent& esd, 
				   UShort_t& nClusters, 
				   UInt_t& triggers) const
{
  nClusters = 0;

  // If this is inel, see if we have a tracklet 
  const AliMultiplicity* spdmult = esd.GetMultiplicity();
  if (!spdmult) {
    AliWarning("No SPD multiplicity");
    return false;
  }

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
  if (nClusters > 0) triggers |= AliAODForwardMult::kNClusterGt0;

  return triggers & AliAODForwardMult::kNClusterGt0;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::CheckNSD(const AliESDEvent& esd, UInt_t& triggers) const
{
  // Analyse some trigger stuff 
  AliTriggerAnalysis ta;
  if (ta.IsOfflineTriggerFired(&esd, AliTriggerAnalysis::kV0AND)) {
    triggers |= AliAODForwardMult::kV0AND;
    if (fUseV0AND) 
      triggers |= AliAODForwardMult::kNSD;
  }
  if (ta.IsOfflineTriggerFired(&esd, AliTriggerAnalysis::kNSD1)) 
    triggers |= AliAODForwardMult::kNSD;
  return triggers & AliAODForwardMult::kNSD;
}
//____________________________________________________________________
Bool_t
AliFMDEventInspector::CheckPileup(const AliESDEvent& esd, 
				  UInt_t& triggers) const
{
  // Check for multiple vertices (pile-up) with at least 3
  // contributors and at least 0.8cm from the primary vertex
  if(fCollisionSystem != AliForwardUtil::kPP) return false;

  Bool_t pileup =  esd.IsPileupFromSPD(fMinPileupContrib,fMinPileupDistance);
  if (pileup) triggers |= AliAODForwardMult::kPileUp;
  return pileup;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::CheckEmpty(const TString& trigStr, UInt_t& triggers) const
{
  if (trigStr.Contains("CBEAMB-ABCE-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kEmpty;
    return true;
  }
  return false;
}
//____________________________________________________________________
Bool_t
AliFMDEventInspector::CheckWords(const AliESDEvent& esd, UInt_t& triggers) const
{
  TObject* word = 0;
  TIter nextColl(&fCollWords);
  while ((word = nextColl())) {
    DMSG(fDebug,10,"Checking if %s trigger %s is fired", 
	 word->GetTitle(), word->GetName());
    if (!esd.IsTriggerClassFired(word->GetName())) continue;

    TString beamSide = word->GetTitle();
    DMSG(fDebug,10,"Found it - this is a %s trigger", beamSide.Data());

    if (!beamSide.EqualTo("B")) continue;
    triggers |= AliAODForwardMult::kB;
    break; // No more to do here
  }
  TIter nextBg(&fBgWords);
  UInt_t all = (AliAODForwardMult::kA | 
		AliAODForwardMult::kC | 
		AliAODForwardMult::kE);
  while ((word = nextBg())) {
    DMSG(fDebug,10,"Checking if %s trigger %s is fired", 
	 word->GetTitle(), word->GetName());
    if (!esd.IsTriggerClassFired(word->GetName())) continue;
    
    TString beamSide = word->GetTitle();
    DMSG(fDebug,10,"Found it - this is a %s trigger", beamSide.Data());

    if (beamSide.Contains("A")) triggers |= AliAODForwardMult::kA;
    if (beamSide.Contains("C")) triggers |= AliAODForwardMult::kC;
    if (beamSide.Contains("E")) triggers |= AliAODForwardMult::kE;

    if ((triggers & all) == all) break; // No more to do
  }
  return true;
}


//____________________________________________________________________
Bool_t
AliFMDEventInspector::ReadVertex(const AliESDEvent& esd, TVector3& ip)
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
  DGUARD(fDebug,2,"Read the vertex in AliFMDEventInspector");
  ip.SetXYZ(1024, 1024, 0);
  
  if(fUseDisplacedVertices) {
    Double_t zvtx = fDisplacedVertex.GetVertexZ();
      
    if(TMath::Abs(zvtx) < 999) {
      ip.SetZ(zvtx);
      return true;
    }
    return false;
  }

  if(fUseFirstPhysicsVertex) return CheckPWGUDVertex(esd, ip);
  
  
  return CheckVertex(esd, ip);
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::CheckPWGUDVertex(const AliESDEvent& esd, 
				       TVector3& ip)  const
{
  // This is the code used by the 1st physics people 
  const AliESDVertex* vertex    = esd.GetPrimaryVertex();
  if (!vertex  || !vertex->GetStatus()) {
    DMSG(fDebug,2,"No primary vertex (%p) or bad status %d", 
	 vertex, (vertex ? vertex->GetStatus() : -1));
    return false;
  }
  const AliESDVertex* vertexSPD = esd.GetPrimaryVertexSPD();
  if (!vertexSPD || !vertexSPD->GetStatus()) {
    DMSG(fDebug,2,"No primary SPD vertex (%p) or bad status %d", 
	 vertexSPD, (vertexSPD ? vertexSPD->GetStatus() : -1));
    return false;
  }
    
  // if vertex is from SPD vertexZ, require more stringent cuts 
  if (vertex->IsFromVertexerZ()) { 
    if (vertex->GetDispersion() > fMaxVzErr || 
	vertex->GetZRes() > 1.25 * fMaxVzErr) {
      DMSG(fDebug,2,"Dispersion %f > %f or resolution %f > %f",
	   vertex->GetDispersion(), fMaxVzErr,
	   vertex->GetZRes(), 1.25 * fMaxVzErr);
      return false;
    }
  }
  ip.SetZ(vertex->GetZ());
  
  if(!vertex->IsFromVertexerZ()) {
    ip.SetX(vertex->GetX());
    ip.SetY(vertex->GetY());
  }
  return true;
}
//____________________________________________________________________
Bool_t
AliFMDEventInspector::CheckVertex(const AliESDEvent& esd, 
				  TVector3& ip) const
{
  // Use standard SPD vertex (perhaps preferable for Pb+Pb)
  // Get the vertex 
  const AliESDVertex* vertex = esd.GetPrimaryVertexSPD();
  if (!vertex) { 
    if (fDebug > 2) {
      AliWarning("No SPD vertex found in ESD"); }
    return false;
  }
    
  // Check that enough tracklets contributed 
  if(vertex->GetNContributors() <= 0) {
    DMSG(fDebug,2,"Number of contributors to vertex is %d<=0",
	 vertex->GetNContributors());
    ip.SetZ(0);
    return false;
  } 
  // Check that the uncertainty isn't too large 
  if (vertex->GetZRes() > fMaxVzErr) { 
    DMSG(fDebug,2,"Uncertaintity in Z of vertex is too large %f > %f", 
	 vertex->GetZRes(), fMaxVzErr);
    return false;
  }
    
  // Get the z coordiante 
  ip.SetZ(vertex->GetZ());
  const AliESDVertex* vertexXY = esd.GetPrimaryVertex();
    
  
  if(!vertexXY->IsFromVertexerZ()) {
    ip.SetX(vertexXY->GetX());
    ip.SetY(vertexXY->GetY());
  }
  return true;
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
  DGUARD(fDebug,2,"Read the run details in AliFMDEventInspector");
  const char* sys  = esd->GetBeamType();
  Float_t     cms  = 2 * esd->GetBeamEnergy();
  Float_t     fld  = esd->GetMagneticField();
  fCollisionSystem = AliForwardUtil::ParseCollisionSystem(sys);
  fEnergy          = AliForwardUtil::ParseCenterOfMassEnergy(fCollisionSystem, 
							     cms);
  fField           = AliForwardUtil::ParseMagneticField(fld);

  StoreInformation(esd->GetRunNumber());
  if (fCollisionSystem   == AliForwardUtil::kUnknown) { 
    AliWarningF("Unknown collision system: %s - please check", sys);
    return false;
  }
  if (fEnergy            <= 0) {
    AliWarningF("Unknown CMS energy: %f (%d) - please check", cms, fEnergy);
    return false;
  }
  if (TMath::Abs(fField) >  10) {
    AliWarningF("Unknown L3 field setting: %f (%d) - please check", fld,fField);
    return false;
  }

  return true;
}


//____________________________________________________________________
const Char_t*
AliFMDEventInspector::CodeString(UInt_t code)
{
  static TString s;
  s = "";
  if (code & kNoEvent)    s.Append("NOEVENT ");
  if (code & kNoTriggers) s.Append("NOTRIGGERS ");
  if (code & kNoSPD)      s.Append("NOSPD ");
  if (code & kNoFMD)      s.Append("NOFMD ");
  if (code & kNoVertex)   s.Append("NOVERTEX ");
  if (code & kBadVertex)  s.Append("BADVERTEX ");
  return s.Data();
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
  
  std::cout << std::boolalpha 
	    << ind << ClassName() << ": " << GetName() << '\n'
	    << ind << " Vertex bins:            " << fVtxAxis.GetNbins() << '\n'
	    << ind << " Vertex range:           [" << fVtxAxis.GetXmin() 
	    << "," << fVtxAxis.GetXmax() << "]\n"
	    << ind << " Low flux cut:           " << fLowFluxCut << '\n'
	    << ind << " Max(delta v_z):         " << fMaxVzErr << " cm\n"
	    << ind << " Min(nContrib_pileup):   " << fMinPileupContrib << '\n'
	    << ind << " Min(v-pileup):          " << fMinPileupDistance << '\n'
	    << ind << " System:                 " 
	    << AliForwardUtil::CollisionSystemString(fCollisionSystem) << '\n'
	    << ind << " CMS energy per nucleon: " << sNN << '\n'
	    << ind << " Field:                  " <<  field << '\n'
	    << ind << " Satellite events:       " << fUseDisplacedVertices<<'\n'
	    << ind << " Centrality method:      " << fCentMethod << '\n'
	    << std::noboolalpha;
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

