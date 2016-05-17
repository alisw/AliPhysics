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
#include "AliProdInfo.h"
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
#include "AliVVZERO.h"
#include <TH1.h>
#include <TList.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TParameter.h>
#include <TMatrixDSym.h>
#include <TPRegexp.h>
#include <iostream>
#include <iomanip>
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliVVZERO.h"

//====================================================================
namespace {
  const UInt_t kTriggerMask = AliAODForwardMult::kOffline;
    // (AliAODForwardMult::kInel      |
    //  AliAODForwardMult::kV0AND     |
    //  AliAODForwardMult::kADOR      |
    //  AliAODForwardMult::kADAND     |
    //  AliAODForwardMult::kSatellite |
    //  AliAODForwardMult::kOffline   |
    //  AliAODForwardMult::kInclusive);


  AliPhysicsSelection* GetPhysicsSelection()
  {
    AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
    if (!am) {
      Warning("GetPhysicsSelection", "No analysis manager");
      return 0;
    }
    
    // Get the input handler - should always be there 
    AliInputEventHandler* ih = 
      dynamic_cast<AliInputEventHandler*>(am->GetInputEventHandler());
    if (!ih) { 
      Warning("GetPhysicsSelection", "No input handler");
      return 0;
    }
    // Get the physics selection - should always be there 
    AliPhysicsSelection* ps = 
      dynamic_cast<AliPhysicsSelection*>(ih->GetEventSelection());
    if (!ps) {
      Warning("GetPhysicsSelection", "No physics selection");
      return 0;
    }
    return ps;
  }
}

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
    fHType(0),
    fHWords(0),
    fHCent(0),
    fHCentVsQual(0),
    fHStatus(0),
    fHVtxStatus(0),
    fHTrgStatus(0), 
    fHPileup(0),
    fLowFluxCut(1000),
    fMaxVzErr(0.2),
    fList(0),
    fEnergy(0),
    fField(999), 
    fCollisionSystem(kUnknown),
    fDebug(0),
    fCentAxis(0),
    fVtxAxis(10,-10,10),
    fVtxMethod(kNormal),
    // fUseFirstPhysicsVertex(false),
    fUseV0AND(false),
    fPileupFlags(kSPD|kTracks|kOutOfBunch|kSPDBins),
    fMinPileupContrib(3), 
    fMinPileupDistance(0.8),
    // fUseDisplacedVertices(false),
    fDisplacedVertex(),
    fCollWords(),
    fBgWords(),
    fCentMethod("default"),
    fMinCent(-1.0),
    fMaxCent(-1.0),
    // fUsepA2012Vertex(false),
    fRunNumber(0),
    fMC(false),
    fProdYear(-1),
    fProdLetter('?'),
    fProdPass(-1),
    fProdSVN(-1),
    fProdMC(false)
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
    fHType(0),
    fHWords(0),
    fHCent(0),
    fHCentVsQual(0),
    fHStatus(0),
    fHVtxStatus(0),
    fHTrgStatus(0),
    fHPileup(0),
    fLowFluxCut(1000),
    fMaxVzErr(0.2),
    fList(0),
    fEnergy(0),
    fField(999), 
    fCollisionSystem(kUnknown),
    fDebug(0),
    fCentAxis(0),
    fVtxAxis(10,-10,10),
    fVtxMethod(kNormal),
// fUseFirstPhysicsVertex(false),
    fUseV0AND(false),
    fPileupFlags(kSPD|kTracks|kOutOfBunch|kSPDBins),
    fMinPileupContrib(3), 
    fMinPileupDistance(0.8),
// fUseDisplacedVertices(false),
  fDisplacedVertex(),
  fCollWords(),
  fBgWords(),
  fCentMethod("default"),
  fMinCent(-1.0),
  fMaxCent(-1.0),
// fUsepA2012Vertex(false),
  fRunNumber(0),
  fMC(false),
  fProdYear(-1),
  fProdLetter('?'),
  fProdPass(-1),
  fProdSVN(-1),
  fProdMC(false)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //   name Name of object
  //
  DGUARD(fDebug,1,"Named CTOR of AliFMDEventInspector: %s", name);
  if (!GetPhysicsSelection()) {
    AliWarning("No physics selection defined in manager\n"
	       "=======================================================\n"
	       " The use of AliFMDEventInspector _suggests_ a valid\n"
	       " AliPhysicsSelection registered with the input handler.\n\n"
	       " Please check our train setup\n"
	       "=======================================================\n");
  }
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
  default:              fCentMethod = "default"; break;
  }
}

//____________________________________________________________________
void 
AliFMDEventInspector::SetMinCentrality(Double_t minCent)
{
  AliWarning("\n"
	     "*******************************************************\n"
	     "* Setting centrality cuts in this stage is deprecated *\n"
	     "*******************************************************");
  fMinCent = minCent;
}
//____________________________________________________________________
void 
AliFMDEventInspector::SetMaxCentrality(Double_t maxCent)
{
  AliWarning("\n"
	     "*******************************************************\n"
	     "* Setting centrality cuts in this stage is deprecated *\n"
	     "*******************************************************");
  fMaxCent = maxCent;
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
  Printf("Inspecting physics selection: %p", o);
  o->Print();
  while ((trigClass = static_cast<TObjString*>(nextClass()))) {
    // Tokenize on space to get positive and negative parts 
    TString     side   = o->GetBeamSide(trigClass->String());
    TObjArray*  parts  = trigClass->String().Tokenize(" ");
    TObjString* part   = 0;
    TIter       nextPart(parts);
    while ((part = static_cast<TObjString*>(nextPart()))) {
      // We only care about the positive ones 
      // if (part->GetName()[0] == '-') continue;
      // part->String().Remove(0,1);
	
      // Tokenize on a comma to get the words 
      TObjArray*  words = part->String().Tokenize(",");
      TObjString* word  = 0;
      TIter       nextWord(words);
      while ((word = static_cast<TObjString*>(nextWord()))) {
	char first = word->String()[0];
	switch (first) {
	case '-': // Negatives 
	case '&': // bunch crossing
	case '*': // ?
	  continue;
	  break;
	case '+':
	  word->String().Remove(0,1);
	  break;
	}
	if (cache.FindObject(word->String())) {
	  DMSG(fDebug, 3, "Word %s already cached", word->GetName());
	  continue;
	}
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

  // Get the physics selection - should always be there 
  AliPhysicsSelection* ps = GetPhysicsSelection();
  if (!ps) {
    AliWarning("No physics selection");
  }
  else {
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
#if 1 
    // This histogram disabled as it causes problems in the merge 
    fHWords = new TH1I("words", "Trigger words seen", 1, 0, 0); 
    fHWords->SetFillColor(kBlue-2);
    fHWords->SetFillStyle(3001);
    fHWords->SetStats(0);
    fHWords->SetDirectory(0);
#endif
    
    CacheConfiguredTriggerClasses(fCollWords, collTriggClasses, oadb);
    CacheConfiguredTriggerClasses(fBgWords,   bgTriggClasses,   oadb);
    fCollWords.Sort();
    fBgWords.Sort();
    DMSG(fDebug, 3, "fHWords=%p", fHWords);
    if (fHWords) {
      Int_t  nColl = fCollWords.GetEntries();
      Int_t  nBg   = fBgWords.GetEntries();
      TAxis* xAxis = fHWords->GetXaxis();

      DMSG(fDebug,3, "Got a total of %d+%d=%d trigger words",
	   nColl, nBg, nColl+nBg); 
      xAxis->Set(nColl+nBg, 0.5, nColl+nBg+.5);
      TIter nextC(&fCollWords);
      TIter nextB(&fBgWords);
      Int_t bin = 1;
      TObject* word = 0;
      while((word = nextC())) {
	DMSG(fDebug, 3, "Defining bin %d to be %s", bin, word->GetName());
	word->SetUniqueID(bin);
	xAxis->SetBinLabel(bin++, word->GetName());
      }
      while((word = nextB())) {
	DMSG(fDebug, 3, "Defining bin %d to be %s", bin, word->GetName());
	word->SetUniqueID(bin);	
	xAxis->SetBinLabel(bin++, word->GetName());
      }
      fHWords->Rebuild();
    }
    if (fDebug >= 3) {
      DMSG(fDebug, 3, "Collision trigger classes");
      fCollWords.ls();
      DMSG(fDebug, 3, "Background trigger classes");
      fBgWords.ls();
    }
  }
  
  TArrayD limits;
  if ((fMinCent < 0 && fMaxCent < 0) || fMaxCent <= fMinCent) {
    // Was:
    // -1.5 -0.5 0.5 1.5 ... 89.5 ... 100.5
    // ----- 92 number --------- ---- 1 ---
    // Now:
    // -0.5 0.5 1.5 ... 89.5 ... 100.5
    // ----- 91 number ---- ---- 1 ---
    limits.Set(92);
    for (Int_t i = 0; i < 91; i++) limits[i] = -.5 + i;
    limits[91] = 100.5;
  }
  else {
    Int_t n = Int_t(fMaxCent-fMinCent+2);
    limits.Set(n);
    for (Int_t i = 0; i < n; i++) { 
      limits[i] = fMinCent + i - .5;
    }
  }

  if (vtxAxis.GetXbins() && vtxAxis.GetXbins()->GetArray()) 
    fVtxAxis.Set(vtxAxis.GetNbins(), vtxAxis.GetXbins()->GetArray());
  else 
    fVtxAxis.Set(vtxAxis.GetNbins(), vtxAxis.GetXmin(), vtxAxis.GetXmax());
  
  fCentAxis  = new TAxis(limits.GetSize()-1, limits.GetArray());
  fHEventsTr = 0;
  if (vtxAxis.GetXbins() && vtxAxis.GetXbins()->GetArray()) 
    fHEventsTr = new TH1I("nEventsTr", "Number of events w/trigger",
			  vtxAxis.GetNbins(), vtxAxis.GetXbins()->GetArray());
  else 
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

  fHEventsAccepted = 0;
  if (vtxAxis.GetXbins() && vtxAxis.GetXbins()->GetArray()) 
    fHEventsAccepted =
      new TH1I("nEventsAccepted", 
	       "Number of events w/trigger and vertex in range",
	       vtxAxis.GetNbins(), vtxAxis.GetXbins()->GetArray());
  else 
    fHEventsAccepted =
      new TH1I("nEventsAccepted", 
	       "Number of events w/trigger and vertex in range",
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

  fHTriggers = AliAODForwardMult::MakeTriggerHistogram("triggers",kTriggerMask);
  fHTriggers->SetStats(0);
  fHTriggers->SetDirectory(0);
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

  if (fHWords) fList->Add(fHWords);

  fHCent = new TH1F("cent", "Centrality", limits.GetSize()-1,limits.GetArray());
  fHCent->SetFillColor(kBlue+1);
  fHCent->SetFillStyle(3001);
  fHCent->SetStats(0);
  fHCent->SetDirectory(0);
  fHCent->SetXTitle("Centrality [%]");
  fHCent->SetYTitle("Events");
  fList->Add(fHCent);

  fHCentVsQual = new TH2F("centVsQuality", "Quality vs Centrality", 
			  9, 196.5, 205.5,
			  limits.GetSize()-1, limits.GetArray());
  fHCentVsQual->SetXTitle("Quality");
  fHCentVsQual->SetYTitle("Centrality [%]");
  fHCentVsQual->SetZTitle("Events");
  TAxis* qAxis = fHCentVsQual->GetXaxis();
  qAxis->SetBinLabel(1, "OK");
  qAxis->SetBinLabel(2, "198: beyond anchor");
  qAxis->SetBinLabel(3, "199: no calib");
  qAxis->SetBinLabel(4, "200: not desired trigger");
  qAxis->SetBinLabel(5, "201: not INEL0 with tracklets");
  qAxis->SetBinLabel(6, "202: vertex Z not within 10cm");
  qAxis->SetBinLabel(7, "203: tagged as pileup (SPD)");
  qAxis->SetBinLabel(8, "204: inconsistent SPD/tracking vertex");
  qAxis->SetBinLabel(9, "205: rejected by tracklets-vs-clusters");

  fHCentVsQual->SetDirectory(0);
  fList->Add(fHCentVsQual);

  fHStatus = new TH1I("status", "Status", 8, 1, 9);
  fHStatus->SetFillColor(kBlue+1);
  fHStatus->SetFillStyle(3001);
  fHStatus->SetStats(0);
  fHStatus->SetDirectory(0);
  TAxis* xAxis = fHStatus->GetXaxis();
  xAxis->SetBinLabel(1, "OK");
  xAxis->SetBinLabel(2, "No event");
  xAxis->SetBinLabel(3, "No triggers");
  xAxis->SetBinLabel(4, "No SPD");
  xAxis->SetBinLabel(5, "No FMD");
  xAxis->SetBinLabel(6, "No vertex");
  xAxis->SetBinLabel(7, "Bad vertex");
  xAxis->SetBinLabel(8, "No centrality");
  fList->Add(fHStatus);

  fHVtxStatus = new TH1I("vtxStatus","Vertex Status",
			 kNotVtxZ,kVtxOK,kNotVtxZ+1);
  fHVtxStatus->SetFillColor(kGreen+1);
  fHVtxStatus->SetFillStyle(3001);
  fHVtxStatus->SetStats(0);
  fHVtxStatus->SetDirectory(0);
  xAxis = fHVtxStatus->GetXaxis();
  xAxis->SetBinLabel(kVtxOK,      "OK");
  xAxis->SetBinLabel(kNoVtx,      "None/bad status");
  xAxis->SetBinLabel(kNoSPDVtx,   "No SPD/bad status");
  xAxis->SetBinLabel(kFewContrib, "N_{contrib} <= 0");
  xAxis->SetBinLabel(kUncertain,  Form("#delta z > %4.2f", fMaxVzErr));
  xAxis->SetBinLabel(kNotVtxZ,    "Not Z vertexer");
  fList->Add(fHVtxStatus);

  fHTrgStatus = new TH1I("trgStatus", "Trigger Status", 
			 kOther, kNoTrgWords, kOther+1);
  fHTrgStatus->SetFillColor(kMagenta+1);
  fHTrgStatus->SetFillStyle(3001);
  fHTrgStatus->SetStats(0);
  fHTrgStatus->SetDirectory(0);
  xAxis = fHTrgStatus->GetXaxis();
  xAxis->SetBinLabel(kNoTrgWords,	"No words");
  xAxis->SetBinLabel(kPP2760Fast,	"FAST in pp@#sqrt{s}=2.76TeV"); 
  xAxis->SetBinLabel(kMUON,		"Muon trigger");
  xAxis->SetBinLabel(kTriggered,	"Triggered");
  xAxis->SetBinLabel(kMinBias,		"CINT1 (V0A||V0C||FASTOR)");
  xAxis->SetBinLabel(kMinBiasNoSPD,	"CINT5 (V0A||V0C)");
  xAxis->SetBinLabel(kV0AndTrg,		"CINT7 (V0A&&V0C)");
  xAxis->SetBinLabel(kMinBiasAD,	"CINT10 (ADA||ADC||V0A||V0C||FASTOR)");
  xAxis->SetBinLabel(kHighMult,		"N>>0");
  xAxis->SetBinLabel(kCentral,	        "Central"); 
  xAxis->SetBinLabel(kSemiCentral,	"Semi-central"); 
  xAxis->SetBinLabel(kDiffractive,	"Diffractive");
  xAxis->SetBinLabel(kADOr,	        "(ADA || ADC)");
  xAxis->SetBinLabel(kADAnd,	        "(ADA && ADC)");
  xAxis->SetBinLabel(kUser,	        "User");
  xAxis->SetBinLabel(kOther,	        "Other");
  fList->Add(fHTrgStatus);

  fHPileup = new TH1I("pileupStatus", "Pile-up status", 4, 0, 4);
  fHPileup->SetFillColor(kYellow+2);
  fHPileup->SetFillStyle(3001);
  fHPileup->SetStats(0);
  fHPileup->SetDirectory(0);
  xAxis = fHPileup->GetXaxis();
  xAxis->SetBinLabel(1,	"SPD tracklets");
  xAxis->SetBinLabel(2,	"Tracks"); 
  xAxis->SetBinLabel(3,	"Out-of-bunch");
  xAxis->SetBinLabel(4,	"SPD tracklets (bins)");
  fList->Add(fHPileup);

  if (AllowDisplaced()) fDisplacedVertex.SetupForData(fList, "", false);
}

//____________________________________________________________________
void
AliFMDEventInspector::StoreInformation()
{
  // Write TNamed objects to output list containing information about
  // the running conditions 
  DGUARD(fDebug,2,"Store information from AliFMDEventInspector");
  if (!fList) return;

  
  fList->Add(AliForwardUtil::MakeParameter("sys", fCollisionSystem));
  fList->Add(AliForwardUtil::MakeParameter("sNN", fEnergy));
  fList->Add(AliForwardUtil::MakeParameter("field", fField));
  fList->Add(AliForwardUtil::MakeParameter("runNo", fRunNumber));
  fList->Add(AliForwardUtil::MakeParameter("lowFlux", fLowFluxCut));
  fList->Add(AliForwardUtil::MakeParameter("ipMethod", fVtxMethod));
  //fList->Add(AliForwardUtil::MakeParameter("fpVtx",fUseFirstPhysicsVertex));
  fList->Add(AliForwardUtil::MakeParameter("v0and",fUseV0AND));
  fList->Add(AliForwardUtil::MakeParameter("pileup",   fPileupFlags));
  fList->Add(AliForwardUtil::MakeParameter("nPileUp", fMinPileupContrib));
  fList->Add(AliForwardUtil::MakeParameter("dPileup", fMinPileupDistance));
  fList->Add(AliForwardUtil::MakeParameter("satellite", AllowDisplaced()));
  fList->Add(AliForwardUtil::MakeParameter("alirootRev", 
					   AliForwardUtil::AliROOTRevision()));
  fList->Add(AliForwardUtil::MakeParameter("alirootBranch", 
					   AliForwardUtil::AliROOTBranch()));
  fList->Add(AliForwardUtil::MakeParameter("mc", fMC));

}

//____________________________________________________________________
void
AliFMDEventInspector::StoreProduction()
{
  if (!fList) return;

  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  if (!inputHandler) {
    AliWarning("Got no input handler");
    return;
  }
  TList *uiList = inputHandler->GetUserInfo();
  if (!uiList) {
    AliWarning("Got no user list from input tree");
    return;
  }

  AliProdInfo p(uiList);
  p.List();
  if (p.GetAlirootSvnVersion() <= 0) return;

  // Make our output list 
  TList* out = new TList;
  out->SetOwner(true);
  out->SetName("production");
  // out->SetTitle("ESD production details");
  fList->Add(out);

  TString period            = p.GetLHCPeriod();
  // TString aliROOTVersion    = p.GetAlirootVersion();
  fProdSVN                  = p.GetAlirootSvnVersion();
  // TString rootVersion       = p.GetRootVersion();
  Int_t   rootSVN           = p.GetRootSvnVersion();
  fProdPass                 = p.GetRecoPass();
  fProdMC                   = p.IsMC();
  
  if (period.Length() > 0) {
    TObjArray* pp = TPRegexp("LHC([0-9]+)([a-z]+)").MatchS(period);
    fProdYear     = static_cast<TObjString*>(pp->At(1))->String().Atoi();
    fProdLetter   = static_cast<TObjString*>(pp->At(2))->String()[0];
    pp->Delete();
  }
  
  out->Add(AliForwardUtil::MakeParameter("year", fProdYear));
  out->Add(AliForwardUtil::MakeParameter("letter", Int_t(fProdLetter)));
  out->Add(AliForwardUtil::MakeParameter("alirootSVN", fProdSVN));
  out->Add(AliForwardUtil::MakeParameter("rootSVN", rootSVN));
  out->Add(AliForwardUtil::MakeParameter("pass", fProdPass));
  out->Add(AliForwardUtil::MakeParameter("mc", fProdMC));
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

  // --- Process satellite event information is requested ------------
  if (AllowDisplaced()) { 
    if (!fDisplacedVertex.Process(event)) 
      AliWarning("Failed to process satellite event");
  }

  // --- Read trigger information from the ESD and store in AOD object
  if (!ReadTriggers(*event, triggers, nClusters)) { 
    if (fDebug > 2) {
      AliWarning("Failed to read triggers from ESD"); }
    fHStatus->Fill(3);
    return kNoTriggers;
  }
  Bool_t trig = (triggers & kTriggerMask);
  if (trig) fHTriggers->Fill(AliAODForwardMult::kWithTrigger);

  // --- Check if this is a high-flux event --------------------------
  const AliMultiplicity* testmult = event->GetMultiplicity();
  if (!testmult) {
    if (fDebug > 3) {
      AliWarning("No central multiplicity object found"); }
  }
  else 
    lowFlux = testmult->GetNumberOfTracklets() < fLowFluxCut;

  fHType->Fill(lowFlux ? 0 : 1);
 
  // --- Get the interaction point -----------------------------------
  Bool_t vzOk = ReadVertex(*event, ip);
  fHEventsTr->Fill(ip.Z());
  if (!vzOk) { 
    if (fDebug > 3) {
      AliWarning("Failed to read vertex from ESD"); }
    fHStatus->Fill(6);
    return kNoVertex;
  }
  if (trig) fHTriggers->Fill(AliAODForwardMult::kWithVertex);
  
   // --- Read centrality information 
  cent          = -10;
  UShort_t qual = 0;
  if (!ReadCentrality(*event, cent, qual)) {
    if (fDebug > 3) 
      AliWarning("Failed to get centrality");
  }
  // --- check centrality cut 
  if (fMinCent > -0.0001 && cent < fMinCent) return  kNoEvent; 
  if (fMaxCent > -0.0001 && cent > fMaxCent) return  kNoEvent; 

  if (cent >= 0) {
    fHCent->Fill(cent);
    fHCentVsQual->Fill((qual == 0 ? 197 : qual), cent);
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
  if (trig) fHTriggers->Fill(AliAODForwardMult::kAccepted);
  
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

  Int_t iqual;
  cent = AliForwardUtil::GetCentrality(esd, "V0M", iqual, (fDebug>1));
  qual = iqual;
  // We overwrite with satellite events, so we can be sure to get the
  // centrality determination from the satellite vertex selection
  if (AllowDisplaced() && fDisplacedVertex.IsSatellite()) {
    cent = fDisplacedVertex.GetCentralityPercentile();
    qual = 0;
  }

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
  AliPhysicsSelection* ps = 
    dynamic_cast<AliPhysicsSelection*>(ih->GetEventSelection());
  UInt_t trgMask = 0;
  AliTriggerAnalysis ta; 
  if (ps) {
    trgMask = ih->IsEventSelected();
  }
  else {
    // Replay HW triggers 
    Int_t v0A = ta.V0Trigger(&esd, AliTriggerAnalysis::kASide, true, false);
    Int_t v0C = ta.V0Trigger(&esd, AliTriggerAnalysis::kCSide, true, false);
    Int_t fso = ta.SPDFiredChips(&esd, 1, false, 0);
    Bool_t v0Abg = (v0A == 2);
    Bool_t v0Cbg = (v0C == 2);
    Bool_t v0Abb = (v0A == 1 || v0A == 3); // Include fake 
    Bool_t v0Cbb = (v0A == 1 || v0A == 3); // Include fake
    Bool_t isBg  = (v0Abg  || v0Cbg);
    if (!isBg && (v0Abb || v0Cbb || fso > 0)) trgMask |= AliVEvent::kMB;
    if (!isBg && (v0Abb || v0Cbb))            trgMask |= AliVEvent::kCINT5; 
    if (!isBg && (v0Abb && v0Cbb))            trgMask |= AliVEvent::kINT7;
    DMSG(fDebug, 1, "HW replay: v0A=%d v0C=%d gfo=%d mask=0x%x",
	 v0A, v0C, fso, trgMask);
  }
    
  // Check if this is a collision candidate (MB)
  ///
  // Historic remark: Note, that we should use the value cached in the
  //   input handler rather than calling IsCollisionCandiate directly
  //   on the AliPhysicsSelection obejct.  If we called the latter
  //   then the AliPhysicsSelection object would overcount by a factor
  //   of 2! :-(
    // UInt_t  trgMask  = ih->IsEventSelected();
  Bool_t  offline  = trgMask;
  Bool_t  fastonly = (trgMask & AliVEvent::kFastOnly);
  TString trigStr  = esd.GetFiredTriggerClasses();

  if (trigStr.IsNull()) fHTrgStatus->Fill(kNoTrgWords);
  
  if(AllowDisplaced()) {
    DMSG(fDebug,3,"Using displaced vertex stuff");
    // if (TMath::Abs(fDisplacedVertex.GetVertexZ()) >= 999) offline = false;
    if (fDisplacedVertex.IsSatellite()) 
      triggers |= AliAODForwardMult::kSatellite;
  }
  
  if (CheckFastPartition(fastonly)) {
    fHTrgStatus->Fill(kPP2760Fast);
    offline = false;
  }
  
  if (offline && CheckCosmics(trigStr)) {
    fHTrgStatus->Fill(kMUON);
    offline = false;
  }
  if (offline) fHTrgStatus->Fill(kTriggered);
  Int_t f = 0;
  if (trgMask & AliVEvent::kMB)          f += fHTrgStatus->Fill(kMinBias);
  if (trgMask & AliVEvent::kCINT5)       f += fHTrgStatus->Fill(kMinBiasNoSPD);
  if (trgMask & AliVEvent::kINT7)        f += fHTrgStatus->Fill(kV0AndTrg);
  if (trgMask & AliVEvent::kHighMult)    f += fHTrgStatus->Fill(kHighMult);
  if (trgMask & AliVEvent::kCentral)     f += fHTrgStatus->Fill(kCentral);
  if (trgMask & AliVEvent::kSemiCentral) f += fHTrgStatus->Fill(kSemiCentral);
  if (trgMask & AliVEvent::kDG5)         f += fHTrgStatus->Fill(kDiffractive);
  if (trgMask & AliVEvent::kUserDefined) f += fHTrgStatus->Fill(kUser);
  if (f <= 0) fHTrgStatus->Fill(kOther);
  
  // if (!CheckpAExtraV0(esd))             offline = false;

  DMSG(fDebug,2,"Event is %striggered by off-line", offline ? "" : "NOT ");

  if (offline)  triggers |= AliAODForwardMult::kOffline;
  // Only flag as in-elastic if a min-bias trigger is here
  if (trgMask & AliVEvent::kMB) {
    triggers |= AliAODForwardMult::kInel;
    CheckINELGT0(esd, nClusters, triggers);
  }
  // Flag as V0-AND if kINT7 is there 
  if (trgMask & AliVEvent::kINT7) {
    triggers |= AliAODForwardMult::kV0AND;
    if (fUseV0AND)
      // If we want to use V0-AND for NSD, flag as NSD 
      triggers |= AliAODForwardMult::kNSD;    
  }
  // Check if NSD trigger fired off-line 
  if (ta.IsOfflineTriggerFired(&esd, AliTriggerAnalysis::kNSD1)) 
    triggers |= AliAODForwardMult::kNSD;
  
  //   AliTriggerAnalysis ta; 
  if (ta.IsSPDClusterVsTrackletBG(&esd, false)) 
    // SPD outlier 
    triggers |= AliAODForwardMult::kSPDOutlier;

  CheckPileup(esd, triggers);
  CheckEmpty(trigStr, triggers);
  // if (CheckPileup(esd, triggers)) fHTriggers->Fill(kPileUp+.5);
  // if (CheckEmpty(trigStr, triggers)) fHTriggers->Fill(kEmpty+.5);

  CheckWords(esd, triggers); 

  AliAODForwardMult::FillTriggerHistogram(kTriggerMask, triggers, fHTriggers);
#if 0
#define TEST_TRIG_BIN(RET,BIN,TRIGGERS) \
  do { switch (BIN) { \
    case kInel:      RET = triggers & AliAODForwardMult::kInel;       break; \
    case kInelGt0:   RET = triggers & AliAODForwardMult::kInelGt0;    break; \
    case kNSD:       RET = triggers & AliAODForwardMult::kNSD;        break; \
    case kV0AND:     RET = triggers & AliAODForwardMult::kV0AND;      break; \
    case kEmpty:     RET = triggers & AliAODForwardMult::kEmpty;      break; \
    case kA:         RET = triggers & AliAODForwardMult::kA;          break; \
    case kB:         RET = triggers & AliAODForwardMult::kB;          break; \
    case kC:         RET = triggers & AliAODForwardMult::kC;          break; \
    case kE:         RET = triggers & AliAODForwardMult::kE;          break; \
    case kPileUp:    RET = triggers & AliAODForwardMult::kPileUp;     break; \
    case kMCNSD:     RET = triggers & AliAODForwardMult::kMCNSD;      break; \
    case kSatellite: RET = triggers & AliAODForwardMult::kSatellite;  break; \
    case kSpdOutlier:RET = triggers & AliAODForwardMult::kSPDOutlier; break; \
    case kOffline:   RET = triggers & AliAODForwardMult::kOffline;    break; \
    default:         RET = false; } } while(false)
      
  if (!fHTriggers) { 
    AliWarning("Histogram of triggers not defined - has init been called");
    return false;
  }
  
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
#endif
  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::CheckFastPartition(bool fastonly) const
{
  // For the 2.76 TeV p+p run (LHC11a), the FMD ran in the slow partition 
  // so it received no triggers from the fast partition. Therefore
  // the fast triggers are removed here but not for MC where all 
  // triggers are fast.
  if (TMath::Abs(fEnergy - 2750.) > 20) return false;
  if (fCollisionSystem != AliForwardUtil::kPP) return false;
  if (fMC) return false; // All MC events for pp @ 2.76TeV are `fast'
  // if (fastonly)
  //   DMSG(fDebug,1,"Fast trigger in pp @ sqrt(s)=2.76TeV removed");

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
  Int_t nTracklets = 0;
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
    if(TMath::Abs(spdmult->GetEta(j)) < 1) nTracklets++;
  }
  // If we at this point had non-zero nClusters, it's INEL>0
  if (nTracklets > 0) triggers |= AliAODForwardMult::kInelGt0;
  nClusters = nTracklets;

  // Loop over single clusters 
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
AliFMDEventInspector::CheckPileup(const AliESDEvent& esd, 
				  UInt_t& triggers) const
{
  // Check for multiple vertices (pile-up) with at least 3
  // contributors and at least 0.8cm from the primary vertex
  // if(fCollisionSystem != AliForwardUtil::kPP) return false;

  UInt_t locTrig = 0;
  // Check for standard SPD pile-up  
  if ((fPileupFlags & kSPD) && 
      esd.IsPileupFromSPD(fMinPileupContrib,fMinPileupDistance)) {
    fHPileup->Fill(0.5, 1);
    locTrig |= AliAODForwardMult::kPileupSPD;
  }

  // Check for multi-vertex pileup 
  if ((fPileupFlags & kTracks) && CheckMultiVertex(esd)) {
    fHPileup->Fill(1.5, 1);
    locTrig |= AliAODForwardMult::kPileupTrack;
  }

  // Check for out-of-bunch pileup
  if ((fPileupFlags & kOutOfBunch) && 
      (esd.GetHeader()->GetIRInt2ClosestInteractionMap() != 0 ||
       esd.GetHeader()->GetIRInt1ClosestInteractionMap() != 0)) {
    fHPileup->Fill(2.5, 1);
    locTrig |= AliAODForwardMult::kPileupBC;
  }
  // check SPD by multiplicity bin
  if ((fPileupFlags & kSPDBins) &&
      (esd.IsPileupFromSPDInMultBins())) {
    // Int_t nTracklets=GetMultiplicity()->GetNumberOfTracklets();
    // if      (nTracklets<20) return IsPileupFromSPD(3,0.8);
    // else if (nTracklets<50) return IsPileupFromSPD(4,0.8);
    // else                    return IsPileupFromSPD(5,0.8);
    
    locTrig |= AliAODForwardMult::kPileupBins;
    fHPileup->Fill(3.5, 1);
  }
    
  
  // Our result
  triggers      |= locTrig;
  Bool_t pileup =  locTrig != 0;
  if (pileup)    triggers |= AliAODForwardMult::kPileUp;
  return pileup;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::CheckMultiVertex(const AliESDEvent& esd,
				       Bool_t             checkOtherBC) const
{
  // Adapted from AliAnalysisUtils 
  // 
  // Parameters 
  const Double_t maxChi2nu    = 5;
  
  // Number of vertices 
  Int_t n = esd.GetNumberOfPileupVerticesTracks();
  if (n < 1) 
    // No additional vertices from tracks -> no pileup 
    return false;
  
  // Primary vertex 
  const AliESDVertex* primary = esd.GetPrimaryVertexTracks();
  if (primary->GetStatus() != 1) 
    // No good primary vertex, but many candidates -> pileup
    return true;
  
  // Bunch crossing number 
  Int_t bcPrimary = primary->GetBC();
  
  for (Int_t i = 0; i < n; i++) { 
    const AliESDVertex* other = esd.GetPileupVertexTracks(i);
    
    if (other->GetNContributors() < fMinPileupContrib) 
      // Not enough contributors to this vertex 
      continue;
    
    if (other->GetChi2perNDF() > maxChi2nu) 
      // Poorly determined vertex 
      continue;
    if (checkOtherBC) {
      Int_t bcOther = other->GetBC();
      if (bcOther != AliVTrack::kTOFBCNA && TMath::Abs(bcOther-bcPrimary) > 2) 
	// Pile-up from other BC 
	return true;      
    }
    // Calculate the distance 
    Double_t d  = -1;
    Double_t dx = primary->GetX() - other->GetX();
    Double_t dy = primary->GetY() - other->GetY();
    Double_t dz = primary->GetZ() - other->GetZ();
    Double_t covPrimary[6], covOther[6];
    primary->GetCovarianceMatrix(covPrimary); 
    other->GetCovarianceMatrix(covOther);
    TMatrixDSym v(3);
    v(0,0) = covPrimary[0] + covOther[0]; // diagonal
    v(1,1) = covPrimary[2] + covOther[2]; // diagonal
    v(2,2) = covPrimary[5] + covOther[5]; // diagonal
    v(1,0) = v(0,1) = covPrimary[1]+covOther[1]; // Band
    v(0,2) = v(1,2) = v(2,0) = v(2,1) = 0; // Off-diagonal+band
    // Try to invert 
    v.InvertFast();
    if (!v.IsValid()) 
      // Question if kStatus bit is every set after InvertFast?
      continue;
    d = (v(0,0) * dx * dx + v(1,1) * dy * dy + v(2,2) * dz * dz + 
	 2 * (v(0,1) * dx * dy + v(0,2) * dx * dz + v(1,2) * dy * dz));
    
    if (d < 0 || TMath::Sqrt(d) < fMinPileupDistance) 
      // Bad distance, or not fare enough from each other 
      continue;
    
    // Well separated vertices -> pileup
    return true;
  }
  return false;
}

//____________________________________________________________________
Bool_t
AliFMDEventInspector::CheckEmpty(const TString& trigStr, UInt_t& triggers) const
{
  if (trigStr.Contains("CBEAMB-ABCE-NOPF-ALL")|| 
      trigStr.Contains("CBEAMB-B-NOPF-ALLNOTRD")) {
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
    DMSG(fDebug,10,"Checking if %s (Coll) trigger %s is fired (%d)", 
	 word->GetTitle(), word->GetName(), word->GetUniqueID());
    if (!esd.IsTriggerClassFired(word->GetName())) continue;
    // if (fHWords) fHWords->Fill(word->GetName(), 1.);
    if (fHWords) fHWords->Fill(word->GetUniqueID());

    TString beamSide = word->GetTitle();
    DMSG(fDebug,10,"Found it - this is a %s trigger", beamSide.Data());

    if (!beamSide.EqualTo("B")) continue;
    triggers |= AliAODForwardMult::kB;

    // break; // No more to do here
  }
  TIter nextBg(&fBgWords);
  UInt_t all = (AliAODForwardMult::kA | 
		AliAODForwardMult::kC | 
		AliAODForwardMult::kE);
  while ((word = nextBg())) {
    DMSG(fDebug,10,"Checking if %s (BG) trigger %s is fired (%d)", 
	 word->GetTitle(), word->GetName(), word->GetUniqueID());
    if (!esd.IsTriggerClassFired(word->GetName())) continue;
    // if (fHWords) fHWords->Fill(word->GetName(), 1.);
    if (fHWords) fHWords->Fill(word->GetUniqueID());
    
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
  
  EVtxStatus s = kNoVtx;
  if (AllowDisplaced() && fDisplacedVertex.IsSatellite()) {
    s = kVtxOK;
    ip.SetZ(fDisplacedVertex.GetVertexZ());
  }
  else if (fVtxMethod == kPWGUD) 
    s = CheckPWGUDVertex(esd, ip);
  else if (fVtxMethod == kpA2012) 
    s = CheckpA2012Vertex(esd,ip);	
  else if (fVtxMethod == kpA2013) 
    s = CheckpA2013Vertex(esd,ip);	
  else 
    s = CheckVertex(esd, ip);
  
  fHVtxStatus->Fill(s);

  return s == kVtxOK;
}

//____________________________________________________________________
AliFMDEventInspector::EVtxStatus
AliFMDEventInspector::CheckPWGUDVertex(const AliESDEvent& esd, 
				       TVector3& ip)  const
{
  // This is the code used by the 1st physics people 
  const AliESDVertex* vertex    = esd.GetPrimaryVertex();
  if (!vertex  || !vertex->GetStatus()) {
    DMSG(fDebug,2,"No primary vertex (%p) or bad status %d", 
	 vertex, (vertex ? vertex->GetStatus() : -1));
    return kNoVtx;
  }
  const AliESDVertex* vertexSPD = esd.GetPrimaryVertexSPD();
  if (!vertexSPD || !vertexSPD->GetStatus()) {
    DMSG(fDebug,2,"No primary SPD vertex (%p) or bad status %d", 
	 vertexSPD, (vertexSPD ? vertexSPD->GetStatus() : -1));
    return kNoSPDVtx;
  }
    
  // if vertex is from SPD vertexZ, require more stringent cuts 
  if (vertex->IsFromVertexerZ()) { 
    if (vertex->GetDispersion() > fMaxVzErr || 
	vertex->GetZRes() > 1.25 * fMaxVzErr) {
      DMSG(fDebug,2,"Dispersion %f > %f or resolution %f > %f",
	   vertex->GetDispersion(), fMaxVzErr,
	   vertex->GetZRes(), 1.25 * fMaxVzErr);
      return kUncertain;
    }
  }
  ip.SetZ(vertex->GetZ());
  
  if(!vertex->IsFromVertexerZ()) {
    ip.SetX(vertex->GetX());
    ip.SetY(vertex->GetY());
  }
  return kVtxOK;
}
//____________________________________________________________________
AliFMDEventInspector::EVtxStatus
AliFMDEventInspector::CheckpA2012Vertex(const AliESDEvent& esd, 
					TVector3& ip)  const
{      
  const Int_t nMinContrib = 0;
  const AliESDVertex *vertex = esd.GetPrimaryVertexSPD();
  if (!vertex) return kNoSPDVtx;
  if (vertex->GetNContributors() <= nMinContrib) return kFewContrib;
  
  if (vertex->IsFromVertexerZ()) return kNotVtxZ;

  if (vertex->GetDispersion() >= 0.04 || vertex->GetZRes()>=0.25) 
    return kUncertain;

  ip.SetX(vertex->GetX());
  ip.SetY(vertex->GetY());
  ip.SetZ(vertex->GetZ());		

  return kVtxOK;
}
//____________________________________________________________________
AliFMDEventInspector::EVtxStatus
AliFMDEventInspector::CheckpA2013Vertex(const AliESDEvent& esd, 
					TVector3& ip)  const
{      
  // This code is adopted from 
  // 
  //   AliAnalysisUtils::IsVertexSelected2013pA(AliVEvent *event)
  // 
  const Int_t nMinContrib = 0;
  const AliESDVertex* primVtx = esd.GetPrimaryVertex();
  if (!primVtx) return kNoVtx;
  if (primVtx->GetNContributors() <= nMinContrib) return kFewContrib;
  
  const AliESDVertex* spdVtx = esd.GetPrimaryVertexSPD();
  if (!spdVtx) return kNoSPDVtx;
  if (spdVtx->GetNContributors() <= nMinContrib) return kFewContrib;
  
  // TString vtxTyp = spdVtx->GetTitle();
  // if (spdVtx->IsFromVertexerZ()) {
  if (spdVtx->IsFromVertexerZ() && spdVtx->GetZRes()>0.25) return kUncertain;

  Bool_t correlateVtx = true;
  if (correlateVtx) { 
    if (TMath::Abs(spdVtx->GetZ() - primVtx->GetZ()) > 0.5) return kUncertain;
  }

  ip.SetX(primVtx->GetX());
  ip.SetY(primVtx->GetY());
  ip.SetZ(primVtx->GetZ());		

  return kVtxOK;
}

//____________________________________________________________________
AliFMDEventInspector::EVtxStatus
AliFMDEventInspector::CheckVertex(const AliESDEvent& esd, 
				  TVector3& ip) const
{
  // Use standard SPD vertex (perhaps preferable for Pb+Pb)
  // Get the vertex 
  const AliESDVertex* vertex = esd.GetPrimaryVertexSPD();
  if (!vertex) { 
    if (fDebug > 2) {
      AliWarning("No SPD vertex found in ESD"); }
    return kNoSPDVtx;
  }
    
  // #if 0 // Check disabled - seem to kill a lot of PbPb events
  // Check that enough tracklets contributed 
  if(vertex->GetNContributors() <= 0) {
    DMSG(fDebug,2,"Number of contributors to vertex is %d<=0",
	 vertex->GetNContributors());
    ip.SetZ(0);
    return kFewContrib;
  } 
  // #endif

  // Check that the uncertainty isn't too large 
  if (vertex->GetZRes() > fMaxVzErr) { 
    DMSG(fDebug,2,"Uncertaintity in Z of vertex is too large %f > %f", 
	 vertex->GetZRes(), fMaxVzErr);
    return kUncertain;
  }
    
  // Get the z coordiante 
  ip.SetZ(vertex->GetZ());
  const AliESDVertex* vertexXY = esd.GetPrimaryVertex();
    
  DMSG(fDebug, 2, "Read vertex: Primary (%f,%f,%f) [%3s] SPD (%f,%f,%f) [%3s]",
       vertexXY->GetX(), vertexXY->GetY(), vertexXY->GetZ(),
       vertexXY->IsFromVertexerZ() ? "Z" : "XYZ",
       vertex  ->GetX(), vertex  ->GetY(), vertex  ->GetZ(),
       vertex  ->IsFromVertexerZ() ? "Z" : "XYZ");

  // If from vertexer Z, we get the mean of X and Y
  ip.SetX(vertexXY->GetX());
  ip.SetY(vertexXY->GetY());
  // if(!vertexXY->IsFromVertexerZ()) {
  //   ip.SetX(vertexXY->GetX());
  //   ip.SetY(vertexXY->GetY());
  // }
  return kVtxOK;
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
  fRunNumber       = esd->GetRunNumber();
  
  Int_t a0 = esd->GetBeamParticleA(0);
  Int_t a1 = esd->GetBeamParticleA(1);
  Int_t z0 = esd->GetBeamParticleZ(0);
  Int_t z1 = esd->GetBeamParticleZ(1);
  if (fCentMethod.EqualTo("default", TString::kIgnoreCase)) {
    if (fCollisionSystem == 1 || fCollisionSystem == 2) 
      SetCentralityMethod("V0M");
    if (fCollisionSystem == 3) { 
      if (a0 != 0 && a1 != 0) 
	SetCentralityMethod(a0 > a1 ? "ZNC" : "ZNA");
      else 
	SetCentralityMethod("V0M");    
    }
  }
  // store the beam species 
  fList->Add(AliForwardUtil::MakeParameter("beam0a", a0));
  fList->Add(AliForwardUtil::MakeParameter("beam0z", z0));
  fList->Add(AliForwardUtil::MakeParameter("beam1a", a1));
  fList->Add(AliForwardUtil::MakeParameter("beam1z", z1));

  StoreProduction();
  StoreInformation();
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
AliFMDEventInspector::Print(Option_t*) const
{
  // 
  // Print information
  // 
  //   option Not used 
  //
  AliForwardUtil::PrintTask(*this);
  TString sNN(AliForwardUtil::CenterOfMassEnergyString(fEnergy));
  sNN.Strip(TString::kBoth, '0');
  sNN.ReplaceAll("GeV", " GeV");
  TString field(AliForwardUtil::MagneticFieldString(fField));
  field.ReplaceAll("p",  "+");
  field.ReplaceAll("m",  "-");
  field.ReplaceAll("kG", " kG");

  gROOT->IncreaseDirLevel();
  PFV("Production", "");
  gROOT->IncreaseDirLevel();
  PF("Period", "LHC%02d%c", fProdYear, fProdLetter);
  PFB("MC anchor", fProdMC);
  PFV("AliROOT Revision", fProdSVN);
  gROOT->DecreaseDirLevel();
  PFV("Vertex bins", fVtxAxis.GetNbins());
  PF("Vertex Range", "[%+6.1f,%+6.1f]", fVtxAxis.GetXmin(), fVtxAxis.GetXmax());
  PFV("Low flux cut",		fLowFluxCut );
  PFV("Max(delta v_z)",		fMaxVzErr );
  PF("Pile-up methods",	        "spd:%s,tracks:%s,out-of-bunch:%s", 
     (fPileupFlags & kSPD        ? "yes" : "no"),
     (fPileupFlags & kTracks     ? "yes" : "no"),
     (fPileupFlags & kOutOfBunch ? "yes" : "no"));
  PFV("Min(nContrib_pileup)",	fMinPileupContrib );
  PFV("Min(v-pileup)",		fMinPileupDistance );
  PFV("System", AliForwardUtil::CollisionSystemString(fCollisionSystem));
  PFV("CMS energy per nucleon",	sNN);
  PFV("Field",			field);
  PFB("Satellite events",	AllowDisplaced());
  // PFB("Use 2012 pA vertex",	fUsepA2012Vertex );
  // PFB("Use PWG-UD vertex",	fUseFirstPhysicsVertex);
  TString vtxMethod("normal");
  switch(fVtxMethod) { 
  case kNormal:    vtxMethod = "normal" ;    break;
  case kpA2012:    vtxMethod = "pA 2012";    break;
  case kpA2013:    vtxMethod = "pA 2013";    break;
  case kPWGUD:     vtxMethod = "PWG-UD";     break;
  case kDisplaced: vtxMethod = "Satellite";  break;
  }
  PFV("IP method",              vtxMethod);
  PFB("Simulation input",	fMC );
  PFV("Centrality method",	fCentMethod);
  PFV("Centrality axis",        (!fCentAxis ? "none" : ""));
  if (!fCentAxis) {return; }
  Int_t nBin = fCentAxis->GetNbins();
  for (Int_t i = 0; i < nBin; i++) { 
    if (i != 0 && (i % 10) == 0) std::cout << '\n';
    if ((i % 10) == 0)           std::cout << "    ";
    std::cout << std::setw(5) << fCentAxis->GetBinLowEdge(i+1) << '-';
  }
  std::cout << std::setw(5) << fCentAxis->GetBinUpEdge(nBin) << std::endl;
  gROOT->DecreaseDirLevel();
}

  
//
// EOF
//

