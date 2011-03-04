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
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliHeader.h"
#include "AliMCEventHandler.h"
//====================================================================
AliFMDEventInspector::AliFMDEventInspector()
  : TNamed(),
    fHEventsTr(0), 
    fHEventsTrVtx(0),
    fHTriggers(0),
    fHType(0),
    fHWords(0),
    fHCent(0),
    fLowFluxCut(1000),
    fMaxVzErr(0.2),
    fList(0),
    fEnergy(0),
    fField(999), 
    fCollisionSystem(kUnknown),
    fDebug(0)
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
    fHTriggers(0),
    fHType(0),
    fHWords(0),
    fHCent(0),
    fLowFluxCut(1000),
    fMaxVzErr(0.2),
    fList(0),
    fEnergy(0),
    fField(999), 
    fCollisionSystem(kUnknown),
    fDebug(0)
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
    fHTriggers(o.fHTriggers),
    fHType(o.fHType),
    fHWords(o.fHWords),
    fHCent(o.fHCent),
    fLowFluxCut(o.fLowFluxCut),
    fMaxVzErr(o.fMaxVzErr),
    fList(o.fList),
    fEnergy(o.fEnergy),
    fField(o.fField), 
    fCollisionSystem(o.fCollisionSystem),
    fDebug(0)
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
  if (fHEventsTr)    delete fHEventsTr;
  if (fHEventsTrVtx) delete fHEventsTrVtx;
  if (fHTriggers)    delete fHTriggers;  
  if (fHType)        delete fHType;
  if (fHWords)       delete fHWords;
  if (fHCent)        delete fHCent;
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
  fHTriggers         = o.fHTriggers;
  fHType             = o.fHType;
  fHWords            = o.fHWords;
  fHCent             = o.fHCent;
  fLowFluxCut        = o.fLowFluxCut;
  fMaxVzErr          = o.fMaxVzErr;
  fDebug             = o.fDebug;
  fList              = (o.fList ? new TList : 0);
  fEnergy            = o.fEnergy;
  fField             = o.fField;
  fCollisionSystem   = o.fCollisionSystem;
  if (fList) { 
    fList->SetName(GetName());
    if (fHEventsTr)    fList->Add(fHEventsTr);
    if (fHEventsTrVtx) fList->Add(fHEventsTrVtx);
    if (fHTriggers)    fList->Add(fHTriggers);
    if (fHType)        fList->Add(fHType);
    if (fHWords)       fList->Add(fHWords);
    if (fHCent)        fList->Add(fHCent);
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
  fHEventsTr = new TH1I("nEventsTr", "Number of events w/trigger", 
			vtxAxis.GetNbins(), 
			vtxAxis.GetXmin(), 
			vtxAxis.GetXmax());
  fHEventsTr->SetXTitle("v_{z} [cm]");
  fHEventsTr->SetYTitle("# of events");
  fHEventsTr->SetFillColor(kRed+1);
  fHEventsTr->SetFillStyle(3001);
  fHEventsTr->SetDirectory(0);
  // fHEventsTr->Sumw2();
  fList->Add(fHEventsTr);

  fHEventsTrVtx = new TH1I("nEventsTrVtx", 
			   "Number of events w/trigger and vertex", 
			   vtxAxis.GetNbins(), 
			   vtxAxis.GetXmin(), 
			   vtxAxis.GetXmax());
  fHEventsTrVtx->SetXTitle("v_{z} [cm]");
  fHEventsTrVtx->SetYTitle("# of events");
  fHEventsTrVtx->SetFillColor(kBlue+1);
  fHEventsTrVtx->SetFillStyle(3001);
  fHEventsTrVtx->SetDirectory(0);
  // fHEventsTrVtx->Sumw2();
  fList->Add(fHEventsTrVtx);

      
  fHTriggers = new TH1I("triggers", "Triggers", 10, 0, 10);
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
  fHTriggers->GetXaxis()->SetBinLabel(kMCNSD  +1,"nsd");
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

  fHCent = new TH1F("cent", "Centrality", 101, -1.5, 100.5);
  fHCent->SetFillColor(kBlue+1);
  fHCent->SetFillStyle(3001);
  fHCent->SetStats(0);
  fHCent->SetDirectory(0);
  fHCent->SetXTitle("Centrality [%]");
  fHCent->SetYTitle("Events");
  fList->Add(fHCent);
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
			      Double_t&          cent)
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

  // Check that we have an event 
  if (!event) { 
    AliWarning("No ESD event found for input event");
    return kNoEvent;
  }

  // Read trigger information from the ESD and store in AOD object
  if (!ReadTriggers(event, triggers)) { 
    if (fDebug > 2) {
      AliWarning("Failed to read triggers from ESD"); }
    return kNoTriggers;
  }

  // Check if this is a high-flux event 
  const AliMultiplicity* testmult = event->GetMultiplicity();
  if (!testmult) {
    if (fDebug > 3) {
      AliWarning("No central multiplicity object found"); }
    return kNoSPD;
  }
  lowFlux = testmult->GetNumberOfTracklets() < fLowFluxCut;

  fHType->Fill(lowFlux ? 0 : 1);
  
  cent = -10;
  AliCentrality* centObj = const_cast<AliESDEvent*>(event)->GetCentrality();
  if (centObj) {
    // AliInfo(Form("Got centrality object %p with quality %d", 
    //              centObj, centObj->GetQuality()));
    // centObj->Print();
    cent = centObj->GetCentralityPercentileUnchecked("V0M");  
  }
  // AliInfo(Form("Centrality is %f", cent));
  fHCent->Fill(cent);

  // Check the FMD ESD data 
  if (!event->GetFMDData()) { 
    if (fDebug > 3) {
      AliWarning("No FMD data found in ESD"); }
    return kNoFMD;
  }

  // Get the vertex information 
  vz          = 0;
  Bool_t vzOk = ReadVertex(event, vz);

  fHEventsTr->Fill(vz);
  if (!vzOk) { 
    if (fDebug > 3) {
      AliWarning("Failed to read vertex from ESD"); }
    return kNoVertex;
  }
  fHEventsTrVtx->Fill(vz);

  // Get the vertex bin 
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
AliFMDEventInspector::ReadTriggers(const AliESDEvent* esd, UInt_t& triggers)
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
  
  AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  AliMCEvent* mcEvent = 0;
  if(mcHandler)
    mcEvent = mcHandler->MCEvent();
  
  if(mcEvent) {
    //Assign MC only triggers : True NSD etc.
    AliHeader* header            = mcEvent->Header();
    AliGenEventHeader* genHeader = header->GenEventHeader();
    
    AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
    AliGenDPMjetEventHeader* dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*>(header->GenEventHeader());
    Bool_t sd = kFALSE;
    if(pythiaGenHeader) {
      Int_t pythiaType = pythiaGenHeader->ProcessType();
      if(pythiaType==92 || pythiaType==93)
	sd = kTRUE;
    }
    if(dpmHeader) {
      Int_t processType = dpmHeader->ProcessType();
      if(processType == 5 || processType == 6)  
	sd = kTRUE;
      
    }
    if(!sd) {
      triggers |= AliAODForwardMult::kMCNSD;
      fHTriggers->Fill(kMCNSD+0.5);
    }
    
  }
    
  // Check if this is a collision candidate (INEL)
  // Note, that we should use the value cached in the input 
  // handler rather than calling IsCollisionCandiate directly 
  // on the AliPhysicsSelection obejct.  If we called the latter
  // then the AliPhysicsSelection object would overcount by a 
  // factor of 2! :-(
  Bool_t inel = ih->IsEventSelected();
  if (inel) { 
    triggers |= AliAODForwardMult::kInel;
    fHTriggers->Fill(kInel+0.5);
  }

  // If this is inel, see if we have a tracklet 
  if (inel) { 
    const AliMultiplicity* spdmult = esd->GetMultiplicity();
    if (!spdmult) {
      AliWarning("No SPD multiplicity");
    }
    else { 
      // Check if we have one or more tracklets 
      // in the range -1 < eta < 1 to set the INEL>0 
      // trigger flag. 
      Int_t n = spdmult->GetNumberOfTracklets();
      for (Int_t j = 0; j < n; j++) { 
	if(TMath::Abs(spdmult->GetEta(j)) < 1) { 
	  triggers |= AliAODForwardMult::kInelGt0;
	  fHTriggers->Fill(kInelGt0+.5);
	  break;
	}
      }
    }
  }

  // Analyse some trigger stuff 
  AliTriggerAnalysis ta;
  if (ta.IsOfflineTriggerFired(esd, AliTriggerAnalysis::kNSD1)) {
    triggers |= AliAODForwardMult::kNSD;
    fHTriggers->Fill(kNSD+.5);
  }
  //Check pileup
  Bool_t pileup =  esd->IsPileupFromSPD(3,0.8);
  if (pileup) {
    triggers |= AliAODForwardMult::kPileUp;
    fHTriggers->Fill(kPileUp+.5);
  }
    
  // Get trigger stuff 
  TString trigStr = esd->GetFiredTriggerClasses();
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

  if (trigStr.Contains("CINT1A-ABCE-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kA;
    fHTriggers->Fill(kA+.5);
  }

  if (trigStr.Contains("CINT1B-ABCE-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kB;
    fHTriggers->Fill(kB+.5);
  }


  if (trigStr.Contains("CINT1C-ABCE-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kC;
    fHTriggers->Fill(kC+.5);
  }

  if (trigStr.Contains("CINT1-E-NOPF-ALL")) {
    triggers |= AliAODForwardMult::kE;
    fHTriggers->Fill(kE+.5);
  }

  return kTRUE;
}
//____________________________________________________________________
Bool_t
AliFMDEventInspector::ReadVertex(const AliESDEvent* esd, Double_t& vz)
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

  // HHD test
  /*
  if (vertex->IsFromVertexerZ()) { 
    return kFALSE;
  }
  if (TMath::Sqrt(TMath::Power(vertex->GetX(),2) + TMath::Power(vertex->GetY(),2)) > 3 ) { 
      return kFALSE;
  }
  */
  
  // Get the z coordiante 
  vz = vertex->GetZ();
  return kTRUE;
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
  
  std::cout << ind << "AliFMDEventInspector: " << GetName() << '\n'
	    << ind << " Low flux cut:           " << fLowFluxCut << '\n'
	    << ind << " Max(delta v_z):         " << fMaxVzErr << " cm\n"
	    << ind << " System:                 " 
	    << AliForwardUtil::CollisionSystemString(fCollisionSystem) << '\n'
	    << ind << " CMS energy per nucleon: " << sNN << '\n'
	    << ind << " Field:                  " <<  field << std::endl;
}

  
//
// EOF
//

