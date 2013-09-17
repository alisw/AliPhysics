// 
// Calculate the multiplicity in the forward regions event-by-event 
// 
// Inputs: 
//   - AliESDEvent 
//   - Kinematics
//   - Track references
//
// Outputs: 
//   - AliAODForwardMult 
// 
// Histograms 
//   
// Corrections used 
// 
#include "AliForwardMCMultiplicityTask.h"
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliAODHandler.h"
#include "AliMultiplicity.h"
#include "AliInputEventHandler.h"
#include "AliForwardCorrectionManager.h"
#include "AliAnalysisManager.h"
#include <TH1.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TROOT.h>

//====================================================================
AliForwardMCMultiplicityTask::AliForwardMCMultiplicityTask()
  : AliForwardMultiplicityBase(),
    fESDFMD(),
    fMCESDFMD(),
    fMCHistos(),
    fMCAODFMD(),
    fMCRingSums(),
    fPrimary(0),
    fEventInspector(),
    fSharingFilter(),
    fDensityCalculator(),
    fCorrections(),
    fHistCollector(),
    fEventPlaneFinder()
{
  // 
  // Constructor
  //
}

//____________________________________________________________________
AliForwardMCMultiplicityTask::AliForwardMCMultiplicityTask(const char* name)
  : AliForwardMultiplicityBase(name), 
    fESDFMD(),
    fMCESDFMD(),
    fMCHistos(),
    fMCAODFMD(kTRUE),
    fMCRingSums(),
    fPrimary(0),
    fEventInspector("event"),
    fSharingFilter("sharing"), 
    fDensityCalculator("density"),
    fCorrections("corrections"),
    fHistCollector("collector"),
    fEventPlaneFinder("eventplane")
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
}

//____________________________________________________________________
AliForwardMCMultiplicityTask::AliForwardMCMultiplicityTask(const AliForwardMCMultiplicityTask& o)
  : AliForwardMultiplicityBase(o),
    fESDFMD(o.fESDFMD),
    fMCESDFMD(o.fMCESDFMD),
    fMCHistos(o.fMCHistos),
    fMCAODFMD(o.fMCAODFMD),
    fMCRingSums(o.fMCRingSums),
    fPrimary(o.fPrimary),
    fEventInspector(o.fEventInspector),
    fSharingFilter(o.fSharingFilter),
    fDensityCalculator(o.fDensityCalculator),
    fCorrections(o.fCorrections),
    fHistCollector(o.fHistCollector),
    fEventPlaneFinder(o.fEventPlaneFinder)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
}

//____________________________________________________________________
AliForwardMCMultiplicityTask&
AliForwardMCMultiplicityTask::operator=(const AliForwardMCMultiplicityTask& o)
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this object 
  //
  if (&o == this) return *this;
  AliForwardMultiplicityBase::operator=(o);

  fEventInspector    = o.fEventInspector;
  fSharingFilter     = o.fSharingFilter;
  fDensityCalculator = o.fDensityCalculator;
  fCorrections       = o.fCorrections;
  fHistCollector     = o.fHistCollector;
  fEventPlaneFinder  = o.fEventPlaneFinder;
  fMCHistos          = o.fMCHistos;
  fMCAODFMD          = o.fMCAODFMD;
  fMCRingSums        = o.fMCRingSums;
  fPrimary           = o.fPrimary;
  return *this;
}

//____________________________________________________________________
void
AliForwardMCMultiplicityTask::SetOnlyPrimary(Bool_t use)
{
  fSharingFilter.GetTrackDensity().SetUseOnlyPrimary(use);
  fCorrections.SetSecondaryForMC(!use);
}

//____________________________________________________________________
void
AliForwardMCMultiplicityTask::CreateBranches(AliAODHandler* ah)
{
  // 
  // Create output objects 
  // 
  //
  AliForwardMultiplicityBase::CreateBranches(ah);

  TObject* mcobj = &fMCAODFMD;
  ah->AddBranch("AliAODForwardMult", &mcobj);

  fPrimary = new TH2D("primary", "MC Primaries", 1,0,1,20,0,TMath::TwoPi());
  fPrimary->SetXTitle("#eta");
  fPrimary->SetYTitle("#varphi [radians]");
  fPrimary->SetZTitle("d^{2}N_{ch}/d#etad#phi");
  fPrimary->Sumw2();
  fPrimary->SetStats(0);
  fPrimary->SetDirectory(0);
    
  ah->AddBranch("TH2D", &fPrimary);
}


//____________________________________________________________________
void
AliForwardMCMultiplicityTask::InitMembers(const TAxis* pe, const TAxis* pv)
{
  // 
  // Initialise the sub objects and stuff.  Called on first event 
  // 
  //
  AliForwardMultiplicityBase::InitMembers(pe, pv);

  fMCHistos.Init(*pe);
  fMCAODFMD.Init(*pe);
  fMCRingSums.Init(*pe);

  AliForwardUtil::Histos::RebinEta(fPrimary, *pe);

  TList* mcRings = new TList;
  mcRings->SetName("mcRingSums");
  mcRings->SetOwner();
  fList->Add(mcRings);

  mcRings->Add(fMCRingSums.Get(1, 'I'));
  mcRings->Add(fMCRingSums.Get(2, 'I'));
  mcRings->Add(fMCRingSums.Get(2, 'O'));
  mcRings->Add(fMCRingSums.Get(3, 'I'));
  mcRings->Add(fMCRingSums.Get(3, 'O'));
  fMCRingSums.Get(1, 'I')->SetMarkerColor(AliForwardUtil::RingColor(1, 'I'));
  fMCRingSums.Get(2, 'I')->SetMarkerColor(AliForwardUtil::RingColor(2, 'I'));
  fMCRingSums.Get(2, 'O')->SetMarkerColor(AliForwardUtil::RingColor(2, 'O'));
  fMCRingSums.Get(3, 'I')->SetMarkerColor(AliForwardUtil::RingColor(3, 'I'));
  fMCRingSums.Get(3, 'O')->SetMarkerColor(AliForwardUtil::RingColor(3, 'O'));
}

//____________________________________________________________________
void
AliForwardMCMultiplicityTask::UserExec(Option_t*)
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  

  // Read production details 
  if (fFirstEvent) 
    fEventInspector.ReadProductionDetails(MCEvent());
    
  // Get the input data 
  AliESDEvent* esd     = GetESDEvent();
  AliMCEvent*  mcEvent = MCEvent();
  if (!esd || !mcEvent) return;

  // Clear stuff 
  fHistos.Clear();
  fESDFMD.Clear();
  fAODFMD.Clear();
  fAODEP.Clear();
  fMCHistos.Clear();
  fMCESDFMD.Clear();
  fMCAODFMD.Clear();
  fPrimary->Reset();

  Bool_t   lowFlux   = kFALSE;
  UInt_t   triggers  = 0;
  UShort_t ivz       = 0;
  TVector3 ip(1024, 1024, 0);
  Double_t cent      = -1;
  UShort_t nClusters = 0;
  UInt_t   found     = fEventInspector.Process(esd, triggers, lowFlux, 
					       ivz, ip, cent, nClusters);
  UShort_t ivzMC    = 0;
  Double_t vzMC     = 0;
  Double_t phiR     = 0;
  Double_t b        = 0;
  Double_t cMC      = 0;
  Int_t    npart    = 0;
  Int_t    nbin     = 0;
  // UInt_t   foundMC  = 
  fEventInspector.ProcessMC(mcEvent, triggers, ivzMC, vzMC, b, cMC,
			    npart, nbin, phiR);
  fEventInspector.CompareResults(ip.Z(), vzMC, cent, cMC, b, npart, nbin);
  
  //Store all events
  MarkEventForStore();
  
  Bool_t isAccepted = true;
  if (found & AliFMDEventInspector::kNoEvent)    isAccepted = false; // return;
  if (found & AliFMDEventInspector::kNoTriggers) isAccepted = false; // return;
  //MarkEventForStore();
  // Always set the B trigger - each MC event _is_ a collision 
  triggers |= AliAODForwardMult::kB;
  // Set trigger bits, and mark this event for storage 
  fAODFMD.SetTriggerBits(triggers);
  fAODFMD.SetSNN(fEventInspector.GetEnergy());
  fAODFMD.SetSystem(fEventInspector.GetCollisionSystem());
  fAODFMD.SetCentrality(cent);
  fAODFMD.SetNClusters(nClusters);

  fMCAODFMD.SetTriggerBits(triggers);
  fMCAODFMD.SetSNN(fEventInspector.GetEnergy());
  fMCAODFMD.SetSystem(fEventInspector.GetCollisionSystem());
  fMCAODFMD.SetCentrality(cent);
  fMCAODFMD.SetNClusters(nClusters);
  
  //All events should be stored - HHD
  //if (isAccepted) MarkEventForStore();

  // Disable this check on SPD - will bias data 
  // if (found & AliFMDEventInspector::kNoSPD)  isAccepted = false; // return;
  if (found & AliFMDEventInspector::kNoFMD)     isAccepted = false; // return;
  if (found & AliFMDEventInspector::kNoVertex)  isAccepted = false; // return;

  if (isAccepted) {
    fAODFMD.SetIpZ(ip.Z());
    fMCAODFMD.SetIpZ(ip.Z());
  }
  if (found & AliFMDEventInspector::kBadVertex) isAccepted = false; // return;

  // We we do not want to use low flux specific code, we disable it here. 
  if (!fEnableLowFlux) lowFlux = false;
  
  

  // Get FMD data 
  AliESDFMD*  esdFMD  = esd->GetFMDData();

  // Apply the sharing filter (or hit merging or clustering if you like)
  if (isAccepted && !fSharingFilter.Filter(*esdFMD, lowFlux, fESDFMD,ip.Z())) { 
    AliWarning("Sharing filter failed!");
    return;
  }
  if (!fSharingFilter.FilterMC(*esdFMD, *mcEvent, ip.Z(),fMCESDFMD,fPrimary)) { 
    AliWarning("MC Sharing filter failed!");
    return;
  }

  // Store some MC parameters in corners of histogram :-)
  fPrimary->SetBinContent(0,                      0,                    vzMC);
  fPrimary->SetBinContent(fPrimary->GetNbinsX()+1,0,                    phiR);
  fPrimary->SetBinContent(fPrimary->GetNbinsX()+1,fPrimary->GetNbinsY(),cMC);
  

  if (!isAccepted) return; // Exit on MC event w/o trigger, vertex, data
  // HHD if (!isAccepted) return; // Exit on MC event w/o trigger, vertex, data
  
  //MarkEventForStore();
  fSharingFilter.CompareResults(fESDFMD, fMCESDFMD);

  // Calculate the inclusive charged particle density 
  if (!fDensityCalculator.Calculate(fESDFMD, fHistos, lowFlux, cent, ip)) { 
    AliWarning("Density calculator failed!");
    return;
  }
  if (!fDensityCalculator.CalculateMC(fMCESDFMD, fMCHistos)) { 
    AliWarning("MC Density calculator failed!");
    return;
  }
  fDensityCalculator.CompareResults(fHistos, fMCHistos);
  
  if (fEventInspector.GetCollisionSystem() == AliFMDEventInspector::kPbPb) {
    if (!fEventPlaneFinder.FindEventplane(esd, fAODEP, 
					  &(fAODFMD.GetHistogram()) , &fHistos))
      AliWarning("Eventplane finder failed!");
  }

  // Do the secondary and other corrections. 
  if (!fCorrections.Correct(fHistos, ivz)) { 
    AliWarning("Corrections failed");
    return;
  }
  if (!fCorrections.CorrectMC(fMCHistos, ivz)) { 
    AliWarning("MC Corrections failed");
    return;
  }
  fCorrections.CompareResults(fHistos, fMCHistos);
    
  if (!fHistCollector.Collect(fHistos, fRingSums, 
			      ivz, fAODFMD.GetHistogram(),
			      fAODFMD.GetCentrality())) {
    AliWarning("Histogram collector failed");
    return;
  }
  if (!fHistCollector.Collect(fMCHistos, fMCRingSums, 
			      ivz, fMCAODFMD.GetHistogram())) {
    AliWarning("MC Histogram collector failed");
    return;
  }
  // Copy underflow bins to overflow bins - always full phi coverage 
  TH2&  hMC  = fMCAODFMD.GetHistogram();
  Int_t nEta = hMC.GetNbinsX();
  Int_t nY   = hMC.GetNbinsY();
  for (Int_t iEta = 1; iEta <= nEta; iEta++) {
    hMC.SetBinContent(iEta, nY+1, hMC.GetBinContent(iEta, 0));
  }

  if (fAODFMD.IsTriggerBits(AliAODForwardMult::kInel))
    fHData->Add(&(fAODFMD.GetHistogram()));

  PostData(1, fList);
}

//____________________________________________________________________
void
AliForwardMCMultiplicityTask::EstimatedNdeta(const TList* input, 
					     TList*       output) const
{
  AliForwardMultiplicityBase::EstimatedNdeta(input, output);
  MakeRingdNdeta(input, "mcRingSums", output, "mcRingResults", 24);
}



//
// EOF
//
