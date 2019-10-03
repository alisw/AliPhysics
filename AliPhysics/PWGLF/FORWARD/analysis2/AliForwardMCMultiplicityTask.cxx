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
#include "AliMCEvent.h"
#include "AliAODHandler.h"
#include "AliMultiplicity.h"
#include "AliInputEventHandler.h"
#include "AliForwardCorrectionManager.h"
#include "AliAnalysisManager.h"
#include <TH1.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TROOT.h>
#define MCAOD_SLOT 4
#define PRIMARY_SLOT 5
#ifdef POST_AOD
# define DEFINE(N,C) DefineOutput(N,C)
# define POST(N,O)   PostData(N,O)
#else
# define DEFINE(N,C) do { } while(false)
# define POST(N,O)   do { } while(false)
#endif
#ifndef ENABLE_TIMING
# define MAKE_SW(NAME) do {} while(false)
# define START_SW(NAME) do {} while(false)
# define FILL_SW(NAME,WHICH) do {} while(false)
#else
# define MAKE_SW(NAME) TStopwatch NAME
# define START_SW(NAME) if (fDoTiming) NAME.Start(true)
# define FILL_SW(NAME,WHICH)				\
  if (fDoTiming) fHTiming->Fill(WHICH,NAME.CpuTime())
#endif

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
    // fMultEventClassifier(),
    fESDFixer(),
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
    // fMultEventClassifier("multClass"),
    fESDFixer("esdFizer"),
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
  fPrimary = new TH2D("primary", "MC Primaries", 1,0,1,20,0,TMath::TwoPi());
  fPrimary->SetXTitle("#eta");
  fPrimary->SetYTitle("#varphi [radians]");
  fPrimary->SetZTitle("d^{2}N_{ch}/d#etad#phi");
  fPrimary->Sumw2();
  fPrimary->SetStats(0);
  fPrimary->SetDirectory(0);
  DEFINE(MCAOD_SLOT,AliAODForwardMult::Class());
  DEFINE(PRIM_SLOT, TH2D::Class());
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
  ah->AddBranch("TH2D", &fPrimary);
}

//____________________________________________________________________
Bool_t
AliForwardMCMultiplicityTask::Book()
{
  // We do this to explicitly disable the noise corrector for MC
  GetESDFixer().SetRecoNoiseFactor(5);

  Bool_t ret = AliForwardMultiplicityBase::Book();
  POST(MCAOD_SLOT, &fMCAODFMD);
  POST(PRIM_SLOT,  fPrimary);
  return ret;
}

//____________________________________________________________________
void
AliForwardMCMultiplicityTask::InitMembers(const TAxis& eta, const TAxis& vertex)
{
  // 
  // Initialise the sub objects and stuff.  Called on first event 
  // 
  //
  AliForwardMultiplicityBase::InitMembers(eta, vertex);

  fMCHistos.Init(eta);
  fMCAODFMD.Init(eta);
  fMCRingSums.Init(eta);

  AliForwardUtil::Histos::RebinEta(fPrimary, eta);
  DMSG(fDebug,0,"Primary histogram rebinned to %d,%f,%f eta axis %d,%f,%f", 
       fPrimary->GetXaxis()->GetNbins(), 
       fPrimary->GetXaxis()->GetXmin(),
       fPrimary->GetXaxis()->GetXmax(),
       eta.GetNbins(), 
       eta.GetXmin(),
       eta.GetXmax());


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
Bool_t
AliForwardMCMultiplicityTask::PreEvent()
{
 if (fFirstEvent) 
    fEventInspector.ReadProductionDetails(MCEvent());
  // Clear stuff 
  fHistos.Clear();
  fESDFMD.Clear();
  fAODFMD.Clear();
  fAODEP.Clear();
  fMCHistos.Clear();
  fMCESDFMD.Clear();
  fMCAODFMD.Clear();
  fPrimary->Reset();
  return true;
}
//____________________________________________________________________
Bool_t
AliForwardMCMultiplicityTask::Event(AliESDEvent& esd)
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  
  MAKE_SW(total);
  MAKE_SW(individual);
  START_SW(total);
  START_SW(individual);

  // Read production details 
    
  // Get the input data 
  AliMCEvent*  mcEvent = MCEvent();
  if (!mcEvent || !mcEvent->Stack()) return false;

  Bool_t   lowFlux   = kFALSE;
  UInt_t   triggers  = 0;
  UShort_t ivz       = 0;
  TVector3 ip(1024, 1024, 0);
  Double_t cent      = -1;
  UShort_t nClusters = 0;
  UInt_t   found     = fEventInspector.Process(&esd, triggers, lowFlux, 
					       ivz, ip, cent, nClusters);
  UShort_t ivzMC    = 0;
  TVector3 ipMC(1024, 1024, 0);  
  Double_t phiR     = 0;
  Double_t b        = 0;
  Double_t cMC      = 0;
  Int_t    npart    = 0;
  Int_t    nbin     = 0;
  // UInt_t   foundMC  = 
  fEventInspector.ProcessMC(mcEvent, triggers, ivzMC, ipMC, b, cMC,
			    npart, nbin, phiR);
  fEventInspector.CompareResults(ip.Z(), ipMC.Z(), cent, cMC, b, npart, nbin);
  // fMultEventClassifier.Process(&esd,&fAODRef);
  FILL_SW(individual,kTimingEventInspector);
  
  //Store all events
  MarkEventForStore();
  
  Bool_t isAccepted = true;
  if (found & AliFMDEventInspector::kNoEvent)    {
    fHStatus->Fill(kStatusNoEvent);
    isAccepted = false;
    // return;
  }
  if (found & AliFMDEventInspector::kNoTriggers) {
    fHStatus->Fill(kStatusNoTrigger);
    isAccepted = false; 
    // return;
  }
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
  
  // Disable this check on SPD - will bias data 
  // if (found & AliFMDEventInspector::kNoSPD)  isAccepted = false; // return;
  if (found & AliFMDEventInspector::kNoFMD)     {
    fHStatus->Fill(kStatusNoFMD);
    isAccepted = false; 
    // return;
  }
  if (found & AliFMDEventInspector::kNoVertex)  {
    fHStatus->Fill(kStatusNoVertex);
    isAccepted = false; 
    // return;
  }

  if (isAccepted) {
    fAODFMD.SetIpZ(ip.Z());
    fMCAODFMD.SetIpZ(ip.Z());
  }
  if (found & AliFMDEventInspector::kBadVertex) isAccepted = false; // return;

  // We we do not want to use low flux specific code, we disable it here. 
  if (!fEnableLowFlux) lowFlux = false;

  // Get FMD data 
  AliESDFMD*  esdFMD  = esd.GetFMDData();
  
  // Fix up the the ESD 
  GetESDFixer().Fix(*esdFMD, ip);

  // Apply the sharing filter (or hit merging or clustering if you like)
  START_SW(individual);
  if (isAccepted && !fSharingFilter.Filter(*esdFMD, lowFlux, fESDFMD,ip.Z())){
    AliWarning("Sharing filter failed!");
    fHStatus->Fill(kStatusFailSharing);
    return false;
  }
  if (!fSharingFilter.FilterMC(*esdFMD, *mcEvent, ip,fMCESDFMD,fPrimary)){
    AliWarning("MC Sharing filter failed!");
    fHStatus->Fill(kStatusFailSharing);
    return false;
  }

  // Store some MC parameters in corners of histogram :-)
  fPrimary->SetBinContent(0,                      0,                ipMC.Z());
  fPrimary->SetBinContent(fPrimary->GetNbinsX()+1,0,                    phiR);
  fPrimary->SetBinContent(fPrimary->GetNbinsX()+1,fPrimary->GetNbinsY(),cMC);
  

  if (!isAccepted) {
    // Exit on MC event w/o trigger, vertex, data - since there's no more 
    // to be done for MC 
    FILL_SW(individual,kTimingSharingFilter);
    return false; 
  }
  
  //MarkEventForStore();
  fSharingFilter.CompareResults(fESDFMD, fMCESDFMD);
  FILL_SW(individual,kTimingSharingFilter);


  // Calculate the inclusive charged particle density 
  START_SW(individual);
  if (!fDensityCalculator.Calculate(fESDFMD, fHistos, lowFlux, cent, ip)) { 
    AliWarning("Density calculator failed!");
    fHStatus->Fill(kStatusFailDensity);
    return false;
  }
  if (!fDensityCalculator.CalculateMC(fMCESDFMD, fMCHistos)) { 
    AliWarning("MC Density calculator failed!");
    fHStatus->Fill(kStatusFailDensity);
    return false;
  }
  fDensityCalculator.CompareResults(fHistos, fMCHistos);
  FILL_SW(individual,kTimingDensityCalculator);
  
  if (fEventInspector.GetCollisionSystem() == AliFMDEventInspector::kPbPb) {
    START_SW(individual);
    if (!fEventPlaneFinder.FindEventplane(&esd, fAODEP, 
					  &(fAODFMD.GetHistogram()), &fHistos)){
      AliWarning("Eventplane finder failed!");
      fHStatus->Fill(kStatusFailEventPlane);
    } 
    FILL_SW(individual,kTimingEventPlaneFinder);   
  }

  // Do the secondary and other corrections. 
  START_SW(individual);
  if (!fCorrections.Correct(fHistos, ivz)) { 
    AliWarning("Corrections failed");
    fHStatus->Fill(kStatusFailCorrector);
    return false;
  }
  if (!fCorrections.CorrectMC(fMCHistos, ivz)) { 
    AliWarning("MC Corrections failed");
    fHStatus->Fill(kStatusFailCorrector);
    return false;
  }
  fCorrections.CompareResults(fHistos, fMCHistos);
  FILL_SW(individual,kTimingCorrections);

  // Collect our 'super' histogram
  Bool_t add = fAODFMD.IsTriggerBits(AliAODForwardMult::kInel);
  START_SW(individual);
  if (!fHistCollector.Collect(fHistos, 
			      fRingSums, 
			      ivz, 
			      fAODFMD.GetHistogram(),
			      fAODFMD.GetCentrality(),
			      false,
			      add)) {
    AliWarning("Histogram collector failed");
    fHStatus->Fill(kStatusFailCollector);
    return false;
  }
  if (!fHistCollector.Collect(fMCHistos, 
			      fMCRingSums, 
			      ivz, 
			      fMCAODFMD.GetHistogram(), 
			      -1, 
			      true,
			      add)) {
    AliWarning("MC Histogram collector failed");
    fHStatus->Fill(kStatusFailCollector);
    return false;
  }
  FILL_SW(individual,kTimingHistCollector);
#if 0
  // Copy underflow bins to overflow bins - always full phi coverage 
  TH2&  hMC  = fMCAODFMD.GetHistogram();
  Int_t nEta = hMC.GetNbinsX();
  Int_t nY   = hMC.GetNbinsY();
  for (Int_t iEta = 1; iEta <= nEta; iEta++) {
    hMC.SetBinContent(iEta, nY+1, hMC.GetBinContent(iEta, 0));
  }
#endif

  if (add) {
    fHData->Add(&(fAODFMD.GetHistogram()));
    fHStatus->Fill(kStatusAllThrough);
  }
  else {
    fHStatus->Fill(kStatusNotAdded);
  }
  FILL_SW(total,kTimingTotal);

  return true;
}

//____________________________________________________________________
Bool_t
AliForwardMCMultiplicityTask::PostEvent()
{
  Bool_t ret = AliForwardMultiplicityBase::PostEvent();
  POST(MCAOD_SLOT, &fMCAODFMD);
  POST(PRIMARY_SLOT, fPrimary);
  return ret;
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
