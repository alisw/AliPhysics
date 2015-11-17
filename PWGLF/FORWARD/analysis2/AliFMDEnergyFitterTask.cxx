// 
// Histogram and fit the energy loss distributions for the FMD
// 
// Inputs: 
//   - AliESDEvent 
//
// Outputs: 
//   - None
// 
// Histograms:
//   
// Corrections used:
//   - None
// 
// 
//
#include "AliFMDEnergyFitterTask.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDFMD.h"
#include "AliMCEvent.h"
#include "AliAODForwardMult.h"
#include "AliAnalysisManager.h"
#include "AliForwardCorrectionManager.h"
#include <TH1.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <iostream>

//====================================================================
AliFMDEnergyFitterTask::AliFMDEnergyFitterTask()
  : AliBaseESDTask(),
    fEventInspector(),				
    fESDFixer(),
    fEnergyFitter(),
    fOnlyMB(false)
{
  // 
  // Constructor
  //
  DGUARD(fDebug, 3,"Default CTOR of AliFMDEnergyFitterTask");
  fCloneList = true;
}

//____________________________________________________________________
AliFMDEnergyFitterTask::AliFMDEnergyFitterTask(const char* name)
  : AliBaseESDTask(name, "AliFMDEnergyFitterTask", 
		   &(AliForwardCorrectionManager::Instance())), 
    fEventInspector("event"),
    fESDFixer("esdFizer"),
    fEnergyFitter("energy"),
    fOnlyMB(false)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
  DGUARD(fDebug, 3,"Named CTOR of AliFMDEnergyFitterTask: %s", name);
  fCloneList = true;
}


//____________________________________________________________________
void
AliFMDEnergyFitterTask::SetDebug(Int_t dbg)
{
  // 
  // Set the debug level 
  // 
  // Parameters:
  //    dbg Debug level
  //
  AliBaseESDTask::SetDebug(dbg);
  fEnergyFitter.SetDebug(dbg);
}
//____________________________________________________________________
TAxis*
AliFMDEnergyFitterTask::DefaultEtaAxis() const
{
  static TAxis* a = new TAxis(0, 0, 0);
  return a;
}
//____________________________________________________________________
TAxis*
AliFMDEnergyFitterTask::DefaultVertexAxis() const
{
  static TAxis* a = new TAxis(10, -10, 10);
  return a;
}

//____________________________________________________________________
Bool_t
AliFMDEnergyFitterTask::Setup()
{
  fEnergyFitter.Init();
  return true;
}

//____________________________________________________________________
Bool_t
AliFMDEnergyFitterTask::Book()
{
  // 
  // Create output objects 
  // 
  //
  DGUARD(fDebug,1,"Create output objects of AliFMDEnergyFitterTask");

  // We don't need any corrections for this task 
  fNeededCorrections = 0; 
  fExtraCorrections  = 0;
  if (fESDFixer.IsUseNoiseCorrection()) 
    fNeededCorrections = AliForwardCorrectionManager::kNoiseGain;

  fESDFixer    .CreateOutputObjects(fList);
  fEnergyFitter.CreateOutputObjects(fList);

  fList->Add(AliForwardUtil::MakeParameter("onlyMB", fOnlyMB));
  return true;
}
//____________________________________________________________________
void
AliFMDEnergyFitterTask::PreCorrections(const AliESDEvent* esd)
{
  if (!esd) return; 
  
  AliESDFMD* esdFMD = esd->GetFMDData();  
  if (!esdFMD) return;

  // TODO: We should always disable this on MC!
  Int_t tgt = fESDFixer.FindTargetNoiseFactor(*esdFMD, false);
  if (tgt <= 0) {
    // If the target noise factor is 0 or less, disable the noise/gain
    // correction.
    fESDFixer.SetRecoNoiseFactor(4);
    fNeededCorrections ^= AliForwardCorrectionManager::kNoiseGain;
  }
  else 
    AliWarning("The noise corrector has been enabled!");
}

//____________________________________________________________________
Bool_t
AliFMDEnergyFitterTask::PreData(const TAxis& /*vertex*/, const TAxis& eta)
{
  // 
  // Initialise the sub objects and stuff.  Called on first event 
  // 
  //
  DGUARD(fDebug,1,"Initialize subs of AliFMDEnergyFitterTask");

  fEnergyFitter.SetupForData(eta,GetEventInspector().GetCollisionSystem());

  Print();
  return true;
}

//____________________________________________________________________
Bool_t
AliFMDEnergyFitterTask::Event(AliESDEvent& esd)
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  

  // static Int_t cnt = 0;
  // cnt++;
  // Get the input data 
  DGUARD(fDebug,3,"Analyse event of AliFMDEnergyFitterTask");
  // --- Read in the data --------------------------------------------
  LoadBranches();

  Bool_t   lowFlux   = kFALSE;
  UInt_t   triggers  = 0;
  UShort_t ivz       = 0;
  TVector3 ip;
  Double_t cent      = 0;
  UShort_t nClusters = 0;
  UInt_t   found     = fEventInspector.Process(&esd, triggers, lowFlux, 
					       ivz, ip, cent, nClusters);
  if (found & AliFMDEventInspector::kNoEvent)    return false;
  if (found & AliFMDEventInspector::kNoTriggers) return false;
  if (found & AliFMDEventInspector::kNoSPD)      return false;
  if (found & AliFMDEventInspector::kNoFMD)      return false;
  if (found & AliFMDEventInspector::kNoVertex)   return false;
  if (found & AliFMDEventInspector::kBadVertex)  return false;

  // do not process pile-up, A, C, and E events 
  if (triggers & AliAODForwardMult::kPileUp)     return false;
  if (triggers & AliAODForwardMult::kA)          return false;
  if (triggers & AliAODForwardMult::kC)          return false;
  if (triggers & AliAODForwardMult::kE)          return false;
  
  // We want only the events found by off-line 
  if (!(triggers & AliAODForwardMult::kOffline)) return false;

  // Perhaps we should also insist on MB only 
  if (fOnlyMB && (!(triggers & AliAODForwardMult::kInel))) return false;

  //  if(cent > 0) {
  //  if( cent < 40 || cent >50 ) return;
  //  else std::cout<<"selecting event with cent "<<cent<<std::endl;
  // }
  
  // Get FMD data 
  AliESDFMD* esdFMD = esd.GetFMDData();

  // Fix up ESD 
  fESDFixer.Fix(*esdFMD, ip);

  // Do the energy stuff 
  if (!fEnergyFitter.Accumulate(*esdFMD, cent, 
				triggers & AliAODForwardMult::kEmpty)){
    AliWarning("Energy fitter failed");
    return false;
  }

  return true;
}

//____________________________________________________________________
Bool_t
AliFMDEnergyFitterTask::Finalize()
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //
  DGUARD(fDebug,1,"Processing merged output of AliFMDEnergyFitterTask");

  AliInfo("Fitting energy loss spectra");
  fEnergyFitter.Fit(fResults);

  return true;
}

#define PFB(N,FLAG)				\
  do {									\
    AliForwardUtil::PrintName(N);					\
    std::cout << std::boolalpha << (FLAG) << std::noboolalpha << std::endl; \
  } while(false)
//____________________________________________________________________
void
AliFMDEnergyFitterTask::Print(Option_t* option) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
  AliBaseESDTask::Print(option);
  gROOT->IncreaseDirLevel();
  PFB("Only MB", fOnlyMB);
  fESDFixer    .Print(option);
  fEnergyFitter.Print(option);
  gROOT->DecreaseDirLevel();
}

//
// EOF
//
