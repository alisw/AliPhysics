//-*- Mode: C++ -*-

#include "TMath.h"
#include "TAxis.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisNetParticleBase.h"

using namespace std;

/**
 * Class for NetParticle Distributions
 * -- Base Class
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

ClassImp(AliAnalysisNetParticleBase)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
AliAnalysisNetParticleBase::AliAnalysisNetParticleBase() :
  TNamed(),
  fHelper(NULL),
  fPdgCode(2212),

  fESD(NULL), 
  fESDTrackCuts(NULL),
  fAOD(NULL),
  fArrayMC(NULL),
  fAODtrackCutBit(1024),
  fIsMC(kFALSE),
  fMCEvent(NULL),
  fStack(NULL),
  
  fCentralityBin(-1.),
  fNTracks(0) {
  // Constructor   

  AliLog::SetClassDebugLevel("AliAnalysisNetParticleBase",10);
}

//________________________________________________________________________
AliAnalysisNetParticleBase::AliAnalysisNetParticleBase(const Char_t* name, const Char_t* title) :
  TNamed(name, title),
  fHelper(NULL),
  fPdgCode(2212),

  fESD(NULL), 
  fESDTrackCuts(NULL),
  fAOD(NULL),
  fArrayMC(NULL),
  fAODtrackCutBit(1024),
  fIsMC(kFALSE),
  fMCEvent(NULL),
  fStack(NULL),
  
  fCentralityBin(-1.),
  fNTracks(0) {
  // Constructor   

  AliLog::SetClassDebugLevel("AliAnalysisNetParticleBase",10);
}

//________________________________________________________________________
AliAnalysisNetParticleBase::~AliAnalysisNetParticleBase() {
  // Destructor
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Public Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisNetParticleBase::Initialize(AliAnalysisNetParticleHelper* helper, AliESDtrackCuts* cuts) {
  // -- Initialize
  //    If esdCuts are provided, they are taken otherwise the efficiency ones are used

  //AliESDtrackCuts *cuts, Bool_t isMC, Int_t trackCutBit

  // -- Get Helper
  fHelper           = helper;

  // -- ESD track cuts
  fESDTrackCuts     = (cuts) ? cuts : helper->GetESDTrackCuts();

  // -- Is MC
  fIsMC             = helper->GetIsMC();

  // -- AOD track filter bit
  fAODtrackCutBit   = helper->GetAODtrackCutBit();

  // -- Get particle species / pdgCode
  if (fIsMC)
    fPdgCode        = AliPID::ParticleCode(fHelper->GetParticleSpecies());

  // -- Init within class
  Init();

  // -- Create histograms
  CreateHistograms();

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Setup/Reset Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Int_t AliAnalysisNetParticleBase::SetupEvent() {
  // -- Setup event

  ResetEvent();

  // -- Get ESD objects
  if (dynamic_cast<AliESDInputHandler*>(fHelper->GetInputEventHandler())) {
    fESD         = dynamic_cast<AliESDEvent*>(fHelper->GetInputEventHandler()->GetEvent());
    fNTracks     = fESD->GetNumberOfTracks();
  }

  // -- Get AOD objects
  else if (dynamic_cast<AliAODInputHandler*>(fHelper->GetInputEventHandler())) {
    fAOD         = dynamic_cast<AliAODEvent*>(fHelper->GetInputEventHandler()->GetEvent());
    fNTracks     = fAOD->GetNumberOfTracks();
    
    if (fIsMC) {
      fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!fArrayMC)
	AliFatal("No array of MC particles found !!!"); // MW  no AliFatal use return values
    }
  }

  // -- Get MC objects
  if (fIsMC) {
    fMCEvent     = fHelper->GetMCEvent();
    if (fMCEvent)
      fStack     = fMCEvent->Stack();
  }

  // -- Get CentralityBin
  fCentralityBin = fHelper->GetCentralityBin();

  // -- Setup in class
  // -------------------
  return Setup();
}

//________________________________________________________________________
void AliAnalysisNetParticleBase::ResetEvent() {
  // -- Reset event

  // -- Reset ESD Event
  fESD       = NULL;

  // -- Reset AOD Event
  fAOD       = NULL;

  // -- Reset MC Event
  if (fIsMC)
    fMCEvent = NULL;

  // -- Reset in class
  Reset();

  return;
}


