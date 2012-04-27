#include "AliAnalysisTaskTransverseEnergy.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"

//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for Et tasks
//  - reconstruction and MonteCarlo output
//
//*-- Authors: Oystein Djuvsland (Bergen)
//_________________________________________________________________________//
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "TH2F.h"
#include <iostream>


ClassImp(AliAnalysisTaskTransverseEnergy)

AliAnalysisTaskTransverseEnergy::AliAnalysisTaskTransverseEnergy(const char* name, Bool_t isMc) :
        AliAnalysisTaskSE(name)
	,fESDEvent(0)
        ,fMCConfigFile("ConfigEtMonteCarlo.C")
        ,fRecoConfigFile("ConfigEtReconstructed.C")
        ,fHistEtRecvsEtMC(0)
	,fHistEtRecOverEtMC(0)
	,fHistDiffEtRecEtMCOverEtMC(0)
        ,fEsdtrackCutsITSTPC(0)
        ,fEsdtrackCutsTPC(0)
        ,fEsdtrackCutsITS(0)
        ,fOutputList(0)
        ,fCentSelTaskName("centralityTask")
	,fIsMc(isMc)
	,fCurrentRunNum(-1)
{
  // Constructor
}

AliAnalysisTaskTransverseEnergy::~AliAnalysisTaskTransverseEnergy()
{    // destructor
  delete fHistEtRecvsEtMC;
  delete fHistEtRecOverEtMC;
  delete fEsdtrackCutsITSTPC;
  delete fEsdtrackCutsTPC;
  delete fEsdtrackCutsITS;
  delete fOutputList;
}


AliCentrality* AliAnalysisTaskTransverseEnergy::GetCentralityObject()
{
  // See header file for class documentation
    if (fESDEvent)return fESDEvent->GetCentrality();
    else return 0;
}


