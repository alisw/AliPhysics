// $Id: AliAnalysisTaskCLQA.cxx 60694 2013-02-04 15:35:56Z morsch $
//
// Constantin's Task
//
// Author: C.Loizides

#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeoParams.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalJet.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliTrackerBase.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliAnalysisTaskCLQA.h"

ClassImp(AliAnalysisTaskCLQA)

//________________________________________________________________________
AliAnalysisTaskCLQA::AliAnalysisTaskCLQA() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskCLQA", kTRUE)
{
  // Default constructor.
}

//________________________________________________________________________
AliAnalysisTaskCLQA::AliAnalysisTaskCLQA(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE)
{
  // Standard constructor.
}

//________________________________________________________________________
AliAnalysisTaskCLQA::~AliAnalysisTaskCLQA()
{
  // Destructor
}

//________________________________________________________________________
void AliAnalysisTaskCLQA::UserCreateOutputObjects()
{
  // Create histograms

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCLQA::RetrieveEventObjects()
{
  // Retrieve event objects.

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskCLQA::FillHistograms()
{
  // Fill histograms.

  return kTRUE;
}
