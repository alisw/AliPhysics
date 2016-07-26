//
// Task to propagate tracks to EMCAL surface.
//
// Author: C.Loizides, S.Aiola

#include "AliEmcalTrackPropagatorTask.h"

#include "AliParticleContainer.h"

#include <TClonesArray.h>

#include <AliVTrack.h>
#include <AliEMCALRecoUtils.h>

ClassImp(AliEmcalTrackPropagatorTask)

//________________________________________________________________________
AliEmcalTrackPropagatorTask::AliEmcalTrackPropagatorTask() : 
  AliAnalysisTaskEmcal("AliEmcalTrackPropagatorTask", kFALSE),
  fDist(440),
  fOnlyIfNotSet(kTRUE),
  fOnlyIfEmcal(kTRUE)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalTrackPropagatorTask::AliEmcalTrackPropagatorTask(const char *name) : 
  AliAnalysisTaskEmcal(name, kFALSE),
  fDist(440),
  fOnlyIfNotSet(kTRUE),
  fOnlyIfEmcal(kTRUE)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalTrackPropagatorTask::~AliEmcalTrackPropagatorTask()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalTrackPropagatorTask::ExecOnce() 
{
  // Executed only once at the beginning of the analysis.

  AliParticleContainer* tracks = GetParticleContainer(0);
  if (tracks) {
    tracks->SetClassName("AliVTrack");
  }

  AliAnalysisTaskEmcal::ExecOnce();
}

//________________________________________________________________________
Bool_t AliEmcalTrackPropagatorTask::Run() 
{
  // Main loop, called for each event.

  AliParticleContainer* tracks = GetParticleContainer(0);

  if (!tracks) return 0;
  
  tracks->ResetCurrentID();
  AliVTrack* track = 0;
  while ((track = static_cast<AliVTrack*>(tracks->GetNextAcceptParticle()))) {
    if (fOnlyIfNotSet && track->IsExtrapolatedToEMCAL()) continue;
    if (fOnlyIfEmcal && !track->IsEMCAL()) continue;
    
    AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(track, fDist);
  }

  return kTRUE;
}
