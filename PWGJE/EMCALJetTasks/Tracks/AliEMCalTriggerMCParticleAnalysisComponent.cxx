/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/*
 * Analysis component for simple Monte-Carlo particles. Loops over all particles, selects
 * those which are physical primary, and fills a THnSparse with them.
 *
 *   Author: Markus Fasel
 */
#include <TAxis.h>
#include <TBinning.h>
#include <TMath.h>
#include <THistManager.h>

#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliVEvent.h"
#include "AliVVertex.h"

#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerKineCuts.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerMCParticleAnalysisComponent.h"
#include "AliEMCalTriggerWeightHandler.h"


ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerMCParticleAnalysisComponent)

using namespace PWGJE::EMCALJetTasks;

//______________________________________________________________________________
AliEMCalTriggerMCParticleAnalysisComponent::AliEMCalTriggerMCParticleAnalysisComponent() :
  AliEMCalTriggerTracksAnalysisComponent()
{
  /*
   * Dummy Constructor for I/O
   */
}

//______________________________________________________________________________
AliEMCalTriggerMCParticleAnalysisComponent::AliEMCalTriggerMCParticleAnalysisComponent(const char* name) :
  AliEMCalTriggerTracksAnalysisComponent(name)
{
  /*
   * Main Constructor, to be called by the users
   */
}

  //______________________________________________________________________________
void AliEMCalTriggerMCParticleAnalysisComponent::CreateHistos() {
  /*
   * Create histograms for the MC truth analysis component
   */

  const TBinning *ptbinning = fBinning->GetBinning("pt"),
      *etabinning = fBinning->GetBinning("eta"),
      *phibinning = fBinning->GetBinning("phi"),
      *vertexbinning = fBinning->GetBinning("zvertex");
  AliEMCalTriggerTracksAnalysisComponent::CreateHistos();
  const TAxis *trackaxes[4] = {
      DefineAxis("pt", *ptbinning),
      DefineAxis("eta", *etabinning),
      DefineAxis("phi", *phibinning),
      DefineAxis("zvertex", *vertexbinning)
  };
  fHistos->CreateTHnSparse("hMCtrueParticles", "Particle-based histogram for MC-true particles", 4, trackaxes, "s");
}

//______________________________________________________________________________
void AliEMCalTriggerMCParticleAnalysisComponent::Process(const AliEMCalTriggerEventData* const data) {
  /*
   * Run event loop
   */
  AliMCEvent *mc = data->GetMCEvent();
  if(!mc) return;
  AliVEvent *rec = data->GetRecEvent();
  Double_t values[4];
  double weight = 1.;
  if(fWeightHandler && data->GetMCEvent()){
    weight = fWeightHandler->GetEventWeight(data->GetMCEvent());
  }
  for(int itrk = 0; itrk < mc->GetNumberOfTracks(); itrk++){
    AliVParticle *track = mc->GetTrack(itrk);
    if(!track->Charge()) continue;
    if(!IsPhysicalPrimary(track, mc)) continue;
    if(!fKineCuts->IsSelected(track)) continue;

    values[0] = TMath::Abs(track->Pt());
    values[1] = track->Eta();
    values[2] = track->Phi();
    values[3] = rec->GetPrimaryVertex()->GetZ();
    fHistos->FillTHnSparse("hMCtrueParticles", values, weight);
  }
}

//______________________________________________________________________________
bool AliEMCalTriggerMCParticleAnalysisComponent::IsPhysicalPrimary(const AliVParticle* const part, const AliMCEvent* const mcev) const {
	/*
	 * Transparent access to Physical primary information
	 */
	const AliAODMCParticle *aodpart = dynamic_cast<const AliAODMCParticle *>(part);
	if(aodpart){
		return aodpart->IsPhysicalPrimary();
	}
	return mcev->IsPhysicalPrimary(part->GetLabel());                 // No abs value since default value is -1
}
