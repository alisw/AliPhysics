/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliClusterContainer.h"
#include "AliEmcalContainer.h"
#include "AliEmcalTriggerEmulation.h"

ClassImp(PWGJE::EMCALJetTasks::AliEmcalTriggerEmulation)

using namespace PWGJE::EMCALJetTasks;

AliEmcalTriggerEmulation::AliEmcalTriggerEmulation():
     TObject(),
     fThresholdEngine()
{
  /*
   * See header file for details
   */
  fThresholdEngine.SetParameter(0, 1);
  fThresholdEngine.SetParameter(1, 0);
  fThresholdEngine.SetParameter(2, 0.5);
}

AliEmcalTriggerEmulation::AliEmcalTriggerEmulation(double energy, double resolution):
      TObject(),
      fThresholdEngine()
{
  /*
   * See header file for details
   */
  fThresholdEngine.SetParameter(0, 1);
  fThresholdEngine.SetParameter(1, energy);
  fThresholdEngine.SetParameter(2, resolution);
}

bool AliEmcalTriggerEmulation::SelectEvent(const AliEmcalContainer * const cont) const {
  /*
   * See header file for details
   */
  TF1 randomengine = fThresholdEngine; // could be done in a more efficient way
  double threshold = randomengine.GetRandom();

  const AliClusterContainer *clustcont = nullptr;
  if((clustcont = dynamic_cast<const AliClusterContainer *>(cont))){
    bool found = false;
    AliVCluster *clust = nullptr;
    for(auto clusten : clustcont->accepted()){
      clust = static_cast<AliVCluster *>(clusten);
      if(clust->E() > threshold){
        found = true;
        break;
      }
    }
    return found;
  }
  return true;
}
