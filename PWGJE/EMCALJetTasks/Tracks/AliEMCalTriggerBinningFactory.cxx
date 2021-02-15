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
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerBinningFactory.h"
#include <map>
#include <vector>
#include <TMath.h>
#include <TArrayD.h>

using namespace PWGJE::EMCALJetTasks;

AliEMCalTriggerBinningFactory::AliEMCalTriggerBinningFactory() {
  /*
   * See header file for details
   */
}

void AliEMCalTriggerBinningFactory::Create(AliEMCalTriggerBinningComponent* const data) {
  /*
   * See header file for details
   */
  if(!data->GetBinning("pt")) data->SetBinning("pt", new DefaultPtBinning);
  if(!data->GetBinning("eta")) data->SetBinning("eta", new DefaultEtaBinning);
  if(!data->GetBinning("phi")) data->SetBinning("phi", new TLinearBinning(100, 0., 2*TMath::Pi()));
  if(!data->GetBinning("zvertex")) data->SetBinning("zvertex", new DefaultZVertexBinning);
  if(!data->GetBinning("centrality")) data->SetBinning("centrality", new TLinearBinning(5, 0., 100.));
}

AliEMCalTriggerBinningFactory::MarkusPtBinning::MarkusPtBinning():
  TCustomBinning()
{
  /*
   * See header file for details
   */
  SetMinimum(0);
  AddStep(2.5, 0.1);
  AddStep(7., 0.25);
  AddStep(10., 0.5);
  AddStep(15., 1.);
  AddStep(20., 2.5);
  AddStep(30., 5.);
  AddStep(100., 10.);
  AddStep(200., 20.);
}

AliEMCalTriggerBinningFactory::DefaultPtBinning::DefaultPtBinning() :
    TCustomBinning()
{
  /*
   * See header file for details
   */
  SetMinimum(0);
  AddStep(1, 0.05);
  AddStep(2, 0.1);
  AddStep(4, 0.2);
  AddStep(7, 0.5);
  AddStep(16, 1);
  AddStep(36, 2);
  AddStep(40, 4);
  AddStep(50, 5);
  AddStep(100, 10);
  AddStep(200, 20);
}
