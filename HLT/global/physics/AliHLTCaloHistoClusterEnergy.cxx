/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Svein Lindal <slindal@fys.uio.no>                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** 
 * @file   AliHLTCaloHistoClusterEnergy
 * @author Svein Lindal
 * @date 
 * @brief  Produces histograms of cluster energy distributions 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTCaloHistoClusterEnergy.h"
#include "AliESDCaloCluster.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRefArray.h"

ClassImp(AliHLTCaloHistoClusterEnergy);

AliHLTCaloHistoClusterEnergy::AliHLTCaloHistoClusterEnergy(TString det) :
  AliHLTCaloHistoProducer(),
  fHistClusterEnergy(NULL),
  fHistClusterEnergyVsNCells(NULL)
{
  // See header file for documentation
  fHistClusterEnergy = new TH1F(Form("%s fHistClusterEnergy", det.Data()), Form("%s Distribution of total energy in clusters", det.Data()), 200, 0, 200);
  fHistClusterEnergy->GetXaxis()->SetTitle("E GeV");
  fHistClusterEnergy->GetYaxis()->SetTitle("Number of counts");
  fHistClusterEnergy->SetMarkerStyle(21);
  fHistArray->AddLast(fHistClusterEnergy);

  fHistClusterEnergyVsNCells = new TH2F(Form("%s fHistClusterEnergyVsNCells", det.Data()), Form("%s Distribution of Energy vs Number of Cells in cluster", det.Data()), 200, 0, 200, 50, 0 , 50);
  fHistClusterEnergyVsNCells->GetXaxis()->SetTitle("Energy in cluster (GeV)");
  fHistClusterEnergyVsNCells->GetYaxis()->SetTitle("Number of Cells in cluster");
  fHistClusterEnergyVsNCells->SetMarkerStyle(21);
  fHistArray->AddLast(fHistClusterEnergyVsNCells);

}

AliHLTCaloHistoClusterEnergy::~AliHLTCaloHistoClusterEnergy()
{
  //destructor
  if(fHistClusterEnergy)
    delete fHistClusterEnergy;
  fHistClusterEnergy = NULL;

  if(fHistClusterEnergyVsNCells)
    delete fHistClusterEnergyVsNCells;
  fHistClusterEnergyVsNCells = NULL;
}

Int_t AliHLTCaloHistoClusterEnergy::FillHistograms(Int_t nc, TRefArray * clusterArray) {
  //See header file for documentation
  for(int ic = 0; ic < nc; ic++) {
    AliESDCaloCluster * cluster = static_cast<AliESDCaloCluster*>(clusterArray->At(ic));
    fHistClusterEnergy->Fill(cluster->E());
    fHistClusterEnergyVsNCells->Fill(cluster->GetNCells(), cluster->E());
  }
  return 0;
}
Int_t AliHLTCaloHistoClusterEnergy::FillHistograms(Int_t nc, vector<AliHLTCaloClusterDataStruct*> &cVec) {
  //See header file for documentation
  HLTInfo("histo");
  for(int ic = 0; ic < nc; ic++) {
    AliHLTCaloClusterDataStruct * cluster = cVec.at(ic);
    fHistClusterEnergy->Fill(cluster->E());
    fHistClusterEnergyVsNCells->Fill(cluster->fNCells, cluster->E());
  }
  return 0;
}
