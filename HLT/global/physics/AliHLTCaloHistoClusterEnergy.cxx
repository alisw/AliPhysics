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
#include "TVector3.h"
#include "TMath.h"

ClassImp(AliHLTCaloHistoClusterEnergy);

AliHLTCaloHistoClusterEnergy::AliHLTCaloHistoClusterEnergy(TString det) :
  AliHLTCaloHistoProducer(),
  fHistClusterEnergy(NULL),
  fHistClusterEnergyVsNCells(NULL),
  fHistClusterEnergyDepositEtaPhi(NULL)
{
  // See header file for documentation
  fHistClusterEnergy = new TH1F(Form("%s fHistClusterEnergy", det.Data()), Form("%s Distribution of total energy in clusters", det.Data()), 5000, 0, 100);
  fHistClusterEnergy->GetXaxis()->SetTitle("E GeV");
  fHistClusterEnergy->GetYaxis()->SetTitle("Number of counts");
  fHistArray->AddLast(fHistClusterEnergy);

  fHistClusterEnergyVsNCells = new TH2F(Form("%s fHistClusterEnergyVsNCells", det.Data()), Form("%s Distribution of Energy vs Number of Cells in cluster", det.Data()), 1000, 0, 100, 50, 0 , 50);
  fHistClusterEnergyVsNCells->GetXaxis()->SetTitle("Energy in cluster (GeV)");
  fHistClusterEnergyVsNCells->GetYaxis()->SetTitle("Number of Cells in cluster");
  fHistArray->AddLast(fHistClusterEnergyVsNCells);
  
  Float_t phiMin = 0.;
  Float_t phiMax = 2*TMath::Pi();
  Float_t etaMin = -1.;
  Float_t etaMax = 1.;
  
  if(det == "PHOS")
  {
     phiMin = 255.0/180.*TMath::Pi();
     phiMax = 325.0/180.*TMath::Pi();
     etaMin = -0.13;
     etaMax = 0.13;
  }
  
  fHistClusterEnergyDepositEtaPhi = new TH2F(Form("%s fHistClusterEnergyDepositedEtaPhi", det.Data()), Form("%s Amount of energy deposited in Phi vs Eta", det.Data()), 200, phiMin, phiMax, 50, etaMin , etaMax);
  fHistClusterEnergyDepositEtaPhi->GetXaxis()->SetTitle("#phi");
  fHistClusterEnergyDepositEtaPhi->GetYaxis()->SetTitle("#eta");
  fHistArray->AddLast(fHistClusterEnergyDepositEtaPhi);

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
    FillClusterEnergyHistos(cluster);
  }
  return 0;
}

Int_t AliHLTCaloHistoClusterEnergy::FillHistograms(Int_t nc, vector<AliHLTCaloClusterDataStruct*> &cVec) {
  //See header file for documentation
  // HLTInfo("histo");
  for(int ic = 0; ic < nc; ic++) {
    AliHLTCaloClusterDataStruct * cluster = cVec.at(ic);
    FillClusterEnergyHistos(cluster);
  }
  return 0;
}

template <class T>
Int_t AliHLTCaloHistoClusterEnergy::FillClusterEnergyHistos(T* cluster) {
  fHistClusterEnergy->Fill(cluster->E());
  fHistClusterEnergyVsNCells->Fill(cluster->E(), cluster->GetNCells());
  
  Float_t pos[3];
  cluster->GetPosition(pos);
  TVector3 vec(pos);
  
  // Stupid hack, too tired to fix
  if(vec.Phi() < 0)  fHistClusterEnergyDepositEtaPhi->Fill(2*TMath::Pi() + vec.Phi(), vec.Eta(), cluster->E());
  else fHistClusterEnergyDepositEtaPhi->Fill(vec.Phi(), vec.Eta(), cluster->E());
  
  return 0;
}
