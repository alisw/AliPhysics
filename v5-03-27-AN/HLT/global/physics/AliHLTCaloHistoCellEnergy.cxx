/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Albin Gaignette, Svein Lindal slindal@fys.uio.no      *
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
 * @file   AliHLTCaloHistoCellEnergy
 * @author Svein Lindal
 * @date 
 * @brief  Produces histograms of cluster energy distributions 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTCaloHistoCellEnergy.h"
#include "TObjArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRefArray.h"
#include "TString.h"
#include "AliESDCaloCluster.h"

ClassImp(AliHLTCaloHistoCellEnergy);

AliHLTCaloHistoCellEnergy::AliHLTCaloHistoCellEnergy(TString det) :
  fHistCellEnergy(NULL),
  fHistCellEnergyVsNCells(NULL),
  fHistArrayPtr(NULL)
{
  // See header file for documentation
  fHistArrayPtr = new TObjArray;

  fHistCellEnergy = new TH1F(Form("%s fHistCellEnergy", det.Data()), Form("%s Distribution of total energy in clusters", det.Data()), 200, 0, 1);
  fHistCellEnergy->GetXaxis()->SetTitle("E GeV");
  fHistCellEnergy->GetYaxis()->SetTitle("Number of counts");
  fHistCellEnergy->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistCellEnergy);

  fHistCellEnergyVsNCells = new TH2F(Form("%s fHistCellEnergyVsNCells", det.Data()), Form("%s Distribution of Energy vs Number of Cells in cluster", det.Data()), 200, 0, 200, 50, 0 , 50);
  fHistCellEnergyVsNCells->GetXaxis()->SetTitle("Energy in cluster (GeV)");
  fHistCellEnergyVsNCells->GetYaxis()->SetTitle("Number of Cells in cluster");
  fHistCellEnergyVsNCells->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistCellEnergyVsNCells);

}

AliHLTCaloHistoCellEnergy::~AliHLTCaloHistoCellEnergy()
{
  //See header file for documentation

  if(fHistCellEnergy)
    delete fHistCellEnergy;
  fHistCellEnergy = NULL;
  
  if(fHistCellEnergyVsNCells)
    delete fHistCellEnergyVsNCells;
  fHistCellEnergyVsNCells = NULL;

  if(fHistArrayPtr)
    delete fHistArrayPtr;
  fHistArrayPtr = NULL;

}

TObjArray* AliHLTCaloHistoCellEnergy::GetHistograms()
{  
  // See header file for documentation
  return fHistArrayPtr;
}


Int_t AliHLTCaloHistoCellEnergy::FillHistograms(Int_t nc, TRefArray * clustersArray) {   
  
  for(int ic = 0; ic < nc; ic++) {
    AliESDCaloCluster * cluster = static_cast<AliESDCaloCluster*>(clustersArray->At(ic));
    for(int i = 0; i < cluster->GetNCells(); i++) {
      fHistCellEnergyVsNCells->Fill(cluster->GetNCells(), cluster->GetCellAmplitudeFraction(i));
      fHistCellEnergy->Fill(cluster->GetCellAmplitudeFraction(i));
    }
  }  
  
  return 0;
}
  
 
 

