
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Albin Gaignette
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
 * @file   AliHLTPHOSHistoProdClusterEnergy
 * @author Albin Gaignette & Svein Lindal
 * @date 
 * @brief  Produces histograms of cluster energy distributions 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSHistoProdClusterEnergy.h"
//#include "AliESDCaloCluster.h"
#include "TMath.h"

#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"
#include "TObjArray.h"
//#include "TClonesArray.h"
//#include <iostream>
#include "TH1F.h"
#include "TH2F.h"

AliHLTPHOSHistoProdClusterEnergy::AliHLTPHOSHistoProdClusterEnergy() :
  fClusterReader(NULL),
  fHistClusterEnergy(NULL),
  fHistClusterEnergyVsNCells(NULL),
  fHistArrayPtr(NULL)
{
  // See header file for documentation
  fHistArrayPtr = new TObjArray;
  fClusterReader = new AliHLTCaloClusterReader();

  fHistClusterEnergy = new TH1F("fHistClusterEnergy", "Distribution of total energy in clusters", 200, 0, 1);
  fHistClusterEnergy->GetXaxis()->SetTitle("E GeV");
  fHistClusterEnergy->GetYaxis()->SetTitle("Number of counts");
  fHistClusterEnergy->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistClusterEnergy);

  fHistClusterEnergyVsNCells = new TH2F("fHistClusterEnergyVsNCells", "Distribution of Energy vs Number of Cells in cluster", 200, 0, 200, 50, 0 , 50);
  fHistClusterEnergyVsNCells->GetXaxis()->SetTitle("Energy in cluster (GeV)");
  fHistClusterEnergyVsNCells->GetYaxis()->SetTitle("Number of Cells in cluster");
  fHistClusterEnergyVsNCells->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistClusterEnergyVsNCells);


}

AliHLTPHOSHistoProdClusterEnergy::~AliHLTPHOSHistoProdClusterEnergy()
{
  if(fHistClusterEnergy){
      delete fHistClusterEnergy;
      fHistClusterEnergy = 0;
    }
}

TObjArray* AliHLTPHOSHistoProdClusterEnergy::GetHistograms()
{  
  // See header file for documentation

  return fHistArrayPtr;
}

Int_t AliHLTPHOSHistoProdClusterEnergy::DoEvent(AliHLTCaloClusterHeaderStruct* cHeader) {   
  
  fClusterReader->SetMemory(cHeader);
  
  int ncls = cHeader->fNClusters;
  Float_t* cPos[ncls];
  Float_t cEnergy[ncls];
  
  AliHLTCaloClusterDataStruct* cluster;
  Int_t icls = 0;
  while ( ( cluster = fClusterReader->NextCluster() ) ) {
    
    cPos[icls] = cluster->fGlobalPos;
    cEnergy[icls] = cluster->fEnergy; 
    
    icls++;
  }  
  
  for(Int_t ipho = 0; ipho<(ncls-1); ipho++) { 
    for(Int_t jpho = ipho+1 ; jpho<ncls ; jpho++) { 
      // Calcul of the theta angle between two photons
      Double_t theta = (2* asin(0.5*TMath::Sqrt((cPos[ipho][0]-cPos[jpho][0])*(cPos[ipho][0]-cPos[jpho][0]) +(cPos[ipho][1]-cPos[jpho][1])*(cPos[ipho][1]-cPos[jpho][1]))/460));
      
      // Calcul of the mass m of the pion 
      Double_t m =(TMath::Sqrt(2 * cEnergy[ipho]* cEnergy[jpho]*(1-TMath::Cos(theta))));
      
      fHistClusterEnergy->Fill(m);
    }
  }
  
  return 0;
}
  
 
 

