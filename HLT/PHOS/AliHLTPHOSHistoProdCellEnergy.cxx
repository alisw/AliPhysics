
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
 * @file   AliHLTPHOSHistoProdCellEnergy
 * @author Albin Gaignette & Svein Lindal
 * @date 
 * @brief  Produces histograms of cluster energy distributions 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSHistoProdCellEnergy.h"
//#include "AliESDCaloCluster.h"
#include "TMath.h"

#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"
#include "TObjArray.h"
//#include "TClonesArray.h"
//#include <iostream>
#include "TH1F.h"
#include "TH2F.h"

AliHLTPHOSHistoProdCellEnergy::AliHLTPHOSHistoProdCellEnergy() :
  fClusterReader(NULL),
  fHistCellEnergy(NULL),
  fHistCellEnergyVsNCells(NULL),
  fHistArrayPtr(NULL)
{
  // See header file for documentation
  fHistArrayPtr = new TObjArray;
  fClusterReader = new AliHLTCaloClusterReader();

  fHistCellEnergy = new TH1F("fHistCellEnergy", "Distribution of total energy in clusters", 200, 0, 1);
  fHistCellEnergy->GetXaxis()->SetTitle("E GeV");
  fHistCellEnergy->GetYaxis()->SetTitle("Number of counts");
  fHistCellEnergy->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistCellEnergy);

  fHistCellEnergyVsNCells = new TH2F("fHistCellEnergyVsNCells", "Distribution of Energy vs Number of Cells in cluster", 200, 0, 200, 50, 0 , 50);
  fHistCellEnergyVsNCells->GetXaxis()->SetTitle("Energy in cluster (GeV)");
  fHistCellEnergyVsNCells->GetYaxis()->SetTitle("Number of Cells in cluster");
  fHistCellEnergyVsNCells->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistCellEnergyVsNCells);


}

AliHLTPHOSHistoProdCellEnergy::~AliHLTPHOSHistoProdCellEnergy()
{
  if(fHistCellEnergy){
      delete fHistCellEnergy;
      fHistCellEnergy = 0;
    }
}

TObjArray* AliHLTPHOSHistoProdCellEnergy::GetHistograms()
{  
  // See header file for documentation

  return fHistArrayPtr;
}

Int_t AliHLTPHOSHistoProdCellEnergy::DoEvent(AliHLTCaloClusterHeaderStruct* cHeader) {   
  
  fClusterReader->SetMemory(cHeader);
  AliHLTCaloClusterDataStruct* cluster;
  while ( ( cluster = fClusterReader->NextCluster() ) ) {
    fHistCellEnergy->Fill(cluster->fNCells);
    fHistCellEnergyVsNCells->Fill(cluster->fNCells, cluster->fCellsAmpFraction);
  }  
  return 0;
}
  
 
 

