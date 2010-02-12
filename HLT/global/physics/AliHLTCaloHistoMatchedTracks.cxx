//-*- Mode: C++ -*-
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Albin Gaignette                                       *
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
 * @file   AliHLTPHOSMatchedclustershistoProducer
 * @author Albin Gaignette
 * @date 
 * @brief  Base Class for the Calo Matched track histograms
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTCaloHistoMatchedTracks.h"
#include "AliESDCaloCluster.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"
#include "TObjArray.h"
#include "TH1F.h"
#include "TObjArray.h"
#include "TRefArray.h"
#include "TString.h"

AliHLTCaloHistoMatchedTracks::AliHLTCaloHistoMatchedTracks(TString det) :
  //  fClusterReader(NULL),
  fHistArrayPtr(NULL),
  fHistMatchDistance(NULL),
  fHistMatchedEnergy(NULL),
  fHistUnMatchedEnergy(NULL)
{
  // See header file for documentation
  fHistArrayPtr = new TObjArray;

  //  fClusterReader = new AliHLTCaloClusterReader();
  fHistMatchDistance = new TH1F( Form("%s fHistMatchDistance", det.Data()), Form("%s Track - Cluster residuals (cm)", det.Data()), 50, 0, 50);
  fHistMatchDistance->GetXaxis()->SetTitle("Distance (cm)");
  fHistMatchDistance->GetYaxis()->SetTitle("Count");
  fHistMatchDistance->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistMatchDistance);

  fHistMatchedEnergy = new TH1F( Form("%s fHistMatchedEnergy", det.Data()), Form("%s Energy distribution of clusters with matching tracks", det.Data()), 200, 0, 200);
  fHistMatchedEnergy->GetXaxis()->SetTitle("Cluster Energy (GeV)");
  fHistMatchedEnergy->GetYaxis()->SetTitle("Number of clusters");
  fHistMatchedEnergy->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistMatchedEnergy);

  fHistUnMatchedEnergy = new TH1F( Form("%s fHistUnMatchedEnergy", det.Data()), Form("%s Energy distribution of clusters with no matching track", det.Data()), 200, 0, 200);
  fHistUnMatchedEnergy->GetXaxis()->SetTitle("Cluster Energy (GeV)");
  fHistUnMatchedEnergy->GetYaxis()->SetTitle("Number of clusters");
  fHistUnMatchedEnergy->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistUnMatchedEnergy);

}


AliHLTCaloHistoMatchedTracks::~AliHLTCaloHistoMatchedTracks()
{

  if(fHistArrayPtr)
    delete fHistArrayPtr;
  fHistArrayPtr = NULL;

  if(fHistMatchDistance){
    delete fHistMatchDistance;
    fHistMatchDistance = NULL;
  }

  if(fHistMatchedEnergy) {
    delete fHistMatchedEnergy;
    fHistMatchedEnergy = NULL;
  }

  if(fHistUnMatchedEnergy) {
    delete fHistUnMatchedEnergy;
    fHistUnMatchedEnergy = NULL;
  }

}


TObjArray* AliHLTCaloHistoMatchedTracks::GetHistograms()
{  
  // See header file for documentation
  return fHistArrayPtr;
}


// Int_t AliHLTCaloHistoMatchedTracks::DoEvent(AliHLTCaloClusterHeaderStruct* cHeader) {   
  
//   fClusterReader->SetMemory(cHeader);
  
//   AliHLTCaloClusterDataStruct* cluster;
//   while ( ( cluster = fClusterReader->NextCluster() ) ) {
    
//     if(cluster->fTracksMatched->GetSize()>0) {
//       fHistMatchedEnergy->Fill(cluster->fEnergy);
//     } else {
//       fHistUnMatchedEnergy->Fill(cluster->fEnergy);
//     }
//   }  
  
  
//   return 0;
// }
  

Int_t AliHLTCaloHistoMatchedTracks::FillHistograms(Int_t nClusters, TRefArray* clustersArray) {   
  
  for(int ic = 0; ic < nClusters; ic++) {
    
    AliESDCaloCluster * cluster = dynamic_cast<AliESDCaloCluster*>(clustersArray->At(ic));
    
    if(cluster->GetNTracksMatched() > 0) {
      fHistMatchedEnergy->Fill(cluster->E());
      fHistMatchDistance->Fill(cluster->GetEmcCpvDistance());
    
    } else {
      fHistUnMatchedEnergy->Fill(cluster->E());
    }
  
  }
  
  return 0;
}
