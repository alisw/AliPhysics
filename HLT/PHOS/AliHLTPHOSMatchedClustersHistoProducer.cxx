
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
 * @file   AliHLTPHOSMatchedclustershistoProducer
 * @author Albin Gaignette
 * @date 
 * @brief  Histogram producer for PHOS HLT 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSMatchedClustersHistoProducer.h"
//#include "AliESDCaloCluster.h"
#include "TMath.h"

#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"
#include "TObjArray.h"
//#include "TClonesArray.h"
//#include <iostream>
#include "TH1F.h"
//#include "TH2F.h"

AliHLTPHOSMatchedClustersHistoProducer::AliHLTPHOSMatchedClustersHistoProducer() :
  fClusterReader(NULL),
  fHistArrayPtr(0),
  fHistMatchQuality(0),
  fHistMatchedEnergy(0),
  fHistUnMatchedEnergy(0)
{
  // See header file for documentation
  fHistArrayPtr = new TObjArray;
  fClusterReader = new AliHLTCaloClusterReader();

  fHistMatchQuality = new TH1F("fHistMatchQuality", "Distance between cluster and track intersection with phos module", 50, 0, 50);
  fHistMatchQuality->GetXaxis()->SetTitle("Distance (cm)");
  fHistMatchQuality->GetYaxis()->SetTitle("Count");
  fHistMatchQuality->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistMatchQuality);

  fHistMatchedEnergy = new  TH1F("fHistMatchedEnergy", "Energy distribution of clusters, negative x is unmatched, positive x is matched", 400, -200, 200);
  fHistMatchedEnergy->GetXaxis()->SetTitle("Cluster Energy (GeV)");
  fHistMatchedEnergy->GetYaxis()->SetTitle("Number of clusters. Negative x direction is unmatched track, positive matched");
  fHistMatchedEnergy->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistMatchedEnergy);

  fHistUnMatchedEnergy = new  TH1F("fHistUnMatchedEnergy", "Energy distribution of clusters, negative x is unmatched, positive x is matched", 400, -200, 200);
  fHistUnMatchedEnergy->GetXaxis()->SetTitle("Cluster Energy (GeV)");
  fHistUnMatchedEnergy->GetYaxis()->SetTitle("Number of clusters. Negative x direction is unmatched track, positive matched");
  fHistUnMatchedEnergy->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistUnMatchedEnergy);
}


AliHLTPHOSMatchedClustersHistoProducer::~AliHLTPHOSMatchedClustersHistoProducer()
{

  if(fHistMatchQuality){
    delete fHistMatchQuality;
    fHistMatchQuality = 0;
  }

  if(fHistMatchedEnergy) {
    delete fHistMatchedEnergy;
    fHistMatchedEnergy = 0;
  }

  if(fHistUnMatchedEnergy) {
    delete fHistUnMatchedEnergy;
    fHistUnMatchedEnergy = 0;
  }

}


TObjArray* AliHLTPHOSMatchedClustersHistoProducer::GetHistograms()
{  
  // See header file for documentation

  return fHistArrayPtr;
}


Int_t AliHLTPHOSMatchedClustersHistoProducer::DoEvent(AliHLTCaloClusterHeaderStruct* cHeader) {   
  
  fClusterReader->SetMemory(cHeader);
  
  AliHLTCaloClusterDataStruct* cluster;
  while ( ( cluster = fClusterReader->NextCluster() ) ) {
    
    fHistMatchQuality->Fill(cluster->fMatchedTrackDistance);
    if(cluster->fTrackDistance > -999) {
      fHistMatchedEnergy->Fill(cluster->fEnergy);
    } else {
      fHistUnMatchedEnergy->Fill(cluster->fEnergy);
    }
  }  
  
  
  return 0;
}
  
 
 

