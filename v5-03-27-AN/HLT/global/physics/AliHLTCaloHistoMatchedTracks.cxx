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
#include "TObjArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "TRefArray.h"
#include "TString.h"

AliHLTCaloHistoMatchedTracks::AliHLTCaloHistoMatchedTracks(TString det) :
  fHistDxy(NULL),
  fHistDz(NULL),
  fHistDxyDz(NULL),
  fHistMatchedEnergy(NULL),
  fHistUnMatchedEnergy(NULL)
{

  fHistMatchedEnergy = new TH1F( Form("%s_fHistMatchedEnergy", det.Data()), Form("%s Energy distribution of clusters with matching tracks", det.Data()), 5000, 0, 100);
  fHistMatchedEnergy->GetXaxis()->SetTitle("Cluster Energy (GeV)");
  fHistMatchedEnergy->GetYaxis()->SetTitle("Number of clusters");
  fHistMatchedEnergy->SetMarkerStyle(21);
  fHistArray->AddLast(fHistMatchedEnergy);

  fHistUnMatchedEnergy = new TH1F( Form("%s_fHistUnMatchedEnergy", det.Data()), Form("%s Energy distribution of clusters with no matching track", det.Data()), 5000, 0, 100);
  fHistUnMatchedEnergy->GetXaxis()->SetTitle("Cluster Energy (GeV)");
  fHistUnMatchedEnergy->GetYaxis()->SetTitle("Number of clusters");
  fHistUnMatchedEnergy->SetMarkerStyle(21);
  fHistArray->AddLast(fHistUnMatchedEnergy);


  fHistDxyDz = new TH2F( Form("%s_fHistdXYdZ", det.Data()), Form("%s dXY - dZ distribution of track - cluster residuals", det.Data()), 50, -50, 50, 50, -50, 50);
  fHistDxyDz->GetXaxis()->SetTitle("sqrt(dx^2 + dy^2)  (cm)");
  fHistDxyDz->GetYaxis()->SetTitle("dz (cm)");
  fHistArray->AddLast(fHistDxyDz);

  fHistDxy = new TH1F( Form("%s_fHistdXY", det.Data()), Form("%s #sqrt(dx^2 + dy^2)", det.Data()), 100, -50, 50);
  fHistDxy->GetXaxis()->SetTitle("sqrt(dx^2 + dy^2)  (cm)");
  fHistArray->AddLast(fHistDxy);
		       
  fHistDz = new TH1F( Form("%s_fHistdZ", det.Data()), Form("%s dZ", det.Data()),100, -50, 50);
  fHistDz->GetXaxis()->SetTitle("dZ (cm)");
  fHistArray->AddLast(fHistDz);

}


AliHLTCaloHistoMatchedTracks::~AliHLTCaloHistoMatchedTracks()
{

  if(fHistMatchedEnergy) 
    delete fHistMatchedEnergy;
  fHistMatchedEnergy = NULL;

  if(fHistUnMatchedEnergy) 
    delete fHistUnMatchedEnergy;
  fHistUnMatchedEnergy = NULL;

  if (fHistDxyDz) 
    delete fHistDxyDz;
  fHistDxyDz = NULL;

  if (fHistDxy) 
    delete fHistDxy;
  fHistDxy = NULL;

  if (fHistDz) 
    delete fHistDz;
  fHistDz = NULL;

}
  

Int_t AliHLTCaloHistoMatchedTracks::FillHistograms(Int_t nc, TRefArray * clusterArray) {
  //See header file for documentation
  for(int ic = 0; ic < nc; ic++) {
    AliESDCaloCluster * cluster = static_cast<AliESDCaloCluster*>(clusterArray->At(ic));
    FillMatchedTracks(cluster);
  }
  return 0;
}

Int_t AliHLTCaloHistoMatchedTracks::FillHistograms(Int_t nc, vector<AliHLTCaloClusterDataStruct*> &cVec) {
  for(int ic = 0; ic < nc; ic++) {
    AliHLTCaloClusterDataStruct * cluster = cVec.at(ic);
    FillMatchedTracks(cluster);
  }
  return 0;
}

template <class T>
Int_t AliHLTCaloHistoMatchedTracks::FillMatchedTracks(T* cluster){
  // HLTInfo("Filling track-matching histograms");

  if(cluster->GetNTracksMatched() > 0) {
    fHistMatchedEnergy->Fill(cluster->E());
    fHistDz->Fill(cluster->GetTrackDz());
    fHistDxy->Fill(cluster->GetTrackDx());
    fHistDxyDz->Fill(cluster->GetTrackDx(), cluster->GetTrackDz());
  } else {
    fHistUnMatchedEnergy->Fill(cluster->E());
  }
  
  return 0;
}

