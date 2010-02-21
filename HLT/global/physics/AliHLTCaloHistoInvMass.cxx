//-*- Mode: C++ -*-
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Svein Lindal                                          *
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
 * @file   AliHLTCaloHistoInvMass
 * @author Svein Lindal <slindal@fys.uio.no>
 * @date 
 * @brief  Produces plots of invariant mass of two clusters
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTCaloHistoInvMass.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "TObjArray.h"
#include "AliESDEvent.h"
#include "TRefArray.h"
#include "TH1F.h"
#include "TString.h"
#include "AliESDCaloCluster.h"

AliHLTCaloHistoInvMass::AliHLTCaloHistoInvMass(TString det) :
  fHistTwoClusterInvMass(NULL)
{
  // See header file for documentation
  fHistTwoClusterInvMass = new TH1F(Form("%s fHistTwoClusterInvMass", det.Data()), Form("%s Invariant mass of two clusters PHOS", det.Data()), 200, 0, 1);
  fHistTwoClusterInvMass->GetXaxis()->SetTitle("m_{#gamma#gamma} GeV");
  fHistTwoClusterInvMass->GetYaxis()->SetTitle("Number of counts");
  fHistTwoClusterInvMass->SetMarkerStyle(21);
  fHistArray->AddLast(fHistTwoClusterInvMass);
}

AliHLTCaloHistoInvMass::~AliHLTCaloHistoInvMass()
{
  if(fHistTwoClusterInvMass)
    delete fHistTwoClusterInvMass;
  fHistTwoClusterInvMass = NULL;
}


Int_t AliHLTCaloHistoInvMass::FillHistograms(Int_t nc, vector<AliHLTCaloClusterDataStruct*> &cVec) {
  //See header file for documentation
  
  Float_t cPos[nc][3];
  Float_t cEnergy[nc];

  for(int ic = 0; ic < nc; ic++) {
    AliHLTCaloClusterDataStruct * cluster = cVec.at(ic);
    cluster->GetPosition(cPos[ic]);
    cEnergy[ic] = cluster->E();
  }

  for(Int_t ipho = 0; ipho<(nc-1); ipho++) { 
    for(Int_t jpho = ipho+1; jpho<nc; jpho++) { 
      
      // Calculate the theta angle between two photons
      Double_t theta = (2* asin(0.5*TMath::Sqrt((cPos[ipho][0]-cPos[jpho][0])*(cPos[ipho][0]-cPos[jpho][0]) +(cPos[ipho][1]-cPos[jpho][1])*(cPos[ipho][1]-cPos[jpho][1]))/460));
      
      // Calculate the mass m of the pion candidate
      Double_t m =(TMath::Sqrt(2 * cEnergy[ipho]* cEnergy[jpho]*(1-TMath::Cos(theta))));
      
      fHistTwoClusterInvMass->Fill(m);
    }
  }

  return 0;
}

Int_t AliHLTCaloHistoInvMass::FillHistograms(Int_t nc, TRefArray * clusterArray) {
  //See header file for documentation
  
  Float_t cPos[nc][3];
  Float_t cEnergy[nc];

  for(int ic = 0; ic < nc; ic++) {
    AliESDCaloCluster * cluster = static_cast<AliESDCaloCluster*>(clusterArray->At(ic));
    cluster->GetPosition(cPos[ic]);
    cEnergy[ic] = cluster->E();
  }

  for(Int_t ipho = 0; ipho<(nc-1); ipho++) { 
    for(Int_t jpho = ipho+1; jpho<nc; jpho++) { 
      
      // Calculate the theta angle between two photons
      Double_t theta = (2* asin(0.5*TMath::Sqrt((cPos[ipho][0]-cPos[jpho][0])*(cPos[ipho][0]-cPos[jpho][0]) +(cPos[ipho][1]-cPos[jpho][1])*(cPos[ipho][1]-cPos[jpho][1]))/460));
      
      // Calculate the mass m of the pion candidate
      Double_t m =(TMath::Sqrt(2 * cEnergy[ipho]* cEnergy[jpho]*(1-TMath::Cos(theta))));
      
      fHistTwoClusterInvMass->Fill(m);
    }
  }

  return 0;
}


// Int_t AliHLTCaloHistoInvMass::FillHistograms(Int_t nc, TRefArray * clustersArray) {
//   //See header file for documentation
  
//   Float_t cPos[nc][3];
//   Float_t cEnergy[nc];

//   for(int ic = 0; ic < nc; ic++) {
//     AliESDCaloCluster * cluster = static_cast<AliESDCaloCluster*>(clustersArray->At(ic));
//     cluster->GetPosition(cPos[ic]);
//     cEnergy[ic] = cluster->E();
//   }

//   for(Int_t ipho = 0; ipho<(nc-1); ipho++) { 
//     for(Int_t jpho = ipho+1; jpho<nc; jpho++) { 
//       // Calculate the theta angle between two photons
//       Double_t theta = (2* asin(0.5*TMath::Sqrt((cPos[ipho][0]-cPos[jpho][0])*(cPos[ipho][0]-cPos[jpho][0]) +(cPos[ipho][1]-cPos[jpho][1])*(cPos[ipho][1]-cPos[jpho][1]))/460));
      
//       // Calculate the mass m of the pion candidate
//       Double_t m =(TMath::Sqrt(2 * cEnergy[ipho]* cEnergy[jpho]*(1-TMath::Cos(theta))));
      
//       fHistTwoClusterInvMass->Fill(m);
//     }
//   }

//   return 0;
// }



// Int_t AliHLTCaloHistoInvMass::DoEvent(AliHLTCaloClusterHeaderStruct* cHeader) {   
  
//   fClusterReader->SetMemory(cHeader);
  
//   int ncls = cHeader->fNClusters;
//   Float_t* cPos[ncls];
//   Float_t cEnergy[ncls];
  
//   AliHLTCaloClusterDataStruct* cluster;
//   Int_t icls = 0;
//   while ( ( cluster = fClusterReader->NextCluster() ) ) {
    
//     cPos[icls] = cluster->fGlobalPos;
//     cEnergy[icls] = cluster->fEnergy; 
    
//     icls++;
//   }  
  
//   for(Int_t ipho = 0; ipho<(ncls-1); ipho++) { 
//     for(Int_t jpho = ipho+1 ; jpho<ncls ; jpho++) { 
//       // Calcul of the theta angle between two photons
//       Double_t theta = (2* asin(0.5*TMath::Sqrt((cPos[ipho][0]-cPos[jpho][0])*(cPos[ipho][0]-cPos[jpho][0]) +(cPos[ipho][1]-cPos[jpho][1])*(cPos[ipho][1]-cPos[jpho][1]))/460));
      
//       // Calcul of the mass m of the pion 
//       Double_t m =(TMath::Sqrt(2 * cEnergy[ipho]* cEnergy[jpho]*(1-TMath::Cos(theta))));
      

//       //BALLE
//       // fHistTwoClusterInvMass->Fill(m);
//     }
//   }
  
//   return 0;
// }
  
 
 

