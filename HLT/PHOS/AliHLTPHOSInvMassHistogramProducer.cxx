
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
 * @file   AliHLTPHOSInvMassHistogramProducer
 * @author Albin Gaignette
 * @date 
 * @brief  Histogram producer for PHOS HLT 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSInvMassHistogramProducer.h"
//#include "AliESDCaloCluster.h"
#include "TMath.h"

#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"
#include "TObjArray.h"
//#include "TClonesArray.h"
//#include <iostream>
#include "TH1F.h"
//#include "TH2F.h"

AliHLTPHOSInvMassHistogramProducer::AliHLTPHOSInvMassHistogramProducer() :
  fClusterReader(NULL),
  fHistTwoClusterInvMass(0),
  fHistArrayPtr(0)
{
  // See header file for documentation
  fHistArrayPtr = new TObjArray;
  fClusterReader = new AliHLTCaloClusterReader();

  fHistTwoClusterInvMass = new TH1F("fHistTwoClusterInvMass", "Invariant mass of two clusters", 200, 0, 1);
  fHistTwoClusterInvMass->GetXaxis()->SetTitle("m_{#gamma#gamma} GeV");
  fHistTwoClusterInvMass->GetYaxis()->SetTitle("Number of counts");
  fHistTwoClusterInvMass->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistTwoClusterInvMass);
}

AliHLTPHOSInvMassHistogramProducer::~AliHLTPHOSInvMassHistogramProducer()
{
  if(fHistTwoClusterInvMass)
    {
      delete fHistTwoClusterInvMass;
      fHistTwoClusterInvMass = 0;
    }
}

TObjArray* AliHLTPHOSInvMassHistogramProducer::GetHistograms()
{  
  // See header file for documentation

  return fHistArrayPtr;
}

Int_t AliHLTPHOSInvMassHistogramProducer::DoEvent(AliHLTCaloClusterHeaderStruct* cHeader) {   
  
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
      
      fHistTwoClusterInvMass->Fill(m);
    }
  }
  
  return 0;
}
  
 
 

