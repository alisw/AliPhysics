/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Francesco Blanco                                      *
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
 * @file   AliHLTEMCALClusterMonitor.cxx
 * @author Francesco Blanco
 * @date 
 * @brief  Online Monitoring Histogram maker for EMCAL  
 */
  

#include "AliHLTEMCALClusterMonitor.h"
#include "AliHLTCaloSharedMemoryInterfacev2.h"
#include "TMath.h"

ClassImp(AliHLTEMCALClusterMonitor);

AliHLTEMCALClusterMonitor::AliHLTEMCALClusterMonitor():
  fClusterReaderPtr(0),
  hList(0),
	hClusterEne(0),
	hClusterEneVsTime(0),
	hClusterCells(0),
	hClusterEneVsCells(0),
	hClusterEtaVsPhi(0)

{
  // See header file for documentation

  fClusterReaderPtr = new AliHLTCaloClusterReader();

  // Booking histograms
  hList = new TObjArray;

	hClusterEne = new TH1F("hClusterEne", "ClusterEnergy (GeV)", 200, 0, 100);
	hList->Add(hClusterEne);

	hClusterEneVsTime = new TH2F("hClusterEneVsTime", "ClusterEnergy vs Time", 200, 0, 100, 40, -0.5, 39.5);
	hList->Add(hClusterEneVsTime);
	
	hClusterCells = new TH1I("hClusterCells", "# of cells per cluster", 50, -0.5, 49.5);
	hList->Add(hClusterCells);

	hClusterEneVsCells = new TH2F("hClusterEneCells", "# of cells per cluster vs cluster energy", 200, 0, 100, 50, -0.5, 49.5);
	hList->Add(hClusterEneVsCells);

	hClusterEtaVsPhi = new TH2F("hClusterEtaVsPhi", "Cluster position in #eta#phi", 100, -0.7, 0.7, 100, 1.38, 3.15);
	hClusterEtaVsPhi->GetXaxis()->SetTitle("#eta");
	hClusterEtaVsPhi->GetYaxis()->SetTitle("#phi");
	hList->Add(hClusterEtaVsPhi);
	
	
}

AliHLTEMCALClusterMonitor::~AliHLTEMCALClusterMonitor() 
{
  //See header file for documentation
}

// Pointer to histograms objects
TObjArray* AliHLTEMCALClusterMonitor::GetHistograms()
{
  return hList;
}


Int_t AliHLTEMCALClusterMonitor::MakeHisto(AliHLTCaloClusterHeaderStruct *caloClusterHeaderPtr)
{

	// Cluster variables
	// Pointer to Cluster struture
	AliHLTCaloClusterDataStruct* caloClusterStructPtr = 0;

	float clusterEne, clusterTime, clusterEta, clusterPhi, clusterX, clusterY, clusterZ;
	int nCells;
	if (caloClusterHeaderPtr) {
	  
	  // stuff to handle clusters here
	  fClusterReaderPtr->SetMemory(caloClusterHeaderPtr);
	  
	  while((caloClusterStructPtr = fClusterReaderPtr->NextCluster()) != 0) {
			clusterX = caloClusterStructPtr->fGlobalPos[0];
			clusterY = caloClusterStructPtr->fGlobalPos[1];
			clusterZ = caloClusterStructPtr->fGlobalPos[2];

			nCells = caloClusterStructPtr->fNCells;
			clusterEne = caloClusterStructPtr->fEnergy;
			clusterTime = caloClusterStructPtr->fTOF;
	  	hClusterEne->Fill(clusterEne);
	  	hClusterEneVsTime->Fill(clusterEne, clusterTime);
	  	hClusterCells->Fill(nCells);
	  	hClusterEneVsCells->Fill(clusterEne, nCells);
			float r = TMath::Sqrt(clusterX*clusterX + clusterY*clusterY + clusterZ*clusterZ);
			clusterEta = 0.5*TMath::Log( (r+clusterZ)/(r-clusterZ) );
			clusterPhi = TMath::ATan2(clusterY, clusterX);
	  	hClusterEtaVsPhi->Fill(clusterEta, clusterPhi);

	    }
	  
	}

return 0; 
}
