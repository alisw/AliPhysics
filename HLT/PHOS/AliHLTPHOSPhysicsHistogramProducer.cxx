
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
 * @file   AliHLTPHOSPhysicsHistogramProducer
 * @author Albin Gaignette
 * @date 
 * @brief  Histogram producer for PHOS HLT 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSPhysicsHistogramProducer.h"
#include "AliESDCaloCluster.h"
#include "TMath.h"
#include "TClonesArray.h"
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"

AliHLTPHOSPhysicsHistogramProducer::AliHLTPHOSPhysicsHistogramProducer() :
  //    AliHLTPHOSBase(),
    fHistNcls(0),
    fHistEnergy(0),
    fHistTotEnergy(0),
    fHistTwoClusterInvMass(0),
    fHistNcells(0),
    fHistNcellsPercentage(0),
    fHistCellsEnergy(0),
    fHistNclusterTotE(0),
    fHistNcellsSumCells(0),
    fHistArrayPtr(0)
{
  // See header file for documentation
  fHistArrayPtr = new TObjArray;

  fHistNcls = new TH1F("fHistNcls", "Number of clusters per event", 20, 0, 20);
  fHistNcls->GetXaxis()->SetTitle("Number of clusters per event");
  fHistNcls->GetYaxis()->SetTitle("# of counts");
  fHistArrayPtr ->AddLast(fHistNcls);
  
  fHistEnergy = new TH1F("fHistEnergy", "Energy in each cluster", 200, 0, 100);
  fHistEnergy->GetXaxis()->SetTitle("E_{T} GeV");
  fHistEnergy->GetYaxis()->SetTitle("# of counts");
  fHistEnergy->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistEnergy);

  fHistTotEnergy = new TH1F("fHistTotEnergy", "Total energy in each event", 200, 0, 200);
  fHistTotEnergy->GetXaxis()->SetTitle("E_{T} GeV");
  fHistTotEnergy->GetYaxis()->SetTitle("# of counts");
  fHistTotEnergy->SetMarkerStyle(21);
  fHistArrayPtr ->AddLast(fHistTotEnergy);

  fHistTwoClusterInvMass = new TH1F("fHistTwoClusterInvMass", "Invariant mass of two clusters", 200, 0, 1);
  fHistTwoClusterInvMass->GetXaxis()->SetTitle("m_{#gamma#gamma} GeV");
  fHistTwoClusterInvMass->GetYaxis()->SetTitle("Number of counts");
  fHistTwoClusterInvMass->SetMarkerStyle(21);
  fHistArrayPtr->AddLast(fHistTwoClusterInvMass);

  fHistNcells = new TH1F("fHistNcells", "Number of cells in each cluster", 200, 0, 100);
  fHistNcells->GetXaxis()->SetTitle("# of cells in the cluster");
  fHistNcells->GetYaxis()->SetTitle("# of counts");
  fHistNcells->SetMarkerStyle(21); 
  fHistArrayPtr ->AddLast(fHistNcells);

  fHistNcellsPercentage = new TH1F("fHistNcellsPercentage", "Percentage of cells in the cluster", 200, 0, 10);
  fHistNcellsPercentage->GetXaxis()->SetTitle("100* #frac{Total # of cells in the cluster}{Total # of cells in PHOS}");
  fHistNcellsPercentage->GetYaxis()->SetTitle("# of counts");
  fHistNcellsPercentage->SetMarkerStyle(21);
  fHistArrayPtr ->AddLast(fHistNcellsPercentage);

  fHistCellsEnergy = new TH2F("fHistCellsEnergy","# of cells in the cluster vs cluster Energy",100,0,60,100,1,60);
  fHistCellsEnergy->GetXaxis()->SetTitle("Cluster Energy [GeV]");
  fHistCellsEnergy->GetYaxis()->SetTitle("# of cells in the cluster");
  fHistArrayPtr ->AddLast(fHistCellsEnergy);

  fHistNclusterTotE = new TH2F("fHistNclusterTotE","# of clusters vs total Energy",100,0,200,100,1,15);
  fHistNclusterTotE->GetXaxis()->SetTitle("total Energy in the event [GeV]");
  fHistNclusterTotE->GetYaxis()->SetTitle("# of clusters in each event");
  fHistArrayPtr ->AddLast(fHistNclusterTotE);

  fHistNcellsSumCells = new TH2F("fHistNcellsSumCells","# of cells in PHOS vs sum cells in the clusters",100,0,300,100,17910,17930);
  fHistNcellsSumCells->GetXaxis()->SetTitle("# of cells in the clusters");
  fHistNcellsSumCells->GetYaxis()->SetTitle("Total # of cells in PHOS");
  fHistArrayPtr->AddLast(fHistNcellsSumCells);
}

AliHLTPHOSPhysicsHistogramProducer::~AliHLTPHOSPhysicsHistogramProducer()
{
  // See header file for documentation
  if(fHistNcls)
    {
      delete fHistNcls;
      fHistNcls = 0;
    }
  if(fHistEnergy)
    {
      delete fHistEnergy;
      fHistEnergy = 0;
    }
  if(fHistTotEnergy)
    {
      delete fHistTotEnergy;
      fHistTotEnergy = 0;
    }
  if(fHistTwoClusterInvMass)
    {
      delete fHistTwoClusterInvMass;
      fHistTwoClusterInvMass = 0;
    }
  if(fHistNcells)
    {
      delete fHistNcells;
      fHistNcells = 0;
    }
  if(fHistNcellsPercentage)
    {
      delete fHistNcellsPercentage;
      fHistNcellsPercentage = 0;
    }
  if(fHistCellsEnergy)
    {
      delete fHistCellsEnergy;
      fHistCellsEnergy = 0;
    }
  if(fHistNclusterTotE)
    {
      delete fHistNclusterTotE;
      fHistNclusterTotE = 0;
    }
  if(fHistNcellsSumCells)
    {
      delete fHistNcellsSumCells;
      fHistNcellsSumCells = 0;
    }
  if(fHistArrayPtr)
    {
      delete fHistArrayPtr;
      fHistArrayPtr = 0;
    }
}

TObjArray* AliHLTPHOSPhysicsHistogramProducer::GetHistograms()
{  
  // See header file for documentation
  return fHistArrayPtr;
}

Int_t AliHLTPHOSPhysicsHistogramProducer::AnalyseClusters(TClonesArray* clusters)
{   
  // See header file for documentation
  Int_t nPHOSModules = 3;
  Int_t totClustersAll=0;
  Int_t totClustersAllwithprob=0;
  Int_t totClustersFile=0; 
  Int_t totClustersFilewithprob=0; 
  Double_t theta=0.0; 
  Double_t m=0.0;
  Double_t TotalCells =0;
  Double_t TotalCellsPHOS= NXCOLUMNSMOD*NZROWSMOD * nPHOSModules;
  // 3584 crystals * 3 modules
  Double_t Nclusters =0;
  Double_t TotEnergy = 0.0;
     
  Int_t ncls = clusters->GetEntriesFast();

  fHistNcls->Fill(ncls);
  totClustersAll+=ncls;
  totClustersFile+=ncls;
  Double_t totE = 0;
   
  Nclusters=0;
  TotEnergy=0.0;
  TotalCells=0;
          
  Float_t** tabposition = new Float_t*[ncls];
  Double_t* tabenergy=new Double_t[ncls];
     
  for(Int_t icls = 0; icls < ncls; icls++)
    {
      AliESDCaloCluster* cluster = (AliESDCaloCluster*)clusters->At(icls);
	  	
      totClustersFilewithprob+=icls;
      totClustersAllwithprob+=icls;
      totE += cluster->E();
      fHistEnergy->Fill(cluster->E());

      tabposition[icls] = new Float_t[3];
      cluster->GetPosition(tabposition[icls]); 
      tabenergy[icls] = cluster->E(); 
  	  
      TotalCells += cluster -> GetNCells();
      fHistNcells->Fill(cluster -> GetNCells());
      fHistCellsEnergy->Fill(cluster -> GetNCells(),cluster->E());
      TotEnergy += cluster->E();
      Nclusters++;
    }
 
  fHistNcellsPercentage->Fill(100*(TotalCells/TotalCellsPHOS));
  fHistNclusterTotE->Fill(TotEnergy,Nclusters);
  fHistNcellsSumCells->Fill(TotalCells,TotalCellsPHOS);	
       
  if(totE > 0)
    fHistTotEnergy->Fill(totE);
  
  for(Int_t ipho = 0; ipho<(ncls-1) ; ipho++)
    { 
      for(Int_t jpho = ipho+1 ; jpho<ncls ; jpho++)
	{ 
	  // Calcul of the theta angle between two photons
	  theta =(2* asin(0.5*TMath::Sqrt((tabposition[ipho][0]-tabposition[jpho][0])*(tabposition[ipho][0]-tabposition[jpho][0]) +(tabposition[ipho][1]-tabposition[jpho][1])*(tabposition[ipho][1]-tabposition[jpho][1]))/460));
 	      
	  // Calcul of the mass m of the pion 
	  m =(TMath::Sqrt(2 * tabenergy[ipho]* tabenergy[jpho]*(1-TMath::Cos(theta))));
	  fHistTwoClusterInvMass->Fill(m);
	}
    }
    
  for(Int_t j=0 ; j<ncls;j++)
    {
      delete[] tabposition[j];
    }
	  
  delete[] tabposition;
  delete[] tabenergy;

  return 0;
}
  
 
 

