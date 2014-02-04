/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "TBits.h"
#include "TRefArray.h"

#include "AliVCluster.h"

// Analysis task to fill histograms with PHOS ESD or AOD clusters and cells
// Authors : Henrik Qvigstad
// Date    : 28.05.2011
/* $Id$ */

#include "AliPHOSClusterSelectionTask.h"
ClassImp(AliPHOSClusterSelectionTask);

AliPHOSClusterSelectionTask::AliPHOSClusterSelectionTask(const char* name = "AliPHOSClusterSelectionTask")
  : AliAnalysisTaskSE(name), 
    fClusters(0x0)
{
  return;
}

AliPHOSClusterSelectionTask::~AliPHOSClusterSelectionTask()
{
  delete fClusters;
}
  
void AliPHOSClusterSelectionTask::UserCreateOutputObjects()
{
  //TODO: implement
  return;
}

void AliPHOSClusterSelectionTask::UserExec(Option_t *option)
{
  AliVEvent* event = InputEvent();
  
  if( 0x0 == fClusters ) fClusters = new TRefArray;
  // Get PHOS Clusters
  event->GetPHOSClusters( fClusters );
  
  // Remove Clusters
  for(int index = 0; index < fClusters->GetEntriesFast(); ++index) {
    AliVCluster* cluster = (AliVCluster*) fClusters->At(iClu);
    
    if( cluster->E() < kMinClusterEnergy // Low Energy Clusters
	|| cluster->GetDistanceToBadChannel() < kMinBCDistance // to close to Bad Channel
	|| cluster->GetNCells() < kMinNCells // to few cells
	|| cluster->GetM02() < kMinM02 
       ) 
      fClusters->RemoveAt(index);
  }
  
  // Compact array after removel of clusters
  fClusters->Compact();
}

TRefArray* AliPHOSClusterSelectionTask::GetPHOSClusters() const
{
  return fClusters;
}

TBits* AliPHOSClusterSelectionTask::GetPHOSClustersSelected(const AliPHOSClusterSelection* selection)
{
  int nClu = fClusters->GetEntriesFast();
  
  static TBits* sBits = 0x0;
  if( 0x0 == sBits ) sBits = new TBits(nClu);

  if( sBits->GetNbits != nClu ){
    delete sBits;
    sBits = new TBits( nClu );
  }

  
  // Determine selection:
  for(int iClu = 0; iClu < nClu; ++iClu) {
    AliVCluster* cluster = (AliVCluster*) fClusters->At(iClu);
    bool sel = selection->IsSelected(cluster);
    sBits->SetBitNumber(iClu, sel);
  }
  
  return sBits;


  // TMap initialization
  if( 0x0 == fSelectionMap ) {
    fSelectionMap = new TMap;
    fSelectionMap->SetOwnerValue(true);
  }

  // Fetch Selection Bits
  TObject* value = fSelectionMap->GetValue(selection);
}


static AliPHOSClusterSelectionTask* AliPHOSClusterSelectionTask::GetTask(const char* name)
{
  // Get AliPHOSClusterSelectionTask from AliAnalysisManager
  
  AliAnalysisManager* analysisManager = dynamic_cast<AliAnalysisManager*>(AliAnalysisManager::GetAnalysisManager());
  AliAnalysisTask* task = analysisManager->GetTask(name);
  if( !task ) AliWarning( Form("No task with name: %s", name) );

  AliPHOSClusterSelectionTask* stask = dynamic_cast<AliAnalysisTask*>(task);
  if( !stask) AliWarning( Form("No AliPHOSClusterSelectionTask with name: %s", name) );
  
  return stask;
}
