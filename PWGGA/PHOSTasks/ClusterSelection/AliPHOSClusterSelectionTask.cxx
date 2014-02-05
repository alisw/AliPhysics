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
    fClusters(0x0),
    fSelectionMap(0x0)
{
  return;
}

AliPHOSClusterSelectionTask::~AliPHOSClusterSelectionTask()
{
  delete fClusters;
}
  
void AliPHOSClusterSelectionTask::UserCreateOutputObjects()
{
  return;
}

void AliPHOSClusterSelectionTask::UserExec(Option_t *option)
{
  AliVEvent* event = InputEvent();
  if( ! event ) 
    AliError("No Event");
  
  // initilise fClusters, array of PHOS Clusters
  if( 0x0 == fClusters ) fClusters = new TRefArray;
  event->GetPHOSClusters( fClusters );
  
  // Remove Clusters
  for(int index = 0; index < fClusters->GetEntriesFast(); ++index) { // TODO: check if array is indexed from 0
    AliVCluster* cluster = (AliVCluster*) fClusters->At(iClu); // TODO: check that fClusters is always compressed
    
    if( cluster->E() < kMinClusterEnergy // Low Energy Clusters
	|| cluster->GetDistanceToBadChannel() < kMinBCDistance // to close to Bad Channel
	|| cluster->GetNCells() < kMinNCells // to few cells
	|| cluster->GetM02() < kMinM02 
	) 
      fClusters->RemoveAt(index);
  }
  
  // Compact array after removel of clusters
  fClusters->Compact();
  
  // initialize fSelectionMap
  if( fSelectionMap )
    fSelectionMap->Clear();
  else {
    fSelectionMap = new TMap;
    fSelectionMap->SetOwnerValue();
  }
}

TRefArray* AliPHOSClusterSelectionTask::GetPHOSClusters() const
{
  if( ! fClusters )
    AliError("fCluster not initialized, do not run this function before ::UserExec");

  return fClusters;
}

TRefArray* AliPHOSClusterSelectionTask::GetPHOSClustersSelected(const AliPHOSClusterSelection* selection, bool useMap, bool addMap )
{
  // useMap - Use The Resulting Array of previous selection (stored in map
  // addMap - Add This Selection to the 'map' for use in future calls of this function.

  if( !fClusters  || !fSelectionMap )
    AliFatal("fCluster not initialized, do not run this function before ::UserExec");
  
  if( useMap ) {
    // Check if Selection is already done
    TRefArray* array = dynamic_selection<TRefArray*> ( fSelectionMap->GetValue(selection) );
    if( array )
      return array;
  }
    
  // if selected clusters not allready determined/in-map, determine and add to map:
  TRefArray* newArray = new TRefArray( * DeterminePHOSClustersSelected(selection) );
  if(addMap)
    fSelectionMap->Add(selection, newArray); // key, value
}

TRefArray* AliPHOSClusterSelectionTask::DeterminePHOSClustersSelected(const AliPHOSClusterSelection* selection)
{
  int nClu = fClusters->GetEntriesFast();
  
  // create/clear array
  static TRefArray* statRefArr = 0x0;
  if( statRefArr )
    statRefArr->Clear();
  else 
    statRefArr = new TRefArray(nClu);
  // array should now exist and be empty, 
  
  
  // fill array with selection:
  for(int iClu = 0; iClu < nClu; ++iClu) {
    AliVCluster* cluster = (AliVCluster*) fClusters->At(iClu);
    if( selection->IsSelected(cluster) )
      statRefArr->AddLast(cluster); // add at end of array
  }
  
  return statRefArr;
}


static AliPHOSClusterSelectionTask* AliPHOSClusterSelectionTask::GetTask(const char* name)
{
  // Get AliPHOSClusterSelectionTask from AliAnalysisManager
  
  AliAnalysisManager* analysisManager = dynamic_cast<AliAnalysisManager*>(AliAnalysisManager::GetAnalysisManager());
  if( !analysisManager ) 
    AliError("No AnalysisManager");
  AliAnalysisTask* task = analysisManager->GetTask(name);
  if( !task ) 
    AliError( Form("No task with name: %s", name) );

  AliPHOSClusterSelectionTask* sTask = dynamic_cast<AliAnalysisTask*>(task);
  if( !sTask) 
    AliError( Form("No AliPHOSClusterSelectionTask with name: %s", name) );
  
  return sTask;
}
