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
#include "TObjArray.h"

#include "AliVCluster.h"

// Analysis task to fill histograms with PHOS ESD or AOD clusters and cells
// Authors : Henrik Qvigstad
// Date    : 28.05.2011
/* $Id$ */

#include "AliPHOSClusterSelectionLogbackTask.h"
ClassImp(AliPHOSClusterSelectionLogbackTask);


AliPHOSClusterSelectionLogbackTask::AliPHOSClusterSelectionLogbackTask(const char* name = "AliPHOSClusterSelectionLogbackTask")
  : AliPHOSClusterSelectionTask(name),
    fClusterArrayList(0x0),
    fSelectionArrayMapLists(0x0)
{
  
}

AliPHOSClusterSelectionLogbackTask::~AliPHOSClusterSelectionLogbackTask()
{
  delete fClusterArrayList;
  delete fSelectionArrayMap;
}
  
void AliPHOSClusterSelectionLogbackTask::UserCreateOutputObjects()
{
  return;
}
void AliPHOSClusterSelectionLogbackTask::UserExec(Option_t *option)
{
  AliPHOSClusterSelectionTask::UserExec(option);
  AliVEvent* event = InputEvent();

  // initialize fClusterArrayLists
  if( !fClusterArrayLists ) {
    fClusterArrayLists = new TMap;
    fClusterArrayLists->SetOwnerValue();
  }
  // initialize fSelectionArrayMapLists
  if( !fSelectionArrayMapLists ) {
    fSelectionArrayMapLists = new TMap;
    fSelectionArrayMapLists->SetOwnerValue();
  }
}

TObjArray* AliPHOSClusterSelectionLogbackTask::GetPHOSClustersLogback(const AliPHOSEventSelection* eventSelection, UInt_t eventBacklogIndex, const AliPHOSClusterSelection* clusterSelection) const
{
  // Gives an array of clusters for an given event and cluster selection, 
  // provided that such an selection has been logged.
  // eventBacklogIndex - How many logged events to go back in time, must be an positive integer 0,1,2... 
  
  if( clusterSelection ) {
    TObject* objEventList = fSelectionArrayMapLists->GetValue(eventSelection);
    TList* eventList = dynamic_cast<TList*> ( objEventList );
    if(!eventList)
      return 0x0;

    TObject* objSelMap = eventList->At(eventBacklogIndex);
    TMap* selMap = dynamic_cast<TMap*> ( objSelMap );
    if(!semMap)
      return 0x0;
    
    TObject* objCluArray = selMap->GetValue(clusterSelection);
    TObjArray* cluArray = dynamic_cast<TObjArray*> ( objCluArray );
    return cluArray;
  }
  else { // i.e. !clusterSelection
    TObject* objEventList = fClusterArrayLists->GetValue(eventSelection);
    TList* eventList = dynamic_cast<TList*> ( objEventList );
    if(!eventList)
      return 0x0;
    
    TObject* objCluArray = eventList->At(eventBacklogIndex);
    TObjArray* cluArray = dynamic_cast<TObjArray*> ( objCluArray );
    return cluArray;
  }
}


void AliPHOSClusterSelectionLogbackTask::LogEvent(const AliPHOSEventSelection* eventSelection, int nEventsToLog)
{
  // Lists:
  TObject* objCluEventList = fClusterArrayLists->GetValue(eventSelection);
  TList* cluEventList = dynamic_cast<TList*> ( objCluEventList );
  TObject* objSelEventList = fSelectionArrayMapLists->GetValue(eventSelection);
  Tlist* selEventList = dynamic_cast<TList*> ( objSelEventList );

  if( !cluEventList || !selEventList) { // if either don't exist
    if( ! (!cluEventList && !selEventList) )
      AliFatal("neither or both should exist!");
    // else, initilize
    cluEventList = new TList;
    cluEventList->SetOwner();
    fClusterArrayLists->Add(eventSelection, cluEventList); // key, value
    
    selEventList = new TList;
    selEventList->SetOwner();
    fSelectionArrayMapLists->Add(eventSelection, selEventList); // key, value
  }

  // Log Clusters
  for(int index = 0; index < fClusters->GetEntriesFast; ++index) {
    //TODO
  }

  // Log Selection Maps
  // TODO
}



AliPHOSClusterSelectionLogbackTask* GetTask(const char* name)
{
  AliPHOSClusterSelectionLogbackTask* task = dynamic_cast<AliPHOSClusterSelectionLogbackTask*>( AliPHOSClusterSelectionTask::GetTask(name) );
  if( !task )
    AliError( Form("No AliPHOSClusterSelectionLogbackTask with name: %s"), name );

  return task;
}
