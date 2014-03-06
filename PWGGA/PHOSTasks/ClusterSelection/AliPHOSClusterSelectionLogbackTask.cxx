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
#include "TObject.h"
#include "TMap.h"

#include "AliVCluster.h"
#include "AliPHOSGeometry.h"

#include "AliPHOSClusterSelection.h"
#include "AliPHOSEventSelection.h"
#include "AliPHOSLogbackCluster.h"

// Analysis task to fill histograms with PHOS ESD or AOD clusters and cells
// Authors : Henrik Qvigstad
// Date    : 28.05.2011
/* $Id$ */

#include "AliPHOSClusterSelectionLogbackTask.h"
ClassImp(AliPHOSClusterSelectionLogbackTask);


AliPHOSClusterSelectionLogbackTask::AliPHOSClusterSelectionLogbackTask(const char* name)
  : AliPHOSClusterSelectionTask(name),
    fMapOfEventLists(0x0)
{
  fMapOfEventLists = new TMap;
  fMapOfEventLists->SetOwnerValue();
}

AliPHOSClusterSelectionLogbackTask::~AliPHOSClusterSelectionLogbackTask()
{
  delete fMapOfEventLists;
}
  
void AliPHOSClusterSelectionLogbackTask::UserCreateOutputObjects()
{
  return;
}

void AliPHOSClusterSelectionLogbackTask::UserExec(Option_t *option)
{
  AliPHOSClusterSelectionTask::UserExec(option);

  // initialize fMapOfEventLists
  if( !fMapOfEventLists ) {
    fMapOfEventLists = new TMap;
    fMapOfEventLists->SetOwnerValue();
  }
}

TObjArray* AliPHOSClusterSelectionLogbackTask::GetPHOSClustersLogback(const AliPHOSEventSelection* eventSelection, UInt_t eventLogbackIndex, const AliPHOSClusterSelection* clusterSelection) const
{
  // Gives an array of clusters for an given event and cluster selection, 
  // provided that such an selection has been logged.
  // eventLogbackIndex - How many logged events to go back in time, must be an positive integer 0,1,2...
  // , an clusterSelection is optional.

  TObject* objEventList = fMapOfEventLists->GetValue(eventSelection);
  TList* eventList = dynamic_cast<TList*> ( objEventList );
  if(!eventList)
    return 0x0;
  
  // eventList is a list of eventArray
  TObject* objEventArray = eventList->At(eventLogbackIndex);
  TObjArray* eventArray = dynamic_cast<TObjArray*> ( objEventArray );
  if( !eventArray )
    return 0x0;
  
  // fMapOfEventLists[eventSelection][eventLogbackIndex][0]
  TObject* objEventClusters = eventArray->At(kCluArray);
  TObjArray* eventClusters = dynamic_cast<TObjArray*> ( objEventClusters );
  if( !eventClusters )
    AliFatal(Form("eventArray should always contain and TObjArray in index %i", kCluArray));
  
  // fMapOfEventLists[eventSelection][eventLogbackIndex][1]
  TObject* objSelectionMap = eventArray->At(kSelMap);
  TMap* selectionMap = dynamic_cast<TMap*> ( objSelectionMap );
  if( !selectionMap )
    AliFatal(Form("eventArray should always contain and TMap in index %i", kSelMap));

  
  // For the given eventSelection and eventBacklog, we should now have 
  // eventClusters and selectionMap !!!
  // which TObjArray we return to the user depends on if an clusterSelection is specified
  if(clusterSelection) {
    TObject* objSelectionArray = selectionMap->GetValue(clusterSelection);
    TObjArray* selectionArray = dynamic_cast<TObjArray*> ( objSelectionArray );
    return selectionArray;
  } 
  else
    return eventClusters;
}


void AliPHOSClusterSelectionLogbackTask::LogEvent(const AliPHOSEventSelection* eventSelection, UInt_t nEventsToLog)
{
  if( nEventsToLog < 1) {
    AliError("nEventsToLog needs to be >0, or logging does not make sense");
    return;
  }

  // Make a copy of the cluster array
  TObjArray* newCluArray = new TObjArray(fClusters->GetEntriesFast());
  newCluArray->SetOwner();
  for(int iclu=0; iclu < fClusters->GetSize(); ++iclu) {
    AliPHOSLogbackCluster* clu = new AliPHOSLogbackCluster((AliVCluster*)fClusters->At(iclu));
    newCluArray->Add(clu);
  }
  
  // Make a copy of the selection map
  // and referance to the new array of clusters
  TMap* newMap = new TMap;
  newMap->SetOwnerValue();
  TMapIter* iter = (TMapIter*) fSelectionMap->MakeIterator();
  TObject* key = 0x0;
  while((key = iter->Next())){ // Loop over Cluster Selections 
    TRefArray* oldSelArray = dynamic_cast<TRefArray*> (fSelectionMap->GetValue(key));
    TObjArray* newSelArray = new TObjArray(oldSelArray->GetEntriesFast());
    newSelArray->SetOwner(false);
    // Referance the new selection array to the new clusters
    for(int selInd=0; selInd<oldSelArray->GetEntriesFast(); ++selInd){
      bool matched = false; // false untill mached
      for(int cluInd=0; cluInd<fClusters->GetEntriesFast(); ++cluInd){
	TObject* oldSelCluster = oldSelArray->At(selInd);
	TObject* oldCluster = fClusters->At(cluInd);
	if(oldSelCluster == oldCluster) { // if old selection matches old cluster,
	  // then new selection should match cluster in the same index as old cluster
	  // at the same old selection index
	  if( newSelArray->At(cluInd) )
	    AliError("should be empty!");
	  newSelArray->AddAt(newCluArray->At(cluInd), selInd);

	  matched = true;
	  break;
	}
      }
      if( ! matched )
	AliError("Should have found a match!");
    }
    
    if( newSelArray->GetEntriesFast() != newSelArray->GetEntries()
	|| newSelArray->GetEntriesFast() != oldSelArray->GetEntriesFast() )
      AliError("Entries should match!");
    
    newMap->Add(key, newSelArray);
  }// endo of creation of new selMap.

  //make an eventArray to hold clusters and maps
  TObjArray* eventArray = new TObjArray(kSize);
  eventArray->SetOwner();
  eventArray->AddAt(newCluArray, kCluArray);
  eventArray->AddAt(newMap, kSelMap);
 
  // at to list of events
  TObject* objEventList = fMapOfEventLists->GetValue(eventSelection);
  TList* eventList = dynamic_cast<TList*> ( objEventList );
  if(!eventList) {
    eventList = new TList;
    eventList->SetOwner();
  }
  eventList->AddFirst(eventArray);
  
  
  // remove old events
  while( eventList->At(nEventsToLog) && nEventsToLog )
    eventList->RemoveLast();
}



AliPHOSClusterSelectionLogbackTask* GetTask(const char* name)
{
  AliPHOSClusterSelectionLogbackTask* task = dynamic_cast<AliPHOSClusterSelectionLogbackTask*>( AliPHOSClusterSelectionTask::GetTask(name) );
  if( !task )
    Printf( Form("No AliPHOSClusterSelectionLogbackTask with name: %s", name) );

  return task;
}
