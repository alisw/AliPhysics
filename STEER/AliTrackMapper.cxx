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

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
// AliTrackMapper.cxx
// description: 
//        create an AliTrackMap - a relation between track label and 
//        it's index in TreeH. Check all branches with hits.
//        Output is in the root file.
//
////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include "TError.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TTree.h"

#include "AliDetector.h"
#include "AliHit.h"
#include "AliLoader.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliTrackMap.h"
#include "AliTrackMapper.h"

ClassImp(AliTrackMapper)

//_______________________________________________________________________
AliTrackMapper::AliTrackMapper(): 
  fDEBUG(0)
{
  //
  // Default ctor
  //
}


//_______________________________________________________________________
void AliTrackMapper::CreateMap(Int_t nEvents, Int_t firstEventNr,
                               const char* fnMap, const char* fnHits)
{
  //
  // method to create a track map for a given number of events starting 
  // from event firstEventNr. This method just opens and closes files and
  // loops over events, the creation of map is delegated to 
  // CreateMap(Int_t eventNr, TFile* fileMap) method
  //
  TStopwatch timer;
  timer.Start();
  
  TFile *fileMap=TFile::Open(fnMap,"new");
  if (!fileMap->IsOpen()) {cerr<<"Can't open output file "<<fnMap<<"!\n"; return;}

  AliRunLoader* rl = AliRunLoader::Open(fnHits);
  if (rl == 0x0) {cerr<<"Can't open input file "<<fnHits<<"!\n"; return;}
  if (rl->LoadgAlice())
    {
      ::Error("AliTrackMapper::CreateMap","Error occured while loading AliRun");
      return;
    }
    
  if (!(gAlice=rl->GetAliRun())) {
    cerr<<"gAlice have not been found on session !\n";
    return;
  }
  
  rl->LoadKinematics();
  
  for (Int_t eventNr = firstEventNr; eventNr < firstEventNr+nEvents;
       eventNr++) {
    CreateMap(eventNr,fileMap,rl);
  } // end loop over events
  
  delete rl;
  fileMap->Close();
  delete fileMap;
  timer.Stop();
  if (fDEBUG > 0) timer.Print();
}

//_______________________________________________________________________
Int_t  AliTrackMapper::CreateMap(Int_t eventNr, TFile* fileMap,AliRunLoader* rl) 
{
  //
  // create an AliTrackMap for a given event
  // correct gAlice must be already present in memory
  //

  rl->GetEvent(eventNr);

  TTree *treeK = rl->TreeK();
  if (!treeK) {
    cerr<<"Error: Event "<<eventNr<<", treeK not found."<<endl;
    return -1;
  }
  Int_t nAllParticles = static_cast<Int_t>(treeK->GetEntries());
  Int_t *trackMap = new Int_t[nAllParticles];
  for (Int_t i = 0; i<nAllParticles; i++) {trackMap[i] = kNoEntry;}


  TObjArray *modules = gAlice->Detectors();
  if (!modules) {
    cerr<<"TObjArray with modules not found."<<endl;
    return -1;
  }
  Int_t nModules = static_cast<Int_t>(modules->GetEntries());
  AliHit* hit;
  for (Int_t iModule = 0; iModule < nModules; iModule++) 
   {
    AliDetector * detector = dynamic_cast<AliDetector*>(modules->At(iModule));
    if (!detector) continue;
    AliLoader* loader = detector->GetLoader();
    if (loader == 0x0)
     {
       ::Warning("AliTrackMapper::CreateMap",
                 "Can not get loader from detector %s.",detector->GetName());
       continue;
     }
    Int_t retval = loader->LoadHits();
    if (retval) {
      ::Error("AliTrackMapper::CreateMap",
            "Event %d: error occured while loading hits for %s",
            eventNr,detector->GetName());
      return -1;
     }
    
    TTree *treeH = loader->TreeH();
    if (!treeH) {
      ::Error("AliTrackMapper::CreateMap","Event %d: Can not get TreeH for %s",
             eventNr,detector->GetName());
      return -1;
     }
    Int_t treeHEntries = static_cast<Int_t>(treeH->GetEntries());
    if (fDEBUG > 1) cout<<"treeHEntries "<<treeHEntries<<endl;
     
    detector->ResetHits();
    
    for (Int_t treeHIndex = 0; treeHIndex < treeHEntries; treeHIndex++)  
     { // process only detectors with shunt = 0
      treeH->GetEvent(treeHIndex);
      if (detector->GetIshunt()) continue; 

      hit=dynamic_cast<AliHit*>(detector->FirstHit(-1));
      Int_t lastLabel=-1, label;
      for( ; hit; hit=dynamic_cast<AliHit*>(detector->NextHit()) ) {
       label=hit->Track();       
       if (lastLabel != label) {
         if (label < 0 || label >= nAllParticles) {
           cerr<<"Error: label out of range. ";
           cerr<<"Event "<<eventNr<<" treeHIndex "<<treeHIndex<<" label = "<<label<<endl;
           return -2;
         }
         if (trackMap[label] >=0 && trackMap[label] != treeHIndex) {
           cerr<<"Error: different treeHIndex for label "<<label
              <<" indeces: "<<trackMap[label]<<" != "<<treeHIndex;
           cerr<<" event "<<eventNr<<" detector "<<detector->GetName()<<endl;
           return -3;
         }
         trackMap[label] = treeHIndex;
         if (fDEBUG > 2) cout<<"treeHIndex, label = "<<treeHIndex<<" "<<label<<endl;
         lastLabel = label;
       }
      }
    }//loop over hits in module
    loader->UnloadHits();
  }//loop over modules

  if (fDEBUG > 2) {
    for (Int_t i = 0; i < nAllParticles; i++) {
      cout<<eventNr<<"\t"<<i<<"\t"<<trackMap[i]<<endl;
    }
  }
  fileMap->cd();
  AliTrackMap* trackMapObject = new AliTrackMap(nAllParticles, trackMap);
  trackMapObject->SetEventNr(eventNr);
  trackMapObject->Write();
  delete trackMapObject;

  delete [] trackMap;
  return 0;
}
  
//_______________________________________________________________________
AliTrackMap* AliTrackMapper::LoadTrackMap(Int_t eventNr, const char* fnMap, TFile* &fileMap) {
  //
  // read an AliTrackMap object for the given event eventNr from 
  // the file fileMap
  //
  fileMap=TFile::Open(fnMap);
  if (!fileMap->IsOpen()) {cerr<<"Can't open file "<<fnMap<<" with map!\n"; return 0;}
  char mapName[20];
  sprintf(mapName,"AliTrackMap_%5.5d",eventNr);
  AliTrackMap* trackMapObject = dynamic_cast<AliTrackMap*>(fileMap->Get(mapName));
  if (!trackMapObject) {
    cerr<<"Error: map named "<<mapName<<" not found."<<endl;
    return 0;
  }
  return trackMapObject;
}

//_______________________________________________________________________
void AliTrackMapper::CheckTrackMap(Int_t eventNr, const char* fnMap) {
  //
  // 
  //
  TFile *fileMap;
  AliTrackMap* trackMapObject = LoadTrackMap(eventNr, fnMap, fileMap);
  if (!trackMapObject) return;
  
  trackMapObject->PrintValues();
  
  delete trackMapObject;
  fileMap->Close();
  delete fileMap;
}

