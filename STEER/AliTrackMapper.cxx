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

/*
$Log$
*/

////////////////////////////////////////////////////////////////////////
//
// AliTrackMapper.cxx
// description: 
//        create an AliTrackMap - a relation between track label and 
//        it's index in TreeH. Check all branches with hits.
//        Output is in the root file.
//
////////////////////////////////////////////////////////////////////////

#include <iostream.h>

#include "TTree.h"
#include "TROOT.h"
#include "TFile.h"
#include "TBenchmark.h"
#include "TStopwatch.h"

#include "AliDetector.h"
#include "AliTrackMapper.h"
#include "AliTrackMap.h"
#include "AliRun.h"
#include "AliHit.h"

ClassImp(AliTrackMapper)

////////////////////////////////////////////////////////////////////////
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

  TFile *fileHits=TFile::Open(fnHits);
  if (!fileHits->IsOpen()) {cerr<<"Can't open input file "<<fnHits<<"!\n"; return;}
  if (!(gAlice=(AliRun*)fileHits->Get("gAlice"))) {
    cerr<<"gAlice have not been found on galice.root !\n";
    return;
  }

  for (Int_t eventNr = firstEventNr; eventNr < firstEventNr+nEvents;
       eventNr++) {
    CreateMap(eventNr,fileMap);
  } // end loop over events

  delete gAlice;
  gAlice = 0;
  fileHits->Close();
  delete fileHits;
  fileMap->Close();
  delete fileMap;
  timer.Stop();
  if (fDEBUG > 0) timer.Print();
}

////////////////////////////////////////////////////////////////////////
Int_t  AliTrackMapper::CreateMap(Int_t eventNr, TFile* fileMap) {
//
// create an AliTrackMap for a given event
// correct gAlice must be already present in memory
//
  Int_t nGenPrimPlusSecParticles = gAlice->GetEvent(eventNr);
  if (fDEBUG > 1) cout<<"nGenPrimPlusSecParticles = "<<nGenPrimPlusSecParticles<<endl;
  if (nGenPrimPlusSecParticles < 1) {
    cerr<<"No primary particles found in event "<<eventNr<<endl;
    return -1;
  }

  TTree *treeK = gAlice->TreeK();
  if (!treeK) {
    cerr<<"Error: Event "<<eventNr<<", treeK not found."<<endl;
    return -1;
  }
  Int_t nAllParticles = static_cast<Int_t>(treeK->GetEntries());
  Int_t *trackMap = new Int_t[nAllParticles];
  for (Int_t i = 0; i<nAllParticles; i++) {trackMap[i] = kNoEntry;}

  TTree *treeH = gAlice->TreeH();
  if (!treeH) {
    cerr<<"Error: Event "<<eventNr<<", treeH not found."<<endl;
    return -1;
  }
  Int_t treeHEntries = static_cast<Int_t>(treeH->GetEntries());
  if (fDEBUG > 1) cout<<"treeHEntries "<<treeHEntries<<endl;


  TObjArray *modules = gAlice->Detectors();
  if (!modules) {
    cerr<<"TObjArray with modules not found."<<endl;
    return -1;
  }
  Int_t nModules = static_cast<Int_t>(modules->GetEntries());
  for (Int_t iModule = 0; iModule < nModules; iModule++) {
//  for (Int_t iModule = nModules-1; iModule >= 0; iModule--) {
    AliDetector * detector = dynamic_cast<AliDetector*>(modules->At(iModule));
    if (!detector) continue;
// process only detectors with shunt = 0
    if (detector->GetIshunt()) continue; 
    if (fDEBUG > 0) cerr<<"Processing detector "<<detector->GetName()<<endl;

    AliHit* hit;
    for (Int_t treeHIndex = 0; treeHIndex < treeHEntries; treeHIndex++) {
      detector->ResetHits();
      treeH->GetEvent(treeHIndex);
      hit=(AliHit*)detector->FirstHit(-1);
      Int_t lastLabel=-1, label;
      for( ; hit; hit=(AliHit*)detector->NextHit() ) {
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
    }
  }

  if (fDEBUG > 0) {
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
  
////////////////////////////////////////////////////////////////////////
AliTrackMap* AliTrackMapper::LoadTrackMap(Int_t eventNr, const char* fnMap, TFile* &fileMap) {
//
// read an AliTrackMap object for the given event eventNr from 
// the file fileMap
//
  fileMap=TFile::Open(fnMap);
  if (!fileMap->IsOpen()) {cerr<<"Can't open file "<<fnMap<<" with map!\n"; return 0;}
  char mapName[20];
  sprintf(mapName,"AliTrackMap_%5.5d",eventNr);
  AliTrackMap* trackMapObject = (AliTrackMap*)fileMap->Get(mapName);
  if (!trackMapObject) {
    cerr<<"Error: map named "<<mapName<<" not found."<<endl;
    return 0;
  }
  return trackMapObject;
}
////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////

