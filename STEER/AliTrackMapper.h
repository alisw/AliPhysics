#ifndef ALITRACKMAPPER_H
#define ALITRACKMAPPER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
// description:
//   create a relation between track label and it's index
//   in TreeH. Check all branches with hits.
//   Output is in the root file.
//   See http://AliSoft.cern.ch/people/chudoba/classes/AliTrackMap.html
//                  
//  Author: Jiri Chudoba (CERN), 2002
//
////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

class TFile;

// --- AliRoot header files ---

class AliTrackMap;
class AliRunLoader;

class AliTrackMapper {

public:
  AliTrackMapper();
  virtual ~AliTrackMapper(){}
  void CreateMap(Int_t nEvents, Int_t firstEventNr, 
                 const char* fnMap = "trackMap.root",
                 const char* fnHits   ="rfio:galice.root");
  Int_t CreateMap(Int_t eventNr, TFile* fileMap,AliRunLoader* rl);
  void SetDebug(Int_t level) {fDEBUG = level;}
  void CheckTrackMap(Int_t eventNr, const char* fnMap = "trackMap.root");
  AliTrackMap* LoadTrackMap(Int_t eventNr, const char* fnMap, TFile* &fileMap);

    
private:

  Int_t fDEBUG;           // Debug flag
  
  ClassDef(AliTrackMapper,0)  // methods to create AliTrackMap
};

#endif // ALITRACKMAPPER_H





