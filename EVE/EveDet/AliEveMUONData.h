// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel & Bogdan Vulpescu: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveMUONData_H
#define AliEveMUONData_H

#include <TEveUtil.h>

#include <TObject.h>

#include <vector>

class TTree;
class TString;

class AliRawReader;

class AliEveMUONChamberData;

class AliEveMUONData : public TObject, public TEveRefCnt
{
public:

  AliEveMUONData();
  virtual ~AliEveMUONData();

  AliEveMUONData(const AliEveMUONData&);
  AliEveMUONData& operator=(const AliEveMUONData&);

  void Reset();

  void LoadDigits(TTree* tree);
  void LoadRecPoints(TTree* tree);
  void LoadRecPointsFromESD(const Char_t *fileName);
  void LoadHits(TTree* tree);
  void LoadRaw(TString fileName);

  void CreateChamber(Int_t chamber);
  void CreateAllChambers();
  void DropAllChambers();
  void DeleteAllChambers();

  void  RegisterTrack(Int_t track);
  Int_t GetNTrackList() const { return fNTrackList; }
  Int_t GetTrack(Int_t index) const;

  AliEveMUONChamberData* GetChamberData(Int_t chamber);

protected:

  std::vector<AliEveMUONChamberData*>   fChambers;   // vector of 14 chambers

  static AliRawReader                  *fgRawReader; // raw reader

  Int_t fNTrackList;      // number of MC tracks which have hits
  Int_t fTrackList[256];  // list of MC tracks which have hits

  ClassDef(AliEveMUONData, 0);           // Manages MUON data for one event

};

#endif
