#ifndef ALIEVE_MUONData_H
#define ALIEVE_MUONData_H

#include <Reve/Reve.h>

#include <TObject.h>

#include <vector>

class TTree;
class TString;

class AliRawReader;

namespace Alieve {

class MUONChamberData;

class MUONData : public TObject, public Reve::ReferenceCount
{

 protected:

  std::vector<MUONChamberData*>   fChambers;           // vector of 14 chambers

  static AliRawReader*            fgRawReader;         // raw reader

  Int_t fNTrackList;      // number of MC tracks which have hits
  Int_t fTrackList[256];  // list of MC tracks which have hits

 public:

  MUONData();
  virtual ~MUONData();

  MUONData(const MUONData&);
  MUONData& operator=(const MUONData&);

  void Reset();

  void LoadDigits(TTree* tree);
  void LoadRecPoints(TTree* tree);
  void LoadHits(TTree* tree);
  void LoadRaw(TString fileName);

  void CreateChamber(Int_t chamber);
  void CreateAllChambers();
  void DropAllChambers();
  void DeleteAllChambers();

  void  RegisterTrack(Int_t track);
  Int_t GetNTrackList() { return fNTrackList; }
  Int_t GetTrack(Int_t index);

  MUONChamberData* GetChamberData(Int_t chamber);

  ClassDef(MUONData,1);           // Manages MUON data for one event

};

}

#endif
