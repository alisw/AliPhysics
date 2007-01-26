#ifndef ALIEVE_MUONData_H
#define ALIEVE_MUONData_H

#include <Reve/Reve.h>

#include <TObject.h>

#include <vector>

class TTree;
class TString;
class TList;

class AliRawReader;
class AliMUONRawStreamTracker;
class AliMUONRawStreamTrigger;
class AliMUONDigit;
class AliMpSegFactory;
class AliMpDDLStore;
class AliMUONLocalTriggerBoard;
class AliMUONLocalStruct;
class AliMUONLocalStruct;

namespace Alieve {

class MUONChamberData;

class MUONData : public TObject, public Reve::ReferenceCount
{

 protected:

  std::vector<MUONChamberData*>   fChambers;           // vector of 14 chambers

  static AliRawReader*            fgRawReader;         // raw reader
  static AliMUONRawStreamTracker* fgRawStreamTracker;  // tracker raw streamer 
  static AliMUONRawStreamTrigger* fgRawStreamTrigger;  // trigger raw streamer
  static AliMpSegFactory*         fgSegFactory;        // segmentation mapping
  static AliMpDDLStore*           fgBusPatchManager;   // bus mapping

  Int_t    fNTracks;              // number of tracks
  Float_t* fTrackPoints;          // array with track points coordinates
  Int_t    fNPoints;              // total number of track points

  Int_t GetTrackerMapping(Int_t buspatchId, UShort_t manuId,
			  UChar_t channelId, AliMUONDigit* digit );

  Int_t GetTriggerMapping(AliMUONLocalTriggerBoard* localBoard, 
			  AliMUONLocalStruct* localStruct,
			  TList& digitList);

  void GetTriggerChamber(AliMUONLocalStruct* localStruct,
                         Int_t& xyPattern, Int_t& iChamber, Int_t& iCath, 
			 Int_t iCase);
 public:

  MUONData();
  virtual ~MUONData();

  MUONData(const MUONData&);
  MUONData& operator=(const MUONData&);

  void LoadDigits(TTree* tree);
  void LoadTracks(TTree* tree);
  void LoadRaw(TString fileName);
  void LoadRawTracker();
  void LoadRawTrigger();

  void CreateChamber(Int_t chamber);
  void CreateAllChambers();
  void DropAllChambers();
  void DeleteAllChambers();

  MUONChamberData* GetChamberData(Int_t chamber);

  Int_t    GetNTracks() { return fNTracks; };
  Int_t    GetNPoints() { return fNPoints; };
  Float_t* GetTrackPoints() { return fTrackPoints; };

  ClassDef(MUONData,1);           // Manages MUON data for one event

};

}

#endif
