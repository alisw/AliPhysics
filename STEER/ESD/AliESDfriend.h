#ifndef ALIESDFRIEND_H
#define ALIESDFRIEND_H

//-------------------------------------------------------------------------
//                     Class AliESDfriend
//               This class contains ESD additions
//       Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TClonesArray.h>

#include "AliESDfriendTrack.h"

class AliESDVZEROfriend;
class AliESDTZEROfriend;

//_____________________________________________________________________________
class AliESDfriend : public TObject {
public:
  AliESDfriend();
  AliESDfriend(const AliESDfriend &);
  AliESDfriend& operator=(const AliESDfriend& esd);  
  virtual ~AliESDfriend();

  Int_t GetNumberOfTracks() const {return fTracks.GetEntriesFast();}
  AliESDfriendTrack *GetTrack(Int_t i) const {
     return (AliESDfriendTrack *)fTracks.At(i);
  }
  Int_t GetEntriesInTracks() const {return fTracks.GetEntries();}
  void AddTrack(const AliESDfriendTrack *t) {
     new(fTracks[fTracks.GetEntriesFast()]) AliESDfriendTrack(*t);
  }

  void AddTrackAt(const AliESDfriendTrack *t, Int_t i) {
     new(fTracks[i]) AliESDfriendTrack(*t);
  }

  void SetVZEROfriend(AliESDVZEROfriend * obj);
  AliESDVZEROfriend *GetVZEROfriend(){ return fESDVZEROfriend; }
  void SetTZEROfriend(AliESDTZEROfriend * obj);
  AliESDTZEROfriend *GetTZEROfriend(){ return fESDTZEROfriend; }

  void Ls(){
	  return fTracks.ls();
  }

  // bit manipulation for filtering
  void SetSkipBit(Bool_t skip){SetBit(23,skip);}
  Bool_t TestSkipBit() {return TestBit(23);}

  //TPC cluster occupancy
  Int_t GetNclustersTPC(UInt_t sector) const {return (sector<72)?fNclustersTPC[sector]:0;}
  Int_t GetNclustersTPCused(UInt_t sector) const {return (sector<72)?fNclustersTPCused[sector]:0;}
  void SetNclustersTPC(UInt_t sector, Int_t occupancy) {if (sector<72) fNclustersTPC[sector]=occupancy;}
  void SetNclustersTPCused(UInt_t sector, Int_t occupancy) {if (sector<72) fNclustersTPCused[sector]=occupancy;}

protected:
  TClonesArray fTracks;    // ESD friend tracks
  AliESDVZEROfriend *fESDVZEROfriend; // VZERO object containing complete raw data
  AliESDTZEROfriend *fESDTZEROfriend; // TZERO calibration object
  
  Int_t fNclustersTPC[72]; //cluster occupancy per sector per sector
  Int_t fNclustersTPCused[72]; //number of clusters used in tracking per sector

  ClassDef(AliESDfriend,4) // ESD friend
};

#endif


