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
#include "AliVfriendEvent.h"

#include "AliESDVZEROfriend.h"

class AliESDTZEROfriend;
class AliESDADfriend;

//_____________________________________________________________________________
class AliESDfriend : public AliVfriendEvent {
public:
  AliESDfriend();
  AliESDfriend(const AliESDfriend &);
  AliESDfriend& operator=(const AliESDfriend& esd);  
  virtual ~AliESDfriend();
  
  // This function will set the ownership
  // needed to read old ESDfriends
  void SetOwner(){
    fTracks.SetOwner();
    Int_t n=fTracks.GetEntriesFast();
    for(;n--;){
      AliESDfriendTrack *t=(AliESDfriendTrack *)fTracks.UncheckedAt(n);
      if(t)t->SetOwner();
    }
  }
  
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

  void SetVZEROfriend(const AliESDVZEROfriend * obj);
  AliESDVZEROfriend *GetVZEROfriend(){ return fESDVZEROfriend; }
  const AliESDVZEROfriend *GetVZEROfriendConst() const { return fESDVZEROfriend; }
  AliVVZEROfriend *GetVVZEROfriend(){ return fESDVZEROfriend; }
  Int_t GetESDVZEROfriend( AliESDVZEROfriend &v ) const {
    if( fESDVZEROfriend ){ v=*fESDVZEROfriend; return 0; }
    return -1;
  }

  void SetTZEROfriend(AliESDTZEROfriend * obj);
  AliESDTZEROfriend *GetTZEROfriend(){ return fESDTZEROfriend; }
  void SetADfriend(AliESDADfriend * obj);
  AliESDADfriend *GetADfriend(){ return fESDADfriend; }

  void Ls() const {
	  return fTracks.ls();
  }

  void Reset();
  // bit manipulation for filtering
  void SetSkipBit(Bool_t skip){SetBit(23,skip);}
  Bool_t TestSkipBit() const { return TestBit(23); }

  //TPC cluster occupancy
  Int_t GetNclustersTPC(UInt_t sector) const {return (sector<72)?fNclustersTPC[sector]:0;}
  Int_t GetNclustersTPCused(UInt_t sector) const {return (sector<72)?fNclustersTPCused[sector]:0;}
  void SetNclustersTPC(UInt_t sector, Int_t occupancy) {if (sector<72) fNclustersTPC[sector]=occupancy;}
  void SetNclustersTPCused(UInt_t sector, Int_t occupancy) {if (sector<72) fNclustersTPCused[sector]=occupancy;}

protected:
  TClonesArray fTracks;    // ESD friend tracks
  AliESDVZEROfriend *fESDVZEROfriend; // VZERO object containing complete raw data
  AliESDTZEROfriend *fESDTZEROfriend; // TZERO calibration object
  AliESDADfriend *fESDADfriend; // AD object containing complete raw data
  
  Int_t fNclustersTPC[72]; //cluster occupancy per sector per sector
  Int_t fNclustersTPCused[72]; //number of clusters used in tracking per sector

  ClassDef(AliESDfriend,5) // ESD friend
};

#endif


