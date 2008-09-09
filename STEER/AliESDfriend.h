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
#include "AliESDVZEROfriend.h"

//_____________________________________________________________________________
class AliESDfriend : public TObject {
public:
  AliESDfriend();
  AliESDfriend(const AliESDfriend &);
  virtual ~AliESDfriend();

  Int_t GetNumberOfTracks() const {return fTracks.GetEntriesFast();}
  AliESDfriendTrack *GetTrack(Int_t i) const {
     return (AliESDfriendTrack *)fTracks.UncheckedAt(i);
  }
  void AddTrack(const AliESDfriendTrack *t) {
     new(fTracks[fTracks.GetEntriesFast()]) AliESDfriendTrack(*t);
  }

  void SetVZEROfriend(AliESDVZEROfriend * obj);
  AliESDVZEROfriend *GetVZEROfriend(){ return fESDVZEROfriend; }

protected:
  TClonesArray fTracks;    // ESD friend tracks
  AliESDVZEROfriend *fESDVZEROfriend; // VZERO object containing complete raw data

  ClassDef(AliESDfriend,2) // ESD friend
};

#endif


