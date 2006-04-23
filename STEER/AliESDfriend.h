#ifndef ALIESDFRIEND_H
#define ALIESDFRIEND_H

//-------------------------------------------------------------------------
//                     Class AliESDfriend
//               This class contains ESD additions
//       Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TClonesArray.h>

class AliESD;
class AliESDfriendTrack;

//_____________________________________________________________________________
class AliESDfriend : public TObject {
public:
  AliESDfriend();
  AliESDfriend(const AliESDfriend &);
  AliESDfriend(const AliESD &);
  virtual ~AliESDfriend();

  Int_t GetNumberOfTracks() const {return fTracks.GetEntriesFast();}
  AliESDfriendTrack *GetTrack(Int_t i) const {
    return (AliESDfriendTrack *)fTracks.UncheckedAt(i);
  }

protected:
  TClonesArray fTracks;    // ESD friend tracks
  ClassDef(AliESDfriend,1) // ESD friend
};

#endif


