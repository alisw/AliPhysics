#ifndef ALIEVE_MUONTracks_H
#define ALIEVE_MUONTracks_H

#include <TObject.h>

class AliMUONTrack;

namespace Reve {
  class Track;
}

namespace Alieve {

class MUONTracks: public TObject
{

  MUONTracks(const MUONTracks&);            // Not implemented
  MUONTracks& operator=(const MUONTracks&); // Not implemented

 public:

  MUONTracks();
  virtual ~MUONTracks();

  void MakeTrack(Int_t label, Reve::Track *rtrack, AliMUONTrack *mtrack);

  ClassDef(MUONTracks, 1);    // Produce Reve:Track from AliMUONTrack

};

}

#endif

