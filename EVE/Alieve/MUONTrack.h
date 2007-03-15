#ifndef ALIEVE_MUONTrack_H
#define ALIEVE_MUONTrack_H

#include <Reve/Track.h>

class AliMUONTrack;
class AliMagF;

namespace Reve {

  class TrackRnrStyle;
  class RecTrack;

}

namespace Alieve {

class MUONTrack: public Reve::Track
{

  MUONTrack(const MUONTrack&);            // Not implemented
  MUONTrack& operator=(const MUONTrack&); // Not implemented

 public:

  MUONTrack(Reve::RecTrack* t, Reve::TrackRnrStyle* rs);
  virtual ~MUONTrack();

  void  MakeTrack(AliMUONTrack *mtrack, AliMagF *fmap);
  void  GetField(Double_t *position, Double_t *field);
  void  Propagate(Float_t *xr, Float_t *yr, Float_t *zr, Int_t i1, Int_t i2);
  void  OneStepRungekutta(Double_t charge, Double_t step, 
			  Double_t* vect, Double_t* vout);
  Int_t ColorIndex(Float_t val);

  void  MUONTrackInfo();          // *MENU*
  void  MUONTriggerInfo();        // *MENU*

 private:

  AliMagF      *fFieldMap;           // pointer to the magnetic field map
  AliMUONTrack *fTrack;              // pointer to the MUON track
  Int_t         fCount;              // track points counter

  ClassDef(MUONTrack, 1);    // Produce Reve:Track from AliMUONTrack

};

}

#endif

