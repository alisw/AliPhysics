#ifndef ALIITSUTRACKHYP_H
#define ALIITSUTRACKHYP_H

#include <TObject.h>
#include <TObjArray.h>
#include "AliITSUSeed.h"


// Container for track hypotheses

class AliITSUTrackHyp: public TObject
{
 public:
  AliITSUTrackHyp(Int_t nlr=0);
  AliITSUTrackHyp(const AliITSUTrackHyp& src);
  AliITSUTrackHyp &operator=(const AliITSUTrackHyp &src);
  virtual ~AliITSUTrackHyp();
  //
  Int_t              GetNLayers()        const {return fNLayers;}
  Int_t              GetNSeeds(Int_t lr) const {return lr<fNLayers ? fLayerSeeds[lr].GetEntriesFast() : 1;}
  AliITSUSeed*       GetSeed(Int_t lr, Int_t id) const;
  AliITSUSeed*       GetESDSeed()        const {return fESDSeed;}
  const TObjArray*   GetLayerSeeds(Int_t lr) const {return lr<fNLayers ? &fLayerSeeds[lr] : 0;}
  void               AddSeed(AliITSUSeed* seed, Int_t lr);
  void               SetESDSeed(AliITSUSeed* seed) {fESDSeed = seed;}
  //
  
  //
  virtual void       Print(Option_t* option = "") const;
  //
 protected:
  UChar_t          fNLayers;               // number of layers
  AliITSUSeed*     fESDSeed;               // bare esd (TPC) seed
  TObjArray*       fLayerSeeds;            // seeds of given layer
  //
  ClassDef(AliITSUTrackHyp,1)
};

//___________________________________________________________________
inline AliITSUSeed* AliITSUTrackHyp::GetSeed(Int_t lr, Int_t id) const
{
  //return requested seed of given layer, no check is done on seed index
  return lr<fNLayers ? (AliITSUSeed*)fLayerSeeds[lr].UncheckedAt(id) : fESDSeed;
}

//___________________________________________________________________
inline void AliITSUTrackHyp::AddSeed(AliITSUSeed* seed, Int_t lr)
{
  //return requested seed of given layer, no check is done
  if (lr<fNLayers) fLayerSeeds[lr].AddLast(seed);
  else             fESDSeed = seed;
}

#endif
