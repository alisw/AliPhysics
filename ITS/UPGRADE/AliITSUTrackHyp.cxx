#include "AliITSUTrackHyp.h"

ClassImp(AliITSUTrackHyp)



//__________________________________________________________________
AliITSUTrackHyp::AliITSUTrackHyp(Int_t nlr) 
: fNLayers(nlr)
  ,fESDSeed(0)
  ,fLayerSeeds(0)
{
  // def. c-tor
  if (fNLayers>0) fLayerSeeds = new TObjArray[fNLayers];
}

//__________________________________________________________________
AliITSUTrackHyp::~AliITSUTrackHyp() 
{
  // d-tor
  delete[] fLayerSeeds;
}

//__________________________________________________________________
AliITSUTrackHyp::AliITSUTrackHyp(const AliITSUTrackHyp &src)
  : TObject(src)
  , fNLayers(src.fNLayers)
  , fESDSeed(src.fESDSeed)
  , fLayerSeeds(0)
{
  // copy c-tor
  if (fNLayers>0) {
    fLayerSeeds = new TObjArray[fNLayers];
    for (int ilr=fNLayers;ilr--;) {
      int ns = src.GetNSeeds(ilr);
      for (int isd=0;isd<ns;isd++) {
	AliITSUSeed* sd = src.GetSeed(ilr,isd);
	if (sd->IsKilled()) continue;
	AddSeed(sd,ilr);
      }      
    }
  }
  //
}

//__________________________________________________________________
AliITSUTrackHyp &AliITSUTrackHyp::operator=(const AliITSUTrackHyp &src)
{
  // copy 
  if (this == &src) return *this;
  this->~AliITSUTrackHyp();
  new(this) AliITSUTrackHyp(src);
  return *this;
  //
}

//__________________________________________________________________
void AliITSUTrackHyp::Print(Option_t* ) const
{
  printf("Track Hyp.#%4d. NSeeds:",GetUniqueID());
  for (int i=0;i<fNLayers;i++) printf(" (%d) %3d",i,GetNSeeds(i)); printf("\n");
}
