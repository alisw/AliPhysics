#ifndef GET_GOOD_PARTICLES
#define GET_GOOD_PARTICLES

struct GoodTrack 
{

  Int_t label;
  Double_t eta;
  Int_t code;
  Double_t px,py,pz;
  Double_t pt;
  Int_t nhits;
};


void GetGoodParticles(Int_t,Int_t,Char_t *,Char_t *);

#endif
