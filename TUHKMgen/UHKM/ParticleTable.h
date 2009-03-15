#ifndef PARTICLETABLE_INCLUDED
#define PARTICLETABLE_INCLUDED

/*                                                                            
                                                                            
        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
                           November. 2, 2005                                

*/

#include <map>

#include <Rtypes.h>

struct ParticleInfo {
  Int_t fBaryonNumber;
  Int_t fStrangeness;
  Int_t fIsospin;
  Int_t fSpin;
  Int_t fCharge;

  ParticleInfo(Int_t bN, Int_t s, Int_t s1, Int_t s2, Int_t c) {
    fBaryonNumber = bN;
    fStrangeness = s;
    fIsospin = s1; //2S
    fSpin = s2; //2I
    fCharge = c; //fCharge = 2 * I3
  }
};

extern const std::map<const Int_t, ParticleInfo> gParticleTable;
typedef std::map<const Int_t, ParticleInfo>::const_iterator MapIt_t;

#endif
