#ifndef ALISTARTHITPHOTON_H
#define ALISTARTHITPHOTON_H
 

#include "AliHit.h"
//------------------------------------------------
// Cerenkov's photon  object
//------------------------------------------------

class AliSTARThitPhoton: public AliHit {
 public:
 
    Int_t	fArray;		// Array number
    Int_t	fPmt;		// PMT number in the array
    Float_t	fTime;		// Time convention in photoelectron
    Float_t	fEtot;		// Energy photon
    Float_t	fMomX;		// Local Momentum
    Float_t	fMomY;		// Local Momentum
    Float_t	fMomZ;		// Local Momentum
    Float_t	fRadius;	// Distance from axis of radiatar
 
 public:
    AliSTARThitPhoton() {}
    AliSTARThitPhoton(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *hit);
    virtual ~AliSTARThitPhoton() {}
    
    ClassDef(AliSTARThitPhoton,1)  //Cerenkov's photons object for set:START
};
#endif
