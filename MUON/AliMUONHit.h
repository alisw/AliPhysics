#ifndef ALIMUONHIT_H
#define ALIMUONHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUONHit
/// \brief MonteCarlo hit
///
/// MUON classe for MonteCarlo Hits, inherited from AliHit for the 
/// In addition to the ALiHit data member fX, fY, fZ and fTrack, AliMUONHit contains some info about the particle crossing the chamber:
/// Impulsion: fPtot, fPx, fPy and fPz
/// Reference position at the center of the chamber (wire plane) fXref, fYref and fZref
/// Cumulated path along the active volume fTlength for spliting of hits for very inclined tracks 
/// Energy loss of the particle inside the gas active volume.
/// Incident fTheta and fPhi angle with respect of the wire plane of the chamber.


#include "AliHit.h"

class AliMUONHit : public AliHit {

 public:
    
    AliMUONHit();
    AliMUONHit(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *hits); 
               // TBR
    AliMUONHit(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *hits,
               Bool_t isNew);

    AliMUONHit(Int_t fIshunt, Int_t track, Int_t iChamber, Int_t idpart, 
               Float_t X, Float_t Y, Float_t Z, Float_t tof, Float_t momentum, 
	       Float_t theta, Float_t phi, Float_t length, Float_t destep);
	       // TBR
    AliMUONHit(Int_t fIshunt, Int_t track, Int_t detElemId, Int_t idpart, 
               Float_t X, Float_t Y, Float_t Z, Float_t tof, Float_t momentum, 
	       Float_t theta, Float_t phi, Float_t length, Float_t destep, Bool_t isNew);

    AliMUONHit(Int_t fIshunt, Int_t track, Int_t iChamber, Int_t idpart, 
               Float_t X, Float_t Y, Float_t Z, Float_t tof, Float_t momentum, 
               Float_t theta, Float_t phi, Float_t length, Float_t destep,
               Float_t Xref, Float_t Yref, Float_t Zref);
	       // TBR
    AliMUONHit(Int_t fIshunt, Int_t track, Int_t detElemId, Int_t idpart, 
               Float_t X, Float_t Y, Float_t Z, Float_t tof, Float_t momentum, 
               Float_t theta, Float_t phi, Float_t length, Float_t destep,
               Float_t Xref, Float_t Yref, Float_t Zref, Bool_t isNew);
    virtual ~AliMUONHit() {}

    Int_t   DetElemId()const;
    Int_t   Chamber()  const;
    Float_t Particle() const {return fParticle;}    
    Float_t Theta()    const {return fTheta;}
    Float_t Phi()      const {return fPhi;}
    Float_t Tlength()  const {return fTlength;}
    Float_t Eloss()    const {return fEloss;}
    Float_t Age()      const {return fAge;}
    Int_t   PHfirst()  const {return fPHfirst;}
    Int_t   PHlast()   const {return fPHlast;}

    Float_t Momentum() const {return fPTot;}
    Float_t Px()       const {return fPx;}
    Float_t Py()       const {return fPy;}
    Float_t Pz()       const {return fPz;}
    Float_t Cx()       const {return fPx/fPTot;} 
    Float_t Cy()       const {return fPy/fPTot;}
    Float_t Cz()       const {return fPz/fPTot;}

    Float_t Xref()     const {return fXref;}
    Float_t Yref()     const {return fYref;}
    Float_t Zref()     const {return fZref;}

 private:  
    Bool_t    fIsDetElemId;   // False if fDetElemId contains Chamber Id (old code) 
    Int_t     fDetElemId;     // Detection element ID
    Float_t   fParticle;      // Geant3 particle type
    Float_t   fTheta ;        // Incident theta angle in degrees      
    Float_t   fPhi   ;        // Incident phi angle in degrees
    Float_t   fTlength;       // Track length inside the chamber
    Float_t   fEloss;         // ionisation energy loss in gas
    Float_t   fAge;           // Particle Age
    Int_t     fPHfirst;       // first padhit
    Int_t     fPHlast;        // last padhit

    Float_t   fPTot;          // Local momentum P of the track when entering in the chamber
    Float_t   fPx;            // Px
    Float_t   fPy;            // Py
    Float_t   fPz;            // Pz
    
    Float_t   fXref;          // X position of hit in the center of the chamber (without angle effect)
    Float_t   fYref;          // Y position of hit in the center of the chamber (without angle effect)
    Float_t   fZref;          // Z position of hit in the center of the chamber (without angle effect)

    
    ClassDef(AliMUONHit,1)    //Hit object for MUON
};
#endif
