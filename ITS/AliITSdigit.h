#ifndef ITSDIGIT_H
#define ITSDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliDigit.h"
#include "AliITS.h"
#include "AliITSgeom.h"

//___________________________________________
class AliITSdigit: public AliDigit  {
////////////////////////////////////////////////////////////////////////
// Version: 0
// Written by Rene Brun, Federico Carminati, and Roberto Barbera
// Minor modifications made and documented by Bjorn S. Nilsen
// July 11 1999
//
// The default ITS digit structure. This should either be replaced
// or added on to later with the proper digit structure defined for
// each detector type. See the proposed Digit structure defined by
// Bjorn S. Nilsen for an example.
//
// Data members:
//
// Int_t fTracks[3]
//     See AliDigit for a full description. The track numbers, up to 3,
// that make up this digit.
//
// Int_t fEvent
//     The event number for this digit. This information is probably
// kept someplace else already (via the TTree structure already in use).
//
// Int_t fLayer
//     The layer number of this digit. This is part of the information
// that determines the detector where this digit is located (layer, ladder,
// and detector numbers).
//
// Int_t fLadder
//     The ladder number of this digit. This is part of the information
// that determines the detector where this digit is located (layer, ladder,
// and detector numbers).
//
// Int_t fDet
//     The detector number of this digit. This is part of the information
// that determines the detector where this digit is located (layer, ladder,
// and detector numbers).
//
// Int_t fNoverl
//     The number of hits that make up this digit.
//
// Member functions:
//
// int *GetTracks()
//     See AliDigit for a full description. Returns a pointer to the
// array fTracks where the tracks number of the tracks that make up
// this digit are stored.
//
// AliITSdigit()
//     The default creator for the AliITSdigit class.
//
// AliITSdigit(Int_t *tracks, Int_t *digits)
//     The creator for the AliITSdigit class. This routine fills the
// AliITSdigit data members from the array digits. The array of track
// numbers are passed to the AliDigit creator. The order of the elements
// in the digits array are fEvent = digits[0], fLayer = digits[1],
// fLadder = digits[2], fDet = digits[3], and fNoverl = digits[4].
// Therefore the array digits is expected to be at least 5 elements long.
//
// ~AliITSdigit()
//     The destructor for the AliITSdigit class. At present the default
// destructor is used.
////////////////////////////////////////////////////////////////////////
 public:
    Int_t fEvent;      // Event number
    Int_t fLayer;      // Layer number
    Int_t fLadder;     // Ladder number
    Int_t fDet;        // Detector number
    Int_t fNoverl;     // Number of overflow

 public:
    AliITSdigit() {}
    AliITSdigit(Int_t *tracks, Int_t *digits);
    virtual ~AliITSdigit() {}

    ClassDef(AliITSdigit,1)  //Digit (Header) object for set:ITS
};
#endif
