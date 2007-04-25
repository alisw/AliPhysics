#ifndef ALIMUONTRACKHITPATTERN_H
#define ALIMUONTRACKHITPATTERN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup rec
/// \class AliMUONTrackHitPattern
/// \brief Class for the MUON track hit pattern
///
//  Author: Diego Stocco

#include <TObject.h>
#include <TList.h>

class AliMUONData;
class AliMUONLocalStruct;
class AliMUONLocalTriggerBoard;
class AliMUONTrackParam;
class AliMUONTriggerCrateStore;

class AliMUONGeometryTransformer;

class TClonesArray;

class AliMUONTrackHitPattern : public TObject {
public:
    AliMUONTrackHitPattern(AliMUONData *MUONData); // Default Constructor
    virtual ~AliMUONTrackHitPattern(); // Destructor
    
    void GetHitPattern(TClonesArray *recTracksArray);
    
    void FindPadMatchingTrack(AliMUONTrackParam *trackParam, Bool_t isMatch[2], Int_t iChamber);
    Float_t MinDistanceFromPad(Float_t xPad, Float_t yPad, Float_t zPad,
			       Float_t dpx, Float_t dpy, AliMUONTrackParam *trackParam);
    void GetPosUncertainty(AliMUONTrackParam *trackParm, Float_t zChamber,
			   Float_t &sigmaX, Float_t &sigmaY, Float_t &sigmaMS);

    Bool_t TriggerDigits();

private:
    AliMUONData *fMUONData; ///< Data container for MUON subsystem
    TList fTriggerDigitsList[4]; ///< List of trigger digits, one per chamber
    AliMUONGeometryTransformer *fTransformer; //!< pointer to transformation
    AliMUONTriggerCrateStore *fCrateManager; ///< pointer to crate manager

    ClassDef(AliMUONTrackHitPattern, 0) // MUON track hit pattern
};

#endif
