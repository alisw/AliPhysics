#ifndef ALIMUONTRACKHITPATTERN_H
#define ALIMUONTRACKHITPATTERN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup rec
/// \class AliMUONTrackHitPattern
/// \brief Class for the MUON track hit pattern
///
//  Author: Diego Stocco

#include <TObject.h>

class AliMUONVTrackStore;
class AliMUONVTriggerStore;
class AliMUONTrackParam;
class AliMUONDigitMaker;
class AliMUONGeometryTransformer;
class AliMUONVDigitStore;

class AliMUONTrackHitPattern : public TObject 
{
public:

  AliMUONTrackHitPattern(const AliMUONGeometryTransformer& transformer,
                         const AliMUONDigitMaker& digitMaker);
  virtual ~AliMUONTrackHitPattern(); // Destructor
    
    void GetHitPattern(AliMUONVTrackStore& trackStore,
                       const AliMUONVTriggerStore& triggerStore) const;
    
    void FindPadMatchingTrack(AliMUONVDigitStore& digitStore,
                              const AliMUONTrackParam& trackParam,
                              Bool_t isMatch[2], Int_t iChamber) const;

    Float_t MinDistanceFromPad(Float_t xPad, Float_t yPad, Float_t zPad,
                               Float_t dpx, Float_t dpy, 
                               const AliMUONTrackParam& trackParam) const;
    
    void ApplyMCSCorrections(AliMUONTrackParam& trackParam) const;

    Bool_t TriggerDigits(const AliMUONVTriggerStore& triggerStore, 
                         AliMUONVDigitStore& digitStore) const;

private:
    /// Not implemented
    AliMUONTrackHitPattern(const AliMUONTrackHitPattern& rhs);
    /// Not implemented
    AliMUONTrackHitPattern& operator = (const AliMUONTrackHitPattern& rhs);

    const AliMUONGeometryTransformer& fTransformer; //!< geometry transformer
    const AliMUONDigitMaker& fDigitMaker; ///< pointer to digit maker

    ClassDef(AliMUONTrackHitPattern, 0) // MUON track hit pattern
};

#endif
