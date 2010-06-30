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
#include <TMatrixD.h>
#include <TArrayI.h>
#include <TArrayF.h>

class AliMUONVTrackStore;
class AliMUONVTriggerStore;
class AliMUONVTriggerTrackStore;
class AliMUONTrackParam;
class AliMUONGeometryTransformer;
class AliMUONVDigitStore;
class AliMUONTriggerTrack;
class AliMUONTrack;
class AliMUONRecoParam;

class AliMUONTrackHitPattern : public TObject 
{
public:

  AliMUONTrackHitPattern(const AliMUONRecoParam* recoParam,
                         const AliMUONGeometryTransformer& transformer,
                         const AliMUONVDigitStore& digitStore);
  virtual ~AliMUONTrackHitPattern(); // Destructor

  void ExecuteValidation(const AliMUONVTrackStore& trackStore,
			 const AliMUONVTriggerTrackStore& triggerTrackStore,
			 const AliMUONVTriggerStore& triggerStore) const;
  
  AliMUONTriggerTrack* MatchTriggerTrack(AliMUONTrack* track,
					 AliMUONTrackParam& trackParam,
					 const AliMUONVTriggerTrackStore& triggerTrackStore,
					 const AliMUONVTriggerStore& triggerStore) const;

  UShort_t GetHitPattern(AliMUONTriggerTrack* matchedTriggerTrack) const;
  
  UShort_t GetHitPattern(AliMUONTrackParam* trackParam) const;

protected:
  void ApplyMCSCorrections(AliMUONTrackParam& trackParam) const;
  
  // Methods for hit pattern from tracker track
  void FindPadMatchingTrack(const AliMUONTrackParam& trackParam,
			    Bool_t isMatch[2], Int_t iChamber) const;

  Float_t MinDistanceFromPad(Float_t xPad, Float_t yPad, Float_t zPad,
			     Float_t dpx, Float_t dpy, 
			     const AliMUONTrackParam& trackParam) const;

  // Methods for hit pattern from matched trigger track
  Bool_t PerformTrigTrackMatch(UShort_t &pattern,
			       const AliMUONTriggerTrack* matchedTrigTrack) const;
  
  Int_t FindPadMatchingTrig(Int_t &detElemId, Float_t coor[2],
			    Bool_t isMatch[2], TArrayI nboard[2],
			    TArrayF &zRealMatch, Float_t y11) const;
  
  Float_t PadMatchTrack(Float_t xPad, Float_t yPad, Float_t dpx, Float_t dpy,
			Float_t xTrackAtPad, Float_t yTrackAtPad) const;
  
  Int_t DetElemIdFromPos(Float_t x, Float_t y, Int_t chamber, Int_t cathode) const;
  
  void LocalBoardFromPos(Float_t x, Float_t y, Int_t detElemId,
			 Int_t cathode, Int_t localBoard[4]) const;

  /// Return reco parameters
  const AliMUONRecoParam* GetRecoParam() const { return fkRecoParam; }
  
private:
  /// Not implemented
  AliMUONTrackHitPattern(const AliMUONTrackHitPattern& rhs);
  /// Not implemented
  AliMUONTrackHitPattern& operator = (const AliMUONTrackHitPattern& rhs);

  void CheckConstants() const;

  const AliMUONRecoParam* fkRecoParam; //!< pointer to reco parameters
  const AliMUONGeometryTransformer& fkTransformer; //!< geometry transformer
  const AliMUONVDigitStore& fkDigitStore; //!< digitStore

  const Float_t fkMaxDistance; //!< Maximum distance for reference
  static const Int_t fgkNcathodes=2; //!<Number of cathodes
  static const Int_t fgkNchambers=4; //!<Number of chambers
  static const Int_t fgkNplanes=8;   //!<Number of planes
  static const Int_t fgkNlocations=4; //!<Number of locations
  
  ClassDef(AliMUONTrackHitPattern, 0) // MUON track hit pattern
};

#endif
