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
class AliMUONDigitMaker;
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
                         const AliMUONDigitMaker& digitMaker);
  virtual ~AliMUONTrackHitPattern(); // Destructor

  void ExecuteValidation(const AliMUONVTrackStore& trackStore,
			 const AliMUONVTriggerTrackStore& triggerTrackStore,
			 const AliMUONVTriggerStore& triggerStore) const;
  
  AliMUONTriggerTrack* MatchTriggerTrack(AliMUONTrack* track,
					 AliMUONTrackParam& trackParam,
					 const AliMUONVTriggerTrackStore& triggerTrackStore,
					 const AliMUONVTriggerStore& triggerStore) const;
    
  UShort_t GetHitPattern(AliMUONTriggerTrack* matchedTriggerTrack,
			 AliMUONVDigitStore& digitStore,
			 AliMUONTrackParam* trackParam=0x0) const;

protected:
  void ApplyMCSCorrections(AliMUONTrackParam& trackParam) const;
  
  void InitMembers();
  
  void SetBit(UShort_t& pattern, Int_t cathode, Int_t chamber) const;
  
  void AddEffInfo(UShort_t& pattern, Int_t slat, Int_t effType) const;
  

  // Methods for hit pattern from tracker track
  void FindPadMatchingTrack(const AliMUONVDigitStore& digitStore,
			    const AliMUONTrackParam& trackParam,
			    Bool_t isMatch[2], Int_t iChamber) const;

  Float_t MinDistanceFromPad(Float_t xPad, Float_t yPad, Float_t zPad,
			     Float_t dpx, Float_t dpy, 
			     const AliMUONTrackParam& trackParam) const;

  // Methods for hit pattern from matched trigger track
  Bool_t PerformTrigTrackMatch(UShort_t &pattern,
			       const AliMUONTriggerTrack* matchedTrigTrack,
			       AliMUONVDigitStore& digitStore) const;
  
  Int_t FindPadMatchingTrig(const AliMUONVDigitStore& digitStore, Int_t &detElemId, Float_t coor[2],
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
  const AliMUONDigitMaker& fkDigitMaker; //!< pointer to digit maker

  Double_t fDeltaZ; //!< distance between stations

  TMatrixD* fTrigCovariance; //!< Covariance matrix 3x3 (X,Y,slopeY) for trigger tracks

  const Float_t fkMaxDistance; //!< Maximum distance for reference
  static const Int_t fgkNcathodes=2; //!<Number of cathodes
  static const Int_t fgkNchambers=4; //!<Number of chambers
  static const Int_t fgkNplanes=8;   //!<Number of planes
  static const Int_t fgkNlocations=4; //!<Number of locations

  enum {
    kNoEff,
    kChEff,
    kSlatEff,
    kBoardEff
  };
  
  ClassDef(AliMUONTrackHitPattern, 0) // MUON track hit pattern
};

#endif
