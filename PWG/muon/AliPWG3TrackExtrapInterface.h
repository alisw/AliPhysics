#ifndef AliPWG3TRACKEXTRAPINTERFACE_H
#define AliPWG3TRACKEXTRAPINTERFACE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

/// \ingroup rec
/// \class AliPWG3TrackExtrapInterface
/// \brief Interface to MUON classes
///
//////////////////////////////////////////
/// Interface for MUON track extrapolation
//////////////////////////////////////////

#include <TObject.h>

class AliMagF;
class AliESDMuonTrack;

class AliPWG3TrackExtrapInterface : public TObject 
{
 public:
	/// Constructor
  AliPWG3TrackExtrapInterface() : TObject(){};
	/// Destructor
  virtual ~AliPWG3TrackExtrapInterface(){};
  
  static void	SetMagField(const AliMagF* magField);
  static Bool_t SetGeometry(const char* fileName = "geometry.root");
  
  // extrapolation without any absorber correction
  static void ExtrapToVertexUncorrected(AliESDMuonTrack* muonTrack, Double_t zVtx);
  // extrapolation corrected for energy loss in absorber only
  static void ExtrapToVertexWithELoss(AliESDMuonTrack* muonTrack, Double_t zVtx);
  // extrapolation corrected for multiple scattering in absorber only
  static void ExtrapToVertexWithBranson(AliESDMuonTrack* muonTrack, Double_t vtx[3]);
  // extrapolation corrected for all absorber effect (multiple scattering and energy loss)
  static void ExtrapToVertex(AliESDMuonTrack* muonTrack, Double_t vtx[3]);
  
  static Double_t TotalMomentumEnergyLoss(AliESDMuonTrack* muonTrack, Double_t zVtx);
  static Double_t TotalMomentumEnergyLoss(AliESDMuonTrack* muonTrack, Double_t vtx[3]);
  
 private:
  // Functions
  AliPWG3TrackExtrapInterface(const AliPWG3TrackExtrapInterface& muonInterface);
  AliPWG3TrackExtrapInterface& operator=(const AliPWG3TrackExtrapInterface& muonInterface);
  
  
  ClassDef(AliPWG3TrackExtrapInterface, 0) // Tools for track extrapolation in ALICE dimuon spectrometer
};
	
#endif
