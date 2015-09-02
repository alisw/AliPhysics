#ifndef AliMFTTrack_H
#define AliMFTTrack_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup MFTrec
/// \class AliMFTTrack
/// \brief Description of an ALICE Standalone MFT track
///
///
/// Description of an ALICE Standalone MFT track
///
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>, IPN-Lyon
/// \date April 27th, 2015

#include "TObject.h"
class AliMFTCATrack;
class AliMFTTrackParam;

//=============================================================================================

class AliMFTTrack : public TObject {
  
  public:
  
  AliMFTTrack();
  AliMFTTrack(AliMFTCATrack *catrack);
  virtual ~AliMFTTrack();

  /// return the minimum value of the function minimized by the fit
  Double_t GetChi2() const {return fChi2;}
  /// set the minimum value of the function minimized by the fit
  void SetChi2(Double_t chi2) { fChi2 = chi2;}

  /// return pointer to track parameters at vertex (can be 0x0)
  AliMFTTrackParam* GetTrackParamAtVertex() const {return fTrackParamAtVertex;}
  void SetTrackParamAtVertex(  AliMFTTrackParam* trackParam){};
 
  /// return pointer to track found by Track Finder (includes clusters)
  AliMFTCATrack* GetCATrack() const {return fCATrack;}
<<<<<<< HEAD
  void SetCATrack(  AliMFTCATrack* track) {};
=======
  void SetCATrack(  AliMFTCATrack* track){};
>>>>>>> Add temporary a simple empty content {} to the definition of the function.

  protected:
  
  Double_t fChi2;                         ///<  Chi2 of the track
  AliMFTTrackParam* fTrackParamAtVertex;  ///<  Track parameters at vertex
  AliMFTCATrack* fCATrack;                ///<  Track found by Track Finder (includes clusters)
  
  /// \cond CLASSIMP
  ClassDef(AliMFTTrack,1);
  /// \endcond
};

//=============================================================================================

#endif
