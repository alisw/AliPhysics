#ifndef ALITRACKRESIDUALS_H
#define ALITRACKRESIDUALS_H

//************************************************************************
// AliTrackResiduals: base class for collecting the track space point    *
// residuals produced by the fast track fitters (AliTrackFitter class).  *
// It provides an interface to the arrays which contain the space points *
// and track extrapolation points within the detector volume to be       *
// aligned. The derived classes should implement method to analyze the   *
// track residuals and minimize their sum in order to get the            *
// AliAlignObj for the given detector volume.                            *
//************************************************************************

#include "TObject.h"

#include "AliAlignObjAngles.h"

class AliTrackPointArray;

class AliTrackResiduals : public TObject {

 public:

  AliTrackResiduals();
  AliTrackResiduals(Int_t ntracks);
  AliTrackResiduals(const AliTrackResiduals &res);
  AliTrackResiduals& operator= (const AliTrackResiduals& res);
  virtual ~AliTrackResiduals();

  void   SetNTracks(Int_t ntracks);
  Bool_t AddTrackPointArrays(AliTrackPointArray *volarray, AliTrackPointArray *trackarray);
  void   SetAlignObj(AliAlignObj *alignobj);

  virtual Bool_t Minimize() = 0;

  Int_t  GetNTracks() const { return fN; }
  Int_t  GetNFilledTracks() const { return fLast; }
  Bool_t GetTrackPointArrays(Int_t i, AliTrackPointArray* &volarray, AliTrackPointArray* &trackarray) const;
  AliAlignObj *GetAlignObj() const { return fAlignObj; }
  Float_t GetChi2() const { return fChi2; }
  Int_t   GetNdf() const  { return fNdf; }

 protected:

  void DeleteTrackPointArrays();

  Int_t              fN;            // Number of tracks
  Int_t              fLast;         // Index of the last filled track arrays
  AliAlignObj        *fAlignObj;    // Pointer to the volume alignment object to be fitted
  AliTrackPointArray **fVolArray;   //! Pointers to the arrays containing space points
  AliTrackPointArray **fTrackArray; //! Pointers to the arrays containing track extrapolation points
  Float_t            fChi2;         // Chi2 (or distance) of residuals minimization
  Int_t              fNdf;          // Number of degrees of freedom
  Bool_t             fIsOwner;      // Track point arrays owned by the object

  ClassDef(AliTrackResiduals,1)

};

#endif
