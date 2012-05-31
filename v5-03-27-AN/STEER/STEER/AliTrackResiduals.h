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

#include "AliAlignObjParams.h"

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
  void   InitAlignObj();
  void   SetMinNPoints(Int_t n) { fMinNPoints = n; }

  virtual Bool_t Minimize() = 0;

  Int_t  GetNTracks() const { return fN; }
  Int_t  GetNFilledTracks() const { return fLast; }
  Bool_t GetTrackPointArrays(Int_t i, AliTrackPointArray* &volarray, AliTrackPointArray* &trackarray) const;
  AliAlignObj *GetAlignObj() const { return fAlignObj; }
  Float_t GetChi2() const { return fChi2; }
  Int_t   GetNdf() const  { return fNdf; }
  Int_t   GetMinNPoints() const  { return fMinNPoints; }
  void    FixParameter(Int_t par,Float_t value=0.) {fBFixed[par]=kTRUE; fFixed[par]= value;}
  Int_t GetNFreeParam();
  void   ReleaseParameter(Int_t par) {fBFixed[par]=kFALSE;}

 protected:

  void DeleteTrackPointArrays();

  Int_t              fN;            // Number of tracks
  Int_t              fLast;         // Index of the last filled track arrays
  AliAlignObj        *fAlignObj;    // Pointer to the volume alignment object to be fitted
  AliTrackPointArray **fVolArray;   //! Pointers to the arrays containing space points
  AliTrackPointArray **fTrackArray; //! Pointers to the arrays containing track extrapolation points
  Float_t            fChi2;         // Chi2 (or distance) of residuals minimization
  Int_t              fNdf;          // Number of degrees of freedom
  Int_t              fMinNPoints;   // Minimum allowed Number of points in the volume which is to be aligned
  Bool_t             fIsOwner;      // Track point arrays owned by the object
  Float_t            fFixed[6];     // The fixed values of parameters 
  Bool_t            fBFixed[6];    // The flag for fixing parameter

  ClassDef(AliTrackResiduals,2)

};

#endif
