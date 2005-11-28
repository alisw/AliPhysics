#ifndef ALITRACKFITTER_H
#define ALITRACKFITTER_H

/*************************************************************************
 * AliTrackFitter: base class for the fast track fitters                 *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *************************************************************************/

#include "TObject.h"

#include "AliTrackPointArray.h"
#include "AliAlignObj.h"

class TMatrixDSym;

class AliTrackFitter : public TObject {

 public:

  AliTrackFitter();
  AliTrackFitter(AliTrackPointArray *array, Bool_t owner = kTRUE);
  AliTrackFitter(const AliTrackFitter &fitter);
  AliTrackFitter& operator= (const AliTrackFitter& fitter);
  virtual ~AliTrackFitter();

  virtual void   Reset();
  virtual void   SetTrackPointArray(AliTrackPointArray *array, Bool_t owner = kTRUE);
  virtual Bool_t Fit(UShort_t volId,
		     AliTrackPointArray *pVolId, AliTrackPointArray *pTrack,
		     AliAlignObj::ELayerID layerRangeMin = AliAlignObj::kFirstLayer,
		     AliAlignObj::ELayerID layerRangeMax = AliAlignObj::kLastLayer) = 0;
  virtual Bool_t GetPCA(const AliTrackPoint &pIn, AliTrackPoint &pOut) const = 0;

  const Float_t* GetX() const {return fPoints->GetX();}
  const Float_t* GetY() const {return fPoints->GetY();}
  const Float_t* GetZ() const {return fPoints->GetZ();}
  const Double_t* GetParam() const {return &fParams[0];}
  const TMatrixDSym &  GetCovariance() const {return *fCov;}

 protected:

  Double_t      fParams[6];    // Track parameters
  TMatrixDSym   *fCov;         // Track cov matrix
  AliTrackPointArray *fPoints; // Pointer to the array with track space points
  Bool_t  fIsOwner;            // Is the object owner of the space points array

 private:

  ClassDef(AliTrackFitter,1) // Abstract class of fast track fitters

};

#endif
