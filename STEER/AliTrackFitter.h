#ifndef ALITRACKFITTER_H
#define ALITRACKFITTER_H

/*************************************************************************
 * AliTrackFitter: base class for the fast track fitters                 *
 *                                                                       *
 *                                                                       *
 *                                                                       *
 *************************************************************************/

#include <TObject.h>
#include <TMatrixDSymfwd.h>

#include "AliTrackPointArray.h"
#include "AliAlignObj.h"

class TArrayI;

class AliTrackFitter : public TObject {

 public:

  AliTrackFitter();
  AliTrackFitter(AliTrackPointArray *array, Bool_t owner = kTRUE);
  AliTrackFitter(const AliTrackFitter &fitter);
  AliTrackFitter& operator= (const AliTrackFitter& fitter);
  virtual ~AliTrackFitter();

  virtual void   Reset();
  virtual void   SetTrackPointArray(AliTrackPointArray *array, Bool_t owner = kTRUE);
  virtual Bool_t Fit(const TArrayI *volIds,const TArrayI *volIdsFit = 0x0,
		     AliGeomManager::ELayerID layerRangeMin = AliGeomManager::kFirstLayer,
		     AliGeomManager::ELayerID layerRangeMax = AliGeomManager::kLastLayer);

  virtual Bool_t Begin(Int_t, Int_t) = 0;
  virtual Bool_t AddPoint(const AliTrackPoint *p) = 0;
  virtual Bool_t Update() = 0;

  virtual Bool_t GetPCA(const AliTrackPoint &pIn, AliTrackPoint &pOut) const = 0;

  Bool_t         FindVolId(const TArrayI *array, UShort_t volid) const;

  void           SetMinNPoints(Int_t n) { fMinNPoints = n;}

  const Float_t* GetX() const {return fPoints->GetX();}
  const Float_t* GetY() const {return fPoints->GetY();}
  const Float_t* GetZ() const {return fPoints->GetZ();}
  const Double_t* GetParam() const {return &fParams[0];}
  const TMatrixDSym &  GetCovariance() const {return *fCov;}
  Float_t        GetChi2() const {return fChi2;}
  Int_t          GetNdf()  const {return fNdf;}
  Int_t          GetMinNPoints()  const {return fMinNPoints;}
  Float_t        GetNormChi2() const { return (fNdf != 0) ? fChi2/fNdf : 0; }
  void           GetTrackResiduals(AliTrackPointArray*& pVolId, AliTrackPointArray*& pTrack) const
    { pVolId = fPVolId; pTrack = fPTrack; }

 protected:

  Double_t      fParams[6];    // Track parameters
  TMatrixDSym  *fCov;          // Track cov matrix
  AliTrackPointArray *fPoints; // Pointer to the array with track space points
  AliTrackPointArray *fPVolId; // Pointer to the array with space-points in volId
  AliTrackPointArray *fPTrack; // Pointer to the array with track extrapolation points in volId
  Float_t       fChi2;         // Chi squared of the fit
  Int_t         fNdf;          // Number of degrees of freedom
  Int_t         fMinNPoints;   // Minimum allowed number of points
  Bool_t  fIsOwner;            // Is the object owner of the space points array

 private:

  ClassDef(AliTrackFitter,1) // Abstract class of fast track fitters

};

#endif
