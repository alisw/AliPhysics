#ifndef ALITRDSEEDV1_H
#define ALITRDSEEDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The TRD track seed                                                    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ROOT_TMath
#include "TMath.h"
#endif

#ifndef ALITRDGEOMETRY_H
#include "AliTRDgeometry.h"
#endif

#ifndef ALIPID_H
#include "AliPID.h"
#endif

#ifndef ALIRIEMAN_H
#include "AliRieman.h"
#endif

class TTreeSRedirector;

class AliRieman;

class AliTRDtrackingChamber;
class AliTRDcluster;
class AliTRDtrackV1;
class AliTRDReconstructor;
class AliTRDseedV1 : public TObject //TODO we should inherit 
{                                   // AliTRDtrackletBase
public:
  enum ETRDtrackletBuffers {    
    kNtb = 32           // max clusters/pad row
   ,kNTimeBins = 2*kNtb // max number of clusters/tracklet
   ,kNSlices = 10       // max dEdx slices
  };

  // bits from 0-13 are reserved by ROOT (see TObject.h)
  enum ETRDtrackletStatus {
    kOwner    = BIT(14) // owner of its clusters
   ,kRowCross = BIT(15) // pad row cross tracklet
   ,kCalib    = BIT(16) // calibrated tracklet
  };

  AliTRDseedV1(Int_t det = -1);
  ~AliTRDseedV1();
  AliTRDseedV1(const AliTRDseedV1 &ref);
  AliTRDseedV1& operator=(const AliTRDseedV1 &ref);

  Bool_t	  AttachClustersIter(
              AliTRDtrackingChamber *chamber, Float_t quality, 
              Bool_t kZcorr = kFALSE, AliTRDcluster *c=0x0);
  Bool_t	  AttachClusters(
              AliTRDtrackingChamber *chamber, Bool_t tilt = kFALSE);
  void      Bootstrap(const AliTRDReconstructor *rec);
  void      Calibrate();
  void      CookdEdx(Int_t nslices);
  void      CookLabels();
  Bool_t    Fit(Bool_t tilt=kTRUE, Int_t errors = 2);
  void      FitMI();
  Bool_t    Init(AliTRDtrackV1 *track);
  inline void      Init(const AliRieman *fit);
  Bool_t    IsEqual(const TObject *inTracklet) const;
  Bool_t    IsCalibrated() const     { return TestBit(kCalib);}
  Bool_t    IsOwner() const          { return TestBit(kOwner);}
  Bool_t    IsOK() const             { return fN2 > 4;}
  Bool_t    IsRowCross() const       { return TestBit(kRowCross);}
  Bool_t    IsUsable(Int_t i) const  { return TESTBIT(fUsable, i);}

  Float_t   GetC() const             { return fC; }
  Float_t   GetChi2() const          { return fChi2; }
  inline Float_t   GetChi2Z() const;
  inline Float_t   GetChi2Y() const;
  static void      GetClusterXY(const AliTRDcluster *c, Double_t &x, Double_t &y);
  void      GetCovAt(Double_t x, Double_t *cov) const;
  void      GetCovXY(Double_t *cov) const { memcpy(cov, &fCov[0], 3*sizeof(Double_t));}
  void      GetCovRef(Double_t *cov) const { memcpy(cov, &fRefCov[0], 3*sizeof(Double_t));}
  Float_t   GetdX() const            { return fdX;}
  Float_t*  GetdEdx()                { return &fdEdx[0];}
  Float_t   GetdQdl(Int_t ic) const;
  Int_t     GetDetector() const      { return fDet;}
  void      GetCalibParam(Float_t &exb, Float_t &vd, Float_t &t0, Float_t &s2, Float_t &dl, Float_t &dt) const    { 
              exb = fExB; vd = fVD; t0 = fT0; s2 = fS2PRF; dl = fDiffL; dt = fDiffT;}
  AliTRDcluster*  GetClusters(Int_t i) const               { return i<0 || i>=kNTimeBins ? 0x0 : fClusters[i];}
  Int_t     GetIndexes(Int_t i) const{ return i<0 || i>=kNTimeBins ? -1 : fIndexes[i];}
  Int_t     GetLabels(Int_t i) const { return fLabels[i];}  
  Double_t  GetMomentum() const      { return fMom;}
  Int_t     GetN() const             { return fN2;}
  Int_t     GetN2() const            { return fN2;}
  Int_t     GetNUsed() const         { return fNUsed;}
  Float_t   GetQuality(Bool_t kZcorr) const;
  Float_t   GetPadLength() const     { return fPadLength;}
  Int_t     GetPlane() const         { return AliTRDgeometry::GetLayer(fDet);    }

  Float_t*  GetProbability();
  Float_t   GetS2Y() const           { return fS2Y;}
  Float_t   GetS2Z() const           { return fS2Z;}
  Float_t   GetSigmaY() const        { return fS2Y > 0. ? TMath::Sqrt(fS2Y) : 0.2;}
  Float_t   GetSnp() const           { return fYref[1]/TMath::Sqrt(1+fYref[1]*fYref[1]);}
  Float_t   GetTgl() const           { return fZref[1];}
  Float_t   GetTilt() const          { return fTilt;}
  Float_t   GetX0() const            { return fX0;}
  Float_t   GetX() const             { return fX0 - fX;}
  Float_t   GetY() const             { return fYfit[0] - fYfit[1] * fX;}
  Double_t  GetYat(Double_t x) const { return fYfit[0] - fYfit[1] * (fX0-x);}
  Float_t   GetYfit(Int_t id) const { return fYfit[id];}
  Float_t   GetYref(Int_t id) const { return fYref[id];}
  Float_t   GetZ() const            { return fZfit[0] - fZfit[1] * fX;}
  Double_t  GetZat(Double_t x) const { return fZfit[0] - fZfit[1] * (fX0-x);}
  Float_t   GetZfit(Int_t id) const { return fZfit[id];}
  Float_t   GetZref(Int_t id) const { return fZref[id];}
  Long_t    GetUsabilityMap() const { return fUsable; }

  inline AliTRDcluster* NextCluster();
  inline AliTRDcluster* PrevCluster();
  void      Print(Option_t *o = "") const;
  inline void ResetClusterIter(Bool_t forward = kTRUE);
  void      Reset();

  void      SetC(Float_t c)         { fC = c;}
  void      SetChi2(Float_t chi2)   { fChi2 = chi2;}
  void      SetCovRef(const Double_t *cov) { memcpy(&fRefCov[0], cov, 3*sizeof(Double_t));}
  void      SetIndexes(Int_t i, Int_t idx) { fIndexes[i]  = idx; }
  void      SetLabels(Int_t *lbls)   { memcpy(fLabels, lbls, 3*sizeof(Int_t)); }
  void      SetMomentum(Double_t mom){ fMom = mom;}
  void      SetOwner();
  void      SetTilt(Float_t tilt)    { fTilt = tilt; }
  void      SetPadLength(Float_t len){ fPadLength = len;}
  void      SetDetector(Int_t d)     { fDet = d;  }
  void      SetDX(Float_t inDX)      { fdX = inDX;}
  void      SetReconstructor(const AliTRDReconstructor *rec) {fReconstructor = rec;}
  void      SetX0(Float_t x0)        { fX0 = x0; }
  void      SetYref(Int_t i, Float_t y) { fYref[i]     = y;}
  void      SetZref(Int_t i, Float_t z) { fZref[i]     = z;}
  void      SetUsabilityMap(Long_t um)  { fUsable = um; }
  void      UpDate(const AliTRDtrackV1* trk);
  void      UpdateUsed();
  void      UseClusters();

protected:
  void Copy(TObject &ref) const;

private:
  const AliTRDReconstructor *fReconstructor;//! local reconstructor
  AliTRDcluster  **fClusterIter;            //! clusters iterator
  Int_t            fIndexes[kNTimeBins];    //! Indexes
  Float_t          fExB;                    //! tg(a_L) @ tracklet location
  Float_t          fVD;                     //! drift velocity @ tracklet location
  Float_t          fT0;                     //! time 0 @ tracklet location
  Float_t          fS2PRF;                  //! sigma^2 PRF for xd->0 and phi=a_L 
  Float_t          fDiffL;                  //! longitudinal diffusion coefficient
  Float_t          fDiffT;                  //! transversal diffusion coefficient
  Char_t           fClusterIdx;             //! clusters iterator
  Long_t           fUsable;                 //! bit map of usable clusters
  UChar_t          fN2;                     // number of clusters attached
  UChar_t          fNUsed;                  // number of used usable clusters
  Short_t          fDet;                    // TRD detector
  Float_t          fTilt;                   // local tg of the tilt angle 
  Float_t          fPadLength;              // local pad length 
  AliTRDcluster   *fClusters[kNTimeBins];   // Clusters
  Float_t          fYref[2];                //  Reference y
  Float_t          fZref[2];                //  Reference z
  Float_t          fYfit[2];                //  Y fit position +derivation
  Float_t          fZfit[2];                //  Z fit position
  Float_t          fMom;                    //  Momentum estimate @ tracklet [GeV/c]
  Float_t          fdX;                     // length of time bin
  Float_t          fX0;                     // anode wire position
  Float_t          fX;                      // radial position of the tracklet
  Float_t          fY;                      // r-phi position of the tracklet
  Float_t          fZ;                      // z position of the tracklet
  Float_t          fS2Y;                    // estimated resolution in the r-phi direction 
  Float_t          fS2Z;                    // estimated resolution in the z direction 
  Float_t          fC;                      // Curvature
  Float_t          fChi2;                   // Global chi2  
  Float_t          fdEdx[kNSlices];         // dE/dx measurements for tracklet
  Float_t          fProb[AliPID::kSPECIES]; //  PID probabilities
  Int_t            fLabels[3];              // most frequent MC labels and total number of different labels
  Double_t         fRefCov[3];              // covariance matrix of the track in the yz plane
  Double_t         fCov[3];                 // covariance matrix of the tracklet in the xy plane

  ClassDef(AliTRDseedV1, 5)                 // The offline TRD tracklet 
};

//____________________________________________________________
inline Float_t AliTRDseedV1::GetChi2Z() const
{
  Double_t dz = fZref[0]-fZfit[0]; dz*=dz;
  Double_t cov[3]; GetCovAt(fX, cov);
  Double_t s2 = fRefCov[2]+cov[2];
  return s2 > 0. ? dz/s2 : 0.; 
}

//____________________________________________________________
inline Float_t AliTRDseedV1::GetChi2Y() const
{
  Double_t dy = fYref[0]-fYfit[0]; dy*=dy;
  Double_t cov[3]; GetCovAt(fX, cov);
  Double_t s2 = fRefCov[0]+cov[0];
  return s2 > 0. ? dy/s2 : 0.; 
}

//____________________________________________________________
inline void AliTRDseedV1::Init(const AliRieman *rieman)
{
  fZref[0] = rieman->GetZat(fX0);
  fZref[1] = rieman->GetDZat(fX0);
  fYref[0] = rieman->GetYat(fX0);
  fYref[1] = rieman->GetDYat(fX0);
  fC       = rieman->GetC(); 
  fChi2    = rieman->GetChi2();
}

//____________________________________________________________
inline AliTRDcluster* AliTRDseedV1::NextCluster()
{
// Mimic the usage of STL iterators.
// Forward iterator

  fClusterIdx++; fClusterIter++;
  while(fClusterIdx < kNTimeBins){
    if(!(*fClusterIter)){ 
      fClusterIdx++; 
      fClusterIter++;
      continue;
    }
    return *fClusterIter;
  }
  return 0x0;
}

//____________________________________________________________
inline AliTRDcluster* AliTRDseedV1::PrevCluster()
{
// Mimic the usage of STL iterators.
// Backward iterator

  fClusterIdx--; fClusterIter--;
  while(fClusterIdx >= 0){
    if(!(*fClusterIter)){ 
      fClusterIdx--; 
      fClusterIter--;
      continue;
    }
    return *fClusterIter;
  }
  return 0x0;
}

//____________________________________________________________
inline void AliTRDseedV1::ResetClusterIter(Bool_t forward) 
{
// Mimic the usage of STL iterators.
// Facilitate the usage of NextCluster for forward like 
// iterator (kTRUE) and PrevCluster for backward like iterator (kFALSE)

  if(forward){
    fClusterIter = &fClusters[0]; fClusterIter--; 
    fClusterIdx=-1;
  } else {
    fClusterIter = &fClusters[kNTimeBins-1]; fClusterIter++; 
    fClusterIdx=kNTimeBins;
  }
}

#endif



