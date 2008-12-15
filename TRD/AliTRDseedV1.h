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

#ifndef ALITRDSEED_H
#include "AliTRDseed.h"
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
class AliTRDseedV1 : public AliTRDseed
{

  public:

  enum {
    knSlices = 10
  };
  // bits from 0-13 are reserved by ROOT (see TObject.h)
  enum AliTRDtrackletStatus {
    kOwner    = BIT(14)
  , kRowCross = BIT(15) 
  };

  AliTRDseedV1(Int_t det = -1);
  ~AliTRDseedV1();
  AliTRDseedV1(const AliTRDseedV1 &ref);
  AliTRDseedV1& operator=(const AliTRDseedV1 &ref);

  Bool_t	  AttachClustersIter(
              AliTRDtrackingChamber *chamber, Float_t quality, 
              Bool_t kZcorr = kFALSE, AliTRDcluster *c=0x0);
  Bool_t	  AttachClusters(
              AliTRDtrackingChamber *chamber, Bool_t kZcorr = kFALSE);
  void      Bootstrap(const AliTRDReconstructor *rec);
  void      CookdEdx(Int_t nslices);
  Bool_t    Fit(Bool_t tilt=kTRUE);

  Bool_t    Init(AliTRDtrackV1 *track);
  inline void      Init(const AliRieman *fit);
  Bool_t    IsEqual(const TObject *inTracklet) const;
  Bool_t    IsOwner() const          { return TestBit(kOwner);}
  Bool_t    IsRowCross() const       { return TestBit(kRowCross);}

  inline Float_t   GetChi2Z(const Float_t z = 999.) const;
  inline Float_t   GetChi2Y(const Float_t y = 999.) const;
  void      GetCovAt(Double_t x, Double_t *cov) const;
  void      GetCovRef(Double_t *cov) const { memcpy(cov, &fRefCov[0], 3*sizeof(Double_t));}
  Double_t* GetCrossXYZ()            { return &fCross[0];}
  Double_t  GetCrossSz2() const      { return fCross[3];}
  Float_t   GetdX() const            { return fdX;}
  Float_t*  GetdEdx()                { return &fdEdx[0];}
  Float_t   GetdQdl(Int_t ic) const;
  Int_t     GetDetector() const      { return fDet;}
  Double_t  GetMomentum() const      { return fMom;}
  Int_t     GetN() const             { return fN2;}
  Float_t   GetQuality(Bool_t kZcorr) const;
  Int_t     GetPlane() const         { return AliTRDgeometry::GetLayer(fDet);    }

  Double_t* GetProbability();
  Double_t  GetSnp() const           { return fSnp;}
  Double_t  GetTgl() const           { return fTgl;}
  Float_t   GetXref() const          { return fXref;}
  Double_t  GetYat(Double_t x) const { return fYfit[0] + fYfit[1] * (fX0-x);}
  Double_t  GetZat(Double_t x) const { return fZfit[0] + fZfit[1] * (fX0-x);}
  
  inline AliTRDcluster* NextCluster();
  inline AliTRDcluster* PrevCluster();
  void      Print(Option_t *o = "") const;
  inline void ResetClusterIter(Bool_t forward = kTRUE);

  void      SetCovRef(const Double_t *cov) { memcpy(&fRefCov[0], cov, 3*sizeof(Double_t));}
  void      SetMomentum(Double_t mom){ fMom = mom;}
  void      SetOwner();
  void      SetDetector(Int_t d)     { fDet = d;  }
  void      SetDX(Float_t inDX)      { fdX = inDX;}
  void      SetSnp(Double_t snp)     { fSnp = snp;}
  void      SetTgl(Double_t tgl)     { fTgl = tgl;}
  void      SetReconstructor(const AliTRDReconstructor *rec) {fReconstructor = rec;}
protected:

  void Copy(TObject &ref) const;

private:
  const AliTRDReconstructor *fReconstructor;//! local reconstructor
  AliTRDcluster    **fClusterIter;          //! clusters iterator
  Char_t           fClusterIdx;             //! clusters iterator
  Int_t            fDet;                    //  TRD detector
  Float_t          fMom;                    //  Momentum estimate for tracklet [GeV/c]
  Float_t          fSnp;                    // sin of track with respect to x direction in XY plane	
  Float_t          fTgl;                    // tg of track with respect to x direction in XZ plane 	
  Float_t          fdX;                     // length of time bin
  Float_t          fXref;                   // average radial position of clusters
  Float_t          fdEdx[knSlices];         // dE/dx measurements for tracklet
  Double_t         fCross[4];               // spatial parameters of the pad row crossing
  Double_t         fRefCov[3];              // covariance matrix of the track in the yz plane
  Double_t         fProb[AliPID::kSPECIES]; //  PID probabilities

  ClassDef(AliTRDseedV1, 3)                 //  New TRD seed 

};

//____________________________________________________________
inline Float_t AliTRDseedV1::GetChi2Z(const Float_t z) const
{
  Float_t z1  = (z == 999.) ? fMeanz : z;
  Float_t chi = fZref[0] - z1;
  return chi*chi;
}

//____________________________________________________________
inline Float_t AliTRDseedV1::GetChi2Y(const Float_t y) const
{
  Float_t y1  = (y == 999.) ? fYfitR[0] : y;
  Float_t chi = fYref[0] - y1;
  return chi*chi;
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
  while(fClusterIdx < AliTRDseed::knTimebins){
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
    fClusterIter = &fClusters[AliTRDseed::knTimebins-1]; fClusterIter++; 
    fClusterIdx=AliTRDseed::knTimebins;
  }
}

#endif



