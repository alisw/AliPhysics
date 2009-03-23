#ifndef ALITRDSEEDV1_H
#define ALITRDSEEDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The TRD offline tracklet                                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDTRACKLETBASE_H
#include "AliTRDtrackletBase.h"
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

#ifndef ALITRDCLUSTER_H 
#include "AliTRDcluster.h"
#endif

class TTreeSRedirector;

class AliRieman;

class AliTRDtrackingChamber;
class AliTRDtrackV1;
class AliTRDReconstructor;
class AliTRDpadPlane;
class AliTRDseedV1 : public AliTRDtrackletBase
{
public:
  enum ETRDtrackletBuffers {    
    kNtb       = 31     // max clusters/pad row
   ,kNclusters = 2*kNtb // max number of clusters/tracklet
   ,kNslices   = 10     // max dEdx slices
  };

  // bits from 0-13 are reserved by ROOT (see TObject.h)
  enum ETRDtrackletStatus {
    kOwner      = BIT(14) // owner of its clusters
   ,kRowCross   = BIT(15) // pad row cross tracklet
   ,kCalib      = BIT(16) // calibrated tracklet
   ,kKink       = BIT(17) // kink prolongation tracklet
   ,kStandAlone = BIT(18)
  };

  AliTRDseedV1(Int_t det = -1);
  ~AliTRDseedV1();
  AliTRDseedV1(const AliTRDseedV1 &ref);
  AliTRDseedV1& operator=(const AliTRDseedV1 &ref);

/*  Bool_t	  AttachClustersIter(
              AliTRDtrackingChamber *chamber, Float_t quality, 
              Bool_t kZcorr = kFALSE, AliTRDcluster *c=0x0);*/
  Bool_t	  AttachClusters(
              AliTRDtrackingChamber *chamber, Bool_t tilt = kFALSE);
  void      Bootstrap(const AliTRDReconstructor *rec);
  void      Calibrate();
  void      CookdEdx(Int_t nslices);
  void      CookLabels();
  Bool_t    CookPID();
  Bool_t    Fit(Bool_t tilt=kTRUE, Int_t errors = 2);
//   void      FitMI();
  Bool_t    Init(AliTRDtrackV1 *track);
  inline void      Init(const AliRieman *fit);
  Bool_t    IsEqual(const TObject *inTracklet) const;
  Bool_t    IsCalibrated() const     { return TestBit(kCalib);}
  Bool_t    IsOwner() const          { return TestBit(kOwner);}
  Bool_t    IsKink() const           { return TestBit(kKink);}
  Bool_t    IsOK() const             { return GetN() > 4 && GetNUsed() < 4;}
  Bool_t    IsRowCross() const       { return TestBit(kRowCross);}
  Bool_t    IsUsable(Int_t i) const  { return fClusters[i] && !fClusters[i]->IsUsed();}
  Bool_t    IsStandAlone() const     { return TestBit(kStandAlone);}

  Float_t   GetC() const             { return fC; }
  Float_t   GetChi2() const          { return fChi2; }
  inline Float_t   GetChi2Z() const;
  inline Float_t   GetChi2Y() const;
  inline Float_t   GetChi2Phi() const;
  static void      GetClusterXY(const AliTRDcluster *c, Double_t &x, Double_t &y);
  void      GetCovAt(Double_t x, Double_t *cov) const;
  void      GetCovXY(Double_t *cov) const { memcpy(cov, &fCov[0], 3*sizeof(Double_t));}
  void      GetCovRef(Double_t *cov) const { memcpy(cov, &fRefCov, 3*sizeof(Double_t));}
  static Double_t GetCovSqrt(Double_t *c, Double_t *d);
  static Double_t GetCovInv(Double_t *c, Double_t *d);
  Float_t   GetdX() const            { return fdX;}
  Float_t*  GetdEdx()                { return &fdEdx[0];}
  Float_t   GetdQdl(Int_t ic) const;
  Float_t   GetdYdX() const          { return fYfit[1]; } 
  Float_t   GetdZdX() const          { return fZref[1]; }
  Int_t     GetdY() const            { return Int_t(GetY()/0.014);}
  Int_t     GetDetector() const      { return fDet;}
  void      GetCalibParam(Float_t &exb, Float_t &vd, Float_t &t0, Float_t &s2, Float_t &dl, Float_t &dt) const    { 
              exb = fExB; vd = fVD; t0 = fT0; s2 = fS2PRF; dl = fDiffL; dt = fDiffT;}
  AliTRDcluster*  GetClusters(Int_t i) const               { return i<0 || i>=kNclusters ? 0x0 : fClusters[i];}
  Int_t     GetIndexes(Int_t i) const{ return i<0 || i>=kNclusters ? -1 : fIndexes[i];}
  Int_t     GetLabels(Int_t i) const { return fLabels[i];}  
  Double_t  GetMomentum() const      { return fMom;}
  Int_t     GetN() const             { return (Int_t)fN&0x1f;}
  Int_t     GetN2() const            { return GetN();}
  Int_t     GetNUsed() const         { return Int_t((fN>>5)&0x1f);}
  Int_t     GetNShared() const       { return Int_t((fN>>10)&0x1f);}
  Float_t   GetQuality(Bool_t kZcorr) const;
  Float_t   GetPadLength() const     { return fPad[0];}
  Float_t   GetPadWidth() const      { return fPad[1];}
  Int_t     GetPlane() const         { return AliTRDgeometry::GetLayer(fDet);    }

  Float_t*  GetProbability(Bool_t force=kFALSE);
  inline Double_t  GetPID(Int_t is=-1) const;
  Float_t   GetS2Y() const           { return fS2Y;}
  Float_t   GetS2Z() const           { return fS2Z;}
  Float_t   GetSigmaY() const        { return fS2Y > 0. ? TMath::Sqrt(fS2Y) : 0.2;}
  Float_t   GetSnp() const           { return fYref[1]/TMath::Sqrt(1+fYref[1]*fYref[1]);}
  Float_t   GetTgl() const           { return fZref[1];}
  Float_t   GetTilt() const          { return fPad[2];}
  UInt_t    GetTrackletWord() const  { return 0;}
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
  Int_t     GetYbin() const         { return Int_t(GetY()/0.016);}
  Int_t     GetZbin() const         { return Int_t(GetZ()/fPad[0]);}

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
  void      SetKink(Bool_t k)        { SetBit(kKink, k);}
  void      SetStandAlone(Bool_t st) { SetBit(kStandAlone, st); }
  void      SetMomentum(Double_t mom){ fMom = mom;}
  void      SetOwner();
  void      SetPadPlane(AliTRDpadPlane *p);
  void      SetPadLength(Float_t l)  { fPad[0] = l;}
  void      SetPadWidth(Float_t w)   { fPad[1] = w;}
  void      SetTilt(Float_t tilt)    { fPad[2] = tilt; }
  void      SetDetector(Int_t d)     { fDet = d;  }
  void      SetDX(Float_t inDX)      { fdX = inDX;}
  void      SetReconstructor(const AliTRDReconstructor *rec) {fReconstructor = rec;}
  void      SetX0(Float_t x0)        { fX0 = x0; }
  void      SetYref(Int_t i, Float_t y) { fYref[i]     = y;}
  void      SetZref(Int_t i, Float_t z) { fZref[i]     = z;}
//   void      SetUsabilityMap(Long_t um)  { fUsable = um; }
  void      UpDate(const AliTRDtrackV1* trk);
  void      UpdateUsed();
  void      UseClusters();

protected:
  void        Copy(TObject &ref) const;

private:
  inline void SetN(Int_t n);
  inline void SetNUsed(Int_t n);
  inline void SetNShared(Int_t n);

  const AliTRDReconstructor *fReconstructor;//! local reconstructor
  AliTRDcluster  **fClusterIter;            //! clusters iterator
  Int_t            fIndexes[kNclusters];    //! Indexes
  Float_t          fExB;                    //! tg(a_L) @ tracklet location
  Float_t          fVD;                     //! drift velocity @ tracklet location
  Float_t          fT0;                     //! time 0 @ tracklet location
  Float_t          fS2PRF;                  //! sigma^2 PRF for xd->0 and phi=a_L 
  Float_t          fDiffL;                  //! longitudinal diffusion coefficient
  Float_t          fDiffT;                  //! transversal diffusion coefficient
  Char_t           fClusterIdx;             //! clusters iterator
  UShort_t         fN;                      // number of clusters attached/used/shared
  Short_t          fDet;                    // TRD detector
  AliTRDcluster   *fClusters[kNclusters];   // Clusters
  Float_t          fPad[3];                 // local pad definition : length/width/tilt 
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
  Float_t          fdEdx[kNslices];         // dE/dx measurements for tracklet
  Float_t          fProb[AliPID::kSPECIES]; //  PID probabilities
  Int_t            fLabels[3];              // most frequent MC labels and total number of different labels
  Double_t         fRefCov[3];              // covariance matrix of the track in the yz plane
  Double_t         fCov[3];                 // covariance matrix of the tracklet in the xy plane

  ClassDef(AliTRDseedV1, 6)                 // The offline TRD tracklet 
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
inline Float_t AliTRDseedV1::GetChi2Phi() const
{
  Double_t dphi = fYref[1]-fYfit[1]; dphi*=dphi;
  Double_t cov[3]; GetCovAt(fX, cov);
  Double_t s2 = fRefCov[2]+cov[2];
  return s2 > 0. ? dphi/s2 : 0.; 
}

//____________________________________________________________
inline Double_t AliTRDseedV1::GetPID(Int_t is) const
{
  if(is<0) return fProb[AliPID::kElectron];
  if(is<AliPID::kSPECIES) return fProb[is];
  return 0.;
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
  while(fClusterIdx < kNclusters){
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
    fClusterIter = &fClusters[kNclusters-1]; fClusterIter++; 
    fClusterIdx=kNclusters;
  }
}

//____________________________________________________________
inline void AliTRDseedV1::SetN(Int_t n)
{
  if(n<0 || n>= (1<<5)) return; 
  fN &= ~0x1f;
  fN |= n;
}

//____________________________________________________________
inline void AliTRDseedV1::SetNUsed(Int_t n)
{
  if(n<0 || n>= (1<<5)) return; 
  fN &= ~(0x1f<<5);
  n <<= 5; fN |= n;
}

//____________________________________________________________
inline void AliTRDseedV1::SetNShared(Int_t n)
{
  if(n<0 || n>= (1<<5)) return; 
  fN &= ~(0x1f<<10);
  n <<= 10; fN |= n;
}


#endif



