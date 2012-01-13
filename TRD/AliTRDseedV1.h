#ifndef ALITRDSEEDV1_H
#define ALITRDSEEDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// \class AliTRDseedV1
// \brief The TRD offline tracklet
// \author Alexandru Bercuci
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


#ifndef ALITRDCLUSTER_H 
#include "AliTRDcluster.h"
#endif


class TTreeSRedirector;
class TLinearFitter;

class AliRieman;

class AliTRDReconstructor;
class AliTRDtrackingChamber;
class AliTRDtrackV1;
class AliTRDpadPlane;
class AliTRDseedV1 : public AliTRDtrackletBase
{
  friend class AliHLTTRDTracklet; // wrapper for HLT

public:
  enum ETRDtrackletBuffers {    
    kNbits     = 6      // bits to store number of clusters
   ,kMask      = 0x3f   // bit mask
   ,kNtb       = 31     // max clusters/pad row
   ,kNclusters = 2*kNtb // max number of clusters/tracklet
   ,kNslices   = 10     // max dEdx slices
  };

  // bits from 0-13 are reserved by ROOT (see TObject.h)
  enum ETRDtrackletStatus {
    kOwner      = BIT(14) // owner of its clusters
   ,kRowCross   = BIT(15) // pad row cross tracklet
   ,kPID        = BIT(16) // PID contributor
   ,kCalib      = BIT(17) // calibrated tracklet
   ,kKink       = BIT(18) // kink prolongation tracklet
   ,kStandAlone = BIT(19) // tracklet build during stand alone track finding
   ,kPrimary    = BIT(20) // tracklet from a primary track candidate
   ,kChmbGood   = BIT(21) // status of the detector from calibration view point
  };

  enum ETRDtrackletError { // up to 8 bits
    kAttachClFound = 0  // not enough clusters found
    ,kAttachRowGap  = 1  // found gap attached rows
    ,kAttachRow     = 2  // found 3 rows
    ,kAttachMultipleCl= 3// multiple clusters attached to time bin
    ,kAttachClAttach= 4  // not enough clusters attached
    ,kFitCl         = 5  // not enough clusters for fit
    ,kFitFailedY    = 6  // fit failed in XY plane failed
    ,kFitFailedZ    = 7  // fit in the QZ plane failed
  };

  AliTRDseedV1(Int_t det = -1);
  ~AliTRDseedV1();
  AliTRDseedV1(const AliTRDseedV1 &ref);
  AliTRDseedV1& operator=(const AliTRDseedV1 &ref);

  Bool_t    AttachClusters(AliTRDtrackingChamber *const chamber, Bool_t tilt = kFALSE, Bool_t ChgPlus=kTRUE, Int_t ev=-1);
  void      Bootstrap(const AliTRDReconstructor *rec);
  void      Calibrate();
  void      CookdEdx(Int_t nslices);
  void      CookLabels();
  Bool_t    CookPID();
  Bool_t    Fit(UChar_t opt=0);
  Bool_t    FitRobust(Bool_t ChgPlus=kTRUE);
  Bool_t    Init(AliTRDtrackV1 *track);
  void      Init(const AliRieman *fit);
  Bool_t    IsEqual(const TObject *inTracklet) const;
  Bool_t    IsCalibrated() const     { return TestBit(kCalib);}
  Bool_t    IsChmbGood() const       { return TestBit(kChmbGood);}
  Bool_t    IsOwner() const          { return TestBit(kOwner);}
  Bool_t    IsKink() const           { return TestBit(kKink);}
  Bool_t    IsPrimary() const        { return TestBit(kPrimary);}
  Bool_t    HasPID() const           { return TestBit(kPID);}
  Bool_t    HasError(ETRDtrackletError err) const
                                     { return TESTBIT(fErrorMsg, err);}
  Bool_t    IsOK() const             { return GetN() > 4 && GetNUsed() < 4;}
  Bool_t    IsRowCross() const       { return TestBit(kRowCross);}
  Bool_t    IsUsable(Int_t i) const  { return fClusters[i] && !fClusters[i]->IsUsed();}
  Bool_t    IsStandAlone() const     { return TestBit(kStandAlone);}

  Float_t   GetAnodeWireOffset(Float_t zt);
  Float_t   GetC(Int_t typ=0) const    { return fC[typ]; }
  Float_t   GetCharge(Bool_t useOutliers=kFALSE) const;
  Float_t   GetChi2() const          { return fChi2; }
  inline Float_t   GetChi2Z() const;
  inline Float_t   GetChi2Y() const;
  inline Float_t   GetChi2Phi() const;
  void      GetCovAt(Double_t x, Double_t *cov) const;
  void      GetCovXY(Double_t *cov) const { memcpy(cov, &fCov[0], 3*sizeof(Double_t));}
  void      GetCovRef(Double_t *cov) const { memcpy(cov, &fRefCov, 7*sizeof(Double_t));}
  static Int_t GetCovSqrt(const Double_t * const c, Double_t *d);
  static Double_t GetCovInv(const Double_t * const c, Double_t *d);
  UChar_t   GetErrorMsg() const      { return fErrorMsg;}
  Float_t   GetdX() const            { return fdX;}
  const Float_t*  GetdEdx() const    { return &fdEdx[0];}
  Float_t   GetQperTB(Int_t tb) const;
  Float_t   GetdQdl() const;
  Float_t   GetdQdl(Int_t ic, Float_t *dx=NULL) const;
  Float_t   GetdYdX() const          { return fYfit[1];}
  Float_t   GetdZdX() const          { return fZfit[1];}
  Int_t     GetdY() const            { return Int_t(GetY()/0.014);}
  Int_t     GetDetector() const      { return fDet;}
  Int_t     GetChargeGaps(Float_t sz[kNtb], Float_t pos[kNtb], Int_t ntb[kNtb]) const;
  void      GetCalibParam(Float_t &exb, Float_t &vd, Float_t &t0, Float_t &s2, Float_t &dl, Float_t &dt) const    { 
              exb = fExB; vd = fVD; t0 = fT0; s2 = fS2PRF; dl = fDiffL; dt = fDiffT;}
  AliTRDcluster*  GetClusters(Int_t i) const               { return i<0 || i>=kNclusters ? NULL: fClusters[i];}
  Bool_t    GetEstimatedCrossPoint(Float_t &x, Float_t &z) const;
  Int_t     GetIndexes(Int_t i) const{ return i<0 || i>=kNclusters ? -1 : fIndexes[i];}
  Int_t     GetLabels(Int_t i) const { return fLabels[i];}  
  Float_t   GetMomentum(Float_t *err = NULL) const;
  Int_t     GetN() const             { return (Int_t)fN&kMask;}
  Int_t     GetN2() const            { return GetN();}
  Int_t     GetNUsed() const         { return Int_t((fN>>kNbits)&kMask);}
  Int_t     GetNShared() const       { return Int_t(((fN>>kNbits)>>kNbits)&kMask);}
  Int_t     GetTBoccupancy() const;
  Int_t     GetTBcross() const;
  Float_t   GetQuality(Bool_t kZcorr) const;
  Float_t   GetPadLength() const     { return fPad[0];}
  Float_t   GetPadWidth() const      { return fPad[1];}
  Int_t     GetPlane() const         { return AliTRDgeometry::GetLayer(fDet);    }

  Float_t*  GetProbability(Bool_t force=kFALSE);
  Float_t   GetPt() const            { return fPt; }
  inline Double_t  GetPID(Int_t is=-1) const;
  Float_t   GetS2Y() const           { return fS2Y;}
  Float_t   GetS2Z() const           { return fS2Z;}
  Float_t   GetSigmaY() const        { return fS2Y > 0. ? TMath::Sqrt(fS2Y) : 0.2;}
  Float_t   GetSnp() const           { return fYref[1]/TMath::Sqrt(1+fYref[1]*fYref[1]);}
  Float_t   GetTgl() const           { return fZref[1]/TMath::Sqrt(1+fYref[1]*fYref[1]);}
  Float_t   GetTilt() const          { return fPad[2];}
  UInt_t    GetTrackletWord() const  { return 0;}
  UShort_t  GetVolumeId() const;
  Float_t   GetX0() const            { return fX0;}
  Float_t   GetX() const             { return fX0 - fX;}
  Float_t   GetY() const             { return fYfit[0] - fYfit[1] * fX;}
  Double_t  GetYat(Double_t x) const { return fYfit[0] - fYfit[1] * (fX0-x);}
  Float_t   GetYfit(Int_t id) const  { return fYfit[id];}
  Float_t   GetYref(Int_t id) const  { return fYref[id];}
  Float_t   GetZ() const             { return fZfit[0] - fZfit[1] * fX;}
  Double_t  GetZat(Double_t x) const { return fZfit[0] - fZfit[1] * (fX0-x);}
  Float_t   GetZfit(Int_t id) const  { return fZfit[id];}
  Float_t   GetZref(Int_t id) const  { return fZref[id];}
  Int_t     GetYbin() const          { return Int_t(GetY()/0.016);}
  Int_t     GetZbin() const          { return Int_t(GetZ()/fPad[0]);}

  inline AliTRDcluster* NextCluster();
  inline AliTRDcluster* PrevCluster();
  void      Print(Option_t *o = "") const;
  inline void ResetClusterIter(Bool_t forward = kTRUE);
  void      Reset(Option_t *opt="");

  void      SetC(Float_t c, Int_t typ=0) { fC[typ] = c;}
  void      SetChmbGood(Bool_t k = kTRUE){ SetBit(kChmbGood, k);}
  void      SetChi2(Float_t chi2)    { fChi2 = chi2;}
  inline void SetCovRef(const Double_t *cov);
  void      SetErrorMsg(ETRDtrackletError err)  { SETBIT(fErrorMsg, err);}
  void      SetIndexes(Int_t i, Int_t idx) { fIndexes[i]  = idx; }
  void      SetLabels(Int_t *lbls)   { memcpy(fLabels, lbls, 3*sizeof(Int_t)); }
  void      SetKink(Bool_t k = kTRUE){ SetBit(kKink, k);}
  void      SetPrimary(Bool_t k = kTRUE){ SetBit(kPrimary, k);}  
  void      SetPID(Bool_t k = kTRUE) { SetBit(kPID, k);}
  void      SetStandAlone(Bool_t st) { SetBit(kStandAlone, st); }
  void      SetPt(Double_t pt)       { fPt = pt;}
  void      SetOwner();
  void      SetPadPlane(AliTRDpadPlane * const p);
  void      SetPadLength(Float_t l)  { fPad[0] = l;}
  void      SetPadWidth(Float_t w)   { fPad[1] = w;}
  void      SetTilt(Float_t tilt)    { fPad[2] = tilt; }
  void      SetDetector(Int_t d)     { fDet = d;  }
  void      SetDX(Float_t inDX)      { fdX = inDX;}
  void      SetReconstructor(const AliTRDReconstructor *rec) {fkReconstructor = rec;}
  void      SetX0(Float_t x0)        { fX0 = x0; }
  void      SetYref(Int_t i, Float_t y) { fYref[i]     = y;}
  void      SetZref(Int_t i, Float_t z) { fZref[i]     = z;}
//   void      SetUsabilityMap(Long_t um)  { fUsable = um; }
  void      Update(const AliTRDtrackV1* trk);
  void      UpdateUsed();
  void      UseClusters();

protected:
  void        Copy(TObject &ref) const;

private:
  inline void SetN(Int_t n);
  inline void SetNUsed(Int_t n);
  inline void SetNShared(Int_t n);
  inline void Swap(Int_t &n1, Int_t &n2) const;
  inline void Swap(Double_t &d1, Double_t &d2) const;

  const AliTRDReconstructor *fkReconstructor;//! local reconstructor
  AliTRDcluster  **fClusterIter;            //! clusters iterator
  Int_t            fIndexes[kNclusters];    //! Indexes
  Float_t          fExB;                    // tg(a_L) @ tracklet location
  Float_t          fVD;                     // drift velocity @ tracklet location
  Float_t          fT0;                     // time 0 @ tracklet location
  Float_t          fS2PRF;                  // sigma^2 PRF for xd->0 and phi=a_L 
  Float_t          fDiffL;                  // longitudinal diffusion coefficient
  Float_t          fDiffT;                  // transversal diffusion coefficient
  Char_t           fClusterIdx;             //! clusters iterator
  UChar_t          fErrorMsg;               // processing error
  UInt_t           fN;                      // number of clusters attached/used/shared
  Short_t          fDet;                    // TRD detector
  AliTRDcluster   *fClusters[kNclusters];   // Clusters
  Float_t          fPad[4];                 // local pad definition : length/width/tilt/anode wire offset 
  Float_t          fYref[2];                //  Reference y, dydx
  Float_t          fZref[2];                //  Reference z, dz/dx
  Float_t          fYfit[2];                //  Fit y, dy/dx
  Float_t          fZfit[2];                //  Fit z
  Float_t          fPt;                     //  Pt estimate @ tracklet [GeV/c]
  Float_t          fdX;                     // length of time bin
  Float_t          fX0;                     // anode wire position
  Float_t          fX;                      // radial position of the tracklet
  Float_t          fY;                      // r-phi position of the tracklet
  Float_t          fZ;                      // z position of the tracklet
  Float_t          fS2Y;                    // estimated resolution in the r-phi direction 
  Float_t          fS2Z;                    // estimated resolution in the z direction 
  Float_t          fC[2];                   // Curvature for standalone [0] rieman [1] vertex constrained 
  Float_t          fChi2;                   // Global chi2  
  Float_t          fdEdx[kNslices];         // dE/dx measurements for tracklet
  Float_t          fProb[AliPID::kSPECIES]; // PID probabilities
  Int_t            fLabels[3];              // most frequent MC labels and total number of different labels
  Double_t         fRefCov[7];              // covariance matrix of the track in the yz plane + the rest of the diagonal elements
  Double_t         fCov[3];                 // covariance matrix of the tracklet in the xy plane

  ClassDef(AliTRDseedV1, 12)                 // The offline TRD tracklet 
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
  return NULL;
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
  return NULL;
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
inline void AliTRDseedV1::SetCovRef(const Double_t *cov)
{ 
// Copy some "important" covariance matrix elements
//  var(y)
// cov(y,z)  var(z)
//                  var(snp)
//                           var(tgl)
//                        cov(tgl, 1/pt)  var(1/pt)

  memcpy(&fRefCov[0], cov, 3*sizeof(Double_t)); // yz full covariance
  fRefCov[3] = cov[ 5];  // snp variance 
  fRefCov[4] = cov[ 9];  // tgl variance
  fRefCov[5] = cov[13];  // cov(tgl, 1/pt)
  fRefCov[6] = cov[14];  // 1/pt variance
}


//____________________________________________________________
inline void AliTRDseedV1::SetN(Int_t n)
{
  if(n<0 || n>kNclusters) return; 
  fN &= ~kMask; 
  fN |= (n&kMask);
}

//____________________________________________________________
inline void AliTRDseedV1::SetNUsed(Int_t n)
{
  if(n<0 || n>kNclusters) return; 
  UInt_t mask(kMask<<kNbits); 
  fN &= ~mask;
  n=n<<kNbits; fN |= (n&mask);
}

//____________________________________________________________
inline void AliTRDseedV1::SetNShared(Int_t n)
{
  if(n<0 || n>kNclusters) return; 
  UInt_t mask((kMask<<kNbits)<<kNbits); 
  fN &= ~mask;
  n = (n<<kNbits)<<kNbits; fN|=(n&mask);
}

//____________________________________________________________
inline void AliTRDseedV1::Swap(Int_t &n1, Int_t &n2) const
{
// swap values of n1 with n2
  Int_t tmp(n1);
  n1=n2; n2=tmp;
}

//____________________________________________________________
inline void AliTRDseedV1::Swap(Double_t &d1, Double_t &d2) const
{
// swap values of d1 with d2
  Double_t tmp(d1);
  d1=d2; d2=tmp;
}


#endif



