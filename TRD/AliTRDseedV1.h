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
  enum AliTRDtrackletStatus {
    kOwner    = BIT(1)
  , kRowCross = BIT(2) 
  };

  AliTRDseedV1(Int_t plane = -1);
  ~AliTRDseedV1();
  AliTRDseedV1(const AliTRDseedV1 &ref);
  AliTRDseedV1& operator=(const AliTRDseedV1 &ref);

  Bool_t	AttachClustersIter(AliTRDtrackingChamber *chamber, Float_t quality, Bool_t kZcorr = kFALSE
                                , AliTRDcluster *c=0x0);
  Bool_t	AttachClusters(AliTRDtrackingChamber *chamber, Bool_t kZcorr = kFALSE);
  void    CookdEdx(Int_t nslices);
  void    Draw(Option_t* o = "");
  Bool_t  Fit();

  Bool_t  Init(AliTRDtrackV1 *track);
  inline void      Init(const AliRieman *fit);
  Bool_t    IsOwner() const          { return TestBit(kOwner);}
  Bool_t    IsRowCross() const       { return TestBit(kRowCross);}

  inline Float_t   GetChi2Z(const Float_t z = 999.) const;
  inline Float_t   GetChi2Y(const Float_t y = 999.) const;
  void      GetCovAt(Double_t x, Double_t *cov) const;
  Double_t* GetCrossXYZ() { return &fCross[0];}
  Double_t  GetCrossSz2() const { return fCross[3];}
  Float_t*  GetdEdx() {return &fdEdx[0];}
  Float_t   GetdQdl(Int_t ic) const;
  Double_t  GetMomentum() const {return fMom;}
  Int_t     GetN() const {return fN2;}
  Float_t   GetQuality(Bool_t kZcorr) const;
  Int_t     GetPlane() const         { return fPlane;    }
  Double_t* GetProbability();
  Double_t  GetSnp() const           { return fSnp;}
  Double_t  GetTgl() const           { return fTgl;}
  Double_t  GetYat(Double_t x) const { return fYfitR[0] + fYfitR[1] * (x - fX0);}
  Double_t  GetZat(Double_t x) const { return fZfitR[0] + fZfitR[1] * (x - fX0);}
  
  void      Print(Option_t *o = "") const;
  
  void      SetMomentum(Double_t mom) {fMom = mom;}
  void      SetOwner(Bool_t own = kTRUE);
  void      SetPlane(Int_t p)                      { fPlane     = p;   }
  void      SetSnp(Double_t snp) {fSnp = snp;}
  void      SetTgl(Double_t tgl) {fTgl = tgl;}
  void      SetReconstructor(const AliTRDReconstructor *rec) {fReconstructor = rec;}
protected:

  void Copy(TObject &ref) const;

private:
  const AliTRDReconstructor *fReconstructor;
  Int_t            fPlane;                  //  TRD plane
  Float_t          fMom;                    //  Momentum estimate for tracklet [GeV/c]
  Float_t          fSnp;                    // sin of track with respect to x direction in XY plane	
  Float_t          fTgl;                    // tg of track with respect to x direction in XZ plane 	
  Float_t          fdX;                     // length of time bin
  Float_t          fdEdx[knSlices];         //  dE/dx measurements for tracklet
  Double_t         fCross[4];            // spatial parameters of the pad row crossing
  Double_t         fProb[AliPID::kSPECIES]; //  PID probabilities

  ClassDef(AliTRDseedV1, 1)                 //  New TRD seed 

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
}

#endif


