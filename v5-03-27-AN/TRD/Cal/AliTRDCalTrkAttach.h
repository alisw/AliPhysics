#ifndef ALITRDCALTRKATTACH_H
#define ALITRDCALTRKATTACH_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////
//
// Container for calibration parameters for AliTRDseedV1::AttachClusters()
// For calibration procedure check AliTRDtrackleOflHelper::CalibrateAttach()
//
// Author Alex Bercuci <A.Bercuci@gsi.de>
//
//////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif
class TObjArray;
class AliTRDCalTrkAttach : public TNamed
{
public:
  enum ETRDCalTrkAttachCalib {
    kResPos = 0    // relative position residual location
   ,kResAng        // angular residual location
   ,kSigma         // relative error location
   ,kNclMean       // mean no. of clusters/tracklet location
   ,kNcalib        // no. of calib objects
  };
  enum ETRDCalTrkAttachParam {
    kNcharge = 2
  };

  AliTRDCalTrkAttach();
  virtual ~AliTRDCalTrkAttach();

  Double_t  CookLikelihood(Bool_t chg, Int_t ly, Float_t pt, Float_t phi, Int_t ncl, Double_t dy, Double_t dphi, Double_t sr) const;
  void      GetNsgmDy(Int_t &n0, Int_t &n1) const {n0 = fNsgmDy[0]; n1 = fNsgmDy[1];}
  void      GetLikeMinRelDecrease(Float_t &p0, Float_t &p1) const {p0 = fLikeMinRelDecrease[0]; p1 = fLikeMinRelDecrease[1];}
  Float_t   GetRClikeLimit() const { return fRClikeLimit;}
  Float_t   GetScaleCov() const { return fScaleCov;}
  Bool_t    LoadReferences(const Char_t *file);
  void      SetNsgmDy(Int_t ns0, Int_t ns1) {fNsgmDy[0] = ns0; fNsgmDy[0] = ns1;}
  void      SetLikeMinRelDecrease(Float_t p0, Float_t p1) {fLikeMinRelDecrease[0]=p0;fLikeMinRelDecrease[1]=p1;}
  void      SetRClikeLimit(Float_t rc) { fRClikeLimit = rc;}
  void      SetScaleCov(Float_t sc) { fScaleCov = sc;}

private:
  Int_t     fNsgmDy[2];             //
  Float_t   fLikeMinRelDecrease[2]; //
  Float_t   fRClikeLimit;           //
  Float_t   fScaleCov;              //
  TObjArray *fLike;                 // array with likelihoods
    
  AliTRDCalTrkAttach(const AliTRDCalTrkAttach& ref);
  AliTRDCalTrkAttach    &operator=(const AliTRDCalTrkAttach &rhs);

  ClassDef(AliTRDCalTrkAttach, 1)    //  Storage for AttachClusters() calibration
};

#endif


