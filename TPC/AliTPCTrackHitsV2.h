#ifndef ALITPCTRACKHITSV2_H
#define ALITPCTRACKHITSV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////
//  Manager class for TPC   hits                   //
////////////////////////////////////////////////
//

#include "TObject.h"

class TClonesArray;
class AliArrayS;
class AliTPChit;
class AliTPCTempHitInfoV2;
class AliTPCCurrentHitV2;
class AliHit;

class AliTrackHitsParamV2 : public TObject {
  friend class   AliTPC;
  friend class   AliTRD;
  friend class   AliTPCTrackHitsV2;
  friend class   AliTPCTempHitInfoV2;
  friend class   AliTRDtrackHits;

public:
  AliTrackHitsParamV2();
  AliTrackHitsParamV2(const AliTrackHitsParamV2 &hit):  TObject(hit)
    {hit.Copy(*this);}
  AliTrackHitsParamV2& operator = (const AliTrackHitsParamV2 &hit)
     {hit.Copy(*this); return (*this);}
  ~AliTrackHitsParamV2();

 private:
  Int_t fTrackID; // ID of the track
  Short_t fVolumeID;// volume ID
  Float_t fR;  //radius
  Float_t fZ;  //z position
  Float_t fFi; //radial angle
  Float_t fAn; //angle with  the radial vector
  Float_t fAd; //derivation of angle
  Float_t fTheta; //theta angle
  Float_t fThetaD; //theta angle derivation
  Int_t   fNHits; //nuber of thits
  Short_t * fHitDistance; //[fNHits] array of hits distances
  Short_t * fCharge; //[fNHits] array of charges
  Short_t * fTime; //[fNHits] array of hits time
  static Int_t fgCounter1; //First internal counter
  static Int_t fgCounter2; // Second internal counter

  void Copy(TObject &) const
  {Error("Copy","Not Implemented");}

  ClassDef(AliTrackHitsParamV2,2)  
};


class AliTPCTrackHitsV2 : public TObject {
  friend class AliTPCTempHitInfoV2;

public:
  AliTPCTrackHitsV2(); 
  ~AliTPCTrackHitsV2();
  AliTPCTrackHitsV2(const AliTPCTrackHitsV2 &hit):  TObject(hit)
    {hit.Copy(*this);}
  AliTPCTrackHitsV2& operator = (const AliTPCTrackHitsV2 &hit)
     {hit.Copy(*this); return (*this);}
  void Clear();
  void AddHitKartez(Int_t volumeID, Int_t trackID, Double_t x, 
		    Double_t y, Double_t z,Int_t q,Float_t time);
  void AddHit(Int_t volumeID, Int_t trackID, Double_t r, 
	      Double_t z, Double_t fi,Int_t q,Float_t time);
 
  Bool_t First(); //set current hit to first hit 
  Bool_t Next();  //set current hit to next
  AliHit * GetHit() const;
  AliTrackHitsParamV2 * GetParam();

  TClonesArray * GetArray(){return fArray;}
  Int_t  GetEntriesFast() const { return fSize;}
  void SetHitPrecision(Double_t prec) {fPrecision=prec;}
  void SetStepPrecision(Double_t prec) {fStep=prec;}
  void SetMaxDistance(UInt_t distance) {fMaxDistance = distance;}
  Bool_t  FlushHitStack(Bool_t force=kTRUE);    //
  Int_t *  GetVolumes(){ return fVolumes;}
  Int_t GetNVolumes() const {return fNVolumes;}

public:
  void AddVolume(Int_t volume); //add volumes to tthe list of volumes
  void FlushHitStack2(Int_t index1, Int_t index2);   //

protected:
  TClonesArray * fArray;  //array of compressed hits
  Int_t fSize;            //total number of hits in track
  Double_t fPrecision;  // required precision
  Double_t fStep;       //unit step size
  UInt_t fMaxDistance;   //maximal distance between two connected hits 
  Int_t fNVolumes;      //number of volumes in track  
  Int_t *  fVolumes;    //[fNVolumes] list of volumes
  AliTPCTempHitInfoV2 * fTempInfo; //!information about track
  AliTPCCurrentHitV2  * fCurrentHit; //!information about current hit 
  AliHit * fHit;                     //! current hit information
  static const Double_t fgkPrecision;  //precision 
  static const Double_t fgkPrecision2;  //precision
  static const Double_t fgkTimePrecision;  //hit time precision 
  static Int_t fgCounter1; // First internal counter
  static Int_t fgCounter2; // Second internal counter

private:
  void Copy(TObject &) const
  {Error("Copy","Not Implemented");}


  ClassDef(AliTPCTrackHitsV2,2) 
};

struct AliTPCCurrentHitV2 {
  Int_t   fParamIndex;//  - current param pointer
  Int_t   fStackIndex; // - current hit stack index
  Double_t fR;   //current Radius
  Bool_t  fStatus; //current status    
};   



#endif //ALITPCTRACKHITSV2_H
