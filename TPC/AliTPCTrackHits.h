#ifndef ALITPCTRACKHITS_H
#define ALITPCTRACKHITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////
//  Manager class for TPC   clusters                   //
////////////////////////////////////////////////

#include "AliCTypes.h"
#include "AliSegmentID.h"
#include "AliArrayS.h"
#include "AliTPC.h"
#include "TVector3.h"
#include "AliObjectArray.h"
#include "TArrayOfArray.h"

class TClonesArray;
class AliArrayS;
class AliTPChit;
class AliTPCTempHitInfo;
class AliTPCCurrentHit;


class AliTrackHitsInfo {
public:
  AliTrackHitsInfo(){fgCounter1++;fgCounter2++;}
  ~AliTrackHitsInfo(){fgCounter1--;}
  Int_t   fTrackID;  //track ID
  Int_t   fVolumeID;   //volume ID
  UInt_t   fHitParamIndex; //corresponding index  
  static Int_t fgCounter1;
  static Int_t fgCounter2;  
  LClassDef(AliTrackHitsInfo,1)  
};


class AliTrackHitsParam {
public:
  AliTrackHitsParam(){fgCounter1++;fgCounter2++;}
  ~AliTrackHitsParam(){fgCounter1--;}
  Float_t fR;  //radius
  Float_t fZ;  //z position
  Float_t fFi; //radial angle
  Float_t fAn; //angle with  the radial vector
  Float_t fAd; //derivation of angle
  Float_t fTheta; //theta angle
  Float_t fThetaD; //theta angle derivation
  static Int_t fgCounter1;
  static Int_t fgCounter2;  
  LClassDef(AliTrackHitsParam,1)  
};


class AliHitInfo {
public:
  AliHitInfo(){fgCounter1++;fgCounter2++;}
  ~AliHitInfo(){fgCounter1--;}
  Short_t fHitDistance; //distance to previous hit
  Short_t fCharge; //deponed charge
  static Int_t fgCounter1;
  static Int_t fgCounter2;  
  LClassDef(AliHitInfo,1)
};



class AliTPCTrackHits : public TObject{
public:
  AliTPCTrackHits(); 
  ~AliTPCTrackHits();
  void Clear();
  void AddHitKartez(Int_t volumeID, Int_t trackID, Double_t x, 
		    Double_t y, Double_t z,Int_t q);
  void AddHit(Int_t volumeID, Int_t trackID, Double_t r, 
	      Double_t z, Double_t fi,Int_t q);
 
  Bool_t First(); //set current hit to first hit 
  Bool_t Next();  //set current hit to next
  AliTPChit * GetHit();
  AliTrackHitsParam * GetParam();
  AliHitInfo * GetHitInfo();
  Int_t  GetEntriesFast() { return fHitsPosAndQ ? fHitsPosAndQ->ArraySize():0;}
  void SetHitPrecision(Double_t prec) {fPrecision=prec;}
  void SetStepPrecision(Double_t prec) {fStep=prec;}
  void SetMaxDistance(UInt_t distance) {fMaxDistance = distance;}
  Bool_t  FlushHitStack(Bool_t force=kTRUE);    //
public:
  void FlushHitStack2(Int_t index1, Int_t index2);   //
  AliObjectArray * fTrackHitsInfo;  //quick information about track
  AliObjectArray * fTrackHitsParam;  //hit information  
  TArrayOfArray_vStack * fHitsPosAndQ;  //position information

  Double_t fPrecision;  // required precision
  Double_t fStep;       //unit step size
  UInt_t fMaxDistance;   //maximal distance between two connected hits 
  AliTPCTempHitInfo * fTempInfo; //!information about track
  AliTPCCurrentHit  * fCurrentHit; //!information about current hit 
  static const Double_t fgkPrecision;  //precision 
  static const Double_t fgkPrecision2;  //precision
  static Int_t fgCounter1;
  static Int_t fgCounter2;  
  ClassDef(AliTPCTrackHits,1) 
};


#endif //ALITPCTRACKHITS_H
