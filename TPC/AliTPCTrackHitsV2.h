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
class AliTPChit;
class AliTPCTempHitInfoV2;
class AliTPCCurrentHitV2;
class AliHit;

class AliTrackHitsParamV2 : public TObject {


public:
  AliTrackHitsParamV2();

   AliTrackHitsParamV2(const AliTrackHitsParamV2 &hit):  TObject(hit),
   fTrackID(0), 
   fVolumeID(0),
   fR(0.),  
   fZ(0.),  
   fFi(0.), 
   fAn(0.), 
   fAd(0.), 
   fTheta(0.), 
   fThetaD(0.), 
   fNHits(0), 
   fHitDistance(0), 
   fCharge(0),
   fTime(0) 
    {hit.Copy(*this);}
  AliTrackHitsParamV2& operator = (const AliTrackHitsParamV2 &hit)
     {hit.Copy(*this); return (*this);}
  ~AliTrackHitsParamV2();

  Int_t   GetTrackID()            const {return fTrackID;}
  Int_t   GetVolumeID()           const {return fVolumeID;}
  Float_t GetR()                  const {return fR;}
  Float_t GetZ()                  const {return fZ;}
  Float_t GetFi()                 const {return fFi;}
  Float_t GetAn()                 const {return fAn;}
  Float_t GetAd()                 const {return fAd;}
  Float_t GetTheta()              const {return fTheta;}
  Float_t GetThetaD()             const {return fThetaD;}
  Int_t   GetNHits()              const {return fNHits;}

  Short_t HitDistance(Int_t i) const {return fHitDistance[i];}
  Short_t Charge(Int_t i)      const {return fCharge[i];}
  Short_t Time(Int_t i)        const {return fTime[i];}

  Short_t& HitDistance(Int_t i) {return fHitDistance[i];}
  Short_t& Charge(Int_t i)      {return fCharge[i];}
  Short_t& Time(Int_t i)        {return fTime[i];}

  void SetHitDistance(Int_t i)
    {Short_t *s=new Short_t[i];
    delete [] fHitDistance; fHitDistance=s;}

  void SetCharge(Int_t i)
    {Short_t *s=new Short_t[i];
    delete [] fCharge; fCharge=s;}

  void SetTime(Int_t i)
    {Short_t *s=new Short_t[i];
    delete [] fTime; fTime=s;}

  void ResizeHitDistance(Int_t i)
    {Short_t *s=new Short_t[i];
    memcpy(s, fHitDistance, sizeof(Short_t)*i); 
    delete [] fHitDistance; fHitDistance=s;}

  void ResizeCharge(Int_t i)
    {Short_t *s=new Short_t[i];
    memcpy(s, fCharge, sizeof(Short_t)*i); 
    delete [] fCharge; fCharge=s;}

  void ResizeTime(Int_t i)
    {Short_t *s=new Short_t[i];
    memcpy(s, fTime, sizeof(Short_t)*i); 
    delete [] fTime; fTime=s;}

  void SetTrackID(Int_t id)    {fTrackID=id;}
  void SetVolumeID(Short_t id) {fVolumeID=id;}
  void SetR(Float_t r)         {fR=r;}
  void SetZ(Float_t z)         {fZ=z;}
  void SetFi(Float_t fi)       {fFi=fi;}
  void SetAn(Float_t an)       {fAn=an;}
  void SetAd(Float_t ad)       {fAd=ad;}
  void SetTheta(Float_t t)     {fTheta=t;}
  void SetThetaD(Float_t t)    {fThetaD=t;}
  void SetNHits(Int_t n)       {fNHits=n;}

  Float_t Eta() const;

 private:
  Int_t fTrackID; // ID of the trac©k
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

public:
  AliTPCTrackHitsV2();
  ~AliTPCTrackHitsV2();
  AliTPCTrackHitsV2(const AliTPCTrackHitsV2 &hit):  TObject(hit),
  fArray(0),  
  fSize(0),           
  fPrecision(0.),  
  fStep(0.),       
  fMaxDistance(0),   
  fNVolumes(0),        
  fVolumes(0),   
  fTempInfo(0), 
  fCurrentHit(0),  
  fHit(0) 
    {hit.Copy(*this);}
  AliTPCTrackHitsV2& operator = (const AliTPCTrackHitsV2 &hit)
     {hit.Copy(*this); return (*this);}
  void Clear(Option_t * /*option*/ ="");
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
  static Double_t GetKPrecision()  {return fgkPrecision;}
  static Double_t GetKPrecision2() {return fgkPrecision2;}

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

class AliTPCCurrentHitV2 {
public:
  Int_t    GetStackIndex() const {return fStackIndex;}
  void     SetStackIndex(Int_t i) {fStackIndex=i;}
  Int_t    GetParamIndex() const {return fParamIndex;}
  void     SetParamIndex(Int_t i) {fParamIndex=i;}
  Double_t GetR() const {return fR;}
  void     SetR(Double_t r) {fR=r;}
  Bool_t   GetStatus() const {return fStatus;}
  void     SetStatus(Bool_t s) {fStatus=s;}
private:
  Int_t   fParamIndex;//  - current param pointer
  Int_t   fStackIndex; // - current hit stack index
  Double_t fR;   //current Radius
  Bool_t  fStatus; //current status    
};   



#endif //ALITPCTRACKHITSV2_H
