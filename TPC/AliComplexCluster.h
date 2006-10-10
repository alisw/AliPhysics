#ifndef ALICOMPLEXCLUSTER_H
#define ALICOMPLEXCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
//
// this is a class
// AliComplexCluster
// by M. Ivanov
//

#include "TObject.h"
#include "TMath.h"
#include "AliTPCclusterMI.h"
//

class AliComplexCluster : public TObject {
public:

  AliComplexCluster() {
    fTracks[0]=fTracks[1]=fTracks[2]=0; 
    fX=fY=fQ=fSigmaX2=fSigmaY2=fSigmaXY=fArea=fMax=0.;
  }
  virtual ~AliComplexCluster() {;}
  Bool_t    IsSortable() const;
  Int_t Compare(const TObject *o) const;
  // the following getters are needed by HLT
  // please dont remove... C. Loizides
  Int_t GetTrack(Int_t i)const {return fTracks[i];} //labels of overlapped tracks
  Float_t GetX()const {return fX;}
  Float_t GetY()const {return fY;}
  Float_t GetQ()const {return fQ;}
  Float_t GetSigmaX2()const {return fSigmaX2;}
  Float_t GetSigmaY2()const {return fSigmaY2;}
  Float_t GetSigmaXY()const {return fSigmaXY;}
  Float_t GetArea()const {return fArea;}
  Float_t GetMax()const {return fMax;}
private:
  Int_t     fTracks[3];//labels of overlapped tracks
  Float_t   fX ;       //Y of cluster
  Float_t   fY ;       //Z of cluster
  Float_t   fQ ;       //Q of cluster (in ADC counts)
  Float_t   fSigmaX2;  //Sigma Y square of cluster
  Float_t   fSigmaY2;  //Sigma Z square of cluster
  Float_t   fSigmaXY;  //      XY moment 
  Float_t   fArea;     //area of cluster
  Float_t   fMax;     //amplitude at maximum 

  ClassDef(AliComplexCluster,1)
    // Cluster manager
};



class AliTPCTrackerPoint  {  
 public:

  AliTPCTrackerPoint():
    fTX(0),
    fTZ(0),
    fTY(0),
    fTAngleZ(0),
    fTAngleY(0),
    fSigmaZ(0),
    fSigmaY(0),
    fErrZ(0),
    fErrY(0),
    fIsShared(0){}
  virtual ~AliTPCTrackerPoint(){}
  Float_t  GetX() const  {return (fTX*0.01);}
  Float_t  GetZ() const {return (fTZ*0.01);}
  Float_t  GetY() const {return (fTY*0.01);}
  Float_t  GetAngleZ() const  {return (Float_t(fTAngleZ)*0.02);}
  Float_t  GetAngleY() const {return (Float_t(fTAngleY)*0.02);}
  //
  void     SetX(Float_t x)  { fTX = Short_t(TMath::Nint(x*100.));} 
  void     SetY(Float_t y)  { fTY = Short_t(TMath::Nint(y*100.));} 
  void     SetZ(Float_t z)  { fTZ = Short_t(TMath::Nint(z*100.));} 
  void     SetAngleZ(Float_t anglez) {fTAngleZ = Char_t(TMath::Nint(anglez*50.));}
  void     SetAngleY(Float_t angley) {fTAngleY = Char_t(TMath::Nint(angley*50.));}
  Float_t  GetSigmaZ() const {return (fSigmaZ*0.02);}
  Float_t  GetSigmaY() const {return (fSigmaY*0.02);}  
  Float_t  GetErrZ()   const {return (fErrZ*0.005);}
  Float_t  GetErrY()   const {return (fErrY*0.005);}
  void     SetErrZ(Float_t errz) {fErrZ = UChar_t(TMath::Nint(errz*200.));}
  void     SetErrY(Float_t erry) {fErrY = UChar_t(TMath::Nint(erry*200.));}

  void     SetSigmaZ(Float_t sigmaz) {fSigmaZ = UChar_t(TMath::Nint(sigmaz*50.));}
  void     SetSigmaY(Float_t sigmay) {fSigmaY = UChar_t(TMath::Nint(sigmay*50.));}

  Char_t   IsShared() const {return fIsShared;}
  void     SetShared(Char_t s) {fIsShared=s;}

  //
 private:
  Short_t   fTX;        // x position of the cluster  in cm - 10 mum prec
  Short_t   fTZ;        // current prolongation in Z  in cm - 10 mum prec.
  Short_t   fTY;        // current prolongation in Y  in cm - 10 mum prec.
  Char_t    fTAngleZ;    // angle 
  Char_t    fTAngleY;    // angle 
  UShort_t  fSigmaZ;     // shape  Z - normalised shape - normaliziation 1 - precision 2 percent
  UShort_t  fSigmaY;     // shape  Y - normalised shape - normaliziation 1 - precision 2 percent
  UShort_t  fErrZ;       // z error estimate - in  mm - 50 mum precision 
  UShort_t  fErrY;       // y error estimate - in  mm - 50 mum precision 
  Char_t   fIsShared;     // indicate sharing of the point between several tracks

  ClassDef(AliTPCTrackerPoint,1)  
};

class AliTPCClusterPoint  {
 public:
  AliTPCClusterPoint(){fCZ=fCY=fSigmaZ=fSigmaY=fQ=fMax=fCType=0;}
  virtual ~AliTPCClusterPoint(){}
  Float_t  GetZ() const    {return (fCZ*0.01);}
  Float_t  GetY() const   {return (fCY*0.01);}
  Float_t  GetSigmaZ() const {return (fSigmaZ*0.02);}
  Float_t  GetSigmaY() const {return (fSigmaY*0.02);}  
  Int_t  GetType() const  {return fCType;}
  Int_t  GetMax()  const {return fMax;}
  Float_t  GetQ()  const {return fQ;}

  //
  void     SetY(Float_t y){ fCY = Short_t(TMath::Nint(y*100.));} 
  void     SetZ(Float_t z){ fCZ = Short_t(TMath::Nint(z*100.));} 
  void     SetSigmaZ(Float_t sigmaz) {fSigmaZ = UChar_t(TMath::Nint(sigmaz*50.));}
  void     SetSigmaY(Float_t sigmay) {fSigmaY = UChar_t(TMath::Nint(sigmay*50.));}
  void     SetQ(Float_t q) {fQ = UShort_t(q);}
  void     SetMax(Float_t max) {fMax = UShort_t(max);}
  void     SetType(Char_t type) {fCType = type;}
 private:
  Short_t  fCZ;       // current cluster position Z in cm - 100 mum precision
  Short_t  fCY;       // current cluster position Y in cm - 100 mum precision
  UChar_t  fSigmaZ;   // shape  Z - normalised shape - normaliziation 1 - precision 2 percent
  UChar_t  fSigmaY;   // shape  Y - normalised shape - normaliziation 1 - precision 2 percent
  UShort_t fQ;        // total charge in cluster 
  UShort_t fMax;      // charge at maximum  
  Char_t   fCType;    // type of the cluster

  //
  ClassDef(AliTPCClusterPoint,1)  
};


class AliTPCExactPoint : public TObject{
 public:
  AliTPCExactPoint(){fEZ=fEY=fEAngleZ=fEAngleY=fEAmp=fEPrim=fTrackID=0;}
 private:
  Float_t fEZ;       // current "exact" position according simulation
  Float_t fEY;       // current "exact" position according simulation
  Float_t fEX;       // x poistion of the cluster
  Float_t fEAngleZ;  // angle Z
  Float_t fEAngleY;  // angle Y
  Float_t fEAmp;     // total charge deposited in row
  Float_t fEPrim;    // primary charge deposited in row
  Int_t   fTrackID;  // id of the track
  Int_t   fRow;      // row
  Int_t   fSec;      //sector
  ClassDef(AliTPCExactPoint,1)  
};


class AliTPCTrackPoint: public TObject{
 public:
  AliTPCTrackPoint(){}
  // AliTPCClusterPoint & GetCPoint(){return fCPoint;}
  AliTPCTrackerPoint & GetTPoint(){return fTPoint;}
  AliTPCclusterMI & GetCPoint(){return fCPoint;}  
 private:
  //  AliTPCClusterPoint fCPoint; 
  //Char_t fIsShared;
  AliTPCTrackerPoint fTPoint;  // track point
  AliTPCclusterMI    fCPoint;  // cluster point
  ClassDef(AliTPCTrackPoint,1)  
};

class AliTPCTrackPoint2: public AliTPCTrackPoint{
 public:
  AliTPCTrackPoint2(){}
 private: 
  Float_t fGX;    //global poition of the point
  Float_t fGY;    //global poition of the point
  Float_t fGZ;    //global poition of the point
  //
  Float_t fDY;    //distortion of the clusters from the global helix (3 point interpolation)
  Float_t fDZ;    //distortion of the clusters from the global helix (3 point interpolation)
  //
  Float_t fDYU;  //derivation in y up
  Float_t fDYD;  //distortion of y down
  //
  Float_t fDZU;  //derivation in y up
  Float_t fDZD;  //distortion of y down
  //
  Float_t fDDY;  //derivation in y,z up-down
  Float_t fDDZ;  //derivation in y,z up-down
  //
  Int_t   fID;            //id of the corresponding track
  Int_t   fLab;           //MC label of the track
  ClassDef(AliTPCTrackPoint2,1)  
};




class AliTPCTrackPointRef: public AliTPCTrackPoint{
 public:
  AliTPCExactPoint & GetExactPoint(){return fEPoint;}
  AliTPCExactPoint & GetNearestPoint(){return fNPoint;}  
 private:
  AliTPCExactPoint fEPoint; //exact point belonging to track
  AliTPCExactPoint fNPoint; //nearest point  
  ClassDef(AliTPCTrackPointRef,1)  
};


#endif //ALICOMPLEXCLUSTER_H
