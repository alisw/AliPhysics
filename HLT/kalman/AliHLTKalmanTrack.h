// @(#) $Id$

#ifndef ALIL3_KALMANTRACK
#define ALIL3_KALMANTRACK

/*
* The state vector is:
*  fP0 (x[0]) : local y-coordinate
*  fP1 (x[1]) : local z-coordinate
*  fP2 (x[2]) : local sine of track momentum dip angle 
*  fP3 (x[3]) : tangent of track momentum dip angle
*  fP4 (x[4]) : 1/pt
*
* The covariance matrix is:
*  fC00                     
*  fC10 fC11
*  fC20 fC21 fC22
*  fC30 fC31 fC32 fC33
*  fC40 fC41 fC42 fC43 fC44
*
* To accsess this use: GetStateVector(Float_t xx[5])
*                      GetCovariance(Float_t xx[15])              
*/

#include "AliHLTRootTypes.h"
#include "AliHLTTrack.h"

// includes for offline comparison, will be removed
#include "AliTPCtrack.h"
// includes for offline comparison, will be removed

class AliHLTSpacePointData;

class AliHLTKalmanTrack : public AliHLTTrack {
//class AliHLTKalmanTrack {

 private:

  Float_t fP0;             // Y-coordinate of a track
  Float_t fP1;             // Z-coordinate of a track
  Float_t fP2;             // C*x0
  Float_t fP3;             // tangent of the track momentum dip angle
  Float_t fP4;             // track curvature

  Float_t fC00;                         // covariance
  Float_t fC10, fC11;                   // matrix
  Float_t fC20, fC21, fC22;             // of the
  Float_t fC30, fC31, fC32, fC33;       // track
  Float_t fC40, fC41, fC42, fC43, fC44; // parameters

  Float_t fChisq; 
  Float_t fMaxChi2;

  Float_t fX;

  Int_t fNHits;

 public:

  AliHLTKalmanTrack();
  virtual ~AliHLTKalmanTrack();
Int_t MakeSeed(AliHLTTrack *track, AliHLTSpacePointData *points0, UInt_t pos0, Int_t slice0, AliHLTSpacePointData *points1, UInt_t pos1, Int_t slice1, AliHLTSpacePointData *points2, UInt_t pos2, Int_t slice2);
  Int_t Init(AliHLTTrack *track, AliHLTSpacePointData *points, UInt_t pos,Int_t slice);
  Int_t Propagate(AliHLTSpacePointData *points, UInt_t pos, Int_t slice);
  Int_t UpdateTrack(AliHLTSpacePointData *points, UInt_t pos, Int_t slice);
  Int_t UpdateTrackII(AliHLTSpacePointData *points, UInt_t pos);
  void AddTrack();
  //  Float_t GetStateVector(Float_t xx[5]) const {
  void GetStateVector(Float_t xx[5]) const {
    xx[0] = fP0;
    xx[1] = fP1;
    xx[2] = fP2;
    xx[3] = fP3;
    xx[4] = fP4;
  }

  Float_t GetX0() {return fP0;}
  Float_t GetX1() {return fP1;}
  Float_t GetX2() {return fP2;}
  Float_t GetX3() {return fP3;}
  Float_t GetX4() {return fP4;}

  Float_t GetC0() {return fC00;}
  Float_t GetC1() {return fC10;}
  Float_t GetC2() {return fC11;}
  Float_t GetC3() {return fC20;}
  Float_t GetC4() {return fC21;}
  Float_t GetC5() {return fC22;}
  Float_t GetC6() {return fC30;}
  Float_t GetC7() {return fC31;}
  Float_t GetC8() {return fC32;}
  Float_t GetC9() {return fC33;}
  Float_t GetC10() {return fC40;}
  Float_t GetC11() {return fC41;}
  Float_t GetC12() {return fC42;}
  Float_t GetC13() {return fC43;}
  Float_t GetC14() {return fC44;}

  void SetX(Float_t f) {fX = f;}
  void SetX0(Float_t f) {fP0 = f;}
  void SetX1(Float_t f) {fP1 = f;}
  void SetX2(Float_t f) {fP2 = f;}
  void SetX3(Float_t f) {fP3 = f;}
  void SetX4(Float_t f) {fP4 = f;}

  void SetC0(Float_t f) {fC00 = f;}
  void SetC1(Float_t f) {fC10 = f;}
  void SetC2(Float_t f) {fC11 = f;}
  void SetC3(Float_t f) {fC20 = f;}
  void SetC4(Float_t f) {fC21 = f;}
  void SetC5(Float_t f) {fC22 = f;}
  void SetC6(Float_t f) {fC30 = f;}
  void SetC7(Float_t f) {fC31 = f;}
  void SetC8(Float_t f) {fC32 = f;}
  void SetC9(Float_t f) {fC33 = f;}
  void SetC10(Float_t f) {fC40 = f;}
  void SetC11(Float_t f) {fC41 = f;}
  void SetC12(Float_t f) {fC42 = f;}
  void SetC13(Float_t f) {fC43 = f;} 
  void SetC14(Float_t f) {fC44 = f;}

  //  Float_t GetCovariance(Float_t cc[15]) const {
  void GetCovariance(Float_t cc[15]) const {
    cc[0 ]=fC00;
    cc[1 ]=fC10;  cc[2 ]=fC11;
    cc[3 ]=fC20;  cc[4 ]=fC21;  cc[5 ]=fC22;
    cc[6 ]=fC40;  cc[7 ]=fC41;  cc[8 ]=fC42;  cc[9 ]=fC44;
    cc[10]=fC30;  cc[11]=fC31;  cc[12]=fC32;  cc[13]=fC43;  cc[14]=fC33;
  }

  Float_t GetChisq() {if(!fChisq) return 0; return fChisq;} 
  Float_t GetX() {return fX;}
  void SetStateVector(Float_t f[5]) {fP0 = f[0]; fP1 = f[1]; fP2 = f[2]; 
                                     fP3 = f[3]; fP4 = f[4];}
  void SetCovariance(Float_t f[15]) {fC00 = f[0];  fC10 = f[1]; fC11 = f[2];
  fC21 = f[3]; fC21 = f[4]; fC22 = f[5]; fC30 = f[6]; fC31 = f[7]; fC32 = f[8];
  fC33 = f[9]; fC40 = f[10]; fC41 = f[11]; fC42 = f[12]; fC43 = f[13]; 
  fC44 = f[14];}
  void SetChisq(Float_t f) {fChisq = f;}
  void SetMaxChi2(Float_t f) {fMaxChi2 = f;}
  //Float_t GetChisqIncrement(AliHLTSpacePointData *points, UInt_t pos);
  Float_t GetChisqIncrement(Float_t y, Float_t error_y, Float_t z, Float_t error_z);
  Float_t f2(Float_t x1,Float_t y1, Float_t x2,Float_t y2, Float_t x3,Float_t y3);
  Float_t f3(Float_t x1,Float_t y1, Float_t x2,Float_t y2, Float_t z1,Float_t z2);
  Float_t f4(Float_t x1,Float_t y1, Float_t x2,Float_t y2, Float_t x3,Float_t y3);
  void Set(AliHLTKalmanTrack *track);

  Int_t GetNHits() const {return fNHits;}
  void SetNHits(Int_t f) {fNHits = f;}

  Int_t PropagateOfflineTrack(Double_t x, Double_t y, Double_t z, Double_t ey, Double_t ez);
  Int_t UpdateOfflineTrack(Double_t x, Double_t y, Double_t z, Double_t ey, Double_t ez);
  Float_t GetChisqIncrementOfflineTrack(Double_t y, Double_t z, Double_t ey, Double_t ez);
  
};

typedef AliHLTKalmanTrack AliL3KalmanTrack; // for backward compatibility

#endif
