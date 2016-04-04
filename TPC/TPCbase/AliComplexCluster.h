#ifndef ALICOMPLEXCLUSTER_H
#define ALICOMPLEXCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliComplexCluster
///
/// \author M. Ivanov

#include "TObject.h"
#include "TMath.h"
#include "AliTPCclusterMI.h"
#include "limits.h"

namespace AliComplexClusterHelper {
  inline Int_t Clamp(Int_t val, Int_t min, Int_t max)
  {
    // limit to range, inspire by boost clamp
    return (val<min)?min:(val>max)?max:val;
  }
}


class AliComplexCluster : public TObject {
public:

  AliComplexCluster();
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
  Int_t     fTracks[3];///< labels of overlapped tracks
  Float_t   fX ;       ///< Y of cluster
  Float_t   fY ;       ///< Z of cluster
  Float_t   fQ ;       ///< Q of cluster (in ADC counts)
  Float_t   fSigmaX2;  ///< Sigma Y square of cluster
  Float_t   fSigmaY2;  ///< Sigma Z square of cluster
  Float_t   fSigmaXY;  ///< XY moment
  Float_t   fArea;     ///< area of cluster
  Float_t   fMax;     ///< amplitude at maximum

  /// \cond CLASSIMP
  ClassDef(AliComplexCluster, 1);
  /// \endcond

    // Cluster manager
};

//RS: this is old, non-economical and virtualized class, for bwd compatibility only
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
  AliTPCTrackerPoint &operator=(const AliTPCTrackerPoint& o);
  Float_t  GetX() const  {return (fTX*0.01);}
  Float_t  GetZ() const {return (fTZ*0.01);}
  Float_t  GetY() const {return (fTY*0.01);}
  Float_t  GetAngleZ() const  {return (Float_t(fTAngleZ)*0.02);}
  Float_t  GetAngleY() const {return (Float_t(fTAngleY)*0.02);}
  //
  void     SetX(Float_t x)  { fTX = Short_t(AliComplexClusterHelper::Clamp(TMath::Nint(x*100.), SHRT_MIN, SHRT_MAX));}
  void     SetY(Float_t y)  { fTY = Short_t(AliComplexClusterHelper::Clamp(TMath::Nint(y*100.), SHRT_MIN, SHRT_MAX));}
  void     SetZ(Float_t z)  { fTZ = Short_t(AliComplexClusterHelper::Clamp(TMath::Nint(z*100.), SHRT_MIN, SHRT_MAX));}
  void     SetAngleZ(Float_t anglez) {fTAngleZ = Char_t(AliComplexClusterHelper::Clamp(TMath::Nint(anglez*50.),SCHAR_MIN,SCHAR_MAX));}
  void     SetAngleY(Float_t angley) {fTAngleY = Char_t(AliComplexClusterHelper::Clamp(TMath::Nint(angley*50.),SCHAR_MIN,SCHAR_MAX));}
  Float_t  GetSigmaZ() const {return (fSigmaZ*0.02);}
  Float_t  GetSigmaY() const {return (fSigmaY*0.02);}  
  Float_t  GetErrZ()   const {return (fErrZ*0.005);}
  Float_t  GetErrY()   const {return (fErrY*0.005);}
  void     SetErrZ(Float_t errz) {fErrZ = UChar_t(AliComplexClusterHelper::Clamp(TMath::Nint(errz*200.),0,UCHAR_MAX));}
  void     SetErrY(Float_t erry) {fErrY = UChar_t(AliComplexClusterHelper::Clamp(TMath::Nint(erry*200.),0,UCHAR_MAX));}

  void     SetSigmaZ(Float_t sigmaz) {fSigmaZ = UChar_t(AliComplexClusterHelper::Clamp(TMath::Nint(sigmaz*50.),0,UCHAR_MAX));}
  void     SetSigmaY(Float_t sigmay) {fSigmaY = UChar_t(AliComplexClusterHelper::Clamp(TMath::Nint(sigmay*50.),0,UCHAR_MAX));}

  Char_t   IsShared() const {return fIsShared;}
  void     SetShared(Char_t s) {fIsShared=s;}
  //
 private:
  Short_t   fTX;        ///< x position of the cluster  in cm - 10 mum prec
  Short_t   fTZ;        ///< current prolongation in Z  in cm - 10 mum prec.
  Short_t   fTY;        ///< current prolongation in Y  in cm - 10 mum prec.
  Char_t    fTAngleZ;    ///< angle
  Char_t    fTAngleY;    ///< angle
  UShort_t  fSigmaZ;     ///< shape  Z - normalised shape - normaliziation 1 - precision 2 percent
  UShort_t  fSigmaY;     ///< shape  Y - normalised shape - normaliziation 1 - precision 2 percent
  UShort_t  fErrZ;       ///< z error estimate - in  mm - 50 mum precision
  UShort_t  fErrY;       ///< y error estimate - in  mm - 50 mum precision
  Char_t   fIsShared;     ///< indicate sharing of the point between several tracks

  /// \cond CLASSIMP
  ClassDef(AliTPCTrackerPoint, 2);
  /// \endcond  
};

//RS: this is new, more economic class for TrackerPoints, avoiding vtable pointer per point
class AliTPCTrackerPoints  {  
  //
 public:
  class Point {
  public:
  Point():fTX(0),fTZ(0),fTY(0),fTAngleZ(0),fTAngleY(0),fSigmaZ(0),fSigmaY(0),fErrZ(0),fErrY(0) {}
    ~Point(){}
    Point &operator=(const Point& o) {
      if(this!=&o) {
	fTX = o.fTX;
	fTY = o.fTY;
	fTZ = o.fTZ;
	fTAngleZ = o.fTAngleZ;
	fTAngleY = o.fTAngleY;
	fSigmaZ = o.fSigmaZ;
	fSigmaY = o.fSigmaY;
	fErrZ   = o.fErrZ;
	fErrY   = o.fErrY;
      }
      return *this;
    }
    //
    Point &operator=(const AliTPCTrackerPoint& o) {
      fTX = o.GetX();
      fTY = o.GetY();
      fTZ = o.GetZ();
      fTAngleZ = o.GetAngleZ();
      fTAngleY = o.GetAngleY();
      fSigmaZ = o.GetSigmaZ();
      fSigmaY = o.GetSigmaY();
      fErrZ   = o.GetErrZ();
      fErrY   = o.GetErrY();
      return *this;
    }
    Float_t  GetX() const  {return (fTX*0.01);}
    Float_t  GetZ() const {return (fTZ*0.01);}
    Float_t  GetY() const {return (fTY*0.01);}
    Float_t  GetAngleZ() const  {return (Float_t(fTAngleZ)*0.02);}
    Float_t  GetAngleY() const {return (Float_t(fTAngleY)*0.02);}
    //
    void     SetX(Float_t x)  { fTX = Short_t(AliComplexClusterHelper::Clamp(TMath::Nint(x*100.), SHRT_MIN, SHRT_MAX));}
    void     SetY(Float_t y)  { fTY = Short_t(AliComplexClusterHelper::Clamp(TMath::Nint(y*100.), SHRT_MIN, SHRT_MAX));}
    void     SetZ(Float_t z)  { fTZ = Short_t(AliComplexClusterHelper::Clamp(TMath::Nint(z*100.), SHRT_MIN, SHRT_MAX));}
    void     SetAngleZ(Float_t anglez) {fTAngleZ = Char_t(AliComplexClusterHelper::Clamp(TMath::Nint(anglez*50.),SCHAR_MIN,SCHAR_MAX));}
    void     SetAngleY(Float_t angley) {fTAngleY = Char_t(AliComplexClusterHelper::Clamp(TMath::Nint(angley*50.),SCHAR_MIN,SCHAR_MAX));}
    //
    Float_t  GetSigmaZ() const {return (fSigmaZ*0.02);}
    Float_t  GetSigmaY() const {return (fSigmaY*0.02);}  
    Float_t  GetErrZ()   const {return (fErrZ*0.005);}
    Float_t  GetErrY()   const {return (fErrY*0.005);}
    //
    void     SetErrZ(Float_t errz) {fErrZ = UChar_t(AliComplexClusterHelper::Clamp(TMath::Nint(errz*200.),0,UCHAR_MAX));}
    void     SetErrY(Float_t erry) {fErrY = UChar_t(AliComplexClusterHelper::Clamp(TMath::Nint(erry*200.),0,UCHAR_MAX));}
    void     SetSigmaZ(Float_t sigmaz) {fSigmaZ = UChar_t(AliComplexClusterHelper::Clamp(TMath::Nint(sigmaz*50.),0,UCHAR_MAX));}
    void     SetSigmaY(Float_t sigmay) {fSigmaY = UChar_t(AliComplexClusterHelper::Clamp(TMath::Nint(sigmay*50.),0,UCHAR_MAX));}

    //
  private:
    Short_t   fTX;        ///< x position of the cluster  in cm - 10 mum prec
    Short_t   fTZ;        ///< current prolongation in Z  in cm - 10 mum prec.
    Short_t   fTY;        ///< current prolongation in Y  in cm - 10 mum prec.
    Char_t    fTAngleZ;    ///< angle
    Char_t    fTAngleY;    ///< angle
    UChar_t   fSigmaZ;     ///< shape  Z - normalised shape - normaliziation 1 - precision 2 percent
    UChar_t   fSigmaY;     ///< shape  Y - normalised shape - normaliziation 1 - precision 2 percent
    UChar_t   fErrZ;       ///< z error estimate - in  mm - 50 mum precision
    UChar_t   fErrY;       ///< y error estimate - in  mm - 50 mum precision
  };
  
 public:
  AliTPCTrackerPoints();
  ~AliTPCTrackerPoints() {}
  AliTPCTrackerPoints &operator=(const AliTPCTrackerPoints& o);
  AliTPCTrackerPoints(const AliTPCTrackerPoints& o);
  void Clear();
  //
  void     SetShared(int i)      {fShared[i/8] |= 0x1<<(i%8);}
  Bool_t   IsShared(int i) const {return fShared[i/8] & 0x1<<(i%8);}
  const Point*   GetPoint(int i) const {return &fPoints[i];}
  void     SetPoint(int i, const AliTPCTrackerPoint* o);
  //
 protected:
  //
  UChar_t  fShared[20];         // shared map
  Point    fPoints[159];        // array of compact non-virtualized points
  //
  ClassDef(AliTPCTrackerPoints, 1);
};

class AliTPCClusterPoint  {
 public:
  AliTPCClusterPoint():
                      fCZ(0),
                      fCY(0), 
                      fSigmaZ(0),
                      fSigmaY(0),
                      fQ(0),
                      fMax(0),
                      fCType(0){}
  virtual ~AliTPCClusterPoint(){}
  Float_t  GetZ() const    {return (fCZ*0.01);}
  Float_t  GetY() const   {return (fCY*0.01);}
  Float_t  GetSigmaZ() const {return (fSigmaZ*0.02);}
  Float_t  GetSigmaY() const {return (fSigmaY*0.02);}  
  Int_t  GetType() const  {return fCType;}
  Int_t  GetMax()  const {return fMax;}
  Float_t  GetQ()  const {return fQ;}

  //
  void     SetY(Float_t y){ fCY = Short_t(AliComplexClusterHelper::Clamp(TMath::Nint(y*100.), SHRT_MIN, SHRT_MAX));}
  void     SetZ(Float_t z){ fCZ = Short_t(AliComplexClusterHelper::Clamp(TMath::Nint(z*100.), SHRT_MIN, SHRT_MAX));}
  void     SetSigmaZ(Float_t sigmaz) {fSigmaZ = UChar_t(AliComplexClusterHelper::Clamp(TMath::Nint(sigmaz*50.),0,UCHAR_MAX));}
  void     SetSigmaY(Float_t sigmay) {fSigmaY = UChar_t(AliComplexClusterHelper::Clamp(TMath::Nint(sigmay*50.),0,UCHAR_MAX));}

  void     SetQ(Float_t q) {fQ = UShort_t(q);}
  void     SetMax(Float_t max) {fMax = UShort_t(max);}
  void     SetType(Char_t type) {fCType = type;}
 private:
  Short_t  fCZ;       ///< current cluster position Z in cm - 100 mum precision
  Short_t  fCY;       ///< current cluster position Y in cm - 100 mum precision
  UChar_t  fSigmaZ;   ///< shape  Z - normalised shape - normaliziation 1 - precision 2 percent
  UChar_t  fSigmaY;   ///< shape  Y - normalised shape - normaliziation 1 - precision 2 percent
  UShort_t fQ;        ///< total charge in cluster
  UShort_t fMax;      ///< charge at maximum
  Char_t   fCType;    ///< type of the cluster

  /// \cond CLASSIMP
  ClassDef(AliTPCClusterPoint, 1);
  /// \endcond  
};


class AliTPCExactPoint : public TObject{
 public:
  AliTPCExactPoint():TObject(),
    fEZ(0.),
    fEY(0.),
    fEX(0.),
    fEAngleZ(0.),
    fEAngleY(0.),
    fEAmp(0.),
    fEPrim(0.),
    fTrackID(0),
    fRow(0),
    fSec(0){}
 private:
  Float_t fEZ;       ///< current "exact" position according simulation
  Float_t fEY;       ///< current "exact" position according simulation
  Float_t fEX;       ///< x poistion of the cluster
  Float_t fEAngleZ;  ///< angle Z
  Float_t fEAngleY;  ///< angle Y
  Float_t fEAmp;     ///< total charge deposited in row
  Float_t fEPrim;    ///< primary charge deposited in row
  Int_t   fTrackID;  ///< id of the track
  Int_t   fRow;      ///< row
  Int_t   fSec;      ///< sector

  /// \cond CLASSIMP
  ClassDef(AliTPCExactPoint, 1);
  /// \endcond 
};


class AliTPCTrackPoint: public TObject{
 public:
  AliTPCTrackPoint():TObject(),
    fTPoint(),
    fCPoint(){}

  // AliTPCClusterPoint & GetCPoint(){return fCPoint;}
  AliTPCTrackerPoint & GetTPoint(){return fTPoint;}
  AliTPCclusterMI & GetCPoint(){return fCPoint;}  
 private:
  //  AliTPCClusterPoint fCPoint; 
  //Char_t fIsShared;
  AliTPCTrackerPoint fTPoint;  ///< track point
  AliTPCclusterMI    fCPoint;  ///< cluster point

  /// \cond CLASSIMP
  ClassDef(AliTPCTrackPoint, 1);
  /// \endcond
};

class AliTPCTrackPoint2: public AliTPCTrackPoint{
 public:
  AliTPCTrackPoint2():AliTPCTrackPoint(),
    fGX(0.),
    fGY(0.),
    fGZ(0.),
    fDY(0.),
    fDZ(0.),
    fDYU(0.),
    fDYD(0),
    fDZU(0.),
    fDZD(0.),
    fDDY(0),
    fDDZ(0.),
    fID(0),
    fLab(0){}
 private: 
  Float_t fGX;    ///< global poition of the point
  Float_t fGY;    ///< global poition of the point
  Float_t fGZ;    ///< global poition of the point
  //
  Float_t fDY;    ///< distortion of the clusters from the global helix (3 point interpolation)
  Float_t fDZ;    ///< distortion of the clusters from the global helix (3 point interpolation)
  //
  Float_t fDYU;  ///< derivation in y up
  Float_t fDYD;  ///< distortion of y down
  //
  Float_t fDZU;  ///< derivation in y up
  Float_t fDZD;  ///< distortion of y down
  //
  Float_t fDDY;  ///< derivation in y,z up-down
  Float_t fDDZ;  ///< derivation in y,z up-down
  //
  Int_t   fID;            ///< id of the corresponding track
  Int_t   fLab;           ///< MC label of the track

  /// \cond CLASSIMP
  ClassDef(AliTPCTrackPoint2, 1);
  /// \endcond
};

#endif //ALICOMPLEXCLUSTER_H
