#ifndef ALITPCTRACK_H
#define ALITPCTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------
//                    TPC Track Class
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

#include <TObject.h>
#include <TMath.h>
#include <iostream.h>

class AliTPCClustersArray;
class AliTPCcluster;

//_____________________________________________________________________________
class AliTPCtrack : public TObject {
public:
   AliTPCtrack();
   AliTPCtrack(UInt_t index, const Double_t xx[5], 
               const Double_t cc[15], Double_t xr, Double_t alpha); 
   AliTPCtrack(const AliTPCtrack& t);
   Int_t Compare(TObject *o);
   Int_t PropagateTo(Double_t xr,
                     Double_t x0=28.94,Double_t rho=0.9e-3,Double_t pm=0.139);
   void PropagateToVertex(
                     Double_t x0=36.66,Double_t rho=1.2e-3,Double_t pm=0.139);
   void Update(const AliTPCcluster* c, Double_t chi2, UInt_t i);
   Int_t Rotate(Double_t angle);
   void CookLabel(AliTPCClustersArray *carray);
   void SetLabel(Int_t lab=0);
   void SetdEdx(Float_t dedx);

   Double_t GetPredictedChi2(const AliTPCcluster *cluster) const;
   Bool_t IsSortable() const;
   Int_t GetLabel() const;
   void GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
   Float_t GetdEdx() const;
   Double_t GetX() const;
   Double_t GetY() const;
   Double_t GetZ() const;
   Double_t GetC() const;
   Double_t GetEta() const;
   Double_t GetTgl() const;
   Double_t GetPt() const;
   Double_t GetP() const;
   void GetCovariance(Double_t cov[15]) const;
   Double_t GetSigmaY2() const;
   Double_t GetSigmaZ2() const;
   Double_t GetSigmaC2() const;
   Double_t GetSigmaTgl2() const;
   Double_t GetAlpha() const;
   Double_t GetChi2() const;
   Int_t GetNumberOfClusters() const;
   void GetCluster(Int_t i, Int_t &sec, Int_t &row, Int_t &ncl) const;

protected: 
   Short_t fN;               // number of clusters 
private:
   Int_t fLab;               // track label
   Double_t fChi2;           // total chi2 value for this track
   Float_t fdEdx;            // dE/dx 

   Double_t fAlpha;          // rotation angle
   Double_t fX;              // X-coordinate of this track (reference plane)

   Double_t fY;              // Y-coordinate of the track
   Double_t fZ;              // Z-coordinate of the track
   Double_t fC;              // track curvature
   Double_t fE;              // C*x0
   Double_t fT;              // tangent of the track dip angle

   Double_t fCyy;                         // covariance
   Double_t fCzy, fCzz;                   // matrix
   Double_t fCcy, fCcz, fCcc;             // of the
   Double_t fCey, fCez, fCec, fCee;       // track
   Double_t fCty, fCtz, fCtc, fCte, fCtt; // parameters

   UInt_t fIndex[200];       // (((row<<8)+sec)<<16)+index

   ClassDef(AliTPCtrack,1)   // Time Projection Chamber reconstructed tracks
};

inline AliTPCtrack::AliTPCtrack() {
  //default contructor
  fN=0;
}

inline void AliTPCtrack::SetLabel(Int_t lab) {
  //just to calm down our rule checker
  fLab=lab;
}

inline void AliTPCtrack::SetdEdx(Float_t dedx) {
  //just to calm down our rule checker
  fdEdx=dedx;
}

inline Bool_t AliTPCtrack::IsSortable() const {
  //just to calm down our rule checker
  return kTRUE;
}

inline Int_t AliTPCtrack::GetLabel() const {
  //just to calm down our rule checker
  return fLab;
}

inline Float_t AliTPCtrack::GetdEdx() const {
  //just to calm down our rule checker
  return fdEdx;
}

inline Double_t AliTPCtrack::GetX() const {
  //just to calm down our rule checker
  return fX;
}

inline Double_t AliTPCtrack::GetY() const {
  //just to calm down our rule checker
  return fY;
}

inline Double_t AliTPCtrack::GetZ() const {
  //just to calm down our rule checker
  return fZ;
}

inline Double_t AliTPCtrack::GetC() const {
  //just to calm down our rule checker
  return fC;
}

inline Double_t AliTPCtrack::GetEta() const {
  //just to calm down our rule checker
  return fE;
}

inline Double_t AliTPCtrack::GetTgl() const {
  //just to calm down our rule checker
  return fT;
}

inline Double_t AliTPCtrack::GetPt() const {
  //just to calm down our rule checker
  return 0.3*0.2/GetC()/100;
}

inline Double_t AliTPCtrack::GetP() const {
  //just to calm down our rule checker
  return TMath::Abs(GetPt())*sqrt(1.+GetTgl()*GetTgl());
}

inline Double_t AliTPCtrack::GetSigmaY2() const {
  //just to calm down our rule checker
  return fCyy;
}

inline Double_t AliTPCtrack::GetSigmaZ2() const {
  //just to calm down our rule checker
  return fCzz;
}

inline Double_t AliTPCtrack::GetSigmaC2() const {
  //just to calm down our rule checker
  return fCcc;
}

inline Double_t AliTPCtrack::GetSigmaTgl2() const {
  //just to calm down our rule checker
  return fCtt;
}

inline Double_t AliTPCtrack::GetAlpha() const {
  //just to calm down our rule checker
  return fAlpha;
}

inline Double_t AliTPCtrack::GetChi2() const {
  //just to calm down our rule checker
  return fChi2;
}

inline Int_t AliTPCtrack::GetNumberOfClusters() const {
  //just to calm down our rule checker
  return fN;
}

inline 
void AliTPCtrack::GetCluster(Int_t i,Int_t &sec,Int_t &row,Int_t &ncl) const {
  //just to calm down our rule checker
  Int_t index=fIndex[i];
  sec=(index&0xff000000)>>24; 
  row=(index&0x00ff0000)>>16; 
  ncl=(index&0x0000ffff)>>00;
}

#endif


