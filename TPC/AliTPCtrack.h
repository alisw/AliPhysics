#ifndef ALITPCTRACK_H
#define ALITPCTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------
//                    TPC Track Class
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

#include <AliKalmanTrack.h>
#include <TMath.h>

class AliTPCClustersArray;
class AliTPCcluster;

//_____________________________________________________________________________
class AliTPCtrack : public AliKalmanTrack {
public:
  AliTPCtrack():AliKalmanTrack(){}
  AliTPCtrack(UInt_t index, const Double_t xx[5], 
              const Double_t cc[15], Double_t xr, Double_t alpha); 
  AliTPCtrack(const AliTPCtrack& t);
  Int_t PropagateToVertex(
                    Double_t x0=36.66,Double_t rho=1.2e-3,Double_t pm=0.139);
  Int_t Rotate(Double_t angle);
  void CookLabel(AliTPCClustersArray *carray);
  void SetdEdx(Float_t dedx) {fdEdx=dedx;}

  Float_t GetdEdx() const {return fdEdx;}
  Double_t GetX()   const {return fX;}
  Double_t GetY()   const {return fP0;}
  Double_t GetZ()   const {return fP1;}
  Double_t GetEta() const {return fP2;}
  Double_t GetC()   const {return fP3;}
  Double_t GetTgl() const {return fP4;}
  Double_t GetSigmaY2() const {return fC00;}
  Double_t GetSigmaZ2() const {return fC11;}
  Double_t GetSigmaC2() const {return fC33;}
  Double_t GetAlpha()   const {return fAlpha;}
  void GetCluster(Int_t i, Int_t &sec, Int_t &row, Int_t &ncl) const;

  Double_t GetPt() const {return 0.299792458*0.2/GetC()/100;}
  Double_t GetP()  const {return TMath::Abs(GetPt())*
    TMath::Sqrt(1.+GetTgl()*GetTgl());}
  void GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const;
  Int_t PropagateTo(Double_t xr,
                    Double_t x0=28.94,Double_t rho=0.9e-3,Double_t pm=0.139);
  void Update(const AliCluster* c, Double_t chi2, UInt_t i);

protected: 
  UInt_t fIndex[200];       // indices of associated clusters 

  Float_t fdEdx;            // dE/dx 

  Double_t fAlpha;          // rotation angle
  Double_t fX;              // X-coordinate of this track (reference plane)

  ClassDef(AliTPCtrack,1)   // Time Projection Chamber reconstructed tracks
};


inline 
void AliTPCtrack::GetCluster(Int_t i,Int_t &sec,Int_t &row,Int_t &ncl) const {
  //return sector, pad row and the index of the i-th cluster of this track 
  Int_t index=fIndex[i];
  sec=(index&0xff000000)>>24; 
  row=(index&0x00ff0000)>>16; 
  ncl=(index&0x0000ffff)>>00;
}

#endif


