#ifndef ALITPCPOLYTRACK_H
#define ALITPCPOLYTRACK_H

//-------------------------------------------------------
//                    TPC Cluster Class
//     track fitting using the polynom
//     allows bending in both direction
//
//   Origin: Marian Ivanov
//-------------------------------------------------------

#include "TObject.h"

//_____________________________________________________________________________
class AliTPCpolyTrack : public TObject {
public:
  AliTPCpolyTrack();
  void Reset();
  void AddPoint(Double_t x, Double_t y, Double_t z, Double_t sy=1, Double_t sz=1);
  inline void GetFitPoint(Double_t x, Double_t &y, Double_t &z) const;
  inline void GetFitDerivation(Double_t x, Double_t &y, Double_t &z) const;
  inline void GetFitDerivation2(Double_t &y, Double_t &z) const;
  void UpdateParameters(Int_t ny, Int_t nz);
  void UpdateParameters();
  Int_t GetN() const {return fNPoints;}
  void GetBoundaries(Double_t &xmin, Double_t &xmax) const 
                     {xmin = fMinX;xmax=fMaxX;}
  void Refit(AliTPCpolyTrack & track, Double_t deltay, Double_t deltaz); 
  void Refit(AliTPCpolyTrack & track, Double_t deltay, Double_t deltaz, Int_t nfirst, Int_t ny, Int_t nz); 
private: 
  void   Fit2( Double_t &a, Double_t &b, Double_t &c);
  void  Fit1(Double_t &a, Double_t &b, Double_t &c);
  //
  Double_t fA; // parameter
  Double_t fB; // parameter
  Double_t fC; // parameter
  Double_t fD; // parameter
  Double_t fE; // parameter
  Double_t fF; // parameter 
  Double_t fMaxX; // X range
  Double_t fMinX; // X range
  //
  Double_t fSumW;   // sum of the weight 

  Double_t fSumX;    // weighted summ of X 
  Double_t fSumX2;   // weighted summ of X**2
  Double_t fSumX3;   // weighted summ of X**3 
  Double_t fSumX4;   // weighted summ of X**4
  Double_t fSumY;    // as X
  Double_t fSumYX;   //  summ of X*Y
  Double_t fSumYX2;  // summ of X**2*Y
  Double_t fSumZ;    // summ of Z 
  Double_t fSumZX;   // summ of Z*X
  Double_t fSumZX2;  // summ of Z*X**2
  
  Double_t fX[200]; // array of coordinates 
  Double_t fY[200]; // array of coordinates
  Double_t fSY[200]; // array of coordinates
  Double_t fZ[200]; // array of coordinates
  Double_t fSZ[200]; // array of coordinates

  Int_t fNPoints; // No of points

  ClassDef(AliTPCpolyTrack,1)  // Time Projection "polynomial track"
};

void AliTPCpolyTrack::GetFitPoint(Double_t x, Double_t &y, Double_t &z) const
{
  y = fA+fB*x+fC*x*x;
  z = fD+fE*x+fF*x*x;
}


void AliTPCpolyTrack::GetFitDerivation(Double_t x, Double_t &y, Double_t &z) const
{

  y = fB+2.*fC*x;
  z = fE+2.*fF*x;
  
}

void AliTPCpolyTrack::GetFitDerivation2(Double_t &y, Double_t &z) const
{

  y = 2.*fC;
  z = 2.*fF;
  
}


#endif


