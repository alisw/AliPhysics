#ifndef ALITPCPOLYTRACK_H
#define ALITPCPOLYTRACK_H

//-------------------------------------------------------
//                    TPC Cluster Class
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
  inline void GetFitPoint(Double_t x, Double_t &y, Double_t &z);
  inline void GetFitDerivation(Double_t x, Double_t &y, Double_t &z);
  inline void GetFitDerivation2(Double_t &y, Double_t &z);
  void UpdateParameters(Int_t ny, Int_t nz);
  void UpdateParameters();
  Int_t GetN(){return fNPoints;}
  void GetBoundaries(Double_t &xmin, Double_t &xmax){xmin = fMinX;xmax=fMaxX;}
  void Refit(AliTPCpolyTrack & track, Double_t deltay, Double_t deltaz); 
  void Refit(AliTPCpolyTrack & track, Double_t deltay, Double_t deltaz, Int_t nfirst, Int_t ny, Int_t nz); 
private: 
  void   Fit2(Double_t fSumY, Double_t fSumYX, Double_t fSumYX2,
	      Double_t fSumX,  Double_t fSumX2, Double_t fSumX3, 
	      Double_t fSumX4, Double_t fSumW, Double_t &a, Double_t &b, Double_t &c);
  void  Fit1(Double_t fSumY, Double_t fSumYX, 
	      Double_t fSumX,  Double_t fSumX2, 
	      Double_t fSumW, Double_t &a, Double_t &b, Double_t &c);
  //
  Double_t fA;
  Double_t fB;
  Double_t fC;
  Double_t fD;
  Double_t fE;
  Double_t fF;
  Double_t fMaxX;
  Double_t fMinX;
  //
  Double_t fSumW;   // sum of the weight 

  Double_t fSumX;    //
  Double_t fSumX2;   //
  Double_t fSumX3;   //  
  Double_t fSumX4;   //
  Double_t fSumY;    //
  Double_t fSumYX;   //  
  Double_t fSumYX2;  //
  Double_t fSumZ;    //
  Double_t fSumZX;   //
  Double_t fSumZX2;  //
  
  Double_t fX[200];
  Double_t fY[200];
  Double_t fSY[200];
  Double_t fZ[200];
  Double_t fSZ[200];

  Int_t fNPoints; 

  ClassDef(AliTPCpolyTrack,1)  // Time Projection "polynomial track"
};

void AliTPCpolyTrack::GetFitPoint(Double_t x, Double_t &y, Double_t &z)
{
  y = fA+fB*x+fC*x*x;
  z = fD+fE*x+fF*x*x;
}


void AliTPCpolyTrack::GetFitDerivation(Double_t x, Double_t &y, Double_t &z)
{

  y = fB+2.*fC*x;
  z = fE+2.*fF*x;
  
}

void AliTPCpolyTrack::GetFitDerivation2(Double_t &y, Double_t &z)
{

  y = 2.*fC;
  z = 2.*fF;
  
}


#endif


