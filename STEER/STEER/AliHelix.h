#ifndef ALIHELIX_H
#define ALIHELIX_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------------------------
//                          Class AliHelix
//
//         Origin: Marian Ivanov marian.ivanov@cern.ch 
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TMath.h>


class AliCluster;
class AliKalmanTrack;
class AliExternalTrackParam;

class AliHelix : public TObject {
public:
  AliHelix();
  AliHelix(const AliHelix &t);
  AliHelix(const AliKalmanTrack &t);
  AliHelix(const AliExternalTrackParam &t);
  AliHelix(Double_t x[3], Double_t p[3], Double_t charge=1, Double_t conversion=0.);
  virtual ~AliHelix(){};
  inline void Evaluate(Double_t t, Double_t r[3]) const;
  void Evaluate(Double_t t, Double_t r[3],  //radius vector
                     Double_t g[3],  //first defivatives
		Double_t gg[3]) const;     //second derivatives
  void GetMomentum(Double_t phase, Double_t p[4], Double_t conversion=0., Double_t *xr=0) const;  // return  momentum  
  void  GetAngle(Double_t t1, const AliHelix &h, Double_t t2, Double_t angle[3]) const;
  inline Double_t GetHelixR(Double_t phase=0) const;
  inline Double_t GetHelixZ(Double_t phase=0) const;

  Double_t  GetPhase(Double_t x0, Double_t y0) const;  //return phase for nearest point
  Double_t  GetPhaseZ(Double_t z0) const;               // return phase for given z
  Int_t     GetPhase(Double_t r0, Double_t t[2]) const;               //return phase for the nearest point
  Int_t     GetRPHIintersections(const AliHelix &h, Double_t phase[2][2], Double_t ri[2], Double_t cut=3.) const;
  Int_t     GetClosestPhases(const AliHelix &h, Double_t phase[2][2]) const;
  Double_t  GetPointAngle(const AliHelix &h, Double_t phase[2],const Float_t *vertex) const;
  Int_t     LinearDCA(const AliHelix &h, Double_t &t1, Double_t &t2,
		      Double_t &R, Double_t &dist) const;
  //
  Int_t     ParabolicDCA(const AliHelix&h,  //helixes
			 Double_t &t1, Double_t &t2, 
			 Double_t &R, Double_t &dist, Int_t iter=1) const;
  Int_t     ParabolicDCA2(const AliHelix&h,  //helixes
			 Double_t &t1, Double_t &t2,
			 Double_t &R, Double_t &dist, Double_t err[3], Int_t iter=1) const;
  Double_t GetHelix(Int_t i) const{return fHelix[i];}
 public:
  Double_t fHelix[9];    //helix parameters
 private:  
  AliHelix &operator=(const AliHelix&helix);

  ClassDef(AliHelix,1)    // AliHelix
};

inline
void AliHelix::Evaluate(Double_t t, Double_t r[3]) const {
  //
  // calculate poitition at given phase t 
  Double_t phase=fHelix[4]*t+fHelix[2];  
  r[0] = fHelix[5] + TMath::Sin(phase)/fHelix[4];
  r[1] = fHelix[0] - TMath::Cos(phase)/fHelix[4];  
  r[2] = fHelix[1] + fHelix[3]*t;
}

inline Double_t  AliHelix::GetHelixR(Double_t phase) const
{
  Double_t x[3];
  Evaluate(phase,x);
  return TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
}

inline Double_t  AliHelix::GetHelixZ(Double_t phase) const
{
  Double_t x[3];
  Evaluate(phase,x);
  return x[2];
}


#endif


