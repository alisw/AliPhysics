// $Id$

#ifndef ALITKETAPHIVECTOR_H
#define ALITKETAPHIVECTOR_H

#include <Riostream.h>

class AliTkEtaPhiVector {
  // simple class to store a vector in eta-phi space
  // some subroutines to calculate distance to another vector
 public:
  AliTkEtaPhiVector() {fX = -9999; fY = -9999;}
  AliTkEtaPhiVector(Double_t x0, Double_t y0) {fX = x0; fY = y0; }
  AliTkEtaPhiVector(const AliTkEtaPhiVector& v)  {fX=v.Eta(); fY=v.Phi();}

  Double_t Eta() const {return fX;}
  Double_t Phi() const {return fY;}

  void setEta(Float_t eta) {fX = eta;}
  void setPhi(Float_t phi) {fY = phi;}
  void setVector(Float_t eta, Float_t phi) {
    fX =eta;
    fY = phi;
  }

  // returns the square difference between the vectors in phi
  Double_t phiDiffSq(AliTkEtaPhiVector& other)
    {
      Double_t phiDiff = TMath::Abs(fY-other.Phi());
      if(phiDiff>TMath::Pi()) phiDiff=2*TMath::Pi()-phiDiff;
      return phiDiff*phiDiff;
    }
  // returns the difference between this vector and the other in phi
  Double_t phiDiff(AliTkEtaPhiVector& other) 
    {
      Double_t phiDiff = TMath::Abs(fY-other.Phi());
      if(phiDiff>TMath::Pi()) phiDiff=2*TMath::Pi()-phiDiff;
      return phiDiff;
    }
  // returns the difference between the vectors in eta
  Double_t etaDiff(AliTkEtaPhiVector& other) 
    {return (fX - other.Eta());}
  // returns the square of the difference between the vectors in eta
  Double_t etaDiffSq(AliTkEtaPhiVector& other) 
    {return (fX - other.Eta())*(fX - other.Eta());}
  // returns the difference between the vectors in eta-phi
  Double_t diff(AliTkEtaPhiVector& other) 
    {return sqrt(diffSq(other));}
  // returns the square of the difference between the vectors in eta-phi
  Double_t diffSq(AliTkEtaPhiVector& other)
    {return (etaDiffSq(other) + phiDiffSq(other));}

 private:
  Double_t fX;
  Double_t fY;
};
  
ostream& operator<<(ostream& s, const AliTkEtaPhiVector& v);
#endif
