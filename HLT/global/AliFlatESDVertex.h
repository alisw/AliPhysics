#ifndef ALIFLATESDVERTEX_H
#define ALIFLATESDVERTEX_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/**
 * >> Flat structure representing a ESD vertex <<
 */

#include "Rtypes.h"
#include "AliVVvertex.h"
#include "AliESDVertex.h"

class AliFlatESDVertex: public AliVVvertex
//class AliFlatESDVertex
{
  public:
  Double32_t fPosition[3];    // vertex position
  Double32_t fCov[6];  // vertex covariance matrix
  Int_t    fNContributors;  // # of tracklets/tracks used for the estimate   
  Double32_t fChi2;  // chi2 of vertex fit
  /*
    Double32_t fSigma;          // track dispersion around found vertex
    Char_t fID;       // ID of this vertex within an ESD event
    Char_t fBCID;     // BC ID assigned to vertex
  */

	AliFlatESDVertex(Bool_t){}
  virtual ~AliFlatESDVertex() {}

  AliFlatESDVertex() :fNContributors(0), fChi2(0){
    for( int i=0; i<3; i++) fPosition[i] = -9999;
    for( int i=0; i<6; i++) fCov[i] = -9999;
  }

  void Set(const AliESDVertex &v );

  Double32_t GetX() const { return fPosition[0]; }
  Double32_t GetY() const { return fPosition[1]; }
  Double32_t GetZ() const { return fPosition[2]; }
  
  void GetXYZ(Double_t pos[3]) const { for(Int_t j=0; j<3; j++) pos[j]=fPosition[j]; }
  void SetXYZ(Double_t pos[3]) { for(Int_t j=0; j<3; j++) fPosition[j]=pos[j]; }

  void   SetX(Double_t xVert) {fPosition[0]=xVert; }
  void   SetY(Double_t yVert) {fPosition[1]=yVert; }
  void   SetZ(Double_t zVert) {fPosition[2]=zVert; } 
  void   SetNContributors(Int_t nContr) {fNContributors=nContr; }


  Int_t    GetNContributors() const { return fNContributors; }

  /*
  void     GetCovarianceMatrix(Double_t covmatrix[6]) const;
  void     SetCovarianceMatrix(const Double_t *) {}
  
  Double_t GetChi2perNDF() const {return -999.;}
  Double_t GetChi2() const {return -999.;}
  void     SetChi2(Double_t ) {}
  Int_t    GetNDF() const {return -999;}

  void     GetSigmaXYZ(Double_t sigma[3]) const;
  void     GetCovMatrix(Double_t covmatrix[6]) const;
  void     GetCovarianceMatrix(Double_t covmatrix[6]) const 
                    {GetCovMatrix(covmatrix);}
  void     GetSNR(Double_t snr[3]) const;
  void     SetCovarianceMatrix(const Double_t *cov);

  Double_t GetXRes() const {return TMath::Sqrt(fCovXX);}
  Double_t GetYRes() const {return TMath::Sqrt(fCovYY);}
  Double_t GetZRes() const {return TMath::Sqrt(fCovZZ);}
  Double_t GetXSNR() const { return fSNR[0]; }
  Double_t GetYSNR() const { return fSNR[1]; }
  Double_t GetZSNR() const { return fSNR[2]; }
  void     SetSNR(double snr, int i) {if (i<3 && i>=0) fSNR[i] = snr;}

  Double_t GetChi2() const { return fChi2; }
  void     SetChi2(Double_t chi) { fChi2 = chi; }
  Double_t GetChi2toNDF() const 
    { return fChi2/(2.*(Double_t)fNContributors-3.); }
  Double_t GetChi2perNDF() const { return GetChi2toNDF();}
  Int_t    GetNDF() const {return (2*fNContributors-3);}

  void     Print(Option_t* option = "") const;
  void     PrintStatus() const {Print();}

  void     Reset() { SetToZero(); SetName("Vertex"); }

  void     SetID(Char_t id) {fID=id;}
  Char_t   GetID() const {return fID;}
  //
  Double_t GetWDist(const AliESDVertex* v) const;
  */


 
};

#endif
