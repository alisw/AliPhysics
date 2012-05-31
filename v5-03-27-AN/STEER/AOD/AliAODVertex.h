#ifndef ALIAODVERTEX_H
#define ALIAODVERTEX_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD vertex base class
//     Author: Markus Oldenburg, CERN
//     Inheritance from AliVVertex: A. Dainese
//-------------------------------------------------------------------------

#include <TRef.h>
#include <TRefArray.h>
#include <TMath.h>

#include "AliVVertex.h"
#include "AliAODRedCov.h"
#include "AliLog.h"

class AliAODVertex : public AliVVertex {

 public :

  enum AODVtx_t {kUndef=-1, kPrimary, kKink, kV0, kCascade, kMulti, kMainSPD, kPileupSPD, kPileupTracks,kMainTPC};

  AliAODVertex();
  AliAODVertex(const Double_t *position, 
	       const Double_t *covMatrix=0x0,
	       Double_t chi2perNDF = -999.,
	       TObject *parent = 0x0,
	       Short_t id=-1,
	       Char_t vtype=kUndef,
	       Int_t nprong = 0);
  AliAODVertex(const Float_t *position, 
	       const Float_t *covMatrix=0x0,
	       Double_t chi2perNDF = -999.,
	       TObject *parent = 0x0,
	       Short_t id=-1,
	       Char_t vtype=kUndef,
	       Int_t nprong = 0);
    AliAODVertex(const Double_t *position, 
		 Double_t chi2perNDF,
		 Char_t vtype=kUndef,
		 Int_t nprong = 0);
    AliAODVertex(const Float_t *position, 
		 Double_t chi2perNDF,
		 Char_t vtype=kUndef,
		 Int_t nprong = 0);

  virtual ~AliAODVertex();
  AliAODVertex(const AliAODVertex& vtx); 
  AliAODVertex& operator=(const AliAODVertex& vtx);

  virtual AliAODVertex* CloneWithoutRefs() const;
  
  void     SetX(Double_t x) { fPosition[0] = x; }
  void     SetY(Double_t y) { fPosition[1] = y; }
  void     SetZ(Double_t z) { fPosition[2] = z; }
  void     SetPosition(Double_t x, Double_t y, Double_t z) { fPosition[0] = x; fPosition[1] = y; fPosition[2] = z; }
  template <class T> void SetPosition(T *pos)
    { fPosition[0] = pos[0]; fPosition[1] = pos[1]; fPosition[2] = pos[2]; }

  void     SetChi2perNDF(Double_t chi2perNDF) { fChi2perNDF = chi2perNDF; }

  void     SetParent(TObject *parent) { fParent = parent; }

  Double_t GetX() const { return fPosition[0]; }
  Double_t GetY() const { return fPosition[1]; }
  Double_t GetZ() const { return fPosition[2]; }
  void     GetXYZ(Double_t position[3]) const 
    {position[0]=fPosition[0]; position[1]=fPosition[1]; position[2]=fPosition[2];}
  template <class T> void GetPosition(T *pos) const
    {pos[0]=fPosition[0]; pos[1]=fPosition[1]; pos[2]=fPosition[2];}

  template <class T> void SetCovMatrix(const T *covMatrix) {
    if(!fCovMatrix) fCovMatrix=new AliAODRedCov<3>();
    fCovMatrix->SetCovMatrix(covMatrix);}

  template <class T> Bool_t GetCovMatrix(T *covMatrix) const {
    if(!fCovMatrix) return kFALSE;
    fCovMatrix->GetCovMatrix(covMatrix); return kTRUE;}

  void GetCovarianceMatrix(Double_t covmatrix[6]) const 
                    {GetCovMatrix(covmatrix);}
  void RemoveCovMatrix() {delete fCovMatrix; fCovMatrix=NULL;}

  template <class T> void     GetSigmaXYZ(T *sigma) const;

  Double_t GetChi2perNDF() const { return fChi2perNDF; }
  Double_t GetChi2() const { return fChi2perNDF*(Double_t)GetNDF(); }
  Int_t    GetNDF() const { return 2*GetNContributors()-3; }

  Short_t  GetID() const { return fID; }
  void     SetID(Short_t id) { fID=id; }

  Char_t   GetType() const { return fType; }
  void     SetType(AODVtx_t vtype) { fType=vtype; }

  TObject* GetParent() const   { return fParent.GetObject(); }
  Bool_t   HasParent(TObject *parent) const { return (fParent.GetObject() == parent) ? kTRUE : kFALSE; }

  void     AddDaughter(TObject *daughter);
  void     RemoveDaughter(TObject *daughter) { fDaughters.Remove(daughter); }
  void     RemoveDaughters() { fDaughters.Clear(); if(fProngs) {delete [] fProngs; fProngs=0; MakeProngs(); fIprong=0;} }
  TObject* GetDaughter(Int_t i); 
  Bool_t   HasDaughter(TObject *daughter) const;
  Int_t    GetNDaughters() const;
  Int_t    GetNContributors() const;
  void     SetNContributors(Int_t nc) {fNContributors = nc;}
  // covariance matrix elements after rotation by phi around z-axis 
  // and, then, by theta around new y-axis
  Double_t  RotatedCovMatrixXX(Double_t phi = 0., Double_t theta = 0.) const;
  Double_t  RotatedCovMatrixXY(Double_t phi = 0., Double_t theta = 0.) const;
  Double_t  RotatedCovMatrixYY(Double_t phi = 0.) const;
  Double_t  RotatedCovMatrixXZ(Double_t phi = 0., Double_t theta = 0.) const;
  Double_t  RotatedCovMatrixYZ(Double_t phi = 0., Double_t theta = 0.) const;
  Double_t  RotatedCovMatrixZZ(Double_t phi = 0., Double_t theta = 0.) const;

  template <class T, class P> void     PhiAndThetaToVertex(AliAODVertex *vtx, P &phi, T &theta) const;
  Double_t  Distance2ToVertex(const AliAODVertex *vtx) const;
  Double_t  DistanceToVertex(AliAODVertex *vtx) const 
     {return TMath::Sqrt(Distance2ToVertex(vtx));}
  Double_t  DistanceXY2ToVertex(const AliAODVertex *vtx) const;
  Double_t  DistanceXYToVertex(AliAODVertex *vtx) const 
     {return TMath::Sqrt(DistanceXY2ToVertex(vtx));}
  Double_t  Error2DistanceToVertex(AliAODVertex *vtx) const;
  Double_t  ErrorDistanceToVertex(AliAODVertex *vtx) const
     {return TMath::Sqrt(Error2DistanceToVertex(vtx));}
  Double_t  Error2DistanceXYToVertex(AliAODVertex *vtx) const;
  Double_t  ErrorDistanceXYToVertex(AliAODVertex *vtx) const
     {return TMath::Sqrt(Error2DistanceXYToVertex(vtx));}
  
  void     PrintIndices() const;
  void     Print(Option_t* option = "") const;

  const char* AsString() const;
  
  static const char* GetTypeName(AODVtx_t type);
  void     SetBC(Int_t bc)               {fBCID = bc;}
  Int_t    GetBC()              const    {return fBCID;}  
private:
  void     MakeProngs() {if (fNprong > 0) {fProngs = new TRef[fNprong]; fIprong=0;}}
	  
 private:

  Double32_t      fPosition[3];   // vertex position
  Double32_t      fChi2perNDF;    // chi2/NDF of vertex fit
  Short_t         fID;            // vertex ID; corresponds to the array index of the appropriate ESD container
  Char_t          fBCID;          // BC ID assigned to vertex
  Char_t          fType;          // vertex type
  Int_t           fNprong;        // number of prongs
  Int_t           fIprong;        //!index  of prong
  Int_t           fNContributors; // Number of contributors for SPD vertex
  AliAODRedCov<3> *fCovMatrix;    // vertex covariance matrix; values of and below the diagonal
  TRef            fParent;        // reference to the parent particle
  TRefArray       fDaughters;     // references to the daughter particles
  TRef            *fProngs;       //[fNprong] alternative daughters for n-prong vertex
  
  ClassDef(AliAODVertex, 8);
};

inline  Int_t AliAODVertex::GetNDaughters() const
{
    if (!fProngs) {
	return fDaughters.GetEntriesFast();
    } else {
	return fNprong;
    }
}

inline TObject* AliAODVertex::GetDaughter(Int_t i)
{
    if (!fProngs) {
	return fDaughters.At(i);
    } else {
	if (i < fNprong) {
	    return fProngs[i].GetObject();
	} else {
	    AliWarning("Daughter index out of range !\n");
	    return 0;
	}
    }
}

#endif
