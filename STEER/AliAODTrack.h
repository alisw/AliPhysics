#ifndef AliAODTrack_H
#define AliAODTrack_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD track base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <TRef.h>

#include "AliVirtualParticle.h"
#include "AliAODVertex.h"

class AliAODTrack : public AliVirtualParticle {

 public:
  
  enum AODTrk_t {kUndef=-1, kPrimary, kSecondary, kOrphan};

  enum AODTrkBits_t {
    kIsDCA=BIT(14)   // set if fPosition is the DCA and not the position of the first point
  };

  enum AODTrkPID_t {
    kUnknown=0, kElectron, kMuon, kPion, kProton, kDeuton, kTriton, kAlpha, kOther};

  AliAODTrack();
  AliAODTrack(Int_t id,
	      Int_t label,
	      Double_t p[3],
	      Bool_t cartesian,
	      Double_t x[3],
	      Bool_t dca,
	      Double_t covMatrix[21],
	      Short_t q,
	      UChar_t itsClusMap,
	      Double_t pid[10],
	      AliAODVertex *prodVertex,
	      AODTrk_t ttype=kUndef);

   AliAODTrack(Int_t id,
	      Int_t label,
	      Float_t p[3],
	      Bool_t cartesian,
	      Float_t x[3],
	      Bool_t dca,
	      Float_t covMatrix[21],
	      Short_t q,
	      UChar_t itsClusMap,
	      Float_t pid[10],
	      AliAODVertex *prodVertex,
	      AODTrk_t ttype=kUndef);

  virtual ~AliAODTrack();
  AliAODTrack(const AliAODTrack& trk); 
  AliAODTrack& operator=(const AliAODTrack& trk);

  // kinematics
  virtual Double_t OneOverPt() const { return fMomentum[0]; }
  virtual Double_t Phi()       const { return fMomentum[1]; }
  virtual Double_t Theta()     const { return fMomentum[2]; }
  
  virtual Double_t Px() const { return TMath::Cos(fMomentum[1])/fMomentum[0]; }
  virtual Double_t Py() const { return TMath::Sin(fMomentum[1])/fMomentum[0]; }
  virtual Double_t Pz() const { return 1./(fMomentum[0] * TMath::Tan(fMomentum[2])); }
  virtual Double_t Pt() const { return 1./fMomentum[0]; }
  virtual Double_t P()  const { return TMath::Sqrt(Pt()*Pt()+Pz()*Pz()); }

          Double_t Chi2() const { return fChi2; }

  virtual Double_t E() const { return -999.; }
  // make a connection to the PID object, here!!!
  virtual Double_t M() const { return -999.; }
  
  virtual Double_t Eta() const { return -TMath::Log(TMath::Tan(0.5 * fMomentum[2])); }
  // make a connection to the PID object, here!!!
  virtual Double_t Y() const { return -999.; }

  virtual Short_t  Charge() const {return fCharge; }

  // PID
  virtual const Double_t *PID() const { return fPID; }

  template <class T> void GetPID(T *pid) const {
    for(Int_t i=0; i<10; ++i) pid[i]=fPID[i];}
 
  template <class T> void SetPID(const T *pid) {
    if(pid) for(Int_t i=0; i<10; ++i) fPID[i]=pid[i];
    else {for(Int_t i=1; i<10; fPID[i++]=0); fPID[0]=1.;}}

  Int_t GetID() const { return fID; }
  Int_t GetLabel() const { return fLabel; } 
  Char_t GetType() const { return fType;}

  template <class T> void GetP(T *p) const {
    p[0]=fMomentum[0]; p[1]=fMomentum[1]; p[2]=fMomentum[2];}

  template <class T> void GetPxPyPz(T *p) const {
    p[0] = Px(); p[1] = Py(); p[2] = Pz();}

  template <class T> Bool_t GetPosition(T *x) const {
    x[0]=fPosition[0]; x[1]=fPosition[1]; x[2]=fPosition[2];
    return TestBit(kIsDCA);}

  template <class T> void SetCovMatrix(const T *covMatrix) {
    if(!fCovMatrix) fCovMatrix=new AliAODTrkCov();
    fCovMatrix->SetCovMatrix(covMatrix);}

  template <class T> Bool_t GetCovMatrix(T *covMatrix) const {
    if(!fCovMatrix) return kFALSE;
    fCovMatrix->GetCovMatrix(covMatrix); return kTRUE;}

  void RemoveCovMatrix() {delete fCovMatrix; fCovMatrix=NULL;}

  UChar_t GetITSClusterMap() const { return fITSClusterMap; }

  AliAODVertex *GetProdVertex() const { return (AliAODVertex*)fProdVertex.GetObject(); }
  
  // print
  void  Print(const Option_t *opt = "") const;

  // setters
  void SetID(const Int_t id) { fID = id; }
  void SetLabel(const Int_t label) {fLabel = label; }

  template <class T> void SetPosition(const T *x, const Bool_t isDCA = kFALSE);
  void SetDCA(Double_t d, Double_t z);

  void SetOneOverPt(const Double_t oneOverPt) { fMomentum[0] = oneOverPt; }
  void SetPt(const Double_t pt) { fMomentum[0] = 1./pt; };
  void SetPhi(const Double_t phi) { fMomentum[1] = phi; }
  void SetTheta(const Double_t theta) { fMomentum[2] = theta; }
  template <class T> void SetP(const T *p, const Bool_t cartesian = kTRUE);
  void SetP() {fMomentum[0]=fMomentum[1]=fMomentum[2]=-999.;}

  void SetCharge(const Short_t q) { fCharge = q; }
  void SetChi2(const Double_t chi2) { fChi2 = chi2; }

  void SetITSClusterMap(const UChar_t itsClusMap) { fITSClusterMap = itsClusMap; }

  void SetProdVertex(TObject *vertex) { fProdVertex = vertex; }

  // name and title
  void SetType(AODTrk_t ttype) { fType=ttype; }


   class AliAODTrkCov {

   //
   //  Class containing the covariance matrix for the track
   //
   //       X          Y          Z         Px        Py        Pz
   //
   // X  fDiag[ 0]  
   //
   // Y  fOdia[ 0]  fDiag[ 1]
   //
   // Z  fOdia[ 1]  fOdia[ 2]  fDiag[ 2]
   //
   // Px fOdia[ 3]  fOdia[ 4]  fOdia[ 5]  fDiag[ 3]
   //
   // Py fOdia[ 6]  fOdia[ 7]  fOdia[ 8]  fOdia[ 9]  fDiag[ 4]
   //
   // Pz fOdia[10]  fOdia[11]  fOdia[12]  fOdia[13]  fOdia[14]  fDiag[ 5]
   //

   public:
   AliAODTrkCov() {}
   virtual ~AliAODTrkCov() {}
   template <class T> void GetCovMatrix(T *cmat) const;
   template <class T> void SetCovMatrix(T *cmat);

   private:
   Double32_t   fDiag[6];  // Diagonal elements
   Double32_t   fODia[15]; // [-1, 1,8] 8 bit precision for off diagonal elements

   ClassDef(AliAODTrack::AliAODTrkCov,1)

 };

 private :

  // Momentum & position
  Double32_t    fMomentum[3];    // momemtum stored in 1/pt, phi, theta
  Double32_t    fPosition[3];    // position of first point on track or dca

  Double32_t    fPID[10];        // [0.,1.,8] pointer to PID object
  Double32_t    fChi2;           // chi2 of mometum fit

  Int_t         fID;             // unique track ID, points back to the ESD track
  Int_t         fLabel;          // track label, points back to MC track
  
  AliAODTrkCov *fCovMatrix;      // covariance matrix (x, y, z, px, py, pz)
  TRef          fProdVertex;     // vertex of origin

  Char_t        fCharge;         // particle charge
  UChar_t       fITSClusterMap;  // map of ITS cluster, one bit per layer
  Char_t        fType;           // Track Type


  ClassDef(AliAODTrack,1);
};

#endif
