#ifndef AliAODNeutral_H
#define AliAODNeutral_H
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
#include "AliAODTrack.h"

class AliAODNeutral : public AliVirtualParticle {

 public:
  
  enum AODTrk_t {kUndef=-1, kPrimary, kSecondary, kOrphan};

  enum AODTrkBits_t {
    kIsDCA=BIT(14)   // set if fPosition is the DCA and not the position of the first point
  };

  enum AODNeuTrkPID_t {
    kUnknown=0, kPhoton, kPi0, kNeutron, kKaon0, kEleCon, kOther};

  AliAODNeutral();
  AliAODNeutral(Int_t id,
		Int_t label,
		Double_t energy,
		Double_t x[3],
		Bool_t dca,
		Double_t covMatrix[21],
		Double_t pid[10],
		AliAODVertex *prodVertex,
		AliAODTrack *primTrack,
		Char_t ttype=kUndef);

   AliAODNeutral(Int_t id,
		 Int_t label,
		 Float_t energy,
		 Float_t x[3],
		 Bool_t dca,
		 Float_t covMatrix[21],
		 Float_t pid[10],
		 AliAODVertex *prodVertex,
		 AliAODTrack *primTrack,
		 Char_t ttype=kUndef);

  virtual ~AliAODNeutral();
  AliAODNeutral(const AliAODNeutral& trk); 
  AliAODNeutral& operator=(const AliAODNeutral& trk);

  Double_t Chi2() const { return fChi2; }

  virtual Double_t E() const { return fEnergy; }
  // make a connection to the PID object, here!!!
  virtual Double_t M() const { return -999.; }
  
  // make a connection to the PID object, here!!!
  virtual Double_t Y() const { return -999.; }

  // PID
  virtual const Double_t *PID() const { return fPID; }

  template <class T> void GetPID(T *pid) const {
    for(Int_t i=0; i<10; ++i) pid[i]=fPID[i];}
 
  template <class T> void SetPID(const T *pid) {
    if(pid) for(Int_t i=0; i<10; ++i) fPID[i]=pid[i];
    else {for(Int_t i=1; i<10; fPID[i++]=0); fPID[0]=1.;}}

  Int_t GetID() const { return fID; }
  Int_t GetLabel() const { return fLabel; } 

  template <class T> Bool_t GetPosition(T *x) const {
    x[0]=fPosition[0]; x[1]=fPosition[1]; x[2]=fPosition[2];
    return TestBit(kIsDCA);}

  template <class T> void SetCovMatrix(const T *covMatrix) {
    if(!fCovMatrix) fCovMatrix=new AliAODNeuCov();
    fCovMatrix->SetCovMatrix(covMatrix);}

  template <class T> Bool_t GetCovMatrix(T *covMatrix) const {
    if(!fCovMatrix) return kFALSE;
    fCovMatrix->GetCovMatrix(covMatrix); return kTRUE;}

  void RemoveCovMatrix() {delete fCovMatrix; fCovMatrix=NULL;}

  AliAODVertex *GetProdVertex() const { return (AliAODVertex*)fProdVertex.GetObject(); }
  AliAODTrack *GetPrimTrack() const { return (AliAODTrack*)fPrimTrack.GetObject(); }
  
  // print
  void  Print(const Option_t *opt = "") const;

  // setters
  void SetID(const Int_t id) { fID = id; }
  void SetLabel(const Int_t label) {fLabel = label; }

  template <class T> void SetPosition(const T *x, const Bool_t isDCA = kFALSE);
  void SetDCA(Double_t d, Double_t z);

  void SetChi2(const Double_t chi2) { fChi2 = chi2; }

  void SetProdVertex(TObject *vertex) { fProdVertex = vertex; }
  void SetPrimTrack(TObject *ptrack) { fPrimTrack = ptrack; }

  virtual Double_t Px() const {return 0.;}
  virtual Double_t Py() const {return 0.;}
  virtual Double_t Pz() const {return 0.;}
  virtual Double_t Pt() const {return 0.;}
  virtual Double_t P() const {return 0.;}
  virtual Double_t OneOverPt() const {return 0.;}
  virtual Double_t Phi() const {return 0.;}
  virtual Double_t Theta() const {return 0.;}
  virtual Double_t Eta() const {return 0.;}
  virtual Short_t Charge() const {return 0.;}

   class AliAODNeuCov {

   //
   //  Class containing the covariance matrix for the cluster
   //
   //       X          Y          Z         E
   //
   // X  fDiag[ 0]  
   //
   // Y  fOdia[ 0]  fDiag[ 1]
   //
   // Z  fOdia[ 1]  fOdia[ 2]  fDiag[ 2]
   //
   // E  fOdia[ 3]  fOdia[ 4]  fOdia[ 5]  fDiag[ 3]
   //

   public:
   AliAODNeuCov() {}
   virtual ~AliAODNeuCov() {}
   template <class T> void GetCovMatrix(T *cmat) const;
   template <class T> void SetCovMatrix(T *cmat);

   private:
   Double32_t   fDiag[4];  // Diagonal elements
   Double32_t   fODia[6];  // [-1, 1,8] 8 bit precision for off diagonal elements

   ClassDef(AliAODNeutral::AliAODNeuCov,1)

 };

 private :

  // Energy & position
  Double32_t    fEnergy;         // energy
  Double32_t    fPosition[3];    // position of first point on track or dca

  Double32_t    fPID[10];        // [0.,1.,8] pointer to PID object
  Double32_t    fChi2;           // chi2 of mometum fit

  Int_t         fID;             // unique track ID, points back to the ESD track
  Int_t         fLabel;          // particle label, points back to MC track
  
  AliAODNeuCov *fCovMatrix;      // covariance matrix (x, y, z, E)
  TRef          fProdVertex;     // vertex of origin
  TRef          fPrimTrack;      // primary track number associated with this cluster

  UChar_t       fType;


  ClassDef(AliAODNeutral,1);
};

#endif
