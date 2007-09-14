#ifndef AliAODCluster_H
#define AliAODCluster_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD cluster base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <TRef.h>

#include "AliVParticle.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"

class AliAODCluster : public AliVParticle {

 public:
  
  enum AODClu_t {kUndef = -1, 
		 kPHOSNeutral, 
		 kPHOSCharged,
		 kEMCALPseudoCluster, 
		 kEMCALClusterv1,
		 kPMDNeutral, 
		 kPMDCharged};

  enum AODCluPID_t {
    kUnknown = 0, 
    kPhoton  = 1, 
    kPi0     = 2, 
    kNeutron = 3, 
    kKaon0   = 4,
    kEleCon  = 5, 
    kCharged = 6, 
    kNeutral = 7 , 
    kOther   = 8};

  AliAODCluster();
  AliAODCluster(Int_t id,
		Int_t label,
		Double_t energy,
		Double_t x[3],
		Double_t covMatrix[10],
		Double_t pid[9],
		AliAODVertex *prodVertex, // not necessary for PMD
		AliAODTrack *primTrack,
		Char_t ttype=kUndef);

   AliAODCluster(Int_t id,
		 Int_t label,
		 Float_t energy,
		 Float_t x[3],
		 Float_t covMatrix[10],
		 Float_t pid[9],
		 AliAODVertex *prodVertex,
		 AliAODTrack *primTrack,
		 Char_t ttype=kUndef);

  virtual ~AliAODCluster();
  AliAODCluster(const AliAODCluster& trk); 
  AliAODCluster& operator=(const AliAODCluster& trk);

  Double_t Chi2() const { return fChi2; }

  virtual Double_t E() const { return fEnergy; }
  // make a connection to the PID object, here!!!
  virtual Double_t M() const { return (fType==AliAODCluster::kPHOSNeutral) ? 0. : -999.; }
  
  // make a connection to the PID object, here!!!
  virtual Double_t Y() const { return (fType==AliAODCluster::kPHOSNeutral) ? Eta() : -999.; }

  // PID
  virtual const Double_t *PID() const { return fPID; }
  AODCluPID_t GetMostProbablePID() const;
 
  template <class T> void GetPID(T *pid) const {
    for(Int_t i=0; i<9; ++i) pid[i]=fPID[i];}
 
  template <class T> void SetPID(const T *pid) {
    if(pid) for(Int_t i=0; i<9; ++i) fPID[i]=pid[i];
    else {for(Int_t i=0; i<9; fPID[i++]=0);} fPID[AliAODCluster::kUnknown]=1.;}

  Int_t GetID() const { return fID; }
  Int_t GetLabel() const { return fLabel; } 

  template <class T> Bool_t GetPosition(T *x) const {
    x[0]=fPosition[0]; x[1]=fPosition[1]; x[2]=fPosition[2];
    return kTRUE;}

  template <class T> void SetCovMatrix(const T *covMatrix) {
    if(!fCovMatrix) fCovMatrix=new AliAODRedCov<4>();
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
  void SetID(Int_t id) { fID = id; }
  void SetLabel(Int_t label) {fLabel = label; }

  template <class T> void SetPosition(const T *x);

  void SetChi2(Double_t chi2) { fChi2 = chi2; }

  void SetProdVertex(TObject *vertex) { fProdVertex = vertex; }
  void SetPrimTrack(TObject *ptrack) { fPrimTrack = ptrack; }

  virtual Double_t Px() const {return (fType==AliAODCluster::kPHOSNeutral) ? fEnergy*fPosition[0] : 0.;}
  virtual Double_t Py() const {return (fType==AliAODCluster::kPHOSNeutral) ? fEnergy*fPosition[1] : 0.;}
  virtual Double_t Pz() const {return (fType==AliAODCluster::kPHOSNeutral) ? fEnergy*fPosition[2] : 0.;}
  virtual Double_t Pt() const {return TMath::Sqrt(Px()*Px() + Py()*Py()); }
  virtual Double_t P() const {return TMath::Sqrt(Px()*Px() + Py()*Py() + Pz()*Pz()); }
  virtual Double_t OneOverPt() const {return Pt() ? 1./Pt(): -999.;}
  virtual Double_t Phi() const {return (fType==AliAODCluster::kPHOSNeutral) ? TMath::Pi()+TMath::ATan2(-fPosition[1], -fPosition[0]) : 0.;}
  virtual Double_t Theta() const {return (fType==AliAODCluster::kPHOSNeutral) ? TMath::ATan2(Pt(), fPosition[2]) : 0.;}
  virtual Double_t Eta() const {return (fType==AliAODCluster::kPHOSNeutral) ? -TMath::Log(TMath::Tan(0.5 * Theta())) : 0.;}
  virtual Short_t  Charge() const {return (fType==AliAODCluster::kPHOSNeutral) ? 0 : -999;}

 private :

  // Energy & position
  Double32_t    fEnergy;         // energy
  Double32_t    fPosition[3];    // position of the cluster

  Double32_t    fChi2;           // chi2 (probably not necessary for PMD)
  Double32_t    fPID[9];         // [0.,1.,8] pointer to PID object

  Int_t         fID;             // unique cluster ID, points back to the ESD cluster
  Int_t         fLabel;          // particle label, points back to MC track
  
  Char_t        fType;           // cluster type

  AliAODRedCov<4> *fCovMatrix;   // covariance matrix (x, y, z, E)
  TRef          fProdVertex;     // vertex of origin (not necessary for PMD)
  TRef          fPrimTrack;      // primary track associated with this cluster (not necessary for PMD)

  // TRef      fAssocCluster;       // for PMD: cluster of other layer associated with this cluster

  ClassDef(AliAODCluster,3);
};

#endif
