#ifndef ALIITSTRACKU_H
#define ALIITSTRACKU_H 
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
////////////////////////////////////////////////////////
//  Stand alone track class UPGRADE                   //
//  Authors: A.Mastroserio                            //
//           C.Terrevoli                              //
//           annalisa.mastroserio@cern.ch             //      
//           cristina.terrevoli@ba.infn.it            //                                                    
////////////////////////////////////////////////////////

#include "AliITStrackV2.h"
#include "AliCluster.h"


class AliESDtrack;

class AliITStrackU : public AliITStrackV2 {


 public:

  AliITStrackU();
  AliITStrackU(Int_t nlay);
  AliITStrackU(AliESDtrack& t,Bool_t c=kFALSE);
  AliITStrackU(const AliITStrackU& t, Bool_t trackMI=kFALSE);
  AliITStrackU(Double_t alpha, Double_t radius,
	       Double_t Ycoor, Double_t Zcoor, Double_t phi, 
	       Double_t tanlambda, Double_t curv, Int_t lab, Int_t nlay);
  Bool_t UpdateMI(const AliCluster *c, Double_t chi2, Int_t i);
  Int_t* ClIndex() {return fClIndex;}
  Float_t GetSigmaY(Int_t i) const {return fSigmaY[i];}
  Float_t GetSigmaZ(Int_t i) const {return fSigmaZ[i];}
  Float_t GetSigmaYZ(Int_t i) const {return fSigmaYZ[i];}
  void SetSigmaY(Int_t i, Float_t s) {fSigmaY[i]=s;}
  void SetSigmaZ(Int_t i, Float_t s) {fSigmaZ[i]=s;}
  void SetSigmaYZ(Int_t i, Float_t s) {fSigmaYZ[i]=s;}
  
  Float_t GetExpQ() const {return fExpQ;}
  void SetExpQ(Float_t f) {fExpQ=f;}

  
  Int_t GetClIndex(Int_t i) const {return fClIndex[i];}
  void SetClIndex(Int_t i, Int_t c) {fClIndex[i]=c;}
  Float_t GetNormChi2(Int_t i) const {return fNormChi2[i];}
  void SetNormChi2(Int_t i, Float_t n) {fNormChi2[i]=n;}
  Float_t GetNormQ(Int_t i) const {return fNormQ[i];}
  void SetNormQ(Int_t i, Float_t q) {fNormQ[i]=q;}
  Float_t GetNy(Int_t i) const {return fNy[i];}
  void SetNy(Int_t i, Float_t f) {fNy[i]=f;}
  Float_t GetNz(Int_t i) const {return fNz[i];}
  void SetNz(Int_t i, Float_t f) {fNz[i]=f;}
  Double_t GetPredictedChi2MI(Double_t cy, Double_t cz, Double_t cerry, Double_t cerrz, Double_t covyz=0.) const;

  Int_t GetClusterIndexU(Int_t i) const {return fSain[i];}
  Int_t GetClusterMark(Int_t layer,Int_t i) const {return fCluMark[layer][i];}
  Int_t GetNumberOfClustersU() const {return fNU;}
  Int_t GetNumberOfMarked(Int_t lay) const {return fNM[lay];}
  static Int_t GetMaxNumberOfClusters() {return kMaxNumberOfClusters;}
  Int_t GetMaxNMarkedPerLayer() const {return kMaxNumberOfClustersL;}
  void  AddClusterU(Int_t layer, Int_t clnumb);
  void  AddClusterV2(Int_t layer,Int_t clnumb);
  void  AddClusterMark(Int_t layer, Int_t clnumb);

  enum {kMaxNumberOfClustersL = 4};// Max. n. of clusters/layer 
  enum {kMaxNumberOfClusters = 15};// Max. number of clusters

 protected: 
  
  void Init(Double_t alpha, Double_t radius,
	    Double_t Ycoor, Double_t Zcoor, Double_t phi, 
	    Double_t tanlambda, Double_t curv, Int_t lab/*,Int_t nlay*/);
  void SetNumberOfClustersU(Int_t n){fNU = n;}
  void SetNumberOfMarked(Int_t lay,Int_t n) {fNM[lay] = n;}
  void ResetIndexU(){for(Int_t k=0; k<kMaxNumberOfClusters; k++) fSain[k]=0;}
  void ResetMarked(); 

  static const Int_t fgMaxNLayer = 8; //max number of layers in ITSUpgrade
  Int_t fNLayers;    // number of layers
  UInt_t  fSain[kMaxNumberOfClusters];   // cluster index (Stand Alone Upgrade)
  Int_t fNU;          // number of clusters Stand Alone Upgrade 

  Int_t fCluMark[fgMaxNLayer][kMaxNumberOfClustersL]; //indices for cluster used
  Int_t fNM[fgMaxNLayer]; //number of marked clusters
  
  Float_t fDy[fgMaxNLayer];           //dy in layer
  Float_t fDz[fgMaxNLayer];           //dz in layer
  Float_t fSigmaY[fgMaxNLayer];       //sigma y 
  Float_t fSigmaZ[fgMaxNLayer];       //sigma z
  Float_t fSigmaYZ[fgMaxNLayer];       //covariance of y and z

  Float_t fNy[fgMaxNLayer];              //expected size of cluster
  Float_t fNz[fgMaxNLayer];              //expected size of cluster
  Float_t fNormQ[fgMaxNLayer];           // normalized Q
  Float_t fExpQ;            // expected Q

  Float_t fNormChi2[fgMaxNLayer];        // normalized chi2
  Int_t    fClIndex[fgMaxNLayer];        //cluster Index
  
  
  ClassDef(AliITStrackU,2)
    };

#endif



