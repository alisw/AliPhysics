#ifndef ALIITSTRACKSA_H
#define ALIITSTRACKSA_H 
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
////////////////////////////////////////////////////
//  Stand alone track class                       //
//  Origin:  Elisabetta Crescio                   //
//  e-mail:  crescio@to.infn.it                   //
//                                                //
////////////////////////////////////////////////////

/* $Id$ */

#include "AliITStrackMI.h"

class AliITStrackSA : public AliITStrackMI {


 public:

  AliITStrackSA();
  AliITStrackSA(const AliITStrackMI& t);
  AliITStrackSA(const AliITStrackSA& t);
  AliITStrackSA(Int_t layer, Int_t ladder, Int_t detector, 
                Double_t Ycoor, Double_t Zcoor, Double_t phi, 
                Double_t tanlambda, Double_t curv, Int_t lab);
  AliITStrackSA(Double_t alpha, Double_t radius,
                Double_t Ycoor, Double_t Zcoor, Double_t phi, 
                Double_t tanlambda, Double_t curv, Int_t lab);

  Int_t GetClusterIndexSA(Int_t i) const {return fSain[i];}
  Int_t GetClusterMark(Int_t layer,Int_t i) const {return fCluMark[layer][i];}
  Int_t GetNumberOfClustersSA() const {return fNSA;}
  Int_t GetNumberOfMarked(Int_t lay) const {return fNM[lay];}
  static Int_t GetMaxNumberOfClusters() {return kMaxNumberOfClusters;}
  Int_t GetMaxNMarkedPerLayer() const {return kMaxNumberOfClustersL;}
  void  AddClusterSA(Int_t layer, Int_t clnumb);
  void  AddClusterV2(Int_t layer,Int_t clnumb);
  void  AddClusterMark(Int_t layer, Int_t clnumb);

  enum {kMaxNumberOfClustersL = 4};// Max. n. of clusters/layer 
  enum {kMaxNumberOfClusters = 15};// Max. number of clusters

 protected: 
  AliITStrackSA& operator=(const AliITStrackSA& t); 
  void Init(Double_t alpha, Double_t radius,
	    Double_t Ycoor, Double_t Zcoor, Double_t phi, 
	    Double_t tanlambda, Double_t curv, Int_t lab);
  void SetNumberOfClustersSA(Int_t n){fNSA = n;}
  void SetNumberOfMarked(Int_t lay,Int_t n) {fNM[lay] = n;}
  void ResetIndexSA(){for(Int_t k=0; k<kMaxNumberOfClusters; k++) fSain[k]=0;}
  void ResetMarked(); 


  UInt_t  fSain[kMaxNumberOfClusters];   // cluster index (SA)
  Int_t fNSA;          // number of clusters SA 

  Int_t fCluMark[AliITSgeomTGeo::kNLayers][kMaxNumberOfClustersL]; //indices for cluster used
  Int_t fNM[AliITSgeomTGeo::kNLayers]; //number of marked clusters

  ClassDef(AliITStrackSA,4)
};

#endif



