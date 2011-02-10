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

#include "AliITStrackMI.h"

class AliITStrackU : public AliITStrackMI {


 public:

  AliITStrackU();
  AliITStrackU(Int_t nlay);
  AliITStrackU(const AliITStrackMI& t);
  AliITStrackU(const AliITStrackU& t);
  AliITStrackU(Double_t alpha, Double_t radius,
                Double_t Ycoor, Double_t Zcoor, Double_t phi, 
                Double_t tanlambda, Double_t curv, Int_t lab, Int_t nlay);

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
	    Double_t tanlambda, Double_t curv, Int_t lab);
  void SetNumberOfClustersU(Int_t n){fNU = n;}
  void SetNumberOfMarked(Int_t lay,Int_t n) {fNM[lay] = n;}
  void ResetIndexU(){for(Int_t k=0; k<kMaxNumberOfClusters; k++) fSain[k]=0;}
  void ResetMarked(); 

  static const Int_t fgMaxNLayer = 8; //max number of layers in ITSUpgrade
  Int_t fNU;          // number of clusters Stand Alone Upgrade 
  Int_t fNLayers;    // number of layers

  UInt_t  fSain[kMaxNumberOfClusters];   // cluster index (Stand Alone Upgrade)

  Int_t fCluMark[fgMaxNLayer][kMaxNumberOfClustersL]; //indices for cluster used
  Int_t fNM[fgMaxNLayer]; //number of marked clusters

  ClassDef(AliITStrackU,1)
};

#endif



