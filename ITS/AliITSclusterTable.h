#ifndef ALIITSCLUSTERTABLE_H
#define ALIITSCLUSTERTABLE_H 
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//////////////////////////////////////////////////////////////////////////
// Class used to simplify some operations with clusters.                 //
// -Function FillArray fills an array wich contains, for each            //
//  ITS module, an array with the indices of all the clusters detected   //
//  by the module. The indices correspond to the cluster indices in class// 
//  AliITSlayer of AliITStrackerV2.                                      //
//  This function is used in AliITStrackerSA::FindTracks.                // 
// -Function FillArrayLabel fills an array wich contains, for each       //
//  particle label, and for each layer, the information on clusters:     //
//  0 if there is no cluster, 1 if there is a cluster with this label.   //
//  This function is used to define trackable tracks.                    //   
///////////////////////////////////////////////////////////////////////////


#include <TArrayI.h>
#include <TArrayD.h>
#include <TObject.h>

class TTree;
class AliITStrackerSA;
class AliITSgeom;
class AliITSclusterV2;

class AliITSclusterTable : public TObject {


 public:

  AliITSclusterTable();
  AliITSclusterTable(AliITSgeom* geom, AliITStrackerSA* tracker, Double_t* primaryVertex);
  void FillArray(TTree* clusterTree);
  void FillArrayLabel(Int_t numberofparticles);
  void FillArrayCoorAngles();
  void GetCoorAngles(AliITSclusterV2* cl,Int_t module,Double_t &phi,Double_t &lambda,Double_t &x,Double_t &y,Double_t &z);
  void GetCoorErrors(AliITSclusterV2* cl, Int_t module,Double_t &sx,Double_t &sy, Double_t &sz);
  virtual ~AliITSclusterTable();


  Int_t      GetNCluster(Int_t mod) const {return fNCl[mod];}
  Int_t      GetClusterIndMod(Int_t mod,Int_t i){return fDet[mod]->At(i);} 
  Int_t      ThereIsClusterOnLayer(Int_t label,Int_t layer)
             {return fLbl[label]->At(layer);}
  Int_t      ThisParticleIsTrackable(Int_t label,Int_t numberofpoints=6);
  Double_t   GetPhiCluster(Int_t layer, Int_t i){return fPhiList[layer]->At(i);}
  Double_t   GetLambdaCluster(Int_t layer, Int_t i) {return fLambdaList[layer]->At(i);}
  Double_t   GetXCluster(Int_t layer, Int_t i){return fXList[layer]->At(i);}
  Double_t   GetYCluster(Int_t layer, Int_t i) {return fYList[layer]->At(i);}
  Double_t   GetZCluster(Int_t layer, Int_t i) {return fZList[layer]->At(i);}
  Double_t   GetXClusterError(Int_t layer, Int_t i) {return fSxList[layer]->At(i);}
  Double_t   GetYClusterError(Int_t layer, Int_t i) {return fSyList[layer]->At(i);}
  Double_t   GetZClusterError(Int_t layer, Int_t i) {return fSzList[layer]->At(i);}

  TArrayI*   GetListOfClusters(Int_t mod) const {return fDet[mod];}
  TArrayI*   GetNClustersSameLabel(Int_t label) const {return fLbl[label];}
  TArrayD*   GetListOfPhi(Int_t layer) const {return fPhiList[layer];}
  TArrayD*   GetListOfLambda(Int_t layer) const {return fLambdaList[layer];}
  TArrayD*   GetListOfX(Int_t layer) const {return fXList[layer];}
  TArrayD*   GetListOfY(Int_t layer) const {return fYList[layer];}
  TArrayD*   GetListOfZ(Int_t layer) const {return fZList[layer];}
  TArrayD*   GetListOfSx(Int_t layer)const {return fSxList[layer];}
  TArrayD*   GetListOfSy(Int_t layer)const {return fSyList[layer];}
  TArrayD*   GetListOfSz(Int_t layer)const {return fSzList[layer];}
 protected: 

  // copy constructor (NO copy allowed: the constructor is protected
  // to avoid misuse)
  AliITSclusterTable(const AliITSclusterTable& tab);
  // assignment operator (NO assignment allowed)
  AliITSclusterTable& operator=(const AliITSclusterTable& /* tab */);

  static Int_t FindIndex(Int_t ndim, Int_t *ptr, Int_t value);

  Int_t        *fNCl;      //number of clusters per module
  Double_t     fPrimaryVertex[3]; //primaryVertex
  TArrayI**    fDet;       //Array of cluster indices for each detector
  TArrayI**    fLbl;       //Array of number of clusters (on each layer) 
                           // with the same label for each label.
  TArrayD**    fPhiList;   //Array of cluster azimuthal angles on each layer
  TArrayD**    fLambdaList;//Array of cluster Lambda angles on each layer
  TArrayD**    fXList;     //Array of cluster x coordinates on each layer
  TArrayD**    fYList;     //Array of cluster y coordinates on each layer
  TArrayD**    fZList;    // Array of cluster z coordinates on each layer
  TArrayD**    fSxList;    //Array of cluster errors on x on each layer
  TArrayD**    fSyList;    //Array of cluster errors on y on each layer
  TArrayD**    fSzList;    //Array of cluster errors on z on each layer
 
  AliITSgeom *fGeom;      //! ITS geometry
  AliITStrackerSA *fTracker; //! SA tracker

  ClassDef(AliITSclusterTable,1)
};

#endif



