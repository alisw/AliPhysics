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
#include <TObject.h>

class TTree;
class AliITStrackerSA;
class AliITSgeom;

class AliITSclusterTable : public TObject {


 public:

  AliITSclusterTable();
  AliITSclusterTable(AliITSgeom* geom, AliITStrackerSA* tracker);
  void FillArray(TTree* clusterTree,Int_t evnumber=0);
  void FillArrayLabel(Int_t numberofparticles,TTree* clusterTree,
		      Int_t evnumber=0);
  virtual ~AliITSclusterTable();


  Int_t      GetNCluster(Int_t mod) const {return fNCl[mod];}
  Int_t      GetClusterIndMod(Int_t mod,Int_t i){return fDet[mod]->At(i);} 
  TArrayI*   GetListOfClusters(Int_t mod) const {return fDet[mod];}
  TArrayI*   GetNClustersSameLabel(Int_t label) const {return fLbl[label];}
  Int_t      ThereIsClusterOnLayer(Int_t label,Int_t layer)
             {return fLbl[label]->At(layer);}
  Int_t      ThisParticleIsTrackable(Int_t label,Int_t numberofpoints=6);

 protected: 

  // copy constructor (NO copy allowed: the constructor is protected
  // to avoid misuse)
  AliITSclusterTable(const AliITSclusterTable& tab);
  // assignment operator (NO assignment allowed)
  AliITSclusterTable& operator=(const AliITSclusterTable& /* tab */);

  static Int_t FindIndex(Int_t ndim, Int_t *ptr, Int_t value);

  Int_t        *fNCl;//number of clusters per module
  TArrayI**    fDet;      //Array of cluster indices for each detector
  TArrayI**    fLbl;      //Array of number of clusters (on each layer) 
                          // with the same label for each label.
  AliITSgeom *fGeom;      //! ITS geometry
  AliITStrackerSA *fTracker; //! SA tracker

  ClassDef(AliITSclusterTable,1)
};

#endif



