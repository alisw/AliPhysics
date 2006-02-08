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


#include <TObject.h>


class AliITSclusterTable : public TObject {


 public:

  AliITSclusterTable();
  AliITSclusterTable(Float_t x, Float_t y, Float_t z, Float_t sx, Float_t sy, Float_t sz, Double_t phi, Double_t lambda, Int_t index);
  virtual ~AliITSclusterTable(){;}

  Int_t   GetOrInd() const {return fOrInd;}
  Float_t GetX() const {return fX;}
  Float_t GetY() const {return fY;}
  Float_t GetZ() const {return fZ;}
  Float_t GetSx() const {return fSx;}
  Float_t GetSy() const {return fSy;}
  Float_t GetSz() const {return fSz;}
  Float_t GetPhi() const {return fPhi;}
  Float_t GetLambda() const {return fLam;}

 protected: 

  Int_t   fOrInd; //! original index in tracker
  Float_t fX;  //!x of cluster 
  Float_t fY;  //!y of cluster
  Float_t fZ;  //!z of cluster
  Float_t fSx; //! error on x
  Float_t fSy; //! error on y
  Float_t fSz; //! error on z
  Double_t fPhi; //! azimuthal angle
  Double_t fLam; //! lambda angle

  ClassDef(AliITSclusterTable,2)
};

#endif



