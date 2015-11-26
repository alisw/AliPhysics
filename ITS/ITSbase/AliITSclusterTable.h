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
  AliITSclusterTable(Double_t x, Double_t y, Double_t z, Double_t sx, Double_t sy, Double_t sz, Double_t phi, Double_t lambda, Int_t index);
  virtual ~AliITSclusterTable(){;}

  Int_t   GetOrInd() const {return int(GetUniqueID());}
  Double_t GetX() const {return fX;}
  Double_t GetY() const {return fY;}
  Double_t GetZ() const {return fZ;}
  Double_t GetSx() const {return fSx;}
  Double_t GetSy() const {return fSy;}
  Double_t GetSz() const {return fSz;}
  Double_t GetPhi() const {return fPhi;}
  Double_t GetLambda() const {return fLam;}

  virtual Bool_t IsEqual(const TObject *obj) const 
    {return fLam == ((AliITSclusterTable*)obj)->fLam;}
  virtual Bool_t      IsSortable() const { return kTRUE; }
  virtual Int_t       Compare(const TObject *obj) const 
    {if(fLam<((AliITSclusterTable*)obj)->fLam) return -1;
    else if(fLam>((AliITSclusterTable*)obj)->fLam) return 1;
    else return 0; }

 protected: 

  Float_t fX;  //!x of cluster 
  Float_t fY;  //!y of cluster
  Float_t fZ;  //!z of cluster
  Float_t fSx; //! error on x
  Float_t fSy; //! error on y
  Float_t fSz; //! error on z
  Float_t fPhi; //! azimuthal angle
  Float_t fLam; //! lambda angle
  
  ClassDef(AliITSclusterTable,4)
};

#endif



