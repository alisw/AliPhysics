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

#include "AliITStrackMI.h"

class AliITStrackSA : public AliITStrackMI {


 public:

  AliITStrackSA();
  AliITStrackSA(const AliITStrackMI& t);
  AliITStrackSA(const AliITStrackSA& t);
  AliITStrackSA(AliITSgeom* geom,Int_t layer, Int_t ladder, Int_t detector, 
                Double_t Ycoor, Double_t Zcoor, Double_t phi, 
                Double_t tanlambda, Double_t curv, Int_t lab);


  Int_t GetClusterIndexSA(Int_t i) const {return fSain[i];}
  Int_t GetNumberOfClustersSA() const {return fNSA;}
  void  AddClusterSA(Int_t layer, Int_t clnumb);
  void  AddClusterV2(Int_t layer,Int_t clnumb);

 protected: 

  void SetNumberOfClustersSA(Int_t n){fNSA = n;}
  void ResetIndexSA(){for(Int_t k=0; k<fgkMaxNumberOfClusters; k++) fSain[k]=0;}
  static const Int_t fgkMaxNumberOfClusters = 20; // Max. number of clusters 
                                            // per trackSA
  UInt_t  fSain[fgkMaxNumberOfClusters];   // cluster index (SA)
  Int_t fNSA;          // number of clusters SA 
 
  ClassDef(AliITStrackSA,2)
};

#endif



