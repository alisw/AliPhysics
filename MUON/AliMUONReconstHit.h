#ifndef ALIMUONRECONSTHIT_H
#define ALIMUONRECONSTHIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>

class AliMUONReconstHit : public TObject {
public:

  // correlation starts from the 1-st cathode  
  // last number in arrays corresponds to cluster on 1-st cathode

   Int_t       fCorrelIndex[4];  // entry number in TreeR for the associated 
                                 // cluster candidates on the 2-nd cathode
   Float_t     fX[4]  ;          // X of clusters on the 2-nd cathode  
   Float_t     fY[4]  ;          // Y of clusters

public:
   AliMUONReconstHit() {
       fCorrelIndex[0]=fCorrelIndex[1]=fCorrelIndex[2]=fCorrelIndex[3]=0;
       fX[0]=fX[1]=fX[2]=fX[3]=0; fY[0]=fY[1]=fY[2]=fY[3]=0; 
   }
   AliMUONReconstHit(Int_t *idx, Float_t *x, Float_t *y);
   virtual ~AliMUONReconstHit() {}
   ClassDef(AliMUONReconstHit,1)  // Reconstructed Hit Object for set:MUON
};
#endif
