#ifndef ALIMUONRECONSTHIT_H
#define ALIMUONRECONSTHIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

#include <TObject.h>

class AliMUONReconstHit : public TObject 
{
  public:
    AliMUONReconstHit();
    AliMUONReconstHit(Int_t *idx, Float_t *x, Float_t *y);
    virtual ~AliMUONReconstHit() {}
   
  private:

    // correlation starts from the 1-st cathode  
    // last number in arrays corresponds to cluster on 1-st cathode

    Int_t       fCorrelIndex[4];  // entry number in TreeR for the associated 
                                 // cluster candidates on the 2-nd cathode
    Float_t     fX[4]  ;          // X of clusters on the 2-nd cathode  
    Float_t     fY[4]  ;          // Y of clusters

  ClassDef(AliMUONReconstHit,1)  // Reconstructed Hit Object for set:MUON
};
#endif
