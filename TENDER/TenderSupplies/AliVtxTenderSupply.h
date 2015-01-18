#ifndef ALIVTXTENDERSUPPLY_H
#define ALIVTXTENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Vertex tender, redo primary vertex on the fly                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <AliTenderSupply.h>

class AliESDVertex;

class AliVtxTenderSupply: public AliTenderSupply {
  
public:
  AliVtxTenderSupply();
  AliVtxTenderSupply(const char *name, const AliTender *tender=NULL);
  
  virtual ~AliVtxTenderSupply(){;}
  
  virtual void              Init(){;}
  virtual void              ProcessEvent();
  //
  Int_t   GetRefitAlgo()              const {return fRefitAlgo;}
  void    SetRefitAlgo(Int_t alg=-1)        {fRefitAlgo = alg;}
  //
private:
  
  AliVtxTenderSupply(const AliVtxTenderSupply&c);
  AliVtxTenderSupply& operator= (const AliVtxTenderSupply&c);

  AliESDVertex *fDiamond;           //!Information about mean vertex  
  Int_t         fRefitAlgo;         //! optional request for vertex refit 

  ClassDef(AliVtxTenderSupply, 1);  // Primary vertex tender task
};


#endif

