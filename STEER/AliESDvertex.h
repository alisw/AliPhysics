#ifndef ALIESDVERTEX_H
#define ALIESDVERTEX_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   Alice ESD vertex object                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TArrayD.h"


class AliESDvertex : public TObject
{
public:
  AliESDvertex();
  virtual ~AliESDvertex() {}
  
protected:
  Int_t        fNPrimary;               // Number of primary tracks
  TArrayF      fCoordinates;            // (3) Vertex coordinates
  TArrayF      fErrorMatrix;            // (6) Error Matrix
  TObjArray    fPrimaryTracks;          // List of primary tracks
  Float_t      fEffectiveMass;          // Effective Mass
  Float_t      fEffectiveMassError;     // Effective Mass Error
private:
  AliESDvertex(const AliESDvertex & esdv);
  AliESDvertex & operator=(const AliESDvertex &) {return (*this);}
  
  ClassDef(AliESDvertex,1)  //ESDvertex 
};


#endif 

