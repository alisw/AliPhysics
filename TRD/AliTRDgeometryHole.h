#ifndef ALITRDGEOMETRYHOLE_H
#define ALITRDGEOMETRYHOLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD geometry with holes                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDgeometry.h"

class AliTRDgeometryHole : public AliTRDgeometry {

 public:

  AliTRDgeometryHole();
  virtual ~AliTRDgeometryHole();

          void    CreateGeometry(Int_t *idtmed); 
          Int_t   IsVersion() const   { return 0; };
          void    Init();

          void    SetPHOShole()       { };
          void    SetRICHhole()       { };
	  Bool_t   IsHole(Int_t iplan, Int_t icham, Int_t /*isec*/) const;
  virtual void    SetOldGeometry();

          Bool_t  GetPHOShole() const { return kTRUE; };
          Bool_t  GetRICHhole() const { return kTRUE; };

 protected:

  ClassDef(AliTRDgeometryHole,2)              // TRD geometry with hole

};

#endif
