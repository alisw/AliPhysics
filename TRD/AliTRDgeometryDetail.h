#ifndef ALITRDGEOMETRYDETAIL_H
#define ALITRDGEOMETRYDETAIL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Detailed TRD geometry for the spaceframe without holes                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDgeometryFull.h"

class AliTRDgeometryDetail : public AliTRDgeometryFull {

 public:

  AliTRDgeometryDetail();
  virtual ~AliTRDgeometryDetail();

          void    CreateGeometry(Int_t *idtmed);
          void    CreateReadout(Int_t *idtmed) const;
          void    CreateCooling(Int_t *idtmed) const;
          void    PositionReadout(Int_t ipla, Int_t icha);
          void    PositionCooling(Int_t ipla, Int_t icha, Int_t idrotm);
          Int_t   IsVersion() const { return 2; };
          void    Init();

 protected:

  ClassDef(AliTRDgeometryDetail,1) // Detailed TRD geometry without hole

};

#endif
