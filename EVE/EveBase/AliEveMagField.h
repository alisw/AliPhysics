// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveMagField_H
#define AliEveMagField_H

#include "TEveTrackPropagator.h"

class AliMagF;

//______________________________________________________________________________
// Short description of AliEveMagField
//

class AliEveMagField : public TEveMagField
{
public:
  AliEveMagField(AliMagF* mf=0);
  virtual ~AliEveMagField() {}

  using TEveMagField::GetField;
  virtual TEveVector GetField(Float_t x, Float_t y, Float_t z) const;

protected:
  AliMagF *fField; //! Pointer to the magnetic field.

private:
  AliEveMagField(const AliEveMagField&);            // Not implemented
  AliEveMagField& operator=(const AliEveMagField&); // Not implemented

  ClassDef(AliEveMagField, 0); // Short description.
};

#endif
