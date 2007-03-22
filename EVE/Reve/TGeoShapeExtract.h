// $Header$

// Copyright (C) 1999-2005, Matevz Tadel. All rights reserved.
// This file is part of GLED, released under GNU General Public License version 2.
// For the licensing terms see $GLEDSYS/LICENSE or http://www.gnu.org/.

#ifndef RootGeo_TGeoShapeExtract_H
#define RootGeo_TGeoShapeExtract_H

#include <TNamed.h>

class TList;
class TGeoShape;

class TGeoShapeExtract : public TNamed
{
  friend class ZGeoRepacker;

  TGeoShapeExtract(const TGeoShapeExtract&);            // Not implemented
  TGeoShapeExtract& operator=(const TGeoShapeExtract&); // Not implemented

protected:
  Double_t    mTrans[16];
  Float_t     mRGBA[4];
  Bool_t      mRnrSelf;
  Bool_t      mRnrElements;
  TGeoShape*  mShape;
  TList*      mElements;

public:
  TGeoShapeExtract(const Text_t* n="TGeoShapeExtract", const Text_t* t=0);
  ~TGeoShapeExtract();

  Bool_t HasElements();

  void SetTrans(const Double_t arr[16]);
  void SetRGBA (const Float_t  arr[4]);

  Double_t*  GetTrans()       { return mTrans; }
  Float_t*   GetRGBA()        { return mRGBA;  }
  Bool_t     GetRnrSelf()     { return mRnrSelf;     }
  Bool_t     GetRnrElements() { return mRnrElements; }
  TGeoShape* GetShape()       { return mShape;    }
  TList*     GetElements()    { return mElements; }

  ClassDef(TGeoShapeExtract, 1)
}; // endclass TGeoShapeExtract

#endif
