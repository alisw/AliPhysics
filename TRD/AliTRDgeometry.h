#ifndef ALITRDGEOMETRY_H
#define ALITRDGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
#include <TMath.h>

#include "AliRun.h"
#include "AliRecPoint.h"
#include "AliGeometry.h"

#include "AliTRDconst.h"

class AliTRDgeometry : public AliGeometry {

 public:

  AliTRDgeometry();
  virtual ~AliTRDgeometry();

  virtual void    CreateGeometry(Int_t *idtmed);
  virtual Int_t   IsVersion() const = 0;
  virtual void    Init();
  virtual Bool_t  Local2Global(Int_t d, Float_t *local, Float_t *global) const;
  virtual Bool_t  Local2Global(Int_t p, Int_t c, Int_t s, Float_t *local, Float_t *global) const;
  virtual Bool_t  Rotate(Int_t d, Float_t *pos, Float_t *rot);
  virtual Bool_t  RotateBack(Int_t d, Float_t *rot, Float_t *pos) const;

  virtual void    SetPHOShole() = 0;
  virtual void    SetRICHhole() = 0;

  virtual void    SetRowPadSize(Float_t size)          { fRowPadSize  = size; };
  virtual void    SetColPadSize(Float_t size)          { fColPadSize  = size; };
  virtual void    SetTimeBinSize(Float_t size)         { fTimeBinSize = size; };

  virtual Bool_t  GetPHOShole() = 0;
  virtual Bool_t  GetRICHhole() = 0;

  virtual Int_t   GetDetector(Int_t p, Int_t c, Int_t s) const;
  virtual Int_t   GetPlane(Int_t d) const;
  virtual Int_t   GetChamber(Int_t d) const;
  virtual Int_t   GetSector(Int_t d) const;

  virtual Int_t   GetRowMax(Int_t p, Int_t c, Int_t s) { return fRowMax[p][c][s]; };
  virtual Int_t   GetColMax(Int_t p)                   { return fColMax[p];       };
  virtual Int_t   GetTimeMax()                         { return fTimeMax;         };
 
  virtual Float_t GetRow0(Int_t p, Int_t c, Int_t s)   const { return fRow0[p][c][s]; };
  virtual Float_t GetCol0(Int_t p)                     const { return fCol0[p];       };
  virtual Float_t GetTime0(Int_t p)                    const { return fTime0[p];      };

  virtual Float_t GetRowPadSize()                      const { return fRowPadSize;  };
  virtual Float_t GetColPadSize()                      const { return fColPadSize;  };
  virtual Float_t GetTimeBinSize()                     const { return fTimeBinSize; };

  virtual void    GetGlobal(const AliRecPoint *p, TVector3 &pos, TMatrix &mat) const; 
  virtual void    GetGlobal(const AliRecPoint *p, TVector3 &pos) const;   

 protected:

  Float_t         fCwidth[kNplan];                 // Width of the chambers

  Int_t           fRowMax[kNplan][kNcham][kNsect]; // Number of pad-rows
  Int_t           fColMax[kNplan];                 // Number of pad-columns
  Int_t           fTimeMax;                        // Number of time buckets

  Float_t         fRow0[kNplan][kNcham][kNsect];   // Row-position of pad 0
  Float_t         fCol0[kNplan];                   // Column-position of pad 0
  Float_t         fTime0[kNplan];                  // Time-position of pad 0

  Float_t         fRowPadSize;                     // Pad size in z-direction
  Float_t         fColPadSize;                     // Pad size in rphi-direction
  Float_t         fTimeBinSize;                    // Size of the time buckets

  ClassDef(AliTRDgeometry,1)                       // TRD geometry base class

};

#endif
