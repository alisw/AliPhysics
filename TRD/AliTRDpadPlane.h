#ifndef ALITRDPADPLANE_H
#define ALITRDPADPLANE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDpadPlane.h,v */

//////////////////////////////////////////////////
//                                              //
//  TRD pad plane class                         //
//                                              //
//  Contains the information on pad postions,   //
//  pad dimensions, tilting angle, etc.         //
//  It also provides methods to identify the    //
//  current pad number from global coordinates. //
//                                              //
//////////////////////////////////////////////////

#include <TObject.h>

class AliTRDgeometry;

//_____________________________________________________________________________
class AliTRDpadPlane : public TObject {

 public:

  AliTRDpadPlane();
  AliTRDpadPlane(Int_t p, Int_t c);
  AliTRDpadPlane(const AliTRDpadPlane &p);
  virtual           ~AliTRDpadPlane();
  AliTRDpadPlane    &operator=(const AliTRDpadPlane &p);
  virtual void       Copy(TObject &p) const;

  Int_t    GetPadRowNumber(const Double_t z);
  Int_t    GetPadColNumber(const Double_t rphi);

  Double_t GetPadRowOffset(const Int_t row, const Double_t z);
  Double_t GetPadColOffset(const Int_t col, const Double_t rphi);
  
 protected:

  AliTRDgeometry *fGeo;       //! TRD geometry       

  Int_t     fPla;             //  Plane number
  Int_t     fCha;             //  Chamber number

  Double_t  fLength;          //  Length of pad plane in z-direction (row)
  Double_t  fWidth;           //  Width of pad plane in rphi-direction (col)

  Double_t  fLengthRim;       //  Length of the rim in z-direction (row)
  Double_t  fWidthRim;        //  Width of the rim in rphi-direction (col)

  Double_t  fLengthOPad;      //  Length of an outer pad in z-direction (row)
  Double_t  fWidthOPad;       //  Width of an outer pad in rphi-direction (col)

  Double_t  fLengthIPad;      //  Length of an inner pad in z-direction (row)
  Double_t  fWidthIPad;       //  Width of an inner pad in rphi-direction (col)

  Double_t  fRowSpacing;      //  Spacing between the pad rows
  Double_t  fColSpacing;      //  Spacing between the pad columns

  Int_t     fNrows;           //  Number of rows
  Int_t     fNcols;           //  Number of columns

  Double_t  fTiltingAngle;    //  Pad tilting angle  

  Double_t *fPadRow;          //! Pad border positions in row direction
  Double_t *fPadCol;          //! Pad border positions in column direction

  ClassDef(AliTRDpadPlane,1)  //  TRD ROC pad plane

};

#endif
