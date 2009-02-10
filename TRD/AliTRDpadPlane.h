#ifndef ALITRDPADPLANE_H
#define ALITRDPADPLANE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD pad plane class                                                   //
//                                                                        //
//  Contains the information on ideal pad positions, pad dimensions,      //
//  tilting angle, etc.                                                   //
//  It also provides methods to identify the current pad number from      //
//  local tracking coordinates.                                           //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

//_____________________________________________________________________________
class AliTRDpadPlane : public TObject {

 public:

  AliTRDpadPlane();
  AliTRDpadPlane(Int_t layer, Int_t stack);
  AliTRDpadPlane(const AliTRDpadPlane &p);
  virtual           ~AliTRDpadPlane();
  AliTRDpadPlane    &operator=(const AliTRDpadPlane &p);
  virtual void       Copy(TObject &p) const;

  void     SetLayer(Int_t l)                   { fLayer          = l; };
  void     SetStack(Int_t s)                   { fStack          = s; };
  void     SetRowSpacing(Double_t s)           { fRowSpacing     = s; };
  void     SetColSpacing(Double_t s)           { fColSpacing     = s; };
  void     SetLengthRim(Double_t l)            { fLengthRim      = l; };
  void     SetWidthRim(Double_t w)             { fWidthRim       = w; };
  void     SetNcols(Int_t n)                   { fNcols          = n;
                                                 if (fPadCol) delete[] fPadCol;
                                                 fPadCol         = new Double_t[fNcols]; };
  void     SetNrows(Int_t n)                   { fNrows          = n;
                                                 if (fPadRow) delete[] fPadRow;
                                                 fPadRow         = new Double_t[fNrows]; };
  void     SetPadCol(Int_t ic, Double_t c)     { if (ic < fNcols) fPadCol[ic] = c;       };
  void     SetPadRow(Int_t ir, Double_t r)     { if (ir < fNrows) fPadRow[ir] = r;       };
  void     SetLength(Double_t l)               { fLength          = l; };
  void     SetWidth(Double_t w)                { fWidth           = w; };
  void     SetLengthOPad(Double_t l)           { fLengthOPad      = l; };
  void     SetWidthOPad(Double_t w)            { fWidthOPad       = w; };
  void     SetLengthIPad(Double_t l)           { fLengthIPad      = l; };
  void     SetWidthIPad(Double_t w)            { fWidthIPad       = w; };
  void     SetPadRowSMOffset(Double_t o)       { fPadRowSMOffset  = o; };
  void     SetAnodeWireOffset(Float_t o)       { fAnodeWireOffset = o; };
  void     SetTiltingAngle(Double_t t);

  Int_t    GetPadRowNumber(Double_t z) const;
  Int_t    GetPadRowNumberROC(Double_t z) const;
  Int_t    GetPadColNumber(Double_t rphi) const;

  Double_t GetTiltOffset(Double_t rowOffset) const 
                                             { return fTiltingTan * (rowOffset - 0.5*fLengthIPad); };

  Double_t GetPadRowOffset(Int_t row, Double_t z) const
                                             { if ((row < 0) || (row >= fNrows))
                                                 return -1.0;
                                                else 
                                                 return fPadRow[row] + fPadRowSMOffset - z;        };
  Double_t GetPadRowOffsetROC(Int_t row, Double_t z) const
                                             { if ((row < 0) || (row >= fNrows))
                                                 return -1.0;
                                               else 
                                                 return fPadRow[row] - z;    };

  Double_t GetPadColOffset(Int_t col, Double_t rphi) const
                                             { if ((col < 0) || (col >= fNcols))
                                                 return -1.0;
                                               else
                                                 return rphi - fPadCol[col]; };

  Double_t GetTiltingAngle() const           { return fTiltingAngle;    };

  Int_t    GetNrows() const                  { return fNrows;           };
  Int_t    GetNcols() const                  { return fNcols;           };

  Double_t GetRow0() const                   { return fPadRow[0] + fPadRowSMOffset;    };
  Double_t GetRow0ROC() const                { return fPadRow[0];       };
  Double_t GetCol0() const                   { return fPadCol[0];       };

  Double_t GetRowEnd() const                 { return fPadRow[fNrows-1] - fLengthOPad + fPadRowSMOffset; };
  Double_t GetRowEndROC() const              { return fPadRow[fNrows-1] - fLengthOPad; };
  Double_t GetColEnd() const                 { return fPadCol[fNcols-1] + fWidthOPad;  };

  Double_t GetRowPos(Int_t row) const        { return fPadRow[row] + fPadRowSMOffset;  };
  Double_t GetRowPosROC(Int_t row) const     { return fPadRow[row];     };
  Double_t GetColPos(Int_t col) const        { return fPadCol[col];     };
  
  Double_t GetRowSize(Int_t row) const       { if ((row == 0) || (row == fNrows-1))
                                                 return fLengthOPad;
                                               else
                                                 return fLengthIPad; };
  Double_t GetColSize(Int_t col) const       { if ((col == 0) || (col == fNcols-1))
                                                 return fWidthOPad;
                                               else
                                                 return fWidthIPad;     };

  Double_t GetLengthRim() const              { return fLengthRim;       };
  Double_t GetWidthRim() const               { return fWidthRim;        };

  Double_t GetRowSpacing() const             { return fRowSpacing;      };
  Double_t GetColSpacing() const             { return fColSpacing;      };

  Double_t GetLengthOPad() const             { return fLengthOPad;      };
  Double_t GetLengthIPad() const             { return fLengthIPad;      };

  Double_t GetWidthOPad() const              { return fWidthOPad;       };
  Double_t GetWidthIPad() const              { return fWidthIPad;       };

  Double_t GetAnodeWireOffset() const        { return fAnodeWireOffset; };

 protected:

  Int_t     fLayer;           //  Layer number
  Int_t     fStack;           //  Stack number

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
  Double_t  fTiltingTan;      //  Tangens of pad tilting angle

  Double_t *fPadRow;          //  Pad border positions in row direction
  Double_t *fPadCol;          //  Pad border positions in column direction

  Double_t  fPadRowSMOffset;  //  To be added to translate local ROC system to local SM system

  Double_t  fAnodeWireOffset; //  Distance of first anode wire from pad edge

  ClassDef(AliTRDpadPlane,6)  //  TRD ROC pad plane

};

#endif
