#ifndef ALICHEB3DCALC_H
#define ALICHEB3DCALC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Author: Ruben Shahoyan
// ruben.shahoyan@cern.ch   09/09/2006
// See Comments in AliCheb3D.h
//
#include <TNamed.h>

class AliCheb3DCalc: public TNamed
{
 public:
  AliCheb3DCalc();
  AliCheb3DCalc(FILE* stream); // read coefs from text file
  AliCheb3DCalc(const AliCheb3DCalc &src);
  AliCheb3DCalc& operator= (const AliCheb3DCalc &rhs);
  ~AliCheb3DCalc() {Clear();}
  //
  void       Print(Option_t* opt="")                                   const;
  void       LoadData(FILE* stream);
  Float_t    Eval(Float_t  *par)                                       const;
  //
#ifdef _INC_CREATION_ALICHEB3D_
  void       SaveData(const char* outfile,Bool_t append=kFALSE)        const;
  void       SaveData(FILE* stream=stdout)                             const;
#endif
  //
  static void ReadLine(TString& str,FILE* stream);
  //
 protected:
  //
  void       Clear(Option_t* option = "");
  void       Init0();
  Float_t    ChebEval1D(Float_t  x, const Float_t * array, int ncf)     const;  
  void       InitRows(int nr);
  void       InitCols(int nc);
  void       InitElemBound2D(int ne);
  void       InitCoefs(int nc);
  Int_t*     GetNColsAtRow()                                            const {return fNColsAtRow;}
  Int_t*     GetColAtRowBg()                                            const {return fColAtRowBg;}
  Int_t*     GetCoefBound2D0()                                          const {return fCoefBound2D0;}
  Int_t*     GetCoefBound2D1()                                          const {return fCoefBound2D1;}
  Float_t *  GetCoefs()                                                 const {return fCoefs;}
  //
 protected:
  Int_t      fNCoefs;            // total number of coeeficients
  Int_t      fNRows;             // number of significant rows in the 3D coeffs matrix
  Int_t      fNCols;             // max number of significant cols in the 3D coeffs matrix
  Int_t      fNElemBound2D;      // number of elements (fNRows*fNCols) to store for the 2D boundary of significant coeffs
  Int_t*     fNColsAtRow;        //[fNRows] number of sighificant columns (2nd dim) at each row of 3D coefs matrix
  Int_t*     fColAtRowBg;        //[fNRows] beginnig of significant columns (2nd dim) for row in the 2D boundary matrix
  Int_t*     fCoefBound2D0;      //[fNElemBound2D] 2D matrix defining the boundary of significance for 3D coeffs.matrix (Ncoefs for col/row)
  Int_t*     fCoefBound2D1;      //[fNElemBound2D] 2D matrix defining the start beginnig of significant coeffs for col/row
  Float_t *  fCoefs;             //[fNCoefs] array of Chebyshev coefficients
  //
  Float_t *  fTmpCf1;            //[fNCols] temp. coeffs for 2d summation
  Float_t *  fTmpCf0;            //[fNRows] temp. coeffs for 1d summation
  //
  ClassDef(AliCheb3DCalc,1)  // Class for interpolation of 3D->1 function by Chebyshev parametrization 
};

inline Float_t AliCheb3DCalc::ChebEval1D(Float_t  x, const Float_t * array, int ncf ) const
{
  // evaluate 1D Chebyshev parameterization. x is the argument mapped to [-1:1] interval
  Float_t b0, b1, b2;
  Float_t x2 = x+x;
  b0 = array[--ncf]; 
  b1 = b2 = 0;
  for (int i=ncf;i--;) {
    b2 = b1;
    b1 = b0;
    b0 = array[i] + x2*b1 -b2;
  }
  return b0 - x*b1;
}

#endif
