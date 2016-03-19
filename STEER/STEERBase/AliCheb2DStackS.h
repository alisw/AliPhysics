#ifndef ALICHEB2DSTACKS_H
#define ALICHEB2DSTACKS_H

/*****************************************************************************
 *         Stack of 2D->ND Chebishev parameterizations with Short_t          *
 *         coefficients representation                                       *
 *                                                                           *
 *  Creation requires user function with signature                           *
 *  void (*fun)(int slice, float* inp2D,float* valND)                        *
 *  and the boundaries of 2D rectangle.                                      *
 *  The precision of interpolation for each dimension is provided in the     *
 *  precD array (max deviation condition). Note, that actual precision       *
 *  will be affected also by the number of nodes evaluated according to      *
 *  requested int np[2] (same number of points for all output dimensions)    *
 *  or int np[dimOut][2] for individual partition for each dimension         *
 *                                                                           *
 *         Author: ruben.shahoyan@cern.ch                                    *
 *****************************************************************************/

#include "AliCheb2DStack.h"

class AliCheb2DStackS : public AliCheb2DStack
{
 public:
  AliCheb2DStackS();
  virtual ~AliCheb2DStackS();
  //
  AliCheb2DStackS(stFun_t fun, int nSlices, int dimOut, const float bmin[2], const float bmax[2], 
		  const int np[2], const float* dead=0, const float *rowXI=0,const float* precD=0);
  AliCheb2DStackS(stFun_t fun, int nSlices, int dimOut, const float bmin[2], const float bmax[2], 
		 const int np[][2], const float* dead=0, const float *rowXI=0, const float* precD=0);
  //
  void          Eval(int sliceID, const float *par, float *res) const;
  Float_t       Eval(int sliceID, int dimOut, const float *par) const;
  void          Print(const Option_t* opt="")            const;
  void          PrintSlice(int isl, const Option_t* opt) const;
  //
 protected:
  //
  void          CreateParams(stFun_t fun, const int *np, const float* prc);
  float         ChebFit(const int np[2], const float* wVals, float* wspace, float prec);
  void          FillFunValues(stFun_t fun, int slice, int dim, const float *grid, const int np[2], float* wVals);
  //
 protected:
  //
  Float_t*      fParScale;          //[fNParams] scaling param. to bring symmetrized variation to +-MaxShort
  Float_t*      fParHVar;           //[fNParams] half of variation within the row
  Short_t*      fCoeffs;            //[fNCoefsTot] coeffs container (all slices)
  //
 private:
  AliCheb2DStackS(const AliCheb2DStackS& src);            // dummy
  AliCheb2DStackS& operator=(const AliCheb2DStackS& rhs); // dummy
  //
  ClassDef(AliCheb2DStackS,2)        // stack of 2D->fDimOut Chebyshev parameterization slices
};


#endif
