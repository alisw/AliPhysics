#ifndef ALICHEB2DSTACKF_H
#define ALICHEB2DSTACKF_H

/*****************************************************************************
 *         Stack of 2D->ND Chebishev parameterizations with Float_t          *
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

class AliCheb2DStackF : public AliCheb2DStack
{
 public:
  AliCheb2DStackF();
  virtual ~AliCheb2DStackF();
  //
  AliCheb2DStackF(stFun_t fun, int nSlices, int dimOut, const float bmin[2], const float bmax[2], 
		 const int np[2], const float* precD=0);
  AliCheb2DStackF(stFun_t fun, int nSlices, int dimOut, const float bmin[2], const float bmax[2], 
		 const int np[][2], const float* precD=0);
  //
  void          Eval(int sliceID, const float *par, float *res) const;
  Float_t       Eval(int sliceID, int dimOut, const float *par) const;
  void          Print(const Option_t* opt="")            const;
  void          PrintSlice(int isl, const Option_t* opt) const;
  //
 protected:
  //
  void          CreateParams(stFun_t fun, const int *np, const float* prc);
  void          ChebFit(const int np[2], float* wspace, float prec);
  void          FillFunValues(stFun_t fun, int slice, int dim, const float *grid, const int np[2]);
  //
 protected:
  //
  Float_t*      fCoeffs;            //[fNCoefsTot] coeffs container (all slices)
  //
 private:
  AliCheb2DStackF(const AliCheb2DStackF& src);            // dummy
  AliCheb2DStackF& operator=(const AliCheb2DStackF& rhs); // dummy
  //
  ClassDef(AliCheb2DStackF,1)        // stack of 2D->fDimOut Chebyshev parameterization slices
};

#endif
