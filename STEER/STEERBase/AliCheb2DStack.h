#ifndef ALICHEB2DSTACK_H
#define ALICHEB2DSTACK_H

#include <TObject.h>

/*****************************************************************************
 *          Base class for stack of 2D->ND Chebishev parameterizations       *
 *                                                                           *
 *          Author: ruben.shahoyan@cern.ch                                   *
 *****************************************************************************/

// when _BRING_TO_BOUNDARY2D_ is defined, the point outside of the fitted folume is assumed
// to be on the surface 
#define _BRING_TO_BOUNDARY2D_
//

typedef void (*stFun_t)(int,float*,float*); // original distribution provided via such function

class AliCheb2DStack : public TObject
{
 public:
  enum {kMaxPoints=255};        // max number of points in input dimension
  //
 public:
  AliCheb2DStack();
  virtual ~AliCheb2DStack();
  //
  AliCheb2DStack(int nSlices, int dimOut, const float bmin[2], const float bmax[2]);
  //
  Int_t           GetNSlices()        const {return fNSlices;}
  Int_t           GetDimOut()         const {return fDimOut;}
  const Float_t*  GetBoundMin()       const {return fBMin;}
  const Float_t*  GetBoundMax()       const {return fBMax;}
  //
  virtual void     Eval(int sliceID, const float *par, float *res) const = 0;
  virtual Float_t  Eval(int sliceID, int dimOut, const float *par) const = 0;
  Bool_t        IsInside(const float *par) const;
  //
  void          Print(const Option_t* opt="")            const;
  virtual void  PrintSlice(int isl, const Option_t* opt) const = 0;
  //
  static float  ChebEval1D(float x, const float* array, int ncf);
  static float  ChebEval1D(float x, const Short_t* array, int ncf);
  static void   SetDefPrecision(float prc=1e-4) {fgkDefPrec = prc>1e-8 ? prc:1e-8;}
  static float  GetDefPrecision() {return fgkDefPrec;}
  //
 protected:
  void          SetSliceDim(Int_t slice, Int_t dim);
  void          MapToInternal(const float* xy, float *xyint) const;
  void          MapToInternal(const float* xy, float &x1, float &x2) const;
  float         MapToExternal(float x,int dim) const;
  float*        DefineGrid(int dim, const int np[2]) const;
  void          CheckDimensions(const int *np) const;
  Int_t         CalcChebCoefs(const float *funval,int np, float *outCoefs, float prec);
  //
 protected:
  //
  Int_t         fDimOut;            // each slice maps 2D to fDimOut values
  Int_t         fNSlices;           // number of slices in the stack
  Int_t         fNParams;           // number of parameterizations = fNSlices*fDimOut
  Int_t         fNCoefsTot;         // dimension of coeffs array for all slices
  Int_t         fNRowsTot;          // total number of 1D Cheb. param rows
  Float_t       fBMin[2];           // min boundaries in each dimension
  Float_t       fBMax[2];           // max boundaries in each dimension  
  Float_t       fBScale[2];         // scale for boundary mapping to [-1:1] interval
  Float_t       fBOffset[2];        // offset for boundary mapping to [-1:1] interval
  //
  UChar_t*      fNRows;             //[fNParams] N of used rows in the 2D coeffs matrix of each param
  UChar_t*      fNCols;             //[fNRowsTot] N of used columns in each row
  Int_t*        fCoeffsEntry;       //[fNSlices] start of the coeffs array in fCoeffs for each slice
  Int_t*        fColEntry;          //[fNSlices] start of the Ncolumns array in fNCols for each slice
  //
  static Float_t fgkDefPrec;           // default precision
  static Float_t fWSpace[kMaxPoints];  // workspace

  //
 private:
  AliCheb2DStack(const AliCheb2DStack& src);            // dummy
  AliCheb2DStack& operator=(const AliCheb2DStack& rhs); // dummy
  //
  ClassDef(AliCheb2DStack,1)        // stack of 2D->fDimOut Chebyshev parameterization slices
};


//_________________________________________________________
inline void AliCheb2DStack::MapToInternal(const float* xy, float args[2]) const
{
  // map xy to [-1:1]
  for (int i=2;i--;) {
    args[i] = (xy[i]-fBOffset[i])*fBScale[i];
#ifdef _BRING_TO_BOUNDARY2D_
    if      (args[i]<-1) args[i]=-1;
    else if (args[i]> 1) args[i]=1;
#endif
  }
}

//_________________________________________________________
inline void AliCheb2DStack::MapToInternal(const float* xy, float &x0, float &x1) const
{
  // map xy to [-1:1]
  x0 = (xy[0]-fBOffset[0])*fBScale[0];
  x1 = (xy[1]-fBOffset[1])*fBScale[1];
#ifdef _BRING_TO_BOUNDARY2D_
  if (x0 < -1) x0 =-1; else if (x0 > 1) x0 = 1;
  if (x1 < -1) x1 =-1; else if (x1 > 1) x1 = 1;
#endif
}

//__________________________________________________________________________________________
inline float AliCheb2DStack::MapToExternal(float x,int dim)  const 
{
  // map from [-1:1] to x for dim-th input dimension
  return x/fBScale[dim]+fBOffset[dim];
}

//__________________________________________________________________________________________
inline Bool_t AliCheb2DStack::IsInside(const float *par) const 
{
  // check if the point is inside of the fitted box
  for (int i=2;i--;) if (fBMin[i]>par[i] || par[i]>fBMax[i]) return kFALSE;
  return kTRUE;
}

//__________________________________________________________________________________________
inline float AliCheb2DStack::ChebEval1D(float x, const float* array, int ncf) 
{
  // evaluate 1D Chebyshev parameterization. x is the argument mapped to [-1:1] interval
  if (!ncf) return 0;
  float b0(array[--ncf]), b1(0), b2(0), x2(x+x);
  for (int i=ncf;i--;) {
    b2 = b1;
    b1 = b0;
    b0 = array[i] + x2*b1 -b2;
  }
  return b0 - x*b1;
  //
}

//__________________________________________________________________________________________
inline float AliCheb2DStack::ChebEval1D(float x, const Short_t* array, int ncf) 
{
  // evaluate 1D Chebyshev parameterization. x is the argument mapped to [-1:1] interval
  if (!ncf) return 0;
  float b0(array[--ncf]), b1(0), b2(0), x2(x+x);
  for (int i=ncf;i--;) {
    b2 = b1;
    b1 = b0;
    b0 = array[i] + x2*b1 -b2;
  }
  return b0 - x*b1;
  //
}



#endif
