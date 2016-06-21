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
  enum {ktgp,kz};
  //
 public:
  AliCheb2DStack();
  virtual ~AliCheb2DStack();
  //
  AliCheb2DStack(int nSlices, int dimOut, const float bmin[2], const float bmax[2],const float* dead=0, const float *rowXI=0);
  //
  Int_t           GetNSlices()        const {return fNSlices;}
  Int_t           GetDimOut()         const {return fDimOut;}
  const Float_t*  GetBoundMin()       const {return fBMin;}
  const Float_t*  GetBoundMax()       const {return fBMax;}
  //
  void            SetXRowInv(const float* xi) {fRowXI = xi;}
  const float*    GetXRowInv() const {return fRowXI;}
  virtual void     Eval(int sliceID, const float *par, float *res) const = 0;
  virtual Float_t  Eval(int sliceID, int dimOut, const float *par) const = 0;
  virtual void     EvalDeriv(int sliceID, int dim, const Float_t  *par, float* res) const = 0;

  Bool_t        IsInside(const float *par) const;
  //
  void          Print(const Option_t* opt="")            const;
  virtual void  PrintSlice(int isl, const Option_t* opt) const = 0;
  //
  static float  ChebEval1D(float x, const float* array, int ncf);
  static float  ChebEval1D(float x, const Short_t* array, int ncf);
  static float  ChebEval1Deriv(float x, const float* array, int ncf);
  static float  ChebEval1Deriv(float x, const Short_t* array, int ncf);
  static void   SetDefPrecision(float prc=1e-4) {fgkDefPrec = prc>1e-8 ? prc:1e-8;}
  static float  GetDefPrecision() {return fgkDefPrec;}
  //
 protected:
  void          MapToInternal(int slice, const float* xy, float *xyint) const;
  void          MapToInternal(int slice, const float* xy, float &x1, float &x2) const;
  float         MapToExternal(int slice, float x,int dim) const;
  float*        DefineGrid(int slice, int dim, const int np[2]) const;
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
  Float_t       fBScaleZ;           // scale for Z boundary mapping to [-1:1] interval
  Float_t       fBOffsetZ;          // offset for Z boundary mapping to [-1:1] interval
  Float_t       fDead[2];           // dead zone in cm to account if X for each row is provided
  const Float_t* fRowXI;            //! optional external!!! set 1/X for each row if dead zones to be accounted
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
  ClassDef(AliCheb2DStack,2)        // stack of 2D->fDimOut Chebyshev parameterization slices
};


//_________________________________________________________
inline void AliCheb2DStack::MapToInternal(int slice, const float* xy, float args[2]) const
{
  // map xy to [-1:1]
  // for Z we don't have dead zones
  args[kz] = (xy[kz]-fBOffsetZ)*fBScaleZ;
  float tmn = fBMin[ktgp], tmx = fBMax[ktgp];
  if (fRowXI) {
    tmn += fDead[0]*fRowXI[slice];
    tmx -= fDead[1]*fRowXI[slice];
  }
  args[ktgp] = 2.f*(xy[ktgp]-tmn)/(tmx-tmn) - 1.f;
#ifdef _BRING_TO_BOUNDARY2D_
  if      (args[kz]<-1.0f) args[kz]=-1.0f;
  else if (args[kz]> 1.0f) args[kz]=1.0f;
  if      (args[ktgp]<-1.0f) args[ktgp]=-1.0f;
  else if (args[ktgp]> 1.0f) args[ktgp]=1.0f;
#endif
}

//_________________________________________________________
inline void AliCheb2DStack::MapToInternal(int slice, const float* xy, float &x0, float &x1) const
{
  // map xy to [-1:1]
  x1 = (xy[kz]-fBOffsetZ)*fBScaleZ;
  float tmn = fBMin[ktgp], tmx = fBMax[ktgp];
  if (fRowXI) {
    tmn += fDead[0]*fRowXI[slice];
    tmx -= fDead[1]*fRowXI[slice];
  }
  x0 = 2.f*(xy[ktgp]-tmn)/(tmx-tmn) - 1.f;
  //
#ifdef _BRING_TO_BOUNDARY2D_
  if (x0 < -1.0f) x0 =-1.0f; else if (x0 > 1.0f) x0 = 1.0f;
  if (x1 < -1.0f) x1 =-1.0f; else if (x1 > 1.0f) x1 = 1.0f;
#endif
}

//__________________________________________________________________________________________
inline float AliCheb2DStack::MapToExternal(int slice, float x,int dim)  const 
{
  // map from [-1:1] to x for dim-th input dimension
  if (dim==kz) return x/fBScaleZ+fBOffsetZ;
  // for y/x need to account for dead zone
  float tmn = fBMin[ktgp], tmx = fBMax[ktgp];
  if (fRowXI) {
    tmn += fDead[0]*fRowXI[slice];
    tmx -= fDead[1]*fRowXI[slice];
  }
  return 0.5*(x+1.0f)*(tmx-tmn)+tmn;
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


//__________________________________________________________________________________________
inline float AliCheb2DStack::ChebEval1Deriv(float x, const float* array, int ncf) 
{
  // evaluate 1D Chebyshev parameterization. x is the argument mapped to [-1:1] interval
  if (--ncf<1) return 0;
  Float_t b0, b1(0), b2(0), x2(x+x);
  float dcf0(0),dcf1,dcf2(0);
  b0 = dcf1 = 2*ncf*array[ncf];
  if (!(--ncf)) return b0/2;

  for (int i=ncf;i--;) {
    b2 = b1;
    b1 = b0;
    dcf0 = dcf2 + 2*(i+1)*array[i+1];
    b0 = dcf0 + x2*b1 -b2;
    dcf2 = dcf1;  
    dcf1 = dcf0;
  }
  return b0 - x*b1 - dcf0/2;
  //
}

//__________________________________________________________________________________________
inline float AliCheb2DStack::ChebEval1Deriv(float x, const Short_t* array, int ncf) 
{
  // evaluate 1D Chebyshev parameterization. x is the argument mapped to [-1:1] interval
  if (--ncf<1) return 0;
  Float_t b0, b1(0), b2(0), x2(x+x);
  float dcf0(0),dcf1,dcf2(0);
  b0 = dcf1 = 2*ncf*array[ncf];
  if (!(--ncf)) return b0/2;

  for (int i=ncf;i--;) {
    b2 = b1;
    b1 = b0;
    dcf0 = dcf2 + 2*(i+1)*array[i+1];
    b0 = dcf0 + x2*b1 -b2;
    dcf2 = dcf1;  
    dcf1 = dcf0;
  }
  return b0 - x*b1 - dcf0/2;
  //
}


#endif
