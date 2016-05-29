/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliCheb2DStackS.h"
#include "AliLog.h"
#include <TMath.h>
#include <limits.h>

ClassImp(AliCheb2DStackS)

//____________________________________________________________________
AliCheb2DStackS::AliCheb2DStackS() 
:  AliCheb2DStack()
  ,fParScale(0)
  ,fParHVar(0)
  ,fCoeffs(0)
{
  // Default constructor
}

//____________________________________________________________________
AliCheb2DStackS::AliCheb2DStackS(stFun_t fun, int nSlices, int dimOut, 
				 const float bmin[2],const float bmax[2], 
				 const int np[2], const float* dead, const float *rowXI,
				 const float* precD)
: AliCheb2DStack(nSlices,dimOut,bmin,bmax,dead,rowXI)
  ,fParScale(new Float_t[nSlices*dimOut])
  ,fParHVar(new Float_t[nSlices*dimOut])
  ,fCoeffs(0)
{
  // create stack of 2D->dimOut Chebyshev parameterizations debined in 2 dimensions between bmin and bmax,
  // and trained with function fun on 2D grid on np points. 
  // Truncate each precision of each output dimension parameterization i to precD[i] if precD!=0, or prec
  //
  int *npd = new int[2*dimOut];
  int maxcoefs=np[0]*np[1]*fNSlices, maxrows=np[0]*fNSlices;
  for (int i=0;i<dimOut;i++) {
    npd[2*i]   = np[0];
    npd[2*i+1] = np[1];
  }
  CheckDimensions(npd);  // basic check
  fNCols  = new UChar_t[maxrows];
  fCoeffs = new Short_t[maxcoefs];
  CreateParams(fun, npd, precD);
  delete[] npd;
  //
  Print();
  //
}

//____________________________________________________________________
AliCheb2DStackS::AliCheb2DStackS(stFun_t fun, int nSlices, int dimOut, 
				 const float bmin[2],const float bmax[2],
				 const int np[][2], const float* dead, const float *rowXI,
				 const float* precD)
:AliCheb2DStack(nSlices,dimOut,bmin,bmax,dead,rowXI)
  ,fParScale(new Float_t[nSlices*dimOut])
  ,fParHVar(new Float_t[nSlices*dimOut])
  ,fCoeffs(0)
{
  // create stack of 2D->dimOut Chebyshev parameterizations debined in 2 dimensions between bmin and bmax,
  // and trained with function fun on 2D grid on np[i][2] points for i-th output dimension. 
  // Truncate each precision of each output dimension parameterization i to precD[i] if precD!=0, or prec   
  //
  int maxcoefs=0, maxrows=0;
  for (int id=fDimOut;id--;) {
    maxcoefs += np[id][0]*np[id][1];
    maxrows  += np[id][0];
  }
  fNCols  = new UChar_t[maxrows*fNSlices];
  fCoeffs = new Short_t[maxcoefs*fNSlices];
  CreateParams(fun, (const int*)np, precD);
  Print();
  //
}

//____________________________________________________________________
AliCheb2DStackS::~AliCheb2DStackS()
{
  // D-tor
  delete[] fCoeffs;
  delete[] fParScale;
  delete[] fParHVar;
}

//____________________________________________________________________
void AliCheb2DStackS::Eval(int sliceID, const float  *par, float *res) const
{
  // evaluate Chebyshev parameterization for 2d->DimOut function at sliceID
  float p0,p1;
  MapToInternal(sliceID,par,p0,p1);
  int pid = sliceID*fDimOut;
  const UChar_t *rows = &fNRows[pid];                      // array of fDimOut rows for current slice
  const UChar_t *cols = &fNCols[fColEntry[sliceID]];       // array of columns per row for current slice
  const Short_t *cfs  = &fCoeffs[fCoeffsEntry[sliceID]];   // array of coefficients for current slice
  const Float_t *scl  = &fParScale[pid];
  const Float_t *hvr  = &fParHVar[pid];
  for (int id=0;id<fDimOut;id++) {
    int nr = *rows++;                            // N rows in the matrix of coeffs for given dimension 
    for (int ir=0;ir<nr;ir++) {
      int nc = *cols++;                          // N of significant colums at this row
      fWSpace[ir] = ChebEval1D(p1,cfs,nc); // interpolation of Cheb. coefs along row
      cfs += nc;                                 // prepare coefs for the next row
    }
    res[id] = ChebEval1D(p0,fWSpace,nr) * (*scl++) + (*hvr++);
  }
  //
}

//____________________________________________________________________
Float_t AliCheb2DStackS::Eval(int sliceID, int dimOut, const float  *par) const
{
  // evaluate Chebyshev parameterization for requested output dimension only at requested sliceID
  float p0,p1;
  MapToInternal(sliceID,par,p0,p1);
  int pid = sliceID*fDimOut;
  const UChar_t *rows = &fNRows[pid];                      // array of fDimOut rows for current slice
  const UChar_t *cols = &fNCols[fColEntry[sliceID]];       // array of columns per row for current slice
  const Short_t *cfs  = &fCoeffs[fCoeffsEntry[sliceID]];   // array of coefficients for current slice
  const Float_t *scl  = &fParScale[pid];
  const Float_t *hvr  = &fParHVar[pid];
  while (dimOut) {
    for (int ir=*rows++;ir--;) cfs += *cols++;  // go to the matrix of needed row
    dimOut--;
    scl++;
    hvr++;
  }
  int nr = *rows++;                             // N rows in the matrix of coeffs for given dimension 
  for (int ir=0;ir<nr;ir++) {
    int nc = *cols++;                          // N of significant colums at this row
    fWSpace[ir] = ChebEval1D(p1,cfs,nc); // interpolation of Cheb. coefs along row
    cfs += nc;                                 // prepare coefs for the next row
  }
  return ChebEval1D(p0,fWSpace,nr) * (*scl) + (*hvr);
   //
}

//____________________________________________________________________
void AliCheb2DStackS::CreateParams(stFun_t fun, const int *np, const float* prc)
{
  // create parameterizations
  //
  // temporary space for max possible coeffs, rows etc
  float **grids = new float*[fDimOut]; // Chebyshev grids for each output dimension
  int maxSpace = 1, totSpace = 1;
  Bool_t sameGrid = kTRUE;
  int ref0=np[0],ref1=np[1];
  for (int id=fDimOut;id--;) {
    int nsp = np[2*id]*np[2*id+1];
    if (maxSpace<nsp) maxSpace = nsp;
    totSpace += nsp;
    if (ref0!=np[2*id] || ref1!=np[2*id+1]) sameGrid = kFALSE;
  }
  //
  // save pointers to recover the beggining of arrays later
  UChar_t* nRows0 = fNRows;
  UChar_t* nCols0 = fNCols;
  Short_t*  coeffs0 = fCoeffs;
  //
  float *tmpCoef2D = new float[maxSpace]; // temporary workspace for coeffs
  float *tmpVals = new float[totSpace]; // temporary workspace for function values
  //
  for (int isl=0;isl<fNSlices;isl++) {
    for (int id=fDimOut;id--;) grids[id] = DefineGrid(isl, id, &np[2*id]);
    fCoeffsEntry[isl] = fCoeffs - coeffs0;   // offset of the given slice coeffs
    fColEntry[isl]    = fNCols  - nCols0;    // offset of the given slice columns dimensions
    //
    if (sameGrid) FillFunValues(fun, isl, grids[0], &np[0], tmpVals);
    else {
      int slot = 0;
      for (int id=0;id<fDimOut;id++) {
	FillFunValues(fun, isl, id, grids[id], &np[2*id], tmpVals+slot); // get values for single dimensions
	slot += np[2*id]*np[2*id+1];
      }
    }
    int slot = 0;
    for (int id=0;id<fDimOut;id++) {
      float prcs = prc ? prc[id] : fgkDefPrec;
      fParScale[isl*fDimOut+id] = ChebFit(&np[2*id], tmpVals+slot, tmpCoef2D, prcs);
      slot += np[2*id]*np[2*id+1];
    }
    for (int id=fDimOut;id--;) delete[] grids[id];
  }
  delete[] grids;
  delete[] tmpCoef2D;
  delete[] tmpVals;
  //
  fNCoefsTot = fCoeffs-coeffs0; // size of final coeffs array
  //
  fNRows = nRows0;
  // redefine arrays in compressed form, clean temp. stuff
  fNCols = new UChar_t[fNRowsTot];
  memcpy(fNCols, nCols0, fNRowsTot*sizeof(UChar_t));
  delete[] nCols0;
  fCoeffs = new Short_t[fNCoefsTot];
  memcpy(fCoeffs, coeffs0, fNCoefsTot*sizeof(Short_t));
  delete[] coeffs0;
  //
}

//____________________________________________________________________
Float_t AliCheb2DStackS::ChebFit(const int np[2], const float* tmpVals, float* tmpCoef2D, float prec)
{
  // prepare Cheb.fit for 2D -> single output dimension
  //
  int maxDim = TMath::Max(np[0],np[1]);
  memset(tmpCoef2D,0,np[0]*np[1]*sizeof(float));
  //
  for (int id1=np[1];id1--;) { // create Cheb.param for each node of 2nd input dimension
    int nc = CalcChebCoefs(&tmpVals[id1*np[0]], np[0], fWSpace, -1);
    for (int id0=nc;id0--;) tmpCoef2D[id1 + id0*np[1]] = fWSpace[id0];
  }
  //
  // once each 1d slice of given 2d slice is parametrized, parametrize the Cheb.coeffs
  float mxAbs = -1;  // we need largest coeff for rescaling to +-MaxShort range
  for (int id0=np[0];id0--;) {
    CalcChebCoefs(&tmpCoef2D[id0*np[1]], np[1], fWSpace, -1);
    for (int id1=np[1];id1--;) {
      tmpCoef2D[id1+id0*np[1]] = fWSpace[id1]; 
      if (TMath::Abs(fWSpace[id1])>mxAbs) mxAbs = TMath::Abs(fWSpace[id1]);
    }
  }  
  //
  float scl = mxAbs>0 ? SHRT_MAX/mxAbs : 1;
  // rescale/truncate to +-MaxShort range
  for (int id0=np[0];id0--;) {
    for (int id1=np[1];id1--;) {
      int id = id1 + id0*np[1];
      Short_t cfs = Short_t(tmpCoef2D[id1+id0*np[1]]*scl);
      tmpCoef2D[id] = float(cfs);
    }
  }
  // 
  float rTiny = 0.5*prec*scl/maxDim; // neglect coefficient below this threshold
  //
  // now find 1D curve which separates significant coefficients of 2D matrix 
  // from nonsignificant ones (up to prec)
  //  double resid = 0;
  int ncNZero=0, nRows = np[0];  // find max significant row
  for (int id0=np[0];id0--;) {
    fNCols[id0]=0;  
    double resid = 0;
    for (int id1=np[1];id1--;) {
      int id = id1 + id0*np[1];     
      float cfa = TMath::Abs(tmpCoef2D[id]);
      if (cfa < rTiny) {tmpCoef2D[id] = 0; continue;} // neglect coefs below the threshold
      resid += cfa;
      //      if (resid<prec) continue; // this coeff is negligible
      if (resid<rTiny) continue; // this coeff is negligible
      //      resid -= cfa;             // otherwise go back 1 step
      fNCols[id0] = id1+1;     // how many coefs to keep
      break;
    }
    if (fNCols[id0]) ncNZero++;
    else if (!ncNZero) nRows--;  // decrease N significant rows untile 1st non-0 column met
  }
  //
  // find max significant column and fill the coeffs in the storage 
  for (int id0=0;id0<nRows;id0++) {
    int nc = *fNCols++;
    for (int id1=0;id1<nc;id1++) {
      *fCoeffs++ = Short_t(tmpCoef2D[id1+id0*np[1]]);
    }
  }
  *fNRows++ = nRows;
  fNRowsTot += nRows;
  //
  return 1./scl;
}

//____________________________________________________________________
void AliCheb2DStackS::FillFunValues(stFun_t fun, int slice, int dim, const float *grid, 
				    const int np[2], float* vals)
{
  // fill function values on the grid
  float args[2];  
  float minv=1e15,maxv=-1e15;
  for (int id1=np[1];id1--;) {
    args[1] = grid[id1];
    for (int id0=np[0];id0--;) { // 
      args[0] = grid[np[1]+id0];
      fun(slice, args, fWSpace);
      vals[id1*np[0] + id0] = fWSpace[dim];
      if (minv>fWSpace[dim]) minv = fWSpace[dim];
      if (maxv<fWSpace[dim]) maxv = fWSpace[dim];
    }
  }
  fParHVar[slice*fDimOut+dim] = (maxv+minv)/2.;
  //
  // map values to +- fWSpace
  float offs = fParHVar[slice*fDimOut+dim];
  for (int id1=np[1];id1--;) for (int id0=np[0];id0--;) vals[id1*np[0] + id0] -= offs;
  //
}

//____________________________________________________________________
void AliCheb2DStackS::FillFunValues(stFun_t fun, int slice, const float *grid, 
				    const int np[2], float* vals)
{
  // fill function values on the grid for all output dimension at same time (same grid)
  float args[2];  
  float *minv = new float[fDimOut], *maxv = new float[fDimOut];
  for (int dim=fDimOut;dim--;) {
    minv[dim] = 1e15;
    maxv[dim] =-1e15;
  }
  int stepDim = np[0]*np[1];
  for (int id1=np[1];id1--;) {
    args[1] = grid[id1];
    for (int id0=np[0];id0--;) { // 
      args[0] = grid[np[1]+id0];
      fun(slice, args, fWSpace);
      for (int dim=fDimOut;dim--;) {
	vals[stepDim*dim + id1*np[0] + id0] = fWSpace[dim];
	if (minv[dim]>fWSpace[dim]) minv[dim] = fWSpace[dim];
	if (maxv[dim]<fWSpace[dim]) maxv[dim] = fWSpace[dim];
      }
    }
  }
  for (int dim=fDimOut;dim--;) {
    float offs = fParHVar[slice*fDimOut+dim] = (maxv[dim]+minv[dim])/2.;
    // map values to +- fWSpace
    float* valsdim = vals+stepDim*dim;
    for (int id1=np[1];id1--;) for (int id0=np[0];id0--;) valsdim[id1*np[0] + id0] -= offs;
  }
  //
  delete[] minv;
  delete[] maxv;
  //
}

//__________________________________________________________________________________________
void AliCheb2DStackS::Print(const Option_t* opt) const
{
  printf("S*"); 
  AliCheb2DStack::Print(opt);
}

//__________________________________________________________________________________________
void AliCheb2DStackS::PrintSlice(int isl, const Option_t* opt) const
{
  // print slice info
  //
  TString opts = opt; opts.ToLower();
  Bool_t showcf = opts.Contains("c");
  int par0 = isl*fDimOut;
  const UChar_t* rows = &fNRows[par0];
  const UChar_t* cols = &fNCols[fColEntry[isl]];
  const Short_t* cfs  = &fCoeffs[fCoeffsEntry[isl]];
  printf("#%3d ",isl);  
  if (showcf) printf("\n");
  for (int id=0;id<fDimOut;id++) {
    int nr = *rows++;
    int ncmax=0,ncf=0;
    for (int ir=0;ir<nr;ir++) {
      ncf += cols[ir];
      if (cols[ir]>ncmax) ncmax = cols[ir];
    }
    printf("D%d: %4d coefs in %3dx%3d (Scl:%.2e Hv:%.2e)| ",
	   id,ncf,nr,ncmax,fParScale[isl*fDimOut+id],fParHVar[isl*fDimOut+id]);
    if (showcf) {
      printf("\n");
      for (int ir=0;ir<nr;ir++) {
	for (int ic=0;ic<cols[ir];ic++) printf("%+6d ",*cfs++); printf("\n");
      }
    }
    cols += nr; // cols entry for next dimension
    //
  }
  if (!showcf) printf("\n");
  //
}
