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


#include "AliCheb2DStackF.h"
#include "AliLog.h"
#include <TMath.h>

ClassImp(AliCheb2DStackF)

//____________________________________________________________________
AliCheb2DStackF::AliCheb2DStackF() 
  :AliCheb2DStack()
  ,fCoeffs(0)
{
  // Default constructor
}

//____________________________________________________________________
AliCheb2DStackF::AliCheb2DStackF(stFun_t fun, int nSlices, int dimOut, const float bmin[2],const float bmax[2], 
			       const int np[2], const float* dead, const float *rowXI, const float* precD)
  :AliCheb2DStack(nSlices,dimOut,bmin,bmax,dead,rowXI)
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
  fCoeffs = new Float_t[maxcoefs];
  CreateParams(fun, npd, precD);
  delete[] npd;
  //
  Print();
  //
}

//____________________________________________________________________
AliCheb2DStackF::AliCheb2DStackF(stFun_t fun, int nSlices, int dimOut, 
				 const float bmin[2],const float bmax[2], 
				 const int np[][2], const float* dead, const float *rowXI, const float* precD)
:AliCheb2DStack(nSlices,dimOut,bmin,bmax,dead,rowXI)
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
  fCoeffs = new Float_t[maxcoefs*fNSlices];
  CreateParams(fun, (const int*)np, precD);
  //
  Print();
}

//____________________________________________________________________
AliCheb2DStackF::~AliCheb2DStackF()
{
  // D-tor
  delete[] fCoeffs;
}

//____________________________________________________________________
void AliCheb2DStackF::Eval(int sliceID, const float  *par, float *res) const
{
  // evaluate Chebyshev parameterization for 2d->DimOut function at sliceID
  float p0,p1;
  MapToInternal(sliceID, par,p0,p1);
  const UChar_t *rows = &fNRows[sliceID*fDimOut];          // array of fDimOut rows for current slice
  const UChar_t *cols = &fNCols[fColEntry[sliceID]];       // array of columns per row for current slice
  const Float_t *cfs  = &fCoeffs[fCoeffsEntry[sliceID]];   // array of coefficients for current slice
  for (int id=0;id<fDimOut;id++) {
    int nr = *rows++;                            // N rows in the matrix of coeffs for given dimension 
    for (int ir=0;ir<nr;ir++) {
      int nc = *cols++;                          // N of significant colums at this row
      fWSpace[ir] = ChebEval1D(p1,cfs,nc); // interpolation of Cheb. coefs along row
      cfs += nc;                                 // prepare coefs for the next row
    }
    res[id] = ChebEval1D(p0,fWSpace,nr);
  }
  //
}

//____________________________________________________________________
Float_t AliCheb2DStackF::Eval(int sliceID, int dimOut, const float *par) const
{
  // evaluate Chebyshev parameterization for requested output dimension only at requested sliceID
  float p0,p1;
  MapToInternal(sliceID,par,p0,p1);
  int pid = sliceID*fDimOut;
  const UChar_t *rows = &fNRows[pid];                      // array of fDimOut rows for current slice
  const UChar_t *cols = &fNCols[fColEntry[sliceID]];       // array of columns per row for current slice
  const Float_t *cfs  = &fCoeffs[fCoeffsEntry[sliceID]];   // array of coefficients for current slice
  while (dimOut) {
    for (int ir=*rows++;ir--;) cfs += *cols++;  // go to the matrix of needed row
    dimOut--;
  }
  int nr = *rows++;                             // N rows in the matrix of coeffs for given dimension 
  for (int ir=0;ir<nr;ir++) {
    int nc = *cols++;                          // N of significant colums at this row
    fWSpace[ir] = ChebEval1D(p1,cfs,nc); // interpolation of Cheb. coefs along row
    cfs += nc;                                 // prepare coefs for the next row
  }
  return ChebEval1D(p0,fWSpace,nr);
  //
}

//____________________________________________________________________
void AliCheb2DStackF::CreateParams(stFun_t fun, const int *np, const float* prc)
{
  // create parameterizations
  //
  // temporary space for max possible coeffs, rows etc
  float **grids = new float*[fDimOut]; // Chebyshev grids for each output dimension
  int maxSpace = 1;
  for (int id=fDimOut;id--;) {
    int nsp = np[2*id]*np[2*id+1];
    if (maxSpace<nsp) maxSpace = nsp;
  }
  //
  // save pointers to recover the beggining of arrays later
  UChar_t* nRows0 = fNRows;
  UChar_t* nCols0 = fNCols;
  float*  coeffs0 = fCoeffs;
  //
  float *tmpCoef2D = new float[maxSpace]; // temporary workspace
  //
  for (int isl=0;isl<fNSlices;isl++) {
    for (int id=fDimOut;id--;) grids[id] = DefineGrid(isl, id, &np[2*id]);
    fCoeffsEntry[isl] = fCoeffs - coeffs0;   // offset of the given slice coeffs
    fColEntry[isl]    = fNCols  - nCols0;    // offset of the given slice columns dimensions
    for (int id=0;id<fDimOut;id++) {
      // fill function values for given output dimension in coeffs array, they will
      // be substited later by coeffs for this dimension
      FillFunValues(fun, isl, id, grids[id], &np[2*id]);
      ChebFit(&np[2*id], tmpCoef2D, prc ? prc[id] : fgkDefPrec);
      for (int id=fDimOut;id--;) delete[] grids[id];
    }
  }
  delete[] grids;
  delete[] tmpCoef2D;
  //
  fNCoefsTot = fCoeffs-coeffs0; // size of final coeffs array
  //
  fNRows = nRows0;
  // redefine arrays in compressed form, clean temp. stuff
  fNCols = new UChar_t[fNRowsTot];
  memcpy(fNCols, nCols0, fNRowsTot*sizeof(UChar_t));
  delete[] nCols0;
  fCoeffs = new Float_t[fNCoefsTot];
  memcpy(fCoeffs, coeffs0, fNCoefsTot*sizeof(Float_t));
  delete[] coeffs0;
  //
}

//____________________________________________________________________
void AliCheb2DStackF::ChebFit(const int np[2], float* tmpCoef2D, float prec)
{
  // prepare Cheb.fit for 2D -> single output dimension
  //
  int ncmax=0, maxDim = TMath::Max(np[0],np[1]);
  memset(tmpCoef2D,0,np[0]*np[1]*sizeof(float));
  //  
  float rTiny = 0.5*prec/maxDim; // neglect coefficient below this threshold
  //
  for (int id1=np[1];id1--;) { // create Cheb.param for each node of 2nd input dimension
    int nc = CalcChebCoefs(&fCoeffs[id1*np[0]], np[0], fWSpace, -1);
    for (int id0=nc;id0--;) tmpCoef2D[id1 + id0*np[1]] = fWSpace[id0];
    if (ncmax<nc) ncmax = nc;              // max coefs to be kept in dim0 to guarantee needed precision
  }
  //
  // once each 1d slice of given 2d slice is parametrized, parametrize the Cheb.coeffs
  for (int id0=np[0];id0--;) {
    CalcChebCoefs(&tmpCoef2D[id0*np[1]], np[1], fWSpace, -1);
    for (int id1=np[1];id1--;) tmpCoef2D[id1+id0*np[1]] = fWSpace[id1];
  }  
  //
  // now find 1D curve which separates significant coefficients of 2D matrix from nonsignificant ones (up to prec)
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
  // find max significant column and fill the storage 
  // for the max sigificant column of each row
  for (int id0=0;id0<nRows;id0++) {
    int nc = *fNCols++;
    for (int id1=0;id1<nc;id1++) *fCoeffs++ = tmpCoef2D[id1+id0*np[1]];
  }
  *fNRows++ = nRows;
  fNRowsTot += nRows;
  //
}

//____________________________________________________________________
void AliCheb2DStackF::FillFunValues(stFun_t fun, int slice, int dim, const float *grid, const int np[2])
{
  // fill function values on the grid
  float args[2];
  for (int id1=np[1];id1--;) {
    args[1] = grid[id1];
    for (int id0=np[0];id0--;) { // 
      args[0] = grid[np[1]+id0];
      fun(slice, args, fWSpace);
      fCoeffs[id1*np[0] + id0] = fWSpace[dim];
    }
  }
}

//__________________________________________________________________________________________
void AliCheb2DStackF::Print(const Option_t* opt) const
{
  printf("F*"); 
  AliCheb2DStack::Print(opt);
}

//__________________________________________________________________________________________
void AliCheb2DStackF::PrintSlice(int isl, const Option_t* opt) const
{
  // print slice info
  //
  TString opts = opt; opts.ToLower();
  Bool_t showcf = opts.Contains("c");
  int par0 = isl*fDimOut;
  const UChar_t* rows = &fNRows[par0];
  const UChar_t* cols = &fNCols[fColEntry[isl]];
  const Float_t* cfs  = &fCoeffs[fCoeffsEntry[isl]];
  printf("#%3d ",isl);  
  if (showcf) printf("\n");
  for (int id=0;id<fDimOut;id++) {
    int nr = *rows++;
    int ncmax=0,ncf=0;
    for (int ir=0;ir<nr;ir++) {
      ncf += cols[ir];
      if (cols[ir]>ncmax) ncmax = cols[ir];
    }
    printf("D%d: %4d coefs in %3dx%3d| ",id,ncf,nr,ncmax);
    if (showcf) {
      printf("\n");
      for (int ir=0;ir<nr;ir++) {
	for (int ic=0;ic<cols[ir];ic++) printf("%+.2e ",*cfs++); printf("\n");
      }
    }
    cols += nr; // cols entry for next dimension
    //
  }
  if (!showcf) printf("\n");
  //
}
