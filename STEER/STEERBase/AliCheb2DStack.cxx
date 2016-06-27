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


#include "AliCheb2DStack.h"
#include "AliLog.h"
#include <TMath.h>

ClassImp(AliCheb2DStack)

float AliCheb2DStack::fgkDefPrec = 1e-4;
float AliCheb2DStack::fWSpace[AliCheb2DStack::kMaxPoints] = {0};
//____________________________________________________________________
AliCheb2DStack::AliCheb2DStack() 
  :fDimOut(0) 
  ,fNSlices(0)
  ,fNParams(0)
  ,fNCoefsTot(0)
  ,fNRowsTot(0)
  ,fRowXI(0)
  ,fNRows(0)
  ,fNCols(0)
  ,fCoeffsEntry(0)
  ,fColEntry(0)
{
  // Default constructor
  for (int i=2;i--;) fBMin[i] = fBMax[i] = 0;
  fBScaleZ = fBOffsetZ = 0;
  fDead[0] = fDead[1] = 0;
}

//____________________________________________________________________
AliCheb2DStack::AliCheb2DStack(int nSlices, int dimOut, const float bmin[2],const float bmax[2], 
			       const float* dead, const float *rowXI)
  :fDimOut(dimOut) 
  ,fNSlices(nSlices)
  ,fNParams(nSlices*dimOut)
  ,fNCoefsTot(0)
  ,fNRowsTot(0)
  ,fRowXI(rowXI)
  ,fNRows(new UChar_t[nSlices*dimOut])
  ,fNCols(0)
  ,fCoeffsEntry(new Int_t[nSlices])
  ,fColEntry(new Int_t[nSlices])
{
  // create stack of 2D->dimOut Chebyshev parameterizations debined in 2 dimensions between bmin and bmax,
  // and trained with function fun on 2D grid on np points. 
  // Truncate each precision of each output dimension parameterization i to precD[i] if precD!=0, or prec
  for (int i=2;i--;) {
    fBMin[i] = bmin[i];
    fBMax[i] = bmax[i];
  }
  fBScaleZ = 2./(fBMax[1] - fBMin[1]); // prepare mapping of boundaries to [-1:1]
  fBOffsetZ = fBMin[1] + 1./fBScaleZ;

  if (dead) {
    fDead[0] = dead[0];
    fDead[1] = dead[1];
    if (!rowXI && (TMath::Abs(fDead[0])>fgkDefPrec || TMath::Abs(fDead[1])>fgkDefPrec)  ) 
      AliErrorF("Dead zones %.2f/%.2f will be ignored since the inverse radii are not provided",fDead[0],fDead[1]);
  }
}

//____________________________________________________________________
AliCheb2DStack::~AliCheb2DStack()
{
  // D-tor
  delete[] fNRows;
  delete[] fNCols;
  delete[] fCoeffsEntry;
  delete[] fColEntry;
}

//____________________________________________________________________
void AliCheb2DStack::CheckDimensions(const int *np) const
{
  // basic consistency check
  if (fDimOut>=kMaxPoints) AliFatalF("N output dimensions=%d > %d",fDimOut,kMaxPoints);
  for (int id=fDimOut;id--;) {
    for (int i=2;i--;) {
      if (np[id*2+i]<1) AliFatalF("N points=%d in input dim. %d is <1 for output dim.%d",
				  np[id*2+i],i,id);
      if (np[id*2+i]<1) AliFatalF("N points=%d in input dim. %d is >%d for output dim.%d",
				  np[id*2+i],i,kMaxPoints,id);
    }
  }
  //  
  for (int i=2;i--;) 
    if (fBMin[i]>=fBMax[i]) AliFatalF("Boundaries for %d-th dim. are not"
				      " increasing: %+.4e %+.4e",i,fBMin[i],fBMax[i]);
  //
}

//____________________________________________________________________
float* AliCheb2DStack::DefineGrid(int slice, int dim, const int np[2]) const
{
  // prepare the grid of Chebyshev roots for dim-th output dimension
  // First np[1] nodes of 2nd input dimesion are stored, then np[0] nodes for 1st dim.
  const int kMinPoints = 1;
  float *grid = new float[np[0]+np[1]];
  int cnt=0;
  for (int id=2;id--;) {
    int npnt = np[id];
    for (int ip=0;ip<npnt;ip++) {
      float x = TMath::Cos( TMath::Pi()*(ip+0.5)/npnt );
      grid[cnt++] = MapToExternal(slice,x,id);
    }
  }
  //
  return grid;
}

//____________________________________________________________________
Int_t AliCheb2DStack::CalcChebCoefs(const float *funval,int np, float *outCoefs, float prec)
{
  // Calculate Chebyshev coeffs using precomputed function values at np roots.
  // If prec>0, estimate the highest coeff number providing the needed precision
  //
  double sm;                 // do summations in double to minimize the roundoff error
  for (int ic=0;ic<np;ic++) { // compute coeffs
    sm = 0;          
    for (int ir=0;ir<np;ir++) sm += funval[ir]*TMath::Cos( ic*(ir+0.5)*TMath::Pi()/np);
    outCoefs[ic] = float( sm * ((ic==0) ? 1./np : 2./np) );
  }
  if (prec<=0) return np;
  //
  sm = 0;
  int cfMax = 0;
  for (cfMax=np;cfMax--;) {
    sm += TMath::Abs(outCoefs[cfMax]);
    if (sm>=prec) break;
  }
  if (++cfMax==0) cfMax=1;
  return cfMax;
  //
}

//__________________________________________________________________________________________
void AliCheb2DStack::Print(const Option_t* opt) const
{
  // print full info
  //
  printf("Cheb.param for %dx2D->%dD in [%+.3e:%+.3e] [%+.3e:%+.3e] | %d coefs in %d rows\n",
	 fNSlices,fDimOut,fBMin[0],fBMax[0],fBMin[1],fBMax[1],fNCoefsTot,fNRowsTot);
  TString opts = opt; opts.ToLower();
  if (opts.Contains("l")) for (int isl=0;isl<fNSlices;isl++) PrintSlice(isl,opt);
  //
}
