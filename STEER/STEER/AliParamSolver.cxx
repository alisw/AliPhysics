/* ----------------------------------------------------------------------------------------
   Class to solve a set of N linearized parametric equations of the type
   Eq(k): sum_i=0^n { g_i G_ik }  + t_k T_k = res_k
   whith n "global" parameters gi and one "local" parameter (per equation) t_k.
   Each measured points provides 3 measured coordinates, with proper covariance matrix.

   Used for Newton-Raphson iteration step in solution of non-linear parametric equations
   F(g,t_k) - res_k = 0, with G_ik = dF(g,t_k)/dg_i and T_k = dF(g,t_k)/dt_k
   Solution is obtained by elimination of local parameters via large (n+N) matrix partitioning 

   Author: ruben.shahoyan@cern.ch
-------------------------------------------------------------------------------------------*/ 
#include "AliParamSolver.h"
#include "AliSymMatrix.h"
#include "AliLog.h"
#include <TMath.h>

ClassImp(AliParamSolver)

//______________________________________________________________________________________
AliParamSolver::AliParamSolver()
: fMatrix(0),fSolGlo(0),fSolLoc(0),fMaxGlobal(0),fNGlobal(0),fNPoints(0),fMaxPoints(0),
  fRHSGlo(0),fRHSLoc(0),fMatGamma(0),fMatG(0),fCovDGl(0)
{ 
  // default constructor
}

//______________________________________________________________________________________
AliParamSolver::AliParamSolver(Int_t maxglo,Int_t locbuff)
: fMatrix(0),fSolGlo(0),fSolLoc(0),fMaxGlobal(maxglo),fNGlobal(maxglo),fNPoints(0),fMaxPoints(0),
  fRHSGlo(0),fRHSLoc(0),fMatGamma(0),fMatG(0),fCovDGl(0)
{ 
  // constructor for nglo globals
  Init(locbuff);
}

//______________________________________________________________________________________
AliParamSolver::AliParamSolver(const AliParamSolver& src)
  : TObject(src),fMatrix(0),fSolGlo(0),fSolLoc(0),fMaxGlobal(src.fMaxGlobal),fNGlobal(src.fNGlobal),
    fNPoints(0),fMaxPoints(0),fRHSGlo(0),fRHSLoc(0),fMatGamma(0),fMatG(0),fCovDGl(0)
{ 
  // copy constructor 
  if (src.fMatrix) {
    Init(src.fMaxPoints);
    (*this) = src;
  }
}

//______________________________________________________________________________________
AliParamSolver& AliParamSolver::operator=(const AliParamSolver& src)
{
  // assignment operator
  if (this==&src) return *this;
  TObject::operator=(src);
  if (src.fMatrix && (fNGlobal!=src.fNGlobal || fMaxPoints<src.fNPoints)) {
    fNGlobal   = src.fNGlobal;
    fMaxGlobal = src.fMaxGlobal;
    if (fMatrix)   delete   fMatrix; fMatrix = 0;
    if (fSolGlo)   delete[] fSolGlo; fSolGlo = 0;
    if (fSolLoc)   delete[] fSolLoc; fSolLoc = 0;
    if (fRHSGlo)   delete[] fRHSGlo; fRHSGlo = 0;
    if (fRHSLoc)   delete[] fRHSLoc; fRHSLoc = 0;
    if (fMatGamma) delete[] fMatGamma; fMatGamma = 0;
    if (fMatG)     delete[] fMatG; fMatG = 0;
    if (fCovDGl)   delete[] fCovDGl; fCovDGl = 0;
    Init(src.fMaxPoints);
  }
  if (src.fMatrix) {
    (*fMatrix) = *(src.fMatrix);
    memcpy(fSolGlo,src.fSolGlo,fNGlobal*sizeof(double));
    memcpy(fSolLoc,src.fSolLoc,fNPoints*sizeof(double));
    memcpy(fRHSGlo,src.fRHSGlo,fNGlobal*sizeof(double));
    memcpy(fRHSLoc,src.fRHSLoc,fNPoints*sizeof(double));
    memcpy(fMatGamma,src.fMatGamma,fNPoints*sizeof(double));
    memcpy(fMatG,src.fMatG,fNPoints*fNGlobal*sizeof(double));
  }
  //
  return *this;
}

//______________________________________________________________________________________
AliParamSolver::~AliParamSolver()
{ 
  // destructor
  delete fMatrix;
  delete[] fSolGlo;
  delete[] fSolLoc;
  delete[] fRHSGlo;
  delete[] fRHSLoc;
  delete[] fMatGamma;
  delete[] fMatG;
  delete[] fCovDGl;
}

//______________________________________________________________________________________
Bool_t AliParamSolver::SolveGlobals(Bool_t obtainCov)
{
  // solve against global vars.
  if (fNPoints<fNGlobal/3) {
    AliError(Form("Number of points: %d is not enough for %d globals",fNPoints,fNGlobal));
    return kFALSE;
  }
  //
  if (!TestBit(kBitGloSol)) {
    if (!fMatrix->SolveChol(fRHSGlo, fSolGlo, obtainCov)) {
      AliDebug(2,"Solution Failed");
      return kFALSE;
    }
    SetBit(kBitGloSol);
    if (obtainCov) SetBit(kBitCInv);
  }
  return kTRUE;
}

//______________________________________________________________________________________
Bool_t AliParamSolver::SolveLocals()
{
  // solve for locals
  const double kTiny = 1e-16;
  if (TestBit(kBitLocSol)) return kTRUE;
  if (!TestBit(kBitGloSol)) {
    AliError("Cannot solve for Locals before SolveGlobals is called");
    return kFALSE;
  }
  for (int i=fNPoints;i--;) {
    if (TMath::Abs(fMatGamma[i])<kTiny) {fSolLoc[i] = 0; continue;}
    double beta = fRHSLoc[i];
    double *mtG = fMatG + i*fNGlobal;  // G_i
    for (int j=fNGlobal;j--;) beta -= mtG[j]*fSolGlo[j];
    fSolLoc[i] = beta/fMatGamma[i];   // Gamma^-1 * (beta - G_i * a)
  }
  SetBit(kBitLocSol);
  return kTRUE;
}

//______________________________________________________________________________________
AliSymMatrix* AliParamSolver::GetCovMatrix()
{
  // obtain cov.mat
  if (!TestBit(kBitGloSol)) {
    AliError("Cannot obtain Cov.Matrix before SolveGlobals is called");
    return 0;
  }
  if (TestBit(kBitCInv)) return fMatrix;
  //
  if (fMatrix->InvertChol()) {
    SetBit(kBitCInv);
    return fMatrix;
  }
  return 0;
}

//______________________________________________________________________________________
Bool_t AliParamSolver::Solve(Bool_t obtainCov)
{
  // solve all
  return (SolveGlobals(obtainCov) && SolveLocals()) ? kTRUE : kFALSE; 
}

//______________________________________________________________________________________
void AliParamSolver::AddEquation(const Double_t* dGl,const Double_t *dLc,const Double_t *res, const Double_t *covI,Double_t sclErrI)
{
  // add the measured point to chi2 normal equations
  // Input: 
  // dGl : NGlo x 3 matrix of derivative for each of the 3 coordinates vs global params
  // dLc : 3-vector of derivative for each coordinate vs local param
  // res : residual of the point (extrapolated - measured)
  // covI: 3 x (3+1)/2 matrix of the (inverted) errors (symmetric)
  // sclErrI: scaling coefficient to apply to inverted errors (used if >0)
  // The contribution of the point to chi2 is  V * covI * V 
  // with component of the vector V(i) = dGl[i]*glob + dLc[i]*loc - res[i] 
  //
  // Instead of the NGlo + NMeasPoints size matrix we create only NGlo size matrix using the 
  // reduction a la millepede : http://www.desy.de/~blobel/millepede1.ps
  const double kTiny = 1e-16;
  //
  if (fNPoints+1 == fMaxPoints) ExpandStorage((fNPoints+1)*2);
  ResetBit(kBitGloSol|kBitLocSol);
  //
  if (TestBit(kBitCInv)) { // solution was obtained and the matrix was inverted for previous points
    fMatrix->InvertChol();
    ResetBit(kBitCInv);
  }
  //
  double *mtG   = fMatG + fNPoints*fNGlobal;  // G_i
  double &beta  = fRHSLoc[fNPoints];
  double &gamma = fMatGamma[fNPoints];
  double cDl[3];
  //
  // cov * dR/dloc
  cDl[kX] = covI[kXX]*dLc[kX] + covI[kXY]*dLc[kY] + covI[kXZ]*dLc[kZ];
  cDl[kY] = covI[kXY]*dLc[kX] + covI[kYY]*dLc[kY] + covI[kYZ]*dLc[kZ];
  cDl[kZ] = covI[kXZ]*dLc[kX] + covI[kYZ]*dLc[kY] + covI[kZZ]*dLc[kZ];
  if (sclErrI>0) { cDl[kX] *= sclErrI; cDl[kY] *= sclErrI; cDl[kZ] *= sclErrI;}
  //
  for (int i=fNGlobal;i--;) {
    const double *dGli = dGl+i*3;  // derivatives of XYZ vs i-th global param
    double *cDgi = fCovDGl+i*3;    // cov * dR/dGl_i
    cDgi[kX] = covI[kXX]*dGli[kX] + covI[kXY]*dGli[kY] + covI[kXZ]*dGli[kZ];
    cDgi[kY] = covI[kXY]*dGli[kX] + covI[kYY]*dGli[kY] + covI[kYZ]*dGli[kZ];
    cDgi[kZ] = covI[kXZ]*dGli[kX] + covI[kYZ]*dGli[kY] + covI[kZZ]*dGli[kZ];
    if (sclErrI>0) { cDgi[kX] *= sclErrI; cDgi[kY] *= sclErrI; cDgi[kZ] *= sclErrI;}
    //
    mtG[i]   = cDl[kX]*dGli[kX] + cDl[kY]*dGli[kY] + cDl[kZ]*dGli[kZ];  // dR/dGl_i * cov * dR/dLoc
  }
  beta    = res[kX]*cDl[kX] + res[kY]*cDl[kY] + res[kZ]*cDl[kZ];  //RHS: res*cov*dR/dloc
  gamma   = dLc[kX]*cDl[kX] + dLc[kY]*cDl[kY] + dLc[kZ]*cDl[kZ];  //RHS: dR/dloc*cov*dR/dloc
  double locSol  = TMath::Abs(gamma)<kTiny ? 0. : beta/gamma;     //local solution: gamma^-1 beta
  //
  AliSymMatrix &matC = *fMatrix;
  for (int i=fNGlobal;i--;) {
    const double *cDgi = fCovDGl+i*3;    // cov * dR/dGl_i
    //
    fRHSGlo[i] += cDgi[kX]*res[kX] + cDgi[kY]*res[kY] + cDgi[kZ]*res[kZ]      // b_i = dR/dGi * cov * res
      -           mtG[i]*locSol;                                              //     - [G gamma^-1 beta ]_i
    //
    for (int j=i+1;j--;) {
      //      const double *cDgj = fCovDGl+j*3;  // cov * dR/dGl_j
      const double *dGlj = dGl+j*3;        // derivatives of XYZ vs i-th global param
      double add = dGlj[kX]*cDgi[kX] + dGlj[kY]*cDgi[kY] + dGlj[kZ]*cDgi[kZ];  // C_ij = dR/dGi * cov * dR/dGj
      if (TMath::Abs(gamma)>kTiny) add -= mtG[i]*mtG[j]/gamma;                 //        - [G gamma^-1 T(G) ]_ij      
      matC(i,j) += add;       
    }
  }
  //
  fNPoints++;
  //
}

//______________________________________________________________________________________
void AliParamSolver::AddConstraint(Int_t parID, Double_t val, Double_t err2inv)
{
  // add gassian constriant to parameter parID
  if (parID>=fNGlobal) {
    AliError(Form("Attempt to constraint non-existing parameter %d",parID));
    return;
  }
  //
  (*fMatrix)(parID,parID) += err2inv;
  fRHSGlo[parID] += val*err2inv;
  //
}

//______________________________________________________________________________________
void AliParamSolver::Init(Int_t npini)
{
  // create storage assuming maximum npini measured points
  fMatrix = new AliSymMatrix(fMaxGlobal);
  fSolGlo = new Double_t[fMaxGlobal];
  fRHSGlo = new Double_t[fMaxGlobal];
  //
  fCovDGl = new Double_t[3*fMaxGlobal];
  ExpandStorage(npini);
  fMatrix->SetSizeUsed(fNGlobal);
  //
}

//______________________________________________________________________________________
void AliParamSolver::ExpandStorage(Int_t newSize)
{
  // increase space to newSize measured points
  newSize = newSize>fMaxPoints ? newSize : fMaxPoints+1;
  double* tmp;
  tmp = new Double_t[newSize];
  if (fSolLoc) delete[] fSolLoc; 
  fSolLoc = tmp;
  //
  tmp = new Double_t[newSize];
  if (fMatGamma) {
    memcpy(tmp,fMatGamma,fNPoints*sizeof(Double_t));
    delete[] fMatGamma;
  }
  fMatGamma = tmp;
  //
  tmp = new Double_t[newSize];
  if (fRHSLoc) {
    memcpy(tmp,fRHSLoc,fNPoints*sizeof(Double_t));
    delete[] fRHSLoc;
  }
  fRHSLoc = tmp;
  //
  tmp = new Double_t[newSize*fMaxGlobal];
  if (fMatG) {
    memcpy(tmp,fMatG,fNPoints*fMaxGlobal*sizeof(Double_t));
    delete[] fMatG;
  }
  fMatG = tmp;
  //
  fMaxPoints = newSize;
  //
}

//______________________________________________________________________________________
void AliParamSolver::Clear(Option_t*)
{
  // reset all
  fNPoints = 0;
  fMatrix->Reset();
  for (int i=fNGlobal;i--;) fRHSGlo[i] = 0;
  ResetBit(kBitGloSol | kBitLocSol | kBitCInv);
}

//______________________________________________________________________________________
void AliParamSolver::Print(Option_t*) const
{
  // print itself
  AliInfo(Form("Solver with %d globals for %d points",fNGlobal,fNPoints));
}

//______________________________________________________________________________________
void AliParamSolver::SetNGlobal(Int_t n) 
{
  // set N global params
  if (n>fMaxGlobal) {
    AliError(Form("Maximum number of globals was set to %d",fMaxGlobal));
    return;
  }
  fNGlobal = n;
  fMatrix->SetSizeUsed(fNGlobal);
}

//______________________________________________________________________________________
void AliParamSolver::SetMaxGlobal(Int_t n) 
{
  // set limit on N glob.
  if (n>0 && n==fMaxGlobal) return;
  fMaxGlobal = n;
  fNGlobal = n;
  if (fMatrix)   delete   fMatrix;   fMatrix = 0;
  if (fSolGlo)   delete[] fSolGlo;   fSolGlo = 0;
  if (fSolLoc)   delete[] fSolLoc;   fSolLoc = 0;
  if (fRHSGlo)   delete[] fRHSGlo;   fRHSGlo = 0;
  if (fRHSLoc)   delete[] fRHSLoc;   fRHSLoc = 0;
  if (fMatGamma) delete[] fMatGamma; fMatGamma = 0;
  if (fMatG)     delete[] fMatG;     fMatG = 0;
  if (fCovDGl)   delete[] fCovDGl;   fCovDGl = 0;
  n = TMath::Max(16,fMaxPoints);
  Init(n);
  //
}

