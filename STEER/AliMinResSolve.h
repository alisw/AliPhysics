#ifndef ALIMINRESSOLVE_H
#define ALIMINRESSOLVE_H

#include <TObject.h>
#include <TVectorD.h>
#include <TMath.h>
#include "AliMatrixSq.h"
#include "AliMatrixSparse.h"
class AliLog;
class TStopwatch;

class AliMinResSolve : public TObject {
  //
 public:
  enum {kPreconDiag=1,kPreconDILU=2,kPreconILU0=100,kPreconILU10=kPreconILU0+10,kPreconsTot};
  enum {kSolMinRes,kSolFGMRes,kNSolvers};
 public:
  AliMinResSolve();
  AliMinResSolve(const AliMatrixSq *mat, const TVectorD* rhs);
  AliMinResSolve(const AliMatrixSq *mat, const double  * rhs);
  AliMinResSolve(const AliMinResSolve& src);
  ~AliMinResSolve();
  AliMinResSolve& operator=(const AliMinResSolve& rhs);
  //
  // ---------  MINRES method (for symmetric matrices)
  Bool_t SolveMinRes(Double_t* VecSol,Int_t precon=0,int itnlim=2000,double rtol=1e-12);
  Bool_t SolveMinRes(TVectorD &VecSol,Int_t precon=0,int itnlim=2000,double rtol=1e-12);
  //
  // ---------  FGMRES method (for general symmetric matrices)
  Bool_t SolveFGMRES(Double_t* VecSol,Int_t precon=0,int itnlim=2000,double rtol=1e-12, int nkrylov=60);
  Bool_t SolveFGMRES(TVectorD &VecSol,Int_t precon=0,int itnlim=2000,double rtol=1e-12, int nkrylov=60);  
  //
  Bool_t InitAuxMinRes();
  Bool_t InitAuxFGMRES(int nkrylov);
  void   ApplyPrecon(const TVectorD& vecRHS, TVectorD& vecOut)     const;
  void   ApplyPrecon(const double*   vecRHS, double*   vecOut)     const;
  //
  Int_t  BuildPrecon(Int_t val=0);
  Int_t  GetPrecon()                                               const    {return fPrecon;} 
  void   ClearAux();
  //
  Int_t  BuildPreconILUK(Int_t lofM);
  Int_t  BuildPreconILUKDense(Int_t lofM);
  Int_t  PreconILUKsymb(Int_t lofM);
  Int_t  PreconILUKsymbDense(Int_t lofM);
  //
 protected:
  //
  Int_t               fSize;                             // dimension of the input matrix
  Int_t               fPrecon;                           // preconditioner type
  const AliMatrixSq*  fMatrix;                           // matrix defining the equations
  const Double_t*     fRHS;                              // right hand side
  //
  Double_t            *fPVecY,*fPVecR1,*fPVecR2,*fPVecV,*fPVecW,*fPVecW1,*fPVecW2;// aux MinRes
  Double_t            **fPvv,**fPvz,**fPhh;
  Double_t            *fPrecDiag,*fPrecAux; // aux space
  AliMatrixSparse     *fMatL,*fMatU;
  //
  ClassDef(AliMinResSolve,0)
};

#endif
