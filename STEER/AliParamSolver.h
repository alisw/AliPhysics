#ifndef ALIPARAMSOLVER_H
#define ALIPARAMSOLVER_H

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

#include <TObject.h>
class AliSymMatrix;

class AliParamSolver: public TObject
{
 public:
  enum {kBitGloSol=BIT(14),kBitLocSol=BIT(15),kBitCInv=BIT(16)};
  enum {kXX=0,kXY=1,kXZ=2,kYX=kXY,kYY=3,kYZ=4,kZX=kXZ,kZY=kYZ,kZZ=5};
  enum {kX,kY,kZ};

  AliParamSolver();
  AliParamSolver(Int_t maxglo,Int_t locsize=16);
  AliParamSolver(AliParamSolver& src);
  AliParamSolver& operator=(const AliParamSolver& src);
  ~AliParamSolver();
  //
  void    AddEquation(const Double_t* dGl,const Double_t *dLc,const Double_t *res,const Double_t *covI,Double_t sclErrI=-1.);
  void    AddConstraint(Int_t parID, Double_t val, Double_t err2inv);
  Bool_t  Solve(Bool_t obtainCov=kFALSE);
  Bool_t  SolveGlobals(Bool_t obtainCov=kFALSE);
  Bool_t  SolveLocals();
  void    SetMaxGlobal(Int_t n);
  void    SetNGlobal(Int_t n);
  void	  Clear(Option_t* = "");
  void	  Print(Option_t* = "") const;
  //
  Int_t   GetNGlobal()          const {return fNGlobal;}
  Int_t   GetMaxGlobal()        const {return fMaxGlobal;}
  AliSymMatrix* GetCovMatrix();
  Double_t*     GetLocals()     const {return (Double_t*)fSolLoc;}
  Double_t*     GetGlobals()    const {return (Double_t*)fSolGlo;}

 protected:
  void    Init(Int_t npini=16);
  void    ExpandStorage(Int_t newSize);

 protected:
  AliSymMatrix*    fMatrix;            // final matrix for global parameters (C in MP)
  Double_t*        fSolGlo;            // solution for globals ( vector a in MP)
  Double_t*        fSolLoc;            // solution for locals  ( vector alpha in MP)
  Int_t            fMaxGlobal;         // max number of globals can process
  Int_t            fNGlobal;           // number of globals
  Int_t            fNPoints;           // number of added points (=number of local parameters)
  //
  // temp storage
  Int_t            fMaxPoints;         // buffer size for storage
  Double_t*        fRHSGlo;            // RHS of globals (vector b in MP)
  Double_t*        fRHSLoc;            // RHS of locals (vector beta in MP)
  Double_t*        fMatGamma;          // diagonals of local partition (Gamma_i in MP)
  Double_t*        fMatG;              // off-diagonals of local partition (G_i in MP)
  Double_t*        fCovDGl;            // temporary matrix of cov*dR/dGlo
  //
  ClassDef(AliParamSolver,0)
};


#endif
