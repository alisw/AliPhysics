#ifndef ALI_MILLEPEDE_H 
#define ALI_MILLEPEDE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \class AliMillepede
/// \brief Class for detector alignment - modified C++ implementation of original
/// millepede package
//
//  Author: Javier Castillo

#include <TObject.h>
#include <TArrayI.h>
#include <TArrayD.h>

class AliMillepede : public TObject {

public: 
  /// Standard constructor
  AliMillepede();

  virtual ~AliMillepede( ); /// Destructor

  /// Initialization
  virtual Int_t InitMille(int nglo, int nloc, int nstd,			  
			  double lResCut, double lResCutInit);
  virtual Int_t GlobalFit(double par[], double error[], double pull[]);
  virtual Int_t SetGlobalParameters(double *param);
  virtual Int_t SetGlobalParameter(int iPar, double vPar);
  virtual Int_t SetParSigma(int iPar, double sigma);
  virtual Int_t SetNonLinear(int index);
  virtual Int_t SetGlobalConstraint(double dercs[], double rhs);
  virtual Int_t SetLocalEquation(double dergb[], double derlc[], double rmeas, double sigma);
  virtual Int_t LocalFit(int n, double localParams[], Bool_t bSingleFit);
                /// Get number of local equations
  virtual Int_t GetNLocalEquations() const {return fNLocalEquations;};
                /// Set number of local equations
  virtual void  SetNLocalEquations(int value) {fNLocalEquations = value;};
  virtual Int_t PrintGlobalParameters() const;
  virtual Int_t SetIterations (double cutfac);
  virtual Double_t GetParError(int iPar) const;

private:

// Max. dimensions

  static const int fgkMaxGlobalPar  = 1000; ///< Max. number of global parameters
  static const int fgkMaxLocalPar   = 20;   ///< Max. number of local parameters
  static const int fgkMaxGloCsts    = 20;   ///< Max. number of constraint equations
  static const int fgkMaxGloPC	    = 1020; ///< fgkMaxGlobalPar+fgkMaxGloCsts


// Private methods 

  // Double_t GetParCorrelation(int i, int j);

  int   SpmInv(double v[][fgkMaxGloPC], double b[], int n);
  int   SpmInv(double v[][fgkMaxLocalPar], double b[], int n);
  int SpAVAt(double v[][fgkMaxLocalPar], double a[][fgkMaxLocalPar], double w[][fgkMaxGlobalPar], int n, int m);
  int SpAX(double a[][fgkMaxLocalPar], double x[], double y[], int n, int m);
  double Chi2DoFLim(int n, int nd);

// Matrices

  double fMatCGlo[fgkMaxGloPC][fgkMaxGloPC];            ///< Matrix C global
  double fMatCLoc[fgkMaxLocalPar][fgkMaxLocalPar];      ///< Matrix C local
  double fMatCGloLoc[fgkMaxGlobalPar][fgkMaxLocalPar];  ///< Rectangular matrix C g*l
  double fMatCGloCorr[fgkMaxGlobalPar][fgkMaxGlobalPar];///< Correction of matrix C global
  double fMatDerConstr[fgkMaxGloCsts][fgkMaxGlobalPar]; ///< Constrained derivatives

// Vectors and useful variables

  double fDiagCGlo[fgkMaxGloPC];        ///< Initial diagonal elements of C global matrix
  double fVecBGlo[fgkMaxGloPC];         ///< Vector B global (parameters) 
  double fVecBGloCorr[fgkMaxGlobalPar]; ///< Correction of vector B global
  double fVecBLoc[fgkMaxLocalPar];      ///< Vector B local (parameters) 

  double fInitPar[fgkMaxGlobalPar];    ///< Initial global parameters
  double fDeltaPar[fgkMaxGlobalPar];   ///< Variation of global parameters 
  double fSigmaPar[fgkMaxGlobalPar];   ///< Sigma of allowed variation of global parameter 

  double fLagMult[fgkMaxGloCsts];   ///< Lagrange multipliers of constrained equations

  Bool_t fIsNonLinear[fgkMaxGlobalPar]; ///< Flag for non linear parameters
  int   fGlo2CGLRow[fgkMaxGlobalPar];   ///< Global parameter to row in "used" g*l matrix
  int   fCGLRow2Glo[fgkMaxGlobalPar];   ///< Row in "used" g*l matrix to global parameter 

  TArrayI fIndexLocEq;  ///< Table of parameter indexes in local equation 
  TArrayD fDerivLocEq;  ///< Table of local equation derivatives wrt. parameter 
  TArrayI fIndexAllEqs; ///< Index in all loc. eq. storage for iterations
  TArrayD fDerivAllEqs; ///< derivative in all loc. eq. storage for iterations
  TArrayI fLocEqPlace;  ///< Loc. Eq. position in AllEqs storage

  Int_t fNIndexLocEq;   ///< Number of entries in fIndexLocEq
  Int_t fNDerivLocEq;   ///< Number of entries in fDerivLocEq
  Int_t fNIndexAllEqs;  ///< Number of entries in fIndexAllEqs
  Int_t fNDerivAllEqs;  ///< Number of entries in fDerivAllEqs
  Int_t fNLocEqPlace;   ///< Number of entries in fLocEqPlace
  

  int    fNLocalEquations;       ///< Number of local equations 
  double fResCutInit;            ///< Cut in residual for first iterartion
  double fResCut;                ///< Cut in residual for other iterartiona

  double fChi2CutFactor;         ///< Cut factor for chi2 cut to accept local fit 
  double fChi2CutRef;            ///< Reference cut for chi2 cut to accept local fit 

  int fIter;                   ///< Current iteration
  int fMaxIter;                ///< Maximum number of iterations
  int fNStdDev;                ///< Number of standard deviations for chi2 cut 
  int fNGlobalConstraints;     ///< Number of constraint equations
  int fNLocalFits;             ///< Number of local fits
  int fNLocalFitsRejected;     ///< Number of local fits rejected
  int fNGlobalPar;             ///< Number of global parameters
  int fNLocalPar;              ///< Number of local parameters

  ClassDef(AliMillepede, 1)  // Millepede Class
};

#endif



