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

/* $Id$ */

//-----------------------------------------------------------------------------
/// \class AliMillepede
/// Detector independent alignment class
///
/// This modified C++ version is based on a C++ translation of Millepede used 
/// for LHCb Vertex Detector alignment (lhcb-2005-101), available here 
/// http://isscvs.cern.ch/cgi-bin/cvsweb.cgi/Alignment/AlignmentTools/src/?cvsroot=lhcb
/// The original millepede fortran package is available at:  
/// http://www.desy.de/~blobel/wwwmille.html  
///
/// \author Javier Castillo
//-----------------------------------------------------------------------------

#include <TArrayI.h>
#include <TArrayD.h>
#include <TMath.h>

#include "AliLog.h"

#include "AliMillepede.h"

//=============================================================================
AliMillepede::AliMillepede()
  : TObject(),
    fIndexLocEq(fgkMaxGlobalPar+fgkMaxLocalPar),
    fDerivLocEq(fgkMaxGlobalPar+fgkMaxLocalPar),
    fIndexAllEqs(1000*fgkMaxGlobalPar+fgkMaxLocalPar),
    fDerivAllEqs(1000*fgkMaxGlobalPar+fgkMaxLocalPar),
    fLocEqPlace(1000), 
    fNIndexLocEq(0),
    fNDerivLocEq(0),
    fNIndexAllEqs(0),
    fNDerivAllEqs(0),
    fNLocEqPlace(0),
    fNLocalEquations(0),
    fResCutInit(0.),
    fResCut(0.),
    fChi2CutFactor(1.0),
    fChi2CutRef(1.0),
    fIter(0),
    fMaxIter(10),
    fNStdDev(3),
    fNGlobalConstraints(0),
    fNLocalFits(0),
    fNLocalFitsRejected(0),
    fNGlobalPar(0),
    fNLocalPar(0)
{
  /// Standard constructor

  AliInfo("                                           ");
  AliInfo("            * o o                   o      ");
  AliInfo("              o o                   o      ");
  AliInfo("   o ooooo  o o o  oo  ooo   oo   ooo  oo  ");
  AliInfo("    o  o  o o o o o  o o  o o  o o  o o  o ");
  AliInfo("    o  o  o o o o oooo o  o oooo o  o oooo ");
  AliInfo("    o  o  o o o o o    ooo  o    o  o o    ");
  AliInfo("    o  o  o o o o  oo  o     oo   ooo  oo  ++ starts");	   
  AliInfo("                                           ");

}

//=============================================================================
AliMillepede::~AliMillepede() {
  /// Destructor

}

//=============================================================================
Int_t AliMillepede::InitMille(int nGlo, int nLoc, int lNStdDev,
			      double lResCut, double lResCutInit)
{
  /// Initialization of millepede
  AliDebug(1,"");
  AliDebug(1,"----------------------------------------------------");
  AliDebug(1,"");
  AliDebug(1,"    Entering InitMille");
  AliDebug(1,"");
  AliDebug(1,"-----------------------------------------------------");
  AliDebug(1,"");

  fNGlobalConstraints = 0;
  fNLocalFits  = 0;                      // Total number of local fits
  fNLocalFitsRejected  = 0;              // Total number of local fits rejected
  fChi2CutRef  = 1.0;                    // Reference value for Chi^2/ndof cut

  AliMillepede::SetNLocalEquations(0);     // Number of local fits (starts at 0)

  fIter    = 0;  // By default iterations are turned off, turned on by use SetIterations
  fMaxIter = 10;

  fResCutInit = lResCutInit; 
  fResCut = lResCut;

  fNGlobalPar = nGlo;     // Number of global derivatives
  fNLocalPar  = nLoc;     // Number of local derivatives
  fNStdDev    = lNStdDev; // Number of StDev for local fit chisquare cut

  AliDebug(1,Form("Number of global parameters   : %d", fNGlobalPar));
  AliDebug(1,Form("Number of local parameters    : %d", fNLocalPar));
  AliDebug(1,Form("Number of standard deviations : %d", fNStdDev));

  if (fNGlobalPar>fgkMaxGlobalPar || fNLocalPar>fgkMaxLocalPar) {
    AliDebug(1,"Two many parameters !!!!!");
    return 0;
  }

  // Global parameters initializations
  for (int i=0; i<fNGlobalPar; i++) {  
    fVecBGlo[i]=0.;
    fInitPar[i]=0.;
    fDeltaPar[i]=0.;
    fSigmaPar[i]=-1.;
    fIsNonLinear[i]=0;
    fGlo2CGLRow[i]=-1;
    fCGLRow2Glo[i]=-1;
    
    for (int j=0; j<fNGlobalPar;j++) {    
      fMatCGlo[i][j]=0.;
    }
  }

  // Local parameters initializations
  
  for (int i=0; i<fNLocalPar; i++) {
    fVecBLoc[i]=0.;
    
    for (int j=0; j<fNLocalPar;j++) {
      fMatCLoc[i][j]=0.;
    }
  }

  // Then we fix all parameters...

  for (int j=0; j<fNGlobalPar; j++) {
    AliMillepede::SetParSigma(j,0.0);
  }

  fDerivLocEq.Reset();  fNDerivLocEq=0;  
  fIndexLocEq.Reset();  fNIndexLocEq=0; 

  fIndexAllEqs.Reset();  fNIndexAllEqs=0;
  fDerivAllEqs.Reset();  fNDerivAllEqs=0;
  fLocEqPlace.Reset();  fNLocEqPlace=0;

  AliDebug(1,"");
  AliDebug(1,"----------------------------------------------------");
  AliDebug(1,"");
  AliDebug(1,"    InitMille has been successfully called!");
  AliDebug(1,"");
  AliDebug(1,"-----------------------------------------------------");
  AliDebug(1,"");
	
  return 1;
}

/*
-----------------------------------------------------------
  PARGLO: initialization of global parameters
-----------------------------------------------------------

  param    = array of starting values

-----------------------------------------------------------
*/
Int_t AliMillepede::SetGlobalParameters(double *param)
{
  /// initialization of global parameters
  for(Int_t iPar=0; iPar<fNGlobalPar; iPar++){
    fInitPar[iPar] = param[iPar];
  }

 return 1;
}

/*
-----------------------------------------------------------
  PARGLO: initialization of global parameters
-----------------------------------------------------------

  iPar    = the index of the global parameter in the 
             result array (equivalent to fDeltaPar[]).

  param    = the starting value

-----------------------------------------------------------
*/
Int_t AliMillepede::SetGlobalParameter(int iPar, double param)
{
  /// initialization of global parameter iPar
  if (iPar<0 || iPar>=fNGlobalPar) {
    return 0;
  }
  else {
    fInitPar[iPar] = param;
  }
  return 1;
}


/*
-----------------------------------------------------------
  PARSIG: define a constraint for a single global param
          param is 'encouraged' to vary within [-sigma;sigma] 
	  range
-----------------------------------------------------------

  iPar    = the index of the global parameter in the 
             result array (equivalent to fDeltaPar[]).

  sigma	   = value of the constraint (sigma <= 0. will 
             mean that parameter is FIXED !!!) 
 
-----------------------------------------------------------
*/ 
Int_t AliMillepede::SetParSigma(int iPar, double sigma)
{
  /// Define a range [-sigma;sigma] where iPar is encourage to vary 
  if (iPar>=fNGlobalPar) {
    return 0;
  }
  else {
    fSigmaPar[iPar] = sigma;
  }

  return 1;
}


/*
-----------------------------------------------------------
  NONLIN: set nonlinear flag for a single global param
          update of param durin iterations will not
	  consider initial starting value
-----------------------------------------------------------

  iPar    = the index of the global parameter in the 
             result array (equivalent to fDeltaPar[]).

-----------------------------------------------------------
*/
Int_t AliMillepede::SetNonLinear(int iPar)
{
  /// Set nonlinear flag for iPar 
  if (iPar<0 || iPar>=fNGlobalPar) {
    return 0;
  }
  else {
    fIsNonLinear[iPar] = 1;
  }
  
  return 1;
}


/*
-----------------------------------------------------------
  INITUN: unit for iteration
-----------------------------------------------------------
  
  lChi2CutFac is used by Fitloc to define the Chi^2/ndof cut value

  A large cutfac value enables to take a wider range of tracks 
  for first iterations, which might be useful if misalignments
  are large.

  As soon as cutfac differs from 0 iteration are requested.
  cutfac is then reduced, from one iteration to the other,
  and iterations are stopped when it reaches the value 1.

  At least one more iteration is often needed in order to remove
  tracks containing outliers.
  
-----------------------------------------------------------
*/
Int_t AliMillepede::SetIterations(double lChi2CutFac)
{
  /// Number of iterations is calculated from lChi2CutFac 
  fChi2CutFactor = TMath::Max(1.0, lChi2CutFac);

  AliInfo(Form("Initial cut factor is %f",fChi2CutFactor));
  fIter = 1; // Initializes the iteration process
  return 1;
}

/*
-----------------------------------------------------------
  CONSTF: define a constraint equation in AliMillepede
-----------------------------------------------------------

  dercs    = the row containing constraint equation 
             derivatives (put into the final matrix)

  lLagMult      = the lagrange multiplier value (sum of equation)	     

-----------------------------------------------------------
*/
Int_t AliMillepede::SetGlobalConstraint(double dercs[], double lLagMult)
{ 
  /// Define a constraint equation
  if (fNGlobalConstraints>=fgkMaxGloCsts) {  
    AliInfo("Too many constraints !!!");
    return 0;
  }
 	
  for (int i=0; i<fNGlobalPar; i++) {
    fMatDerConstr[fNGlobalConstraints][i] = dercs[i];
  }
 	
  fLagMult[fNGlobalConstraints] = lLagMult;
  fNGlobalConstraints++ ;
  AliInfo(Form("Number of constraints increased to %d",fNGlobalConstraints));
  return 1;
}


/*
-----------------------------------------------------------
  EQULOC: write ONE equation in the matrices
-----------------------------------------------------------

  dergb[1..fNGlobalPar]	= global parameters derivatives
  derlc[1..fNLocalPar] 	= local parameters derivatives
  rmeas  		= measured value
  sigma 		= error on measured value (nothing to do with SetParSigma!!!)

-----------------------------------------------------------
*/
Int_t AliMillepede::SetLocalEquation(double dergb[], double derlc[], double lMeas, double lSigma)
{	
  /// Write one local equation
  if (lSigma<=0.0) { // If parameter is fixed, then no equation
    for (int i=0; i<fNLocalPar; i++) {
      derlc[i] = 0.0;
    }
    for (int i=0; i<fNGlobalPar; i++) {
      dergb[i] = 0.0;
    }
    return 1;
  }
    
  // Serious equation, initialize parameters   	
  double lWeight =  1.0/(lSigma*lSigma);
  int iLocFirst  = -1;
  int iLocLast   = -1;
  int iGlobFirst = -1;
  int iGlobLast  = -1;
 
  for (int i=0; i<fNLocalPar; i++) { // Retrieve local param interesting indices
    if (derlc[i]!=0.0) {
      if (iLocFirst == -1) { 
	iLocFirst = i;	// first index
      }
      iLocLast = i;     // last index
    }
  }
  AliDebug(2,Form("%d / %d",iLocFirst, iLocLast));
	
  for (int i=0; i<fNGlobalPar; i++) { // Idem for global parameters
    if (dergb[i]!=0.0) {
      if (iGlobFirst == -1) {
	iGlobFirst = i;  // first index
      }
      iGlobLast = i; 	 // last index
    }
  }
  AliDebug(2,Form("%d / %d",iGlobFirst,iGlobLast));

  if (fNIndexLocEq==fIndexLocEq.GetSize()) fIndexLocEq.Set(2*fNIndexLocEq);
  fIndexLocEq.AddAt(-1,fNIndexLocEq++);    
  if (fNDerivLocEq==fDerivLocEq.GetSize()) fDerivLocEq.Set(2*fNDerivLocEq);
  fDerivLocEq.AddAt(lMeas,fNDerivLocEq++);    


  for (int i=iLocFirst; i<=iLocLast; i++) { // Store interesting local parameters
    if (derlc[i]!=0.0) {
      if (fNIndexLocEq==fIndexLocEq.GetSize()) fIndexLocEq.Set(2*fNIndexLocEq);
      fIndexLocEq.AddAt(i,fNIndexLocEq++);    
      if (fNDerivLocEq==fDerivLocEq.GetSize()) fDerivLocEq.Set(2*fNDerivLocEq);
      fDerivLocEq.AddAt(derlc[i],fNDerivLocEq++);    
      derlc[i] = 0.0;
    }
  }

  if (fNIndexLocEq==fIndexLocEq.GetSize()) fIndexLocEq.Set(2*fNIndexLocEq);
  fIndexLocEq.AddAt(-1,fNIndexLocEq++);    
  if (fNDerivLocEq==fDerivLocEq.GetSize()) fDerivLocEq.Set(2*fNDerivLocEq);
  fDerivLocEq.AddAt(lWeight,fNDerivLocEq++);    

  for (int i=iGlobFirst; i<=iGlobLast; i++) { // Store interesting global parameters
    if (dergb[i]!=0.0) {
      if (fNIndexLocEq==fIndexLocEq.GetSize()) fIndexLocEq.Set(2*fNIndexLocEq);
      fIndexLocEq.AddAt(i,fNIndexLocEq++);    
      if (fNDerivLocEq==fDerivLocEq.GetSize()) fDerivLocEq.Set(2*fNDerivLocEq);
      fDerivLocEq.AddAt(dergb[i],fNDerivLocEq++);    
      dergb[i] = 0.0;
    }
  }
  
  AliDebug(2,Form("Out Equloc --  NST = %d",fNDerivLocEq));
  return 1; 	
}

/*
-----------------------------------------------------------
  FITLOC:  perform local params fit, once all the equations
           have been written by EquLoc
-----------------------------------------------------------

  iFit        = number of the fit, it is used to store 
                fit parameters and then retrieve them 
		for iterations (via FINDEXALLEQS and FDERIVALLEQS)

  localParams = contains the fitted track parameters and
                related errors

  bSingleFit  = is an option, if it is set to 1, we don't 
                perform the last loop. It is used to update 
		the track parameters without modifying global
		matrices

-----------------------------------------------------------
*/
Int_t AliMillepede::LocalFit(int iFit, double localParams[], bool bSingleFit)
{
  /// Perform local parameters fit once all the local equations have been set

  // Few initializations
  int iEqTerm = 0;
  int i, j, k;
  int iIdx, jIdx, kIdx;
  int iIdxIdx, jIdxIdx;
  int iMeas   = -1;
  int iWeight = 0;
  int iLocFirst = 0;
  int iLocLast  = 0;
  int iGloFirst = 0;
  int iGloLast  = 0;
  int nGloInFit = 0;
	
  double lMeas    = 0.0;
  double lWeight  = 0.0;

  double lChi2    = 0.0;
  double lRedChi2 = 0.0;
  double lChi2Cut = 0.0;
  int nEq  = 0;
  int nDoF = 0;
  int nEqTerms = fNDerivLocEq; // Number of terms (local + global derivatives, 		 
			       // measurement and weight) involved in this local equation


  // Fill the track store at first pass

  if (fIter < 2 && !bSingleFit) { // Do it only once 
    AliDebug(1,Form("Store equation no: %d", iFit)); 

    for (i=0; i<nEqTerms; i++) { // Store the track parameters
      if (fNIndexAllEqs==fIndexAllEqs.GetSize()) fIndexAllEqs.Set(2*fNIndexAllEqs);
      fIndexAllEqs.AddAt(fIndexLocEq[i],fNIndexAllEqs++);
      if (fNDerivAllEqs==fDerivAllEqs.GetSize()) fDerivAllEqs.Set(2*fNDerivAllEqs);
      fDerivAllEqs.AddAt(fDerivLocEq[i],fNDerivAllEqs++);
    }
    if (fNLocEqPlace==fLocEqPlace.GetSize()) fLocEqPlace.Set(2*fNLocEqPlace);
    fLocEqPlace.AddAt(fNIndexAllEqs,fNLocEqPlace++);

    
    AliDebug(2,Form("FLocEqPlace size = %d",fLocEqPlace[iFit])); 
    AliDebug(2,Form("FIndexAllEqs size   = %d",fNIndexAllEqs)); 
  }


  for (i=0; i<fNLocalPar; i++) { // reset local params
    fVecBLoc[i] = 0.0;

    for (j=0; j<fNLocalPar; j++) {
      fMatCLoc[i][j] = 0.0;
    }
  }
	
  for (i=0; i<fNGlobalPar; i++) {
    fGlo2CGLRow[i] = -1;  // reset mixed params
  }


/*

  LOOPS : HOW DOES IT WORKS ?	

  Now start by reading the informations stored with EquLoc.
  Those informations are in vector FINDEXSTORE and FDERIVSTORE.
  Each -1 in FINDEXSTORE delimits the equation parameters:
  
  First -1  ---> lMeas in FDERIVSTORE 
  Then we have indices of local eq in FINDEXSTORE, and derivatives in FDERIVSTORE
  Second -1 ---> weight in FDERIVSTORE
  Then follows indices and derivatives of global eq.
  ....
  
  We took them and store them into matrices.
  
  As we want ONLY local params, we substract the part of the estimated value
  due to global params. Indeed we could have already an idea of these params,
  with previous alignment constants for example (set with PARGLO). Also if there
  are more than one iteration (FITLOC could be called by FITGLO)

*/

    
//
// FIRST LOOP : local track fit
//
	
  iEqTerm = 0;
//   fIndexLocEq.push_back(-1);
  if (fNIndexLocEq==fIndexLocEq.GetSize()) fIndexLocEq.Set(2*fNIndexLocEq);
  fIndexLocEq.AddAt(-1,fNIndexLocEq++);    
  
  while (iEqTerm <= nEqTerms) {
    if (fIndexLocEq[iEqTerm] == -1) {
      if (iMeas == -1) {        // First  -1 : lMeas
	iMeas = iEqTerm;
	iLocFirst = iEqTerm+1;
      }  
      else if (iWeight == 0) {  // Second -1 : weight 
	iWeight = iEqTerm;
	iLocLast = iEqTerm-1;
	iGloFirst = iEqTerm+1;
      } 
      else {                    // Third  -1 : end of equation; start of next  
	iGloLast = iEqTerm-1;
	lMeas	= fDerivLocEq[iMeas];
	lWeight 	= fDerivLocEq[iWeight];
// 	AliDebug(1,Form("lMeas = %f", lMeas));
// 	AliDebug(1,Form("lWeight = %f", lWeight));
        
	// Now suppress the global part (only relevant with iterations)
	// 
	for (i=iGloFirst; i<=iGloLast; i++) {
	  iIdx = fIndexLocEq[i];              // Global param indice
// 	  AliDebug(2,Form("fDeltaPar[%d] = %f", iIdx, fDeltaPar[iIdx]));        
// 	  AliDebug(2,Form("Starting misalignment = %f",fInitPar[iIdx]));        
	  if (fIsNonLinear[iIdx] == 0)
	    lMeas -= fDerivLocEq[i]*(fInitPar[iIdx]+fDeltaPar[iIdx]); // linear parameter
	  else
	    lMeas -= fDerivLocEq[i]*(fDeltaPar[iIdx]); // nonlinear parameter
	}
// 	AliDebug(2,Form("lMeas after global stuff removal = %f", lMeas));
				
	for (i=iLocFirst; i<=iLocLast; i++) { // Finally fill local matrix and vector
	  iIdx = fIndexLocEq[i];   // Local param indice (the matrix line) 
	  fVecBLoc[iIdx] += lWeight*lMeas*fDerivLocEq[i];  
// 	  AliDebug(2,Form("fVecBLoc[%d] = %f", iIdx, fVecBLoc[iIdx]));
					
	  for (j=iLocFirst; j<=i ; j++) { // Symmetric matrix, don't bother j>i coeffs
	    jIdx = fIndexLocEq[j];						
	    fMatCLoc[iIdx][jIdx] += lWeight*fDerivLocEq[i]*fDerivLocEq[j];	    
// 	    AliDebug(2,Form("fMatCLoc[%d][%d] = ", iIdx, jIdx, fMatCLoc[iIdx][jIdx]));
	  }
	}
	iMeas   = -1;
	iWeight = 0;
	iEqTerm--;   // end of one equation is the beginning of next
      } // End of "end of equation" operations
    } // End of loop on equation
    iEqTerm++;
  } // End of loop on all equations used in the fit


//
// Local params matrix is completed, now invert to solve...
//
	
  Int_t nRank = AliMillepede::SpmInv(fMatCLoc, fVecBLoc, fNLocalPar);
  // nRank is the number of nonzero diagonal elements 
     	
  AliDebug(1,"");
  AliDebug(1," __________________________________________________");
  AliDebug(1,Form(" Printout of local fit  (FITLOC)  with rank= %d", nRank));
  AliDebug(1," Result of local fit :      (index/parameter/error)");
  
  for (i=0; i<fNLocalPar; i++) {
    AliDebug(1,Form("%d   /   %.6f   /   %.6f", i, fVecBLoc[i], TMath::Sqrt(fMatCLoc[i][i])));	
  }
  

// Store the track params and errors

  for (i=0; i<fNLocalPar; i++) {
    localParams[2*i] = fVecBLoc[i];
    localParams[2*i+1] = TMath::Sqrt(TMath::Abs(fMatCLoc[i][i]));
  }

    
//
// SECOND LOOP : residual calculation
//
  
  iEqTerm = 0;
  iMeas = -1;
  iWeight = 0;

  while (iEqTerm <= nEqTerms) {
    if (fIndexLocEq[iEqTerm] == -1) {
      if (iMeas == -1) {        // First  -1 : lMeas
	iMeas = iEqTerm;
	iLocFirst = iEqTerm+1;
      }  
      else if (iWeight == 0) {  // Second -1 : weight 
	iWeight = iEqTerm;
	iLocLast = iEqTerm-1;
	iGloFirst = iEqTerm+1;
      } 
      else {                    // Third  -1 : end of equation; start of next  
	iGloLast = iEqTerm-1;
	lMeas	= fDerivLocEq[iMeas];
	lWeight 	= fDerivLocEq[iWeight];
	
	// Print all (for debugging purposes)

// 	int nDerLoc = iLocLast-iLocFirst+1;   // Number of local derivatives involved
// 	int nDerGlo = iGloLast-iGloFirst+1;   // Number of global derivatives involved

// 	AliDebug(2,"");
// 	AliDebug(2,Form(". equation:  measured value %.6f +/- %.6f", lMeas, 1.0/TMath::Sqrt(lWeight)));
// 	AliDebug(2,Form("Number of derivatives (global, local): %d, %d",nDerGlo,nDerLoc));
// 	AliDebug(2,"Global derivatives are: (index/derivative/parvalue) ");
	
// 	for (i=iGloFirst; i<=iGloLast; i++) {
// 	  AliDebug(2,Form("%d / %.6f / %.6f",fIndexLocEq[i],fDerivLocEq[i],fInitPar[fIndexLocEq[i]]));
//      } 

// 	AliDebug(2,"Local derivatives are: (index/derivative) ");
	
// 	for (i=(ja+1); i<jb; i++) {AliDebug(2,Form("%d / %.6f",fIndexLocEq[i], fDerivLocEq[i]));}	  

	// Now suppress local and global parts to LMEAS;
	//
	// First the local part 
	for (i=iLocFirst; i<=iLocLast; i++) { 
	  iIdx = fIndexLocEq[i];
	  lMeas -= fDerivLocEq[i]*fVecBLoc[iIdx];
	}
	// Then the global part
	for (i=iGloFirst; i<=iGloLast; i++) {
	  iIdx = fIndexLocEq[i];
	  if (fIsNonLinear[iIdx] == 0)
	    lMeas -= fDerivLocEq[i]*(fInitPar[iIdx]+fDeltaPar[iIdx]); // linear parameter
	  else
	    lMeas -= fDerivLocEq[i]*(fDeltaPar[iIdx]); // nonlinear parameter
	}

	// lMeas contains now the residual value
	AliDebug(2,Form("Residual value : %.6f", lMeas));

	// reject the track if lMeas is too important (outlier)
	if (TMath::Abs(lMeas) >= fResCutInit && fIter <= 1) {
	  AliDebug(2,"Rejected track !!!!!");
    	  fNLocalFitsRejected++;      
	  fIndexLocEq.Reset();  fNIndexLocEq=0; // reset stores and go to the next track 
	  fDerivLocEq.Reset();  fNDerivLocEq=0;	  
	  return 0;
	}

	if (TMath::Abs(lMeas) >= fResCut && fIter > 1) {
	  AliDebug(2,"Rejected track !!!!!");
    	  fNLocalFitsRejected++;      
	  fIndexLocEq.Reset();  fNIndexLocEq=0; // reset stores and go to the next track 
	  fDerivLocEq.Reset();  fNDerivLocEq=0;	  
	  return 0;
	}

	lChi2 += lWeight*lMeas*lMeas ; // total chi^2
	nEq++;                    // number of equations			
	iMeas   = -1;
	iWeight = 0;
	iEqTerm--;
      } // End of "end of equation" operations
    }   // End of loop on equation
    iEqTerm++;
  } // End of loop on all equations used in the fit

  nDoF = nEq-nRank;	
  lRedChi2 = 0.0;

  AliDebug(1,Form("Final chi square / degrees of freedom %.2f / %d",lChi2, nDoF));
  
  if (nDoF > 0) lRedChi2 = lChi2/float(nDoF);  // Chi^2/dof
	
  fNLocalFits++;

  if (fNStdDev != 0 && nDoF > 0 && !bSingleFit) // Chisquare cut
  {
    lChi2Cut = AliMillepede::Chi2DoFLim(fNStdDev, nDoF)*fChi2CutFactor;
    
    AliDebug(1,Form("Reject if Chisq/Ndf = %.4f > %.4f",lRedChi2,lChi2Cut));
 
    if (lRedChi2 > lChi2Cut) // Reject the track if too much...
    {
      AliDebug(2,"Rejected track !!!!!");
      fNLocalFitsRejected++;      
      fIndexLocEq.Reset();  fNIndexLocEq=0; // reset stores and go to the next track 
      fDerivLocEq.Reset();  fNDerivLocEq=0;
      return 0;
    }
  }

  if (bSingleFit) // Stop here if just updating the track parameters
  {
    fIndexLocEq.Reset();  fNIndexLocEq=0; // Reset store for the next track 
    fDerivLocEq.Reset();  fNDerivLocEq=0;
    return 1;
  }

//  
// THIRD LOOP: local operations are finished, track is accepted 
// We now update the global parameters (other matrices)
//

  iEqTerm = 0;
  iMeas = -1;
  iWeight = 0;

  while (iEqTerm <= nEqTerms)
  {
    if (fIndexLocEq[iEqTerm] == -1)
    {
      if (iMeas == -1) {        // First  -1 : lMeas
	iMeas = iEqTerm;
	iLocFirst = iEqTerm+1;
      }  
      else if (iWeight == 0) {  // Second -1 : weight 
	iWeight = iEqTerm;
	iLocLast = iEqTerm-1;
	iGloFirst = iEqTerm+1;
      } 
      else {                    // Third  -1 : end of equation; start of next  
	iGloLast = iEqTerm-1;
	lMeas	= fDerivLocEq[iMeas];
	lWeight 	= fDerivLocEq[iWeight];

        // Now suppress the global part
	for (i=iGloFirst; i<=iGloLast; i++) {
	  iIdx = fIndexLocEq[i];   // Global param indice
	  if (fIsNonLinear[iIdx] == 0)
	    lMeas -= fDerivLocEq[i]*(fInitPar[iIdx]+fDeltaPar[iIdx]); // linear parameter
	  else
	    lMeas -= fDerivLocEq[i]*(fDeltaPar[iIdx]); // nonlinear parameter
	}
        
	for (i=iGloFirst; i<=iGloLast; i++) {
	  iIdx = fIndexLocEq[i];   // Global param indice (the matrix line)          
        
	  fVecBGlo[iIdx] += lWeight*lMeas*fDerivLocEq[i];  
// 	  AliDebug(2,Form("fVecBGlo[%d] = %.6f", j, fVecBGlo[j] ));

	  // First of all, the global/global terms (exactly like local matrix)
	  //	  
	  for (j=iGloFirst; j<=iGloLast; j++) {	  
	    jIdx = fIndexLocEq[j];			
            fMatCGlo[iIdx][jIdx] += lWeight*fDerivLocEq[i]*fDerivLocEq[j];
// 	    AliDebug(2,Form("fMatCGlo[%d][%d] = %.6f",iIdx,jIdx,fMatCGlo[iIdx][jIdx]));
	  } 

	  // Now we have also rectangular matrices containing global/local terms.
	  //
	  iIdxIdx = fGlo2CGLRow[iIdx];  // Index of index          
	  if (iIdxIdx == -1) {	  // New global variable	 
	    for (k=0; k<fNLocalPar; k++) {
	      fMatCGloLoc[nGloInFit][k] = 0.0;  // Initialize the row
	    }
	    fGlo2CGLRow[iIdx] = nGloInFit;
	    fCGLRow2Glo[nGloInFit] = iIdx;
	    iIdxIdx = nGloInFit;
	    nGloInFit++;
	  }

	  // Now fill the rectangular matrix
	  for (k=iLocFirst; k<=iLocLast ; k++) {
	    kIdx = fIndexLocEq[k];						
	    fMatCGloLoc[iIdxIdx][kIdx] += lWeight*fDerivLocEq[i]*fDerivLocEq[k];
// 	    AliDebug(2,Form("fMatCGloLoc[%d][%d] = %.6f",iIdxIdx,kIdx,fMatCGloLoc[iIdxIdx][kIdx]));
	  } 
	}
	iMeas   = -1;
	iWeight =  0;
	iEqTerm--;
      } // End of "end of equation" operations
    }   // End of loop on equation
    iEqTerm++;
  } // End of loop on all equations used in the fit
	
  // Third loop is finished, now we update the correction matrices
  AliMillepede::SpAVAt(fMatCLoc, fMatCGloLoc, fMatCGloCorr, fNLocalPar, nGloInFit);
  AliMillepede::SpAX(fMatCGloLoc, fVecBLoc, fVecBGloCorr, fNLocalPar, nGloInFit);

  for (iIdxIdx=0; iIdxIdx<nGloInFit; iIdxIdx++) {
    iIdx = fCGLRow2Glo[iIdxIdx];
    fVecBGlo[iIdx] -= fVecBGloCorr[iIdxIdx];
    
    for (jIdxIdx=0; jIdxIdx<=iIdxIdx; jIdxIdx++) {    
      int jIdx = fCGLRow2Glo[jIdxIdx];
      fMatCGlo[iIdx][jIdx] -= fMatCGloCorr[iIdxIdx][jIdxIdx];
      fMatCGlo[jIdx][iIdx] = fMatCGlo[iIdx][jIdx];
    }
  }
	
  fIndexLocEq.Reset();  fNIndexLocEq=0; // Reset store for the next track 
  fDerivLocEq.Reset();  fNDerivLocEq=0;

  return 1;
}
 

/*
-----------------------------------------------------------
  MAKEGLOBALFIT:  perform global params fit, once all the 'tracks'
                  have been fitted by FitLoc
-----------------------------------------------------------

  par[]        = array containing the computed global 
                 parameters (the misalignment constants)

  error[]      = array containing the error on global 
                 parameters (estimated by AliMillepede)

  pull[]        = array containing the corresponding pulls 

-----------------------------------------------------------
*/
Int_t AliMillepede::GlobalFit(double par[], double error[], double pull[])
{
  /// perform global parameters fit once all the local equations have been fitted
  int i, j;
  int nVar    = 0;
  int nGloFix = 0;
  double lConstraint;

  double step[1010];

  double localPars[2*fgkMaxLocalPar];

  int nLocFitsGood = 0;
  int nLocFitsTot  = 0;
  int nLocFits     = 0;

  AliInfo("..... Making global fit .....");

  nLocFitsTot = AliMillepede::GetNLocalEquations();
	
  while (fIter < fMaxIter)  // Iteration for the final loop
  {

    nLocFits = AliMillepede::GetNLocalEquations();
    AliInfo(Form("...using %d local fits...",nLocFits));

// Start by saving the diagonal elements
    
    for (i=0; i<fNGlobalPar; i++) {
      fDiagCGlo[i] = fMatCGlo[i][i];
    }

//  Then we retrieve the different constraints: fixed parameter or global equation

    nGloFix = 0; // First look at the fixed global params
    
    for (i=0; i<fNGlobalPar; i++) {    
      if (fSigmaPar[i] <= 0.0) {  // fixed global param
	nGloFix++;
	for (j=0; j<fNGlobalPar; j++) {
	  fMatCGlo[i][j] = 0.0;  // Reset row and column
	  fMatCGlo[j][i] = 0.0;
	}
      }
      else {
	fMatCGlo[i][i] += 1.0/(fSigmaPar[i]*fSigmaPar[i]);
      }
    }
        
    nVar = fNGlobalPar;  // Current number of equations	
    AliDebug(1,Form("Number of constraint equations : %d", fNGlobalConstraints));
    
    for (i=0; i<fNGlobalConstraints; i++) { // Then the constraint equation    
      lConstraint = fLagMult[i];
      for (j=0; j<fNGlobalPar; j++) {	
	fMatCGlo[nVar][j] = float(nLocFits)*fMatDerConstr[i][j];
	fMatCGlo[j][nVar] = float(nLocFits)*fMatDerConstr[i][j];          
	lConstraint -= fMatDerConstr[i][j]*(fInitPar[j]+fDeltaPar[j]);
      }
	
      fMatCGlo[nVar][nVar] = 0.0;
      fVecBGlo[nVar] = float(nLocFits)*lConstraint;
      nVar++;
    }


    // Intended to compute the final global chisquare

    double lFinalCor = 0.0;

    if (fIter > 1) {    
      for (i=0; i<fNGlobalPar; i++) {	
	for (j=0; j<fNGlobalPar; j++) {
// 	  printf("%d, %d, %.6f  %.6f  %.6f\n",i,j,step[i],fMatCGlo[i][j],step[j]);
	  lFinalCor += step[i]*fMatCGlo[i][j]*step[j]; 
	  if (i == j && fSigmaPar[i] != 0) {
	    lFinalCor -= step[i]*step[i]/(fSigmaPar[i]*fSigmaPar[i]);
	  }
	}
      }
    }

    AliInfo(Form(" Final coeff is %.6f",lFinalCor));		
    AliInfo(Form(" Final NDOFs = %d", fNGlobalPar));

    //  The final matrix inversion

    Int_t nRank = AliMillepede::SpmInv(fMatCGlo, fVecBGlo, nVar);

    for (i=0; i<fNGlobalPar; i++) {    
      fDeltaPar[i] += fVecBGlo[i];    // Update global parameters values (for iterations)
      AliDebug(1,Form("fDeltaPar[%d] = %.6f", i, fDeltaPar[i]));
      AliDebug(1,Form("fMatCGlo[%d][%d] = %.6f", i, i, fMatCGlo[i][i]));
      AliDebug(1,Form("err = %.6f", TMath::Sqrt(TMath::Abs(fMatCGlo[i][i]))));

      step[i] = fVecBGlo[i];

      if (fIter == 1) error[i] = fMatCGlo[i][i]; // Unfitted error
    }
    AliInfo("");
    AliInfo(Form("The rank defect of the symmetric %d by %d matrix is %d (bad if non 0)",
	    nVar, nVar, nVar-nGloFix-nRank));
		
    AliInfo("");
    AliInfo(Form("Total : %d local fits, %d rejected.", fNLocalFits, fNLocalFitsRejected));
    if (fIter == 0)  break;  // No iterations set     
    fIter++;
        
    if (fIter == fMaxIter)  break;  // End of story         

    // Reinitialize parameters for iteration
    //
    fNLocalFits = 0;
    fNLocalFitsRejected = 0;

    if (fChi2CutFactor != fChi2CutRef) {    
      fChi2CutFactor = TMath::Sqrt(fChi2CutFactor);
      if (fChi2CutFactor < 1.2*fChi2CutRef) {
	fChi2CutFactor = fChi2CutRef;
	fIter = fMaxIter - 1;     // Last iteration
      }
    }
    AliInfo(Form("Iteration %d with cut factor %.2f", fIter, fChi2CutFactor));
    
    // Reset global variables
    //    
    for (i=0; i<nVar; i++) {    
      fVecBGlo[i] = 0.0;
      for (j=0; j<nVar; j++) {     
	fMatCGlo[i][j] = 0.0;
      }
    }

    //
    // We start a new iteration
    //

    // First we read the stores for retrieving the local params
    //
    nLocFitsGood = 0;

    for (i=0; i<nLocFitsTot; i++) {
      int iEqFirst = 0;
      int iEqLast = 0;

      (i>0) ? iEqFirst = fLocEqPlace[i-1] : iEqFirst = 0;
      iEqLast = fLocEqPlace[i];

      AliDebug(2,Form("Track %d : ",i));
      AliDebug(2,Form("Starts at %d", iEqFirst));
      AliDebug(2,Form("Ends at %d",iEqLast));

      if (fIndexAllEqs[iEqFirst] != -999) { // Fit is still OK      
	fIndexLocEq.Reset();  fNIndexLocEq=0;
	fDerivLocEq.Reset();  fNDerivLocEq=0;
	
	for (j=iEqFirst; j<iEqLast; j++) {
	  if (fNIndexLocEq==fIndexLocEq.GetSize()) fIndexLocEq.Set(2*fNIndexLocEq);
	  fIndexLocEq.AddAt(fIndexAllEqs[j],fNIndexLocEq++);    
	  if (fNDerivLocEq==fDerivLocEq.GetSize()) fDerivLocEq.Set(2*fNDerivLocEq);
	  fDerivLocEq.AddAt(fDerivAllEqs[j],fNDerivLocEq++);    
	}

	for (j=0; j<2*fNLocalPar; j++) {
	  localPars[j] = 0.;
	}	

//  	Int_t sc = AliMillepede::LocalFit(i,localPars,0);
// 	(sc) ? nLocFitsGood++ : fIndexAllEqs[iEqFirst] = -999;   
	AliMillepede::LocalFit(i,localPars,0);
	nLocFitsGood++;
      }
    } // End of loop on fits

    AliMillepede::SetNLocalEquations(nLocFitsGood);

  } // End of iteration loop
	
  AliMillepede::PrintGlobalParameters(); // Print the final results

  for (j=0; j<fNGlobalPar; j++) {  
    par[j]   = fInitPar[j]+fDeltaPar[j];
    pull[j]  = (fSigmaPar[j] <= 0.) ? 0. : fDeltaPar[j]/TMath::Sqrt(fSigmaPar[j]*fSigmaPar[j]-fMatCGlo[j][j]);
    error[j] = TMath::Sqrt(TMath::Abs(fMatCGlo[j][j]));
  }

  AliInfo(" ");
  AliInfo("            * o o                   o      ");
  AliInfo("              o o                   o      ");
  AliInfo("   o ooooo  o o o  oo  ooo   oo   ooo  oo  ");
  AliInfo("    o  o  o o o o o  o o  o o  o o  o o  o ");
  AliInfo("    o  o  o o o o oooo o  o oooo o  o oooo ");
  AliInfo("    o  o  o o o o o    ooo  o    o  o o    ");
  AliInfo("    o  o  o o o o  oo  o     oo   ooo  oo ++ ends.");
  AliInfo("                       o                   ");	  

  return 1;
}

/*
-----------------------------------------------------------
  ERRPAR: return error for parameter iPar
-----------------------------------------------------------

  iPar     = the index of the global parameter in the 
             result array (equivalent to fDeltaPar[]).
 
-----------------------------------------------------------
*/ 
Double_t AliMillepede::GetParError(Int_t iPar) const
{
  /// return error for parameter iPar
  Double_t lErr = -1.;
  if (iPar>=0 && iPar<fNGlobalPar) {
    lErr = TMath::Sqrt(TMath::Abs(fMatCGlo[iPar][iPar]));
  }
  return lErr;
}


/*
-----------------------------------------------------------
  SPMINV:  obtain solution of a system of linear equations with symmetric matrix 
 	   and the inverse (using 'singular-value friendly' GAUSS pivot)
-----------------------------------------------------------

	Solve the equation :  V * X = B
	
	V is replaced by inverse matrix and B by X, the solution vector
-----------------------------------------------------------
*/
int AliMillepede::SpmInv(double matV[][fgkMaxGloPC], double vecB[], int nGlo)
{
  ///  Obtain solution of a system of linear equations with symmetric matrix 
  ///  and the inverse (using 'singular-value friendly' GAUSS pivot)
    
  Int_t nRank = 0;
  int iPivot;
  double vPivot = 0.;
  double eps = 0.00000000000001;

  bool *bUnUsed = new bool[nGlo];
  double *diagV = new double[nGlo];
  double *rowMax = new double[nGlo];
  double *colMax = new double[nGlo];

  double *temp = new double[nGlo];

  for (Int_t i=0; i<nGlo; i++) {
    rowMax[i] = 0.0;
    colMax[i] = 0.0;
    bUnUsed[i] = true;

    for (Int_t j=0; j<i; j++) {
      if (matV[j][i] == 0) {
	matV[j][i] = matV[i][j];
      }
    }
  }
  
  // Small loop for matrix equilibration (gives a better conditioning) 

  for (Int_t i=0; i<nGlo; i++) {
    for (Int_t j=0; j<nGlo; j++) { 
      if (TMath::Abs(matV[i][j]) >= rowMax[i]) rowMax[i] = TMath::Abs(matV[i][j]); // Max elemt of row i
      if (TMath::Abs(matV[j][i]) >= colMax[i]) colMax[i] = TMath::Abs(matV[j][i]); // Max elemt of column i
    }
  }

  for (Int_t i=0; i<nGlo; i++) {
    if (0.0 != rowMax[i]) rowMax[i] = 1./rowMax[i]; // Max elemt of row i
    if (0.0 != colMax[i]) colMax[i] = 1./colMax[i]; // Max elemt of column i
  }

  for (Int_t i=0; i<nGlo; i++) {
    for (Int_t j=0; j<nGlo; j++) {
      matV[i][j] = TMath::Sqrt(rowMax[i])*matV[i][j]*TMath::Sqrt(colMax[j]); // Equilibrate the V matrix
    }
    diagV[i] = TMath::Abs(matV[i][i]); // save diagonal elem absolute values 	
  }


  for (Int_t i=0; i<nGlo; i++) {
    vPivot = 0.0;
    iPivot = -1;
    
    for (Int_t j=0; j<nGlo; j++) { // First look for the pivot, ie max unused diagonal element       
      if (bUnUsed[j] && (TMath::Abs(matV[j][j])>std::max(TMath::Abs(vPivot),eps*diagV[j]))) {    
	vPivot = matV[j][j];
	iPivot = j;
      }
    }
    
    if (iPivot >= 0) {   // pivot found          
      nRank++;
      bUnUsed[iPivot] = false; // This value is used
      vPivot = 1.0/vPivot;
      matV[iPivot][iPivot] = -vPivot; // Replace pivot by its inverse
      
      for (Int_t j=0; j<nGlo; j++) {      
	for (Int_t jj=0; jj<nGlo; jj++) {  
	  if (j != iPivot && jj != iPivot) {// Other elements (!!! do them first as you use old matV[k][j]'s !!!)	  
	    matV[j][jj] = matV[j][jj] - vPivot*matV[j][iPivot]*matV[iPivot][jj];
	  }
	}
      }
      
      for (Int_t j=0; j<nGlo; j++) {      
	if (j != iPivot) { // Pivot row or column elements 
	  matV[j][iPivot] = matV[j][iPivot]*vPivot;	// Column
	  matV[iPivot][j] = matV[iPivot][j]*vPivot;	// Line
	}
      }
    }
    else {  // No more pivot value (clear those elements)
      for (Int_t j=0; j<nGlo; j++) {
	if (bUnUsed[j]) {
	  vecB[j] = 0.0;

	  for (Int_t k=0; k<nGlo; k++) {
	    matV[j][k] = 0.0;
	    matV[k][j] = 0.0;
	  }
	}
      }
      break;  // No more pivots anyway, stop here
    }
  }
  
  for (Int_t i=0; i<nGlo; i++) {
    for (Int_t j=0; j<nGlo; j++) {
      matV[i][j] = TMath::Sqrt(colMax[i])*matV[i][j]*TMath::Sqrt(rowMax[j]); // Correct matrix V
    }
  }
  
  for (Int_t j=0; j<nGlo; j++) {
    temp[j] = 0.0;
    
    for (Int_t jj=0; jj<nGlo; jj++) { // Reverse matrix elements
      matV[j][jj] = -matV[j][jj];
      temp[j] += matV[j][jj]*vecB[jj];
    }		
  }

  for (Int_t j=0; j<nGlo; j++) {
    vecB[j] = temp[j]; // The final result
  }
  
  delete [] temp;
  delete [] bUnUsed;
  delete [] diagV;
  delete [] rowMax;
  delete [] colMax;

  return nRank;
}

//
// Same method but for local fit, so heavily simplified
//
int AliMillepede::SpmInv(double matV[][fgkMaxLocalPar], double vecB[], int nLoc)
{
  ///  Obtain solution of a system of linear equations with symmetric matrix 
  ///  and the inverse (using 'singular-value friendly' GAUSS pivot)

  Int_t nRank = 0;
  Int_t iPivot = -1;
  double vPivot = 0.;
  double eps = 0.0000000000001;

  bool *bUnUsed = new bool[nLoc];
  double *diagV  = new double[nLoc];
  double *temp  = new double[nLoc];
  
  for (Int_t i=0; i<nLoc; i++) {
    bUnUsed[i] = true;
    diagV[i] = TMath::Abs(matV[i][i]);     // save diagonal elem absolute values
    for (Int_t j=0; j<i; j++) {
      matV[j][i] = matV[i][j] ;
    }
  }
	
  for (Int_t i=0; i<nLoc; i++) {
    vPivot = 0.0;
    iPivot = -1;
		
    for (Int_t j=0; j<nLoc; j++) { // First look for the pivot, ie max unused diagonal element 
      if (bUnUsed[j] && (TMath::Abs(matV[j][j])>TMath::Max(TMath::Abs(vPivot),eps*diagV[j]))) {
	vPivot = matV[j][j];
	iPivot = j;
      }
    }
		
    if (iPivot >= 0) {   // pivot found
      nRank++;
      bUnUsed[iPivot] = false;
      vPivot = 1.0/vPivot;
      matV[iPivot][iPivot] = -vPivot; // Replace pivot by its inverse
      
      for (Int_t j=0; j<nLoc; j++) {
	if (j != iPivot) {
	  for (Int_t jj=0; jj<=j; jj++) {	
	    if (jj != iPivot) {// Other elements (!!! do them first as you use old matV[k][j]'s !!!)
	      matV[j][jj] = matV[j][jj] - vPivot*matV[j][iPivot]*matV[iPivot][jj];
	      matV[jj][j] = matV[j][jj];
	    }
	  }
	}
      }

      for (Int_t j=0; j<nLoc; j++) {      
	if (j != iPivot) {      // Pivot row or column elements 	
	  matV[j][iPivot] = matV[j][iPivot]*vPivot; // Column
	  matV[iPivot][j] = matV[j][iPivot];
	}
      }
    }
    else { // No more pivot value (clear those elements)
      for (Int_t j=0; j<nLoc; j++) {
	if (bUnUsed[j]) {
	  vecB[j] = 0.0;
	  matV[j][j] = 0.0;
	  for (Int_t k=0; k<j; k++) {	  
	    matV[j][k] = 0.0;
	    matV[k][j] = 0.0;
	  }
	}
      }

      break;  // No more pivots anyway, stop here
    }
  }

  for (Int_t j=0; j<nLoc; j++) {  
    temp[j] = 0.0;    
    for (Int_t jj=0; jj<nLoc; jj++) { // Reverse matrix elements
      matV[j][jj] = -matV[j][jj];
      temp[j] += matV[j][jj]*vecB[jj];
    }			
  }

  for (Int_t j=0; j<nLoc; j++) {
    vecB[j] = temp[j];
  }

  delete [] bUnUsed;
  delete [] diagV;
  delete [] temp;
  
  return nRank;
}


/*
-----------------------------------------------------------
  SPAVAT
-----------------------------------------------------------

  multiply symmetric N-by-N matrix from the left with general M-by-N
  matrix and from the right with the transposed of the same  general
  matrix  to  form  symmetric  M-by-M   matrix.
  
                                                       T
  CALL SPAVAT(V,A,W,N,M)      W   =   A   *   V   *   A
   		           M*M     M*N     N*N     N*M
  
  where V = symmetric N-by-N matrix
        A = general N-by-M matrix
        W = symmetric M-by-M matrix
-----------------------------------------------------------
*/
Int_t AliMillepede::SpAVAt(double matV[][fgkMaxLocalPar], double matA[][fgkMaxLocalPar], double matW[][fgkMaxGlobalPar], int nLoc, int nGlo)
{
  ///  multiply symmetric N-by-N matrix from the left with general M-by-N
  ///  matrix and from the right with the transposed of the same general
  ///  matrix to form symmetric M-by-M matrix.

  for (Int_t i=0; i<nGlo; i++) {
    for (Int_t j=0; j<=i; j++) {  // Matrix w is symmetric
      matW[i][j] = 0.0;      // Reset final matrix			
      for (Int_t k=0; k<nLoc; k++) {	
	matW[i][j] += matA[i][k]*matV[k][k]*matA[j][k];    // diagonal terms of v
	for (Int_t l=0; l<k; l++) {	  
	  matW[i][j] += matA[i][k]*matV[k][l]*matA[j][l];  // Use symmetric properties
	  matW[i][j] += matA[i][l]*matV[k][l]*matA[j][k];  // of matrix v
	}
      }
      if (i!=j){
	matW[j][i] = matW[i][j]; // Matrix w is symmetric
      }
    }
  }
	
  return 1;
}


/*
-----------------------------------------------------------
  SPAX
-----------------------------------------------------------

  multiply general M-by-N matrix A and N-vector X
 
  CALL  SPAX(A,X,Y,M,N)          Y =  A * X
                                 M   M*N  N
 
  where A = general M-by-N matrix (A11 A12 ... A1N  A21 A22 ...)
        X = N vector
        Y = M vector
-----------------------------------------------------------
*/
Int_t AliMillepede::SpAX(double matA[][fgkMaxLocalPar], double vecX[], double vecY[], int nCol, int nRow)
{
  ///   multiply general M-by-N matrix A and N-vector X
  for (Int_t i=0; i<nRow; i++) {
    vecY[i] = 0.0;	    // Reset final vector			
    for (Int_t j=0; j<nCol; j++) {
      vecY[i] += matA[i][j]*vecX[j];  // fill the vector
    }
  }
	
  return 1;
}


/*
-----------------------------------------------------------
  PRTGLO
-----------------------------------------------------------

  Print the final results into the logfile

-----------------------------------------------------------
*/
Int_t AliMillepede::PrintGlobalParameters() const
{
  ///  Print the final results into the logfile
  double lError = 0.;
  double lGlobalCor =0.;
	
  AliInfo("");
  AliInfo("   Result of fit for global parameters");
  AliInfo("   ===================================");
  AliInfo("    I       initial       final       differ        lastcor        error       gcor");
  AliInfo("-----------------------------------------------------------------------------------");
	
  for (int i=0; i<fNGlobalPar; i++) {
    lError = TMath::Sqrt(TMath::Abs(fMatCGlo[i][i]));
    if (fMatCGlo[i][i] < 0.0) lError = -lError;
    lGlobalCor = 0.0;
		
    if (TMath::Abs(fMatCGlo[i][i]*fDiagCGlo[i]) > 0) {    
      lGlobalCor = TMath::Sqrt(TMath::Abs(1.0-1.0/(fMatCGlo[i][i]*fDiagCGlo[i])));
      AliInfo(Form("%d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f",
		   i,fInitPar[i],fInitPar[i]+fDeltaPar[i],fDeltaPar[i],fVecBGlo[i],lError,lGlobalCor));
    }
    else {    
      AliInfo(Form("%d\t %.6f\t %.6f\t %.6f\t %.6f\t OFF\t OFF",i,fInitPar[i],fInitPar[i]+fDeltaPar[i],fDeltaPar[i],fVecBGlo[i]));
    }
  }
  return 1;
}


/*
----------------------------------------------------------------
  CHI2DOFLIM:  return the limit in chi^2/nd for n sigmas stdev authorized
----------------------------------------------------------------

  Only n=1, 2, and 3 are expected in input
----------------------------------------------------------------
*/
double AliMillepede::Chi2DoFLim(int nSig, int nDoF)
{
  /// return the limit in chi^2/nd for n sigmas stdev authorized
  int lNSig;
  double sn[3]        =	{0.47523, 1.690140, 2.782170};
  double table[3][30] = {{1.0000, 1.1479, 1.1753, 1.1798, 1.1775, 1.1730, 1.1680, 1.1630,
			  1.1581, 1.1536, 1.1493, 1.1454, 1.1417, 1.1383, 1.1351, 1.1321,
			  1.1293, 1.1266, 1.1242, 1.1218, 1.1196, 1.1175, 1.1155, 1.1136,
			  1.1119, 1.1101, 1.1085, 1.1070, 1.1055, 1.1040},
			 {4.0000, 3.0900, 2.6750, 2.4290, 2.2628, 2.1415, 2.0481, 1.9736,
			  1.9124, 1.8610, 1.8171, 1.7791, 1.7457, 1.7161, 1.6897, 1.6658,
			  1.6442, 1.6246, 1.6065, 1.5899, 1.5745, 1.5603, 1.5470, 1.5346,
			  1.5230, 1.5120, 1.5017, 1.4920, 1.4829, 1.4742},
			 {9.0000, 5.9146, 4.7184, 4.0628, 3.6410, 3.3436, 3.1209, 2.9468,
			  2.8063, 2.6902, 2.5922, 2.5082, 2.4352, 2.3711, 2.3143, 2.2635,
			  2.2178, 2.1764, 2.1386, 2.1040, 2.0722, 2.0428, 2.0155, 1.9901,
			  1.9665, 1.9443, 1.9235, 1.9040, 1.8855, 1.8681}};

  if (nDoF < 1) {
    return 0.0;
  }
  else {  
    lNSig = TMath::Max(1,TMath::Min(nSig,3));

    if (nDoF <= 30) {    
      return table[lNSig-1][nDoF-1];
    }
    else { // approximation
      return ((sn[lNSig-1]+TMath::Sqrt(float(2*nDoF-3)))*
	      (sn[lNSig-1]+TMath::Sqrt(float(2*nDoF-3))))/float(2*nDoF-2);
    }
  }
}
