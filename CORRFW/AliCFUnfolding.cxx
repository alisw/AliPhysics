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

//---------------------------------------------------------------------//
//                                                                     //
// AliCFUnfolding Class                                                //
// Class to handle general unfolding procedure                         // 
// For the moment only bayesian unfolding is supported                 //
// The next steps are to add chi2 minimisation and weighting methods   //
//                                                                     //
//                                                                     //
//                                                                     //
// Use :                                                               //
// -------                                                             //
// The Bayesian unfolding consists of several iterations.              //
// At each iteration, an inverse response matrix is calculated, given  //
// the measured spectrum, the a priori (guessed) spectrum,             //
// the efficiency spectrum and the response matrix.                    //
//                                                                     //
// Then at each iteration, the unfolded spectrum is calculated using   //
// the inverse response : the goal is to get an unfolded spectrum      //
// similar (according to some criterion) to the a priori one.          //
// If the difference is too big, another iteration is performed :      //
// the a priori spectrum is updated to the unfolded one from the       //
// previous iteration, and so on so forth, until the maximum number    //
// of iterations or the similarity criterion is reached.               //
//                                                                     //
// Chi2 calculation became absolute with the correlated error          //
// calculation.                                                        //
// Errors on the unfolded distribution are not known until the end     //
// Use the convergence criterion instead                               //
//                                                                     //
// Currently the user has to define the max. number of iterations      //
// (::SetMaxNumberOfIterations)                                        //
// and                                                                 //
//     - the chi2 below which the procedure will stop                  //
// (::SetMaxChi2 or ::SetMaxChi2PerDOF)   (OBSOLETE)                   //
//     - the convergence criterion below which the procedure will stop //
// SetMaxConvergencePerDOF(Double_t val);                              //
//                                                                     //
// Correlated error calculation can be activated by using:             //
// SetUseCorrelatedErrors(Bool_t b) in combination with convergence    //
// criterion                                                           //
// Documentation about correlated error calculation method can be      //
// found in AliCFUnfolding::CalculateCorrelatedErrors()                //
// Author: marta.verweij@cern.ch                                       //
//                                                                     //
// An optional possibility is to smooth the unfolded spectrum at the   //
// end of each iteration, either using a fit function                  //
// (only if #dimensions <=3)                                           //
// or a simple averaging using the neighbouring bins values.           //
// This is possible calling the function ::UseSmoothing                //
// If no argument is passed to this function, then the second option   //
// is used.                                                            //
//                                                                     //
// IMPORTANT:                                                          //
//-----------                                                          //
// With this approach, the efficiency map must be calculated           //
// with *simulated* values only, otherwise the method won't work.      //
//                                                                     //
// ex: efficiency(bin_pt) = number_rec(bin_pt) / number_sim(bin_pt)    //
//                                                                     //
// the pt bin "bin_pt" must always be the same in both the efficiency  //
// numerator and denominator.                                          //
// This is why the efficiency map has to be created by a method        //
// from which both reconstructed and simulated values are accessible   //
// simultaneously.                                                     //
//                                                                     //
//                                                                     //
//---------------------------------------------------------------------//
// Author : renaud.vernet@cern.ch                                      //
//---------------------------------------------------------------------//


#include "AliCFUnfolding.h"
#include "TMath.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TRandom3.h"


ClassImp(AliCFUnfolding)

//______________________________________________________________

AliCFUnfolding::AliCFUnfolding() :
  TNamed(),
  fResponseOrig(0x0),
  fPriorOrig(0x0),
  fEfficiencyOrig(0x0),
  fMeasuredOrig(0x0),
  fMaxNumIterations(0),
  fNVariables(0),
  fUseSmoothing(kFALSE),
  fSmoothFunction(0x0),
  fSmoothOption("iremn"),
  fMaxConvergence(0),
  fNRandomIterations(0),
  fResponse(0x0),
  fPrior(0x0),
  fEfficiency(0x0),
  fMeasured(0x0),
  fInverseResponse(0x0),
  fMeasuredEstimate(0x0),
  fConditional(0x0),
  fUnfolded(0x0),
  fUnfoldedFinal(0x0),
  fCoordinates2N(0x0),
  fCoordinatesN_M(0x0),
  fCoordinatesN_T(0x0),
  fRandomResponse(0x0),
  fRandomEfficiency(0x0),
  fRandomMeasured(0x0),
  fRandom3(0x0),
  fDeltaUnfoldedP(0x0),
  fDeltaUnfoldedN(0x0),
  fNCalcCorrErrors(0),
  fRandomSeed(0)
{
  //
  // default constructor
  //
}

//______________________________________________________________

AliCFUnfolding::AliCFUnfolding(const Char_t* name, const Char_t* title, const Int_t nVar, 
			       const THnSparse* response, const THnSparse* efficiency, const THnSparse* measured, const THnSparse* prior ,
			       Double_t maxConvergencePerDOF, UInt_t randomSeed, Int_t maxNumIterations
			       ) :
  TNamed(name,title),
  fResponseOrig((THnSparse*)response->Clone()),
  fPriorOrig(0x0),
  fEfficiencyOrig((THnSparse*)efficiency->Clone()),
  fMeasuredOrig((THnSparse*)measured->Clone()),
  fMaxNumIterations(maxNumIterations),
  fNVariables(nVar),
  fUseSmoothing(kFALSE),
  fSmoothFunction(0x0),
  fSmoothOption("iremn"),
  fMaxConvergence(0),
  fNRandomIterations(maxNumIterations),
  fResponse((THnSparse*)response->Clone()),
  fPrior(0x0),
  fEfficiency((THnSparse*)efficiency->Clone()),
  fMeasured((THnSparse*)measured->Clone()),
  fInverseResponse(0x0),
  fMeasuredEstimate(0x0),
  fConditional(0x0),
  fUnfolded(0x0),
  fUnfoldedFinal(0x0),
  fCoordinates2N(0x0),
  fCoordinatesN_M(0x0),
  fCoordinatesN_T(0x0),
  fRandomResponse((THnSparse*)response->Clone()),
  fRandomEfficiency((THnSparse*)efficiency->Clone()),
  fRandomMeasured((THnSparse*)measured->Clone()),
  fRandom3(0x0),
  fDeltaUnfoldedP(0x0),
  fDeltaUnfoldedN(0x0),
  fNCalcCorrErrors(0),
  fRandomSeed(randomSeed)
{
  //
  // named constructor
  //

  AliInfo(Form("\n\n--------------------------\nCreating an unfolder :\n--------------------------\nresponse matrix has %d dimension(s)",fResponse->GetNdimensions()));

  if (!prior) CreateFlatPrior(); // if no prior distribution declared, simply use a flat distribution
  else {
    fPrior = (THnSparse*) prior->Clone();
    fPriorOrig = (THnSparse*)fPrior->Clone();
    if (fPrior->GetNdimensions() != fNVariables) 
      AliFatal(Form("The prior matrix should have %d dimensions, and it has actually %d",fNVariables,fPrior->GetNdimensions()));
  }

  if (fEfficiency->GetNdimensions() != fNVariables) 
    AliFatal(Form("The efficiency matrix should have %d dimensions, and it has actually %d",fNVariables,fEfficiency->GetNdimensions()));
  if (fMeasured->GetNdimensions() != fNVariables) 
    AliFatal(Form("The measured matrix should have %d dimensions, and it has actually %d",fNVariables,fMeasured->GetNdimensions()));
  if (fResponse->GetNdimensions() != 2*fNVariables) 
    AliFatal(Form("The response matrix should have %d dimensions, and it has actually %d",2*fNVariables,fResponse->GetNdimensions()));
  

  for (Int_t iVar=0; iVar<fNVariables; iVar++) {
    AliInfo(Form("prior      matrix has %d bins in dimension %d",fPrior     ->GetAxis(iVar)->GetNbins(),iVar));
    AliInfo(Form("efficiency matrix has %d bins in dimension %d",fEfficiency->GetAxis(iVar)->GetNbins(),iVar));
    AliInfo(Form("measured   matrix has %d bins in dimension %d",fMeasured  ->GetAxis(iVar)->GetNbins(),iVar));
  }

  fRandomResponse  ->SetTitle("Randomized response matrix");
  fRandomEfficiency->SetTitle("Randomized efficiency");
  fRandomMeasured  ->SetTitle("Randomized measured");
  SetMaxConvergencePerDOF(maxConvergencePerDOF)  ;
  Init();
}

//______________________________________________________________

AliCFUnfolding::~AliCFUnfolding() {
  //
  // destructor
  //
  if (fResponse)           delete fResponse;
  if (fPrior)              delete fPrior;
  if (fEfficiency)         delete fEfficiency;
  if (fEfficiencyOrig)     delete fEfficiencyOrig;
  if (fMeasured)           delete fMeasured;
  if (fMeasuredOrig)       delete fMeasuredOrig;
  if (fSmoothFunction)     delete fSmoothFunction;
  if (fPriorOrig)          delete fPriorOrig;
  if (fInverseResponse)    delete fInverseResponse;
  if (fMeasuredEstimate)   delete fMeasuredEstimate;
  if (fConditional)        delete fConditional;
  if (fUnfolded)           delete fUnfolded;
  if (fUnfoldedFinal)      delete fUnfoldedFinal;
  if (fCoordinates2N)      delete [] fCoordinates2N; 
  if (fCoordinatesN_M)     delete [] fCoordinatesN_M; 
  if (fCoordinatesN_T)     delete [] fCoordinatesN_T; 
  if (fRandomResponse)     delete fRandomResponse;
  if (fRandomEfficiency)   delete fRandomEfficiency;
  if (fRandomMeasured)     delete fRandomMeasured;
  if (fRandom3)            delete fRandom3;
  if (fDeltaUnfoldedP)     delete fDeltaUnfoldedP;
  if (fDeltaUnfoldedN)     delete fDeltaUnfoldedN;
 
}

//______________________________________________________________

void AliCFUnfolding::Init() {
  //
  // initialisation function : creates internal settings
  //

  fRandom3 = new TRandom3(fRandomSeed);

  fCoordinates2N  = new Int_t[2*fNVariables];
  fCoordinatesN_M = new Int_t[fNVariables];
  fCoordinatesN_T = new Int_t[fNVariables];

  // create the matrix of conditional probabilities P(M|T)
  CreateConditional(); //done only once at initialization
  
  // create the frame of the inverse response matrix
  fInverseResponse  = (THnSparse*) fResponse->Clone();
  // create the frame of the unfolded spectrum
  fUnfolded = (THnSparse*) fPrior->Clone();
  fUnfolded->SetTitle("Unfolded");
  // create the frame of the measurement estimate spectrum
  fMeasuredEstimate = (THnSparse*) fMeasured->Clone();
  
  // create the frame of the delta profiles
  fDeltaUnfoldedP = (THnSparse*)fPrior->Clone();
  fDeltaUnfoldedP->SetTitle("#Delta unfolded");
  fDeltaUnfoldedP->Reset();
  fDeltaUnfoldedN = (THnSparse*)fPrior->Clone();
  fDeltaUnfoldedN->SetTitle("");
  fDeltaUnfoldedN->Reset();


}


//______________________________________________________________

void AliCFUnfolding::CreateEstMeasured() {
  //
  // This function creates a estimate (M) of the reconstructed spectrum 
  // given the a priori distribution (T), the efficiency (E) and the conditional matrix (COND)
  //
  // --> P(M) = SUM   { P(M|T)    * P(T) }
  // --> M(i) = SUM_k { COND(i,k) * T(k) * E (k)}
  //
  // This is needed to calculate the inverse response matrix
  //


  // clean the measured estimate spectrum
  fMeasuredEstimate->Reset();

  THnSparse* priorTimesEff = (THnSparse*) fPrior->Clone();
  priorTimesEff->Multiply(fEfficiency);

  // fill it
  for (Long_t iBin=0; iBin<fConditional->GetNbins(); iBin++) {
    Double_t conditionalValue = fConditional->GetBinContent(iBin,fCoordinates2N);
    GetCoordinates();
    Double_t priorTimesEffValue = priorTimesEff->GetBinContent(fCoordinatesN_T);
    Double_t fill = conditionalValue * priorTimesEffValue ;
    
    if (fill>0.) {
      fMeasuredEstimate->AddBinContent(fCoordinatesN_M,fill);
      fMeasuredEstimate->SetBinError(fCoordinatesN_M,0.);
    }
  }
  delete priorTimesEff ;
}

//______________________________________________________________

void AliCFUnfolding::CreateInvResponse() {
  //
  // Creates the inverse response matrix (INV) with Bayesian method
  //  : uses the conditional matrix (COND), the prior probabilities (T) and the efficiency map (E)
  //
  // --> P(T|M)   = P(M|T)    * P(T) * eff(T) / SUM   { P(M|T)    * P(T) }
  // --> INV(i,j) = COND(i,j) * T(j) * E(j)   / SUM_k { COND(i,k) * T(k) }
  //

  THnSparse* priorTimesEff = (THnSparse*) fPrior->Clone();
  priorTimesEff->Multiply(fEfficiency);

  for (Long_t iBin=0; iBin<fConditional->GetNbins(); iBin++) {
    Double_t conditionalValue = fConditional->GetBinContent(iBin,fCoordinates2N);
    GetCoordinates();
    Double_t estMeasuredValue   = fMeasuredEstimate->GetBinContent(fCoordinatesN_M);
    Double_t priorTimesEffValue = priorTimesEff    ->GetBinContent(fCoordinatesN_T);
    Double_t fill = (estMeasuredValue>0. ? conditionalValue * priorTimesEffValue / estMeasuredValue : 0. ) ;
    if (fill>0. || fInverseResponse->GetBinContent(fCoordinates2N)>0.) {
      fInverseResponse->SetBinContent(fCoordinates2N,fill);
      fInverseResponse->SetBinError  (fCoordinates2N,0.);
    }
  } 
  delete priorTimesEff ;
}

//______________________________________________________________

void AliCFUnfolding::Unfold() {
  //
  // Main routine called by the user : 
  // it calculates the unfolded spectrum from the response matrix, measured spectrum and efficiency
  // several iterations are performed until a reasonable chi2 or convergence criterion is reached
  //

  Int_t iIterBayes     = 0 ;
  Double_t convergence = 0.;

  for (iIterBayes=0; iIterBayes<fMaxNumIterations; iIterBayes++) { // bayes iterations

    CreateEstMeasured(); // create measured estimate from prior
    CreateInvResponse(); // create inverse response  from prior
    CreateUnfolded();    // create unfoled spectrum  from measured and inverse response

    convergence = GetConvergence();
    AliDebug(0,Form("convergence at iteration %d is %e",iIterBayes,convergence));

    if (fMaxConvergence>0. && convergence<fMaxConvergence && fNCalcCorrErrors == 0) {
      fNRandomIterations = iIterBayes;
      AliDebug(0,Form("convergence is met at iteration %d",iIterBayes));
      break;
    }

    if (fUseSmoothing) {
      if (Smooth()) {
	AliError("Couldn't smooth the unfolded spectrum!!");
	if (fNCalcCorrErrors>0) {
	  AliInfo(Form("=======================\nUnfold of randomized distribution finished at iteration %d with convergence %e \n",iIterBayes,convergence));
	}
	else {
	  AliInfo(Form("\n\n=======================\nFinish at iteration %d : convergence is %e and you required it to be < %e\n=======================\n\n",iIterBayes,convergence,fMaxConvergence));
	}
	return;
      }
    }

    // update the prior distribution
    if (fPrior) delete fPrior ;
    fPrior = (THnSparse*)fUnfolded->Clone() ;
    fPrior->SetTitle("Prior");

  } // end bayes iteration

  if (fNCalcCorrErrors==0) fUnfoldedFinal = (THnSparse*) fUnfolded->Clone() ;

  //
  //for (Long_t iBin=0; iBin<fUnfoldedFinal->GetNbins(); iBin++) AliDebug(2,Form("%e\n",fUnfoldedFinal->GetBinError(iBin)));
  //

  if (fNCalcCorrErrors == 0) {
    AliInfo("\n================================================\nFinished bayes iteration, now calculating errors...\n================================================\n");
    fNCalcCorrErrors = 1;
    CalculateCorrelatedErrors();
  }

  if (fNCalcCorrErrors >1 ) {
    AliInfo(Form("\n\n=======================\nFinished at iteration %d : convergence is %e and you required it to be < %e\n=======================\n\n",iIterBayes,convergence,fMaxConvergence));
  }
  else if(fNCalcCorrErrors>0) {
    AliInfo(Form("=======================\nUnfolding of randomized distribution finished at iteration %d with convergence %e \n",iIterBayes,convergence));
  }
}

//______________________________________________________________

void AliCFUnfolding::CreateUnfolded() {
  //
  // Creates the unfolded (T) spectrum from the measured spectrum (M) and the inverse response matrix (INV)
  // We have P(T) = SUM   { P(T|M)   * P(M) } 
  //   -->   T(i) = SUM_k { INV(i,k) * M(k) }
  //


  // clear the unfolded spectrum
  // if in the process of error calculation, the random unfolded spectrum is created
  // otherwise the normal unfolded spectrum is created

  fUnfolded->Reset();
  
  for (Long_t iBin=0; iBin<fInverseResponse->GetNbins(); iBin++) {
    Double_t invResponseValue = fInverseResponse->GetBinContent(iBin,fCoordinates2N);
    GetCoordinates();
    Double_t effValue      = fEfficiency->GetBinContent(fCoordinatesN_T);
    Double_t measuredValue = fMeasured  ->GetBinContent(fCoordinatesN_M);
    Double_t fill = (effValue>0. ? invResponseValue * measuredValue / effValue : 0.) ;

    if (fill>0.) {
      // set errors to zero
      // true errors will be filled afterwards
      Double_t err = 0.;
      fUnfolded->SetBinError  (fCoordinatesN_T,err);
      fUnfolded->AddBinContent(fCoordinatesN_T,fill);
    }
  }
}

//______________________________________________________________

void AliCFUnfolding::CalculateCorrelatedErrors() {

  // Step 1: Create randomized distribution (fRandomXXXX) of each bin of 
  //         the measured spectrum to calculate correlated errors. 
  //         Poisson statistics: mean = measured value of bin
  // Step 2: Unfold randomized distribution
  // Step 3: Store difference of unfolded spectrum from measured distribution and 
  //         unfolded distribution from randomized distribution 
  //         -> fDeltaUnfoldedP (TProfile with option "S")
  // Step 4: Repeat Step 1-3 several times (fNRandomIterations)
  // Step 5: The spread of fDeltaUnfoldedP for each bin is the error on the unfolded spectrum of that specific bin


  //Do fNRandomIterations = bayes iterations performed
  for (int i=0; i<fNRandomIterations; i++) {
    
    // reset prior to original one
    if (fPrior) delete fPrior ;
    fPrior = (THnSparse*) fPriorOrig->Clone();

    // create randomized distribution and stick measured spectrum to it
    CreateRandomizedDist();

    if (fResponse) delete fResponse ;
    fResponse = (THnSparse*) fRandomResponse->Clone();
    fResponse->SetTitle("Response");

    if (fEfficiency) delete fEfficiency ;
    fEfficiency = (THnSparse*) fRandomEfficiency->Clone();
    fEfficiency->SetTitle("Efficiency");

    if (fMeasured)   delete fMeasured   ;
    fMeasured = (THnSparse*) fRandomMeasured->Clone();
    fMeasured->SetTitle("Measured");

    //unfold with randomized distributions
    Unfold();
    FillDeltaUnfoldedProfile();
  }

  // Get statistical errors for final unfolded spectrum
  // ie. spread of each pt bin in fDeltaUnfoldedP
  Double_t meanx2 = 0.;
  Double_t mean = 0.;
  Double_t checksigma = 0.;
  Double_t entriesInBin = 0.;
  for (Long_t iBin=0; iBin<fUnfoldedFinal->GetNbins(); iBin++) {
    fUnfoldedFinal->GetBinContent(iBin,fCoordinatesN_M);
    mean = fDeltaUnfoldedP->GetBinContent(fCoordinatesN_M);
    meanx2 = fDeltaUnfoldedP->GetBinError(fCoordinatesN_M);
    entriesInBin = fDeltaUnfoldedN->GetBinContent(fCoordinatesN_M);
    if(entriesInBin > 1.) checksigma = TMath::Sqrt((entriesInBin/(entriesInBin-1.))*TMath::Abs(meanx2-mean*mean));
    //printf("mean %f, meanx2 %f, sigmacheck %f, nentries %f\n",mean, meanx2, checksigma,entriesInBin);
    //AliDebug(2,Form("filling error %e\n",sigma));
    fUnfoldedFinal->SetBinError(fCoordinatesN_M,checksigma);
  }

  // now errors are calculated
  fNCalcCorrErrors = 2;
}

//______________________________________________________________
void AliCFUnfolding::CreateRandomizedDist() {
  //
  // Create randomized dist from original measured distribution
  // This distribution is created several times, each time with a different random number
  //

  for (Long_t iBin=0; iBin<fResponseOrig->GetNbins(); iBin++) {
    Double_t val = fResponseOrig->GetBinContent(iBin,fCoordinatesN_M); //used as mean
    Double_t err = fResponseOrig->GetBinError(fCoordinatesN_M);        //used as sigma
    Double_t ran = fRandom3->Gaus(val,err);
    // random        = fRandom3->PoissonD(measuredValue); //doesn't work for normalized spectra, use Gaus (assuming raw counts in bin is large >10)
    fRandomResponse->SetBinContent(iBin,ran);
  }
  for (Long_t iBin=0; iBin<fEfficiencyOrig->GetNbins(); iBin++) {
    Double_t val = fEfficiencyOrig->GetBinContent(iBin,fCoordinatesN_M); //used as mean
    Double_t err = fEfficiencyOrig->GetBinError(fCoordinatesN_M);        //used as sigma
    Double_t ran = fRandom3->Gaus(val,err);
    // random        = fRandom3->PoissonD(measuredValue); //doesn't work for normalized spectra, use Gaus (assuming raw counts in bin is large >10)
    fRandomEfficiency->SetBinContent(iBin,ran);
  }
  for (Long_t iBin=0; iBin<fMeasuredOrig->GetNbins(); iBin++) {
    Double_t val = fMeasuredOrig->GetBinContent(iBin,fCoordinatesN_M); //used as mean
    Double_t err = fMeasuredOrig->GetBinError(fCoordinatesN_M);        //used as sigma
    Double_t ran = fRandom3->Gaus(val,err);
    // random        = fRandom3->PoissonD(measuredValue); //doesn't work for normalized spectra, use Gaus (assuming raw counts in bin is large >10)
    fRandomMeasured->SetBinContent(iBin,ran);
  }
}

//______________________________________________________________
void AliCFUnfolding::FillDeltaUnfoldedProfile() {
  //
  // Store difference of unfolded spectrum from measured distribution and unfolded spectrum from randomized distribution
  // The delta profile has been set to a THnSparse to handle N dimension
  // The THnSparse contains in each bin the mean value and spread of the difference 
  // This function updates the profile wrt to its previous mean and error
  // The relation between iterations (n+1) and n is as follows :
  //  mean_{n+1} = (n*mean_n + value_{n+1}) / (n+1)
  // sigma_{n+1} = sqrt { 1/(n+1) * [ n*sigma_n^2 + (n^2+n)*(mean_{n+1}-mean_n)^2 ] }    (can this be optimized?)

  for (Long_t iBin=0; iBin<fUnfoldedFinal->GetNbins(); iBin++) {
    Double_t deltaInBin   = fUnfoldedFinal->GetBinContent(iBin,fCoordinatesN_M) - fUnfolded->GetBinContent(fCoordinatesN_M);
    Double_t entriesInBin = fDeltaUnfoldedN->GetBinContent(fCoordinatesN_M);
    //AliDebug(2,Form("%e %e ==> delta = %e\n",fUnfoldedFinal->GetBinContent(iBin,fCoordinatesN_M),fUnfolded->GetBinContent(iBin),deltaInBin));

    //printf("deltaInBin %f\n",deltaInBin);
    //printf("pt %f\n",ptaxis->GetBinCenter(iBin+1));
    
    Double_t mean_n = fDeltaUnfoldedP->GetBinContent(fCoordinatesN_M) ;
    Double_t mean_nplus1 = mean_n ;
    mean_nplus1 *= entriesInBin ;
    mean_nplus1 += deltaInBin ;
    mean_nplus1 /= (entriesInBin+1) ;

    Double_t meanx2_n = fDeltaUnfoldedP->GetBinError(fCoordinatesN_M) ;
    Double_t meanx2_nplus1 = meanx2_n ;
    meanx2_nplus1 *= entriesInBin ;
    meanx2_nplus1 += (deltaInBin*deltaInBin) ;
    meanx2_nplus1 /= (entriesInBin+1) ;

    //AliDebug(2,Form("sigma = %e\n",sigma));

    fDeltaUnfoldedP->SetBinError(fCoordinatesN_M,meanx2_nplus1) ;
    fDeltaUnfoldedP->SetBinContent(fCoordinatesN_M,mean_nplus1) ;
    fDeltaUnfoldedN->SetBinContent(fCoordinatesN_M,entriesInBin+1);
  }
}

//______________________________________________________________

void AliCFUnfolding::GetCoordinates() {
  //
  // assign coordinates in Measured and True spaces (dim=N) from coordinates in global space (dim=2N)
  //
  for (Int_t i = 0; i<fNVariables ; i++) {
    fCoordinatesN_M[i] = fCoordinates2N[i];
    fCoordinatesN_T[i] = fCoordinates2N[i+fNVariables];
  }
}

//______________________________________________________________

void AliCFUnfolding::CreateConditional() {
  //
  // creates the conditional probability matrix (R*) holding the P(M|T), given the reponse matrix R
  //
  //  --> R*(i,j) = R(i,j) / SUM_k{ R(k,j) }
  //

  fConditional = (THnSparse*) fResponse->Clone();  // output of this function

  Int_t* dim = new Int_t [fNVariables];
  for (Int_t iDim=0; iDim<fNVariables; iDim++) dim[iDim] = fNVariables+iDim ; //dimensions corresponding to TRUE values (i.e. from N to 2N-1)

  THnSparse* responseInT = fConditional->Projection(fNVariables,dim,"E");     // output denominator : 
                                                                              // projection of the response matrix on the TRUE axis
  delete [] dim; 

  // fill the conditional probability matrix
  for (Long_t iBin=0; iBin<fResponse->GetNbins(); iBin++) {
    Double_t responseValue = fResponse->GetBinContent(iBin,fCoordinates2N);
    GetCoordinates();
    Double_t projValue = responseInT->GetBinContent(fCoordinatesN_T);
   
    Double_t fill = responseValue / projValue ;
    if (fill>0. || fConditional->GetBinContent(fCoordinates2N)>0.) {
      fConditional->SetBinContent(fCoordinates2N,fill);
      Double_t err = 0.;
      fConditional->SetBinError  (fCoordinates2N,err);
    }
  }
  delete responseInT ;
}
//______________________________________________________________

Int_t AliCFUnfolding::GetDOF() {
  //
  // number of dof = number of bins
  //

  Int_t nDOF = 1 ;
  for (Int_t iDim=0; iDim<fNVariables; iDim++) {
    nDOF *= fPrior->GetAxis(iDim)->GetNbins();
  }
  AliDebug(0,Form("Number of degrees of freedom = %d",nDOF));
  return nDOF;
}

//______________________________________________________________

Double_t AliCFUnfolding::GetChi2() {
  //
  // Returns the chi2 between unfolded and a priori spectrum
  // This function became absolute with the correlated error calculation.
  // Errors on the unfolded distribution are not known until the end
  // Use the convergence criterion instead
  //

  Double_t chi2      = 0. ;
  Double_t error_unf = 0.;
  for (Long_t iBin=0; iBin<fPrior->GetNbins(); iBin++) {
    Double_t priorValue = fPrior->GetBinContent(iBin,fCoordinatesN_T);
    error_unf  = fUnfolded->GetBinError(fCoordinatesN_T);
    chi2 += (error_unf > 0. ? TMath::Power((fUnfolded->GetBinContent(fCoordinatesN_T) - priorValue)/error_unf,2) / priorValue : 0.) ;
  }
  return chi2;
}

//______________________________________________________________

Double_t AliCFUnfolding::GetConvergence() {
  //
  // Returns convergence criterion = \sum_t ((U_t^{n-1}-U_t^n)/U_t^{n-1})^2
  // U is unfolded spectrum, t is the bin, n = current, n-1 = previous
  //
  Double_t convergence = 0.;
  Double_t priorValue  = 0.;
  Double_t currentValue = 0.;
  for (Long_t iBin=0; iBin < fPrior->GetNbins(); iBin++) {
    priorValue = fPrior->GetBinContent(iBin,fCoordinatesN_T);
    currentValue = fUnfolded->GetBinContent(fCoordinatesN_T);

    if (priorValue > 0.)
      convergence += ((priorValue-currentValue)/priorValue)*((priorValue-currentValue)/priorValue);
    else 
      AliWarning(Form("priorValue = %f. Adding 0 to convergence criterion.",priorValue)); 
  }
  return convergence;
}

//______________________________________________________________

void AliCFUnfolding::SetMaxConvergencePerDOF(Double_t val) {
  //
  // Max. convergence criterion per degree of freedom : user setting
  // convergence criterion = DOF*val; DOF = number of bins
  // In Jan-Fiete's multiplicity note: Convergence criterion = DOF*0.001^2
  //

  Int_t nDOF = GetDOF() ;
  fMaxConvergence = val * nDOF ;
  AliInfo(Form("MaxConvergence = %e. Number of degrees of freedom = %d",fMaxConvergence,nDOF));
}

//______________________________________________________________

Short_t AliCFUnfolding::Smooth() {
  //
  // Smoothes the unfolded spectrum
  //
  // By default each cell content is replaced by the average with the neighbouring bins (but not diagonally-neighbouring bins)
  // However, if a specific function fcn has been defined in UseSmoothing(fcn), the unfolded will be fit and updated using fcn 
  //
  
  if (fSmoothFunction) {
    AliDebug(2,Form("Smoothing spectrum with fit function %p",fSmoothFunction));
    return SmoothUsingFunction();
  }
  else return SmoothUsingNeighbours(fUnfolded);
}

//______________________________________________________________

Short_t AliCFUnfolding::SmoothUsingNeighbours(THnSparse* hist) {
  //
  // Smoothes the unfolded spectrum using neighouring bins
  //

  Int_t  const nDimensions = hist->GetNdimensions() ;
  Int_t* coordinates = new Int_t[nDimensions];

  Int_t* numBins = new Int_t[nDimensions];
  for (Int_t iVar=0; iVar<nDimensions; iVar++) numBins[iVar] = hist->GetAxis(iVar)->GetNbins();
  
  //need a copy because hist will be updated during the loop, and this creates problems
  THnSparse* copy = (THnSparse*)hist->Clone();

  for (Long_t iBin=0; iBin<copy->GetNbins(); iBin++) { //loop on non-empty bins
    Double_t content = copy->GetBinContent(iBin,coordinates);
    Double_t error2  = TMath::Power(copy->GetBinError(iBin),2);

    // skip the under/overflow bins...
    Bool_t isOutside = kFALSE ;
    for (Int_t iVar=0; iVar<nDimensions; iVar++) {
      if (coordinates[iVar]<1 || coordinates[iVar]>numBins[iVar]) {
	isOutside=kTRUE;
	break;
      }
    }
    if (isOutside) continue;
    
    Int_t neighbours = 0; // number of neighbours to average with

    for (Int_t iVar=0; iVar<nDimensions; iVar++) {
      if (coordinates[iVar] > 1) { // must not be on low edge border
	coordinates[iVar]-- ; //get lower neighbouring bin 
	content += copy->GetBinContent(coordinates);
	error2  += TMath::Power(copy->GetBinError(coordinates),2);
	neighbours++;
	coordinates[iVar]++ ; //back to initial coordinate
      }
      if (coordinates[iVar] < numBins[iVar]) { // must not be on up edge border
	coordinates[iVar]++ ; //get upper neighbouring bin
	content += copy->GetBinContent(coordinates);
	error2  += TMath::Power(copy->GetBinError(coordinates),2);
	neighbours++;
	coordinates[iVar]-- ; //back to initial coordinate
      }
    }
    // make an average
    hist->SetBinContent(coordinates,content/(1.+neighbours));
    hist->SetBinError  (coordinates,TMath::Sqrt(error2)/(1.+neighbours));
  }
  delete [] numBins;
  delete [] coordinates ;
  delete copy;
  return 0;
}

//______________________________________________________________

Short_t AliCFUnfolding::SmoothUsingFunction() {
  //
  // Fits the unfolded spectrum using the function fSmoothFunction
  //

  AliDebug(0,Form("Smooth function is a %s with option \"%s\" and has %d parameters : ",fSmoothFunction->ClassName(),fSmoothOption,fSmoothFunction->GetNpar()));

  for (Int_t iPar=0; iPar<fSmoothFunction->GetNpar(); iPar++) AliDebug(0,Form("par[%d]=%e",iPar,fSmoothFunction->GetParameter(iPar)));

  Int_t fitResult = 0;

  switch (fNVariables) {
  case 1 : fitResult = fUnfolded->Projection(0)    ->Fit(fSmoothFunction,fSmoothOption); break;
  case 2 : fitResult = fUnfolded->Projection(1,0)  ->Fit(fSmoothFunction,fSmoothOption); break; // (1,0) instead of (0,1) -> TAxis issue
  case 3 : fitResult = fUnfolded->Projection(0,1,2)->Fit(fSmoothFunction,fSmoothOption); break;
  default: AliFatal(Form("Cannot handle such fit in %d dimensions",fNVariables)) ; return 1;
  }

  if (fitResult != 0) {
    AliWarning(Form("Fit failed with status %d, stopping the loop",fitResult));
    return 1;
  }

  Int_t nDim = fNVariables;
  Int_t* bins = new Int_t[nDim]; // number of bins for each variable
  Long_t nBins = 1; // used to calculate the total number of bins in the THnSparse

  for (Int_t iVar=0; iVar<nDim; iVar++) {
    bins[iVar] = fUnfolded->GetAxis(iVar)->GetNbins();
    nBins *= bins[iVar];
  }

  Int_t *bin  = new Int_t[nDim];    // bin to fill the THnSparse (holding the bin coordinates)
  Double_t x[3] = {0,0,0} ;         // value in bin center (max dimension is 3 (TF3))

  // loop on the bins and update of fUnfolded
  // THnSparse::Multiply(TF1*) doesn't exist, so let's do it bin by bin
  for (Long_t iBin=0; iBin<nBins; iBin++) {
    Long_t bin_tmp = iBin ;
    for (Int_t iVar=0; iVar<nDim; iVar++) {
      bin[iVar] = 1 + bin_tmp % bins[iVar] ;
      bin_tmp /= bins[iVar] ;
      x[iVar] = fUnfolded->GetAxis(iVar)->GetBinCenter(bin[iVar]);
    }
    Double_t functionValue = fSmoothFunction->Eval(x[0],x[1],x[2]) ;
    fUnfolded->SetBinError  (bin,fUnfolded->GetBinError(bin)*functionValue/fUnfolded->GetBinContent(bin));
    fUnfolded->SetBinContent(bin,functionValue);
  }
  delete [] bins;
  delete [] bin ;
  return 0;
}

//______________________________________________________________

void AliCFUnfolding::CreateFlatPrior() {
  //
  // Creates a flat prior distribution
  // 

  AliInfo("Creating a flat a priori distribution");
  
  // create the frame of the THnSparse given (for example) the one from the efficiency map
  fPrior = (THnSparse*) fEfficiency->Clone();
  fPrior->SetTitle("Prior");

  if (fNVariables != fPrior->GetNdimensions()) 
    AliFatal(Form("The prior matrix should have %d dimensions, and it has actually %d",fNVariables,fPrior->GetNdimensions()));

  Int_t nDim = fNVariables;
  Int_t* bins = new Int_t[nDim]; // number of bins for each variable
  Long_t nBins = 1; // used to calculate the total number of bins in the THnSparse

  for (Int_t iVar=0; iVar<nDim; iVar++) {
    bins[iVar] = fPrior->GetAxis(iVar)->GetNbins();
    nBins *= bins[iVar];
  }

  Int_t *bin = new Int_t[nDim]; // bin to fill the THnSparse (holding the bin coordinates)

  // loop that sets 1 in each bin
  for (Long_t iBin=0; iBin<nBins; iBin++) {
    Long_t bin_tmp = iBin ;
    for (Int_t iVar=0; iVar<nDim; iVar++) {
      bin[iVar] = 1 + bin_tmp % bins[iVar] ;
      bin_tmp /= bins[iVar] ;
    }
    fPrior->SetBinContent(bin,1.); // put 1 everywhere
    fPrior->SetBinError  (bin,0.); // put 0 everywhere
  }
  
  fPriorOrig = (THnSparse*)fPrior->Clone();

  delete [] bin;
  delete [] bins;
}
