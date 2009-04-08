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

//--------------------------------------------------------------------//
//                                                                    //
// AliCFUnfolding Class                                               //
// Class to handle general unfolding procedure                        // 
// For the moment only bayesian unfolding is supported                //
// The next steps are to add chi2 minimisation and weighting methods  //
//                                                                    //
// Author : renaud.vernet@cern.ch                                     //
//--------------------------------------------------------------------//


#include "AliCFUnfolding.h"
#include "TMath.h"
#include "TAxis.h"
#include "AliLog.h"

ClassImp(AliCFUnfolding)

//______________________________________________________________

AliCFUnfolding::AliCFUnfolding() :
  TNamed(),
  fResponse(0x0),
  fPrior(0x0),
  fOriginalPrior(0x0),
  fEfficiency(0x0),
  fMeasured(0x0),
  fMaxNumIterations(0),
  fNVariables(0),
  fMaxChi2(0),
  fUseSmoothing(kFALSE),
  fInverseResponse(0x0),
  fMeasuredEstimate(0x0),
  fConditional(0x0),
  fProjResponseInT(0x0),
  fUnfolded(0x0),
  fCoordinates2N(0x0),
  fCoordinatesN_M(0x0),
  fCoordinatesN_T(0x0)
{
  //
  // default constructor
  //
}

//______________________________________________________________

AliCFUnfolding::AliCFUnfolding(const Char_t* name, const Char_t* title, const Int_t nVar, 
			       const THnSparse* response, const THnSparse* efficiency, const THnSparse* measured, const THnSparse* prior) :
  TNamed(name,title),
  fResponse((THnSparse*)response->Clone()),
  fPrior(0x0),
  fOriginalPrior(0x0),
  fEfficiency((THnSparse*)efficiency->Clone()),
  fMeasured((THnSparse*)measured->Clone()),
  fMaxNumIterations(0),
  fNVariables(nVar),
  fMaxChi2(0),
  fUseSmoothing(kFALSE),
  fInverseResponse(0x0),
  fMeasuredEstimate(0x0),
  fConditional(0x0),
  fProjResponseInT(0x0),
  fUnfolded(0x0),
  fCoordinates2N(0x0),
  fCoordinatesN_M(0x0),
  fCoordinatesN_T(0x0)
{
  //
  // named constructor
  //

  AliInfo(Form("\n\n--------------------------\nCreating an unfolder :\n--------------------------\nresponse matrix has %d dimension(s)",fResponse->GetNdimensions()));
  
  if (!prior) CreateFlatPrior(); // if no prior distribution declared, simply use a flat distribution
  else {
    fPrior = (THnSparse*) prior->Clone();
    fOriginalPrior = (THnSparse*)fPrior->Clone();
  }
  
  for (Int_t iVar=0; iVar<fNVariables; iVar++) {
    AliInfo(Form("prior      matrix has %d bins in dimension %d",fPrior     ->GetAxis(iVar)->GetNbins(),iVar));
    AliInfo(Form("efficiency matrix has %d bins in dimension %d",fEfficiency->GetAxis(iVar)->GetNbins(),iVar));
    AliInfo(Form("measured   matrix has %d bins in dimension %d",fMeasured  ->GetAxis(iVar)->GetNbins(),iVar));
  }
  Init();
}


//______________________________________________________________

AliCFUnfolding::AliCFUnfolding(const AliCFUnfolding& c) :
  TNamed(c),
  fResponse((THnSparse*)c.fResponse->Clone()),
  fPrior((THnSparse*)c.fPrior->Clone()),
  fOriginalPrior((THnSparse*)c.fOriginalPrior->Clone()),
  fEfficiency((THnSparse*)c.fEfficiency->Clone()),
  fMeasured((THnSparse*)c.fMeasured->Clone()),
  fMaxNumIterations(c.fMaxNumIterations),
  fNVariables(c.fNVariables),
  fMaxChi2(c.fMaxChi2),
  fUseSmoothing(c.fUseSmoothing),
  fInverseResponse((THnSparse*)c.fInverseResponse->Clone()),
  fMeasuredEstimate((THnSparse*)fMeasuredEstimate->Clone()),
  fConditional((THnSparse*)c.fConditional->Clone()),
  fProjResponseInT((THnSparse*)c.fProjResponseInT->Clone()),
  fUnfolded((THnSparse*)c.fUnfolded->Clone()),
  fCoordinates2N(new Int_t(*c.fCoordinates2N)),
  fCoordinatesN_M(new Int_t(*c.fCoordinatesN_M)),
  fCoordinatesN_T(new Int_t(*c.fCoordinatesN_T))
{
  //
  // copy constructor
  //
}

//______________________________________________________________

AliCFUnfolding& AliCFUnfolding::operator=(const AliCFUnfolding& c) {
  //
  // assignment operator
  //
  
  if (this!=&c) {
    TNamed::operator=(c);
    fResponse = (THnSparse*)c.fResponse->Clone() ;
    fPrior = (THnSparse*)c.fPrior->Clone() ;
    fOriginalPrior = (THnSparse*)c.fOriginalPrior->Clone() ;
    fEfficiency = (THnSparse*)c.fEfficiency->Clone() ;
    fMeasured = (THnSparse*)c.fMeasured->Clone() ;
    fMaxNumIterations = c.fMaxNumIterations ;
    fNVariables = c.fNVariables ;
    fMaxChi2 = c.fMaxChi2 ;
    fUseSmoothing = c.fUseSmoothing ;
    fInverseResponse = (THnSparse*)c.fInverseResponse->Clone() ;
    fMeasuredEstimate = (THnSparse*)fMeasuredEstimate->Clone() ;
    fConditional = (THnSparse*)c.fConditional->Clone() ;
    fProjResponseInT = (THnSparse*)c.fProjResponseInT->Clone() ;
    fUnfolded = (THnSparse*)c.fUnfolded->Clone() ;
    fCoordinates2N  = new Int_t(*c.fCoordinates2N)  ;
    fCoordinatesN_M = new Int_t(*c.fCoordinatesN_M) ;
    fCoordinatesN_T = new Int_t(*c.fCoordinatesN_T) ;
  }
  return *this;
}

//______________________________________________________________

AliCFUnfolding::~AliCFUnfolding() {
  //
  // destructor
  //
  if (fResponse)           delete fResponse;
  if (fPrior)              delete fPrior;
  if (fOriginalPrior)      delete fOriginalPrior;
  if (fEfficiency)         delete fEfficiency;
  if (fMeasured)           delete fMeasured;
  if (fInverseResponse)    delete fInverseResponse;
  if (fMeasuredEstimate)   delete fMeasuredEstimate;
  if (fConditional)        delete fConditional;
  if (fProjResponseInT)    delete fProjResponseInT;
  if (fCoordinates2N)      delete [] fCoordinates2N; 
  if (fCoordinatesN_M)     delete [] fCoordinatesN_M; 
  if (fCoordinatesN_T)     delete [] fCoordinatesN_T; 
}

//______________________________________________________________

void AliCFUnfolding::Init() {
  //
  // initialisation function : creates internal settings
  //

  fCoordinates2N  = new Int_t[2*fNVariables];
  fCoordinatesN_M = new Int_t[fNVariables];
  fCoordinatesN_T = new Int_t[fNVariables];

  // create the matrix of conditional probabilities P(M|T)
  CreateConditional();
  
  // create the frame of the inverse response matrix
  fInverseResponse  = (THnSparse*) fResponse->Clone();
  // create the frame of the unfolded spectrum
  fUnfolded = (THnSparse*) fPrior->Clone();
  // create the frame of the measurement estimate spectrum
  fMeasuredEstimate = (THnSparse*) fMeasured->Clone();
}

//______________________________________________________________

void AliCFUnfolding::CreateEstMeasured() {
  //
  // This function creates a estimate (M) of the reconstructed spectrum 
  // given the a priori distribution (T) and the conditional matrix (COND)
  //
  // --> P(M) = SUM   { P(M|T)    * P(T) }
  // --> M(i) = SUM_k { COND(i,k) * T(k) }
  //
  // This is needed to calculate the inverse response matrix
  //


  // clean the measured estimate spectrum
  for (Long64_t i=0; i<fMeasuredEstimate->GetNbins(); i++) {
    fMeasuredEstimate->GetBinContent(i,fCoordinatesN_M);
    fMeasuredEstimate->SetBinContent(fCoordinatesN_M,0.);
  }
  
  // fill it
  for (Int_t iBin=0; iBin<fConditional->GetNbins(); iBin++) {
    Double_t conditionalValue = fConditional->GetBinContent(iBin,fCoordinates2N);
    GetCoordinates();
    Double_t priorValue = fPrior->GetBinContent(fCoordinatesN_T);
    Double_t fill = fMeasuredEstimate->GetBinContent(fCoordinatesN_M) + conditionalValue * priorValue * fEfficiency->GetBinContent(fCoordinatesN_T);
    if (fill>0.) fMeasuredEstimate->SetBinContent(fCoordinatesN_M,fill);
  }
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

  for (Int_t iBin=0; iBin<fConditional->GetNbins(); iBin++) {
    Double_t conditionalValue = fConditional->GetBinContent(iBin,fCoordinates2N);
    GetCoordinates();
    Double_t priorValue = fPrior->GetBinContent(fCoordinatesN_T);
    Double_t estimatedMeasured = fMeasuredEstimate->GetBinContent(fCoordinatesN_M);
    Double_t fill = (estimatedMeasured>0. ? conditionalValue * priorValue * fEfficiency->GetBinContent(fCoordinatesN_T) / estimatedMeasured : 0. ) ;
    if (fill>0. || fInverseResponse->GetBinContent(fCoordinates2N)>0.) fInverseResponse->SetBinContent(fCoordinates2N,fill);
  }
}

//______________________________________________________________

void AliCFUnfolding::Unfold() {
  //
  // Main routine called by the user : 
  // it calculates the unfolded spectrum from the response matrix and the measured spectrum
  // several iterations are performed until a reasonable chi2 is reached
  //

  Int_t iIterBayes=0 ;
  Double_t chi2=0 ;

  for (iIterBayes=0; iIterBayes<fMaxNumIterations; iIterBayes++) { // bayes iterations
    CreateEstMeasured();
    CreateInvResponse();
    CreateUnfolded();
    chi2 = GetChi2();
    //printf("chi2 = %e\n",chi2);
    if (fMaxChi2>0. && chi2<fMaxChi2) {
      break;
    }
    // update the prior distribution
    if (fUseSmoothing) Smooth();
    fPrior = (THnSparse*)fUnfolded->Clone() ; // this should be changed (memory)
  }
  AliInfo(Form("Finished at iteration %d : Chi2 is %e and you required it to be < %e",iIterBayes,chi2,fMaxChi2));
}

//______________________________________________________________

void AliCFUnfolding::CreateUnfolded() {
  //
  // Creates the unfolded (T) spectrum from the measured spectrum (M) and the inverse response matrix (INV)
  // We have P(T) = SUM   { P(T|M)   * P(M) } 
  //   -->   T(i) = SUM_k { INV(i,k) * M(k) }
  //


  // clear the unfolded spectrum
  for (Long64_t i=0; i<fUnfolded->GetNbins(); i++) {
    fUnfolded->GetBinContent(i,fCoordinatesN_T);
    fUnfolded->SetBinContent(fCoordinatesN_T,0.);
  }
  
  for (Int_t iBin=0; iBin<fInverseResponse->GetNbins(); iBin++) {
    Double_t invResponseValue = fInverseResponse->GetBinContent(iBin,fCoordinates2N);
    GetCoordinates();
    Double_t effValue = fEfficiency->GetBinContent(fCoordinatesN_T);
    Double_t fill = fUnfolded->GetBinContent(fCoordinatesN_T) + (effValue>0. ? invResponseValue*fMeasured->GetBinContent(fCoordinatesN_M)/effValue : 0.) ;
    if (fill>0.) fUnfolded->SetBinContent(fCoordinatesN_T,fill);
  }
}
    
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

  fConditional     = (THnSparse*) fResponse->Clone(); // output of this function
  fProjResponseInT = (THnSparse*) fPrior->Clone();    // output denominator : 
                                                      // projection of the response matrix on the TRUE axis
  
  // set in fProjResponseInT zero everywhere
  for (Int_t iBin=0; iBin<fProjResponseInT->GetNbins(); iBin++) {
    fProjResponseInT->GetBinContent(iBin,fCoordinatesN_T);
    fProjResponseInT->SetBinContent(fCoordinatesN_T,0.);
  }

  // calculate the response projection on T axis
  for (Int_t iBin=0; iBin<fResponse->GetNbins(); iBin++) {
    Double_t responseValue = fResponse->GetBinContent(iBin,fCoordinates2N);
    GetCoordinates();
    Double_t fill = fProjResponseInT->GetBinContent(fCoordinatesN_T) + responseValue ;
    if (fill>0.) fProjResponseInT->SetBinContent(fCoordinatesN_T,fill);
  }
  
  // fill the conditional probability matrix
  for (Int_t iBin=0; iBin<fResponse->GetNbins(); iBin++) {
    Double_t responseValue = fResponse->GetBinContent(iBin,fCoordinates2N);
    GetCoordinates();
    Double_t fill = responseValue / fProjResponseInT->GetBinContent(fCoordinatesN_T) ;
    if (fill>0. || fConditional->GetBinContent(fCoordinates2N)) fConditional->SetBinContent(fCoordinates2N,fill);
  }
}

//______________________________________________________________

Double_t AliCFUnfolding::GetChi2() {
  //
  // Returns the chi2 between unfolded and a priori spectrum
  //

  Double_t chi2 = 0. ;
  for (Int_t iBin=0; iBin<fPrior->GetNbins(); iBin++) {
    Double_t priorValue = fPrior->GetBinContent(iBin);
    chi2 += (priorValue>0. ? TMath::Power(fUnfolded->GetBinContent(iBin) - priorValue,2) / priorValue : 0.) ;
  }
  return chi2;
}

//______________________________________________________________

void AliCFUnfolding::SetMaxChi2PerDOF(Double_t val) {
  //
  // Max. chi2 per degree of freedom : user setting
  //

  Int_t nDOF = 1 ;
  for (Int_t iDim=0; iDim<fNVariables; iDim++) {
    nDOF *= fPrior->GetAxis(iDim)->GetNbins();
  }
  AliInfo(Form("Number of degrees of freedom = %d",nDOF));
  fMaxChi2 = val * nDOF ;
}

//______________________________________________________________

void AliCFUnfolding::Smooth() {
  //
  // Smoothes the unfolded spectrum
  // Each cell content is replaced by the average with the neighbouring bins (but not diagonally-neighbouring bins)
  //
  
  Int_t* numBins = new Int_t[fNVariables];
  for (Int_t iVar=0; iVar<fNVariables; iVar++) numBins[iVar]=fUnfolded->GetAxis(iVar)->GetNbins();
  
  //need a copy because fUnfolded will be updated during the loop, and this creates problems
  THnSparse* copy = (THnSparse*)fUnfolded->Clone();

  for (Int_t iBin=0; iBin<copy->GetNbins(); iBin++) { //loop on non-empty bins
    Double_t content = copy->GetBinContent(iBin,fCoordinatesN_T);

    // skip the under/overflow bins...
    Bool_t isOutside = kFALSE ;
    for (Int_t iVar=0; iVar<fNVariables; iVar++) {
      if (fCoordinatesN_T[iVar]<1 || fCoordinatesN_T[iVar]>numBins[iVar]) {
	isOutside=kTRUE;
	break;
      }
    }
    if (isOutside) continue;
    
    Int_t neighbours = 0; // number of neighbours to average with

    for (Int_t iVar=0; iVar<fNVariables; iVar++) {
      if (fCoordinatesN_T[iVar] > 1) { // must not be on low edge border
	fCoordinatesN_T[iVar]-- ; //get lower neighbouring bin 
	Double_t contentNeighbour = copy->GetBinContent(fCoordinatesN_T);
	content += contentNeighbour;
	neighbours++;
	fCoordinatesN_T[iVar]++ ; //back to initial coordinate
      }
      if (fCoordinatesN_T[iVar] < numBins[iVar]) { // must not be on up edge border
	fCoordinatesN_T[iVar]++ ; //get upper neighbouring bin
	Double_t contentNeighbour = copy->GetBinContent(fCoordinatesN_T);
	content += contentNeighbour ;
	neighbours++;
	fCoordinatesN_T[iVar]-- ; //back to initial coordinate
      }
    }
    content /= (1+neighbours) ; // make an average
    fUnfolded->SetBinContent(fCoordinatesN_T,content);
  }
  delete [] numBins;
  delete copy;
}


//______________________________________________________________

void AliCFUnfolding::CreateFlatPrior() {
  //
  // Creates a flat prior distribution
  // 

  AliInfo("Creating a flat a priori distribution");
  
  // create the frame of the THnSparse given (for example) the one from the efficiency map
  fPrior = (THnSparse*) fEfficiency->Clone();

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
  }
  
  fOriginalPrior = (THnSparse*)fPrior->Clone();

  delete [] bin;
  delete [] bins;
}
