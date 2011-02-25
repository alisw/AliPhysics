/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class to calculate statistical          //
// significance from AliMultiVeector objects with signal and     //
// background counts vs. cut values                              //
// Origin: Francesco Prino (prino@to.infn.it)                    //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliSignificanceCalculator.h"
#include "TMath.h"
#include "AliLog.h"

ClassImp(AliSignificanceCalculator)
//___________________________________________________________________________
AliSignificanceCalculator::AliSignificanceCalculator():
fSignal(0),
fErrSquareSignal(0),
fBackground(0),
fErrSquareBackground(0),
fSignificance(0),
fErrSignificance(0),
fNormSig(1.),  
fNormBkg(1.)
{
  // default constructor
  fSignal=new AliMultiDimVector();
  fBackground=new AliMultiDimVector();
}
//___________________________________________________________________________
AliSignificanceCalculator::AliSignificanceCalculator(AliMultiDimVector* sig, AliMultiDimVector* bkg, Float_t normsig, Float_t normbkg):
fSignal(sig),
fErrSquareSignal(0),
fBackground(bkg),
fErrSquareBackground(0),
fSignificance(0),
fErrSignificance(0),
fNormSig(normsig),
fNormBkg(normbkg)
{
  // standard constructor
  if(fSignal && fBackground) CalculateSignificance();
}
//___________________________________________________________________________
AliSignificanceCalculator::AliSignificanceCalculator(AliMultiDimVector* sig, AliMultiDimVector* bkg, AliMultiDimVector* err2sig, AliMultiDimVector* err2bkg, Float_t normsig, Float_t normbkg):
fSignal(sig),
fErrSquareSignal(err2sig),
fBackground(bkg),
fErrSquareBackground(err2bkg),
fSignificance(0),
fErrSignificance(0),
fNormSig(normsig),
fNormBkg(normbkg)
{
  // standard constructor
  if(fSignal && fBackground) CalculateSignificance();
}
//___________________________________________________________________________
AliSignificanceCalculator::~AliSignificanceCalculator(){
  // destructor
  if(fSignal) delete fSignal;
  if(fBackground) delete fBackground;
  if(fSignificance) delete fSignificance;
  if(fErrSignificance) delete fErrSignificance;
}
//___________________________________________________________________________
Bool_t AliSignificanceCalculator::Check() const {
  // checks AliMultiDimVector dimension and normalization
  if(fSignal==0 || fBackground==0) return kFALSE;
  if(fNormSig==0. || fNormBkg==0.) return kFALSE;
  if(fSignal->GetNTotCells() != fBackground->GetNTotCells()) return kFALSE;
  return kTRUE;
}
//___________________________________________________________________________
void AliSignificanceCalculator::CalculateSignificance(){
  // calculates significance and its error
  if(!Check()) AliFatal("Signal and Background AliMultiDimVector dimensions do not match!");
  
  if(fSignificance) delete fSignificance;
  if(fErrSignificance) delete fErrSignificance;
  fSignificance=new AliMultiDimVector();
  fSignificance->CopyStructure(fSignal);
  fErrSignificance=new AliMultiDimVector();
  fErrSignificance->CopyStructure(fSignal);

  for(ULong64_t i=0;i<fSignal->GetNTotCells();i++) {
    if(fSignal->GetElement(i)!=-1 && fBackground->GetElement(i)!=-1){
      Float_t s=fSignal->GetElement(i)*fNormSig;
      Float_t b=fBackground->GetElement(i)*fNormBkg;
      Float_t signif=0.;
      Float_t errsig=0.;
      if((s+b)>0){ 
	signif=s/TMath::Sqrt(s+b);
	Float_t errs,errb;
	if(fErrSquareSignal) errs=TMath::Sqrt(fErrSquareSignal->GetElement(i))*fNormSig;
	else errs=TMath::Sqrt(fSignal->GetElement(i))*fNormSig; // Poisson statistics
	if(fErrSquareBackground) errb=TMath::Sqrt(fErrSquareBackground->GetElement(i))*fNormBkg;
	else errb=TMath::Sqrt(fBackground->GetElement(i))*fNormBkg; // Poisson
	Float_t dsigds=(s+2*b)/2./(s+b)/TMath::Sqrt(s+b);
	Float_t dsigdb=-s/2./(s+b)/TMath::Sqrt(s+b);
	errsig=TMath::Sqrt(dsigds*dsigds*errs*errs+dsigdb*dsigdb*errb*errb);
      }
      fSignificance->SetElement(i,signif);
      fErrSignificance->SetElement(i,errsig);
    }
  }
  fSignificance->SetNameTitle("Significance","Significance");
  fErrSignificance->SetNameTitle("ErrorOnSignificance","ErrorOnSignificance");
}
//___________________________________________________________________________
AliMultiDimVector* AliSignificanceCalculator::CalculatePurity() const {
  // calculates purity 
  if(!Check()) AliFatal("Signal and Background AliMultiDimVector dimensions do not match!");
  
  AliMultiDimVector* purity=new AliMultiDimVector();
  purity->CopyStructure(fSignal);
  for(ULong64_t i=0;i<fSignal->GetNTotCells();i++) {
    if(fSignal->GetElement(i)!=-1 && fBackground->GetElement(i)!=-1){
      Float_t s=fSignal->GetElement(i)*fNormSig;
      Float_t b=fBackground->GetElement(i)*fNormBkg;
      Float_t pur=0.;
      if((s+b)>0) pur=s/(s+b);
      purity->SetElement(i,pur);
    }
  }
  purity->SetNameTitle("Purity","Purity");
  return purity;
}
//___________________________________________________________________________
AliMultiDimVector* AliSignificanceCalculator::CalculatePurityError() const {
  // calculates error on purity 
  if(!Check()) AliFatal("Signal and Background AliMultiDimVector dimensions do not match!");
  
  AliMultiDimVector* epurity=new AliMultiDimVector();
  epurity->CopyStructure(fSignal);
  for(ULong64_t i=0;i<fSignal->GetNTotCells();i++) {
    if(fSignal->GetElement(i)!=-1 && fBackground->GetElement(i)!=-1){
      Float_t s=fSignal->GetElement(i)*fNormSig;
      Float_t b=fBackground->GetElement(i)*fNormBkg;
      Float_t epur=0.;
      if((s+b)>0){
	Float_t errs,errb;
	if(fErrSquareSignal) errs=TMath::Sqrt(fErrSquareSignal->GetElement(i))*fNormSig;
	else errs=TMath::Sqrt(fSignal->GetElement(i))*fNormSig; // Poisson statistics
	if(fErrSquareBackground) errb=TMath::Sqrt(fErrSquareBackground->GetElement(i))*fNormBkg;
	else errb=TMath::Sqrt(fBackground->GetElement(i))*fNormBkg; // Poisson
	Float_t dpurds=b/(s+b)/(s+b);
	Float_t dpurdb=-s/(s+b)/(s+b);
	epur=TMath::Sqrt(dpurds*dpurds*errs*errs+dpurdb*dpurdb*errb*errb);
      }
      epurity->SetElement(i,epur);
    }
  }
  epurity->SetNameTitle("ErrorOnPurity","ErrorOnPurity");
  return epurity;
}
//___________________________________________________________________________
AliMultiDimVector* AliSignificanceCalculator::CalculateSOverB() const {
  // Signal over Background
  if(!Check()) AliFatal("Signal and Background AliMultiDimVector dimensions do not match!");
  
  AliMultiDimVector* sob=new AliMultiDimVector();
  sob->CopyStructure(fSignal);
  for(ULong64_t i=0;i<fSignal->GetNTotCells();i++) {
    if(fSignal->GetElement(i)!=-1 && fBackground->GetElement(i)!=-1){
      Float_t s=fSignal->GetElement(i)*fNormSig;
      Float_t b=fBackground->GetElement(i)*fNormBkg;
      Float_t soverb=0.;
      if(b>0) soverb=s/b;
      sob->SetElement(i,soverb);
    }
  }
  sob->SetNameTitle("SoverB","SoverB");
  return sob;
}
//___________________________________________________________________________
AliMultiDimVector* AliSignificanceCalculator::CalculateSOverBError() const {
  // Error on Signal over Background
  if(!Check()) AliFatal("Signal and Background AliMultiDimVector dimensions do not match!");
  
  AliMultiDimVector* esob=new AliMultiDimVector();
  esob->CopyStructure(fSignal);
  for(ULong64_t i=0;i<fSignal->GetNTotCells();i++) {
    if(fSignal->GetElement(i)!=-1 && fBackground->GetElement(i)!=-1){
      Float_t s=fSignal->GetElement(i)*fNormSig;
      Float_t b=fBackground->GetElement(i)*fNormBkg;
      Float_t esoverb=0.;
      if(b>0){
	Float_t soverb=s/b;
	Float_t errs,errb;
	if(fErrSquareSignal) errs=TMath::Sqrt(fErrSquareSignal->GetElement(i))*fNormSig;
	else errs=TMath::Sqrt(fSignal->GetElement(i))*fNormSig; // Poisson statistics
	if(fErrSquareBackground) errb=TMath::Sqrt(fErrSquareBackground->GetElement(i))*fNormBkg;
	else errb=TMath::Sqrt(fBackground->GetElement(i))*fNormBkg; // Poisson
	esoverb=soverb*TMath::Sqrt(errs*errs/s/s+errb*errb/b/b);
      }
      esob->SetElement(i,esoverb);
    }
  }
  esob->SetNameTitle("ErrorOnSoverB","ErrorOnSoverB");
  return esob;
}
