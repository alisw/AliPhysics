/*
***********************************************************
  Implementation of the AliResonanceFits
  Contact: i.c.arsene@cern.ch
  2014/10/09
  *********************************************************
*/

#ifndef ALIRESONANCEFITS_H
#include "AliResonanceFits.h"
#endif

#include <iostream>
#include <algorithm>
#include <iomanip>
using std::cout;
using std::endl;
using std::flush;
using std::swap;
using std::setw;

#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THn.h>
#include <TMath.h>
#include <TRandom.h>
#include <TMinuit.h>

#include "AliReducedVarManager.h"


ClassImp(AliResonanceFits)

// initialization of static members needed by the Minuit fitter
TH1* AliResonanceFits::fgTempSignal = 0x0;
TH1* AliResonanceFits::fgTempBkg = 0x0;
Bool_t AliResonanceFits::fgOptionUse2DMatching = kFALSE;
Double_t AliResonanceFits::fgPtFitRange[2] = {0.0, 100.};
Double_t AliResonanceFits::fgMassFitRange[2] = {0.0, 15.};
Int_t AliResonanceFits::fgOptionMEMatching = AliResonanceFits::kMatchSEOS;
Double_t AliResonanceFits::fgMassExclusionRange[2] = {2.5+1.0e-6, 3.2-1.0e-6};
Bool_t AliResonanceFits::fgOptionUseSignificantZero = kFALSE;

//_______________________________________________________________________________
AliResonanceFits::AliResonanceFits() :
  fSEOS(0x0),
  fSELSleg1(0x0),
  fSELSleg2(0x0),
  fMEOS(0x0),
  fMELSleg1(0x0),
  fMELSleg2(0x0),
  fNVariables(0),
  fVariables(),
  fVarLimits(),
  fVarBinLimits(),
  fVarIndices(),
  fMassVariable(AliReducedVarManager::kMass),
  fPtVariable(AliReducedVarManager::kPt),
  fNLoopingVariables(0),
  fCurrentVariable(0),
  fIter(),
  fOptionBkgMethod(kBkgMixedEvent),
  fOptionUseRfactorCorrection(kFALSE),
  fOptionScale(kScaleEntries),
  fOptionLSmethod(kLSGeometricMean),
  fWeightedAveragePower(2.0),
  fOptionMinuit(kMinuitMethodChi2),
  fOptionScaleSummedBkg(kFALSE),
  fUserEnabledMassFitRange(kFALSE),
  fUserEnabledPtFitRange(kFALSE),
  fUserEnabledMassExclusionRange(kFALSE),
  fSplusB(0x0),
  fBkg(0x0),
  fSig(0x0),
  fSoverB(0x0),
  fFitValues(),
  fMatchingIsDone(kFALSE),
  fMinuitFitter(0x0)
{
  //
  // Default constructor
  //
   for(Int_t i=0; i<kNMaxVariables; ++i) {
      fVariables[i] = -1; fVarLimits[i][0] = 0; fVarLimits[i][1] = 0;
      fVarIndices[i] = -1;
      fIter[i] = -1; fVarBinLimits[i][0] = -1; fVarBinLimits[i][1] = -1;
   }
   for(Int_t i=0; i<kNFitValues; ++i) fFitValues[i] = 0.;
}

//_______________________________________________________________________________
AliResonanceFits::~AliResonanceFits()
{
  //
  // De-constructor
  //
  if(fSEOS) delete fSEOS;
  if(fMEOS) delete fMEOS;
  if(fSELSleg1) delete fSELSleg1;
  if(fSELSleg2) delete fSELSleg2;
  if(fMELSleg1) delete fMELSleg1;
  if(fMELSleg2) delete fMELSleg2;
  if(fMinuitFitter) delete fMinuitFitter;
}

//_______________________________________________________________________________
void AliResonanceFits::SetHistograms(THnF* seos, THnF* meos /*=0x0*/, 
		                     THnF* selsLeg1 /*=0x0*/, THnF* selsLeg2 /*=0x0*/, 
		                     THnF* melsLeg1 /*=0x0*/, THnF* melsLeg2 /*=0x0*/) 
{
  //
  // Set the multi-dim histograms
  //
  fSEOS = seos; fMEOS = meos;
  fSELSleg1 = selsLeg1; fSELSleg2 = selsLeg2;
  fMELSleg1 = melsLeg1; fMELSleg2 = melsLeg2;
  fMatchingIsDone = kFALSE;
}


//_______________________________________________________________________________
void AliResonanceFits::AddVariables(Int_t nVars, Int_t* vars, Int_t* indices)
{
  //
  // initialize variable types and mapping of the dimensions in the THn histograms
  //
  for(Int_t i=0; i<nVars; ++i) {
     if(fNVariables>=kNMaxVariables) return;
     fVariables[fNVariables] = vars[i];
     fVarIndices[fNVariables] = indices[i];
     if(fSEOS) {
        fVarLimits[fNVariables][0] = fSEOS->GetAxis(indices[i])->GetXmin()+1.0e-6;
        fVarLimits[fNVariables][1] = fSEOS->GetAxis(indices[i])->GetXmax()-1.0e-6;
     }
     fNVariables++;
  }
  fMatchingIsDone = kFALSE;
}


//_______________________________________________________________________________
void AliResonanceFits::AddVariable(Int_t var, Int_t index)
{
   //
   // initialize a variable being used in the THnF
   //
   if(fNVariables>=kNMaxVariables) return;
   fVariables[fNVariables] = var;
   fVarIndices[fNVariables] = index;
   if(fSEOS) {
      fVarLimits[fNVariables][0] = fSEOS->GetAxis(index)->GetXmin()+1.0e-6;
      fVarLimits[fNVariables][1] = fSEOS->GetAxis(index)->GetXmax()-1.0e-6;
   }
   fNVariables++;
   fMatchingIsDone = kFALSE;
}


//_______________________________________________________________________________
void AliResonanceFits::SetVarRange(Int_t var, Double_t* lims) {
   //
   // set the user range for variable var 
   // NOTE: Variable var is encoded using the AliReducedVarManager::Variables enum
   // NOTE: If the var is not found in fVariables, nothing happens
   // NOTE: The user provided limits are slightly modified to avoid bin edge problems
   Int_t idx = -1;
   for(Int_t i=0; i<fNVariables; ++i)
      if(fVariables[i]==var) 
         idx = fVarIndices[i];
   if(idx==-1) return;
   fVarLimits[idx][0] = lims[0]+1.0e-6;
   fVarLimits[idx][1] = lims[1]-1.0e-6;
   fMatchingIsDone = kFALSE;
}


//_______________________________________________________________________________
void AliResonanceFits::ApplyUserRanges(THnF* h) {
   //
   // apply user ranges to the THnF histogram
   // NOTE: user ranges are applied only to the specified fNVariables variables
   //             The THnF may contain more than fNVariables, but those unspecified will be automatically integrtaed over
   // 
   for(Int_t i=0; i<fNVariables; ++i) 
      h->GetAxis(fVarIndices[i])->SetRangeUser(fVarLimits[i][0], fVarLimits[i][1]);
}


//_______________________________________________________________________________
Bool_t AliResonanceFits::Initialize() {
   //
   // make sure all prerequisites for signal extraction are met
   //   
   AliReducedVarManager::SetDefaultVarNames();
   
   // clean up the output histograms
   if(fSplusB) {delete fSplusB; fSplusB = 0;}
   if(fSig) {delete fSig; fSig = 0;}
   if(fBkg) {delete fBkg; fBkg = 0;}
   if(fSoverB) {delete fSoverB; fSoverB = 0;}
   
   // Check the needed user histograms
   if(!fSEOS) {
      cout << "AliResonanceFits::Initialize() Fatal: No SE-OS histogram provided! This is always needed" << endl;
      return kFALSE;
   }
   if(fOptionBkgMethod==kBkgMixedEvent || (fOptionBkgMethod==kBkgLikeSign && fOptionUseRfactorCorrection)) {
      if(!fMEOS) {
         cout << "AliResonanceFits::Initialize() Fatal: No ME-OS histogram provided! This is needed with the current matching options" << endl;
         return kFALSE;
      }
   }
   if(fOptionBkgMethod==kBkgLikeSign || (fOptionBkgMethod==kBkgMixedEvent && fgOptionMEMatching==kMatchSELS)) {
      if(!fSELSleg1) {
         cout << "AliResonanceFits::Initialize() Fatal: No SE-LS leg1 histogram provided! This is needed with the current matching options" << endl;
         return kFALSE;
      }
      if(!fSELSleg2) {
         cout << "AliResonanceFits::Initialize() Fatal: No SE-LS leg2 histogram provided! This is needed with the current matching options" << endl;
         return kFALSE;
      }
   }
   if((fOptionBkgMethod==kBkgMixedEvent && fgOptionMEMatching==kMatchSELS && fOptionUseRfactorCorrection) ||
      (fOptionBkgMethod==kBkgLikeSign && fOptionUseRfactorCorrection)) {
      if(!fMELSleg1) {
         cout << "AliResonanceFits::Initialize() Fatal: No ME-LS leg1 histogram provided! This is needed with the current matching options" << endl;
         return kFALSE;
      }
      if(!fMELSleg2) {
         cout << "AliResonanceFits::Initialize() Fatal: No ME-LS leg2 histogram provided! This is needed with the current matching options" << endl;
         return kFALSE;
      }
   }
   
   // at least one variable should be enabled (the invariant mass)
   if(fNVariables<1) {
      cout << "AliResonanceFits::Initialize() Fatal: No dimensions have been defined. Use AddVariable(s)  !" << endl;
      return kFALSE;      
   }
   // The mass variable should be among the user enabled variables.
   // Move the mass variable to be the last in the variable arrays
   Bool_t massVarFound = kFALSE;
   for(Int_t i=0; i<fNVariables; ++i) {
      if(fVariables[i]==fMassVariable) {
         swap(fVariables[fNVariables-1], fVariables[i]);
         swap(fVarIndices[fNVariables-1], fVarIndices[i]);
         swap(fVarLimits[fNVariables-1][0], fVarLimits[i][0]);
         swap(fVarLimits[fNVariables-1][1], fVarLimits[i][1]);
         massVarFound = kTRUE;
      }
   }
   if(!massVarFound) {
      cout << "AliResonanceFits::Initialize() Fatal: Mass dimension was not found in the list of user defined dimensions. " << endl;
      cout << "                      The current needed mass variable is " << AliReducedVarManager::fgVariableNames[fMassVariable] << endl;
      cout << "                 Use SetMassVariable() to change the desired mass variable, or define the mass dimension via AddVariable(s)()" << endl;
      return kFALSE;
   }
   if(!fUserEnabledMassFitRange) {
      cout << "AliResonanceFits::Initialize() Info: Mass fit range used from default values or those set via SetVarRange()" << endl;
      cout << "                    Use SetMassFitRange() if you want a different mass fit range" << endl;
      fgMassFitRange[0] = fVarLimits[fNVariables-1][0];
      fgMassFitRange[1] = fVarLimits[fNVariables-1][1];
   }
   
   // Look for the pt variable
   // Move the pt variable to be the next to last in the variable arrays
   Bool_t ptVarFound = kFALSE;
   for(Int_t i=0; i<fNVariables; ++i) {
      if(fVariables[i]==fPtVariable) {
         swap(fVariables[fNVariables-2], fVariables[i]);
         swap(fVarIndices[fNVariables-2], fVarIndices[i]);
         swap(fVarLimits[fNVariables-2][0], fVarLimits[i][0]);
         swap(fVarLimits[fNVariables-2][1], fVarLimits[i][1]);
         ptVarFound = kTRUE;
         break;
      }
   }
   if(!ptVarFound) {
      fPtVariable = -1;
      cout << "AliResonanceFits::Initialize() Info: No pt dimension indicated. It is fine with the current options." << endl;
   }
   if(!ptVarFound && (fgOptionUse2DMatching || fUserEnabledPtFitRange)) {
      cout << "AliResonanceFits::Initialize() Fatal: Pt dimension needed with the current options but was not found in the list of user defined dimensions  !" << endl;
      cout << "                      The current needed pt variable is " << AliReducedVarManager::fgVariableNames[fPtVariable] << endl;
      cout << "                 Use SetPtVariable() to change the desired pt variable, or define the pt dimension via AddVariable(s)()" << endl;
      return kFALSE;       // pt dimension needed for 2D matching or for selected fit pt range
   }
   if(ptVarFound && !fUserEnabledPtFitRange) {
      cout << "AliResonanceFits::Initialize() Info: Pt fit range used from default values or those set via SetVarRange()" << endl;
      cout << "                    Use SetPtFitRange() if you want a pt fit range different from the one used for signal counting" << endl;
      fgPtFitRange[0] = fVarLimits[fNVariables-2][0];
      fgPtFitRange[1] = fVarLimits[fNVariables-2][1];
   } 
   
   for(Int_t i=0; i<fNVariables; ++i) {
      fVarBinLimits[i][0] = fSEOS->GetAxis(fVarIndices[i])->FindBin(fVarLimits[i][0]);
      fVarBinLimits[i][1] = fSEOS->GetAxis(fVarIndices[i])->FindBin(fVarLimits[i][1]);
   }
   
   // Initialize looping variables =======================================
   // NOTE: These variables are used in the Slice() function, please make sure you know what you're doing
   //        The "looping variables" should be defined event-wise variables only (all added variables except mass and pt)
   fNLoopingVariables = (fPtVariable == -1 ? fNVariables-1 : fNVariables-2);
   fCurrentVariable = 0;
   for(Int_t i=0; i<fNLoopingVariables; ++i) fIter[i] = 0;
   // =======================================================
   
   ApplyUserRanges(fSEOS);
   if(fMEOS) ApplyUserRanges(fMEOS);
   if(fSELSleg1) ApplyUserRanges(fSELSleg1);
   if(fSELSleg2) ApplyUserRanges(fSELSleg2);
   if(fMELSleg1) ApplyUserRanges(fMELSleg1);
   if(fMELSleg2) ApplyUserRanges(fMELSleg2);
   
   // Print summary of all the options
   Print();
   
   return kTRUE;
}

//_______________________________________________________________________________
void AliResonanceFits::Slice() {
   //
   // Extract signal
   //
   
   if(fCurrentVariable==fNLoopingVariables) {
      // NOTE: The code in this if statement is run in the innermost loop
      
      // ************************************************************************************
      // NOTE: Needed for the proper handling of the recursive calling of the Slice function
      fCurrentVariable--;
      //************************************************************************************
      
      AddSlice();
   }
   // ************************************************************************************
   // NOTE: Needed for the proper handling of the recursive calling of the Slice function
   else {
      for(fIter[fCurrentVariable]=fVarBinLimits[fCurrentVariable][0]; fIter[fCurrentVariable]<=fVarBinLimits[fCurrentVariable][1]; fIter[fCurrentVariable]++) {
   // ************************************************************************************
         
         // outer loops -> set limits to the current variable axes (select individual bins)         
         fSEOS->GetAxis(fVarIndices[fCurrentVariable])->SetRange(fIter[fCurrentVariable], fIter[fCurrentVariable]);
         if(fOptionBkgMethod==kBkgMixedEvent || (fOptionBkgMethod==kBkgLikeSign && fOptionUseRfactorCorrection))
            fMEOS->GetAxis(fVarIndices[fCurrentVariable])->SetRange(fIter[fCurrentVariable], fIter[fCurrentVariable]);
         if(fOptionBkgMethod==kBkgLikeSign || (fOptionBkgMethod==kBkgMixedEvent && fgOptionMEMatching==kMatchSELS)) {
            fSELSleg1->GetAxis(fVarIndices[fCurrentVariable])->SetRange(fIter[fCurrentVariable], fIter[fCurrentVariable]);
            fSELSleg2->GetAxis(fVarIndices[fCurrentVariable])->SetRange(fIter[fCurrentVariable], fIter[fCurrentVariable]);
         }
         if((fOptionBkgMethod==kBkgMixedEvent && fgOptionMEMatching==kMatchSELS && fOptionUseRfactorCorrection) ||
            (fOptionBkgMethod==kBkgLikeSign && fOptionUseRfactorCorrection)) {
            fMELSleg1->GetAxis(fVarIndices[fCurrentVariable])->SetRange(fIter[fCurrentVariable], fIter[fCurrentVariable]);
            fMELSleg2->GetAxis(fVarIndices[fCurrentVariable])->SetRange(fIter[fCurrentVariable], fIter[fCurrentVariable]);
         }
         
   // ************************************************************************************
   // NOTE: Needed for the proper handling of the recursive calling of the Slice function
         fCurrentVariable++;
         Slice();
         while(fCurrentVariable>0) {
            if(fIter[fCurrentVariable]==fVarBinLimits[fCurrentVariable][1]) 
               fCurrentVariable--;
            else 
               break;
         }
      }
   }
   //************************************************************************************   
}

//____________________________________________________________________________________
void AliResonanceFits::AddSlice() {
   //
   // Take mass or (mass,pt) projections in the innermost loop, add them to the permanent
   //    histograms, construct backgrounds
   //
   TH1* projSEOS = 0x0;
   TH1* projMEOS = 0x0;
   TH1* projSELSleg1 = 0x0; TH1* projSELSleg2 = 0x0;
   TH1* projMELSleg1 = 0x0; TH1* projMELSleg2 = 0x0;
   // below are projections needed if a user pt range for the matching has been specified
   TH1* projSEOS_ptRange = 0x0;
   TH1* projMEOS_ptRange = 0x0;
   TH1* projSELSleg1_ptRange = 0x0; TH1* projSELSleg2_ptRange = 0x0;
   TH1* projMELSleg1_ptRange = 0x0; TH1* projMELSleg2_ptRange = 0x0;
   
   // NOTE: Remember, last variable in array is mass (fNVariables-1), and next to last, if its the case, is pt (fNVariables-2) 
   if(fgOptionUse2DMatching || fUserEnabledPtFitRange) {
      // SE-OS slices
      projSEOS = (TH2D*)fSEOS->Projection(fVarIndices[fNVariables-2], fVarIndices[fNVariables-1]);
      projSEOS->SetName(Form("projSEOS_%.6f", gRandom->Rndm()));
      if(!fgOptionUse2DMatching) 
         projSEOS_ptRange = ((TH2D*)projSEOS)->ProjectionX(Form("projSEOS_ptRange_%.6f", gRandom->Rndm()), 
                                                                                projSEOS->GetYaxis()->FindBin(fgPtFitRange[0]), 
                                                                                projSEOS->GetYaxis()->FindBin(fgPtFitRange[1]), "eo");
      // ME-OS slices
      if(fOptionBkgMethod==kBkgMixedEvent || (fOptionBkgMethod==kBkgLikeSign && fOptionUseRfactorCorrection)) {
         projMEOS = (TH2D*)fMEOS->Projection(fVarIndices[fNVariables-2], fVarIndices[fNVariables-1]);
         projMEOS->SetName(Form("projMEOS_%.6f", gRandom->Rndm()));
         if(!fgOptionUse2DMatching)
            projMEOS_ptRange = ((TH2D*)projMEOS)->ProjectionX(Form("projMEOS_ptRange_%.6f", gRandom->Rndm()), 
                                                                                    projMEOS->GetYaxis()->FindBin(fgPtFitRange[0]), 
                                                                                    projMEOS->GetYaxis()->FindBin(fgPtFitRange[1]), "eo");
      }
      // SE-LS slices
      if(fOptionBkgMethod==kBkgLikeSign || (fOptionBkgMethod==kBkgMixedEvent && fgOptionMEMatching==kMatchSELS)) {
         projSELSleg1 = (TH2D*)fSELSleg1->Projection(fVarIndices[fNVariables-2], fVarIndices[fNVariables-1]);
         projSELSleg1->SetName(Form("projSELSleg1_%.6f", gRandom->Rndm()));
         projSELSleg2 = (TH2D*)fSELSleg2->Projection(fVarIndices[fNVariables-2], fVarIndices[fNVariables-1]);
         projSELSleg2->SetName(Form("projSELSleg2_%.6f", gRandom->Rndm()));
         if(!fgOptionUse2DMatching) {
            projSELSleg1_ptRange = ((TH2D*)projSELSleg1)->ProjectionX(Form("projSELSleg1_ptRange_%.6f", gRandom->Rndm()), 
                                                                                                projSELSleg1->GetYaxis()->FindBin(fgPtFitRange[0]), 
                                                                                                projSELSleg1->GetYaxis()->FindBin(fgPtFitRange[1]), "eo");
            projSELSleg2_ptRange = ((TH2D*)projSELSleg2)->ProjectionX(Form("projSELSleg2_ptRange_%.6f", gRandom->Rndm()), 
                                                                                                projSELSleg2->GetYaxis()->FindBin(fgPtFitRange[0]), 
                                                                                                projSELSleg2->GetYaxis()->FindBin(fgPtFitRange[1]), "eo");
         }
      }
      // ME-LS slices
      if((fOptionBkgMethod==kBkgMixedEvent && fgOptionMEMatching==kMatchSELS && fOptionUseRfactorCorrection) ||
         (fOptionBkgMethod==kBkgLikeSign && fOptionUseRfactorCorrection)) {
         projMELSleg1 = (TH2D*)fMELSleg1->Projection(fVarIndices[fNVariables-2], fVarIndices[fNVariables-1]);
         projMELSleg1->SetName(Form("projMELSleg1_%.6f", gRandom->Rndm()));
         projMELSleg2 = (TH2D*)fMELSleg2->Projection(fVarIndices[fNVariables-2], fVarIndices[fNVariables-1]);
         projMELSleg2->SetName(Form("projMELSleg2_%.6f", gRandom->Rndm()));
         if(!fgOptionUse2DMatching) {
            projMELSleg1_ptRange = ((TH2D*)projMELSleg1)->ProjectionX(Form("projMELSleg1_ptRange_%.6f", gRandom->Rndm()), 
                                                                                                projMELSleg1->GetYaxis()->FindBin(fgPtFitRange[0]), 
                                                                                                projMELSleg1->GetYaxis()->FindBin(fgPtFitRange[1]), "eo");
            projMELSleg2_ptRange = ((TH2D*)projMELSleg2)->ProjectionX(Form("projMELSleg2_ptRange_%.6f", gRandom->Rndm()), 
                                                                                                projMELSleg2->GetYaxis()->FindBin(fgPtFitRange[0]), 
                                                                                                projMELSleg2->GetYaxis()->FindBin(fgPtFitRange[1]), "eo");
         }
      }
   }        // end if fgOptionUse2DMatching || fUserEnabledPtFitRange
   if(!fgOptionUse2DMatching) {        // use just 1D matching
      projSEOS = (TH1D*)fSEOS->Projection(fVarIndices[fNVariables-1]);
      projSEOS->SetName(Form("projSEOS_%.6f", gRandom->Rndm()));
      
      if(fOptionBkgMethod==kBkgMixedEvent || (fOptionBkgMethod==kBkgLikeSign && fOptionUseRfactorCorrection)) {
         projMEOS = (TH1D*)fMEOS->Projection(fVarIndices[fNVariables-1]);
         projMEOS->SetName(Form("projMEOS_%.6f", gRandom->Rndm()));
      }
      if(fOptionBkgMethod==kBkgLikeSign || (fOptionBkgMethod==kBkgMixedEvent && fgOptionMEMatching==kMatchSELS)) {
         projSELSleg1 = (TH1D*)fSELSleg1->Projection(fVarIndices[fNVariables-1]);
         projSELSleg1->SetName(Form("projSELSleg1_%.6f", gRandom->Rndm()));
         projSELSleg2 = (TH1D*)fSELSleg2->Projection(fVarIndices[fNVariables-1]);
         projSELSleg2->SetName(Form("projSELSleg2_%.6f", gRandom->Rndm()));
      }
      if((fOptionBkgMethod==kBkgMixedEvent && fgOptionMEMatching==kMatchSELS && fOptionUseRfactorCorrection) ||
         (fOptionBkgMethod==kBkgLikeSign && fOptionUseRfactorCorrection)) {
         projMELSleg1 = (TH1D*)fMELSleg1->Projection(fVarIndices[fNVariables-1]);
         projMELSleg1->SetName(Form("projMELSleg1_%.6f", gRandom->Rndm()));
         projMELSleg2 = (TH1D*)fMELSleg2->Projection(fVarIndices[fNVariables-1]);
         projMELSleg2->SetName(Form("projMELSleg2_%.6f", gRandom->Rndm()));
      }
   }          // end else
   
   // Add the temporary SEOS slice to the S+B histogram(s)
   if(!fSplusB) {
      if(fgOptionUse2DMatching) 
         fSplusB = (TH2D*)projSEOS->Clone(Form("fSplusB_%.6f", gRandom->Rndm()));
      else {
         if(!fUserEnabledPtFitRange) fSplusB = (TH1D*)projSEOS->Clone(Form("fSplusB_%.6f", gRandom->Rndm()));
         else fSplusB = ((TH2D*)projSEOS)->ProjectionX(Form("fSplusB_%.6f", gRandom->Rndm()), 0, -1, "eo");
      }
      fSplusB->SetDirectory(0x0);
   }
   else {
      if(fgOptionUse2DMatching || (!fgOptionUse2DMatching && !fUserEnabledPtFitRange))
         fSplusB->Add(projSEOS);
      else {
         TH1D* tempHist = ((TH2D*)projSEOS)->ProjectionX(Form("fSplusB_%.6f", gRandom->Rndm()), 0, -1, "eo");
         fSplusB->Add(tempHist);
         delete tempHist;
      }
   }
   
   TH1* bkgSlice = 0x0;
   // Construct the mixed event background, if its the case
   if(fOptionBkgMethod==kBkgMixedEvent) {
      TH1* scaleHist = 0x0;
      // Scale to the SE-OS in the side bands
      if(fgOptionMEMatching==kMatchSEOS) scaleHist = (!fgOptionUse2DMatching && fUserEnabledPtFitRange ? projSEOS_ptRange : projSEOS);
      // Scale to the SE-LS in the full allowed mass range
      if(fgOptionMEMatching==kMatchSELS) { 
         if(!fgOptionUse2DMatching && fUserEnabledPtFitRange)
            scaleHist = BuildLSbkg(projSELSleg1_ptRange, projSELSleg2_ptRange, projMEOS_ptRange, projMELSleg1_ptRange, projMELSleg2_ptRange);
         else
            scaleHist = BuildLSbkg(projSELSleg1, projSELSleg2, projMEOS, projMELSleg1, projMELSleg2);
      }
      TH1* bkgHist = projMEOS;
      if(!fgOptionUse2DMatching && fUserEnabledPtFitRange) bkgHist = projMEOS_ptRange;
      
      if(fOptionScale==kScaleEntries)                  ComputeEntryScale(scaleHist, bkgHist);
      if(fOptionScale==kScaleWeightedAverage) ComputeWeightedScale(scaleHist, bkgHist);
      if(fOptionScale==kScaleFit)                         FitScale(scaleHist, bkgHist);
      
      // scale the current bkg projection with the scale factor obtained above 
      projMEOS->Scale(fFitValues[kBkgScale]);
      bkgSlice = projMEOS;
   }     // end if mixed event bkg
   
   // Construct the like-sign background
   if(fOptionBkgMethod==kBkgLikeSign) {
      if(!fgOptionUse2DMatching && fUserEnabledPtFitRange)
         bkgSlice = BuildLSbkg(projSELSleg1_ptRange, projSELSleg2_ptRange, projMEOS_ptRange, projMELSleg1_ptRange, projMELSleg2_ptRange);
      else
         bkgSlice = BuildLSbkg(projSELSleg1, projSELSleg2, projMEOS, projMELSleg1, projMELSleg2);
   }
   
   // Add the current bkg slice to the total
   if(bkgSlice) {
      if(!fBkg) {
         if(fgOptionUse2DMatching) 
            fBkg = (TH2D*)bkgSlice->Clone(Form("fBkg_%.6f", gRandom->Rndm()));
         else 
            fBkg = (TH1D*)bkgSlice->Clone(Form("fBkg_%.6f", gRandom->Rndm()));
      }
      else
         fBkg->Add(bkgSlice);
   }
   
   delete projSEOS;
   if(projMEOS) delete projMEOS;
   if(projSELSleg1) delete projSELSleg1; 
   if(projSELSleg2) delete projSELSleg2;
   if(projMELSleg1) delete projMELSleg1; 
   if(projMELSleg2) delete projMELSleg2;
   if(projSEOS_ptRange) delete projSEOS_ptRange;
   if(projMEOS_ptRange) delete projMEOS_ptRange;
   if(projSELSleg1_ptRange) delete projSELSleg1_ptRange;
   if(projSELSleg2_ptRange) delete projSELSleg2_ptRange;
   if(projMELSleg1_ptRange) delete projMELSleg1_ptRange;
   if(projMELSleg2_ptRange) delete projMELSleg2_ptRange;
}

//_____________________________________________________________________________________________
TH1* AliResonanceFits::BuildLSbkg(TH1* selsLeg1, TH1* selsLeg2, TH1* meos /*=0x0*/, TH1* melsLeg1 /*=0x0*/, TH1* melsLeg2 /*=0x0*/) {
   //
   // Compute the SE-LS projection with R-factor correction if its the case 
   //
   TH1* mels = 0x0; TH1* sels = 0x0;
   if(fgOptionUse2DMatching)
      sels = (TH2D*)selsLeg1->Clone(Form("sels%.6f", gRandom->Rndm()));
   else
      sels = (TH1D*)selsLeg1->Clone(Form("sels%.6f", gRandom->Rndm()));
   
   if(fOptionUseRfactorCorrection) {
      if(fgOptionUse2DMatching)
         mels = (TH2D*)melsLeg1->Clone(Form("mels%.6f", gRandom->Rndm()));
      else
         mels = (TH1D*)melsLeg1->Clone(Form("mels%.6f", gRandom->Rndm()));
   }
   
   sels->Sumw2();
   selsLeg2->Sumw2();
   if(fOptionUseRfactorCorrection) {
      mels->Sumw2(); meos->Sumw2(); melsLeg2->Sumw2();
   }
   
   // Arithmetic mean
   if(fOptionLSmethod==kLSArithmeticMean) {
      sels->Add(selsLeg2);

      if(fOptionUseRfactorCorrection) {       // run R-factor correction
         mels->Add(melsLeg2);
         mels->Scale(0.5);
            
         TH1* rFactor = 0x0;
         if(fgOptionUse2DMatching) 
            rFactor = (TH2D*)meos->Clone(Form("rFactor%.6f", gRandom->Rndm()));
         else 
            rFactor = (TH1D*)meos->Clone(Form("rFactor%.6f", gRandom->Rndm())); 
         rFactor->Sumw2();
         rFactor->Divide(mels);
         sels->Multiply(rFactor);
         
         delete mels;
         delete rFactor;
      }
      return sels;
   }   // end if arithmetic mean
   
   // Geometric mean
   // Unfortunately there are no TH1 ROOT functions to compute the square root of a histogram
   //     so we need to do it by hand
   if(fOptionLSmethod==kLSGeometricMean) {
      sels->Multiply(selsLeg2);
      if(fOptionUseRfactorCorrection) mels->Multiply(melsLeg2);
      
      // compute the sqrt for the sels and mels histograms
      for(Int_t ipt=1; ipt<=(fgOptionUse2DMatching ? sels->GetYaxis()->GetNbins() : 1); ++ipt) {
         for(Int_t im=1; im<=sels->GetXaxis()->GetNbins(); ++im) {
            Double_t counts = (fgOptionUse2DMatching ? sels->GetBinContent(im,ipt) : sels->GetBinContent(im));
            Double_t countsErr = (fgOptionUse2DMatching ? sels->GetBinError(im,ipt) : sels->GetBinError(im));
            // NOTE: Relative error of sqrt(x) is 0.5 * relative error of x
            countsErr = (counts>1.0e-6 ? 0.5 * countsErr / counts : 0.0);
            if(fgOptionUse2DMatching) {
               sels->SetBinContent(im,ipt, (counts>1.0e-6 ? TMath::Sqrt(counts) : 0.0));
               sels->SetBinError(im,ipt, (counts>1.0e-6 ? countsErr*TMath::Sqrt(counts) : 0.0));
            }
            else {
               sels->SetBinContent(im, (counts>1.0e-6 ? TMath::Sqrt(counts) : 0.0));
               sels->SetBinError(im, (counts>1.0e-6 ? countsErr*TMath::Sqrt(counts) : 0.0));
            }
            
            if(!fOptionUseRfactorCorrection) continue;
            counts = (fgOptionUse2DMatching ? mels->GetBinContent(im,ipt) : mels->GetBinContent(im));
            countsErr = (fgOptionUse2DMatching ? mels->GetBinError(im,ipt) : mels->GetBinError(im));
            // NOTE: Relative error of sqrt(x) is 0.5 * relative error of x
            countsErr = (counts>1.0e-6 ? 0.5 * countsErr / counts : 0.0);
            if(fgOptionUse2DMatching) {
               mels->SetBinContent(im,ipt, (counts>1.0e-6 ? TMath::Sqrt(counts) : 0.0));
               mels->SetBinError(im,ipt, (counts>1.0e-6 ? countsErr*TMath::Sqrt(counts) : 0.0));
            }
            else {
               mels->SetBinContent(im, (counts>1.0e-6 ? TMath::Sqrt(counts) : 0.0));
               mels->SetBinError(im, (counts>1.0e-6 ? countsErr*TMath::Sqrt(counts) : 0.0));
            }
         }  // end loop over mass bins
      }  // end loop over pt bins
      sels->Scale(2.0);
      if(!fOptionUseRfactorCorrection)
         return sels;
      
      TH1* rFactor = 0x0;
      if(fgOptionUse2DMatching) 
         rFactor = (TH2D*)meos->Clone(Form("rFactor%.6f", gRandom->Rndm()));
      else 
         rFactor = (TH1D*)meos->Clone(Form("rFactor%.6f", gRandom->Rndm())); 
      rFactor->Sumw2();
      rFactor->Divide(mels);
      sels->Multiply(rFactor);
      
      delete mels;
      delete rFactor;
      return sels;
   }  // end if(geometric mean)
   
   return 0x0;
}

//____________________________________________________________________________________
void AliResonanceFits::ComputeEntryScale(TH1* sig, TH1* bkg) {
   //
   // Compute the scale of the bkg histogram to sig based on the number of entries in the fitting range fgMassFitRange
   //    and excluding fgMassExclusionRange
   // NOTE: The exclusion range is considered to be a sub-interval of the fgMassFitRange
   //
   Double_t entriesSig = 0.; Double_t entriesSigErr = 0.;
   Double_t entriesSigExcl = 0.; Double_t entriesSigExclErr = 0.;
   Double_t entriesBkg = 0.; Double_t entriesBkgErr = 0.;
   Double_t entriesBkgExcl = 0.; Double_t entriesBkgExclErr = 0.;
   if(fgOptionUse2DMatching) {
      entriesSig = ((TH2*)sig)->IntegralAndError(sig->GetXaxis()->FindBin(fgMassFitRange[0]), sig->GetXaxis()->FindBin(fgMassFitRange[1]),
                                                              sig->GetXaxis()->FindBin(fgPtFitRange[0]), sig->GetXaxis()->FindBin(fgPtFitRange[1]), entriesSigErr);
      if(fgOptionMEMatching != kMatchSELS) 
         entriesSigExcl = ((TH2*)sig)->IntegralAndError(sig->GetXaxis()->FindBin(fgMassExclusionRange[0]), sig->GetXaxis()->FindBin(fgMassExclusionRange[1]),
                                                                       sig->GetXaxis()->FindBin(fgPtFitRange[0]), sig->GetXaxis()->FindBin(fgPtFitRange[1]), entriesSigExclErr);
         entriesBkg = ((TH2*)bkg)->IntegralAndError(bkg->GetXaxis()->FindBin(fgMassFitRange[0]), bkg->GetXaxis()->FindBin(fgMassFitRange[1]),
                                                                bkg->GetXaxis()->FindBin(fgPtFitRange[0]), bkg->GetXaxis()->FindBin(fgPtFitRange[1]), entriesBkgErr);
      if(fgOptionMEMatching != kMatchSELS)
         entriesBkgExcl = ((TH2*)bkg)->IntegralAndError(bkg->GetXaxis()->FindBin(fgMassExclusionRange[0]), bkg->GetXaxis()->FindBin(fgMassExclusionRange[1]),
                                                                          bkg->GetXaxis()->FindBin(fgPtFitRange[0]), bkg->GetXaxis()->FindBin(fgPtFitRange[1]), entriesBkgExclErr);
   }
   else {
      entriesSig = sig->IntegralAndError(sig->GetXaxis()->FindBin(fgMassFitRange[0]), 
                                                              sig->GetXaxis()->FindBin(fgMassFitRange[1]), entriesSigErr);
      if(fgOptionMEMatching != kMatchSELS)
         entriesSigExcl = sig->IntegralAndError(sig->GetXaxis()->FindBin(fgMassExclusionRange[0]), 
                                                                       sig->GetXaxis()->FindBin(fgMassExclusionRange[1]), entriesSigExclErr);
      entriesBkg = bkg->IntegralAndError(bkg->GetXaxis()->FindBin(fgMassFitRange[0]), 
                                                                 bkg->GetXaxis()->FindBin(fgMassFitRange[1]), entriesBkgErr);
      if(fgOptionMEMatching != kMatchSELS)
         entriesBkgExcl = bkg->IntegralAndError(bkg->GetXaxis()->FindBin(fgMassExclusionRange[0]), 
                                                                          bkg->GetXaxis()->FindBin(fgMassExclusionRange[1]), entriesBkgExclErr);
   }
   // If not matching to the SE-LS subtract the yield in the exclusion range from the total; recompute uncertainty using quadrature
   if(fgOptionMEMatching != kMatchSELS) {
      entriesSig -= entriesSigExcl;
      entriesBkg -= entriesBkgExcl;
      entriesSigErr = TMath::Sqrt(entriesSigErr*entriesSigErr + entriesSigExclErr*entriesSigExclErr);
      entriesBkgErr = TMath::Sqrt(entriesBkgErr*entriesBkgErr + entriesBkgExclErr*entriesBkgExclErr);
   }
   
   fFitValues[kBkgScale] = (entriesSig>1.0e-6 && entriesBkg>1.0e-6 ? entriesSig/entriesBkg : 0.0);
   fFitValues[kBkgScaleErr] = 0.0;
   if(entriesSig>1.0e-6 && entriesBkg>1.0e-6)
      fFitValues[kBkgScaleErr] = fFitValues[kBkgScale]*TMath::Sqrt(entriesSigErr*entriesSigErr/entriesSig/entriesSig + 
                                                                                                         entriesBkgErr*entriesBkgErr/entriesBkg/entriesBkg);
}

//____________________________________________________________________________________
void AliResonanceFits::ComputeWeightedScale(TH1* sig, TH1* bkg) {
   //
   // Compute the scale of the bkg histogram to sig based on the number of entries in the fitting range fgMassFitRange
   //    and excluding fgMassExclusionRange
   // NOTE: The exclusion range is considered to be a sub-interval of the fgMassFitRange
   //
   
   // obtain the S/B histogram
   TH1* soverb = 0x0;
   if(fgOptionUse2DMatching) 
      soverb = (TH2D*)sig->Clone(Form("soverb_%.6f", gRandom->Rndm()));
   else
      soverb = (TH1D*)sig->Clone(Form("soverb_%.6f", gRandom->Rndm()));
   soverb->Divide(bkg);   
   
   // loop to compute the weighted average
   Double_t sweights = 0.0; Double_t avWeights = 0.0; Double_t nMassBins=0;
   Double_t serror = 0.0;
   
   for(Int_t ipt=1; ipt<=(fgOptionUse2DMatching ? soverb->GetYaxis()->GetNbins() : 1); ++ipt) {
      Float_t pt = (fgOptionUse2DMatching ? soverb->GetYaxis()->GetBinCenter(ipt) : 0.0);
      if(fgOptionUse2DMatching && (pt<fgPtFitRange[0] || pt>fgPtFitRange[1])) continue;      // only the selected pt fit range 
      
      for(Int_t im=1; im<=sig->GetXaxis()->GetNbins(); ++im) {
         Double_t m = soverb->GetXaxis()->GetBinCenter(im);
         if(m<fgMassFitRange[0] || m>fgMassFitRange[1]) continue;    // only the selected mass fit range
         if(fgOptionMEMatching != kMatchSELS)    // exclude the region around the peak only if not matching to the SE-LS
            if(m>fgMassExclusionRange[0] && m<fgMassExclusionRange[1]) continue;     // exclude the fgMassExclusionRange
         
         Double_t s = (fgOptionUse2DMatching ? soverb->GetBinContent(im,ipt) : soverb->GetBinContent(im)); 
         Double_t sErr = (fgOptionUse2DMatching ? soverb->GetBinError(im,ipt) : soverb->GetBinError(im));
         if(sErr<1.0e-5) continue;
         
         // weighting using S/B error
         sweights  += 1.0/TMath::Power(sErr, fWeightedAveragePower);
         avWeights += s/TMath::Power(sErr, fWeightedAveragePower);
         serror += TMath::Power(sErr, 2.0-2.0*fWeightedAveragePower);
         
         nMassBins += 1.0;
      }   // end loop over mass bins
   }   // end loop over pt bins
   
   delete soverb;
   
   if(sweights>0.0) avWeights /= sweights;
   fFitValues[kBkgScale] = avWeights;  
   
   if(sweights>0.0) fFitValues[kBkgScaleErr] = TMath::Sqrt(serror)/sweights;
   else fFitValues[kBkgScaleErr] = TMath::Sqrt(serror);
}

//____________________________________________________________________________________
void AliResonanceFits::FitScale(TH1* sig, TH1* bkg, Bool_t fixScale /*=kFALSE*/) {
   //
   // Compute the bkg scaling by fitting
   //
   fgTempSignal = sig;
   fgTempBkg = bkg;
   
   Double_t arglist[2];
   Int_t ierflg=0;
   // NOTE: Initialize Minuit to fit one single parameter (the scale)
   if(!fMinuitFitter) {
      fMinuitFitter = new TMinuit(1);
      fMinuitFitter->SetFCN(Fcn);
      
      if(fOptionMinuit==kMinuitMethodChi2) arglist[0] = 1.0;
      if(fOptionMinuit==kMinuitMethodLikelihood) arglist[0] = 0.5;
      fMinuitFitter->mnexcm("SET ERR", arglist, 1, ierflg);
   }
   
   // Set starting values and step sizes for parameters
   Double_t scaleStart=1.0;
   Double_t step=0.01;
     
   fMinuitFitter->mnparm(0, "scale", scaleStart, step, 0.0, 1000000.0, ierflg);
   // NOTE: use the fixScale method just to compute the chi2
   if(fixScale) fMinuitFitter->FixParameter(fFitValues[kBkgScale]);
                      
   arglist[0] = 500;
   arglist[1] = 1.;
   
   fMinuitFitter->SetMaxIterations(10000);
   fMinuitFitter->mnexcm("SIMPLEX", arglist ,2,ierflg);
   fMinuitFitter->mnexcm("MIGRAD", arglist ,2,ierflg);
   
   fMinuitFitter->GetParameter(0, fFitValues[kBkgScale], fFitValues[kBkgScaleErr]);
}

//____________________________________________________________________________________
void AliResonanceFits::Fcn(Int_t&, Double_t*, Double_t &f, Double_t *par, Int_t) {
   //
   // chi2 function used as interface with minuit
   // par[0] - scale factor for the background histogram
   //
   f = AliResonanceFits::Chi2(fgTempSignal, fgTempBkg, par[0]);
}

//____________________________________________________________________________________
Double_t AliResonanceFits::Chi2(TH1* sig, TH1* bkg, Double_t scale, Double_t scaleError /*=0.0*/) {
   //
   // Compute the chi2 for the difference between the signal and scaled background
   // Assume signal and background uncertainties are uncorrelated
   //
   // NOTE:
   // In case of 2D matching, the signal and bkg are TH2D histograms with mass on X-axis and pt on Y-axis
   //
   Float_t chi2 = 0.0;
   Int_t ndf = 0;
   
   for(Int_t ipt=1; ipt<=(fgOptionUse2DMatching ? sig->GetYaxis()->GetNbins() : 1); ++ipt) {
      Float_t pt = (fgOptionUse2DMatching ? sig->GetYaxis()->GetBinCenter(ipt) : 0.0);
      if(fgOptionUse2DMatching && (pt<fgPtFitRange[0] || pt>fgPtFitRange[1])) continue;      // only the selected pt fit range 
      
      for(Int_t im=1; im<=sig->GetXaxis()->GetNbins(); ++im) {
         Double_t m = sig->GetXaxis()->GetBinCenter(im);
         if(m<fgMassFitRange[0] || m>fgMassFitRange[1]) continue;    // only the selected mass fit range
         if(fgOptionMEMatching != kMatchSELS)         //  exclude the region around the peak only if not matching to the SE-LS
            if(m>fgMassExclusionRange[0] && m<fgMassExclusionRange[1]) continue;     // exclude the fgMassExclusionRange
            
         Double_t sigVal = (fgOptionUse2DMatching ? sig->GetBinContent(im,ipt) : sig->GetBinContent(im));
         if(!fgOptionUseSignificantZero && sigVal<0.0001) continue;
         Double_t bkgVal = (fgOptionUse2DMatching ? bkg->GetBinContent(im,ipt) : bkg->GetBinContent(im));
         if(bkgVal<=0.0001) continue;
         Double_t sigErr = (fgOptionUse2DMatching ? sig->GetBinError(im,ipt) : sig->GetBinError(im));
         if(sigVal<0.0001)   // when considering zero entry bins as significant, assume error to be 1
            sigErr = 1.0;
         Double_t bkgErr = (fgOptionUse2DMatching ? bkg->GetBinError(im,ipt) : bkg->GetBinError(im));
         
         Float_t err = 0.0;
         if(scale>0.0)
            err = bkgVal*scale*TMath::Sqrt(scaleError*scaleError/scale/scale + bkgErr*bkgErr/bkgVal/bkgVal);
         err = sigErr*sigErr+err*err;
         
         chi2 += (err>0.0 ? (sigVal-scale*bkgVal)*(sigVal-scale*bkgVal)/err : 0.0);
         
         ndf += 1;
      }   // end loop over mass bins
   }         // end loop over pt bins
   
   return (ndf > 0 ? chi2/Double_t(ndf) : 1000.);
}

//_______________________________________________________________________________
Bool_t AliResonanceFits::Process() {
   //
   // Main function handling the signal extraction
   //
   fMatchingIsDone = kFALSE;
   
   // Initialize and make sure all prerequisites are met
   Bool_t initState = Initialize();
   if(!initState) return kFALSE;
                      
   // Loop over the defined dimensions and build the SplusB and bkg histograms
   if(fOptionBkgMethod==kBkgMixedEvent || fOptionBkgMethod==kBkgLikeSign)
   Slice();
   if(fOptionBkgMethod==kBkgFunction) {
      // TODO
   }

   if(fOptionScaleSummedBkg) {
      if(fOptionScale==kScaleEntries)                  ComputeEntryScale(fSplusB, fBkg);
      if(fOptionScale==kScaleWeightedAverage) ComputeWeightedScale(fSplusB, fBkg);
      if(fOptionScale==kScaleFit)                         FitScale(fSplusB, fBkg);
      fBkg->Scale(fFitValues[kBkgScale]);
   }
   
   // build the signal projection
   if(fSig) {delete fSig; fSig = 0;}
   if(fgOptionUse2DMatching)
      fSig = (TH2D*)fSplusB->Clone(Form("fSig_%.6f", gRandom->Rndm()));
   else
      fSig = (TH1D*)fSplusB->Clone(Form("fSig_%.6f", gRandom->Rndm()));
   fSig->Add(fBkg, -1.0);
   
   // build the S/B projection
   if(fSoverB) {delete fSoverB; fSoverB = 0;}
   if(fgOptionUse2DMatching)
      fSoverB = (TH2D*)fSig->Clone(Form("fSoverB_%.6f", gRandom->Rndm()));
   else
      fSoverB = (TH1D*)fSig->Clone(Form("fSoverB_%.6f", gRandom->Rndm()));
   fSoverB->Divide(fBkg);
   
   fMatchingIsDone = kTRUE;
   
   return kTRUE;
}

//____________________________________________________________________________________
Double_t* AliResonanceFits::ComputeOutputValues(Double_t minMass, Double_t maxMass, Double_t minPt /*=-1.*/, Double_t maxPt /*=-1.*/) {
   //
   // compute signal, bkg, S/B, significance, etc.
   // Signal is integrated in the specified mass and optionally pt range specified
   // NOTE: 
   // the pt interval used for signal integration may be different wrt pt interval used for matching
   //            The pt range at this step is available only in the case of 2D matching
   
   if(!fMatchingIsDone) {
      cout  << "AliResonanceFits::ComputeOutputValues(): Matching / fitting procedure was not performed, so no values can be computed" << endl;
      cout  <<  "      Please run AliResonanceFits::Process() succesfully first !" << endl;
      return 0x0;
   }
   
   Int_t minMassBin = fSig->GetXaxis()->FindBin(minMass+1.0e-6);
   Int_t maxMassBin = fSig->GetXaxis()->FindBin(maxMass-1.0e-6);
   if(fgOptionUse2DMatching) {
      // if min and max pt are not specified, then integrate over the full available pt range
      Int_t minPtBin = (minPt<0. ? 1 : fSig->GetYaxis()->FindBin(minPt+1.0e-6));
      Int_t maxPtBin = (maxPt<0. ? fSig->GetYaxis()->GetNbins() : fSig->GetYaxis()->FindBin(maxPt-1.0e-6));
      
      fFitValues[kSig] = ((TH2*)fSig)->IntegralAndError(minMassBin, maxMassBin, minPtBin, maxPtBin, fFitValues[kSigErr]);
      fFitValues[kBkg] = ((TH2*)fBkg)->IntegralAndError(minMassBin, maxMassBin, minPtBin, maxPtBin, fFitValues[kBkgErr]);
      fFitValues[kSplusB] = ((TH2*)fSplusB)->IntegralAndError(minMassBin, maxMassBin, minPtBin, maxPtBin, fFitValues[kSplusBerr]);
   }
   else {
      fFitValues[kSig] = fSig->IntegralAndError(minMassBin, maxMassBin, fFitValues[kSigErr]);
      fFitValues[kBkg] = fBkg->IntegralAndError(minMassBin, maxMassBin, fFitValues[kBkgErr]);
      fFitValues[kSplusB] = fSplusB->IntegralAndError(minMassBin, maxMassBin, fFitValues[kSplusBerr]);
   }
   if(fOptionScaleSummedBkg) {
      // Update the signal error with the error from the bkg scaling
      fFitValues[kSigErr] = TMath::Sqrt(fFitValues[kBkg]*fFitValues[kBkg]*fFitValues[kBkgScaleErr]*fFitValues[kBkgScaleErr]/fFitValues[kBkgScale]/fFitValues[kBkgScale] + fFitValues[kSigErr]*fFitValues[kSigErr]);
   }
   
   fFitValues[kSoverB] = (fFitValues[kBkg]>0.001 ? fFitValues[kSig] / fFitValues[kBkg] : 0.0);
   fFitValues[kSoverBerr] = TMath::Sqrt(fFitValues[kSigErr]*fFitValues[kSigErr]/fFitValues[kSig]/fFitValues[kSig] + 
   fFitValues[kBkgErr]*fFitValues[kBkgErr]/fFitValues[kBkg]/fFitValues[kBkg]);
   fFitValues[kSoverBerr] *= fFitValues[kSoverB];
   
   if(fOptionBkgMethod==kBkgMixedEvent)
      fFitValues[kSignif] = ((fFitValues[kSig]+fFitValues[kBkg]>0.001) ? (fFitValues[kSig]/(TMath::Sqrt(fFitValues[kSig]+fFitValues[kBkg]))) : 0.0);
   if(fOptionBkgMethod==kBkgLikeSign)
      fFitValues[kSignif] = ((fFitValues[kSig]+2.0*fFitValues[kBkg]>0.001) ? (fFitValues[kSig]/(TMath::Sqrt(fFitValues[kSig]+2.0*fFitValues[kBkg]))) : 0.0);
   fFitValues[kChisqSideBands] = Chi2(fSig, fBkg, 1.0, 0.0);
   
   return fFitValues;
}

//_____________________________________________________________________________________________
void AliResonanceFits::PrintFitValues() {
   //
   // print fit values
   //
   if(!fMatchingIsDone) {
      cout << "AliResonanceFits::PrintFitValues()  Matching/fitting was not performed or the current fit values are not consistent with the current user options" << endl;
      cout << "          Run Process() first" << endl;
      return;
   }
   cout << endl;
   cout << "AliResonanceFits summary of matching/fit result ======================================" << endl;
   cout << std::left << setw(20) << "Signal counts" << " :: " << fFitValues[kSig] << " +/- " << fFitValues[kSigErr] << endl;
   cout << setw(20) << "Bkg counts" << " :: " << fFitValues[kBkg] << " +/- " << fFitValues[kBkgErr] << endl;
   cout << setw(20) << "Signal + bkg counts" << " :: " << fFitValues[kSplusB] << " +/- " << fFitValues[kSplusBerr] << endl;
   cout << setw(20) << "S/B" << " :: " << fFitValues[kSoverB] << " +/- " << fFitValues[kSoverBerr] << endl;
   cout << setw(20) << (fOptionBkgMethod==kBkgMixedEvent ? "S/sqrt(S+B)" : "S/sqrt(S+2B)") << " :: " << fFitValues[kSignif] << endl;
   cout << setw(20) << "Chi2" << " :: " << fFitValues[kChisqSideBands] << endl;
   if(fOptionScaleSummedBkg)
      cout << setw(20) << "Bkg scale" << " :: " << fFitValues[kBkgScale] << " +/- " << fFitValues[kBkgScaleErr] << endl;
   cout << "===================================================================" << endl;
}

//_____________________________________________________________________________________________
void AliResonanceFits::Print() {
   //
   // Print all user options
   //
   cout << endl;
   cout << "AliResonanceFits summary of all user options ======================================" << endl;
   cout << "fSEOS ::\t" << fSEOS << endl;
   if(fMEOS) cout << "fMEOS ::\t" << fMEOS << endl;
   if(fSELSleg1) cout << "fSELSleg1 ::\t" << fSELSleg1 << endl;
   if(fSELSleg2) cout << "fSELSleg2 ::\t" << fSELSleg2 << endl;
   if(fMELSleg1) cout << "fMELSleg1 ::\t" << fMELSleg1 << endl;
   if(fMELSleg2) cout << "fMELSleg2 ::\t" << fMELSleg2 << endl;
   cout << "fNVariables ::\t" << fNVariables << endl << endl;
   for(Int_t i=0; i<fNVariables; ++i) {
      cout << "Dim #" << i << "   Name: " << std::left << setw(25) << AliReducedVarManager::fgVariableNames[fVariables[i]] << "  ";
      cout << "index #" << fVarIndices[i] << "   ";
      cout << "user range: [" << std::right << setw(7) << fVarLimits[i][0] << ", " << setw(7) << fVarLimits[i][1] << "]   ";
      cout << "bin range: [" << std::right << setw(5) << fVarBinLimits[i][0] << ", " << setw(5) << fVarBinLimits[i][1] << "]" << endl;
   }
   cout << endl;
   cout << std::left << setw(43) << "fMassVariable" << " :: " << AliReducedVarManager::fgVariableNames[fMassVariable] << endl;
   cout << std::left << setw(43) << "fPtVariable" << " :: " << AliReducedVarManager::fgVariableNames[fPtVariable] << endl;
   cout << std::left << setw(43) << "fNLoopingVariables" << " :: " << fNLoopingVariables << endl;
   cout << std::left << setw(43) << "Use 2D matching" << " :: " << fgOptionUse2DMatching << endl;
   cout << std::left << setw(43) << "Bkg method" << " :: " << (fOptionBkgMethod==kBkgMixedEvent ? "Mixed event" : (fOptionBkgMethod==kBkgLikeSign ? "Like-sign" : "Fit function")) << endl;
   if(fOptionBkgMethod==kBkgMixedEvent)
      cout << std::left << setw(43) << "Mixed event bkg matching" << " :: " << (fgOptionMEMatching==kMatchSEOS ? "Sidebands of SE-OS" : "Full range of SE-LS") << endl;
   if(fOptionBkgMethod==kBkgLikeSign || fgOptionMEMatching==kMatchSELS)
      cout << std::left << setw(43) << "Use R-factor correction for LS bkg" << " :: " << fOptionUseRfactorCorrection << endl;
   if(fOptionBkgMethod==kBkgMixedEvent)
      cout << std::left << setw(43) << "Bkg scaling option" << " :: " << (fOptionScale==kScaleEntries ? "Counts" : (fOptionScale==kScaleWeightedAverage ? "Weighted average" : "Fit")) << endl;
   if(fOptionBkgMethod==kBkgMixedEvent && fOptionScale==kScaleWeightedAverage)
      cout << std::left << setw(43) << "Stat error weights inverse power" << " :: " << fWeightedAveragePower << endl;
   if(fOptionBkgMethod==kBkgMixedEvent && fOptionScale==kScaleFit)
      cout << std::left << setw(43) << "Minuit fit method" << " :: " << (fOptionMinuit==kMinuitMethodChi2 ? "Chi2 minimization" : "Likelihood maximization") << endl;
   if(fOptionBkgMethod==kBkgMixedEvent && fOptionScale==kScaleFit)
      cout << std::left << setw(43) << "Set zero entry bins as significant in Minuit fit" << " :: " << fgOptionUseSignificantZero << endl;
   if(fOptionBkgMethod==kBkgLikeSign || fgOptionMEMatching==kMatchSELS)
      cout << std::left << setw(43) << "Like-sign bkg method" << " :: " << (fOptionLSmethod==kLSGeometricMean ? "Geometric mean" : "Arithmetic mean") << endl;
   if(fOptionBkgMethod==kBkgMixedEvent || fOptionBkgMethod==kBkgLikeSign)
      cout << std::left << setw(43) << "Run additional scaling after all summations" << " :: " << fOptionScaleSummedBkg << endl;
   
   cout << "Mass fit range " << (fUserEnabledMassFitRange ? "(Set via SetMassFitRange()) :: [" : " :: [") 
   << fgMassFitRange[0] << " - " << fgMassFitRange[1] << "] GeV/c^2" << endl;
   cout << "Pt fit range " << (fUserEnabledPtFitRange ? "(Set via SetPtFitRange()) :: [" : " :: [") 
   << fgPtFitRange[0] << " - " << fgPtFitRange[1] << "] GeV/c" << endl;
   cout << "Mass exclusion range " << (fUserEnabledMassExclusionRange ? "(Set via SetMassExclusionRange()) :: [" : " :: [") 
   << fgMassExclusionRange[0] << " - " << fgMassExclusionRange[1] << "] GeV/c^2" << endl;
   
   cout << "==========================================================" << endl;
   cout << endl;
}
