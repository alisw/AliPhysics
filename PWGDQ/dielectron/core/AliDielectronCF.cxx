
/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//       Dielectron Correction framework manager                         //
//                                                                       //
/*









*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TList.h>
#include <TObjArray.h>
#include <TVectorD.h>
#include <TString.h>
#include <TObjString.h>

#include <AliCFContainer.h>
#include <AliAnalysisFilter.h>
#include <AliAnalysisCuts.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliDielectronCF.h"
#include "AliDielectronMC.h"
#include "AliDielectronPair.h"
#include "AliDielectronSignalMC.h"

ClassImp(AliDielectronCF)

AliDielectronCF::AliDielectronCF() :
  TNamed("DielectronCF","DielectronCF"),
  fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
  fNSteps(0),
  fNVars(0),
  fVarBinLimits(0x0),
  fNVarsLeg(0),
  fVarBinLimitsLeg(0x0),
  fNCuts(0),
  fValues(0x0),
  fIsMCTruth(0x0),
  fStepForMCtruth(kFALSE),
  fStepForNoCutsMCmotherPid(kFALSE),
  fStepForAfterAllCuts(kTRUE),
  fStepForPreFilter(kFALSE),
  fStepsForEachCut(kFALSE),
  fStepsForCutsIncreasing(kFALSE),
  fStepsForSignal(kTRUE),
  fStepsForBackground(kFALSE),
  fStepsForMCtruthOnly(kFALSE),
  fNStepMasks(0),
  fPdgMother(-1),
  fSignalsMC(0x0),
  fCfContainer(0x0),
  fHasMC(kFALSE),
  fNAddSteps(0)
{
  //
  // Default constructor
  //
  for (Int_t i=0; i<AliDielectronVarManager::kNMaxValues; ++i){
    fVariables[i]=0;
    fVariablesLeg[i]=0;
  }

  for (Int_t i=0; i<kNmaxAddSteps; ++i){
    fStepMasks[i]=0xFFFFFF;
  }
}

//________________________________________________________________
AliDielectronCF::AliDielectronCF(const char* name, const char* title) :
  TNamed(name, title),
  fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
  fNSteps(0),
  fNVars(0),
  fVarBinLimits(0x0),
  fNVarsLeg(0),
  fVarBinLimitsLeg(0x0),
  fNCuts(0),
  fValues(0x0),
  fIsMCTruth(0x0),
  fStepForMCtruth(kFALSE),
  fStepForNoCutsMCmotherPid(kFALSE),
  fStepForAfterAllCuts(kTRUE),
  fStepForPreFilter(kFALSE),
  fStepsForEachCut(kFALSE),
  fStepsForCutsIncreasing(kFALSE),
  fStepsForSignal(kTRUE),
  fStepsForBackground(kFALSE),
  fStepsForMCtruthOnly(kFALSE),
  fNStepMasks(0),
  fPdgMother(-1),
  fSignalsMC(0x0),
  fCfContainer(0x0),
  fHasMC(kFALSE),
  fNAddSteps(0)
{
  //
  // Named constructor
  //
  for (Int_t i=0; i<AliDielectronVarManager::kNMaxValues; ++i){
    fVariables[i]=0;
    fVariablesLeg[i]=0;
  }
  
  for (Int_t i=0; i<kNmaxAddSteps; ++i){
    fStepMasks[i]=0xFFFFFF;
  }
}

//________________________________________________________________
AliDielectronCF::~AliDielectronCF()
{
  //
  // Destructor
  //
  if (fUsedVars) delete fUsedVars;
  if (fValues) delete [] fValues;
  if (fIsMCTruth) delete [] fIsMCTruth;
  if (fVarBinLimits) delete fVarBinLimits;
  if (fVarBinLimitsLeg) delete fVarBinLimitsLeg;
}

//________________________________________________________________
void AliDielectronCF::AddVariable(AliDielectronVarManager::ValueTypes type, Int_t nbins,
                                  Double_t min, Double_t max, Bool_t leg, Bool_t log)
{
  //
  // Add a variable to the CF configuration
  // if leg is true it will add the variables of the leg
  // if log is true log binning will be created
  //

  TVectorD *binLimits=0x0;
  if (!log) binLimits=MakeLinBinning(nbins,min,max);
  else binLimits=MakeLogBinning(nbins,min,max);
  AddVariable(type,binLimits,leg);
}

//________________________________________________________________
void AliDielectronCF::AddVariable(AliDielectronVarManager::ValueTypes type, const char* binLimitStr, Bool_t leg/*=kFALSE*/)
{
  //
  // Add a variable to the CF configuration
  // specify arbitrary binning in a string.
  // Bin limits need to be separated by a ","
  //
  TString limits(binLimitStr);
  if (limits.IsNull()){
    AliError(Form("Bin Limit string is empty, cannot add the variable '%s'",AliDielectronVarManager::GetValueName(type)));
    return;
  }

  TObjArray *arr=limits.Tokenize(",");
  Int_t nLimits=arr->GetEntries();
  if (nLimits<2){
    AliError(Form("Need at leas 2 bin limits, cannot add the variable '%s'",AliDielectronVarManager::GetValueName(type)));
    delete arr;
    return;
  }

  TVectorD *binLimits=new TVectorD(nLimits);
  for (Int_t iLim=0; iLim<nLimits; ++iLim){
    (*binLimits)[iLim]=(static_cast<TObjString*>(arr->At(iLim)))->GetString().Atof();
  }

  delete arr;
  AddVariable(type,binLimits,leg);
}

//________________________________________________________________
void AliDielectronCF::AddVariable(AliDielectronVarManager::ValueTypes type, TVectorD *binLimits, Bool_t leg/*=kFALSE*/)
{
  //
  // Add variable with the binning given in the TVectorD
  //
  fUsedVars->SetBitNumber(type,kTRUE);
  if (!leg){
    if (!fVarBinLimits){
      fVarBinLimits=new TObjArray;
      fVarBinLimits->SetOwner();
    }
    fVarBinLimits->Add(binLimits);
    fVariables[fNVars]  = (UInt_t)type;
    ++fNVars;
  } else {
    if (!fVarBinLimitsLeg){
      fVarBinLimitsLeg=new TObjArray;
      fVarBinLimitsLeg->SetOwner();
    }
    fVarBinLimitsLeg->Add(binLimits);
    fVariablesLeg[fNVarsLeg]  = (UInt_t)type;
    ++fNVarsLeg;
  }
}

//________________________________________________________________
void AliDielectronCF::InitialiseContainer(const AliAnalysisFilter& filter)
{
  //
  // Initialise container based on the cuts in the analysis filter
  //

  fNCuts=filter.GetCuts()->GetEntries();

  fHasMC=AliDielectronMC::Instance()->HasMC();
  fNAddSteps=1;
  if (fHasMC){
    if (fStepsForSignal && fSignalsMC) fNAddSteps+=fSignalsMC->GetEntries();
    if (fStepsForBackground) ++fNAddSteps;
    if (fStepsForMCtruthOnly) --fNAddSteps; // No Step for Pair information
  } else {
    //if 
    fStepForMCtruth=kFALSE;
    fStepForNoCutsMCmotherPid=kFALSE;
    fStepsForSignal=kFALSE;
    fStepsForBackground=kFALSE;
  }
  // consitency checks to not duplicate steps
  if (fStepsForCutsIncreasing)     fStepForAfterAllCuts=kFALSE;
  if (fStepsForEachCut&&fNCuts==1) fStepForAfterAllCuts=kFALSE;
  
  fNSteps=0;
  if (fStepForMCtruth && fSignalsMC)           fNSteps+=fSignalsMC->GetEntries();
  if (fStepForNoCutsMCmotherPid && fSignalsMC) fNSteps+=fSignalsMC->GetEntries();
  if (fStepForAfterAllCuts)                    fNSteps+=fNAddSteps;
  
  if (fStepsForEachCut)         fNSteps+=(fNAddSteps*fNCuts); //one step for each cut + Signal (MC)
  if (fStepsForCutsIncreasing)  fNSteps+=(fNAddSteps*fNCuts); //one step for the increasing cuts + Signal (MC)
  // e.g. cut1, cut1&cut2, cut1&cut2&cut3, ...
  
  fNSteps+=(fNAddSteps*fNStepMasks);                              // cuts for the additional cut masks

  if (fStepForPreFilter) fNSteps+=fNAddSteps; //Add at the end for Prefilter (maxcutmask+1)
  
  // create the container
  
  Int_t *nbins=new Int_t[fNVars+2*fNVarsLeg];
  for (Int_t i=0;i<fNVars;++i) {
    Int_t nBins=(static_cast<TVectorD*>(fVarBinLimits->At(i)))->GetNrows()-1;
    nbins[i]=nBins;
  }
  for (Int_t i=0;i<fNVarsLeg;++i){
    Int_t nBins=(static_cast<TVectorD*>(fVarBinLimitsLeg->At(i)))->GetNrows()-1;
    nbins[i+fNVars]=nBins;
    nbins[i+fNVars+fNVarsLeg]=nBins;
  }
  
  fCfContainer = new AliCFContainer(GetName(), GetTitle(), fNSteps, fNVars+2*fNVarsLeg, nbins);
  delete [] nbins;
  
  // initialize the variables and their bin limits
  for (Int_t iVar=0; iVar<fNVars; iVar++) {
    UInt_t type=fVariables[iVar];
    Double_t *binLim = (static_cast<TVectorD*>(fVarBinLimits->At(iVar)))->GetMatrixArray();

    fCfContainer->SetBinLimits(iVar, binLim);
    fCfContainer->SetVarTitle(iVar, AliDielectronVarManager::GetValueName(type));
  }
  
  // initialize the variables and their bin limits for the Legs
  for (Int_t iVar=0; iVar<fNVarsLeg; iVar++) {
    UInt_t type=fVariablesLeg[iVar];
    Double_t *binLim=(static_cast<TVectorD*>(fVarBinLimitsLeg->At(iVar)))->GetMatrixArray();

    //Leg1
    fCfContainer->SetBinLimits(iVar+fNVars, binLim);
    fCfContainer->SetVarTitle(iVar+fNVars, Form("Leg1_%s",AliDielectronVarManager::GetValueName(type)));
    
    //Leg2
    fCfContainer->SetBinLimits(iVar+fNVars+fNVarsLeg, binLim);
    fCfContainer->SetVarTitle(iVar+fNVars+fNVarsLeg, Form("Leg2_%s",AliDielectronVarManager::GetValueName(type)));
  }

  // array for storing values
  fValues = new Double_t[fNVars+2*fNVarsLeg];

  // array for storing MC info
  if (fHasMC && fSignalsMC && fSignalsMC->GetEntries()>0) fIsMCTruth=new Bool_t[fSignalsMC->GetEntries()];
  //=================//
  // Set step titles //
  //=================//
  Int_t step=0;

  //Pure MC truth
  if(fStepForMCtruth && fSignalsMC) {
    for(Int_t i=0; i<fSignalsMC->GetEntries(); i++)
      fCfContainer->SetStepTitle(step++, Form("MC truth (Signal: %s)", fSignalsMC->At(i)->GetTitle()));
  }

  //before cuts (MC truth)
  if (fStepForNoCutsMCmotherPid && fSignalsMC){
    for(Int_t i=0; i<fSignalsMC->GetEntries(); i++)
      fCfContainer->SetStepTitle(step++,Form("No cuts (Signal: %s)",fSignalsMC->At(i)->GetTitle()));
  }
  
  TString cutName;
  //Steps for each of the cuts
  if (fStepsForEachCut){
    for (Int_t iCut=0; iCut<fNCuts;++iCut) {
      cutName=filter.GetCuts()->At(iCut)->GetName(); //TODO: User GetTitle???
      if (!fStepsForMCtruthOnly) {
        fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut
      }
      if (fHasMC){
        if (fStepsForSignal && fSignalsMC) {
          for(Int_t i=0; i<fSignalsMC->GetEntries(); i++) {
            fCfContainer->SetStepTitle(step++, Form("%s (Signal: %s)", cutName.Data(), fSignalsMC->At(i)->GetTitle())); //Step for  the cut with MC truth
          } 
        }
        if (fStepsForBackground)
          fCfContainer->SetStepTitle(step++, (cutName+" (Background)").Data()); //Step for the cut with MC truth
      }
    }
  }

  //Steps for increasing cut match
  if (fStepsForCutsIncreasing){
    cutName=""; //TODO: User GetTitle???
    for (Int_t iCut=0; iCut<fNCuts;++iCut) {
      if (!cutName.IsNull()) cutName+="&";
      cutName+=filter.GetCuts()->At(iCut)->GetName();
      if (!fStepsForMCtruthOnly) {
        fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut
      }
      if (fHasMC){
        if (fStepsForSignal && fSignalsMC)
          for(Int_t i=0; i<fSignalsMC->GetEntries(); i++) {
            fCfContainer->SetStepTitle(step++, Form("%s (Signal: %s)", cutName.Data(), fSignalsMC->At(i)->GetTitle())); //Step for the cut with MC truth
          }
          if (fStepsForBackground)
            fCfContainer->SetStepTitle(step++, (cutName+" (Background)").Data()); //Step for the cut with MC truth
      }
    }
  }

  //Steps of user defined cut combinations
  for (UInt_t iComb=0; iComb<fNStepMasks; ++iComb){
    cutName="";
    UInt_t mask=fStepMasks[iComb];
    for (Int_t iCut=0; iCut<fNCuts;++iCut) {
      if (mask&(1<<iCut)){
        if (cutName.IsNull()){
          cutName=filter.GetCuts()->At(iCut)->GetName();
        }else{
          cutName+="&";
          cutName+=filter.GetCuts()->At(iCut)->GetName();
        }
      }
    }
    
    fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut
    
    if (fHasMC){
      if (fStepsForSignal && fSignalsMC)
        for(Int_t i=0; i<fSignalsMC->GetEntries(); i++) {
          fCfContainer->SetStepTitle(step++, Form("%s (Signal: %s)", cutName.Data(), fSignalsMC->At(i)->GetTitle())); //Step for the cut with MC truth
        }
      if (fStepsForBackground)
        fCfContainer->SetStepTitle(step++, (cutName+" (Background)").Data()); //Step for the cut with MC truth
    }
  }

  //After All cuts
  if (fStepForAfterAllCuts){
    cutName="No pair cuts";
    if (filter.GetCuts()->At(0)){
      cutName=filter.GetCuts()->At(0)->GetName(); //TODO: User GetTitle???
      for (Int_t iCut=1; iCut<fNCuts;++iCut) {
        cutName+="&";
        cutName+=filter.GetCuts()->At(iCut)->GetName();
      }
    }
    if (!fStepsForMCtruthOnly) {
      fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut
    }
    if (fHasMC){
      if (fStepsForSignal && fSignalsMC)
        for(Int_t i=0; i<fSignalsMC->GetEntries(); i++) {
          fCfContainer->SetStepTitle(step++, Form("%s (Signal: %s)", cutName.Data(), fSignalsMC->At(i)->GetTitle())); //Step for the cut with MC truth
        }
        if (fStepsForBackground)
          fCfContainer->SetStepTitle(step++, (cutName+" (Background)").Data()); //Step for the cut with MC truth
    }
  }

  //Additional Step for result after PreFilter
  if (fStepForPreFilter){
    cutName="PreFilter";
    fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut
    if (fHasMC){
      if (fStepsForSignal && fSignalsMC)
        for(Int_t i=0; i<fSignalsMC->GetEntries(); i++) {
          fCfContainer->SetStepTitle(step++, Form("%s (Signal %s)", cutName.Data(), fSignalsMC->At(i)->GetTitle())); //Step for the cut with MC truth
        }
        if (fStepsForBackground)
          fCfContainer->SetStepTitle(step++, (cutName+" (Background)").Data()); //Step for the cut with MC truth
    }
  }



  if (step!=fNSteps) {
    AliError(Form("Something went wrong in the naming of the steps!!! (%d != %d)",step,fNSteps));
  }
}

//________________________________________________________________
void AliDielectronCF::Fill(UInt_t mask, const AliDielectronPair *particle)
{
  //
  // Fill the containers
  //

  // Check the MC truths
  if(fIsMCTruth) {
    for(Int_t i=0; i<fSignalsMC->GetEntries(); i++) fIsMCTruth[i]=kFALSE;
  }

  //TODO: for the moment don't fill truth information for mixed event paris. No valid MC info is available
  //      in the mixing handler
  Bool_t isMixedPair=(particle->GetType()>2&&particle->GetType()<10);
  
  Bool_t isBackground = kFALSE;
  if(fIsMCTruth && !isMixedPair) {
    for(Int_t i=0; i<fSignalsMC->GetEntries(); ++i) { 
      fIsMCTruth[i] = AliDielectronMC::Instance()->IsMCTruth(particle, (AliDielectronSignalMC*)fSignalsMC->At(i));
      isBackground = (isBackground || fIsMCTruth[i]);
    }
    // background is considered that pair which does not fulfill any of the signals
    isBackground = !isBackground;
  }
  
  Double_t valuesPair[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::SetFillMap(fUsedVars);
  AliDielectronVarManager::Fill(particle,valuesPair);

  // make MC truth information avaiable also in non-MCtruth steps
  AliDielectronMC* dieMC = AliDielectronMC::Instance();
  if(dieMC && dieMC->GetMCEvent()){
    AliDielectronVarManager::Fill(dieMC->GetMCEvent(), valuesPair);
  }
  
  for (Int_t iVar=0; iVar<fNVars; ++iVar){
    Int_t var=fVariables[iVar];
    fValues[iVar]=valuesPair[var];
  }

  if (fNVarsLeg>0){
    Double_t valuesLeg1[AliDielectronVarManager::kNMaxValues]={0};
    AliDielectronVarManager::Fill(particle->GetFirstDaughterP(),valuesLeg1);
    Double_t valuesLeg2[AliDielectronVarManager::kNMaxValues]={0};
    AliDielectronVarManager::Fill(particle->GetSecondDaughterP(),valuesLeg2);

    for (Int_t iVar=0; iVar<fNVarsLeg; ++iVar){
      Int_t var=fVariablesLeg[iVar];
      fValues[iVar+fNVars]=valuesLeg1[var];
      fValues[iVar+fNVars+fNVarsLeg]=valuesLeg2[var];
    }
  }
  
  UInt_t selectedMask=(1<<fNCuts)-1;

  //============//
  // Fill steps //
  //============//
  // Pure MC steps are handled in FillMC
  Int_t step=0;
  if (fStepForMCtruth && fIsMCTruth) step+=fSignalsMC->GetEntries();
  
  //No cuts (MC truth)
  if (fStepForNoCutsMCmotherPid && fIsMCTruth){
    for(Int_t i=0; i<fSignalsMC->GetEntries(); ++i) {
      if(fIsMCTruth[i]) {
        fCfContainer->Fill(fValues,step);
      }
      ++step;
    }
  }
  
  //Steps for each of the cuts
  if (fStepsForEachCut){
    for (Int_t iCut=0; iCut<fNCuts;++iCut) {
      UInt_t cutMask=1<<iCut;
      if ((mask&cutMask)==cutMask) {
        if(!fStepsForMCtruthOnly) {
          fCfContainer->Fill(fValues,step);
          ++step;
        }
        if (fHasMC){
          if ( fStepsForSignal && fIsMCTruth){
            for(Int_t i=0; i<fSignalsMC->GetEntries(); ++i) {
              if(fIsMCTruth[i]) {
                fCfContainer->Fill(fValues,step);
              }
              ++step;
            }
          }
          if ( fStepsForBackground ){
            if (isBackground) fCfContainer->Fill(fValues,step);
            ++step;
          }
        }
      } else {
        step+=fNAddSteps;
      }
    }
  }
  

  //Steps for increasing cut match
  if (fStepsForCutsIncreasing&&fNCuts>2){
    for (Int_t iCut=0; iCut<fNCuts;++iCut) {
      UInt_t cutMask=(1<<(iCut+1))-1;
      if ((mask&cutMask)==cutMask) {
        if(!fStepsForMCtruthOnly) {
          fCfContainer->Fill(fValues,step);
          ++step;
        }

        if (fHasMC){
          if ( fStepsForSignal && fIsMCTruth){
            for(Int_t i=0; i<fSignalsMC->GetEntries(); ++i) {
              if(fIsMCTruth[i]) {
                fCfContainer->Fill(fValues,step);
              }
              ++step;
            }
          }
          if ( fStepsForBackground ){
            if (isBackground) fCfContainer->Fill(fValues,step);
            ++step;
          }
        }
      } else {
        step+=fNAddSteps;
      }
    }
  }

  //Steps of user defined cut combinations
  for (UInt_t iComb=0; iComb<fNStepMasks; ++iComb){
    UInt_t userMask=fStepMasks[iComb];
    if ((mask&userMask)==userMask) {
      if(!fStepsForMCtruthOnly) {
        fCfContainer->Fill(fValues,step);
        ++step;
      }
      if (fHasMC){
        if ( fStepsForSignal && fIsMCTruth){
          for(Int_t i=0; i<fSignalsMC->GetEntries(); ++i) {
            if(fIsMCTruth[i]) {
              fCfContainer->Fill(fValues,step);
            }
            ++step;
          }
        }
        if ( fStepsForBackground ){
          if (isBackground) fCfContainer->Fill(fValues,step);
          ++step;
        }
      }
    } else {
      step+=fNAddSteps;
    }
  }
  
  //All cuts
  if (fStepForAfterAllCuts){
    if (mask == selectedMask){
      if(!fStepsForMCtruthOnly) {
        fCfContainer->Fill(fValues,step);
        ++step;
      }

      if (fHasMC){
        if ( fStepsForSignal && fIsMCTruth){
          for(Int_t i=0; i<fSignalsMC->GetEntries(); ++i) {
            if(fIsMCTruth[i]) {
              fCfContainer->Fill(fValues,step);
            }
            ++step;
          }
        }
        if ( fStepsForBackground ){
          if (isBackground) fCfContainer->Fill(fValues,step);
          ++step;
        }
      }
    } else {
      step+=fNAddSteps;
    }
  }

  //prefilter
  if (fStepForPreFilter) {
    if (mask&(1<<fNCuts)) {
      if(!fStepsForMCtruthOnly) {
        fCfContainer->Fill(fValues,step);
        ++step;
      }
      if (fHasMC){
        if ( fStepsForSignal && fIsMCTruth){
          for(Int_t i=0; i<fSignalsMC->GetEntries(); ++i) {
            if(fIsMCTruth[i]) {
              fCfContainer->Fill(fValues,step);
            }
            ++step;
          }
        }
        if ( fStepsForBackground ){
          if (isBackground) fCfContainer->Fill(fValues,step);
          ++step;
        }
      }
    }
    else {
      step+=fNAddSteps;
    }
  }
  
  if (step!=fNSteps) {
    AliError("Something went wrong in the step filling!!!");
  }
}

//________________________________________________________________
void AliDielectronCF::FillMC(const TObject *particle)
{
  //
  // fill MC part of the Container
  //
  if (!fStepForMCtruth) return;
  
  Double_t valuesPair[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::SetFillMap(fUsedVars);
  AliDielectronVarManager::Fill(particle,valuesPair);

  AliVParticle *d1=0x0;
  AliVParticle *d2=0x0;
  AliDielectronMC::Instance()->GetDaughters(particle,d1,d2);
  
  //TODO: temporary solution, set manually the pair type to 1: unlikesign SE
  valuesPair[AliDielectronVarManager::kPairType]=1;
  
  for (Int_t iVar=0; iVar<fNVars; ++iVar){
    Int_t var=fVariables[iVar];
    fValues[iVar]=valuesPair[var];
  }
  
  if (fNVarsLeg>0){
    Double_t valuesLeg1[AliDielectronVarManager::kNMaxValues];
    Double_t valuesLeg2[AliDielectronVarManager::kNMaxValues];
    if (d1->Pt()>d2->Pt()){
      AliDielectronVarManager::Fill(d1,valuesLeg1);
      AliDielectronVarManager::Fill(d2,valuesLeg2);
    } else {
      AliDielectronVarManager::Fill(d2,valuesLeg1);
      AliDielectronVarManager::Fill(d1,valuesLeg2);
    }
    
    for (Int_t iVar=0; iVar<fNVarsLeg; ++iVar){
      Int_t var=fVariablesLeg[iVar];
      fValues[iVar+fNVars]=valuesLeg1[var];
      fValues[iVar+fNVars+fNVarsLeg]=valuesLeg2[var];
    }
  }
  
  fCfContainer->Fill(fValues,0);
}


//________________________________________________________________
void AliDielectronCF::FillMC(Int_t label1, Int_t label2, Int_t nSignal) {
  //
  // fill the pure MC part of the container starting from a pair of 2 particles (part1 and part2 are legs)
  //
  if (!fStepForMCtruth) return;

  AliVParticle* part1 = AliDielectronMC::Instance()->GetMCTrackFromMCEvent(label1);
  AliVParticle* part2 = AliDielectronMC::Instance()->GetMCTrackFromMCEvent(label2);
  if(!part1 || !part2) return;

  AliDielectronMC* dieMC = AliDielectronMC::Instance();

  Int_t mLabel1 = dieMC->GetMothersLabel(label1);    // should work for both ESD and AOD
  Int_t mLabel2 = dieMC->GetMothersLabel(label2);
  // check the same mother option
  AliDielectronSignalMC* sigMC = (AliDielectronSignalMC*)fSignalsMC->At(nSignal);
  if(sigMC->GetMothersRelation()==AliDielectronSignalMC::kSame && mLabel1!=mLabel2) return;
  if(sigMC->GetMothersRelation()==AliDielectronSignalMC::kDifferent && mLabel1==mLabel2) return;

  AliDielectronVarManager::SetFillMap(fUsedVars);
  // fill the leg variables
  if (fNVarsLeg>0){
    Double_t valuesLeg1[AliDielectronVarManager::kNMaxValues];
    Double_t valuesLeg2[AliDielectronVarManager::kNMaxValues];
    if (part1->Pt()>part2->Pt()){
      AliDielectronVarManager::Fill(part1,valuesLeg1);
      AliDielectronVarManager::Fill(part2,valuesLeg2);
    } else {
      AliDielectronVarManager::Fill(part2,valuesLeg1);
      AliDielectronVarManager::Fill(part1,valuesLeg2);
    }
    
    for (Int_t iVar=0; iVar<fNVarsLeg; ++iVar){
      Int_t var=fVariablesLeg[iVar];
      fValues[iVar+fNVars]=valuesLeg1[var];
      fValues[iVar+fNVars+fNVarsLeg]=valuesLeg2[var];
    }
  }

  Double_t valuesPair[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::FillVarMCParticle2(part1,part2,valuesPair);
  AliDielectronVarManager::Fill(dieMC->GetMCEvent(), valuesPair);

  if(part1->Charge()*part2->Charge()<0)
    valuesPair[AliDielectronVarManager::kPairType]=1;
  else if(part1->Charge()>0)
    valuesPair[AliDielectronVarManager::kPairType]=0;
  else
    valuesPair[AliDielectronVarManager::kPairType]=2; // if one of the two particles is neutral, the pair will go here

  for(Int_t iVar=0; iVar<fNVars; ++iVar){
    Int_t var=fVariables[iVar];
    fValues[iVar]=valuesPair[var];
  }

  fCfContainer->Fill(fValues,nSignal);
}


//_____________________________________________________________________________
TVectorD* AliDielectronCF::MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax) const
{
  //
  // Make logarithmic binning
  // the user has to delete the array afterwards!!!
  //
  
  //check limits
  if (xmin<1e-20 || xmax<1e-20){
    AliError("For Log binning xmin and xmax must be > 1e-20. Using linear binning instead!");
    return MakeLinBinning(nbinsX, xmin, xmax);
  }
  if (xmax<xmin){
    Double_t tmp=xmin;
    xmin=xmax;
    xmax=tmp;
  }
  TVectorD *binLim=new TVectorD(nbinsX+1);
  Double_t first=xmin;
  Double_t last=xmax;
  Double_t expMax=TMath::Log(last/first);
  for (Int_t i=0; i<nbinsX+1; ++i){
    (*binLim)[i]=first*TMath::Exp(expMax/nbinsX*(Double_t)i);
  }
  return binLim;
}

//_____________________________________________________________________________
TVectorD* AliDielectronCF::MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax) const
{
  //
  // Make logarithmic binning
  // the user has to delete the array afterwards!!!
  //
  if (xmax<xmin){
    Double_t tmp=xmin;
    xmin=xmax;
    xmax=tmp;
  }
  TVectorD *binLim=new TVectorD(nbinsX+1);
  Double_t first=xmin;
  Double_t last=xmax;
  Double_t binWidth=(last-first)/nbinsX;
  for (Int_t i=0; i<nbinsX+1; ++i){
    (*binLim)[i]=first+binWidth*(Double_t)i;
  }
  return binLim;
}

