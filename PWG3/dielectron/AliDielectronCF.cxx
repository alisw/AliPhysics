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

#include <AliCFContainer.h>
#include <AliAnalysisFilter.h>
#include <AliAnalysisCuts.h>
#include <AliLog.h>

#include "AliDielectronCF.h"
#include "AliDielectronMC.h"

ClassImp(AliDielectronCF)

AliDielectronCF::AliDielectronCF() :
  TNamed("DielectronCF","DielectronCF"),
  fNSteps(0),
  fNVars(0),
  fNCuts(0),
  fStepForMCtruth(kFALSE),
  fStepForNoCutsMCmotherPid(kFALSE),
  fStepForAfterAllCuts(kTRUE),
  fStepsForEachCut(kFALSE),
  fStepsForCutsIncreasing(kFALSE),
  fNStepMasks(0),
  fPdgMother(-1),
  fCfContainer(0x0)
{
  //
  // Default constructor
  //
  for (Int_t i=0; i<AliDielectronVarManager::kNMaxValues; ++i){
    fVariables[i]=0;
  }

  for (Int_t i=0; i<kNmaxAddSteps; ++i){
    fNBins[i]=0;
    fVarLoLimit[i]=0.;
    fVarUpLimit[i]=0.;
    fStepMasks[i]=0xFFFFFF;
  }
}

//________________________________________________________________
AliDielectronCF::AliDielectronCF(const char* name, const char* title) :
  TNamed(name, title),
  fNSteps(0),
  fNVars(0),
  fNCuts(0),
  fStepForMCtruth(kFALSE),
  fStepForNoCutsMCmotherPid(kFALSE),
  fStepForAfterAllCuts(kTRUE),
  fStepsForEachCut(kFALSE),
  fStepsForCutsIncreasing(kFALSE),
  fNStepMasks(0),
  fPdgMother(-1),
  fCfContainer(0x0)
{
  //
  // Named constructor
  //
  for (Int_t i=0; i<AliDielectronVarManager::kNMaxValues; ++i){
    fVariables[i]=0;
  }
  
  for (Int_t i=0; i<kNmaxAddSteps; ++i){
    fNBins[i]=0;
    fVarLoLimit[i]=0.;
    fVarUpLimit[i]=0.;
    fStepMasks[i]=0xFFFFFF;
  }
}

//________________________________________________________________
AliDielectronCF::~AliDielectronCF()
{
  //
  // Destructor
  //
  
}

//________________________________________________________________
void AliDielectronCF::AddVariable(AliDielectronVarManager::ValueTypes type, Int_t nbins, Double_t min, Double_t max)
{
  //
  // Add a variable to the CF configuration
  //
  fVariables[fNVars]  = (UInt_t)type;
  fVarLoLimit[fNVars] = min;
  fVarUpLimit[fNVars] = max;
  fNBins[fNVars]      = nbins;
  ++fNVars;
}

//________________________________________________________________
void AliDielectronCF::InitialiseContainer(const AliAnalysisFilter& filter)
{
  //
  // Initialise container based on the cuts in the analysis filter
  //

  fNCuts=filter.GetCuts()->GetEntries();

  fNSteps=0;
  if (fStepForMCtruth) ++fNSteps;
  if (fStepForNoCutsMCmotherPid) ++fNSteps;
  if (fStepForAfterAllCuts) fNSteps+=2;
  
  if (fStepsForEachCut&&fNCuts>1)        fNSteps+=(2*fNCuts);     //one step for each cut + MC truth
  if (fStepsForCutsIncreasing&&fNCuts>2) fNSteps+=(2*(fNCuts-2)); //one step for the increasing cuts + MC truth
                                                      // e.g. cut2&cut3, cut2&cut3&cut4
  fNSteps+=(2*fNStepMasks);                            // cuts for the additional cut masks
  // create the container
  fCfContainer = new AliCFContainer(GetName(), GetTitle(), fNSteps, fNVars, fNBins);
  
  // initialize the variables and their bin limits
  for (Int_t iVar=0; iVar<fNVars; iVar++) {
    UInt_t type=fVariables[iVar];
    Int_t    nBins = fNBins[iVar];
    Double_t loLim = fVarLoLimit[iVar];
    Double_t upLim = fVarUpLimit[iVar];
    Double_t *binLim = new Double_t[nBins+1];
    
    // set the bin limits
    for(Int_t iBin=0; iBin<=nBins; iBin++) binLim[iBin] = loLim + (upLim-loLim) / nBins*(Double_t)iBin;
    
    fCfContainer->SetBinLimits(iVar, binLim);
    fCfContainer->SetVarTitle(iVar, AliDielectronVarManager::GetValueName(type));
    delete binLim;
  }
  
  //=================//
  // Set step titles //
  //=================//
  Int_t step=0;

  //Pure MC truth
  if (fStepForMCtruth){
      fCfContainer->SetStepTitle(step++,"MC truth");
  }

  //before cuts (MC truth)
  if (fStepForNoCutsMCmotherPid){
    fCfContainer->SetStepTitle(step++,"No cuts (MC mother)");
  }
  
  //After All cuts
  TString cutName;
  if (fStepForAfterAllCuts){
    cutName="All Cuts"; //TODO: User GetTitle???
    fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut
    cutName+=" (MC truth)";
    fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut with MC truth
  }

  //Steps for each of the cuts
  if (fStepsForEachCut&&fNCuts>1){
    for (Int_t iCut=0; iCut<fNCuts;++iCut) {
      cutName=filter.GetCuts()->At(iCut)->GetName(); //TODO: User GetTitle???
      fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut
      cutName+=" (MC mother)";
      fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut with MC truth
    }
  }

  //Steps for increasing cut match
  if (fStepsForCutsIncreasing&&fNCuts>2){
    cutName=filter.GetCuts()->At(0)->GetName(); //TODO: User GetTitle???
    for (Int_t iCut=1; iCut<fNCuts-1;++iCut) {
      cutName+="&";
      cutName+=filter.GetCuts()->At(iCut)->GetName();
      fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut
      cutName+=" (MC mother)";
      fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut with MC truth
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
    cutName+=" (MC mother)";
    fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut with MC truth
  }

  if (step!=fNSteps) {
    AliError("Something went wrong in the naming of the steps!!!");
  }
}

//________________________________________________________________
void AliDielectronCF::Fill(UInt_t mask, const TObject *particle)
{
  //
  // Fill the containers
  //

  Bool_t isMCTruth=kFALSE;
  if (fPdgMother>=0) isMCTruth=AliDielectronMC::Instance()->IsMotherPdg(particle,fPdgMother);

  Double_t valuesAll[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::Fill(particle,valuesAll);

  Double_t values[AliDielectronVarManager::kNMaxValues];
  for (Int_t iVar=0; iVar<fNVars; ++iVar){
    Int_t var=fVariables[iVar];
    values[iVar]=valuesAll[var];
  }

  UInt_t selectedMask=(1<<fNCuts)-1;

  //============//
  // Fill steps //
  //============//
  // step 0 would be full MC truth and is handled in FillMC
  Int_t step=0;
  if (fStepForMCtruth) ++step;

  //No cuts (MC truth)
  if (fStepForNoCutsMCmotherPid){
    if (isMCTruth) fCfContainer->Fill(values,step);
    ++step;
  }
  
  //All cuts
  if (fStepForAfterAllCuts){
    if (mask == selectedMask){
      fCfContainer->Fill(values,step);
      ++step;
      if (isMCTruth) fCfContainer->Fill(values,step);
      ++step;
    } else {
      step+=2;
    }
  }
  
  //Steps for each of the cuts
  if (fStepsForEachCut&&fNCuts>1){
    for (Int_t iCut=0; iCut<fNCuts;++iCut) {
      if (mask&(1<<iCut)) {
        fCfContainer->Fill(values,step);
        ++step;
        if (isMCTruth) fCfContainer->Fill(values,step);
        ++step;
      } else {
        step+=2;
      }
    }
  }

  //Steps for increasing cut match
  if (fStepsForCutsIncreasing&&fNCuts>2){
    for (Int_t iCut=1; iCut<fNCuts-1;++iCut) {
      if (mask&(1<<((iCut+1)-1))) {
        fCfContainer->Fill(values,step);
        ++step;
        if (isMCTruth) fCfContainer->Fill(values,step);
        ++step;
      } else {
        step+=2;
      }
    }
  }

  //Steps of user defined cut combinations
  for (UInt_t iComb=0; iComb<fNStepMasks; ++iComb){
    UInt_t userMask=fStepMasks[iComb];
    if (mask&userMask) {
      fCfContainer->Fill(values,step);
      ++step;
      if (isMCTruth) fCfContainer->Fill(values,step);
      ++step;
    } else {
      step+=2;
    }
  }
  
}

//________________________________________________________________
void AliDielectronCF::FillMC(const TObject *particle)
{
  //
  // fill MC part of the Container
  //
  if (!fStepForMCtruth) return;
  
  Double_t valuesAll[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::Fill(particle,valuesAll);
  
  Double_t values[AliDielectronVarManager::kNMaxValues];
  for (Int_t iVar=0; iVar<fNVars; ++iVar){
    Int_t var=fVariables[iVar];
    values[iVar]=valuesAll[var];
  }
  //TODO: temporary solution, set manually the pair type to 1: mixed e+-
  values[AliDielectronVarManager::kPairType]=1;
  fCfContainer->Fill(values,0);
}

