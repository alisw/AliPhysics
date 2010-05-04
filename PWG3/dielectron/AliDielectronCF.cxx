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
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliDielectronCF.h"
#include "AliDielectronMC.h"
#include "AliDielectronPair.h"

ClassImp(AliDielectronCF)

AliDielectronCF::AliDielectronCF() :
  TNamed("DielectronCF","DielectronCF"),
  fNSteps(0),
  fNVars(0),
  fNVarsLeg(0),
  fNCuts(0),
  fValues(0x0),
  fStepForMCtruth(kFALSE),
  fStepForNoCutsMCmotherPid(kFALSE),
  fStepForAfterAllCuts(kTRUE),
  fStepsForEachCut(kFALSE),
  fStepsForCutsIncreasing(kFALSE),
  fStepsForSignal(kTRUE),
  fStepsForBackground(kFALSE),
  fNStepMasks(0),
  fPdgMother(-1),
  fCfContainer(0x0),
  fHasMC(kFALSE),
  fNAddSteps(0)
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
    fNBinsLeg[i]=0;
    fVarLoLimitLeg[i]=0.;
    fVarUpLimitLeg[i]=0.;
    fStepMasks[i]=0xFFFFFF;
  }
}

//________________________________________________________________
AliDielectronCF::AliDielectronCF(const char* name, const char* title) :
  TNamed(name, title),
  fNSteps(0),
  fNVars(0),
  fNVarsLeg(0),
  fNCuts(0),
  fValues(0x0),
  fStepForMCtruth(kFALSE),
  fStepForNoCutsMCmotherPid(kFALSE),
  fStepForAfterAllCuts(kTRUE),
  fStepsForEachCut(kFALSE),
  fStepsForCutsIncreasing(kFALSE),
  fStepsForSignal(kTRUE),
  fStepsForBackground(kFALSE),
  fNStepMasks(0),
  fPdgMother(-1),
  fCfContainer(0x0),
  fHasMC(kFALSE),
  fNAddSteps(0)
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
    fNBinsLeg[i]=0;
    fVarLoLimitLeg[i]=0.;
    fVarUpLimitLeg[i]=0.;
    fStepMasks[i]=0xFFFFFF;
  }
}

//________________________________________________________________
AliDielectronCF::~AliDielectronCF()
{
  //
  // Destructor
  //
  if (fValues) delete [] fValues;
}

//________________________________________________________________
void AliDielectronCF::AddVariable(AliDielectronVarManager::ValueTypes type, Int_t nbins, Double_t min, Double_t max, Bool_t leg)
{
  //
  // Add a variable to the CF configuration
  //

  if (!leg){
    fVariables[fNVars]  = (UInt_t)type;
    fVarLoLimit[fNVars] = min;
    fVarUpLimit[fNVars] = max;
    fNBins[fNVars]      = nbins;
    ++fNVars;
  } else {
    fVariablesLeg[fNVarsLeg]  = (UInt_t)type;
    fVarLoLimitLeg[fNVarsLeg] = min;
    fVarUpLimitLeg[fNVarsLeg] = max;
    fNBinsLeg[fNVarsLeg]      = nbins;
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
    if (fStepsForSignal) ++fNAddSteps;
    if (fStepsForBackground) ++fNAddSteps;
  } else {
    //if 
    fStepForMCtruth=kFALSE;
    fStepForNoCutsMCmotherPid=kFALSE;
    fStepsForSignal=kFALSE;
    fStepsForBackground=kFALSE;
  }
    
  fNSteps=0;
  if (fStepForMCtruth) ++fNSteps;
  if (fStepForNoCutsMCmotherPid) ++fNSteps;
  if (fStepForAfterAllCuts) fNSteps+=fNAddSteps;

  if (fStepsForEachCut&&fNCuts>1)        fNSteps+=(fNAddSteps*fNCuts);     //one step for each cut + Signal (MC)
  if (fStepsForCutsIncreasing&&fNCuts>2) fNSteps+=(fNAddSteps*(fNCuts-2)); //one step for the increasing cuts + Signal (MC)
                                                      // e.g. cut2&cut3, cut2&cut3&cut4
  fNSteps+=(fNAddSteps*fNStepMasks);                            // cuts for the additional cut masks
  // create the container
  Int_t *nbins=new Int_t[fNVars+2*fNVarsLeg];
  for (Int_t i=0;i<fNVars;++i)    nbins[i]=fNBins[i];
  for (Int_t i=0;i<fNVarsLeg;++i) nbins[i+fNVars]=fNBinsLeg[i];
  for (Int_t i=0;i<fNVarsLeg;++i) nbins[i+fNVars+fNVarsLeg]=fNBinsLeg[i];
  
  fCfContainer = new AliCFContainer(GetName(), GetTitle(), fNSteps, fNVars+2*fNVarsLeg, nbins);
  delete [] nbins;
  
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
    delete [] binLim;
  }
  
  // initialize the variables and their bin limits for the Legs
  for (Int_t iVar=0; iVar<fNVarsLeg; iVar++) {
    UInt_t type=fVariablesLeg[iVar];
    Int_t    nBins = fNBinsLeg[iVar];
    Double_t loLim = fVarLoLimitLeg[iVar];
    Double_t upLim = fVarUpLimitLeg[iVar];
    Double_t *binLim = new Double_t[nBins+1];
    
    // set the bin limits
    for(Int_t iBin=0; iBin<=nBins; iBin++) binLim[iBin] = loLim + (upLim-loLim) / nBins*(Double_t)iBin;

    //Leg1
    fCfContainer->SetBinLimits(iVar+fNVars, binLim);
    fCfContainer->SetVarTitle(iVar+fNVars, Form("Leg1_%s",AliDielectronVarManager::GetValueName(type)));
    
    //Leg2
    fCfContainer->SetBinLimits(iVar+fNVars+fNVarsLeg, binLim);
    fCfContainer->SetVarTitle(iVar+fNVars+fNVarsLeg, Form("Leg2_%s",AliDielectronVarManager::GetValueName(type)));
    
    delete [] binLim;
  }

  // array for storing values
  fValues = new Double_t[fNVars+2*fNVarsLeg];
  
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
    fCfContainer->SetStepTitle(step++,"No cuts (Signal)");
  }
  
  TString cutName;
  //Steps for each of the cuts
  if (fStepsForEachCut&&fNCuts>1){
    for (Int_t iCut=0; iCut<fNCuts;++iCut) {
      cutName=filter.GetCuts()->At(iCut)->GetName(); //TODO: User GetTitle???
      
      fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut
      
      if (fHasMC){
        if (fStepsForSignal)
          fCfContainer->SetStepTitle(step++, (cutName+" (Signal)").Data()); //Step for the cut with MC truth
        if (fStepsForBackground)
          fCfContainer->SetStepTitle(step++, (cutName+" (Background)").Data()); //Step for the cut with MC truth
      }
    }
  }

  //Steps for increasing cut match
  if (fStepsForCutsIncreasing&&fNCuts>2){
    cutName=filter.GetCuts()->At(0)->GetName(); //TODO: User GetTitle???
    for (Int_t iCut=1; iCut<fNCuts-1;++iCut) {
      cutName+="&";
      cutName+=filter.GetCuts()->At(iCut)->GetName();
      
      fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut
      
      if (fHasMC){
        if (fStepsForSignal)
          fCfContainer->SetStepTitle(step++, (cutName+" (Signal)").Data()); //Step for the cut with MC truth
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
      if (fStepsForSignal)
        fCfContainer->SetStepTitle(step++, (cutName+" (Signal)").Data()); //Step for the cut with MC truth
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
      fCfContainer->SetStepTitle(step++, cutName.Data()); //Step for the cut
    }
    if (fHasMC){
      if (fStepsForSignal)
        fCfContainer->SetStepTitle(step++, (cutName+" (Signal)").Data()); //Step for the cut with MC truth
      if (fStepsForBackground)
        fCfContainer->SetStepTitle(step++, (cutName+" (Background)").Data()); //Step for the cut with MC truth
    }
  }

  if (step!=fNSteps) {
    AliError("Something went wrong in the naming of the steps!!!");
  }
}

//________________________________________________________________
void AliDielectronCF::Fill(UInt_t mask, const AliDielectronPair *particle)
{
  //
  // Fill the containers
  //

  Bool_t isMCTruth=kFALSE;
  if (fHasMC) isMCTruth=AliDielectronMC::Instance()->IsMotherPdg(particle,fPdgMother);

  Double_t valuesPair[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::Fill(particle,valuesPair);

  for (Int_t iVar=0; iVar<fNVars; ++iVar){
    Int_t var=fVariables[iVar];
    fValues[iVar]=valuesPair[var];
  }

  if (fNVarsLeg>0){
    Double_t valuesLeg1[AliDielectronVarManager::kNMaxValues];
    AliDielectronVarManager::Fill(particle->GetFirstDaughter(),valuesLeg1);
    Double_t valuesLeg2[AliDielectronVarManager::kNMaxValues];
    AliDielectronVarManager::Fill(particle->GetSecondDaughter(),valuesLeg2);

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
  // step 0 would be full MC truth and is handled in FillMC
  Int_t step=0;
  if (fStepForMCtruth) ++step;

  //No cuts (MC truth)
  if (fStepForNoCutsMCmotherPid){
    if (isMCTruth) fCfContainer->Fill(fValues,step);
    ++step;
  }
  
  //Steps for each of the cuts
  if (fStepsForEachCut&&fNCuts>1){
    for (Int_t iCut=0; iCut<fNCuts;++iCut) {
      if (mask&(1<<iCut)) {
        fCfContainer->Fill(fValues,step);
        ++step;

        if (fHasMC){
          if ( fStepsForSignal){
            if (isMCTruth) fCfContainer->Fill(fValues,step);
            ++step;
          }
          if ( fStepsForBackground ){
            if (!isMCTruth) fCfContainer->Fill(fValues,step);
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
    for (Int_t iCut=1; iCut<fNCuts-1;++iCut) {
      if (mask&(1<<((iCut+1)-1))) {
        fCfContainer->Fill(fValues,step);
        ++step;
        
        if (fHasMC){
          if ( fStepsForSignal){
            if (isMCTruth) fCfContainer->Fill(fValues,step);
            ++step;
          }
          if ( fStepsForBackground ){
            if (!isMCTruth) fCfContainer->Fill(fValues,step);
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
    if (mask&userMask) {
      fCfContainer->Fill(fValues,step);
      ++step;
      
      if (fHasMC){
        if ( fStepsForSignal){
          if (isMCTruth) fCfContainer->Fill(fValues,step);
          ++step;
        }
        if ( fStepsForBackground ){
          if (!isMCTruth) fCfContainer->Fill(fValues,step);
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
      fCfContainer->Fill(fValues,step);
      ++step;
      
      if (fHasMC){
        if ( fStepsForSignal){
          if (isMCTruth) fCfContainer->Fill(fValues,step);
          ++step;
        }
        if ( fStepsForBackground ){
          if (!isMCTruth) fCfContainer->Fill(fValues,step);
          ++step;
        }
      }
    } else {
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
  AliDielectronVarManager::Fill(particle,valuesPair);
  
  //TODO: temporary solution, set manually the pair type to 1: unlikesign SE
  valuesPair[AliDielectronVarManager::kPairType]=1;
  
  for (Int_t iVar=0; iVar<fNVars; ++iVar){
    Int_t var=fVariables[iVar];
    fValues[iVar]=valuesPair[var];
  }
  
  if (fNVarsLeg>0){
    AliVParticle *d1=0x0;
    AliVParticle *d2=0x0;
    AliDielectronMC::Instance()->GetDaughters(particle,d1,d2);
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

