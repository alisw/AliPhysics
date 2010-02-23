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
//   Cut class providing cuts to all infomation                          //
//     available for the AliVParticle interface                          //                                                     //
//                                                                       //
// Authors:                                                              //
//   Jens Wiechula <Jens.Wiechula@cern.ch>                               //
/*



*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#include "AliDielectronVarCuts.h"


ClassImp(AliDielectronVarCuts)


AliDielectronVarCuts::AliDielectronVarCuts() :
  AliAnalysisCuts(),
  fNActiveCuts(0),
  fActiveCutsMask(0),
  fSelectedCutsMask(0)
{
  //
  // Default costructor
  //
  for (Int_t i=0; i<AliDielectronVarManager::kNMaxValues; ++i){
    fActiveCuts[i]=0;
    fCutMin[i]=0;
    fCutMax[i]=0;
  }
}

//________________________________________________________________________
AliDielectronVarCuts::AliDielectronVarCuts(const char* name, const char* title) :
  AliAnalysisCuts(name,title),
  fNActiveCuts(0),
  fActiveCutsMask(0),
  fSelectedCutsMask(0)
{
  //
  // Named contructor
  //
  for (Int_t i=0; i<AliDielectronVarManager::kNMaxValues; ++i){
    fActiveCuts[i]=0;
    fCutMin[i]=0;
    fCutMax[i]=0;
  }
}

//________________________________________________________________________
AliDielectronVarCuts::~AliDielectronVarCuts()
{
  //
  // Destructor
  //
}

//________________________________________________________________________
Bool_t AliDielectronVarCuts::IsSelected(TObject* track)
{
  //
  // Make cut decision
  //

  Double_t values[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::Fill(track,values);
  fSelectedCutsMask=0;
  SetSelected(kFALSE);
  
  for (Int_t iCut=0; iCut<fNActiveCuts; ++iCut){
    Int_t cut=fActiveCuts[iCut];
    SETBIT(fSelectedCutsMask,iCut);
    if ( (values[cut]<fCutMin[cut]) || (values[cut]>fCutMax[cut]) ) CLRBIT(fSelectedCutsMask,iCut);
  }
  Bool_t isSelected=(fSelectedCutsMask==fActiveCutsMask);
  SetSelected(isSelected);
  return isSelected;
}

//________________________________________________________________________
void AliDielectronVarCuts::AddCut(Double_t min, Double_t max, AliDielectronVarManager::ValueTypes type)
{
  //
  // Set cut range and activate it
  //
  if (min>max){
    Double_t tmp=min;
    min=max;
    max=tmp;
  }
  fCutMin[type]=min;
  fCutMax[type]=max;
  ActivateCut(type);
}

//________________________________________________________________________
void AliDielectronVarCuts::ActivateCut(AliDielectronVarManager::ValueTypes cutName)
{
  //
  // Add the cut to the list of active cuts
  //

  if (IsCutActive(cutName)) return;
  SETBIT(fActiveCutsMask,fNActiveCuts);
  fActiveCuts[fNActiveCuts++]=(UChar_t)cutName;
}

//________________________________________________________________________
Bool_t AliDielectronVarCuts::IsCutActive(AliDielectronVarManager::ValueTypes cut)
{
  //
  // Check if this cut is already activated
  //
  for (Int_t iCut=0; iCut<fNActiveCuts; ++iCut){
    if (fActiveCuts[iCut]==(UChar_t)cut) return kTRUE;
  }
  
  return kFALSE;
}

//________________________________________________________________________
void AliDielectronVarCuts::Print(const Option_t* /*option*/) const
{
  //
  // Print cuts and the range
  //
  printf("cut ranges for '%s'\n",GetTitle());
  for (Int_t iCut=0; iCut<fNActiveCuts; ++iCut){
    UChar_t cut=fActiveCuts[iCut];
    printf("Cut %02d: %f < %s < %f\n", iCut,
           fCutMin[cut], AliDielectronVarManager::GetValueName((Int_t)cut), fCutMax[cut]);
  }
}
