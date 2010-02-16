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
    SETBIT(fSelectedCutsMask,cut);
    if ( (values[cut]<fCutMin[cut]) || (values[cut]>fCutMax[cut]) ) CLRBIT(fSelectedCutsMask,cut);
  }
  Bool_t isSelected=(fSelectedCutsMask==fActiveCutsMask);
  SetSelected(isSelected);
  return isSelected;
}
//________________________________________________________________________
void AliDielectronVarCuts::ActivateCut(AliDielectronVarManager::ValueTypes cutName)
{
  //
  // Add the cut to the list of active cuts
  //

  if (IsCutActive(cutName)) return;
  fActiveCuts[fNActiveCuts++]=(UChar_t)cutName;
  SETBIT(fActiveCutsMask,cutName);
}


