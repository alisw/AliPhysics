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
#include "AliDielectronMC.h"

ClassImp(AliDielectronVarCuts)


AliDielectronVarCuts::AliDielectronVarCuts() :
  AliAnalysisCuts(),
  fNActiveCuts(0),
  fActiveCutsMask(0),
  fSelectedCutsMask(0),
  fCutOnMCtruth(kFALSE),
  fCutType(kAll)
{
  //
  // Default costructor
  //
  for (Int_t i=0; i<AliDielectronVarManager::kNMaxValues; ++i){
    fActiveCuts[i]=0;
    fCutMin[i]=0;
    fCutMax[i]=0;
    fCutExclude[i]=kFALSE;
  }
}

//________________________________________________________________________
AliDielectronVarCuts::AliDielectronVarCuts(const char* name, const char* title) :
  AliAnalysisCuts(name,title),
  fNActiveCuts(0),
  fActiveCutsMask(0),
  fSelectedCutsMask(0),
  fCutOnMCtruth(kFALSE),
  fCutType(kAll)
{
  //
  // Named contructor
  //
  for (Int_t i=0; i<AliDielectronVarManager::kNMaxValues; ++i){
    fActiveCuts[i]=0;
    fCutMin[i]=0;
    fCutMax[i]=0;
    fCutExclude[i]=kFALSE;
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

  //reset
  fSelectedCutsMask=0;
  SetSelected(kFALSE);

  if (!track) return kFALSE;
  
  //If MC cut, get MC truth
  if (fCutOnMCtruth){
    AliVParticle *part=static_cast<AliVParticle*>(track);
    track=AliDielectronMC::Instance()->GetMCTrackFromMCEvent(part->GetLabel());
    if (!track) return kFALSE;
  }

  //Fill values
  Double_t values[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::Fill(track,values);
  
  for (Int_t iCut=0; iCut<fNActiveCuts; ++iCut){
    Int_t cut=fActiveCuts[iCut];
    SETBIT(fSelectedCutsMask,iCut);
    if ( ((values[cut]<fCutMin[iCut]) || (values[cut]>fCutMax[iCut]))^fCutExclude[iCut] ) CLRBIT(fSelectedCutsMask,iCut);
    if ( fCutType==kAll && !TESTBIT(fSelectedCutsMask,iCut) ) return kFALSE; // option to (minor) speed improvement
  }
  
  Bool_t isSelected=(fSelectedCutsMask==fActiveCutsMask);
  if ( fCutType==kAny ) isSelected=(fSelectedCutsMask>0);
  SetSelected(isSelected);
  return isSelected;
}

//________________________________________________________________________
void AliDielectronVarCuts::AddCut(AliDielectronVarManager::ValueTypes type, Double_t min, Double_t max, Bool_t excludeRange)
{
  //
  // Set cut range and activate it
  //
  if (min>max){
    Double_t tmp=min;
    min=max;
    max=tmp;
  }
  fCutMin[fNActiveCuts]=min;
  fCutMax[fNActiveCuts]=max;
  fCutExclude[fNActiveCuts]=excludeRange;
  SETBIT(fActiveCutsMask,fNActiveCuts);
  fActiveCuts[fNActiveCuts]=(UShort_t)type;
  ++fNActiveCuts;
}

//________________________________________________________________________
void AliDielectronVarCuts::Print(const Option_t* /*option*/) const
{
  //
  // Print cuts and the range
  //
  printf("cut ranges for '%s'\n",GetTitle());
  if (fCutType==kAll){
    printf("All Cuts have to be fulfilled\n");
  } else {
    printf("Any Cut has to be fulfilled\n");
  }
  for (Int_t iCut=0; iCut<fNActiveCuts; ++iCut){
    Int_t cut=(Int_t)fActiveCuts[iCut];
    Bool_t inverse=fCutExclude[iCut];

    if (!inverse){
      printf("Cut %02d: %f < %s < %f\n", iCut,
             fCutMin[iCut], AliDielectronVarManager::GetValueName((Int_t)cut), fCutMax[iCut]);
    } else {
      printf("Cut %02d: !(%f < %s < %f)\n", iCut,
             fCutMin[iCut], AliDielectronVarManager::GetValueName((Int_t)cut), fCutMax[iCut]);
    }
  }
}
