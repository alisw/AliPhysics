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


#include <TH1.h>

#include "AliDielectronVarCuts.h"
#include "AliDielectronMC.h"

ClassImp(AliDielectronVarCuts)


AliDielectronVarCuts::AliDielectronVarCuts() :
  AliAnalysisCuts(),
  fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
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
    fUpperCut[i]=0x0;
  }
}

//________________________________________________________________________
AliDielectronVarCuts::AliDielectronVarCuts(const char* name, const char* title) :
  AliAnalysisCuts(name,title),
  fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
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
    fUpperCut[i]=0x0;
  }
}

//________________________________________________________________________
AliDielectronVarCuts::~AliDielectronVarCuts()
{
  //
  // Destructor
  //
  if (fUsedVars) delete fUsedVars;
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
  AliDielectronVarManager::SetFillMap(fUsedVars);
  AliDielectronVarManager::Fill(track,values);

  for (Int_t iCut=0; iCut<fNActiveCuts; ++iCut){
    Int_t cut=fActiveCuts[iCut];
    SETBIT(fSelectedCutsMask,iCut);
    if ( !fUpperCut[iCut] && ((values[cut]<fCutMin[iCut]) || (values[cut]>fCutMax[iCut]))^fCutExclude[iCut] ) CLRBIT(fSelectedCutsMask,iCut);
    else if ( fUpperCut[iCut]) {
      // use a TH1 inherited cut object
      Float_t x=0.,y=0.,z=0.;
      switch (fUpperCut[iCut]->GetDimension()) {
      case 3: z=values[fUpperCut[iCut]->GetZaxis()->GetUniqueID()];
      case 2: y=values[fUpperCut[iCut]->GetYaxis()->GetUniqueID()];
      case 1: x=values[fUpperCut[iCut]->GetXaxis()->GetUniqueID()];
      }
      Int_t bin = fUpperCut[iCut]->FindBin(x,y,z);
      if ( ((values[cut]<fCutMin[iCut]) || (values[cut]>fUpperCut[iCut]->GetBinContent(bin)))^fCutExclude[iCut] ) CLRBIT(fSelectedCutsMask,iCut);
    }
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
  fUsedVars->SetBitNumber(type,kTRUE);
  ++fNActiveCuts;
}

//________________________________________________________________________
void AliDielectronVarCuts::AddCut(AliDielectronVarManager::ValueTypes type, Double_t min, TH1 * const max,  Bool_t excludeRange)
{
  //
  // Set cut range and activate it
  //
  fCutMin[fNActiveCuts]=min;
  fCutMax[fNActiveCuts]=0.0;
  fCutExclude[fNActiveCuts]=excludeRange;
  SETBIT(fActiveCutsMask,fNActiveCuts);
  fActiveCuts[fNActiveCuts]=(UShort_t)type;
  fUsedVars->SetBitNumber(type,kTRUE);
  // cut dependencies
  UInt_t var =0;
  switch(max->GetDimension()) {
  case 3:
  case 2:
    var=AliDielectronVarManager::GetValueType(max->GetZaxis()->GetName());
    fUsedVars->SetBitNumber(var,kTRUE);
    max->GetZaxis()->SetUniqueID(var);
  case 1:
    var=AliDielectronVarManager::GetValueType(max->GetYaxis()->GetName());
    fUsedVars->SetBitNumber(var,kTRUE);
    max->GetYaxis()->SetUniqueID(var);
  /*case 1:*/
    var=AliDielectronVarManager::GetValueType(max->GetXaxis()->GetName());
    fUsedVars->SetBitNumber(var,kTRUE);
    max->GetXaxis()->SetUniqueID(var);
  }
  fUpperCut[fNActiveCuts]=(TH1*)max->Clone("histCut");
  fUpperCut[fNActiveCuts]->SetDirectory(0x0);
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
