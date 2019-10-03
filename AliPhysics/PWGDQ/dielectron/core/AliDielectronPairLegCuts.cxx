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
//   Cut class providing cuts for both legs in the AliDielectronPair     //
//                                                                       //
//                                                                       //
/*
Add any number of leg cuts using e.g. for leg 1
GetFilterLeg1().AddCuts(mycut)
where mycut has to inherit from AliAnalysisCuts

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TList.h>

#include "AliDielectronPair.h"
#include "AliVParticle.h"

#include "AliDielectronPairLegCuts.h"

ClassImp(AliDielectronPairLegCuts)


AliDielectronPairLegCuts::AliDielectronPairLegCuts() :
  AliAnalysisCuts(),
  fFilterLeg1("PairFilterLeg1","PairFilterLeg1"),
  fFilterLeg2("PairFilterLeg2","PairFilterLeg2"),
  fCutType(kBothLegs)
{
  //
  // Default contructor
  //
}

//________________________________________________________________________
AliDielectronPairLegCuts::AliDielectronPairLegCuts(const char* name, const char* title) :
  AliAnalysisCuts(name,title),
  fFilterLeg1("PairFilterLeg1","PairFilterLeg1"),
  fFilterLeg2("PairFilterLeg2","PairFilterLeg2"),
  fCutType(kBothLegs)
{
  //
  // Named contructor
  //
}

AliDielectronPairLegCuts::~AliDielectronPairLegCuts()
{
  //
  // Dtor
  //

  // remove all cuts in the leg2 filter which are also in the leg1 filter
  // to avoid doule deletion
  TIter nextCut(fFilterLeg1.GetCuts());
  TObject *o=0x0;
  while ( (o=nextCut()) ) fFilterLeg2.GetCuts()->Remove(o);
}

//________________________________________________________________________
Bool_t AliDielectronPairLegCuts::IsSelected(TObject* track)
{
  //
  // check if cuts are fulfilled
  //
  
  //check if we have a AliDielectronPair
  AliDielectronPair *pair=dynamic_cast<AliDielectronPair*>(track);
  if (!pair) return kFALSE;

  //get both legs
  AliVParticle *leg1=pair->GetFirstDaughterP();
  AliVParticle *leg2=pair->GetSecondDaughterP();

  // In case of the internal train, there seems to be an issue in the event mixing (storage of track pointers?),
  // which sometimes causes missing legs. Histograms of such pair types show much less entries compared to 
  // running without the internal train, while same-events histograms are identical.
  if (!leg1) { Printf("AliDielectronPairLegCuts::IsSelected(): leg 1 not found! pair type: %i", pair->GetType());
    return kFALSE;
  }
  if (!leg2) { Printf("AliDielectronPairLegCuts::IsSelected(): leg 2 not found! pair type: %i", pair->GetType());
    return kFALSE;
  }
  
  //mask used to require that all cuts are fulfilled
  UInt_t selectedMaskLeg1=(1<<fFilterLeg1.GetCuts()->GetEntries())-1;
  UInt_t selectedMaskLeg2=(1<<fFilterLeg2.GetCuts()->GetEntries())-1;
  
  //test cuts
  Bool_t isLeg1selected=(fFilterLeg1.IsSelected(leg1)==selectedMaskLeg1);
  if(fCutType==kBothLegs && !isLeg1selected) {
    SetSelected(isLeg1selected);
    return isLeg1selected;
  }
  Bool_t isLeg2selected=(fFilterLeg2.IsSelected(leg2)==selectedMaskLeg2);
  
  Bool_t isLeg1selectedMirror=(fFilterLeg1.IsSelected(leg2)==selectedMaskLeg1);
  Bool_t isLeg2selectedMirror=(fFilterLeg2.IsSelected(leg1)==selectedMaskLeg2);
  
  Bool_t isSelected=isLeg1selected&&isLeg2selected;
  if (fCutType==kAnyLeg)
    isSelected=isLeg1selected||isLeg2selected;
  
  if (fCutType==kMixLegs)
    isSelected=(isLeg1selected&&isLeg2selected)||(isLeg1selectedMirror&&isLeg2selectedMirror);
  
  SetSelected(isSelected);
  return isSelected;
}



