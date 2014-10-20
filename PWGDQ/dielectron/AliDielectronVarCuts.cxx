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
//   Julian Book   <Julian.Book@cern.ch>                                 //
/*



*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#include <THnBase.h>

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
    fBitCut[i]=kFALSE;
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
    fBitCut[i]=kFALSE;
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

    // apply 'bit cut'
    if(fBitCut[iCut]) {
      if ( (TESTBIT((UInt_t)values[cut],(UInt_t)fCutMin[iCut]))^(!fCutExclude[iCut]) )  CLRBIT(fSelectedCutsMask,iCut);
    }
    else {
      // standard var cuts
      if ( !fUpperCut[iCut] && ((values[cut]<fCutMin[iCut]) || (values[cut]>fCutMax[iCut]))^fCutExclude[iCut] ) {
	CLRBIT(fSelectedCutsMask,iCut);
      }
      else if ( fUpperCut[iCut]) {
	/// use a THnBase inherited cut object //
	Double_t *vals = new Double_t[fUpperCut[iCut]->GetNdimensions()];//={-1};
	// get array of values for the corresponding dimensions using axis names
	for(Int_t idim=0; idim<fUpperCut[iCut]->GetNdimensions(); idim++) {
	  vals[idim] = values[AliDielectronVarManager::GetValueType(fUpperCut[iCut]->GetAxis(idim)->GetName())];
	  // printf(" \t %s %.3f ",fUpperCut[iCut]->GetAxis(idim)->GetName(),vals[idim]);
	}
	// find bin for values (w/o creating it in case it is not filled)
	Long_t bin = fUpperCut[iCut]->GetBin(vals,kFALSE);
	Double_t cutMax = (bin>0 ? fUpperCut[iCut]->GetBinContent(bin) : -999. );
	if ( ((values[cut]<fCutMin[iCut]) || (values[cut]>cutMax))^fCutExclude[iCut] ) CLRBIT(fSelectedCutsMask,iCut);
	delete [] vals;
      }
    }
    // cut type and decision
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
void AliDielectronVarCuts::AddBitCut(AliDielectronVarManager::ValueTypes type, UInt_t bit, Bool_t excludeRange)
{
  //
  // Set cut range and activate it
  //
  fCutMin[fNActiveCuts]=bit;
  fCutExclude[fNActiveCuts]=excludeRange;
  fBitCut[fNActiveCuts]=kTRUE;
  SETBIT(fActiveCutsMask,fNActiveCuts);
  fActiveCuts[fNActiveCuts]=(UShort_t)type;
  fUsedVars->SetBitNumber(type,kTRUE);
  ++fNActiveCuts;
}

//________________________________________________________________________
void AliDielectronVarCuts::AddCut(AliDielectronVarManager::ValueTypes type, Double_t min, THnBase * const max,  Bool_t excludeRange)
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

  fUpperCut[fNActiveCuts]=max;
  // fill used variables into map
  for(Int_t idim=0; idim<fUpperCut[fNActiveCuts]->GetNdimensions(); idim++) {
    TString var(fUpperCut[fNActiveCuts]->GetAxis(idim)->GetName());
    fUsedVars->SetBitNumber(AliDielectronVarManager::GetValueType(var.Data()), kTRUE);
  }
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
    Bool_t bitcut=fBitCut[iCut];
    Bool_t objcut=fUpperCut[iCut];

    if(!bitcut && !objcut) {
      // standard cut
      if (!inverse){
	printf("Cut %02d: %f < %s < %f\n", iCut,
	       fCutMin[iCut], AliDielectronVarManager::GetValueName((Int_t)cut), fCutMax[iCut]);
      } else {
	printf("Cut %02d: !(%f < %s < %f)\n", iCut,
	       fCutMin[iCut], AliDielectronVarManager::GetValueName((Int_t)cut), fCutMax[iCut]);
      }
    }
    else if(bitcut) {
      // bit cut
      if (!inverse){
	printf("Cut %02d: %s & (1ULL<<%d) \n", iCut,
	       AliDielectronVarManager::GetValueName((Int_t)cut), (UInt_t)fCutMin[iCut]);
      } else {
	printf("Cut %02d: !(%s & (1ULL<<%d)) \n", iCut,
	       AliDielectronVarManager::GetValueName((Int_t)cut), (UInt_t)fCutMin[iCut]);
      }
    }
    else if(objcut) {
      // upper cut limit provided by object depending on 'dep'
      TString dep="";
      for(Int_t idim=0; idim<fUpperCut[iCut]->GetNdimensions(); idim++)
	dep+=Form("%s%s",(idim?",":""),fUpperCut[iCut]->GetAxis(idim)->GetName());

      if (!inverse){
	printf("Cut %02d: %f < %s < obj(%s)\n", iCut,
	       fCutMin[iCut], AliDielectronVarManager::GetValueName((Int_t)cut), dep.Data());
      } else {
	printf("Cut %02d: !(%f < %s < obj(%s))\n", iCut,
	       fCutMin[iCut], AliDielectronVarManager::GetValueName((Int_t)cut), dep.Data());
      }
    }
  } //loop over cuts
}
