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
    fVarOperation[i]=AliDielectronVarCuts::kNone;
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
    fVarOperation[i]=AliDielectronVarCuts::kNone;
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
  Double_t opResultValue = 0.;

  for (Int_t iCut=0; iCut<fNActiveCuts; ++iCut){
    Int_t cut=fActiveCuts[iCut];
    SETBIT(fSelectedCutsMask,iCut);

    // apply 'bit cut'
    if(fBitCut[iCut]) {
      if ( (TESTBIT((UInt_t)values[cut],(UInt_t)fCutMin[iCut]))^(!fCutExclude[iCut]) )  CLRBIT(fSelectedCutsMask,iCut);
    }
    else {
      if(!(fVarOperation[iCut] == AliDielectronVarCuts::kNone)){
        Int_t cutB = fActiveCuts[iCut+1];
        switch (fVarOperation[iCut]) {
          case AliDielectronVarCuts::kNone :
            AliFatal("AliDielectronVarCuts: You should not be here check the code!!!");
          break;
          case AliDielectronVarCuts::kAdd :
            opResultValue = values[cut]+values[cutB];
          break;
          case AliDielectronVarCuts::kSubtract :
            opResultValue = values[cut]-values[cutB];
          break;
          case AliDielectronVarCuts::kMutliply :
            opResultValue = values[cut]*values[cutB];
          break;
          case AliDielectronVarCuts::kDivide :
            opResultValue = values[cut]/values[cutB];
          break;
        }
        if( ((opResultValue<fCutMin[iCut]) || (opResultValue>fCutMax[iCut]))^fCutExclude[iCut] )          CLRBIT(fSelectedCutsMask,iCut);

        iCut++;
        SETBIT(fSelectedCutsMask,iCut); // Set bit to true for cut only containing second variable
      }
      else{
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
Bool_t AliDielectronVarCuts::IsSelected(Double_t* values)
{
  //
  // Make cut decision
  //

  //reset
  fSelectedCutsMask=0;
  SetSelected(kFALSE);

  Double_t opResultValue = 0.;
  for (Int_t iCut=0; iCut<fNActiveCuts; ++iCut){
    Int_t cut=fActiveCuts[iCut];
    SETBIT(fSelectedCutsMask,iCut);

    // apply 'bit cut'
    if(fBitCut[iCut]) {
      if ( (TESTBIT((UInt_t)values[cut],(UInt_t)fCutMin[iCut]))^(!fCutExclude[iCut]) )  CLRBIT(fSelectedCutsMask,iCut);
    }
    else {
      if(!(fVarOperation[iCut] == AliDielectronVarCuts::kNone)){
        Int_t cutB = fActiveCuts[iCut+1];
        switch (fVarOperation[iCut]) {
          case AliDielectronVarCuts::kNone :
            AliFatal("AliDielectronVarCuts: You should not be here check the code!!!");
            break;
          case AliDielectronVarCuts::kAdd :
            opResultValue = values[cut]+values[cutB];
            break;
          case AliDielectronVarCuts::kSubtract :
            opResultValue = values[cut]-values[cutB];
            break;
          case AliDielectronVarCuts::kMutliply :
            opResultValue = values[cut]*values[cutB];
            break;
          case AliDielectronVarCuts::kDivide :
            opResultValue = values[cut]/values[cutB];
            break;
        }
        if( ((opResultValue<fCutMin[iCut]) || (opResultValue>fCutMax[iCut]))^fCutExclude[iCut] )          CLRBIT(fSelectedCutsMask,iCut);

        iCut++;
        SETBIT(fSelectedCutsMask,iCut); // Set bit to true for cut only containing second variable
      }
      else{
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

//_________________________________________________________________________
void AliDielectronVarCuts::AddCut(AliDielectronVarManager::ValueTypes typeA, AliDielectronVarManager::ValueTypes typeB, Double_t min, Double_t max, EVarCutsOperation operation, Bool_t excludeRange) {
  //
  // Set cut range, define operation between the two VarManager vars and activate cut
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

  fActiveCuts[fNActiveCuts]=(UShort_t)typeA;
  fUsedVars->SetBitNumber(typeA,kTRUE);
  fUsedVars->SetBitNumber(typeB,kTRUE);

  fVarOperation[fNActiveCuts] = operation;
  ++fNActiveCuts;
  SETBIT(fActiveCutsMask,fNActiveCuts);
  fActiveCuts[fNActiveCuts]=(UShort_t)typeB;
  ++fNActiveCuts;

}

//________________________________________________________________________
void AliDielectronVarCuts::InvertCuts()
{
  //
  // Invert the cut logic of this AliDielectronVarCuts object
  //
  for (Int_t iCut=0; iCut<fNActiveCuts; ++iCut){
    fCutExclude[iCut] = !fCutExclude[iCut];
  }
  if (fCutType==kAll) fCutType=kAny;
  else                fCutType=kAll;
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

//________________________________________________________________________
const char* AliDielectronVarCuts::GetCutName(Int_t iCut) const
{
  //
  // Return name of the cut at position iCut
  //
  if (iCut > fNActiveCuts-1) {
    printf("AliDielectronVarCuts::GetCutName(iCut): out of range! iCut=%i \n", iCut);
    return "ERROR";
  }
  Int_t cut=(Int_t)fActiveCuts[iCut];
  return AliDielectronVarManager::GetValueName((Int_t)cut);
}

//________________________________________________________________________
Bool_t AliDielectronVarCuts::IsCutOnVariableX(Int_t iCut, Int_t varNumber) const
{
  //
  // Use it to check if the cut at position iCut is for the variable varNumber (e.g. varNumber = AliDielectronVarManager::kPt)
  //
  if (iCut > fNActiveCuts-1) {
    printf("AliDielectronVarCuts::IsCutOnVariableX(iCut, ...): out of range! iCut=%i \n", iCut);
    return kFALSE;
  }
  Int_t cut=(Int_t)fActiveCuts[iCut];
  return (varNumber == (Int_t)cut);
}

//________________________________________________________________________
Int_t AliDielectronVarCuts::GetCutLimits(Int_t iCut, Double_t &cutMin, Double_t &cutMax) const
{
  //
  // Return min and max of the cut at position iCut
  //
  if (iCut > fNActiveCuts-1) {
    printf("AliDielectronVarCuts::GetCutLimits(iCut, ...): out of range! iCut=%i \n", iCut);
    return -1;
  }
  cutMin = -999;
  cutMax = -999;

  //Bool_t inverse=fCutExclude[iCut];
  Bool_t bitcut=fBitCut[iCut];
  Bool_t objcut=fUpperCut[iCut];

  if(!bitcut && !objcut) {
    // standard cut
    cutMin = fCutMin[iCut];
    cutMax = fCutMax[iCut];
    //if (inverse){}
  }
  else { // TODO: add the other possibilities.
  	printf("Cut %02d: not a standard cut. Not supported by this getter.\n", iCut);
    return -1;
  }

  return iCut;
}
