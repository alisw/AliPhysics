/*
***********************************************************
  Implementation of AliReducedCompositeCut class.
  Contact: iarsene@cern.ch
  2018/02/06
  *********************************************************
*/

#include <iostream>
using std::cout;
using std::endl;

#ifndef ALIREDUCEDCOMPOSITECUT_H
#include "AliReducedCompositeCut.h"
#endif

ClassImp(AliReducedCompositeCut)

//____________________________________________________________________________
AliReducedCompositeCut::AliReducedCompositeCut(Bool_t useAND /* = kTRUE */) :
  AliReducedInfoCut(),
  fOptionUseAND(useAND)
{
  //
  // default constructor
  //
  fCuts.SetOwner(kTRUE);
}

//____________________________________________________________________________
AliReducedCompositeCut::AliReducedCompositeCut(const Char_t* name, const Char_t* title, Bool_t useAND /* = kTRUE */) :
  AliReducedInfoCut(name, title),
  fOptionUseAND(useAND)
{
  //
  // named constructor
  //
   fCuts.SetOwner(kTRUE);
}

//____________________________________________________________________________
AliReducedCompositeCut::~AliReducedCompositeCut() {
  //
  // destructor
  //
}

//____________________________________________________________________________
Bool_t AliReducedCompositeCut::IsSelected(TObject* obj, Float_t* values) {
  //
  // apply cuts
  //
  for(Int_t iCut=0; iCut<fCuts.GetEntries(); ++iCut) {
     AliReducedInfoCut* cut = (AliReducedInfoCut*)fCuts.At(iCut);
     
     if(fOptionUseAND && !cut->IsSelected(obj, values)) return kFALSE;
     if(!fOptionUseAND && cut->IsSelected(obj, values)) return kTRUE;
  }
  
  if(fOptionUseAND) return kTRUE;
  else return kFALSE;
}

//____________________________________________________________________________
Bool_t AliReducedCompositeCut::IsSelected(Float_t* values) {
  //
  // apply cuts
  //
  for(Int_t iCut=0; iCut<fCuts.GetEntries(); ++iCut) {
     AliReducedInfoCut* cut = (AliReducedInfoCut*)fCuts.At(iCut);

     if(fOptionUseAND && !cut->IsSelected(values)) return kFALSE;
     if(!fOptionUseAND && cut->IsSelected(values)) return kTRUE;
  }

  if(fOptionUseAND) return kTRUE;
  else return kFALSE;
}

//____________________________________________________________________________
Bool_t AliReducedCompositeCut::IsSelected(TObject* obj) {
   //
   // apply cuts
   //
   for(Int_t iCut=0; iCut<fCuts.GetEntries(); ++iCut) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fCuts.At(iCut);
      
      if(fOptionUseAND && !cut->IsSelected(obj)) return kFALSE;
      if(!fOptionUseAND && cut->IsSelected(obj)) return kTRUE;
   }
   
   if(fOptionUseAND) return kTRUE;
   else return kFALSE;
}
