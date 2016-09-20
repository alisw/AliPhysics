/*
***********************************************************
  Implementation of AliReducedVarCut class.
  Contact: iarsene@cern.ch
  2016/09/07
  *********************************************************
*/

#include <iostream>
using std::cout;
using std::endl;

#ifndef ALIREDUCEDVARCUT_H
#include "AliReducedVarCut.h"
#endif

#include "AliReducedBaseTrack.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedPairInfo.h"
#include "AliReducedVarManager.h"

ClassImp(AliReducedVarCut)

//____________________________________________________________________________
AliReducedVarCut::AliReducedVarCut() :
  AliReducedInfoCut(),
  fNCuts(0)
{
  //
  // default constructor
  //
  for(Int_t i=0; i<AliReducedVarManager::kNVars; ++i) { 
     fCutVariables[i] = AliReducedVarManager::kNothing;
     fDependentVariable[i] = AliReducedVarManager::kNothing;
     fCutLow[i] = 0.;
     fCutHigh[i] = 0.;
     fCutExclude[i] = kFALSE;
     fCutHasDependentVariable[i] = kFALSE;
     fDependentVariableCutLow[i] = 0.;
     fDependentVariableCutHigh[i] = 0.;
     fDependentVariableExclude[i] = kFALSE;
  }
}

//____________________________________________________________________________
AliReducedVarCut::AliReducedVarCut(const Char_t* name, const Char_t* title) :
  AliReducedInfoCut(name, title),
  fNCuts(0)
{
  //
  // named constructor
  //
   for(Int_t i=0; i<AliReducedVarManager::kNVars; ++i) { 
      fCutVariables[i] = AliReducedVarManager::kNothing;
      fDependentVariable[i] = AliReducedVarManager::kNothing;
      fCutLow[i] = 0.;
      fCutHigh[i] = 0.;
      fCutExclude[i] = kFALSE;
      fCutHasDependentVariable[i] = kFALSE;
      fDependentVariableCutLow[i] = 0.;
      fDependentVariableCutHigh[i] = 0.;
      fDependentVariableExclude[i] = kFALSE;
   }
}

//____________________________________________________________________________
AliReducedVarCut::~AliReducedVarCut() {
  //
  // destructor
  //
}


//____________________________________________________________________________
void AliReducedVarCut::AddCut(AliReducedVarManager::Variables var, Float_t cutLow, Float_t cutHigh, Bool_t exclude /*= kFALSE*/, 
                                                   AliReducedVarManager::Variables dependentVar /*=AliReducedVarManager::kNothing*/, 
                                                   Float_t depCutLow /*=0.*/, Float_t depCutHigh /*=0.*/, Bool_t depCutExclude /*=kFALSE*/) {
   //
   //  Add a cut
   //
   fCutVariables[fNCuts] = var; fCutLow[fNCuts] = cutLow; fCutHigh[fNCuts] = cutHigh; fCutExclude[fNCuts] = exclude;
   AliReducedVarManager::SetUseVariable(var);
   if(dependentVar != AliReducedVarManager::kNothing) {
      fCutHasDependentVariable[fNCuts] = kTRUE;
      fDependentVariable[fNCuts] = dependentVar; 
      fDependentVariableCutLow[fNCuts] = depCutLow; fDependentVariableCutHigh[fNCuts] = depCutHigh; 
      fDependentVariableExclude[fNCuts] = depCutExclude;
      AliReducedVarManager::SetUseVariable(var);
   }
   fNCuts++;
}


//____________________________________________________________________________
Bool_t AliReducedVarCut::IsSelected(TObject* obj) {
  //
  // apply cuts
  //
  if(!obj->InheritsFrom(AliReducedBaseTrack::Class())) return kFALSE;
  
  //Fill values
  Float_t values[AliReducedVarManager::kNVars];
  if(obj->InheritsFrom(AliReducedBaseEvent::Class())) AliReducedVarManager::FillEventInfo((AliReducedBaseEvent*)obj, values);
  if(obj->InheritsFrom(AliReducedBaseTrack::Class())) AliReducedVarManager::FillTrackInfo((AliReducedBaseTrack*)obj, values);
  if(obj->InheritsFrom(AliReducedPairInfo::Class())) AliReducedVarManager::FillPairInfo((AliReducedPairInfo*)obj, values);
  
  return IsSelected(values);
}


//____________________________________________________________________________
Bool_t AliReducedVarCut::IsSelected(Float_t* values) {
   //
   // apply cuts
   //         
   for(Int_t i=0; i<fNCuts; ++i) {
      //cout << "AliReducedVarCut::IsSelected() icut " << i << endl;
      if(fCutHasDependentVariable[i]) {
         //cout << "AliReducedVarCut::IsSelected() has dependent var " << endl;
         Bool_t inRangeDep = (values[fDependentVariable[i]]>=fDependentVariableCutLow[i] && values[fDependentVariable[i]]<=fDependentVariableCutHigh[i]);
         //cout << "AliReducedVarCut::IsSelected() inRangeDep/fDepVar/val/cutLow/cutHigh: " << inRangeDep 
         //        << "/" << fDependentVariable[i] << "/" << values[fDependentVariable[i]] << "/"
         //        << fDependentVariableCutLow[i] << "/" << fDependentVariableCutHigh[i] << endl;
         //cout << "AliReducedVarCut::IsSelected() fDependentVariableExclude: " << fDependentVariableExclude[i] << endl;        
         // do not apply this cut if outside of the applicability range
         if(!inRangeDep && !fDependentVariableExclude[i]) continue;
         if(inRangeDep && fDependentVariableExclude[i]) continue;
      }
      Bool_t inRange = (values[fCutVariables[i]]>=fCutLow[i] && values[fCutVariables[i]]<=fCutHigh[i]);
      //cout << "AliReducedVarCut::IsSelected() inRange/fCutVariables/val/cutLow/cutHigh: " << inRange << "/" << fCutVariables[i]
      //        << "/" << values[fCutVariables[i]] << "/" << fCutLow[i] << "/" << fCutHigh[i] << "/" << fCutExclude[i] << endl;
      if(!inRange && !fCutExclude[i]) return kFALSE;
      if(inRange && fCutExclude[i]) return kFALSE;
   }
   
   return kTRUE;
}
