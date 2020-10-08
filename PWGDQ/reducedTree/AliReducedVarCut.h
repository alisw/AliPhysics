// Class for cutting on ALICE AliReducedVarManager information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
//   07/09/2016

#ifndef ALIREDUCEDVARCUT_H
#define ALIREDUCEDVARCUT_H

#include <TF1.h>
#include "AliReducedInfoCut.h"
#include "AliReducedVarManager.h"

//_________________________________________________________________________
class AliReducedVarCut : public AliReducedInfoCut {
   
 public:
  AliReducedVarCut();
  AliReducedVarCut(const Char_t* name, const Char_t* title);
  virtual ~AliReducedVarCut();

  enum Constants {
     kNMaxCuts=200       // maximum number of cuts
  }; 
  
  // NOTE: Apply a selection on variable "var" to be in the range [cutLow,cutHigh] or outside this range if "exclude" is set to true
  // NOTE: If a dependent variable is specified, then the selection is applied only if the dependent variable is in the range [depCutLow,depCutHigh]
  // NOTE:       or outside if "depCutExclude" is true
  void AddCut(Int_t var, Float_t cutLow, Float_t cutHigh, Bool_t exclude = kFALSE, 
              Int_t dependentVar=AliReducedVarManager::kNothing, Float_t depCutLow=0., Float_t depCutHigh=0., Bool_t depCutExclude=kFALSE,
              Int_t dependentVar2=AliReducedVarManager::kNothing, Float_t depCut2Low=0., Float_t depCut2High=0., Bool_t depCut2Exclude=kFALSE);
  // NOTE: Define cuts which use functions of a defined variable instead of a constant cut; the logic of the arguments is the same as for the above function
  // NOTE: The use case is for cuts on correlations between 2 variables
  void AddCut(Int_t var, Float_t cutLow, TF1* funcCutHigh, Bool_t exclude = kFALSE,
              Int_t dependentVar=AliReducedVarManager::kNothing, Float_t depCutLow=0., Float_t depCutHigh=0., Bool_t depCutExclude=kFALSE,
              Int_t dependentVar2=AliReducedVarManager::kNothing, Float_t depCut2Low=0., Float_t depCut2High=0., Bool_t depCut2Exclude=kFALSE);
  void AddCut(Int_t var, TF1* funcCutLow, Float_t cutHigh, Bool_t exclude = kFALSE,
              Int_t dependentVar=AliReducedVarManager::kNothing, Float_t depCutLow=0., Float_t depCutHigh=0., Bool_t depCutExclude=kFALSE,
              Int_t dependentVar2=AliReducedVarManager::kNothing, Float_t depCut2Low=0., Float_t depCut2High=0., Bool_t depCut2Exclude=kFALSE);
  void AddCut(Int_t var, TF1* funcCutLow, TF1* funcCutHigh, Bool_t exclude = kFALSE,
              Int_t dependentVar=AliReducedVarManager::kNothing, Float_t depCutLow=0., Float_t depCutHigh=0., Bool_t depCutExclude=kFALSE,
              Int_t dependentVar2=AliReducedVarManager::kNothing, Float_t depCut2Low=0., Float_t depCut2High=0., Bool_t depCut2Exclude=kFALSE);

  virtual Bool_t IsSelected(TObject* obj);
  virtual Bool_t IsSelected(Float_t* values);
  virtual Bool_t IsSelected(TObject* obj, Float_t* values);
  
 protected: 
  
   Int_t       fNCuts;                                    // number of enabled cuts
   Short_t   fCutVariables[kNMaxCuts];    // list of variables enabled to cut on
   Short_t   fDependentVariable[kNMaxCuts];      // one can specify that a cut is applied only if a specified variable is within a given range (or outside)   
   Short_t   fDependentVariable2[kNMaxCuts];     // one can specify that a cut is applied only if a specified variable is within a given range (or outside)
   Float_t   fCutLow[kNMaxCuts];          // lower limit for all quantities defined in the AliReducedVarManager
   Float_t   fCutHigh[kNMaxCuts];         // higher limit for all quantities defined in the AliReducedVarManager 
   Bool_t    fCutExclude[kNMaxCuts];    // if true, then use the selection range for exclusion
   Bool_t    fCutHasDependentVariable[kNMaxCuts];    // if true, check the requirements of the dependent variable
   Bool_t    fCutHasDependentVariable2[kNMaxCuts];   // if true, check the requirements of the dependent variable
   Float_t   fDependentVariableCutLow[kNMaxCuts];    // lower limit for the dependent variable
   Float_t   fDependentVariable2CutLow[kNMaxCuts];   // lower limit for the dependent variable
   Float_t   fDependentVariableCutHigh[kNMaxCuts];   // higher limit for the dependent variable
   Float_t   fDependentVariable2CutHigh[kNMaxCuts];  // higher limit for the dependent variable
   Bool_t    fDependentVariableExclude[kNMaxCuts];  // if true, then use the dependent variable range as exclusion
   Bool_t    fDependentVariable2Exclude[kNMaxCuts];  // if true, then use the dependent variable range as exclusion
   TF1*      fFuncCutLow[kNMaxCuts];     // low cut functions
   TF1*      fFuncCutHigh[kNMaxCuts];    // high cut functions
      
  AliReducedVarCut(const AliReducedVarCut &c);
  AliReducedVarCut& operator= (const AliReducedVarCut &c);
  
  ClassDef(AliReducedVarCut,4);
};

#endif
