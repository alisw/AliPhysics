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

  // NOTE: Apply a selection on variable "var" to be in the range [cutLow,cutHigh] or outside this range if "exclude" is set to true
  // NOTE: If a dependent variable is specified, then the selection is applied only if the dependent variable is in the range [depCutLow,depCutHigh]
  // NOTE:       or outside if "depCutExclude" is true
  void AddCut(AliReducedVarManager::Variables var, Float_t cutLow, Float_t cutHigh, Bool_t exclude = kFALSE, 
              AliReducedVarManager::Variables dependentVar=AliReducedVarManager::kNothing, Float_t depCutLow=0., Float_t depCutHigh=0., Bool_t depCutExclude=kFALSE);
  // NOTE: Define cuts which use functions of a defined variable instead of a constant cut; the logic of the arguments is the same as for the above function
  // NOTE: The use case is for cuts on correlations between 2 variables
  void AddCut(AliReducedVarManager::Variables var, Float_t cutLow, TF1* funcCutHigh, Bool_t exclude = kFALSE, AliReducedVarManager::Variables  
              dependentVar=AliReducedVarManager::kNothing, Float_t depCutLow=0., Float_t depCutHigh=0., Bool_t depCutExclude=kFALSE);
  void AddCut(AliReducedVarManager::Variables var, TF1* funcCutLow, Float_t cutHigh, Bool_t exclude = kFALSE, AliReducedVarManager::Variables 
              dependentVar=AliReducedVarManager::kNothing, Float_t depCutLow=0., Float_t depCutHigh=0., Bool_t depCutExclude=kFALSE);
  void AddCut(AliReducedVarManager::Variables var, TF1* funcCutLow, TF1* funcCutHigh, Bool_t exclude = kFALSE, AliReducedVarManager::Variables  
              dependentVar=AliReducedVarManager::kNothing, Float_t depCutLow=0., Float_t depCutHigh=0., Bool_t depCutExclude=kFALSE);
  
  virtual Bool_t IsSelected(TObject* obj);
  virtual Bool_t IsSelected(Float_t* values);
  virtual Bool_t IsSelected(TObject* obj, Float_t* values);
  
 protected: 
  
   Int_t       fNCuts;                                                                 // number of enabled cuts
   Short_t   fCutVariables[AliReducedVarManager::kNVars];    // list of variables enabled to cut on
   Short_t   fDependentVariable[AliReducedVarManager::kNVars];      // one can specify that a cut is applied only if a specified variable is within a given range (or outside)   
   Float_t   fCutLow[AliReducedVarManager::kNVars];          // lower limit for all quantities defined in the AliReducedVarManager 
   Float_t   fCutHigh[AliReducedVarManager::kNVars];         // higher limit for all quantities defined in the AliReducedVarManager 
   Bool_t    fCutExclude[AliReducedVarManager::kNVars];    // if true, then use the selection range for exclusion
   Bool_t    fCutHasDependentVariable[AliReducedVarManager::kNVars];    // if true, check the requirements of the dependent variable
   Float_t   fDependentVariableCutLow[AliReducedVarManager::kNVars];    // lower limit for the dependent variable
   Float_t   fDependentVariableCutHigh[AliReducedVarManager::kNVars];   // higher limit for the dependent variable
   Bool_t    fDependentVariableExclude[AliReducedVarManager::kNVars];  // if true, then use the dependent variable range as exclusion
   TF1*      fFuncCutLow[AliReducedVarManager::kNVars];     // low cut functions
   TF1*      fFuncCutHigh[AliReducedVarManager::kNVars];    // high cut functions
      
  AliReducedVarCut(const AliReducedVarCut &c);
  AliReducedVarCut& operator= (const AliReducedVarCut &c);
  
  ClassDef(AliReducedVarCut,1);
};

#endif
