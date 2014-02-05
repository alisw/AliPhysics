#ifndef ALIDIELECTRONHFHELPER_H
#define ALIDIELECTRONHFHELPER_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           #
//#             Class AliDielectronHFhelper                   #
//#       Dielectron Histogram Framework helper               #
//#                                                           #
//#  Authors:                                                 #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#                                                           #
//#############################################################



#include <TNamed.h>
#include <TVectorD.h>

#include "AliDielectronVarManager.h"

class AliDielectronHFhelper : public TNamed {
public:
  enum { kMaxCuts=20 };

  //AliDielectronHFhelper();
  AliDielectronHFhelper(const char* filename, const char* container);

  virtual ~AliDielectronHFhelper();
  void SetHFArray(const char* filename, const char* container);

  void SetRangeUser(const char *varname, Double_t min, Double_t max, Bool_t leg=kFALSE);
  void SetRangeUser(AliDielectronVarManager::ValueTypes type, Double_t min, Double_t max, Bool_t leg=kFALSE);
  void UnsetRangeUser(const char* varname, Bool_t leg=kFALSE);
  void UnsetRangeUser(AliDielectronVarManager::ValueTypes type, Bool_t leg=kFALSE);

  // getter functions
  Int_t GetNSteps() const {return fMainArr->GetEntries(); }

  TObjArray* CollectHistos(AliDielectronVarManager::ValueTypes varx,
			   AliDielectronVarManager::ValueTypes vary=AliDielectronVarManager::kNMaxValues,
			   AliDielectronVarManager::ValueTypes varz=AliDielectronVarManager::kNMaxValues,
			   AliDielectronVarManager::ValueTypes varw=AliDielectronVarManager::kNMaxValues)
  { return CollectProfiles("hist",varx,vary,varz,varw); }
  TObjArray* CollectHistos(TString option,
			   AliDielectronVarManager::ValueTypes varx,
			   AliDielectronVarManager::ValueTypes vary=AliDielectronVarManager::kNMaxValues,
			   AliDielectronVarManager::ValueTypes varz=AliDielectronVarManager::kNMaxValues,
			   AliDielectronVarManager::ValueTypes varw=AliDielectronVarManager::kNMaxValues)
  { return CollectProfiles(Form("%s:hist",option.Data()),varx,vary,varz,varw); }
  TObjArray* CollectProfiles(TString option,
			     AliDielectronVarManager::ValueTypes varx,
			     AliDielectronVarManager::ValueTypes vary=AliDielectronVarManager::kNMaxValues,
			     AliDielectronVarManager::ValueTypes varz=AliDielectronVarManager::kNMaxValues,
			     AliDielectronVarManager::ValueTypes vart=AliDielectronVarManager::kNMaxValues);

  TObjArray* FindObjects(TObjArray *histos);
  TObjArray* Merge(TObjArray *arr);

  void CheckCuts(TObjArray *arr);
  virtual void Print(const Option_t* option ="") const ;
  void PrintCuts();

private:
  TObjArray *fMainArr;         // main array of pair types or sources
  TObjArray *fCutVars;         // array for cut variables
  TVectorD fCutLowLimits;      // vector to store the lower cut limits
  TVectorD fCutUpLimits;       // vector to store the upper cut limits

  AliDielectronHFhelper(const AliDielectronHFhelper &c);
  AliDielectronHFhelper &operator=(const AliDielectronHFhelper &c);

  ClassDef(AliDielectronHFhelper,1)                   // HF  helper class
};

//
// Inline functions
//
#endif

