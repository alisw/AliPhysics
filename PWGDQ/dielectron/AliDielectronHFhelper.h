#ifndef ALIDIELECTRONHFHELPER_H
#define ALIDIELECTRONHFHELPER_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           #
//#             Class AliDielectronHF                         #
//#       Dielectron Histogram Framework helper         #
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
  enum ECollectType { kSE=0, kME, kMEOS, kROT, kAll };
  enum { kMaxCuts=20 };

  //AliDielectronHFhelper();
  AliDielectronHFhelper(const char* filename, const char* container);

  virtual ~AliDielectronHFhelper();
  void SetHFArray(const char* filename, const char* container);

  void SetRangeUser(const char *varname, Double_t min, Double_t max, Bool_t leg=kFALSE);
  void SetRangeUser(AliDielectronVarManager::ValueTypes type, Double_t min, Double_t max, Bool_t leg=kFALSE);
  void UnsetRangeUser(const char* varname, Bool_t leg=kFALSE);
  void UnsetRangeUser(AliDielectronVarManager::ValueTypes type, Bool_t leg=kFALSE);

  TObjArray* CollectHistos();

  TH1* GetHistogram(const char *step, TObjArray *histArr=0x0);
  TH1* FindHistograms(TObjArray *histos);
  TH1* MergeHistos(TObjArray *arr);

  void CheckCuts(TObjArray *arr);
  virtual void Print(const Option_t* option ="") const ;
  void PrintCuts();

private:
  TObjArray *fArrPairType;         // array of pair types, sources or steps
  TObjArray *fCutVars;             // array for cut variables
  TVectorD fCutLowLimits;
  TVectorD fCutUpLimits;

  AliDielectronHFhelper(const AliDielectronHFhelper &c);
  AliDielectronHFhelper &operator=(const AliDielectronHFhelper &c);

  ClassDef(AliDielectronHFhelper,0)                   // HF  helper class
};

//
// Inline functions
//
#endif

