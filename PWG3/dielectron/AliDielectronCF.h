#ifndef ALIDIELECTRONCF_H
#define ALIDIELECTRONCF_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//#############################################################
//#                                                           #
//#             Class AliDielectronCF                         #
//#       Dielectron Correction Framework Manager             #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Frederick Kramer,   Uni Ffm, / Frederick.Kramer@cern.ch #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################



#include <TNamed.h>
#include <TVectorDfwd.h>
#include "AliDielectronVarManager.h"

class AliAnalysisCuts;
class AliAnalysisFilter;
class AliCFContainer;
class AliDielectronPair;
class TObjArray;

class AliDielectronCF : public TNamed {
public:
  enum {kNmaxAddSteps=50};
  
  AliDielectronCF();
  AliDielectronCF(const char* name, const char* title);
  virtual ~AliDielectronCF();

  void SetStepForMCtruth(Bool_t steps=kTRUE)           { fStepForMCtruth=steps;           }
  void SetStepForNoCutsMCmotherPid(Bool_t steps=kTRUE) { fStepForNoCutsMCmotherPid=steps; }
  void SetStepForAfterAllCuts(Bool_t steps=kTRUE)      { fStepForAfterAllCuts=steps;      }
  void SetStepForPreFilter(Bool_t steps=kTRUE)         { fStepForPreFilter=steps;         }
  void SetStepsForEachCut(Bool_t steps=kTRUE)          { fStepsForEachCut=steps;          }
  void SetStepsForCutsIncreasing(Bool_t steps=kTRUE)   { fStepsForCutsIncreasing=steps;   }
  void SetStepsForSignal(Bool_t steps=kTRUE)           { fStepsForSignal=steps;           }
  void SetStepsForBackground(Bool_t steps=kTRUE)       { fStepsForBackground=steps;       }
  
  void SetPdgMother(Int_t pdg) { fPdgMother=pdg; }
  
  void AddStepMask(UInt_t mask)                  { fStepMasks[fNStepMasks++]=mask; }
  
  void AddVariable(AliDielectronVarManager::ValueTypes type, Int_t nbins,
                   Double_t min, Double_t max, Bool_t leg=kFALSE, Bool_t log=kFALSE);
  void AddVariable(AliDielectronVarManager::ValueTypes type, const char* binLimitStr, Bool_t leg=kFALSE);
  void AddVariable(AliDielectronVarManager::ValueTypes type, TVectorD *binLimits, Bool_t leg=kFALSE);
  
  void InitialiseContainer(const AliAnalysisFilter& filter);
  
//   void Fill(UInt_t mask, const TObject *particle);
  void Fill(UInt_t mask, const AliDielectronPair *particle);
  void FillMC(const TObject *particle);
  
  AliCFContainer* GetContainer() const { return fCfContainer; }
  
private:
  UInt_t          fVariables[AliDielectronVarManager::kNMaxValues]; //configured variables
  UInt_t          fVariablesLeg[AliDielectronVarManager::kNMaxValues]; //configured variables for the legs
  
  Int_t           fNSteps;                     // number of selection steps
  
  Int_t           fNVars;                      // number of variables
  TObjArray      *fVarBinLimits;               // array of bin limits
  
  Int_t           fNVarsLeg;                   // number of variables for the legs
  TObjArray      *fVarBinLimitsLeg;            //  array of bin limits of the legs
  
  Int_t           fNCuts;                         // Number of cuts in the filter concerned

  Double_t        *fValues;                       //! Value array for filling the container
  
  Bool_t fStepForMCtruth;               //create a step for the MC truth
  Bool_t fStepForNoCutsMCmotherPid;     //create a step for before cuts, but with MC truth of the mother
  Bool_t fStepForAfterAllCuts;          //create a step for before cuts, but with MC truth of the mother
  Bool_t fStepForPreFilter;             //create a step after PreFilter application
  Bool_t fStepsForEachCut;              //create steps for each cut?
  Bool_t fStepsForCutsIncreasing;       //create steps for increasing cut combinatons?
                                        //e.g. cut1&cut2, cut1&cut2&cut3 ...
  Bool_t fStepsForSignal;               //steps for pure signal
  Bool_t fStepsForBackground;           //steps for pure background
  
  UInt_t fStepMasks[kNmaxAddSteps];      //steps for additional cut combinatons
  UInt_t fNStepMasks;                    //number of configured step masks

  Int_t fPdgMother;                      //Pdg code of MCtruth validation
  AliCFContainer* fCfContainer;          //the CF container

  Bool_t fHasMC;                         //if MC info is available
  Int_t  fNAddSteps;                     //number of additional MC related steps per cut step

  TVectorD* MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax) const;
  TVectorD* MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax) const;
  
  AliDielectronCF(const AliDielectronCF &c);
  AliDielectronCF &operator=(const AliDielectronCF &c);
  
  ClassDef(AliDielectronCF,4)  //Dielectron Correction Framework handler
};

#endif
