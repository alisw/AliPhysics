#ifndef ALIDIELECTRONVARCUTS_H
#define ALIDIELECTRONVARCUTS_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           #
//#         Class AliDielectronVarCuts                        #
//#         Provide cuts for all variables handled in         #
//#           AliDielectronVarManager                         #
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

#include <Rtypes.h>
#include <TBits.h>

#include <AliAnalysisCuts.h>
#include "AliDielectronVarManager.h"

class THnBase;
class AliDielectronVarCuts : public AliAnalysisCuts {
public:
  // Whether all cut criteria have to be fulfilled of just any
  enum CutType { kAll=0, kAny };
  enum EVarCutsOperation{ kNone = 0,
                          kAdd,
                          kSubtract,
                          kMutliply,
                          kDivide
                        };

  AliDielectronVarCuts();
  AliDielectronVarCuts(const char* name, const char* title);
  virtual ~AliDielectronVarCuts();
  //TODO: make copy constructor and assignment operator public
  void AddCut(AliDielectronVarManager::ValueTypes type, Double_t min, Double_t max, Bool_t excludeRange=kFALSE);
  void AddCut(AliDielectronVarManager::ValueTypes type, Double_t value, Bool_t excludeRange=kFALSE);
  void AddBitCut(AliDielectronVarManager::ValueTypes type, UInt_t bit, Bool_t excludeRange=kFALSE);
  void AddCut(AliDielectronVarManager::ValueTypes type, Double_t min, THnBase * const max,  Bool_t excludeRange=kFALSE);
  void AddCut(AliDielectronVarManager::ValueTypes typeA, AliDielectronVarManager::ValueTypes typeB, Double_t min, Double_t max, EVarCutsOperation operation, Bool_t excludeRange=kFALSE);
  void InvertCuts();

  // setters
  void    SetCutOnMCtruth(Bool_t mc=kTRUE) { fCutOnMCtruth=mc; }
  void    SetCutType(CutType type)         { fCutType=type;    }

  // getters
  Bool_t  GetCutOnMCtruth() const { return fCutOnMCtruth; }
  CutType GetCutType()      const { return fCutType;      }

  Int_t GetNCuts() { return fNActiveCuts; }

  //
  //Analysis cuts interface
  //
  virtual Bool_t IsSelected(TObject* track);
  virtual Bool_t IsSelected(Double_t* values);
  virtual Bool_t IsSelected(TList*   /* list */ ) {return kFALSE;}

//   virtual Bool_t IsSelected(TObject* track, TObject */*event*/=0);
//   virtual Long64_t Merge(TCollection* /* list */)      { return 0; }

  //
  // Cut information
  //
  virtual UInt_t GetSelectedCutsMask() const { return fSelectedCutsMask; }

  virtual void Print(const Option_t* option = "") const;
  const char*  GetCutName(Int_t iCut) const;
  Bool_t       IsCutOnVariableX(Int_t iCut, Int_t varNumber) const;
  Int_t        GetCutLimits(Int_t iCut, Double_t &cutMin, Double_t &cutMax) const;


 private:

  TBits     *fUsedVars;            // list of used variables
  UShort_t  fActiveCuts[AliDielectronVarManager::kNMaxValues];       // list of activated cuts
  UShort_t  fNActiveCuts;                      // number of acive cuts
  UInt_t    fActiveCutsMask;                   // mask of active cuts

  UInt_t   fSelectedCutsMask;                 // Maks of selected cuts, is available after calling IsSelected

  Bool_t   fCutOnMCtruth;                     // whether to cut on the MC truth of the particle

  CutType  fCutType;                          // type of the cut: any, all

  Double_t fCutMin[AliDielectronVarManager::kNMaxValues];           // minimum values for the cuts
  Double_t fCutMax[AliDielectronVarManager::kNMaxValues];           // maximum values for the cuts
  Bool_t fCutExclude[AliDielectronVarManager::kNMaxValues];         // inverse cut logic?
  Bool_t fBitCut[AliDielectronVarManager::kNMaxValues];             // bit cut
  THnBase  *fUpperCut[AliDielectronVarManager::kNMaxValues];        // use object as upper cut
  EVarCutsOperation fVarOperation[AliDielectronVarManager::kNMaxValues]; // operation between two vars, attention in principle kNMaxValues could be exceeded by the cut logic use with care

  AliDielectronVarCuts(const AliDielectronVarCuts &c);
  AliDielectronVarCuts &operator=(const AliDielectronVarCuts &c);

  ClassDef(AliDielectronVarCuts,7)         //Cut class providing cuts to all infomation available for the AliVParticle interface
};


//
//Inline functions
//
inline void AliDielectronVarCuts::AddCut(AliDielectronVarManager::ValueTypes type, Double_t value, Bool_t excludeRange)
{
  //
  // Set cut in a small delta around value
  //
  const Double_t kDelta=1e-20;
  AddCut(type,value-kDelta,value+kDelta, excludeRange);
}

#endif
