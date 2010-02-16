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
#include <AliAnalysisCuts.h>
#include <AliDielectronVarManager.h>

class AliDielectronVarCuts : public AliAnalysisCuts {
public:
  AliDielectronVarCuts();
  AliDielectronVarCuts(const char* name, const char* title);
  virtual ~AliDielectronVarCuts();
  //TODO: make copy constructor and assignment operator public
  void AddCut(Double_t min, Double_t max, AliDielectronVarManager::ValueTypes type);

  //
  //Analysis cuts interface
  //
  virtual Bool_t IsSelected(TObject* track);
  virtual Bool_t IsSelected(TList*   /* list */ ) {return kFALSE;}
  
//   virtual Bool_t IsSelected(TObject* track, TObject */*event*/=0);
//   virtual Long64_t Merge(TCollection* /* list */)      { return 0; }
  
  //
  // Cut information
  //
  virtual UInt_t GetSelectedCutsMask() const { return fSelectedCutsMask; }
  Bool_t IsCutActive(AliDielectronVarManager::ValueTypes cut) { return TESTBIT(fActiveCutsMask,cut); }

private:

  UChar_t  fActiveCuts[AliDielectronVarManager::kNMaxValues];       // list of activated cuts
  UChar_t  fNActiveCuts;                      // number of acive cuts
  UInt_t   fActiveCutsMask;                   // maks of active cuts
  
  UInt_t   fSelectedCutsMask;                 // Maks of selected cuts, is available after calling IsSelected
  
  Double_t fCutMin[AliDielectronVarManager::kNMaxValues];           // minimum values for the cuts
  Double_t fCutMax[AliDielectronVarManager::kNMaxValues];           // maximum values for the cuts

  void ActivateCut(AliDielectronVarManager::ValueTypes cutName);
  
  AliDielectronVarCuts(const AliDielectronVarCuts &c);
  AliDielectronVarCuts &operator=(const AliDielectronVarCuts &c);
  
  ClassDef(AliDielectronVarCuts,1)         //Cut class providing cuts to all infomation available for the AliVParticle interface
};


//
//Inline functions
//
inline void AliDielectronVarCuts::AddCut(Double_t min, Double_t max, AliDielectronVarManager::ValueTypes type)
{
  //
  // Set cut range and activate it
  //
  fCutMin[type]=min;
  fCutMax[type]=max;
  ActivateCut(type);
}

#endif

