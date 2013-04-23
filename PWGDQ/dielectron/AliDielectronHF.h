#ifndef ALIDIELECTRONHF_H
#define ALIDIELECTRONHF_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#         Class AliDielectronHF                     #
//#                                                           #
//#  Authors:                                                 #
//#   Julian Book, Uni-Frankfurt / Julian.Book@cern.ch        #
//#                                                           #
//#############################################################

#include <TNamed.h>
#include <TObjArray.h>

#include "AliDielectronVarManager.h"

class AliDielectronHF : public TNamed {
public:
  enum { kMaxCuts=20 };
  enum EBinType {
    kStdBin=0,
    kBinToMax,
    kBinFromMin,
    kSymBin
  };
  enum EPairType {
    kOSonly=0,  kOSandLS, kOSandROT, kOSandMIX, kAll
  };

  AliDielectronHF();
  AliDielectronHF(const char*name, const char* title);

  virtual ~AliDielectronHF();

  void Init();
  void SetSignalsMC(TObjArray* array)    {fSignalsMC = array;}
  void SetStepForMCGenerated(Bool_t switcher=kTRUE)    {fStepGenerated = switcher;}
  void SetPairTypes(EPairType ptype=kOSonly) { fPairType=ptype; }
  void SetRefHist(TH1 *obj, UInt_t vars[4]);

  void AddCutVariable(AliDielectronVarManager::ValueTypes type, Int_t nbins,
		      Double_t min, Double_t max, Bool_t log=kFALSE, Bool_t leg=kFALSE, EBinType btype=kStdBin);
  void AddCutVariable(AliDielectronVarManager::ValueTypes type, 
		      const char* binLimitStr, Bool_t leg=kFALSE, EBinType btype=kStdBin);
  void AddCutVariable(AliDielectronVarManager::ValueTypes type, 
		      TVectorD * binLimits, Bool_t leg=kFALSE, EBinType btype=kStdBin);

  void Fill(Int_t pairIndex, const AliDielectronPair *particle);
  void Fill(Int_t label1, Int_t label2, Int_t nSignal);
  void Fill(Int_t Index, Double_t * const valuesPair, Double_t * const valuesLeg1, Double_t * const valuesLeg2);

  Bool_t IsPairTypeSelected(Int_t itype);

  Int_t GetNumberOfBins() const;
  const TObjArray * GetHistArray() const { return &fArrPairType; }
  Bool_t GetStepForMCGenerated() const   { return fStepGenerated; }

private:
  TObjArray fArrPairType;           //-> array of pair types
  EPairType fPairType;              // which pair combinations to include
  TObjArray* fSignalsMC;            //! array of MC signals to be studied

  UShort_t  fVarCuts[kMaxCuts];     // cut variables
  Bool_t    fVarCutType[kMaxCuts];  // array to store leg booleans
  TObjArray fAxes;                  // Axis descriptions of the cut binning
  UShort_t  fBinType[kMaxCuts];     // binning type of the axes
  
  Bool_t    fHasMC;
  Bool_t    fStepGenerated;         // switcher for generated particles

  TH1       *fRefObj;               // reference object

  AliDielectronHF(const AliDielectronHF &c);
  AliDielectronHF &operator=(const AliDielectronHF &c);

  
  ClassDef(AliDielectronHF,2)         // Dielectron HF
};



#endif
