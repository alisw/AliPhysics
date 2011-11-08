/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Class AliHFEvarManager
// Common place for definiton of variables to be filled into the 
// correction framework container
// More information can be found inside the implementation file
//
#ifndef ALIHFEVARMANAGER_H
#define ALIHFEVARMANAGER_H

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

template <class X>
class THnSparseT;
class TArrayF;
typedef class THnSparseT<TArrayF> THnSparseF;
class TH3F;
class TF3;
class AliCFContainer;
class AliHFEcontainer;
class AliHFEsignalCuts;
class AliVParticle;


class AliHFEvarManager : public TNamed{
public:
  enum EVarCode_t{
    kPt = 1,
    kEta,
    kPhi,
    kCharge,
    kSource,
    kCentrality,
    kSpecies,
  };
	AliHFEvarManager();
	AliHFEvarManager(const Char_t *name);
  AliHFEvarManager(const AliHFEvarManager &ref);
  AliHFEvarManager &operator=(const AliHFEvarManager &ref);
  void Copy(TObject &o) const;
	~AliHFEvarManager();

  TObjArray *GetVariables() const { return fVariables; }
  
  void SetOwner(Bool_t owner = kTRUE) { SetBit(kOwner, owner); }
  Bool_t IsOwner() const { return TestBit(kOwner); }

  void AddVariable(TString name);
  void AddVariable(TString name, Int_t nBins, Double_t min, Double_t max, Bool_t isLogarithmic = kFALSE);
  void AddVariable(TString name, Int_t nBins, const Double_t *binning);
  Bool_t IsVariableDefined(TString name);
  void DefineVariables(AliHFEcontainer *cont);
  void NewTrack(AliVParticle *track, AliVParticle *mcTrack = NULL, Float_t centrality = 99.0, Int_t aprioriPID = -1, Bool_t signal = kTRUE);
  Bool_t IsSignalTrack() const { return fSignalTrack; }
  void FillContainer(AliCFContainer *const cont, Int_t step, Bool_t useMC = kFALSE) const;
  void FillContainer(const AliHFEcontainer *const cont, const Char_t *contname, UInt_t step, Bool_t useMC = kFALSE, Double_t externalWeight = 1.) const;
  void FillContainerStepname(const AliHFEcontainer *const cont, const Char_t *contname, const Char_t *step, Bool_t useMC = kFALSE, Double_t externalWeight = 1.) const;
  void FillCorrelationMatrix(THnSparseF *matrix) const;
  
  void SetSignalCuts(AliHFEsignalCuts *signal) { fSignal = signal; }
  void SetWeightFactors(TH3F *weightFactors);
  void SetWeightFactorsFunction(TF3*weightFactorsFunction);
  
  struct AliHFEvariable : public TNamed{
    public:
      AliHFEvariable();
      AliHFEvariable(const Char_t *name, const Char_t *title, UInt_t fCode, UInt_t nBins, Double_t min, Double_t max, Bool_t isLogarithmic = kFALSE);
      AliHFEvariable(const Char_t *name, const Char_t *title, UInt_t fCode, UInt_t nBins, const Double_t *binning);
      AliHFEvariable(const AliHFEvariable &ref);
      AliHFEvariable &operator=(const AliHFEvariable &ref);
      ~AliHFEvariable();

      UInt_t GetVarCode() const { return fCode; }
      UInt_t GetNumberOfBins() const { return fNBins; }
      const Double_t* GetBinning() const { return fBinning; }
      Bool_t HasUserDefinedBinning() const { return fBinning != NULL; }
      Double_t GetMinimum() const { return fMin; }
      Double_t GetMaximum() const { return fMax; } 
      Int_t IsLogarithmic() const { return fIsLogarithmic; }
    private:
      UInt_t    fCode;              // Unique variable code
      UInt_t    fNBins;             // Number of bins
      Double_t  fMin;               // Minimum
      Double_t  fMax;               // Maximum
      Double_t *fBinning;           // User defined binning
      Bool_t    fIsLogarithmic;     // Logarithmic binning

      ClassDef(AliHFEvarManager::AliHFEvariable, 1) // HFE variable definition
  };
 
protected:
  Double_t GetValue(AliVParticle *track, UInt_t code, Float_t centrality = 99.0, Int_t aprioriPID = -1) const;
  void FillArray(AliVParticle *track, Double_t *container, Float_t centrality = 99.0, Int_t aprioriPID = -1) const;
  Double_t FindWeight(Double_t pt, Double_t eta, Double_t phi) const;

private:
  enum{
    kOwner = BIT(14)
  };
  TObjArray *fVariables;                // Variables to process
  Double_t *fContent;                   //! Cache values for track in classmember 
  Double_t *fContentMC;                 //! Cache content of the asssociated MC track in class member
  Double_t fWeightFactor;               // Cache weighting factor
  Bool_t  fSignalTrack;                 // Signal Track
	Bool_t  fWeighting;                   // Weighting or not for the efficiency maps
  AliHFEsignalCuts *fSignal;            // MC Signal Definition
	TH3F *fWeightFactors;                 // Weight factors
	TF3  *fWeightFactorsFunction;         // Weight factors

	ClassDef(AliHFEvarManager, 1)         // The variable Manager for the HFE Analysis
};

#endif /* ALIHFEVARMANAGER_H */
