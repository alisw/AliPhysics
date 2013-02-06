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
// Class for spectrum correction
// Subtraction of hadronic background, Unfolding of the data and
// Renormalization done here
// For more information see the implementation file
//
#ifndef ALIHFECORRECTSPECTRUMBASE_H
#define ALIHFECORRECTSPECTRUMBASE_H

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class TGraphErrors;
class TObject;
class TH1;
class TF1;
class TList;
class TObjArray;
class AliCFContainer;
class AliHFEcontainer;
class AliCFDataGrid;
class AliCFEffGrid;

class AliHFECorrectSpectrumBase : public TNamed{
  public:
    enum CFContainer_t{
      kDataContainer = 0,
      kBackgroundData = 1,
      kMCContainerMC = 2,
      kMCContainerESD = 3,
      kMCContainerCharmMC = 4,
      kMCWeightedContainerNonHFEESD = 5,
      kMCWeightedContainerConversionESD = 6,
      kDataContainerV0 = 7,
      kMCWeightedContainerNonHFEESDSig = 8,
      kMCWeightedContainerConversionESDSig = 9,
      kPhotonicBackground = 10,
      kNbCFContainers = 11
    };

    enum Chargetype_t{
      kNegCharge = -1,
      kPosCharge = 1,
      kAllCharge = 0
    };
   
    AliHFECorrectSpectrumBase(const char* name);
    ~AliHFECorrectSpectrumBase();
    

    virtual Bool_t Init(const AliHFEcontainer */*datahfecontainer*/, const AliHFEcontainer */*mchfecontainer*/, const AliHFEcontainer */*bghfecontainer*/, const AliHFEcontainer */*v0hfecontainer*/,AliCFContainer */*photoniccontainerD*/) { return kTRUE;};
    virtual Bool_t Correct(Bool_t /*subtractcontamination*/, Bool_t /*subtractphotonic*/) { return kTRUE;};
   
    TGraphErrors *Normalize(THnSparse * const spectrum) const;
    TGraphErrors *Normalize(AliCFDataGrid * const spectrum) const;
    TGraphErrors *NormalizeTH1(TH1 *input) const;
    void CorrectFromTheWidth(TH1D *h1) const;
    void CorrectStatErr(AliCFDataGrid *backgroundGrid) const;
    
    void SetCorrelation(THnSparseF * const correlation) {fCorrelation = correlation; };
    void SetContainer(AliCFContainer *cont, AliHFECorrectSpectrumBase::CFContainer_t type);
    void SetEfficiencyFunction(TF1 *efficiencyFunction) { fEfficiencyFunction = efficiencyFunction; };
    
    void SetNumberOfEvents(Int_t nEvents) { fNEvents = nEvents; };
    void SetMCEffStep(Int_t step) { fStepMC = step; };
    void SetMCTruthStep(Int_t step) { fStepTrue = step; };
    void SetStepToCorrect(Int_t step) { fStepData = step; };
    void SetStepBeforeCutsV0(Int_t step) { fStepBeforeCutsV0 = step; };
    void SetStepAfterCutsV0(Int_t step) { fStepAfterCutsV0 = step; };
    void SetNbDimensions(Int_t nbDimensions) { fNbDimensions = nbDimensions; };
    void SetChargeChoosen(Chargetype_t chargechoosen) {fChargeChoosen = chargechoosen; };
    void SetEtaRange(Double_t etamin, Double_t etamax) { fEtaRange[0] = etamin; fEtaRange[1] = etamax; fEtaSelected = kTRUE; }
    void SetSmoothing(Bool_t setSmoothing) {fSetSmoothing = setSmoothing;};
    void SetTestOneBinCentrality(Double_t centralitymin, Double_t centralitymax) { fTestCentralityLow = centralitymin; fTestCentralityHigh = centralitymax;}
    void SetStepGuessedUnfolding(Int_t stepGuessedUnfolding) { fStepGuessedUnfolding = stepGuessedUnfolding; };
    void SetNumberOfIteration(Int_t numberOfIteration) { fNumberOfIterations = numberOfIteration; };
    
    

 protected:
    AliHFECorrectSpectrumBase(const AliHFECorrectSpectrumBase &ref);
    AliHFECorrectSpectrumBase &operator=(const AliHFECorrectSpectrumBase &ref);
    virtual void Copy(TObject &o) const;
    AliCFContainer *GetContainer(AliHFECorrectSpectrumBase::CFContainer_t contt);
    AliCFContainer *GetSlicedContainer(AliCFContainer *cont, Int_t ndim, Int_t *dimensions,Int_t source=-1,Chargetype_t charge=kAllCharge,Int_t centralitylow=-1, Int_t centralityhigh=-1);
    THnSparseF *GetSlicedCorrelation(THnSparseF *correlationmatrix,Int_t nDim, Int_t *dimensions,Chargetype_t charge=kAllCharge,Int_t centralitylow=-1, Int_t centralityhigh=-1) const;
    TObject* GetSpectrum(const AliCFContainer * const c, Int_t step);
    TObject* GetEfficiency(const AliCFContainer * const c, Int_t step, Int_t step0);

    TObjArray *fCFContainers;     // List of Correction Framework Containers
    THnSparseF *fCorrelation;     // Correlation Matrices
    TF1 *fEfficiencyFunction;     // Efficiency Function
   
    Bool_t fEtaSelected;              // Switch for eta selection
    Bool_t fSetSmoothing;             // Set smoothing

    Int_t fNbDimensions;          // Number of dimensions for the correction
    Int_t fNEvents;               // Number of Events
    Int_t fStepMC;                // MC step (for unfolding)
    Int_t fStepTrue;              // MC step of the final spectrum
    Int_t fStepData;              // Data Step (various applications)
    Int_t fStepBeforeCutsV0;      // Before cuts V0
    Int_t fStepAfterCutsV0;       // After cuts V0
    Int_t fStepGuessedUnfolding;  // Step for first guessed unfolding
    Int_t fNumberOfIterations;    // Number of iterations
    Chargetype_t fChargeChoosen;         // Select positive or negative electrons

    Double_t fEtaRange[2];        // Eta range 
    Double_t fEtaRangeNorm[2];    // Eta range used in the normalization

    Int_t fTestCentralityLow;     // To test one bin in centrality only
    Int_t fTestCentralityHigh;    // To test one bin in centrality only
      


  private:
   
 
   
    ClassDef(AliHFECorrectSpectrumBase, 1) 
};
#endif

