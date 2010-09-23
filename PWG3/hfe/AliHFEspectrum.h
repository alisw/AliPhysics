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
#ifndef ALIHFESPECTRUM_H
#define ALIHFESPECTRUM_H

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class TGraphErrors;
class TObject;
class TH1;
class TList;
class AliCFContainer;
class AliCFDataGrid;
class AliCFEffGrid;

class AliHFEspectrum : public TNamed{
  public:
    enum CFContainer_t{
      kDataContainer = 0,
      kMCContainer = 1,
      kBackgroundData = 2,
      kBackgroundMC = 3
    };
    enum BackgroundSource_t{
      kMCbackground = 0,
      kDataBackground = 1
    };
    AliHFEspectrum(const char* name);
    ~AliHFEspectrum();

    void Correct(AliCFContainer *datacontainer,AliCFContainer *mccontainer,THnSparseF *mccorrelation,AliCFContainer *contaminationcontainer=0x0);
    AliCFDataGrid *SubtractBackground(Int_t dimensions, Bool_t setBackground = kFALSE);
    AliCFDataGrid *TakeBackgroundFromData(Int_t nDim);    

    TList *Unfold(Int_t dimensions, AliCFDataGrid* const bgsubpectrum = 0x0);
    AliCFDataGrid *CorrectForEfficiency(Int_t dimensions, AliCFDataGrid* const bgsubpectrum = 0x0);
 
    TGraphErrors *Normalize(THnSparse * const spectrum) const;
    TGraphErrors *Normalize(AliCFDataGrid * const spectrum) const;
 
    void SetCorrelation(THnSparseF * const correlation) {fCorrelation = correlation; };
    void SetContainer(AliCFContainer *cont, AliHFEspectrum::CFContainer_t type);
    void SetNumberOfEvents(Int_t nEvents) { fNEvents = nEvents; }
    void SetBackgroundSource(BackgroundSource_t source) { fBackgroundSource = source; };
    void SetMCEffStep(Int_t step) { fStepMC = step; };
    void SetMCTruthStep(Int_t step) { fStepTrue = step; };
    void SetStepToCorrect(Int_t step) { fStepData = step; };
    void SetStepGuessedUnfolding(Int_t stepGuessedUnfolding) { fStepGuessedUnfolding = stepGuessedUnfolding; };
    void SetNumberOfIteration(Int_t numberOfIteration) { fNumberOfIterations = numberOfIteration; };
    void SetDumpToFile(Bool_t dumpToFile) { fDumpToFile=dumpToFile; }; 
  

  protected:
    AliCFDataGrid *MakeBackgroundEstimateFromMC(Int_t nDimensions);
    
    AliCFContainer *GetContainer(AliHFEspectrum::CFContainer_t contt);
    AliCFContainer *GetSlicedContainer(AliCFContainer *cont, Int_t ndim, Int_t *dimensions);
    THnSparse *GetSlicedCorrelation(Int_t nDim, Int_t *dimensions) const;
    TObject* GetSpectrum(AliCFContainer * const c, Int_t step);
    TObject* GetEfficiency(AliCFContainer * const c, Int_t step, Int_t step0);
 
    void AddTemporaryObject(TObject *cont);
    void ClearObject(TObject *o);
    
    void CorrectFromTheWidth(TH1D *h1) const;
    TGraphErrors *NormalizeTH1(TH1 *input) const;


  private:
    AliHFEspectrum(const AliHFEspectrum &);
    AliHFEspectrum &operator=(const AliHFEspectrum &);
 
    TList *fCFContainers;         // List of Correction Framework Containers
    TList *fTemporaryObjects;     // Emulate garbage collection
    THnSparseF *fCorrelation;     // Correlation Matrices
    AliCFDataGrid *fBackground;   // Background Grid

    BackgroundSource_t fBackgroundSource;   // Source for the background estimate

    Bool_t fDumpToFile;           // Write Result in a file

    Int_t fNEvents;               // Number of Events
    Int_t fStepMC;                // MC step (for unfolding)
    Int_t fStepTrue;              // MC step of the final spectrum
    Int_t fStepData;              // Data Step (various applications)
    Int_t fStepGuessedUnfolding;  // Step for first guessed unfolding
    Int_t fNumberOfIterations;    // Number of iterations


    ClassDef(AliHFEspectrum, 1) 
};
#endif

