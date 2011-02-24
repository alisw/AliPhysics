#ifndef ALIDIELECTRONSPECTRUM_H
#define ALIDIELECTRONSPECTRUM_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//#############################################################
//#                                                           # 
//#         Class AliDielectronSpectrum                       #
//#         Manage Cuts on the legs of the pair               #
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

#include <TObjArray.h>
#include <TString.h>

#include <TNamed.h>

class AliCFGridSparse;
class AliDielectronSignalBase;
class AliCFGridSparse;
class AliCFContainer;

class AliDielectronSpectrum : public TNamed {
public:
  AliDielectronSpectrum();
  AliDielectronSpectrum(const char*name, const char* title);

  virtual ~AliDielectronSpectrum();
  
  void AddMethod(AliDielectronSignalBase * const method) {fSignalMethods.Add(method);}

  void SetCorrectionContainer(AliCFContainer * const container, Int_t nominator, Int_t denominator);
  void SetSignalContainer(AliCFContainer * const container, Int_t step);
  
  void SetVariables(const char* vars)   { fVariables=vars; }

  void SetStepForSignal(Bool_t step=kTRUE)       { fStepSignal=step;       }
  void SetStepForSignificance(Bool_t step=kTRUE) { fStepSignificance=step; }
  void SetStepForSignalOverBackground(Bool_t step=kTRUE) { fStepSOB=step;  }
  void SetStepForMass(Bool_t step=kTRUE)         { fStepMass=step;         }
  void SetStepForMassWidth(Bool_t step=kTRUE)    { fStepMassWidth=step;    }
  
  void SetNoOwnerSpectrum(Bool_t noOwner=kTRUE) { fOwnerSpectrum=!noOwner; }

  void SetVisualDebug(Bool_t vis=kTRUE) { fVisualDebug=vis; }
  
  AliCFContainer* GetSpectrumContainer() const { return fCFSpectrum;   }
  AliCFGridSparse* GetCorrectionMatrix() const { return fCFCorrMatrix; }
  
  void Process();


private:
  AliCFContainer  *fCFSignal;               // CF container with from which to extract the Signal
  AliCFContainer  *fCFCorrection;           // CF container from which to extract the correction matrix
  AliCFContainer  *fCFSpectrum;             // CF container with extracted signal
  AliCFGridSparse *fCFCorrMatrix;           // correction matrix

  Bool_t fStepSignal;                       // if to create a step for the signal
  Bool_t fStepSignificance;                 // if to create a step for the significance
  Bool_t fStepSOB;                          // if to create a step for signal over background
  Bool_t fStepMass;                         // if to create a setp for the mass
  Bool_t fStepMassWidth;                    // if to create a setp for the mass width
  
  Int_t fSignalStep;                        // step to use from the signal container
  
  Int_t fCorrNom;                           // Nominator to use from corr matrix container
  Int_t fCorrDenom;                         // Deominator to use from corr matrix container
  
  TObjArray    fSignalMethods;              // array with signal extraction methods
  TString      fVariables;                  // variable names as a function of which to extract the signal

  Bool_t fOwnerSpectrum;                    // if we own the creted spectrum

  Bool_t fVisualDebug;                      // if we want to show the fit and print it
  
  Int_t        fNvars;                      //! number of variables
  Int_t        *fVars;                      //! variable numbers translated from fVariables
  Int_t        *fNbins;                     //! number of bins for each variable
  Int_t        *fCurrentBins;               //! bin currently selected for each variable
  Double_t     *fCurrentPositions;          //! variables values currently selected
  
  void Fill(Int_t *bin, Int_t step, Double_t value, Double_t error);
  Bool_t SetupVariables();
  void CreateCorrectionMatrix();
  void CreateCFSpectrum();
  void ExtractSignalInBins(Int_t variable=0);
  
  AliDielectronSpectrum(const AliDielectronSpectrum &c);
  AliDielectronSpectrum &operator=(const AliDielectronSpectrum &c);
  
  ClassDef(AliDielectronSpectrum,1)         //Cut class providing cuts for both legs of a pair
    
};

//
// inline functions
//
inline void AliDielectronSpectrum::SetCorrectionContainer(AliCFContainer * const container, Int_t nominator, Int_t denominator)
{
  fCFCorrection=container;
  fCorrNom=nominator;
  fCorrDenom=denominator;
}

inline void AliDielectronSpectrum::SetSignalContainer(AliCFContainer * const container, Int_t step)
{
  fCFSignal=container;
  fSignalStep=step;
}

#endif
