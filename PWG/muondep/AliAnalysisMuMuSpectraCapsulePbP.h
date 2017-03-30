
#ifndef ALINANALYSISMUMUSPECTRACAPSULEPBP_H
#define ALINANALYSISMUMUSPECTRACAPSULEPBP_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// @ingroup pwg_muondep_mumu
/// @class AliAnalysisMuMuSpectraCapsulePbP
/// @brief helper class to deal with results stored in a spectra with pPb methods.
///
/// author : Benjamin Audurier (Subatech)

#include "TNamed.h"
#include "TMath.h"
#include <TString.h>
#include "TGraphErrors.h"
#include "AliAnalysisMuMuSpectra.h"
#include "AliCounterCollection.h"
#include "AliMergeableCollection.h"
#include "AliAnalysisMuMuSpectraCapsule.h"

class TGraphErrors;
class AliAnalysisMuMuSpectra;
class AliAnalysisMuMuSpectraCapsulePbP : public AliAnalysisMuMuSpectraCapsule
{

public:
  //ctor
  AliAnalysisMuMuSpectraCapsulePbP(
                             const AliAnalysisMuMuSpectra            *  spectra=0x0,
                             const TString                           spectraPath ="",
                             const char                              * externFile="",
                             const char                              * externFile2="");
  // dtor
  virtual ~AliAnalysisMuMuSpectraCapsulePbP();
  // Compute Yield
  TGraphErrors* ComputeYield(const char* what="", const TH1* histo=0x0, const char* sResName="");
  // Draw fit results and save them if wanted
  void DrawResults(const char* particle="PSI",const char* subresults="")const;
  // Print some data members
  void Print(Option_t* opt="") const;
  // Print constants used
  void PrintConst() const;
  // Compute R_pA
  TGraphErrors* RpAAsGraphic(Double_t MUL) const;


  // Return some data member. Double "const" on purpose to avoid leverage on data members
  const Double_t              * GetConstArray() const {return fConstArray;};
  const AliAnalysisMuMuSpectra* GetSpectra()    const {return fSpectra;};
  const TString                GetSpectraName() const {return fSpectraName;};

private:
  // Read and compute values from extern file
  void GetValuesFromExternFile(TString sbin, Double_t numArray[],Double_t MUL) const;
  // Set global constants according to centrality
  Bool_t SetConstantFromExternFile(const char* file);
  // Equality operator
  AliAnalysisMuMuSpectraCapsulePbP(const AliAnalysisMuMuSpectraCapsulePbP& rhs);// not implemented on purpose
  AliAnalysisMuMuSpectraCapsulePbP& operator=(const AliAnalysisMuMuSpectraCapsulePbP& rhs);// not implemented on purpose


private:

  Double_t fConstArray[13]; // Array to store constant according to centrality bins

 /// \cond CLASSIMP
 ClassDef(AliAnalysisMuMuSpectraCapsulePbP,2);
 /// \endcond
};







#endif
