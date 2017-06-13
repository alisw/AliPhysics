
#ifndef ALINANALYSISMUMUSPECTRAPROCESSORPBPB_H
#define ALINANALYSISMUMUSPECTRAPROCESSORPBPB_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// @ingroup pwg_muondep_mumu
/// @class AliAnalysisMuMuSpectraCapsulePbPb
/// @brief helper class to deal with results stored in a spectra with PbPb methods.
///
/// author : Benjamin Audurier (Subatech)



#include "TNamed.h"
#include "TMath.h"
#include "TCanvas.h"
#include <TString.h>
#include "TGraphErrors.h"
#include "AliAnalysisMuMuSpectra.h"
#include "AliCounterCollection.h"
#include "AliMergeableCollection.h"
#include "AliAnalysisMuMuSpectraProcessor.h"

class TGraphErrors;
class AliAnalysisMuMuSpectra;
class AliAnalysisMuMuSpectraProcessorPbPb : public AliAnalysisMuMuSpectraProcessor
{

public:
  //ctor
  AliAnalysisMuMuSpectraProcessorPbPb(
                             const AliAnalysisMuMuSpectra            *  spectra=0x0,
                             const TString                           spectraPath ="",
                             const char                              * externFile="",
                             const char                              * externFile2="");
  // dtor
  virtual ~AliAnalysisMuMuSpectraProcessorPbPb();
  // Compute Yield
  TGraphErrors* ComputeYield(const char* what="", const TH1* histo=0x0, const char* sResName="",Double_t MUL=0.);
  // Draw fit results and save them if wanted
  void DrawResults(const char* what="NofJPsi", const char* particle="PSI",const char* subresults="")const;
  // Print Flag
  void SetPrintFlag(){fPrintFlag=kTRUE;};
  // Style for canvas
  void SetCanvasStyle(TCanvas *can) const ;
  // Print some data members
  void Print(Option_t* opt="") const;
  // Print constants used
  void PrintConst() const;
  // Compute quantities linked to RAA
  TList* RAAasGraphic(Double_t MUL) const;


  // Return some data member. Double "const" on purpose to avoid leverage on data members
  const Double_t              * GetConstArray()     const {return fConstArray;};
  const AliAnalysisMuMuSpectra* GetSpectra()       const {return fSpectra;};
  const TString                GetSpectraName()    const {return fSpectraName;};

private:
  // Read and compute values from extern file
  Bool_t ComputeRAA(TString sbin, Double_t numArray[],Double_t MUL, Double_t binwidth) const;
  // Read extern file for Pt and Y case
  Bool_t ReadFromFile(TString sbin, float valueArray[]) const;
  // Set global constants according to centrality
  Bool_t SetConstantFromExternFile(const char* file);
  // Equality operator
  AliAnalysisMuMuSpectraProcessorPbPb(const AliAnalysisMuMuSpectraProcessorPbPb& rhs);// not implemented on purpose
  AliAnalysisMuMuSpectraProcessorPbPb& operator=(const AliAnalysisMuMuSpectraProcessorPbPb& rhs);// not implemented on purpose


private:
  Double_t fConstArray[13]; // Array to store constant according to centrality bins

  /// \cond CLASSIMP
  ClassDef(AliAnalysisMuMuSpectraProcessorPbPb,1);
 /// \endcond
};







#endif
