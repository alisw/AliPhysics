
#ifndef ALINANALYSISMUMUSPECTRACAPSULEPP_H
#define ALINANALYSISMUMUSPECTRACAPSULEPP_H


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
#include <TString.h>
#include "TGraphErrors.h"
#include "AliAnalysisMuMuSpectra.h"
#include "AliCounterCollection.h"
#include "AliMergeableCollection.h"
#include "AliAnalysisMuMuSpectraCapsule.h"

class TGraphErrors;
class AliAnalysisMuMuSpectra;
class AliAnalysisMuMuSpectraCapsulePP : public AliAnalysisMuMuSpectraCapsule
{

public:
  //ctor
  AliAnalysisMuMuSpectraCapsulePP(
                             const AliAnalysisMuMuSpectra            *  spectra=0x0,
                             const TString                           spectraPath ="",
                             const char                              * externFile="",
                             const char                              * externFile2="");
  // dtor
  virtual ~AliAnalysisMuMuSpectraCapsulePP();
  // Compute Yield
  TGraphErrors* ComputeYield(const char* what="", const TH1* histo=0x0, const char* sResName="",Double_t MUL=0.);
  // Compute Cross-Section
  TList* ComputeJpsiPPCrossSection(const char* what ="CorrNofJPsi") const ;
  // Draw fit results and save them if wanted
  void DrawResults(const char* particle="PSI",const char* subresults="")const;
  // Print Flag
  void SetPrintFlag(){fPrintFlag=kTRUE;};
  // Print some data members
  void Print(Option_t* opt="") const;
  // Print constants used
  void PrintConst() const;


  // Return some data member. Double "const" on purpose to avoid leverage on data members
  const AliAnalysisMuMuSpectra* GetSpectra()       const {return fSpectra;};
  const TString                GetSpectraName()    const {return fSpectraName;};
   const Double_t              * GetConstArray()     const {return fConstArray;};

private:
  // Read extern file for Pt and Y case
  Bool_t ReadFromFile(TString sbin, float valueArray[]) const;
  // Equality operator
  AliAnalysisMuMuSpectraCapsulePP(const AliAnalysisMuMuSpectraCapsulePP& rhs);// not implemented on purpose
  AliAnalysisMuMuSpectraCapsulePP& operator=(const AliAnalysisMuMuSpectraCapsulePP& rhs);// not implemented on purpose


private:

  Double_t      fConstArray[13]; // Array to store constant according to centrality bins

  /// \cond CLASSIMP
  ClassDef(AliAnalysisMuMuSpectraCapsulePP,3);
  /// \endcond
};







#endif
