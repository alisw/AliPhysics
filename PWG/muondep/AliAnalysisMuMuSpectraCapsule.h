
#ifndef ALINANALYSISMUMUSPECTRACAPSULE_H
#define ALINANALYSISMUMUSPECTRACAPSULE_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// @ingroup pwg_muondep_mumu
/// @class AliAnalysisMuMuSpectraCapsule
/// @bried Mother class to all capsule class.
///
/// author : Benjamin Audurier (Subatech)


#include <TString.h>
#include "TObject.h"
#include "TGraphErrors.h"
#include "AliAnalysisMuMuSpectra.h"
#include "AliCounterCollection.h"
#include "AliMergeableCollection.h"

class TGraphErrors;
class AliAnalysisMuMuSpectraCapsule : public TObject
{

public:
  //ctor
  AliAnalysisMuMuSpectraCapsule(
    const AliAnalysisMuMuSpectra            *  spectra=0x0,
    const TString                           spectraPath ="",
    const char                              * externFile="",
    const char                              * externFile2="");
  //dctor
  virtual ~AliAnalysisMuMuSpectraCapsule();
  // Compute Yield
  virtual TGraphErrors* ComputeYield(const char* what, const TH1* histo, const char* sResName, Double_t MUL) = 0;
  // Print some data members
  virtual void Print(Option_t* opt) const = 0;
  // Print constants used
  virtual void PrintConst() const = 0;
  // Number of "what" for all subresults
  void PrintNofWhat(const char* what="") const;
  // Set global constants according to centrality
  Bool_t SetConstantFromExternFile(const char* file, Double_t* constantArray, const TString* spectraName);

  // Protection to be sure each daughter class returns data members
  const virtual Double_t              * GetConstArray()   const =0;

  const virtual AliAnalysisMuMuSpectra* GetSpectra()      const =0;

  const virtual TString                 GetSpectraName() const  =0;

  private:
  // Equality operator
  AliAnalysisMuMuSpectraCapsule(const AliAnalysisMuMuSpectraCapsule& rhs);// not implemented on purpose
  AliAnalysisMuMuSpectraCapsule& operator=(const AliAnalysisMuMuSpectraCapsule& rhs);// not implemented on purpose

  protected:
  const AliAnalysisMuMuSpectra* fSpectra;     // Spectra with result and subresults
  TString                       fExternFile;  // name of spectra selected
  TString                       fExternFile2; // name of spectra selected
  Bool_t                        fPrintFlag;
  const TString                 fSpectraName; // SpectraName




/// \cond CLASSIMP
ClassDef(AliAnalysisMuMuSpectraCapsule,2)
/// \endcond;
};




#endif
