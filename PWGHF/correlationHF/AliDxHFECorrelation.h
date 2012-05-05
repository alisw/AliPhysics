//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliDxHFECorrelation.h
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-04-25
/// @brief  Worker class for D0-HF electron correlation
///

#ifndef ALIDXHFECORRELATION_H
#define ALIDXHFECORRELATION_H

#include "TNamed.h"

class TH1F;
class TH2F;

class AliDxHFECorrelation : public TNamed {
 public:
  /// default constructor
  AliDxHFECorrelation(const char* name=NULL);
  /// destructor
  virtual ~AliDxHFECorrelation();

  // init
  int Init();

  /// fill histograms from particles
  int Fill(const TObjArray* candidatesD0, const TObjArray* candidatesElectron);

  /// overloaded from TObject: cleanup
  virtual void Clear(Option_t * option ="");
  /// overloaded from TObject: print info
  virtual void Print(Option_t *option="") const;
  /// overloaded from TObject: draw histograms
  virtual void Draw(Option_t *option="");
  /// overloaded from TObject: find object by name
  virtual TObject* FindObject(const char *name) const;
  /// overloaded from TObject: find object by pointer
  virtual TObject* FindObject(const TObject *obj) const;
  /// overloaded from TObject: save to file
  virtual void     SaveAs(const char *filename="",Option_t *option="") const; // *MENU*

  AliDxHFECorrelation& operator+=(const AliDxHFECorrelation& other);

  enum {
    khD0pT,         // TH1F
    khD0Phi,        // TH1F
    khD0Eta,        // TH1F
    khElectronpT,   // TH1F
    khElectronPhi,  // TH1F
    khElectronEta,  // TH1F
    khDeltaPhi,        // TH1F
    kNofHistograms
  };

 protected:
 private:
  /// copy constructor
  AliDxHFECorrelation(const AliDxHFECorrelation& other);
  /// assignment operator
  AliDxHFECorrelation& operator=(const AliDxHFECorrelation& other);

  TObjArray* fHistograms; // the histograms

  ClassDef(AliDxHFECorrelation, 1)
};
#endif
