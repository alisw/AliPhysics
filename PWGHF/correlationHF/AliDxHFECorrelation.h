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

class AliRDHFCutsD0toKpi;
class TH1;
class THnSparse;
class TObject;
class TList;

class AliDxHFECorrelation : public TNamed {
 public:
  /// default constructor
  AliDxHFECorrelation(const char* name=NULL);
  /// destructor
  virtual ~AliDxHFECorrelation();

  // event control histogram
  enum {
    kEventsAll = 0, // all events
    kEventsSel,     // selected events
    kEventsD0 ,     // events with D0s
    kEventsD0e,     // events with correlated D0s
    kNEventControlLabels
  };

  // init
  int Init();

  /// fill histograms from particles
  int Fill(const TObjArray* candidatesD0, const TObjArray* candidatesElectron);

  /// histogram event properties
  virtual int HistogramEventProperties(int bin);

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

  virtual void SetCuts(AliRDHFCutsD0toKpi* cuts) {fCuts=cuts;}
  virtual void SetUseMC(Bool_t useMC){fUseMC=useMC;}

  Bool_t GetUseMC() const {return fUseMC;}
  const TList* GetControlObjects() const {return fControlObjects;}


  AliDxHFECorrelation& operator+=(const AliDxHFECorrelation& other);


  // Probably not needed anymore, since code was changed to THnSparse
  // but keep here in case we need it later
  enum {
    khD0pT,         // TH1F
    khD0Phi,        // TH1F
    khD0Eta,        // TH1F
    khElectronpT,   // TH1F
    khElectronPhi,  // TH1F
    khElectronEta,  // TH1F
    kNofHistograms
  };

 protected:
  /// add control object to list, the base class becomes owner of the object
  int AddControlObject(TObject* pObj);

 private:
  /// copy constructor
  AliDxHFECorrelation(const AliDxHFECorrelation& other);
  /// assignment operator
  AliDxHFECorrelation& operator=(const AliDxHFECorrelation& other);

  TObjArray* fHistograms;     // the histograms - for the moment not in use. 
  TList* fControlObjects;     // list of control objects
  THnSparse* fCorrProperties; // the Correlation properties of selected particles
  TH1* fhEventControlCorr;    // event control histogram
  AliRDHFCutsD0toKpi *fCuts;  //  Cuts 
  Bool_t fUseMC;              // use MC info

  static const char* fgkEventControlBinNames[];
  static const char* fgkCorrControlBinNames[];

  ClassDef(AliDxHFECorrelation, 2)
};
#endif
