/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef ALIFLOWCUMUCONSTANTS_H
#define ALIFLOWCUMUCONSTANTS_H

#include <TROOT.h>

// Description: constants for the LYZ flow makers
// Author: Ante Bilandzic, Nikhef
// modified: Mikolaj Krzewicki, Nikhef (mikolaj.krzewicki@cern.ch)

#include <TNamed.h>

class AliFlowCumuConstants : public TNamed {

public:
  AliFlowCumuConstants();
  virtual ~AliFlowCumuConstants();
  static AliFlowCumuConstants* GetMaster();

  //static consts for compatibility only!
  static const Int_t fgkQmax = 11; //
  static const Int_t fgkPmax = 5; //
  static const Int_t fgkQmax4 = 5; //
  static const Int_t fgkPmax4 = 2; //
  static const Int_t fgkQmax6 = 7; //
  static const Int_t fgkPmax6 = 3; //
  static const Int_t fgkQmax8 = 9; //
  static const Int_t fgkPmax8 = 4; //
  static const Int_t fgkQmax16 = 17; //
  static const Int_t fgkPmax16 = 8; //
  static const Int_t fgkFlow = 2; //
  static const Int_t fgkMltpl = 1; //
  //static const Double_t  fgkBinWidth = 0.1; //
  //static const Double_t  fgkR0 = 2.2; //
  //static const Double_t  fgkPtMax = 3.1; //
  //static const Double_t  fgkPtMin = 0.0; //
  //static const Bool_t  fgkOtherEquations = kFALSE; //
  
  Int_t GetQmax() {return fQmax;}
  Int_t GetPmax() {return fPmax;}
  Int_t GetQmax4() {return fQmax4;}
  Int_t GetPmax4() {return fPmax4;}
  Int_t GetQmax6() {return fQmax6;}
  Int_t GetPmax6() {return fPmax6;}
  Int_t GetQmax8() {return fQmax8;}
  Int_t GetPmax8() {return fPmax8;}
  Int_t GetQmax16() {return fQmax16;}
  Int_t GetPmax16() {return fPmax16;}
  Int_t GetFlow() {return fFlow;}
  Int_t GetMltpl() {return fMltpl;}
  Double_t GetBinWidth() {return fBinWidth;}
  Double_t GetR0() {return fR0;}
  Double_t GetPtMax() {return fPtMax;}
  Double_t GetPtMin() {return fPtMin;}
  Bool_t GetOtherEquations() {return fOtherEquations;}
  
  void SetQmax( Int_t i ) {fQmax = i;}
  void SetPmax( Int_t i ) {fPmax = i;}
  void SetQmax4( Int_t i ) {fQmax4 = i;}
  void SetPmax4( Int_t i ) {fPmax4 = i;}
  void SetQmax6( Int_t i ) {fQmax6 = i;}
  void SetPmax6( Int_t i ) {fPmax6 = i;}
  void SetQmax8( Int_t i ) {fQmax8 = i;}
  void SetPmax8( Int_t i ) {fPmax8 = i;}
  void SetQmax16( Int_t i ) {fQmax16 = i;}
  void SetPmax16( Int_t i ) {fPmax16 = i;}
  void SetFlow( Int_t i ) {fFlow = i;}
  void SetMltpl( Int_t i ) {fMltpl = i;}
  void SetBinWidth( Double_t i ) {fBinWidth = i;}
  void SetR0( Double_t i ) {fR0 = i;}
  void SetPtMax( Double_t i ) {fPtMax = i;}
  void SetPtMin( Double_t i ) {fPtMin = i;}
  void SetOtherEquations( Bool_t i ) {fOtherEquations = i;}
  
private:
  AliFlowCumuConstants& operator= (const AliFlowCumuConstants& c);
  AliFlowCumuConstants(const AliFlowCumuConstants& a);

  Int_t fQmax; //
  Int_t fPmax; //
  Int_t fQmax4; //
  Int_t fPmax4; //
  Int_t fQmax6; //
  Int_t fPmax6; //
  Int_t fQmax8; //
  Int_t fPmax8; //
  Int_t fQmax16; //
  Int_t fPmax16; //
  Int_t fFlow; //
  Int_t fMltpl; //

  // Histograms limits
  Double_t  fBinWidth; //
  Double_t  fR0; //
  Double_t  fPtMax; //
  Double_t  fPtMin; //
  
  // Other numerical equations for cumulants 
  Bool_t  fOtherEquations; //

  static AliFlowCumuConstants* fgPMasterConfig;

  ClassDef(AliFlowCumuConstants,1)  // macro for rootcint
};

#endif

