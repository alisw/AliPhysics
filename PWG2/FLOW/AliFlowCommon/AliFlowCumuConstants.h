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
  static const Int_t fgkQmax = 11; // Qmax 
  static const Int_t fgkPmax = 5; // Pmax 
  static const Int_t fgkQmax4 = 5; // Qmax4
  static const Int_t fgkPmax4 = 2; // Pmax4
  static const Int_t fgkQmax6 = 7; // Qmax6
  static const Int_t fgkPmax6 = 3; // Pmax6
  static const Int_t fgkQmax8 = 9; // Qmax8
  static const Int_t fgkPmax8 = 4; // Pmax8
  static const Int_t fgkQmax16 = 17; // Qmax16
  static const Int_t fgkPmax16 = 8; // Pmax16
  static const Int_t fgkFlow = 2; // flow estimate
  static const Int_t fgkMltpl = 1; // multiple
  //static const Double_t  fgkBinWidth = 0.1; //
  //static const Double_t  fgkR0 = 2.2; //
  //static const Double_t  fgkPtMax = 3.1; //
  //static const Double_t  fgkPtMin = 0.0; //
  //static const Bool_t  fgkOtherEquations = kFALSE; //
  
  Int_t GetQmax() const {return fQmax;}
  Int_t GetPmax() const {return fPmax;}
  Int_t GetQmax4() const {return fQmax4;}
  Int_t GetPmax4() const {return fPmax4;}
  Int_t GetQmax6() const {return fQmax6;}
  Int_t GetPmax6() const {return fPmax6;}
  Int_t GetQmax8() const {return fQmax8;}
  Int_t GetPmax8() const {return fPmax8;}
  Int_t GetQmax16() const {return fQmax16;}
  Int_t GetPmax16() const {return fPmax16;}
  Int_t GetFlow() const {return fFlow;}
  Int_t GetMltpl() const {return fMltpl;}
  Double_t GetBinWidth() const {return fBinWidth;}
  Double_t GetR0() const {return fR0;}
  Double_t GetPtMax() const {return fPtMax;}
  Double_t GetPtMin() const {return fPtMin;}
  Bool_t GetOtherEquations() const {return fOtherEquations;}
  
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

  Int_t fQmax; // Qmax
  Int_t fPmax; // Pmax
  Int_t fQmax4; // Qmax4
  Int_t fPmax4; // Pmax4
  Int_t fQmax6; // Qmax6
  Int_t fPmax6; // Pmax6
  Int_t fQmax8; // Qmax8
  Int_t fPmax8; // Pmax8
  Int_t fQmax16; // Qmax16
  Int_t fPmax16; // Pmax16
  Int_t fFlow; // flow estimate
  Int_t fMltpl; // multiple

  // Histograms limits
  Double_t  fBinWidth; // bin width
  Double_t  fR0; // R0
  Double_t  fPtMax; // Pt max
  Double_t  fPtMin; // Pt min
  
  // Other numerical equations for cumulants 
  Bool_t  fOtherEquations; // othe equations

  static AliFlowCumuConstants* fgPMasterConfig; // master config

  ClassDef(AliFlowCumuConstants,1)  // macro for rootcint
};

#endif

