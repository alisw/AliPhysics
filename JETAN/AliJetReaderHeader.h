#ifndef ALIJETREADERHEADER_H
#define ALIJETREADERHEADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// base class for Jet Reader Header 
//
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------
  
#include <TNamed.h>
#include <TString.h>
 
class AliJetReaderHeader : public TNamed
{

 public:
  AliJetReaderHeader(const char* name);
  AliJetReaderHeader();
  virtual ~AliJetReaderHeader();
  
  // Getters
  virtual const TString GetComment() {return fComment;}
  virtual const char* GetDirectory() {return fDir.Data();}
  virtual const char* GetPattern() {return fPattern.Data();}
  virtual Float_t     GetFiducialEtaMin() const {return fFiducialEtaMin;}
  virtual Float_t     GetFiducialEtaMax() const {return fFiducialEtaMax;}  
  virtual Float_t     GetPtCut() const {return fPtCut;}
  Int_t   GetNEvents()    const {return fLast-fFirst;}
  Int_t   GetFirstEvent() const {return fFirst;}
  Int_t   GetLastEvent()  const {return fLast;}
  Int_t   GetDetector()   const {return fOption;}
  Int_t   GetDebug()      const {return fDebug;}
  // Setters
  virtual void SetComment(const char* s) {fComment=TString(s);}
  virtual void SetPattern(const char* s) {fPattern=TString(s);}
  virtual void SetDirectory(const char* s) {fDir=TString(s);}
  virtual void SetFirstEvent(Int_t i=0) {fFirst=i;}
  virtual void SetLastEvent(Int_t i=-1) {fLast=i;}
  virtual void SetFiducialEta(Float_t etamin, Float_t etamax) 
      { fFiducialEtaMin = etamin; fFiducialEtaMax = etamax;}
  virtual void SetPtCut(Float_t par = 2.0) {fPtCut = par;}

  virtual void SetDetector(Int_t option = 0) {fOption = option;}
  virtual void SetDebug(Int_t debug = 0) {fDebug = debug;}
 protected:
  Int_t   fFirst;          // First and last events analyzed
  Int_t   fLast;           // in current set of files
  Int_t   fOption;         // detector used for jet reconstruction   
  Int_t   fDebug;          // debug option
  Float_t fFiducialEtaMin; // Fiducial minimum eta
  Float_t fFiducialEtaMax; // Fiducial maximum eta
  Float_t fPtCut;          // pt cut
  TString fComment;        // a comment
  TString fDir;            // directory with input files
  TString fPattern;        // pattern to look for input files
  
  ClassDef(AliJetReaderHeader,2);
};
 
#endif
