#ifndef ALIJETREADERHEADER_H
#define ALIJETREADERHEADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// base class for Jet Reader Header 
//
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------
#include <Riostream.h>  
#include <TNamed.h>
#include <TString.h>
 
class AliJetReaderHeader : public TNamed
{

 public:
  AliJetReaderHeader(const char* name);
  AliJetReaderHeader();
  virtual ~AliJetReaderHeader();
  
  // Getters
  virtual const TString GetComment()  {return fComment;}
  virtual const char* GetDirectory()  {return fDir.Data();}
  virtual const char* GetBgDirectory(){return fBgDir.Data();}
  virtual const char* GetPattern() {return fPattern.Data();}
  virtual const TString GetEMCALmatrices2bLoad() {return fMatricesEMCAL;}
  virtual const TString GetEMCALgeo2bLoad() {return fGeomEMCAL;}
  virtual const TString GetMyOADBfile() {return fMyOADBfile;} 
  
  virtual Float_t     GetFiducialEtaMin() const {return fFiducialEtaMin;}
  virtual Float_t     GetFiducialEtaMax() const {return fFiducialEtaMax;} 
  virtual Float_t     GetFiducialPhiMin() const {return fFiducialPhiMin;}
  virtual Float_t     GetFiducialPhiMax() const {return fFiducialPhiMax;}  
  virtual Float_t     GetPtCut() const {return fPtCut;}
  Int_t   GetNEvents()     const {return fLast-fFirst;}
  Int_t   GetFirstEvent()  const {return fFirst;}
  Int_t   GetLastEvent()   const {return fLast;}
  Int_t   GetDetector()    const {return fOption;}
  Int_t   GetCluster()     const {return fCluster;}
  Bool_t  GetDZ()          const {return fDZ;}
  Int_t   GetDebug()       const {return fDebug;}
  Int_t   GetSignalPerBg() const {return fSignalPerBg;}
  	  
  // Setters
  virtual void SetComment(const char* s)     {fComment=TString(s);}
  virtual void SetPattern(const char* s)     {fPattern=TString(s);}
  virtual void SetDirectory(const char* s)   {fDir=TString(s);}
  virtual void SetBgDirectory(const char* s, Int_t n = 1)
      {fBgDir=TString(s); fSignalPerBg = n;}
  virtual void SetEMCALgeo2bLoad(const char* s)   {fGeomEMCAL=TString(s);} 
  virtual void SetEMCALmatrices2bLoad(const char* s)   {fMatricesEMCAL=TString(s);}
  virtual void SetMyOADBfile(const char* s)   {fMyOADBfile=TString(s);}
  virtual void SetFirstEvent(Int_t i=0) {fFirst=i;}
  virtual void SetLastEvent(Int_t i=-1) {fLast=i;}
  virtual void SetFiducialEta(Float_t etamin, Float_t etamax) 
      { fFiducialEtaMin = etamin; fFiducialEtaMax = etamax;}
  virtual void SetFiducialPhi(Float_t phimin, Float_t phimax) 
      { fFiducialPhiMin = phimin; fFiducialPhiMax = phimax;}
  virtual void SetPtCut(Float_t par = 2.0) {fPtCut = par;}
  virtual void SetDZ(Bool_t deadzone = 0) {fDZ = deadzone;}
  virtual void SetDetector(Int_t option = 0) {fOption = option;}
  virtual void SetCluster(Int_t option = 0) {fCluster = option;}
  virtual void SetDebug(Int_t debug = 0) {fDebug = debug;}

 protected:
  Int_t   fFirst;          // First and last events analyzed
  Int_t   fLast;           // in current set of files
  Int_t   fOption;         // detector used for jet reconstruction  
  Int_t   fCluster;        // cluster type
  Int_t   fDebug;          // debug option
  Bool_t  fDZ;             // include dead zones or not 
  Int_t   fSignalPerBg;
  Float_t fFiducialEtaMin; // Fiducial minimum eta
  Float_t fFiducialEtaMax; // Fiducial maximum eta
  Float_t fFiducialPhiMin; // Fiducial minimum phi
  Float_t fFiducialPhiMax; // Fiducial maximum phi
  Float_t fPtCut;          // pt cut
  TString fComment;        // a comment
  TString fDir;            // directory with input files for signal
  TString fBgDir;          // directory with input files for background
  TString fPattern;        // pattern to look for input files
  TString fMatricesEMCAL;		// survey/matrices version for EMCAL
  TString fGeomEMCAL;        // geometry version for EMCAL
  TString fMyOADBfile;      //  private version of the OADB file with EMCAL matrices 
  
  ClassDef(AliJetReaderHeader,3);
};
 
#endif
