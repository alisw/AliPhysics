#ifndef ALIJETREADERHEADER_H
#define ALIJETREADERHEADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
 
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
  AliJetReaderHeader(Int_t det);
  virtual ~AliJetReaderHeader();
  
  // Getters
  virtual const TString GetComment()                     {return fComment;}
  virtual const char*   GetDirectory()                   {return fDir.Data();}
  virtual const TString GetEMCALmatrices2bLoad()         {return fMatricesEMCAL;}
  virtual const TString GetEMCALgeo2bLoad()              {return fGeomEMCAL;}
  virtual const TString GetMyOADBfile()                  {return fMyOADBfile;}
  virtual Float_t       GetFiducialEtaMin() const        {return fFiducialEtaMin;}
  virtual Float_t       GetFiducialEtaMax() const        {return fFiducialEtaMax;} 
  virtual Float_t       GetFiducialPhiMin() const        {return fFiducialPhiMin;}
  virtual Float_t       GetFiducialPhiMax() const        {return fFiducialPhiMax;}  
  virtual Float_t       GetPtCut() const                 {return fPtCut;}
  virtual Float_t       GetEtCellCut() const             {return fEtCellCut;}
  Int_t                 GetDetector() const              {return fOption;}
  Int_t                 GetCluster() const               {return fCluster;}
  Int_t                 GetDebug() const                 {return fDebug;}
  UInt_t                GetTestFilterMask() const        {return fTestFilterMask;}
  UInt_t                GetFilterType() const            {return fFilterType;}
  TString               GetDataType() const              {return fDataType;}
  Bool_t                GetIsHighMult() const            {return fIsHighMult;}

  // Setters
  virtual void SetComment(const char* s)                 {fComment=TString(s);}
  virtual void SetDirectory(const char* s)               {fDir=TString(s);}
  virtual void SetEMCALgeo2bLoad(const char* s)          {fGeomEMCAL=TString(s);}
  virtual void SetEMCALmatrices2bLoad(const char* s)     {fMatricesEMCAL=TString(s);}
  virtual void SetMyOADBfile(const char* s)              {fMyOADBfile=TString(s);}
  virtual void SetFiducialEta(Float_t etamin, Float_t etamax) 
      { fFiducialEtaMin = etamin; fFiducialEtaMax = etamax;}
  virtual void SetFiducialPhi(Float_t phimin, Float_t phimax) 
      { fFiducialPhiMin = phimin; fFiducialPhiMax = phimax;}
  virtual void SetPtCut(Float_t par = 2.0)               {fPtCut = par;}
  virtual void SetEtCellCut(Float_t par = 0.0)           {fEtCellCut = par;}
  virtual void SetDetector(Int_t option = 0)             {fOption = option;}
  virtual void SetCluster(Int_t option = 0)              {fCluster = option;}
  virtual void SetDebug(Int_t debug = 0)                 {fDebug = debug;}

  virtual void SetDataType(const char* type /*= "AOD"*/) {fDataType = TString(type);}
  virtual void SetTestFilterMask(UInt_t i)               {fTestFilterMask = i;}
  virtual void SetFilterType(UInt_t i)                   {fFilterType = i;}
  virtual void SetReadSignalOnly(Bool_t flag = kTRUE)    {fReadSignalOnly = flag;}
  virtual void SetReadBkgdOnly(Bool_t flag = kTRUE)      {fReadBkgdOnly = flag;}
  virtual void SetIsHighMult(Bool_t mult = kFALSE)       {fIsHighMult = mult;}

  // Other
  Bool_t       ReadSignalOnly() const                    {return fReadSignalOnly;}
  Bool_t       ReadBkgdOnly() const                      {return fReadBkgdOnly;}

  // we have different cases
  // AOD reading -> MC from AOD
  // ESD reading -> MC from Kinematics
  // this has to match with our selection of input events
  enum {kTrackUndef = 0, kTrackESD, kTrackAOD, kTrackKineAll,kTrackKineCharged, kTrackAODMCAll, kTrackAODMCCharged, kTrackAODMCChargedAcceptance, kTrackAODextra, kTrackAODextraonly};

 protected:
  Int_t   fOption;                 // detector used for jet reconstruction  
  Int_t   fCluster;                // cluster type
  Int_t   fDebug;                  // debug option
  Float_t fFiducialEtaMin;         // Fiducial minimum eta
  Float_t fFiducialEtaMax;         // Fiducial maximum eta
  Float_t fFiducialPhiMin;         // Fiducial minimum phi
  Float_t fFiducialPhiMax;         // Fiducial maximum phi
  Float_t fPtCut;                  // pt cut
  Float_t fEtCellCut;              // et cell cut
  TString fComment;                // a comment
  TString fDir;                    // directory with input files for signal
  TString fMatricesEMCAL;          // survey/matrices version for EMCAL
  TString fGeomEMCAL;              // geometry version for EMCAL
  TString fMyOADBfile;             // private version of the OADB file with EMCAL matrices
  UInt_t  fTestFilterMask;         // Filter Mask for jets, not tested if 0
  UInt_t  fFilterType;             // Filter type: 0 = all, 1 = ITSTPC, 2 = TPC
  Bool_t  fReadSignalOnly;         // read particles from signal event only
  Bool_t  fReadBkgdOnly;           // read particles from bkgd event only
  TString fDataType;               // Input data type
  Bool_t  fIsHighMult;             // High multiplicity flag

  ClassDef(AliJetReaderHeader,4)   // jet reader header base class

};
 
#endif
