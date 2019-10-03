/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

// $Id: AliJEventHeader.h,v 1.1 2008/05/02 11:56:23 djkim Exp $
////////////////////////////////////////////////////
/*!
  \file AliJEventHeader.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.1 $
  \date $Date: 2008/05/02 11:56:39 $

*/
////////////////////////////////////////////////////

#ifndef ALIJEVENTHEADER_H
#define ALIJEVENTHEADER_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "AliJBaseEventHeader.h"

class AliJEventHeader : public AliJBaseEventHeader {
 public:

  enum { kcV0M, kcFMD, kcTRK, kcTKL, kcCL0, kcCL1, kcV0MvsFMD, kcTKLvsV0, kcZEMvsZDC, kcV0A, kcV0C, kcNTYPE };

  /*
   * V0M = V0 multiplicity
   * FMD = FMD raw multiplicity
   * TRK = N. of tracks
   * TKL = N. of tracklets
   * CL0 = N. of clusters in layer 0
   * CL1 = N. of clusters in layer 1
   * V0MvsFMD = correlation between V0 and FMD
   * TKLvsV0 = correlation between tracklets and V0
   * ZEMvsZDC = correlation between ZEM and ZDC 
   */
  AliJEventHeader();              // default constructor
  AliJEventHeader(int eventid,
                  float cent,
                  float vrtz,
                  ULong64_t triggmaskAli,
                  UInt_t triggmaskJC,
                  Int_t  refmult,
                  Int_t  refmult1,
                  Int_t  refmult2,
                  Float_t  v0mult,
                  Float_t  v0Amult,
                  Float_t  v0Cmult,
                  UInt_t eventType
                 );

  AliJEventHeader(const AliJEventHeader& a);                           

  virtual ~AliJEventHeader(){;}     // destructor

  ULong64_t  GetTriggerMaskAlice()   const {return fTriggerMaskAlice;}  
  UInt_t     GetTriggerMaskJCorran() const {return fTriggerMaskJCorran;}  
  Int_t      GetSPDTrackletMult()    const {return fSPDTrackletMult;}
  Int_t      GetITSSATrackletMult()    const {return fTrackletsITSSA;}
  Int_t      GetITSTPCTrackletMult()    const {return fTrackletsITSTPC;}
  UInt_t     GetEventType()          const {return fEventType;}
  Float_t    GetV0Mult()             const {return fV0Mult;}
  Float_t    GetV0AMult()             const {return fV0AMult;}
  Float_t    GetV0CMult()             const {return fV0CMult;}
  Int_t      GetVtxMult()            const { return fVtxMult; }
  UShort_t   GetBunchCrossNumber()   const { return fBunchCrossNumber; }

  Float_t GetCentralityArray( UInt_t it ) const { return it<kcNTYPE ? fCentralityArray[it] : -1; }

  void SetTriggerMaskAlice(ULong64_t mask) {fTriggerMaskAlice = mask;}
  void SetTriggerMaskJCorran(UInt_t mask) {fTriggerMaskJCorran = mask;}
  void SetSPDTrackletMult(Int_t ref) { fSPDTrackletMult = ref;}
  void SetITSSATrackletMult(Int_t ref)  { fTrackletsITSSA = ref;}
  void SetITSTPCTrackletMult(Int_t ref)  { fTrackletsITSTPC = ref;}
  void SetEventType(UInt_t eventype) {fEventType = eventype;}
  void SetV0Mult(Float_t multV0) {fV0Mult = multV0;}
  void SetV0AMult(Float_t multV0) {fV0AMult = multV0;}
  void SetV0CMult(Float_t multV0) {fV0CMult = multV0;}
  void SetVtxMult(Int_t m){ fVtxMult = m; };
  void SetCentralityArray(UInt_t it, Float_t cen ){ if( it < kcNTYPE ) fCentralityArray[it]=cen; }
  void SetBunchCrossNumber( UShort_t n ){ fBunchCrossNumber = n; }

  TString GetFiredTriggers() const { return fFiredTriggers; }
  void SetFiredTriggers(TString s){ fFiredTriggers=s; }

  AliJEventHeader&  operator=(const AliJEventHeader& header);

  TString GetESDFileName() const { return fESDFileName; }
  void SetESDFileName(TString s){ fESDFileName=s; }
  Int_t GetEventNumberESDFile() const { return fEventNumberESDFile; }
  void SetEventNumberESDFile(Int_t s){ fEventNumberESDFile=s; }

  void SetL0TriggerInputs(UInt_t n) {fL0TriggerInputs=n;}
  UInt_t      GetL0TriggerInputs() const {return fL0TriggerInputs;}  


 private:

  ULong64_t   fTriggerMaskAlice;           //Alice Trigger MASK
  UInt_t      fTriggerMaskJCorran;         // JCorran Trigger MASK
  Int_t       fSPDTrackletMult;             //SPD tracklet multiplicity
  Int_t       fTrackletsITSTPC;           // Multiplicity ITS and TPC
  Int_t       fTrackletsITSSA;            // Multiplicity ITS standalone + ITS
  Double32_t   fV0Mult;                   // VZERO multiplicity
  Double32_t   fV0AMult;                   // VZERO multiplicity
  Double32_t   fV0CMult;                   // VZERO multiplicity
  UInt_t      fEventType;                 // Type of Event
  TString     fFiredTriggers;       // String with fired triggers from AOD
  Int_t       fVtxMult;                   //FK// EFF number of vertex contributors 
  Double32_t  fCentralityArray[kcNTYPE];  //?//
  UShort_t    fBunchCrossNumber;   // bunch crossing identifier
  TString   fESDFileName;          // file name for the ESD file
  Int_t     fEventNumberESDFile;   // Number of event in the ESD file
  UInt_t    fL0TriggerInputs;    //L0 Trigger Inputs (mask)

  ClassDef(AliJEventHeader,3)

};

#endif
