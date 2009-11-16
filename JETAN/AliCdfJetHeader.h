#ifndef ALICDFJETHEADER_H
#define ALICDFJETHEADER_H

/*
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice
 *
*/

// settings for jet finder process

#include "AliJetHeader.h"

class AliCdfJetHeader : public AliJetHeader
  {
  public:

  AliCdfJetHeader();
  virtual ~AliCdfJetHeader() { }

  // Getters
  Double_t GetRadius   () const { return fRadius; }
  Double_t GetJetPtCut () const { return fJetPtCut ; }
  Int_t GetMinPartJet  () const { return fMinPartJet ; }

  // Setters
  void SetRadius         ( Double_t radius )          { fRadius = radius; }
  void SetJetPtCut       ( Double_t jetptcut )        { fJetPtCut = jetptcut; }
  void SetAODwrite       ( Bool_t aodwrite )          { fAODwrite = aodwrite ; }
  void SetAODtracksWrite ( Bool_t aodtrackswrite )    { fAODtracksWrite = aodtrackswrite ; }
  void SetMinPartJet     ( Int_t npart )              { fMinPartJet = npart ; }

//  void SetCDFJetHeader   () { fCDFheader = (AliCdfJetHeader*)fHeader; }

  Bool_t IsAODwrite() const { return fAODwrite ; }
  Bool_t IsAODtracksWrite() const { return fAODtracksWrite ; }

//     void PrintParameters() const ;

  protected:

  AliCdfJetHeader(const AliCdfJetHeader &jh);
  AliCdfJetHeader& operator=(const AliCdfJetHeader &jh);

  // parameters of algorithm
  Double_t fRadius ;      //  Cone radius
  Int_t  fMinPartJet ;       // minimum number of particles in jet

  // JET Pt cut
  Double_t fJetPtCut ;  // pt cut of jets

  Bool_t fAODwrite ;         // flag for writing to AOD
  Bool_t fAODtracksWrite ;   // flag for writing tracks to AOD

//  AliCdfJetHeader* fCDFheader ; // local pointer to CDF Jet Header

  ClassDef ( AliCdfJetHeader, 1 )

  };
#endif
