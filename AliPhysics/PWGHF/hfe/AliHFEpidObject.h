/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Object used in the electron identification
// For more information see the implementation file
//
#ifndef ALIHFEPIDOBJECT_H
#define ALIHFEPIDOBJECT_H

#include <Rtypes.h>

class AliVTrack;
class AliVParticle;

class AliHFEpidObject{
  public:
    enum AnalysisType_t{ 
      kESDanalysis,
      kAODanalysis
    };
    AliHFEpidObject():
      fkRecTrack(NULL), 
      fAnalysisType(kESDanalysis),
      fAbInitioPID(-1),
      fCentrality(99),
      fMultiplicity(0),
      fCorrTPCnSigma(0),
      fIsPbPb(kFALSE),         // Default: pp
      fIspPb(kFALSE),         // Default: pp
      fHasCorrTPCnSigma(kFALSE)
      {
      }
    AliHFEpidObject(const AliHFEpidObject &ref):
      fkRecTrack(ref.fkRecTrack), 
      fAnalysisType(ref.fAnalysisType),
      fAbInitioPID(ref.fAbInitioPID),
      fCentrality(ref.fCentrality),
      fMultiplicity(ref.fMultiplicity),
      fCorrTPCnSigma(ref.fCorrTPCnSigma),
      fIsPbPb(ref.fIsPbPb),
      fIspPb(ref.fIspPb),
      fHasCorrTPCnSigma(ref.fHasCorrTPCnSigma)
      {
      }
    AliHFEpidObject &operator=(const AliHFEpidObject &ref);
    ~AliHFEpidObject(){};

    void SetRecTrack(const AliVTrack * recTrack) {fkRecTrack = recTrack; }
    void SetMCTrack(const AliVParticle * mcTrack);
    void SetAnalysisType(AnalysisType_t type) { fAnalysisType = type; }
    void SetAbInitioPID(Int_t abInitioPID) { fAbInitioPID = abInitioPID; }
    void SetCentrality(Int_t centrality) { fCentrality = centrality; }
    void SetMulitplicity(Double_t mult) { fMultiplicity = mult; }
    void SetCorrectedTPCnSigma(Double_t sigm) { fCorrTPCnSigma=sigm; fHasCorrTPCnSigma=true; }
    void SetPbPb() { fIsPbPb = kTRUE; }
    void SetpPb() { fIsPbPb = kFALSE; fIspPb=kTRUE; }
    void SetPP() { fIsPbPb = kFALSE; }

    const AliVTrack *GetRecTrack() const { return fkRecTrack; }
    Int_t GetAbInitioPID() const { return fAbInitioPID; }
    Int_t GetCentrality() const { return fCentrality; }
    Double_t GetMultiplicity() const { return fMultiplicity; }
    Double_t GetCorrectedTPCnSigma()  const { return fCorrTPCnSigma; }
    Bool_t HasCorrectedTPCnSigma()  const { return fHasCorrTPCnSigma; }
    Bool_t IsAODanalysis() const { return fAnalysisType == static_cast<UChar_t>(kAODanalysis); }
    Bool_t IsESDanalysis() const { return fAnalysisType == static_cast<UChar_t>(kESDanalysis); }
    Bool_t IsPbPb() const { return fIsPbPb; }
    Bool_t IspPb() const { return fIspPb; }

  private:
    const AliVTrack *fkRecTrack;        // Reconstructed track
    UChar_t fAnalysisType;              // Analysis Mode (ESD or AOD)
    Int_t fAbInitioPID;                 // AbInitio PID
    Int_t fCentrality;                  // Centrality Information
    Double_t fMultiplicity;             // Multiplicity information
    Double_t fCorrTPCnSigma;            // Corrected TPC n sigma
    Bool_t fIsPbPb;                     // Collision type PbPb
    Bool_t fIspPb;                      // Collision type pPb
    Bool_t fHasCorrTPCnSigma;           // whether corrected TPC n sigma is set
};
#endif

