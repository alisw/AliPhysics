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

class AliVParticle;

class AliHFEpidObject{
  public:
    enum AnalysisType_t { 
      kESDanalysis,
      kAODanalysis
    };
    AliHFEpidObject():
      fkRecTrack(NULL), 
      fAnalysisType(kESDanalysis),
      fAbInitioPID(-1),
      fCentrality(99),
      fIsPbPb(kFALSE)         // Default: pp
      {
      }
    AliHFEpidObject(const AliHFEpidObject &ref):
      fkRecTrack(ref.fkRecTrack), 
      fAnalysisType(ref.fAnalysisType),
      fAbInitioPID(ref.fAbInitioPID),
      fCentrality(ref.fCentrality),
      fIsPbPb(ref.fIsPbPb)
      {
      }
    AliHFEpidObject &operator=(const AliHFEpidObject &ref);
    ~AliHFEpidObject(){};

    void SetRecTrack(const AliVParticle * recTrack) {fkRecTrack = recTrack; }
    void SetMCTrack(const AliVParticle * mcTrack);
    void SetAnalysisType(AnalysisType_t type) { fAnalysisType = type; }
    void SetAbInitioPID(Int_t abInitioPID) { fAbInitioPID = abInitioPID; }
    void SetCentrality(Int_t centrality) { fCentrality = centrality; }
    void SetPbPb() { fIsPbPb = kTRUE; }
    void SetPP() { fIsPbPb = kFALSE; }

    const AliVParticle *GetRecTrack() const { return fkRecTrack; }
    Int_t GetAbInitioPID() const { return fAbInitioPID; }
    Int_t GetCentrality() const { return fCentrality; }
    Bool_t IsAODanalysis() const { return fAnalysisType == static_cast<UChar_t>(kAODanalysis); }
    Bool_t IsESDanalysis() const { return fAnalysisType == static_cast<UChar_t>(kESDanalysis); }
    Bool_t IsPbPb() const { return fIsPbPb; }

  private:
    const AliVParticle *fkRecTrack;     // Reconstructed track
    UChar_t fAnalysisType;              // Analysis Mode (ESD or AOD)
    Int_t fAbInitioPID;                 // AbInitio PID
    Int_t fCentrality;                  // Centrality Information
    Bool_t fIsPbPb;                     // Collision type
};
#endif

