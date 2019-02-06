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
// Debug tree to look at the distribution of the variable we are cutting on
//
//
#ifndef ALIHFEREDUCEDEVENTCREATORESD_H
#define ALIHFEREDUCEDEVENTCREATORESD_H

#include "AliAnalysisTaskSE.h"

class TString;
class TTree;
class AliAnalysisUtils;
class AliPIDResponse;
class AliHFEcuts;
class AliHFEextraCuts;
class AliHFEpidTPC;
class AliTRDTriggerAnalysis;
class AliHFEsignalCuts;
class AliHFEreducedEvent;
class AliHFEV0taginfo;

class AliHFEreducedEventCreatorESD : public AliAnalysisTaskSE{
  public:
    AliHFEreducedEventCreatorESD();
    AliHFEreducedEventCreatorESD(const char *name);
    virtual ~AliHFEreducedEventCreatorESD();

    enum EPweightType { kNone, kChi, kSigmaSquared}; // event plane weight type
    
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);

    // Setters for cuts
    void SetMinNclustersTPC(Int_t mincl) { fNclustersTPC = mincl; };
    void SetMinNclustersTPCPID(Int_t mincl) { fNclustersTPCPID = mincl; };
    void SetMinNclustersITS(Int_t mincl) { fNclustersITS = mincl; };
    void SetRemoveFirstEventFromChunk() { fRemoveFirstEvent = kTRUE; }
    void SetFlagPileupEvents() { fFlagPileupEvents = kTRUE; }
    void SetSelectSignalOnly(Bool_t select = kTRUE) { fSelectSignalOnly = select; } 
    AliHFEpidTPC *GetTPCResponse() { return fTPCpid; }

    Bool_t IsTOFmismatch(const AliVTrack *const track, const AliPIDResponse *const pid) const;
    
    void ReadVZEROCalibration2010h();
    void CalculateQvectorVZERO(Double_t Qa2[2], Double_t Qc2[2], Double_t Qa3[2], Double_t Qc3[2]) const;
    void CalculateQvectorCombinedVZERO(Double_t Q2[2], Double_t Q3[2]) const;
    Int_t GetVZEROCentralityBin() const;
    void CalculateEventPlaneVZERO(Double_t vzero[2][2]) const;
    void CalculateEventPlaneCombinedVZERO(Double_t* comb) const;
    
  private:
    AliHFEreducedEventCreatorESD(const AliHFEreducedEventCreatorESD &);
    AliHFEreducedEventCreatorESD &operator=(const AliHFEreducedEventCreatorESD &);
    
    TTree *fHFEtree;                  // HFE tree 
    AliAnalysisUtils *fAnalysisUtils; // Analysis Utils
    AliHFEreducedEvent *fHFEevent;    // hfe event
    AliHFEcuts *fTrackCuts;           // Track
    AliHFEextraCuts *fExtraCuts;      // HFE IP info
    AliHFEsignalCuts *fSignalCuts;    // Signal Cuts
    AliHFEpidTPC *fTPCpid;            // TPC PID
    AliHFEV0taginfo *fV0Tagger;       // Tags v0 tracks per Event
    AliTRDTriggerAnalysis *fTRDTriggerAnalysis; //! TRD Trigger Analysis
    Int_t fEventNumber;               // Event Number
    Int_t fNclustersTPC;              // Min Number of clusters in TPC
    Int_t fNclustersTPCPID;           // Min Number of clusters for TPC PID
    Int_t fNclustersITS;              // Min Number of clusters in ITS
    Bool_t fRemoveFirstEvent;         // Remove first event from chunk
    Bool_t fFlagPileupEvents;         // Flag pileup events
    Bool_t fSelectSignalOnly;         // Select signal-only tracks
    
    // EP correction variables
    //vzero event plane calibration cache for 10h data ,  from Redmers task
    Int_t fRunNumber; //! current runnumber (for QA and jet, track selection)
    Int_t fRunNumberCaliInfo; //! runnumber of the cached calibration info
    Float_t                 fMeanQ[9][2][2];                //! recentering
    Float_t                 fWidthQ[9][2][2];               //! recentering
    Float_t                 fMeanQv3[9][2][2];              //! recentering
    Float_t                 fWidthQv3[9][2][2];             //! recentering
    TH1*                    fVZEROgainEqualization;         //! equalization histo
    Float_t                 fVZEROApol;                     //! calibration info per disc
    Float_t                 fVZEROCpol;                     //! calibration info per disc
    TArrayD*                fChi2A;                         // chi vs cent for vzero A ep_2
    TArrayD*                fChi2C;                         // chi vs cent for vzero C ep_2
    TArrayD*                fChi3A;                         // chi vs cent for vzero A ep_3
    TArrayD*                fChi3C;                         // chi vs cent for vzero C ep_3
    TArrayD*                fSigma2A;                       // chi vs cent for vzero A ep_2
    TArrayD*                fSigma2C;                       // chi vs cent for vzero C ep_2
    TArrayD*                fSigma3A;                       // chi vs cent for vzero A ep_3
    TArrayD*                fSigma3C;                       // chi vs cent for vzero C ep_3
    EPweightType            fWeightForVZERO;                // use chi weight for vzero
    TFile* fOADB; //! fOADB
    TArrayD* fCentralityClasses; //-> centrality classes (maximum 10)
    Int_t fInCentralitySelection; //! centrality bin
        
    ClassDef(AliHFEreducedEventCreatorESD, 2)
};
#endif


