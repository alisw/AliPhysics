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
#ifndef ALIHFEDEBUGTREETASK_H
#define ALIHFEDEBUGTREETASK_H

#include "AliAnalysisTaskSE.h"

class AliHFEcuts;
class AliHFEsignalCuts;
class AliTrackReference;
class TString;
class TTreeSRedirector;
class AliHFEpidTRD;
class AliHFEpidTPC;

class AliHFEdebugTreeTask : public AliAnalysisTaskSE{
  public:
    AliHFEdebugTreeTask();
    AliHFEdebugTreeTask(const char *name);
    virtual ~AliHFEdebugTreeTask();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *) {};

    // Setters for cuts
    void SetFileName(const char *filename);
    void SetMinNclustersTPC(Int_t mincl) { fNclustersTPC = mincl; };
    void SetMinNclustersTPCPID(Int_t mincl) { fNclustersTPCPID = mincl; };
    void SetMinNclustersITS(Int_t mincl) { fNclustersITS = mincl; };
    Int_t GetElecSourceMC(TParticle * const mcpart); // return electron source id
    AliHFEpidTPC *GetTPCResponse() { return fTPCpid; }
    
  private:
    AliHFEdebugTreeTask(const AliHFEdebugTreeTask &);
    AliHFEdebugTreeTask &operator=(const AliHFEdebugTreeTask &);

    AliTrackReference *FindTrackReference(AliMCParticle *track, Float_t minRadius, Float_t maxRadius, Int_t detectorID);
    
    AliHFEcuts *fTrackCuts;           // Track
    AliHFEsignalCuts *fSignalCuts;    // Signal Cuts
    AliHFEpidTRD *fTRDpid;            // TRD PID
    AliHFEpidTPC *fTPCpid;            // TPC PID
    AliHFEextraCuts *fExtraCuts;      // HFE IP info
    Int_t fNclustersTPC;              // Min Number of clusters in TPC
    Int_t fNclustersTPCPID;           // Min Number of clusters for TPC PID
    Int_t fNclustersITS;              // Min Number of clusters in ITS
    TString fFilename;                // file name for the debug tree
    TTreeSRedirector *fDebugTree;     // Debug Tree
    Int_t fNparents; // number of heavy hadrons to be considered
    Int_t fParentSelect[2][7]; // heavy hadron species
    static const Int_t fgkGluon=21; // gluon pdg code
    static const Int_t fgkMaxGener=10; // ancester level wanted to be checked 
    static const Int_t fgkMaxIter=100; // number of iteration to find out matching particle 
    static const Int_t fgkqType=7; // number of particle type to be checked
    static const Int_t fgkEtaRanges=3; // cuts for different eta ranges
    ClassDef(AliHFEdebugTreeTask, 1)
};
#endif

