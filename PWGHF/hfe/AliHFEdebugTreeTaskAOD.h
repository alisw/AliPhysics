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
#ifndef ALIHFEDEBUGTREETASKAOD_H
#define ALIHFEDEBUGTREETASKAOD_H

#include "AliAnalysisTaskSE.h"

class AliHFEcuts;
class TString;
class TTreeSRedirector;
class AliHFEpidTPC;
class AliAODMCHeader;
class TClonesArray;
class AliHFEsignalCuts;
class AliHFEextraCuts;

class AliHFEdebugTreeTaskAOD : public AliAnalysisTaskSE{
  public:
    AliHFEdebugTreeTaskAOD();
    AliHFEdebugTreeTaskAOD(const char *name);
    virtual ~AliHFEdebugTreeTaskAOD();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);

    // Setters for cuts
    void SetFileName(const char *filename);
    void SetMinNclustersTPC(Int_t mincl) { fNclustersTPC = mincl; };
    void SetMinNclustersTPCPID(Int_t mincl) { fNclustersTPCPID = mincl; };
    void SetMinNclustersITS(Int_t mincl) { fNclustersITS = mincl; };
    AliHFEpidTPC *GetTPCResponse() { return fTPCpid; }
    
  private:
    AliHFEdebugTreeTaskAOD(const AliHFEdebugTreeTaskAOD &);
    AliHFEdebugTreeTaskAOD &operator=(const AliHFEdebugTreeTaskAOD &);

    AliAODMCHeader *fAODMCHeader;     // ! MC info AOD
    TClonesArray *fAODArrayMCInfo;    // ! MC info particle AOD
    AliHFEcuts *fTrackCuts;           // Track
    AliHFEextraCuts *fExtraCuts;      // HFE IP info
    AliHFEsignalCuts *fSignalCuts;    // Signal Cuts
    AliHFEpidTPC *fTPCpid;            // TPC PID
    Int_t fEventNumber;               // Event Number
    Int_t fNclustersTPC;              // Min Number of clusters in TPC
    Int_t fNclustersTPCPID;           // Min Number of clusters for TPC PID
    Int_t fNclustersITS;              // Min Number of clusters in ITS
    TString fFilename;                // file name for the debug tree
    TTreeSRedirector *fDebugTree;     // Debug Tree
  
    ClassDef(AliHFEdebugTreeTaskAOD, 1)
};
#endif

