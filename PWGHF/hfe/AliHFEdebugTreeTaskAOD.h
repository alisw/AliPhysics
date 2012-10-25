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
class TTree;

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
    void SetDebugStream(Bool_t on) { fDebugstream = on; };
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
    Bool_t  fDebugstream;             // to choose the way
    TTreeSRedirector *fDebugTree;     // Debug Tree
    TTree *fDebugTreee;               // Debug Tree

    Float_t fCentrality;              // variable
    Int_t   fRun;                     // run
    Int_t   fDoublec;                 // double counted
    Float_t fMomentum;                // Momentum
    Float_t fMomentumTPC;             // Momentum TPC
    Float_t fTransverseMomentum;      // Transverse Momentum
    Float_t fEta;                     // Eta
    Float_t fPhi;                     // Phi
    Int_t   fCharge;                  // charge
    Int_t   fNClustersTPCall;         // Nb of TPC clusters TPC all
    Int_t   fNClustersTPCPID;         // Nb of TPC clusters TPC PID
    Int_t   fNClustersTPCshared;      // Nb of TPC clusters shared
    Int_t   fNCrossedRowsTPC;         // Nb of cross row TPC
    Float_t fClusterRatioTPCall;      // cls ratio TPC all
    Int_t   fNClustersITS;            // Nb of ITS clusters
    Int_t   fStatusL0;                // status L0
    Int_t   fStatusL1;                // status L1
    Float_t fSigmaTOF;                // Sigma TOF
    Float_t fSigmaTPC;                // Sigma TPC
    Float_t fDcaxy;                   // Dcaxy
    Float_t fDcaz;                    // Dcaz
    Int_t   fFilter2;                 // filter2
    Int_t   fFilter4;                 // filter4
    Int_t   fSource;                  // source
    Float_t fEr;                      // er
    Float_t fSignal;                  // signal

    
  
    ClassDef(AliHFEdebugTreeTaskAOD, 2)
};
#endif

