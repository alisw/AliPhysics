#ifndef ALIHFCUTOPTTREEHANDLER_H
#define ALIHFCUTOPTTREEHANDLER_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// Class AliHFCutOptTreeHandler
// helper class to handle a tree for cut optimisation and MVA analyses
// Authors:
// F. Catalano, fabio.catalano@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
// F. Prino, prino@to.infn.it
// L. Vermunt, luuk.vermunt@cern.ch
/////////////////////////////////////////////////////////////

#include <TTree.h>
#include <TClonesArray.h>
#include "AliAODRecoDecayHF.h"
#include "AliAODPidHF.h"

class AliHFCutOptTreeHandler : public TObject
{
  public:
    enum decaychannel {
      kD0toKpi,
      kDplustoKpipi,
      kDstoKKpi
    };

    enum candtype {
      kBkg,
      kPromptSig,
      kFDSig,
      kPromptRefl,
      kFDRefl
    };

    enum optpid {
      kNoPID,
      kNsigmaPID,
      kNsigmaPIDchar,
      kNsigmaPIDfloatandchar, //--> to test
      kNsigmaCombPID,
      kNsigmaCombPIDchar,
      kNsigmaCombPIDfloatandchar //--> to test
    };

    AliHFCutOptTreeHandler();
    AliHFCutOptTreeHandler(int decay, int PIDopt, bool isMC);
    virtual ~AliHFCutOptTreeHandler();

    bool SetVariables(AliAODRecoDecayHF* d, int masshypo, AliAODPidHF* pidHF=0x0, TClonesArray* arrayMC=0x0);
    TTree* BuildTree(TString name="tree", TString title="tree");
    bool FillTree() {
      if(!fTreeTopolVar) return false;
      if(!(fCandType==kBkg && fFillOnlySignal)) fTreeTopolVar->Fill();
      return true;
    }

    void SetDecayChannel(int decay=kD0toKpi) {fDecayChannel=decay; SetPdgCodes();}
    void SetOptPID(int PIDopt) {fPidOpt=PIDopt;}
    void SetCentrality(char centrality) {fCentrality=centrality;}
    void SetUseCentrality(bool usecent=true) {fUseCentrality=usecent;}
    void SetIsSelectedStd(bool isselected=true) {fIsSelStd=isselected;}
    void SetUseSelectedStdFlag(bool useselflag=true) {fUseSelFlag=useselflag;}
    void SetFillOnlySignal(bool fillopt=true) {fFillOnlySignal=fillopt;}
    void SetIsMC(bool isMC=true) {fIsMC=isMC;}
    void SetCandidateType(bool issignal, bool isprompt, bool isreflected) {
      if(issignal) fIsSignal=1;
      else fIsSignal=0;
      if(isprompt) fIsPrompt=1;
      else fIsPrompt=0;
      if(fDecayChannel!=kDplustoKpipi) {
        if(isreflected) fIsRefl=1;
        else fIsRefl=0;
      }
      else fIsRefl=0; //D+ -> Kpipi signal never reflected
    }

  private:
    void SetPidVars(AliAODRecoDecayHF* d, AliAODPidHF* pidHF);
    void SetPdgCodes();

    enum {knMaxProngs=3,knPidVars=12,knTopolVars=15,knTopolVarsCommon=10,knTopolVarsDzero=4,knTopolVarsDs=5,knTopolVarsDplus=2};

    TTree* fTreeTopolVar; /// tree for cut optimisation
    int fDecayChannel; /// decay channel
    int fPdgCode; /// absolute value of pdg code of the particle of interest
    int fPdgCodeProngs[knMaxProngs]; ///absolute values of pdg codes of the daughters
    float fTopolVarVector[knTopolVars]; /// array with topological variables 
    float fPIDnSigmaVector[knPidVars]; /// array with nSigma PID
    float fPIDnSigmaCharVector[knPidVars]; /// array with nSigma PID (char)
    int fPidOpt; /// option for PID variables
    char fCandType; ///flag for candidate type (bkg, prompt signal, FD signal, prompt refl, FD refl)
    bool fUseCentrality; ///flag to enable centrality
    char fCentrality; ///centrality in case of p-Pb or Pb-Pb
    bool fFillOnlySignal; ///flag to enable only signal filling
    bool fIsMC; ///flag to enable checks on MC truth 
    int fIsSignal; ///flag for signal=1 (including prompt, FD, reflected), bkg=0 
    int fIsPrompt; ///flag for prompt=1 (inluding reflected), FD=0
    int fIsRefl; ///flag for reflected signal=1, non-reflected signal=0 
    char fIsSelStd; ///flag to tag selected candidates by "standard" cuts
    bool fUseSelFlag; ///flag to enable branch on candidates selected by "standard" cuts

  /// \cond CLASSIMP
  ClassDef(AliHFCutOptTreeHandler,1); /// 
  /// \endcond
};

#endif
