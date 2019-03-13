/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandler
// \brief helper class to handle a tree for cut optimisation and MVA analyses
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
// F. Prino, prino@to.infn.it
// L. Vermunt, luuk.vermunt@cern.ch
// L. van Doremalen, lennart.van.doremalen@cern.ch
// J. Norman, jaime.norman@cern.ch
// G. Luparello, grazia.luparello@cern.ch
/////////////////////////////////////////////////////////////

#include <cmath>
#include <limits>
#include "AliHFTreeHandler.h"
#include "AliPID.h"
#include "AliAODRecoDecayHF.h"
#include "AliPIDResponse.h"
#include "AliESDtrack.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandler);
/// \endcond

//________________________________________________________________
AliHFTreeHandler::AliHFTreeHandler():
  TObject(),
  fTreeVar(0x0),
  fNProngs(),
  fNCandidates(0),
  fCandTypeMap(0),
  fCandType(),
  fInvMass(),
  fPt(),
  fY(),
  fEta(),
  fPhi(),
  fDecayLength(),
  fDecayLengthXY(),
  fNormDecayLengthXY(),
  fCosP(),
  fCosPXY(),
  fImpParXY(),
  fDCA(),
  fPProng(),
  fTPCPProng(),
  fTOFPProng(),
  fPtProng(),
  fEtaProng(),
  fPhiProng(),
  fNTPCclsProng(),
  fNTPCclsPidProng(),
  fNTPCCrossedRowProng(),
  fChi2perNDFProng(),
  fNITSclsProng(),
  fITSclsMapProng(),
  fTrackIntegratedLengthProng(),
  fStartTimeResProng(),
  fPIDNsigmaVector(),
  fPIDNsigmaIntVector(),
  fPIDrawVector(),
  fPidOpt(kNsigmaPID),
  fSingleTrackOpt(kRedSingleTrackVars),
  fFillOnlySignal(false),
  fIsMCGenTree(false),
  fDauInAccFlag(false),
  fDauInAcceptance(),
  fEvID(),
  fRunNumber()
{
  //
  // Default constructor
  //
}

//________________________________________________________________
AliHFTreeHandler::AliHFTreeHandler(int PIDopt):
  TObject(),
  fTreeVar(0x0),
  fNProngs(),
  fNCandidates(0),
  fCandTypeMap(0),
  fCandType(),
  fInvMass(),
  fPt(),
  fY(),
  fEta(),
  fPhi(),
  fDecayLength(),
  fDecayLengthXY(),
  fNormDecayLengthXY(),
  fCosP(),
  fCosPXY(),
  fImpParXY(),
  fDCA(),
  fPProng(),
  fTPCPProng(),
  fTOFPProng(),
  fPtProng(),
  fEtaProng(),
  fPhiProng(),
  fNTPCclsProng(),
  fNTPCclsPidProng(),
  fNTPCCrossedRowProng(),
  fChi2perNDFProng(),
  fNITSclsProng(),
  fITSclsMapProng(),
  fTrackIntegratedLengthProng(),
  fStartTimeResProng(),
  fPIDNsigmaVector(),
  fPIDNsigmaIntVector(),
  fPIDrawVector(),
  fPidOpt(PIDopt),
  fSingleTrackOpt(kRedSingleTrackVars),
  fFillOnlySignal(false),
  fIsMCGenTree(false),
  fDauInAccFlag(false),
  fDauInAcceptance(),
  fEvID(),
  fRunNumber()
{
  //
  // Standard constructor
  //
}

//________________________________________________________________
AliHFTreeHandler::~AliHFTreeHandler()
{
  //
  // Destructor
  //

  if(fTreeVar) delete fTreeVar;
}

//________________________________________________________________
TTree* AliHFTreeHandler::BuildTreeMCGen(TString name, TString title) {

  fIsMCGenTree = true;

  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=0x0;
  }
  fTreeVar = new TTree(name.Data(),title.Data());
  fTreeVar->Branch("run_number",&fRunNumber);
  fTreeVar->Branch("ev_id",&fEvID);
  fTreeVar->Branch("n_cand",&fNCandidates);
  fTreeVar->Branch("cand_type",&fCandType);
  fTreeVar->Branch("pt_cand",&fPt);
  fTreeVar->Branch("y_cand",&fY);
  fTreeVar->Branch("eta_cand",&fEta);
  fTreeVar->Branch("phi_cand",&fPhi);
  fTreeVar->Branch("dau_in_acc",&fDauInAcceptance);

  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandler::SetMCGenVariables(int runnumber, unsigned int eventID, AliAODMCParticle* mcpart) {

  if(!mcpart) return false;
  if(!(fCandTypeMap&kSignal)) return true; // fill only signal in the generated

  fNCandidates++;

  fRunNumber.push_back(runnumber);
  fEvID.push_back(eventID);
  fCandType.push_back(fCandTypeMap);
  fPt.push_back(mcpart->Pt());
  fY.push_back(mcpart->Y());
  fEta.push_back(mcpart->Eta());
  fPhi.push_back(mcpart->Phi());
  fDauInAcceptance.push_back(fDauInAccFlag);
  
  return true;
}

//________________________________________________________________
void AliHFTreeHandler::SetCandidateType(bool issignal, bool isbkg, bool isprompt, bool isFD, bool isreflected) 
{  
  if(issignal) fCandTypeMap |= kSignal;
  else fCandTypeMap &= ~kSignal;
  if(isbkg && !fIsMCGenTree) fCandTypeMap |= kBkg;
  else fCandTypeMap &= ~kBkg;
  if(isprompt) fCandTypeMap |= kPrompt;
  else fCandTypeMap &= ~kPrompt;
  if(isFD) fCandTypeMap |= kFD;
  else fCandTypeMap &= ~kFD;
  if(isreflected && !fIsMCGenTree) fCandTypeMap |= kRefl;
  else fCandTypeMap &= ~kRefl;
}

//________________________________________________________________
void AliHFTreeHandler::AddCommonDmesonVarBranches() {

  fTreeVar->Branch("run_number",&fRunNumber);
  fTreeVar->Branch("ev_id",&fEvID);
  fTreeVar->Branch("n_cand",&fNCandidates);
  fTreeVar->Branch("cand_type",&fCandType);
  fTreeVar->Branch("inv_mass",&fInvMass);
  fTreeVar->Branch("pt_cand",&fPt);
  fTreeVar->Branch("y_cand",&fY);
  fTreeVar->Branch("eta_cand",&fEta);
  fTreeVar->Branch("phi_cand",&fPhi);
  fTreeVar->Branch("d_len",&fDecayLength);
  fTreeVar->Branch("d_len_xy",&fDecayLengthXY);
  fTreeVar->Branch("norm_dl_xy",&fNormDecayLengthXY);
  fTreeVar->Branch("cos_p",&fCosP);
  fTreeVar->Branch("cos_p_xy",&fCosPXY);
  fTreeVar->Branch("imp_par_xy",&fImpParXY);
  fTreeVar->Branch("dca",&fDCA);
} 

//________________________________________________________________
void AliHFTreeHandler::AddSingleTrackBranches() {

  if(fSingleTrackOpt==kNoSingleTrackVars) return;

  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {

    if(fSingleTrackOpt==kRedSingleTrackVars) {
      fTreeVar->Branch(Form("pt_prong%d",iProng),&fPtProng[iProng]);
      fTreeVar->Branch(Form("eta_prong%d",iProng),&fEtaProng[iProng]);
      fTreeVar->Branch(Form("phi_prong%d",iProng),&fPhiProng[iProng]);
    }
    else if(fSingleTrackOpt==kAllSingleTrackVars) {
      fTreeVar->Branch(Form("pt_prong%d",iProng),&fPtProng[iProng]);
      fTreeVar->Branch(Form("eta_prong%d",iProng),&fEtaProng[iProng]);
      fTreeVar->Branch(Form("phi_prong%d",iProng),&fPhiProng[iProng]);
      fTreeVar->Branch(Form("p_prong%d",iProng),&fPProng[iProng]);
      fTreeVar->Branch(Form("nTPCcls_prong%d",iProng),&fNTPCclsProng[iProng]);
      fTreeVar->Branch(Form("nTPCclspid_prong%d",iProng),&fNTPCclsPidProng[iProng]);
      fTreeVar->Branch(Form("nTPCcrossrow_prong%d",iProng),&fNTPCCrossedRowProng[iProng]);
      fTreeVar->Branch(Form("chi2perndf_prong%d",iProng),&fChi2perNDFProng[iProng]);
      fTreeVar->Branch(Form("nITScls_prong%d",iProng),&fNITSclsProng[iProng]);
      fTreeVar->Branch(Form("ITSclsmap_prong%d",iProng),&fITSclsMapProng[iProng]);
    }
  }
}

//________________________________________________________________
void AliHFTreeHandler::AddPidBranches(bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF) 
{

  if(fPidOpt==kNoPID) return;
  if(fPidOpt>kRawAndNsigmaPID) {
    AliWarning("Wrong PID setting!");
    return;
  }

  bool useHypo[knMaxHypo4Pid] = {usePionHypo,useKaonHypo,useProtonHypo};
  bool useDet[knMaxDet4Pid] = {useTPC,useTOF};
  TString partHypoName[knMaxHypo4Pid] = {"Pi","K","Pr"};
  TString detName[knMaxDet4Pid] = {"TPC","TOF"};
  TString rawPidName[knMaxDet4Pid] = {"dEdxTPC","ToF"};

  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    if((fPidOpt>=kNsigmaPID && fPidOpt<=kNsigmaPIDfloatandint) || fPidOpt==kRawAndNsigmaPID) {
      for(unsigned int iDet=0; iDet<knMaxDet4Pid; iDet++) {
        if(!useDet[iDet]) continue;
        for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
          if(!useHypo[iPartHypo]) continue;
          if(fPidOpt==kNsigmaPID || fPidOpt==kNsigmaPIDfloatandint || fPidOpt==kRawAndNsigmaPID) fTreeVar->Branch(Form("nsig%s_%s_%d",detName[iDet].Data(),partHypoName[iPartHypo].Data(),iProng),&fPIDNsigmaVector[iProng][iDet][iPartHypo]);
          if(fPidOpt==kNsigmaPIDint || fPidOpt==kNsigmaPIDfloatandint) fTreeVar->Branch(Form("nsig%s_%s_%d",detName[iDet].Data(),partHypoName[iPartHypo].Data(),iProng),&fPIDNsigmaIntVector[iProng][iDet][iPartHypo]);
        }
      }
    }
    if(fPidOpt>=kNsigmaCombPID && fPidOpt<=kNsigmaCombPIDfloatandint) {
      for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
        if(!useHypo[iPartHypo]) continue;
        if(fPidOpt==kNsigmaCombPID || fPidOpt==kNsigmaCombPIDfloatandint) fTreeVar->Branch(Form("nsigComb_%s_%d",partHypoName[iPartHypo].Data(),iProng),&fPIDNsigmaVector[iProng][0][iPartHypo]);
        if(fPidOpt==kNsigmaCombPIDint || fPidOpt==kNsigmaCombPIDfloatandint) fTreeVar->Branch(Form("int_nsigComb_%s_%d",partHypoName[iPartHypo].Data(),iProng),&fPIDNsigmaIntVector[iProng][0][iPartHypo]);
      }
    }
    if(fPidOpt==kRawPID || fPidOpt==kRawAndNsigmaPID) {
      for(unsigned int iDet=0; iDet<knMaxDet4Pid; iDet++) {
        if(!useDet[iDet]) continue;
        fTreeVar->Branch(Form("%s_%d",rawPidName[iDet].Data(),iProng),&fPIDrawVector[iProng][iDet]);
      }
      if(useTPC) fTreeVar->Branch(Form("pTPC_prong%d",iProng),&fTPCPProng[iProng]);
      if(useTOF) {
        fTreeVar->Branch(Form("pTOF_prong%d",iProng),&fTOFPProng[iProng]);
        fTreeVar->Branch(Form("trlen_prong%d",iProng),&fTrackIntegratedLengthProng[iProng]);
        fTreeVar->Branch(Form("start_time_res_prong%d",iProng),&fStartTimeResProng[iProng]);
      }
    }
  }
}

//________________________________________________________________
bool AliHFTreeHandler::SetSingleTrackVars(AliAODTrack* prongtracks[]) {

  //Impact parameters of the prongs are defined as a species dependent variable because the prongs 
  //cannot be obtained in similar way for the different AliAODRecoDecay objects (AliAODTrack cannot
  //be used because of recomputation PV)

  if(fSingleTrackOpt==kNoSingleTrackVars) return true;

  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    if(!prongtracks[iProng]) {
      AliWarning("Prong track not found!");
      return false;
    }
  }

  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {

    if(fSingleTrackOpt==kRedSingleTrackVars) {
      fPtProng[iProng].push_back(prongtracks[iProng]->Pt());
      fEtaProng[iProng].push_back(prongtracks[iProng]->Eta());
      fPhiProng[iProng].push_back(prongtracks[iProng]->Phi());
    }
    else if(fSingleTrackOpt==kAllSingleTrackVars) {
      fPtProng[iProng].push_back(prongtracks[iProng]->Pt());
      fEtaProng[iProng].push_back(prongtracks[iProng]->Eta());
      fPhiProng[iProng].push_back(prongtracks[iProng]->Phi());
      fPProng[iProng].push_back(prongtracks[iProng]->P());
      fNTPCclsProng[iProng].push_back(prongtracks[iProng]->GetTPCNcls());
      fNTPCclsPidProng[iProng].push_back(prongtracks[iProng]->GetTPCsignalN());
      fNTPCCrossedRowProng[iProng].push_back(prongtracks[iProng]->GetTPCNCrossedRows());
      fChi2perNDFProng[iProng].push_back(prongtracks[iProng]->Chi2perNDF());
      fNITSclsProng[iProng].push_back(prongtracks[iProng]->GetITSNcls());
      fITSclsMapProng[iProng].push_back(static_cast<int>(prongtracks[iProng]->GetITSClusterMap()));
    }
  }

  return true;
}

//________________________________________________________________
bool AliHFTreeHandler::SetPidVars(AliAODTrack* prongtracks[], AliPIDResponse* pidrespo, bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF) 
{
  if(!pidrespo) return false;
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    if(!prongtracks[iProng]) {
      AliWarning("Prong track not found!");
      return false;
    }
  }

  //PID variables
  double sig[knMaxProngs][knMaxDet4Pid][knMaxHypo4Pid];
  double sigComb[knMaxProngs][knMaxHypo4Pid];
  double rawPID[knMaxProngs][knMaxDet4Pid];
  bool useHypo[knMaxHypo4Pid] = {usePionHypo,useKaonHypo,useProtonHypo};
  bool useDet[knMaxDet4Pid] = {useTPC,useTOF};
  AliPID::EParticleType parthypo[knMaxHypo4Pid] = {AliPID::kPion,AliPID::kKaon,AliPID::kProton};
  
  //compute PID variables for different options
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    if((fPidOpt>=kNsigmaPID && fPidOpt<=kNsigmaCombPIDfloatandint) || fPidOpt==kRawAndNsigmaPID) {
      for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
        if(useHypo[iPartHypo]) {
          if(useTPC) sig[iProng][kTPC][iPartHypo] = pidrespo->NumberOfSigmasTPC(prongtracks[iProng],parthypo[iPartHypo]);
          if(useTOF) sig[iProng][kTOF][iPartHypo] = pidrespo->NumberOfSigmasTOF(prongtracks[iProng],parthypo[iPartHypo]);
          if(fPidOpt>=kNsigmaCombPID && fPidOpt<=kNsigmaCombPIDfloatandint) {
            sigComb[iProng][iPartHypo] = CombineNsigmaDiffDet(sig[iProng][kTPC][iPartHypo],sig[iProng][kTOF][iPartHypo]);
          }
        }
      }
    }
    if(fPidOpt==kRawPID || fPidOpt==kRawAndNsigmaPID) {
      if(useTPC) rawPID[iProng][kTPC] = prongtracks[iProng]->GetTPCsignal();
      if(useTOF) {
        if (!(prongtracks[iProng]->GetStatus() & AliESDtrack::kTOFout) || !(prongtracks[iProng]->GetStatus() & AliESDtrack::kTIME)) {
          rawPID[iProng][kTOF] = -1.;
        }
        else {
          float len = prongtracks[iProng]->GetIntegratedLength();
          if (len < 350.f) {
            rawPID[iProng][kTOF] = -1.;
          }
          else {
            rawPID[iProng][kTOF] = prongtracks[iProng]->GetTOFsignal();
            float time0 = pidrespo->GetTOFResponse().GetStartTime(prongtracks[iProng]->P());
            rawPID[iProng][kTOF] -= time0;
          }
        }
      }
    }
  }

  //fill PID arrays for different options
  switch(fPidOpt) {
    case 1: //kNsigmaPID
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(int iDet=kTPC; iDet<=kTOF; iDet++) {
          if(!useDet[iDet]) continue;
          for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
            if(!useHypo[iPartHypo]) continue;
            fPIDNsigmaVector[iProng][iDet][iPartHypo].push_back(sig[iProng][iDet][iPartHypo]);
          }
        }
      }
    break;
    case 2: //kNsigmaPIDint
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(int iDet=kTPC; iDet<=kTOF; iDet++) {
          if(!useDet[iDet]) continue;
          for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
            if(!useHypo[iPartHypo]) continue;
            fPIDNsigmaIntVector[iProng][iDet][iPartHypo].push_back(RoundFloatToInt(sig[iProng][iDet][iPartHypo]*100));
          }
        }
      }
    break;
    case 3: //kNsigmaPIDfloatandint
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(int iDet=kTPC; iDet<=kTOF; iDet++) {
          if(!useDet[iDet]) continue;
          for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
            if(!useHypo[iPartHypo]) continue;
            fPIDNsigmaVector[iProng][iDet][iPartHypo].push_back(sig[iProng][iDet][iPartHypo]*100);
            fPIDNsigmaIntVector[iProng][iDet][iPartHypo].push_back(RoundFloatToInt(sig[iProng][iDet][iPartHypo]*100));
          }
        }
      }
    break;
    case 4: //kNsigmaCombPID
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
          if(!useHypo[iPartHypo]) continue;
          fPIDNsigmaVector[iProng][0][iPartHypo].push_back(sigComb[iProng][iPartHypo]);
        }
      }
    break;
    case 5: //kNsigmaCombPIDint
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
          if(!useHypo[iPartHypo]) continue;
          fPIDNsigmaIntVector[iProng][0][iPartHypo].push_back(RoundFloatToInt(sigComb[iProng][iPartHypo]*100));
        }
      }
    break;
    case 6: //kNsigmaCombPIDfloatandint
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
          if(!useHypo[iPartHypo]) continue;
          fPIDNsigmaVector[iProng][0][iPartHypo].push_back(sigComb[iProng][iPartHypo]*100);
          fPIDNsigmaIntVector[iProng][0][iPartHypo].push_back(RoundFloatToInt(sigComb[iProng][iPartHypo]*100));
        }
      }
    break;
    case 7: //kRawPID
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(int iDet=kTPC; iDet<=kTOF; iDet++) {
          if(!useDet[iDet]) continue;
          fPIDrawVector[iProng][iDet].push_back(rawPID[iProng][iDet]);
        }
        if(useTPC) fTPCPProng[iProng].push_back(prongtracks[iProng]->GetTPCmomentum());
        if(useTOF) {
          fTOFPProng[iProng].push_back(GetTOFmomentum(prongtracks[iProng],pidrespo));
          fTrackIntegratedLengthProng[iProng].push_back(prongtracks[iProng]->GetIntegratedLength());
          fStartTimeResProng[iProng].push_back(pidrespo->GetTOFResponse().GetStartTimeRes(prongtracks[iProng]->P()));
        }
      }
    break;
    case 8: //kRawAndNsigmaPID
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(int iDet=kTPC; iDet<=kTOF; iDet++) {
          if(!useDet[iDet]) continue;
          for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
            if(!useHypo[iPartHypo]) continue;
            fPIDNsigmaVector[iProng][iDet][iPartHypo].push_back(sig[iProng][iDet][iPartHypo]);
          }
          fPIDrawVector[iProng][iDet].push_back(rawPID[iProng][iDet]);
        }
        if(useTPC) fTPCPProng[iProng].push_back(prongtracks[iProng]->GetTPCmomentum());
        if(useTOF) {
          fTOFPProng[iProng].push_back(GetTOFmomentum(prongtracks[iProng],pidrespo));
          fTrackIntegratedLengthProng[iProng].push_back(prongtracks[iProng]->GetIntegratedLength());
          fStartTimeResProng[iProng].push_back(pidrespo->GetTOFResponse().GetStartTimeRes(prongtracks[iProng]->P()));
        }
      }
    break;
    default:
      AliWarning("Wrong PID setting!");
      return false;
    break;
  }

  return true;
}

//________________________________________________________________
void AliHFTreeHandler::ResetDmesonCommonVarVectors() {
  
  fRunNumber.clear();
  fEvID.clear();
  fCandType.clear();
  fInvMass.clear();
  fPt.clear();
  fY.clear();
  fEta.clear();
  fPhi.clear();
  fDecayLength.clear();
  fDecayLengthXY.clear();
  fNormDecayLengthXY.clear();
  fCosP.clear();
  fCosPXY.clear();
  fImpParXY.clear();
  fDCA.clear();
}

//________________________________________________________________
void AliHFTreeHandler::ResetMCGenVectors() {
  
  fRunNumber.clear();
  fEvID.clear();
  fCandType.clear();
  fInvMass.clear();
  fPt.clear();
  fY.clear();
  fEta.clear();
  fPhi.clear();
  fDauInAcceptance.clear();
}

//________________________________________________________________
void AliHFTreeHandler::ResetSingleTrackVarVectors() {
  
  if(fSingleTrackOpt==kNoSingleTrackVars) return;
  
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {

    if(fSingleTrackOpt==kRedSingleTrackVars) {
      fPtProng[iProng].clear();
      fEtaProng[iProng].clear();
      fPhiProng[iProng].clear();
    }
    else if(fSingleTrackOpt==kAllSingleTrackVars) {
      fPtProng[iProng].clear();
      fEtaProng[iProng].clear();
      fPhiProng[iProng].clear();
      fPProng[iProng].clear();
      fNTPCclsProng[iProng].clear();
      fNTPCclsPidProng[iProng].clear();
      fNTPCCrossedRowProng[iProng].clear();
      fChi2perNDFProng[iProng].clear();
      fNITSclsProng[iProng].clear();
      fITSclsMapProng[iProng].clear();
    }
  }
}

//________________________________________________________________
void AliHFTreeHandler::ResetPidVarVectors() {

  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    for(unsigned int iDet=0; iDet<knMaxDet4Pid; iDet++) {
      for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
        fPIDNsigmaVector[iProng][iDet][iPartHypo].clear();
        fPIDNsigmaIntVector[iProng][iDet][iPartHypo].clear();
      }
      fPIDrawVector[iProng][iDet].clear();
    }
    fTPCPProng[iProng].clear();
    fTOFPProng[iProng].clear();
    fTrackIntegratedLengthProng[iProng].clear();
    fStartTimeResProng[iProng].clear();
  }
}

//________________________________________________________________
double AliHFTreeHandler::CombineNsigmaDiffDet(double nsigmaTPC, double nsigmaTOF)
{
  if(nsigmaTPC > -998. && nsigmaTOF > -998.) return TMath::Sqrt((nsigmaTPC*nsigmaTPC+nsigmaTOF*nsigmaTOF)/2);
  else if(nsigmaTPC > -998. && nsigmaTOF < -998.) return nsigmaTPC;
  else if(nsigmaTPC > -998. && nsigmaTOF < -998.) return nsigmaTPC;
  else return -999.;
}

//________________________________________________________________
int AliHFTreeHandler::RoundFloatToInt(double num) 
{
  if(num>=static_cast<double>(std::numeric_limits<int>::max())) return std::numeric_limits<int>::max();
  else if(num<=static_cast<double>(std::numeric_limits<int>::min())) return std::numeric_limits<int>::min();
 
  return static_cast<int>(num);
}

//________________________________________________________________
float AliHFTreeHandler::ComputeMaxd0MeasMinusExp(AliAODRecoDecayHF* cand, float bfield)
{
  float dd0max=0;
  unsigned int fNProngs_cand = (unsigned int)cand->GetNProngs();
  for(unsigned int iProng=0; iProng<fNProngs_cand; iProng++) {
    double d0diff, errd0diff;
    cand->Getd0MeasMinusExpProng(iProng,bfield,d0diff,errd0diff);
    float normdd0 = d0diff/errd0diff;
    if(iProng==0) dd0max=normdd0;
    else if(TMath::Abs(normdd0)>TMath::Abs(dd0max)) dd0max=normdd0;
  }
  return dd0max;
}

//________________________________________________________________
float AliHFTreeHandler::GetTOFmomentum(AliAODTrack* track, AliPIDResponse* pidrespo)
{
  float t_d = pidrespo->GetTOFResponse().GetExpectedSignal(track, AliPID::kTriton); //largest mass possible with Z=1
  float len = track->GetIntegratedLength();
  float beta_d = len / (t_d * kCSPEED);
  float mass = AliPID::ParticleMassZ(AliPID::kTriton); //largest mass possible with Z=1

  if(TMath::Abs(beta_d-1.) < 1.e-12) return track->GetTPCmomentum();
  else return mass*beta_d/sqrt(1.-(beta_d*beta_d));
}
