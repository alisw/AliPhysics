/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerApply
// \brief helper class to handle a tree for cut optimisation and MVA analyses
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
/////////////////////////////////////////////////////////////

#include <cmath>
#include <limits>
#include "AliHFTreeHandlerApply.h"
#include "AliPID.h"
#include "AliAODRecoDecayHF.h"
#include "AliPIDResponse.h"
#include "AliESDtrack.h"
#include "TMath.h"

/// \cond CLASSIMP
ClassImp(AliHFTreeHandlerApply);
/// \endcond

//________________________________________________________________
AliHFTreeHandlerApply::AliHFTreeHandlerApply():
AliHFTreeHandlerApply(kNsigmaPID)
{
  //
  // Default constructor
  //
}

//________________________________________________________________
AliHFTreeHandlerApply::AliHFTreeHandlerApply(int PIDopt):
TObject(),
fTreeVar(nullptr),
fNProngs(-1),
fCandType(0),
fInvMass(-9999.),
fPt(-9999.),
fPtGen(-9999.),
fY(-9999.),
fEta(-9999.),
fPhi(-9999.),
fMLProb(-9999.),
fDecayLength(-9999.),
fDecayLengthXY(-9999.),
fNormDecayLengthXY(-9999.),
fCosP(-9999.),
fCosPXY(-9999.),
fImpParXY(-9999.),
fDCA(-9999.),
fPProng{},
fSPDhitsProng{},
fTPCPProng{},
fTOFPProng{},
fPtProng{},
fEtaProng{},
fPhiProng{},
fNTPCclsProng{},
fNTPCclsPidProng{},
fNTPCCrossedRowProng{},
fChi2perNDFProng{},
fNITSclsProng{},
fITSclsMapProng{},
fTrackIntegratedLengthProng{},
fStartTimeResProng{},
fPIDNsigmaVector{},
fPIDNsigmaIntVector{},
fPIDrawVector{},
fPidOpt(PIDopt),
fSingleTrackOpt(kRedSingleTrackVars),
fFillOnlySignal(false),
fIsMCGenTree(false),
fDauInAcceptance(false),
fEvID(9999),
fRunNumber(9999),
fRunNumberPrevCand(9999),
fApplyNsigmaTPCDataCorr(false),
fSystNsigmaTPCDataCorr(AliAODPidHF::kNone),
fMeanNsigmaTPCPionData{},
fMeanNsigmaTPCKaonData{},
fMeanNsigmaTPCProtonData{},
fSigmaNsigmaTPCPionData{},
fSigmaNsigmaTPCKaonData{},
fSigmaNsigmaTPCProtonData{},
fPlimitsNsigmaTPCDataCorr{},
fNPbinsNsigmaTPCDataCorr(0),
fEtalimitsNsigmaTPCDataCorr{},
fNEtabinsNsigmaTPCDataCorr(0)
{
  //
  // Default constructor
  //
  for(unsigned int iProng=0; iProng<knMaxProngs; iProng++) {
    fPProng[iProng] = -9999.;
    fTPCPProng[iProng] = -9999.;
    fTOFPProng[iProng] = -9999.;
    fPtProng[iProng] = -9999.;
    fEtaProng[iProng] = -9999.;
    fPhiProng[iProng] = -9999.;
    fNTPCclsProng[iProng] = -9999;
    fNTPCclsPidProng[iProng] = -9999;
    fNTPCCrossedRowProng[iProng] = -9999.;
    fChi2perNDFProng[iProng] = -9999.;
    fNITSclsProng[iProng] = -9999;
    fITSclsMapProng[iProng] = -9999;
    fTrackIntegratedLengthProng[iProng] = -9999.;
    fStartTimeResProng[iProng] = -9999.;
    for(unsigned int iDet=0; iDet<knMaxDet4Pid; iDet++)
      fPIDrawVector[iProng][iDet] = -999.;
    for(unsigned int iDet=0; iDet<knMaxDet4Pid+1; iDet++) {
      for(unsigned int iHypo=0; iHypo<knMaxHypo4Pid; iHypo++) {
        fPIDNsigmaVector[iProng][iDet][iHypo] = -999.;
        fPIDNsigmaIntVector[iProng][iDet][iHypo] = -999;
      }
    }
  }
  
  for(int iP=0; iP<=AliAODPidHF::kMaxPBins; iP++) {
    fPlimitsNsigmaTPCDataCorr[iP] = 0.;
  }
  for(int iEta=0; iEta<=AliAODPidHF::kMaxEtaBins; iEta++) {
    fEtalimitsNsigmaTPCDataCorr[iEta] = 0.;
  }
}

//________________________________________________________________
AliHFTreeHandlerApply::~AliHFTreeHandlerApply()
{
  //
  // Destructor
  //
  
  delete fTreeVar;
}

//________________________________________________________________
TTree* AliHFTreeHandlerApply::BuildTreeMCGen(TString name, TString title) {
  
  fIsMCGenTree = true;
  
  if(fTreeVar) {
    delete fTreeVar;
    fTreeVar=nullptr;
  }
  fTreeVar = new TTree(name.Data(),title.Data());
  fTreeVar->Branch("run_number",&fRunNumber);
  fTreeVar->Branch("ev_id",&fEvID);
  fTreeVar->Branch("cand_type",&fCandType);
  fTreeVar->Branch("pt_cand",&fPt);
  fTreeVar->Branch("y_cand",&fY);
  fTreeVar->Branch("eta_cand",&fEta);
  fTreeVar->Branch("phi_cand",&fPhi);
  fTreeVar->Branch("dau_in_acc",&fDauInAcceptance);
  
  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandlerApply::SetMCGenVariables(int runnumber, unsigned int eventID, AliAODMCParticle* mcpart) {
  
  if(!mcpart) return false;
  if(!(fCandType&kSignal)) return true; // fill only signal in the generated
  
  fRunNumber = runnumber;
  fEvID = eventID;
  fPt = mcpart->Pt();
  fY = mcpart->Y();
  fEta = mcpart->Eta();
  fPhi = mcpart->Phi();
  
  return true;
}

//________________________________________________________________
void AliHFTreeHandlerApply::SetCandidateType(bool issignal, bool isbkg, bool isprompt, bool isFD, bool isreflected)
{  
  if(issignal) fCandType |= kSignal;
  else fCandType &= ~kSignal;
  if(isbkg && !fIsMCGenTree) fCandType |= kBkg;
  else fCandType &= ~kBkg;
  if(isprompt) fCandType |= kPrompt;
  else fCandType &= ~kPrompt;
  if(isFD) fCandType |= kFD;
  else fCandType &= ~kFD;
  if(isreflected && !fIsMCGenTree) fCandType |= kRefl;
  else fCandType &= ~kRefl;
}

//________________________________________________________________
void AliHFTreeHandlerApply::AddCommonDmesonVarBranches(Bool_t HasSecVtx) {
  
  fTreeVar->Branch("run_number",&fRunNumber);
  fTreeVar->Branch("ev_id",&fEvID);
  fTreeVar->Branch("cand_type",&fCandType);
  fTreeVar->Branch("inv_mass",&fInvMass);
  fTreeVar->Branch("pt_cand",&fPt);
  fTreeVar->Branch("pt_gen_cand",&fPtGen);
  fTreeVar->Branch("y_cand",&fY);
  fTreeVar->Branch("eta_cand",&fEta);
  fTreeVar->Branch("phi_cand",&fPhi);
  fTreeVar->Branch("ml_prob",&fMLProb);
  if(HasSecVtx){
    fTreeVar->Branch("d_len",&fDecayLength);
    fTreeVar->Branch("d_len_xy",&fDecayLengthXY);
    fTreeVar->Branch("norm_dl_xy",&fNormDecayLengthXY);
    fTreeVar->Branch("cos_p",&fCosP);
    fTreeVar->Branch("cos_p_xy",&fCosPXY);
    fTreeVar->Branch("imp_par_xy",&fImpParXY);
    fTreeVar->Branch("dca",&fDCA);
  }
} 

//________________________________________________________________
void AliHFTreeHandlerApply::AddSingleTrackBranches() {
  
  if(fSingleTrackOpt==kNoSingleTrackVars) return;
  
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    
    if(fSingleTrackOpt==kRedSingleTrackVars) {
      fTreeVar->Branch(Form("pt_prong%d",iProng),&fPtProng[iProng]);
      fTreeVar->Branch(Form("eta_prong%d",iProng),&fEtaProng[iProng]);
      fTreeVar->Branch(Form("phi_prong%d",iProng),&fPhiProng[iProng]);
      fTreeVar->Branch(Form("p_prong%d",iProng),&fPProng[iProng]);
      fTreeVar->Branch(Form("spdhits_prong%d",iProng),&fSPDhitsProng[iProng]);
    }
    else if(fSingleTrackOpt==kAllSingleTrackVars) {
      fTreeVar->Branch(Form("pt_prong%d",iProng),&fPtProng[iProng]);
      fTreeVar->Branch(Form("eta_prong%d",iProng),&fEtaProng[iProng]);
      fTreeVar->Branch(Form("phi_prong%d",iProng),&fPhiProng[iProng]);
      fTreeVar->Branch(Form("p_prong%d",iProng),&fPProng[iProng]);
      fTreeVar->Branch(Form("spdhits_prong%d",iProng),&fSPDhitsProng[iProng]);
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
void AliHFTreeHandlerApply::AddPidBranches(bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF)
{
  
  if(fPidOpt==kNoPID) return;
  if(fPidOpt>kNsigmaDetAndCombPID) {
    AliWarning("Wrong PID setting!");
    return;
  }
  
  bool useHypo[knMaxHypo4Pid] = {usePionHypo,useKaonHypo,useProtonHypo};
  bool useDet[knMaxDet4Pid] = {useTPC,useTOF};
  TString partHypoName[knMaxHypo4Pid] = {"Pi","K","Pr"};
  TString detName[knMaxDet4Pid] = {"TPC","TOF"};
  TString rawPidName[knMaxDet4Pid] = {"dEdxTPC","ToF"};
  
  for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
    if((fPidOpt>=kNsigmaPID && fPidOpt<=kNsigmaPIDfloatandint) || fPidOpt>=kRawAndNsigmaPID) {
      for(unsigned int iDet=0; iDet<knMaxDet4Pid; iDet++) {
        if(!useDet[iDet]) continue;
        for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
          if(!useHypo[iPartHypo]) continue;
          if(fPidOpt==kNsigmaPID || fPidOpt==kNsigmaPIDfloatandint || fPidOpt>=kRawAndNsigmaPID)
            fTreeVar->Branch(Form("nsig%s_%s_%d",detName[iDet].Data(),partHypoName[iPartHypo].Data(),iProng),&fPIDNsigmaVector[iProng][iDet][iPartHypo]);
          if(fPidOpt==kNsigmaPIDint || fPidOpt==kNsigmaPIDfloatandint)
            fTreeVar->Branch(Form("int_nsig%s_%s_%d",detName[iDet].Data(),partHypoName[iPartHypo].Data(),iProng),&fPIDNsigmaIntVector[iProng][iDet][iPartHypo]);
        }
      }
    }
    if((fPidOpt>=kNsigmaCombPID && fPidOpt<=kNsigmaCombPIDfloatandint) || fPidOpt==kNsigmaDetAndCombPID) {
      for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
        if(!useHypo[iPartHypo]) continue;
        if(fPidOpt==kNsigmaCombPID || fPidOpt==kNsigmaCombPIDfloatandint || fPidOpt==kNsigmaDetAndCombPID)
          fTreeVar->Branch(Form("nsigComb_%s_%d",partHypoName[iPartHypo].Data(),iProng),&fPIDNsigmaVector[iProng][kCombTPCTOF][iPartHypo]);
        if(fPidOpt==kNsigmaCombPIDint || fPidOpt==kNsigmaCombPIDfloatandint)
          fTreeVar->Branch(Form("int_nsigComb_%s_%d",partHypoName[iPartHypo].Data(),iProng),&fPIDNsigmaIntVector[iProng][kCombTPCTOF][iPartHypo]);
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
bool AliHFTreeHandlerApply::SetSingleTrackVars(AliAODTrack* prongtracks[]) {
  
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
      fPtProng[iProng]=prongtracks[iProng]->Pt();
      fEtaProng[iProng]=prongtracks[iProng]->Eta();
      fPhiProng[iProng]=prongtracks[iProng]->Phi();
      fPProng[iProng]=prongtracks[iProng]->P();
      fSPDhitsProng[iProng] = prongtracks[iProng]->GetITSClusterMap() & 0x3;
    }
    else if(fSingleTrackOpt==kAllSingleTrackVars) {
      fPtProng[iProng]=prongtracks[iProng]->Pt();
      fEtaProng[iProng]=prongtracks[iProng]->Eta();
      fPhiProng[iProng]=prongtracks[iProng]->Phi();
      fPProng[iProng]=prongtracks[iProng]->P();
      fSPDhitsProng[iProng] = prongtracks[iProng]->GetITSClusterMap() & 0x3;
      fNTPCclsProng[iProng]=prongtracks[iProng]->GetTPCNcls();
      fNTPCclsPidProng[iProng]=prongtracks[iProng]->GetTPCsignalN();
      fNTPCCrossedRowProng[iProng]=prongtracks[iProng]->GetTPCNCrossedRows();
      fChi2perNDFProng[iProng]=prongtracks[iProng]->Chi2perNDF();
      fNITSclsProng[iProng]=prongtracks[iProng]->GetITSNcls();
      fITSclsMapProng[iProng]=static_cast<int>(prongtracks[iProng]->GetITSClusterMap());
    }
  }
  
  return true;
}

//________________________________________________________________
bool AliHFTreeHandlerApply::SetPidVars(AliAODTrack* prongtracks[], AliPIDResponse* pidrespo, bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF)
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
    if((fPidOpt>=kNsigmaPID && fPidOpt<=kNsigmaCombPIDfloatandint) || fPidOpt>=kRawAndNsigmaPID) {
      for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
        if(useHypo[iPartHypo]) {
          if(useTPC) {
            float nSigmaTPC = pidrespo->NumberOfSigmasTPC(prongtracks[iProng],parthypo[iPartHypo]);
            if(fApplyNsigmaTPCDataCorr && nSigmaTPC>-990.) {
              float sigma=1., mean=0.;
              GetNsigmaTPCMeanSigmaData(mean, sigma, parthypo[iPartHypo], prongtracks[iProng]->GetTPCmomentum(), prongtracks[iProng]->Eta());
              nSigmaTPC = (nSigmaTPC-mean)/sigma;
            }
            sig[iProng][kTPC][iPartHypo] = nSigmaTPC;
          }
          if(useTOF) sig[iProng][kTOF][iPartHypo] = pidrespo->NumberOfSigmasTOF(prongtracks[iProng],parthypo[iPartHypo]);
          if(((fPidOpt>=kNsigmaCombPID && fPidOpt<=kNsigmaCombPIDfloatandint) || fPidOpt==kNsigmaDetAndCombPID) && useTPC && useTOF) {
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
            fPIDNsigmaVector[iProng][iDet][iPartHypo]=sig[iProng][iDet][iPartHypo];
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
            fPIDNsigmaIntVector[iProng][iDet][iPartHypo]=RoundFloatToInt(sig[iProng][iDet][iPartHypo]*100);
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
            fPIDNsigmaVector[iProng][iDet][iPartHypo]=sig[iProng][iDet][iPartHypo]*100;
            fPIDNsigmaIntVector[iProng][iDet][iPartHypo]=RoundFloatToInt(sig[iProng][iDet][iPartHypo]*100);
          }
        }
      }
      break;
    case 4: //kNsigmaCombPID
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
          if(!useHypo[iPartHypo]) continue;
          fPIDNsigmaVector[iProng][kCombTPCTOF][iPartHypo]=sigComb[iProng][iPartHypo];
        }
      }
      break;
    case 5: //kNsigmaCombPIDint
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
          if(!useHypo[iPartHypo]) continue;
          fPIDNsigmaIntVector[iProng][kCombTPCTOF][iPartHypo]=RoundFloatToInt(sigComb[iProng][iPartHypo]*100);
        }
      }
      break;
    case 6: //kNsigmaCombPIDfloatandint
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
          if(!useHypo[iPartHypo]) continue;
          fPIDNsigmaVector[iProng][kCombTPCTOF][iPartHypo]=sigComb[iProng][iPartHypo]*100;
          fPIDNsigmaIntVector[iProng][kCombTPCTOF][iPartHypo]=RoundFloatToInt(sigComb[iProng][iPartHypo]*100);
        }
      }
      break;
    case 7: //kRawPID
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(int iDet=kTPC; iDet<=kTOF; iDet++) {
          if(!useDet[iDet]) continue;
          fPIDrawVector[iProng][iDet]=rawPID[iProng][iDet];
        }
        if(useTPC) fTPCPProng[iProng]=prongtracks[iProng]->GetTPCmomentum();
        if(useTOF) {
          fTOFPProng[iProng]=GetTOFmomentum(prongtracks[iProng],pidrespo);
          fTrackIntegratedLengthProng[iProng]=prongtracks[iProng]->GetIntegratedLength();
          fStartTimeResProng[iProng]=pidrespo->GetTOFResponse().GetStartTimeRes(prongtracks[iProng]->P());
        }
      }
      break;
    case 8: //kRawAndNsigmaPID
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(int iDet=kTPC; iDet<=kTOF; iDet++) {
          if(!useDet[iDet]) continue;
          for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
            if(!useHypo[iPartHypo]) continue;
            fPIDNsigmaVector[iProng][iDet][iPartHypo]=sig[iProng][iDet][iPartHypo];
          }
          fPIDrawVector[iProng][iDet]=rawPID[iProng][iDet];
        }
        if(useTPC) fTPCPProng[iProng]=prongtracks[iProng]->GetTPCmomentum();
        if(useTOF) {
          fTOFPProng[iProng]=GetTOFmomentum(prongtracks[iProng],pidrespo);
          fTrackIntegratedLengthProng[iProng]=prongtracks[iProng]->GetIntegratedLength();
          fStartTimeResProng[iProng]=pidrespo->GetTOFResponse().GetStartTimeRes(prongtracks[iProng]->P());
        }
      }
      break;
    case 9: //kNsigmaDetAndCombPID
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
          if(!useHypo[iPartHypo]) continue;
          fPIDNsigmaVector[iProng][kCombTPCTOF][iPartHypo]=sigComb[iProng][iPartHypo];
          for(int iDet=kTPC; iDet<=kTOF; iDet++) {
            if(!useDet[iDet]) continue;
            fPIDNsigmaVector[iProng][iDet][iPartHypo]=sig[iProng][iDet][iPartHypo];
          }
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
double AliHFTreeHandlerApply::CombineNsigmaDiffDet(double nsigmaTPC, double nsigmaTOF)
{
  if(nsigmaTPC > -998. && nsigmaTOF > -998.) return TMath::Sqrt((nsigmaTPC*nsigmaTPC+nsigmaTOF*nsigmaTOF)/2);
  else if(nsigmaTPC > -998. && nsigmaTOF < -998.) return TMath::Abs(nsigmaTPC);
  else if(nsigmaTPC < -998. && nsigmaTOF > -998.) return TMath::Abs(nsigmaTOF);
  else return -999.;
}

//________________________________________________________________
int AliHFTreeHandlerApply::RoundFloatToInt(double num)
{
  if(num>=static_cast<double>(std::numeric_limits<int>::max())) return std::numeric_limits<int>::max();
  else if(num<=static_cast<double>(std::numeric_limits<int>::min())) return std::numeric_limits<int>::min();
  
  return static_cast<int>(num);
}

//________________________________________________________________
float AliHFTreeHandlerApply::ComputeMaxd0MeasMinusExp(AliAODRecoDecayHF* cand, float bfield)
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
float AliHFTreeHandlerApply::GetTOFmomentum(AliAODTrack* track, AliPIDResponse* pidrespo)
{
  float t_d = pidrespo->GetTOFResponse().GetExpectedSignal(track, AliPID::kTriton); //largest mass possible with Z=1
  float len = track->GetIntegratedLength();
  float beta_d = len / (t_d * kCSPEED);
  float mass = AliPID::ParticleMassZ(AliPID::kTriton); //largest mass possible with Z=1
  
  if(TMath::Abs(beta_d-1.) < 1.e-12) return track->GetTPCmomentum();
  else return mass*beta_d/sqrt(1.-(beta_d*beta_d));
}

//________________________________________________________________
void AliHFTreeHandlerApply::GetNsigmaTPCMeanSigmaData(float &mean, float &sigma, AliPID::EParticleType species, float pTPC, float eta) {
  
  if(fRunNumber!=fRunNumberPrevCand)
    AliAODPidHF::SetNsigmaTPCDataDrivenCorrection(fRunNumber, fSystNsigmaTPCDataCorr, fNPbinsNsigmaTPCDataCorr, fPlimitsNsigmaTPCDataCorr,
                                                  fNEtabinsNsigmaTPCDataCorr, fEtalimitsNsigmaTPCDataCorr, fMeanNsigmaTPCPionData, fMeanNsigmaTPCKaonData,
                                                  fMeanNsigmaTPCProtonData, fSigmaNsigmaTPCPionData, fSigmaNsigmaTPCKaonData, fSigmaNsigmaTPCProtonData);
  
  int bin = TMath::BinarySearch(fNPbinsNsigmaTPCDataCorr,fPlimitsNsigmaTPCDataCorr,pTPC);
  if(bin<0) bin=0; //underflow --> equal to min value
  else if(bin>fNPbinsNsigmaTPCDataCorr-1) bin=fNPbinsNsigmaTPCDataCorr-1; //overflow --> equal to max value
  
  int etabin = TMath::BinarySearch(fNEtabinsNsigmaTPCDataCorr,fEtalimitsNsigmaTPCDataCorr,TMath::Abs(eta));
  if(etabin<0) etabin=0; //underflow --> equal to min value
  else if(etabin>fNEtabinsNsigmaTPCDataCorr-1) etabin=fNEtabinsNsigmaTPCDataCorr-1; //overflow --> equal to max value
  
  switch(species) {
    case AliPID::kPion:
    {
      mean = fMeanNsigmaTPCPionData[etabin][bin];
      sigma = fSigmaNsigmaTPCPionData[etabin][bin];
      break;
    }
    case AliPID::kKaon:
    {
      mean = fMeanNsigmaTPCKaonData[etabin][bin];
      sigma = fSigmaNsigmaTPCKaonData[etabin][bin];
      break;
    }
    case AliPID::kProton:
    {
      mean = fMeanNsigmaTPCProtonData[etabin][bin];
      sigma = fSigmaNsigmaTPCProtonData[etabin][bin];
      break;
    }
    default:
    {
      mean = 0.;
      sigma = 1.;
      break;
    }
  }
}
