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
#include <array>
#include "AliHFTreeHandler.h"
#include "AliPID.h"
#include "AliAODRecoDecayHF.h"
#include "AliPIDResponse.h"
#include "AliESDtrack.h"
#include "TMath.h"

using std::array;

/// \cond CLASSIMP
ClassImp(AliHFTreeHandler);
/// \endcond

//________________________________________________________________
AliHFTreeHandler::AliHFTreeHandler():
  TObject(),
  fTreeVar(nullptr),
  fNProngs(-1),
  fCandType(0),
  fInvMass(-9999.),
  fPt(-9999.),
  fY(-9999.),
  fEta(-9999.),
  fPhi(-9999.),
  fDecayLength(-9999.),
  fDecayLengthXY(-9999.),
  fNormDecayLengthXY(-9999.),
  fCosP(-9999.),
  fCosPXY(-9999.),
  fImpParXY(-9999.),
  fDCA(-9999.),
  fPidOpt(kNsigmaPID),
  fSingleTrackOpt(kRedSingleTrackVars),
  fFillOnlySignal(false),
  fIsMCGenTree(false),
  fDauInAcceptance(false),
  fEvID(9999),
  fRunNumber(9999),
  fRunNumberPrevCand(9999),
  fApplyNsigmaTPCDataCorr(false),
  fSystNsigmaTPCDataCorr(kNone),
  fMeanNsigmaTPCPionData{},
  fMeanNsigmaTPCKaonData{},
  fMeanNsigmaTPCProtonData{},
  fSigmaNsigmaTPCPionData{},
  fSigmaNsigmaTPCKaonData{},
  fSigmaNsigmaTPCProtonData{},
  fPlimitsNsigmaTPCDataCorr{},
  fNPbinsNsigmaTPCDataCorr(0)
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
    for(unsigned int iDet=0; iDet<knMaxDet4Pid; iDet++) {
      fPIDrawVector[iProng][iDet] = -9999.;
      for(unsigned int iHypo=0; iHypo<knMaxHypo4Pid; iHypo++) {
        fPIDNsigmaVector[iProng][iDet][iHypo] = -9999.;
        fPIDNsigmaIntVector[iProng][iDet][iHypo] = -999999;      
      }
    }
  }
}

//________________________________________________________________
AliHFTreeHandler::AliHFTreeHandler(int PIDopt):
    TObject(),
  fTreeVar(nullptr),
  fNProngs(-1),
  fCandType(0),
  fInvMass(-9999.),
  fPt(-9999.),
  fY(-9999.),
  fEta(-9999.),
  fPhi(-9999.),
  fDecayLength(-9999.),
  fDecayLengthXY(-9999.),
  fNormDecayLengthXY(-9999.),
  fCosP(-9999.),
  fCosPXY(-9999.),
  fImpParXY(-9999.),
  fDCA(-9999.),
  fPidOpt(PIDopt),
  fSingleTrackOpt(kRedSingleTrackVars),
  fFillOnlySignal(false),
  fIsMCGenTree(false),
  fDauInAcceptance(false),
  fEvID(9999),
  fRunNumber(9999),
  fRunNumberPrevCand(9999),
  fApplyNsigmaTPCDataCorr(false),
  fSystNsigmaTPCDataCorr(kNone),
  fMeanNsigmaTPCPionData{},
  fMeanNsigmaTPCKaonData{},
  fMeanNsigmaTPCProtonData{},
  fSigmaNsigmaTPCPionData{},
  fSigmaNsigmaTPCKaonData{},
  fSigmaNsigmaTPCProtonData{},
  fPlimitsNsigmaTPCDataCorr{},
  fNPbinsNsigmaTPCDataCorr(0)
{
  //
  // Standard constructor
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
    for(unsigned int iDet=0; iDet<knMaxDet4Pid; iDet++) {
      fPIDrawVector[iProng][iDet] = -9999.;
      for(unsigned int iHypo=0; iHypo<knMaxHypo4Pid; iHypo++) {
        fPIDNsigmaVector[iProng][iDet][iHypo] = -9999.;
        fPIDNsigmaIntVector[iProng][iDet][iHypo] = -999999;      
      }
    }
  }
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
bool AliHFTreeHandler::SetMCGenVariables(int runnumber, unsigned int eventID, AliAODMCParticle* mcpart) {

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
void AliHFTreeHandler::SetCandidateType(bool issignal, bool isbkg, bool isprompt, bool isFD, bool isreflected) 
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
void AliHFTreeHandler::AddCommonDmesonVarBranches() {

  fTreeVar->Branch("run_number",&fRunNumber);
  fTreeVar->Branch("ev_id",&fEvID);
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
      fPtProng[iProng]=prongtracks[iProng]->Pt();
      fEtaProng[iProng]=prongtracks[iProng]->Eta();
      fPhiProng[iProng]=prongtracks[iProng]->Phi();
    }
    else if(fSingleTrackOpt==kAllSingleTrackVars) {
      fPtProng[iProng]=prongtracks[iProng]->Pt();
      fEtaProng[iProng]=prongtracks[iProng]->Eta();
      fPhiProng[iProng]=prongtracks[iProng]->Phi();
      fPProng[iProng]=prongtracks[iProng]->P();
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
          if(useTPC) {
          float nSigmaTPC = pidrespo->NumberOfSigmasTPC(prongtracks[iProng],parthypo[iPartHypo]);
            if(fApplyNsigmaTPCDataCorr) {
              float sigma=1., mean=0.;
              GetNsigmaTPCMeanSigmaData(mean, sigma, parthypo[iPartHypo], prongtracks[iProng]->GetTPCmomentum());
              nSigmaTPC = (nSigmaTPC-mean)/sigma;
            }
            sig[iProng][kTPC][iPartHypo] = nSigmaTPC;
          }
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
          fPIDNsigmaVector[iProng][0][iPartHypo]=sigComb[iProng][iPartHypo];
        }
      }
    break;
    case 5: //kNsigmaCombPIDint
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
          if(!useHypo[iPartHypo]) continue;
          fPIDNsigmaIntVector[iProng][0][iPartHypo]=RoundFloatToInt(sigComb[iProng][iPartHypo]*100);
        }
      }
    break;
    case 6: //kNsigmaCombPIDfloatandint
      for(unsigned int iProng=0; iProng<fNProngs; iProng++) {
        for(unsigned int iPartHypo=0; iPartHypo<knMaxHypo4Pid; iPartHypo++) {
          if(!useHypo[iPartHypo]) continue;
          fPIDNsigmaVector[iProng][0][iPartHypo]=sigComb[iProng][iPartHypo]*100;
          fPIDNsigmaIntVector[iProng][0][iPartHypo]=RoundFloatToInt(sigComb[iProng][iPartHypo]*100);
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
    default:
      AliWarning("Wrong PID setting!");
      return false;
    break;
  }

  return true;
}

//________________________________________________________________
double AliHFTreeHandler::CombineNsigmaDiffDet(double nsigmaTPC, double nsigmaTOF)
{
  if(nsigmaTPC > -998. && nsigmaTOF > -998.) return TMath::Sqrt((nsigmaTPC*nsigmaTPC+nsigmaTOF*nsigmaTOF)/2);
  else if(nsigmaTPC > -998. && nsigmaTOF < -998.) return nsigmaTPC;
  else if(nsigmaTPC > -998. && nsigmaTOF < -998.) return nsigmaTPC;
  else return -9999.;
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

//________________________________________________________________
void AliHFTreeHandler::SetNsigmaTPCDataCorr() {
  
    if(fRunNumber>=295585 && fRunNumber<=296623 && fSystNsigmaTPCDataCorr==kPbPb010) { //LHC18q 0-10%
    fNPbinsNsigmaTPCDataCorr = 8;
    array<float,9> pTPClims = {0.3,0.5,0.75,1.,1.5,2.,3.,5.,10.};    
    array<float,8> meanPion = {-0.476642, -0.611512, -0.70491, -0.785863, -0.858335, -0.913384, -0.926733, -1.03424};
    array<float,8> meanKaon = {-0.376284, -0.689586, -0.752243, -0.922438, -0.95792, -0.958785, -1.00629, -1.10473};
    array<float,8> meanProton = {-0.162057, -0.222369, -0.517459, -0.874908, -0.961924, -1.01193, -0.839815, -0.691694};
    array<float,8> sigmaPion = {0.98579, 0.962247, 0.945548, 0.920657, 0.909255, 0.957158, 0.907777, 0.954516};
    array<float,8> sigmaKaon = {0.851531, 0.909522, 0.96582, 0.900314, 0.887377, 0.880861, 0.848008, 0.916044};
    array<float,8> sigmaProton = {0.748482, 0.79806, 0.852967, 0.979616, 0.997911, 0.860067, 0.883535, 0.929892};
    std::copy(pTPClims.begin(),pTPClims.end(),fPlimitsNsigmaTPCDataCorr);
    std::copy(meanPion.begin(),meanPion.end(),fMeanNsigmaTPCPionData);
    std::copy(meanKaon.begin(),meanKaon.end(),fMeanNsigmaTPCKaonData);
    std::copy(meanProton.begin(),meanProton.end(),fMeanNsigmaTPCProtonData);
    std::copy(sigmaPion.begin(),sigmaPion.end(),fSigmaNsigmaTPCPionData);
    std::copy(sigmaKaon.begin(),sigmaKaon.end(),fSigmaNsigmaTPCKaonData);
    std::copy(sigmaProton.begin(),sigmaProton.end(),fSigmaNsigmaTPCProtonData);
  }
  else if(fRunNumber>=295585 && fRunNumber<=296623 && fSystNsigmaTPCDataCorr==kPbPb3050) { //LHC18q 30-50%
    fNPbinsNsigmaTPCDataCorr = 8;
    array<float,9> pTPClims = {0.3,0.5,0.75,1.,1.5,2.,3.,5.,10.};    
    array<float,8> meanPion = {-0.282783, -0.351074, -0.370549, -0.368398, -0.37546, -0.332551, -0.304001, -0.329724};
    array<float,8> meanKaon = {-0.147986, -0.26169, -0.339263, -0.566137, -0.619671, -0.586758, -0.430222, -0.243858};
    array<float,8> meanProton = {-0.116342, -0.204619, -0.18317, -0.288015, -0.431383, -0.496598, -0.476154, -0.464085};
    array<float,8> sigmaPion = {1.15354, 1.1191, 1.11418, 1.11474, 1.11121, 1.10621, 1.05918, 1.05666};
    array<float,8> sigmaKaon = {1.07676, 1.11978, 1.14182, 1.09804, 1.09674, 1.08182, 1.07091, 1.00419};
    array<float,8> sigmaProton = {1.07493, 1.11953, 1.14044, 1.14904, 1.12158, 1.09405, 1.06465, 1.11152};
    std::copy(pTPClims.begin(),pTPClims.end(),fPlimitsNsigmaTPCDataCorr);
    std::copy(meanPion.begin(),meanPion.end(),fMeanNsigmaTPCPionData);
    std::copy(meanKaon.begin(),meanKaon.end(),fMeanNsigmaTPCKaonData);
    std::copy(meanProton.begin(),meanProton.end(),fMeanNsigmaTPCProtonData);
    std::copy(sigmaPion.begin(),sigmaPion.end(),fSigmaNsigmaTPCPionData);
    std::copy(sigmaKaon.begin(),sigmaKaon.end(),fSigmaNsigmaTPCKaonData);
    std::copy(sigmaProton.begin(),sigmaProton.end(),fSigmaNsigmaTPCProtonData);
  }
  else if(fRunNumber>296690 && fRunNumber<297595 && fSystNsigmaTPCDataCorr==kPbPb010) { //LHC18r 0-10%
    fNPbinsNsigmaTPCDataCorr = 8;
    array<float,9> pTPClims = {0.3,0.5,0.75,1.,1.5,2.,3.,5.,10.};    
    array<float,8> meanPion = {-0.4669, -0.651889, -0.731293, -0.750849, -0.767918, -0.769948, -0.729383, -0.7741};
    array<float,8> meanKaon = {-0.420412, -0.656824, -0.728482, -0.99377, -1.11258, -1.04111, -1.05214, -0.778762};
    array<float,8> meanProton = {-0.346431, -0.445263, -0.504456, -0.80259, -0.971442, -1.00859, -0.853291, -0.595747};
    array<float,8> sigmaPion = {1.31686, 1.24606, 1.21786, 1.21274, 1.21565, 1.29167, 1.26293, 1.27201};
    array<float,8> sigmaKaon = {1.1904, 1.27156, 1.27005, 1.15127, 1.09914, 1.12193, 1.07542, 1.27068};
    array<float,8> sigmaProton = {1.10662, 1.18216, 1.25083, 1.3166, 1.25666, 1.12755, 1.12149, 1.20881};
    std::copy(pTPClims.begin(),pTPClims.end(),fPlimitsNsigmaTPCDataCorr);
    std::copy(meanPion.begin(),meanPion.end(),fMeanNsigmaTPCPionData);
    std::copy(meanKaon.begin(),meanKaon.end(),fMeanNsigmaTPCKaonData);
    std::copy(meanProton.begin(),meanProton.end(),fMeanNsigmaTPCProtonData);
    std::copy(sigmaPion.begin(),sigmaPion.end(),fSigmaNsigmaTPCPionData);
    std::copy(sigmaKaon.begin(),sigmaKaon.end(),fSigmaNsigmaTPCKaonData);
    std::copy(sigmaProton.begin(),sigmaProton.end(),fSigmaNsigmaTPCProtonData);
  }
  else if(fRunNumber>296690 && fRunNumber<297595 && fSystNsigmaTPCDataCorr==kPbPb3050) { //LHC18r 30-50%
    fNPbinsNsigmaTPCDataCorr = 8;
    array<float,9> pTPClims = {0.3,0.5,0.75,1.,1.5,2.,3.,5.,10.};    
    array<float,8> meanPion = {-0.298388, -0.342664, -0.396873, -0.462451, -0.538392, -0.58011, -0.612805, -0.642393};
    array<float,8> meanKaon = {-0.152114, -0.00548698, -0.425063, -0.546696, -0.530307, -0.599359, -0.605992, -0.774291};
    array<float,8> meanProton = {-0.0368995, -0.0849347, -0.30911, -0.479177, -0.584153, -0.578737, -0.568087, -0.489077};
    array<float,8> sigmaPion = {0.878966, 0.88176, 0.883022, 0.873402, 0.859262, 0.847664, 0.828672, 0.800194};
    array<float,8> sigmaKaon = {0.775762, 1.00624, 0.886505, 0.885096, 0.91232, 0.867032, 0.833827, 0.987321};
    array<float,8> sigmaProton = {0.751132, 0.788041, 0.81738, 0.862206, 0.876342, 0.863609, 0.825849, 0.937087};
    std::copy(pTPClims.begin(),pTPClims.end(),fPlimitsNsigmaTPCDataCorr);
    std::copy(meanPion.begin(),meanPion.end(),fMeanNsigmaTPCPionData);
    std::copy(meanKaon.begin(),meanKaon.end(),fMeanNsigmaTPCKaonData);
    std::copy(meanProton.begin(),meanProton.end(),fMeanNsigmaTPCProtonData);
    std::copy(sigmaPion.begin(),sigmaPion.end(),fSigmaNsigmaTPCPionData);
    std::copy(sigmaKaon.begin(),sigmaKaon.end(),fSigmaNsigmaTPCKaonData);
    std::copy(sigmaProton.begin(),sigmaProton.end(),fSigmaNsigmaTPCProtonData);
  }
  else { //default: no correction applied
    fNPbinsNsigmaTPCDataCorr = 1;
    array<float,2> pTPClims = {0.,1000.};    
    array<float,1> meanPion = {0.};
    array<float,1> meanKaon = {0.};
    array<float,1> meanProton = {0.};
    array<float,1> sigmaPion = {1.};
    array<float,1> sigmaKaon = {1.};
    array<float,1> sigmaProton = {1.};
    std::copy(pTPClims.begin(),pTPClims.end(),fPlimitsNsigmaTPCDataCorr);
    std::copy(meanPion.begin(),meanPion.end(),fMeanNsigmaTPCPionData);
    std::copy(meanKaon.begin(),meanKaon.end(),fMeanNsigmaTPCKaonData);
    std::copy(meanProton.begin(),meanProton.end(),fMeanNsigmaTPCProtonData);
    std::copy(sigmaPion.begin(),sigmaPion.end(),fSigmaNsigmaTPCPionData);
    std::copy(sigmaKaon.begin(),sigmaKaon.end(),fSigmaNsigmaTPCKaonData);
    std::copy(sigmaProton.begin(),sigmaProton.end(),fSigmaNsigmaTPCProtonData);
  }
}

//________________________________________________________________
void AliHFTreeHandler::GetNsigmaTPCMeanSigmaData(float &mean, float &sigma, AliPID::EParticleType species, float pTPC) {
    
  if(fRunNumber!=fRunNumberPrevCand)
    SetNsigmaTPCDataCorr();

  int bin = TMath::BinarySearch(fNPbinsNsigmaTPCDataCorr,fPlimitsNsigmaTPCDataCorr,pTPC);
  if(bin<0) bin=0; //underflow --> equal to min value
  else if(bin>fNPbinsNsigmaTPCDataCorr-1) bin=fNPbinsNsigmaTPCDataCorr-1; //overflow --> equal to max value

  switch(species) {
    case AliPID::kPion: 
    {
      mean = fMeanNsigmaTPCPionData[bin];
      sigma = fSigmaNsigmaTPCPionData[bin];
      break;
    }
    case AliPID::kKaon: 
    {
      mean = fMeanNsigmaTPCKaonData[bin];
      sigma = fSigmaNsigmaTPCKaonData[bin];
      break;
    }
    case AliPID::kProton: 
    {
      mean = fMeanNsigmaTPCProtonData[bin];
      sigma = fSigmaNsigmaTPCProtonData[bin];
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
