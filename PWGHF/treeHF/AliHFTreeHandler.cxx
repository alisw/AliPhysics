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
// N. Zardoshti, nima.zardoshti@cern.ch
/////////////////////////////////////////////////////////////

#include <cmath>
#include <limits>
#include "AliHFTreeHandler.h"
#include "AliPID.h"
#include "AliAODRecoDecayHF.h"
#include "AliPIDResponse.h"
#include "AliESDtrack.h"
#include "TMath.h"

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
  fPtGen(-9999.),
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
  fEvIDExt(9999),
  fEvIDLong(9999),
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
  fNEtabinsNsigmaTPCDataCorr(0),
  fPtJet(-99.),
  fPtGenJet(-99.),
  fEtaJet(-99.),
  fEtaGenJet(-99.),
  fPhiJet(-99.),
  fPhiGenJet(-99.),
  fLeadingPtJet(-99.),
  fLeadingPtGenJet(-99.),
  fDeltaEtaJetHadron(-99.),
  fDeltaEtaGenJetHadron(-99.),
  fDeltaPhiJetHadron(-99.),
  fDeltaPhiGenJetHadron(-99.),
  fDeltaRJetHadron(-99.),
  fDeltaRGenJetHadron(-99.),
  fNTracksJet(-99.),
  fNTracksGenJet(-99.),
  fZgJet(-99.),
  fZgGenJet(-99.),
  fRgJet(-99.),
  fRgGenJet(-99.),
  fNsdJet(-99.),
  fNsdGenJet(-99.),
  fPt_splittingJet(-99.),
  fPt_splittingGenJet(-99.),
  fk0Jet(-99.),
  fk0GenJet(-99.),
  fZk0Jet(-99.),
  fZk0GenJet(-99.),
  fRk0Jet(-99.),
  fRk0GenJet(-99.),
  fk1Jet(-99.),
  fk1GenJet(-99.),
  fZk1Jet(-99.),
  fZk1GenJet(-99.),
  fRk1Jet(-99.),
  fRk1GenJet(-99.),
  fk2Jet(-99.),
  fk2GenJet(-99.),
  fZk2Jet(-99.),
  fZk2GenJet(-99.),
  fRk2Jet(-99.),
  fRk2GenJet(-99.),
  fkTJet(-99.),
  fkTGenJet(-99.),
  fZkTJet(-99.),
  fZkTGenJet(-99.),
  fRkTJet(-99.),
  fRkTGenJet(-99.),
  fFillJets(false),
  fDoJetSubstructure(false), 
  fJetRadius(0.4),
  fSubJetRadius(0.2),
  fJetAlgorithm(0),
  fSubJetAlgorithm(2),
  fMinJetPt(0.0),
  fSoftDropZCut(0.1),
  fSoftDropBeta(0.0),
  fTrackingEfficiency(1.0)
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
AliHFTreeHandler::AliHFTreeHandler(int PIDopt):
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
  fEvIDExt(9999),
  fEvIDLong(9999),
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
  fNEtabinsNsigmaTPCDataCorr(0),
  fPtJet(-99.),
  fPtGenJet(-99.),
  fEtaJet(-99.),
  fEtaGenJet(-99.),
  fPhiJet(-99.),
  fPhiGenJet(-99.),
  fLeadingPtJet(-99.),
  fLeadingPtGenJet(-99.),
  fDeltaEtaJetHadron(-99.),
  fDeltaEtaGenJetHadron(-99.),
  fDeltaPhiJetHadron(-99.),
  fDeltaPhiGenJetHadron(-99.),
  fDeltaRJetHadron(-99.),
  fDeltaRGenJetHadron(-99.),
  fNTracksJet(-99.),
  fNTracksGenJet(-99.),
  fZgJet(-99.),
  fZgGenJet(-99.),
  fRgJet(-99.),
  fRgGenJet(-99.),
  fNsdJet(-99.),
  fNsdGenJet(-99.),
  fPt_splittingJet(-99.),
  fPt_splittingGenJet(-99.),
  fk0Jet(-99.),
  fk0GenJet(-99.),
  fZk0Jet(-99.),
  fZk0GenJet(-99.),
  fRk0Jet(-99.),
  fRk0GenJet(-99.),
  fk1Jet(-99.),
  fk1GenJet(-99.),
  fZk1Jet(-99.),
  fZk1GenJet(-99.),
  fRk1Jet(-99.),
  fRk1GenJet(-99.),
  fk2Jet(-99.),
  fk2GenJet(-99.),
  fZk2Jet(-99.),
  fZk2GenJet(-99.),
  fRk2Jet(-99.),
  fRk2GenJet(-99.),
  fkTJet(-99.),
  fkTGenJet(-99.),
  fZkTJet(-99.),
  fZkTGenJet(-99.),
  fRkTJet(-99.),
  fRkTGenJet(-99.),
  fFillJets(false),
  fDoJetSubstructure(false), 
  fJetRadius(0.4),
  fSubJetRadius(0.2),
  fJetAlgorithm(0),
  fSubJetAlgorithm(2),
  fMinJetPt(0.0),
  fSoftDropZCut(0.1),
  fSoftDropBeta(0.0),
  fTrackingEfficiency(1.0)
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
  fTreeVar->Branch("ev_id_ext",&fEvIDExt);
  fTreeVar->Branch("ev_id_long",&fEvIDLong);
  fTreeVar->Branch("cand_type",&fCandType);
  fTreeVar->Branch("pt_cand",&fPt);
  fTreeVar->Branch("y_cand",&fY);
  fTreeVar->Branch("eta_cand",&fEta);
  fTreeVar->Branch("phi_cand",&fPhi);
  fTreeVar->Branch("dau_in_acc",&fDauInAcceptance);

  if (fFillJets) AddGenJetBranches(); //Gen Jet Branches added here
  
  return fTreeVar;
}

//________________________________________________________________
bool AliHFTreeHandler::SetMCGenVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, AliAODMCParticle* mcpart) {

  fRunNumber = runnumber;
  fEvID = eventID;
  fEvIDExt = eventID_Ext;
  fEvIDLong = eventID_Long;
  
  if(!mcpart) return false;
  if(!(fCandType&kSignal)) return true; // fill only signal in the generated


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
void AliHFTreeHandler::AddCommonDmesonVarBranches(Bool_t HasSecVtx) {

  fTreeVar->Branch("run_number",&fRunNumber);
  fTreeVar->Branch("ev_id",&fEvID);
  fTreeVar->Branch("ev_id_ext",&fEvIDExt);
  fTreeVar->Branch("ev_id_long",&fEvIDLong);
  fTreeVar->Branch("cand_type",&fCandType);
  fTreeVar->Branch("inv_mass",&fInvMass);
  fTreeVar->Branch("pt_cand",&fPt);
  fTreeVar->Branch("pt_gen_cand",&fPtGen);
  fTreeVar->Branch("y_cand",&fY);
  fTreeVar->Branch("eta_cand",&fEta);
  fTreeVar->Branch("phi_cand",&fPhi);
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
void AliHFTreeHandler::AddSingleTrackBranches() {

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
void AliHFTreeHandler::AddJetBranches() { //Jet branches added

  fTreeVar->Branch("pt_jet",&fPtJet);
  fTreeVar->Branch("pt_gen_jet",&fPtGenJet);
  fTreeVar->Branch("eta_jet",&fEtaJet);
  fTreeVar->Branch("eta_gen_jet",&fEtaGenJet);
  fTreeVar->Branch("phi_jet",&fPhiJet);
  fTreeVar->Branch("phi_gen_jet",&fPhiGenJet);
  fTreeVar->Branch("ptleadingtrack_jet",&fLeadingPtJet);
  fTreeVar->Branch("ptleadingtrack_genjet",&fLeadingPtGenJet);
  fTreeVar->Branch("delta_eta_jet",&fDeltaEtaJetHadron);
  fTreeVar->Branch("delta_eta_gen_jet",&fDeltaEtaGenJetHadron);
  fTreeVar->Branch("delta_phi_jet",&fDeltaPhiJetHadron);
  fTreeVar->Branch("delta_phi_gen_jet",&fDeltaPhiGenJetHadron);
  fTreeVar->Branch("delta_r_jet",&fDeltaRJetHadron);
  fTreeVar->Branch("delta_r_gen_jet",&fDeltaRGenJetHadron);
  fTreeVar->Branch("ntracks_jet",&fNTracksJet);
  fTreeVar->Branch("ntracks_gen_jet",&fNTracksGenJet);
  fTreeVar->Branch("zg_jet",&fZgJet);
  fTreeVar->Branch("zg_gen_jet",&fZgGenJet);
  fTreeVar->Branch("rg_jet",&fRgJet);
  fTreeVar->Branch("rg_gen_jet",&fRgGenJet);
  fTreeVar->Branch("nsd_jet",&fNsdJet);
  fTreeVar->Branch("nsd_gen_jet",&fNsdGenJet);
  fTreeVar->Branch("Pt_splitting_jet",&fPt_splittingJet);
  fTreeVar->Branch("Pt_splitting_gen_jet",&fPt_splittingGenJet);
  fTreeVar->Branch("k0_jet",&fk0Jet);
  fTreeVar->Branch("k0_gen_jet",&fk0GenJet);
  fTreeVar->Branch("zk0_jet",&fZk0Jet);
  fTreeVar->Branch("zk0_gen_jet",&fZk0GenJet);
  fTreeVar->Branch("rk0_jet",&fRk0Jet);
  fTreeVar->Branch("rk0_gen_jet",&fRk0GenJet);
  fTreeVar->Branch("k1_jet",&fk1Jet);
  fTreeVar->Branch("k1_gen_jet",&fk1GenJet);
  fTreeVar->Branch("zk1_jet",&fZk1Jet);
  fTreeVar->Branch("zk1_gen_jet",&fZk1GenJet);
  fTreeVar->Branch("rk1_jet",&fRk1Jet);
  fTreeVar->Branch("rk1_gen_jet",&fRk1GenJet);
  fTreeVar->Branch("k2_jet",&fk2Jet);
  fTreeVar->Branch("k2_gen_jet",&fk2GenJet);
  fTreeVar->Branch("zk2_jet",&fZk2Jet);
  fTreeVar->Branch("zk2_gen_jet",&fZk2GenJet);
  fTreeVar->Branch("rk2_jet",&fRk2Jet);
  fTreeVar->Branch("rk2_gen_jet",&fRk2GenJet);
  fTreeVar->Branch("kT_jet",&fkTJet);
  fTreeVar->Branch("kT_gen_jet",&fkTGenJet);
  fTreeVar->Branch("zkT_jet",&fZkTJet);
  fTreeVar->Branch("zkT_gen_jet",&fZkTGenJet);
  fTreeVar->Branch("rkT_jet",&fRkTJet);
  fTreeVar->Branch("rkT_gen_jet",&fRkTGenJet);
  

    
}

//________________________________________________________________
void AliHFTreeHandler::AddGenJetBranches() { //Gen jet branches added

  fTreeVar->Branch("pt_jet",&fPtGenJet);
  fTreeVar->Branch("eta_jet",&fEtaGenJet);
  fTreeVar->Branch("phi_jet",&fPhiGenJet);
  fTreeVar->Branch("ptleadingtrack_jet",&fLeadingPtGenJet);
  fTreeVar->Branch("delta_eta_jet",&fDeltaEtaGenJetHadron);
  fTreeVar->Branch("delta_phi_jet",&fDeltaPhiGenJetHadron);
  fTreeVar->Branch("delta_r_jet",&fDeltaRGenJetHadron);
  fTreeVar->Branch("ntracks_jet",&fNTracksGenJet);
  fTreeVar->Branch("zg_jet",&fZgGenJet);
  fTreeVar->Branch("rg_jet",&fRgGenJet);
  fTreeVar->Branch("nsd_jet",&fNsdGenJet);
  fTreeVar->Branch("Pt_splitting_jet",&fPt_splittingGenJet);
  fTreeVar->Branch("k0_jet",&fk0GenJet);
  fTreeVar->Branch("zk0_jet",&fZk0GenJet);
  fTreeVar->Branch("rk0_jet",&fRk0GenJet);
  fTreeVar->Branch("k1_jet",&fk1GenJet);
  fTreeVar->Branch("zk1_jet",&fZk1GenJet);
  fTreeVar->Branch("rk1_jet",&fRk1GenJet);
  fTreeVar->Branch("k2_jet",&fk2GenJet);
  fTreeVar->Branch("zk2_jet",&fZk2GenJet);
  fTreeVar->Branch("rk2_jet",&fRk2GenJet);
  fTreeVar->Branch("kT_jet",&fkTGenJet);
  fTreeVar->Branch("zkT_jet",&fZkTGenJet);
  fTreeVar->Branch("rkT_jet",&fRkTGenJet);

    
}

//________________________________________________________________
void AliHFTreeHandler::AddPidBranches(bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF) 
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
void AliHFTreeHandler::SetJetVars(TClonesArray *array, AliAODRecoDecayHF* cand, Double_t invmass, TClonesArray *mcarray, AliAODMCParticle* mcPart) {
#ifdef HAVE_FASTJET
  AliHFJetFinder hfjetfinder;
  SetJetParameters(hfjetfinder); 
  AliHFJet hfjet(hfjetfinder.GetHFJet(array,cand,invmass));
  SetJetTreeVars(hfjet);
  
  AliHFJet hfgenjet;
  if (mcarray && mcPart){
    AliHFJetFinder hfgenjetfinder;
    SetJetParameters(hfgenjetfinder); 
    hfgenjet=hfgenjetfinder.GetHFMCJet(mcarray,mcPart);
  }
  SetGenJetTreeVars(hfgenjet);
#else
  std::cout << "You need to have fastjet installed to get meaningful results" <<std::endl;
#endif 

}

//________________________________________________________________
void AliHFTreeHandler::SetAndFillInclusiveJetVars(TClonesArray *array, TClonesArray *mcarray) {
#ifdef HAVE_FASTJET
  AliHFJetFinder hfjetfinder;
  SetJetParameters(hfjetfinder); 
  std::vector<AliHFJet> jets(hfjetfinder.GetJets(array));

  std::vector<AliHFJet> genjets;
  if (mcarray){
    AliHFJetFinder hfgenjetfinder;
    SetJetParameters(hfgenjetfinder); 
    genjets=hfgenjetfinder.GetMCJets(mcarray);
  }
    
  AliHFJet jet;
  AliHFJet genjet;
  for (Int_t i=0; i<jets.size(); i++){
    jet=jets[i];
    SetJetTreeVars(jet);
    if (mcarray){
      for (Int_t j=0; j<genjets.size(); j++){
	genjet=genjets[j];
	if (TMath::Sqrt(((jet.fEta-genjet.fEta)*(jet.fEta-genjet.fEta))+((jet.fPhi-genjet.fPhi)*(jet.fPhi-genjet.fPhi))) < 0.2) break;
      }
    }
    SetGenJetTreeVars(genjet);
    AliHFTreeHandler::FillTree();
  }
#else
  std::cout << "You need to have fastjet installed to get meaningful results" <<std::endl;
#endif 

}


//________________________________________________________________
void AliHFTreeHandler::SetGenJetVars(TClonesArray *array, AliAODMCParticle* mcPart) {
#ifdef HAVE_FASTJET
  AliHFJetFinder hfjetfinder;
  SetJetParameters(hfjetfinder);
  AliHFJet hfjet(hfjetfinder.GetHFMCJet(array,mcPart));
	  
  SetGenJetTreeVars(hfjet);
#else
  std::cout << "You need to have fastjet installed to get meaningful results" <<std::endl;
#endif 
}


//________________________________________________________________
void AliHFTreeHandler::SetAndFillInclusiveGenJetVars(TClonesArray *array) {
#ifdef HAVE_FASTJET
  AliHFJetFinder hfjetfinder;
  SetJetParameters(hfjetfinder);
  std::vector<AliHFJet> jets(hfjetfinder.GetMCJets(array));
  AliHFJet jet;
  for (Int_t i=0; i<jets.size(); i++){
    jet=jets[i];
    SetGenJetTreeVars(jet);
    AliHFTreeHandler::FillTree();
  }
#else
  std::cout << "You need to have fastjet installed to get meaningful results" <<std::endl;
#endif 
}
#ifdef HAVE_FASTJET

//________________________________________________________________
void AliHFTreeHandler::SetJetParameters(AliHFJetFinder& hfjetfinder){

  hfjetfinder.SetJetRadius(fJetRadius);
  hfjetfinder.SetJetAlgorithm(fJetAlgorithm);
  hfjetfinder.SetMinJetPt(fMinJetPt);
  hfjetfinder.SetSubJetRadius(fSubJetRadius);
  hfjetfinder.SetSubJetAlgorithm(fSubJetAlgorithm);
  hfjetfinder.SetDoJetSubstructure(fDoJetSubstructure);
  hfjetfinder.SetSoftDropParams(fSoftDropZCut,fSoftDropBeta);
  hfjetfinder.SetTrackingEfficiency(fTrackingEfficiency);

}
#endif
//________________________________________________________________
void AliHFTreeHandler::SetJetTreeVars(AliHFJet hfjet){

  fPtJet=hfjet.GetPt();
  fEtaJet=hfjet.GetEta();
  fPhiJet=hfjet.GetPhi();
  fLeadingPtJet=hfjet.GetLeadingPt();
  fDeltaEtaJetHadron=hfjet.GetDeltaEta();
  fDeltaPhiJetHadron=hfjet.GetDeltaPhi();
  fDeltaRJetHadron=hfjet.GetDeltaR();
  fNTracksJet=hfjet.GetN();
  if (fDoJetSubstructure){
    fZgJet=hfjet.GetZg();
    fRgJet=hfjet.GetRg();
    fNsdJet=hfjet.GetNsd();
    fPt_splittingJet=hfjet.GetPt_splitting();
    fk0Jet=hfjet.Getk0();
    fZk0Jet=hfjet.GetZk0();
    fRk0Jet=hfjet.GetRk0();
    fk1Jet=hfjet.Getk1();
    fZk1Jet=hfjet.GetZk1();
    fRk1Jet=hfjet.GetRk1();
    fk2Jet=hfjet.Getk2();
    fZk2Jet=hfjet.GetZk2();
    fRk2Jet=hfjet.GetRk2();
    fkTJet=hfjet.GetkT();
    fZkTJet=hfjet.GetZkT();
    fRkTJet=hfjet.GetRkT();
  }

}

//________________________________________________________________
void AliHFTreeHandler::SetGenJetTreeVars(AliHFJet hfjet){

  fPtGenJet=hfjet.GetPt();
  fEtaGenJet=hfjet.GetEta();
  fPhiGenJet=hfjet.GetPhi();
  fLeadingPtGenJet=hfjet.GetLeadingPt();
  fDeltaEtaGenJetHadron=hfjet.GetDeltaEta();
  fDeltaPhiGenJetHadron=hfjet.GetDeltaPhi();
  fDeltaRGenJetHadron=hfjet.GetDeltaR();
  fNTracksGenJet=hfjet.GetN();
  if (fDoJetSubstructure){
    fZgGenJet=hfjet.GetZg();
    fRgGenJet=hfjet.GetRg();
    fNsdGenJet=hfjet.GetNsd();
    fPt_splittingGenJet=hfjet.GetPt_splitting();
    fk0GenJet=hfjet.Getk0();
    fZk0GenJet=hfjet.GetZk0();
    fRk0GenJet=hfjet.GetRk0();
    fk1GenJet=hfjet.Getk1();
    fZk1GenJet=hfjet.GetZk1();
    fRk1GenJet=hfjet.GetRk1();
    fk2GenJet=hfjet.Getk2();
    fZk2GenJet=hfjet.GetZk2();
    fRk2GenJet=hfjet.GetRk2();
    fkTGenJet=hfjet.GetkT();
    fZkTGenJet=hfjet.GetZkT();
    fRkTGenJet=hfjet.GetRkT();
  }

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
double AliHFTreeHandler::CombineNsigmaDiffDet(double nsigmaTPC, double nsigmaTOF)
{
  if(nsigmaTPC > -998. && nsigmaTOF > -998.) return TMath::Sqrt((nsigmaTPC*nsigmaTPC+nsigmaTOF*nsigmaTOF)/2);
  else if(nsigmaTPC > -998. && nsigmaTOF < -998.) return TMath::Abs(nsigmaTPC);
  else if(nsigmaTPC < -998. && nsigmaTOF > -998.) return TMath::Abs(nsigmaTOF);
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

//________________________________________________________________
void AliHFTreeHandler::GetNsigmaTPCMeanSigmaData(float &mean, float &sigma, AliPID::EParticleType species, float pTPC, float eta) {
    
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
