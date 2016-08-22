/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/* AliAnalysisTaskPSHFE.h
 *
 * 
 */
#ifndef ALIANALYSISTASKEID_H
#define ALIANALYSISTASKEID_H

class TH1F;
class TList;
class AliESDtrackCuts;
class AliPIDResponse;
class AliESDEvent;
class AliESDtrack;
class TGraph;

#ifndef ALIANALYSISTASKSE_H
#include <AliAnalysisTaskSE.h>
#endif

class AliAnalysisTaskPSHFE : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskPSHFE();
    AliAnalysisTaskPSHFE(const char *name);
    virtual ~AliAnalysisTaskPSHFE();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     FillRegionHistos(AliESDEvent *esd, Int_t *elecIDs, Int_t elecCnt);
    virtual void     FillPIDHistos(AliESDEvent *esd, AliESDtrack *esdtrack, AliPIDResponse *fPIDResponse);
    virtual void     FillDPhiHistos(AliESDEvent *esd, AliESDtrack *esdtrack, Int_t i);
    void             FillPhotoElecHistos(AliESDEvent *esd, AliESDtrack *esdtrack, AliPIDResponse *fPIDResponse, Int_t i);
    void             SetTrackCuts(AliESDtrackCuts *gtrkCuts, AliESDtrackCuts *ctrkCuts);
    void             SetElectronTrackCuts(Bool_t trkCutBool);
    virtual void     Terminate(Option_t *);
    
 private:
    TList           *fOutputMB;        //! 
    TList           *fOutputEMC7; //!
    TList           *fOutputEMC8; //!
    TList           *fOutputEMCJet; //!
    AliESDtrackCuts *fTrackCutsStrong;     //
    AliESDtrackCuts *fTrackCutsWeak;     //
    AliESDtrackCuts *globaltrackCuts; 
    AliESDtrackCuts *comptrackCuts; 
    
    //Physics selection booleans
    Bool_t          MBtrg;
    Bool_t          EMC7trg;
    Bool_t          EMC8trg;
    Bool_t          EMCJettrg;
    
    //tag bools
    Bool_t          tagStrong;
    Bool_t          tagPhot;
    
    //elec track cut bool
    Bool_t          trackCutsStrong=kTRUE;
    
    //MB Histos
    //Track cut QA histos
    TH1F            *fHistTPCNClus_MB;
    TH1F            *fHistITSNClus_MB;
    TH1F            *fHistTPCSig_MB;
    TH1F            *fHistTPCSigCut_MB;
    TH1F            *fHistImpPar_MB;
    TH1F            *fHistImpParTag_MB;
    //PID QA histos
    //Pt bins of (1-2GeV, 2-3GeV, 3-4GeV, 4-5GeV, 5-6GeV, >6GeV)
    //TPC nSigma Plots
    TH2F            *fHistTPC_TOF_MB[6];
    TH2F            *fHistTPC_EMC_MB[6];
    TH2F            *fHistTPC_TRD_MB[6];
    TH2F            *fHistTPC_TOFEMC_MB[6];
    TH2F            *fHistTPC_TOFTRD_MB[6];
    TH2F            *fHistTPC_EMCTRD_MB[6];
    TH2F            *fHistTPC_TOFEMCTRD_MB[6];
    //TOF nSigma plots
    TH2F            *fHistTOF_TPC_MB[6];
    TH2F            *fHistTOF_EMC_MB[6];
    TH2F            *fHistTOF_TRD_MB[6];
    TH2F            *fHistTOF_TPCEMC_MB[6];
    TH2F            *fHistTOF_TPCTRD_MB[6];
    TH2F            *fHistTOF_EMCTRD_MB[6];
    TH2F            *fHistTOF_TPCEMCTRD_MB[6];
    //E/P Plots
    TH1F            *fHistEMC_TPC_MB[6];
    TH1F            *fHistEMC_TOF_MB[6];
    TH1F            *fHistEMC_TRD_MB[6];
    TH1F            *fHistEMC_TPCTOF_MB[6];
    TH1F            *fHistEMC_TPCTRD_MB[6];
    TH1F            *fHistEMC_TOFTRD_MB[6];
    TH1F            *fHistEMC_TPCTOFTRD_MB[6];
    TH1F            *fHistEMC_Had_MB_1Gev;
    //TRD nSigma plots
    TH2F            *fHistTRD_TPC_MB[6];
    TH2F            *fHistTRD_TOF_MB[6];
    TH2F            *fHistTRD_EMC_MB[6];
    TH2F            *fHistTRD_TPCTOF_MB[6];
    TH2F            *fHistTRD_TPCEMC_MB[6];
    TH2F            *fHistTRD_TOFEMC_MB[6];
    TH2F            *fHistTRD_TPCTOFEMC_MB[6];
    //EMCal Shower SHape plots
    TH2F            *fHistM02_All_MB[6];
    TH2F            *fHistM20_All_MB[6];
    TH2F            *fHistM02_Elec_MB[6];
    TH2F            *fHistM20_Elec_MB[6];
    //General Event histos
    TH1F            *fHistPtSum_MB;
    TH1F            *fHistPtSumTag_MB;
    TH2F            *fHistEtaPhi_MB;
    TH2F            *fHistEtaPhiTag_MB;
    TH2F            *fHistEtaPhiTPCOnly_MB;
    TH1F            *fHistDPhi300_500_MB[3];
    TH1F            *fHistDPhi500_800_MB[3];
    TH1F            *fHistDPhi800_1_MB[3];
    TH1F            *fHistDPhi1_2_MB[3];
    TH1F            *fHistDPhi2_3_MB[3];
    TH1F            *fHistDPhi3_4_MB[3];
    TH1F            *fHistDPhi4_MB[3];
    TH1F            *fHistDPhi28_MB;
    TH2F            *fHistDPhiDEta28_MB;
    TH1F            *fHistNevents_MB;
    TH1F            *fHistInvMassElecLike_MB;
    TH1F            *fHistOpAngElecLike_MB;
    TH1F            *fHistInvMassElecUnLike_MB;
    TH1F            *fHistOpAngElecUnLike_MB;
    TH1F            *fHistPtAssoc_MB;
    TH1F            *fHistPtTag_MB;
    
    
    //EMC7 Histos
    //Track Cut QA histos
    TH1F            *fHistTPCNClus_EMC7;
    TH1F            *fHistITSNClus_EMC7;
    TH1F            *fHistTPCSig_EMC7;
    TH1F            *fHistTPCSigCut_EMC7;
    TH1F            *fHistImpPar_EMC7;
    TH1F            *fHistImpParTag_EMC7;
    //PID QA histos
    //Pt bins of (1-2GeV, 2-3GeV, 3-4GeV, 4-5GeV, 5-6GeV, >6GeV)
    //DeDx Plots
    TH2F            *fHistTPC_TOF_EMC7[6];
    TH2F            *fHistTPC_EMC_EMC7[6];
    TH2F            *fHistTPC_TRD_EMC7[6];
    TH2F            *fHistTPC_TOFEMC_EMC7[6];
    TH2F            *fHistTPC_TOFTRD_EMC7[6];
    TH2F            *fHistTPC_EMCTRD_EMC7[6];
    TH2F            *fHistTPC_TOFEMCTRD_EMC7[6];
    //TOF nSigma plots
    TH2F            *fHistTOF_TPC_EMC7[6];
    TH2F            *fHistTOF_EMC_EMC7[6];
    TH2F            *fHistTOF_TRD_EMC7[6];
    TH2F            *fHistTOF_TPCEMC_EMC7[6];
    TH2F            *fHistTOF_TPCTRD_EMC7[6];
    TH2F            *fHistTOF_EMCTRD_EMC7[6];
    TH2F            *fHistTOF_TPCEMCTRD_EMC7[6];
    //E/P Plots
    TH1F            *fHistEMC_TPC_EMC7[6];
    TH1F            *fHistEMC_TOF_EMC7[6];
    TH1F            *fHistEMC_TRD_EMC7[6];
    TH1F            *fHistEMC_TPCTOF_EMC7[6];
    TH1F            *fHistEMC_TPCTRD_EMC7[6];
    TH1F            *fHistEMC_TOFTRD_EMC7[6];
    TH1F            *fHistEMC_TPCTOFTRD_EMC7[6];
    //TRD Liklihood plots
    TH2F            *fHistTRD_TPC_EMC7[6];
    TH2F            *fHistTRD_TOF_EMC7[6];
    TH2F            *fHistTRD_EMC_EMC7[6];
    TH2F            *fHistTRD_TPCTOF_EMC7[6];
    TH2F            *fHistTRD_TPCEMC_EMC7[6];
    TH2F            *fHistTRD_TOFEMC_EMC7[6];
    TH2F            *fHistTRD_TPCTOFEMC_EMC7[6];
    //EMCal Shower SHape plots
    TH2F            *fHistM02_All_EMC7[6];
    TH2F            *fHistM20_All_EMC7[6];
    TH2F            *fHistM02_Elec_EMC7[6];
    TH2F            *fHistM20_Elec_EMC7[6];
    //General Event histos 
    TH1F            *fHistPtSum_EMC7;
    TH1F            *fHistPtSumTag_EMC7;
    TH2F            *fHistEtaPhi_EMC7;
    TH2F            *fHistEtaPhiTag_EMC7;
    TH1F            *fHistDPhi300_500_EMC7[3];
    TH1F            *fHistDPhi500_800_EMC7[3];
    TH1F            *fHistDPhi800_1_EMC7[3];
    TH1F            *fHistDPhi1_2_EMC7[3];
    TH1F            *fHistDPhi2_3_EMC7[3];
    TH1F            *fHistDPhi3_4_EMC7[3];
    TH1F            *fHistDPhi4_EMC7[3];
    TH1F            *fHistDPhi28_EMC7;
    TH2F            *fHistDPhiDEta28_EMC7;
    TH1F            *fHistNevents_EMC7;
    TH1F            *fHistInvMassElecLike_EMC7;
    TH1F            *fHistOpAngElecLike_EMC7;
    TH1F            *fHistInvMassElecUnLike_EMC7;
    TH1F            *fHistOpAngElecUnLike_EMC7;
    TH1F            *fHistPtAssoc_EMC7;
    TH1F            *fHistPtTag_EMC7;
    
    
    //EMC8 Histos
    //Track Cut histos
    TH1F            *fHistTPCNClus_EMC8;
    TH1F            *fHistITSNClus_EMC8;
    TH1F            *fHistTPCSig_EMC8;
    TH1F            *fHistTPCSigCut_EMC8;
    TH1F            *fHistImpPar_EMC8;
    TH1F            *fHistImpParTag_EMC8;
    //PID QA histos
    //Pt bins of (1-2GeV, 2-3GeV, 3-4GeV, 4-5GeV, 5-6GeV, >6GeV)
    //DeDx Plots
    TH2F            *fHistTPC_TOF_EMC8[6];
    TH2F            *fHistTPC_EMC_EMC8[6];
    TH2F            *fHistTPC_TRD_EMC8[6];
    TH2F            *fHistTPC_TOFEMC_EMC8[6];
    TH2F            *fHistTPC_TOFTRD_EMC8[6];
    TH2F            *fHistTPC_EMCTRD_EMC8[6];
    TH2F            *fHistTPC_TOFEMCTRD_EMC8[6];
    //TOF nSigma plots
    TH2F            *fHistTOF_TPC_EMC8[6];
    TH2F            *fHistTOF_EMC_EMC8[6];
    TH2F            *fHistTOF_TRD_EMC8[6];
    TH2F            *fHistTOF_TPCEMC_EMC8[6];
    TH2F            *fHistTOF_TPCTRD_EMC8[6];
    TH2F            *fHistTOF_EMCTRD_EMC8[6];
    TH2F            *fHistTOF_TPCEMCTRD_EMC8[6];
    //E/P Plots
    TH1F            *fHistEMC_TPC_EMC8[6];
    TH1F            *fHistEMC_TOF_EMC8[6];
    TH1F            *fHistEMC_TRD_EMC8[6];
    TH1F            *fHistEMC_TPCTOF_EMC8[6];
    TH1F            *fHistEMC_TPCTRD_EMC8[6];
    TH1F            *fHistEMC_TOFTRD_EMC8[6];
    TH1F            *fHistEMC_TPCTOFTRD_EMC8[6];
    //TRD Liklihood plots
    TH2F            *fHistTRD_TPC_EMC8[6];
    TH2F            *fHistTRD_TOF_EMC8[6];
    TH2F            *fHistTRD_EMC_EMC8[6];
    TH2F            *fHistTRD_TPCTOF_EMC8[6];
    TH2F            *fHistTRD_TPCEMC_EMC8[6];
    TH2F            *fHistTRD_TOFEMC_EMC8[6];
    TH2F            *fHistTRD_TPCTOFEMC_EMC8[6];
    //EMCal Shower SHape plots
    TH2F            *fHistM02_All_EMC8[6];
    TH2F            *fHistM20_All_EMC8[6];
    TH2F            *fHistM02_Elec_EMC8[6];
    TH2F            *fHistM20_Elec_EMC8[6];
    //General Event histos
    TH1F            *fHistPtSum_EMC8;
    TH1F            *fHistPtSumTag_EMC8;
    TH2F            *fHistEtaPhi_EMC8;
    TH2F            *fHistEtaPhiTag_EMC8;
    TH1F            *fHistDPhi300_500_EMC8[3];
    TH1F            *fHistDPhi500_800_EMC8[3];
    TH1F            *fHistDPhi800_1_EMC8[3];
    TH1F            *fHistDPhi1_2_EMC8[3];
    TH1F            *fHistDPhi2_3_EMC8[3];
    TH1F            *fHistDPhi3_4_EMC8[3];
    TH1F            *fHistDPhi4_EMC8[3];
    TH1F            *fHistDPhi28_EMC8;
    TH2F            *fHistDPhiDEta28_EMC8;
    TH1F            *fHistNevents_EMC8;
    TH1F            *fHistInvMassElecLike_EMC8;
    TH1F            *fHistOpAngElecLike_EMC8;
    TH1F            *fHistInvMassElecUnLike_EMC8;
    TH1F            *fHistOpAngElecUnLike_EMC8;
    TH1F            *fHistPtAssoc_EMC8;
    TH1F            *fHistPtTag_EMC8;
    
    
    //EMCJet Histos
    //Track Cut histos
    TH1F            *fHistTPCNClus_EMCJet;
    TH1F            *fHistITSNClus_EMCJet;
    TH1F            *fHistTPCSig_EMCJet;
    TH1F            *fHistTPCSigCut_EMCJet;
    TH1F            *fHistImpPar_EMCJet;
    TH1F            *fHistImpParTag_EMCJet;
    //PID QA histos
    //Pt bins of (1-2GeV, 2-3GeV, 3-4GeV, 4-5GeV, 5-6GeV, >6GeV)
    //DeDx Plots
    TH2F            *fHistTPC_TOF_EMCJet[6];
    TH2F            *fHistTPC_EMC_EMCJet[6];
    TH2F            *fHistTPC_TRD_EMCJet[6];
    TH2F            *fHistTPC_TOFEMC_EMCJet[6];
    TH2F            *fHistTPC_TOFTRD_EMCJet[6];
    TH2F            *fHistTPC_EMCTRD_EMCJet[6];
    TH2F            *fHistTPC_TOFEMCTRD_EMCJet[6];
    //TOF nSigma plots
    TH2F            *fHistTOF_TPC_EMCJet[6];
    TH2F            *fHistTOF_EMC_EMCJet[6];
    TH2F            *fHistTOF_TRD_EMCJet[6];
    TH2F            *fHistTOF_TPCEMC_EMCJet[6];
    TH2F            *fHistTOF_TPCTRD_EMCJet[6];
    TH2F            *fHistTOF_EMCTRD_EMCJet[6];
    TH2F            *fHistTOF_TPCEMCTRD_EMCJet[6];
    //E/P Plots
    TH1F            *fHistEMC_TPC_EMCJet[6];
    TH1F            *fHistEMC_TOF_EMCJet[6];
    TH1F            *fHistEMC_TRD_EMCJet[6];
    TH1F            *fHistEMC_TPCTOF_EMCJet[6];
    TH1F            *fHistEMC_TPCTRD_EMCJet[6];
    TH1F            *fHistEMC_TOFTRD_EMCJet[6];
    TH1F            *fHistEMC_TPCTOFTRD_EMCJet[6];
    //TRD Liklihood plots
    TH2F            *fHistTRD_TPC_EMCJet[6];
    TH2F            *fHistTRD_TOF_EMCJet[6];
    TH2F            *fHistTRD_EMC_EMCJet[6];
    TH2F            *fHistTRD_TPCTOF_EMCJet[6];
    TH2F            *fHistTRD_TPCEMC_EMCJet[6];
    TH2F            *fHistTRD_TOFEMC_EMCJet[6];
    TH2F            *fHistTRD_TPCTOFEMC_EMCJet[6];
    //EMCal Shower SHape plots
    TH2F            *fHistM02_All_EMCJet[6];
    TH2F            *fHistM20_All_EMCJet[6];
    TH2F            *fHistM02_Elec_EMCJet[6];
    TH2F            *fHistM20_Elec_EMCJet[6];
    //General Event histos
    TH1F            *fHistPtSum_EMCJet;
    TH1F            *fHistPtSumTag_EMCJet;
    TH2F            *fHistEtaPhi_EMCJet;
    TH2F            *fHistEtaPhiTag_EMCJet;
    TH1F            *fHistDPhi300_500_EMCJet[3];
    TH1F            *fHistDPhi500_800_EMCJet[3];
    TH1F            *fHistDPhi800_1_EMCJet[3];
    TH1F            *fHistDPhi1_2_EMCJet[3];
    TH1F            *fHistDPhi2_3_EMCJet[3];
    TH1F            *fHistDPhi3_4_EMCJet[3];
    TH1F            *fHistDPhi4_EMCJet[3];
    TH1F            *fHistDPhi28_EMCJet;
    TH2F            *fHistDPhiDEta28_EMCJet;
    TH1F            *fHistNevents_EMCJet;
    TH1F            *fHistInvMassElecLike_EMCJet;
    TH1F            *fHistOpAngElecLike_EMCJet;
    TH1F            *fHistInvMassElecUnLike_EMCJet;
    TH1F            *fHistOpAngElecUnLike_EMCJet;
    TH1F            *fHistPtAssoc_EMCJet;
    TH1F            *fHistPtTag_EMCJet;
    
    //Region Histos
    
    //Tag-Side Histos
    
    //Track Multiplicity Histos
    TH1F            *fHistTrkMultTag_MB[4];
    TH1F            *fHistTrkMultTag_EMC7[4];
    TH1F            *fHistTrkMultTag_EMC8[4];
    TH1F            *fHistTrkMultTag_EMCJet[4];
    
    //Pt Distribution Histos
    TH1F            *fHistTrkPtTag_MB[4];
    TH1F            *fHistTrkPtTag_EMC7[4];
    TH1F            *fHistTrkPtTag_EMC8[4];
    TH1F            *fHistTrkPtTag_EMCJet[4];
    
    //dEdx by Pt Histos
    TH2F            *fHistDeDxPtTag_MB[4];
    TH2F            *fHistDeDxPtTag_EMC7[4];
    TH2F            *fHistDeDxPtTag_EMC8[4];
    TH2F            *fHistDeDxPtTag_EMCJet[4];
    
    //Away-Side Histos
    
    //Track Multiplicity Histos
    TH1F            *fHistTrkMultAway_MB[4];
    TH1F            *fHistTrkMultAway_EMC7[4];
    TH1F            *fHistTrkMultAway_EMC8[4];
    TH1F            *fHistTrkMultAway_EMCJet[4];
    
    //Pt Distribution Histos
    TH1F            *fHistTrkPtAway_MB[4];
    TH1F            *fHistTrkPtAway_EMC7[4];
    TH1F            *fHistTrkPtAway_EMC8[4];
    TH1F            *fHistTrkPtAway_EMCJet[4];
    
    //dEdx by Pt Histos
    TH2F            *fHistDeDxPtAway_MB[4];
    TH2F            *fHistDeDxPtAway_EMC7[4];
    TH2F            *fHistDeDxPtAway_EMC8[4];
    TH2F            *fHistDeDxPtAway_EMCJet[4];
    
    //TransMax-Side Histos
    
    //Track Multiplicity Histos
    TH1F            *fHistTrkMultTransMax_MB[4];
    TH1F            *fHistTrkMultTransMax_EMC7[4];
    TH1F            *fHistTrkMultTransMax_EMC8[4];
    TH1F            *fHistTrkMultTransMax_EMCJet[4];
    
    //Pt Distribution Histos
    TH1F            *fHistTrkPtTransMax_MB[4];
    TH1F            *fHistTrkPtTransMax_EMC7[4];
    TH1F            *fHistTrkPtTransMax_EMC8[4];
    TH1F            *fHistTrkPtTransMax_EMCJet[4];
    
    //dEdx by Pt Histos
    TH2F            *fHistDeDxPtTransMax_MB[4];
    TH2F            *fHistDeDxPtTransMax_EMC7[4];
    TH2F            *fHistDeDxPtTransMax_EMC8[4];
    TH2F            *fHistDeDxPtTransMax_EMCJet[4];
    
    //TransMin-Side Histos
    
    //Track Multiplicity Histos
    TH1F            *fHistTrkMultTransMin_MB[4];
    TH1F            *fHistTrkMultTransMin_EMC7[4];
    TH1F            *fHistTrkMultTransMin_EMC8[4];
    TH1F            *fHistTrkMultTransMin_EMCJet[4];
    
    //Pt Distribution Histos
    TH1F            *fHistTrkPtTransMin_MB[4];
    TH1F            *fHistTrkPtTransMin_EMC7[4];
    TH1F            *fHistTrkPtTransMin_EMC8[4];
    TH1F            *fHistTrkPtTransMin_EMCJet[4];
    
    //dEdx by Pt Histos
    TH2F            *fHistDeDxPtTransMin_MB[4];
    TH2F            *fHistDeDxPtTransMin_EMC7[4];
    TH2F            *fHistDeDxPtTransMin_EMC8[4];
    TH2F            *fHistDeDxPtTransMin_EMCJet[4];
    
    //Rejection Histos
    //Conscious decision not to split them into trigger or pt bins. Ask klay about it.
    TH1F            *fHistPIDRejection;
    TH1F            *fHistBadEMCclusID;
    
    //Number of tagged electrons per event
    TH1F            *fHistNElecPerEvent;
    
    //Mismatched photonic electrons
    TH1F            *fHistPhotoMismatch;
    
    //PtSum Histos
    TH2F            *fHistPtSumTransMaxB2B;
    TH2F            *fHistPtSumTransMinB2B;
    TH2F            *fHistPtSumTransMaxLead;
    TH2F            *fHistPtSumTransMinLead;
    
    
    AliAnalysisTaskPSHFE(const AliAnalysisTaskPSHFE&); // not implemented
    AliAnalysisTaskPSHFE& operator=(const AliAnalysisTaskPSHFE&); // not implemented
    
    ClassDef(AliAnalysisTaskPSHFE, 1); // example of analysis
};

#endif

