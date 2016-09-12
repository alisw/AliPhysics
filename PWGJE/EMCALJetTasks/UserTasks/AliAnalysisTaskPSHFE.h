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
class AliAODEvent;
class AliESDtrack;
class AliAODTrack;
class AliEventPoolManager;
class AliEventPool;

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
    void             FillPIDHistos(AliESDEvent *esd, AliESDtrack *esdtrack, AliPIDResponse *fPIDResponse);
    void             FillDPhiHistos(AliESDEvent *esd, AliESDtrack *esdtrack, Int_t i);
    void             FillPhotoElecHistos(AliESDEvent *esd, AliESDtrack *esdtrack, AliPIDResponse *fPIDResponse, Int_t i);
    void             FillPIDHistos(AliAODEvent *aod, AliAODTrack *aodtrack, AliPIDResponse *fPIDResponse);
    void             FillDPhiHistos(AliAODEvent *esd, AliAODTrack *aodtrack, Int_t i);
    void             FillMEDPhiHistos(AliAODTrack *aodtrack);
    void             FillPhotoElecHistos(AliAODEvent *aod, AliAODTrack *aodtrack, AliPIDResponse *fPIDResponse, Int_t i);
    void             SetTrackCuts(AliESDtrackCuts *gtrkCuts, AliESDtrackCuts *ctrkCuts);
    void             SetElectronTrackCuts(Bool_t trkCutBool);
    void             SetSSCutBool(Bool_t SSCutBool);
    void             SetAODEvent(Bool_t AODq){aodEv=AODq;}
    void             SetUseNonSignalEvents(Bool_t use){UseNonSignalEvents=use;}
    TObjArray*       MakeTrkArr(AliAODEvent *aod);
    virtual void     Terminate(Option_t *);
    
 private:
    TList           *fOutputMB;        //! 
    TList           *fOutputEMC7; //!
    TList           *fOutputEMCEGA; //!
    TList           *fOutputEMCJet; //!
    AliESDtrackCuts *fTrackCutsStrong;     //!
    AliESDtrackCuts *fTrackCutsWeak;     //!
    AliESDtrackCuts *globaltrackCuts; //!
    AliESDtrackCuts *comptrackCuts; //!
    AliEventPoolManager *fPoolMan;//!
    AliEventPool    *fPool;//!
    TObjArray       *trkArr;//!
    
    //Boolean to keep track of whether we are using aod
    Bool_t          aodEv;
    Bool_t          UseNonSignalEvents;
    
    //Physics selection booleans
    Bool_t          MBtrg;//!
    Bool_t          EMC7trg;//!
    Bool_t          EMCEGAtrg;//!
    Bool_t          EMCJettrg;//!
    
    //tag bools
    Bool_t          tagStrong;//!
    Bool_t          tagPhot;//!
    
    //elec track cut bool
    Bool_t          trackCutsStrong=kFALSE;//!
    
    //Shower Shape Bool
    Bool_t          applySSCuts=kFALSE;//!
    
    //MB Histos
    //Track cut QA histos
    TH1F            *fHistTPCNClus_MB;//!
    TH1F            *fHistITSNClus_MB;//!
    TH1F            *fHistTPCSig_MB;//!
    TH1F            *fHistTPCSigCut_MB;//!
    TH1F            *fHistImpPar_MB;//!
    TH1F            *fHistImpParTag_MB;//!
    //PID QA histos
    //Pt bins of (1-2GeV, 2-3GeV, 3-4GeV, 4-5GeV, 5-6GeV, >6GeV)
    //TPC nSigma Plots
    TH2F            *fHistTPC_EMCTRD_MB[6];//!
    //E/P Plots
    TH1F            *fHistEMC_TPCTRD_MB[6];//!
    TH1F            *fHistEMC_Had_MB_1Gev;//!
    //TRD nSigma plots
    TH2F            *fHistTRD_TPCEMC_MB[6];//!
    //EMCal Shower SHape plots
    TH2F            *fHistM02_All_MB[6];//!
    TH2F            *fHistM20_All_MB[6];//!
    TH2F            *fHistM02_Elec_MB[6];//!
    TH2F            *fHistM20_Elec_MB[6];//!
    //General Event histos
    TH1F            *fHistPtSum_MB;//!
    TH1F            *fHistPtSumTag_MB;//!
    TH1F            *fHistPtSumEMC_MB;//!
    TH2F            *fHistEtaPhi_MB;//!
    TH2F            *fHistEtaPhiTag_MB;//!
    TH2F            *fHistEtaPhiTPCOnly_MB;//!
    TH1F            *fHistDPhi300_500_MB[3];//!
    TH1F            *fHistDPhi500_800_MB[3];//!
    TH1F            *fHistDPhi800_1_MB[3];//!
    TH1F            *fHistDPhi1_2_MB[3];//!
    TH1F            *fHistDPhi2_3_MB[3];//!
    TH1F            *fHistDPhi3_4_MB[3];//!
    TH1F            *fHistDPhi4_MB[3];//!
    TH1F            *fHistDPhi28_MB;//!
    TH2F            *fHistDPhiDEta28_MB;//!
    TH1F            *fHistNevents_MB;//!
    TH1F            *fHistInvMassElecLike_MB;//!
    TH1F            *fHistOpAngElecLike_MB;//!
    TH1F            *fHistInvMassElecUnLike_MB;//!
    TH1F            *fHistOpAngElecUnLike_MB;//!
    TH1F            *fHistPtAssoc_MB;//!
    TH1F            *fHistPtAssocMix_MB;//!
    TH1F            *fHistPtTag_MB;//!
    TH1F            *fHistPhotoMismatch_MB;//!
    TH2F            *fHistDPhi18Spe_MB;//!
    
    //Mixed Event DPhi histos
    TH1F            *fHistDPhiMix300_500_MB[3];//!
    TH1F            *fHistDPhiMix500_800_MB[3];//!
    TH1F            *fHistDPhiMix800_1_MB[3];//!
    TH1F            *fHistDPhiMix1_2_MB[3];//!
    TH1F            *fHistDPhiMix2_3_MB[3];//!
    TH1F            *fHistDPhiMix3_4_MB[3];//!
    TH1F            *fHistDPhiMix4_MB[3];//!
    TH1F            *fHistDPhiMix28_MB;//!
    TH2F            *fHistDPhiDEtaMix28_MB;//!
    
    
    //EMC7 Histos
    //Track Cut QA histos
    TH1F            *fHistTPCNClus_EMC7;//!
    TH1F            *fHistITSNClus_EMC7;//!
    TH1F            *fHistTPCSig_EMC7;//!
    TH1F            *fHistTPCSigCut_EMC7;//!
    TH1F            *fHistImpPar_EMC7;//!
    TH1F            *fHistImpParTag_EMC7;//!
    //PID QA histos
    //Pt bins of (1-2GeV, 2-3GeV, 3-4GeV, 4-5GeV, 5-6GeV, >6GeV)
    //DeDx Plots
    TH2F            *fHistTPC_EMCTRD_EMC7[6];//!
    //E/P Plots
    TH1F            *fHistEMC_TPCTRD_EMC7[6];//!
    //TRD Liklihood plots
    TH2F            *fHistTRD_TPCEMC_EMC7[6];//!
    //EMCal Shower SHape plots
    TH2F            *fHistM02_All_EMC7[6];//!
    TH2F            *fHistM20_All_EMC7[6];//!
    TH2F            *fHistM02_Elec_EMC7[6];//!
    TH2F            *fHistM20_Elec_EMC7[6];//!
    //General Event histos 
    TH1F            *fHistPtSum_EMC7;//!
    TH1F            *fHistPtSumTag_EMC7;//!
    TH1F            *fHistPtSumEMC_EMC7;//!
    TH2F            *fHistEtaPhi_EMC7;//!
    TH2F            *fHistEtaPhiTag_EMC7;//!
    TH1F            *fHistDPhi300_500_EMC7[3];//!
    TH1F            *fHistDPhi500_800_EMC7[3];//!
    TH1F            *fHistDPhi800_1_EMC7[3];//!
    TH1F            *fHistDPhi1_2_EMC7[3];//!
    TH1F            *fHistDPhi2_3_EMC7[3];//!
    TH1F            *fHistDPhi3_4_EMC7[3];//!
    TH1F            *fHistDPhi4_EMC7[3];//!
    TH1F            *fHistDPhi28_EMC7;//!
    TH2F            *fHistDPhiDEta28_EMC7;//!
    TH1F            *fHistNevents_EMC7;//!
    TH1F            *fHistInvMassElecLike_EMC7;//!
    TH1F            *fHistOpAngElecLike_EMC7;//!
    TH1F            *fHistInvMassElecUnLike_EMC7;//!
    TH1F            *fHistOpAngElecUnLike_EMC7;//!
    TH1F            *fHistPtAssoc_EMC7;//!
    TH1F            *fHistPtAssocMix_EMC7;//!
    TH1F            *fHistPtTag_EMC7;//!
    TH1F            *fHistPhotoMismatch_EMC7;//!
    TH2F            *fHistDPhi18Spe_EMC7;//!
    
    //Mixed Event DPhi histos
    TH1F            *fHistDPhiMix300_500_EMC7[3];//!
    TH1F            *fHistDPhiMix500_800_EMC7[3];//!
    TH1F            *fHistDPhiMix800_1_EMC7[3];//!
    TH1F            *fHistDPhiMix1_2_EMC7[3];//!
    TH1F            *fHistDPhiMix2_3_EMC7[3];//!
    TH1F            *fHistDPhiMix3_4_EMC7[3];//!
    TH1F            *fHistDPhiMix4_EMC7[3];//!
    TH1F            *fHistDPhiMix28_EMC7;//!
    TH2F            *fHistDPhiDEtaMix28_EMC7;//!
    
    //EMCEGA Histos
    //Track Cut QA histos
    TH1F            *fHistTPCNClus_EMCEGA;//!
    TH1F            *fHistITSNClus_EMCEGA;//!
    TH1F            *fHistTPCSig_EMCEGA;//!
    TH1F            *fHistTPCSigCut_EMCEGA;//!
    TH1F            *fHistImpPar_EMCEGA;//!
    TH1F            *fHistImpParTag_EMCEGA;//!
    //PID QA histos
    //Pt bins of (1-2GeV, 2-3GeV, 3-4GeV, 4-5GeV, 5-6GeV, >6GeV)
    //DeDx Plots
    TH2F            *fHistTPC_EMCTRD_EMCEGA[6];//!
    //E/P Plots
    TH1F            *fHistEMC_TPCTRD_EMCEGA[6];//!
    //TRD Liklihood plots
    TH2F            *fHistTRD_TPCEMC_EMCEGA[6];//!
    //EMCal Shower SHape plots
    TH2F            *fHistM02_All_EMCEGA[6];//!
    TH2F            *fHistM20_All_EMCEGA[6];//!
    TH2F            *fHistM02_Elec_EMCEGA[6];//!
    TH2F            *fHistM20_Elec_EMCEGA[6];//!
    //General Event histos
    TH1F            *fHistPtSum_EMCEGA;//!
    TH1F            *fHistPtSumTag_EMCEGA;//!
    TH1F            *fHistPtSumEMC_EMCEGA;//!
    TH2F            *fHistEtaPhi_EMCEGA;//!
    TH2F            *fHistEtaPhiTag_EMCEGA;//!
    TH1F            *fHistDPhi300_500_EMCEGA[3];//!
    TH1F            *fHistDPhi500_800_EMCEGA[3];//!
    TH1F            *fHistDPhi800_1_EMCEGA[3];//!
    TH1F            *fHistDPhi1_2_EMCEGA[3];//!
    TH1F            *fHistDPhi2_3_EMCEGA[3];//!
    TH1F            *fHistDPhi3_4_EMCEGA[3];//!
    TH1F            *fHistDPhi4_EMCEGA[3];//!
    TH1F            *fHistDPhi28_EMCEGA;//!
    TH2F            *fHistDPhiDEta28_EMCEGA;//!
    TH1F            *fHistNevents_EMCEGA;//!
    TH1F            *fHistInvMassElecLike_EMCEGA;//!
    TH1F            *fHistOpAngElecLike_EMCEGA;//!
    TH1F            *fHistInvMassElecUnLike_EMCEGA;//!
    TH1F            *fHistOpAngElecUnLike_EMCEGA;//!
    TH1F            *fHistPtAssoc_EMCEGA;//!
    TH1F            *fHistPtAssocMix_EMCEGA;//!
    TH1F            *fHistPtTag_EMCEGA;//!
    TH1F            *fHistPhotoMismatch_EMCEGA;//!
    TH2F            *fHistDPhi18Spe_EMCEGA;//!
    //Mixed Event DPhi histos
    TH1F            *fHistDPhiMix300_500_EMCEGA[3];//!
    TH1F            *fHistDPhiMix500_800_EMCEGA[3];//!
    TH1F            *fHistDPhiMix800_1_EMCEGA[3];//!
    TH1F            *fHistDPhiMix1_2_EMCEGA[3];//!
    TH1F            *fHistDPhiMix2_3_EMCEGA[3];//!
    TH1F            *fHistDPhiMix3_4_EMCEGA[3];//!
    TH1F            *fHistDPhiMix4_EMCEGA[3];//!
    TH1F            *fHistDPhiMix28_EMCEGA;//!
    TH2F            *fHistDPhiDEtaMix28_EMCEGA;//!
    
    //EMCJet Histos
    //Track Cut histos
    TH1F            *fHistTPCNClus_EMCJet;//!
    TH1F            *fHistITSNClus_EMCJet;//!
    TH1F            *fHistTPCSig_EMCJet;//!
    TH1F            *fHistTPCSigCut_EMCJet;//!
    TH1F            *fHistImpPar_EMCJet;//!
    TH1F            *fHistImpParTag_EMCJet;//!
    //PID QA histos
    //Pt bins of (1-2GeV, 2-3GeV, 3-4GeV, 4-5GeV, 5-6GeV, >6GeV)
    //DeDx Plots
    TH2F            *fHistTPC_EMCTRD_EMCJet[6];//!
    //E/P Plots
    TH1F            *fHistEMC_TPCTRD_EMCJet[6];//!
    //TRD Liklihood plots
    TH2F            *fHistTRD_TPCEMC_EMCJet[6];//!
    //EMCal Shower SHape plots
    TH2F            *fHistM02_All_EMCJet[6];//!
    TH2F            *fHistM20_All_EMCJet[6];//!
    TH2F            *fHistM02_Elec_EMCJet[6];//!
    TH2F            *fHistM20_Elec_EMCJet[6];//!
    //General Event histos
    TH1F            *fHistPtSum_EMCJet;//!
    TH1F            *fHistPtSumTag_EMCJet;//!
    TH1F            *fHistPtSumEMC_EMCJet;//!
    TH2F            *fHistEtaPhi_EMCJet;//!
    TH2F            *fHistEtaPhiTag_EMCJet;//!
    TH1F            *fHistDPhi300_500_EMCJet[3];//!
    TH1F            *fHistDPhi500_800_EMCJet[3];//!
    TH1F            *fHistDPhi800_1_EMCJet[3];//!
    TH1F            *fHistDPhi1_2_EMCJet[3];//!
    TH1F            *fHistDPhi2_3_EMCJet[3];//!
    TH1F            *fHistDPhi3_4_EMCJet[3];//!
    TH1F            *fHistDPhi4_EMCJet[3];//!
    TH1F            *fHistDPhi28_EMCJet;//!
    TH2F            *fHistDPhiDEta28_EMCJet;//!
    TH1F            *fHistNevents_EMCJet;//!
    TH1F            *fHistInvMassElecLike_EMCJet;//!
    TH1F            *fHistOpAngElecLike_EMCJet;//!
    TH1F            *fHistInvMassElecUnLike_EMCJet;//!
    TH1F            *fHistOpAngElecUnLike_EMCJet;//!
    TH1F            *fHistPtAssoc_EMCJet;//!
    TH1F            *fHistPtAssocMix_EMCJet;//!
    TH1F            *fHistPtTag_EMCJet;//!
    TH1F            *fHistPhotoMismatch_EMCJet;//!
    TH2F            *fHistDPhi18Spe_EMCJet;//!
    
    //Mixed Event DPhi histos
    TH1F            *fHistDPhiMix300_500_EMCJet[3];//!
    TH1F            *fHistDPhiMix500_800_EMCJet[3];//!
    TH1F            *fHistDPhiMix800_1_EMCJet[3];//!
    TH1F            *fHistDPhiMix1_2_EMCJet[3];//!
    TH1F            *fHistDPhiMix2_3_EMCJet[3];//!
    TH1F            *fHistDPhiMix3_4_EMCJet[3];//!
    TH1F            *fHistDPhiMix4_EMCJet[3];//!
    TH1F            *fHistDPhiMix28_EMCJet;//!
    TH2F            *fHistDPhiDEtaMix28_EMCJet;//!
    
    //Rejection Histos
    //Conscious decision not to split them into trigger or pt bins. Ask klay about it.
    TH1F            *fHistPIDRejection;//!
    
    //Number of tagged electrons per event
    TH1F            *fHistNElecPerEvent;//!
    
    AliAnalysisTaskPSHFE(const AliAnalysisTaskPSHFE&); // not implemented
    AliAnalysisTaskPSHFE& operator=(const AliAnalysisTaskPSHFE&); // not implemented
    
    ClassDef(AliAnalysisTaskPSHFE, 1); // example of analysis
};

#endif

