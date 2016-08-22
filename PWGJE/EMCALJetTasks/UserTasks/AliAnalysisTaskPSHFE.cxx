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

/* $Id$ */

/* AliAnalysisTaskPSHFE.cxx
 *
 * Macro for analysis of heavy flavour electrons
 * -Patrick Steffanic
 */
#include "AliAnalysisTaskPSHFE.h"

#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TMath.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TF1.h"
#include "TGraph.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliPID.h"
#include "AliTOFcluster.h"
#include "AliPIDResponse.h"
#include "AliTRDPIDResponse.h"

using std::cout;
using std::endl;


ClassImp(AliAnalysisTaskPSHFE)

//________________________________________________________________________
AliAnalysisTaskPSHFE::AliAnalysisTaskPSHFE() // All data members should be initialised here
   :AliAnalysisTaskSE(),
    fOutputMB(0),
    fOutputEMC7(0),
    fOutputEMC8(0),
    fOutputEMCJet(0),
    fTrackCutsStrong(0),
    fTrackCutsWeak(0),
    globaltrackCuts(0),
    comptrackCuts(0),
    EMC7trg(0),
    EMC8trg(0),
    EMCJettrg(0),
    MBtrg(0),
    tagStrong(0),

    fHistPIDRejection(0),
    fHistBadEMCclusID(0),
    fHistNElecPerEvent(0),
    fHistPtSumTransMaxB2B(0),
    fHistPtSumTransMinB2B(0),
    fHistPtSumTransMaxLead(0),
    fHistPtSumTransMinLead(0),
    fHistPhotoMismatch(0),

    fHistTPCNClus_MB(0),
    fHistITSNClus_MB(0),
    fHistTPCSig_MB(0),
    fHistTPCSigCut_MB(0),
    fHistImpPar_MB(0),
    fHistImpParTag_MB(0),
    fHistNevents_MB(0),
    fHistPtSum_MB(0),
    fHistPtSumTag_MB(0),
    fHistEtaPhi_MB(0),
    fHistEtaPhiTag_MB(0),
    fHistEtaPhiTPCOnly_MB(0),
    fHistDPhi28_MB(0),
    fHistDPhiDEta28_MB(0),
    fHistEMC_Had_MB_1Gev(0),
    fHistInvMassElecLike_MB(0),
    fHistInvMassElecUnLike_MB(0),
    fHistOpAngElecLike_MB(0),
    fHistOpAngElecUnLike_MB(0),
    fHistPtAssoc_MB(0),
    fHistPtTag_MB(0),

    fHistTPCNClus_EMC7(0),
    fHistITSNClus_EMC7(0),
    fHistTPCSig_EMC7(0),
    fHistTPCSigCut_EMC7(0),
    fHistImpPar_EMC7(0),
    fHistImpParTag_EMC7(0),
    fHistNevents_EMC7(0),
    fHistPtSum_EMC7(0),
    fHistPtSumTag_EMC7(0),
    fHistEtaPhi_EMC7(0),
    fHistEtaPhiTag_EMC7(0),
    fHistDPhi28_EMC7(0),
    fHistDPhiDEta28_EMC7(0),
    fHistInvMassElecLike_EMC7(0),
    fHistInvMassElecUnLike_EMC7(0),
    fHistOpAngElecLike_EMC7(0),
    fHistOpAngElecUnLike_EMC7(0),
    fHistPtAssoc_EMC7(0),
    fHistPtTag_EMC7(0),

    fHistTPCNClus_EMC8(0),
    fHistITSNClus_EMC8(0),
    fHistTPCSig_EMC8(0),
    fHistTPCSigCut_EMC8(0),
    fHistImpPar_EMC8(0),
    fHistImpParTag_EMC8(0),
    fHistNevents_EMC8(0),
    fHistPtSum_EMC8(0),
    fHistPtSumTag_EMC8(0),
    fHistEtaPhi_EMC8(0),
    fHistEtaPhiTag_EMC8(0),
    fHistDPhi28_EMC8(0),
    fHistDPhiDEta28_EMC8(0),
    fHistInvMassElecLike_EMC8(0),
    fHistInvMassElecUnLike_EMC8(0),
    fHistOpAngElecLike_EMC8(0),
    fHistOpAngElecUnLike_EMC8(0),
    fHistPtAssoc_EMC8(0),
    fHistPtTag_EMC8(0),

    fHistTPCNClus_EMCJet(0),
    fHistITSNClus_EMCJet(0),
    fHistTPCSig_EMCJet(0),
    fHistTPCSigCut_EMCJet(0),
    fHistImpPar_EMCJet(0),
    fHistImpParTag_EMCJet(0),
    fHistNevents_EMCJet(0),
    fHistPtSum_EMCJet(0),
    fHistPtSumTag_EMCJet(0),
    fHistEtaPhi_EMCJet(0),
    fHistEtaPhiTag_EMCJet(0),
    fHistDPhi28_EMCJet(0),
    fHistDPhiDEta28_EMCJet(0),
    fHistInvMassElecLike_EMCJet(0),
    fHistInvMassElecUnLike_EMCJet(0),
    fHistOpAngElecLike_EMCJet(0),
    fHistOpAngElecUnLike_EMCJet(0),
    fHistPtAssoc_EMCJet(0),
    fHistPtTag_EMCJet(0)
        
        // The last in the above list should not have a comma after it
{
        
        //Init the DPhi Plots here because they are stored in arrays
    for(Int_t i=0; i<3; i++){
        fHistDPhi300_500_MB[i]=0;
        fHistDPhi500_800_MB[i]=0;
        fHistDPhi800_1_MB[i]=0;
        fHistDPhi1_2_MB[i]=0;
        fHistDPhi2_3_MB[i]=0;
        fHistDPhi3_4_MB[i]=0;
        fHistDPhi4_MB[i]=0;
            
        fHistDPhi300_500_EMC7[i]=0;
        fHistDPhi500_800_EMC7[i]=0;
        fHistDPhi800_1_EMC7[i]=0;
        fHistDPhi1_2_EMC7[i]=0;
        fHistDPhi2_3_EMC7[i]=0;
        fHistDPhi3_4_EMC7[i]=0;
        fHistDPhi4_EMC7[i]=0;
            
        fHistDPhi300_500_EMC8[i]=0;
        fHistDPhi500_800_EMC8[i]=0;
        fHistDPhi800_1_EMC8[i]=0;
        fHistDPhi1_2_EMC8[i]=0;
        fHistDPhi2_3_EMC8[i]=0;
        fHistDPhi3_4_EMC8[i]=0;
        fHistDPhi4_EMC8[i]=0;
            
        fHistDPhi300_500_EMCJet[i]=0;
        fHistDPhi500_800_EMCJet[i]=0;
        fHistDPhi800_1_EMCJet[i]=0;
        fHistDPhi1_2_EMCJet[i]=0;
        fHistDPhi2_3_EMCJet[i]=0;
        fHistDPhi3_4_EMCJet[i]=0;
        fHistDPhi4_EMCJet[i]=0;
    }
        
        //Init PID Plots here since they are stored in Arrays
    for(int i=0;i<6;i++){
        //MB Plots
        fHistTPC_TOF_MB[i]=0;
        fHistTPC_EMC_MB[i]=0;
        fHistTPC_TRD_MB[i]=0;
        fHistTPC_TOFEMC_MB[i]=0;
        fHistTPC_TOFTRD_MB[i]=0;
        fHistTPC_EMCTRD_MB[i]=0;
        fHistTPC_TOFEMCTRD_MB[i]=0;
        
        fHistTOF_TPC_MB[i]=0;
        fHistTOF_EMC_MB[i]=0;
        fHistTOF_TRD_MB[i]=0;
        fHistTOF_TPCEMC_MB[i]=0;
        fHistTOF_TPCTRD_MB[i]=0;
        fHistTOF_EMCTRD_MB[i]=0;
        fHistTOF_TPCEMCTRD_MB[i]=0;
        
        fHistEMC_TPC_MB[i]=0;
        fHistEMC_TOF_MB[i]=0;
        fHistEMC_TRD_MB[i]=0;
        fHistEMC_TPCTOF_MB[i]=0;
        fHistEMC_TPCTRD_MB[i]=0;
        fHistEMC_TOFTRD_MB[i]=0;
        fHistEMC_TPCTOFTRD_MB[i]=0;
        
        fHistTRD_TPC_MB[i]=0;
        fHistTRD_TOF_MB[i]=0;
        fHistTRD_EMC_MB[i]=0;
        fHistTRD_TPCTOF_MB[i]=0;
        fHistTRD_TPCEMC_MB[i]=0;
        fHistTRD_TOFEMC_MB[i]=0;
        fHistTRD_TPCTOFEMC_MB[i]=0;
        
        //EMC7 Plots
        fHistTPC_TOF_EMC7[i]=0;
        fHistTPC_EMC_EMC7[i]=0;
        fHistTPC_TRD_EMC7[i]=0;
        fHistTPC_TOFEMC_EMC7[i]=0;
        fHistTPC_TOFTRD_EMC7[i]=0;
        fHistTPC_EMCTRD_EMC7[i]=0;
        fHistTPC_TOFEMCTRD_EMC7[i]=0;
        
        fHistTOF_TPC_EMC7[i]=0;
        fHistTOF_EMC_EMC7[i]=0;
        fHistTOF_TRD_EMC7[i]=0;
        fHistTOF_TPCEMC_EMC7[i]=0;
        fHistTOF_TPCTRD_EMC7[i]=0;
        fHistTOF_EMCTRD_EMC7[i]=0;
        fHistTOF_TPCEMCTRD_EMC7[i]=0;
        
        fHistEMC_TPC_EMC7[i]=0;
        fHistEMC_TOF_EMC7[i]=0;
        fHistEMC_TRD_EMC7[i]=0;
        fHistEMC_TPCTOF_EMC7[i]=0;
        fHistEMC_TPCTRD_EMC7[i]=0;
        fHistEMC_TOFTRD_EMC7[i]=0;
        fHistEMC_TPCTOFTRD_EMC7[i]=0;
        
        fHistTRD_TPC_EMC7[i]=0;
        fHistTRD_TOF_EMC7[i]=0;
        fHistTRD_EMC_EMC7[i]=0;
        fHistTRD_TPCTOF_EMC7[i]=0;
        fHistTRD_TPCEMC_EMC7[i]=0;
        fHistTRD_TOFEMC_EMC7[i]=0;
        fHistTRD_TPCTOFEMC_EMC7[i]=0;
        
        //EMC8 Plots
        fHistTPC_TOF_EMC8[i]=0;
        fHistTPC_EMC_EMC8[i]=0;
        fHistTPC_TRD_EMC8[i]=0;
        fHistTPC_TOFEMC_EMC8[i]=0;
        fHistTPC_TOFTRD_EMC8[i]=0;
        fHistTPC_EMCTRD_EMC8[i]=0;
        fHistTPC_TOFEMCTRD_EMC8[i]=0;
        
        fHistTOF_TPC_EMC8[i]=0;
        fHistTOF_EMC_EMC8[i]=0;
        fHistTOF_TRD_EMC8[i]=0;
        fHistTOF_TPCEMC_EMC8[i]=0;
        fHistTOF_TPCTRD_EMC8[i]=0;
        fHistTOF_EMCTRD_EMC8[i]=0;
        fHistTOF_TPCEMCTRD_EMC8[i]=0;
        
        fHistEMC_TPC_EMC8[i]=0;
        fHistEMC_TOF_EMC8[i]=0;
        fHistEMC_TRD_EMC8[i]=0;
        fHistEMC_TPCTOF_EMC8[i]=0;
        fHistEMC_TPCTRD_EMC8[i]=0;
        fHistEMC_TOFTRD_EMC8[i]=0;
        fHistEMC_TPCTOFTRD_EMC8[i]=0;
        
        fHistTRD_TPC_EMC8[i]=0;
        fHistTRD_TOF_EMC8[i]=0;
        fHistTRD_EMC_EMC8[i]=0;
        fHistTRD_TPCTOF_EMC8[i]=0;
        fHistTRD_TPCEMC_EMC8[i]=0;
        fHistTRD_TOFEMC_EMC8[i]=0;
        fHistTRD_TPCTOFEMC_EMC8[i]=0;
        
        //EMCJet Plots
        fHistTPC_TOF_EMCJet[i]=0;
        fHistTPC_EMC_EMCJet[i]=0;
        fHistTPC_TRD_EMCJet[i]=0;
        fHistTPC_TOFEMC_EMCJet[i]=0;
        fHistTPC_TOFTRD_EMCJet[i]=0;
        fHistTPC_EMCTRD_EMCJet[i]=0;
        fHistTPC_TOFEMCTRD_EMCJet[i]=0;
        
        fHistTOF_TPC_EMCJet[i]=0;
        fHistTOF_EMC_EMCJet[i]=0;
        fHistTOF_TRD_EMCJet[i]=0;
        fHistTOF_TPCEMC_EMCJet[i]=0;
        fHistTOF_TPCTRD_EMCJet[i]=0;
        fHistTOF_EMCTRD_EMCJet[i]=0;
        fHistTOF_TPCEMCTRD_EMCJet[i]=0;
        
        fHistEMC_TPC_EMCJet[i]=0;
        fHistEMC_TOF_EMCJet[i]=0;
        fHistEMC_TRD_EMCJet[i]=0;
        fHistEMC_TPCTOF_EMCJet[i]=0;
        fHistEMC_TPCTRD_EMCJet[i]=0;
        fHistEMC_TOFTRD_EMCJet[i]=0;
        fHistEMC_TPCTOFTRD_EMCJet[i]=0;
        
        fHistTRD_TPC_EMCJet[i]=0;
        fHistTRD_TOF_EMCJet[i]=0;
        fHistTRD_EMC_EMCJet[i]=0;
        fHistTRD_TPCTOF_EMCJet[i]=0;
        fHistTRD_TPCEMC_EMCJet[i]=0;
        fHistTRD_TOFEMC_EMCJet[i]=0;
        fHistTRD_TPCTOFEMC_EMCJet[i]=0;
        
        fHistM02_All_MB[i]=0;
        fHistM02_Elec_MB[i]=0;
        fHistM20_All_MB[i]=0;
        fHistM20_Elec_MB[i]=0;
        
        fHistM02_All_EMC7[i]=0;
        fHistM02_Elec_EMC7[i]=0;
        fHistM20_All_EMC7[i]=0;
        fHistM20_Elec_EMC7[i]=0;
        
        fHistM02_All_EMC8[i]=0;
        fHistM02_Elec_EMC8[i]=0;
        fHistM20_All_EMC8[i]=0;
        fHistM20_Elec_EMC8[i]=0;
        
        fHistM02_All_EMCJet[i]=0;
        fHistM02_Elec_EMCJet[i]=0;
        fHistM20_All_EMCJet[i]=0;
        fHistM20_Elec_EMCJet[i]=0;
    }
        
    //Region Histos
        
    for(Int_t i=0;i<4;i++){
        //Tag Side Histos
        fHistTrkMultTag_MB[i]=0;
        fHistTrkMultTag_EMC7[i]=0;
        fHistTrkMultTag_EMC8[i]=0;
        fHistTrkMultTag_EMCJet[i]=0;
            
        fHistTrkPtTag_MB[i]=0;
        fHistTrkPtTag_EMC7[i]=0;
        fHistTrkPtTag_EMC8[i]=0;
        fHistTrkPtTag_EMCJet[i]=0;
            
        fHistDeDxPtTag_MB[i]=0;
        fHistDeDxPtTag_EMC7[i]=0;
        fHistDeDxPtTag_EMC8[i]=0;
        fHistDeDxPtTag_EMCJet[i]=0;
            
        //Away Side Histos
        fHistTrkMultAway_MB[i]=0;
        fHistTrkMultAway_EMC7[i]=0;
        fHistTrkMultAway_EMC8[i]=0;
        fHistTrkMultAway_EMCJet[i]=0;
            
        fHistTrkPtAway_MB[i]=0;
        fHistTrkPtAway_EMC7[i]=0;
        fHistTrkPtAway_EMC8[i]=0;
        fHistTrkPtAway_EMCJet[i]=0;
            
        fHistDeDxPtAway_MB[i]=0;
        fHistDeDxPtAway_EMC7[i]=0;
        fHistDeDxPtAway_EMC8[i]=0;
        fHistDeDxPtAway_EMCJet[i]=0;
            
        //TransMax Side Histos
        fHistTrkMultTransMax_MB[i]=0;
        fHistTrkMultTransMax_EMC7[i]=0;
        fHistTrkMultTransMax_EMC8[i]=0;
        fHistTrkMultTransMax_EMCJet[i]=0;
            
        fHistTrkPtTransMax_MB[i]=0;
        fHistTrkPtTransMax_EMC7[i]=0;
        fHistTrkPtTransMax_EMC8[i]=0;
        fHistTrkPtTransMax_EMCJet[i]=0;
            
        fHistDeDxPtTransMax_MB[i]=0;
        fHistDeDxPtTransMax_EMC7[i]=0;
        fHistDeDxPtTransMax_EMC8[i]=0;
        fHistDeDxPtTransMax_EMCJet[i]=0;
            
        //TransMin Side Histos
        fHistTrkMultTransMin_MB[i]=0;
        fHistTrkMultTransMin_EMC7[i]=0;
        fHistTrkMultTransMin_EMC8[i]=0;
        fHistTrkMultTransMin_EMCJet[i]=0;
            
        fHistTrkPtTransMin_MB[i]=0;
        fHistTrkPtTransMin_EMC7[i]=0;
        fHistTrkPtTransMin_EMC8[i]=0;
        fHistTrkPtTransMin_EMCJet[i]=0;
           
        fHistDeDxPtTransMin_MB[i]=0;
        fHistDeDxPtTransMin_EMC7[i]=0;
        fHistDeDxPtTransMin_EMC8[i]=0;
        fHistDeDxPtTransMin_EMCJet[i]=0;
    }
    
    // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskPSHFE::AliAnalysisTaskPSHFE(const char *name) // All data members should be initialised here
   :AliAnalysisTaskSE(name),
    fOutputMB(0),
    fOutputEMC7(0),
    fOutputEMC8(0),
    fOutputEMCJet(0),
    fTrackCutsStrong(0),
    fTrackCutsWeak(0),
    globaltrackCuts(0),
    comptrackCuts(0),
    EMC7trg(0),
    EMC8trg(0),
    EMCJettrg(0),
    MBtrg(0),
    tagStrong(0),

    fHistPIDRejection(0),
    fHistBadEMCclusID(0),
    fHistNElecPerEvent(0),
    fHistPtSumTransMaxB2B(0),
    fHistPtSumTransMinB2B(0),
    fHistPtSumTransMaxLead(0),
    fHistPtSumTransMinLead(0),
    fHistPhotoMismatch(0),

    fHistTPCNClus_MB(0),
    fHistITSNClus_MB(0),
    fHistTPCSig_MB(0),
    fHistTPCSigCut_MB(0),
    fHistImpPar_MB(0),
    fHistImpParTag_MB(0),
    fHistNevents_MB(0),
    fHistPtSum_MB(0),
    fHistPtSumTag_MB(0),
    fHistEtaPhi_MB(0),
    fHistEtaPhiTag_MB(0),
    fHistEtaPhiTPCOnly_MB(0),
    fHistDPhi28_MB(0),
    fHistDPhiDEta28_MB(0),
    fHistEMC_Had_MB_1Gev(0),
    fHistInvMassElecLike_MB(0),
    fHistInvMassElecUnLike_MB(0),
    fHistOpAngElecLike_MB(0),
    fHistOpAngElecUnLike_MB(0),
    fHistPtAssoc_MB(0),
    fHistPtTag_MB(0),

    fHistTPCNClus_EMC7(0),
    fHistITSNClus_EMC7(0),
    fHistTPCSig_EMC7(0),
    fHistTPCSigCut_EMC7(0),
    fHistImpPar_EMC7(0),
    fHistImpParTag_EMC7(0),
    fHistNevents_EMC7(0),
    fHistPtSum_EMC7(0),
    fHistPtSumTag_EMC7(0),
    fHistEtaPhi_EMC7(0),
    fHistEtaPhiTag_EMC7(0),
    fHistDPhi28_EMC7(0),
    fHistDPhiDEta28_EMC7(0),
    fHistInvMassElecLike_EMC7(0),
    fHistInvMassElecUnLike_EMC7(0),
    fHistOpAngElecLike_EMC7(0),
    fHistOpAngElecUnLike_EMC7(0),
    fHistPtAssoc_EMC7(0),
    fHistPtTag_EMC7(0),

    fHistTPCNClus_EMC8(0),
    fHistITSNClus_EMC8(0),
    fHistTPCSig_EMC8(0),
    fHistTPCSigCut_EMC8(0),
    fHistImpPar_EMC8(0),
    fHistImpParTag_EMC8(0),
    fHistNevents_EMC8(0),
    fHistPtSum_EMC8(0),
    fHistPtSumTag_EMC8(0),
    fHistEtaPhi_EMC8(0),
    fHistEtaPhiTag_EMC8(0),
    fHistDPhi28_EMC8(0),
    fHistDPhiDEta28_EMC8(0),
    fHistInvMassElecLike_EMC8(0),
    fHistInvMassElecUnLike_EMC8(0),
    fHistOpAngElecLike_EMC8(0),
    fHistOpAngElecUnLike_EMC8(0),
    fHistPtAssoc_EMC8(0),
    fHistPtTag_EMC8(0),

    fHistTPCNClus_EMCJet(0),
    fHistITSNClus_EMCJet(0),
    fHistTPCSig_EMCJet(0),
    fHistTPCSigCut_EMCJet(0),
    fHistImpPar_EMCJet(0),
    fHistImpParTag_EMCJet(0),
    fHistNevents_EMCJet(0),
    fHistPtSum_EMCJet(0),
    fHistPtSumTag_EMCJet(0),
    fHistEtaPhi_EMCJet(0),
    fHistEtaPhiTag_EMCJet(0),
    fHistDPhi28_EMCJet(0),
    fHistDPhiDEta28_EMCJet(0),
    fHistInvMassElecLike_EMCJet(0),
    fHistInvMassElecUnLike_EMCJet(0),
    fHistOpAngElecLike_EMCJet(0),
    fHistOpAngElecUnLike_EMCJet(0),
    fHistPtAssoc_EMCJet(0),
    fHistPtTag_EMCJet(0)
        
        // The last in the above list should not have a comma after it
{
    // Constructor
    // Define input and output slots here (never in the dummy constructor)
    // Input slot #0 works with a TChain - it is connected to the default input container
    // Output slot #1 writes into a TH1 container
        
    for(Int_t i=0; i<3; i++){
        fHistDPhi300_500_MB[i]=0;
        fHistDPhi500_800_MB[i]=0;
        fHistDPhi800_1_MB[i]=0;
        fHistDPhi1_2_MB[i]=0;
        fHistDPhi2_3_MB[i]=0;
        fHistDPhi3_4_MB[i]=0;
        fHistDPhi4_MB[i]=0;
            
        fHistDPhi300_500_EMC7[i]=0;
        fHistDPhi500_800_EMC7[i]=0;
        fHistDPhi800_1_EMC7[i]=0;
        fHistDPhi1_2_EMC7[i]=0;
        fHistDPhi2_3_EMC7[i]=0;
        fHistDPhi3_4_EMC7[i]=0;
        fHistDPhi4_EMC7[i]=0;
            
        fHistDPhi300_500_EMC8[i]=0;
        fHistDPhi500_800_EMC8[i]=0;
        fHistDPhi800_1_EMC8[i]=0;
        fHistDPhi1_2_EMC8[i]=0;
        fHistDPhi2_3_EMC8[i]=0;
        fHistDPhi3_4_EMC8[i]=0;
        fHistDPhi4_EMC8[i]=0;
            
        fHistDPhi300_500_EMCJet[i]=0;
        fHistDPhi500_800_EMCJet[i]=0;
        fHistDPhi800_1_EMCJet[i]=0;
        fHistDPhi1_2_EMCJet[i]=0;
        fHistDPhi2_3_EMCJet[i]=0;
        fHistDPhi3_4_EMCJet[i]=0;
        fHistDPhi4_EMCJet[i]=0;
    }
        
            //Init PID Plots here since they are stored in Arrays
    for(Int_t i=0;i<6;i++){
        //MB Plots
        fHistTPC_TOF_MB[i]=0;
        fHistTPC_EMC_MB[i]=0;
        fHistTPC_TRD_MB[i]=0;
        fHistTPC_TOFEMC_MB[i]=0;
        fHistTPC_TOFTRD_MB[i]=0;
        fHistTPC_EMCTRD_MB[i]=0;
        fHistTPC_TOFEMCTRD_MB[i]=0;
        
        fHistTOF_TPC_MB[i]=0;
        fHistTOF_EMC_MB[i]=0;
        fHistTOF_TRD_MB[i]=0;
        fHistTOF_TPCEMC_MB[i]=0;
        fHistTOF_TPCTRD_MB[i]=0;
        fHistTOF_EMCTRD_MB[i]=0;
        fHistTOF_TPCEMCTRD_MB[i]=0;
        
        fHistEMC_TPC_MB[i]=0;
        fHistEMC_TOF_MB[i]=0;
        fHistEMC_TRD_MB[i]=0;
        fHistEMC_TPCTOF_MB[i]=0;
        fHistEMC_TPCTRD_MB[i]=0;
        fHistEMC_TOFTRD_MB[i]=0;
        fHistEMC_TPCTOFTRD_MB[i]=0;
        
        fHistTRD_TPC_MB[i]=0;
        fHistTRD_TOF_MB[i]=0;
        fHistTRD_EMC_MB[i]=0;
        fHistTRD_TPCTOF_MB[i]=0;
        fHistTRD_TPCEMC_MB[i]=0;
        fHistTRD_TOFEMC_MB[i]=0;
        fHistTRD_TPCTOFEMC_MB[i]=0;
        
        //EMC7 Plots
        fHistTPC_TOF_EMC7[i]=0;
        fHistTPC_EMC_EMC7[i]=0;
        fHistTPC_TRD_EMC7[i]=0;
        fHistTPC_TOFEMC_EMC7[i]=0;
        fHistTPC_TOFTRD_EMC7[i]=0;
        fHistTPC_EMCTRD_EMC7[i]=0;
        fHistTPC_TOFEMCTRD_EMC7[i]=0;
        
        fHistTOF_TPC_EMC7[i]=0;
        fHistTOF_EMC_EMC7[i]=0;
        fHistTOF_TRD_EMC7[i]=0;
        fHistTOF_TPCEMC_EMC7[i]=0;
        fHistTOF_TPCTRD_EMC7[i]=0;
        fHistTOF_EMCTRD_EMC7[i]=0;
        fHistTOF_TPCEMCTRD_EMC7[i]=0;
        
        fHistEMC_TPC_EMC7[i]=0;
        fHistEMC_TOF_EMC7[i]=0;
        fHistEMC_TRD_EMC7[i]=0;
        fHistEMC_TPCTOF_EMC7[i]=0;
        fHistEMC_TPCTRD_EMC7[i]=0;
        fHistEMC_TOFTRD_EMC7[i]=0;
        fHistEMC_TPCTOFTRD_EMC7[i]=0;
        
        fHistTRD_TPC_EMC7[i]=0;
        fHistTRD_TOF_EMC7[i]=0;
        fHistTRD_EMC_EMC7[i]=0;
        fHistTRD_TPCTOF_EMC7[i]=0;
        fHistTRD_TPCEMC_EMC7[i]=0;
        fHistTRD_TOFEMC_EMC7[i]=0;
        fHistTRD_TPCTOFEMC_EMC7[i]=0;
        
        //EMC8 Plots
        fHistTPC_TOF_EMC8[i]=0;
        fHistTPC_EMC_EMC8[i]=0;
        fHistTPC_TRD_EMC8[i]=0;
        fHistTPC_TOFEMC_EMC8[i]=0;
        fHistTPC_TOFTRD_EMC8[i]=0;
        fHistTPC_EMCTRD_EMC8[i]=0;
        fHistTPC_TOFEMCTRD_EMC8[i]=0;
        
        fHistTOF_TPC_EMC8[i]=0;
        fHistTOF_EMC_EMC8[i]=0;
        fHistTOF_TRD_EMC8[i]=0;
        fHistTOF_TPCEMC_EMC8[i]=0;
        fHistTOF_TPCTRD_EMC8[i]=0;
        fHistTOF_EMCTRD_EMC8[i]=0;
        fHistTOF_TPCEMCTRD_EMC8[i]=0;
        
        fHistEMC_TPC_EMC8[i]=0;
        fHistEMC_TOF_EMC8[i]=0;
        fHistEMC_TRD_EMC8[i]=0;
        fHistEMC_TPCTOF_EMC8[i]=0;
        fHistEMC_TPCTRD_EMC8[i]=0;
        fHistEMC_TOFTRD_EMC8[i]=0;
        fHistEMC_TPCTOFTRD_EMC8[i]=0;
        
        fHistTRD_TPC_EMC8[i]=0;
        fHistTRD_TOF_EMC8[i]=0;
        fHistTRD_EMC_EMC8[i]=0;
        fHistTRD_TPCTOF_EMC8[i]=0;
        fHistTRD_TPCEMC_EMC8[i]=0;
        fHistTRD_TOFEMC_EMC8[i]=0;
        fHistTRD_TPCTOFEMC_EMC8[i]=0;
        
        //EMCJet Plots
        fHistTPC_TOF_EMCJet[i]=0;
        fHistTPC_EMC_EMCJet[i]=0;
        fHistTPC_TRD_EMCJet[i]=0;
        fHistTPC_TOFEMC_EMCJet[i]=0;
        fHistTPC_TOFTRD_EMCJet[i]=0;
        fHistTPC_EMCTRD_EMCJet[i]=0;
        fHistTPC_TOFEMCTRD_EMCJet[i]=0;
        
        fHistTOF_TPC_EMCJet[i]=0;
        fHistTOF_EMC_EMCJet[i]=0;
        fHistTOF_TRD_EMCJet[i]=0;
        fHistTOF_TPCEMC_EMCJet[i]=0;
        fHistTOF_TPCTRD_EMCJet[i]=0;
        fHistTOF_EMCTRD_EMCJet[i]=0;
        fHistTOF_TPCEMCTRD_EMCJet[i]=0;
        
        fHistEMC_TPC_EMCJet[i]=0;
        fHistEMC_TOF_EMCJet[i]=0;
        fHistEMC_TRD_EMCJet[i]=0;
        fHistEMC_TPCTOF_EMCJet[i]=0;
        fHistEMC_TPCTRD_EMCJet[i]=0;
        fHistEMC_TOFTRD_EMCJet[i]=0;
        fHistEMC_TPCTOFTRD_EMCJet[i]=0;
        
        fHistTRD_TPC_EMCJet[i]=0;
        fHistTRD_TOF_EMCJet[i]=0;
        fHistTRD_EMC_EMCJet[i]=0;
        fHistTRD_TPCTOF_EMCJet[i]=0;
        fHistTRD_TPCEMC_EMCJet[i]=0;
        fHistTRD_TOFEMC_EMCJet[i]=0;
        fHistTRD_TPCTOFEMC_EMCJet[i]=0;
        
        fHistM02_All_MB[i]=0;
        fHistM02_Elec_MB[i]=0;
        fHistM20_All_MB[i]=0;
        fHistM20_Elec_MB[i]=0;
        
        fHistM02_All_EMC7[i]=0;
        fHistM02_Elec_EMC7[i]=0;
        fHistM20_All_EMC7[i]=0;
        fHistM20_Elec_EMC7[i]=0;
        
        fHistM02_All_EMC8[i]=0;
        fHistM02_Elec_EMC8[i]=0;
        fHistM20_All_EMC8[i]=0;
        fHistM20_Elec_EMC8[i]=0;
        
        fHistM02_All_EMCJet[i]=0;
        fHistM02_Elec_EMCJet[i]=0;
        fHistM20_All_EMCJet[i]=0;
        fHistM20_Elec_EMCJet[i]=0;
    }
        
        
    //Region Histos
   
    for(Int_t i=0;i<4;i++){
        //Tag Side Histos
        fHistTrkMultTag_MB[i]=0;
        fHistTrkMultTag_EMC7[i]=0;
        fHistTrkMultTag_EMC8[i]=0;
        fHistTrkMultTag_EMCJet[i]=0;
            
        fHistTrkPtTag_MB[i]=0;
        fHistTrkPtTag_EMC7[i]=0;
        fHistTrkPtTag_EMC8[i]=0;
        fHistTrkPtTag_EMCJet[i]=0;
            
        fHistDeDxPtTag_MB[i]=0;
        fHistDeDxPtTag_EMC7[i]=0;
        fHistDeDxPtTag_EMC8[i]=0;
        fHistDeDxPtTag_EMCJet[i]=0;
            
        //Away Side Histos
        fHistTrkMultAway_MB[i]=0;
        fHistTrkMultAway_EMC7[i]=0;
        fHistTrkMultAway_EMC8[i]=0;
        fHistTrkMultAway_EMCJet[i]=0;
            
        fHistTrkPtAway_MB[i]=0;
        fHistTrkPtAway_EMC7[i]=0;
        fHistTrkPtAway_EMC8[i]=0;
        fHistTrkPtAway_EMCJet[i]=0;
          
        fHistDeDxPtAway_MB[i]=0;
        fHistDeDxPtAway_EMC7[i]=0;
        fHistDeDxPtAway_EMC8[i]=0;
        fHistDeDxPtAway_EMCJet[i]=0;
            
        //TransMax Side Histos
        fHistTrkMultTransMax_MB[i]=0;
        fHistTrkMultTransMax_EMC7[i]=0;
        fHistTrkMultTransMax_EMC8[i]=0;
        fHistTrkMultTransMax_EMCJet[i]=0;
           
        fHistTrkPtTransMax_MB[i]=0;
        fHistTrkPtTransMax_EMC7[i]=0;
        fHistTrkPtTransMax_EMC8[i]=0;
        fHistTrkPtTransMax_EMCJet[i]=0;
            
        fHistDeDxPtTransMax_MB[i]=0;
        fHistDeDxPtTransMax_EMC7[i]=0;
        fHistDeDxPtTransMax_EMC8[i]=0;
        fHistDeDxPtTransMax_EMCJet[i]=0;
            
        //TransMin Side Histos
        fHistTrkMultTransMin_MB[i]=0;
        fHistTrkMultTransMin_EMC7[i]=0;
        fHistTrkMultTransMin_EMC8[i]=0;
        fHistTrkMultTransMin_EMCJet[i]=0;
           
        fHistTrkPtTransMin_MB[i]=0;
        fHistTrkPtTransMin_EMC7[i]=0;
        fHistTrkPtTransMin_EMC8[i]=0;
        fHistTrkPtTransMin_EMCJet[i]=0;
        
        fHistDeDxPtTransMin_MB[i]=0;
        fHistDeDxPtTransMin_EMC7[i]=0;
        fHistDeDxPtTransMin_EMC8[i]=0;
        fHistDeDxPtTransMin_EMCJet[i]=0;
    }
        
        
    DefineOutput(1, TList::Class());//MB
    DefineOutput(2, TList::Class());//EMC7
    DefineOutput(3, TList::Class());//EMC8
    DefineOutput(4, TList::Class());//EMCJet
}

//________________________________________________________________________
AliAnalysisTaskPSHFE::~AliAnalysisTaskPSHFE()
{
    // Destructor. Clean-up the output list, but not the histograms that are put inside
    // (the list is owner and will clean-up these histograms). Protect in PROOF case.
    for(Int_t i=0;i<3;i++){
        delete fHistDPhi300_500_MB[i];
        delete fHistDPhi500_800_MB[i];
        delete fHistDPhi800_1_MB[i];
        delete fHistDPhi1_2_MB[i];
        delete fHistDPhi2_3_MB[i];
        delete fHistDPhi3_4_MB[i];
        delete fHistDPhi4_MB[i];
        
        delete fHistDPhi300_500_EMC7[i];
        delete fHistDPhi500_800_EMC7[i];
        delete fHistDPhi800_1_EMC7[i];
        delete fHistDPhi1_2_EMC7[i];
        delete fHistDPhi2_3_EMC7[i];
        delete fHistDPhi3_4_EMC7[i];
        delete fHistDPhi4_EMC7[i];
        
        delete fHistDPhi300_500_EMC8[i];
        delete fHistDPhi500_800_EMC8[i];
        delete fHistDPhi800_1_EMC8[i];
        delete fHistDPhi1_2_EMC8[i];
        delete fHistDPhi2_3_EMC8[i];
        delete fHistDPhi3_4_EMC8[i];
        delete fHistDPhi4_EMC8[i];
        
        delete fHistDPhi300_500_EMCJet[i];
        delete fHistDPhi500_800_EMCJet[i];
        delete fHistDPhi800_1_EMCJet[i];
        delete fHistDPhi1_2_EMCJet[i];
        delete fHistDPhi2_3_EMCJet[i];
        delete fHistDPhi3_4_EMCJet[i];
        delete fHistDPhi4_EMCJet[i];
    }
    
        //Init PID Plots here since they are stored in Arrays
    for(Int_t i=0;i<6;i++){
        //MB Plots
        delete fHistTPC_TOF_MB[i];
        delete fHistTPC_EMC_MB[i];
        delete fHistTPC_TRD_MB[i];
        delete fHistTPC_TOFEMC_MB[i];
        delete fHistTPC_TOFTRD_MB[i];
        delete fHistTPC_EMCTRD_MB[i];
        delete fHistTPC_TOFEMCTRD_MB[i];
        
        delete fHistTOF_TPC_MB[i];
        delete fHistTOF_EMC_MB[i];
        delete fHistTOF_TRD_MB[i];
        delete fHistTOF_TPCEMC_MB[i];
        delete fHistTOF_TPCTRD_MB[i];
        delete fHistTOF_EMCTRD_MB[i];
        delete fHistTOF_TPCEMCTRD_MB[i];
        
        delete fHistEMC_TPC_MB[i];
        delete fHistEMC_TOF_MB[i];
        delete fHistEMC_TRD_MB[i];
        delete fHistEMC_TPCTOF_MB[i];
        delete fHistEMC_TPCTRD_MB[i];
        delete fHistEMC_TOFTRD_MB[i];
        delete fHistEMC_TPCTOFTRD_MB[i];
        
        delete fHistTRD_TPC_MB[i];
        delete fHistTRD_TOF_MB[i];
        delete fHistTRD_EMC_MB[i];
        delete fHistTRD_TPCTOF_MB[i];
        delete fHistTRD_TPCEMC_MB[i];
        delete fHistTRD_TOFEMC_MB[i];
        delete fHistTRD_TPCTOFEMC_MB[i];
        
        //EMC7 Plots
        delete fHistTPC_TOF_EMC7[i];
        delete fHistTPC_EMC_EMC7[i];
        delete fHistTPC_TRD_EMC7[i];
        delete fHistTPC_TOFEMC_EMC7[i];
        delete fHistTPC_TOFTRD_EMC7[i];
        delete fHistTPC_EMCTRD_EMC7[i];
        delete fHistTPC_TOFEMCTRD_EMC7[i];
        
        delete fHistTOF_TPC_EMC7[i];
        delete fHistTOF_EMC_EMC7[i];
        delete fHistTOF_TRD_EMC7[i];
        delete fHistTOF_TPCEMC_EMC7[i];
        delete fHistTOF_TPCTRD_EMC7[i];
        delete fHistTOF_EMCTRD_EMC7[i];
        delete fHistTOF_TPCEMCTRD_EMC7[i];
        
        delete fHistEMC_TPC_EMC7[i];
        delete fHistEMC_TOF_EMC7[i];
        delete fHistEMC_TRD_EMC7[i];
        delete fHistEMC_TPCTOF_EMC7[i];
        delete fHistEMC_TPCTRD_EMC7[i];
        delete fHistEMC_TOFTRD_EMC7[i];
        delete fHistEMC_TPCTOFTRD_EMC7[i];
        
        delete fHistTRD_TPC_EMC7[i];
        delete fHistTRD_TOF_EMC7[i];
        delete fHistTRD_EMC_EMC7[i];
        delete fHistTRD_TPCTOF_EMC7[i];
        delete fHistTRD_TPCEMC_EMC7[i];
        delete fHistTRD_TOFEMC_EMC7[i];
        delete fHistTRD_TPCTOFEMC_EMC7[i];
        
        //EMC8 Plots
        delete fHistTPC_TOF_EMC8[i];
        delete fHistTPC_EMC_EMC8[i];
        delete fHistTPC_TRD_EMC8[i];
        delete fHistTPC_TOFEMC_EMC8[i];
        delete fHistTPC_TOFTRD_EMC8[i];
        delete fHistTPC_EMCTRD_EMC8[i];
        delete fHistTPC_TOFEMCTRD_EMC8[i];
        
        delete fHistTOF_TPC_EMC8[i];
        delete fHistTOF_EMC_EMC8[i];
        delete fHistTOF_TRD_EMC8[i];
        delete fHistTOF_TPCEMC_EMC8[i];
        delete fHistTOF_TPCTRD_EMC8[i];
        delete fHistTOF_EMCTRD_EMC8[i];
        delete fHistTOF_TPCEMCTRD_EMC8[i];
        
        delete fHistEMC_TPC_EMC8[i];
        delete fHistEMC_TOF_EMC8[i];
        delete fHistEMC_TRD_EMC8[i];
        delete fHistEMC_TPCTOF_EMC8[i];
        delete fHistEMC_TPCTRD_EMC8[i];
        delete fHistEMC_TOFTRD_EMC8[i];
        delete fHistEMC_TPCTOFTRD_EMC8[i];
        
        delete fHistTRD_TPC_EMC8[i];
        delete fHistTRD_TOF_EMC8[i];
        delete fHistTRD_EMC_EMC8[i];
        delete fHistTRD_TPCTOF_EMC8[i];
        delete fHistTRD_TPCEMC_EMC8[i];
        delete fHistTRD_TOFEMC_EMC8[i];
        delete fHistTRD_TPCTOFEMC_EMC8[i];
        
        //EMCJet Plots
        delete fHistTPC_TOF_EMCJet[i];
        delete fHistTPC_EMC_EMCJet[i];
        delete fHistTPC_TRD_EMCJet[i];
        delete fHistTPC_TOFEMC_EMCJet[i];
        delete fHistTPC_TOFTRD_EMCJet[i];
        delete fHistTPC_EMCTRD_EMCJet[i];
        delete fHistTPC_TOFEMCTRD_EMCJet[i];
        
        delete fHistTOF_TPC_EMCJet[i];
        delete fHistTOF_EMC_EMCJet[i];
        delete fHistTOF_TRD_EMCJet[i];
        delete fHistTOF_TPCEMC_EMCJet[i];
        delete fHistTOF_TPCTRD_EMCJet[i];
        delete fHistTOF_EMCTRD_EMCJet[i];
        delete fHistTOF_TPCEMCTRD_EMCJet[i];
        
        delete fHistEMC_TPC_EMCJet[i];
        delete fHistEMC_TOF_EMCJet[i];
        delete fHistEMC_TRD_EMCJet[i];
        delete fHistEMC_TPCTOF_EMCJet[i];
        delete fHistEMC_TPCTRD_EMCJet[i];
        delete fHistEMC_TOFTRD_EMCJet[i];
        delete fHistEMC_TPCTOFTRD_EMCJet[i];
        
        delete fHistTRD_TPC_EMCJet[i];
        delete fHistTRD_TOF_EMCJet[i];
        delete fHistTRD_EMC_EMCJet[i];
        delete fHistTRD_TPCTOF_EMCJet[i];
        delete fHistTRD_TPCEMC_EMCJet[i];
        delete fHistTRD_TOFEMC_EMCJet[i];
        delete fHistTRD_TPCTOFEMC_EMCJet[i];
        
        delete fHistM02_All_MB[i];
        delete fHistM02_Elec_MB[i];
        delete fHistM20_All_MB[i];
        delete fHistM20_Elec_MB[i];
        
        delete fHistM02_All_EMC7[i];
        delete fHistM02_Elec_EMC7[i];
        delete fHistM20_All_EMC7[i];
        delete fHistM20_Elec_EMC7[i];
        
        delete fHistM02_All_EMC8[i];
        delete fHistM02_Elec_EMC8[i];
        delete fHistM20_All_EMC8[i];
        delete fHistM20_Elec_EMC8[i];
        
        delete fHistM02_All_EMCJet[i];
        delete fHistM02_Elec_EMCJet[i];
        delete fHistM20_All_EMCJet[i];
        delete fHistM20_Elec_EMCJet[i];
    }
    
    //Region Histos
    
    for(Int_t i=0;i<4;i++){
        //Tag Side Histos
        delete fHistTrkMultTag_MB[i];
        delete fHistTrkMultTag_EMC7[i];
        delete fHistTrkMultTag_EMC8[i];
        delete fHistTrkMultTag_EMCJet[i];
        
        delete fHistTrkPtTag_MB[i];
        delete fHistTrkPtTag_EMC7[i];
        delete fHistTrkPtTag_EMC8[i];
        delete fHistTrkPtTag_EMCJet[i];
          
        delete fHistDeDxPtTag_MB[i];
        delete fHistDeDxPtTag_EMC7[i];
        delete fHistDeDxPtTag_EMC8[i];
        delete fHistDeDxPtTag_EMCJet[i];
          
        //Away Side Histos
        delete fHistTrkMultAway_MB[i];
        delete fHistTrkMultAway_EMC7[i];
        delete fHistTrkMultAway_EMC8[i];
        delete fHistTrkMultAway_EMCJet[i];
        
        delete fHistTrkPtAway_MB[i];
        delete fHistTrkPtAway_EMC7[i];
        delete fHistTrkPtAway_EMC8[i];
        delete fHistTrkPtAway_EMCJet[i];
         
        delete fHistDeDxPtAway_MB[i];
        delete fHistDeDxPtAway_EMC7[i];
        delete fHistDeDxPtAway_EMC8[i];
        delete fHistDeDxPtAway_EMCJet[i];
            
        //TransMax Side Histos
        delete fHistTrkMultTransMax_MB[i];
        delete fHistTrkMultTransMax_EMC7[i];
        delete fHistTrkMultTransMax_EMC8[i];
        delete fHistTrkMultTransMax_EMCJet[i];
        
        delete fHistTrkPtTransMax_MB[i];
        delete fHistTrkPtTransMax_EMC7[i];
        delete fHistTrkPtTransMax_EMC8[i];
        delete fHistTrkPtTransMax_EMCJet[i];
        
        delete fHistDeDxPtTransMax_MB[i];
        delete fHistDeDxPtTransMax_EMC7[i];
        delete fHistDeDxPtTransMax_EMC8[i];
        delete fHistDeDxPtTransMax_EMCJet[i];
        
        //TransMin Side Histos
        delete fHistTrkMultTransMin_MB[i];
        delete fHistTrkMultTransMin_EMC7[i];
        delete fHistTrkMultTransMin_EMC8[i];
        delete fHistTrkMultTransMin_EMCJet[i];
            
        delete fHistTrkPtTransMin_MB[i];
        delete fHistTrkPtTransMin_EMC7[i];
        delete fHistTrkPtTransMin_EMC8[i];
        delete fHistTrkPtTransMin_EMCJet[i];
            
        delete fHistDeDxPtTransMin_MB[i];
        delete fHistDeDxPtTransMin_EMC7[i];
        delete fHistDeDxPtTransMin_EMC8[i];
        delete fHistDeDxPtTransMin_EMCJet[i];
    }
    
    delete fTrackCutsStrong;
    delete fTrackCutsWeak;
}

//________________________________________________________________________
void AliAnalysisTaskPSHFE::UserCreateOutputObjects(){
    // Create histograms
    // Called once (on the worker node)
        
    fOutputMB = new TList();
    OpenFile(1);
    fOutputMB->SetOwner();  // IMPORTANT!
    fOutputEMC7 = new TList();
    OpenFile(2);
    fOutputEMC7->SetOwner();  // IMPORTANT!
    fOutputEMC8 = new TList();
    OpenFile(3);
    fOutputEMC8->SetOwner();  // IMPORTANT!
    fOutputEMCJet = new TList();
    OpenFile(4);
    fOutputEMCJet->SetOwner();  // IMPORTANT!
    
    //Some strings for histograms
    TString ptRangesDPhi[3] = {"1-2Gev", "2-4Gev", "4-8Gev"};
    TString ptRangesPID[6] = {"1-2GeV", "2-3GeV", "3-4GeV", "4-5GeV", "5-6GeV", ">6GeV"};
    TString ptRangesRegion[4] = {"1-2Gev", "2-4Gev", "4-6Gev", ">6Gev"};
    
    //Strong cuts for heavy flavour
    
    fTrackCutsStrong = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
    fTrackCutsStrong->SetMinNClustersTPC(120);
    fTrackCutsStrong->SetMinNClustersITS(4);
    fTrackCutsStrong->SetPtRange(1,10e10);
    fTrackCutsStrong->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kBoth);
    fTrackCutsStrong->SetMaxChi2PerClusterTPC(2);
    fTrackCutsStrong->SetMaxDCAToVertexXY(1);
    fTrackCutsStrong->SetMaxDCAToVertexZ(2);
    fTrackCutsStrong->SetEtaRange(-.6, .6);
    fTrackCutsStrong->SetAcceptKinkDaughters(kFALSE);
    
    //Weak cuts for heavy flavour
    
    fTrackCutsWeak = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
    fTrackCutsWeak->SetMinNClustersTPC(80);
    fTrackCutsWeak->SetMinNClustersITS(3);
    fTrackCutsWeak->SetPtRange(1,10e10);
    fTrackCutsWeak->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kBoth);
    fTrackCutsWeak->SetMaxChi2PerClusterTPC(4);
    //fTrackCutsWeak->SetMaxDCAToVertexXY(2);
    //fTrackCutsWeak->SetMaxDCAToVertexZ(3);
    fTrackCutsWeak->SetEtaRange(-.7, .7);
    fTrackCutsWeak->SetAcceptKinkDaughters(kFALSE);
    
    // Create histograms
    
    //Photonic e mismatch histo
    
    fHistPhotoMismatch = new TH1F("fHistPhotoMismatch", "Electrons identified as 'heavy flavour' that fall in photonic invariant mass and opening angle cuts", 2, 0, 1);
    fHistPhotoMismatch->GetXaxis()->SetTitle("Electrons");
    fHistPhotoMismatch->GetYaxis()->SetTitle("Cts");
    
    //Invariant mass histos
    
    fHistInvMassElecLike_MB = new TH1F("fHistInvMassElecLike_MB", "Invariant mass for all like-signed electron pairs", 100, 0, .5);
    fHistInvMassElecLike_MB->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecLike_MB->GetYaxis()->SetTitle("Cts");
    
    fHistInvMassElecLike_EMC7 = new TH1F("fHistInvMassElecLike_EMC7", "Invariant mass for all like-signed electron pairs", 100, 0, .5);
    fHistInvMassElecLike_EMC7->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecLike_EMC7->GetYaxis()->SetTitle("Cts");
    
    fHistInvMassElecLike_EMC8 = new TH1F("fHistInvMassElecLike_EMC8", "Invariant mass for all like-signed electron pairs", 100, 0, .5);
    fHistInvMassElecLike_EMC8->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecLike_EMC8->GetYaxis()->SetTitle("Cts");
    
    fHistInvMassElecLike_EMCJet = new TH1F("fHistInvMassElecLike_EMCJet", "Invariant mass for all like-signed electron pairs", 100, 0, .5);
    fHistInvMassElecLike_EMCJet->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecLike_EMCJet->GetYaxis()->SetTitle("Cts");
    
    fHistInvMassElecUnLike_MB = new TH1F("fHistInvMassElecUnLike_MB", "Invariant mass for all unlike-signed electron pairs", 100, 0, .5);
    fHistInvMassElecUnLike_MB->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecUnLike_MB->GetYaxis()->SetTitle("Cts");
    
    fHistInvMassElecUnLike_EMC7 = new TH1F("fHistInvMassElecUnLike_EMC7", "Invariant mass for all unlike-signed electron pairs", 100, 0, .5);
    fHistInvMassElecUnLike_EMC7->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecUnLike_EMC7->GetYaxis()->SetTitle("Cts");
    
    fHistInvMassElecUnLike_EMC8 = new TH1F("fHistInvMassElecUnLike_EMC8", "Invariant mass for all unlike-signed electron pairs", 100, 0, .5);
    fHistInvMassElecUnLike_EMC8->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecUnLike_EMC8->GetYaxis()->SetTitle("Cts");
    
    fHistInvMassElecUnLike_EMCJet = new TH1F("fHistInvMassElecUnLike_EMCJet", "Invariant mass for all unlike-signed electron pairs", 100, 0, .5);
    fHistInvMassElecUnLike_EMCJet->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecUnLike_EMCJet->GetYaxis()->SetTitle("Cts");
    
    //Opening Angle Histos
    
    fHistOpAngElecLike_MB = new TH1F("fHistOpAngElecLike_MB", "Opening angle for all like-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecLike_MB->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecLike_MB->GetYaxis()->SetTitle("Cts");
    
    fHistOpAngElecLike_EMC7 = new TH1F("fHistOpAngElecLike_EMC7", "Opening angle for all like-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecLike_EMC7->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecLike_EMC7->GetYaxis()->SetTitle("Cts");
    
    fHistOpAngElecLike_EMC8 = new TH1F("fHistOpAngElecLike_EMC8", "Opening angle for all like-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecLike_EMC8->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecLike_EMC8->GetYaxis()->SetTitle("Cts");
    
    fHistOpAngElecLike_EMCJet = new TH1F("fHistOpAngElecLike_EMCJet", "Opening angle for all like-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecLike_EMCJet->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecLike_EMCJet->GetYaxis()->SetTitle("Cts");
    
    fHistOpAngElecUnLike_MB = new TH1F("fHistOpAngElecUnLike_MB", "Opening angle for all unlike-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecUnLike_MB->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecUnLike_MB->GetYaxis()->SetTitle("Cts");
    
    fHistOpAngElecUnLike_EMC7 = new TH1F("fHistOpAngElecUnLike_EMC7", "Opening angle for all unlike-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecUnLike_EMC7->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecUnLike_EMC7->GetYaxis()->SetTitle("Cts");
    
    fHistOpAngElecUnLike_EMC8 = new TH1F("fHistOpAngElecUnLike_EMC8", "Opening angle for all unlike-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecUnLike_EMC8->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecUnLike_EMC8->GetYaxis()->SetTitle("Cts");
    
    fHistOpAngElecUnLike_EMCJet = new TH1F("fHistOpAngElecUnLike_EMCJet", "Opening angle for all unlike-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecUnLike_EMCJet->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecUnLike_EMCJet->GetYaxis()->SetTitle("Cts");
    
    //PtSum Histos
    
    fHistPtSumTransMaxB2B = new TH2F("fHistPtSumTransMaxB2B", "PtSum density vs candidate electron Pt in the TransMax region for back-to-back events", 800, 0, 8, 1000, 0, 300);
    fHistPtSumTransMaxB2B->GetXaxis()->SetTitle("Candidate Pt");
    fHistPtSumTransMaxB2B->GetYaxis()->SetTitle("PtSum Density");
    
    fHistPtSumTransMinB2B = new TH2F("fHistPtSumTransMinB2B", "PtSum density vs candidate electron Pt in the TransMin region for back-to-back events", 800, 0, 8, 1000, 0, 300);
    fHistPtSumTransMinB2B->GetXaxis()->SetTitle("Candidate Pt");
    fHistPtSumTransMinB2B->GetYaxis()->SetTitle("PtSum Density");
    
    fHistPtSumTransMaxLead = new TH2F("fHistPtSumTransMaxLead", "PtSum density vs candidate electron Pt in the TransMax region for 'leading jet' events", 800, 0, 8, 1000, 0, 300);
    fHistPtSumTransMaxLead->GetXaxis()->SetTitle("Candidate Pt");
    fHistPtSumTransMaxLead->GetYaxis()->SetTitle("PtSum Density");
    
    fHistPtSumTransMinLead = new TH2F("fHistPtSumTransMinLead", "PtSum density vs candidate electron Pt in the TransMin region for 'leading jet' events", 800, 0, 8, 1000, 0, 300);
    fHistPtSumTransMinLead->GetXaxis()->SetTitle("Candidate Pt");
    fHistPtSumTransMinLead->GetYaxis()->SetTitle("PtSum Density");
    
    //Rejection Histos
    
    fHistPIDRejection = new TH1F("fHistPIDRejection", "PID rejection counts for each detector.", 4, 1, 4);
    fHistPIDRejection->GetXaxis()->SetTitle("Detector");
    fHistPIDRejection->GetYaxis()->SetTitle("Cts");
    fHistPIDRejection->GetXaxis()->SetBinLabel(1, "TPC");
    fHistPIDRejection->GetXaxis()->SetBinLabel(2, "TOF");
    fHistPIDRejection->GetXaxis()->SetBinLabel(3, "TRD");
    fHistPIDRejection->GetXaxis()->SetBinLabel(4, "EMC");
    
    //Make emcal bad cluster id histo
    fHistBadEMCclusID = new TH1F("fHistBadEMCclusID", "Number of EMCal clusters with ID -99999", 2, 1, 2);
    fHistBadEMCclusID->GetXaxis()->SetBinLabel(1, "Bad Clusters");
    fHistBadEMCclusID->GetYaxis()->SetTitle("Cts");
    
    //Number of electrons per event histo
    fHistNElecPerEvent = new TH1F("fHistNElecPerEvent", "Number of tagged electrons per event", 5, 1, 5);
    fHistNElecPerEvent->GetXaxis()->SetTitle("Num. of Electrons");
    fHistNElecPerEvent->GetYaxis()->SetTitle("Cts");
    
    //Region Histos
    
    for(Int_t i=0;i<4;i++){
           
           //Tag Side Histos
       
           //Multiplicity Histos
           fHistTrkMultTag_MB[i] = new TH1F(TString::Format("fHistTrkMultTag_MB_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the Tag Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultTag_MB[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultTag_MB[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkMultTag_EMC7[i] = new TH1F(TString::Format("fHistTrkMultTag_EMC7_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the Tag Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultTag_EMC7[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultTag_EMC7[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkMultTag_EMC8[i] = new TH1F(TString::Format("fHistTrkMultTag_EMC8_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the Tag Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultTag_EMC8[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultTag_EMC8[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkMultTag_EMCJet[i] = new TH1F(TString::Format("fHistTrkMultTag_EMCJet_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the Tag Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultTag_EMCJet[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultTag_EMCJet[i]->GetYaxis()->SetTitle("Cts");
            
           //Pt Histos
           fHistTrkPtTag_MB[i] = new TH1F(TString::Format("fHistTrkPtTag_MB_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the Tag Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtTag_MB[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtTag_MB[i]->GetYaxis()->SetTitle("Cts");
      
           fHistTrkPtTag_EMC7[i] = new TH1F(TString::Format("fHistTrkPtTag_EMC7_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the Tag Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtTag_EMC7[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtTag_EMC7[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkPtTag_EMC8[i] = new TH1F(TString::Format("fHistTrkPtTag_EMC8_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the Tag Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtTag_EMC8[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtTag_EMC8[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkPtTag_EMCJet[i] = new TH1F(TString::Format("fHistTrkPtTag_EMCJet_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the Tag Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtTag_EMCJet[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtTag_EMCJet[i]->GetYaxis()->SetTitle("Cts");
            
           //DeDx by Pt Histos
           fHistDeDxPtTag_MB[i] = new TH2F(TString::Format("fHistDeDxPtTag_MB_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the Tag Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtTag_MB[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtTag_MB[i]->GetYaxis()->SetTitle("TPC dE/dx");
       
           fHistDeDxPtTag_EMC7[i] = new TH2F(TString::Format("fHistDeDxPtTag_EMC7_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the Tag Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtTag_EMC7[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtTag_EMC7[i]->GetYaxis()->SetTitle("TPC dE/dx");
       
           fHistDeDxPtTag_EMC8[i] = new TH2F(TString::Format("fHistDeDxPtTag_EMC8_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the Tag Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtTag_EMC8[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtTag_EMC8[i]->GetYaxis()->SetTitle("TPC dE/dx");
       
           fHistDeDxPtTag_EMCJet[i] = new TH2F(TString::Format("fHistDeDxPtTag_EMCJet_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the Tag Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtTag_EMCJet[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtTag_EMCJet[i]->GetYaxis()->SetTitle("TPC dE/dx");
        
           //Away Side Histos
       
           //Multiplicity Histos
           fHistTrkMultAway_MB[i] = new TH1F(TString::Format("fHistTrkMultAway_MB_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the Away Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultAway_MB[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultAway_MB[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkMultAway_EMC7[i] = new TH1F(TString::Format("fHistTrkMultAway_EMC7_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the Away Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultAway_EMC7[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultAway_EMC7[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkMultAway_EMC8[i] = new TH1F(TString::Format("fHistTrkMultAway_EMC8_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the Away Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultAway_EMC8[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultAway_EMC8[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkMultAway_EMCJet[i] = new TH1F(TString::Format("fHistTrkMultAway_EMCJet_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the Away Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultAway_EMCJet[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultAway_EMCJet[i]->GetYaxis()->SetTitle("Cts");
            
           //Pt Histos
           fHistTrkPtAway_MB[i] = new TH1F(TString::Format("fHistTrkPtAway_MB_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the Away Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtAway_MB[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtAway_MB[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkPtAway_EMC7[i] = new TH1F(TString::Format("fHistTrkPtAway_EMC7_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the Away Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtAway_EMC7[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtAway_EMC7[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkPtAway_EMC8[i] = new TH1F(TString::Format("fHistTrkPtAway_EMC8_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the Away Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtAway_EMC8[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtAway_EMC8[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkPtAway_EMCJet[i] = new TH1F(TString::Format("fHistTrkPtAway_EMCJet_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the Away Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtAway_EMCJet[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtAway_EMCJet[i]->GetYaxis()->SetTitle("Cts");
            
           //DeDx by Pt Histos
           fHistDeDxPtAway_MB[i] = new TH2F(TString::Format("fHistDeDxPtAway_MB_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the Away Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtAway_MB[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtAway_MB[i]->GetYaxis()->SetTitle("TPC dE/dx");
       
           fHistDeDxPtAway_EMC7[i] = new TH2F(TString::Format("fHistDeDxPtAway_EMC7_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the Away Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtAway_EMC7[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtAway_EMC7[i]->GetYaxis()->SetTitle("TPC dE/dx");
       
           fHistDeDxPtAway_EMC8[i] = new TH2F(TString::Format("fHistDeDxPtAway_EMC8_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the Away Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtAway_EMC8[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtAway_EMC8[i]->GetYaxis()->SetTitle("TPC dE/dx");
       
           fHistDeDxPtAway_EMCJet[i] = new TH2F(TString::Format("fHistDeDxPtAway_EMCJet_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the Away Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtAway_EMCJet[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtAway_EMCJet[i]->GetYaxis()->SetTitle("TPC dE/dx");
        
           //transMin Side Histos
       
           //Multiplicity Histos
           fHistTrkMultTransMin_MB[i] = new TH1F(TString::Format("fHistTrkMultTransMin_MB_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the TransMin Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultTransMin_MB[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultTransMin_MB[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkMultTransMin_EMC7[i] = new TH1F(TString::Format("fHistTrkMultTransMin_EMC7_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the TransMin Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultTransMin_EMC7[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultTransMin_EMC7[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkMultTransMin_EMC8[i] = new TH1F(TString::Format("fHistTrkMultTransMin_EMC8_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the TransMin Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultTransMin_EMC8[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultTransMin_EMC8[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkMultTransMin_EMCJet[i] = new TH1F(TString::Format("fHistTrkMultTransMin_EMCJet_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the TransMin Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultTransMin_EMCJet[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultTransMin_EMCJet[i]->GetYaxis()->SetTitle("Cts");
            
           //Pt Histos
           fHistTrkPtTransMin_MB[i] = new TH1F(TString::Format("fHistTrkPtTransMin_MB_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the TransMin Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtTransMin_MB[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtTransMin_MB[i]->GetYaxis()->SetTitle("Cts");
      
           fHistTrkPtTransMin_EMC7[i] = new TH1F(TString::Format("fHistTrkPtTransMin_EMC7_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the TransMin Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtTransMin_EMC7[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtTransMin_EMC7[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkPtTransMin_EMC8[i] = new TH1F(TString::Format("fHistTrkPtTransMin_EMC8_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the TransMin Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtTransMin_EMC8[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtTransMin_EMC8[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkPtTransMin_EMCJet[i] = new TH1F(TString::Format("fHistTrkPtTransMin_EMCJet_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the TransMin Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtTransMin_EMCJet[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtTransMin_EMCJet[i]->GetYaxis()->SetTitle("Cts");
            
           //DeDx by Pt Histos
           fHistDeDxPtTransMin_MB[i] = new TH2F(TString::Format("fHistDeDxPtTransMin_MB_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the TransMin Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtTransMin_MB[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtTransMin_MB[i]->GetYaxis()->SetTitle("TPC dE/dx");
       
           fHistDeDxPtTransMin_EMC7[i] = new TH2F(TString::Format("fHistDeDxPtTransMin_EMC7_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the TransMin Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtTransMin_EMC7[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtTransMin_EMC7[i]->GetYaxis()->SetTitle("TPC dE/dx");
       
           fHistDeDxPtTransMin_EMC8[i] = new TH2F(TString::Format("fHistDeDxPtTransMin_EMC8_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the TransMin Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtTransMin_EMC8[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtTransMin_EMC8[i]->GetYaxis()->SetTitle("TPC dE/dx");
       
           fHistDeDxPtTransMin_EMCJet[i] = new TH2F(TString::Format("fHistDeDxPtTransMin_EMCJet_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the TransMin Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtTransMin_EMCJet[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtTransMin_EMCJet[i]->GetYaxis()->SetTitle("TPC dE/dx");
            
           //TransMax Side Histos
       
           //Multiplicity Histos
           fHistTrkMultTransMax_MB[i] = new TH1F(TString::Format("fHistTrkMultTransMax_MB_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the TransMax Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultTransMax_MB[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultTransMax_MB[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkMultTransMax_EMC7[i] = new TH1F(TString::Format("fHistTrkMultTransMax_EMC7_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the TransMax Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultTransMax_EMC7[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultTransMax_EMC7[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkMultTransMax_EMC8[i] = new TH1F(TString::Format("fHistTrkMultTransMax_EMC8_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the TransMax Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultTransMax_EMC8[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultTransMax_EMC8[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkMultTransMax_EMCJet[i] = new TH1F(TString::Format("fHistTrkMultTransMax_EMCJet_%s", ptRangesRegion[i].Data()), TString::Format("Multiplicity of Tracks on the TransMax Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 400, 0, 400);
           fHistTrkMultTransMax_EMCJet[i]->GetXaxis()->SetTitle("Track Multiplicity");
           fHistTrkMultTransMax_EMCJet[i]->GetYaxis()->SetTitle("Cts");
            
           //Pt Histos
           fHistTrkPtTransMax_MB[i] = new TH1F(TString::Format("fHistTrkPtTransMax_MB_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the TransMax Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtTransMax_MB[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtTransMax_MB[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkPtTransMax_EMC7[i] = new TH1F(TString::Format("fHistTrkPtTransMax_EMC7_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the TransMax Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtTransMax_EMC7[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtTransMax_EMC7[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkPtTransMax_EMC8[i] = new TH1F(TString::Format("fHistTrkPtTransMax_EMC8_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the TransMax Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtTransMax_EMC8[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtTransMax_EMC8[i]->GetYaxis()->SetTitle("Cts");
       
           fHistTrkPtTransMax_EMCJet[i] = new TH1F(TString::Format("fHistTrkPtTransMax_EMCJet_%s", ptRangesRegion[i].Data()), TString::Format("Pt Distribution of Tracks on the TransMax Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8);
           fHistTrkPtTransMax_EMCJet[i]->GetXaxis()->SetTitle("Track Pt");
           fHistTrkPtTransMax_EMCJet[i]->GetYaxis()->SetTitle("Cts");
            
           //DeDx by Pt Histos
           fHistDeDxPtTransMax_MB[i] = new TH2F(TString::Format("fHistDeDxPtTransMax_MB_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the TransMax Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtTransMax_MB[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtTransMax_MB[i]->GetYaxis()->SetTitle("TPC dE/dx");
       
           fHistDeDxPtTransMax_EMC7[i] = new TH2F(TString::Format("fHistDeDxPtTransMax_EMC7_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the TransMax Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtTransMax_EMC7[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtTransMax_EMC7[i]->GetYaxis()->SetTitle("TPC dE/dx");
       
           fHistDeDxPtTransMax_EMC8[i] = new TH2F(TString::Format("fHistDeDxPtTransMax_EMC8_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the TransMax Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtTransMax_EMC8[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtTransMax_EMC8[i]->GetYaxis()->SetTitle("TPC dE/dx");
       
           fHistDeDxPtTransMax_EMCJet[i] = new TH2F(TString::Format("fHistDeDxPtTransMax_EMCJet_%s", ptRangesRegion[i].Data()), TString::Format("DeDx by Pt for Tracks on the TransMax Side of Event with a %s pt electron", ptRangesRegion[i].Data()), 800, 0, 8, 300, 0, 12);
           fHistDeDxPtTransMax_EMCJet[i]->GetXaxis()->SetTitle("Track Pt");
           fHistDeDxPtTransMax_EMCJet[i]->GetYaxis()->SetTitle("TPC dE/dx");
    }
    
    //PID Plots
    
    //TPC PID Plots
    for(Int_t i=0; i<6; i++){
        //MB
        fHistTPC_TOF_MB[i] = new TH2F(TString::Format("fHistTPC_TOF_MB_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOF_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOF_MB[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_EMC_MB[i] = new TH2F(TString::Format("fHistTPC_EMC_MB_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_EMC_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_EMC_MB[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TRD_MB[i] = new TH2F(TString::Format("fHistTPC_TRD_MB_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TRD_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TRD_MB[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TOFEMC_MB[i] = new TH2F(TString::Format("fHistTPC_TOFEMC_MB_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOFEMC_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOFEMC_MB[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TOFTRD_MB[i] = new TH2F(TString::Format("fHistTPC_TOFTRD_MB_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOFTRD_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOFTRD_MB[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_EMCTRD_MB[i] = new TH2F(TString::Format("fHistTPC_EMCTRD_MB_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_EMCTRD_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_EMCTRD_MB[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TOFEMCTRD_MB[i] = new TH2F(TString::Format("fHistTPC_TOFEMCTRD_MB_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF, EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOFEMCTRD_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOFEMCTRD_MB[i]->GetYaxis()->SetTitle("nSigma");
    
        //EMC7
        fHistTPC_TOF_EMC7[i] = new TH2F(TString::Format("fHistTPC_TOF_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOF_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOF_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_EMC_EMC7[i] = new TH2F(TString::Format("fHistTPC_EMC_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_EMC_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_EMC_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TRD_EMC7[i] = new TH2F(TString::Format("fHistTPC_TRD_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TRD_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TRD_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TOFEMC_EMC7[i] = new TH2F(TString::Format("fHistTPC_TOFEMC_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOFEMC_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOFEMC_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TOFTRD_EMC7[i] = new TH2F(TString::Format("fHistTPC_TOFTRD_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOFTRD_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOFTRD_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_EMCTRD_EMC7[i] = new TH2F(TString::Format("fHistTPC_EMCTRD_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_EMCTRD_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_EMCTRD_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TOFEMCTRD_EMC7[i] = new TH2F(TString::Format("fHistTPC_TOFEMCTRD_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF, EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOFEMCTRD_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOFEMCTRD_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        //EMC8
        fHistTPC_TOF_EMC8[i] = new TH2F(TString::Format("fHistTPC_TOF_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOF_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOF_EMC8[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_EMC_EMC8[i] = new TH2F(TString::Format("fHistTPC_EMC_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_EMC_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_EMC_EMC8[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TRD_EMC8[i] = new TH2F(TString::Format("fHistTPC_TRD_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TRD_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TRD_EMC8[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TOFEMC_EMC8[i] = new TH2F(TString::Format("fHistTPC_TOFEMC_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOFEMC_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOFEMC_EMC8[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TOFTRD_EMC8[i] = new TH2F(TString::Format("fHistTPC_TOFTRD_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOFTRD_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOFTRD_EMC8[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_EMCTRD_EMC8[i] = new TH2F(TString::Format("fHistTPC_EMCTRD_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_EMCTRD_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_EMCTRD_EMC8[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TOFEMCTRD_EMC8[i] = new TH2F(TString::Format("fHistTPC_TOFEMCTRD_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF, EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOFEMCTRD_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOFEMCTRD_EMC8[i]->GetYaxis()->SetTitle("nSigma");
        
        //EMCJet
        fHistTPC_TOF_EMCJet[i] = new TH2F(TString::Format("fHistTPC_TOF_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOF_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOF_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_EMC_EMCJet[i] = new TH2F(TString::Format("fHistTPC_EMC_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_EMC_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_EMC_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TRD_EMCJet[i] = new TH2F(TString::Format("fHistTPC_TRD_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TRD_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TRD_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TOFEMC_EMCJet[i] = new TH2F(TString::Format("fHistTPC_TOFEMC_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOFEMC_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOFEMC_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TOFTRD_EMCJet[i] = new TH2F(TString::Format("fHistTPC_TOFTRD_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOFTRD_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOFTRD_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_EMCTRD_EMCJet[i] = new TH2F(TString::Format("fHistTPC_EMCTRD_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_EMCTRD_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_EMCTRD_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTPC_TOFEMCTRD_EMCJet[i] = new TH2F(TString::Format("fHistTPC_TOFEMCTRD_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after TOF, EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_TOFEMCTRD_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_TOFEMCTRD_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
    }
    
    //TOF PID Plots
    for(Int_t i=0; i<6; i++){
        //MB
        fHistTOF_TPC_MB[i] = new TH2F(TString::Format("fHistTOF_TPC_MB_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPC_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPC_MB[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_EMC_MB[i] = new TH2F(TString::Format("fHistTOF_EMC_MB_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_EMC_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_EMC_MB[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TRD_MB[i] = new TH2F(TString::Format("fHistTOF_TRD_MB_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TRD_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TRD_MB[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TPCEMC_MB[i] = new TH2F(TString::Format("fHistTOF_TPCEMC_MB_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPCEMC_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPCEMC_MB[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TPCTRD_MB[i] = new TH2F(TString::Format("fHistTOF_TPCTRD_MB_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPCTRD_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPCTRD_MB[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_EMCTRD_MB[i] = new TH2F(TString::Format("fHistTOF_EMCTRD_MB_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_EMCTRD_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_EMCTRD_MB[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TPCEMCTRD_MB[i] = new TH2F(TString::Format("fHistTOF_TPCEMCTRD_MB_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC, EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPCEMCTRD_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPCEMCTRD_MB[i]->GetYaxis()->SetTitle("nSigma");
    
        //EMC7
        fHistTOF_TPC_EMC7[i] = new TH2F(TString::Format("fHistTOF_TPC_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPC_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPC_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_EMC_EMC7[i] = new TH2F(TString::Format("fHistTOF_EMC_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_EMC_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_EMC_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TRD_EMC7[i] = new TH2F(TString::Format("fHistTOF_TRD_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TRD_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TRD_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TPCEMC_EMC7[i] = new TH2F(TString::Format("fHistTOF_TPCEMC_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPCEMC_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPCEMC_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TPCTRD_EMC7[i] = new TH2F(TString::Format("fHistTOF_TPCTRD_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPCTRD_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPCTRD_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_EMCTRD_EMC7[i] = new TH2F(TString::Format("fHistTOF_EMCTRD_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_EMCTRD_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_EMCTRD_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TPCEMCTRD_EMC7[i] = new TH2F(TString::Format("fHistTOF_TPCEMCTRD_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC, EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPCEMCTRD_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPCEMCTRD_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        //EMC8
        fHistTOF_TPC_EMC8[i] = new TH2F(TString::Format("fHistTOF_TPC_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPC_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPC_EMC8[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_EMC_EMC8[i] = new TH2F(TString::Format("fHistTOF_EMC_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_EMC_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_EMC_EMC8[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TRD_EMC8[i] = new TH2F(TString::Format("fHistTOF_TRD_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TRD_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TRD_EMC8[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TPCEMC_EMC8[i] = new TH2F(TString::Format("fHistTOF_TPCEMC_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPCEMC_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPCEMC_EMC8[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TPCTRD_EMC8[i] = new TH2F(TString::Format("fHistTOF_TPCTRD_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPCTRD_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPCTRD_EMC8[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_EMCTRD_EMC8[i] = new TH2F(TString::Format("fHistTOF_EMCTRD_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_EMCTRD_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_EMCTRD_EMC8[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TPCEMCTRD_EMC8[i] = new TH2F(TString::Format("fHistTOF_TPCEMCTRD_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC, EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPCEMCTRD_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPCEMCTRD_EMC8[i]->GetYaxis()->SetTitle("nSigma");
        
        //EMCJet
        fHistTOF_TPC_EMCJet[i] = new TH2F(TString::Format("fHistTOF_TPC_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPC_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPC_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_EMC_EMCJet[i] = new TH2F(TString::Format("fHistTOF_EMC_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_EMC_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_EMC_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TRD_EMCJet[i] = new TH2F(TString::Format("fHistTOF_TRD_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TRD_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TRD_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TPCEMC_EMCJet[i] = new TH2F(TString::Format("fHistTOF_TPCEMC_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPCEMC_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPCEMC_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TPCTRD_EMCJet[i] = new TH2F(TString::Format("fHistTOF_TPCTRD_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TPC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPCTRD_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPCTRD_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_EMCTRD_EMCJet[i] = new TH2F(TString::Format("fHistTOF_EMCTRD_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_EMCTRD_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_EMCTRD_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
        
        fHistTOF_TPCEMCTRD_EMCJet[i] = new TH2F(TString::Format("fHistTOF_TPCEMCTRD_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TOF nSigma for tracks with Pt between %s after TOF, EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTOF_TPCEMCTRD_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTOF_TPCEMCTRD_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
    }
    
    //EMC PID Plots
    for(Int_t i=0; i<6; i++){
        //MB
        fHistEMC_TOF_MB[i] = new TH1F(TString::Format("fHistEMC_TOF_MB_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TOF_MB[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TOF_MB[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPC_MB[i] = new TH1F(TString::Format("fHistEMC_TPC_MB_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TPC cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPC_MB[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPC_MB[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TRD_MB[i] = new TH1F(TString::Format("fHistEMC_TRD_MB_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TRD_MB[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TRD_MB[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPCTOF_MB[i] = new TH1F(TString::Format("fHistEMC_TPCTOF_MB_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF and TPC cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTOF_MB[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTOF_MB[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TOFTRD_MB[i] = new TH1F(TString::Format("fHistEMC_TOFTRD_MB_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TOFTRD_MB[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TOFTRD_MB[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPCTRD_MB[i] = new TH1F(TString::Format("fHistEMC_TPCTRD_MB_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TPC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTRD_MB[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTRD_MB[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPCTOFTRD_MB[i] = new TH1F(TString::Format("fHistEMC_TPCTOFTRD_MB_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF, TPC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTOFTRD_MB[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTOFTRD_MB[i]->GetYaxis()->SetTitle("Cts");
    
        //EMC7
        fHistEMC_TOF_EMC7[i] = new TH1F(TString::Format("fHistEMC_TOF_EMC7_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TOF_EMC7[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TOF_EMC7[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPC_EMC7[i] = new TH1F(TString::Format("fHistEMC_TPC_EMC7_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TPC cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPC_EMC7[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPC_EMC7[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TRD_EMC7[i] = new TH1F(TString::Format("fHistEMC_TRD_EMC7_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TRD_EMC7[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TRD_EMC7[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPCTOF_EMC7[i] = new TH1F(TString::Format("fHistEMC_TPCTOF_EMC7_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF and TPC cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTOF_EMC7[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTOF_EMC7[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TOFTRD_EMC7[i] = new TH1F(TString::Format("fHistEMC_TOFTRD_EMC7_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TOFTRD_EMC7[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TOFTRD_EMC7[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPCTRD_EMC7[i] = new TH1F(TString::Format("fHistEMC_TPCTRD_EMC7_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TPC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTRD_EMC7[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTRD_EMC7[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPCTOFTRD_EMC7[i] = new TH1F(TString::Format("fHistEMC_TPCTOFTRD_EMC7_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF, TPC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTOFTRD_EMC7[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTOFTRD_EMC7[i]->GetYaxis()->SetTitle("Cts");
        
        //EMC8
        fHistEMC_TOF_EMC8[i] = new TH1F(TString::Format("fHistEMC_TOF_EMC8_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TOF_EMC8[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TOF_EMC8[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPC_EMC8[i] = new TH1F(TString::Format("fHistEMC_TPC_EMC8_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TPC cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPC_EMC8[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPC_EMC8[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TRD_EMC8[i] = new TH1F(TString::Format("fHistEMC_TRD_EMC8_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TRD_EMC8[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TRD_EMC8[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPCTOF_EMC8[i] = new TH1F(TString::Format("fHistEMC_TPCTOF_EMC8_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF and TPC cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTOF_EMC8[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTOF_EMC8[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TOFTRD_EMC8[i] = new TH1F(TString::Format("fHistEMC_TOFTRD_EMC8_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TOFTRD_EMC8[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TOFTRD_EMC8[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPCTRD_EMC8[i] = new TH1F(TString::Format("fHistEMC_TPCTRD_EMC8_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TPC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTRD_EMC8[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTRD_EMC8[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPCTOFTRD_EMC8[i] = new TH1F(TString::Format("fHistEMC_TPCTOFTRD_EMC8_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF, TPC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTOFTRD_EMC8[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTOFTRD_EMC8[i]->GetYaxis()->SetTitle("Cts");
        
        //EMCJet
        fHistEMC_TOF_EMCJet[i] = new TH1F(TString::Format("fHistEMC_TOF_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TOF_EMCJet[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TOF_EMCJet[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPC_EMCJet[i] = new TH1F(TString::Format("fHistEMC_TPC_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TPC cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPC_EMCJet[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPC_EMCJet[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TRD_EMCJet[i] = new TH1F(TString::Format("fHistEMC_TRD_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TRD_EMCJet[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TRD_EMCJet[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPCTOF_EMCJet[i] = new TH1F(TString::Format("fHistEMC_TPCTOF_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF and TPC cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTOF_EMCJet[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTOF_EMCJet[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TOFTRD_EMCJet[i] = new TH1F(TString::Format("fHistEMC_TOFTRD_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TOFTRD_EMCJet[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TOFTRD_EMCJet[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPCTRD_EMCJet[i] = new TH1F(TString::Format("fHistEMC_TPCTRD_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTRD_EMCJet[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTRD_EMCJet[i]->GetYaxis()->SetTitle("Cts");
        
        fHistEMC_TPCTOFTRD_EMCJet[i] = new TH1F(TString::Format("fHistEMC_TPCTOFTRD_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TOF, TPC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTOFTRD_EMCJet[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTOFTRD_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    //TRD PID Plots
    for(Int_t i=0; i<6; i++){
        //MB
        fHistTRD_TOF_MB[i] = new TH2F(TString::Format("fHistTRD_TOF_MB_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TOF_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TOF_MB[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPC_MB[i] = new TH2F(TString::Format("fHistTRD_TPC_MB_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TPC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPC_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPC_MB[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_EMC_MB[i] = new TH2F(TString::Format("fHistTRD_EMC_MB_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_EMC_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_EMC_MB[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPCTOF_MB[i] = new TH2F(TString::Format("fHistTRD_TPCTOF_MB_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF and TPC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCTOF_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCTOF_MB[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TOFEMC_MB[i] = new TH2F(TString::Format("fHistTRD_TOFEMC_MB_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TOFEMC_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TOFEMC_MB[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPCEMC_MB[i] = new TH2F(TString::Format("fHistTRD_TPCEMC_MB_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TPC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCEMC_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCEMC_MB[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPCTOFEMC_MB[i] = new TH2F(TString::Format("fHistTRD_TPCTOFEMC_MB_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF, TPC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCTOFEMC_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCTOFEMC_MB[i]->GetYaxis()->SetTitle("electron Likelihood");
    
        //EMC7
        fHistTRD_TOF_EMC7[i] = new TH2F(TString::Format("fHistTRD_TOF_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TOF_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TOF_EMC7[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPC_EMC7[i] = new TH2F(TString::Format("fHistTRD_TPC_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TPC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPC_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPC_EMC7[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_EMC_EMC7[i] = new TH2F(TString::Format("fHistTRD_EMC_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_EMC_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_EMC_EMC7[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPCTOF_EMC7[i] = new TH2F(TString::Format("fHistTRD_TPCTOF_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF and TPC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCTOF_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCTOF_EMC7[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TOFEMC_EMC7[i] = new TH2F(TString::Format("fHistTRD_TOFEMC_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TOFEMC_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TOFEMC_EMC7[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPCEMC_EMC7[i] = new TH2F(TString::Format("fHistTRD_TPCEMC_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TPC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCEMC_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCEMC_EMC7[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPCTOFEMC_EMC7[i] = new TH2F(TString::Format("fHistTRD_TPCTOFEMC_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF, TPC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCTOFEMC_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCTOFEMC_EMC7[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        //EMC8
        fHistTRD_TOF_EMC8[i] = new TH2F(TString::Format("fHistTRD_TOF_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TOF_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TOF_EMC8[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPC_EMC8[i] = new TH2F(TString::Format("fHistTRD_TPC_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TPC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPC_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPC_EMC8[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_EMC_EMC8[i] = new TH2F(TString::Format("fHistTRD_EMC_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_EMC_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_EMC_EMC8[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPCTOF_EMC8[i] = new TH2F(TString::Format("fHistTRD_TPCTOF_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF and TPC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCTOF_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCTOF_EMC8[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TOFEMC_EMC8[i] = new TH2F(TString::Format("fHistTRD_TOFEMC_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TOFEMC_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TOFEMC_EMC8[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPCEMC_EMC8[i] = new TH2F(TString::Format("fHistTRD_TPCEMC_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TPC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCEMC_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCEMC_EMC8[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPCTOFEMC_EMC8[i] = new TH2F(TString::Format("fHistTRD_TPCTOFEMC_EMC8_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF, TPC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCTOFEMC_EMC8[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCTOFEMC_EMC8[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        //EMCJet
        fHistTRD_TOF_EMCJet[i] = new TH2F(TString::Format("fHistTRD_TOF_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TOF_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TOF_EMCJet[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPC_EMCJet[i] = new TH2F(TString::Format("fHistTRD_TPC_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TPC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPC_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPC_EMCJet[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_EMC_EMCJet[i] = new TH2F(TString::Format("fHistTRD_EMC_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_EMC_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_EMC_EMCJet[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPCTOF_EMCJet[i] = new TH2F(TString::Format("fHistTRD_TPCTOF_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF and TPC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCTOF_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCTOF_EMCJet[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TOFEMC_EMCJet[i] = new TH2F(TString::Format("fHistTRD_TOFEMC_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TOFEMC_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TOFEMC_EMCJet[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPCEMC_EMCJet[i] = new TH2F(TString::Format("fHistTRD_TPCEMC_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after EMC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCEMC_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCEMC_EMCJet[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        fHistTRD_TPCTOFEMC_EMCJet[i] = new TH2F(TString::Format("fHistTRD_TPCTOFEMC_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TOF, TPC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCTOFEMC_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCTOFEMC_EMCJet[i]->GetYaxis()->SetTitle("electron Likelihood");
    }
    
    //EMC Shower Shape PID Plots
    for(Int_t i=0; i<6; i++){
        //MB
        fHistM02_All_MB[i] = new TH2F(TString::Format("fHistM02_All_MB_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M02(semi-major axis) for tracks with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM02_All_MB[i]->GetXaxis()->SetTitle("M02");
        fHistM02_All_MB[i]->GetYaxis()->SetTitle("E/p");
        
        fHistM20_All_MB[i] = new TH2F(TString::Format("fHistM20_All_MB_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M20(semi-minor axis) for tracks with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM20_All_MB[i]->GetXaxis()->SetTitle("M20");
        fHistM20_All_MB[i]->GetYaxis()->SetTitle("E/p");
        
        fHistM02_Elec_MB[i] = new TH2F(TString::Format("fHistM02_Elec_MB_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M02(semi-major axis) for electron candidates with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM02_Elec_MB[i]->GetXaxis()->SetTitle("M02");
        fHistM02_Elec_MB[i]->GetYaxis()->SetTitle("E/p");
        
        fHistM20_Elec_MB[i] = new TH2F(TString::Format("fHistM20_Elec_MB_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M20(semi-minor axis) for electron candidates with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM20_Elec_MB[i]->GetXaxis()->SetTitle("M20");
        fHistM20_Elec_MB[i]->GetYaxis()->SetTitle("E/p");
        
        //EMC7
        fHistM02_All_EMC7[i] = new TH2F(TString::Format("fHistM02_All_EMC7_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M02(semi-major axis) for tracks with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM02_All_EMC7[i]->GetXaxis()->SetTitle("M02");
        fHistM02_All_EMC7[i]->GetYaxis()->SetTitle("E/p");
        
        fHistM20_All_EMC7[i] = new TH2F(TString::Format("fHistM20_All_EMC7_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M20(semi-minor axis) for tracks with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM20_All_EMC7[i]->GetXaxis()->SetTitle("M20");
        fHistM20_All_EMC7[i]->GetYaxis()->SetTitle("E/p");
        
        fHistM02_Elec_EMC7[i] = new TH2F(TString::Format("fHistM02_Elec_EMC7_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M02(semi-major axis) for electron candidates with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM02_Elec_EMC7[i]->GetXaxis()->SetTitle("M02");
        fHistM02_Elec_EMC7[i]->GetYaxis()->SetTitle("E/p");
        
        fHistM20_Elec_EMC7[i] = new TH2F(TString::Format("fHistM20_Elec_EMC7_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M20(semi-minor axis) for electron candidates with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM20_Elec_EMC7[i]->GetXaxis()->SetTitle("M20");
        fHistM20_Elec_EMC7[i]->GetYaxis()->SetTitle("E/p");
        
        //EMC8
        fHistM02_All_EMC8[i] = new TH2F(TString::Format("fHistM02_All_EMC8_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M02(semi-major axis) for tracks with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM02_All_EMC8[i]->GetXaxis()->SetTitle("M02");
        fHistM02_All_EMC8[i]->GetYaxis()->SetTitle("E/p");
        
        fHistM20_All_EMC8[i] = new TH2F(TString::Format("fHistM20_All_EMC8_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M20(semi-minor axis) for tracks with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM20_All_EMC8[i]->GetXaxis()->SetTitle("M20");
        fHistM20_All_EMC8[i]->GetYaxis()->SetTitle("E/p");
        
        fHistM02_Elec_EMC8[i] = new TH2F(TString::Format("fHistM02_Elec_EMC8_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M02(semi-major axis) for electron candidates with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM02_Elec_EMC8[i]->GetXaxis()->SetTitle("M02");
        fHistM02_Elec_EMC8[i]->GetYaxis()->SetTitle("E/p");
        
        fHistM20_Elec_EMC8[i] = new TH2F(TString::Format("fHistM20_Elec_EMC8_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M20(semi-minor axis) for electron candidates with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM20_Elec_EMC8[i]->GetXaxis()->SetTitle("M20");
        fHistM20_Elec_EMC8[i]->GetYaxis()->SetTitle("E/p");
        
        //EMCJet
        fHistM02_All_EMCJet[i] = new TH2F(TString::Format("fHistM02_All_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M02(semi-major axis) for tracks with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM02_All_EMCJet[i]->GetXaxis()->SetTitle("M02");
        fHistM02_All_EMCJet[i]->GetYaxis()->SetTitle("E/p");
        
        fHistM20_All_EMCJet[i] = new TH2F(TString::Format("fHistM20_All_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M20(semi-minor axis) for tracks with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM20_All_EMCJet[i]->GetXaxis()->SetTitle("M20");
        fHistM20_All_EMCJet[i]->GetYaxis()->SetTitle("E/p");
        
        fHistM02_Elec_EMCJet[i] = new TH2F(TString::Format("fHistM02_Elec_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M02(semi-major axis) for electron candidates with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM02_Elec_EMCJet[i]->GetXaxis()->SetTitle("M02");
        fHistM02_Elec_EMCJet[i]->GetYaxis()->SetTitle("E/p");
        
        fHistM20_Elec_EMCJet[i] = new TH2F(TString::Format("fHistM20_Elec_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("E/p vs M20(semi-minor axis) for electron candidates with Pt between %s",ptRangesPID[i].Data()), 100, 0, 1.5, 300, 0, 2);
        fHistM20_Elec_EMCJet[i]->GetXaxis()->SetTitle("M20");
        fHistM20_Elec_EMCJet[i]->GetYaxis()->SetTitle("E/p");
    }
    
    //DPhi for candidate electrons 2-8 gev and assoc. particles >3gev
    fHistDPhi28_MB = new TH1F("fHistDPhi28_MB", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhi28_MB->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi28_MB->GetYaxis()->SetTitle("Cts");
    
    fHistDPhi28_EMC7 = new TH1F("fHistDPhi28_EMC7", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhi28_EMC7->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi28_EMC7->GetYaxis()->SetTitle("Cts");
    
    fHistDPhi28_EMC8 = new TH1F("fHistDPhi28_EMC8", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhi28_EMC8->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi28_EMC8->GetYaxis()->SetTitle("Cts");
    
    fHistDPhi28_EMCJet = new TH1F("fHistDPhi28_EMCJet", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhi28_EMCJet->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi28_EMCJet->GetYaxis()->SetTitle("Cts");
    
    //DPhi by Eta for triggered particles 2-8 gev and assoc. particles >3gev
    fHistDPhiDEta28_MB = new TH2F("fHistDPhiDEta28_MB", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2, 100, -0.9, 0.9);
    fHistDPhiDEta28_MB->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEta28_MB->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEta28_MB->GetZaxis()->SetTitle("Cts");
    
    fHistDPhiDEta28_EMC7 = new TH2F("fHistDPhiDEta28_EMC7", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2, 100, -0.9, 0.9);
    fHistDPhiDEta28_EMC7->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEta28_EMC7->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEta28_EMC7->GetZaxis()->SetTitle("Cts");
    
    fHistDPhiDEta28_EMC8 = new TH2F("fHistDPhiDEta28_EMC8", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2, 100, -0.9, 0.9);
    fHistDPhiDEta28_EMC8->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEta28_EMC8->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEta28_EMC8->GetZaxis()->SetTitle("Cts");
    
    fHistDPhiDEta28_EMCJet = new TH2F("fHistDPhiDEta28_EMCJet", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2, 100, -0.9, 0.9);
    fHistDPhiDEta28_EMCJet->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEta28_EMCJet->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEta28_EMCJet->GetZaxis()->SetTitle("Cts");
        
    // Delta Phi for tracks > 300MeV
    for(Int_t i=0; i<3; i++){
        fHistDPhi300_500_MB[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.3-.5Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<.5Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi300_500_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi300_500_MB[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi300_500_EMC7[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.3-.5Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<.5Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi300_500_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi300_500_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi300_500_EMC8[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.3-.5Gev_EMC8",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<.5Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi300_500_EMC8[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi300_500_EMC8[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi300_500_EMCJet[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.3-.5Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<.5Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi300_500_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi300_500_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    // Delta Phi for tracks > 500MeV
    for(Int_t i=0; i<3; i++){
        fHistDPhi500_800_MB[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.5-.8Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .5Gev<pt<.8Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi500_800_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi500_800_MB[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi500_800_EMC7[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.5-.8Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .5Gev<pt<.8Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi500_800_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi500_800_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi500_800_EMC8[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.5-.8Gev_EMC8",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .5Gev<pt<.8Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi500_800_EMC8[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi500_800_EMC8[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi500_800_EMCJet[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.5-.8Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .5Gev<pt<.8Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi500_800_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi500_800_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    // Delta Phi for tracks > 800MeV
    for(Int_t i=0; i<3; i++){
        fHistDPhi800_1_MB[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.8-1Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .8Gev<pt<1Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi800_1_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi800_1_MB[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi800_1_EMC7[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.8-1Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .8Gev<pt<1Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi800_1_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi800_1_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi800_1_EMC8[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.8-1Gev_EMC8",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .8Gev<pt<1Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi800_1_EMC8[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi800_1_EMC8[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi800_1_EMCJet[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.8-1Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .8Gev<pt<1Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi800_1_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi800_1_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    // Delta Phi for tracks > 1GeV
    for(Int_t i=0; i<3; i++){
        fHistDPhi1_2_MB[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_1-2Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi1_2_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi1_2_MB[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi1_2_EMC7[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_1-2Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi1_2_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi1_2_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi1_2_EMC8[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_1-2Gev_EMC8",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi1_2_EMC8[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi1_2_EMC8[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi1_2_EMCJet[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_1-2Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi1_2_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi1_2_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    // Delta Phi for tracks > 2GeV
    for(Int_t i=0; i<3; i++){
        fHistDPhi2_3_MB[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_2-3Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<3Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi2_3_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi2_3_MB[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi2_3_EMC7[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_2-3Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<3Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi2_3_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi2_3_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi2_3_EMC8[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_2-3Gev_EMC8",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<3Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi2_3_EMC8[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi2_3_EMC8[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi2_3_EMCJet[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_2-3Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<3Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi2_3_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi2_3_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    // Delta Phi for tracks > 3GeV
    for(Int_t i=0; i<3; i++){
        fHistDPhi3_4_MB[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_3-4Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 3Gev<pt<4Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi3_4_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi3_4_MB[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi3_4_EMC7[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_3-4Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 3Gev<pt<4Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi3_4_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi3_4_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi3_4_EMC8[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_3-4Gev_EMC8",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 3Gev<pt<4Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi3_4_EMC8[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi3_4_EMC8[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi3_4_EMCJet[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_3-4Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 3Gev<pt<4Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi3_4_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi3_4_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    // Delta Phi for tracks > 5GeV
    for(Int_t i=0; i<3; i++){
        fHistDPhi4_MB[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_4Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated pt>4Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi4_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi4_MB[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi4_EMC7[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_4Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated pt>4Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi4_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi4_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi4_EMC8[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_4Gev_EMC8",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated pt>4Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi4_EMC8[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi4_EMC8[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhi4_EMCJet[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_4Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated pt>4Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi4_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi4_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    //Hadron e/p plot
    fHistEMC_Had_MB_1Gev = new TH1F("fHistEMC_Had_MB_1Gev", "E/p for hadrons with Pt between 1-2Gev", 100, 0, 1.5);
    fHistEMC_Had_MB_1Gev->GetXaxis()->SetTitle("E/p");
    fHistEMC_Had_MB_1Gev->GetYaxis()->SetTitle("Cts");
    
    // Eta-Phi distribution for tagged events
    fHistEtaPhiTag_MB = new TH2F("fHistEtaPhiTag_MB", "Eta-Phi distribution of tracks in tagged events", 100, -.9,.9,100,0,2*TMath::Pi());
    fHistEtaPhiTag_MB->GetXaxis()->SetTitle("Eta");
    fHistEtaPhiTag_MB->GetYaxis()->SetTitle("Phi");
    
    fHistEtaPhiTag_EMC7 = new TH2F("fHistEtaPhiTag_EMC7", "Eta-Phi distribution of tracks in tagged events", 100, -.9,.9,100,0,2*TMath::Pi());
    fHistEtaPhiTag_EMC7->GetXaxis()->SetTitle("Eta");
    fHistEtaPhiTag_EMC7->GetYaxis()->SetTitle("Phi");
    
    fHistEtaPhiTag_EMC8 = new TH2F("fHistEtaPhiTag_EMC8", "Eta-Phi distribution of tracks in tagged events", 100, -.9,.9,100,0,2*TMath::Pi());
    fHistEtaPhiTag_EMC8->GetXaxis()->SetTitle("Eta");
    fHistEtaPhiTag_EMC8->GetYaxis()->SetTitle("Phi");
    
    fHistEtaPhiTag_EMCJet = new TH2F("fHistEtaPhiTag_EMCJet", "Eta-Phi distribution of tracks in tagged events", 100, -.9,.9,100,0,2*TMath::Pi());
    fHistEtaPhiTag_EMCJet->GetXaxis()->SetTitle("Eta");
    fHistEtaPhiTag_EMCJet->GetYaxis()->SetTitle("Phi");
    
    // Eta-Phi distribution
    fHistEtaPhi_MB = new TH2F("fHistEtaPhi_MB", "Eta-Phi distribution of tracks", 100, -.9,.9,100,0,2*TMath::Pi());
    fHistEtaPhi_MB->GetXaxis()->SetTitle("Eta");
    fHistEtaPhi_MB->GetYaxis()->SetTitle("Phi");
    
    fHistEtaPhi_EMC7 = new TH2F("fHistEtaPhi_EMC7", "Eta-Phi distribution of tracks", 100, -.9,.9,100,0,2*TMath::Pi());
    fHistEtaPhi_EMC7->GetXaxis()->SetTitle("Eta");
    fHistEtaPhi_EMC7->GetYaxis()->SetTitle("Phi");
    
    fHistEtaPhi_EMC8 = new TH2F("fHistEtaPhi_EMC8", "Eta-Phi distribution of tracks", 100, -.9,.9,100,0,2*TMath::Pi());
    fHistEtaPhi_EMC8->GetXaxis()->SetTitle("Eta");
    fHistEtaPhi_EMC8->GetYaxis()->SetTitle("Phi");
    
    fHistEtaPhi_EMCJet = new TH2F("fHistEtaPhi_EMCJet", "Eta-Phi distribution of tracks", 100, -.9,.9,100,0,2*TMath::Pi());
    fHistEtaPhi_EMCJet->GetXaxis()->SetTitle("Eta");
    fHistEtaPhi_EMCJet->GetYaxis()->SetTitle("Phi");
    
    fHistEtaPhiTPCOnly_MB = new TH2F("fHistEtaPhiTPCOnly_MB", "Eta-Phi distribution of TPC only tracks", 100, -.9,.9,100,0,2*TMath::Pi());
    fHistEtaPhiTPCOnly_MB->GetXaxis()->SetTitle("Eta");
    fHistEtaPhiTPCOnly_MB->GetYaxis()->SetTitle("Phi");
    
    // Energy per event histos
    fHistPtSum_MB = new TH1F("fHistPtSum_MB", "Pt sum for events w/o an electron candidate", 500, 0, 500);
    fHistPtSum_MB->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSum_MB->GetYaxis()->SetTitle("Cts");
    
    fHistPtSum_EMC7 = new TH1F("fHistPtSum_EMC7", "Pt sum for events w/o an electron candidate", 500, 0, 500);
    fHistPtSum_EMC7->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSum_EMC7->GetYaxis()->SetTitle("Cts");
    
    fHistPtSum_EMC8 = new TH1F("fHistPtSum_EMC8", "Pt sum for events w/o an electron candidate", 500, 0, 500);
    fHistPtSum_EMC8->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSum_EMC8->GetYaxis()->SetTitle("Cts");
    
    fHistPtSum_EMCJet = new TH1F("fHistPtSum_EMCJet", "Pt sum for events w/o an electron candidate", 500, 0, 500);
    fHistPtSum_EMCJet->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSum_EMCJet->GetYaxis()->SetTitle("Cts");
    
    // Energy per tagged event histos
    fHistPtSumTag_MB = new TH1F("fHistPtSumTag_MB", "Pt sum for events w/ an electron candidate", 500, 0, 500);
    fHistPtSumTag_MB->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSumTag_MB->GetYaxis()->SetTitle("Cts");
    
    fHistPtSumTag_EMC7 = new TH1F("fHistPtSumTag_EMC7", "Pt sum for events w/ an electron candidate", 500, 0, 500);
    fHistPtSumTag_EMC7->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSumTag_EMC7->GetYaxis()->SetTitle("Cts");
    
    fHistPtSumTag_EMC8 = new TH1F("fHistPtSumTag_EMC8", "Pt sum for events w/ an electron candidate", 500, 0, 500);
    fHistPtSumTag_EMC8->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSumTag_EMC8->GetYaxis()->SetTitle("Cts");
    
    fHistPtSumTag_EMCJet = new TH1F("fHistPtSumTag_EMCJet", "Pt sum for events w/ an electron candidate", 500, 0, 500);
    fHistPtSumTag_EMCJet->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSumTag_EMCJet->GetYaxis()->SetTitle("Cts");
    
    // Numbers of events
    fHistNevents_MB = new TH1F("fHistNevents_MB", "Number of events that have an 'electron'", 2,0,1);
    fHistNevents_MB->GetXaxis()->SetBinLabel(1,"Events");
    fHistNevents_MB->GetXaxis()->SetBinLabel(2,"Events containing candidates");
    fHistNevents_MB->GetYaxis()->SetTitle("Cts");
    
    fHistNevents_EMC7 = new TH1F("fHistNevents_EMC7", "Number of events that have an 'electron'", 2,0,1);
    fHistNevents_EMC7->GetXaxis()->SetBinLabel(1,"Events");
    fHistNevents_EMC7->GetXaxis()->SetBinLabel(2,"Events containing candidates");
    fHistNevents_EMC7->GetYaxis()->SetTitle("Cts");
    
    fHistNevents_EMC8 = new TH1F("fHistNevents_EMC8", "Number of events that have an 'electron'", 2,0,1);
    fHistNevents_EMC8->GetXaxis()->SetBinLabel(1,"Events");
    fHistNevents_EMC8->GetXaxis()->SetBinLabel(2,"Events containing candidates");
    fHistNevents_EMC8->GetYaxis()->SetTitle("Cts");
    
    fHistNevents_EMCJet = new TH1F("fHistNevents_EMCJet", "Number of events that have an 'electron'", 2,0,1);
    fHistNevents_EMCJet->GetXaxis()->SetBinLabel(1,"Events");
    fHistNevents_EMCJet->GetXaxis()->SetBinLabel(2,"Events containing candidates");
    fHistNevents_EMCJet->GetYaxis()->SetTitle("Cts");
    
    //TPC Sigma
    fHistTPCSig_MB = new TH1F("fHistTPCSig_MB", "dEdx Resolution in TPC", 100,0,100);
    fHistTPCSig_MB->GetXaxis()->SetTitle("Resolution");
    fHistTPCSig_MB->GetYaxis()->SetTitle("Count");
    
    fHistTPCSig_EMC7 = new TH1F("fHistTPCSig_EMC7", "dEdx Resolution in TPC", 100,0,100);
    fHistTPCSig_EMC7->GetXaxis()->SetTitle("Resolution");
    fHistTPCSig_EMC7->GetYaxis()->SetTitle("Count");
    
    fHistTPCSig_EMC8 = new TH1F("fHistTPCSig_EMC8", "dEdx Resolution in TPC", 100,0,100);
    fHistTPCSig_EMC8->GetXaxis()->SetTitle("Resolution");
    fHistTPCSig_EMC8->GetYaxis()->SetTitle("Count");
    
    fHistTPCSig_EMCJet = new TH1F("fHistTPCSig_EMCJet", "dEdx Resolution in TPC", 100,0,100);
    fHistTPCSig_EMCJet->GetXaxis()->SetTitle("Resolution");
    fHistTPCSig_EMCJet->GetYaxis()->SetTitle("Count");
    
    //TPC Sigma after general cuts
    fHistTPCSigCut_MB = new TH1F("fHistTPCSigCut_MB", "dEdx Resolution in TPC for tracks that pass the basic track cuts", 100,0,100);
    fHistTPCSigCut_MB->GetXaxis()->SetTitle("Resolution");
    fHistTPCSigCut_MB->GetYaxis()->SetTitle("Count");
    
    fHistTPCSigCut_EMC7 = new TH1F("fHistTPCSigCut_EMC7", "dEdx Resolution in TPC for tracks that pass the basic track cuts", 100,0,100);
    fHistTPCSigCut_EMC7->GetXaxis()->SetTitle("Resolution");
    fHistTPCSigCut_EMC7->GetYaxis()->SetTitle("Count");
    
    fHistTPCSigCut_EMC8 = new TH1F("fHistTPCSigCut_EMC8", "dEdx Resolution in TPC for tracks that pass the basic track cuts", 100,0,100);
    fHistTPCSigCut_EMC8->GetXaxis()->SetTitle("Resolution");
    fHistTPCSigCut_EMC8->GetYaxis()->SetTitle("Count");
    
    fHistTPCSigCut_EMCJet = new TH1F("fHistTPCSigCut_EMCJet", "dEdx Resolution in TPC for tracks that pass the basic track cuts", 100,0,100);
    fHistTPCSigCut_EMCJet->GetXaxis()->SetTitle("Resolution");
    fHistTPCSigCut_EMCJet->GetYaxis()->SetTitle("Count");
    
    //Impact Parameter histos
    fHistImpPar_MB = new TH1F("fHistImpPar_MB", "Impact Parameter distribution in xy plane for all tracks", 100,-.5, .5);
    fHistImpPar_MB->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpPar_MB->GetYaxis()->SetTitle("Count");
    
    fHistImpPar_EMC7 = new TH1F("fHistImpPar_EMC7", "Impact Parameter distribution in xy plane for all tracks", 100,-.5, .5);
    fHistImpPar_EMC7->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpPar_EMC7->GetYaxis()->SetTitle("Count");
    
    fHistImpPar_EMC8 = new TH1F("fHistImpPar_EMC8", "Impact Parameter distribution in xy plane for all tracks", 100,-.5, .5);
    fHistImpPar_EMC8->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpPar_EMC8->GetYaxis()->SetTitle("Count");
    
    fHistImpPar_EMCJet = new TH1F("fHistImpPar_EMCJet", "Impact Parameter distribution in xy plane for all tracks", 100,-.5, .5);
    fHistImpPar_EMCJet->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpPar_EMCJet->GetYaxis()->SetTitle("Count");
    
    //Impact Parameter for tagged electrons histos
    fHistImpParTag_MB = new TH1F("fHistImpParTag_MB", "Impact Parameter distribution in xy plane for electron candidates", 100,-.5, .5);
    fHistImpParTag_MB->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpParTag_MB->GetYaxis()->SetTitle("Count");
    
    fHistImpParTag_EMC7 = new TH1F("fHistImpParTag_EMC7", "Impact Parameter distribution in xy plane for electron candidates", 100,-.5, .5);
    fHistImpParTag_EMC7->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpParTag_EMC7->GetYaxis()->SetTitle("Count");
    
    fHistImpParTag_EMC8 = new TH1F("fHistImpParTag_EMC8", "Impact Parameter distribution in xy plane for electron candidates", 100,-.5, .5);
    fHistImpParTag_EMC8->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpParTag_EMC8->GetYaxis()->SetTitle("Count");
    
    fHistImpParTag_EMCJet = new TH1F("fHistImpParTag_EMCJet", "Impact Parameter distribution in xy plane for electron candidates", 100,-.5, .5);
    fHistImpParTag_EMCJet->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpParTag_EMCJet->GetYaxis()->SetTitle("Count");
    
    //Number of clusters per track in TPC
    fHistTPCNClus_MB = new TH1F("fHistTPCNClus_MB", "Number of Clusters per track in TPC", 159, 0, 159);
    fHistTPCNClus_MB->GetXaxis()->SetTitle("Number of TPC Clusters");
    fHistTPCNClus_MB->GetYaxis()->SetTitle("Number of Tracks");
    
    fHistTPCNClus_EMC7 = new TH1F("fHistTPCNClus_EMC7", "Number of Clusters per track in TPC", 159, 0, 159);
    fHistTPCNClus_EMC7->GetXaxis()->SetTitle("Number of TPC Clusters");
    fHistTPCNClus_EMC7->GetYaxis()->SetTitle("Number of Tracks");
    
    fHistTPCNClus_EMC8 = new TH1F("fHistTPCNClus_EMC8", "Number of Clusters per track in TPC", 159, 0, 159);
    fHistTPCNClus_EMC8->GetXaxis()->SetTitle("Number of TPC Clusters");
    fHistTPCNClus_EMC8->GetYaxis()->SetTitle("Number of Tracks");
    
    fHistTPCNClus_EMCJet = new TH1F("fHistTPCNClus_EMCJet", "Number of Clusters per track in TPC", 159, 0, 159);
    fHistTPCNClus_EMCJet->GetXaxis()->SetTitle("Number of TPC Clusters");
    fHistTPCNClus_EMCJet->GetYaxis()->SetTitle("Number of Tracks");
    
    //Number of clusters per track in ITS
    fHistITSNClus_MB = new TH1F("fHistITSNClus_MB", "Number of Clusters per Track in ITS", 10, 0, 10);
    fHistITSNClus_MB->GetXaxis()->SetTitle("Number of ITS Clusters");
    fHistITSNClus_MB->GetYaxis()->SetTitle("Number of Tracks");
    
    fHistITSNClus_EMC7 = new TH1F("fHistITSNClus_EMC7", "Number of Clusters per Track in ITS", 10, 0, 10);
    fHistITSNClus_EMC7->GetXaxis()->SetTitle("Number of ITS Clusters");
    fHistITSNClus_EMC7->GetYaxis()->SetTitle("Number of Tracks");
    
    fHistITSNClus_EMC8 = new TH1F("fHistITSNClus_EMC8", "Number of Clusters per Track in ITS", 10, 0, 10);
    fHistITSNClus_EMC8->GetXaxis()->SetTitle("Number of ITS Clusters");
    fHistITSNClus_EMC8->GetYaxis()->SetTitle("Number of Tracks");
    
    fHistITSNClus_EMCJet = new TH1F("fHistITSNClus_EMCJet", "Number of Clusters per Track in ITS", 10, 0, 10);
    fHistITSNClus_EMCJet->GetXaxis()->SetTitle("Number of ITS Clusters");
    fHistITSNClus_EMCJet->GetYaxis()->SetTitle("Number of Tracks");
        
    // Pt distribution of all tracks and tagged tracks
    fHistPtAssoc_MB = new TH1F("fHistPtAssoc_MB", "Pt distribution for associated tracks", 100,0, 15);
    fHistPtAssoc_MB->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtAssoc_MB->GetYaxis()->SetTitle("Count");
    
    fHistPtAssoc_EMC7 = new TH1F("fHistPtAssoc_EMC7", "Pt distribution for associated tracks", 100,0, 15);
    fHistPtAssoc_EMC7->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtAssoc_EMC7->GetYaxis()->SetTitle("Count");
    
    fHistPtAssoc_EMC8 = new TH1F("fHistPtAssoc_EMC8", "Pt distribution for associated tracks", 100,0, 15);
    fHistPtAssoc_EMC8->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtAssoc_EMC8->GetYaxis()->SetTitle("Count");
    
    fHistPtAssoc_EMCJet = new TH1F("fHistPtAssoc_EMCJet", "Pt distribution for associated tracks", 100,0, 15);
    fHistPtAssoc_EMCJet->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtAssoc_EMCJet->GetYaxis()->SetTitle("Count");
    
    //Impact Parameter for tagged electrons histos
    fHistPtTag_MB = new TH1F("fHistPtTag_MB", "Pt distribution for electron candidates", 100,0, 15);
    fHistPtTag_MB->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtTag_MB->GetYaxis()->SetTitle("Count");
    
    fHistPtTag_EMC7 = new TH1F("fHistPtTag_EMC7", "Pt distribution for electron candidates", 100,0, 15);
    fHistPtTag_EMC7->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtTag_EMC7->GetYaxis()->SetTitle("Count");
    
    fHistPtTag_EMC8 = new TH1F("fHistPtTag_EMC8", "Pt distribution for electron candidates", 100,0, 15);
    fHistPtTag_EMC8->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtTag_EMC8->GetYaxis()->SetTitle("Count");
    
    fHistPtTag_EMCJet = new TH1F("fHistPtTag_EMCJet", "Pt distribution for electron candidates", 100,0, 15);
    fHistPtTag_EMCJet->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtTag_EMCJet->GetYaxis()->SetTitle("Count");
        
    //Add rejection plots to MB plots since it is the easiest place
    fOutputMB->Add(fHistPIDRejection);
    fOutputMB->Add(fHistBadEMCclusID);
    fOutputMB->Add(fHistNElecPerEvent);
    fOutputMB->Add(fHistPhotoMismatch);
    //ditto for the pt sum plots
    fOutputMB->Add(fHistPtSumTransMaxB2B);
    fOutputMB->Add(fHistPtSumTransMinB2B);
    fOutputMB->Add(fHistPtSumTransMaxLead);
    fOutputMB->Add(fHistPtSumTransMinLead);
    
    fOutputMB->Add(fHistPtAssoc_MB);
    fOutputMB->Add(fHistPtTag_MB);
    fOutputMB->Add(fHistTPCNClus_MB);
    fOutputMB->Add(fHistITSNClus_MB);
    fOutputMB->Add(fHistTPCSig_MB);
    fOutputMB->Add(fHistTPCSigCut_MB);
    fOutputMB->Add(fHistImpPar_MB);
    fOutputMB->Add(fHistImpParTag_MB);
    fOutputMB->Add(fHistNevents_MB);
    fOutputMB->Add(fHistPtSum_MB);
    fOutputMB->Add(fHistPtSumTag_MB);
    fOutputMB->Add(fHistEtaPhi_MB);
    fOutputMB->Add(fHistEtaPhiTag_MB);
    fOutputMB->Add(fHistEtaPhiTPCOnly_MB);
    fOutputMB->Add(fHistEMC_Had_MB_1Gev);
    fOutputMB->Add(fHistInvMassElecLike_MB);
    fOutputMB->Add(fHistInvMassElecUnLike_MB);
    fOutputMB->Add(fHistOpAngElecLike_MB);
    fOutputMB->Add(fHistOpAngElecUnLike_MB);
    for(Int_t i=0; i<6;i++){
        fOutputMB->Add(fHistTPC_TOF_MB[i]);
        fOutputMB->Add(fHistTPC_EMC_MB[i]);
        fOutputMB->Add(fHistTPC_TRD_MB[i]);
        fOutputMB->Add(fHistTPC_TOFEMC_MB[i]);
        fOutputMB->Add(fHistTPC_TOFTRD_MB[i]);
        fOutputMB->Add(fHistTPC_EMCTRD_MB[i]);
        fOutputMB->Add(fHistTPC_TOFEMCTRD_MB[i]);
        
        fOutputMB->Add(fHistTOF_TPC_MB[i]);
        fOutputMB->Add(fHistTOF_EMC_MB[i]);
        fOutputMB->Add(fHistTOF_TRD_MB[i]);
        fOutputMB->Add(fHistTOF_TPCEMC_MB[i]);
        fOutputMB->Add(fHistTOF_TPCTRD_MB[i]);
        fOutputMB->Add(fHistTOF_EMCTRD_MB[i]);
        fOutputMB->Add(fHistTOF_TPCEMCTRD_MB[i]);
        
        fOutputMB->Add(fHistEMC_TPC_MB[i]);
        fOutputMB->Add(fHistEMC_TOF_MB[i]);
        fOutputMB->Add(fHistEMC_TRD_MB[i]);
        fOutputMB->Add(fHistEMC_TPCTOF_MB[i]);
        fOutputMB->Add(fHistEMC_TPCTRD_MB[i]);
        fOutputMB->Add(fHistEMC_TOFTRD_MB[i]);
        fOutputMB->Add(fHistEMC_TPCTOFTRD_MB[i]);
        
        fOutputMB->Add(fHistTRD_TPC_MB[i]);
        fOutputMB->Add(fHistTRD_TOF_MB[i]);
        fOutputMB->Add(fHistTRD_EMC_MB[i]);
        fOutputMB->Add(fHistTRD_TPCTOF_MB[i]);
        fOutputMB->Add(fHistTRD_TPCEMC_MB[i]);
        fOutputMB->Add(fHistTRD_TOFEMC_MB[i]);
        fOutputMB->Add(fHistTRD_TPCTOFEMC_MB[i]);
        
        fOutputMB->Add(fHistM02_All_MB[i]);
        fOutputMB->Add(fHistM20_All_MB[i]);
        fOutputMB->Add(fHistM02_Elec_MB[i]);
        fOutputMB->Add(fHistM20_Elec_MB[i]);
    }
    for(Int_t i=0; i<3;i++){
    fOutputMB->Add(fHistDPhi300_500_MB[i]);
    fOutputMB->Add(fHistDPhi500_800_MB[i]);
    fOutputMB->Add(fHistDPhi800_1_MB[i]);
    fOutputMB->Add(fHistDPhi1_2_MB[i]);
    fOutputMB->Add(fHistDPhi2_3_MB[i]);
    fOutputMB->Add(fHistDPhi3_4_MB[i]);
    fOutputMB->Add(fHistDPhi4_MB[i]);
    }
    fOutputMB->Add(fHistDPhi28_MB);
    fOutputMB->Add(fHistDPhiDEta28_MB);
    
    fOutputEMC7->Add(fHistPtAssoc_EMC7);
    fOutputEMC7->Add(fHistPtTag_EMC7);
    fOutputEMC7->Add(fHistTPCNClus_EMC7);
    fOutputEMC7->Add(fHistITSNClus_EMC7);
    fOutputEMC7->Add(fHistTPCSig_EMC7);
    fOutputEMC7->Add(fHistTPCSigCut_EMC7);
    fOutputEMC7->Add(fHistImpPar_EMC7);
    fOutputEMC7->Add(fHistImpParTag_EMC7);
    fOutputEMC7->Add(fHistNevents_EMC7);
    fOutputEMC7->Add(fHistPtSum_EMC7);
    fOutputEMC7->Add(fHistPtSumTag_EMC7);
    fOutputEMC7->Add(fHistEtaPhi_EMC7);
    fOutputEMC7->Add(fHistEtaPhiTag_EMC7);
    fOutputEMC7->Add(fHistInvMassElecLike_EMC7);
    fOutputEMC7->Add(fHistInvMassElecUnLike_EMC7);
    fOutputEMC7->Add(fHistOpAngElecLike_EMC7);
    fOutputEMC7->Add(fHistOpAngElecUnLike_EMC7);
    for(Int_t i=0; i<6;i++){
        fOutputEMC7->Add(fHistTPC_TOF_EMC7[i]);
        fOutputEMC7->Add(fHistTPC_EMC_EMC7[i]);
        fOutputEMC7->Add(fHistTPC_TRD_EMC7[i]);
        fOutputEMC7->Add(fHistTPC_TOFEMC_EMC7[i]);
        fOutputEMC7->Add(fHistTPC_TOFTRD_EMC7[i]);
        fOutputEMC7->Add(fHistTPC_EMCTRD_EMC7[i]);
        fOutputEMC7->Add(fHistTPC_TOFEMCTRD_EMC7[i]);
        
        fOutputEMC7->Add(fHistTOF_TPC_EMC7[i]);
        fOutputEMC7->Add(fHistTOF_EMC_EMC7[i]);
        fOutputEMC7->Add(fHistTOF_TRD_EMC7[i]);
        fOutputEMC7->Add(fHistTOF_TPCEMC_EMC7[i]);
        fOutputEMC7->Add(fHistTOF_TPCTRD_EMC7[i]);
        fOutputEMC7->Add(fHistTOF_EMCTRD_EMC7[i]);
        fOutputEMC7->Add(fHistTOF_TPCEMCTRD_EMC7[i]);
        
        fOutputEMC7->Add(fHistEMC_TPC_EMC7[i]);
        fOutputEMC7->Add(fHistEMC_TOF_EMC7[i]);
        fOutputEMC7->Add(fHistEMC_TRD_EMC7[i]);
        fOutputEMC7->Add(fHistEMC_TPCTOF_EMC7[i]);
        fOutputEMC7->Add(fHistEMC_TPCTRD_EMC7[i]);
        fOutputEMC7->Add(fHistEMC_TOFTRD_EMC7[i]);
        fOutputEMC7->Add(fHistEMC_TPCTOFTRD_EMC7[i]);
        
        fOutputEMC7->Add(fHistTRD_TPC_EMC7[i]);
        fOutputEMC7->Add(fHistTRD_TOF_EMC7[i]);
        fOutputEMC7->Add(fHistTRD_EMC_EMC7[i]);
        fOutputEMC7->Add(fHistTRD_TPCTOF_EMC7[i]);
        fOutputEMC7->Add(fHistTRD_TPCEMC_EMC7[i]);
        fOutputEMC7->Add(fHistTRD_TOFEMC_EMC7[i]);
        fOutputEMC7->Add(fHistTRD_TPCTOFEMC_EMC7[i]);
        
        fOutputEMC7->Add(fHistM02_All_EMC7[i]);
        fOutputEMC7->Add(fHistM20_All_EMC7[i]);
        fOutputEMC7->Add(fHistM02_Elec_EMC7[i]);
        fOutputEMC7->Add(fHistM20_Elec_EMC7[i]);
    }
    for(Int_t i=0; i<3;i++){
    fOutputEMC7->Add(fHistDPhi300_500_EMC7[i]);
    fOutputEMC7->Add(fHistDPhi500_800_EMC7[i]);
    fOutputEMC7->Add(fHistDPhi800_1_EMC7[i]);
    fOutputEMC7->Add(fHistDPhi1_2_EMC7[i]);
    fOutputEMC7->Add(fHistDPhi2_3_EMC7[i]);
    fOutputEMC7->Add(fHistDPhi3_4_EMC7[i]);
    fOutputEMC7->Add(fHistDPhi4_EMC7[i]);
    }
    fOutputEMC7->Add(fHistDPhi28_EMC7);
    fOutputEMC7->Add(fHistDPhiDEta28_EMC7);
    
    fOutputEMC8->Add(fHistPtAssoc_EMC8);
    fOutputEMC8->Add(fHistPtTag_EMC8);
    fOutputEMC8->Add(fHistTPCNClus_EMC8);
    fOutputEMC8->Add(fHistITSNClus_EMC8);
    fOutputEMC8->Add(fHistTPCSig_EMC8);
    fOutputEMC8->Add(fHistTPCSigCut_EMC8);
    fOutputEMC8->Add(fHistImpPar_EMC8);
    fOutputEMC8->Add(fHistImpParTag_EMC8);
    fOutputEMC8->Add(fHistNevents_EMC8);
    fOutputEMC8->Add(fHistPtSum_EMC8);
    fOutputEMC8->Add(fHistPtSumTag_EMC8);
    fOutputEMC8->Add(fHistEtaPhi_EMC8);
    fOutputEMC8->Add(fHistEtaPhiTag_EMC8);
    fOutputEMC8->Add(fHistInvMassElecLike_EMC8);
    fOutputEMC8->Add(fHistInvMassElecUnLike_EMC8);
    fOutputEMC8->Add(fHistOpAngElecLike_EMC8);
    fOutputEMC8->Add(fHistOpAngElecUnLike_EMC8);
    for(Int_t i=0; i<6;i++){
        fOutputEMC8->Add(fHistTPC_TOF_EMC8[i]);
        fOutputEMC8->Add(fHistTPC_EMC_EMC8[i]);
        fOutputEMC8->Add(fHistTPC_TRD_EMC8[i]);
        fOutputEMC8->Add(fHistTPC_TOFEMC_EMC8[i]);
        fOutputEMC8->Add(fHistTPC_TOFTRD_EMC8[i]);
        fOutputEMC8->Add(fHistTPC_EMCTRD_EMC8[i]);
        fOutputEMC8->Add(fHistTPC_TOFEMCTRD_EMC8[i]);
        
        fOutputEMC8->Add(fHistTOF_TPC_EMC8[i]);
        fOutputEMC8->Add(fHistTOF_EMC_EMC8[i]);
        fOutputEMC8->Add(fHistTOF_TRD_EMC8[i]);
        fOutputEMC8->Add(fHistTOF_TPCEMC_EMC8[i]);
        fOutputEMC8->Add(fHistTOF_TPCTRD_EMC8[i]);
        fOutputEMC8->Add(fHistTOF_EMCTRD_EMC8[i]);
        fOutputEMC8->Add(fHistTOF_TPCEMCTRD_EMC8[i]);
        
        fOutputEMC8->Add(fHistEMC_TPC_EMC8[i]);
        fOutputEMC8->Add(fHistEMC_TOF_EMC8[i]);
        fOutputEMC8->Add(fHistEMC_TRD_EMC8[i]);
        fOutputEMC8->Add(fHistEMC_TPCTOF_EMC8[i]);
        fOutputEMC8->Add(fHistEMC_TPCTRD_EMC8[i]);
        fOutputEMC8->Add(fHistEMC_TOFTRD_EMC8[i]);
        fOutputEMC8->Add(fHistEMC_TPCTOFTRD_EMC8[i]);
        
        fOutputEMC8->Add(fHistTRD_TPC_EMC8[i]);
        fOutputEMC8->Add(fHistTRD_TOF_EMC8[i]);
        fOutputEMC8->Add(fHistTRD_EMC_EMC8[i]);
        fOutputEMC8->Add(fHistTRD_TPCTOF_EMC8[i]);
        fOutputEMC8->Add(fHistTRD_TPCEMC_EMC8[i]);
        fOutputEMC8->Add(fHistTRD_TOFEMC_EMC8[i]);
        fOutputEMC8->Add(fHistTRD_TPCTOFEMC_EMC8[i]);
        
        fOutputEMC8->Add(fHistM02_All_EMC8[i]);
        fOutputEMC8->Add(fHistM20_All_EMC8[i]);
        fOutputEMC8->Add(fHistM02_Elec_EMC8[i]);
        fOutputEMC8->Add(fHistM20_Elec_EMC8[i]);
    }
    for(Int_t i=0; i<3;i++){
    fOutputEMC8->Add(fHistDPhi300_500_EMC8[i]);
    fOutputEMC8->Add(fHistDPhi500_800_EMC8[i]);
    fOutputEMC8->Add(fHistDPhi800_1_EMC8[i]);
    fOutputEMC8->Add(fHistDPhi1_2_EMC8[i]);
    fOutputEMC8->Add(fHistDPhi2_3_EMC8[i]);
    fOutputEMC8->Add(fHistDPhi3_4_EMC8[i]);
    fOutputEMC8->Add(fHistDPhi4_EMC8[i]);
    }
    fOutputEMC8->Add(fHistDPhi28_EMC8);
    fOutputEMC8->Add(fHistDPhiDEta28_EMC8);
    
    fOutputEMCJet->Add(fHistPtAssoc_EMCJet);
    fOutputEMCJet->Add(fHistPtTag_EMCJet);
    fOutputEMCJet->Add(fHistTPCNClus_EMCJet);
    fOutputEMCJet->Add(fHistITSNClus_EMCJet);
    fOutputEMCJet->Add(fHistTPCSig_EMCJet);
    fOutputEMCJet->Add(fHistTPCSigCut_EMCJet);
    fOutputEMCJet->Add(fHistImpPar_EMCJet);
    fOutputEMCJet->Add(fHistImpParTag_EMCJet);
    fOutputEMCJet->Add(fHistNevents_EMCJet);
    fOutputEMCJet->Add(fHistPtSum_EMCJet);
    fOutputEMCJet->Add(fHistPtSumTag_EMCJet);
    fOutputEMCJet->Add(fHistEtaPhi_EMCJet);
    fOutputEMCJet->Add(fHistEtaPhiTag_EMCJet);
    fOutputEMCJet->Add(fHistInvMassElecLike_EMCJet);
    fOutputEMCJet->Add(fHistInvMassElecUnLike_EMCJet);
    fOutputEMCJet->Add(fHistOpAngElecLike_EMCJet);
    fOutputEMCJet->Add(fHistOpAngElecUnLike_EMCJet);
    for(Int_t i=0; i<6;i++){
        fOutputEMCJet->Add(fHistTPC_TOF_EMCJet[i]);
        fOutputEMCJet->Add(fHistTPC_EMC_EMCJet[i]);
        fOutputEMCJet->Add(fHistTPC_TRD_EMCJet[i]);
        fOutputEMCJet->Add(fHistTPC_TOFEMC_EMCJet[i]);
        fOutputEMCJet->Add(fHistTPC_TOFTRD_EMCJet[i]);
        fOutputEMCJet->Add(fHistTPC_EMCTRD_EMCJet[i]);
        fOutputEMCJet->Add(fHistTPC_TOFEMCTRD_EMCJet[i]);
        
        fOutputEMCJet->Add(fHistTOF_TPC_EMCJet[i]);
        fOutputEMCJet->Add(fHistTOF_EMC_EMCJet[i]);
        fOutputEMCJet->Add(fHistTOF_TRD_EMCJet[i]);
        fOutputEMCJet->Add(fHistTOF_TPCEMC_EMCJet[i]);
        fOutputEMCJet->Add(fHistTOF_TPCTRD_EMCJet[i]);
        fOutputEMCJet->Add(fHistTOF_EMCTRD_EMCJet[i]);
        fOutputEMCJet->Add(fHistTOF_TPCEMCTRD_EMCJet[i]);
        
        fOutputEMCJet->Add(fHistEMC_TPC_EMCJet[i]);
        fOutputEMCJet->Add(fHistEMC_TOF_EMCJet[i]);
        fOutputEMCJet->Add(fHistEMC_TRD_EMCJet[i]);
        fOutputEMCJet->Add(fHistEMC_TPCTOF_EMCJet[i]);
        fOutputEMCJet->Add(fHistEMC_TPCTRD_EMCJet[i]);
        fOutputEMCJet->Add(fHistEMC_TOFTRD_EMCJet[i]);
        fOutputEMCJet->Add(fHistEMC_TPCTOFTRD_EMCJet[i]);
        
        fOutputEMCJet->Add(fHistTRD_TPC_EMCJet[i]);
        fOutputEMCJet->Add(fHistTRD_TOF_EMCJet[i]);
        fOutputEMCJet->Add(fHistTRD_EMC_EMCJet[i]);
        fOutputEMCJet->Add(fHistTRD_TPCTOF_EMCJet[i]);
        fOutputEMCJet->Add(fHistTRD_TPCEMC_EMCJet[i]);
        fOutputEMCJet->Add(fHistTRD_TOFEMC_EMCJet[i]);
        fOutputEMCJet->Add(fHistTRD_TPCTOFEMC_EMCJet[i]);
        
        fOutputEMCJet->Add(fHistM02_All_EMCJet[i]);
        fOutputEMCJet->Add(fHistM20_All_EMCJet[i]);
        fOutputEMCJet->Add(fHistM02_Elec_EMCJet[i]);
        fOutputEMCJet->Add(fHistM20_Elec_EMCJet[i]);
    }
    for(Int_t i=0; i<3;i++){
    fOutputEMCJet->Add(fHistDPhi300_500_EMCJet[i]);
    fOutputEMCJet->Add(fHistDPhi500_800_EMCJet[i]);
    fOutputEMCJet->Add(fHistDPhi800_1_EMCJet[i]);
    fOutputEMCJet->Add(fHistDPhi1_2_EMCJet[i]);
    fOutputEMCJet->Add(fHistDPhi2_3_EMCJet[i]);
    fOutputEMCJet->Add(fHistDPhi3_4_EMCJet[i]);
    fOutputEMCJet->Add(fHistDPhi4_EMCJet[i]);
    }
    fOutputEMCJet->Add(fHistDPhi28_EMCJet);
    fOutputEMCJet->Add(fHistDPhiDEta28_EMCJet);
    
    //Add Region Histos
    
    for(Int_t i=0;i<4;i++){
        
        //Tag Side Histos
        
        //Multiplicity Histos
        fOutputMB->Add(fHistTrkMultTag_MB[i]);
        fOutputEMC7->Add(fHistTrkMultTag_EMC7[i]);
        fOutputEMC8->Add(fHistTrkMultTag_EMC8[i]);
        fOutputEMCJet->Add(fHistTrkMultTag_EMCJet[i]);
        
        //Pt Histos
        fOutputMB->Add(fHistTrkPtTag_MB[i]);
        fOutputEMC7->Add(fHistTrkPtTag_EMC7[i]);
        fOutputEMC8->Add(fHistTrkPtTag_EMC8[i]);
        fOutputEMCJet->Add(fHistTrkPtTag_EMCJet[i]);
        
        //DeDx by Pt Histos
        fOutputMB->Add(fHistDeDxPtTag_MB[i]);
        fOutputEMC7->Add(fHistDeDxPtTag_EMC7[i]);
        fOutputEMC8->Add(fHistDeDxPtTag_EMC8[i]);
        fOutputEMCJet->Add(fHistDeDxPtTag_EMCJet[i]);
        
        //Away Side Histos
        
        //Multiplicity Histos
        fOutputMB->Add(fHistTrkMultAway_MB[i]);
        fOutputEMC7->Add(fHistTrkMultAway_EMC7[i]);
        fOutputEMC8->Add(fHistTrkMultAway_EMC8[i]);
        fOutputEMCJet->Add(fHistTrkMultAway_EMCJet[i]);
        
        //Pt Histos
        fOutputMB->Add(fHistTrkPtAway_MB[i]);
        fOutputEMC7->Add(fHistTrkPtAway_EMC7[i]);
        fOutputEMC8->Add(fHistTrkPtAway_EMC8[i]);
        fOutputEMCJet->Add(fHistTrkPtAway_EMCJet[i]);
        
        //DeDx by Pt Histos
        fOutputMB->Add(fHistDeDxPtAway_MB[i]);
        fOutputEMC7->Add(fHistDeDxPtAway_EMC7[i]);
        fOutputEMC8->Add(fHistDeDxPtAway_EMC8[i]);
        fOutputEMCJet->Add(fHistDeDxPtAway_EMCJet[i]);
        
        //TransMin Side Histos
        
        //Multiplicity Histos
        fOutputMB->Add(fHistTrkMultTransMin_MB[i]);
        fOutputEMC7->Add(fHistTrkMultTransMin_EMC7[i]);
        fOutputEMC8->Add(fHistTrkMultTransMin_EMC8[i]);
        fOutputEMCJet->Add(fHistTrkMultTransMin_EMCJet[i]);
        
        //Pt Histos
        fOutputMB->Add(fHistTrkPtTransMin_MB[i]);
        fOutputEMC7->Add(fHistTrkPtTransMin_EMC7[i]);
        fOutputEMC8->Add(fHistTrkPtTransMin_EMC8[i]);
        fOutputEMCJet->Add(fHistTrkPtTransMin_EMCJet[i]);
        
        //DeDx by Pt Histos
        fOutputMB->Add(fHistDeDxPtTransMin_MB[i]);
        fOutputEMC7->Add(fHistDeDxPtTransMin_EMC7[i]);
        fOutputEMC8->Add(fHistDeDxPtTransMin_EMC8[i]);
        fOutputEMCJet->Add(fHistDeDxPtTransMin_EMCJet[i]);
        
        //TransMax Side Histos
        
        //Multiplicity Histos
        fOutputMB->Add(fHistTrkMultTransMax_MB[i]);
        fOutputEMC7->Add(fHistTrkMultTransMax_EMC7[i]);
        fOutputEMC8->Add(fHistTrkMultTransMax_EMC8[i]);
        fOutputEMCJet->Add(fHistTrkMultTransMax_EMCJet[i]);
        
        //Pt Histos
        fOutputMB->Add(fHistTrkPtTransMax_MB[i]);
        fOutputEMC7->Add(fHistTrkPtTransMax_EMC7[i]);
        fOutputEMC8->Add(fHistTrkPtTransMax_EMC8[i]);
        fOutputEMCJet->Add(fHistTrkPtTransMax_EMCJet[i]);
        
        //DeDx by Pt Histos
        fOutputMB->Add(fHistDeDxPtTransMax_MB[i]);
        fOutputEMC7->Add(fHistDeDxPtTransMax_EMC7[i]);
        fOutputEMC8->Add(fHistDeDxPtTransMax_EMC8[i]);
        fOutputEMCJet->Add(fHistDeDxPtTransMax_EMCJet[i]);
    }

    // NEW HISTO added to fOutput here
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMC8);
    PostData(4, fOutputEMCJet);// Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskPSHFE::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
        
    //_______________________________Major event-level stuff____________________________________
    // Create pointer to reconstructed event
    AliVEvent *event = InputEvent();
    if (!event) { Printf("ERROR: Could not retrieve event"); return; }
        
    // create pointer to event
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
    if (!esd) {
        AliError("Cannot get the ESD event");
        return;
    }  
    
    // input handler
  const AliAnalysisManager* man(AliAnalysisManager::GetAnalysisManager());
  if (NULL == man) {
    AliWarning("AliAnalysisManager is not available");
    return;
  }

  AliESDInputHandler* inputHandler(dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler()));  
  if (NULL == inputHandler) {
    AliWarning("AliESDInputHandler is not available");
    return;
  }

  UInt_t fSelectMask = inputHandler->IsEventSelected();
  Bool_t isSelected = fSelectMask & (AliVEvent::kEMC7 | AliVEvent::kEMCEJE | AliVEvent::kEMC8);
  if(!isSelected){
        AliWarning("This is not an EMCal triggered event");
  }
    
  MBtrg = fSelectMask & AliVEvent::kAnyINT;
  EMC7trg = fSelectMask & AliVEvent::kEMC7;
  EMC8trg = fSelectMask & AliVEvent::kEMC8;
  EMCJettrg = fSelectMask & AliVEvent::kEMCEJE;
    
  Int_t elecIDs[1000];
  Int_t elecCnt=0;
    
    if(!globaltrackCuts||!comptrackCuts){
        AliWarning("The hybrid track cuts are null");
        return;
    }
  AliPIDResponse* fPIDResponse = inputHandler->GetPIDResponse();
      if(!fPIDResponse){
        AliWarning("NULL PIDResponse");}
    
    //__________________________End major event stuff_____________________________
    
    
    //Fill the histogram cataloguing # of events vs. events tagged
    if(MBtrg){
        fHistNevents_MB->Fill("Events",1);
    }
    if(EMC7trg){
            fHistNevents_EMC7->Fill("Events",1);
    }
    
    if(EMC8trg){
            fHistNevents_EMC8->Fill("Events",1);
    }
    
    if(EMCJettrg){
            fHistNevents_EMCJet->Fill("Events",1);
    }
    
    //Initialize energy variable and tagging flags
    Double_t PtSum=0;
    tagStrong=kFALSE;
    Bool_t tagEvt=kFALSE;

    //Initialize the # of tracks variable and the Eta Phi arrays
    Int_t ntracks = esd->GetNumberOfTracks();
    std::vector<Double_t> Eta;
    std::vector<Double_t> Phi;
    
    // Track loop for reconstructed event
    for(Int_t i = 0; i < ntracks; i++) {
        
        tagStrong=kFALSE;
        tagPhot=kFALSE;
        
        AliESDtrack* esdtrack = esd->GetTrack(i); // pointer to reconstructed to track      
        
        if(!esdtrack) { 
            AliError(Form("ERROR: Could not retrieve esdtrack %d",i)); 
            continue; 
        }

        //Fill TPCOnly track eta-phi
        if(AliESDtrackCuts::GetStandardTPCOnlyTrackCuts()->AcceptTrack(esdtrack)){
            fHistEtaPhiTPCOnly_MB->Fill(esdtrack->Eta(),esdtrack->Phi());
        }
        
        //Do hybrid track cuts
        if(!globaltrackCuts->AcceptTrack(esdtrack)&&!comptrackCuts->AcceptTrack(esdtrack)){continue;}
        
        //Add this tracks energy to the running total
        PtSum=PtSum+esdtrack->Pt();
        
        //Fill the Eta Phi arrays with this tracks Eta and Phi
        Eta.push_back(esdtrack->Eta());
        Phi.push_back(esdtrack->Phi());
        
        if(MBtrg){
            fHistEtaPhi_MB->Fill(esdtrack->Eta(),esdtrack->Phi());
        }
        if(EMC7trg){
                fHistEtaPhi_EMC7->Fill(esdtrack->Eta(),esdtrack->Phi());
        }
        if(EMC8trg){
                fHistEtaPhi_EMC8->Fill(esdtrack->Eta(),esdtrack->Phi());
        }
        if(EMCJettrg){
                fHistEtaPhi_EMCJet->Fill(esdtrack->Eta(),esdtrack->Phi());
        }
        
        //do Cut level histograms
        if(MBtrg){
            if(esdtrack->GetTPCncls()>0){
            fHistTPCNClus_MB->Fill(esdtrack->GetTPCncls());}
            fHistITSNClus_MB->Fill(esdtrack->GetNcls(0));
        }
        if(EMC7trg){
                if(esdtrack->GetTPCncls()>0){
                    fHistTPCNClus_EMC7->Fill(esdtrack->GetTPCncls());}
                fHistITSNClus_EMC7->Fill(esdtrack->GetNcls(0));
        }
        
        if(EMC8trg){
                if(esdtrack->GetTPCncls()>0){
                    fHistTPCNClus_EMC8->Fill(esdtrack->GetTPCncls());}
                fHistITSNClus_EMC8->Fill(esdtrack->GetNcls(0));
        }
        
        if(EMCJettrg){
                if(esdtrack->GetTPCncls()>0){
                    fHistTPCNClus_EMCJet->Fill(esdtrack->GetTPCncls());}
                fHistITSNClus_EMCJet->Fill(esdtrack->GetNcls(0)); 
        }
        
        //Fill histogram for TPC resolution
        if(MBtrg){
            fHistTPCSig_MB->Fill(esdtrack->GetTPCsignalSigma());
        }
        if(EMC7trg){
                fHistTPCSig_EMC7->Fill(esdtrack->GetTPCsignalSigma());
        }
        
        if(EMC8trg){
                fHistTPCSig_EMC8->Fill(esdtrack->GetTPCsignalSigma());
        }
        
        if(EMCJettrg){
                fHistTPCSig_EMCJet->Fill(esdtrack->GetTPCsignalSigma());
        }
        
        //Impact parameter
        Float_t xy;
        Float_t z;
        
        esdtrack->GetImpactParameters(xy, z);
        
        if(MBtrg){
            fHistImpPar_MB->Fill(xy);
        }
        if(EMC7trg){
                fHistImpPar_EMC7->Fill(xy);
        }
        
        if(EMC8trg){
                fHistImpPar_EMC8->Fill(xy);
        }
        
        if(EMCJettrg){
                fHistImpPar_EMCJet->Fill(xy);
        }
        
        FillPhotoElecHistos(esd, esdtrack, fPIDResponse, i);
        
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //If the track doesn't pass the cuts, move on to the next one
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if(trackCutsStrong){
            if(!fTrackCutsStrong->AcceptTrack(esdtrack) || esdtrack->GetTPCsignalN()<80){continue;}
        }else{
            if(!fTrackCutsWeak->AcceptTrack(esdtrack) || esdtrack->GetTPCsignalN()<80){continue;}
        }
        
        //Fill histogram for TPC resolution
        if(MBtrg){
            fHistTPCSigCut_MB->Fill(esdtrack->GetTPCsignalSigma());
        }
        if(EMC7trg){
                fHistTPCSigCut_EMC7->Fill(esdtrack->GetTPCsignalSigma());
        }
        
        if(EMC8trg){
                fHistTPCSigCut_EMC8->Fill(esdtrack->GetTPCsignalSigma());
        }
        
        if(EMCJettrg){
                fHistTPCSigCut_EMCJet->Fill(esdtrack->GetTPCsignalSigma());
        }
        
        FillPIDHistos(esd, esdtrack, fPIDResponse);//Fill PID histos and set "tagStrong" boolean if this track satisfies cuts
        
        //If the track made it through any 3 detector cuts
        if(tagStrong){
            tagEvt=kTRUE;//This event contains a candidate electron
            //Increment the candidate count and put track ID of candidate into array
            
            elecIDs[elecCnt]=i;
            elecCnt+=1;
            
            //Fill impact parameter plots
            if(MBtrg){
                fHistImpParTag_MB->Fill(xy);
            }
            if(EMC7trg){
                    fHistImpParTag_EMC7->Fill(xy);
            }
        
            if(EMC8trg){
                    fHistImpParTag_EMC8->Fill(xy);
            }
        
            if(EMCJettrg){
                    fHistImpParTag_EMCJet->Fill(xy);
            }
            
            //Pt distribution
            if(MBtrg){
                fHistPtTag_MB->Fill(esdtrack->Pt());
            }
            if(EMC7trg){
                    fHistPtTag_EMC7->Fill(esdtrack->Pt());
            }

            if(EMC8trg){
                    fHistPtTag_EMC8->Fill(esdtrack->Pt());
            }

            if(EMCJettrg){
                    fHistPtTag_EMCJet->Fill(esdtrack->Pt());
            }
            
            FillDPhiHistos(esd, esdtrack, i);//Fill DPhi histos
            
            if(tagPhot){fHistPhotoMismatch->Fill(1);}
            
        }//end if(tagStrong)
        
    }//end main track loop
    
    //Call function to fill Region histos and pass it int array of IDs for identified electron tracks
    Int_t elecIDsSparse[elecCnt];
    for(Int_t i=0;i<elecCnt;i++){
        elecIDsSparse[i]=elecIDs[i];
    }
    
    FillRegionHistos(esd, elecIDsSparse, elecCnt);
    
    //Fill Number of electrons plot
    fHistNElecPerEvent->Fill(elecCnt);
    
    //Fill the total pt sum histogram
    if(tagEvt){
        if(MBtrg){
            fHistPtSumTag_MB->Fill(PtSum);
        }
        if(EMC7trg){
            fHistPtSumTag_EMC7->Fill(PtSum);
        }
        if(EMC8trg){
            fHistPtSumTag_EMC8->Fill(PtSum);
        }
        if(EMCJettrg){
            fHistPtSumTag_EMCJet->Fill(PtSum);
        }
    }else{
        if(MBtrg){
                fHistPtSum_MB->Fill(PtSum);
        }
        if(EMC7trg){
                fHistPtSum_EMC7->Fill(PtSum);
        }
        if(EMC8trg){
                fHistPtSum_EMC8->Fill(PtSum);
        }
        if(EMCJettrg){
                fHistPtSum_EMCJet->Fill(PtSum);
        }
    }
    
    //Fill Nevent histos
    if(tagEvt){
        if(MBtrg){
            fHistNevents_MB->Fill("Events containing candidates",1);
        }  
                if(EMC7trg){
                        fHistNevents_EMC7->Fill("Events containing candidates",1);
                }
                if(EMC8trg){
                        fHistNevents_EMC8->Fill("Events containing candidates",1);
                }
                if(EMCJettrg){
                        fHistNevents_EMCJet->Fill("Events containing candidates",1);
                }  
    }
    
    //Fill the Eta Phi histograms
if(tagEvt){
    for(Int_t i=0;i<Eta.size();i++){
        if(MBtrg){
            fHistEtaPhiTag_MB->Fill(Eta[i],Phi[i]);
        }
        if(EMC7trg){
                fHistEtaPhiTag_EMC7->Fill(Eta[i],Phi[i]);
        }
        if(EMC8trg){
                fHistEtaPhiTag_EMC8->Fill(Eta[i],Phi[i]);
        }
        if(EMCJettrg){
                fHistEtaPhiTag_EMCJet->Fill(Eta[i],Phi[i]);
        }
    }
}
    
    
    
    // NEW HISTO should be filled before this point, as PostData puts the
    // information for this iteration of the UserExec in the container
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMC8);
    PostData(4, fOutputEMCJet);
}

//________________________________________________________________________
void AliAnalysisTaskPSHFE::FillRegionHistos(AliESDEvent *esd, Int_t *elecIDs, Int_t elecCnt)
{
    //First check if this event has any candidate electrons
    if(elecCnt==0){
        return;
    }
    if(!esd){
        AliError("The esd event passed to FillRegionHistos is null...how did this even happen");
        return;
    }
    if(!elecIDs){
        AliError("The list electron candidate track IDs is null");
        return;
    }
    
    Int_t ntracks = esd->GetNumberOfTracks();
    
    Double_t EtaMax = .8;
    Double_t PtSumLeft = 0;
    Double_t EtSumLeft = 0;
    Double_t PtSumRight = 0;
    Double_t EtSumRight = 0;
    
    for(Int_t i=0;i<elecCnt;i++){
        AliESDtrack *elecTrk = esd->GetTrack(elecIDs[i]);
        
        if(!elecTrk){
            AliWarning("The candidate electron track in null, something is wrong.");
            continue;
        }
        
        //Initialize variables for Multiplicity
        Int_t tag_Mult = 0;
        Int_t away_Mult = 0;
        Int_t left_Mult = 0;
        Int_t right_Mult = 0;
        
        //Fix the phi of our candidate electron
        Double_t Phi_0 = elecTrk->Phi();
        
        /*//Set the boundaries for our wedges based on the Phi_0
        Double_t awaySidePhiUpper = Phi_0-3.0*TMath::Pi()/4.0;
        Double_t awaySidePhiLower = Phi_0-5.0*TMath::Pi()/4.0;
        Double_t tagSidePhiUpper = Phi_0+TMath::Pi()/4.0;
        Double_t tagSidePhiLower = Phi_0-TMath::Pi()/4.0;
        
        //Away side Phi corrections
        if(awaySidePhiUpper>2.0*TMath::Pi()){
            awaySidePhiUpper = awaySidePhiUpper-2.0*TMath::Pi();
        }
        
        if(awaySidePhiLower>2.0*TMath::Pi()){
            awaySidePhiLower = awaySidePhiLower-2.0*TMath::Pi();
        }
        
        if(awaySidePhiUpper<0){
            awaySidePhiUpper = 2.0*TMath::Pi()+awaySidePhiUpper;
        }
        
        if(awaySidePhiLower<0){
            awaySidePhiLower = 2.0*TMath::Pi()+awaySidePhiLower;
        }
        
        //Tag side Phi Corrections
        if(tagSidePhiUpper>2.0*TMath::Pi()){
            tagSidePhiUpper = tagSidePhiUpper-2.0*TMath::Pi();
        }
        
        if(tagSidePhiLower>2.0*TMath::Pi()){
            tagSidePhiLower = tagSidePhiLower-2.0*TMath::Pi();
        }
        
        if(tagSidePhiUpper<0){
            tagSidePhiUpper = 2.0*TMath::Pi()+tagSidePhiUpper;
        }
        
        if(tagSidePhiLower<0){
            tagSidePhiLower = 2.0*TMath::Pi()+tagSidePhiLower;
        }*/
        
        //Now cycle through the rest of the tracks and fill the appropriate histograms
        for(Int_t j=0; j<ntracks; j++){
            
            if(i==j){continue;}//Don't double count our candidate
            
            AliESDtrack *esdtrack = esd->GetTrack(j);
            if(!esdtrack){continue;}
           
            //Do hybrid track cuts
            if(!globaltrackCuts->AcceptTrack(esdtrack)&&!comptrackCuts->AcceptTrack(esdtrack)){continue;}
            
            //Create dphi variable
            Double_t dphi = esdtrack->Phi()-Phi_0;
            if(dphi<0){dphi = (2*TMath::Pi()-Phi_0)+esdtrack->Phi();}
            
            //Tag side check
            if(((dphi>0&&dphi<TMath::Pi()/4)||(dphi>7*TMath::Pi()/4&&dphi<2*TMath::Pi()))&&TMath::Abs(esdtrack->Eta())<EtaMax){
                //Add it to the multiplicity count
                tag_Mult = tag_Mult + 1;
                
                //Fill the Pt histos with the tracks Pt
                
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    
                    if(MBtrg){
                        fHistTrkPtTag_MB[0]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTag_EMC7[0]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTag_EMC8[0]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTag_EMCJet[0]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistTrkPtTag_MB[1]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTag_EMC7[1]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTag_EMC8[1]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTag_EMCJet[1]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistTrkPtTag_MB[2]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTag_EMC7[2]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTag_EMC8[2]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTag_EMCJet[2]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistTrkPtTag_MB[3]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTag_EMC7[3]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTag_EMC8[3]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTag_EMCJet[3]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //Fill the DeDx by Pt plots
                
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistDeDxPtTag_MB[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTag_EMC7[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTag_EMC8[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTag_EMCJet[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistDeDxPtTag_MB[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTag_EMC7[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTag_EMC8[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTag_EMCJet[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistDeDxPtTag_MB[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTag_EMC7[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTag_EMC8[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTag_EMCJet[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistDeDxPtTag_MB[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTag_EMC7[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTag_EMC8[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTag_EMCJet[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
            }//end tag side check
            
            //Away side check
            if((dphi>3*TMath::Pi()/4&&dphi<5*TMath::Pi()/4)&&TMath::Abs(esdtrack->Eta())<EtaMax){
                //Add it to the multiplicity count
                away_Mult = away_Mult + 1;
                
                //Fill the Pt histos with the tracks Pt
                
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistTrkPtAway_MB[0]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtAway_EMC7[0]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtAway_EMC8[0]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtAway_EMCJet[0]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistTrkPtAway_MB[1]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtAway_EMC7[1]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtAway_EMC8[1]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtAway_EMCJet[1]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistTrkPtAway_MB[2]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtAway_EMC7[2]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtAway_EMC8[2]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtAway_EMCJet[2]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistTrkPtAway_MB[3]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtAway_EMC7[3]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtAway_EMC8[3]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtAway_EMCJet[3]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //Fill the DeDx by Pt plots
                
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistDeDxPtAway_MB[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtAway_EMC7[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtAway_EMC8[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtAway_EMCJet[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistDeDxPtAway_MB[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtAway_EMC7[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtAway_EMC8[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtAway_EMCJet[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistDeDxPtAway_MB[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtAway_EMC7[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtAway_EMC8[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtAway_EMCJet[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistDeDxPtAway_MB[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtAway_EMC7[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtAway_EMC8[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtAway_EMCJet[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
            }//end away side check
            
            //Left side check
            if((dphi>TMath::Pi()/4&&dphi<3*TMath::Pi()/4)&&TMath::Abs(esdtrack->Eta())<EtaMax){
                //Add it to the multiplicity count
                left_Mult = left_Mult + 1;
                
                //Add the Pt to the PtSum
                PtSumLeft = PtSumLeft + esdtrack->Pt();
            }
            
            //Right side check
            if((dphi>5*TMath::Pi()/4&&dphi<7*TMath::Pi()/4)&&TMath::Abs(esdtrack->Eta())<EtaMax){
                //Add it to the multiplicity count
                right_Mult = right_Mult + 1;
                
                //Add Pt to PtSum
                PtSumRight = PtSumRight + esdtrack->Pt();
                
                
            }
            
        }//End loop over all tracks
        
        //Fill multiplicity histos
        //tag side
        //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistTrkMultTag_MB[0]->Fill(tag_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTag_EMC7[0]->Fill(tag_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTag_EMC8[0]->Fill(tag_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTag_EMCJet[0]->Fill(tag_Mult);
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistTrkMultTag_MB[1]->Fill(tag_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTag_EMC7[1]->Fill(tag_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTag_EMC8[1]->Fill(tag_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTag_EMCJet[1]->Fill(tag_Mult);
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistTrkMultTag_MB[2]->Fill(tag_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTag_EMC7[2]->Fill(tag_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTag_EMC8[2]->Fill(tag_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTag_EMCJet[2]->Fill(tag_Mult);
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistTrkMultTag_MB[3]->Fill(tag_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTag_EMC7[3]->Fill(tag_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTag_EMC8[3]->Fill(tag_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTag_EMCJet[3]->Fill(tag_Mult);
                    }
                    
                }
        
        //away side
        //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistTrkMultAway_MB[0]->Fill(away_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultAway_EMC7[0]->Fill(away_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultAway_EMC8[0]->Fill(away_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultAway_EMCJet[0]->Fill(away_Mult);
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistTrkMultAway_MB[1]->Fill(away_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultAway_EMC7[1]->Fill(away_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultAway_EMC8[1]->Fill(away_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultAway_EMCJet[1]->Fill(away_Mult);
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistTrkMultAway_MB[2]->Fill(away_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultAway_EMC7[2]->Fill(away_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultAway_EMC8[2]->Fill(away_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultAway_EMCJet[2]->Fill(away_Mult);
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistTrkMultAway_MB[3]->Fill(away_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultAway_EMC7[3]->Fill(away_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultAway_EMC8[3]->Fill(away_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultAway_EMCJet[3]->Fill(away_Mult);
                    }
                    
                }
        
        
        
        
        
        //Fill transMin and transMax histos------------------------------------------------------------------------------------------
        if(PtSumLeft>PtSumRight){
            //Left is transMax
            
            //fill multiplicity histos
            //left side
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistTrkMultTransMax_MB[0]->Fill(left_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMax_EMC7[0]->Fill(left_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMax_EMC8[0]->Fill(left_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMax_EMCJet[0]->Fill(left_Mult);
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistTrkMultTransMax_MB[1]->Fill(left_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMax_EMC7[1]->Fill(left_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMax_EMC8[1]->Fill(left_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMax_EMCJet[1]->Fill(left_Mult);
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistTrkMultTransMax_MB[2]->Fill(left_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMax_EMC7[2]->Fill(left_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMax_EMC8[2]->Fill(left_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMax_EMCJet[2]->Fill(left_Mult);
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistTrkMultTransMax_MB[3]->Fill(left_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMax_EMC7[3]->Fill(left_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMax_EMC8[3]->Fill(left_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMax_EMCJet[3]->Fill(left_Mult);
                    }
                    
                }
            
                //right side
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistTrkMultTransMin_MB[0]->Fill(right_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMin_EMC7[0]->Fill(right_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMin_EMC8[0]->Fill(right_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMin_EMCJet[0]->Fill(right_Mult);
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistTrkMultTransMin_MB[1]->Fill(right_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMin_EMC7[1]->Fill(right_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMin_EMC8[1]->Fill(right_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMin_EMCJet[1]->Fill(right_Mult);
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistTrkMultTransMin_MB[2]->Fill(right_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMin_EMC7[2]->Fill(right_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMin_EMC8[2]->Fill(right_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMin_EMCJet[2]->Fill(right_Mult);
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistTrkMultTransMin_MB[3]->Fill(right_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMin_EMC7[3]->Fill(right_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMin_EMC8[3]->Fill(right_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMin_EMCJet[3]->Fill(right_Mult);
                    }
                    
                }
            
            for(Int_t j=0; j<ntracks; j++){
            
                if(i==j){continue;}//Don't double count our candidate
            
                AliESDtrack *esdtrack = esd->GetTrack(j);
                if(!esdtrack){continue;}
                
                Double_t dphi = esdtrack->Phi()-Phi_0;
                if(dphi<0){
                    dphi = (2*TMath::Pi()-Phi_0)+esdtrack->Phi();
                }
                
            //Left side check
            if((dphi>TMath::Pi()/4&&dphi<3*TMath::Pi()/4)&&TMath::Abs(esdtrack->Eta())<EtaMax){                
                //Fill the Pt histos with the tracks Pt
                
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistTrkPtTransMax_MB[0]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMax_EMC7[0]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMax_EMC8[0]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMax_EMCJet[0]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistTrkPtTransMax_MB[1]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMax_EMC7[1]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMax_EMC8[1]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMax_EMCJet[1]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistTrkPtTransMax_MB[2]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMax_EMC7[2]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMax_EMC8[2]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMax_EMCJet[2]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistTrkPtTransMax_MB[3]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMax_EMC7[3]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMax_EMC8[3]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMax_EMCJet[3]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //Fill the DeDx by Pt plots
                
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistDeDxPtTransMax_MB[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMax_EMC7[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMax_EMC8[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMax_EMCJet[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistDeDxPtTransMax_MB[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMax_EMC7[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMax_EMC8[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMax_EMCJet[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistDeDxPtTransMax_MB[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMax_EMC7[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMax_EMC8[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMax_EMCJet[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistDeDxPtTransMax_MB[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMax_EMC7[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMax_EMC8[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMax_EMCJet[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
            }//end left side check
                
                //transmin is da right
                if((dphi>5*TMath::Pi()/4&&dphi<7*TMath::Pi()/4)&&TMath::Abs(esdtrack->Eta())<EtaMax){
                //Fill the Pt histos with the tracks Pt
                
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistTrkPtTransMin_MB[0]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMin_EMC7[0]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMin_EMC8[0]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMin_EMCJet[0]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistTrkPtTransMin_MB[1]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMin_EMC7[1]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMin_EMC8[1]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMin_EMCJet[1]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistTrkPtTransMin_MB[2]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMin_EMC7[2]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMin_EMC8[2]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMin_EMCJet[2]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistTrkPtTransMin_MB[3]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMin_EMC7[3]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMin_EMC8[3]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMin_EMCJet[3]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //Fill the DeDx by Pt plots
                
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistDeDxPtTransMin_MB[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMin_EMC7[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMin_EMC8[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMin_EMCJet[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistDeDxPtTransMin_MB[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMin_EMC7[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMin_EMC8[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMin_EMCJet[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistDeDxPtTransMin_MB[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMin_EMC7[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMin_EMC8[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMin_EMCJet[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistDeDxPtTransMin_MB[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMin_EMC7[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMin_EMC8[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMin_EMCJet[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                }//end of transmin right side
                
            }//End for(regular tracks)
            }//End If left>right
            
            else{
                //Right is transMax
            
            //fill multiplicity histos
            //left side
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistTrkMultTransMin_MB[0]->Fill(left_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMin_EMC7[0]->Fill(left_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMin_EMC8[0]->Fill(left_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMin_EMCJet[0]->Fill(left_Mult);
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistTrkMultTransMin_MB[1]->Fill(left_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMin_EMC7[1]->Fill(left_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMin_EMC8[1]->Fill(left_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMin_EMCJet[1]->Fill(left_Mult);
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistTrkMultTransMin_MB[2]->Fill(left_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMin_EMC7[2]->Fill(left_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMin_EMC8[2]->Fill(left_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMin_EMCJet[2]->Fill(left_Mult);
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistTrkMultTransMin_MB[3]->Fill(left_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMin_EMC7[3]->Fill(left_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMin_EMC8[3]->Fill(left_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMin_EMCJet[3]->Fill(left_Mult);
                    }
                    
                }
            
                //right side
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistTrkMultTransMax_MB[0]->Fill(right_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMax_EMC7[0]->Fill(right_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMax_EMC8[0]->Fill(right_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMax_EMCJet[0]->Fill(right_Mult);
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistTrkMultTransMax_MB[1]->Fill(right_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMax_EMC7[1]->Fill(right_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMax_EMC8[1]->Fill(right_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMax_EMCJet[1]->Fill(right_Mult);
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistTrkMultTransMax_MB[2]->Fill(right_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMax_EMC7[2]->Fill(right_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMax_EMC8[2]->Fill(right_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMax_EMCJet[2]->Fill(right_Mult);
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistTrkMultTransMax_MB[3]->Fill(right_Mult);
                    }
                    if(EMC7trg){
                        fHistTrkMultTransMax_EMC7[3]->Fill(right_Mult);
                    }
                    
                    if(EMC8trg){
                        fHistTrkMultTransMax_EMC8[3]->Fill(right_Mult);
                    }
                    
                    if(EMCJettrg){
                        fHistTrkMultTransMax_EMCJet[3]->Fill(right_Mult);
                    }
                    
                }
            
            for(Int_t j=0; j<ntracks; j++){
            
                if(i==j){continue;}//Don't double count our candidate
            
                AliESDtrack *esdtrack = esd->GetTrack(j);
                if(!esdtrack){continue;}
                
                //Do hybrid track cuts
                if(!globaltrackCuts->AcceptTrack(esdtrack)&&!comptrackCuts->AcceptTrack(esdtrack)){continue;}
                
                Double_t dphi = esdtrack->Phi()-Phi_0;
                if(dphi<0){
                    dphi = (2*TMath::Pi()-Phi_0)+esdtrack->Phi();
                }
                
            //Left side check
            if((dphi>TMath::Pi()/4&&dphi<3*TMath::Pi()/4)&&TMath::Abs(esdtrack->Eta())<EtaMax){                
                //Fill the Pt histos with the tracks Pt
                
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistTrkPtTransMin_MB[0]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMin_EMC7[0]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMin_EMC8[0]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMin_EMCJet[0]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistTrkPtTransMin_MB[1]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMin_EMC7[1]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMin_EMC8[1]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMin_EMCJet[1]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistTrkPtTransMin_MB[2]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMin_EMC7[2]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMin_EMC8[2]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMin_EMCJet[2]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistTrkPtTransMin_MB[3]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMin_EMC7[3]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMin_EMC8[3]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMin_EMCJet[3]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //Fill the DeDx by Pt plots
                
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistDeDxPtTransMin_MB[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMin_EMC7[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMin_EMC8[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMin_EMCJet[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistDeDxPtTransMin_MB[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMin_EMC7[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMin_EMC8[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMin_EMCJet[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistDeDxPtTransMin_MB[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMin_EMC7[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMin_EMC8[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMin_EMCJet[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistDeDxPtTransMin_MB[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMin_EMC7[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMin_EMC8[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMin_EMCJet[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
            }
                
                //transmax is da right
                if((dphi>5*TMath::Pi()/4&&dphi<7*TMath::Pi()/4)&&TMath::Abs(esdtrack->Eta())<EtaMax){
                //Fill the Pt histos with the tracks Pt
                
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistTrkPtTransMax_MB[0]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMax_EMC7[0]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMax_EMC8[0]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMax_EMCJet[0]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistTrkPtTransMax_MB[1]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMax_EMC7[1]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMax_EMC8[1]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMax_EMCJet[1]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistTrkPtTransMax_MB[2]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMax_EMC7[2]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMax_EMC8[2]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMax_EMCJet[2]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistTrkPtTransMax_MB[3]->Fill(esdtrack->Pt());
                    }
                    if(EMC7trg){
                        fHistTrkPtTransMax_EMC7[3]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMC8trg){
                        fHistTrkPtTransMax_EMC8[3]->Fill(esdtrack->Pt());
                    }
                    
                    if(EMCJettrg){
                        fHistTrkPtTransMax_EMCJet[3]->Fill(esdtrack->Pt());
                    }
                    
                }
                
                //Fill the DeDx by Pt plots
                
                //tagged track 1<Pt<2
                if(elecTrk->Pt()>1&&elecTrk->Pt()<2){
                    if(MBtrg){
                        fHistDeDxPtTransMax_MB[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMax_EMC7[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMax_EMC8[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMax_EMCJet[0]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track 2<Pt<4
                if(elecTrk->Pt()>2&&elecTrk->Pt()<4){
                    if(MBtrg){
                        fHistDeDxPtTransMax_MB[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMax_EMC7[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMax_EMC8[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMax_EMCJet[1]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track 4<Pt<6
                if(elecTrk->Pt()>4&&elecTrk->Pt()<6){
                    if(MBtrg){
                        fHistDeDxPtTransMax_MB[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMax_EMC7[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMax_EMC8[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMax_EMCJet[2]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                
                //tagged track Pt>6
                if(elecTrk->Pt()>6){
                    if(MBtrg){
                        fHistDeDxPtTransMax_MB[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    if(EMC7trg){
                        fHistDeDxPtTransMax_EMC7[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMC8trg){
                        fHistDeDxPtTransMax_EMC8[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                    if(EMCJettrg){
                        fHistDeDxPtTransMax_EMCJet[3]->Fill(esdtrack->Pt(), esdtrack->GetTPCsignal()/esdtrack->GetTPCsignalSigma());
                    }
                    
                }
                }//end of transmin right side
                
            }//End for(regular tracks)
            }//End of else left<right
        
        
        
        //Check for another candidate in the away side, if yes call this a back-to-back and fill a different set of histograms
        
        Bool_t B2B=kFALSE;
        
        for(Int_t k=i+1;k<elecCnt;k++){
            AliESDtrack *elecTrk2 = esd->GetTrack(elecIDs[k]);
  
            if(!elecTrk2){
                AliWarning("The candidate electron track in null, something is wrong.");
                continue;
            }
            
            Double_t dphi = elecTrk2->Phi()-Phi_0;
            if(dphi<0){
                dphi = (2*TMath::Pi()-Phi_0)+elecTrk2->Phi();
            }
            
            if((dphi>3*TMath::Pi()/4&&dphi<5*TMath::Pi()/4)){
                B2B=kTRUE;
            }
        }
        //Compute which transverse region had the most Pt and fill the appropriate histogram
        if(B2B){
            if(PtSumLeft>PtSumRight){
                fHistPtSumTransMaxB2B->Fill(elecTrk->Pt(), PtSumLeft/(.4*TMath::Pi()));
                fHistPtSumTransMinB2B->Fill(elecTrk->Pt(), PtSumRight/(.4*TMath::Pi()));
            }else{
                fHistPtSumTransMaxB2B->Fill(elecTrk->Pt(), PtSumRight/(.4*TMath::Pi()));
                fHistPtSumTransMinB2B->Fill(elecTrk->Pt(), PtSumLeft/(.4*TMath::Pi()));
            }
        }else{
            if(PtSumLeft>PtSumRight){
                fHistPtSumTransMaxLead->Fill(elecTrk->Pt(), PtSumLeft/(.4*TMath::Pi()));
                fHistPtSumTransMinLead->Fill(elecTrk->Pt(), PtSumRight/(.4*TMath::Pi()));
            }else{
                fHistPtSumTransMaxLead->Fill(elecTrk->Pt(), PtSumRight/(.4*TMath::Pi()));
                fHistPtSumTransMinLead->Fill(elecTrk->Pt(), PtSumLeft/(.4*TMath::Pi()));
            }
        }
        
    }//End loop over candidate tracks
    
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMC8);
    PostData(4, fOutputEMCJet);
    return;
}

void AliAnalysisTaskPSHFE::FillPIDHistos(AliESDEvent *esd, AliESDtrack *esdtrack, AliPIDResponse *fPIDResponse){
        
    if(!esdtrack){
        AliWarning("esdtrack is null, no point in doing PID");
        return;
    }
    
    
        Bool_t isPIDRej = kFALSE;
    
        //Fill TOF and TPC status variables
        AliPIDResponse::EDetPidStatus TOFStatus=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, esdtrack);
        
        AliPIDResponse::EDetPidStatus TPCStatus=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, esdtrack);
        
        AliPIDResponse::EDetPidStatus TRDStatus=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTRD, esdtrack);
        
        AliPIDResponse::EDetPidStatus EMCStatus=fPIDResponse->CheckPIDStatus(AliPIDResponse::kEMCAL, esdtrack);
        
        //Check validity of PID, TODO: Add a rejection histogram
        if(TOFStatus!=AliPIDResponse::kDetPidOk){
            fHistPIDRejection->Fill(2);
            //isPIDRej=kTRUE;
        }
        
        if(TPCStatus!=AliPIDResponse::kDetPidOk){
            fHistPIDRejection->Fill(1);
            isPIDRej=kTRUE;
        }
        
        if(TRDStatus!=AliPIDResponse::kDetPidOk){
            fHistPIDRejection->Fill(3);
            isPIDRej=kTRUE;
        }
        
        if(EMCStatus!=AliPIDResponse::kDetPidOk){
            fHistPIDRejection->Fill(4);
            isPIDRej=kTRUE;
        }
        
        
        //Get the # of sigmas around an electron hypothesis for TOF and TPC
        Double_t nSigmaTOF;
        nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(esdtrack,AliPID::kElectron);
        
        Double_t nSigmaTPC;
        nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(esdtrack,AliPID::kElectron);
        
        Double_t elecLikeTRD[1];
        if(fPIDResponse->ComputeTRDProbability(esdtrack, AliPID::kElectron, elecLikeTRD, AliTRDPIDResponse::kLQ2D) != AliPIDResponse::kDetPidOk || esdtrack->GetTRDntrackletsPID()<4){
            fHistPIDRejection->Fill(3);
            isPIDRej=kTRUE;
        }
        
        if(isPIDRej){return;}
       
        //declare emcal cluster PID variables
        Double_t EOP=-1;
        Double_t M02=-1;
        Double_t M20=-1;
        Int_t caloId=esdtrack->GetEMCALcluster();
        
        if(caloId==-99999){
            fHistBadEMCclusID->Fill(1);
            return;
        }
        
        AliESDCaloCluster* tagEMCclus=esd->GetCaloCluster(caloId);
        
        if(tagEMCclus->E()>.5){
        EOP = tagEMCclus->E()/esdtrack->Pt();
        M02 = tagEMCclus->GetM02();
        M20 = tagEMCclus->GetM20();
        }
        else{
            return;
        }
      
        
        
        //=========================================================================================================================
        //PID cuts and histogram filling
        //=========================================================================================================================
        
        //Some Double arrays for pt ranges
        Double_t ptUpper[6] = {2, 3, 4, 5, 6, 1000000};
        Double_t ptLower[6] = {1, 2, 3, 4, 5, 6};
        Double_t TPCcut = 2;
        Double_t TOFcut = 2;
        Double_t TRDcut = .9;
        Double_t EMCcutLower[6] = {.85,.85,.85,.85,.85,.95};
        Double_t EMCcutHigher[6] = {1.15,1.15,1.15,1.15,1.15,1.25};
        
        for(Int_t i=0; i<6; i++){
        if(esdtrack->Pt()>ptLower[i]&&esdtrack->Pt()<ptUpper[i]){
            
            //Fill general shower shape plots
            if(MBtrg){
                    fHistM02_All_MB[i]->Fill(M02, EOP);
                    fHistM20_All_MB[i]->Fill(M20, EOP);
                }
                if(EMC7trg){
                    fHistM02_All_EMC7[i]->Fill(M02, EOP);
                    fHistM20_All_EMC7[i]->Fill(M20, EOP);
                }
                
                if(EMC8trg){
                    fHistM02_All_EMC8[i]->Fill(M02, EOP);
                    fHistM20_All_EMC8[i]->Fill(M20, EOP);
                }
                
                if(EMCJettrg){
                    fHistM02_All_EMCJet[i]->Fill(M02, EOP);
                    fHistM20_All_EMCJet[i]->Fill(M20, EOP);
                }
            
        //TPC Plots
            
            //TOF cuts
            if(nSigmaTOF<TOFcut&&nSigmaTOF>-TOFcut){
                if(MBtrg){
                    fHistTPC_TOF_MB[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                if(EMC7trg){
                    fHistTPC_TOF_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                
                if(EMC8trg){
                    fHistTPC_TOF_EMC8[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                
                if(EMCJettrg){
                    fHistTPC_TOF_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
            }
            
            //EMC cuts
            if(EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]){
                
                if(MBtrg){
                    fHistTPC_EMC_MB[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                if(EMC7trg){
                    fHistTPC_EMC_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                
                if(EMC8trg){
                    fHistTPC_EMC_EMC8[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                
                if(EMCJettrg){
                    fHistTPC_EMC_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
            }
            
            //TRD cuts
            if(elecLikeTRD[0]>TRDcut){
                
                if(MBtrg){
                    fHistTPC_TRD_MB[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                if(EMC7trg){
                    fHistTPC_TRD_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                
                if(EMC8trg){
                    fHistTPC_TRD_EMC8[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                
                if(EMCJettrg){
                    fHistTPC_TRD_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
            }
            
            //TOF+EMC cuts
            if(nSigmaTOF<TOFcut&&nSigmaTOF>-TOFcut&&EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]){
                
                if(MBtrg){
                    fHistTPC_TOFEMC_MB[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                if(EMC7trg){
                    fHistTPC_TOFEMC_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                
                if(EMC8trg){
                    fHistTPC_TOFEMC_EMC8[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                
                if(EMCJettrg){
                    fHistTPC_TOFEMC_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
            }
            
            //TOF+TRD cuts
            if(nSigmaTOF<TOFcut&&nSigmaTOF>-TOFcut&&elecLikeTRD[0]>TRDcut){
                
                if(MBtrg){
                    fHistTPC_TOFTRD_MB[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                if(EMC7trg){
                    fHistTPC_TOFTRD_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                
                if(EMC8trg){
                    fHistTPC_TOFTRD_EMC8[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                
                if(EMCJettrg){
                    fHistTPC_TOFTRD_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
            }
            
            //EMC+TRD cuts
            if(EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]&&elecLikeTRD[0]>TRDcut){
                
                if(MBtrg){
                    fHistTPC_EMCTRD_MB[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                if(EMC7trg){
                    fHistTPC_EMCTRD_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                
                if(EMC8trg){
                    fHistTPC_EMCTRD_EMC8[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                
                if(EMCJettrg){
                    fHistTPC_EMCTRD_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
            }
            
            //TOF+EMC+TRD cuts
            if(nSigmaTOF<TOFcut&&nSigmaTOF>-TOFcut&&EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]&&elecLikeTRD[0]>TRDcut){
                //tagStrong=kTRUE;
                if(MBtrg){
                    fHistTPC_TOFEMCTRD_MB[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                if(EMC7trg){
                    fHistTPC_TOFEMCTRD_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                
                if(EMC8trg){
                    fHistTPC_TOFEMCTRD_EMC8[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                
                if(EMCJettrg){
                    fHistTPC_TOFEMCTRD_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
            }
            
        //TOF Plots
            
            //TPC cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut){
                
                if(MBtrg){
                    fHistTOF_TPC_MB[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                if(EMC7trg){
                    fHistTOF_TPC_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                
                if(EMC8trg){
                    fHistTOF_TPC_EMC8[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                
                if(EMCJettrg){
                    fHistTOF_TPC_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
            }
            
            //EMC cuts
            if(EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]){
                
                if(MBtrg){
                    fHistTOF_EMC_MB[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                if(EMC7trg){
                    fHistTOF_EMC_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                
                if(EMC8trg){
                    fHistTOF_EMC_EMC8[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                
                if(EMCJettrg){
                    fHistTOF_EMC_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
            }
            
            //TRD cuts
            if(elecLikeTRD[0]>TRDcut){
                
                if(MBtrg){
                    fHistTOF_TRD_MB[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                if(EMC7trg){
                    fHistTOF_TRD_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                
                if(EMC8trg){
                    fHistTOF_TRD_EMC8[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                
                if(EMCJettrg){
                    fHistTOF_TRD_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
            }
            
            //TPC+EMC cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]){
                
                if(MBtrg){
                    fHistTOF_TPCEMC_MB[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                if(EMC7trg){
                    fHistTOF_TPCEMC_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                
                if(EMC8trg){
                    fHistTOF_TPCEMC_EMC8[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                
                if(EMCJettrg){
                    fHistTOF_TPCEMC_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
            }
            
            //TPC+TRD cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut){
                
                if(MBtrg){
                    fHistTOF_TPCTRD_MB[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                if(EMC7trg){
                    fHistTOF_TPCTRD_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                
                if(EMC8trg){
                    fHistTOF_TPCTRD_EMC8[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                
                if(EMCJettrg){
                    fHistTOF_TPCTRD_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
            }
            
            //EMC+TRD cuts
            if(EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]&&elecLikeTRD[0]>TRDcut){
                
                if(MBtrg){
                    fHistTOF_EMCTRD_MB[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                if(EMC7trg){
                    fHistTOF_EMCTRD_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                
                if(EMC8trg){
                    fHistTOF_EMCTRD_EMC8[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                
                if(EMCJettrg){
                    fHistTOF_EMCTRD_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
            }
            
            //TPC+EMC+TRD cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]&&elecLikeTRD[0]>TRDcut){
                //tagStrong=kTRUE;
                if(MBtrg){
                    fHistTOF_TPCEMCTRD_MB[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                if(EMC7trg){
                    fHistTOF_TPCEMCTRD_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                
                if(EMC8trg){
                    fHistTOF_TPCEMCTRD_EMC8[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
                
                if(EMCJettrg){
                    fHistTOF_TPCEMCTRD_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTOF);
                }
            }
            
        //EMC Plots
            
            //TPC cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut){
                
                if(MBtrg){
                    fHistEMC_TPC_MB[i]->Fill(EOP);
                }
                if(EMC7trg){
                    fHistEMC_TPC_EMC7[i]->Fill(EOP);
                }
                
                if(EMC8trg){
                    fHistEMC_TPC_EMC8[i]->Fill(EOP);
                }
                
                if(EMCJettrg){
                    fHistEMC_TPC_EMCJet[i]->Fill(EOP);
                }
            }
            
            //TOF cuts
            if(nSigmaTOF<TOFcut&&nSigmaTOF>-TOFcut){
                
                if(MBtrg){
                    fHistEMC_TOF_MB[i]->Fill(EOP);
                }
                if(EMC7trg){
                    fHistEMC_TOF_EMC7[i]->Fill(EOP);
                }
                
                if(EMC8trg){
                    fHistEMC_TOF_EMC8[i]->Fill(EOP);
                }
                
                if(EMCJettrg){
                    fHistEMC_TOF_EMCJet[i]->Fill(EOP);
                }
            }
            
            //TRD cuts
            if(elecLikeTRD[0]>TRDcut){
                
                if(MBtrg){
                    fHistEMC_TRD_MB[i]->Fill(EOP);
                }
                if(EMC7trg){
                    fHistEMC_TRD_EMC7[i]->Fill(EOP);
                }
                
                if(EMC8trg){
                    fHistEMC_TRD_EMC8[i]->Fill(EOP);
                }
                
                if(EMCJettrg){
                    fHistEMC_TRD_EMCJet[i]->Fill(EOP);
                }
            }
            
            //TPC+TOF cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&nSigmaTOF<TOFcut&&nSigmaTOF>-TOFcut){
                
                if(MBtrg){
                    fHistEMC_TPCTOF_MB[i]->Fill(EOP);
                }
                if(EMC7trg){
                    fHistEMC_TPCTOF_EMC7[i]->Fill(EOP);
                }
                
                if(EMC8trg){
                    fHistEMC_TPCTOF_EMC8[i]->Fill(EOP);
                }
                
                if(EMCJettrg){
                    fHistEMC_TPCTOF_EMCJet[i]->Fill(EOP);
                }
            }
            
            //TPC+TRD cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut){
                
                if(MBtrg){
                    fHistEMC_TPCTRD_MB[i]->Fill(EOP);
                }
                if(EMC7trg){
                    fHistEMC_TPCTRD_EMC7[i]->Fill(EOP);
                }
                
                if(EMC8trg){
                    fHistEMC_TPCTRD_EMC8[i]->Fill(EOP);
                }
                
                if(EMCJettrg){
                    fHistEMC_TPCTRD_EMCJet[i]->Fill(EOP);
                }
                
                //Fill the shower shape plots here also
                if(MBtrg){
                    fHistM02_Elec_MB[i]->Fill(M02, EOP);
                    fHistM20_Elec_MB[i]->Fill(M20, EOP);
                }
                if(EMC7trg){
                    fHistM02_Elec_EMC7[i]->Fill(M02, EOP);
                    fHistM20_Elec_EMC7[i]->Fill(M20, EOP);
                }
                
                if(EMC8trg){
                    fHistM02_Elec_EMC8[i]->Fill(M02, EOP);
                    fHistM20_Elec_EMC8[i]->Fill(M20, EOP);
                }
                
                if(EMCJettrg){
                    fHistM02_Elec_EMCJet[i]->Fill(M02, EOP);
                    fHistM20_Elec_EMCJet[i]->Fill(M20, EOP);
                }
            }
            
            //TOF+TRD cuts
            if(nSigmaTOF<TOFcut&&nSigmaTOF>-TOFcut&&elecLikeTRD[0]>TRDcut){
                
                if(MBtrg){
                    fHistEMC_TOFTRD_MB[i]->Fill(EOP);
                }
                if(EMC7trg){
                    fHistEMC_TOFTRD_EMC7[i]->Fill(EOP);
                }
                
                if(EMC8trg){
                    fHistEMC_TOFTRD_EMC8[i]->Fill(EOP);
                }
                
                if(EMCJettrg){
                    fHistEMC_TOFTRD_EMCJet[i]->Fill(EOP);
                }
            }
            
            //TPC+TOF+TRD cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&nSigmaTOF<TOFcut&&nSigmaTOF>-TOFcut&&elecLikeTRD[0]>TRDcut){
                //tagStrong=kTRUE;
                if(MBtrg){
                    fHistEMC_TPCTOFTRD_MB[i]->Fill(EOP);
                }
                if(EMC7trg){
                    fHistEMC_TPCTOFTRD_EMC7[i]->Fill(EOP);
                }
                
                if(EMC8trg){
                    fHistEMC_TPCTOFTRD_EMC8[i]->Fill(EOP);
                }
                
                if(EMCJettrg){
                    fHistEMC_TPCTOFTRD_EMCJet[i]->Fill(EOP);
                }
            }
            
        //TRD Plots
            
            //TPC cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut){
                
                if(MBtrg){
                    fHistTRD_TPC_MB[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                if(EMC7trg){
                    fHistTRD_TPC_EMC7[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                
                if(EMC8trg){
                    fHistTRD_TPC_EMC8[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                
                if(EMCJettrg){
                    fHistTRD_TPC_EMCJet[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
            }
            
            //TOF cuts
            if(nSigmaTOF<TOFcut&&nSigmaTOF>-TOFcut){
                
                if(MBtrg){
                    fHistTRD_TOF_MB[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                if(EMC7trg){
                    fHistTRD_TOF_EMC7[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                
                if(EMC8trg){
                    fHistTRD_TOF_EMC8[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                
                if(EMCJettrg){
                    fHistTRD_TOF_EMCJet[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
            }
            
            //EMC cuts
            if(EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]){
                
                if(MBtrg){
                    fHistTRD_EMC_MB[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                if(EMC7trg){
                    fHistTRD_EMC_EMC7[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                
                if(EMC8trg){
                    fHistTRD_EMC_EMC8[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                
                if(EMCJettrg){
                    fHistTRD_EMC_EMCJet[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
            }
            
            //TPC+TOF cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&nSigmaTOF<TOFcut&&nSigmaTOF>-TOFcut){
                
                if(MBtrg){
                    fHistTRD_TPCTOF_MB[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                if(EMC7trg){
                    fHistTRD_TPCTOF_EMC7[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                
                if(EMC8trg){
                    fHistTRD_TPCTOF_EMC8[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                
                if(EMCJettrg){
                    fHistTRD_TPCTOF_EMCJet[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
            }
            
            //TPC+EMC cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]){
                
                if(MBtrg){
                    fHistTRD_TPCEMC_MB[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                if(EMC7trg){
                    fHistTRD_TPCEMC_EMC7[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                
                if(EMC8trg){
                    fHistTRD_TPCEMC_EMC8[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                
                if(EMCJettrg){
                    fHistTRD_TPCEMC_EMCJet[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
            }
            
            //TOF+EMC cuts
            if(nSigmaTOF<TOFcut&&nSigmaTOF>-TOFcut&&EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]){
                
                if(MBtrg){
                    fHistTRD_TOFEMC_MB[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                if(EMC7trg){
                    fHistTRD_TOFEMC_EMC7[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                
                if(EMC8trg){
                    fHistTRD_TOFEMC_EMC8[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                
                if(EMCJettrg){
                    fHistTRD_TOFEMC_EMCJet[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
            }
            
            //TPC+TOF+EMC cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&nSigmaTOF<TOFcut&&nSigmaTOF>-TOFcut&&EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]){
                //tagStrong=kTRUE;
                if(MBtrg){
                    fHistTRD_TPCTOFEMC_MB[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                if(EMC7trg){
                    fHistTRD_TPCTOFEMC_EMC7[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                
                if(EMC8trg){
                    fHistTRD_TPCTOFEMC_EMC8[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                
                if(EMCJettrg){
                    fHistTRD_TPCTOFEMC_EMCJet[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
            }
        }
        }
    if(MBtrg){
    if(esdtrack->Pt()<2&&esdtrack->Pt()>1){
        if(nSigmaTPC<-2&&nSigmaTPC>-8){
            fHistEMC_Had_MB_1Gev->Fill(EOP);
        }
    }
    }
        //An electron candidate is one that passes TPC +-2Sig, TRD>.9, 0.85<E/p<1.15
        if(esdtrack->Pt()<6){
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut&&EOP<EMCcutHigher[0]&&EOP>EMCcutLower[0]){
                tagStrong=kTRUE;
            }
        }
        else{
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut&&EOP<EMCcutHigher[5]&&EOP>EMCcutLower[5]){
                tagStrong=kTRUE;
            }
        }
        
        /*//Check if any tracks pass all cuts
        if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&nSigmaTOF<TOFcut&&nSigmaTOF>-TOFcut&&elecLikeTRD[0]>TRDcut&&EOP<EMCcutHigher&&EOP>EMCcutLower){
            tagAll=kTRUE;
        }*/
        
        //=========================================================================
        //                     End of PID stuff
        //=========================================================================
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMC8);
    PostData(4, fOutputEMCJet);
    return;
}

void AliAnalysisTaskPSHFE::FillDPhiHistos(AliESDEvent *esd, AliESDtrack *esdtrack, Int_t i){
//Run through all tracks in this 
        Int_t ntracks = esd->GetNumberOfTracks();
            for(Int_t j = 0; j < ntracks; j++) {
                
                //Don't double count the tagged tracks
                if(i==j){continue;}
                
                
                AliESDtrack* esdtrackassoc = esd->GetTrack(j); // pointer to reconstructed to track          
                
                if(!esdtrackassoc) { 
                        AliError(Form("ERROR: Could not retrieve esdtrack %d",j)); 
                        continue; 
                }
                
                //Do hybrid track cuts
                if(!globaltrackCuts->AcceptTrack(esdtrackassoc)&&!comptrackCuts->AcceptTrack(esdtrackassoc)){continue;}
            
                //Pt distribution
                if(MBtrg){
                        fHistPtAssoc_MB->Fill(esdtrack->Pt());
                }
                if(EMC7trg){
                        fHistPtAssoc_EMC7->Fill(esdtrack->Pt());
                }
        
                if(EMC8trg){
                        fHistPtAssoc_EMC8->Fill(esdtrack->Pt());
                }
        
                if(EMCJettrg){
                        fHistPtAssoc_EMCJet->Fill(esdtrack->Pt());
                }
                
                //Fill Delta Phi variable and correct for periodicity
                Double_t DPhi=esdtrackassoc->Phi()-esdtrack->Phi();
                
                if(DPhi<-TMath::Pi()/2){DPhi=TMath::Abs(2*TMath::Pi()+DPhi);}
                
                if(DPhi>3*TMath::Pi()/2){DPhi=-TMath::Abs(2*TMath::Pi()-DPhi);}
                
                Double_t DEta=esdtrackassoc->Eta()-esdtrack->Eta();
                
                //candidate 1<pt<2
                if(esdtrack->Pt()>1&&esdtrack->Pt()<2){
                    
                //300-500MeV
                if(esdtrackassoc->Pt()>.3&&esdtrackassoc->Pt()<.5){
                    if(MBtrg){
                        fHistDPhi300_500_MB[0]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi300_500_EMC7[0]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi300_500_EMC8[0]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi300_500_EMCJet[0]->Fill(DPhi);
                    }
                }
                
                //500MeV
                if(esdtrackassoc->Pt()>.5&&esdtrackassoc->Pt()<.8){
                    if(MBtrg){
                        fHistDPhi500_800_MB[0]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi500_800_EMC7[0]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi500_800_EMC8[0]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi500_800_EMCJet[0]->Fill(DPhi);
                    }
                }
                
                //800MeV
                if(esdtrackassoc->Pt()>.8&&esdtrackassoc->Pt()<1){
                    if(MBtrg){
                        fHistDPhi800_1_MB[0]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi800_1_EMC7[0]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi800_1_EMC8[0]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi800_1_EMCJet[0]->Fill(DPhi);
                    }
                }
                
                //1GeV
                if(esdtrackassoc->Pt()>1&&esdtrackassoc->Pt()<2){
                    if(MBtrg){
                        fHistDPhi1_2_MB[0]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi1_2_EMC7[0]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi1_2_EMC8[0]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi1_2_EMCJet[0]->Fill(DPhi);
                    }
                }
                
                //2GeV
                if(esdtrackassoc->Pt()>2&&esdtrackassoc->Pt()<3){
                    if(MBtrg){
                        fHistDPhi2_3_MB[0]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi2_3_EMC7[0]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi2_3_EMC8[0]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi2_3_EMCJet[0]->Fill(DPhi);
                    }
                }
                
                //3GeV
                if(esdtrackassoc->Pt()>3&&esdtrackassoc->Pt()<4){
                    if(MBtrg){
                        fHistDPhi3_4_MB[0]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi3_4_EMC7[0]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi3_4_EMC8[0]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi3_4_EMCJet[0]->Fill(DPhi);
                    }
                }
                
                //5GeV
                if(esdtrackassoc->Pt()>4){
                    if(MBtrg){
                        fHistDPhi4_MB[0]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi4_EMC7[0]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi4_EMC8[0]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi4_EMCJet[0]->Fill(DPhi);
                    }
                }
                }
                
                //candidate 2<pt<4
                if(esdtrack->Pt()>2&&esdtrack->Pt()<4){
                    
                //300MeV
                if(esdtrackassoc->Pt()>.3&&esdtrackassoc->Pt()<.5){
                    if(MBtrg){
                        fHistDPhi300_500_MB[1]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi300_500_EMC7[1]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi300_500_EMC8[1]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi300_500_EMCJet[1]->Fill(DPhi);
                    }
                }
                
                //500MeV
                if(esdtrackassoc->Pt()>.5&&esdtrackassoc->Pt()<.8){
                    if(MBtrg){
                        fHistDPhi500_800_MB[1]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi500_800_EMC7[1]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi500_800_EMC8[1]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi500_800_EMCJet[1]->Fill(DPhi);
                    }
                }
                
                //800MeV
                if(esdtrackassoc->Pt()>.8&&esdtrackassoc->Pt()<1){
                    if(MBtrg){
                        fHistDPhi800_1_MB[1]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi800_1_EMC7[1]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi800_1_EMC8[1]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi800_1_EMCJet[1]->Fill(DPhi);
                    }
                }
                
                //1GeV
                if(esdtrackassoc->Pt()>1&&esdtrackassoc->Pt()<2){
                    if(MBtrg){
                        fHistDPhi1_2_MB[1]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi1_2_EMC7[1]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi1_2_EMC8[1]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi1_2_EMCJet[1]->Fill(DPhi);
                    }
                }
                
                //2GeV
                if(esdtrackassoc->Pt()>2&&esdtrackassoc->Pt()<3){
                    if(MBtrg){
                        fHistDPhi2_3_MB[1]->Fill(DPhi);
                        fHistDPhi28_MB->Fill(DPhi);
                        fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    }
                    if(EMC7trg){
                            fHistDPhi2_3_EMC7[1]->Fill(DPhi);
                            fHistDPhi28_EMC7->Fill(DPhi);
                            fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    }
                    if(EMC8trg){
                            fHistDPhi2_3_EMC8[1]->Fill(DPhi);
                            fHistDPhi28_EMC8->Fill(DPhi);
                            fHistDPhiDEta28_EMC8->Fill(DPhi, DEta);
                    }
                    if(EMCJettrg){
                            fHistDPhi2_3_EMCJet[1]->Fill(DPhi);
                            fHistDPhi28_EMCJet->Fill(DPhi);
                            fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    }
                }
                
                //3GeV
                if(esdtrackassoc->Pt()>3&&esdtrackassoc->Pt()<4){
                    if(MBtrg){
                        fHistDPhi3_4_MB[1]->Fill(DPhi);
                        fHistDPhi28_MB->Fill(DPhi);
                        fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    }
                    if(EMC7trg){
                            fHistDPhi3_4_EMC7[1]->Fill(DPhi);
                            fHistDPhi28_EMC7->Fill(DPhi);
                            fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    }
                    if(EMC8trg){
                            fHistDPhi3_4_EMC8[1]->Fill(DPhi);
                            fHistDPhi28_EMC8->Fill(DPhi);
                            fHistDPhiDEta28_EMC8->Fill(DPhi, DEta);
                    }
                    if(EMCJettrg){
                            fHistDPhi3_4_EMCJet[1]->Fill(DPhi);
                            fHistDPhi28_EMCJet->Fill(DPhi);
                            fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    }
                }
                
                //5GeV
                if(esdtrackassoc->Pt()>4){
                    if(MBtrg){
                        fHistDPhi4_MB[1]->Fill(DPhi);
                        fHistDPhi28_MB->Fill(DPhi);
                        fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    }
                    if(EMC7trg){
                            fHistDPhi4_EMC7[1]->Fill(DPhi);
                            fHistDPhi28_EMC7->Fill(DPhi);
                            fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    }
                    if(EMC8trg){
                            fHistDPhi4_EMC8[1]->Fill(DPhi);
                            fHistDPhi28_EMC8->Fill(DPhi);
                            fHistDPhiDEta28_EMC8->Fill(DPhi, DEta);
                    }
                    if(EMCJettrg){
                            fHistDPhi4_EMCJet[1]->Fill(DPhi);
                            fHistDPhi28_EMCJet->Fill(DPhi);
                            fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    }
                }                    
                }
            
                //candidate 4<pt<8
                if(esdtrack->Pt()>4&&esdtrack->Pt()<8){
                    
                //300MeV
                if(esdtrackassoc->Pt()>.3&&esdtrackassoc->Pt()<.5){
                    if(MBtrg){
                        fHistDPhi300_500_MB[2]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi300_500_EMC7[2]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi300_500_EMC8[2]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi300_500_EMCJet[2]->Fill(DPhi);
                    }
                }
                
                //500MeV
                if(esdtrackassoc->Pt()>.5&&esdtrackassoc->Pt()<.8){
                    if(MBtrg){
                           fHistDPhi500_800_MB[2]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi500_800_EMC7[2]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi500_800_EMC8[2]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi500_800_EMCJet[2]->Fill(DPhi);
                    }
                }
                
                //800MeV
                if(esdtrackassoc->Pt()>.8&&esdtrackassoc->Pt()<1){
                    if(MBtrg){
                        fHistDPhi800_1_MB[2]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi800_1_EMC7[2]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi800_1_EMC8[2]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi800_1_EMCJet[2]->Fill(DPhi);
                    }
                }
                
                //1GeV
                if(esdtrackassoc->Pt()>1&&esdtrackassoc->Pt()<2){
                    if(MBtrg){
                        fHistDPhi1_2_MB[2]->Fill(DPhi);
                    }
                    if(EMC7trg){
                            fHistDPhi1_2_EMC7[2]->Fill(DPhi);
                    }
                    if(EMC8trg){
                            fHistDPhi1_2_EMC8[2]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                            fHistDPhi1_2_EMCJet[2]->Fill(DPhi);
                    }
                }
                
                //2GeV
                if(esdtrackassoc->Pt()>2&&esdtrackassoc->Pt()<3){
                    if(MBtrg){
                        fHistDPhi2_3_MB[2]->Fill(DPhi);
                        fHistDPhi28_MB->Fill(DPhi);
                        fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    }
                    if(EMC7trg){
                            fHistDPhi2_3_EMC7[2]->Fill(DPhi);
                            fHistDPhi28_EMC7->Fill(DPhi);
                            fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    }
                    if(EMC8trg){
                            fHistDPhi2_3_EMC8[2]->Fill(DPhi);
                            fHistDPhi28_EMC8->Fill(DPhi);
                            fHistDPhiDEta28_EMC8->Fill(DPhi, DEta);
                    }
                    if(EMCJettrg){
                            fHistDPhi2_3_EMCJet[2]->Fill(DPhi);
                            fHistDPhi28_EMCJet->Fill(DPhi);
                            fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    }
                }
                
                //3GeV
                if(esdtrackassoc->Pt()>3&&esdtrackassoc->Pt()<4){
                    if(MBtrg){
                        fHistDPhi3_4_MB[2]->Fill(DPhi);
                        fHistDPhi28_MB->Fill(DPhi);
                        fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    }
                    if(EMC7trg){
                            fHistDPhi3_4_EMC7[2]->Fill(DPhi);
                            fHistDPhi28_EMC7->Fill(DPhi);
                            fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    }
                    if(EMC8trg){
                            fHistDPhi3_4_EMC8[2]->Fill(DPhi);
                            fHistDPhi28_EMC8->Fill(DPhi);
                            fHistDPhiDEta28_EMC8->Fill(DPhi, DEta);
                    }
                    if(EMCJettrg){
                            fHistDPhi3_4_EMCJet[2]->Fill(DPhi);
                            fHistDPhi28_EMCJet->Fill(DPhi);
                            fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    }
                }
                
                //5GeV
                if(esdtrackassoc->Pt()>4){
                    if(MBtrg){
                        fHistDPhi4_MB[2]->Fill(DPhi);
                        fHistDPhi28_MB->Fill(DPhi);
                        fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    }
                    if(EMC7trg){
                            fHistDPhi4_EMC7[2]->Fill(DPhi);
                            fHistDPhi28_EMC7->Fill(DPhi);
                            fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    }
                    if(EMC8trg){
                            fHistDPhi4_EMC8[2]->Fill(DPhi);
                            fHistDPhi28_EMC8->Fill(DPhi);
                            fHistDPhiDEta28_EMC8->Fill(DPhi, DEta);
                    }
                    if(EMCJettrg){
                            fHistDPhi4_EMCJet[2]->Fill(DPhi);
                            fHistDPhi28_EMCJet->Fill(DPhi);
                            fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    }
                }
                }
                
             }
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMC8);
    PostData(4, fOutputEMCJet);
    return;
}

void AliAnalysisTaskPSHFE::SetTrackCuts(AliESDtrackCuts *gtrkCuts, AliESDtrackCuts *ctrkCuts){
    cout<<"assigning hybrid track cuts \n";
    globaltrackCuts=gtrkCuts;
    comptrackCuts=ctrkCuts;
}

void AliAnalysisTaskPSHFE::SetElectronTrackCuts(Bool_t trkCutBool){
    trackCutsStrong=trkCutBool;
}

void AliAnalysisTaskPSHFE::FillPhotoElecHistos(AliESDEvent *esd, AliESDtrack *esdtrack, AliPIDResponse *fPIDResponse, Int_t i){
   
    if(!esdtrack){
        AliWarning("esdtrack is null, no point in doing Photonic Electron stuff");
        return;
    }
    
    Bool_t isElec=kFALSE;
    Double_t ElecMass=.0005109989;    
    //Fill TOF and TPC status variables  
    AliPIDResponse::EDetPidStatus TPCStatus=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, esdtrack);
        
    AliPIDResponse::EDetPidStatus TRDStatus=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTRD, esdtrack);
        
    AliPIDResponse::EDetPidStatus EMCStatus=fPIDResponse->CheckPIDStatus(AliPIDResponse::kEMCAL, esdtrack);
        
    //Check validity of PID, TODO: Add a rejection histogram
    if(TPCStatus!=AliPIDResponse::kDetPidOk){
        return;
    }
        
    if(TRDStatus!=AliPIDResponse::kDetPidOk){
        return;
    }
        
    if(EMCStatus!=AliPIDResponse::kDetPidOk){
        return;
    }
        
        
    //Get the # of sigmas around an electron hypothesis for TOF and TPC
    Double_t nSigmaTPC;
    nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(esdtrack,AliPID::kElectron);
        
    Double_t elecLikeTRD[1];
    if(fPIDResponse->ComputeTRDProbability(esdtrack, AliPID::kElectron, elecLikeTRD, AliTRDPIDResponse::kLQ2D) != AliPIDResponse::kDetPidOk || esdtrack->GetTRDntrackletsPID()<4){
        return;
    }
       
    //declare emcal cluster PID variables
    Double_t EOP=-1;
    Int_t caloId=esdtrack->GetEMCALcluster();
        
    if(caloId==-99999){
        return;
    }
        
    AliESDCaloCluster* tagEMCclus=esd->GetCaloCluster(caloId);
        
    EOP = tagEMCclus->E()/esdtrack->Pt();
    
    if((nSigmaTPC<2&&nSigmaTPC>-2)||(EOP<1.4&&EOP>.8)||(elecLikeTRD[0]>.8)){
        isElec=kTRUE;
    }
    
    if(!isElec){return;}
    
    Int_t ntracks=esd->GetNumberOfTracks();
    
    for(Int_t j = 0; j < ntracks; j++) {
                
        //Don't double count the tagged tracks
        if(i==j){continue;}
                
        AliESDtrack* esdtrackassoc = esd->GetTrack(j); // pointer to reconstructed to track          
                
        if(!esdtrackassoc) { 
            AliError(Form("ERROR: Could not retrieve esdtrack %d",j)); 
            continue; 
        }
        
        Bool_t isElecToo=kFALSE;
        
        //Fill TOF and TPC status variables  
        AliPIDResponse::EDetPidStatus TPCStatusassoc=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, esdtrackassoc);

        AliPIDResponse::EDetPidStatus TRDStatusassoc=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTRD, esdtrackassoc);

        //Check validity of PID, TODO: Add a rejection histogram
        if(TPCStatusassoc!=AliPIDResponse::kDetPidOk){
            continue;
        }

        if(TRDStatusassoc!=AliPIDResponse::kDetPidOk){
            continue;
        }


        //Get the # of sigmas around an electron hypothesis for TOF and TPC
        Double_t nSigmaTPCassoc;
        nSigmaTPCassoc = fPIDResponse->NumberOfSigmasTPC(esdtrackassoc,AliPID::kElectron);

        Double_t elecLikeTRDassoc[1];
        if(fPIDResponse->ComputeTRDProbability(esdtrackassoc, AliPID::kElectron, elecLikeTRDassoc, AliTRDPIDResponse::kLQ2D) != AliPIDResponse::kDetPidOk || esdtrackassoc->GetTRDntrackletsPID()<4){
            return;
        }

        if((nSigmaTPCassoc<2&&nSigmaTPCassoc>-2)||(elecLikeTRDassoc[0]>.8)){
            isElecToo=kTRUE;
        }
        
        if(!isElecToo){continue;}
        
        Double_t elecE1=TMath::Sqrt(esdtrack->P()*esdtrack->P()+ElecMass*ElecMass);
        Double_t elecE2=TMath::Sqrt(esdtrackassoc->P()*esdtrackassoc->P()+ElecMass*ElecMass);
        
        TLorentzVector elec1(esdtrack->Px(), esdtrack->Py(), esdtrack->Pz(), elecE1);
        TLorentzVector elec2(esdtrackassoc->Px(), esdtrackassoc->Py(), esdtrackassoc->Pz(), elecE2);
        
        Double_t InvMass=(elec1+elec2).M();
        Double_t OpAng=elec1.Angle(elec2.Vect());
        
        if(esdtrack->GetSign()==esdtrackassoc->GetSign()){
            
            if(MBtrg){
                fHistInvMassElecLike_MB->Fill(InvMass);
                fHistOpAngElecLike_MB->Fill(OpAng);
            }
            
            if(EMC7trg){
                fHistInvMassElecLike_EMC7->Fill(InvMass);
                fHistOpAngElecLike_EMC7->Fill(OpAng);
            }
            
            if(EMC8trg){
                fHistInvMassElecLike_EMC8->Fill(InvMass);
                fHistOpAngElecLike_EMC8->Fill(OpAng);
            }
            
            if(EMCJettrg){
                fHistInvMassElecLike_EMCJet->Fill(InvMass);
                fHistOpAngElecLike_EMCJet->Fill(OpAng);
            }
            
        }else{
            if(InvMass<0.1&&OpAng<1){tagPhot=kTRUE;}
            if(MBtrg){
                fHistInvMassElecUnLike_MB->Fill(InvMass);
                fHistOpAngElecUnLike_MB->Fill(OpAng);
            }
            
            if(EMC7trg){
                fHistInvMassElecUnLike_EMC7->Fill(InvMass);
                fHistOpAngElecUnLike_EMC7->Fill(OpAng);
            }
            
            if(EMC8trg){
                fHistInvMassElecUnLike_EMC8->Fill(InvMass);
                fHistOpAngElecUnLike_EMC8->Fill(OpAng);
            }
            
            if(EMCJettrg){
                fHistInvMassElecUnLike_EMCJet->Fill(InvMass);
                fHistOpAngElecUnLike_EMCJet->Fill(OpAng);
            }
        }
    }//end loop over tracks
    
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMC8);
    PostData(4, fOutputEMCJet);
    return;
}

//________________________________________________________________________
void AliAnalysisTaskPSHFE::Terminate(Option_t *) 
{
    
}
