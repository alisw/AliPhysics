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
#include "AliAODTrack.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliTRDPIDResponse.h"
#include "AliEventPoolManager.h"

using std::cout;
using std::endl;


ClassImp(AliAnalysisTaskPSHFE)

//________________________________________________________________________
AliAnalysisTaskPSHFE::AliAnalysisTaskPSHFE() // All data members should be initialised here
   :AliAnalysisTaskSE(),
    fOutputMB(0),
    fOutputEMC7(0),
    fOutputEMCJet(0),
    fTrackCutsStrong(0),
    fTrackCutsWeak(0),
    globaltrackCuts(0),
    comptrackCuts(0),
    fPoolMan(0),
    fPool(0),
    aodEv(0),
    trkArr(0),
    EMC7trg(0),
    EMCJettrg(0),
    MBtrg(0),
    tagStrong(0),
    tagPhot(0),

    fHistPIDRejection(0),
    fHistNElecPerEvent(0),

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
    fHistDPhiMix28_MB(0),
    fHistDPhiDEtaMix28_MB(0),
    fHistEMC_Had_MB_1Gev(0),
    fHistInvMassElecLike_MB(0),
    fHistInvMassElecUnLike_MB(0),
    fHistOpAngElecLike_MB(0),
    fHistOpAngElecUnLike_MB(0),
    fHistPtAssoc_MB(0),
    fHistPtAssocMix_MB(0),
    fHistPtTag_MB(0),
    fHistPhotoMismatch_MB(0),
    fHistDPhi28dEdx_MB(0),

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
    fHistDPhiMix28_EMC7(0),
    fHistDPhiDEtaMix28_EMC7(0),
    fHistInvMassElecLike_EMC7(0),
    fHistInvMassElecUnLike_EMC7(0),
    fHistOpAngElecLike_EMC7(0),
    fHistOpAngElecUnLike_EMC7(0),
    fHistPtAssoc_EMC7(0),
    fHistPtAssocMix_EMC7(0),
    fHistPtTag_EMC7(0),
    fHistPhotoMismatch_EMC7(0),
    fHistDPhi28dEdx_EMC7(0),

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
    fHistDPhiMix28_EMCJet(0),
    fHistDPhiDEtaMix28_EMCJet(0),
    fHistInvMassElecLike_EMCJet(0),
    fHistInvMassElecUnLike_EMCJet(0),
    fHistOpAngElecLike_EMCJet(0),
    fHistOpAngElecUnLike_EMCJet(0),
    fHistPtAssoc_EMCJet(0),
    fHistPtAssocMix_EMCJet(0),
    fHistPtTag_EMCJet(0),
    fHistPhotoMismatch_EMCJet(0),
    fHistDPhi28dEdx_EMCJet(0)
        
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

        fHistDPhi300_500_EMCJet[i]=0;
        fHistDPhi500_800_EMCJet[i]=0;
        fHistDPhi800_1_EMCJet[i]=0;
        fHistDPhi1_2_EMCJet[i]=0;
        fHistDPhi2_3_EMCJet[i]=0;
        fHistDPhi3_4_EMCJet[i]=0;
        fHistDPhi4_EMCJet[i]=0;
        
        //ME histos
        fHistDPhiMix300_500_MB[i]=0;
        fHistDPhiMix500_800_MB[i]=0;
        fHistDPhiMix800_1_MB[i]=0;
        fHistDPhiMix1_2_MB[i]=0;
        fHistDPhiMix2_3_MB[i]=0;
        fHistDPhiMix3_4_MB[i]=0;
        fHistDPhiMix4_MB[i]=0;
            
        fHistDPhiMix300_500_EMC7[i]=0;
        fHistDPhiMix500_800_EMC7[i]=0;
        fHistDPhiMix800_1_EMC7[i]=0;
        fHistDPhiMix1_2_EMC7[i]=0;
        fHistDPhiMix2_3_EMC7[i]=0;
        fHistDPhiMix3_4_EMC7[i]=0;
        fHistDPhiMix4_EMC7[i]=0;

        fHistDPhiMix300_500_EMCJet[i]=0;
        fHistDPhiMix500_800_EMCJet[i]=0;
        fHistDPhiMix800_1_EMCJet[i]=0;
        fHistDPhiMix1_2_EMCJet[i]=0;
        fHistDPhiMix2_3_EMCJet[i]=0;
        fHistDPhiMix3_4_EMCJet[i]=0;
        fHistDPhiMix4_EMCJet[i]=0;
    }
        
        //Init PID Plots here since they are stored in Arrays
    for(int i=0;i<6;i++){
        //MB Plots
        fHistTPC_EMCTRD_MB[i]=0;
        
        fHistEMC_TPCTRD_MB[i]=0;
        
        fHistTRD_TPCEMC_MB[i]=0;
        
        //EMC7 Plots
        fHistTPC_EMCTRD_EMC7[i]=0;

        fHistEMC_TPCTRD_EMC7[i]=0;
        
        fHistTRD_TPCEMC_EMC7[i]=0;
        
        //EMCJet Plots
        fHistTPC_EMCTRD_EMCJet[i]=0;

        fHistEMC_TPCTRD_EMCJet[i]=0;
        
        fHistTRD_TPCEMC_EMCJet[i]=0;

        fHistM02_All_MB[i]=0;
        fHistM02_Elec_MB[i]=0;
        fHistM20_All_MB[i]=0;
        fHistM20_Elec_MB[i]=0;
        
        fHistM02_All_EMC7[i]=0;
        fHistM02_Elec_EMC7[i]=0;
        fHistM20_All_EMC7[i]=0;
        fHistM20_Elec_EMC7[i]=0;
        
        fHistM02_All_EMCJet[i]=0;
        fHistM02_Elec_EMCJet[i]=0;
        fHistM20_All_EMCJet[i]=0;
        fHistM20_Elec_EMCJet[i]=0;
    }
    
    // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskPSHFE::AliAnalysisTaskPSHFE(const char *name) // All data members should be initialised here
   :AliAnalysisTaskSE(name),
    fOutputMB(0),
    fOutputEMC7(0),
    fOutputEMCJet(0),
    fTrackCutsStrong(0),
    fTrackCutsWeak(0),
    globaltrackCuts(0),
    comptrackCuts(0),
    fPoolMan(0),
    fPool(0),
    aodEv(0),
    trkArr(0),
    EMC7trg(0),
    EMCJettrg(0),
    MBtrg(0),
    tagStrong(0),
    tagPhot(0),

    fHistPIDRejection(0),
    fHistNElecPerEvent(0),

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
    fHistDPhiMix28_MB(0),
    fHistDPhiDEtaMix28_MB(0),
    fHistEMC_Had_MB_1Gev(0),
    fHistInvMassElecLike_MB(0),
    fHistInvMassElecUnLike_MB(0),
    fHistOpAngElecLike_MB(0),
    fHistOpAngElecUnLike_MB(0),
    fHistPtAssoc_MB(0),
    fHistPtAssocMix_MB(0),
    fHistPtTag_MB(0),
    fHistPhotoMismatch_MB(0),
    fHistDPhi28dEdx_MB(0),

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
    fHistDPhiMix28_EMC7(0),
    fHistDPhiDEtaMix28_EMC7(0),
    fHistInvMassElecLike_EMC7(0),
    fHistInvMassElecUnLike_EMC7(0),
    fHistOpAngElecLike_EMC7(0),
    fHistOpAngElecUnLike_EMC7(0),
    fHistPtAssoc_EMC7(0),
    fHistPtAssocMix_EMC7(0),
    fHistPtTag_EMC7(0),
    fHistPhotoMismatch_EMC7(0),
    fHistDPhi28dEdx_EMC7(0),

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
    fHistDPhiMix28_EMCJet(0),
    fHistDPhiDEtaMix28_EMCJet(0),
    fHistInvMassElecLike_EMCJet(0),
    fHistInvMassElecUnLike_EMCJet(0),
    fHistOpAngElecLike_EMCJet(0),
    fHistOpAngElecUnLike_EMCJet(0),
    fHistPtAssoc_EMCJet(0),
    fHistPtAssocMix_EMCJet(0),
    fHistPtTag_EMCJet(0),
    fHistPhotoMismatch_EMCJet(0),
    fHistDPhi28dEdx_EMCJet(0)
        
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
        
        fHistDPhi300_500_EMCJet[i]=0;
        fHistDPhi500_800_EMCJet[i]=0;
        fHistDPhi800_1_EMCJet[i]=0;
        fHistDPhi1_2_EMCJet[i]=0;
        fHistDPhi2_3_EMCJet[i]=0;
        fHistDPhi3_4_EMCJet[i]=0;
        fHistDPhi4_EMCJet[i]=0;
        
        //ME histos
        fHistDPhiMix300_500_MB[i]=0;
        fHistDPhiMix500_800_MB[i]=0;
        fHistDPhiMix800_1_MB[i]=0;
        fHistDPhiMix1_2_MB[i]=0;
        fHistDPhiMix2_3_MB[i]=0;
        fHistDPhiMix3_4_MB[i]=0;
        fHistDPhiMix4_MB[i]=0;
            
        fHistDPhiMix300_500_EMC7[i]=0;
        fHistDPhiMix500_800_EMC7[i]=0;
        fHistDPhiMix800_1_EMC7[i]=0;
        fHistDPhiMix1_2_EMC7[i]=0;
        fHistDPhiMix2_3_EMC7[i]=0;
        fHistDPhiMix3_4_EMC7[i]=0;
        fHistDPhiMix4_EMC7[i]=0;

        fHistDPhiMix300_500_EMCJet[i]=0;
        fHistDPhiMix500_800_EMCJet[i]=0;
        fHistDPhiMix800_1_EMCJet[i]=0;
        fHistDPhiMix1_2_EMCJet[i]=0;
        fHistDPhiMix2_3_EMCJet[i]=0;
        fHistDPhiMix3_4_EMCJet[i]=0;
        fHistDPhiMix4_EMCJet[i]=0;
    }
        
            //Init PID Plots here since they are stored in Arrays
    for(Int_t i=0;i<6;i++){
        //MB Plots
        fHistTPC_EMCTRD_MB[i]=0;

        fHistEMC_TPCTRD_MB[i]=0;

        fHistTRD_TPCEMC_MB[i]=0;
        
        //EMC7 Plots
        fHistTPC_EMCTRD_EMC7[i]=0;
        
        fHistEMC_TPCTRD_EMC7[i]=0;
        
        fHistTRD_TPCEMC_EMC7[i]=0;

        //EMCJet Plots
        fHistTPC_EMCTRD_EMCJet[i]=0;

        fHistEMC_TPCTRD_EMCJet[i]=0;
        
        fHistTRD_TPCEMC_EMCJet[i]=0;
        
        fHistM02_All_MB[i]=0;
        fHistM02_Elec_MB[i]=0;
        fHistM20_All_MB[i]=0;
        fHistM20_Elec_MB[i]=0;
        
        fHistM02_All_EMC7[i]=0;
        fHistM02_Elec_EMC7[i]=0;
        fHistM20_All_EMC7[i]=0;
        fHistM20_Elec_EMC7[i]=0;
        
        fHistM02_All_EMCJet[i]=0;
        fHistM02_Elec_EMCJet[i]=0;
        fHistM20_All_EMCJet[i]=0;
        fHistM20_Elec_EMCJet[i]=0;
    }
        
        
    DefineOutput(1, TList::Class());//MB
    DefineOutput(2, TList::Class());//EMC7
    DefineOutput(3, TList::Class());//EMCJet
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
        
        delete fHistDPhi300_500_EMCJet[i];
        delete fHistDPhi500_800_EMCJet[i];
        delete fHistDPhi800_1_EMCJet[i];
        delete fHistDPhi1_2_EMCJet[i];
        delete fHistDPhi2_3_EMCJet[i];
        delete fHistDPhi3_4_EMCJet[i];
        delete fHistDPhi4_EMCJet[i];
        
        //ME histos
        delete fHistDPhiMix300_500_MB[i];
        delete fHistDPhiMix500_800_MB[i];
        delete fHistDPhiMix800_1_MB[i];
        delete fHistDPhiMix1_2_MB[i];
        delete fHistDPhiMix2_3_MB[i];
        delete fHistDPhiMix3_4_MB[i];
        delete fHistDPhiMix4_MB[i];
        
        delete fHistDPhiMix300_500_EMC7[i];
        delete fHistDPhiMix500_800_EMC7[i];
        delete fHistDPhiMix800_1_EMC7[i];
        delete fHistDPhiMix1_2_EMC7[i];
        delete fHistDPhiMix2_3_EMC7[i];
        delete fHistDPhiMix3_4_EMC7[i];
        delete fHistDPhiMix4_EMC7[i];

        delete fHistDPhiMix300_500_EMCJet[i];
        delete fHistDPhiMix500_800_EMCJet[i];
        delete fHistDPhiMix800_1_EMCJet[i];
        delete fHistDPhiMix1_2_EMCJet[i];
        delete fHistDPhiMix2_3_EMCJet[i];
        delete fHistDPhiMix3_4_EMCJet[i];
        delete fHistDPhiMix4_EMCJet[i];
    }
    
        //Init PID Plots here since they are stored in Arrays
    for(Int_t i=0;i<6;i++){
        //MB Plots
        delete fHistTPC_EMCTRD_MB[i];

        delete fHistEMC_TPCTRD_MB[i];
        
        delete fHistTRD_TPCEMC_MB[i];
        
        //EMC7 Plots
        delete fHistTPC_EMCTRD_EMC7[i];

        delete fHistEMC_TPCTRD_EMC7[i];

        delete fHistTRD_TPCEMC_EMC7[i];
        
        //EMCJet Plots
        delete fHistTPC_EMCTRD_EMCJet[i];

        delete fHistEMC_TPCTRD_EMCJet[i];
        
        delete fHistTRD_TPCEMC_EMCJet[i];

        delete fHistM02_All_MB[i];
        delete fHistM02_Elec_MB[i];
        delete fHistM20_All_MB[i];
        delete fHistM20_Elec_MB[i];
        
        delete fHistM02_All_EMC7[i];
        delete fHistM02_Elec_EMC7[i];
        delete fHistM20_All_EMC7[i];
        delete fHistM20_Elec_EMC7[i];
        
        delete fHistM02_All_EMCJet[i];
        delete fHistM02_Elec_EMCJet[i];
        delete fHistM20_All_EMCJet[i];
        delete fHistM20_Elec_EMCJet[i];
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
    fOutputEMCJet = new TList();
    OpenFile(3);
    fOutputEMCJet->SetOwner();  // IMPORTANT!
    
    //Initialize event pool stuff
    Double_t vertexBins[5] = { -10, -4,  0, 4, 10 };
    Int_t nZvtxBins  = 4;
    Double_t multBins[4] = {0, 200, 500, 1000};
    Int_t nMultBins = 3;
    
    fPoolMan = new AliEventPoolManager(50, 50, nMultBins, multBins, nZvtxBins, vertexBins);
    fPoolMan->Validate();
    
    //Some strings for histograms
    TString ptRangesDPhi[3] = {"1-2Gev", "2-4Gev", "4-8Gev"};
    TString ptRangesPID[6] = {"1-2GeV", "2-3GeV", "3-4GeV", "4-5GeV", "5-6GeV", ">6GeV"};
    TString ptRangesRegion[4] = {"1-2Gev", "2-4Gev", "4-6Gev", ">6Gev"};
    
    //Strong cuts for heavy flavour
    fTrackCutsStrong = new AliESDtrackCuts();
    // TPC
    fTrackCutsStrong->SetRequireTPCRefit(kTRUE);
    // ITS
    fTrackCutsStrong->SetRequireITSRefit(kTRUE);
    // 7*(0.0026+0.0050/pt^1.01)
    fTrackCutsStrong->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    fTrackCutsStrong->SetMaxDCAToVertexZ(2);
    fTrackCutsStrong->SetDCAToVertex2D(kFALSE);
    fTrackCutsStrong->SetRequireSigmaToVertex(kFALSE);
    fTrackCutsStrong->SetMaxChi2PerClusterITS(36);
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
    fTrackCutsWeak = new AliESDtrackCuts();
    // TPC
    fTrackCutsWeak->SetRequireTPCRefit(kTRUE);
    // ITS
    fTrackCutsWeak->SetRequireITSRefit(kTRUE);
    // 7*(0.0026+0.0050/pt^1.01)
    fTrackCutsWeak->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    fTrackCutsWeak->SetMaxDCAToVertexZ(2);
    fTrackCutsWeak->SetDCAToVertex2D(kFALSE);
    fTrackCutsWeak->SetRequireSigmaToVertex(kFALSE);
    fTrackCutsWeak->SetMaxChi2PerClusterITS(36);
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
    
    fHistPhotoMismatch_MB = new TH1F("fHistPhotoMismatch_MB", "Electrons identified as 'heavy flavour' that fall in photonic invariant mass and opening angle cuts", 2, 0, 1);
    fHistPhotoMismatch_MB->GetXaxis()->SetTitle("Electrons");
    fHistPhotoMismatch_MB->GetYaxis()->SetTitle("Cts");

    fHistPhotoMismatch_EMC7 = new TH1F("fHistPhotoMismatch_EMC7", "Electrons identified as 'heavy flavour' that fall in photonic invariant mass and opening angle cuts", 2, 0, 1);
    fHistPhotoMismatch_EMC7->GetXaxis()->SetTitle("Electrons");
    fHistPhotoMismatch_EMC7->GetYaxis()->SetTitle("Cts");
    
    fHistPhotoMismatch_EMCJet = new TH1F("fHistPhotoMismatch_EMCJet", "Electrons identified as 'heavy flavour' that fall in photonic invariant mass and opening angle cuts", 2, 0, 1);
    fHistPhotoMismatch_EMCJet->GetXaxis()->SetTitle("Electrons");
    fHistPhotoMismatch_EMCJet->GetYaxis()->SetTitle("Cts");
    
    //Invariant mass histos
    
    fHistInvMassElecLike_MB = new TH1F("fHistInvMassElecLike_MB", "Invariant mass for all like-signed electron pairs", 100, 0, .5);
    fHistInvMassElecLike_MB->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecLike_MB->GetYaxis()->SetTitle("Cts");
    
    fHistInvMassElecLike_EMC7 = new TH1F("fHistInvMassElecLike_EMC7", "Invariant mass for all like-signed electron pairs", 100, 0, .5);
    fHistInvMassElecLike_EMC7->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecLike_EMC7->GetYaxis()->SetTitle("Cts");
    
    fHistInvMassElecLike_EMCJet = new TH1F("fHistInvMassElecLike_EMCJet", "Invariant mass for all like-signed electron pairs", 100, 0, .5);
    fHistInvMassElecLike_EMCJet->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecLike_EMCJet->GetYaxis()->SetTitle("Cts");
    
    fHistInvMassElecUnLike_MB = new TH1F("fHistInvMassElecUnLike_MB", "Invariant mass for all unlike-signed electron pairs", 100, 0, .5);
    fHistInvMassElecUnLike_MB->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecUnLike_MB->GetYaxis()->SetTitle("Cts");
    
    fHistInvMassElecUnLike_EMC7 = new TH1F("fHistInvMassElecUnLike_EMC7", "Invariant mass for all unlike-signed electron pairs", 100, 0, .5);
    fHistInvMassElecUnLike_EMC7->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecUnLike_EMC7->GetYaxis()->SetTitle("Cts");
    
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
    
    fHistOpAngElecLike_EMCJet = new TH1F("fHistOpAngElecLike_EMCJet", "Opening angle for all like-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecLike_EMCJet->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecLike_EMCJet->GetYaxis()->SetTitle("Cts");
    
    fHistOpAngElecUnLike_MB = new TH1F("fHistOpAngElecUnLike_MB", "Opening angle for all unlike-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecUnLike_MB->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecUnLike_MB->GetYaxis()->SetTitle("Cts");
    
    fHistOpAngElecUnLike_EMC7 = new TH1F("fHistOpAngElecUnLike_EMC7", "Opening angle for all unlike-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecUnLike_EMC7->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecUnLike_EMC7->GetYaxis()->SetTitle("Cts");
    
    fHistOpAngElecUnLike_EMCJet = new TH1F("fHistOpAngElecUnLike_EMCJet", "Opening angle for all unlike-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecUnLike_EMCJet->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecUnLike_EMCJet->GetYaxis()->SetTitle("Cts");
    
    //Rejection Histos
    
    fHistPIDRejection = new TH1F("fHistPIDRejection", "PID rejection counts for each detector.", 4, 1, 4);
    fHistPIDRejection->GetXaxis()->SetTitle("Detector");
    fHistPIDRejection->GetYaxis()->SetTitle("Cts");
    fHistPIDRejection->GetXaxis()->SetBinLabel(1, "TPC");
    fHistPIDRejection->GetXaxis()->SetBinLabel(2, "TOF");
    fHistPIDRejection->GetXaxis()->SetBinLabel(3, "TRD");
    fHistPIDRejection->GetXaxis()->SetBinLabel(4, "EMC");
    
    //Number of electrons per event histo
    fHistNElecPerEvent = new TH1F("fHistNElecPerEvent", "Number of tagged electrons per event", 5, 1, 5);
    fHistNElecPerEvent->GetXaxis()->SetTitle("Num. of Electrons");
    fHistNElecPerEvent->GetYaxis()->SetTitle("Cts");
    
    //PID Plots
    
    //TPC PID Plots
    for(Int_t i=0; i<6; i++){
        //MB        
        fHistTPC_EMCTRD_MB[i] = new TH2F(TString::Format("fHistTPC_EMCTRD_MB_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_EMCTRD_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_EMCTRD_MB[i]->GetYaxis()->SetTitle("nSigma");
        
        //EMC7
        fHistTPC_EMCTRD_EMC7[i] = new TH2F(TString::Format("fHistTPC_EMCTRD_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_EMCTRD_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_EMCTRD_EMC7[i]->GetYaxis()->SetTitle("nSigma");
        
        //EMCJet
        fHistTPC_EMCTRD_EMCJet[i] = new TH2F(TString::Format("fHistTPC_EMCTRD_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_EMCTRD_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_EMCTRD_EMCJet[i]->GetYaxis()->SetTitle("nSigma");
    }
    
    //EMC PID Plots
    for(Int_t i=0; i<6; i++){
        //MB
        fHistEMC_TPCTRD_MB[i] = new TH1F(TString::Format("fHistEMC_TPCTRD_MB_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TPC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTRD_MB[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTRD_MB[i]->GetYaxis()->SetTitle("Cts");
        
        //EMC7
        fHistEMC_TPCTRD_EMC7[i] = new TH1F(TString::Format("fHistEMC_TPCTRD_EMC7_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TPC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTRD_EMC7[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTRD_EMC7[i]->GetYaxis()->SetTitle("Cts");
        
        //EMCJet
        fHistEMC_TPCTRD_EMCJet[i] = new TH1F(TString::Format("fHistEMC_TPCTRD_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTRD_EMCJet[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTRD_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    //TRD PID Plots
    for(Int_t i=0; i<6; i++){
        //MB
        fHistTRD_TPCEMC_MB[i] = new TH2F(TString::Format("fHistTRD_TPCEMC_MB_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TPC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCEMC_MB[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCEMC_MB[i]->GetYaxis()->SetTitle("electron Likelihood");
    
        //EMC7
        fHistTRD_TPCEMC_EMC7[i] = new TH2F(TString::Format("fHistTRD_TPCEMC_EMC7_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TPC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCEMC_EMC7[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCEMC_EMC7[i]->GetYaxis()->SetTitle("electron Likelihood");
        
        //EMCJet
        fHistTRD_TPCEMC_EMCJet[i] = new TH2F(TString::Format("fHistTRD_TPCEMC_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after EMC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCEMC_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCEMC_EMCJet[i]->GetYaxis()->SetTitle("electron Likelihood");
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
    
    fHistDPhiDEta28_EMCJet = new TH2F("fHistDPhiDEta28_EMCJet", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2, 100, -0.9, 0.9);
    fHistDPhiDEta28_EMCJet->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEta28_EMCJet->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEta28_EMCJet->GetZaxis()->SetTitle("Cts");
    
    //DPhi for candidate electrons 2-8 gev and assoc. particles >2Gev for Mixed Events
    fHistDPhiMix28_MB = new TH1F("fHistDPhiMix28_MB", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev for Mixed Events", 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhiMix28_MB->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiMix28_MB->GetYaxis()->SetTitle("Cts");
    
    fHistDPhiMix28_EMC7 = new TH1F("fHistDPhiMix28_EMC7", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev for Mixed Events", 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhiMix28_EMC7->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiMix28_EMC7->GetYaxis()->SetTitle("Cts");
    
    fHistDPhiMix28_EMCJet = new TH1F("fHistDPhiMix28_EMCJet", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev for Mixed Events", 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhiMix28_EMCJet->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiMix28_EMCJet->GetYaxis()->SetTitle("Cts");
    
    //DPhi by Eta for triggered particles 2-8 gev and assoc. particles >3gev
    fHistDPhiDEtaMix28_MB = new TH2F("fHistDPhiDEtaMix28_MB", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev for Mixed Events", 100, -TMath::Pi()/2, 3*TMath::Pi()/2, 100, -0.9, 0.9);
    fHistDPhiDEtaMix28_MB->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEtaMix28_MB->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEtaMix28_MB->GetZaxis()->SetTitle("Cts");
    
    fHistDPhiDEtaMix28_EMC7 = new TH2F("fHistDPhiDEtaMix28_EMC7", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev for Mixed Events", 100, -TMath::Pi()/2, 3*TMath::Pi()/2, 100, -0.9, 0.9);
    fHistDPhiDEtaMix28_EMC7->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEtaMix28_EMC7->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEtaMix28_EMC7->GetZaxis()->SetTitle("Cts");
    
    fHistDPhiDEtaMix28_EMCJet = new TH2F("fHistDPhiDEtaMix28_EMCJet", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev for Mixed Events", 100, -TMath::Pi()/2, 3*TMath::Pi()/2, 100, -0.9, 0.9);
    fHistDPhiDEtaMix28_EMCJet->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEtaMix28_EMCJet->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEtaMix28_EMCJet->GetZaxis()->SetTitle("Cts");
    
    //DPhi by dEdx for triggered particles 2-8 gev and assoc. particles >2gev
    fHistDPhi28dEdx_MB = new TH2F("fHistDPhi28dEdx_MB", "Delta-Phi by dE/dx for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2, 300, -30, 130);
    fHistDPhi28dEdx_MB->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi28dEdx_MB->GetYaxis()->SetTitle("dE/dx");
    fHistDPhi28dEdx_MB->GetZaxis()->SetTitle("Cts");
    
    fHistDPhi28dEdx_EMC7 = new TH2F("fHistDPhi28dEdx_EMC7", "Delta-Phi by dE/dx for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2, 300, -30, 130);
    fHistDPhi28dEdx_EMC7->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi28dEdx_EMC7->GetYaxis()->SetTitle("dE/dx");
    fHistDPhi28dEdx_EMC7->GetZaxis()->SetTitle("Cts");
    
    fHistDPhi28dEdx_EMCJet = new TH2F("fHistDPhiDEta28_EMCJet", "Delta-Phi by dE/dx for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2, 300, -30, 130);
    fHistDPhi28dEdx_EMCJet->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi28dEdx_EMCJet->GetYaxis()->SetTitle("dE/dx");
    fHistDPhi28dEdx_EMCJet->GetZaxis()->SetTitle("Cts");
        
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
        fHistDPhi4_EMCJet[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_4Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated pt>4Gev",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi4_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi4_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    //Mixed Event DPhi plots
    
    // Delta Phi for tracks > 300MeV
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix300_500_MB[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_.3-.5Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<.5Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix300_500_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix300_500_MB[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix300_500_EMC7[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_.3-.5Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<.5Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix300_500_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix300_500_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix300_500_EMCJet[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_.3-.5Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<.5Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix300_500_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix300_500_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    // Delta Phi for tracks > 500MeV
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix500_800_MB[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_.5-.8Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .5Gev<pt<.8Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix500_800_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix500_800_MB[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix500_800_EMC7[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_.5-.8Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .5Gev<pt<.8Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix500_800_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix500_800_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix500_800_EMCJet[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_.5-.8Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .5Gev<pt<.8Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix500_800_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix500_800_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    // Delta Phi for tracks > 800MeV
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix800_1_MB[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_.8-1Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .8Gev<pt<1Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix800_1_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix800_1_MB[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix800_1_EMC7[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_.8-1Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .8Gev<pt<1Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix800_1_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix800_1_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix800_1_EMCJet[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_.8-1Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .8Gev<pt<1Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix800_1_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix800_1_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    // Delta Phi for tracks > 1GeV
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix1_2_MB[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_1-2Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix1_2_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix1_2_MB[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix1_2_EMC7[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_1-2Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix1_2_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix1_2_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix1_2_EMCJet[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_1-2Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix1_2_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix1_2_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    // Delta Phi for tracks > 2GeV
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix2_3_MB[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_2-3Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<3Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix2_3_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix2_3_MB[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix2_3_EMC7[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_2-3Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<3Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix2_3_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix2_3_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix2_3_EMCJet[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_2-3Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<3Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix2_3_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix2_3_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    // Delta Phi for tracks > 3GeV
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix3_4_MB[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_3-4Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 3Gev<pt<4Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix3_4_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix3_4_MB[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix3_4_EMC7[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_3-4Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 3Gev<pt<4Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix3_4_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix3_4_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix3_4_EMCJet[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_3-4Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 3Gev<pt<4Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix3_4_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix3_4_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }
    
    // Delta Phi for tracks > 5GeV
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix4_MB[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_4Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated pt>4Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix4_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix4_MB[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix4_EMC7[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_4Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated pt>4Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix4_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix4_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }
    
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix4_EMCJet[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_4Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated pt>4Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix4_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix4_EMCJet[i]->GetYaxis()->SetTitle("Cts");
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
    
    fHistTPCSigCut_EMCJet = new TH1F("fHistTPCSigCut_EMCJet", "dEdx Resolution in TPC for tracks that pass the basic track cuts", 100,0,100);
    fHistTPCSigCut_EMCJet->GetXaxis()->SetTitle("Resolution");
    fHistTPCSigCut_EMCJet->GetYaxis()->SetTitle("Count");
    
    //Impact Parameter histos
    fHistImpPar_MB = new TH1F("fHistImpPar_MB", "Impact Parameter distribution in xy plane for all tracks", 100,-1, 1);
    fHistImpPar_MB->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpPar_MB->GetYaxis()->SetTitle("Count");
    
    fHistImpPar_EMC7 = new TH1F("fHistImpPar_EMC7", "Impact Parameter distribution in xy plane for all tracks", 100,-1, 1);
    fHistImpPar_EMC7->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpPar_EMC7->GetYaxis()->SetTitle("Count");
    
    fHistImpPar_EMCJet = new TH1F("fHistImpPar_EMCJet", "Impact Parameter distribution in xy plane for all tracks", 100,-1, 1);
    fHistImpPar_EMCJet->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpPar_EMCJet->GetYaxis()->SetTitle("Count");
    
    //Impact Parameter for tagged electrons histos
    fHistImpParTag_MB = new TH1F("fHistImpParTag_MB", "Impact Parameter distribution in xy plane for electron candidates", 100,-1, 1);
    fHistImpParTag_MB->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpParTag_MB->GetYaxis()->SetTitle("Count");
    
    fHistImpParTag_EMC7 = new TH1F("fHistImpParTag_EMC7", "Impact Parameter distribution in xy plane for electron candidates", 100,-1, 1);
    fHistImpParTag_EMC7->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpParTag_EMC7->GetYaxis()->SetTitle("Count");
    
    fHistImpParTag_EMCJet = new TH1F("fHistImpParTag_EMCJet", "Impact Parameter distribution in xy plane for electron candidates", 100,-1, 1);
    fHistImpParTag_EMCJet->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpParTag_EMCJet->GetYaxis()->SetTitle("Count");
    
    //Number of clusters per track in TPC
    fHistTPCNClus_MB = new TH1F("fHistTPCNClus_MB", "Number of Clusters per track in TPC", 159, 0, 159);
    fHistTPCNClus_MB->GetXaxis()->SetTitle("Number of TPC Clusters");
    fHistTPCNClus_MB->GetYaxis()->SetTitle("Number of Tracks");
    
    fHistTPCNClus_EMC7 = new TH1F("fHistTPCNClus_EMC7", "Number of Clusters per track in TPC", 159, 0, 159);
    fHistTPCNClus_EMC7->GetXaxis()->SetTitle("Number of TPC Clusters");
    fHistTPCNClus_EMC7->GetYaxis()->SetTitle("Number of Tracks");
    
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
    
    fHistPtAssoc_EMCJet = new TH1F("fHistPtAssoc_EMCJet", "Pt distribution for associated tracks", 100,0, 15);
    fHistPtAssoc_EMCJet->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtAssoc_EMCJet->GetYaxis()->SetTitle("Count");
    
    // Pt distribution of all tracks for mixed events
    fHistPtAssocMix_MB = new TH1F("fHistPtAssocMix_MB", "Pt distribution for associated tracks in mixed events", 100,0, 15);
    fHistPtAssocMix_MB->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtAssocMix_MB->GetYaxis()->SetTitle("Count");
    
    fHistPtAssocMix_EMC7 = new TH1F("fHistPtAssocMix_EMC7", "Pt distribution for associated tracks in mixed events", 100,0, 15);
    fHistPtAssocMix_EMC7->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtAssocMix_EMC7->GetYaxis()->SetTitle("Count");
    
    fHistPtAssocMix_EMCJet = new TH1F("fHistPtAssocMix_EMCJet", "Pt distribution for associated tracks in mixed events", 100,0, 15);
    fHistPtAssocMix_EMCJet->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtAssocMix_EMCJet->GetYaxis()->SetTitle("Count");
    
    //Impact Parameter for tagged electrons histos
    fHistPtTag_MB = new TH1F("fHistPtTag_MB", "Pt distribution for electron candidates", 100,0, 15);
    fHistPtTag_MB->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtTag_MB->GetYaxis()->SetTitle("Count");
    
    fHistPtTag_EMC7 = new TH1F("fHistPtTag_EMC7", "Pt distribution for electron candidates", 100,0, 15);
    fHistPtTag_EMC7->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtTag_EMC7->GetYaxis()->SetTitle("Count");
    
    fHistPtTag_EMCJet = new TH1F("fHistPtTag_EMCJet", "Pt distribution for electron candidates", 100,0, 15);
    fHistPtTag_EMCJet->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtTag_EMCJet->GetYaxis()->SetTitle("Count");
        
    //Add rejection plots to MB plots since it is the easiest place
    fOutputMB->Add(fHistPIDRejection);
    fOutputMB->Add(fHistNElecPerEvent);
    
    fOutputMB->Add(fHistPhotoMismatch_MB);
    fOutputMB->Add(fHistPtAssoc_MB);
    fOutputMB->Add(fHistPtAssocMix_MB);
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
        fOutputMB->Add(fHistTPC_EMCTRD_MB[i]);
        
        fOutputMB->Add(fHistEMC_TPCTRD_MB[i]);
        
        fOutputMB->Add(fHistTRD_TPCEMC_MB[i]);
        
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
    fOutputMB->Add(fHistDPhiMix300_500_MB[i]);
    fOutputMB->Add(fHistDPhiMix500_800_MB[i]);
    fOutputMB->Add(fHistDPhiMix800_1_MB[i]);
    fOutputMB->Add(fHistDPhiMix1_2_MB[i]);
    fOutputMB->Add(fHistDPhiMix2_3_MB[i]);
    fOutputMB->Add(fHistDPhiMix3_4_MB[i]);
    fOutputMB->Add(fHistDPhiMix4_MB[i]);
    }
    fOutputMB->Add(fHistDPhi28_MB);
    fOutputMB->Add(fHistDPhiDEta28_MB);
    fOutputMB->Add(fHistDPhiMix28_MB);
    fOutputMB->Add(fHistDPhiDEtaMix28_MB);
    fOutputMB->Add(fHistDPhi28dEdx_MB);
    
    fOutputEMC7->Add(fHistPhotoMismatch_EMC7);
    fOutputEMC7->Add(fHistPtAssoc_EMC7);
    fOutputEMC7->Add(fHistPtAssocMix_EMC7);
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
        fOutputEMC7->Add(fHistTPC_EMCTRD_EMC7[i]);
        
        fOutputEMC7->Add(fHistEMC_TPCTRD_EMC7[i]);
        
        fOutputEMC7->Add(fHistTRD_TPCEMC_EMC7[i]);
        
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
    fOutputEMC7->Add(fHistDPhiMix300_500_EMC7[i]);
    fOutputEMC7->Add(fHistDPhiMix500_800_EMC7[i]);
    fOutputEMC7->Add(fHistDPhiMix800_1_EMC7[i]);
    fOutputEMC7->Add(fHistDPhiMix1_2_EMC7[i]);
    fOutputEMC7->Add(fHistDPhiMix2_3_EMC7[i]);
    fOutputEMC7->Add(fHistDPhiMix3_4_EMC7[i]);
    fOutputEMC7->Add(fHistDPhiMix4_EMC7[i]);
    }
    fOutputEMC7->Add(fHistDPhi28_EMC7);
    fOutputEMC7->Add(fHistDPhiDEta28_EMC7);
    fOutputEMC7->Add(fHistDPhiMix28_EMC7);
    fOutputEMC7->Add(fHistDPhiDEtaMix28_EMC7);
    fOutputEMC7->Add(fHistDPhi28dEdx_EMC7);
    
    fOutputEMCJet->Add(fHistPhotoMismatch_EMCJet);
    fOutputEMCJet->Add(fHistPtAssoc_EMCJet);
    fOutputEMCJet->Add(fHistPtAssocMix_EMCJet);
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
        fOutputEMCJet->Add(fHistTPC_EMCTRD_EMCJet[i]);
 
        fOutputEMCJet->Add(fHistEMC_TPCTRD_EMCJet[i]);
       
        fOutputEMCJet->Add(fHistTRD_TPCEMC_EMCJet[i]);
       
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
    fOutputEMCJet->Add(fHistDPhiMix300_500_EMCJet[i]);
    fOutputEMCJet->Add(fHistDPhiMix500_800_EMCJet[i]);
    fOutputEMCJet->Add(fHistDPhiMix800_1_EMCJet[i]);
    fOutputEMCJet->Add(fHistDPhiMix1_2_EMCJet[i]);
    fOutputEMCJet->Add(fHistDPhiMix2_3_EMCJet[i]);
    fOutputEMCJet->Add(fHistDPhiMix3_4_EMCJet[i]);
    fOutputEMCJet->Add(fHistDPhiMix4_EMCJet[i]);
    }
    fOutputEMCJet->Add(fHistDPhi28_EMCJet);
    fOutputEMCJet->Add(fHistDPhiDEta28_EMCJet);
    fOutputEMCJet->Add(fHistDPhiMix28_EMCJet);
    fOutputEMCJet->Add(fHistDPhiDEtaMix28_EMCJet);
    fOutputEMCJet->Add(fHistDPhi28dEdx_EMCJet);

    // NEW HISTO added to fOutput here
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMCJet);// Post data for ALL output slots >0 here, to get at least an empty histogram
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
    
    if(aodEv){cout<<"This is an AOD event\n";}
    
    AliESDEvent* esd;
    // create pointer to event
    AliAODEvent* aod = dynamic_cast<AliAODEvent*>(event);
    if (!aod) {
        AliError("Cannot get the AOD event");
        return;
    }  
    if(!aodEv){
        esd = dynamic_cast<AliESDEvent*>(event);
        cout<<"Made ESD Event\n";
        if (!esd) {
            AliError("Cannot get the ESD event");
            return;
        }  
    }
    
    // input handler
    const AliAnalysisManager* man(AliAnalysisManager::GetAnalysisManager());
    if (NULL == man) {
        AliWarning("AliAnalysisManager is not available");
        return;
    }

    AliInputEventHandler* inputHandler = (AliInputEventHandler*)man->GetInputEventHandler();
    if (!inputHandler) {
            AliWarning("AliInputEventHandler is not available");
            return;
    }
    /*if(!aodEv){
        AliESDInputHandler* inputHandlerESD = dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler());  
        if (!inputHandlerESD) {
            AliWarning("AliESDInputHandler is not available");
            return;
        }
    }else{
        AliAODInputHandler* inputHandlerAOD = (AliAODInputHandler*)(man->GetInputEventHandler());  
        if (!inputHandlerAOD) {
            AliWarning("AliAODInputHandler is not available");
            return;
        }
    }*/
    
    UInt_t fSelectMask = inputHandler->IsEventSelected();
    /*if(!aodEv){
        UInt_t fSelectMask = inputHandler->IsEventSelected();
    }else{
        UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    }       */ 
    
    Bool_t isSelected = fSelectMask & (AliVEvent::kEMC7 | AliVEvent::kEMCEJE | AliVEvent::kEMC8);
    if(!isSelected){
        AliWarning("This is not an EMCal triggered event");
    }

    MBtrg = fSelectMask & AliVEvent::kAnyINT;
    EMC7trg = fSelectMask & AliVEvent::kEMC7;
    EMCJettrg = fSelectMask & AliVEvent::kEMCEJE;

    Int_t elecIDs[1000];
    Int_t elecCnt=0;
    
    if(!aodEv){
        if(!globaltrackCuts||!comptrackCuts){
            AliWarning("The hybrid track cuts are null");
            return;
        }
    }
    
    AliPIDResponse* fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
    /*if(!aodEv){
        AliPIDResponse* fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
    }else{
        AliPIDResponse* fPIDResponse = (inputHandlerAOD)->GetPIDResponse();
    }*/
    if(!fPIDResponse){
        AliWarning("NULL PIDResponse");
    }

    if(aodEv){
        trkArr = MakeTrkArr(aod);
    }
    //__________________________End major event stuff_____________________________
    
    
    //Fill the histogram cataloguing # of events vs. events tagged
    if(MBtrg){
        fHistNevents_MB->Fill("Events",1);
    }
    if(EMC7trg){
        fHistNevents_EMC7->Fill("Events",1);
    }
    
    if(EMCJettrg){
        fHistNevents_EMCJet->Fill("Events",1);
    }
    
    //Initialize energy variable and tagging flags
    Double_t PtSum=0;
    tagStrong=kFALSE;
    Bool_t tagEvt=kFALSE;

    //Initialize the # of tracks variable and the Eta Phi arrays
    Int_t ntracks=0;
    ntracks = aod->GetNumberOfTracks();
    if(!aodEv){
        ntracks = esd->GetNumberOfTracks();
    }
    
    std::vector<Double_t> Eta;
    std::vector<Double_t> Phi;
    
    // Track loop for reconstructed event
    for(Int_t i = 0; i < ntracks; i++) {

        tagStrong=kFALSE;
        tagPhot=kFALSE;
        AliESDtrack* esdtrack = 0;
        AliAODTrack* aodtrack = (AliAODTrack*)aod->GetTrack(i); // pointer to reconstructed to track       
        if(!aodEv){
            esdtrack = (AliESDtrack*)esd->GetTrack(i); // pointer to reconstructed to track 
        }

        if(!aodtrack) { 
            AliError(Form("ERROR: Could not retrieve track %d",i)); 
            continue; 
        }
        
        if(!aodEv){
            if(!esdtrack) { 
                AliError(Form("ERROR: Could not retrieve track %d",i)); 
                continue; 
            }
        }

        //Fill TPCOnly track eta-phi

        if(aodtrack->IsTPCOnly()){
            fHistEtaPhiTPCOnly_MB->Fill(aodtrack->Eta(),aodtrack->Phi());
        }

        if(!aodEv){
            if(AliESDtrackCuts::GetStandardTPCOnlyTrackCuts()->AcceptTrack(esdtrack)){
                fHistEtaPhiTPCOnly_MB->Fill(esdtrack->Eta(),esdtrack->Phi());
            }
        }


        //Do hybrid track cuts
        if(!aodtrack->IsHybridGlobalConstrainedGlobal()){continue;}
        if(!aodEv){
            if(!globaltrackCuts->AcceptTrack(esdtrack)&&!comptrackCuts->AcceptTrack(esdtrack)){continue;}
        }

        //Add this tracks energy to the running total
        PtSum=PtSum+aodtrack->Pt();

        //Fill the Eta Phi arrays with this tracks Eta and Phi
        Eta.push_back(aodtrack->Eta());
        Phi.push_back(aodtrack->Phi());

        if(MBtrg){
            fHistEtaPhi_MB->Fill(aodtrack->Eta(),aodtrack->Phi());
        }
        if(EMC7trg){
            fHistEtaPhi_EMC7->Fill(aodtrack->Eta(),aodtrack->Phi());
        }
        if(EMCJettrg){
            fHistEtaPhi_EMCJet->Fill(aodtrack->Eta(),aodtrack->Phi());
        }

        //do Cut level histograms
        if(MBtrg){
            if(aodtrack->GetTPCncls()>0){
                fHistTPCNClus_MB->Fill(aodtrack->GetTPCncls());
            }
            fHistITSNClus_MB->Fill(aodtrack->GetNcls(0));
        }
        if(EMC7trg){
            if(aodtrack->GetTPCncls()>0){
                fHistTPCNClus_EMC7->Fill(aodtrack->GetTPCncls());
            }
            fHistITSNClus_EMC7->Fill(aodtrack->GetNcls(0));
        }

        if(EMCJettrg){
            if(aodtrack->GetTPCncls()>0){
                fHistTPCNClus_EMCJet->Fill(aodtrack->GetTPCncls());
            }
            fHistITSNClus_EMCJet->Fill(aodtrack->GetNcls(0)); 
        }

        if(!aodEv){
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
            if(EMCJettrg){
                fHistEtaPhi_EMCJet->Fill(esdtrack->Eta(),esdtrack->Phi());
            }

            //do Cut level histograms
            if(MBtrg){
                if(esdtrack->GetTPCncls()>0){
                    fHistTPCNClus_MB->Fill(esdtrack->GetTPCncls());
                }
                fHistITSNClus_MB->Fill(esdtrack->GetNcls(0));
            }
            if(EMC7trg){
                if(esdtrack->GetTPCncls()>0){
                    fHistTPCNClus_EMC7->Fill(esdtrack->GetTPCncls());
                }
                fHistITSNClus_EMC7->Fill(esdtrack->GetNcls(0));
            }

            if(EMCJettrg){
                if(esdtrack->GetTPCncls()>0){
                    fHistTPCNClus_EMCJet->Fill(esdtrack->GetTPCncls());
                }
                fHistITSNClus_EMCJet->Fill(esdtrack->GetNcls(0)); 
            }
        }


        //Fill histogram for TPC resolution
        if(!aodEv){
            if(MBtrg){
                fHistTPCSig_MB->Fill(esdtrack->GetTPCsignalSigma());
            }
            if(EMC7trg){
                fHistTPCSig_EMC7->Fill(esdtrack->GetTPCsignalSigma());
            }

            if(EMCJettrg){
                fHistTPCSig_EMCJet->Fill(esdtrack->GetTPCsignalSigma());
            }
        }
        
        //Impact parameter
        Float_t xy;
        Float_t z;
        
        xy=TMath::Sqrt(aodtrack->XAtDCA()*aodtrack->XAtDCA()+aodtrack->YAtDCA()*aodtrack->YAtDCA());
        
        if(!aodEv){
            esdtrack->GetImpactParameters(xy, z);
        }
            
        if(MBtrg){
            fHistImpPar_MB->Fill(xy);
        }
        if(EMC7trg){
            fHistImpPar_EMC7->Fill(xy);
        }

        if(EMCJettrg){
            fHistImpPar_EMCJet->Fill(xy);
        }
        
        FillPhotoElecHistos(aod, aodtrack, fPIDResponse, i);
        
        if(!aodEv){
            FillPhotoElecHistos(esd, esdtrack, fPIDResponse, i);
        }
            
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //If the track doesn't pass the cuts, move on to the next one
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if(trackCutsStrong){
            if(!fTrackCutsStrong->AcceptVTrack((AliVTrack*)aodtrack) || aodtrack->GetTPCsignalN()<80){continue;}
        }else{
            if(!fTrackCutsWeak->AcceptVTrack((AliVTrack*)aodtrack) || aodtrack->GetTPCsignalN()<80){continue;}
        }

        if(!aodEv){
            if(trackCutsStrong){
                if(!fTrackCutsStrong->AcceptTrack(esdtrack) || esdtrack->GetTPCsignalN()<80){continue;}
            }else{
                if(!fTrackCutsWeak->AcceptTrack(esdtrack) || esdtrack->GetTPCsignalN()<80){continue;}
            }
        }
            

        
        //Fill histogram for TPC resolution
        if(!aodEv){
            if(MBtrg){
                fHistTPCSigCut_MB->Fill(esdtrack->GetTPCsignalSigma());
            }
            if(EMC7trg){
                fHistTPCSigCut_EMC7->Fill(esdtrack->GetTPCsignalSigma());
            }

            if(EMCJettrg){
                fHistTPCSigCut_EMCJet->Fill(esdtrack->GetTPCsignalSigma());
            }
        }
        
        FillPIDHistos(aod, aodtrack, fPIDResponse);//Fill PID histos and set "tagStrong" boolean if this track satisfies cuts

        if(!aodEv){
            FillPIDHistos(esd, esdtrack, fPIDResponse);//Fill PID histos and set "tagStrong" boolean if this track satisfies cuts
        }

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

            if(EMCJettrg){
                fHistImpParTag_EMCJet->Fill(xy);
            }

            //Pt distribution
            
            if(MBtrg){
                fHistPtTag_MB->Fill(aodtrack->Pt());
            }
            if(EMC7trg){
                fHistPtTag_EMC7->Fill(aodtrack->Pt());
            }

            if(EMCJettrg){
                fHistPtTag_EMCJet->Fill(aodtrack->Pt());
            }

            if(!aodEv){
                if(MBtrg){
                    fHistPtTag_MB->Fill(esdtrack->Pt());
                }
                if(EMC7trg){
                    fHistPtTag_EMC7->Fill(esdtrack->Pt());
                }

                if(EMCJettrg){
                    fHistPtTag_EMCJet->Fill(esdtrack->Pt());
                }
            }
                
            FillDPhiHistos(aod, aodtrack, i);//Fill DPhi histos
            
            if(!aodEv){
                FillDPhiHistos(esd, esdtrack, i);//Fill DPhi histos
            }   
         

            fPool = fPoolMan->GetEventPool(ntracks, aod->GetPrimaryVertex()->GetZ());
            fPool->PrintInfo();

            if(!fPool){cout<<"No Pool for this event man\n"; continue;}

            if(fPool->IsReady() ){
                FillMEDPhiHistos(aodtrack);
            }
            else{
                cout<<"Pool wasn't ready\n";
            }


            if(tagPhot){
                if(MBtrg){
                    fHistPhotoMismatch_MB->Fill(0);
                }
                if(EMC7trg){
                    fHistPhotoMismatch_EMC7->Fill(0);
                }
                if(EMCJettrg){
                    fHistPhotoMismatch_EMCJet->Fill(0);
                }
            }


        }//end if(tagStrong)

    }//end main track loop
    if(tagEvt){
        fPool->UpdatePool(trkArr);
    }
    //Call function to fill Region histos and pass it int array of IDs for identified electron tracks
    Int_t elecIDsSparse[elecCnt];
    for(Int_t i=0;i<elecCnt;i++){
        elecIDsSparse[i]=elecIDs[i];
    }

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
            if(EMCJettrg){
                fHistEtaPhiTag_EMCJet->Fill(Eta[i],Phi[i]);
            }
        }
    }


    
    // NEW HISTO should be filled before this point, as PostData puts the
    // information for this iteration of the UserExec in the container
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMCJet);
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

            if(EMCJettrg){
                fHistM02_All_EMCJet[i]->Fill(M02, EOP);
                fHistM20_All_EMCJet[i]->Fill(M20, EOP);
            }

            //TPC Plots
            //EMC+TRD cuts
            if(EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]&&elecLikeTRD[0]>TRDcut){

                if(MBtrg){
                    fHistTPC_EMCTRD_MB[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
                if(EMC7trg){
                    fHistTPC_EMCTRD_EMC7[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }

                if(EMCJettrg){
                    fHistTPC_EMCTRD_EMCJet[i]->Fill(esdtrack->Pt(), nSigmaTPC);
                }
            }

            //EMC Plots

            //TPC+TRD cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut){

                if(!applySSCuts){
                    if(MBtrg){
                        fHistEMC_TPCTRD_MB[i]->Fill(EOP);
                    }
                    if(EMC7trg){
                        fHistEMC_TPCTRD_EMC7[i]->Fill(EOP);
                    }

                    if(EMCJettrg){
                        fHistEMC_TPCTRD_EMCJet[i]->Fill(EOP);
                    }
                }else{
                    if(M02<.5&&M02>0&&M20<.3&&M20>0){
                        if(MBtrg){
                            fHistEMC_TPCTRD_MB[i]->Fill(EOP);
                        }
                        if(EMC7trg){
                            fHistEMC_TPCTRD_EMC7[i]->Fill(EOP);
                        }

                        if(EMCJettrg){
                            fHistEMC_TPCTRD_EMCJet[i]->Fill(EOP);
                        }
                    }
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

                if(EMCJettrg){
                    fHistM02_Elec_EMCJet[i]->Fill(M02, EOP);
                    fHistM20_Elec_EMCJet[i]->Fill(M20, EOP);
                }
            }

            //TRD Plots

            //TPC+EMC cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]){

                if(MBtrg){
                    fHistTRD_TPCEMC_MB[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }
                if(EMC7trg){
                    fHistTRD_TPCEMC_EMC7[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
                }

                if(EMCJettrg){
                    fHistTRD_TPCEMC_EMCJet[i]->Fill(esdtrack->Pt(), elecLikeTRD[0]);
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

    if(applySSCuts){

        //An electron candidate is one that passes TPC +-2Sig, TRD>.9, 0.85<E/p<1.15, M02=(0,.5), M20=(0,.3)

        if(esdtrack->Pt()<6){
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut&&EOP<EMCcutHigher[0]&&EOP>EMCcutLower[0]&&M02<.5&&M02>0&&M20<.3&&M20>0){

                tagStrong=kTRUE;

            }

        }

        else{

            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut&&EOP<EMCcutHigher[5]&&EOP>EMCcutLower[5]&&M02<.5&&M02>0&&M20<.3&&M20>0){

                tagStrong=kTRUE;

            }

        }

    }else{//no sscuts

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
    PostData(3, fOutputEMCJet);
    return;
}

void AliAnalysisTaskPSHFE::FillPIDHistos(AliAODEvent *aod, AliAODTrack *aodtrack, AliPIDResponse *fPIDResponse){

    if(!aodtrack){
        AliWarning("aodtrack is null, no point in doing PID");
        return;
    }


    Bool_t isPIDRej = kFALSE;

    //Fill TOF and TPC status variables
    AliPIDResponse::EDetPidStatus TOFStatus=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, aodtrack);

    AliPIDResponse::EDetPidStatus TPCStatus=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, aodtrack);

    AliPIDResponse::EDetPidStatus TRDStatus=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTRD, aodtrack);

    AliPIDResponse::EDetPidStatus EMCStatus=fPIDResponse->CheckPIDStatus(AliPIDResponse::kEMCAL, aodtrack);

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
    nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(aodtrack,AliPID::kElectron);

    Double_t nSigmaTPC;
    nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kElectron);

    Double_t elecLikeTRD[1];
    if(fPIDResponse->ComputeTRDProbability(aodtrack, AliPID::kElectron, elecLikeTRD, AliTRDPIDResponse::kLQ2D) != AliPIDResponse::kDetPidOk || aodtrack->GetTRDntrackletsPID()<4){
        fHistPIDRejection->Fill(3);
        isPIDRej=kTRUE;
    }

    if(isPIDRej){return;}

    //declare emcal cluster PID variables
    Double_t EOP=-1;
    Double_t M02=-1;
    Double_t M20=-1;
    Int_t caloId=aodtrack->GetEMCALcluster();

    if(caloId==-99999){
        return;
    }

    AliAODCaloCluster* tagEMCclus=aod->GetCaloCluster(caloId);

    if(tagEMCclus->E()>.5){
        EOP = tagEMCclus->E()/aodtrack->Pt();
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
        if(aodtrack->Pt()>ptLower[i]&&aodtrack->Pt()<ptUpper[i]){

            //Fill general shower shape plots
            if(MBtrg){
                fHistM02_All_MB[i]->Fill(M02, EOP);
                fHistM20_All_MB[i]->Fill(M20, EOP);
            }
            if(EMC7trg){
                fHistM02_All_EMC7[i]->Fill(M02, EOP);
                fHistM20_All_EMC7[i]->Fill(M20, EOP);
            }

            if(EMCJettrg){
                fHistM02_All_EMCJet[i]->Fill(M02, EOP);
                fHistM20_All_EMCJet[i]->Fill(M20, EOP);
            }

            //TPC Plots
            //EMC+TRD cuts
            if(EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]&&elecLikeTRD[0]>TRDcut){

                if(MBtrg){
                    fHistTPC_EMCTRD_MB[i]->Fill(aodtrack->Pt(), nSigmaTPC);
                }
                if(EMC7trg){
                    fHistTPC_EMCTRD_EMC7[i]->Fill(aodtrack->Pt(), nSigmaTPC);
                }

                if(EMCJettrg){
                    fHistTPC_EMCTRD_EMCJet[i]->Fill(aodtrack->Pt(), nSigmaTPC);
                }
            }

            //EMC Plots

            //TPC+TRD cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut){

                if(!applySSCuts){
                    if(MBtrg){
                        fHistEMC_TPCTRD_MB[i]->Fill(EOP);
                    }
                    if(EMC7trg){
                        fHistEMC_TPCTRD_EMC7[i]->Fill(EOP);
                    }

                    if(EMCJettrg){
                        fHistEMC_TPCTRD_EMCJet[i]->Fill(EOP);
                    }
                }else{
                    if(M02<.5&&M02>0&&M20<.3&&M20>0){
                        if(MBtrg){
                            fHistEMC_TPCTRD_MB[i]->Fill(EOP);
                        }
                        if(EMC7trg){
                            fHistEMC_TPCTRD_EMC7[i]->Fill(EOP);
                        }

                        if(EMCJettrg){
                            fHistEMC_TPCTRD_EMCJet[i]->Fill(EOP);
                        }
                    }
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

                if(EMCJettrg){
                    fHistM02_Elec_EMCJet[i]->Fill(M02, EOP);
                    fHistM20_Elec_EMCJet[i]->Fill(M20, EOP);
                }
            }

            //TRD Plots

            //TPC+EMC cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]){

                if(MBtrg){
                    fHistTRD_TPCEMC_MB[i]->Fill(aodtrack->Pt(), elecLikeTRD[0]);
                }
                if(EMC7trg){
                    fHistTRD_TPCEMC_EMC7[i]->Fill(aodtrack->Pt(), elecLikeTRD[0]);
                }

                if(EMCJettrg){
                    fHistTRD_TPCEMC_EMCJet[i]->Fill(aodtrack->Pt(), elecLikeTRD[0]);
                }
            }
        }
    }
    if(MBtrg){
        if(aodtrack->Pt()<2&&aodtrack->Pt()>1){
            if(nSigmaTPC<-2&&nSigmaTPC>-8){
                fHistEMC_Had_MB_1Gev->Fill(EOP);
            }
        }
    }

    if(applySSCuts){

        //An electron candidate is one that passes TPC +-2Sig, TRD>.9, 0.85<E/p<1.15, M02=(0,.5), M20=(0,.3)

        if(aodtrack->Pt()<6){
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut&&EOP<EMCcutHigher[0]&&EOP>EMCcutLower[0]&&M02<.5&&M02>0&&M20<.3&&M20>0){

                tagStrong=kTRUE;

            }

        }

        else{

            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut&&EOP<EMCcutHigher[5]&&EOP>EMCcutLower[5]&&M02<.5&&M02>0&&M20<.3&&M20>0){

                tagStrong=kTRUE;

            }

        }

    }else{//no sscuts

        //An electron candidate is one that passes TPC +-2Sig, TRD>.9, 0.85<E/p<1.15

        if(aodtrack->Pt()<6){
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut&&EOP<EMCcutHigher[0]&&EOP>EMCcutLower[0]){

                tagStrong=kTRUE;

            }

        }

        else{

            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut&&EOP<EMCcutHigher[5]&&EOP>EMCcutLower[5]){

                tagStrong=kTRUE;

            }

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
    PostData(3, fOutputEMCJet);
    return;
}

void AliAnalysisTaskPSHFE::FillDPhiHistos(AliESDEvent *esd, AliESDtrack *esdtrack, Int_t i){
    //Run through all tracks in this 
    Int_t ntracks = esd->GetNumberOfTracks();
    for(Int_t j = 0; j < ntracks; j++) {

        //Don't double count the tagged tracks
        if(i==j){continue;}


        AliESDtrack* esdtrackassoc = (AliESDtrack*)esd->GetTrack(j); // pointer to reconstructed to track          

        if(!esdtrackassoc) { 
            AliError(Form("ERROR: Could not retrieve esdtrack %d",j)); 
            continue; 
        }

        //Do hybrid track cuts
        if(!globaltrackCuts->AcceptTrack(esdtrackassoc)&&!comptrackCuts->AcceptTrack(esdtrackassoc)){continue;}

        //Pt distribution
        if(MBtrg){
            fHistPtAssoc_MB->Fill(esdtrackassoc->Pt());
        }
        if(EMC7trg){
            fHistPtAssoc_EMC7->Fill(esdtrackassoc->Pt());
        }

        if(EMCJettrg){
            fHistPtAssoc_EMCJet->Fill(esdtrackassoc->Pt());
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
                    fHistDPhi28dEdx_MB->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
                if(EMC7trg){
                    fHistDPhi2_3_EMC7[1]->Fill(DPhi);
                    fHistDPhi28_EMC7->Fill(DPhi);
                    fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMC7->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
                if(EMCJettrg){
                    fHistDPhi2_3_EMCJet[1]->Fill(DPhi);
                    fHistDPhi28_EMCJet->Fill(DPhi);
                    fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMCJet->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
            }

            //3GeV
            if(esdtrackassoc->Pt()>3&&esdtrackassoc->Pt()<4){
                if(MBtrg){
                    fHistDPhi3_4_MB[1]->Fill(DPhi);
                    fHistDPhi28_MB->Fill(DPhi);
                    fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_MB->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
                if(EMC7trg){
                    fHistDPhi3_4_EMC7[1]->Fill(DPhi);
                    fHistDPhi28_EMC7->Fill(DPhi);
                    fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMC7->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
                if(EMCJettrg){
                    fHistDPhi3_4_EMCJet[1]->Fill(DPhi);
                    fHistDPhi28_EMCJet->Fill(DPhi);
                    fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMCJet->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
            }

            //5GeV
            if(esdtrackassoc->Pt()>4){
                if(MBtrg){
                    fHistDPhi4_MB[1]->Fill(DPhi);
                    fHistDPhi28_MB->Fill(DPhi);
                    fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_MB->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
                if(EMC7trg){
                    fHistDPhi4_EMC7[1]->Fill(DPhi);
                    fHistDPhi28_EMC7->Fill(DPhi);
                    fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMC7->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
                if(EMCJettrg){
                    fHistDPhi4_EMCJet[1]->Fill(DPhi);
                    fHistDPhi28_EMCJet->Fill(DPhi);
                    fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMCJet->Fill(DPhi, esdtrackassoc->GetTPCsignal());
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
                    fHistDPhi28dEdx_MB->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
                if(EMC7trg){
                    fHistDPhi2_3_EMC7[2]->Fill(DPhi);
                    fHistDPhi28_EMC7->Fill(DPhi);
                    fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMC7->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
                if(EMCJettrg){
                    fHistDPhi2_3_EMCJet[2]->Fill(DPhi);
                    fHistDPhi28_EMCJet->Fill(DPhi);
                    fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMCJet->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
            }

            //3GeV
            if(esdtrackassoc->Pt()>3&&esdtrackassoc->Pt()<4){
                if(MBtrg){
                    fHistDPhi3_4_MB[2]->Fill(DPhi);
                    fHistDPhi28_MB->Fill(DPhi);
                    fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_MB->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
                if(EMC7trg){
                    fHistDPhi3_4_EMC7[2]->Fill(DPhi);
                    fHistDPhi28_EMC7->Fill(DPhi);
                    fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMC7->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
                if(EMCJettrg){
                    fHistDPhi3_4_EMCJet[2]->Fill(DPhi);
                    fHistDPhi28_EMCJet->Fill(DPhi);
                    fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMCJet->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
            }

            //5GeV
            if(esdtrackassoc->Pt()>4){
                if(MBtrg){
                    fHistDPhi4_MB[2]->Fill(DPhi);
                    fHistDPhi28_MB->Fill(DPhi);
                    fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_MB->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
                if(EMC7trg){
                    fHistDPhi4_EMC7[2]->Fill(DPhi);
                    fHistDPhi28_EMC7->Fill(DPhi);
                    fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMC7->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
                if(EMCJettrg){
                    fHistDPhi4_EMCJet[2]->Fill(DPhi);
                    fHistDPhi28_EMCJet->Fill(DPhi);
                    fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMCJet->Fill(DPhi, esdtrackassoc->GetTPCsignal());
                }
            }
        }

    }
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMCJet);
    return;
}

void AliAnalysisTaskPSHFE::FillDPhiHistos(AliAODEvent *aod, AliAODTrack *aodtrack, Int_t i){
    //Run through all tracks in this 
    Int_t ntracks = aod->GetNumberOfTracks();
    for(Int_t j = 0; j < ntracks; j++) {

        //Don't double count the tagged tracks
        if(i==j){continue;}


        AliAODTrack* aodtrackassoc = (AliAODTrack*)aod->GetTrack(j); // pointer to reconstructed to track          

        if(!aodtrackassoc) { 
            AliError(Form("ERROR: Could not retrieve aodtrack %d",j)); 
            continue; 
        }

        //Do hybrid track cuts
        if(aodtrackassoc->IsHybridGlobalConstrainedGlobal()){continue;}

        //Pt distribution
        if(MBtrg){
            fHistPtAssoc_MB->Fill(aodtrackassoc->Pt());
        }
        if(EMC7trg){
            fHistPtAssoc_EMC7->Fill(aodtrackassoc->Pt());
        }

        if(EMCJettrg){
            fHistPtAssoc_EMCJet->Fill(aodtrackassoc->Pt());
        }

        //Fill Delta Phi variable and correct for periodicity
        Double_t DPhi=aodtrackassoc->Phi()-aodtrack->Phi();

        if(DPhi<-TMath::Pi()/2){DPhi=TMath::Abs(2*TMath::Pi()+DPhi);}

        if(DPhi>3*TMath::Pi()/2){DPhi=-TMath::Abs(2*TMath::Pi()-DPhi);}

        Double_t DEta=aodtrackassoc->Eta()-aodtrack->Eta();

        //candidate 1<pt<2
        if(aodtrack->Pt()>1&&aodtrack->Pt()<2){

            //300-500MeV
            if(aodtrackassoc->Pt()>.3&&aodtrackassoc->Pt()<.5){
                if(MBtrg){
                    fHistDPhi300_500_MB[0]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi300_500_EMC7[0]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi300_500_EMCJet[0]->Fill(DPhi);
                }
            }

            //500MeV
            if(aodtrackassoc->Pt()>.5&&aodtrackassoc->Pt()<.8){
                if(MBtrg){
                    fHistDPhi500_800_MB[0]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi500_800_EMC7[0]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi500_800_EMCJet[0]->Fill(DPhi);
                }
            }

            //800MeV
            if(aodtrackassoc->Pt()>.8&&aodtrackassoc->Pt()<1){
                if(MBtrg){
                    fHistDPhi800_1_MB[0]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi800_1_EMC7[0]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi800_1_EMCJet[0]->Fill(DPhi);
                }
            }

            //1GeV
            if(aodtrackassoc->Pt()>1&&aodtrackassoc->Pt()<2){
                if(MBtrg){
                    fHistDPhi1_2_MB[0]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi1_2_EMC7[0]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi1_2_EMCJet[0]->Fill(DPhi);
                }
            }

            //2GeV
            if(aodtrackassoc->Pt()>2&&aodtrackassoc->Pt()<3){
                if(MBtrg){
                    fHistDPhi2_3_MB[0]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi2_3_EMC7[0]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi2_3_EMCJet[0]->Fill(DPhi);
                }
            }

            //3GeV
            if(aodtrackassoc->Pt()>3&&aodtrackassoc->Pt()<4){
                if(MBtrg){
                    fHistDPhi3_4_MB[0]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi3_4_EMC7[0]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi3_4_EMCJet[0]->Fill(DPhi);
                }
            }

            //5GeV
            if(aodtrackassoc->Pt()>4){
                if(MBtrg){
                    fHistDPhi4_MB[0]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi4_EMC7[0]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi4_EMCJet[0]->Fill(DPhi);
                }
            }
        }

        //candidate 2<pt<4
        if(aodtrack->Pt()>2&&aodtrack->Pt()<4){

            //300MeV
            if(aodtrackassoc->Pt()>.3&&aodtrackassoc->Pt()<.5){
                if(MBtrg){
                    fHistDPhi300_500_MB[1]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi300_500_EMC7[1]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi300_500_EMCJet[1]->Fill(DPhi);
                }
            }

            //500MeV
            if(aodtrackassoc->Pt()>.5&&aodtrackassoc->Pt()<.8){
                if(MBtrg){
                    fHistDPhi500_800_MB[1]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi500_800_EMC7[1]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi500_800_EMCJet[1]->Fill(DPhi);
                }
            }

            //800MeV
            if(aodtrackassoc->Pt()>.8&&aodtrackassoc->Pt()<1){
                if(MBtrg){
                    fHistDPhi800_1_MB[1]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi800_1_EMC7[1]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi800_1_EMCJet[1]->Fill(DPhi);
                }
            }

            //1GeV
            if(aodtrackassoc->Pt()>1&&aodtrackassoc->Pt()<2){
                if(MBtrg){
                    fHistDPhi1_2_MB[1]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi1_2_EMC7[1]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi1_2_EMCJet[1]->Fill(DPhi);
                }
            }

            //2GeV
            if(aodtrackassoc->Pt()>2&&aodtrackassoc->Pt()<3){
                if(MBtrg){
                    fHistDPhi2_3_MB[1]->Fill(DPhi);
                    fHistDPhi28_MB->Fill(DPhi);
                    fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_MB->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
                if(EMC7trg){
                    fHistDPhi2_3_EMC7[1]->Fill(DPhi);
                    fHistDPhi28_EMC7->Fill(DPhi);
                    fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMC7->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
                if(EMCJettrg){
                    fHistDPhi2_3_EMCJet[1]->Fill(DPhi);
                    fHistDPhi28_EMCJet->Fill(DPhi);
                    fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMCJet->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
            }

            //3GeV
            if(aodtrackassoc->Pt()>3&&aodtrackassoc->Pt()<4){
                if(MBtrg){
                    fHistDPhi3_4_MB[1]->Fill(DPhi);
                    fHistDPhi28_MB->Fill(DPhi);
                    fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_MB->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
                if(EMC7trg){
                    fHistDPhi3_4_EMC7[1]->Fill(DPhi);
                    fHistDPhi28_EMC7->Fill(DPhi);
                    fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMC7->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
                if(EMCJettrg){
                    fHistDPhi3_4_EMCJet[1]->Fill(DPhi);
                    fHistDPhi28_EMCJet->Fill(DPhi);
                    fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMCJet->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
            }

            //5GeV
            if(aodtrackassoc->Pt()>4){
                if(MBtrg){
                    fHistDPhi4_MB[1]->Fill(DPhi);
                    fHistDPhi28_MB->Fill(DPhi);
                    fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_MB->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
                if(EMC7trg){
                    fHistDPhi4_EMC7[1]->Fill(DPhi);
                    fHistDPhi28_EMC7->Fill(DPhi);
                    fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMC7->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
                if(EMCJettrg){
                    fHistDPhi4_EMCJet[1]->Fill(DPhi);
                    fHistDPhi28_EMCJet->Fill(DPhi);
                    fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMCJet->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
            }                    
        }

        //candidate 4<pt<8
        if(aodtrack->Pt()>4&&aodtrack->Pt()<8){

            //300MeV
            if(aodtrackassoc->Pt()>.3&&aodtrackassoc->Pt()<.5){
                if(MBtrg){
                    fHistDPhi300_500_MB[2]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi300_500_EMC7[2]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi300_500_EMCJet[2]->Fill(DPhi);
                }
            }

            //500MeV
            if(aodtrackassoc->Pt()>.5&&aodtrackassoc->Pt()<.8){
                if(MBtrg){
                    fHistDPhi500_800_MB[2]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi500_800_EMC7[2]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi500_800_EMCJet[2]->Fill(DPhi);
                }
            }

            //800MeV
            if(aodtrackassoc->Pt()>.8&&aodtrackassoc->Pt()<1){
                if(MBtrg){
                    fHistDPhi800_1_MB[2]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi800_1_EMC7[2]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi800_1_EMCJet[2]->Fill(DPhi);
                }
            }

            //1GeV
            if(aodtrackassoc->Pt()>1&&aodtrackassoc->Pt()<2){
                if(MBtrg){
                    fHistDPhi1_2_MB[2]->Fill(DPhi);
                }
                if(EMC7trg){
                    fHistDPhi1_2_EMC7[2]->Fill(DPhi);
                }
                if(EMCJettrg){
                    fHistDPhi1_2_EMCJet[2]->Fill(DPhi);
                }
            }

            //2GeV
            if(aodtrackassoc->Pt()>2&&aodtrackassoc->Pt()<3){
                if(MBtrg){
                    fHistDPhi2_3_MB[2]->Fill(DPhi);
                    fHistDPhi28_MB->Fill(DPhi);
                    fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_MB->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
                if(EMC7trg){
                    fHistDPhi2_3_EMC7[2]->Fill(DPhi);
                    fHistDPhi28_EMC7->Fill(DPhi);
                    fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMC7->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
                if(EMCJettrg){
                    fHistDPhi2_3_EMCJet[2]->Fill(DPhi);
                    fHistDPhi28_EMCJet->Fill(DPhi);
                    fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMCJet->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
            }

            //3GeV
            if(aodtrackassoc->Pt()>3&&aodtrackassoc->Pt()<4){
                if(MBtrg){
                    fHistDPhi3_4_MB[2]->Fill(DPhi);
                    fHistDPhi28_MB->Fill(DPhi);
                    fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_MB->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
                if(EMC7trg){
                    fHistDPhi3_4_EMC7[2]->Fill(DPhi);
                    fHistDPhi28_EMC7->Fill(DPhi);
                    fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMC7->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
                if(EMCJettrg){
                    fHistDPhi3_4_EMCJet[2]->Fill(DPhi);
                    fHistDPhi28_EMCJet->Fill(DPhi);
                    fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMCJet->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
            }

            //5GeV
            if(aodtrackassoc->Pt()>4){
                if(MBtrg){
                    fHistDPhi4_MB[2]->Fill(DPhi);
                    fHistDPhi28_MB->Fill(DPhi);
                    fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_MB->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
                if(EMC7trg){
                    fHistDPhi4_EMC7[2]->Fill(DPhi);
                    fHistDPhi28_EMC7->Fill(DPhi);
                    fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMC7->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
                if(EMCJettrg){
                    fHistDPhi4_EMCJet[2]->Fill(DPhi);
                    fHistDPhi28_EMCJet->Fill(DPhi);
                    fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                    fHistDPhi28dEdx_EMCJet->Fill(DPhi, aodtrackassoc->GetTPCsignal());
                }
            }
        }

    }
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMCJet);
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

void AliAnalysisTaskPSHFE::SetSSCutBool(Bool_t SSCutBool){
    applySSCuts=SSCutBool;
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

        AliESDtrack* esdtrackassoc = (AliESDtrack*)esd->GetTrack(j); // pointer to reconstructed to track          

        if(!esdtrackassoc) { 
            AliError(Form("ERROR: Could not retrieve esdtrack %d",j)); 
            continue; 
        }

        //Do hybrid track cuts
        if(!globaltrackCuts->AcceptTrack(esdtrackassoc)&&!comptrackCuts->AcceptTrack(esdtrackassoc)){continue;}

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

            if(EMCJettrg){
                fHistInvMassElecLike_EMCJet->Fill(InvMass);
                fHistOpAngElecLike_EMCJet->Fill(OpAng);
            }

        }else{
            if(InvMass<0.1&&OpAng<0.1){tagPhot=kTRUE;}
            if(MBtrg){
                fHistInvMassElecUnLike_MB->Fill(InvMass);
                fHistOpAngElecUnLike_MB->Fill(OpAng);
            }

            if(EMC7trg){
                fHistInvMassElecUnLike_EMC7->Fill(InvMass);
                fHistOpAngElecUnLike_EMC7->Fill(OpAng);
            }

            if(EMCJettrg){
                fHistInvMassElecUnLike_EMCJet->Fill(InvMass);
                fHistOpAngElecUnLike_EMCJet->Fill(OpAng);
            }
        }
    }//end loop over tracks

    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMCJet);
    return;
}

void AliAnalysisTaskPSHFE::FillPhotoElecHistos(AliAODEvent *aod, AliAODTrack *aodtrack, AliPIDResponse *fPIDResponse, Int_t i){

    if(!aodtrack){
        AliWarning("aodtrack is null, no point in doing Photonic Electron stuff");
        return;
    }

    Bool_t isElec=kFALSE;
    Double_t ElecMass=.0005109989;    
    //Fill TOF and TPC status variables  
    AliPIDResponse::EDetPidStatus TPCStatus=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, aodtrack);

    AliPIDResponse::EDetPidStatus TRDStatus=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTRD, aodtrack);

    AliPIDResponse::EDetPidStatus EMCStatus=fPIDResponse->CheckPIDStatus(AliPIDResponse::kEMCAL, aodtrack);

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
    nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kElectron);

    Double_t elecLikeTRD[1];
    if(fPIDResponse->ComputeTRDProbability(aodtrack, AliPID::kElectron, elecLikeTRD, AliTRDPIDResponse::kLQ2D) != AliPIDResponse::kDetPidOk || aodtrack->GetTRDntrackletsPID()<4){
        return;
    }

    //declare emcal cluster PID variables
    Double_t EOP=-1;
    Int_t caloId=aodtrack->GetEMCALcluster();

    if(caloId==-99999){
        return;
    }

    AliAODCaloCluster* tagEMCclus=aod->GetCaloCluster(caloId);

    EOP = tagEMCclus->E()/aodtrack->Pt();

    if((nSigmaTPC<2&&nSigmaTPC>-2)||(EOP<1.4&&EOP>.8)||(elecLikeTRD[0]>.8)){
        isElec=kTRUE;
    }

    if(!isElec){return;}

    Int_t ntracks=aod->GetNumberOfTracks();

    for(Int_t j = 0; j < ntracks; j++) {

        //Don't double count the tagged tracks
        if(i==j){continue;}

        AliAODTrack* aodtrackassoc = (AliAODTrack*)aod->GetTrack(j); // pointer to reconstructed to track          

        if(!aodtrackassoc) { 
            AliError(Form("ERROR: Could not retrieve aodtrack %d",j)); 
            continue; 
        }

        //Do hybrid track cuts
        if(aodtrackassoc->IsHybridGlobalConstrainedGlobal()){continue;}

        Bool_t isElecToo=kFALSE;

        //Fill TOF and TPC status variables  
        AliPIDResponse::EDetPidStatus TPCStatusassoc=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, aodtrackassoc);

        AliPIDResponse::EDetPidStatus TRDStatusassoc=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTRD, aodtrackassoc);

        //Check validity of PID, TODO: Add a rejection histogram
        if(TPCStatusassoc!=AliPIDResponse::kDetPidOk){
            continue;
        }

        if(TRDStatusassoc!=AliPIDResponse::kDetPidOk){
            continue;
        }


        //Get the # of sigmas around an electron hypothesis for TOF and TPC
        Double_t nSigmaTPCassoc;
        nSigmaTPCassoc = fPIDResponse->NumberOfSigmasTPC(aodtrackassoc,AliPID::kElectron);

        Double_t elecLikeTRDassoc[1];
        if(fPIDResponse->ComputeTRDProbability(aodtrackassoc, AliPID::kElectron, elecLikeTRDassoc, AliTRDPIDResponse::kLQ2D) != AliPIDResponse::kDetPidOk || aodtrackassoc->GetTRDntrackletsPID()<4){
            return;
        }

        if((nSigmaTPCassoc<2&&nSigmaTPCassoc>-2)||(elecLikeTRDassoc[0]>.8)){
            isElecToo=kTRUE;
        }

        if(!isElecToo){continue;}

        Double_t elecE1=TMath::Sqrt(aodtrack->P()*aodtrack->P()+ElecMass*ElecMass);
        Double_t elecE2=TMath::Sqrt(aodtrackassoc->P()*aodtrackassoc->P()+ElecMass*ElecMass);

        TLorentzVector elec1(aodtrack->Px(), aodtrack->Py(), aodtrack->Pz(), elecE1);
        TLorentzVector elec2(aodtrackassoc->Px(), aodtrackassoc->Py(), aodtrackassoc->Pz(), elecE2);

        Double_t InvMass=(elec1+elec2).M();
        Double_t OpAng=elec1.Angle(elec2.Vect());

        if(aodtrack->GetSign()==aodtrackassoc->GetSign()){

            if(MBtrg){
                fHistInvMassElecLike_MB->Fill(InvMass);
                fHistOpAngElecLike_MB->Fill(OpAng);
            }

            if(EMC7trg){
                fHistInvMassElecLike_EMC7->Fill(InvMass);
                fHistOpAngElecLike_EMC7->Fill(OpAng);
            }

            if(EMCJettrg){
                fHistInvMassElecLike_EMCJet->Fill(InvMass);
                fHistOpAngElecLike_EMCJet->Fill(OpAng);
            }

        }else{
            if(InvMass<0.1&&OpAng<0.1){tagPhot=kTRUE;}
            if(MBtrg){
                fHistInvMassElecUnLike_MB->Fill(InvMass);
                fHistOpAngElecUnLike_MB->Fill(OpAng);
            }

            if(EMC7trg){
                fHistInvMassElecUnLike_EMC7->Fill(InvMass);
                fHistOpAngElecUnLike_EMC7->Fill(OpAng);
            }

            if(EMCJettrg){
                fHistInvMassElecUnLike_EMCJet->Fill(InvMass);
                fHistOpAngElecUnLike_EMCJet->Fill(OpAng);
            }
        }
    }//end loop over tracks

    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMCJet);
    return;
}

TObjArray* AliAnalysisTaskPSHFE::MakeTrkArr(AliAODEvent *aod)
{
    if(!aod){AliWarning("Invalid AOD Event");}
    Int_t nTracks = aod->GetNumberOfTracks();
    TObjArray* accTracks = new TObjArray;

    for(Int_t i=0;i<nTracks;i++){
        AliAODTrack *aodtrack = (AliAODTrack*)aod->GetTrack(i);

        if(!aodtrack){
            continue;
        }

        if(aodtrack->IsHybridGlobalConstrainedGlobal()){
            continue;
        }

        accTracks->Add(new AliAODTrack(*aodtrack));
    }

    return accTracks;
}

void AliAnalysisTaskPSHFE::FillMEDPhiHistos(AliAODTrack *aodtrack)
{
    Int_t nEvents = fPool->GetCurrentNEvents();

    for(Int_t Ev=0;Ev<nEvents;Ev++){

        TObjArray *mixedTracks = fPool->GetEvent(Ev);

        Int_t nMixedTracks = mixedTracks->GetEntriesFast();

        for(Int_t j=0;j<nMixedTracks;j++){
            AliAODTrack* aodtrackassoc = (AliAODTrack*)mixedTracks->At(j); // pointer to reconstructed to track          

            if(!aodtrackassoc) { 
                AliError(Form("ERROR: Could not retrieve aodtrack %d",j)); 
                continue; 
            }

            //Do hybrid track cuts
            if(aodtrackassoc->IsHybridGlobalConstrainedGlobal()){continue;}

            //Pt distribution
            if(MBtrg){
                fHistPtAssocMix_MB->Fill(aodtrack->Pt());
            }
            if(EMC7trg){
                fHistPtAssocMix_EMC7->Fill(aodtrack->Pt());
            }

            if(EMCJettrg){
                fHistPtAssocMix_EMCJet->Fill(aodtrack->Pt());
            }

            //Fill Delta Phi variable and correct for periodicity
            Double_t DPhi=aodtrackassoc->Phi()-aodtrack->Phi();

            if(DPhi<-TMath::Pi()/2){DPhi=TMath::Abs(2*TMath::Pi()+DPhi);}

            if(DPhi>3*TMath::Pi()/2){DPhi=-TMath::Abs(2*TMath::Pi()-DPhi);}

            Double_t DEta=aodtrackassoc->Eta()-aodtrack->Eta();

            //candidate 1<pt<2
            if(aodtrack->Pt()>1&&aodtrack->Pt()<2){

                //300-500MeV
                if(aodtrackassoc->Pt()>.3&&aodtrackassoc->Pt()<.5){
                    if(MBtrg){
                        fHistDPhiMix300_500_MB[0]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix300_500_EMC7[0]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix300_500_EMCJet[0]->Fill(DPhi);
                    }
                }

                //500MeV
                if(aodtrackassoc->Pt()>.5&&aodtrackassoc->Pt()<.8){
                    if(MBtrg){
                        fHistDPhiMix500_800_MB[0]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix500_800_EMC7[0]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix500_800_EMCJet[0]->Fill(DPhi);
                    }
                }

                //800MeV
                if(aodtrackassoc->Pt()>.8&&aodtrackassoc->Pt()<1){
                    if(MBtrg){
                        fHistDPhiMix800_1_MB[0]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix800_1_EMC7[0]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix800_1_EMCJet[0]->Fill(DPhi);
                    }
                }

                //1GeV
                if(aodtrackassoc->Pt()>1&&aodtrackassoc->Pt()<2){
                    if(MBtrg){
                        fHistDPhiMix1_2_MB[0]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix1_2_EMC7[0]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix1_2_EMCJet[0]->Fill(DPhi);
                    }
                }

                //2GeV
                if(aodtrackassoc->Pt()>2&&aodtrackassoc->Pt()<3){
                    if(MBtrg){
                        fHistDPhiMix2_3_MB[0]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix2_3_EMC7[0]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix2_3_EMCJet[0]->Fill(DPhi);
                    }
                }

                //3GeV
                if(aodtrackassoc->Pt()>3&&aodtrackassoc->Pt()<4){
                    if(MBtrg){
                        fHistDPhiMix3_4_MB[0]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix3_4_EMC7[0]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix3_4_EMCJet[0]->Fill(DPhi);
                    }
                }

                //5GeV
                if(aodtrackassoc->Pt()>4){
                    if(MBtrg){
                        fHistDPhiMix4_MB[0]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix4_EMC7[0]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix4_EMCJet[0]->Fill(DPhi);
                    }
                }
            }

            //candidate 2<pt<4
            if(aodtrack->Pt()>2&&aodtrack->Pt()<4){

                //300MeV
                if(aodtrackassoc->Pt()>.3&&aodtrackassoc->Pt()<.5){
                    if(MBtrg){
                        fHistDPhiMix300_500_MB[1]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix300_500_EMC7[1]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix300_500_EMCJet[1]->Fill(DPhi);
                    }
                }

                //500MeV
                if(aodtrackassoc->Pt()>.5&&aodtrackassoc->Pt()<.8){
                    if(MBtrg){
                        fHistDPhiMix500_800_MB[1]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix500_800_EMC7[1]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix500_800_EMCJet[1]->Fill(DPhi);
                    }
                }

                //800MeV
                if(aodtrackassoc->Pt()>.8&&aodtrackassoc->Pt()<1){
                    if(MBtrg){
                        fHistDPhiMix800_1_MB[1]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix800_1_EMC7[1]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix800_1_EMCJet[1]->Fill(DPhi);
                    }
                }

                //1GeV
                if(aodtrackassoc->Pt()>1&&aodtrackassoc->Pt()<2){
                    if(MBtrg){
                        fHistDPhiMix1_2_MB[1]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix1_2_EMC7[1]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix1_2_EMCJet[1]->Fill(DPhi);
                    }
                }

                //2GeV
                if(aodtrackassoc->Pt()>2&&aodtrackassoc->Pt()<3){
                    if(MBtrg){
                        fHistDPhiMix2_3_MB[1]->Fill(DPhi);
                        fHistDPhiMix28_MB->Fill(DPhi);
                        fHistDPhiDEtaMix28_MB->Fill(DPhi, DEta);
                    }
                    if(EMC7trg){
                        fHistDPhiMix2_3_EMC7[1]->Fill(DPhi);
                        fHistDPhiMix28_EMC7->Fill(DPhi);
                        fHistDPhiDEtaMix28_EMC7->Fill(DPhi, DEta);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix2_3_EMCJet[1]->Fill(DPhi);
                        fHistDPhiMix28_EMCJet->Fill(DPhi);
                        fHistDPhiDEtaMix28_EMCJet->Fill(DPhi, DEta);
                    }
                }

                //3GeV
                if(aodtrackassoc->Pt()>3&&aodtrackassoc->Pt()<4){
                    if(MBtrg){
                        fHistDPhiMix3_4_MB[1]->Fill(DPhi);
                        fHistDPhiMix28_MB->Fill(DPhi);
                        fHistDPhiDEtaMix28_MB->Fill(DPhi, DEta);
                    }
                    if(EMC7trg){
                        fHistDPhiMix3_4_EMC7[1]->Fill(DPhi);
                        fHistDPhiMix28_EMC7->Fill(DPhi);
                        fHistDPhiDEtaMix28_EMC7->Fill(DPhi, DEta);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix3_4_EMCJet[1]->Fill(DPhi);
                        fHistDPhiMix28_EMCJet->Fill(DPhi);
                        fHistDPhiDEtaMix28_EMCJet->Fill(DPhi, DEta);
                    }
                }

                //5GeV
                if(aodtrackassoc->Pt()>4){
                    if(MBtrg){
                        fHistDPhiMix4_MB[1]->Fill(DPhi);
                        fHistDPhiMix28_MB->Fill(DPhi);
                        fHistDPhiDEtaMix28_MB->Fill(DPhi, DEta);
                    }
                    if(EMC7trg){
                        fHistDPhiMix4_EMC7[1]->Fill(DPhi);
                        fHistDPhiMix28_EMC7->Fill(DPhi);
                        fHistDPhiDEtaMix28_EMC7->Fill(DPhi, DEta);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix4_EMCJet[1]->Fill(DPhi);
                        fHistDPhiMix28_EMCJet->Fill(DPhi);
                        fHistDPhiDEtaMix28_EMCJet->Fill(DPhi, DEta);
                    }
                }                    
            }

            //candidate 4<pt<8
            if(aodtrack->Pt()>4&&aodtrack->Pt()<8){

                //300MeV
                if(aodtrackassoc->Pt()>.3&&aodtrackassoc->Pt()<.5){
                    if(MBtrg){
                        fHistDPhiMix300_500_MB[2]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix300_500_EMC7[2]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix300_500_EMCJet[2]->Fill(DPhi);
                    }
                }

                //500MeV
                if(aodtrackassoc->Pt()>.5&&aodtrackassoc->Pt()<.8){
                    if(MBtrg){
                        fHistDPhiMix500_800_MB[2]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix500_800_EMC7[2]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix500_800_EMCJet[2]->Fill(DPhi);
                    }
                }

                //800MeV
                if(aodtrackassoc->Pt()>.8&&aodtrackassoc->Pt()<1){
                    if(MBtrg){
                        fHistDPhiMix800_1_MB[2]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix800_1_EMC7[2]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix800_1_EMCJet[2]->Fill(DPhi);
                    }
                }

                //1GeV
                if(aodtrackassoc->Pt()>1&&aodtrackassoc->Pt()<2){
                    if(MBtrg){
                        fHistDPhiMix1_2_MB[2]->Fill(DPhi);
                    }
                    if(EMC7trg){
                        fHistDPhiMix1_2_EMC7[2]->Fill(DPhi);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix1_2_EMCJet[2]->Fill(DPhi);
                    }
                }

                //2GeV
                if(aodtrackassoc->Pt()>2&&aodtrackassoc->Pt()<3){
                    if(MBtrg){
                        fHistDPhiMix2_3_MB[2]->Fill(DPhi);
                        fHistDPhiMix28_MB->Fill(DPhi);
                        fHistDPhiDEtaMix28_MB->Fill(DPhi, DEta);
                    }
                    if(EMC7trg){
                        fHistDPhiMix2_3_EMC7[2]->Fill(DPhi);
                        fHistDPhiMix28_EMC7->Fill(DPhi);
                        fHistDPhiDEtaMix28_EMC7->Fill(DPhi, DEta);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix2_3_EMCJet[2]->Fill(DPhi);
                        fHistDPhiMix28_EMCJet->Fill(DPhi);
                        fHistDPhiDEtaMix28_EMCJet->Fill(DPhi, DEta);
                    }
                }

                //3GeV
                if(aodtrackassoc->Pt()>3&&aodtrackassoc->Pt()<4){
                    if(MBtrg){
                        fHistDPhiMix3_4_MB[2]->Fill(DPhi);
                        fHistDPhiMix28_MB->Fill(DPhi);
                        fHistDPhiDEtaMix28_MB->Fill(DPhi, DEta);
                    }
                    if(EMC7trg){
                        fHistDPhiMix3_4_EMC7[2]->Fill(DPhi);
                        fHistDPhiMix28_EMC7->Fill(DPhi);
                        fHistDPhiDEtaMix28_EMC7->Fill(DPhi, DEta);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix3_4_EMCJet[2]->Fill(DPhi);
                        fHistDPhiMix28_EMCJet->Fill(DPhi);
                        fHistDPhiDEtaMix28_EMCJet->Fill(DPhi, DEta);
                    }
                }

                //5GeV
                if(aodtrackassoc->Pt()>4){
                    if(MBtrg){
                        fHistDPhiMix4_MB[2]->Fill(DPhi);
                        fHistDPhiMix28_MB->Fill(DPhi);
                        fHistDPhiDEtaMix28_MB->Fill(DPhi, DEta);
                    }
                    if(EMC7trg){
                        fHistDPhiMix4_EMC7[2]->Fill(DPhi);
                        fHistDPhiMix28_EMC7->Fill(DPhi);
                        fHistDPhiDEtaMix28_EMC7->Fill(DPhi, DEta);
                    }
                    if(EMCJettrg){
                        fHistDPhiMix4_EMCJet[2]->Fill(DPhi);
                        fHistDPhiMix28_EMCJet->Fill(DPhi);
                        fHistDPhiDEtaMix28_EMCJet->Fill(DPhi, DEta);
                    }
                }
            }
        }
    }
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMCJet);
    return;
}

//________________________________________________________________________
void AliAnalysisTaskPSHFE::Terminate(Option_t *) 
{
    
}
