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
    #include "TSystem.h"

    #include "AliAnalysisTaskSE.h"
    #include "AliAnalysisManager.h"
    #include "AliStack.h"
    #include "AliESDtrackCuts.h"
    #include "AliAODTrack.h"
    #include "AliAODCaloCluster.h"
    #include "AliAODEvent.h"
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
fOutputEMCEGA(0),
fOutputEMCJet(0),
fTrackCutsStrong(0),
fTrackCutsWeak(0),
fPoolMan(0),
fPool(0),
EMC7trg(0),
EMCEGAtrg(0),
EMCJettrg(0),
MBtrg(0),
tagStrong(0),
tagPhot(0),
UseNonSignalEvents(0),

fHistPIDRejection(0),
fHistNElecPerEvent(0),
fHistTestDCA(0),
fHistTestEMCEnergy(0),
fHistTestTPCdEdx(0),
fHistTestEOP(0),
fHistTestOGDPhi(0),
fHistTestPt(0),
fHistTestInvMassElecLike(0),
fHistTestInvMassElecUnLike(0),
fHistTestInvMassPionLike(0),
fHistTestInvMassPionUnLike(0),
fHistTestDPhiSpeNoSec(0),
fHistTestDPhi18Sec(0),
fHistTestDPhi18NoSec(0),
fHistTestDPhiType(0),

fHistTPCNClus_MB(0),
fHistITSNClus_MB(0),
fHistImpPar_MB(0),
fHistImpParTag_MB(0),
fHistNevents_MB(0),
fHistPtSum_MB(0),
fHistPtSumTag_MB(0),
fHistPtSumEMC_MB(0),
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
fHistDPhi18Spe_MB(0),

fHistTPCNClus_EMC7(0),
fHistITSNClus_EMC7(0),
fHistImpPar_EMC7(0),
fHistImpParTag_EMC7(0),
fHistNevents_EMC7(0),
fHistPtSum_EMC7(0),
fHistPtSumTag_EMC7(0),
fHistPtSumEMC_EMC7(0),
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
fHistDPhi18Spe_EMC7(0),

fHistTPCNClus_EMCEGA(0),
fHistITSNClus_EMCEGA(0),
fHistImpPar_EMCEGA(0),
fHistImpParTag_EMCEGA(0),
fHistNevents_EMCEGA(0),
fHistPtSum_EMCEGA(0),
fHistPtSumTag_EMCEGA(0),
fHistPtSumEMC_EMCEGA(0),
fHistEtaPhi_EMCEGA(0),
fHistEtaPhiTag_EMCEGA(0),
fHistDPhi28_EMCEGA(0),
fHistDPhiDEta28_EMCEGA(0),
fHistDPhiMix28_EMCEGA(0),
fHistDPhiDEtaMix28_EMCEGA(0),
fHistInvMassElecLike_EMCEGA(0),
fHistInvMassElecUnLike_EMCEGA(0),
fHistOpAngElecLike_EMCEGA(0),
fHistOpAngElecUnLike_EMCEGA(0),
fHistPtAssoc_EMCEGA(0),
fHistPtAssocMix_EMCEGA(0),
fHistPtTag_EMCEGA(0),
fHistPhotoMismatch_EMCEGA(0),
fHistDPhi18Spe_EMCEGA(0),

fHistTPCNClus_EMCJet(0),
fHistITSNClus_EMCJet(0),
fHistImpPar_EMCJet(0),
fHistImpParTag_EMCJet(0),
fHistNevents_EMCJet(0),
fHistPtSum_EMCJet(0),
fHistPtSumTag_EMCJet(0),
fHistPtSumEMC_EMCJet(0),
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
fHistDPhi18Spe_EMCJet(0)

    // The last in the above list should not have a comma after it
{

    //Init the DPhi Plots here because they are stored in arrays
    for(Int_t i=0; i<3; i++){
        fHistDPhi300_1_MB[i]=0;
        fHistDPhi1_2_MB[i]=0;
        fHistDPhi2_4_MB[i]=0;
        fHistDPhi4_8_MB[i]=0;


        fHistDPhi300_1_EMC7[i]=0;
        fHistDPhi1_2_EMC7[i]=0;
        fHistDPhi2_4_EMC7[i]=0;
        fHistDPhi4_8_EMC7[i]=0;


        fHistDPhi300_1_EMCEGA[i]=0;
        fHistDPhi1_2_EMCEGA[i]=0;
        fHistDPhi2_4_EMCEGA[i]=0;
        fHistDPhi4_8_EMCEGA[i]=0;

        fHistDPhi300_1_EMCJet[i]=0;
        fHistDPhi1_2_EMCJet[i]=0;
        fHistDPhi2_4_EMCJet[i]=0;
        fHistDPhi4_8_EMCJet[i]=0;

        //ME histos

        fHistDPhiMix300_1_MB[i]=0;
        fHistDPhiMix1_2_MB[i]=0;
        fHistDPhiMix2_4_MB[i]=0;
        fHistDPhiMix4_8_MB[i]=0;

        fHistDPhiMix300_1_EMC7[i]=0;
        fHistDPhiMix1_2_EMC7[i]=0;
        fHistDPhiMix2_4_EMC7[i]=0;
        fHistDPhiMix4_8_EMC7[i]=0;

        fHistDPhiMix300_1_EMCEGA[i]=0;
        fHistDPhiMix1_2_EMCEGA[i]=0;
        fHistDPhiMix2_4_EMCEGA[i]=0;
        fHistDPhiMix4_8_EMCEGA[i]=0;

        fHistDPhiMix300_1_EMCJet[i]=0;
        fHistDPhiMix1_2_EMCJet[i]=0;
        fHistDPhiMix2_4_EMCJet[i]=0;
        fHistDPhiMix4_8_EMCJet[i]=0;
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

        //EMCEGA Plots
        fHistTPC_EMCTRD_EMCEGA[i]=0;

        fHistEMC_TPCTRD_EMCEGA[i]=0;

        fHistTRD_TPCEMC_EMCEGA[i]=0;

        //EMCJet Plots
        fHistTPC_EMCTRD_EMCJet[i]=0;

        fHistEMC_TPCTRD_EMCJet[i]=0;

        fHistTRD_TPCEMC_EMCJet[i]=0;

    }

    // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskPSHFE::AliAnalysisTaskPSHFE(const char *name) // All data members should be initialised here
    :AliAnalysisTaskSE(name),
fOutputMB(0),
fOutputEMC7(0),
fOutputEMCEGA(0),
fOutputEMCJet(0),
fTrackCutsStrong(0),
fTrackCutsWeak(0),
fPoolMan(0),
fPool(0),
EMC7trg(0),
EMCEGAtrg(0),
EMCJettrg(0),
MBtrg(0),
tagStrong(0),
tagPhot(0),
UseNonSignalEvents(0),

fHistPIDRejection(0),
fHistNElecPerEvent(0),
fHistTestDCA(0),
fHistTestEMCEnergy(0),
fHistTestTPCdEdx(0),
fHistTestEOP(0),
fHistTestOGDPhi(0),
fHistTestPt(0),
fHistTestInvMassElecLike(0),
fHistTestInvMassElecUnLike(0),
fHistTestInvMassPionLike(0),
fHistTestInvMassPionUnLike(0),
fHistTestDPhiSpeNoSec(0),
fHistTestDPhi18Sec(0),
fHistTestDPhi18NoSec(0),
fHistTestDPhiType(0),

fHistTPCNClus_MB(0),
fHistITSNClus_MB(0),
fHistImpPar_MB(0),
fHistImpParTag_MB(0),
fHistNevents_MB(0),
fHistPtSum_MB(0),
fHistPtSumTag_MB(0),
fHistPtSumEMC_MB(0),
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
fHistDPhi18Spe_MB(0),

fHistTPCNClus_EMC7(0),
fHistITSNClus_EMC7(0),
fHistImpPar_EMC7(0),
fHistImpParTag_EMC7(0),
fHistNevents_EMC7(0),
fHistPtSum_EMC7(0),
fHistPtSumTag_EMC7(0),
fHistPtSumEMC_EMC7(0),
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
fHistDPhi18Spe_EMC7(0),

fHistTPCNClus_EMCEGA(0),
fHistITSNClus_EMCEGA(0),
fHistImpPar_EMCEGA(0),
fHistImpParTag_EMCEGA(0),
fHistNevents_EMCEGA(0),
fHistPtSum_EMCEGA(0),
fHistPtSumTag_EMCEGA(0),
fHistPtSumEMC_EMCEGA(0),
fHistEtaPhi_EMCEGA(0),
fHistEtaPhiTag_EMCEGA(0),
fHistDPhi28_EMCEGA(0),
fHistDPhiDEta28_EMCEGA(0),
fHistDPhiMix28_EMCEGA(0),
fHistDPhiDEtaMix28_EMCEGA(0),
fHistInvMassElecLike_EMCEGA(0),
fHistInvMassElecUnLike_EMCEGA(0),
fHistOpAngElecLike_EMCEGA(0),
fHistOpAngElecUnLike_EMCEGA(0),
fHistPtAssoc_EMCEGA(0),
fHistPtAssocMix_EMCEGA(0),
fHistPtTag_EMCEGA(0),
fHistPhotoMismatch_EMCEGA(0),
fHistDPhi18Spe_EMCEGA(0),

fHistTPCNClus_EMCJet(0),
fHistITSNClus_EMCJet(0),
fHistImpPar_EMCJet(0),
fHistImpParTag_EMCJet(0),
fHistNevents_EMCJet(0),
fHistPtSum_EMCJet(0),
fHistPtSumTag_EMCJet(0),
fHistPtSumEMC_EMCJet(0),
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
fHistDPhi18Spe_EMCJet(0)

    // The last in the above list should not have a comma after it
{
    // Constructor
    // Define input and output slots here (never in the dummy constructor)
    // Input slot #0 works with a TChain - it is connected to the default input container
    // Output slot #1 writes into a TH1 container

    for(Int_t i=0; i<3; i++){

        fHistDPhi300_1_MB[i]=0;
        fHistDPhi1_2_MB[i]=0;
        fHistDPhi2_4_MB[i]=0;
        fHistDPhi4_8_MB[i]=0;


        fHistDPhi300_1_EMC7[i]=0;
        fHistDPhi1_2_EMC7[i]=0;
        fHistDPhi2_4_EMC7[i]=0;
        fHistDPhi4_8_EMC7[i]=0;


        fHistDPhi300_1_EMCEGA[i]=0;
        fHistDPhi1_2_EMCEGA[i]=0;
        fHistDPhi2_4_EMCEGA[i]=0;
        fHistDPhi4_8_EMCEGA[i]=0;


        fHistDPhi300_1_EMCJet[i]=0;
        fHistDPhi1_2_EMCJet[i]=0;
        fHistDPhi2_4_EMCJet[i]=0;
        fHistDPhi4_8_EMCJet[i]=0;


        //ME histos

        fHistDPhiMix300_1_MB[i]=0;
        fHistDPhiMix1_2_MB[i]=0;
        fHistDPhiMix2_4_MB[i]=0;
        fHistDPhiMix4_8_MB[i]=0;


        fHistDPhiMix300_1_EMC7[i]=0;
        fHistDPhiMix1_2_EMC7[i]=0;
        fHistDPhiMix2_4_EMC7[i]=0;
        fHistDPhiMix4_8_EMC7[i]=0;


        fHistDPhiMix300_1_EMCEGA[i]=0;
        fHistDPhiMix1_2_EMCEGA[i]=0;
        fHistDPhiMix2_4_EMCEGA[i]=0;
        fHistDPhiMix4_8_EMCEGA[i]=0;


        fHistDPhiMix300_1_EMCJet[i]=0;
        fHistDPhiMix1_2_EMCJet[i]=0;
        fHistDPhiMix2_4_EMCJet[i]=0;
        fHistDPhiMix4_8_EMCJet[i]=0;

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

        //EMCEGA Plots
        fHistTPC_EMCTRD_EMCEGA[i]=0;

        fHistEMC_TPCTRD_EMCEGA[i]=0;

        fHistTRD_TPCEMC_EMCEGA[i]=0;

        //EMCJet Plots
        fHistTPC_EMCTRD_EMCJet[i]=0;

        fHistEMC_TPCTRD_EMCJet[i]=0;

        fHistTRD_TPCEMC_EMCJet[i]=0;

    }


    DefineOutput(1, TList::Class());//MB
    DefineOutput(2, TList::Class());//EMC7
    DefineOutput(3, TList::Class());//EMCEGA
    DefineOutput(4, TList::Class());//EMCJet
}

//________________________________________________________________________
AliAnalysisTaskPSHFE::~AliAnalysisTaskPSHFE()
{
    // Destructor. Clean-up the output list, but not the histograms that are put inside
    // (the list is owner and will clean-up these histograms). Protect in PROOF case.
    for(Int_t i=0;i<3;i++){

        delete fHistDPhi300_1_MB[i];
        delete fHistDPhi1_2_MB[i];
        delete fHistDPhi2_4_MB[i];
        delete fHistDPhi4_8_MB[i];

        delete fHistDPhi300_1_EMC7[i];
        delete fHistDPhi1_2_EMC7[i];
        delete fHistDPhi2_4_EMC7[i];
        delete fHistDPhi4_8_EMC7[i];


        delete fHistDPhi300_1_EMCEGA[i];
        delete fHistDPhi1_2_EMCEGA[i];
        delete fHistDPhi2_4_EMCEGA[i];
        delete fHistDPhi4_8_EMCEGA[i];


        delete fHistDPhi300_1_EMCJet[i];
        delete fHistDPhi1_2_EMCJet[i];
        delete fHistDPhi2_4_EMCJet[i];
        delete fHistDPhi4_8_EMCJet[i];


        //ME histos

        delete fHistDPhiMix300_1_MB[i];
        delete fHistDPhiMix1_2_MB[i];
        delete fHistDPhiMix2_4_MB[i];
        delete fHistDPhiMix4_8_MB[i];


        delete fHistDPhiMix300_1_EMC7[i];
        delete fHistDPhiMix1_2_EMC7[i];
        delete fHistDPhiMix2_4_EMC7[i];
        delete fHistDPhiMix4_8_EMC7[i];


        delete fHistDPhiMix300_1_EMCEGA[i];
        delete fHistDPhiMix1_2_EMCEGA[i];
        delete fHistDPhiMix2_4_EMCEGA[i];
        delete fHistDPhiMix4_8_EMCEGA[i];


        delete fHistDPhiMix300_1_EMCJet[i];
        delete fHistDPhiMix1_2_EMCJet[i];
        delete fHistDPhiMix2_4_EMCJet[i];
        delete fHistDPhiMix4_8_EMCJet[i];

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

        //EMCEGA Plots
        delete fHistTPC_EMCTRD_EMCEGA[i];

        delete fHistEMC_TPCTRD_EMCEGA[i];

        delete fHistTRD_TPCEMC_EMCEGA[i];

        //EMCJet Plots
        delete fHistTPC_EMCTRD_EMCJet[i];

        delete fHistEMC_TPCTRD_EMCJet[i];

        delete fHistTRD_TPCEMC_EMCJet[i];

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
    fOutputEMCEGA = new TList();
    OpenFile(3);
    fOutputEMCEGA->SetOwner();  // IMPORTANT!
    fOutputEMCJet = new TList();
    OpenFile(4);
    fOutputEMCJet->SetOwner();  // IMPORTANT!

    //Initialize event pool stuff
    Double_t vertexBins[5] = { -10, -4,  0, 4, 10 };
    Int_t nZvtxBins  = 4;
    Double_t multBins[4] = {0, 100, 300, 500};
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

    fHistPhotoMismatch_EMCEGA = new TH1F("fHistPhotoMismatch_EMCEGA", "Electrons identified as 'heavy flavour' that fall in photonic invariant mass and opening angle cuts", 2, 0, 1);
    fHistPhotoMismatch_EMCEGA->GetXaxis()->SetTitle("Electrons");
    fHistPhotoMismatch_EMCEGA->GetYaxis()->SetTitle("Cts");

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

    fHistInvMassElecLike_EMCEGA = new TH1F("fHistInvMassElecLike_EMCEGA", "Invariant mass for all like-signed electron pairs", 100, 0, .5);
    fHistInvMassElecLike_EMCEGA->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecLike_EMCEGA->GetYaxis()->SetTitle("Cts");

    fHistInvMassElecLike_EMCJet = new TH1F("fHistInvMassElecLike_EMCJet", "Invariant mass for all like-signed electron pairs", 100, 0, .5);
    fHistInvMassElecLike_EMCJet->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecLike_EMCJet->GetYaxis()->SetTitle("Cts");

    fHistInvMassElecUnLike_MB = new TH1F("fHistInvMassElecUnLike_MB", "Invariant mass for all unlike-signed electron pairs", 100, 0, .5);
    fHistInvMassElecUnLike_MB->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecUnLike_MB->GetYaxis()->SetTitle("Cts");

    fHistInvMassElecUnLike_EMC7 = new TH1F("fHistInvMassElecUnLike_EMC7", "Invariant mass for all unlike-signed electron pairs", 100, 0, .5);
    fHistInvMassElecUnLike_EMC7->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecUnLike_EMC7->GetYaxis()->SetTitle("Cts");

    fHistInvMassElecUnLike_EMCEGA = new TH1F("fHistInvMassElecUnLike_EMCEGA", "Invariant mass for all unlike-signed electron pairs", 100, 0, .5);
    fHistInvMassElecUnLike_EMCEGA->GetXaxis()->SetTitle("Invariant Mass(Gev/c^2)");
    fHistInvMassElecUnLike_EMCEGA->GetYaxis()->SetTitle("Cts");

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

    fHistOpAngElecLike_EMCEGA = new TH1F("fHistOpAngElecLike_EMCEGA", "Opening angle for all like-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecLike_EMCEGA->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecLike_EMCEGA->GetYaxis()->SetTitle("Cts");

    fHistOpAngElecLike_EMCJet = new TH1F("fHistOpAngElecLike_EMCJet", "Opening angle for all like-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecLike_EMCJet->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecLike_EMCJet->GetYaxis()->SetTitle("Cts");

    fHistOpAngElecUnLike_MB = new TH1F("fHistOpAngElecUnLike_MB", "Opening angle for all unlike-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecUnLike_MB->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecUnLike_MB->GetYaxis()->SetTitle("Cts");

    fHistOpAngElecUnLike_EMC7 = new TH1F("fHistOpAngElecUnLike_EMC7", "Opening angle for all unlike-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecUnLike_EMC7->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecUnLike_EMC7->GetYaxis()->SetTitle("Cts");

    fHistOpAngElecUnLike_EMCEGA = new TH1F("fHistOpAngElecUnLike_EMCEGA", "Opening angle for all unlike-signed electron pairs", 100, 0, TMath::Pi());
    fHistOpAngElecUnLike_EMCEGA->GetXaxis()->SetTitle("Opening Angle(rad)");
    fHistOpAngElecUnLike_EMCEGA->GetYaxis()->SetTitle("Cts");

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

        //EMCEGA
        fHistTPC_EMCTRD_EMCEGA[i] = new TH2F(TString::Format("fHistTPC_EMCTRD_EMCEGA_%s",ptRangesPID[i].Data()), TString::Format("TPC nSigma for tracks with Pt between %s after EMC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, -10, 10);
        fHistTPC_EMCTRD_EMCEGA[i]->GetXaxis()->SetTitle("Pt");
        fHistTPC_EMCTRD_EMCEGA[i]->GetYaxis()->SetTitle("nSigma");

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

        //EMCEGA
        fHistEMC_TPCTRD_EMCEGA[i] = new TH1F(TString::Format("fHistEMC_TPCTRD_EMCEGA_%s",ptRangesPID[i].Data()), TString::Format("E/p for tracks with Pt between %s after TPC and TRD cuts",ptRangesPID[i].Data()), 100, 0, 1.5);
        fHistEMC_TPCTRD_EMCEGA[i]->GetXaxis()->SetTitle("E/p");
        fHistEMC_TPCTRD_EMCEGA[i]->GetYaxis()->SetTitle("Cts");

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

        //EMCEGA
        fHistTRD_TPCEMC_EMCEGA[i] = new TH2F(TString::Format("fHistTRD_TPCEMC_EMCEGA_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after TPC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCEMC_EMCEGA[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCEMC_EMCEGA[i]->GetYaxis()->SetTitle("electron Likelihood");

        //EMCJet
        fHistTRD_TPCEMC_EMCJet[i] = new TH2F(TString::Format("fHistTRD_TPCEMC_EMCJet_%s",ptRangesPID[i].Data()), TString::Format("TRD electron Likelihood for tracks with Pt between %s after EMC and EMC cuts",ptRangesPID[i].Data()), 100, 0, 10, 800, 0, 1);
        fHistTRD_TPCEMC_EMCJet[i]->GetXaxis()->SetTitle("Pt");
        fHistTRD_TPCEMC_EMCJet[i]->GetYaxis()->SetTitle("electron Likelihood");
    }

    //DPhi for candidate electrons 2-8 gev and assoc. particles >3gev
    fHistDPhi28_MB = new TH1F("fHistDPhi28_MB", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhi28_MB->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi28_MB->GetYaxis()->SetTitle("Cts");

    fHistDPhi28_EMC7 = new TH1F("fHistDPhi28_EMC7", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhi28_EMC7->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi28_EMC7->GetYaxis()->SetTitle("Cts");

    fHistDPhi28_EMCEGA = new TH1F("fHistDPhi28_EMCEGA", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhi28_EMCEGA->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi28_EMCEGA->GetYaxis()->SetTitle("Cts");

    fHistDPhi28_EMCJet = new TH1F("fHistDPhi28_EMCJet", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhi28_EMCJet->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi28_EMCJet->GetYaxis()->SetTitle("Cts");

    //DPhi by Eta for triggered particles 2-8 gev and assoc. particles >3gev
    fHistDPhiDEta28_MB = new TH2F("fHistDPhiDEta28_MB", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 50, -TMath::Pi()/2, 3*TMath::Pi()/2, 50, -0.9, 0.9);
    fHistDPhiDEta28_MB->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEta28_MB->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEta28_MB->GetZaxis()->SetTitle("Cts");

    fHistDPhiDEta28_EMC7 = new TH2F("fHistDPhiDEta28_EMC7", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 550, -TMath::Pi()/2, 3*TMath::Pi()/2, 50, -0.9, 0.9);
    fHistDPhiDEta28_EMC7->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEta28_EMC7->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEta28_EMC7->GetZaxis()->SetTitle("Cts");

    fHistDPhiDEta28_EMCEGA = new TH2F("fHistDPhiDEta28_EMCEGA", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 50, -TMath::Pi()/2, 3*TMath::Pi()/2, 50, -0.9, 0.9);
    fHistDPhiDEta28_EMCEGA->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEta28_EMCEGA->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEta28_EMCEGA->GetZaxis()->SetTitle("Cts");

    fHistDPhiDEta28_EMCJet = new TH2F("fHistDPhiDEta28_EMCJet", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev", 50, -TMath::Pi()/2, 3*TMath::Pi()/2, 50, -0.9, 0.9);
    fHistDPhiDEta28_EMCJet->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEta28_EMCJet->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEta28_EMCJet->GetZaxis()->SetTitle("Cts");

    //DPhi for candidate electrons 2-8 gev and assoc. particles >2Gev for Mixed Events
    fHistDPhiMix28_MB = new TH1F("fHistDPhiMix28_MB", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev for Mixed Events", 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhiMix28_MB->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiMix28_MB->GetYaxis()->SetTitle("Cts");

    fHistDPhiMix28_EMC7 = new TH1F("fHistDPhiMix28_EMC7", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev for Mixed Events", 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhiMix28_EMC7->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiMix28_EMC7->GetYaxis()->SetTitle("Cts");

    fHistDPhiMix28_EMCEGA = new TH1F("fHistDPhiMix28_EMCEGA", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev for Mixed Events", 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhiMix28_EMCEGA->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiMix28_EMCEGA->GetYaxis()->SetTitle("Cts");

    fHistDPhiMix28_EMCJet = new TH1F("fHistDPhiMix28_EMCJet", "Delta-Phi for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev for Mixed Events", 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistDPhiMix28_EMCJet->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiMix28_EMCJet->GetYaxis()->SetTitle("Cts");

    //DPhi by Eta for triggered particles 2-8 gev and assoc. particles >3gev
    fHistDPhiDEtaMix28_MB = new TH2F("fHistDPhiDEtaMix28_MB", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev for Mixed Events", 50, -TMath::Pi()/2, 3*TMath::Pi()/2, 50, -0.9, 0.9);
    fHistDPhiDEtaMix28_MB->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEtaMix28_MB->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEtaMix28_MB->GetZaxis()->SetTitle("Cts");

    fHistDPhiDEtaMix28_EMC7 = new TH2F("fHistDPhiDEtaMix28_EMC7", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev for Mixed Events", 50, -TMath::Pi()/2, 3*TMath::Pi()/2, 50, -0.9, 0.9);
    fHistDPhiDEtaMix28_EMC7->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEtaMix28_EMC7->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEtaMix28_EMC7->GetZaxis()->SetTitle("Cts");

    fHistDPhiDEtaMix28_EMCEGA = new TH2F("fHistDPhiDEtaMix28_EMCEGA", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev for Mixed Events", 50, -TMath::Pi()/2, 3*TMath::Pi()/2, 50, -0.9, 0.9);
    fHistDPhiDEtaMix28_EMCEGA->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEtaMix28_EMCEGA->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEtaMix28_EMCEGA->GetZaxis()->SetTitle("Cts");

    fHistDPhiDEtaMix28_EMCJet = new TH2F("fHistDPhiDEtaMix28_EMCJet", "Delta-Phi by Delta-Eta for candidate electrons with 2<pt<8Gev and assoc. with pt>2Gev for Mixed Events", 50, -TMath::Pi()/2, 3*TMath::Pi()/2, 50, -0.9, 0.9);
    fHistDPhiDEtaMix28_EMCJet->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhiDEtaMix28_EMCJet->GetYaxis()->SetTitle("Delta-Eta");
    fHistDPhiDEtaMix28_EMCJet->GetZaxis()->SetTitle("Cts");

    //DPhi by dEdx for triggered particles 2-8 gev and assoc. particles >2gev
    fHistDPhi18Spe_MB = new TH2F("fHistDPhi18Spe_MB", "Delta-Phi by most probable species for candidate electrons with 1<pt<8Gev and assoc. with pt>.3Gev", 50, -TMath::Pi()/2, 3*TMath::Pi()/2, 10, 0, 10);
    fHistDPhi18Spe_MB->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi18Spe_MB->GetYaxis()->SetTitle("Species");
    fHistDPhi18Spe_MB->GetYaxis()->SetBinLabel(1, "Unkown");
    fHistDPhi18Spe_MB->GetYaxis()->SetBinLabel(2, "Electron");
    fHistDPhi18Spe_MB->GetYaxis()->SetBinLabel(3, "Muon");
    fHistDPhi18Spe_MB->GetYaxis()->SetBinLabel(4, "Pion");
    fHistDPhi18Spe_MB->GetYaxis()->SetBinLabel(5, "Kaon");
    fHistDPhi18Spe_MB->GetYaxis()->SetBinLabel(6, "Proton");
    fHistDPhi18Spe_MB->GetYaxis()->SetBinLabel(7, "Deuteron");
    fHistDPhi18Spe_MB->GetZaxis()->SetTitle("Cts");

    fHistDPhi18Spe_EMC7 = new TH2F("fHistDPhi18Spe_EMC7", "Delta-Phi by most probable species for candidate electrons with 1<pt<8Gev and assoc. with pt>.3Gev", 50, -TMath::Pi()/2, 3*TMath::Pi()/2, 10, 0, 10);
    fHistDPhi18Spe_EMC7->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi18Spe_EMC7->GetYaxis()->SetTitle("Species");
    fHistDPhi18Spe_EMC7->GetYaxis()->SetBinLabel(1, "Unkown");
    fHistDPhi18Spe_EMC7->GetYaxis()->SetBinLabel(2, "Electron");
    fHistDPhi18Spe_EMC7->GetYaxis()->SetBinLabel(3, "Muon");
    fHistDPhi18Spe_EMC7->GetYaxis()->SetBinLabel(4, "Pion");
    fHistDPhi18Spe_EMC7->GetYaxis()->SetBinLabel(5, "Kaon");
    fHistDPhi18Spe_EMC7->GetYaxis()->SetBinLabel(6, "Proton");
    fHistDPhi18Spe_EMC7->GetYaxis()->SetBinLabel(7, "Deuteron");
    fHistDPhi18Spe_EMC7->GetZaxis()->SetTitle("Cts");

    fHistDPhi18Spe_EMCEGA = new TH2F("fHistDPhi18Spe_EMCEGA", "Delta-Phi by most probable species for candidate electrons with 1<pt<8Gev and assoc. with pt>.3Gev", 50, -TMath::Pi()/2, 3*TMath::Pi()/2, 10, 0, 10);
    fHistDPhi18Spe_EMCEGA->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi18Spe_EMCEGA->GetYaxis()->SetTitle("Species");
    fHistDPhi18Spe_EMCEGA->GetYaxis()->SetBinLabel(1, "Unkown");
    fHistDPhi18Spe_EMCEGA->GetYaxis()->SetBinLabel(2, "Electron");
    fHistDPhi18Spe_EMCEGA->GetYaxis()->SetBinLabel(3, "Muon");
    fHistDPhi18Spe_EMCEGA->GetYaxis()->SetBinLabel(4, "Pion");
    fHistDPhi18Spe_EMCEGA->GetYaxis()->SetBinLabel(5, "Kaon");
    fHistDPhi18Spe_EMCEGA->GetYaxis()->SetBinLabel(6, "Proton");
    fHistDPhi18Spe_EMCEGA->GetYaxis()->SetBinLabel(7, "Deuteron");
    fHistDPhi18Spe_EMCEGA->GetZaxis()->SetTitle("Cts");

    fHistDPhi18Spe_EMCJet = new TH2F("fHistDPhi18Spe_EMCJet", "Delta-Phi by most probable species for candidate electrons with 1<pt<8Gev and assoc. with pt>.3Gev", 50, -TMath::Pi()/2, 3*TMath::Pi()/2, 10, 0, 10);
    fHistDPhi18Spe_EMCJet->GetXaxis()->SetTitle("Delta-Phi");
    fHistDPhi18Spe_EMCJet->GetYaxis()->SetTitle("Species");
    fHistDPhi18Spe_EMCJet->GetYaxis()->SetBinLabel(1, "Unkown");
    fHistDPhi18Spe_EMCJet->GetYaxis()->SetBinLabel(2, "Electron");
    fHistDPhi18Spe_EMCJet->GetYaxis()->SetBinLabel(3, "Muon");
    fHistDPhi18Spe_EMCJet->GetYaxis()->SetBinLabel(4, "Pion");
    fHistDPhi18Spe_EMCJet->GetYaxis()->SetBinLabel(5, "Kaon");
    fHistDPhi18Spe_EMCJet->GetYaxis()->SetBinLabel(6, "Proton");
    fHistDPhi18Spe_EMCJet->GetYaxis()->SetBinLabel(7, "Deuteron");
    fHistDPhi18Spe_EMCJet->GetZaxis()->SetTitle("Cts");

    // Delta Phi for tracks > 300MeV
    for(Int_t i=0; i<3; i++){
        fHistDPhi300_1_MB[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.3-1Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<1Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi300_1_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi300_1_MB[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhi300_1_EMC7[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.3-1Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<1Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi300_1_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi300_1_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhi300_1_EMCEGA[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.3-1Gev_EMCEGA",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<1Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi300_1_EMCEGA[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi300_1_EMCEGA[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhi300_1_EMCJet[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_.3-1Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<1Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi300_1_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi300_1_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }

    // Delta Phi for tracks > 1GeV
    for(Int_t i=0; i<3; i++){
        fHistDPhi1_2_MB[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_1-2Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi1_2_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi1_2_MB[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhi1_2_EMC7[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_1-2Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev",ptRangesDPhi[i].Data()).Data(),50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi1_2_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi1_2_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhi1_2_EMCEGA[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_1-2Gev_EMCEGA",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi1_2_EMCEGA[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi1_2_EMCEGA[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhi1_2_EMCJet[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_1-2Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi1_2_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi1_2_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }

    // Delta Phi for tracks > 2GeV
    for(Int_t i=0; i<3; i++){
        fHistDPhi2_4_MB[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_2-4Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<4Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi2_4_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi2_4_MB[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhi2_4_EMC7[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_2-4Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<4Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi2_4_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi2_4_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhi2_4_EMCEGA[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_2-4Gev_EMCEGA",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<4Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi2_4_EMCEGA[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi2_4_EMCEGA[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhi2_4_EMCJet[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_2-4Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<4Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi2_4_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi2_4_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }

    // Delta Phi for tracks > 3GeV
    for(Int_t i=0; i<3; i++){
        fHistDPhi4_8_MB[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_4-8Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 4Gev<pt<8Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi4_8_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi4_8_MB[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhi4_8_EMC7[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_4-8Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 4Gev<pt<8Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi4_8_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi4_8_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhi4_8_EMCEGA[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_4-8Gev_EMCEGA",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 4Gev<pt<8Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi4_8_EMCEGA[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi4_8_EMCEGA[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhi4_8_EMCJet[i] = new TH1F(TString::Format("fHistDPhi_trig_%s_assoc_4-8Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 4Gev<pt<8Gev",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhi4_8_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhi4_8_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }

    

    //Mixed Event DPhi plots

    // Delta Phi for tracks > 300MeV

    for(Int_t i=0; i<3; i++){
        fHistDPhiMix300_1_MB[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_.3-1Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<1Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix300_1_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix300_1_MB[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhiMix300_1_EMC7[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_.3-1Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<1Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix300_1_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix300_1_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhiMix300_1_EMCEGA[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_.3-1Gev_EMCEGA",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<1Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix300_1_EMCEGA[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix300_1_EMCEGA[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhiMix300_1_EMCJet[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_.3-1Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated .3Gev<pt<1Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix300_1_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix300_1_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }

    // Delta Phi for tracks > 1GeV
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix1_2_MB[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_1-2Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix1_2_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix1_2_MB[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhiMix1_2_EMC7[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_1-2Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix1_2_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix1_2_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhiMix1_2_EMCEGA[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_1-2Gev_EMCEGA",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix1_2_EMCEGA[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix1_2_EMCEGA[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhiMix1_2_EMCJet[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_1-2Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 1Gev<pt<2Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix1_2_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix1_2_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }

    // Delta Phi for tracks > 2GeV
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix2_4_MB[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_2-4Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<4Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix2_4_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix2_4_MB[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhiMix2_4_EMC7[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_2-4Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<4Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix2_4_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix2_4_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhiMix2_4_EMCEGA[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_2-4Gev_EMCEGA",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<4Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix2_4_EMCEGA[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix2_4_EMCEGA[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhiMix2_4_EMCJet[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_2-4Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 2Gev<pt<4Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix2_4_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix2_4_EMCJet[i]->GetYaxis()->SetTitle("Cts");
    }

    // Delta Phi for tracks > 3GeV
    for(Int_t i=0; i<3; i++){
        fHistDPhiMix4_8_MB[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_4-8Gev_MB",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 4Gev<pt<8Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix4_8_MB[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix4_8_MB[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhiMix4_8_EMC7[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_4-8Gev_EMC7",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 4Gev<pt<8Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix4_8_EMC7[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix4_8_EMC7[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhiMix4_8_EMCEGA[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_4-8Gev_EMCEGA",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 4Gev<pt<8Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix4_8_EMCEGA[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix4_8_EMCEGA[i]->GetYaxis()->SetTitle("Cts");
    }

    for(Int_t i=0; i<3; i++){
        fHistDPhiMix4_8_EMCJet[i] = new TH1F(TString::Format("fHistDPhiMix_trig_%s_assoc_4-8Gev_EMCJet",ptRangesDPhi[i].Data()).Data(), TString::Format("Delta-Phi for candidate electrons w/ pt=%s and Associated 4Gev<pt<8Gev for Mixed Events",ptRangesDPhi[i].Data()).Data(), 50, -TMath::Pi()/2, 3*TMath::Pi()/2);
        fHistDPhiMix4_8_EMCJet[i]->GetXaxis()->SetTitle("Delta-Phi");
        fHistDPhiMix4_8_EMCJet[i]->GetYaxis()->SetTitle("Cts");
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

    fHistEtaPhiTag_EMCEGA = new TH2F("fHistEtaPhiTag_EMCEGA", "Eta-Phi distribution of tracks in tagged events", 100, -.9,.9,100,0,2*TMath::Pi());
    fHistEtaPhiTag_EMCEGA->GetXaxis()->SetTitle("Eta");
    fHistEtaPhiTag_EMCEGA->GetYaxis()->SetTitle("Phi");

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

    fHistEtaPhi_EMCEGA = new TH2F("fHistEtaPhi_EMCEGA", "Eta-Phi distribution of tracks", 100, -.9,.9,100,0,2*TMath::Pi());
    fHistEtaPhi_EMCEGA->GetXaxis()->SetTitle("Eta");
    fHistEtaPhi_EMCEGA->GetYaxis()->SetTitle("Phi");

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

    fHistPtSum_EMCEGA = new TH1F("fHistPtSum_EMCEGA", "Pt sum for events w/o an electron candidate", 500, 0, 500);
    fHistPtSum_EMCEGA->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSum_EMCEGA->GetYaxis()->SetTitle("Cts");

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

    fHistPtSumTag_EMCEGA = new TH1F("fHistPtSumTag_EMCEGA", "Pt sum for events w/ an electron candidate", 500, 0, 500);
    fHistPtSumTag_EMCEGA->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSumTag_EMCEGA->GetYaxis()->SetTitle("Cts");

    fHistPtSumTag_EMCJet = new TH1F("fHistPtSumTag_EMCJet", "Pt sum for events w/ an electron candidate", 500, 0, 500);
    fHistPtSumTag_EMCJet->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSumTag_EMCJet->GetYaxis()->SetTitle("Cts");

    // Energy per event in EMC acceptance histos
    fHistPtSumEMC_MB = new TH1F("fHistPtSumEMC_MB", "Pt sum for events in EMCal acceptance", 500, 0, 500);
    fHistPtSumEMC_MB->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSumEMC_MB->GetYaxis()->SetTitle("Cts");

    fHistPtSumEMC_EMC7 = new TH1F("fHistPtSumEMC_EMC7", "Pt sum for events in EMCal acceptance", 500, 0, 500);
    fHistPtSumEMC_EMC7->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSumEMC_EMC7->GetYaxis()->SetTitle("Cts");

    fHistPtSumEMC_EMCEGA = new TH1F("fHistPtSumEMC_EMCEGA", "Pt sum for events in EMCal acceptance", 500, 0, 500);
    fHistPtSumEMC_EMCEGA->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSumEMC_EMCEGA->GetYaxis()->SetTitle("Cts");

    fHistPtSumEMC_EMCJet = new TH1F("fHistPtSumEMC_EMCJet", "Pt sum for events in EMCal acceptance", 500, 0, 500);
    fHistPtSumEMC_EMCJet->GetXaxis()->SetTitle("Pt Sum");
    fHistPtSumEMC_EMCJet->GetYaxis()->SetTitle("Cts");

    // Numbers of events
    fHistNevents_MB = new TH1F("fHistNevents_MB", "Number of events that have an 'electron'", 2,0,1);
    fHistNevents_MB->GetXaxis()->SetBinLabel(1,"Events");
    fHistNevents_MB->GetXaxis()->SetBinLabel(2,"Events containing candidates");
    fHistNevents_MB->GetYaxis()->SetTitle("Cts");

    fHistNevents_EMC7 = new TH1F("fHistNevents_EMC7", "Number of events that have an 'electron'", 2,0,1);
    fHistNevents_EMC7->GetXaxis()->SetBinLabel(1,"Events");
    fHistNevents_EMC7->GetXaxis()->SetBinLabel(2,"Events containing candidates");
    fHistNevents_EMC7->GetYaxis()->SetTitle("Cts");

    fHistNevents_EMCEGA = new TH1F("fHistNevents_EMCEGA", "Number of events that have an 'electron'", 2,0,1);
    fHistNevents_EMCEGA->GetXaxis()->SetBinLabel(1,"Events");
    fHistNevents_EMCEGA->GetXaxis()->SetBinLabel(2,"Events containing candidates");
    fHistNevents_EMCEGA->GetYaxis()->SetTitle("Cts");

    fHistNevents_EMCJet = new TH1F("fHistNevents_EMCJet", "Number of events that have an 'electron'", 2,0,1);
    fHistNevents_EMCJet->GetXaxis()->SetBinLabel(1,"Events");
    fHistNevents_EMCJet->GetXaxis()->SetBinLabel(2,"Events containing candidates");
    fHistNevents_EMCJet->GetYaxis()->SetTitle("Cts");

    //Impact Parameter histos
    fHistImpPar_MB = new TH1F("fHistImpPar_MB", "Impact Parameter distribution in xy plane for all tracks", 100,-1, 1);
    fHistImpPar_MB->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpPar_MB->GetYaxis()->SetTitle("Count");

    fHistImpPar_EMC7 = new TH1F("fHistImpPar_EMC7", "Impact Parameter distribution in xy plane for all tracks", 100,-1, 1);
    fHistImpPar_EMC7->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpPar_EMC7->GetYaxis()->SetTitle("Count");

    fHistImpPar_EMCEGA = new TH1F("fHistImpPar_EMCEGA", "Impact Parameter distribution in xy plane for all tracks", 100,-1, 1);
    fHistImpPar_EMCEGA->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpPar_EMCEGA->GetYaxis()->SetTitle("Count");

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

    fHistImpParTag_EMCEGA = new TH1F("fHistImpParTag_EMCEGA", "Impact Parameter distribution in xy plane for electron candidates", 100,-1, 1);
    fHistImpParTag_EMCEGA->GetXaxis()->SetTitle("Impact Parameter(cm)");
    fHistImpParTag_EMCEGA->GetYaxis()->SetTitle("Count");

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

    fHistTPCNClus_EMCEGA = new TH1F("fHistTPCNClus_EMCEGA", "Number of Clusters per track in TPC", 159, 0, 159);
    fHistTPCNClus_EMCEGA->GetXaxis()->SetTitle("Number of TPC Clusters");
    fHistTPCNClus_EMCEGA->GetYaxis()->SetTitle("Number of Tracks");

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

    fHistITSNClus_EMCEGA = new TH1F("fHistITSNClus_EMCEGA", "Number of Clusters per Track in ITS", 10, 0, 10);
    fHistITSNClus_EMCEGA->GetXaxis()->SetTitle("Number of ITS Clusters");
    fHistITSNClus_EMCEGA->GetYaxis()->SetTitle("Number of Tracks");

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

    fHistPtAssoc_EMCEGA = new TH1F("fHistPtAssoc_EMCEGA", "Pt distribution for associated tracks", 100,0, 15);
    fHistPtAssoc_EMCEGA->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtAssoc_EMCEGA->GetYaxis()->SetTitle("Count");

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

    fHistPtAssocMix_EMCEGA = new TH1F("fHistPtAssocMix_EMCEGA", "Pt distribution for associated tracks in mixed events", 100,0, 15);
    fHistPtAssocMix_EMCEGA->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtAssocMix_EMCEGA->GetYaxis()->SetTitle("Count");

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

    fHistPtTag_EMCEGA = new TH1F("fHistPtTag_EMCEGA", "Pt distribution for electron candidates", 100,0, 15);
    fHistPtTag_EMCEGA->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtTag_EMCEGA->GetYaxis()->SetTitle("Count");

    fHistPtTag_EMCJet = new TH1F("fHistPtTag_EMCJet", "Pt distribution for electron candidates", 100,0, 15);
    fHistPtTag_EMCJet->GetXaxis()->SetTitle("Pt(Gev)");
    fHistPtTag_EMCJet->GetYaxis()->SetTitle("Count");

    //test histos
    fHistTestDCA = new TH1F("fHistTestDCA", "DCA distribution for all tracks with DPhi to any candidate electron <0.1 rad", 100,-15, 15);
    fHistTestDCA->GetXaxis()->SetTitle("DCA(cm)");
    fHistTestDCA->GetYaxis()->SetTitle("Count");

    fHistTestEMCEnergy = new TH1F("fHistTestEMCEnergy", "Energy from EMCal for all tracks with DPhi to any candidate electron <0.1 rad", 100, 0, 10);
    fHistTestEMCEnergy->GetXaxis()->SetTitle("EMC Energy[GeV]");
    fHistTestEMCEnergy->GetYaxis()->SetTitle("Counts");

    fHistTestTPCdEdx = new TH2F("fHistTestTPCdEdx", "TPC dE/dx for all tracks with DPhi to any candidate electron <0.1 rad", 100, 0,8, 300, -30, 180);
    fHistTestTPCdEdx->GetYaxis()->SetTitle("TPC dE/dx[a.u.]");
    fHistTestTPCdEdx->GetXaxis()->SetTitle("pT[GeV/c]");

    fHistTestEOP = new TH1F("fHistTestEOP", "E/p for all tracks with DPhi to any candidate electron <0.1 rad", 30, 0, 1.5);
    fHistTestEOP->GetXaxis()->SetTitle("E/p[c]");
    fHistTestEOP->GetYaxis()->SetTitle("Counts");

    fHistTestOGDPhi = new TH1F("fHistTestOGDPhi", "Original DPhi before periodicity correction", 100, -2*TMath::Pi(), 2*TMath::Pi());
    fHistTestOGDPhi->GetXaxis()->SetTitle("DPhi[rad]");
    fHistTestOGDPhi->GetYaxis()->SetTitle("Counts");

    fHistTestPt = new TH1F("fHistTestPt", "Pt distribution for associated particles nearly on top of tagged particle", 30, 0, 8);
    fHistTestPt->GetXaxis()->SetTitle("Pt[Gev/c]");
    fHistTestPt->GetYaxis()->SetTitle("Cts");

    fHistTestInvMassElecLike = new TH1F("fHistTestInvMassElecLike", "Invariant Mass distribution for associated electrons of like sign in |DPhi|<0.1rad", 30, 0, 8);
    fHistTestInvMassElecLike->GetXaxis()->SetTitle("Mass[Gev/c^2]");
    fHistTestInvMassElecLike->GetYaxis()->SetTitle("Cts");

    fHistTestInvMassElecUnLike = new TH1F("fHistTestInvMassElecUnLike", "Invariant Mass distribution for associated electrons of unlike sign with |DPhi|<0.1rad", 30, 0, 8);
    fHistTestInvMassElecUnLike->GetXaxis()->SetTitle("Mass[Gev/c^2]");
    fHistTestInvMassElecUnLike->GetYaxis()->SetTitle("Cts");

    fHistTestInvMassPionLike = new TH1F("fHistTestInvMassPionLike", "Invariant Mass distribution for associated pions with |DPhi|<0.1rad", 30, 0, 8);
    fHistTestInvMassPionLike->GetXaxis()->SetTitle("Mass[Gev/c^2]");
    fHistTestInvMassPionLike->GetYaxis()->SetTitle("Cts");

    fHistTestInvMassPionUnLike = new TH1F("fHistTestInvMassPionUnLike", "Invariant Mass distribution for associated pions with |DPhi|<0.1rad", 30, 0, 8);
    fHistTestInvMassPionUnLike->GetXaxis()->SetTitle("Mass[Gev/c^2]");
    fHistTestInvMassPionUnLike->GetYaxis()->SetTitle("Cts");

    //DPhi by dEdx for triggered particles 1-8 gev and assoc. particles >.3gev
    fHistTestDPhiSpeNoSec = new TH2F("fHistTestDPhiSpeNoSec", "Delta-Phi by most probable species for candidate electrons with 1<pt<8Gev and assoc. with pt>.3Gev with no secondary tracks", 100, -TMath::Pi()/2, 3*TMath::Pi()/2, 10, 0, 10);
    fHistTestDPhiSpeNoSec->GetXaxis()->SetTitle("Delta-Phi");
    fHistTestDPhiSpeNoSec->GetYaxis()->SetTitle("Species");
    fHistTestDPhiSpeNoSec->GetYaxis()->SetBinLabel(1, "Unkown");
    fHistTestDPhiSpeNoSec->GetYaxis()->SetBinLabel(2, "Electron");
    fHistTestDPhiSpeNoSec->GetYaxis()->SetBinLabel(3, "Muon");
    fHistTestDPhiSpeNoSec->GetYaxis()->SetBinLabel(4, "Pion");
    fHistTestDPhiSpeNoSec->GetYaxis()->SetBinLabel(5, "Kaon");
    fHistTestDPhiSpeNoSec->GetYaxis()->SetBinLabel(6, "Proton");
    fHistTestDPhiSpeNoSec->GetYaxis()->SetBinLabel(7, "Deuteron");
    fHistTestDPhiSpeNoSec->GetZaxis()->SetTitle("Cts");

    //DPhi by track type to test for tpconly tracks for triggered particles 1-8 gev and assoc. particles >.3gev
    fHistTestDPhiType = new TH2F("fHistTestDPhiType", "Delta-Phi by track type for candidate electrons with 1<pt<8Gev and assoc. with pt>.3Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2, 4, 0, 4);
    fHistTestDPhiType->GetXaxis()->SetTitle("Delta-Phi");
    fHistTestDPhiType->GetYaxis()->SetTitle("Track Type");
    fHistTestDPhiType->GetYaxis()->SetBinLabel(1, "Hybrid");
    fHistTestDPhiType->GetYaxis()->SetBinLabel(2, "Complementary");
    fHistTestDPhiType->GetYaxis()->SetBinLabel(3, "TPC Only");
    fHistTestDPhiType->GetZaxis()->SetTitle("Cts");

    //DPhi for candidate electrons 1-8 gev and assoc. particles >.3gev with secondaries
    fHistTestDPhi18Sec = new TH1F("fHistTestDPhi18Sec", "Delta-Phi for candidate electrons with 1<pt<8Gev and assoc. with pt>.3Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistTestDPhi18Sec->GetXaxis()->SetTitle("Delta-Phi");
    fHistTestDPhi18Sec->GetYaxis()->SetTitle("Cts");

    //DPhi for candidate electrons 1-8 gev and assoc. particles >.3gev without secondaries
    fHistTestDPhi18NoSec = new TH1F("fHistTestDPhi18NoSec", "Delta-Phi for candidate electrons with 1<pt<8Gev and assoc. with pt>.3Gev", 100, -TMath::Pi()/2, 3*TMath::Pi()/2);
    fHistTestDPhi18NoSec->GetXaxis()->SetTitle("Delta-Phi");
    fHistTestDPhi18NoSec->GetYaxis()->SetTitle("Cts");

    //Add rejection plots to MB plots since it is the easiest place
    fOutputMB->Add(fHistPIDRejection);
    fOutputMB->Add(fHistNElecPerEvent);
    //samesies for the test plots
    fOutputMB->Add(fHistTestDCA);
    fOutputMB->Add(fHistTestEMCEnergy);
    fOutputMB->Add(fHistTestTPCdEdx);
    fOutputMB->Add(fHistTestEOP);
    fOutputMB->Add(fHistTestOGDPhi);
    fOutputMB->Add(fHistTestPt);
    fOutputMB->Add(fHistTestInvMassElecLike);
    fOutputMB->Add(fHistTestInvMassElecUnLike);
    fOutputMB->Add(fHistTestInvMassPionLike);
    fOutputMB->Add(fHistTestInvMassPionUnLike);
    fOutputMB->Add(fHistTestDPhiSpeNoSec);
    fOutputMB->Add(fHistTestDPhiType);
    fOutputMB->Add(fHistTestDPhi18Sec);
    fOutputMB->Add(fHistTestDPhi18NoSec);

    fOutputMB->Add(fHistPhotoMismatch_MB);
    fOutputMB->Add(fHistPtAssoc_MB);
    fOutputMB->Add(fHistPtAssocMix_MB);
    fOutputMB->Add(fHistPtTag_MB);
    fOutputMB->Add(fHistTPCNClus_MB);
    fOutputMB->Add(fHistITSNClus_MB);
    fOutputMB->Add(fHistImpPar_MB);
    fOutputMB->Add(fHistImpParTag_MB);
    fOutputMB->Add(fHistNevents_MB);
    fOutputMB->Add(fHistPtSum_MB);
    fOutputMB->Add(fHistPtSumTag_MB);
    fOutputMB->Add(fHistPtSumEMC_MB);
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

    }
    for(Int_t i=0; i<3;i++){
        fOutputMB->Add(fHistDPhi300_1_MB[i]);
        fOutputMB->Add(fHistDPhi1_2_MB[i]);
        fOutputMB->Add(fHistDPhi2_4_MB[i]);
        fOutputMB->Add(fHistDPhi4_8_MB[i]);
        fOutputMB->Add(fHistDPhiMix300_1_MB[i]);
        fOutputMB->Add(fHistDPhiMix1_2_MB[i]);
        fOutputMB->Add(fHistDPhiMix2_4_MB[i]);
        fOutputMB->Add(fHistDPhiMix4_8_MB[i]);
    }
    fOutputMB->Add(fHistDPhi28_MB);
    fOutputMB->Add(fHistDPhiDEta28_MB);
    fOutputMB->Add(fHistDPhiMix28_MB);
    fOutputMB->Add(fHistDPhiDEtaMix28_MB);
    fOutputMB->Add(fHistDPhi18Spe_MB);

    fOutputEMC7->Add(fHistPhotoMismatch_EMC7);
    fOutputEMC7->Add(fHistPtAssoc_EMC7);
    fOutputEMC7->Add(fHistPtAssocMix_EMC7);
    fOutputEMC7->Add(fHistPtTag_EMC7);
    fOutputEMC7->Add(fHistTPCNClus_EMC7);
    fOutputEMC7->Add(fHistITSNClus_EMC7);
    fOutputEMC7->Add(fHistImpPar_EMC7);
    fOutputEMC7->Add(fHistImpParTag_EMC7);
    fOutputEMC7->Add(fHistNevents_EMC7);
    fOutputEMC7->Add(fHistPtSum_EMC7);
    fOutputEMC7->Add(fHistPtSumTag_EMC7);
    fOutputEMC7->Add(fHistPtSumEMC_EMC7);
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

    }
    for(Int_t i=0; i<3;i++){
        fOutputEMC7->Add(fHistDPhi300_1_EMC7[i]);
        fOutputEMC7->Add(fHistDPhi1_2_EMC7[i]);
        fOutputEMC7->Add(fHistDPhi2_4_EMC7[i]);
        fOutputEMC7->Add(fHistDPhi4_8_EMC7[i]);
        fOutputEMC7->Add(fHistDPhiMix300_1_EMC7[i]);
        fOutputEMC7->Add(fHistDPhiMix1_2_EMC7[i]);
        fOutputEMC7->Add(fHistDPhiMix2_4_EMC7[i]);
        fOutputEMC7->Add(fHistDPhiMix4_8_EMC7[i]);
    }
    fOutputEMC7->Add(fHistDPhi28_EMC7);
    fOutputEMC7->Add(fHistDPhiDEta28_EMC7);
    fOutputEMC7->Add(fHistDPhiMix28_EMC7);
    fOutputEMC7->Add(fHistDPhiDEtaMix28_EMC7);
    fOutputEMC7->Add(fHistDPhi18Spe_EMC7);

    fOutputEMCEGA->Add(fHistPhotoMismatch_EMCEGA);
    fOutputEMCEGA->Add(fHistPtAssoc_EMCEGA);
    fOutputEMCEGA->Add(fHistPtAssocMix_EMCEGA);
    fOutputEMCEGA->Add(fHistPtTag_EMCEGA);
    fOutputEMCEGA->Add(fHistTPCNClus_EMCEGA);
    fOutputEMCEGA->Add(fHistITSNClus_EMCEGA);
    fOutputEMCEGA->Add(fHistImpPar_EMCEGA);
    fOutputEMCEGA->Add(fHistImpParTag_EMCEGA);
    fOutputEMCEGA->Add(fHistNevents_EMCEGA);
    fOutputEMCEGA->Add(fHistPtSum_EMCEGA);
    fOutputEMCEGA->Add(fHistPtSumTag_EMCEGA);
    fOutputEMCEGA->Add(fHistPtSumEMC_EMCEGA);
    fOutputEMCEGA->Add(fHistEtaPhi_EMCEGA);
    fOutputEMCEGA->Add(fHistEtaPhiTag_EMCEGA);
    fOutputEMCEGA->Add(fHistInvMassElecLike_EMCEGA);
    fOutputEMCEGA->Add(fHistInvMassElecUnLike_EMCEGA);
    fOutputEMCEGA->Add(fHistOpAngElecLike_EMCEGA);
    fOutputEMCEGA->Add(fHistOpAngElecUnLike_EMCEGA);
    for(Int_t i=0; i<6;i++){
        fOutputEMCEGA->Add(fHistTPC_EMCTRD_EMCEGA[i]);

        fOutputEMCEGA->Add(fHistEMC_TPCTRD_EMCEGA[i]);

        fOutputEMCEGA->Add(fHistTRD_TPCEMC_EMCEGA[i]);

    }
    for(Int_t i=0; i<3;i++){
        fOutputEMCEGA->Add(fHistDPhi300_1_EMCEGA[i]);
        fOutputEMCEGA->Add(fHistDPhi1_2_EMCEGA[i]);
        fOutputEMCEGA->Add(fHistDPhi2_4_EMCEGA[i]);
        fOutputEMCEGA->Add(fHistDPhi4_8_EMCEGA[i]);
        fOutputEMCEGA->Add(fHistDPhiMix300_1_EMCEGA[i]);
        fOutputEMCEGA->Add(fHistDPhiMix1_2_EMCEGA[i]);
        fOutputEMCEGA->Add(fHistDPhiMix2_4_EMCEGA[i]);
        fOutputEMCEGA->Add(fHistDPhiMix4_8_EMCEGA[i]);
    }
    fOutputEMCEGA->Add(fHistDPhi28_EMCEGA);
    fOutputEMCEGA->Add(fHistDPhiDEta28_EMCEGA);
    fOutputEMCEGA->Add(fHistDPhiMix28_EMCEGA);
    fOutputEMCEGA->Add(fHistDPhiDEtaMix28_EMCEGA);
    fOutputEMCEGA->Add(fHistDPhi18Spe_EMCEGA);

    fOutputEMCJet->Add(fHistPhotoMismatch_EMCJet);
    fOutputEMCJet->Add(fHistPtAssoc_EMCJet);
    fOutputEMCJet->Add(fHistPtAssocMix_EMCJet);
    fOutputEMCJet->Add(fHistPtTag_EMCJet);
    fOutputEMCJet->Add(fHistTPCNClus_EMCJet);
    fOutputEMCJet->Add(fHistITSNClus_EMCJet);
    fOutputEMCJet->Add(fHistImpPar_EMCJet);
    fOutputEMCJet->Add(fHistImpParTag_EMCJet);
    fOutputEMCJet->Add(fHistNevents_EMCJet);
    fOutputEMCJet->Add(fHistPtSum_EMCJet);
    fOutputEMCJet->Add(fHistPtSumTag_EMCJet);
    fOutputEMCJet->Add(fHistPtSumEMC_EMCJet);
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

    }
    for(Int_t i=0; i<3;i++){
        fOutputEMCJet->Add(fHistDPhi300_1_EMCJet[i]);
        fOutputEMCJet->Add(fHistDPhi1_2_EMCJet[i]);
        fOutputEMCJet->Add(fHistDPhi2_4_EMCJet[i]);
        fOutputEMCJet->Add(fHistDPhi4_8_EMCJet[i]);
        fOutputEMCJet->Add(fHistDPhiMix300_1_EMCJet[i]);
        fOutputEMCJet->Add(fHistDPhiMix1_2_EMCJet[i]);
        fOutputEMCJet->Add(fHistDPhiMix2_4_EMCJet[i]);
        fOutputEMCJet->Add(fHistDPhiMix4_8_EMCJet[i]);
    }
    fOutputEMCJet->Add(fHistDPhi28_EMCJet);
    fOutputEMCJet->Add(fHistDPhiDEta28_EMCJet);
    fOutputEMCJet->Add(fHistDPhiMix28_EMCJet);
    fOutputEMCJet->Add(fHistDPhiDEtaMix28_EMCJet);
    fOutputEMCJet->Add(fHistDPhi18Spe_EMCJet);

    // NEW HISTO added to fOutput here
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMCEGA);
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
    AliAODEvent* aod = dynamic_cast<AliAODEvent*>(event);
    if (!aod) {
        AliError("Cannot get the AOD event");
        return;
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

    UInt_t fSelectMask = inputHandler->IsEventSelected();

    Bool_t isSelected = fSelectMask & (AliVEvent::kEMC7 | AliVEvent::kEMCEJE | AliVEvent::kEMCEGA);
    if(!isSelected){
        AliWarning("This is not an EMCal triggered event");
    }

    MBtrg = fSelectMask & AliVEvent::kAnyINT;
    EMC7trg = fSelectMask & AliVEvent::kEMC7;
    EMCEGAtrg = fSelectMask & AliVEvent::kEMCEGA;
    EMCJettrg = fSelectMask & AliVEvent::kEMCEJE;

    trigVal= EMCJettrg?EMCJE:EMC7trg?EMC7:EMCEGAtrg?EMCEGA:NONE;

    Int_t elecIDs[1000];
    Int_t elecCnt=0;

    AliPIDResponse* fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();

    if(!fPIDResponse){
        AliWarning("NULL PIDResponse");
    }

    TObjArray* trkArr = MakeTrkArr(aod);
    Bool_t delTrk = kTRUE;

    //__________________________End major event stuff_____________________________


    //Fill the histogram cataloguing # of events vs. events tagged
    if(MBtrg){
        fHistNevents_MB->Fill("Events",1);
    }
    switch(trigVal){
        case(EMC7):
            fHistNevents_EMC7->Fill("Events",1);
            break;

        case(EMCEGA):
            fHistNevents_EMCEGA->Fill("Events",1);
            break;

        case(EMCJE):
            fHistNevents_EMCJet->Fill("Events",1);
            break;
    }

    //Initialize energy variable and tagging flags
    Double_t PtSum=0;
    Double_t PtSumEMC=0;
    tagStrong=kFALSE;
    Bool_t tagEvt=kFALSE;

    //Initialize the # of tracks variable and the Eta Phi arrays
    Int_t ntracks=0;
    ntracks = aod->GetNumberOfTracks();

    std::vector<Double_t> Eta;
    std::vector<Double_t> Phi;

    fPool = fPoolMan->GetEventPool(ntracks, aod->GetPrimaryVertex()->GetZ());

    // Track loop for reconstructed event
    for(Int_t i = 0; i < ntracks; i++) {

        tagStrong=kFALSE;
        tagPhot=kFALSE;
        AliAODTrack* aodtrack = (AliAODTrack*)aod->GetTrack(i); // pointer to reconstructed to track       

        if(!aodtrack) { 
            AliError(Form("ERROR: Could not retrieve track %d",i)); 
            continue; 
        }

        //Fill TPCOnly track eta-phi

        if(aodtrack->IsTPCOnly()){
            fHistEtaPhiTPCOnly_MB->Fill(aodtrack->Eta(),aodtrack->Phi());
        }

        //Do hybrid track cuts
        if(!aodtrack->IsHybridGlobalConstrainedGlobal()){continue;}

        //Add this tracks energy to the running total
        PtSum=PtSum+aodtrack->Pt();
        if(aodtrack->Eta()<.7&&aodtrack->Eta()>-.7&&aodtrack->Phi()>80&&aodtrack->Phi()<180){
            PtSumEMC=PtSumEMC+aodtrack->Pt();
        }

        //Fill the Eta Phi arrays with this tracks Eta and Phi
        Eta.push_back(aodtrack->Eta());
        Phi.push_back(aodtrack->Phi());

        if(MBtrg){
            fHistEtaPhi_MB->Fill(aodtrack->Eta(),aodtrack->Phi());
        }
        switch(trigVal){
            case(EMC7):
                fHistEtaPhi_EMC7->Fill(aodtrack->Eta(),aodtrack->Phi());
                break;
            case(EMCEGA):
                fHistEtaPhi_EMCEGA->Fill(aodtrack->Eta(),aodtrack->Phi());
                break;
            case(EMCJE):
                fHistEtaPhi_EMCJet->Fill(aodtrack->Eta(),aodtrack->Phi());
                break;
        }
        //do Cut level histograms
        if(MBtrg){
            if(aodtrack->GetTPCncls()>0){
                fHistTPCNClus_MB->Fill(aodtrack->GetTPCncls());
            }
            fHistITSNClus_MB->Fill(aodtrack->GetNcls(0));
        }
        switch(trigVal){
            case(EMC7):
                if(aodtrack->GetTPCncls()>0){
                    fHistTPCNClus_EMC7->Fill(aodtrack->GetTPCncls());
                }
                fHistITSNClus_EMC7->Fill(aodtrack->GetNcls(0));
                break;

            case(EMCEGA):
                if(aodtrack->GetTPCncls()>0){
                    fHistTPCNClus_EMCEGA->Fill(aodtrack->GetTPCncls());
                }
                fHistITSNClus_EMCEGA->Fill(aodtrack->GetNcls(0));
                break;

            case(EMCJE):
                if(aodtrack->GetTPCncls()>0){
                    fHistTPCNClus_EMCJet->Fill(aodtrack->GetTPCncls());
                }
                fHistITSNClus_EMCJet->Fill(aodtrack->GetNcls(0)); 
                break;
        }

        //Impact parameter
        Float_t xy;
        Float_t z;

        xy=TMath::Sqrt(aodtrack->XAtDCA()*aodtrack->XAtDCA()+aodtrack->YAtDCA()*aodtrack->YAtDCA());

        if(MBtrg){
            fHistImpPar_MB->Fill(xy);
        }
        switch(trigVal){
            case(EMC7):
                fHistImpPar_EMC7->Fill(xy);
                break;

            case(EMCEGA):
                fHistImpPar_EMCEGA->Fill(xy);
                break;

            case(EMCJE):
                fHistImpPar_EMCJet->Fill(xy);
                break;
        }

        FillPhotoElecHistos(aod, aodtrack, fPIDResponse, i);

        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //If the track doesn't pass the cuts, move on to the next one
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if(trackCutsStrong){
            if(!fTrackCutsStrong->AcceptVTrack((AliVTrack*)aodtrack) || aodtrack->GetTPCsignalN()<80){continue;}
        }else{
            if(!fTrackCutsWeak->AcceptVTrack((AliVTrack*)aodtrack) || aodtrack->GetTPCsignalN()<80){continue;}
        }

        FillPIDHistos(aod, aodtrack, fPIDResponse);//Fill PID histos and set "tagStrong" boolean if this track satisfies cuts

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
            switch(trigVal){
                case(EMC7):
                    fHistImpParTag_EMC7->Fill(xy);
                    break;

                case(EMCEGA):
                    fHistImpParTag_EMCEGA->Fill(xy);
                    break;

                case(EMCJE):
                    fHistImpParTag_EMCJet->Fill(xy);
                    break;
            }

            //Pt distribution

            if(MBtrg){
                fHistPtTag_MB->Fill(aodtrack->Pt());
            }
            switch(trigVal){
                case(EMC7):
                    fHistPtTag_EMC7->Fill(aodtrack->Pt());
                    break;

                case(EMCEGA):
                    fHistPtTag_EMCEGA->Fill(aodtrack->Pt());
                    break;

                case(EMCJE):
                    fHistPtTag_EMCJet->Fill(aodtrack->Pt());
                    break;
            }
            FillDPhiHistos(aod, aodtrack, i);//Fill DPhi histos

            if(tagPhot){
                if(MBtrg){
                    fHistPhotoMismatch_MB->Fill(0);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistPhotoMismatch_EMC7->Fill(0);
                        break;
                    case(EMCEGA):
                        fHistPhotoMismatch_EMCEGA->Fill(0);
                        break;
                    case(EMCJE):
                        fHistPhotoMismatch_EMCJet->Fill(0);
                        break;
                }
            }




            if(!fPool){cout<<"No Pool for this event man\n"; continue;}

            fPool->PrintInfo();
            if(fPool->IsReady() ){
                FillMEDPhiHistos(aodtrack);
            }
            else{
                cout<<"Pool wasn't ready\n";
            }
        }//end if(tagStrong)

    }//end main track loop
    if(!fPool){
        cout<<"No pool exists, can't update it"<<'\n';
    }
    else{
        if(UseNonSignalEvents){
            if(trkArr){
                delTrk=kFALSE;
                fPool->UpdatePool(trkArr);
            }
        }else{
            if(tagEvt){
                if(trkArr){
                    delTrk=kFALSE;
                    fPool->UpdatePool(trkArr);
                }
            }
        }
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
        switch(trigVal){
            case(EMC7):
                fHistPtSumTag_EMC7->Fill(PtSum);
                break;
            case(EMCEGA):
                fHistPtSumTag_EMCEGA->Fill(PtSum);
                break;
            case(EMCJE):
                fHistPtSumTag_EMCJet->Fill(PtSum);
                break;
        }
    }else{
        if(MBtrg){
            fHistPtSum_MB->Fill(PtSum);
        }
        switch(trigVal){
            case(EMC7):
                fHistPtSum_EMC7->Fill(PtSum);
                break;
            case(EMCEGA):
                fHistPtSum_EMCEGA->Fill(PtSum);
                break;
            case(EMCJE):
                fHistPtSum_EMCJet->Fill(PtSum);
                break;
        }
    }

    if(MBtrg){
        fHistPtSumEMC_MB->Fill(PtSumEMC);
    }
    switch(trigVal){
        case(EMC7):
            fHistPtSumEMC_EMC7->Fill(PtSumEMC);
            break;
        case(EMCEGA):
            fHistPtSumEMC_EMCEGA->Fill(PtSumEMC);
            break;
        case(EMCJE):
            fHistPtSumEMC_EMCJet->Fill(PtSumEMC);
            break;
    }

    //Fill Nevent histos
    if(tagEvt){
        if(MBtrg){
            fHistNevents_MB->Fill("Events containing candidates",1);
        }  
        switch(trigVal){
            case(EMC7):
                fHistNevents_EMC7->Fill("Events containing candidates",1);
                break;
            case(EMCEGA):
                fHistNevents_EMCEGA->Fill("Events containing candidates",1);
                break;
            case(EMCJE):
                fHistNevents_EMCJet->Fill("Events containing candidates",1);
                break;
        }
    }

    //Fill the Eta Phi histograms
    if(tagEvt){
        for(Int_t i=0;i<Eta.size();i++){
            if(MBtrg){
                fHistEtaPhiTag_MB->Fill(Eta[i],Phi[i]);
            }
            switch(trigVal){
                case(EMC7):
                    fHistEtaPhiTag_EMC7->Fill(Eta[i],Phi[i]);
                    break;
                case(EMCEGA):
                    fHistEtaPhiTag_EMCEGA->Fill(Eta[i],Phi[i]);
                    break;
                case(EMCJE):
                    fHistEtaPhiTag_EMCJet->Fill(Eta[i],Phi[i]);
                    break;
            }
        }
    }

    if(trkArr&&delTrk){
        delete trkArr;
    }



    // NEW HISTO should be filled before this point, as PostData puts the
    // information for this iteration of the UserExec in the container
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMCEGA);
    PostData(4, fOutputEMCJet);
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
    Int_t caloId=aodtrack->GetEMCALcluster();

    if(caloId==-99999){
        return;
    }

    AliAODCaloCluster* tagEMCclus=aod->GetCaloCluster(caloId);

    if(tagEMCclus->E()>.5){
        EOP = tagEMCclus->E()/aodtrack->Pt();
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
    Double_t EMCcutLower[6] = {.85,.85,.85,.85,.85,.85};
    Double_t EMCcutHigher[6] = {1.3,1.3,1.3,1.3,1.3,1.3};

    for(Int_t i=0; i<6; i++){
        if(aodtrack->Pt()>ptLower[i]&&aodtrack->Pt()<ptUpper[i]){

            //TPC Plots
            //EMC+TRD cuts
            if(EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]&&elecLikeTRD[0]>TRDcut){

                if(MBtrg){
                    fHistTPC_EMCTRD_MB[i]->Fill(aodtrack->Pt(), nSigmaTPC);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistTPC_EMCTRD_EMC7[i]->Fill(aodtrack->Pt(), nSigmaTPC);
                        break;

                    case(EMCEGA):
                        fHistTPC_EMCTRD_EMCEGA[i]->Fill(aodtrack->Pt(), nSigmaTPC);
                        break;

                    case(EMCJE):
                        fHistTPC_EMCTRD_EMCJet[i]->Fill(aodtrack->Pt(), nSigmaTPC);
                        break;
                }
            }

            //EMC Plots

            //TPC+TRD cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut){

                if(!applySSCuts){
                    if(MBtrg){
                        fHistEMC_TPCTRD_MB[i]->Fill(EOP);
                    }
                    switch(trigVal){
                        case(EMC7):
                            fHistEMC_TPCTRD_EMC7[i]->Fill(EOP);
                            break;

                        case(EMCEGA):
                            fHistEMC_TPCTRD_EMCEGA[i]->Fill(EOP);
                            break;

                        case(EMCJE):
                            fHistEMC_TPCTRD_EMCJet[i]->Fill(EOP);
                            break;
                    }
                }else{
                    if(MBtrg){
                        fHistEMC_TPCTRD_MB[i]->Fill(EOP);
                    }
                    switch(trigVal){
                        case(EMC7):
                            fHistEMC_TPCTRD_EMC7[i]->Fill(EOP);
                            break;

                        case(EMCEGA):
                            fHistEMC_TPCTRD_EMCEGA[i]->Fill(EOP);
                            break;

                        case(EMCJE):
                            fHistEMC_TPCTRD_EMCJet[i]->Fill(EOP);
                            break;
                    }
                }
            }

            //TRD Plots

            //TPC+EMC cuts
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&EOP<EMCcutHigher[i]&&EOP>EMCcutLower[i]){

                if(MBtrg){
                    fHistTRD_TPCEMC_MB[i]->Fill(aodtrack->Pt(), elecLikeTRD[0]);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistTRD_TPCEMC_EMC7[i]->Fill(aodtrack->Pt(), elecLikeTRD[0]);
                        break;

                    case(EMCEGA):
                        fHistTRD_TPCEMC_EMCEGA[i]->Fill(aodtrack->Pt(), elecLikeTRD[0]);
                        break;

                    case(EMCJE):
                        fHistTRD_TPCEMC_EMCJet[i]->Fill(aodtrack->Pt(), elecLikeTRD[0]);
                        break;
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
            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut&&EOP<EMCcutHigher[0]&&EOP>EMCcutLower[0]){

                tagStrong=kTRUE;

            }

        }

        else{

            if(nSigmaTPC<TPCcut&&nSigmaTPC>-TPCcut&&elecLikeTRD[0]>TRDcut&&EOP<EMCcutHigher[5]&&EOP>EMCcutLower[5]){

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
    PostData(3, fOutputEMCEGA);
    PostData(4, fOutputEMCJet);
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
        if(!aodtrackassoc->IsHybridGlobalConstrainedGlobal()&&!aodtrackassoc->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)){continue;}

        //Pt distribution
        if(MBtrg){
            fHistPtAssoc_MB->Fill(aodtrackassoc->Pt());
        }
        switch(trigVal){
            case(EMC7):
                fHistPtAssoc_EMC7->Fill(aodtrackassoc->Pt());
                break;

            case(EMCEGA):
                fHistPtAssoc_EMCEGA->Fill(aodtrackassoc->Pt());
                break;

            case(EMCJE):
                fHistPtAssoc_EMCJet->Fill(aodtrackassoc->Pt());
                break;
        }

        //Fill Delta Phi variable and correct for periodicity
        Double_t DPhi=aodtrackassoc->Phi()-aodtrack->Phi();

        if(DPhi<-TMath::Pi()/2){DPhi=TMath::Abs(2*TMath::Pi()+DPhi);}

        if(DPhi>3*TMath::Pi()/2){DPhi=-TMath::Abs(2*TMath::Pi()-DPhi);}

        Double_t DEta=aodtrackassoc->Eta()-aodtrack->Eta();


        Int_t PID=0;
        cout<<"most probPID"<<AliAODTrack::kElectron<<":"<<aodtrackassoc->GetMostProbablePID()<<'\n';
        switch(aodtrackassoc->GetMostProbablePID()){
            case AliAODTrack::kElectron:
                PID=1;
                break;
            case AliAODTrack::kMuon:
                PID=2;
                break;
            case AliAODTrack::kPion:
                PID=3;
                break;
            case AliAODTrack::kKaon:
                PID=4;
                break;
            case AliAODTrack::kProton:
                PID=5;
                break;
            case AliAODTrack::kDeuteron:
                PID=6;
                break;
            case AliAODTrack::kUnknown:
                PID=0;
                break;
        }

        if(MBtrg){
            fHistDPhi18Spe_MB->Fill(DPhi, PID);
        }
        switch(trigVal){
            case(EMC7):
                fHistDPhi18Spe_EMC7->Fill(DPhi, PID);
                break;
            case(EMCEGA):
                fHistDPhi18Spe_EMCEGA->Fill(DPhi, PID);
                break;
            case(EMCJE):
                fHistDPhi18Spe_EMCJet->Fill(DPhi, PID);
                break;
        }

        if(aodtrackassoc->GetType()==AliAODTrack::kPrimary){

            if(aodtrack->Pt()<8&&aodtrack->Pt()>1&&aodtrackassoc->Pt()>0.3)
            {
                fHistTestDPhi18NoSec->Fill(DPhi);
            }

            fHistTestDPhiSpeNoSec->Fill(DPhi, PID);

        }else{
            if(aodtrack->Pt()<8&&aodtrack->Pt()>1&&aodtrackassoc->Pt()>0.3)
            {
                fHistTestDPhi18Sec->Fill(DPhi);
            }
        }



        //Fill DPhi by Type here
        if(aodtrackassoc->IsHybridGlobalConstrainedGlobal()) fHistTestDPhiType->Fill(DPhi, 1);

        if((aodtrackassoc->IsGlobalConstrained())) fHistTestDPhiType->Fill(DPhi, 2);

        if((aodtrackassoc->IsTPCOnly())) fHistTestDPhiType->Fill(DPhi, 3);


        if(DPhi<0.1&&DPhi>-0.1&&DEta<0.1&&DEta>-0.1){

            fHistTestPt->Fill(aodtrackassoc->Pt());
            fHistTestOGDPhi->Fill(aodtrackassoc->Phi()-aodtrack->Phi());
            fHistTestDCA->Fill(aodtrackassoc->DCA());
            fHistTestTPCdEdx->Fill(aodtrackassoc->Pt(), aodtrackassoc->GetTPCsignal());


            Int_t partOneID = 0;
            Int_t partTwoID = 0;

            if(aodtrackassoc->GetMostProbablePID()==AliAODTrack::kElectron) partOneID=1;

            if(aodtrackassoc->GetMostProbablePID()==AliAODTrack::kPion) partOneID=2;


            for(Int_t k=0;k<ntracks;k++){
                if(i==k || j==k){continue;}

                AliAODTrack* aodtrackassoc2 = (AliAODTrack*)aod->GetTrack(k);

                if(aodtrackassoc2->GetMostProbablePID()==AliAODTrack::kElectron) partTwoID=1;

                if(aodtrackassoc2->GetMostProbablePID()==AliAODTrack::kPion) partTwoID=2;

                Double_t DPhi=aodtrackassoc2->Phi()-aodtrack->Phi();

                if(DPhi<-TMath::Pi()/2){DPhi=TMath::Abs(2*TMath::Pi()+DPhi);}

                if(DPhi>3*TMath::Pi()/2){DPhi=-TMath::Abs(2*TMath::Pi()-DPhi);}

                Double_t DEta=aodtrackassoc2->Eta()-aodtrack->Eta();

                Double_t ElecMass=.000511;

                Double_t PionMass=.139;
                //fill inv mass plot

                if(partOneID==1&&partTwoID==1&&aodtrackassoc->Charge()==aodtrackassoc2->Charge()&&aodtrackassoc->Charge()!=0) 
                {
                    Double_t assocE1=TMath::Sqrt(aodtrackassoc->P()*aodtrackassoc->P()+ElecMass*ElecMass);
                    Double_t assocE2=TMath::Sqrt(aodtrackassoc2->P()*aodtrackassoc2->P()+ElecMass*ElecMass);

                    TLorentzVector assoc1(aodtrackassoc->Px(), aodtrackassoc->Py(), aodtrackassoc->Pz(), assocE1);
                    TLorentzVector assoc2(aodtrackassoc2->Px(), aodtrackassoc2->Py(), aodtrackassoc2->Pz(), assocE2);

                    Double_t InvMass=(assoc1+assoc2).M();

                    fHistTestInvMassElecLike->Fill(InvMass);   
                }

                if(partOneID==1&&partTwoID==1&&aodtrackassoc->Charge()!=aodtrackassoc2->Charge()&&aodtrackassoc->Charge()!=0&&aodtrackassoc2->Charge()!=0) 
                {
                    Double_t assocE1=TMath::Sqrt(aodtrackassoc->P()*aodtrackassoc->P()+ElecMass*ElecMass);
                    Double_t assocE2=TMath::Sqrt(aodtrackassoc2->P()*aodtrackassoc2->P()+ElecMass*ElecMass);

                    TLorentzVector assoc1(aodtrackassoc->Px(), aodtrackassoc->Py(), aodtrackassoc->Pz(), assocE1);
                    TLorentzVector assoc2(aodtrackassoc2->Px(), aodtrackassoc2->Py(), aodtrackassoc2->Pz(), assocE2);

                    Double_t InvMass=(assoc1+assoc2).M();

                    fHistTestInvMassElecUnLike->Fill(InvMass);   
                }

                if(partOneID==2&&partTwoID==2&&aodtrackassoc->Charge()==aodtrackassoc2->Charge()&&aodtrackassoc->Charge()!=0) 
                {
                    Double_t assocE1=TMath::Sqrt(aodtrackassoc->P()*aodtrackassoc->P()+PionMass*PionMass);
                    Double_t assocE2=TMath::Sqrt(aodtrackassoc2->P()*aodtrackassoc2->P()+PionMass*PionMass);

                    TLorentzVector assoc1(aodtrackassoc->Px(), aodtrackassoc->Py(), aodtrackassoc->Pz(), assocE1);
                    TLorentzVector assoc2(aodtrackassoc2->Px(), aodtrackassoc2->Py(), aodtrackassoc2->Pz(), assocE2);

                    Double_t InvMass=(assoc1+assoc2).M();

                    fHistTestInvMassPionLike->Fill(InvMass);   
                }

                if(partOneID==2&&partTwoID==2&&aodtrackassoc->Charge()!=aodtrackassoc2->Charge()&&aodtrackassoc->Charge()!=0&&aodtrackassoc2->Charge()!=0) 
                {
                    Double_t assocE1=TMath::Sqrt(aodtrackassoc->P()*aodtrackassoc->P()+PionMass*PionMass);
                    Double_t assocE2=TMath::Sqrt(aodtrackassoc2->P()*aodtrackassoc2->P()+PionMass*PionMass);

                    TLorentzVector assoc1(aodtrackassoc->Px(), aodtrackassoc->Py(), aodtrackassoc->Pz(), assocE1);
                    TLorentzVector assoc2(aodtrackassoc2->Px(), aodtrackassoc2->Py(), aodtrackassoc2->Pz(), assocE2);

                    Double_t InvMass=(assoc1+assoc2).M();

                    fHistTestInvMassPionUnLike->Fill(InvMass);   
                }



            }
            //Fill other plots
            Int_t cid = aodtrackassoc->GetEMCALcluster();
            if(cid > 0){
                AliAODCaloCluster *aodcl = aod->GetCaloCluster(cid);

                fHistTestEMCEnergy->Fill(aodcl->E());
                fHistTestEOP->Fill(aodcl->E()/aodtrackassoc->Pt());
            }
            else{cout<<"No EMCal cluster for this anamolous Peak\n";}




        }

        if(PID==1||PID==2||PID==0){continue;}

        //candidate 1<pt<2
        if(aodtrack->Pt()>1&&aodtrack->Pt()<2){

            

            //300MeV
            if(aodtrackassoc->Pt()>.3&&aodtrackassoc->Pt()<1){
                if(MBtrg){
                    fHistDPhi300_1_MB[0]->Fill(DPhi);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistDPhi300_1_EMC7[0]->Fill(DPhi);
                        break;
                    case(EMCEGA):
                        fHistDPhi300_1_EMCEGA[0]->Fill(DPhi);
                        break;
                    case(EMCJE):
                        fHistDPhi300_1_EMCJet[0]->Fill(DPhi);
                        break;
                }
            }

            //1GeV
            if(aodtrackassoc->Pt()>1&&aodtrackassoc->Pt()<2){
                if(MBtrg){
                    fHistDPhi1_2_MB[0]->Fill(DPhi);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistDPhi1_2_EMC7[0]->Fill(DPhi);
                        break;
                    case(EMCEGA):
                        fHistDPhi1_2_EMCEGA[0]->Fill(DPhi);
                        break;
                    case(EMCJE):
                        fHistDPhi1_2_EMCJet[0]->Fill(DPhi);
                        break;
                }
            }

            //2GeV
            if(aodtrackassoc->Pt()>2&&aodtrackassoc->Pt()<4){
                if(MBtrg){
                    fHistDPhi2_4_MB[0]->Fill(DPhi);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistDPhi2_4_EMC7[0]->Fill(DPhi);
                        break;
                    case(EMCEGA):
                        fHistDPhi2_4_EMCEGA[0]->Fill(DPhi);
                        break;
                    case(EMCJE):
                        fHistDPhi2_4_EMCJet[0]->Fill(DPhi);
                        break;
                }
            }

            //4GeV
            if(aodtrackassoc->Pt()>4&&aodtrackassoc->Pt()<8){
                if(MBtrg){
                    fHistDPhi4_8_MB[0]->Fill(DPhi);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistDPhi4_8_EMC7[0]->Fill(DPhi);
                        break;
                    case(EMCEGA):
                        fHistDPhi4_8_EMCEGA[0]->Fill(DPhi);
                        break;
                    case(EMCJE):
                        fHistDPhi4_8_EMCJet[0]->Fill(DPhi);
                        break;
                }
            }

           
        }

        //candidate 2<pt<4
        if(aodtrack->Pt()>2&&aodtrack->Pt()<4){

            

            //300MeV
            if(aodtrackassoc->Pt()>.3&&aodtrackassoc->Pt()<1){
                if(MBtrg){
                    fHistDPhi300_1_MB[1]->Fill(DPhi);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistDPhi300_1_EMC7[1]->Fill(DPhi);
                        break;
                    case(EMCEGA):
                        fHistDPhi300_1_EMCEGA[1]->Fill(DPhi);
                        break;
                    case(EMCJE):
                        fHistDPhi300_1_EMCJet[1]->Fill(DPhi);
                        break;
                }
            }

            //1GeV
            if(aodtrackassoc->Pt()>1&&aodtrackassoc->Pt()<2){
                if(MBtrg){
                    fHistDPhi1_2_MB[1]->Fill(DPhi);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistDPhi1_2_EMC7[1]->Fill(DPhi);
                        break;
                    case(EMCEGA):
                        fHistDPhi1_2_EMCEGA[1]->Fill(DPhi);
                        break;
                    case(EMCJE):
                        fHistDPhi1_2_EMCJet[1]->Fill(DPhi);
                        break;
                }
            }

            //2GeV
            if(aodtrackassoc->Pt()>2&&aodtrackassoc->Pt()<4){
                if(MBtrg){
                    fHistDPhi2_4_MB[1]->Fill(DPhi);
                    fHistDPhi28_MB->Fill(DPhi);
                    fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistDPhi2_4_EMC7[1]->Fill(DPhi);
                        fHistDPhi28_EMC7->Fill(DPhi);
                        fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                        break;
                    case(EMCEGA):
                        fHistDPhi2_4_EMCEGA[1]->Fill(DPhi);
                        fHistDPhi28_EMCEGA->Fill(DPhi);
                        fHistDPhiDEta28_EMCEGA->Fill(DPhi, DEta);
                        break;
                    case(EMCJE):
                        fHistDPhi2_4_EMCJet[1]->Fill(DPhi);
                        fHistDPhi28_EMCJet->Fill(DPhi);
                        fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                        break;
                }
            }

            //4GeV
            if(aodtrackassoc->Pt()>4&&aodtrackassoc->Pt()<8){
                if(MBtrg){
                    fHistDPhi4_8_MB[1]->Fill(DPhi);
                    fHistDPhi28_MB->Fill(DPhi);
                    fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistDPhi4_8_EMC7[1]->Fill(DPhi);
                        fHistDPhi28_EMC7->Fill(DPhi);
                        fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                        break;
                    case(EMCEGA):
                        fHistDPhi4_8_EMCEGA[1]->Fill(DPhi);
                        fHistDPhi28_EMCEGA->Fill(DPhi);
                        fHistDPhiDEta28_EMCEGA->Fill(DPhi, DEta);
                        break;
                    case(EMCJE):
                        fHistDPhi4_8_EMCJet[1]->Fill(DPhi);
                        fHistDPhi28_EMCJet->Fill(DPhi);
                        fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                        break;
                }
            }

            
        }

        //candidate 4<pt<8
        if(aodtrack->Pt()>4&&aodtrack->Pt()<8){

            

            //300MeV
            if(aodtrackassoc->Pt()>.3&&aodtrackassoc->Pt()<1){
                if(MBtrg){
                    fHistDPhi300_1_MB[2]->Fill(DPhi);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistDPhi300_1_EMC7[2]->Fill(DPhi);
                        break;
                    case(EMCEGA):
                        fHistDPhi300_1_EMCEGA[2]->Fill(DPhi);
                        break;
                    case(EMCJE):
                        fHistDPhi300_1_EMCJet[2]->Fill(DPhi);
                        break;
                }
            }

            //1GeV
            if(aodtrackassoc->Pt()>1&&aodtrackassoc->Pt()<2){
                if(MBtrg){
                    fHistDPhi1_2_MB[2]->Fill(DPhi);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistDPhi1_2_EMC7[2]->Fill(DPhi);
                        break;
                    case(EMCEGA):
                        fHistDPhi1_2_EMCEGA[2]->Fill(DPhi);
                        break;
                    case(EMCJE):
                        fHistDPhi1_2_EMCJet[2]->Fill(DPhi);
                        break;
                }
            }

            //2GeV
            if(aodtrackassoc->Pt()>2&&aodtrackassoc->Pt()<4){
                if(MBtrg){
                    fHistDPhi2_4_MB[2]->Fill(DPhi);
                    fHistDPhi28_MB->Fill(DPhi);
                    fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistDPhi2_4_EMC7[2]->Fill(DPhi);
                        fHistDPhi28_EMC7->Fill(DPhi);
                        fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                        break;
                    case(EMCEGA):
                        fHistDPhi2_4_EMCEGA[2]->Fill(DPhi);
                        fHistDPhi28_EMCEGA->Fill(DPhi);
                        fHistDPhiDEta28_EMCEGA->Fill(DPhi, DEta);
                        break;
                    case(EMCJE):
                        fHistDPhi2_4_EMCJet[2]->Fill(DPhi);
                        fHistDPhi28_EMCJet->Fill(DPhi);
                        fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                        break;
                }
            }

            //4GeV
            if(aodtrackassoc->Pt()>4&&aodtrackassoc->Pt()<8){
                if(MBtrg){
                    fHistDPhi4_8_MB[2]->Fill(DPhi);
                    fHistDPhi28_MB->Fill(DPhi);
                    fHistDPhiDEta28_MB->Fill(DPhi, DEta);
                }
                switch(trigVal){
                    case(EMC7):
                        fHistDPhi4_8_EMC7[2]->Fill(DPhi);
                        fHistDPhi28_EMC7->Fill(DPhi);
                        fHistDPhiDEta28_EMC7->Fill(DPhi, DEta);
                        break;
                    case(EMCEGA):
                        fHistDPhi4_8_EMCEGA[2]->Fill(DPhi);
                        fHistDPhi28_EMCEGA->Fill(DPhi);
                        fHistDPhiDEta28_EMCEGA->Fill(DPhi, DEta);
                        break;
                    case(EMCJE):
                        fHistDPhi4_8_EMCJet[2]->Fill(DPhi);
                        fHistDPhi28_EMCJet->Fill(DPhi);
                        fHistDPhiDEta28_EMCJet->Fill(DPhi, DEta);
                        break;
                }
            }

            
        }

    }
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMCEGA);
    PostData(4, fOutputEMCJet);
    return;
}

void AliAnalysisTaskPSHFE::SetElectronTrackCuts(Bool_t trkCutBool){
    trackCutsStrong=trkCutBool;
}

void AliAnalysisTaskPSHFE::SetSSCutBool(Bool_t SSCutBool){
    applySSCuts=SSCutBool;
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

        if(aodtrack->Charge()==aodtrackassoc->Charge()&&aodtrack->Charge()!=0){

            if(MBtrg){
                fHistInvMassElecLike_MB->Fill(InvMass);
                fHistOpAngElecLike_MB->Fill(OpAng);
            }
            switch(trigVal){
                case(EMC7):
                    fHistInvMassElecLike_EMC7->Fill(InvMass);
                    fHistOpAngElecLike_EMC7->Fill(OpAng);
                    break;

                case(EMCEGA):
                    fHistInvMassElecLike_EMCEGA->Fill(InvMass);
                    fHistOpAngElecLike_EMCEGA->Fill(OpAng);
                    break;

                case(EMCJE):
                    fHistInvMassElecLike_EMCJet->Fill(InvMass);
                    fHistOpAngElecLike_EMCJet->Fill(OpAng);
                    break;
            }

        }else if(aodtrack->Charge()!=aodtrackassoc->Charge()&&aodtrack->Charge()!=0&&aodtrackassoc->Charge()!=0){
            if(InvMass<0.1&&OpAng<0.1){tagPhot=kTRUE;}
            if(MBtrg){
                fHistInvMassElecUnLike_MB->Fill(InvMass);
                fHistOpAngElecUnLike_MB->Fill(OpAng);
            }
            switch(trigVal){
                case(EMC7):
                    fHistInvMassElecUnLike_EMC7->Fill(InvMass);
                    fHistOpAngElecUnLike_EMC7->Fill(OpAng);
                    break;

                case(EMCEGA):
                    fHistInvMassElecUnLike_EMCEGA->Fill(InvMass);
                    fHistOpAngElecUnLike_EMCEGA->Fill(OpAng);
                    break;

                case(EMCJE):
                    fHistInvMassElecUnLike_EMCJet->Fill(InvMass);
                    fHistOpAngElecUnLike_EMCJet->Fill(OpAng);
                    break;
            }
        }
    }//end loop over tracks

    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMCEGA);
    PostData(4, fOutputEMCJet);
    return;
}

TObjArray* AliAnalysisTaskPSHFE::MakeTrkArr(AliAODEvent *aod)
{
    if(!aod){AliWarning("Invalid AOD Event");}
    Int_t nTracks = aod->GetNumberOfTracks();
    TObjArray* accTracks = new TObjArray;
    accTracks->SetOwner();

    for(Int_t i=0;i<nTracks;i++){
        AliAODTrack *aodtrack = (AliAODTrack*)aod->GetTrack(i);

        if(!aodtrack){
            continue;
        }

        if(aodtrack->IsHybridGlobalConstrainedGlobal()){
            continue;
        }
        AliAODTrack* temp = new AliAODTrack(*aodtrack);
        accTracks->Add(temp);
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
                fHistPtAssocMix_MB->Fill(aodtrackassoc->Pt());
            }
            switch(trigVal){
                case(EMC7):
                    fHistPtAssocMix_EMC7->Fill(aodtrackassoc->Pt());
                    break;

                case(EMCEGA):
                    fHistPtAssocMix_EMCEGA->Fill(aodtrackassoc->Pt());
                    break;

                case(EMCJE):
                    fHistPtAssocMix_EMCJet->Fill(aodtrackassoc->Pt());
                    break;
            }

            //Fill Delta Phi variable and correct for periodicity
            Double_t DPhi=aodtrackassoc->Phi()-aodtrack->Phi();

            if(DPhi<-TMath::Pi()/2){DPhi=TMath::Abs(2*TMath::Pi()+DPhi);}

            if(DPhi>3*TMath::Pi()/2){DPhi=-TMath::Abs(2*TMath::Pi()-DPhi);}

            Double_t DEta=aodtrackassoc->Eta()-aodtrack->Eta();

            //candidate 1<pt<2
            if(aodtrack->Pt()>1&&aodtrack->Pt()<2){

                

                //300MeV
                if(aodtrackassoc->Pt()>.3&&aodtrackassoc->Pt()<1){
                    if(MBtrg){
                        fHistDPhiMix300_1_MB[0]->Fill(DPhi);
                    }
                    switch(trigVal){
                        case(EMC7):
                            fHistDPhiMix300_1_EMC7[0]->Fill(DPhi);
                            break;
                        case(EMCEGA):
                            fHistDPhiMix300_1_EMCEGA[0]->Fill(DPhi);
                            break;
                        case(EMCJE):
                            fHistDPhiMix300_1_EMCJet[0]->Fill(DPhi);
                            break;
                    }
                }

                //1GeV
                if(aodtrackassoc->Pt()>1&&aodtrackassoc->Pt()<2){
                    if(MBtrg){
                        fHistDPhiMix1_2_MB[0]->Fill(DPhi);
                    }
                    switch(trigVal){
                        case(EMC7):
                            fHistDPhiMix1_2_EMC7[0]->Fill(DPhi);
                            break;
                        case(EMCEGA):
                            fHistDPhiMix1_2_EMCEGA[0]->Fill(DPhi);
                            break;
                        case(EMCJE):
                            fHistDPhiMix1_2_EMCJet[0]->Fill(DPhi);
                            break;
                    }
                }

                //2GeV
                if(aodtrackassoc->Pt()>2&&aodtrackassoc->Pt()<4){
                    if(MBtrg){
                        fHistDPhiMix2_4_MB[0]->Fill(DPhi);
                    }
                    switch(trigVal){                       
                        case(EMC7):
                            fHistDPhiMix2_4_EMC7[0]->Fill(DPhi);
                            break;
                        case(EMCEGA):
                            fHistDPhiMix2_4_EMCEGA[0]->Fill(DPhi);
                            break;
                        case(EMCJE):
                            fHistDPhiMix2_4_EMCJet[0]->Fill(DPhi);
                            break;
                    }
                }

                //4GeV
                if(aodtrackassoc->Pt()>4&&aodtrackassoc->Pt()<8){
                    if(MBtrg){
                        fHistDPhiMix4_8_MB[0]->Fill(DPhi);
                    }
                    switch(trigVal){                        
                        case(EMC7):
                            fHistDPhiMix4_8_EMC7[0]->Fill(DPhi);
                            break;
                        case(EMCEGA):
                            fHistDPhiMix4_8_EMCEGA[0]->Fill(DPhi);
                            break;
                        case(EMCJE):
                            fHistDPhiMix4_8_EMCJet[0]->Fill(DPhi);
                            break;
                    }
                }

                
            }

            //candidate 2<pt<4
            if(aodtrack->Pt()>2&&aodtrack->Pt()<4){

                

                //300MeV
                if(aodtrackassoc->Pt()>.3&&aodtrackassoc->Pt()<1){
                    if(MBtrg){
                        fHistDPhiMix300_1_MB[1]->Fill(DPhi);
                    }
                    switch(trigVal){
                        case(EMC7):
                            fHistDPhiMix300_1_EMC7[1]->Fill(DPhi);
                            break;
                        case(EMCEGA):
                            fHistDPhiMix300_1_EMCEGA[1]->Fill(DPhi);
                            break;
                        case(EMCJE):
                            fHistDPhiMix300_1_EMCJet[1]->Fill(DPhi);
                            break;
                    }
                }

                //1GeV
                if(aodtrackassoc->Pt()>1&&aodtrackassoc->Pt()<2){
                    if(MBtrg){
                        fHistDPhiMix1_2_MB[1]->Fill(DPhi);
                    }
                    switch(trigVal){
                        case(EMC7):
                            fHistDPhiMix1_2_EMC7[1]->Fill(DPhi);
                            break;
                        case(EMCEGA):
                            fHistDPhiMix1_2_EMCEGA[1]->Fill(DPhi);
                            break;
                        case(EMCJE):
                            fHistDPhiMix1_2_EMCJet[1]->Fill(DPhi);
                            break;
                    }
                }

                //2GeV
                if(aodtrackassoc->Pt()>2&&aodtrackassoc->Pt()<4){
                    if(MBtrg){
                        fHistDPhiMix2_4_MB[1]->Fill(DPhi);
                        fHistDPhiMix28_MB->Fill(DPhi);
                        fHistDPhiDEtaMix28_MB->Fill(DPhi, DEta);
                    }
                    switch(trigVal){
                        case(EMC7):
                            fHistDPhiMix2_4_EMC7[1]->Fill(DPhi);
                            fHistDPhiMix28_EMC7->Fill(DPhi);
                            fHistDPhiDEtaMix28_EMC7->Fill(DPhi, DEta);
                            break;
                        case(EMCEGA):
                            fHistDPhiMix2_4_EMCEGA[1]->Fill(DPhi);
                            fHistDPhiMix28_EMCEGA->Fill(DPhi);
                            fHistDPhiDEtaMix28_EMCEGA->Fill(DPhi, DEta);
                            break;
                        case(EMCJE):
                            fHistDPhiMix2_4_EMCJet[1]->Fill(DPhi);
                            fHistDPhiMix28_EMCJet->Fill(DPhi);
                            fHistDPhiDEtaMix28_EMCJet->Fill(DPhi, DEta);
                            break;
                    }
                }

                //4GeV
                if(aodtrackassoc->Pt()>4&&aodtrackassoc->Pt()<8){
                    if(MBtrg){
                        fHistDPhiMix4_8_MB[1]->Fill(DPhi);
                        fHistDPhiMix28_MB->Fill(DPhi);
                        fHistDPhiDEtaMix28_MB->Fill(DPhi, DEta);
                    }
                    switch(trigVal){
                        case(EMC7):
                            fHistDPhiMix4_8_EMC7[1]->Fill(DPhi);
                            fHistDPhiMix28_EMC7->Fill(DPhi);
                            fHistDPhiDEtaMix28_EMC7->Fill(DPhi, DEta);
                            break;
                        case(EMCEGA):
                            fHistDPhiMix4_8_EMCEGA[1]->Fill(DPhi);
                            fHistDPhiMix28_EMCEGA->Fill(DPhi);
                            fHistDPhiDEtaMix28_EMCEGA->Fill(DPhi, DEta);
                            break;
                        case(EMCJE):
                            fHistDPhiMix4_8_EMCJet[1]->Fill(DPhi);
                            fHistDPhiMix28_EMCJet->Fill(DPhi);
                            fHistDPhiDEtaMix28_EMCJet->Fill(DPhi, DEta);
                            break;
                    }
                }

                           
            }

            //candidate 4<pt<8
            if(aodtrack->Pt()>4&&aodtrack->Pt()<8){

                

                //300MeV
                if(aodtrackassoc->Pt()>.3&&aodtrackassoc->Pt()<1){
                    if(MBtrg){
                        fHistDPhiMix300_1_MB[2]->Fill(DPhi);
                    }
                    switch(trigVal){
                        case(EMC7):
                            fHistDPhiMix300_1_EMC7[2]->Fill(DPhi);
                            break;
                        case(EMCEGA):
                            fHistDPhiMix300_1_EMCEGA[2]->Fill(DPhi);
                            break;
                        case(EMCJE):
                            fHistDPhiMix300_1_EMCJet[2]->Fill(DPhi);
                            break;
                    }
                }

                //1GeV
                if(aodtrackassoc->Pt()>1&&aodtrackassoc->Pt()<2){
                    if(MBtrg){
                        fHistDPhiMix1_2_MB[2]->Fill(DPhi);
                    }
                    switch(trigVal){
                        case(EMC7):
                            fHistDPhiMix1_2_EMC7[2]->Fill(DPhi);
                            break;
                        case(EMCEGA):
                            fHistDPhiMix1_2_EMCEGA[2]->Fill(DPhi);
                            break;
                        case(EMCJE):
                            fHistDPhiMix1_2_EMCJet[2]->Fill(DPhi);
                            break;
                    }
                }

                //2GeV
                if(aodtrackassoc->Pt()>2&&aodtrackassoc->Pt()<4){
                    if(MBtrg){
                        fHistDPhiMix2_4_MB[2]->Fill(DPhi);
                        fHistDPhiMix28_MB->Fill(DPhi);
                        fHistDPhiDEtaMix28_MB->Fill(DPhi, DEta);
                    }
                    switch(trigVal){
                        case(EMC7):
                            fHistDPhiMix2_4_EMC7[2]->Fill(DPhi);
                            fHistDPhiMix28_EMC7->Fill(DPhi);
                            fHistDPhiDEtaMix28_EMC7->Fill(DPhi, DEta);
                            break;
                        case(EMCEGA):
                            fHistDPhiMix2_4_EMCEGA[2]->Fill(DPhi);
                            fHistDPhiMix28_EMCEGA->Fill(DPhi);
                            fHistDPhiDEtaMix28_EMCEGA->Fill(DPhi, DEta);
                            break;
                        case(EMCJE):
                            fHistDPhiMix2_4_EMCJet[2]->Fill(DPhi);
                            fHistDPhiMix28_EMCJet->Fill(DPhi);
                            fHistDPhiDEtaMix28_EMCJet->Fill(DPhi, DEta);
                            break;
                    }
                }

                //4GeV
                if(aodtrackassoc->Pt()>4&&aodtrackassoc->Pt()<8){
                    if(MBtrg){
                        fHistDPhiMix4_8_MB[2]->Fill(DPhi);
                        fHistDPhiMix28_MB->Fill(DPhi);
                        fHistDPhiDEtaMix28_MB->Fill(DPhi, DEta);
                    }
                    switch(trigVal){
                        case(EMC7):
                            fHistDPhiMix4_8_EMC7[2]->Fill(DPhi);
                            fHistDPhiMix28_EMC7->Fill(DPhi);
                            fHistDPhiDEtaMix28_EMC7->Fill(DPhi, DEta);
                            break;
                        case(EMCEGA):
                            fHistDPhiMix4_8_EMCEGA[2]->Fill(DPhi);
                            fHistDPhiMix28_EMCEGA->Fill(DPhi);
                            fHistDPhiDEtaMix28_EMCEGA->Fill(DPhi, DEta);
                            break;
                        case(EMCJE):
                            fHistDPhiMix4_8_EMCJet[2]->Fill(DPhi);
                            fHistDPhiMix28_EMCJet->Fill(DPhi);
                            fHistDPhiDEtaMix28_EMCJet->Fill(DPhi, DEta);
                            break;
                    }
                }

                
            }
        }
    }
    PostData(1, fOutputMB);
    PostData(2, fOutputEMC7);
    PostData(3, fOutputEMCEGA);
    PostData(4, fOutputEMCJet);
    return;
}

//________________________________________________________________________
void AliAnalysisTaskPSHFE::Terminate(Option_t *) 
{

}
