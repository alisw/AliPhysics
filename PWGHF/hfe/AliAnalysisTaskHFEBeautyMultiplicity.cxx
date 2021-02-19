/* ************************************************************************
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
 ************************************************************************ */

/*===============================================*
 *                                               *
 * AliAnalysisTaskHFEBeautyMultiplicity          *
 * Auther : Shunya Chiba, University of Tsukuba  *
 *                                               *
 *===============================================*/

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TRandom.h"
#include "TMath.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskHFEBeautyMultiplicity.h"
#include "AliVVertex.h"
#include "AliVCluster.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerPatchFinder.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliVertexingHFUtils.h"
#include "AliKFParticle.h"

//---- Header for Monte Carlo
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"



class AliAnalysisTaskHFEBeautyMultiplicity;    // analysis class

using namespace std;

ClassImp(AliAnalysisTaskHFEBeautyMultiplicity) // classimp: necessary for root

AliAnalysisTaskHFEBeautyMultiplicity::AliAnalysisTaskHFEBeautyMultiplicity() : AliAnalysisTaskSE(),
    fAOD(0),
    fOutputList(0),
    fpidResponse(0),
    fVevent(0),
    fMultSelection(0),
    fCentralityEstimator("V0M"),

    //---- Tender ----//
    fTracks_tender(0),
    fUseTender(kTRUE),
    fCaloClusters_tender(0),

    fNevents(0),
    fCent(0),               // Centrality
    fMult(0),               // Multiplicity
    fVtxZ(0),               // Zvertex
    fVtxZ_2(0),             // Zvertex after exent cut
    fVtxX(0),               // Xvertex
    fVtxY(0),               // Yvertex

    fVtxCorrelation(0),     // Primary Zvertex vs. SPD Zvertex
    fNcont(0),		    // NcontV vs. NcontVSPD

    fZvtx_Ntrklet(0),	    // Z vertex vs N tracklets
    fZvtx_Ntrklet_Corr(0),  // Z vertex vs N tracklets (Corrected)

    fEMCClsEtaPhi(0),       // EMCal Cluster Eta vs. Phi
    fHistNCells(0),         // No. of EMCal Cells in a cluster
    fHistCalCells(0),       // EMCal cells in a cluster
    fHistClustE(0),         // EMCal cluster energy distributioon
    fHistNClsE(0),          // No. of EMCal clusters in the event (ALL)
    fHistNClsE1(0),         // No. of EMCal clusters in the event (ClsE > 0.1 GeV)
    fHistNClsE2(0),         // No. of EMCal clusters in the event (ClsE > 0.2 GeV)
    fHistNClsE3(0),         // No. of EMCal clusters in the event (ClsE > 0.5 GeV)
    fAllTrkPt(0),           // pT distribution (ALL)
    fEMCTrkPt(0),           // pT distribution (matched EMCal)
    fAllTrkEta(0),          // eta distribution (All)
    fEMCTrkEta(0),          // eta distribution (matched EMCal)
    fAllTrkPhi(0),          // phi distribution (All)
    fEMCTrkPhi(0),          // phi distribution (matched EMCal)
    fPhiEta(0),		    // Trk Phi vs. Eta
    fTPCCrossedRow(0),      // TPC CrossedRows

    fdEdx(0),               // All track dE/dx distribution
    fTPCnsig(0),            // ALL track TPC Nsigma distribution (electron)
    fTPCnsig_Pi(0),         // All track TPC Nsigma distribution (pion)
    fTOFnsig(0),            // All track TOF Nsigma distribution
    fITSnsig(0),            // All track ITS Nsigma distribution
    fTPCnsigEta0(0),        // TPC Nsigma vs. Eta (pT > 2 GeV/c)
    fTPCnsigEta1(0),        // TPC Nsigma vs. Eta (pT > 3 GeV/c)
    fTPCnsigEta2(0),        // TPC Nsigma vs. Eta (pT > 5 GeV/c)
    fClsEtaPhiAftMatch(0),          //EMCal cluster Eta vs. Phi (after track matching)
    fClsEtaPhiAftMatchEMCin(0),     //EMCal cluster Eta vs. Phi (after track matching inside EMCal)
    fClsEtaPhiAftMatchEMCout(0),    //EMCal cluster Eta vs. Phi (after track matching outside EMCal)
    fHistEMCTrkMatch_Eta(0),        // Distance of EMCal cluster to its closest track (#Delta#eta)
    fHistEMCTrkMatch_Phi(0),        // Distance of EMCal cluster to its closest track (#Delta#phi)
    fEMCTrkMatch_EtaPhi(0),         // deltaEta vs deltaPhi
    fEMCTrkMatch_EtaPhi_AfterCut(0),// deltaEta vs deltaPhi after cut

    fTrkPt_2(0),	    // track pT (after track cut)
    fTrkEta_2(0),	    // track Eta (after track cut)
    fTrkPhi_2(0),	    // track Phi (after track cut)
    fdEdx_2(0),		    // dE/dx (after track cut)
    fTPCnsig_2(0),	    // TPC Nsigma (after track cut)
    fTOFnsig_2(0),	    // TOF Nsigma (after track cut)
    fITSnsig_2(0),	    // ITS Nsigma (after track cut)
    fTPCCrossedRow_2(0),    // TPC CrossedRows (after track cut)

    fM02_1(0),              // long axis
    fM20_1(0),              // short axis
    fM02_2(0),              // long axis (after PID)
    fM20_2(0),              // short axis (after PID)
    fNtracks(0),            // track selection
    fTrkEtaPhi_AfterCut(0), // Track Eta vs. Phi (after cut)_

    fHistEopAll(0),         // E/p Histgram (ALL)
    fEopElectron1(0),       // pT vs electron E/p
    fEopHadron1(0),         // pT vs Hadron E/p

    fInvmassLS(0),          // Invariant mass vs. pT (Like-sign)
    fInvmassULS(0),         // Invariant mass vs. pT (Unlike-sign)

    fDCAxy_Ele_1(0),        // DCA Electron (DCA*charge*Bsign)
    fDCAxy_Had_1(0),        // DCA Hadron (DCA*charge*Bsign)
    fDCAxy_LS_1(0),         // DCA Like-sign (DCA*charge*Bsign)
    fDCAxy_ULS_1(0),        // DCA Unlike-sign (DCA*charge*Bsign)

    fDCAxy_Ele_2(0),        // DCA Electron (DCA*charge)
    fDCAxy_Had_2(0),        // DCA Hadron (DCA*charge)
    fDCAxy_LS_2(0),         // DCA Like-sign (DCA*charge)
    fDCAxy_ULS_2(0),        // DCA Unlike-sign (DCA*charge)

    fDCAxy_Ele_3(0),        // DCA Electron (DCA)
    fDCAxy_Had_3(0),        // DCA Hadron (DCA)
    fDCAxy_LS_3(0),         // DCA Like-sign (DCA)
    fDCAxy_ULS_3(0),        // DCA Unlike-sign (DCA)

    fDCAxy_Ele_4(0),

    fEopElectron2(0),       // electron pT (tight)
    fEopHadron2(0),         // hadron pT (tight)
    fEopElectron3(0),       // electron except photonic(invariant mass)

    fHistConv_R(0),	    // conversion R
    fElectronEtaPhi(0),	    // eta vs. phi (electron)
    fHadronEtaPhi(0),	    // eta vs. phi (hadron)


    fMultiEstimatorAvg(0),  // TProfile
    Nref(0),


    //---- MC data ----//
    fMCcheckMother(0),
    fMCarray(0),
    fMCparticle(0),
    fMCTrackpart(0),
    fMCheader(0),

    fHistPho_Reco0(0),
    fHistPho_Reco1(0),
    fHistPho_Reco2(0),
    NembMCpi0(0),
    NembMCeta(0),
    NpureMCproc(0),
    NpureMC(0),
    Nch(0),
    iBevt(kFALSE),
    fNoB(0),
    fNoD(0),

    fCheckEtaMC(0),
    fHistMCorg_Pi0(0),
    fHistMCorg_Eta(0),
    fHistMCorg_D(0),
    fHistMCorg_BD(0),
    fHistMCorg_B(0),
    fHistMCorg_Lc(0),
    fPt_Btoe(0),
    fHistPt_HFE_MC_B(0),
    fHistPt_HFE_MC_D(0),
    fHistPt_HFE_MC_Lc(0),

    fDCAxy_MC_B(0),	// DCA from B
    fDCAxy_MC_D(0),	// DCA from D
    fDCAxy_MC_Dpm(0),	// DCA from D+,D*+
    fDCAxy_MC_D0(0),	// DCA from D0,D*0
    fDCAxy_MC_Ds(0),	// DCA from Ds+,D*+s
    fDCAxy_MC_Lc(0),	// DCA from Lambda

    fDCAxy_MC_ele(0),
    fDCAxy_MC_Phot(0),

    fHistPt_B_TrkCut0(0),
    fHistPt_B_TrkCut1(0),
    fHistPt_B_TrkCut2(0),
    fHistPt_B_TrkCut3(0),
    fHistPt_B_TrkCut4(0),
    fHistPt_B_TrkCut5(0),
    fHistPt_B_TrkCut6(0),
    fHistPt_B_TrkCut7(0),
    fHistPt_B_TrkCut8(0),
    fHistPt_B_TrkCut9(0),

    fHistPt_D_TrkCut0(0),
    fHistPt_D_TrkCut1(0),
    fHistPt_D_TrkCut2(0),
    fHistPt_D_TrkCut3(0),
    fHistPt_D_TrkCut4(0),
    fHistPt_D_TrkCut5(0),
    fHistPt_D_TrkCut6(0),
    fHistPt_D_TrkCut7(0),
    fHistPt_D_TrkCut8(0),
    fHistPt_D_TrkCut9(0)







{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskHFEBeautyMultiplicity::AliAnalysisTaskHFEBeautyMultiplicity(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0),
    fOutputList(0),
    fpidResponse(0),
    fVevent(0),
    fMultSelection(0),
    fCentralityEstimator("V0M"),

    //---- Tender ----//
    fTracks_tender(0),
    fUseTender(kTRUE),
    fCaloClusters_tender(0),

    fNevents(0),
    fCent(0),               // Centrality
    fMult(0),               // Multiplicity
    fVtxZ(0),               // Zvertex
    fVtxZ_2(0),             // Zvertex after exent cut
    fVtxX(0),               // Xvertex
    fVtxY(0),               // Yvertex

    fVtxCorrelation(0),     // Primary Zvertex vs. SPD Zvertex
    fNcont(0),		    // NcontV vs. NcontVSPD
    
    fZvtx_Ntrklet(0),	    // Z vertex vs N tracklets
    fZvtx_Ntrklet_Corr(0),  // Z vertex vs N tracklets (Corrected)

    fEMCClsEtaPhi(0),       // EMCal Cluster Eta vs. Phi
    fHistNCells(0),         // No. of EMCal Cells in a cluster
    fHistCalCells(0),       // EMCal cells in a cluster
    fHistClustE(0),         // EMCal cluster energy distributioon
    fHistNClsE(0),          // No. of EMCal clusters in the event (ALL)
    fHistNClsE1(0),         // No. of EMCal clusters in the event (ClsE > 0.1 GeV)
    fHistNClsE2(0),         // No. of EMCal clusters in the event (ClsE > 0.2 GeV)
    fHistNClsE3(0),         // No. of EMCal clusters in the event (ClsE > 0.5 GeV)
    fAllTrkPt(0),           // pT distribution (ALL)
    fEMCTrkPt(0),           // pT distribution (matched EMCal)
    fAllTrkEta(0),          // eta distribution (All)
    fEMCTrkEta(0),          // eta distribution (matched EMCal)
    fAllTrkPhi(0),          // phi distribution (All)
    fEMCTrkPhi(0),          // phi distribution (matched EMCal)
    fPhiEta(0),		    // Trk Phi vs. Eta
    fTPCCrossedRow(0),      // TPC CrossedRows

    fdEdx(0),               // All track dE/dx distribution
    fTPCnsig(0),            // ALL track TPC Nsigma distribution (electron)
    fTPCnsig_Pi(0),         // All track TPC Nsigma distribution (pion)
    fTOFnsig(0),            // All track TOF Nsigma distribution
    fITSnsig(0),            // All track ITS Nsigma distribution
    fTPCnsigEta0(0),        // TPC Nsigma vs. Eta (pT > 2 GeV/c)
    fTPCnsigEta1(0),        // TPC Nsigma vs. Eta (pT > 3 GeV/c)
    fTPCnsigEta2(0),        // TPC Nsigma vs. Eta (pT > 5 GeV/c)
    fClsEtaPhiAftMatch(0),          //EMCal cluster Eta vs. Phi (after track matching)
    fClsEtaPhiAftMatchEMCin(0),     //EMCal cluster Eta vs. Phi (after track matching inside EMCal)
    fClsEtaPhiAftMatchEMCout(0),    //EMCal cluster Eta vs. Phi (after track matching outside EMCal)
    fHistEMCTrkMatch_Eta(0),        // Distance of EMCal cluster to its closest track (#Delta#eta)
    fHistEMCTrkMatch_Phi(0),        // Distance of EMCal cluster to its closest track (#Delta#phi)
    fEMCTrkMatch_EtaPhi(0),         // deltaEta vs deltaPhi
    fEMCTrkMatch_EtaPhi_AfterCut(0),// deltaEta vs deltaPhi after cut

    fTrkPt_2(0),	    // track pT (after track cut)
    fTrkEta_2(0),	    // track Eta (after track cut)
    fTrkPhi_2(0),	    // track Phi (after track cut)
    fdEdx_2(0),		    // dE/dx (after track cut)
    fTPCnsig_2(0),	    // TPC Nsigma (after track cut)
    fTOFnsig_2(0),	    // TOF Nsigma (after track cut)
    fITSnsig_2(0),	    // ITS Nsigma (after track cut)
    fTPCCrossedRow_2(0),    // TPC CrossedRows (after track cut)

    fM02_1(0),              // long axis
    fM20_1(0),              // short axis
    fM02_2(0),              // long axis (after PID)
    fM20_2(0),              // short axis (after PID)
    fNtracks(0),            // track selection
    fTrkEtaPhi_AfterCut(0), // Track Eta vs. Phi (after cut)_

    fHistEopAll(0),         // E/p Histgram (ALL)
    fEopElectron1(0),       // pT vs electron E/p
    fEopHadron1(0),         // pT vs Hadron E/p

    fInvmassLS(0),          // Invariant mass vs. pT (Like-sign)
    fInvmassULS(0),         // Invariant mass vs. pT (Unlike-sign)

    fDCAxy_Ele_1(0),        // DCA Electron (DCA*charge*Bsign)
    fDCAxy_Had_1(0),        // DCA Hadron (DCA*charge*Bsign)
    fDCAxy_LS_1(0),         // DCA Like-sign (DCA*charge*Bsign)
    fDCAxy_ULS_1(0),        // DCA Unlike-sign (DCA*charge*Bsign)

    fDCAxy_Ele_2(0),        // DCA Electron (DCA*charge)
    fDCAxy_Had_2(0),        // DCA Hadron (DCA*charge)
    fDCAxy_LS_2(0),         // DCA Like-sign (DCA*charge)
    fDCAxy_ULS_2(0),        // DCA Unlike-sign (DCA*charge)

    fDCAxy_Ele_3(0),        // DCA Electron (DCA)
    fDCAxy_Had_3(0),        // DCA Hadron (DCA)
    fDCAxy_LS_3(0),         // DCA Like-sign (DCA)
    fDCAxy_ULS_3(0),        // DCA Unlike-sign (DCA)

    fDCAxy_Ele_4(0),

    fEopElectron2(0),       // electron pT (tight)
    fEopHadron2(0),         // hadron pT (tight)
    fEopElectron3(0),       // electron except photonic(invariant mass)

    fHistConv_R(0),	    // conversion R
    fElectronEtaPhi(0),	    // eta vs. phi (electron)
    fHadronEtaPhi(0),	    // eta vs. phi (hadron)

    fMultiEstimatorAvg(0),
    Nref(0),



    //---- MC data ----//
    fMCcheckMother(0),
    fMCarray(0),
    fMCparticle(0),
    fMCTrackpart(0),
    fMCheader(0),

    fHistPho_Reco0(0),
    fHistPho_Reco1(0),
    fHistPho_Reco2(0),
    NembMCpi0(0),
    NembMCeta(0),
    NpureMCproc(0),
    NpureMC(0),
    Nch(0),
    iBevt(kFALSE),
    fNoB(0),
    fNoD(0),

    fCheckEtaMC(0),
    fHistMCorg_Pi0(0),
    fHistMCorg_Eta(0),
    fHistMCorg_D(0),
    fHistMCorg_BD(0),
    fHistMCorg_B(0),
    fHistMCorg_Lc(0),
    fPt_Btoe(0),
    fHistPt_HFE_MC_B(0),
    fHistPt_HFE_MC_D(0),
    fHistPt_HFE_MC_Lc(0),

    fDCAxy_MC_B(0),	// DCA from B
    fDCAxy_MC_D(0),	// DCA from D
    fDCAxy_MC_Dpm(0),	// DCA from D+,D*+
    fDCAxy_MC_D0(0),	// DCA from D0,D*0
    fDCAxy_MC_Ds(0),	// DCA from Ds+,D*+s
    fDCAxy_MC_Lc(0),	// DCA from Lambda

    fDCAxy_MC_ele(0),
    fDCAxy_MC_Phot(0),    

    fHistPt_B_TrkCut0(0),
    fHistPt_B_TrkCut1(0),
    fHistPt_B_TrkCut2(0),
    fHistPt_B_TrkCut3(0),
    fHistPt_B_TrkCut4(0),
    fHistPt_B_TrkCut5(0),
    fHistPt_B_TrkCut6(0),
    fHistPt_B_TrkCut7(0),
    fHistPt_B_TrkCut8(0),
    fHistPt_B_TrkCut9(0),

    fHistPt_D_TrkCut0(0),
    fHistPt_D_TrkCut1(0),
    fHistPt_D_TrkCut2(0),
    fHistPt_D_TrkCut3(0),
    fHistPt_D_TrkCut4(0),
    fHistPt_D_TrkCut5(0),
    fHistPt_D_TrkCut6(0),
    fHistPt_D_TrkCut7(0),
    fHistPt_D_TrkCut8(0),
    fHistPt_D_TrkCut9(0)



{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskHFEBeautyMultiplicity::~AliAnalysisTaskHFEBeautyMultiplicity()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
        delete fTracks_tender;
	}
}


//_____________________________________________________________________________
void AliAnalysisTaskHFEBeautyMultiplicity::UserCreateOutputObjects()
{


//**********************************************//
// Automatic determination of the analysis mode //
//**********************************************//
	AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
		SetAODAnalysis();
	} else {
		SetESDAnalysis();
	}


    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)


//********************//
// create a histogram //
//********************//

  //Nevents
    fNevents = new TH1F("fNevents","Number of events",7, -0.5, 6.5);
    fOutputList->Add(fNevents);
    fNevents->GetYaxis()->SetTitle("counts");
    fNevents->GetXaxis()->SetBinLabel(1,"ALL events");
    fNevents->GetXaxis()->SetBinLabel(2,"pile-up cut");
    fNevents->GetXaxis()->SetBinLabel(3,"Global track cont cut");
    fNevents->GetXaxis()->SetBinLabel(4,"SPD track cont cut");
    fNevents->GetXaxis()->SetBinLabel(5,"SPD&Global mutch");
    fNevents->GetXaxis()->SetBinLabel(6,"SPD resolusion cut");
    fNevents->GetXaxis()->SetBinLabel(7,"|Zvertex| < 10cm");
    
  //Primary Zvertex vs. SPD Zvertex
    fVtxCorrelation = new TH2F("fVtxCorrelation",";Z_{vertex}^{Primary} (cm);Z_{vertex}^{SPD} (cm)",1200,-30,30,1200,-30,30);
    fOutputList->Add(fVtxCorrelation);

  //Number of contribution
    fNcont = new TH2F("fNcont","Number of contribution;Ncont (primary);Ncont (SPD)",500,0,500,500,0,500);
    fOutputList->Add(fNcont);
    
  //centrality
    fCent = new TH1F("fCent","Centrality;centrality(%);counts",100,0,100);
    fOutputList->Add(fCent);

  //Multiplicity
    fMult = new TH2F("fMult","Track multiplicity;centrality(%);",100,0,100,20000,0,20000);
    fOutputList->Add(fMult);

  //Zvertex  
    fVtxZ = new TH1F("fVtxZ","Z vertex position; Z_{vertex} [cm]; counts", 500, -25, 25);
    fOutputList->Add(fVtxZ);

  //Zvertex after event cut  
    fVtxZ_2 = new TH1F("fVtxZ_2","Z vertex position (after event cut); Z_{vertex} [cm]; counts", 500, -25, 25);
    fOutputList->Add(fVtxZ_2);

  //Xvertex  
    fVtxX = new TH1F("fVtxX","X vertex position; Vtx_{x} [cm]; counts", 500, -25, 25);
    fOutputList->Add(fVtxX);
    
  //Yvertex  
    fVtxY = new TH1F("fVtxY","Y vertex position; Vtx_{y} [cm]; counts", 500, -25, 25);
    fOutputList->Add(fVtxY);

  //Z vertex vs N tracklets
    fZvtx_Ntrklet = new TH2F("fZvtx_Ntrklet","Zvertex vs N tracklets;Z_{vertex} [cm];N^{SPD}_{tracklets}",400,-20,20,4001,-0.5,4000.5);
    fOutputList->Add(fZvtx_Ntrklet);

  //Z vertex vs N tracklets (Corrected)
    fZvtx_Ntrklet_Corr = new TH2F("fZvtx_Ntrklet_Corr","Zvertex vs N tracklets (Corrected);Z_{vertex} [cm];N^{SPD}_{tracklets}",400,-20,20,4001,-0.5,4000.5);
    fOutputList->Add(fZvtx_Ntrklet_Corr);

  //EMCal Cluster Eta and Phi
    fEMCClsEtaPhi = new TH2F("fEMCClsEtaPhi","EMCal&DCal Cluster #eta and #phi distribution; #eta; #phi",150,-0.75,0.75,63,0,6.3);
    fOutputList->Add(fEMCClsEtaPhi);

  //No of EMCal Cells in a cluster
    fHistNCells = new TH2F("fHistNCells","No. of EMCal cells in a cluster; Cluster E (GeV); N^{EMC}_{Cells}",500,0,50,30,0,30);
    fOutputList->Add(fHistNCells);

  //EMCal cells in a cluster
    fHistCalCells = new TH2F("fHistCalCells","EMCal cells in a cluster;cell ID;E (GeV)",15000,-0.5,14999.5,300,0,30);
    fOutputList->Add(fHistCalCells);

  //EMCal cluster energy distributioon
    fHistClustE = new TH1F("fHistClustE","EMCal cluster energy distribution; Cluster E (GeV); counts",500,0,50);
    fOutputList->Add(fHistClustE);

  //NCls ALL
    fHistNClsE = new TH1F("fHistNClsE","No. of EMCal cluster in the event; N^{EMC}_{cls}; counts",2000,0,2000);
    fOutputList->Add(fHistNClsE);

  //NCls  E > 0.1 GeV
    fHistNClsE1 = new TH1F("fHistNClsE1","No. of EMCal cluster in the event (E > 0,1 GeV); N^{EMC}_{cls}; counts",2000,0,2000);
    fOutputList->Add(fHistNClsE1);

  //NCls  E > 0.2 GeV
    fHistNClsE2 = new TH1F("fHistNClsE2","No. of EMCal cluster in the event (E > 0.2 GeV); N^{EMC}_{cls}; counts",2000,0,2000);
    fOutputList->Add(fHistNClsE2);

  //NCls  E > 0.5 GeV
    fHistNClsE3 = new TH1F("fHistNClsE3","No. of EMCal cluster in the event (E > 0.5 GeV); N^{EMC}_{cls}; counts",2000,0,2000);
    fOutputList->Add(fHistNClsE3);

  //pT distribution All
    fAllTrkPt = new TH1F("fAllTrkPt","p_{T} distribution ; p_{T} (GeV/c); counts",1000,0,50.0);
    fOutputList->Add(fAllTrkPt);

  //pT distribution matched to EMCal
    fEMCTrkPt = new TH1F("fEMCTrkPt","p_{T} distribution of tracks matched to EMCal; p_{T} (GeV/c); counts",1000,0,50.0);
    fOutputList->Add(fEMCTrkPt);

  //Eta distribution (All)
    fAllTrkEta = new TH1F("fAllTrkEta","Track #eta distribution; #eta; counts",80,-2.0,2.0);
    fOutputList->Add(fAllTrkEta);

  //Eta distribution (EMCal)
    fEMCTrkEta = new TH1F("fEMCTrkEta","#eta distribution of  tracks matched to EMCal; #eta; counts",80,-2.0,2.0);
    fOutputList->Add(fEMCTrkEta);

  //Phi distribution (All)
    fAllTrkPhi = new TH1F("fAllTrkPhi","Track #phi distribution; #phi; counts",63,0,6.3);
    fOutputList->Add(fAllTrkPhi);

  //Phi distribution (EMCal)
    fEMCTrkPhi = new TH1F("fEMCTrkPhi","#phi distribution of  tracks matched to EMCal; #phi; counts",63,0,6.3);
    fOutputList->Add(fEMCTrkPhi);

  //TPC CrossedRows
    fTPCCrossedRow = new TH1F("fTPCCrossedRow","No of TPC CrossedRow; N^{ITS}_{CrossedRow}; counts",200,0.,200.); 
    fOutputList->Add(fTPCCrossedRow);

  //Phi vs. Eta
    fPhiEta = new TH2F("fPhiEta","#phi vs. #eta",630,0,6.3,400,-2.0,2.0);
    fOutputList->Add(fPhiEta);


  //dE/dx distribution (electron)
    fdEdx = new TH2F("fdEdx","All track dE/dx distribution;p (GeV/c);dE/dx",300,0,15,1600,0,160);
    fOutputList->Add(fdEdx);

  //TPC Nsigma distribution (electron)
    fTPCnsig = new TH2F("fTPCnsig","All track TPC Nsigma distribution (electron);p (GeV/c);n^{TPC}_{#sigma_{electron}}",300,0,15,200,-10,10);
    fOutputList->Add(fTPCnsig);

  //TPC Nsigma distribution (pion)
    fTPCnsig_Pi = new TH2F("fTPCnsig_Pi","All track TPC Nsigma distribution (pion);p (GeV/c);n^{TPC}_{#sigma_{pion}}",300,0,15,200,-10,10);
    fOutputList->Add(fTPCnsig_Pi);
    
  //TOF Nsigma distribution (electron)
    fTOFnsig = new TH2F("fTOFnsig","All track TOF Nsigma distribution (electron);p (GeV/c);n^{TOF}_{#sigma_{electron}}",300,0,15,200,-10,10);
    fOutputList->Add(fTOFnsig);

  //ITS Nsigma distribution (electron)
    fITSnsig = new TH2F("fITSnsig","All track ITS Nsigma distribution (electron);p (GeV/c);n^{ITS}_{#sigma_{electron}}",300,0,15,200,-10,10);
    fOutputList->Add(fITSnsig);

  //TPC Nsigma vs. Eta (pT > 2 GeV/c)
    fTPCnsigEta0 = new TH2F("fTPCnsigEta0","TPC Nsigma (electron) vs. Eta;#eta;n^{TPC}_{#sigma_{electron}}",40,-1,1,200,-10,10);
    fOutputList->Add(fTPCnsigEta0);

  //TPC Nsigma vs. Eta (pT > 3 GeV/c)
    fTPCnsigEta1 = new TH2F("fTPCnsigEta1","TPC Nsigma (electron) vs. Eta;#eta;n^{TPC}_{#sigma_{electron}}",40,-1,1,200,-10,10);
    fOutputList->Add(fTPCnsigEta1);

  //TPC Nsigma vs. Eta (pT > 5 GeV/c)
    fTPCnsigEta2 = new TH2F("fTPCnsigEta2","TPC Nsigma (electron) vs. Eta;#eta;n^{TPC}_{#sigma_{electron}}",40,-1,1,200,-10,10);
    fOutputList->Add(fTPCnsigEta2);

  //EMCal cluster Eta and Phi (after track matching)
    fClsEtaPhiAftMatch = new TH2F("fClsEtaPhiAftMatch","EMCal cluster #eta and #phi distribution after track matching;#eta;#phi",180,-0.9,0.9,630,0,6.3);
    fOutputList->Add(fClsEtaPhiAftMatch);

  //EMCal cluster Eta and Phi (after track matching inside EMCal)
    fClsEtaPhiAftMatchEMCin = new TH2F("fClsEtaPhiAftMatchEMCin","EMCal cluster #eta and #phi distribution after track matching (inside EMCal #phi acceptance);#eta;#phi",160,-0.8,0.8,630,0,6.3);
    fOutputList->Add(fClsEtaPhiAftMatchEMCin);

  //EMCal cluster Eta and Phi (after track matching outside EMCal)
    fClsEtaPhiAftMatchEMCout = new TH2F("fClsEtaPhiAftMatchEMCout","EMCal cluster #eta and #phi distribution after track matching (outside EMCal #phi acceptance);#eta;#phi",160,-0.8,0.8,630,0,6.3);
    fOutputList->Add(fClsEtaPhiAftMatchEMCout);

  //Distance of EMCal cluster to its closest track (#Delta#eta)
    fHistEMCTrkMatch_Eta = new TH1F("fHistEMCTrkMatch_Eta","Distance of EMCal cluster to its closest track (#Delta#eta);#Delta#eta;counts",600,-0.3,0.3);
    fOutputList->Add(fHistEMCTrkMatch_Eta);

  //Distance of EMCal cluster to its closest track (#Delta#phi)
    fHistEMCTrkMatch_Phi = new TH1F("fHistEMCTrkMatch_Phi","Distance of EMCal cluster to its closest track (#Delta#phi);#Delta#phi;counts",600,-0.3,0.3);
    fOutputList->Add(fHistEMCTrkMatch_Phi);

  //Distance of EMCal cluster (#Eta vs #Phi)
    fEMCTrkMatch_EtaPhi = new TH2F("fEMCTrkMatch_EtaPhi","Distance of EMCal Cluster (#Delta#eta vs #Delta#phi);#Delta#eta;#Delta#phi (rad)",600,-0.3,0.3,600,-0.3,0.3);
    fOutputList->Add(fEMCTrkMatch_EtaPhi);
  
  //Distance of EMCal cluster (#Eta vs #Phi) after cut
    fEMCTrkMatch_EtaPhi_AfterCut = new TH2F("fEMCTrkMatch_EtaPhi_AfterCut","Distance of EMCal Cluster (#Delta#eta vs #Delta#phi) after cut;#Delta#eta;#Delta#phi (rad)",200,-0.1,0.1,200,-0.1,0.1);
    fOutputList->Add(fEMCTrkMatch_EtaPhi_AfterCut);

  //pT distribution (after track cut)
    fTrkPt_2 = new TH1F("fTrkPt_2","p_{T} distribution (after track cut) ; p_{T} (GeV/c); counts",1000,0,50.0);
    fOutputList->Add(fTrkPt_2);

  //Eta distribution (after track cut)
    fTrkEta_2 = new TH1F("fTrkEta_2","Track #eta distribution (after track cut); #eta; counts",80,-2.0,2.0);
    fOutputList->Add(fTrkEta_2);

  //Phi distribution (after track cut)
    fTrkPhi_2 = new TH1F("fTrkPhi_2","Track #phi distribution (after track cut); #phi (rad); counts",63,0,6.3);
    fOutputList->Add(fTrkPhi_2);

  //dE/dx (after track cut)
    fdEdx_2 = new TH2F("fdEdx_2","dE/dx distribution (after track cut);p (GeV/c);dE/dx",300,0,15,1600,0,160);
    fOutputList->Add(fdEdx_2);

  //TPC Nsigma distribution (after track cut)
    fTPCnsig_2 = new TH2F("fTPCnsig_2","TPC Nsigma distribution (after track cut);p (GeV/c);n^{TPC}_{#sigma_{electron}}",300,0,15,200,-10,10);
    fOutputList->Add(fTPCnsig_2);
    
  //TOF Nsigma distribution (after track cut)
    fTOFnsig_2 = new TH2F("fTOFnsig_2","All track TOF Nsigma distribution (after track cut);p (GeV/c);n^{TOF}_{#sigma_{electron}}",300,0,15,200,-10,10);
    fOutputList->Add(fTOFnsig_2);

  //ITS Nsigma distribution (after track cut)
    fITSnsig_2 = new TH2F("fITSnsig_2","ITS Nsigma distribution (after track cut);p (GeV/c);n^{ITS}_{#sigma_{electron}}",300,0,15,200,-10,10);
    fOutputList->Add(fITSnsig_2);

  //TPC CrossedRows (after track cut)
    fTPCCrossedRow_2 = new TH1F("fTPCCrossedRow_2","No of TPC CrossedRow; N^{ITS}_{CrossedRow}; counts",200,0.,200.); 
    fOutputList->Add(fTPCCrossedRow_2);
    
  //M02
    fM02_1 = new TH2F("fM02_1","M02 vs p_{T} distribution;p_{T} [GeV/c];long axis of ellipse : M02 [cm]",400,0,20,800,0,4);
    fOutputList -> Add(fM02_1);

  //M20
    fM20_1 = new TH2F("fM20_1","M20 vs p_{T} distribution;p_{T} [GeV/c];short axis of ellipse : M20 [cm]",400,0,20,800,0,4);
    fOutputList -> Add(fM20_1);
    
  //M02 (after PID)
    fM02_2 = new TH2F("fM02_2","M02 vs p_{T} distribution;p_{T} [GeV/c];long axis of ellipse (after PID) : M02 [cm]",400,0,20,800,0,4);
    fOutputList -> Add(fM02_2);

  //M20 (after PID)
    fM20_2 = new TH2F("fM20_2","M20 vs p_{T} distribution;p_{T} [GeV/c];short axis of ellipse (after PID) : M20 [cm]",400,0,20,800,0,4);
    fOutputList -> Add(fM20_2);

  //Track Eta vs. Phi (after Track cut)
    fTrkEtaPhi_AfterCut = new TH2F("fTrkEtaPhi_AfterCut","Track Eta vs. Phi (after cut);#eta,#phi (rad)",180,-0.9,0.9,630,0,6.3);
    fOutputList->Add(fTrkEtaPhi_AfterCut);

  //E/p (all)
    fHistEopAll = new TH1F("fHistEopAll","E/p;E/p;counts",60,0,3.0);
    fOutputList->Add(fHistEopAll);

  //Ntracks
    fNtracks = new TH1F("fNtracks","Number of tracks",11, -0.5, 10.5);
    fOutputList->Add(fNtracks);
    fNtracks->GetYaxis()->SetTitle("counts");
    fNtracks->GetXaxis()->SetBinLabel(1,"matching tracks");
    fNtracks->GetXaxis()->SetBinLabel(2,"AOD standard");
    fNtracks->GetXaxis()->SetBinLabel(3,"TPC and ITS refit");
    fNtracks->GetXaxis()->SetBinLabel(4,"TPCCrossedRow cut");
    //fNtracks->GetXaxis()->SetBinLabel(4,"TPC cluster cut");
    fNtracks->GetXaxis()->SetBinLabel(5,"ITS cluster cut");
    fNtracks->GetXaxis()->SetBinLabel(6,"dE/dx calculation");
    fNtracks->GetXaxis()->SetBinLabel(7,"SPD hit cut");
    fNtracks->GetXaxis()->SetBinLabel(8,"chi2 cut");
    fNtracks->GetXaxis()->SetBinLabel(9,"DCA cut");
    fNtracks->GetXaxis()->SetBinLabel(10,"eta & phi diff cut");
    fNtracks->GetXaxis()->SetBinLabel(11,"Eta cut");
    
  //pT vs E/p (electron)
    fEopElectron1 = new TH2F("fEopElectron1","Electron;p_{T} [GeV/c];E/p",600,0,30,150,0,3.0);
    fOutputList->Add(fEopElectron1);
    
  //pT vs E/p (hadron)
    fEopHadron1 = new TH2F("fEopHadron1","Hadron;p_{T} [GeV/c];E/p",600,0,30,150,0,3.0);
    fOutputList->Add(fEopHadron1);

  //Invariant mass vs. pT (like-sign)
    fInvmassLS = new TH2F("fInvmassLS", "Invariant mass vs. p_{T} (LS);p_{T} [GeV/c];mass [GeV/c^{2}]",600,0,30,200,0,0.5);
    fOutputList->Add(fInvmassLS);
    
  //Invariant mass vs. pT (unlike-sign)
    fInvmassULS = new TH2F("fInvmassULS", "Invariant mass vs. p_{T} (ULS);p_{T} [GeV/c];mass [GeV/c^{2}]",600,0,30,200,0,0.5);
    fOutputList->Add(fInvmassULS);
  

  //---- DCA * charge * Bsign ----//
  //pT vs. DCA (electron)
    fDCAxy_Ele_1 = new TH2F("fDCAxy_Ele_1","p_{T} vs. DCA_{xy} (Electron);p_{T} [GeV/c];DCA_{xy} #times charge #times Bsign [cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_Ele_1);
    
  //pT vs. DCA (Hadron)
    fDCAxy_Had_1 = new TH2F("fDCAxy_Had_1","p_{T} vs. DCA_{xy} (Hadron);p_{T} [GeV/c];DCA_{xy} #times charge #times Bsign [cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_Had_1);
    
  //pT vs. DCA (Like-sign)
    fDCAxy_LS_1 = new TH2F("fDCAxy_LS_1","p_{T} vs. DCA_{xy} (Like-sign pairs);p_{T} [GeV/c];DCA_{xy} #times charge #times Bsign [cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_LS_1);
      
  //pT vs. DCA (Unlike-sign)
    fDCAxy_ULS_1 = new TH2F("fDCAxy_ULS_1","p_{T} vs. DCA_{xy} (Unlike-sign pairs);p_{T} [GeV/c];DCA_{xy} #times charge #times Bsign [cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_ULS_1);
  

  //---- DCA * charge ----//
  //pT vs. DCA (electron)
    fDCAxy_Ele_2 = new TH2F("fDCAxy_Ele_2","p_{T} vs. DCA_{xy} (Electron);p_{T} [GeV/c];DCA_{xy} #times charge [cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_Ele_2);
    
  //pT vs. DCA (Hadron)
    fDCAxy_Had_2 = new TH2F("fDCAxy_Had_2","p_{T} vs. DCA_{xy} (Hadron);p_{T} [GeV/c];DCA_{xy} #times charge [cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_Had_2);
    
  //pT vs. DCA (Like-sign)
    fDCAxy_LS_2 = new TH2F("fDCAxy_LS_2","p_{T} vs. DCA_{xy} (Like-sign pairs);p_{T} [GeV/c];DCA_{xy} #times charge [cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_LS_2);
      
  //pT vs. DCA (Unlike-sign)
    fDCAxy_ULS_2 = new TH2F("fDCAxy_ULS_2","p_{T} vs. DCA_{xy} (Unlike-sign pairs);p_{T} [GeV/c];DCA_{xy} #times charge [cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_ULS_2);
    
  
  //---- DCA ----//
  //pT vs. DCA (electron)
    fDCAxy_Ele_3 = new TH2F("fDCAxy_Ele_3","p_{T} vs. DCA_{xy} (Electron);p_{T} [GeV/c];DCA_{xy} [cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_Ele_3);
    
  //pT vs. DCA (Hadron)
    fDCAxy_Had_3 = new TH2F("fDCAxy_Had_3","p_{T} vs. DCA_{xy} (Hadron);p_{T} [GeV/c];DCA_{xy} [cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_Had_3);
    
  //pT vs. DCA (Like-sign)
    fDCAxy_LS_3 = new TH2F("fDCAxy_LS_3","p_{T} vs. DCA_{xy} (Like-sign pairs);p_{T} [GeV/c];DCA_{xy} [cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_LS_3);
      
  //pT vs. DCA (Unlike-sign)
    fDCAxy_ULS_3 = new TH2F("fDCAxy_ULS_3","p_{T} vs. DCA_{xy} (Unlike-sign pairs);p_{T} [GeV/c];DCA_{xy} [cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_ULS_3);


  //pT vs. DCA (Semi  electron)
    fDCAxy_Ele_4 = new TH2F("fDCAxy_Ele_4","p_{T} vs. DCA_{xy} (Electron);p_{T} [GeV/c];DCA_{xy} #times charge #times Bsign [cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_Ele_4);

  
  //electron pT (tight)
    fEopElectron2 = new TH1F("fEopElectron2","Electron;p_{T} [GeV/c];",600,0,30);
    fOutputList->Add(fEopElectron2);
    
  //hadron pT (tight)
    fEopHadron2 = new TH1F("fEopHadron2","Hadron;p_{T} [GeV/c];",600,0,30);
    fOutputList->Add(fEopHadron2);

  //electron pT (Except photonic(invariant mass))
    fEopElectron3 = new TH1F("fEopElectron3","Electron;p_{T} [GeV/c];",600,0,30);
    fOutputList->Add(fEopElectron3);

  //conversion R
    fHistConv_R = new TH2F("fHistConv_R","conversion R;p_{T} [GeV/c];R [cm]",600,0,30,500,0,50);
    fOutputList->Add(fHistConv_R);

  //Electron Eta vs. Phi
    fElectronEtaPhi = new TH2F("fElectronEtaPhi","Eta vs. Phi (electron)",180,-0.9,0.9,180,-0.9,0.9);
    fOutputList->Add(fElectronEtaPhi);
    
  //Hadron Eta vs. Phi
    fHadronEtaPhi = new TH2F("fHadronEtaPhi","Eta vs. Phi (hadron)",180,-0.9,0.9,180,-0.9,0.9);
    fOutputList->Add(fHadronEtaPhi);
    




    
//************************************ MC data ************************************//
  //Number of B
    fNoB = new TH1F("fNoB","Number of B event",13,-0.5,12.5);
    fOutputList->Add(fNoB);
    fNoB->GetYaxis()->SetTitle("counts");
    fNoB->GetXaxis()->SetBinLabel(1,"B->e");
    fNoB->GetXaxis()->SetBinLabel(2,"correct");
    fNoB->GetXaxis()->SetBinLabel(3,"Track mutch");
    fNoB->GetXaxis()->SetBinLabel(4,"AOD standard");
    fNoB->GetXaxis()->SetBinLabel(5,"TPC & ITS refit");
    fNoB->GetXaxis()->SetBinLabel(6,"TPC CrossedRow");
    fNoB->GetXaxis()->SetBinLabel(7,"ITS cluster");
    fNoB->GetXaxis()->SetBinLabel(8,"dE/dx calculation");
    fNoB->GetXaxis()->SetBinLabel(9,"SPD hit");
    fNoB->GetXaxis()->SetBinLabel(10,"chi2 cut");
    fNoB->GetXaxis()->SetBinLabel(11,"DCA cut");
    fNoB->GetXaxis()->SetBinLabel(12,"eta cut");
    fNoB->GetXaxis()->SetBinLabel(13,"after PID ");
  
  //Number of D
    fNoD = new TH1F("fNoD","Number of D event",13,-0.5,12.5);
    fOutputList->Add(fNoD);
    fNoD->GetYaxis()->SetTitle("counts");
    fNoD->GetXaxis()->SetBinLabel(1,"D->e & B->D->e");
    fNoD->GetXaxis()->SetBinLabel(2,"correct");
    fNoD->GetXaxis()->SetBinLabel(3,"Track mutch");
    fNoD->GetXaxis()->SetBinLabel(4,"AOD standard");
    fNoD->GetXaxis()->SetBinLabel(5,"TPC & ITS refit");
    fNoD->GetXaxis()->SetBinLabel(6,"TPC CrossedRow");
    fNoD->GetXaxis()->SetBinLabel(7,"ITS cluster");
    fNoD->GetXaxis()->SetBinLabel(8,"dE/dx calculation");
    fNoD->GetXaxis()->SetBinLabel(9,"SPD hit");
    fNoD->GetXaxis()->SetBinLabel(10,"chi2 cut");
    fNoD->GetXaxis()->SetBinLabel(11,"DCA cut");
    fNoD->GetXaxis()->SetBinLabel(12,"eta cut");
    fNoD->GetXaxis()->SetBinLabel(13,"after PID");


  //Total photonic electron(MC)
    fHistPho_Reco0 = new TH1F("fHistPho_Reco0", "Total photonic electron (MC); p_{T} [GeV/c];",1200,0,60);
    fOutputList->Add(fHistPho_Reco0);
    
  //Reconstructed photonic electron in data (MC)
    fHistPho_Reco1 = new TH1F("fHistPho_Reco1", "Reconstructed photonic electron in data (MC); p_{T} [GeV/c];",1200,0,60);
    fOutputList->Add(fHistPho_Reco1);

  //Non-Reconstructed photonic electron in data (MC)
    fHistPho_Reco2 = new TH1F("fHistPho_Reco2", "Non-Reconstructed photonic electron in data (MC); p_{T} [GeV/c];",1200,0,60);
    fOutputList->Add(fHistPho_Reco2);
    
  //check Eta range cut in MC
    fCheckEtaMC = new TH1F("fCheckEtaMC","check Eta range cut in MC", 160, -0.8, 0.8);
    fOutputList->Add(fCheckEtaMC);
    
  //MC org Pi0
    fHistMCorg_Pi0 = new TH2F("fHistMCorg_Pi0","MC org Pi0", 2, -0.5, 1.5, 1200, 0, 60);
    fOutputList->Add(fHistMCorg_Pi0);
    
  //MC org Eta
    fHistMCorg_Eta = new TH2F("fHistMCorg_Eta","MC org Eta", 2, -0.5, 1.5, 1200, 0, 60);
    fOutputList->Add(fHistMCorg_Eta);
    
  //MC org D (D->e)
    fHistMCorg_D = new TH1F("fHistMCorg_D","MC org D;p_{T} [GeV/c];", 1200, 0, 60);
    fOutputList->Add(fHistMCorg_D);
    
  //MC org B->D (B->D->e)
    fHistMCorg_BD = new TH1F("fHistMCorg_BD","MC org B->D;p_{T} [GeV/c];", 1200, 0, 60);
    fOutputList->Add(fHistMCorg_BD);
    
  //MC org B (B->e)
    fHistMCorg_B = new TH1F("fHistMCorg_B","MC org B;p_{T} [GeV/c];", 1200, 0, 60);
    fOutputList->Add(fHistMCorg_B);
    
  //MC org Lc (Lc->e)
    fHistMCorg_Lc = new TH1F("fHistMCorg_Lc","MC org Lc;p_{T} [GeV/c];", 1200, 0, 60);
    fOutputList->Add(fHistMCorg_Lc);
    
  //B meson(GM) pt vs electron pt
    fPt_Btoe = new TH2F("fPt_Btoe","B meson&baryon (M&GM) vs electron;electron p_{T} [GeV/c];B-meson p_{T} [GeV/c]",1200,0,60,1200,0,60);
    fOutputList->Add(fPt_Btoe);
    
  //HFE from B meson
    fHistPt_HFE_MC_B = new TH1F("fHistPt_HFE_MC_B","HFE from B;p_{T} [GeV/c];", 1200,0,60);
    fOutputList->Add(fHistPt_HFE_MC_B);
    
  //HFE from D meson
    fHistPt_HFE_MC_D = new TH1F("fHistPt_HFE_MC_D","HFE from D;p_{T} [GeV/c];", 1200,0,60);
    fOutputList->Add(fHistPt_HFE_MC_D);
    
  //HFE from Lc
    fHistPt_HFE_MC_Lc = new TH1F("fHistPt_HFE_MC_Lc","HFE from Lc;p_{T} [GeV/c];", 1200,0,60);
    fOutputList->Add(fHistPt_HFE_MC_Lc);
  
  //DCAxy from B
    fDCAxy_MC_B = new TH2F("fDCAxy_MC_B","p_{T} vs DCA_{xy} (MC : B-meson);p_{T} [GeV/c];DCA_{xy} #times charge #times Bsign[cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_MC_B);
  
  //DCAxy from D
    fDCAxy_MC_D = new TH2F("fDCAxy_MC_D","p_{T} vs DCA_{xy} (MC : D-meson);p_{T} [GeV/c];DCA_{xy} #times charge #times Bsign[cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_MC_D);

  //DCAxy from Dpm
    fDCAxy_MC_Dpm = new TH2F("fDCAxy_MC_Dpm","p_{T} vs DCA_{xy} (MC : D^{+},D^{*+});p_{T} [GeV/c];DCA_{xy} #times charge #times Bsign[cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_MC_Dpm);
  
  //DCAxy from D0
    fDCAxy_MC_D0 = new TH2F("fDCAxy_MC_D0","p_{T} vs DCA_{xy} (MC : D^{0},D^{*0});p_{T} [GeV/c];DCA_{xy} #times charge #times Bsign[cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_MC_D0);
  
  //DCAxy from Ds
    fDCAxy_MC_Ds = new TH2F("fDCAxy_MC_Ds","p_{T} vs DCA_{xy} (MC : D^{+}_{s},D^{*+}_{s});p_{T} [GeV/c];DCA_{xy} #times charge #times Bsign[cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_MC_Ds);
  
  //DCAxy from Lambda c 
    fDCAxy_MC_Lc = new TH2F("fDCAxy_MC_Lc","p_{T} vs DCA_{xy} (MC : #Lambda_{c});p_{T} [GeV/c];DCA_{xy} #times charge #times Bsign[cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_MC_Lc);

  //DCAxy elrctron
    fDCAxy_MC_ele = new TH2F("fDCAxy_MC_ele","p_{T} vs DCA_{xy} (MC : electron);p_{T} [GeV/c];DCA_{xy} #times charge #times Bsign[cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_MC_ele);
    
  //DCAxy photonic elrctron
    fDCAxy_MC_Phot = new TH2F("fDCAxy_MC_Phot","p_{T} vs DCA_{xy} (MC : photonic electron);p_{T} [GeV/c];DCA_{xy} #times charge #times Bsign[cm]",600,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_MC_Phot);

  //B Hist
    fHistPt_B_TrkCut0 = new TH1F("fHistPt_B_TrkCut0","B (matching track);p_{T} [GeV/c];",1200,0,60);	fOutputList->Add(fHistPt_B_TrkCut0);
    fHistPt_B_TrkCut1 = new TH1F("fHistPt_B_TrkCut1","B (AOD standard);p_{T} [GeV/c];",1200,0,60);	fOutputList->Add(fHistPt_B_TrkCut1);
    fHistPt_B_TrkCut2 = new TH1F("fHistPt_B_TrkCut2","B (TPC & ITS refit);p_{T} [GeV/c];",1200,0,60);	fOutputList->Add(fHistPt_B_TrkCut2);
    fHistPt_B_TrkCut3 = new TH1F("fHistPt_B_TrkCut3","B (TPCCrossedRow);p_{T} [GeV/c];",1200,0,60);	fOutputList->Add(fHistPt_B_TrkCut3);
    fHistPt_B_TrkCut4 = new TH1F("fHistPt_B_TrkCut4","B (ITS cluster cut);p_{T} [GeV/c];",1200,0,60);	fOutputList->Add(fHistPt_B_TrkCut4);
    fHistPt_B_TrkCut5 = new TH1F("fHistPt_B_TrkCut5","B (dE/dx calculation);p_{T} [GeV/c];",1200,0,60);	fOutputList->Add(fHistPt_B_TrkCut5);
    fHistPt_B_TrkCut6 = new TH1F("fHistPt_B_TrkCut6","B (SPD hit cut);p_{T} [GeV/c];",1200,0,60);	fOutputList->Add(fHistPt_B_TrkCut6);
    fHistPt_B_TrkCut7 = new TH1F("fHistPt_B_TrkCut7","B (chi2 cut);p_{T} [GeV/c];",1200,0,60);		fOutputList->Add(fHistPt_B_TrkCut7);
    fHistPt_B_TrkCut8 = new TH1F("fHistPt_B_TrkCut8","B (DCA cut);p_{T} [GeV/c];",1200,0,60);		fOutputList->Add(fHistPt_B_TrkCut8);
    fHistPt_B_TrkCut9 = new TH1F("fHistPt_B_TrkCut9","B (eta cut);p_{T} [GeV/c];",1200,0,60);		fOutputList->Add(fHistPt_B_TrkCut9);
    
  //D Hist
    fHistPt_D_TrkCut0 = new TH1F("fHistPt_D_TrkCut0","D (matching track);p_{T} [GeV/c];",1200,0,60);	fOutputList->Add(fHistPt_D_TrkCut0);
    fHistPt_D_TrkCut1 = new TH1F("fHistPt_D_TrkCut1","D (AOD standard);p_{T} [GeV/c];",1200,0,60);	fOutputList->Add(fHistPt_D_TrkCut1);
    fHistPt_D_TrkCut2 = new TH1F("fHistPt_D_TrkCut2","D (TPC & ITS refit);p_{T} [GeV/c];",1200,0,60);	fOutputList->Add(fHistPt_D_TrkCut2);
    fHistPt_D_TrkCut3 = new TH1F("fHistPt_D_TrkCut3","D (TPCCrossedRow);p_{T} [GeV/c];",1200,0,60);	fOutputList->Add(fHistPt_D_TrkCut3);
    fHistPt_D_TrkCut4 = new TH1F("fHistPt_D_TrkCut4","D (ITS cluster cut);p_{T} [GeV/c];",1200,0,60);	fOutputList->Add(fHistPt_D_TrkCut4);
    fHistPt_D_TrkCut5 = new TH1F("fHistPt_D_TrkCut5","D (dE/dx calculation);p_{T} [GeV/c];",1200,0,60);	fOutputList->Add(fHistPt_D_TrkCut5);
    fHistPt_D_TrkCut6 = new TH1F("fHistPt_D_TrkCut6","D (SPD hit cut);p_{T} [GeV/c];",1200,0,60);	fOutputList->Add(fHistPt_D_TrkCut6);
    fHistPt_D_TrkCut7 = new TH1F("fHistPt_D_TrkCut7","D (chi2 cut);p_{T} [GeV/c];",1200,0,60);		fOutputList->Add(fHistPt_D_TrkCut7);
    fHistPt_D_TrkCut8 = new TH1F("fHistPt_D_TrkCut8","D (DCA cut);p_{T} [GeV/c];",1200,0,60);		fOutputList->Add(fHistPt_D_TrkCut8);
    fHistPt_D_TrkCut9 = new TH1F("fHistPt_D_TrkCut9","D (eta cut);p_{T} [GeV/c];",1200,0,60);		fOutputList->Add(fHistPt_D_TrkCut9);

    

    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}

//_____________________________________________________________________________
void AliAnalysisTaskHFEBeautyMultiplicity::UserExec(Option_t *)
{
    // user exec
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you 
    // have access to the current event. 
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    // get an event (called fAOD) from the input file
                                                        // there's another event format (ESD) which works in a similar wya
                                                        // but is more cpu/memory unfriendly. for now, we'll stick with aod's
    if(!fAOD) return;                                   // if the pointer to the event is empty (getting it failed) skip this event
        // example part: i'll show how to loop over the tracks in an event 
        // and extract some information from them which we'll store in a histogram
    
    fMCarray  = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    
     
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if(!fVevent) {
        cout << "ERROR: fVEvent not available" << endl;
        return;
    }
    
    
    //__________ Track Cut __________
    Double_t CutTrackEta[2] = {-0.6, 0.6};
    Double_t CutTPCNCls = 100;
    Double_t CutITSNCls = 3;
    Double_t CutTPCdEdx = 80;
    Double_t CutDCAxy   = 2.4;
    Double_t CutDCAz    = 3.2;
    Int_t CutTPCNCrossedRow = 100;
    Double_t CutDiff	= 0.05;
    
    //___________ PID Cut __________
    Double_t CutTPCNsigma[2] = {-1.0, 3.0};
    Double_t CutM20[2] = {0.015, 0.3};
    Double_t CutEop[2] = {0.8, 1.2};
    Double_t CutHadNsigma = -3.5;
        


//************************//
//        if tender       //
//************************//
    if(fUseTender) {
        //new branches with calibrated tracks and clusters
        if(IsAODanalysis()) fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("tracks"));
    }
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());



//**********************//
//    PID initialized   //
//**********************//
    fpidResponse = fInputHandler -> GetPIDResponse();


//___________________________ Centrality ___________________________
    Double_t centrality = -1;
    AliCentrality *fCentrality = (AliCentrality*)fAOD->GetCentrality();
    //centrality = fCentrality->GetCentralityPercentile("V0M");
    
    Float_t lPercentile = 300;
    AliMultSelection *fMultSelection = 0x0;
    fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
    if(!fMultSelection) {
       //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
        centrality = fCentrality -> GetCentralityPercentile(fCentralityEstimator.Data());
    }else{
        lPercentile = fMultSelection->GetMultiplicityPercentile("V0M");
        centrality  = fMultSelection->GetMultiplicityPercentile(fCentralityEstimator.Data(), false);
    }




//***************************//
//     Event information     //
//***************************//
    Int_t ntracks = -999;
    if(!fUseTender) ntracks = fVevent -> GetNumberOfTracks();
    if( fUseTender) ntracks = fTracks_tender -> GetEntries();
    fMult -> Fill(centrality, ntracks);
    
    Int_t Bsign = 0;
    if(fAOD->GetMagneticField() < 0) Bsign = -1;
    if(fAOD->GetMagneticField() > 0) Bsign = 1;
    //cout << Bsign << endl;
    
    
//___________________________  events selection ___________________________
    //==== Global Vertex ====
    Double_t Xvertex = -100, Yvertex = -100, Zvertex = -100;
    const AliVVertex *pVtx = fVevent -> GetPrimaryVertex();
    Double_t NcontV = pVtx -> GetNContributors();
    Xvertex = pVtx -> GetX();   // X vertex
    Yvertex = pVtx -> GetY();   // Y vertex
    Zvertex = pVtx -> GetZ();   // Z vertex
    fVtxX->Fill(Xvertex);
    fVtxY->Fill(Yvertex);
    fVtxZ->Fill(Zvertex);
    //==== SPD Vertex ====
    const AliVVertex *pVtxSPD = fVevent -> GetPrimaryVertexSPD();
    Double_t ZvertexSPD = pVtxSPD -> GetZ();
    Double_t NcontVSPD  = pVtxSPD -> GetNContributors();
    Double_t cov[6]={0};
    pVtxSPD -> GetCovarianceMatrix(cov);

    fVtxCorrelation->Fill(Zvertex,ZvertexSPD);
    fNcont->Fill(NcontV, NcontVSPD);
    
    //----1.All events
    fNevents->Fill(0);

    //----2.pile-up event cut
    if(fVevent->IsPileupFromSPDInMultBins()) return;
    fNevents->Fill(1);
    
    //----3.Global track contribution cut
    if(NcontV < 2) return;
    fNevents->Fill(2);
 
    //----4.SPD track contribution cut
    if(NcontVSPD < 2) return;
    fNevents->Fill(3);

    //----5.SPD vertex & Global vertex match
    if(TMath::Abs(ZvertexSPD - Zvertex) > 0.5) return;
    fNevents->Fill(4);
    
    //----6.SPD vertex resolution cut
    if(TMath::Sqrt(cov[5]) > 0.25) return;
    fNevents->Fill(5);

    //----7.Zvertex cut
    if(TMath::Abs(Zvertex) > 10.0) return;
    fNevents->Fill(6);

    fVtxZ_2->Fill(Zvertex);	//Z vertex after event cut

//__________________________________________________________________________
    
    fCent -> Fill(centrality);  // centrality
    //cout << "centrality = " << centrality << " %" << endl;
    
    iBevt = kFALSE;	// b,bbar identificetion
    if(fMCarray) CheckMCgen(fMCheader, CutTrackEta[1]);   // True production of HFE



    //---------- SPD tracklets ----------//
    Int_t nTracklets = 0;
    Int_t nAcc = 0;
    Double_t etaRange = 1.0;

    AliAODTracklets *tracklets = static_cast<const AliAODEvent*>(fAOD)->GetTracklets();
    nTracklets = tracklets->GetNumberOfTracklets();

    for(Int_t nn=0; nn<nTracklets; nn++) {
	Double_t theta = tracklets->GetTheta(nn);
	Double_t eta = tracklets->GetEta(nn);
	if(TMath::Abs(eta) < etaRange) nAcc++;	// No. of tracklet in |eta|<1.0 (TPC coverage)
    }

    fZvtx_Ntrklet->Fill(Zvertex, nAcc);



    //-------- SPD tracklets correction --------//
    Double_t correctednAcc = nAcc;
    Double_t fRefMult = Nref;
    //Double_t WeightNtrklet = -1.;
    //Double_t WeightZvtx = -1.;
    TProfile* estimatorAvg;
    if(!fMCarray) estimatorAvg = GetEstimatorHistogram(fAOD);
    if(fMCarray)  estimatorAvg = GetEstimatorHistogramMC(fAOD);

    if(estimatorAvg) {
	    correctednAcc = static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg, nAcc, Zvertex,fRefMult));
    }

    fZvtx_Ntrklet_Corr->Fill(Zvertex, correctednAcc);





    
//*********************//
//        EMCAL        //
//*********************//
    
    fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("tracks"));
    fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("caloClusters"));   //<- Get cluster&track information updated by correction framework
    
    
    Int_t Nclust = fVevent->GetNumberOfCaloClusters(); //Get number of Cluster

    int NclustAll = 0;
    int NclustE1  = 0;	// of clust E > 0.1
    int NclustE2  = 0;	// of clust E > 0.2
    int NclustE3  = 0;	// of clust E > 0.5

    Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;
    
//___________________________  EMCal cluster loop ___________________________
    for(Int_t icl=0; icl<Nclust; icl++)
       {
        AliVCluster *clust = 0x0;     
        //clust = (AliVCluster*)fVevent->GetCaloCluster(icl); // address cluster matched to track(cluster numbering)
        clust = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(icl));        //EMCal correction framework

        if(clust && clust->IsEMCAL())
          {
              Double_t clustE = clust->E();       //cluster energy
              // if(clustE < 0.001) continue;

              Float_t emcx[3];                                  //cluster position
              clust -> GetPosition(emcx);                       //Get EMCal position
              TVector3 clustpos(emcx[0], emcx[1], emcx[2]);
              Double_t emcphi = clustpos.Phi();                 //phi information (-pi < phi < pi)
              if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi());   //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
              
              Double_t emceta = clustpos.Eta();                 //eta information

              fHistClustE -> Fill(clustE);                      //EMCal cluster energy distribution

              fEMCClsEtaPhi -> Fill(emceta, emcphi);
              fHistNCells->Fill(clustE,clust->GetNCells());

              NclustAll++;
              if(clustE > 0.1) NclustE1++;
              if(clustE > 0.2) NclustE2++;
              if(clustE > 0.5) NclustE3++;
              
          }
           
       }
    
    fHistNClsE -> Fill(NclustAll);
    fHistNClsE1-> Fill(NclustE1);
    fHistNClsE2-> Fill(NclustE2);
    fHistNClsE3-> Fill(NclustE3);
// _________________________________________________________________________________


    //---- cell Information ----//
    AliVCaloCells *fCaloCells = fVevent -> GetEMCALCells();


    Short_t cellAddr, nSACell;
    Int_t mclabel;
    Short_t iSACell;
    Double_t cellAmp=-1., cellTimeT=-1., clusterTime=-1., efrac=-1.;

    nSACell = fCaloCells->GetNumberOfCells();
    for(iSACell = 0; iSACell < nSACell; iSACell++)
    {
	    Bool_t haveCell = fCaloCells->GetCell(iSACell, cellAddr, cellAmp, cellTimeT, mclabel, efrac);
	    if(haveCell)fHistCalCells->Fill(cellAddr,cellAmp);
    }



     
//*********************//
//        track        //
//*********************//

//___________________________  track loop ___________________________
   
    //Int_t iTracks = fAOD->GetNumberOfTracks();     // see how many tracks there are in the event
    Int_t iTracks(fTracks_tender->GetEntries());     // EMCal correction framework

    for(Int_t i(0); i < iTracks; i++)
    {                                                                            // loop overall these tracks
        //AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));     // get a track (type AliAODTrack) from the event
        
        AliAODTrack* track = dynamic_cast<AliAODTrack*>(fTracks_tender->At(i));  // EMCal correction framework

        if(!track) continue;                                                     // if we failed, skip this track
        
        
    //********************//
    // Get MC Information //
    //********************//
        Int_t ilabel = TMath::Abs(track->GetLabel());
        Int_t pdg = -999;
        Double_t pid_ele = 0.0;
        Double_t pTMom = -1.0;
        Int_t pidM = -1;
        Int_t ilabelM = -1;
        Double_t pTGMom = -1.0;
        Int_t pidGM = -1;
        Int_t ilabelGM = -1;
        Bool_t pid_eleD = kFALSE;
        Bool_t pid_eleB = kFALSE;
        Bool_t pid_eleP = kFALSE;
        Bool_t iEmbPi0 = kFALSE;
        Bool_t iEmbEta = kFALSE;
        
        if(ilabel>0 && fMCarray)
        {
            fMCTrackpart = (AliAODMCParticle*) fMCarray->At(ilabel);
            pdg = fMCTrackpart->GetPdgCode();
            if(TMath::Abs(pdg)==11) pid_ele = 1.0;
            if(pid_ele==1.0) FindMother(fMCTrackpart, ilabelM, pidM, pTMom);
            
            pid_eleD = IsDdecay(pidM);  // D->e, B->D->e
            pid_eleB = IsBdecay(pidM);  // B->e
            pid_eleP = IsPdecay(pidM);  // photon->e

	    if(pid_eleB) fNoB -> Fill(0);
	    if(pid_eleD) fNoD -> Fill(0);
            

            if(pid_eleD && iBevt)    // mother is D-meson, but GM is B-meson&baryon
            {
                AliAODMCParticle* fMCTrackpartMom = (AliAODMCParticle*) fMCarray->At(ilabelM);
                FindMother(fMCTrackpartMom, ilabelGM, pidGM, pTGMom);
                
                pid_eleB = kTRUE;
                pid_eleD = kFALSE;
		pidM = pidGM;
		pTMom = pTGMom;
            }
            
	    if(pid_eleB) fNoB -> Fill(1);
	    if(pid_eleD) fNoD -> Fill(1);

            
            if(pidM==111)   //pi0
            {
                if(ilabelM >= NembMCpi0 && ilabelM < NembMCeta) iEmbPi0 = kTRUE;
                if(ilabelM >= NembMCeta && ilabelM < NpureMCproc) iEmbEta = kTRUE;
            }
            
            if(pidM==221)   //eta
            {
                if(ilabelM >= NembMCeta && ilabelM < NpureMCproc) iEmbEta = kTRUE;
            }
 
            if(pidM==22)    //electron from photon
            {
                AliAODMCParticle* fMCparticleM = (AliAODMCParticle*) fMCarray->At(ilabelM);
                FindMother(fMCparticleM, ilabelM, pidM, pTMom);
                
                if(pidM==111)
                {
                    if(ilabelM >= NembMCpi0 && ilabelM < NembMCeta) iEmbPi0 = kTRUE;
                    if(ilabelM >= NembMCeta && ilabelM < NpureMCproc) iEmbEta = kTRUE;
                }
                if(pidM==221)
                {
                    if(ilabelM >= NembMCeta && ilabelM < NpureMCproc) iEmbEta = kTRUE;
                }
            }
        }
        
        if(pidM==443)continue; // remove enhanced J/psi in MC !
        if(pidM==-99)continue; // remove e from no mother !
        
        
        
        
    //******************//
    // Track properties //
    //******************//
    //___________________________  All track information ___________________________
        Double_t TrkPt = -999, TrkEta = -999, TrkPhi = -999, TrkP = -999, charge = -999;
        Double_t dEdx = -999, TPCnSigma = -999, TPCnSigma_Pi = -999;
        Double_t TOFnSigma = -999, ITSnSigma = -999;
        Double_t ITSchi2 = -999, TPCchi2NDF = -999, TPCCrossedRows = -999;


        TrkP   = track -> P();
        TrkPt  = track -> Pt();
        TrkEta = track -> Eta();
        TrkPhi = track -> Phi();
        dEdx   = track -> GetTPCsignal();
        charge = track -> Charge();
        TPCchi2NDF = track -> Chi2perNDF();
        ITSchi2 = track -> GetITSchi2();
	TPCCrossedRows = track -> GetTPCCrossedRows();
        
        fAllTrkPt  -> Fill(TrkPt);      //All track (Pt)
        fAllTrkEta -> Fill(TrkEta);     //All track (Eta)
        fAllTrkPhi -> Fill(TrkPhi);     //All track (Phi)
        fdEdx      -> Fill(TrkP,dEdx);  //All track (P vs dE/dx)
	fPhiEta	   -> Fill(TrkPhi,TrkEta);
	fTPCCrossedRow -> Fill(TPCCrossedRows);

        
        //---- electron ----//
        TPCnSigma = fpidResponse -> NumberOfSigmasTPC(track, AliPID::kElectron);
        TOFnSigma = fpidResponse -> NumberOfSigmasTOF(track, AliPID::kElectron);
        ITSnSigma = fpidResponse -> NumberOfSigmasITS(track, AliPID::kElectron);
        fTPCnsig -> Fill(TrkP,TPCnSigma);
        fTOFnsig -> Fill(TrkP,TOFnSigma);
        fITSnsig -> Fill(TrkP,ITSnSigma);
        
        //---- pion ----//
        TPCnSigma_Pi = fpidResponse -> NumberOfSigmasTPC(track, AliPID::kPion);
        fTPCnsig_Pi -> Fill(TrkP,TPCnSigma_Pi);


        
        if(TrkPt>2.0)fTPCnsigEta0 -> Fill(TrkEta,TPCnSigma);    //track pT > 2 GeV/c
        if(TrkPt>3.0)fTPCnsigEta1 -> Fill(TrkEta,TPCnSigma);    //track pT > 3 GeV/c
        if(TrkPt>5.0)fTPCnsigEta2 -> Fill(TrkEta,TPCnSigma);    //track pT > 5 GeV/c
        
        
        
        
    //************************//
    //Track matching to EMCAL //
    //************************//
        if(!track->IsEMCAL()) continue;
        Int_t EMCalIndex = -1;
        EMCalIndex = track -> GetEMCALcluster();  // get index of EMCal cluster which matched to track

        
        //___________________________  EMCal matched track ___________________________
        if(EMCalIndex < 0) continue;

        fEMCTrkPt->Fill(TrkPt);     //track matched EMCal (Pt)
        
        AliVCluster *clustMatch=0x0;
        clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex); // address cluster matched to track

        
        Double_t emceta = -999, emcphi = -999;
        if(clustMatch && clustMatch->IsEMCAL())
        {
        //___________________________  track selection ___________________________
          //---- 0.matching tracks ----
            fNtracks->Fill(0);
	    	if(pid_eleB){ fNoB -> Fill(2); fHistPt_B_TrkCut0 -> Fill(TrkPt);}
	    	if(pid_eleD){ fNoD -> Fill(2); fHistPt_D_TrkCut0 -> Fill(TrkPt);}
            
          //---- 1.AOD standard track cut ----
            if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;
            fNtracks->Fill(1);
	    	if(pid_eleB){ fNoB -> Fill(3); fHistPt_B_TrkCut1 -> Fill(TrkPt);}
	    	if(pid_eleD){ fNoD -> Fill(3); fHistPt_D_TrkCut1 -> Fill(TrkPt);}
            
          //---- 2.TPC and ITS refit cut ----
            if((!(track->GetStatus()&AliESDtrack::kITSrefit) || (!(track->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
            fNtracks->Fill(2);
	    	if(pid_eleB){ fNoB -> Fill(4); fHistPt_B_TrkCut2 -> Fill(TrkPt);}
	    	if(pid_eleD){ fNoD -> Fill(4); fHistPt_D_TrkCut2 -> Fill(TrkPt);}
            
          //---- 3.TPC cluster cut ----
            //if(track->GetTPCNcls() < CutTPCNCls) continue;
            //fNtracks->Fill(3);

	  //---- 3.TPC CrossedRow cut ----
            if(TPCCrossedRows < CutTPCNCrossedRow) continue;
            fNtracks->Fill(3);
	    	if(pid_eleB){ fNoB -> Fill(5); fHistPt_B_TrkCut3 -> Fill(TrkPt);}
	    	if(pid_eleD){ fNoD -> Fill(5); fHistPt_D_TrkCut3 -> Fill(TrkPt);}
            
          //---- 4.ITS cluster cut ----
            if(track->GetITSNcls() < CutITSNCls) continue;
            fNtracks->Fill(4);
	    	if(pid_eleB){ fNoB -> Fill(6); fHistPt_B_TrkCut4 -> Fill(TrkPt);}
	    	if(pid_eleD){ fNoD -> Fill(6); fHistPt_D_TrkCut4 -> Fill(TrkPt);}
            
          //---- 5.TPC cluster cut for dE/dx calculation ----
            if(track->GetTPCsignalN() < CutTPCdEdx) continue;
            fNtracks->Fill(5);
	    	if(pid_eleB){ fNoB -> Fill(7); fHistPt_B_TrkCut5 -> Fill(TrkPt);}
	    	if(pid_eleD){ fNoD -> Fill(7); fHistPt_D_TrkCut5 -> Fill(TrkPt);}
            
          //---- 6.SPD hit cut ----
            if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) continue;
            fNtracks->Fill(6);
	    	if(pid_eleB){ fNoB -> Fill(8); fHistPt_B_TrkCut6 -> Fill(TrkPt);}
	    	if(pid_eleD){ fNoD -> Fill(8); fHistPt_D_TrkCut6 -> Fill(TrkPt);}

	  //---- 7.chi2 cut ----
	    if((ITSchi2 >= 25) || (TPCchi2NDF >= 4)) continue;
            fNtracks->Fill(7);
	    	if(pid_eleB){ fNoB -> Fill(9); fHistPt_B_TrkCut7 -> Fill(TrkPt);}
	    	if(pid_eleD){ fNoD -> Fill(9); fHistPt_D_TrkCut7 -> Fill(TrkPt);}
            
            
          //---- 8.DCA cut ----
            Double_t DCA[2], CovarianceMatrix[3];
            if(track->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., DCA, CovarianceMatrix))
            if(TMath::Abs(DCA[0]) > CutDCAxy || TMath::Abs(DCA[1]) > CutDCAz) continue;
            fNtracks->Fill(8);
	    	if(pid_eleB){ fNoB -> Fill(10); fHistPt_B_TrkCut8 -> Fill(TrkPt);}
	    	if(pid_eleD){ fNoD -> Fill(10); fHistPt_D_TrkCut8 -> Fill(TrkPt);}
            

	  //---- 9.eta & phi diff cut ----

            // calculate phi and eta difference between a track and a cluster
            Double_t fPhiDiff = -999, fEtaDiff = -999;
            GetTrkClsEtaPhiDiff(track, clustMatch, fPhiDiff, fEtaDiff);
            fHistEMCTrkMatch_Eta -> Fill(fEtaDiff);         //delta Eta
            fHistEMCTrkMatch_Phi -> Fill(fPhiDiff);         //delta Phi
            fEMCTrkMatch_EtaPhi  -> Fill(fEtaDiff,fPhiDiff);//delta Eta vs delta Phi

	    if(TMath::Abs(fPhiDiff) > CutDiff || TMath::Abs(fEtaDiff) > CutDiff) continue;
	    fNtracks->Fill(9);
            fEMCTrkMatch_EtaPhi_AfterCut -> Fill(fEtaDiff,fPhiDiff);






        //*******************************//
        // Select EMCal or DCal clusters //
        //*******************************//
            Float_t emcx[3]; // cluster position
            clustMatch -> GetPosition(emcx);
            TVector3 clustpos(emcx[0], emcx[1], emcx[2]);
            emceta = clustpos.Eta();
            emcphi = clustpos.Phi();


        //*******************************************//
        // Properties of tracks matched to the EMCal //
        //*******************************************//
	
            fTrkPt_2   -> Fill(TrkPt);      //after cut track (Pt)
            fTrkEta_2  -> Fill(TrkEta);     //after cut track (Eta)
            fTrkPhi_2  -> Fill(TrkPhi);     //after cut track (Phi)
            fdEdx_2    -> Fill(TrkP,dEdx);  //after cut track (P vs dE/dx)
            fTPCnsig_2 -> Fill(TrkP,TPCnSigma);
            fTOFnsig_2 -> Fill(TrkP,TOFnSigma);
            fITSnsig_2 -> Fill(TrkP,ITSnSigma);
	    fTPCCrossedRow_2 -> Fill(TPCCrossedRows);

	
	
            if(TrkPt>1.0)   //reject low-pt track
            {
                fEMCTrkEta->Fill(TrkEta);   //track matched EMCal (Eta)
                fEMCTrkPhi->Fill(TrkPhi);   //track matched EMCal (Phi)
            }

            fClsEtaPhiAftMatch->Fill(emceta,emcphi);

            if(TrkPhi > 1.396 && TrkPhi < 3.264){    //EMCal acceptance (80 to 187 degrees)
                fClsEtaPhiAftMatchEMCin -> Fill(emceta,emcphi);} // inside
            else{
                fClsEtaPhiAftMatchEMCout -> Fill(emceta,emcphi);}// outside


        //---- EMCal Electron IDentification info ----//
            Double_t eop = -1.0;
            Double_t m02 = -99999.0, m20 = -99999.0, sqm02m20 = -99999.0;
            Double_t clustMatchE = clustMatch->E();
            if(TrkP>0) eop = clustMatchE/TrkP;
            m02 = clustMatch -> GetM02();   // long axis
            m20 = clustMatch -> GetM20();   // short axis
            sqm02m20 = sqrt(pow(m02,2) + pow(m20,2));
            
            fM02_1 -> Fill(TrkPt,m02);
            fM20_1 -> Fill(TrkPt,m20);


        //*******************************//
        //    Electron Identification    //
        //*******************************//

            //---- 10.Eta cut ----
            if(TrkEta > CutTrackEta[1] && TrkEta < CutTrackEta[0]) continue;
            fNtracks->Fill(10);
	    	if(pid_eleB){ fNoB -> Fill(11); fHistPt_B_TrkCut9 -> Fill(TrkPt);}
	    	if(pid_eleD){ fNoD -> Fill(11); fHistPt_D_TrkCut9 -> Fill(TrkPt);}

            fHistEopAll -> Fill(eop);
	    fTrkEtaPhi_AfterCut->Fill(TrkEta, TrkPhi);
            
            Bool_t fFlagNonHFE = kFALSE;    // photonic electron identification
            
            //========== Electron E/p ==========//
            if((TPCnSigma >= CutTPCNsigma[0] && TPCnSigma <= CutTPCNsigma[1]) && (m20 >= CutM20[0] && m20 <= CutM20[1]))  // TPC nsigma & shower shape cut
            {
                fEopElectron1 -> Fill(TrkPt,eop);
                
                if(eop >= CutEop[0] && eop <= CutEop[1]) // E/p cut (with TPC nsigma & shower shape cut)
                {
                    fEopElectron2 -> Fill(TrkPt); 
		    fElectronEtaPhi -> Fill(TrkEta, TrkPhi);

                    fM02_2 -> Fill(TrkPt,m02);
                    fM20_2 -> Fill(TrkPt,m20);
                    fDCAxy_Ele_1 -> Fill(TrkPt, DCA[0]*charge*Bsign);
                    fDCAxy_Ele_2 -> Fill(TrkPt, DCA[0]*charge);
                    fDCAxy_Ele_3 -> Fill(TrkPt, DCA[0]);
                    //cout << "pT = " << TrkPt << " ,  DCAxy = " << DCA[0] << " ,  charge = " << charge << endl;
                    
                    
                    //---- Photonic electron ----
                    SelectPhotonicElectron(i, track, fFlagNonHFE, pidM, TrkPt, DCA[0], Bsign);
                    if(!fFlagNonHFE){
                        fEopElectron3 -> Fill(TrkPt);
                        fDCAxy_Ele_4 -> Fill(TrkPt, DCA[0]*charge*Bsign);
                    }
                    
                    
		    if(pid_ele == 1.0)
		    {
			fDCAxy_MC_ele -> Fill(TrkPt, DCA[0]*charge*Bsign);	//DCA (all electron)

                    	//-------- Photonic electron --------//
                    	if(pid_eleP)    // electron from photon & pi0 & eta (<- TMath::Abs(pdg)==11 && (pidM==22 || pidM==111 || pidM==221))
                    	{
                        	fHistPho_Reco0->Fill(TrkPt);     // all information of photonic electron
				fDCAxy_MC_Phot -> Fill(TrkPt, DCA[0]*charge*Bsign);	//DCA (total photonic electron)
                       
                       		if(fFlagNonHFE)
                       		{
                          		fHistPho_Reco1->Fill(TrkPt); // Reconstructed by EMCal & TPC & InvMass
                        	}
                        	else
                       		{
                           		fHistPho_Reco2->Fill(TrkPt); // Non-Reconstructed by EMCal & TPC & InvMass
                        	}
                    	}


		    	//-------- Heavy Flavour electron --------//
                    	if(pid_eleB)
			{	
				fNoB -> Fill(12);
			    	fHistPt_HFE_MC_B -> Fill(track->Pt()); // HFE from B meson&baryon (MC)
			    	fDCAxy_MC_B -> Fill(TrkPt, DCA[0]*charge*Bsign);
		    	}

                    	if(pid_eleD)
			{
				fNoD -> Fill(12);
			    	fHistPt_HFE_MC_D -> Fill(track->Pt()); // HFE from D meson (MC)
			    	fDCAxy_MC_D -> Fill(TrkPt, DCA[0]*charge*Bsign);

			    	if(TMath::Abs(pidM)==411 || TMath::Abs(pidM)==413) fDCAxy_MC_Dpm -> Fill(TrkPt, DCA[0]*charge*Bsign);
			    	if(TMath::Abs(pidM)==421 || TMath::Abs(pidM)==423) fDCAxy_MC_D0  -> Fill(TrkPt, DCA[0]*charge*Bsign);
			    	if(TMath::Abs(pidM)==431 || TMath::Abs(pidM)==433) fDCAxy_MC_Ds  -> Fill(TrkPt, DCA[0]*charge*Bsign);

		    		if(TMath::Abs(pidM)==4122)
				{			   // HFE from Lambda c (MC)
			    		fHistPt_HFE_MC_Lc -> Fill(track->Pt());
			    		fDCAxy_MC_Lc -> Fill(TrkPt, DCA[0]*charge*Bsign);
		    		}
		    	}
		    }
                }
            }
            
            
            //========== Hadron E/p ==========//
            if((TPCnSigma <= CutHadNsigma) && (m20 >= CutM20[0] && m20 <= CutM20[1]))
            {
                fEopHadron1 -> Fill(TrkPt,eop);

                if(eop >= CutEop[0] && eop <= CutEop[1])  // Hadron DCA in 0.8<E/p<1.2 
                {
                    fEopHadron2 -> Fill(TrkPt);
		    fHadronEtaPhi -> Fill(TrkEta, TrkPhi);

                    fDCAxy_Had_1 -> Fill(TrkPt, DCA[0]*charge*Bsign); 
                    fDCAxy_Had_2 -> Fill(TrkPt, DCA[0]*charge);
                    fDCAxy_Had_3 -> Fill(TrkPt, DCA[0]);
                }

            }
            

        }
  } // continue until all the tracks are processed (track roop)
//_________________________________________________________________________________________________


PostData(1, fOutputList);               // stream the results the analysis of this event to
                                        // the output manager which will take care of writing
                                        // it to a file
}




//_________________________________________________________________________________________________
void AliAnalysisTaskHFEBeautyMultiplicity::GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff)
{
	//Calculate phi and eta difference between a track and a cluster.
	//The position of the track is obtained on the EMCal surface.

	phidiff = 999;
	etadiff = 999;

	if(!t||!v) return;

	Double_t veta = t -> GetTrackEtaOnEMCal();
	Double_t vphi = t -> GetTrackPhiOnEMCal();

	Float_t pos[3] = {0., 0., 0.};
	v -> GetPosition(pos);
	TVector3 cpos(pos);
	Double_t ceta = cpos.Eta();
	Double_t cphi = cpos.Phi();
	etadiff = veta-ceta;
	phidiff = TVector2::Phi_mpi_pi(vphi-cphi);

}

//_________________________________________________________________________________________________
void AliAnalysisTaskHFEBeautyMultiplicity::SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC, Double_t TrkPt, Double_t DCAxy, Int_t Bsign)
{
    //*******************************//
    // Non-HFE Invariant mass method //
    //*******************************//
    
   /* AliESDtrackCuts* esdTrackCutsAsso = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    esdTrackCutsAsso->SetAcceptKinkDaughters(kFALSE);
    esdTrackCutsAsso->SetRequireTPCRefit(kTRUE);
    esdTrackCutsAsso->SetRequireITSRefit(kTRUE);
    esdTrackCutsAsso->SetEtaRange(-0.9,0.9);
    esdTrackCutsAsso->SetMaxChi2PerClusterTPC(4);
    esdTrackCutsAsso->SetMinNClustersTPC(70);
    esdTrackCutsAsso->SetMaxDCAToVertexZ(3.2);
    esdTrackCutsAsso->SetMaxDCAToVertexXY(2.4);
    esdTrackCutsAsso->SetDCAToVertex2D(kTRUE);
    */
    
    Bool_t flagPhotonicElec = kFALSE;
    
    int ntracks = -999;
    if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
    if(fUseTender) ntracks = fTracks_tender->GetEntries();
    
    
    for (Int_t jtrack = 0; jtrack < ntracks; jtrack++)
    {
        AliVParticle* VAssotrack = 0x0;
        if(!fUseTender) VAssotrack  = fVevent->GetTrack(jtrack);
        if(fUseTender) VAssotrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jtrack)); //take tracks from Tender list
        
        if (!VAssotrack) {
            printf("ERROR: Could not receive track %d\n", jtrack);
            continue;
        }
        
        AliVTrack *Assotrack = dynamic_cast<AliVTrack*>(VAssotrack);
        AliESDtrack *eAssotrack = dynamic_cast<AliESDtrack*>(VAssotrack);
        AliAODTrack *aAssotrack = dynamic_cast<AliAODTrack*>(VAssotrack);
        
        
        //----reject same track
        if(jtrack == itrack) continue;
        if(aAssotrack->Px()==track->Px() && aAssotrack->Py()==track->Py() && aAssotrack->Pz()==track->Pz()) continue;
        
        
        Bool_t fFlagLS = kFALSE, fFlagULS = kFALSE;
        Double_t ptAsso = -999., AssoTrackNsigma = -999.0, mass = -999., width = -999, R = -999, R_error = -999;
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
        
        AssoTrackNsigma = fpidResponse -> NumberOfSigmasTPC(Assotrack, AliPID::kElectron);
        ptAsso = Assotrack -> Pt();
        Int_t chargeAsso = Assotrack -> Charge();
        Int_t charge = track -> Charge();
        Double_t AssoTPCchi2perNDF = aAssotrack -> Chi2perNDF();
        if(charge > 0) fPDGe1 = -11;
        if(chargeAsso > 0) fPDGe2 = -11;
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        //----track cuts applied
        if(fAOD) {
            if(!aAssotrack -> TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            if(aAssotrack -> GetTPCNcls() < 80) continue;
            if(aAssotrack -> GetITSNcls() < 3 ) continue;
            if((!(aAssotrack -> GetStatus()&AliESDtrack::kITSrefit)|| (!(aAssotrack -> GetStatus()&AliESDtrack::kTPCrefit)))) continue;
        }
        else{
            //if(!esdTrackCutsAsso -> AcceptTrack(eAssotrack)) continue;
        }
        
        //-----loose cut on partner electron
        if(ptAsso < 0.1) continue;
        if(aAssotrack->Eta()<-0.9 || aAssotrack->Eta()>0.9) continue;
        if(AssoTrackNsigma < -3 || AssoTrackNsigma > 3) continue;
        //if(AssoTPCchi2perNDF >= 4) continue;
        
        
        //----define KFParticle to get mass
        AliKFParticle::SetField(fVevent->GetMagneticField());
        AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
        AliKFParticle ge2 = AliKFParticle(*Assotrack, fPDGe2);
        AliKFParticle recg(ge1, ge2);
    
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg)) > 3.) continue;
        
        //-------Get mass
        Int_t MassCorrect, RCorrect;
        MassCorrect = recg.GetMass(mass,width);
	RCorrect = recg.GetR(R,R_error);

        
        if(fFlagLS){    // Like-sign
            if(TrkPt >= 1.0){
                fInvmassLS -> Fill(TrkPt,mass);
                if(mass <= 0.15){
                    fDCAxy_LS_1 -> Fill(TrkPt, DCAxy*charge*Bsign);
                    fDCAxy_LS_2 -> Fill(TrkPt, DCAxy*charge);
                    fDCAxy_LS_3 -> Fill(TrkPt, DCAxy);
                }
            }
        }
        
        if(fFlagULS){   // Unlike-sign
            if(TrkPt >= 1.0){
                fInvmassULS -> Fill(TrkPt,mass);
                if(mass <= 0.15){
                    fDCAxy_ULS_1 -> Fill(TrkPt, DCAxy*charge*Bsign);
                    fDCAxy_ULS_2 -> Fill(TrkPt, DCAxy*charge);
                    fDCAxy_ULS_3 -> Fill(TrkPt, DCAxy);

		    fHistConv_R -> Fill(TrkPt,R);
                }
            }
        }
       
        if(mass <= 0.15 && fFlagULS && !flagPhotonicElec) flagPhotonicElec = kTRUE; // Tag Non-HFE (photonic electron by Invariant-mass method)

    }
    fFlagPhotonicElec = flagPhotonicElec;
}

//_________________________________________________________________________________________________
void AliAnalysisTaskHFEBeautyMultiplicity::FindMother(AliAODMCParticle *part, int &label, int &pid, double &ptMom)
{
    if(part -> GetMother() > -1)
    {
        label = part -> GetMother();
        AliAODMCParticle *partM = (AliAODMCParticle*)fMCarray -> At(label);
        pid = partM -> GetPdgCode();
        ptMom = partM -> Pt();
    }
    else
    {
        pid = -99;
    }
}

//_________________________________________________________________________________________________
Bool_t AliAnalysisTaskHFEBeautyMultiplicity::IsDdecay(int mpid)      // D mason
{
    int abmpid = TMath::Abs(mpid);
    // D+ : 411,        D0 : 421,        D*+ : 413,       D*0 : 423,       Ds+ : 431,       Ds*+ : 433       Lc : 4122
    if(abmpid == 411 || abmpid == 421 || abmpid == 413 || abmpid == 423 || abmpid == 431 || abmpid == 433 || abmpid == 4122)
      {
          return kTRUE;
      }
    else
      {
          return kFALSE;
      }
}

//_________________________________________________________________________________________________
Bool_t AliAnalysisTaskHFEBeautyMultiplicity::IsBdecay(int mpid)      // B mason
{
    int abmpid = TMath::Abs(mpid);
    // B0 : 511,        B+ : 521,        B*0 : 513,       B*+ : 523,       Bs0 : 531,       Bs*0 : 533
    if(abmpid == 511 || abmpid == 521 || abmpid == 513 || abmpid == 523 || abmpid == 531 || abmpid == 533)
      {
          return kTRUE;
      }
    else
      {
          return kFALSE;
      }
}

//_________________________________________________________________________________________________
Bool_t AliAnalysisTaskHFEBeautyMultiplicity::IsPdecay(int mpid)      // photon
{
    int abmpid = TMath::Abs(mpid);
    // photon : 22,    pi0 : 111,       eta : 221
    if(abmpid == 22 || abmpid == 111 || abmpid == 221)
      {
          return kTRUE;
      }
    else
      {
          return kFALSE;
      }
}

//_________________________________________________________________________________________________
void AliAnalysisTaskHFEBeautyMultiplicity::CheckMCgen(AliAODMCHeader* fMCheader, Double_t CutEta)
{
    TList *lh=fMCheader->GetCocktailHeaders();
    NpureMC = 0;
    NpureMCproc = 0;
    NembMCpi0 = 0;
    NembMCeta = 0;
    Nch = 0;
    TString MCgen;
    TString embpi0("pi");
    TString embeta("eta");
    TString embbeauty("bele");
    TString embcharm("cele");


    if(lh)
    {
        for(int igene=0; igene<lh->GetEntries(); igene++)
        {
            AliGenEventHeader* gh = (AliGenEventHeader*)lh->At(igene);
            if(gh)
            {
                MCgen = gh -> GetName();    // Get particle name
                if(igene == 0) NpureMC = gh -> NProduced();
                
                if(MCgen.Contains(embpi0))NembMCpi0 = NpureMCproc;  //if "pi" contains
                if(MCgen.Contains(embeta))NembMCeta = NpureMCproc;  //if "eta" contains

		if(MCgen.Contains(embbeauty))
		{
			iBevt = kTRUE;		 // b,bbar
		}
                
                NpureMCproc += gh->NProduced();  // generate by PYTHIA or HIJING
            }
        }
    }

    
    for(int imc=0; imc < fMCarray->GetEntries(); imc++)
    {
        Bool_t iEnhance = kFALSE;
        if(imc >= NpureMC) iEnhance = kTRUE;
        Int_t iHijing = 1;  // select particle from Hijing or PYTHIA
        
        fMCparticle = (AliAODMCParticle*) fMCarray -> At(imc);
        Int_t pdgGen = TMath::Abs(fMCparticle->GetPdgCode());
        
        Double_t pdgEta = fMCparticle->Eta();                   // eta
        Double_t pTtrue = fMCparticle->Pt();                    // Pt
        Int_t chargetrue = fMCparticle->Charge();               // charge
        Bool_t isPhysPrim = fMCparticle->IsPhysicalPrimary();   // ?
        
        //Get N Charge
        //if(chargetrue!=0 && TMath::Abs(pdgEta)<1.0 && isPhysPrim) Nch++;    // ?
        
        
        if(TMath::Abs(pdgEta)>CutEta) continue;
        fCheckEtaMC -> Fill(pdgEta);
        
        Int_t pdgMom = -99;
        Int_t labelMom = -1;
        Double_t pTMom = -1.0;
        
        FindMother(fMCparticle, labelMom, pdgMom, pTMom);
        if(pdgMom == -99 && iEnhance) iHijing = 0;  // particles from enhance
        if(pdgMom > 0 && iEnhance) iHijing = -1;    // particles from enhance but feeddown
        
        if(iHijing > -1)
        {
            if(pdgGen==111) fHistMCorg_Pi0 -> Fill(iHijing, pTtrue); // pi0
            if(pdgGen==221) fHistMCorg_Eta -> Fill(iHijing, pTtrue); // eta
        }
        
        if(TMath::Abs(pdgGen)!=11) continue;    //except Non-electrons
        if(pTMom < 2.0) continue;
        
        Int_t pdgGM = -99;
        Int_t labelGM = -1;
        Double_t pTGMom = -1.0;
        
        if(pdgMom!=0)
        {
            AliAODMCParticle* fMCparticleMom = (AliAODMCParticle*) fMCarray -> At(labelMom);
            if(IsDdecay(pdgMom))
            {
                fHistMCorg_D -> Fill(fMCparticle->Pt());             // D->e
                FindMother(fMCparticleMom, labelGM, pdgGM, pTGMom);
                
                //if(IsBdecay(pdgGM))
		if(iBevt)
                {
                    fHistMCorg_BD -> Fill(fMCparticle->Pt());        // B->D->e
                    fPt_Btoe -> Fill(fMCparticle->Pt(), pTGMom);
                }
                
            }
            
            if(IsBdecay(pdgMom))
            {
                fHistMCorg_B->Fill(fMCparticle->Pt());               // B->e
                fPt_Btoe->Fill(fMCparticle->Pt(), pTMom);
            }

	    if(TMath::Abs(pdgMom)==4122)			     // Lambda_c -> e
	    {
		fHistMCorg_Lc -> Fill(fMCparticle->Pt());
	    }
	   
            
        }
    }
    
    return;
}



//_________________________________________________________________________________________________
TProfile* AliAnalysisTaskHFEBeautyMultiplicity::GetEstimatorHistogram(const AliAODEvent* fAOD)
{
	//Int_t runNo = fAOD->GetRunNumber();
	//Int_t period = 1;

	return fMultiEstimatorAvg;

}


//_________________________________________________________________________________________________
TProfile* AliAnalysisTaskHFEBeautyMultiplicity::GetEstimatorHistogramMC(const AliAODEvent* fAOD)
{
	//Int_t runNo = fAOD->GetRunNumber();
	//Int_t period = 2;

	return fMultiEstimatorAvg;
}


//_________________________________________________________________________________________________
Double_t AliAnalysisTaskHFEBeautyMultiplicity::GetCorrectedNtrackletsD(TProfile* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult)
{
	if(TMath::Abs(vtxZ)>10.0){
		return uncorrectedNacc;	// out of range for correction of multiplicity
	}

	if(!estimatorAvg){
		cout << "ERROR:: Missing TProfile for correction of multiplicity" << endl;
		return uncorrectedNacc;	// missing TProfile
	}

	Double_t localAvg = estimatorAvg->GetBinContent(estimatorAvg->FindBin(vtxZ));
	Double_t deltaM = 0;
	deltaM = uncorrectedNacc*(refMult/localAvg -1);


	Double_t correctedNacc = -1;
	correctedNacc = uncorrectedNacc + (deltaM>0 ? 1 : -1) * gRandom->PoissonD(TMath::Abs(deltaM));

	if(correctedNacc < 0) correctedNacc = 0;

	return correctedNacc;

}


//_________________________________________________________________________________________________
void AliAnalysisTaskHFEBeautyMultiplicity::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_________________________________________________________________________________________________
