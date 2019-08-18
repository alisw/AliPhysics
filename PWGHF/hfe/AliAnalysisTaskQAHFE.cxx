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
//
// QA task
// 
// Authors:
//   Shreyasi Acharya <shreyasi.acharya@cern.ch>
//

#include <Riostream.h>
#include "TChain.h"
#include "TTree.h"
 
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliPPVsMultUtils.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliVEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskQAHFE.h"
#include "AliNormalizationCounter.h"
#include <TFile.h>
#include <THnSparse.h>
#include "AliEventPoolManager.h"
#include "TParticle.h"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include "AliHFEtools.h"
#include "AliCFManager.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include  "AliTOFPIDResponse.h"
#include  "AliTPCPIDResponse.h"

#include "TChain.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliEMCALTrack.h"
#include "AliExternalTrackParam.h"
#include "AliPhysicsSelection.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"
#include "AliLog.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliSelectNonHFE.h"
#include "AliHFEpidTPC.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"
#include "AliESDHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliEventPoolManager.h"
#include "TObjArray.h"
#include "AliMultSelection.h"

#include "AliSelectNonHFE.h"
#include "TClonesArray.h"
#include "AliAODMCParticle.h"

#include <TProfile.h>
#include "AliKFParticle.h"
#include "AliVertexingHFUtils.h"
#include "AliAnalysisUtils.h"
#include "AliESDUtils.h"
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskQAHFE)

//________________________________________________________________________//
AliAnalysisTaskQAHFE::AliAnalysisTaskQAHFE(): 
AliAnalysisTaskSE(),

//------------General Variables
fAOD(0), 
fPID(0), 
fpidResponse(0),
fAODVertex(0), 
fOutputList(0), 
fcount(0), 

//-----------

ftrigger(AliVEvent::kINT7),
fTPCNclus(100),
fITSNclus(3),
fTPCNclusPID(80),
fSPDBoth(kTRUE),
fSPDAny(kFALSE),
fSPDFirst(kFALSE),
	
fDCAxy(0),
fDCAz(0),
fDCAxyCut(1),
fDCAzCut(2),
	
fpTMin(0.3),
fEtarange(0.8),
fFilterbit(AliAODTrack::kTrkGlobalNoDCA),
fTPCnsigmin(-1),
fTPCnsigmax(3),
fTOFnsig(3),
	
fInvmassCut(0.14),	
fAssoTPCCluster(60),
fAssoITSRefit(kTRUE),
fAssopTMin(0.1),
fAssoEtarange(0.9),
fAssoTPCnsig(3.5),




//-------------Analysis
fHistPt(0), fHistMult(0), fTPCSignal(0), fNentries(0), fNentries2(0),
fVtxZ(0),fVtxZ_corr(0),fVtxX(0),fVtxY(0),
fDCAZvsPt(0),fDCAXYvsPt(0),
fdEdxVsP_TPC(0),
fnSigmaVsP_TPC(0),
fdEdxVsP_TPC_cut(0),
fnSigmaVsP_TPC_cut(0),
fnSigmaVsP_TOF(0), 
fnBetaVsP_TOF(0),
fnSigmaVsEta_TPC(0), 
fnSigmaVsEta_TPC_cut(0), 

fdEdxVsPt_TPC(0),
fnSigmaVsPt_TPC(0),
fnSigmaVsPt_TOF(0), 
fnBetaVsPt_TOF(0), 

fLandau(0),fErr(0),fHadCot_Landau(0),fHadCot_Err(0),fPt_incl_e_Landau(0),fPt_incl_e_Err(0),

//Multiplicity dep

fPteV0M_Lan(0),fPteV0M_Err(0), fPteSPD_Lan(0),fPteSPD_Err(0),
fPteV0M(0),fPteSPD(0),
fHadV0MLan(0),fHadSPDLan(0),fHadV0MErr(0),fHadSPDErr(0),
/*fnSigmaVsP_TPC_cut_multV0M(0),fnSigmaVsP_TPC_cut_multSPD(0),*/
/*fnSigmaVsPt_TPC_cut_multV0M(0),fnSigmaVsPt_TPC_cut_multSPD(0),*/

//------------Ntracklets Correction and Multiplicity------------------

fMultiPercentile(0),fMultiPercentileSPD(0),
fV0Mult(0),fV0MultCorr(0),fV0AMult(0),fV0CMult(0),fSPD_tracklet_NoEtaCut(0),fSPD_tracklet(0), 
fCent(0),fCentSPD(0), 
fMultV0M_vs_alimult(0),fMultSPD_vs_alimult(0),
fileEstimator(0),
fSPDCorrMultDistNoEtaCut_min(0),fSPDCorrMultDistNoEtaCut_max(0),fSPDCorrMultDist_min(0),fSPDCorrMultDist_max(0),
fHistNtrVsZvtxNoEtaCut(0),fHistNtrCorrVsZvtxNoEtaCut_min(0),fHistNtrCorrVsZvtxNoEtaCut_max(0),fHistNtrVsZvtx(0),fHistNtrCorrVsZvtx_min(0),fHistNtrCorrVsZvtx_max(0),fSPDtr_V0Mmult(0),
fRefMult(13.079),
fV0MultVsZvtx(0),fV0MultCorrVsZvtx(0),
fSPDCorrMultDist_min_vs_AliSPDMult(0),
fNchVsZvtx(0),fSPDNtrCorrVsNch(0),	fV0MCorrVsNch(0),/*fSPDNtrCorrVsV0MCorrVsNch(0),*/
 
//------------------------Non-Hfe
fNonHFE(new AliSelectNonHFE()),
fInvmassLS1(0),fInvmassULS1(0),
fPte_ULS(0),fPte_LS(0),
fPte_ULS_KF(0),fPte_LS_KF(0),
fPte_ULS_multV0M(0),fPte_LS_multV0M(0),
fPte_ULS_multSPD(0),fPte_LS_multSPD(0),

//------------------------------MC 
fIsMC(kFALSE),
fMCarray(0),fMCparticle(0),fMCmother(0),

fPtHFEMCtot(0),fPtHFEMCtot_SPD(0),fPtHFEMCtot_V0M(0),
fPtHFEMC(0),
fPtHFEMC_SPD(0),
fPtHFEMC_V0M(0),

fPtHFEMC_afterfilterbit(0),
fPtHFEMC_aftertrackcuts(0),  
fPtHFEMC_aftertrackcuts_SPD(0),  
fPtHFEMC_aftertrackcuts_V0M(0), 
fPtHFEMC_aftertracktofcuts(0),
fPtHFEMC_aftertracktofcuts_SPD(0),
fPtHFEMC_aftertracktofcuts_V0M(0),

fPt_elec_phot(0),
fPt_elec_phot_multV0M(0),
fPt_elec_phot_multSPD(0),

fPt_elec_phot1(0),
fPt_elec_phot1_multV0M(0),
fPt_elec_phot1_multSPD(0),

fpt_e_nonphot_MC(0),
fElecNos(0),

pi0MC(0),
etaMC(0),
gammaMC(0),

pi0MC_1(0),
etaMC_1(0),
gammaMC_1(0),


fPT_elec_MCtrue(0),
fpt_incl_elec_MC(0),
fPi0_Pt_MC(0),fEta_Pt_MC(0),fGamma_Pt_MC(0),
fPte_LS_invMss_MC(0),fPte_ULS_invMss_MC(0),
fPte_LS_MC(0),fPte_ULS_MC(0),
fPte_LS_MC_multV0M(0),fPte_ULS_MC_multV0M(0),
fPte_LS_MC_multSPD(0),fPte_ULS_MC_multSPD(0),

fInvMULS(0),
//fInvMULS_multV0M(0),
//fInvMULS_multSPD(0),

pdg(0),

fPt_elec_pi0_MB(0),
fPt_elec_pi0_MB_ULS(0),
fPt_elec_pi0_MB_LS(0),
fPt_elec_eta_MB(0),
fPt_elec_eta_MB_ULS(0),
fPt_elec_eta_MB_LS(0),
fPt_elec_gamma_MB(0),
fPt_elec_gamma_MB_ULS(0),
fPt_elec_gamma_MB_LS(0),

fInvMULSnSp(0),
fInvMULSpi0(0),
fInvMULSeta(0),
fInvMULSgamma(0),

fPtHFEMC_reco(0),
fPtHFEMC_reco_true(0),
fPtHFEMC_reco_SPD(0),
fPtHFEMC_reco_V0M(0),

fPi0EtaSpectra(0),
fPi0EtaSpectraSp(0)


{
  fPID = new AliHFEpid("hfePid");
  for(Int_t i=0; i<2; i++) fMultEstimatorAvg[i]=0;
}


AliAnalysisTaskQAHFE::AliAnalysisTaskQAHFE(const char *name) 
  : AliAnalysisTaskSE(name), 
 
//------------General Variables
fAOD(0), 
fPID(0), 
fpidResponse(0),
fAODVertex(0), 
fOutputList(0), 
fcount(0), 

//-----------
ftrigger(AliVEvent::kINT7),
fTPCNclus(100),
fITSNclus(3),
fTPCNclusPID(80),
fSPDBoth(kTRUE),
fSPDAny(kFALSE),
fSPDFirst(kFALSE),
	
fDCAxy(0),
fDCAz(0),
fDCAxyCut(1),
fDCAzCut(2),
	
fpTMin(0.3),
fEtarange(0.8),
fFilterbit(AliAODTrack::kTrkGlobalNoDCA),
fTPCnsigmin(-1),
fTPCnsigmax(3),
fTOFnsig(3),

fInvmassCut(0.14),
fAssoTPCCluster(60),
fAssoITSRefit(kTRUE),
fAssopTMin(0.1),
fAssoEtarange(0.9),
fAssoTPCnsig(3.5),

//-------------Analysis
fHistPt(0), fHistMult(0), fTPCSignal(0), fNentries(0), fNentries2(0),
fVtxZ(0),fVtxZ_corr(0),fVtxX(0),fVtxY(0),
fDCAZvsPt(0),fDCAXYvsPt(0),
fdEdxVsP_TPC(0),
fnSigmaVsP_TPC(0),
fdEdxVsP_TPC_cut(0),
fnSigmaVsP_TPC_cut(0),
fnSigmaVsP_TOF(0), 
fnBetaVsP_TOF(0),
fnSigmaVsEta_TPC(0), 
fnSigmaVsEta_TPC_cut(0),  

fdEdxVsPt_TPC(0),
fnSigmaVsPt_TPC(0),
fnSigmaVsPt_TOF(0), 
fnBetaVsPt_TOF(0),

fLandau(0),fErr(0),fHadCot_Landau(0),fHadCot_Err(0),fPt_incl_e_Landau(0),fPt_incl_e_Err(0),

//Multiplicity dep

fPteV0M_Lan(0),fPteV0M_Err(0), fPteSPD_Lan(0),fPteSPD_Err(0),
fPteV0M(0),fPteSPD(0),
fHadV0MLan(0),fHadSPDLan(0),fHadV0MErr(0),fHadSPDErr(0),
/*fnSigmaVsP_TPC_cut_multV0M(0),fnSigmaVsP_TPC_cut_multSPD(0),*/
/*fnSigmaVsPt_TPC_cut_multV0M(0),fnSigmaVsPt_TPC_cut_multSPD(0),*/

//------------Ntracklets Correction and Multiplicity------------------

fMultiPercentile(0),fMultiPercentileSPD(0),
fV0Mult(0),fV0MultCorr(0),fV0AMult(0),fV0CMult(0),fSPD_tracklet_NoEtaCut(0),fSPD_tracklet(0), 
fCent(0),fCentSPD(0), 
fMultV0M_vs_alimult(0),fMultSPD_vs_alimult(0),
fileEstimator(0),
fSPDCorrMultDistNoEtaCut_min(0),fSPDCorrMultDistNoEtaCut_max(0),fSPDCorrMultDist_min(0),fSPDCorrMultDist_max(0),
fHistNtrVsZvtxNoEtaCut(0),fHistNtrCorrVsZvtxNoEtaCut_min(0),fHistNtrCorrVsZvtxNoEtaCut_max(0),fHistNtrVsZvtx(0),fHistNtrCorrVsZvtx_min(0),fHistNtrCorrVsZvtx_max(0),fSPDtr_V0Mmult(0),
fRefMult(13.079),
fV0MultVsZvtx(0),fV0MultCorrVsZvtx(0),
fSPDCorrMultDist_min_vs_AliSPDMult(0),
fNchVsZvtx(0),fSPDNtrCorrVsNch(0),	fV0MCorrVsNch(0),/*fSPDNtrCorrVsV0MCorrVsNch(0),*/
 
//------------------------Non-Hfe
fNonHFE(new AliSelectNonHFE()),
fInvmassLS1(0),fInvmassULS1(0),
fPte_ULS(0),fPte_LS(0),
fPte_ULS_KF(0),fPte_LS_KF(0),
fPte_ULS_multV0M(0),fPte_LS_multV0M(0),
fPte_ULS_multSPD(0),fPte_LS_multSPD(0),

//------------------------------MC 
fIsMC(kFALSE),
fMCarray(0),fMCparticle(0),fMCmother(0),

fPtHFEMCtot(0),fPtHFEMCtot_SPD(0),fPtHFEMCtot_V0M(0),
fPtHFEMC(0),
fPtHFEMC_SPD(0),
fPtHFEMC_V0M(0),

fPtHFEMC_afterfilterbit(0),
fPtHFEMC_aftertrackcuts(0),  
fPtHFEMC_aftertrackcuts_SPD(0),  
fPtHFEMC_aftertrackcuts_V0M(0), 
fPtHFEMC_aftertracktofcuts(0),
fPtHFEMC_aftertracktofcuts_SPD(0),
fPtHFEMC_aftertracktofcuts_V0M(0),

fPt_elec_phot(0),
fPt_elec_phot_multV0M(0),
fPt_elec_phot_multSPD(0),

fPt_elec_phot1(0),
fPt_elec_phot1_multV0M(0),
fPt_elec_phot1_multSPD(0),

fpt_e_nonphot_MC(0),
fElecNos(0),

pi0MC(0),
etaMC(0),
gammaMC(0),

pi0MC_1(0),
etaMC_1(0),
gammaMC_1(0),

fPT_elec_MCtrue(0),
fpt_incl_elec_MC(0),
fPi0_Pt_MC(0),fEta_Pt_MC(0),fGamma_Pt_MC(0),
fPte_LS_invMss_MC(0),fPte_ULS_invMss_MC(0),
fPte_LS_MC(0),fPte_ULS_MC(0),
fPte_LS_MC_multV0M(0),fPte_ULS_MC_multV0M(0),
fPte_LS_MC_multSPD(0),fPte_ULS_MC_multSPD(0),

fInvMULS(0),
//fInvMULS_multV0M(0),
//fInvMULS_multSPD(0),

pdg(0),

fPt_elec_pi0_MB(0),
fPt_elec_pi0_MB_ULS(0),
fPt_elec_pi0_MB_LS(0),
fPt_elec_eta_MB(0),
fPt_elec_eta_MB_ULS(0),
fPt_elec_eta_MB_LS(0),
fPt_elec_gamma_MB(0),
fPt_elec_gamma_MB_ULS(0),
fPt_elec_gamma_MB_LS(0),

fInvMULSnSp(0),
fInvMULSpi0(0),
fInvMULSeta(0),
fInvMULSgamma(0),

fPtHFEMC_reco(0),
fPtHFEMC_reco_true(0),
fPtHFEMC_reco_SPD(0),
fPtHFEMC_reco_V0M(0),

fPi0EtaSpectra(0),
fPi0EtaSpectraSp(0)




{
  // Constructor
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container

  fPID = new AliHFEpid("hfePid");
  DefineInput(0, TChain::Class());  
  DefineOutput(1, TList::Class());
  DefineOutput(2, TH1F::Class());
  for(Int_t i=0; i<2; i++) fMultEstimatorAvg[i]=0;
}

//_________________________Destructer_____________________________________
AliAnalysisTaskQAHFE::~AliAnalysisTaskQAHFE()
{
  delete fPID;
  if (fOutputList) { delete fOutputList; fOutputList = 0;}
  if (fNentries){ delete fNentries; fNentries = 0;}
  if (fNentries2){ delete fNentries2; fNentries2 = 0;}
  
     for(Int_t i=0; i<2; i++) {
      if (fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i];
  }
	//delete fPIDqa;
}

//_______________________________________________________________________
void AliAnalysisTaskQAHFE::Init()
{
  // Initialization
  	
	
  if(fDebug > 1) printf("AliAnalysisTaskQAHFE::Init() \n");
  return;
}




//________________________________________________________________________
void AliAnalysisTaskQAHFE::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
	//Initialize PID

	fErr=new TF1("fErr","[0]+[1]*TMath::Erf([2]*x-[3])",0,4.5);
	fErr->SetParameters(0.489999,0.490001,0.501613,3.58221);
	
	fLandau=new TF1("fLandau","[0]*TMath::Landau(x,[1],[2])",0,4.5);
	fLandau->SetParameters(1.15545,10.5563,2.64482);
	
	if(fIsMC)
	{
	
		fErr->SetParameters(2.43937,2.43937,0.313236,3.79046);
		fLandau->SetParameters(0.0269315,9.17024,2.38605);

	}


  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
	
  Double_t pi=TMath::Pi();

  fcount = new TH1D("fcount", "fcount", 10, 0.0, 10.0);
  fcount->Sumw2();  
  
  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 15, 0.1, 15.1);
  fOutputList->Add(fHistPt);
  fHistMult = new TH1F("fHistMult", "Multiplicity distribution", 500, 0, 500);
  fOutputList->Add(fHistMult);
  fTPCSignal = new TH1F("fTPCSignal", "fTPCSignal distribution", 100, 0, 100);
  fOutputList->Add(fTPCSignal);
	
	//--------------------------Vertex Dist------------------------------------
	fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
  fOutputList->Add(fVtxZ);

	fVtxZ_corr=new TH1F("fVtxZ_corr","Z vertex position;Vtx_{z};counts",1000,-50,50);
  fOutputList->Add(fVtxZ_corr);

  fVtxY = new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",1000,-50,50);
  fOutputList->Add(fVtxY);

  fVtxX = new TH1F("fVtxX","X vertex position;Vtx_{x};counts",1000,-50,50);
  fOutputList->Add(fVtxX);
  
 
  //--------------------------DCA Dist------------------------------------
	fDCAZvsPt = new TH2F("fDCAZvsPt","DCAZ ;p_{T};DCA_{z}",300,0,100,100,-5,5);
  fOutputList->Add(fDCAZvsPt);

  fDCAXYvsPt = new TH2F("fDCAXYvsPt","DCAXY ;p_{T};DCA_{xy}",300,0,100,100,-5,5);
  fOutputList->Add(fDCAXYvsPt);

	
	//--------------------------TPC dedx nsig vs p ----------------------------------------------
  fdEdxVsP_TPC = new TH2F("fdEdxVsP_TPC", "fdEdxVsP_TPC distribution",300,0,15,750,10,160);
  fdEdxVsP_TPC->Sumw2();
  fOutputList->Add(fdEdxVsP_TPC);
  
  fnSigmaVsP_TPC= new TH2F("fnSigmaVsP_TPC", "fnSigmaVsP_TPC distribution",300,0.,15.,750,-15.,15.);
  fnSigmaVsP_TPC->Sumw2(); 
  fOutputList->Add(fnSigmaVsP_TPC);
  
  fdEdxVsP_TPC_cut = new TH2F("fdEdxVsP_TPC_cut", "fdEdxVsP_TPC_cut distribution",300,0,15,750,10,160);
  fdEdxVsP_TPC_cut->Sumw2();
  fOutputList->Add(fdEdxVsP_TPC_cut);
  
 	fnSigmaVsP_TPC_cut = new TH2F("fnSigmaVsP_TPC_cut", "fnSigmaVsP_TPC_cut distribution",300,0,15,250,-15,10); 
  fnSigmaVsP_TPC_cut->Sumw2(); 
  fOutputList->Add(fnSigmaVsP_TPC_cut);

  
  fnSigmaVsP_TOF= new TH2F("fnSigmaVsP_TOF", "fnSigmaVsP_TOF distribution",300,0.,15.,1000,-10.,20.);
  fnSigmaVsP_TOF->Sumw2(); 
  fOutputList->Add(fnSigmaVsP_TOF);
  
  fnBetaVsP_TOF= new TH2F("fnBetaVsP_TOF", "fnBetaVsP_TOF distribution",300,0.,15.,1000,-0.5,1.5);
	fnBetaVsP_TOF->Sumw2(); 
	fOutputList->Add(fnBetaVsP_TOF);


 //--------------------------TPC dedx nsig vs pt ----------------------------------------------
  fdEdxVsPt_TPC = new TH2F("fdEdxVsPt_TPC", "fdEdxVsPt_TPC distribution",300,0,15,750,10,160);
  fdEdxVsPt_TPC->Sumw2();
  fOutputList->Add(fdEdxVsPt_TPC);
  
  fnSigmaVsPt_TPC= new TH2F("fnSigmaVsPt_TPC", "fnSigmaVsPt_TPC distribution",300,0.,15.,750,-15.,15.);
  fnSigmaVsPt_TPC->Sumw2(); 
  fOutputList->Add(fnSigmaVsPt_TPC);
  
  fnSigmaVsPt_TOF= new TH2F("fnSigmaVsPt_TOF", "fnSigmaVsPt_TOF distribution",300,0.,15.,1000,-10.,20.);
  fnSigmaVsPt_TOF->Sumw2(); 
  fOutputList->Add(fnSigmaVsPt_TOF);
  
  fnBetaVsPt_TOF= new TH2F("fnBetaVsPt_TOF", "fnBetaVsPt_TOF distribution",300,0.,15.,1000,-0.5,1.5);
	fnBetaVsPt_TOF->Sumw2(); 
	fOutputList->Add(fnBetaVsPt_TOF);


 
 //--------------------------TPC  nsig vs eta ----------------------------------------------
  fnSigmaVsEta_TPC = new TH2F("fnSigmaVsEta_TPC", "fnSigmaVsEta_TPC distribution",1600,-0.8,0.8,750,-15.,15.);
  fnSigmaVsEta_TPC->Sumw2();
  fOutputList->Add(fnSigmaVsEta_TPC);
  
  fnSigmaVsEta_TPC_cut = new TH2F("fnSigmaVsEta_TPC_cut", "fnSigmaVsEta_TPC_cut distribution",1600,-0.8,0.8,750,-15.,15.);
  fnSigmaVsEta_TPC_cut->Sumw2(); 
  fOutputList->Add(fnSigmaVsEta_TPC_cut);
  
 
  //-----------------Incl e Spectrum and Hadron Contamination---------------------------------------------
  fHadCot_Landau = new TH1D("fHadCot_Landau", "fHadCot_Landau Distribution",300,0,15);
  fHadCot_Landau->Sumw2();
  fOutputList->Add(fHadCot_Landau); 
  
  fHadCot_Err = new TH1D("fHadCot_Err", "fHadCot_Err Distribution",300,0,15);
  fHadCot_Err->Sumw2();
  fOutputList->Add(fHadCot_Err);
  
  fPt_incl_e_Landau = new TH1D("fPt_incl_e_Landau", "fPt_incl_e_Landau Distribution",300,0,15);
  fPt_incl_e_Landau->Sumw2();
  fOutputList->Add(fPt_incl_e_Landau); 
  
  fPt_incl_e_Err = new TH1D("fPt_incl_e_Err", "fPt_incl_e_Err Distribution",300,0,15);
  fPt_incl_e_Err->Sumw2();
  fOutputList->Add(fPt_incl_e_Err);
  
  //------------------Centrality and Multiplicity------------------------------------------------- 
 
  fV0Mult = new TH1F("fV0Mult", "fV0Mult Distribution",1000,0,1000);
  fV0Mult->Sumw2();
 // fOutputList->Add(fV0Mult);
  
  fV0MultCorr = new TH1F("fV0MultCorr", "fV0MultCorr Distribution",1000,0,1000);
  fV0MultCorr->Sumw2();
 // fOutputList->Add(fV0MultCorr);
  
  fV0AMult = new TH1F("fV0AMult", "fV0AMult Distribution",1000,0,1000);
  fV0AMult->Sumw2();
  //fOutputList->Add(fV0AMult);
  
  fV0CMult = new TH1F("fV0CMult", "fV0CMult Distribution",1000,0,1000);
  fV0CMult->Sumw2();
 // fOutputList->Add(fV0CMult); 
  
  
 	//-----------------AliMultSelection Class------------------------------------------------------
 	fCent = new TH1F("fCent","Centrality",100,0,100);
 	fCent->Sumw2();
  fOutputList->Add(fCent);
  
  fCentSPD = new TH1F("fCentSPD","Centrality",100,0,100);
 //	fCentSPD->Sumw2();
  fOutputList->Add(fCentSPD);
  
	fMultV0M_vs_alimult = new TH2F("fMultV0M_vs_alimult","Track multiplicity",1000,0,1000.,100,0,100);
//  fOutputList->Add(fMultV0M_vs_alimult);
  
  fMultSPD_vs_alimult = new TH2F("fMultSPD_vs_alimult","Track multiplicity after cuts",200,0,200,100,0,100);
//  fOutputList->Add(fMultSPD_vs_alimult);
  
  //................................
  
  fPteV0M_Lan = new TH2F("fPteV0M_Lan","Track multiplicity vs pt",300,0,15.,1000,0,1000.);
	fPteV0M_Lan->Sumw2();
  //fOutputList->Add(fPteV0M_Lan);
  
  fPteV0M_Err = new TH2F("fPteV0M_Err","Track multiplicity vs pt",300,0,15.,1000,0,1000.);
	fPteV0M_Err->Sumw2();
  //fOutputList->Add(fPteV0M_Err);
  
  fPteSPD_Lan = new TH2F("fPteSPD_Lan","Track multiplicity vs pt",300,0,15.,200,0,200.);
  fPteSPD_Lan->Sumw2();
 // fOutputList->Add(fPteSPD_Lan);
  
  fPteSPD_Err = new TH2F("fPteSPD_Err","Track multiplicity vs pt",300,0,15.,200,0,200.);
  fPteSPD_Err->Sumw2();
  //fOutputList->Add(fPteSPD_Err);
 
 //................................
 
 	fPteV0M= new TH2F("fPteV0M","Track multiplicity vs pt",300,0,15.,1000,0,1000.);
	fPteV0M->Sumw2();
 // fOutputList->Add(fPteV0M);
  
  fPteSPD= new TH2F("fPteSPD","Track multiplicity vs pt",300,0,15.,200,0,200.);
	fPteSPD->Sumw2();
  //fOutputList->Add(fPteSPD);
  	
  //................................
  
  fHadV0MLan= new TH2F("fHadV0MLan","Track multiplicity vs pt",300,0,15.,1000,0,1000.);
	fHadV0MLan->Sumw2();
  //fOutputList->Add(fHadV0MLan);
  
  fHadSPDLan= new TH2F("fHadSPDLan","Track multiplicity vs pt",300,0,15.,200,0,200.);
	fHadSPDLan->Sumw2();
  //fOutputList->Add(fHadSPDLan);
  
  fHadV0MErr= new TH2F("fHadV0MErr","Track multiplicity vs pt",300,0,15.,1000,0,1000.);
	fHadV0MErr->Sumw2();
  //fOutputList->Add(fHadV0MErr);
  
  fHadSPDErr= new TH2F("fHadSPDErr","Track multiplicity vs pt",300,0,15.,200,0,200.);
	fHadSPDErr->Sumw2();
  //fOutputList->Add(fHadSPDErr);
  
  /*fnSigmaVsP_TPC_cut_multV0M= new TH3F("fnSigmaVsP_TPC_cut_multV0M","",100,0,5.,250,-15,10,800,0,800.);
	fnSigmaVsP_TPC_cut_multV0M->Sumw2();
  fOutputList->Add(fnSigmaVsP_TPC_cut_multV0M);
  
  fnSigmaVsP_TPC_cut_multSPD= new TH3F("fnSigmaVsP_TPC_cut_multSPD","",300,0,15.,250,-15,10,200,0,200.);
	fnSigmaVsP_TPC_cut_multSPD->Sumw2();
  fOutputList->Add(fnSigmaVsP_TPC_cut_multSPD);*/
/*  
  fnSigmaVsPt_TPC_cut_multV0M= new TH3F("fnSigmaVsPt_TPC_cut_multV0M","",100,0,5.,250,-15,10,1000,0,1000.);
	fnSigmaVsPt_TPC_cut_multV0M->Sumw2();
  fOutputList->Add(fnSigmaVsPt_TPC_cut_multV0M);
  
  fnSigmaVsPt_TPC_cut_multSPD= new TH3F("fnSigmaVsPt_TPC_cut_multSPD","",300,0,15.,250,-15,10,200,0,200.);
	fnSigmaVsPt_TPC_cut_multSPD->Sumw2();
  fOutputList->Add(fnSigmaVsPt_TPC_cut_multSPD);*/
  
  //-----------------NonHFE-----------------------------------------------------------------------------------------
  fPte_ULS = new TH1F("fPte_ULS", "ULS electron pt",300,0,15);
  fPte_ULS->Sumw2();
  fOutputList->Add(fPte_ULS);
  
  fPte_LS = new TH1F("fPte_LS", "LS electron pt",300,0,15);
  fPte_LS->Sumw2();
  fOutputList->Add(fPte_LS);
  
  fPte_ULS_KF = new TH1F("fPte_ULS_KF", "ULS electron_KF pt",300,0,15);
  fPte_ULS_KF->Sumw2();
  fOutputList->Add(fPte_ULS_KF);
  
  fPte_LS_KF = new TH1F("fPte_LS_KF", "LS electron_KF pt",300,0,15);
  fPte_LS_KF->Sumw2();
  fOutputList->Add(fPte_LS_KF);
  
  fInvmassLS1 = new TH1F("fInvmassLS1","Inv mass of LS (e,e) for pt^{e}; mass(GeV/c^2); counts;",1000,0,1.0);
	fOutputList->Add(fInvmassLS1);
	
	fInvmassULS1 = new TH1F("fInvmassULS1","Inv mass of ULS (e,e) for pt^{e}; mass(GeV/c^2); counts;",1000,0,1.0);
	fOutputList->Add(fInvmassULS1);
	
	//....................................................................................
	
	fPte_ULS_multV0M = new TH2F("fPte_ULS_multV0M", "ULS electron pt in percentile",300,0,15,1000,0,1000);
	fPte_ULS_multV0M->Sumw2();
  //fOutputList->Add(fPte_ULS_multV0M);
  
  fPte_LS_multV0M = new TH2F("fPte_LS_multV0M", "LS electron pt in percentile",300,0,15,1000,0,1000);
  fPte_LS_multV0M->Sumw2();
  //fOutputList->Add(fPte_LS_multV0M);
  
 
  fPte_ULS_multSPD = new TH2F("fPte_ULS_multSPD", "ULS electron pt in percentile",300,0,15,200,0,200);
  fPte_ULS_multSPD->Sumw2();
  //fOutputList->Add(fPte_ULS_multSPD);
  
  fPte_LS_multSPD = new TH2F("fPte_LS_multSPD", "LS electron pt in percentile",300,0,15,200,0,200);
  fPte_LS_multSPD->Sumw2();
  //fOutputList->Add(fPte_LS_multSPD);
  
 
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fPtHFEMCtot = new TH1F("fPtHFEMCtot",";p_{t} (GeV/c)",300,0,15);
  fPtHFEMCtot->Sumw2();
  fOutputList->Add(fPtHFEMCtot);
  
  fPtHFEMCtot_SPD = new TH2F("fPtHFEMCtot_SPD",";p_{t} (GeV/c)",300,0,15,200,0,200);
  fPtHFEMCtot_SPD->Sumw2();
 // fOutputList->Add(fPtHFEMCtot_SPD);
  
  fPtHFEMCtot_V0M = new TH2F("fPtHFEMCtot_V0M",";p_{t} (GeV/c)",300,0,15,1000,0,1000);
  fPtHFEMCtot_V0M->Sumw2();
 // fOutputList->Add(fPtHFEMCtot_V0M);
  
  fPtHFEMC = new TH1F("fPtHFEMC",";p_{t} (GeV/c)",300,0,15);
  fPtHFEMC->Sumw2();
  fOutputList->Add(fPtHFEMC);
  
  fPtHFEMC_SPD = new TH2F("fPtHFEMC_SPD",";p_{t} (GeV/c)",300,0,15,200,0,200);
  fPtHFEMC_SPD->Sumw2();
 // fOutputList->Add(fPtHFEMC_SPD);
  
  fPtHFEMC_V0M = new TH2F("fPtHFEMC_V0M",";p_{t} (GeV/c)",300,0,15,1000,0,1000);
  fPtHFEMC_V0M->Sumw2();
  //fOutputList->Add(fPtHFEMC_V0M);

  fPtHFEMC_afterfilterbit= new TH1F("fPtHFEMC_afterfilterbit",";p_{t} (GeV/c)",300,0,15.);
  fPtHFEMC_afterfilterbit->Sumw2();
  fOutputList->Add(fPtHFEMC_afterfilterbit);
  
  fPtHFEMC_aftertrackcuts= new TH1F("fPtHFEMC_aftertrackcuts",";p_{t} (GeV/c)",300,0,15.);
  fPtHFEMC_aftertrackcuts->Sumw2();
  fOutputList->Add(fPtHFEMC_aftertrackcuts);
  
  fPtHFEMC_aftertrackcuts_SPD= new TH2F("fPtHFEMC_aftertrackcuts_SPD",";p_{t} (GeV/c)",300,0,15.,200,0,200);
  fPtHFEMC_aftertrackcuts_SPD->Sumw2();
  //fOutputList->Add(fPtHFEMC_aftertrackcuts_SPD);
  
  fPtHFEMC_aftertrackcuts_V0M= new TH2F("fPtHFEMC_aftertrackcuts_V0M",";p_{t} (GeV/c)",300,0,15.,1000,0,1000);
  fPtHFEMC_aftertrackcuts_V0M->Sumw2();
  //fOutputList->Add(fPtHFEMC_aftertrackcuts_V0M);
  
  fPtHFEMC_aftertracktofcuts= new TH1F("fPtHFEMC_aftertracktofcuts",";p_{t} (GeV/c)",300,0,15.);
  fPtHFEMC_aftertracktofcuts->Sumw2();
  fOutputList->Add(fPtHFEMC_aftertracktofcuts);
  
  fPtHFEMC_aftertracktofcuts_SPD= new TH2F("fPtHFEMC_aftertracktofcuts_SPD",";p_{t} (GeV/c)",300,0,15.,200,0,200);
  fPtHFEMC_aftertracktofcuts_SPD->Sumw2();
  //fOutputList->Add(fPtHFEMC_aftertracktofcuts_SPD);
  
  fPtHFEMC_aftertracktofcuts_V0M= new TH2F("fPtHFEMC_aftertracktofcuts_V0M",";p_{t} (GeV/c)",300,0,15.,1000,0,1000);
  fPtHFEMC_aftertracktofcuts_V0M->Sumw2();
  //fOutputList->Add(fPtHFEMC_aftertracktofcuts_V0M);
  
  pi0MC= new TH1F("pi0MC",";p_{t} (GeV/c)",300,0,15.);
  pi0MC->Sumw2();
  fOutputList->Add(pi0MC);
  
  etaMC= new TH1F("etaMC",";p_{t} (GeV/c)",300,0,15.);
  etaMC->Sumw2();
  fOutputList->Add(etaMC);
  
  gammaMC= new TH1F("gammaMC",";p_{t} (GeV/c)",300,0,15.);
  gammaMC->Sumw2();
  fOutputList->Add(gammaMC);
  
  
  pi0MC_1= new TH1F("pi0MC_1",";p_{t} (GeV/c)",300,0,15.);
  pi0MC_1->Sumw2();
  fOutputList->Add(pi0MC_1);
  
  etaMC_1= new TH1F("etaMC_1",";p_{t} (GeV/c)",300,0,15.);
  etaMC_1->Sumw2();
  fOutputList->Add(etaMC_1);
  
  gammaMC_1= new TH1F("gammaMC_1",";p_{t} (GeV/c)",300,0,15.);
  gammaMC_1->Sumw2();
  fOutputList->Add(gammaMC_1);
  
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  fPt_elec_phot= new TH1F("fPt_elec_phot", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_phot->Sumw2();
  fOutputList->Add(fPt_elec_phot);
  
  fPt_elec_phot1= new TH1F("fPt_elec_phot1", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_phot1->Sumw2();
  fOutputList->Add(fPt_elec_phot1);
  
  fpt_e_nonphot_MC= new TH1F("fpt_e_nonphot_MC", "; p_{T}(GeV/c); counts;",300,0,15.);
  fpt_e_nonphot_MC->Sumw2();
  fOutputList->Add(fpt_e_nonphot_MC);
  
  fElecNos = new TH1F("fElecNos","no of electrons",10,-0.5,9.5);
  fElecNos->GetXaxis()->SetBinLabel(1,"Total no. of Particles");
  fElecNos->GetXaxis()->SetBinLabel(2,"Inclusive Electrons");
  fElecNos->GetXaxis()->SetBinLabel(3,"Non-hfe Electrons");
  fOutputList->Add(fElecNos);
 
  fPT_elec_MCtrue = new TH1F("fPT_elec_MCtrue","Pt of e from MC_True",300,0.,15.);
  fPT_elec_MCtrue->Sumw2();
  fOutputList->Add(fPT_elec_MCtrue);
  
  fpt_incl_elec_MC= new TH1F("fpt_incl_elec_MC","Pt of inclusive e from MC",300,0.,15.);
  fpt_incl_elec_MC->Sumw2();
  fOutputList->Add(fpt_incl_elec_MC);
  
  fPi0_Pt_MC = new TH1F("fPi0_Pt_MC","Pt of e mother Pi0 from MC",300,0.,15.);
  fPi0_Pt_MC->Sumw2();
  fOutputList->Add(fPi0_Pt_MC);
  
  fEta_Pt_MC = new TH1F("fEta_Pt_MC","Pt of e mother Eta from MC",300,0.,15.);
  fEta_Pt_MC->Sumw2();
  fOutputList->Add(fEta_Pt_MC);
  
  fGamma_Pt_MC = new TH1F("fGamma_Pt_MC","Pt of e mother Gamma from MC",300,0.,15.);
  fGamma_Pt_MC->Sumw2();
  fOutputList->Add(fGamma_Pt_MC);
  
  fPte_ULS_MC = new TH1F("fPte_ULS_MC","Pt of ULS from MC",300,0.,15.);
  fPte_ULS_MC->Sumw2(); 
  fOutputList->Add(fPte_ULS_MC);
  
  fPte_LS_MC = new TH1F("fPte_LS_MC","Pt of LS from MC",300,0.,15.);
  fPte_LS_MC->Sumw2();
  fOutputList->Add(fPte_LS_MC);
  
  fPte_ULS_MC_multV0M = new TH2F("fPte_ULS_MC_multV0M","Pt of ULS from MC",300,0.,15.,1000,0.,1000); 
  fPte_ULS_MC_multV0M->Sumw2();
  //fOutputList->Add(fPte_ULS_MC_multV0M);
  
  fPte_LS_MC_multV0M = new TH2F("fPte_LS_MC_multV0M","Pt of LS from MC",300,0.,15.,1000,0.,1000);
  fPte_LS_MC_multV0M->Sumw2();
  //fOutputList->Add(fPte_LS_MC_multV0M);
  
  fPte_ULS_MC_multSPD = new TH2F("fPte_ULS_MC_multSPD","Pt of ULS from MC",300,0.,15.,200,0.,200.);
  fPte_ULS_MC_multSPD->Sumw2(); 
  //fOutputList->Add(fPte_ULS_MC_multSPD);
  
  fPte_LS_MC_multSPD = new TH2F("fPte_LS_MC_multSPD","Pt of LS from MC",300,0.,15.,200,0.,200.);
  fPte_LS_MC_multSPD->Sumw2();
  //fOutputList->Add(fPte_LS_MC_multSPD);
  
  fPte_ULS_invMss_MC = new TH1F("fPte_ULS_invMss_MC","Pt of ULS invmass from MC",300,0.,15.);
  fPte_ULS_invMss_MC->Sumw2();
  fOutputList->Add(fPte_ULS_invMss_MC);
  
  fPte_LS_invMss_MC = new TH1F("fPte_LS_invMss_MC","Pt of LS invmss from MC",300,0.,15.);
  fPte_LS_invMss_MC->Sumw2();
  fOutputList->Add(fPte_LS_invMss_MC);
  
  fInvMULS = new TH2F("fInvMULS","Pt and InvMass from MC",300,0.,15.,1000,0,1.0);
  fInvMULS->Sumw2();
  fOutputList->Add(fInvMULS);
  
 /* fInvMULS_multV0M = new TH3F("fInvMULS_multV0M","Pt and InvMass from MC",300,0.,15.,100,0,1.0,1000,0.,1000.);
  fInvMULS_multV0M->Sumw2();
  fOutputList->Add(fInvMULS_multV0M);
  
  fInvMULS_multSPD = new TH3F("fInvMULS_multSPD","Pt and InvMass from MC",300,0.,15.,100,0,1.0,200,0.,200.);
  fInvMULS_multSPD->Sumw2();
  fOutputList->Add(fInvMULS_multSPD);*/
  
  
  fPt_elec_phot_multV0M= new TH2F("fPt_elec_phot_multV0M", "; p_{T}(GeV/c); counts;",300,0,15.,1000,0,1000.);
  fPt_elec_phot_multV0M->Sumw2();
 // fOutputList->Add(fPt_elec_phot_multV0M);
  
  fPt_elec_phot_multSPD= new TH2F("fPt_elec_phot_multSPD", "; p_{T}(GeV/c); counts;",300,0,15.,200,0,200.);
  fPt_elec_phot_multSPD->Sumw2();
///  fOutputList->Add(fPt_elec_phot_multSPD);
  
  fPt_elec_phot1_multV0M= new TH2F("fPt_elec_phot1_multV0M", "; p_{T}(GeV/c); counts;",300,0,15.,1000,0,1000.);
  fPt_elec_phot1_multV0M->Sumw2();
 // fOutputList->Add(fPt_elec_phot1_multV0M);
  
  fPt_elec_phot1_multSPD= new TH2F("fPt_elec_phot1_multSPD", "; p_{T}(GeV/c); counts;",300,0,15.,200,0,200.);
  fPt_elec_phot1_multSPD->Sumw2();
 // fOutputList->Add(fPt_elec_phot1_multSPD);
  
  //----pi0---------------
  fPt_elec_pi0_MB= new TH1F("fPt_elec_pi0_MB", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_pi0_MB->Sumw2();
  fOutputList->Add(fPt_elec_pi0_MB);
  
  fPt_elec_pi0_MB_LS= new TH1F("fPt_elec_pi0_MB_LS", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_pi0_MB_LS->Sumw2();
  fOutputList->Add(fPt_elec_pi0_MB_LS);
  
  fPt_elec_pi0_MB_ULS= new TH1F("fPt_elec_pi0_MB_ULS", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_pi0_MB_ULS->Sumw2();
  fOutputList->Add(fPt_elec_pi0_MB_ULS);
  
  //----eta---------------
  fPt_elec_eta_MB= new TH1F("fPt_elec_eta_MB", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_eta_MB->Sumw2();
  fOutputList->Add(fPt_elec_eta_MB);
  
  fPt_elec_eta_MB_LS= new TH1F("fPt_elec_eta_MB_LS", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_eta_MB_LS->Sumw2();
  fOutputList->Add(fPt_elec_eta_MB_LS);
  
  fPt_elec_eta_MB_ULS= new TH1F("fPt_elec_eta_MB_ULS", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_eta_MB_ULS->Sumw2();
  fOutputList->Add(fPt_elec_eta_MB_ULS);
  
  
  //----gamma---------------
  fPt_elec_gamma_MB= new TH1F("fPt_elec_gamma_MB", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_gamma_MB->Sumw2();
  fOutputList->Add(fPt_elec_gamma_MB);
  
  fPt_elec_gamma_MB_LS= new TH1F("fPt_elec_gamma_MB_LS", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_gamma_MB_LS->Sumw2();
  fOutputList->Add(fPt_elec_gamma_MB_LS);
  
  fPt_elec_gamma_MB_ULS= new TH1F("fPt_elec_gamma_MB_ULS", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_gamma_MB_ULS->Sumw2();
  fOutputList->Add(fPt_elec_gamma_MB_ULS); 
  
  //-------------------------------------------------------------------------------------------------------------
  
  Int_t nbinsInvMULS[4] = {100, 6, 100, 200};
  Double_t binlowInvMULS[4] = {0., 0, 0., 0.};
  Double_t binupInvMULS[4] = {10, 6, 1., 200};

	fInvMULSnSp = new THnSparseF("fInvMULSnSp", "fInvMULSnSp;pt;source;mass;multdep", 4, nbinsInvMULS,  binlowInvMULS,binupInvMULS);
	fInvMULSnSp->Sumw2();
	fOutputList->Add(fInvMULSnSp);
	
  fInvMULSpi0 = new TH2F("fInvMULSpi0","Pt and InvMass from MC",300,0.,15.,100,0,1.0);
  fInvMULSpi0->Sumw2();
  fOutputList->Add(fInvMULSpi0);
  
  fInvMULSeta = new TH2F("fInvMULSeta","Pt and InvMass from MC",300,0.,15.,100,0,1.0);
  fInvMULSeta->Sumw2();
  fOutputList->Add(fInvMULSeta);
  
  fInvMULSgamma = new TH2F("fInvMULSgamma","Pt and InvMass from MC",300,0.,15.,100,0,1.0);
  fInvMULSgamma->Sumw2();
  fOutputList->Add(fInvMULSgamma);
  
  

  //--------------------------------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------------------------------------
  fPtHFEMC_reco = new TH1F("fPtHFEMC_reco",";p_{t} (GeV/c)",300,0,15);
  fPtHFEMC_reco->Sumw2();
  fOutputList->Add(fPtHFEMC_reco);
  
   fPtHFEMC_reco_true = new TH1F("fPtHFEMC_reco_true",";p_{t} (GeV/c)",300,0,15);
  fPtHFEMC_reco_true->Sumw2();
  fOutputList->Add(fPtHFEMC_reco_true);
  
  fPtHFEMC_reco_SPD = new TH2F("fPtHFEMC_reco_SPD",";p_{t} (GeV/c)",300,0,15,200,0,200);
  fPtHFEMC_reco_SPD->Sumw2();
//  fOutputList->Add(fPtHFEMC_reco_SPD);
  
  fPtHFEMC_reco_V0M = new TH2F("fPtHFEMC_reco_V0M",";p_{t} (GeV/c)",300,0,15,1000,0,1000);
  fPtHFEMC_reco_V0M->Sumw2();
//  fOutputList->Add(fPtHFEMC_reco_V0M);
 
 //--------------------------------------------------------------------------------------------------
 //--------SPD correction---------------------------------
  
  const char *estimatorName="tr"; //SPD
  
  fHistNtrVsZvtxNoEtaCut = new TH2F("hNtrVsZvtxNoEtaCut",Form("N%s vs VtxZ; Z_{vtx};N_{%s};",estimatorName,estimatorName),300,-15,15,1000,0,1000); //
  //fOutputList->Add(fHistNtrVsZvtxNoEtaCut);
  
  fHistNtrCorrVsZvtxNoEtaCut_min = new TH2F("hNtrCorrVsZvtxNoEtaCut_min",Form("N%s vs VtxZ; Z_{vtx};N_{%s};",estimatorName,estimatorName),300,-15,15,1000,0,1000); //
  //fOutputList->Add(fHistNtrCorrVsZvtxNoEtaCut_min);
  
  fHistNtrCorrVsZvtxNoEtaCut_max = new TH2F("hNtrCorrVsZvtxNoEtaCut_max",Form("N%s vs VtxZ; Z_{vtx};N_{%s};",estimatorName,estimatorName),300,-15,15,1000,0,1000); //
  //fOutputList->Add(fHistNtrCorrVsZvtxNoEtaCut_max);
  
  fSPDCorrMultDistNoEtaCut_min = new TH1F("fSPDCorrMultDistNoEtaCut_min",Form("Corrected Mult Dist; N_{%s};Counts;",estimatorName),1000,0,1000); //
  //fOutputList->Add(fSPDCorrMultDistNoEtaCut_min);
  
 	fSPDCorrMultDistNoEtaCut_max = new TH1F("fSPDCorrMultDistNoEtaCut_max",Form("Corrected Mult Dist;  N_{%s};Counts;",estimatorName),1000,0,1000); //
  //fOutputList->Add(fSPDCorrMultDistNoEtaCut_max);
  
  fSPD_tracklet_NoEtaCut = new TH1F("fSPD_tracklet_NoEtaCut", "fSPD_tracklet_NoEtaCut Distribution",300,0,300);
  fSPD_tracklet_NoEtaCut->Sumw2(); 
  //fOutputList->Add(fSPD_tracklet_NoEtaCut); 
  
  fSPD_tracklet = new TH1F("fSPD_tracklet", "fSPD_tracklet Distribution",300,0,300);
  fSPD_tracklet->Sumw2(); 
  //fOutputList->Add(fSPD_tracklet); 
  
  fSPDCorrMultDist_min = new TH1F("fSPDCorrMultDist_min",Form("Corrected Mult Dist;  N_{%s};Counts;",estimatorName),300,0,300); //
  //fOutputList->Add(fSPDCorrMultDist_min);
  
 	fSPDCorrMultDist_max = new TH1F("fSPDCorrMultDist_max",Form("Corrected Mult Dist; N_{%s};Counts;",estimatorName),300,0,300); //
  //fOutputList->Add(fSPDCorrMultDist_max);
  
  fHistNtrVsZvtx = new TH2F("hNtrVsZvtx",Form("N%s vs VtxZ; Z_{vtx};N_{%s};",estimatorName,estimatorName),300,-15,15,300,0,300); //
  //fOutputList->Add(fHistNtrVsZvtx);
  
  fHistNtrCorrVsZvtx_min = new TH2F("hNtrCorrVsZvtx_min",Form("N%s vs VtxZ; Z_{vtx};N_{%s};",estimatorName,estimatorName),300,-15,15,300,0,300); //
  //fOutputList->Add(fHistNtrCorrVsZvtx_min);
  					
  fHistNtrCorrVsZvtx_max = new TH2F("hNtrCorrVsZvtx_max",Form("N%s vs VtxZ; Z_{vtx};N_{%s};",estimatorName,estimatorName),300,-15,15,300,0,300); //
 // fOutputList->Add(fHistNtrCorrVsZvtx_max);
   
  fV0MultVsZvtx= new TH2F("fV0MultVsZvtx","UnCorrected Mult Dist; N_VOM;Counts;",300,-15,15,1000,0,1000); //
  //fOutputList->Add(fV0MultVsZvtx);
  
  fV0MultCorrVsZvtx= new TH2F("fV0MultCorrVsZvtx","Corrected Mult Dist; N_VOM;Counts;",300,-15,15,1000,0,1000); //
  //fOutputList->Add(fV0MultCorrVsZvtx);	
 
  fNchVsZvtx= new TH2F("fNchVsZvtx","Corrected Mult Dist;Nch;Counts;",300,-15,15,1000,0,1000); //
 /// fOutputList->Add(fNchVsZvtx);
  
  fSPDNtrCorrVsNch = new TH2F("fSPDNtrCorrVsNch","",200,0,200.,400,0,400.); //
  fSPDNtrCorrVsNch->Sumw2();
 // fOutputList->Add(fSPDNtrCorrVsNch);	
  
	fV0MCorrVsNch= new TH2F("fV0MCorrVsNch","",1000,0,1000.,400,0,400.); //
	fV0MCorrVsNch->Sumw2();
  //fOutputList->Add(fV0MCorrVsNch);
  
    
  fSPDtr_V0Mmult= new TH2F("fSPDtr_V0Mmult","",200,0,200.,1000,0,1000);
	fSPDtr_V0Mmult->Sumw2();
  //fOutputList->Add(fSPDtr_V0Mmult);
  
  //fSPDNtrCorrVsV0MCorrVsNch= new TH3F("fSPDNtrCorrVsV0MCorrVsNch","SPDtr V0M N_{ch};SPD_{tr}; N_{V0M};Counts;",100,0,200,500,0,1000,150,0,300.); //
  /*fSPDNtrCorrVsV0MCorrVsNch= new TH3F("fSPDNtrCorrVsV0MCorrVsNch","SPDtr V0M N_{ch};SPD_{tr}; N_{V0M};Counts;",50,0,200,500,0,500,50,0,300.);
  fSPDNtrCorrVsV0MCorrVsNch->Sumw2();
  fOutputList->Add(fSPDNtrCorrVsV0MCorrVsNch);
  */
  /*
  fnSigmaVsPt_TPC_cut_multSPD= new TH3F("fnSigmaVsPt_TPC_cut_multSPD","",300,0,15.,250,-15,10,200,0,200.);
	fnSigmaVsPt_TPC_cut_multSPD->Sumw2();
  fOutputList->Add(fnSigmaVsPt_TPC_cut_multSPD);*/
  	
	// ----- weights for tagging efficiency -----
    // QA weights
    Int_t nBinspdg = 3;
    Double_t minpdg = 0.;
    Double_t maxpdg = 3.;
    Double_t binLimpdg[nBinspdg+1];
    for(Int_t i=0; i<=nBinspdg; i++) binLimpdg[i]=(Double_t)minpdg + (maxpdg-minpdg)/nBinspdg*(Double_t)i ;
    
    Int_t nBinsg = 2;
    Double_t ming = 0.;
    Double_t maxg = 2.;
    Double_t binLimg[nBinsg+1];
    for(Int_t i=0; i<=nBinsg; i++) binLimg[i]=(Double_t)ming + (maxg-ming)/nBinsg*(Double_t)i ;
    
    Int_t nBinstype = 7;
    Double_t mintype = -1.;
    Double_t maxtype = 6.;
    Double_t binLimtype[nBinstype+1];
    for(Int_t i=0; i<=nBinstype; i++) binLimtype[i]=(Double_t)mintype + (maxtype-mintype)/nBinstype*(Double_t)i ;
    
    Int_t nBinspt = 44;
    Double_t binLimpt[45]= {0.1,0.112797,0.127231,0.143512,0.161877,0.182592,0.205957,0.232313,0.262041,0.295573,0.333397,
0.37606,0.424183,0.478465,0.539692,0.608754,0.686654,0.774523,0.873636,0.985432,1.11153,1.25377,1.41421,1.59519,1.79932,2.02957,
2.28928,2.58223,2.91267,3.2854,3.70582,4.18004,4.71494,5.3183,5.99886,6.76651,7.6324,8.60909,9.71076,10.9534,12.3551,13.9361,15.7195,17.731,20};//bin limits from the measured pi0 spectrum
    //for(Int_t i=0; i<=nBinspt; i++) binLimpt[i]=(Double_t)minpt + (maxpt-minpt)/nBinspt*(Double_t)i ;
    
    const Int_t nDima=4;
    Int_t nBina[nDima] = {44,nBinspdg,nBinsg,nBinstype};
    fPi0EtaSpectra = new THnSparseF("fPi0EtaSpectra","fPi0EtaSpectra",nDima,nBina);
    fPi0EtaSpectra->SetBinEdges(0,binLimpt);
    fPi0EtaSpectra->SetBinEdges(1,binLimpdg);
    fPi0EtaSpectra->SetBinEdges(2,binLimg);
    fPi0EtaSpectra->SetBinEdges(3,binLimtype);
    fPi0EtaSpectra->Sumw2();
    //fOutputList->Add(fPi0EtaSpectra); ---> put at the end of the list
    // ------------------------------------------
    fOutputList->Add(fPi0EtaSpectra);
  
  
  
  
  Int_t nbinspt[3] = {1000, 3, 7};
  Double_t binlow[3] = {0., 0, -1.};
  Double_t binup[3] = {10, 3, 6.};
	fPi0EtaSpectraSp = new THnSparseF("fPi0EtaSpectraSp", "fPi0EtaSpectraSp;pt;source;type", 3, nbinspt, binlow,binup);
	fPi0EtaSpectraSp->Sumw2();
   fOutputList->Add(fPi0EtaSpectraSp);
   
 //===================================================================================================//
 
 
  fNentries2=new TH1F("CutSet", "", 20,-0.5,19.5);
  fNentries2->GetXaxis()->SetBinLabel(1,"trigger");
  fNentries2->GetXaxis()->SetBinLabel(2,"TPCNclus");
  fNentries2->GetXaxis()->SetBinLabel(3,"ITSNclus");
  fNentries2->GetXaxis()->SetBinLabel(4,"TPCNclusPID");
  fNentries2->GetXaxis()->SetBinLabel(5,"SPDBoth");
  fNentries2->GetXaxis()->SetBinLabel(6,"SPDAny");
  fNentries2->GetXaxis()->SetBinLabel(7,"SPDFirst");
  fNentries2->GetXaxis()->SetBinLabel(8,"DCAxyCut");
  fNentries2->GetXaxis()->SetBinLabel(9,"DCAzCut");
  fNentries2->GetXaxis()->SetBinLabel(10,"Etarange");  
  fNentries2->GetXaxis()->SetBinLabel(11,"TPCnsigmin");
  fNentries2->GetXaxis()->SetBinLabel(12,"TPCnsigmax");
  fNentries2->GetXaxis()->SetBinLabel(13,"TOFnsig");
  fNentries2->GetXaxis()->SetBinLabel(14,"InvmassCut");
  fNentries2->GetXaxis()->SetBinLabel(15,"AssoTPCCluster");  
  fNentries2->GetXaxis()->SetBinLabel(16,"AssoITSRefit");
  fNentries2->GetXaxis()->SetBinLabel(17,"AssopTMin");
  fNentries2->GetXaxis()->SetBinLabel(18,"AssoEtarange");
  fNentries2->GetXaxis()->SetBinLabel(19,"AssoTPCnsig");
  fNentries2->GetXaxis()->SetBinLabel(19,"Filterbit");
   fOutputList->Add(fNentries2);
  
  fNentries2->SetBinContent(1, ftrigger);
  fNentries2->SetBinContent(2,fTPCNclus);
  fNentries2->SetBinContent(3,fITSNclus);
  fNentries2->SetBinContent(4,fTPCNclusPID);
  fNentries2->SetBinContent(5,fSPDBoth);
  fNentries2->SetBinContent(6,fSPDAny);
  fNentries2->SetBinContent(7,fSPDFirst);
  fNentries2->SetBinContent(8,fDCAxyCut);
  fNentries2->SetBinContent(9,fDCAzCut);
  fNentries2->SetBinContent(10,fEtarange);  
  fNentries2->SetBinContent(11,fTPCnsigmin);
  fNentries2->SetBinContent(12,fTPCnsigmax);
  fNentries2->SetBinContent(13,fTOFnsig);
  fNentries2->SetBinContent(14,fInvmassCut);
  fNentries2->SetBinContent(15,fAssoTPCCluster);  
  fNentries2->SetBinContent(16,fAssoITSRefit);
  fNentries2->SetBinContent(17,fAssopTMin);
  fNentries2->SetBinContent(18,fAssoEtarange);
  fNentries2->SetBinContent(19,fAssoTPCnsig);
  fNentries2->SetBinContent(20,fFilterbit);

   
  //==================================================================================================//
  //==================================================================================================//
  					
  const char* nameoutput=GetOutputSlot(2)->GetContainer()->GetName();  // gets the output slot from mgr-> in AddTask
  fNentries=new TH1F(nameoutput, "", 5,-0.5,4.5);
  fNentries->GetXaxis()->SetBinLabel(1,"No. of Events Analyzed");
  fNentries->GetXaxis()->SetBinLabel(2,"No. of Events Accepted");
  fNentries->GetXaxis()->SetBinLabel(3,"After Pile Up cuts");
  fNentries->GetXaxis()->SetBinLabel(4,"After Z vtx cut");
  fNentries->GetXaxis()->SetBinLabel(5,"After N contributors");
  //fNentries->GetXaxis()->SetBinLabel(6," ");
  //fNentries->GetXaxis()->SetBinLabel(6,"Pile-up Rej");
  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);
  fNentries->Sumw2();
  fNentries->SetMinimum(0);
  
   
 // const char* nameoutput2=GetOutputSlot(3)->GetContainer()->GetName();  // gets the output slot from mgr-> in AddTask
  
 
  
  fNentries2->GetXaxis()->SetNdivisions(1,kFALSE);
  fNentries2->Sumw2();
  fNentries2->SetMinimum(0);
 
  PostData(1,fOutputList);
  PostData(2,fNentries);

}

//////////////=================================================///////////////////////

void AliAnalysisTaskQAHFE::UserExec(Option_t *) 
{
	fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
	if (!fAOD) { return;}
  
	fNentries->Fill(0);                //No. of Events Analyzed 
	Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
	//Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kHighMultV0);
	//Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & ftrigger);
  
	if(!isSelected) return;  
	fNentries->Fill(1);                 //No. of Events Accepted
	
	/*//----------------------------For 15f --------------------------------------------------------------
	TString firedTriggerClasses = static_cast<const AliAODEvent*>(InputEvent())->GetFiredTriggerClasses();
	if (!(firedTriggerClasses.Contains("ALLNOTRD"))) return; // events without SDD
	fNentries->Fill(2);                //After SDD cut
	//-----------------------------------------------------------------------------------------------------*/
	
	const AliAODVertex *pVtx= fAOD->GetPrimaryVertex();
	if(!pVtx){return;}
	TString vtxTtl = pVtx->GetTitle();
	//if(!vtxTtl.Contains("VertexerTracks")) {printf("No PrimaryVertexTracks\n");}
	
	Double_t Zvertex=-999,Xvertex=-999,Yvertex=-999;
  
	Zvertex = pVtx->GetZ();
	Yvertex = pVtx->GetY();
	Xvertex = pVtx->GetX();
	fVtxZ->Fill(Zvertex);
	fVtxX->Fill(Xvertex);
	fVtxY->Fill(Yvertex);
  	
	if(pVtx->GetNContributors()<2) return;
	
	fNentries->Fill(2);
	
	Bool_t isPileupfromSPDmulbins=fAOD->IsPileupFromSPDInMultBins(); //This function checks if there was a pile up reconstructed with SPD
	
	if(isPileupfromSPDmulbins) return;
	
	Int_t minContributors=5;    //minimum contributors to the pilup vertices, multi-vertex
	Double_t minChi2=5.; 
	Double_t minWeiZDiff=15;   //minimum of the sqrt of weighted distance between the primary and the pilup vertex, multi-vertex
	Bool_t checkPlpFromDifferentBC=kFALSE;
	
	AliAnalysisUtils utils;
	utils.SetMinPlpContribMV(minContributors); //Multi Vertex pileup selection
	utils.SetMaxPlpChi2MV(minChi2);   //max value of Chi2perNDF of the pileup vertex, multi-vertex
	utils.SetMinWDistMV(minWeiZDiff);
	utils.SetCheckPlpFromDifferentBCMV(checkPlpFromDifferentBC); //SPD Pileup slection
	Bool_t isPileupFromMV = utils.IsPileUpMV(fAOD);      //check for multi-vertexer pile-up
	
	if(isPileupFromMV) return;
	fNentries->Fill(3);
/*
	//For Pb-Pb ---------------
   ///Correlation cuts between vertexes:
	AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
	if (pVtx->GetNContributors()<2 || vtSPD->GetNContributors()<1) return; // one of vertices is missing

	double covTrc[6],covSPD[6];
	pVtx->GetCovarianceMatrix(covTrc);
	vtSPD->GetCovarianceMatrix(covSPD);
	double dz =pVtx->GetZ()-vtSPD->GetZ();
	double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
	double errTrc = TMath::Sqrt(covTrc[5]);
	double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
	if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return; // bad vertexing 
	fNentries->Fill(3);
	///Number of events rejected by the correlation cuts between SPD and track vertexes
	///Track vertex cut:
*/		
	if(TMath::Abs(Zvertex)>10){return;} 
	fNentries->Fill(4);
  
	fVtxZ_corr->Fill(Zvertex);
	
	Int_t runNumber = fAOD->GetRunNumber();
	Int_t num=0;
	Int_t nTracks = fAOD->GetNumberOfTracks();

//Look for kink mother for AOD
	Double_t *fListOfmotherkink = 0;
	Int_t fNumberOfVertices = 0; 
	Int_t fNumberOfMotherkink = 0;

	fNumberOfVertices = fAOD->GetNumberOfVertices();		
	fListOfmotherkink = new Double_t[fNumberOfVertices];
			
	for(Int_t ivertex=0; ivertex < fNumberOfVertices; ivertex++) 
	{
		AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
		if(!aodvertex) continue;
		if(aodvertex->GetType()==AliAODVertex::kKink) 
		{
			AliAODTrack *mother1 = (AliAODTrack *) aodvertex->GetParent();
			if(!mother1) continue;
			Int_t idmother = mother1->GetID();
			fListOfmotherkink[fNumberOfMotherkink] = idmother;
			fNumberOfMotherkink++;
		}
	}

	
	AliMultSelection *fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
	//if(!fMultSelection-> IsEventSelected()){return;}
	if(!fMultSelection) 
	{
		//If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
		AliWarning("AliMultSelection object not found!");
	}
	else{
		fMultiPercentile = fMultSelection->GetMultiplicityPercentile("V0M", false); //or use "SPDTracklets" V0M, i used CL0 also
		fMultiPercentileSPD = fMultSelection->GetMultiplicityPercentile("SPDTracklets", false);
	}

	fCent->Fill(fMultiPercentile); //centrality dist.
	fCentSPD->Fill(fMultiPercentileSPD);

	//cout<<" fMultiPercentileSPD "<<fMultiPercentileSPD<<" fMultiPercentile "<< fMultiPercentile <<endl;
	
	fpidResponse = fInputHandler->GetPIDResponse();
	if(!fpidResponse)
	{
		AliDebug(1, "Using default PID Response");
		fpidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
		fPID->SetPIDResponse(fpidResponse);
	}
	fPID->SetPIDResponse(fpidResponse);
	
	//===============V0 Multiplicity===========================
	AliAODVZERO *vzeroAOD = dynamic_cast<AliAODVZERO *>( dynamic_cast<AliAODEvent *>(fAOD)->GetVZEROData());
	Int_t V0AMult = static_cast<Int_t>(vzeroAOD->GetMTotV0A());
	Int_t V0CMult = static_cast<Int_t>(vzeroAOD->GetMTotV0C());
	fV0AMult->Fill(V0AMult);
	fV0CMult->Fill(V0CMult);
	Int_t V0Mult=V0AMult+V0CMult;
  
	fV0Mult->Fill(V0Mult);
	fV0MultVsZvtx->Fill(Zvertex,V0Mult);
	//V0M Correction
  Int_t V0AMultCorr=V0AMult, V0CMultCorr=V0CMult, V0MultCorr=V0Mult;
  V0AMultCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(V0AMult,Zvertex));
  V0CMultCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0C(V0CMult,Zvertex));
 	 V0MultCorr = V0AMultCorr + V0CMultCorr;

 	fV0MultCorr->Fill(V0MultCorr);
 	fV0MultCorrVsZvtx->Fill(Zvertex,V0MultCorr);
	Double_t V0Mmult = V0MultCorr;
 
	//===============SPD Tracklet==============================
	Int_t nAcceta=0.;
	AliAODTracklets *tracklets = ((AliAODEvent*)fAOD)->GetTracklets();
	Int_t nTracklets = tracklets->GetNumberOfTracklets();
	for (Int_t nn = 0; nn < nTracklets; nn++)
	{
		Double_t theta = tracklets->GetTheta(nn);
		Double_t eta = -TMath::Log(TMath::Tan(theta/2.0));
		if (TMath::Abs(eta) < 1.) nAcceta++;
	}
	fSPD_tracklet_NoEtaCut->Fill(nTracklets); 
	fSPD_tracklet->Fill(nAcceta);
 	

  	
	Double_t counteta1Corr_min=nAcceta;
	Double_t N_corr_tr_eta1=nAcceta;
  
  
 	// Double_t fRefMult_min=14.227,fRefMult_max=17.633, fRefMult1_min=9.003, fRefMult1_max=13.077; // June15
 	//Double_t fRefMult_min=14.276,fRefMult_max=17.632,
 	//Double_t fRefMult1_max=13.0793; // June22
  Double_t fRefMult1_min=9.016;
  
	//Double_t countCorr_min=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(fMultEstimatorAvg[0],nTracklets,Zvertex,fRefMult_min));
	//Double_t countCorr_max=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(fMultEstimatorAvg[0],nTracklets,Zvertex,fRefMult_max));
	//fSPDCorrMultDistNoEtaCut_min->Fill(countCorr_min);
  //fSPDCorrMultDistNoEtaCut_max->Fill(countCorr_max);	
  //fHistNtrVsZvtxNoEtaCut->Fill(Zvertex,nTracklets);
  //fHistNtrCorrVsZvtxNoEtaCut_min->Fill(Zvertex,countCorr_min);
  //fHistNtrCorrVsZvtxNoEtaCut_max->Fill(Zvertex,countCorr_max);
  
  
	counteta1Corr_min=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(fMultEstimatorAvg[1],nAcceta,Zvertex,fRefMult1_min));
	N_corr_tr_eta1=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(fMultEstimatorAvg[1],nAcceta,Zvertex,fRefMult)); //48.75 fRefMult
	
  Double_t SPDntr = N_corr_tr_eta1;

  fSPDCorrMultDist_min->Fill(counteta1Corr_min);
  fSPDCorrMultDist_max->Fill(SPDntr);
  	
  fHistNtrVsZvtx->Fill(Zvertex,nAcceta);
  fHistNtrCorrVsZvtx_min->Fill(Zvertex,counteta1Corr_min);
  fHistNtrCorrVsZvtx_max->Fill(Zvertex,SPDntr);
  	
 	fMultV0M_vs_alimult->Fill(fMultiPercentile,V0Mmult);
	fMultSPD_vs_alimult->Fill(fMultiPercentileSPD,SPDntr);
	
	fSPDtr_V0Mmult->Fill(SPDntr,V0Mmult);
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////========================================MC LOOP===================================================/////
//////==================================================================================================/////  
	if(fIsMC) 
	{    
		Double_t qaweights[5];
		Double_t pi0etaweights[3];
		
		fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
		if(!fMCarray)
		{
			AliError("Array of MC particles not found");
		 	return;
		}
		for(Int_t iMC = 0; iMC < fMCarray->GetEntries(); iMC++)
		{

			fMCparticle = (AliAODMCParticle*) fMCarray->At(iMC);
			pdg = fMCparticle->GetPdgCode();
			
			if(TMath::Abs(pdg) == 11){
				Int_t IsElecHf=GetHFE(fMCparticle,fMCarray);
				if((IsElecHf==kBeauty) || (IsElecHf==kCharm)){
			 		fPtHFEMCtot->Fill(fMCparticle->Pt());
			 		fPtHFEMCtot_SPD->Fill(fMCparticle->Pt(),SPDntr);
					fPtHFEMCtot_V0M->Fill(fMCparticle->Pt(),V0Mmult);
			 	}
			}
			///For the reconstruction efficiency of HFE:-----------
			if(TMath::Abs(fMCparticle->Eta()) <=fEtarange  )
			{					
				if(TMath::Abs(pdg) == 11)
				{		
					Int_t IsElecHf=GetHFE(fMCparticle,fMCarray);
					if((IsElecHf==kBeauty) || (IsElecHf==kCharm)){
						fPtHFEMC->Fill(fMCparticle->Pt());
						fPtHFEMC_SPD->Fill(fMCparticle->Pt(),SPDntr);
						fPtHFEMC_V0M->Fill(fMCparticle->Pt(),V0Mmult);
					}
				}
			} //eta <0.8 condition
			
			//pt spectra for pi0 and eta
			if(TMath::Abs(fMCparticle->Y()>0.5 )) continue;
			
			if (TMath::Abs(fMCparticle->GetPdgCode())==111) pi0MC->Fill(fMCparticle->Pt());  // pdg=111 pi0
			if (TMath::Abs(fMCparticle->GetPdgCode())==221) etaMC->Fill(fMCparticle->Pt());  // pdg=221 eta	
			if (TMath::Abs(fMCparticle->GetPdgCode())==22) gammaMC->Fill(fMCparticle->Pt());   // pdg=22  gamma
			
			Int_t type = GetPi0EtaType(fMCparticle,fMCarray);
				
			
			
			///Using thnSparse--------------------------------------
				//Pt
				qaweights[0] = fMCparticle->Pt();
				pi0etaweights[0] = fMCparticle->Pt();
						
				// What pdg
				qaweights[1]=-1.;
				if (TMath::Abs(fMCparticle->GetPdgCode())==111) qaweights[1]=0.2;  // pdg=111 pi0
				if (TMath::Abs(fMCparticle->GetPdgCode())==221) qaweights[1]=1.2;  // pdg=221 eta
				if (TMath::Abs(fMCparticle->GetPdgCode())==22) qaweights[1]=2.2;   // pdg=22  gamma
				
				pi0etaweights[1]=qaweights[1];
				
				// What type
				//Int_t type = GetPi0EtaType(fMCparticle,fMCarray);
				qaweights[3]=type;
				pi0etaweights[2]=type;
					
				// Fill
				if(qaweights[1]>0.) fPi0EtaSpectra->Fill(qaweights);
				if(pi0etaweights[1]>0.) fPi0EtaSpectraSp->Fill(pi0etaweights);
				
				if(type==kNoFeedDown)
				{
						if (TMath::Abs(fMCparticle->GetPdgCode())==111) pi0MC_1->Fill(fMCparticle->Pt());  // pdg=111 pi0
						if (TMath::Abs(fMCparticle->GetPdgCode())==221) etaMC_1->Fill(fMCparticle->Pt());  // pdg=221 eta	
						if (TMath::Abs(fMCparticle->GetPdgCode())==22) gammaMC_1->Fill(fMCparticle->Pt());   // pdg=22  gamma
				
				}
				///-----------------------------------------------------
		
		       			
					
					
		////====================================================
		} //For loop
	}



	if(fIsMC) 
	{ 
		Int_t nch = GetNcharged();
		fSPDNtrCorrVsNch->Fill(SPDntr,nch);
		fV0MCorrVsNch->Fill(V0Mmult,nch);
		fNchVsZvtx->Fill(Zvertex,nch);
		//fSPDNtrCorrVsV0MCorrVsNch->Fill(SPDntr,V0Mmult,nch);
	}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////======================================TRACK LOOP==================================================/////
//////==================================================================================================/////

	for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) 
  {
		AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
   	
    if(!track) AliFatal("Not a standard AOD");
  			Double_t pt=0., p=0., eta=-999, dEdx_TPC=-999., Signal_TOF=-999., TOFbeta=-999., fTPCnSigma=-999.0, fTOFnSigma=-999.0;    
    
		pt = track->Pt();
		p = track->P();  
		eta=track->Eta();
	
	if (TMath::Abs(eta)>=fEtarange) continue; 
	if(!track->TestFilterMask(fFilterbit)) continue; 
   
	if(fIsMC && track->GetLabel()>=0)
   		{
   			fMCparticle=(AliAODMCParticle*)fMCarray->At(track->GetLabel());
				pdg = fMCparticle->GetPdgCode();
				if(TMath::Abs(pdg) == 11)
				{		
					Int_t IsElecHf=GetHFE(fMCparticle,fMCarray);
					if((IsElecHf==kBeauty) || (IsElecHf==kCharm)){
					fPtHFEMC_afterfilterbit->Fill(track->Pt());
					}
				}
   		}
    Int_t tracktypeTrig=ClassifyTrack(track,pVtx);  //track classify
    
    fDCAZvsPt->Fill(track->Pt(),fDCAz);
	 fDCAXYvsPt->Fill(track->Pt(),fDCAxy);
    	
    if(tracktypeTrig!=1) continue;  //if(tracktype==0) continue; if(tracktype==1) //tracks "not" passed AliAODTrack::kPrimary at  	
     				  													//reconstructed level & have proper TPC PID response(?)
	  // Reject kink mother
 
		Bool_t kinkmotherpass = kTRUE;
		for(Int_t kinkmother = 0; kinkmother < fNumberOfMotherkink; kinkmother++) 
		{
			if(track->GetID() == fListOfmotherkink[kinkmother]) 
			{
				kinkmotherpass = kFALSE;
				continue;
			}
		}
		if(!kinkmotherpass) continue;
		
		if(fIsMC && track->GetLabel()>=0)
   		{
   		fMCparticle=(AliAODMCParticle*)fMCarray->At(track->GetLabel());
			pdg = fMCparticle->GetPdgCode();
			if(TMath::Abs(pdg) == 11)
			{		
				Int_t IsElecHf=GetHFE(fMCparticle,fMCarray);
				if((IsElecHf==kBeauty) || (IsElecHf==kCharm)){
					fPtHFEMC_aftertrackcuts->Fill(track->Pt());
					fPtHFEMC_aftertrackcuts_SPD->Fill(track->Pt(),SPDntr);
					fPtHFEMC_aftertrackcuts_V0M->Fill(track->Pt(),V0Mmult);
				}
			}
   	}
   	
		num++;
    	
		dEdx_TPC = track->GetTPCsignal();
     	
		TOFbeta=Beta(track);
		Double_t weight_lan=0.,weight_lan_inv=0.;
		Double_t weight_err=0.,weight_err_inv=0.;

		//Signal_TOF = track->GetTOFsignal();     // returns: TOF signal - t0 (T0 interaction time)
		//TOFbeta = fHelPID->TOFBetaCalc(track);

		
		fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
		fTOFnSigma = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
						
		fHistPt->Fill(pt);
		fTPCSignal->Fill(dEdx_TPC);

		fdEdxVsP_TPC->Fill(p,dEdx_TPC);
		fnSigmaVsP_TPC->Fill(p,fTPCnSigma);
		fnSigmaVsP_TOF->Fill(p,fTOFnSigma);
		fnBetaVsP_TOF->Fill(p,TOFbeta);
		
		fnSigmaVsEta_TPC->Fill(eta,fTPCnSigma);
	
		fdEdxVsPt_TPC->Fill(pt,dEdx_TPC);
		fnSigmaVsPt_TPC->Fill(pt,fTPCnSigma);
		fnSigmaVsPt_TOF->Fill(pt,fTOFnSigma);
		fnBetaVsPt_TOF->Fill(pt,TOFbeta);
		  
		if(TMath::Abs(fTOFnSigma) < fTOFnsig) 
		{
			fdEdxVsP_TPC_cut->Fill(p,dEdx_TPC);
			fnSigmaVsP_TPC_cut->Fill(p,fTPCnSigma);
			fnSigmaVsEta_TPC_cut->Fill(eta,fTPCnSigma);
       		
			//fnSigmaVsP_TPC_cut_multV0M->Fill(p,fTPCnSigma,V0Mmult);
			//fnSigmaVsP_TPC_cut_multSPD->Fill(p,fTPCnSigma,SPDntr);
       	
      if(fIsMC && track->GetLabel()>=0)
   		{
   			fMCparticle=(AliAODMCParticle*)fMCarray->At(track->GetLabel());
				pdg = fMCparticle->GetPdgCode();
				if(TMath::Abs(pdg) == 11)
				{		
					Int_t IsElecHf=GetHFE(fMCparticle,fMCarray);
					if((IsElecHf==kBeauty) || (IsElecHf==kCharm)){
					fPtHFEMC_aftertracktofcuts->Fill(pt);
					fPtHFEMC_aftertracktofcuts_SPD->Fill(pt,SPDntr);
					fPtHFEMC_aftertracktofcuts_V0M->Fill(pt,V0Mmult);
					}
				}
   		}
       	
       		
			weight_lan=fLandau->Eval(p);
			weight_err=fErr->Eval(p);
			weight_lan_inv=1-weight_lan;
			weight_err_inv=1-weight_err;
       		       		
			if(fTPCnSigma>fTPCnsigmin && fTPCnSigma<fTPCnsigmax) 
			{ 		
				fHadCot_Landau->Fill(pt,weight_lan);
				fHadCot_Err->Fill(pt,weight_err);
				fPt_incl_e_Landau->Fill(pt,weight_lan_inv);
				fPt_incl_e_Err->Fill(pt,weight_err_inv);
  						
				//---------Multiplicity----------------------
				fPteV0M_Lan->Fill(pt,V0Mmult,weight_lan_inv);
				fPteSPD_Lan->Fill(pt,SPDntr,weight_lan_inv);
				fPteV0M_Err->Fill(pt,V0Mmult,weight_err_inv);
				fPteSPD_Err->Fill(pt,SPDntr,weight_err_inv);
  				
				fPteV0M->Fill(pt,V0Mmult);
				fPteSPD->Fill(pt,SPDntr);
  				
				fHadV0MLan->Fill(pt,V0Mmult,weight_lan);
				fHadSPDLan->Fill(pt,SPDntr,weight_lan);
				fHadV0MErr->Fill(pt,V0Mmult,weight_err);
 				fHadSPDErr->Fill(pt,SPDntr,weight_err);	
  				
				//SelectPhotonicElectron(iTracks,track,fpidResponse,V0Mult,nAcceta,fIsMC);
  				
				fNonHFE = new AliSelectNonHFE();
				fNonHFE->SetAODanalysis(kTRUE);
				fNonHFE->SetInvariantMassCut(fInvmassCut);
				fNonHFE->SetAlgorithm("DCA"); //KF,DCA
				fNonHFE->SetPIDresponse(fpidResponse);
				fNonHFE->SetTrackCuts(-1*fAssoTPCnsig,fAssoTPCnsig); //TPCnsigma cuts
				//fNonHFE->SetChi2OverNDFCut(4.0);
				//fNonHFE->SetDCACut(fDCAcut);
				//fNonHFE-> SetDCAPartnerCuts(1,2);   //SetDCAPartnerCuts(Double_t xy, Double_t z)
				fNonHFE->SetAdditionalCuts(fAssopTMin,fAssoTPCCluster);  //SetAdditionalCuts(fPtMinAsso,fTpcNclsAsso);
				fNonHFE->SetHistMassBack(fInvmassLS1);
				fNonHFE->SetHistMass(fInvmassULS1);
				fNonHFE->FindNonHFE(iTracks,track,fAOD);
			
				Int_t fNULS = fNonHFE->GetNULS();
				Int_t fNLS = fNonHFE->GetNLS();
			 
				if(fNonHFE->IsULS())
				{
					fPte_ULS->Fill(track->Pt(),fNULS);
					fPte_ULS_multV0M->Fill(track->Pt(),V0Mmult,fNULS);
					fPte_ULS_multSPD->Fill(track->Pt(),SPDntr,fNULS);
				}
				if(fNonHFE->IsLS())
				{
					fPte_LS->Fill(track->Pt(),fNLS);
					fPte_LS_multV0M->Fill(track->Pt(),V0Mmult,fNLS);
					fPte_LS_multSPD->Fill(track->Pt(),SPDntr,fNLS);
				}
  
				if(fIsMC && track->GetLabel()>=0)
   			{
					//cout<<"track->GetLabel()  "<<track->GetLabel()<<endl;
					fElecNos->Fill(0);   						
					fMCparticle=(AliAODMCParticle*)fMCarray->At(track->GetLabel());
					pdg = fMCparticle->GetPdgCode();
					Float_t ptMC= fMCparticle->Pt();
					fPT_elec_MCtrue->Fill(ptMC);
					
					if(TMath::Abs(pdg)!=11) continue ;      //Is electron:
					
					fpt_incl_elec_MC->Fill(pt);
					fElecNos->Fill(1);
								
					Int_t fMCmotherindex=fMCparticle->GetMother(); //Getting Electron Mother
					if(fMCmotherindex<0) continue ;
					
					fMCmother= (AliAODMCParticle*)fMCarray->At(fMCparticle->GetMother());
					Int_t pdg_mother = fMCmother->GetPdgCode();
								
					//----If mother is Pi0 , eta or gamma --------------------------------------------------------------
					if(TMath::Abs(pdg_mother) == 111 || TMath::Abs(pdg_mother) == 22 || TMath::Abs(pdg_mother) == 221)
					{
						fElecNos->Fill(2);
						fPt_elec_phot->Fill(pt);
						fPt_elec_phot_multSPD->Fill(pt,SPDntr);
						fPt_elec_phot_multV0M->Fill(pt,V0Mmult);
							
						if(TMath::Abs(pdg_mother) == 111) fPi0_Pt_MC->Fill(fMCmother->Pt());
						if(TMath::Abs(pdg_mother) == 221) fEta_Pt_MC->Fill(fMCmother->Pt());
						if(TMath::Abs(pdg_mother) == 22) fGamma_Pt_MC->Fill(fMCmother->Pt());
							
						Double_t ptmotherw = -1.;	
						Int_t elec_source = GetElecSourceType(fMCparticle,fMCarray,ptmotherw);
						//cout<<" elec_source ="<<elec_source<<endl;
										
						if((elec_source==kPi0NoFeedDown) || (elec_source==kGPi0NoFeedDown) || (elec_source==kEtaNoFeedDown) || (elec_source==kGEtaNoFeedDown))
						{						
							fPt_elec_phot1->Fill(pt);
							fPt_elec_phot1_multSPD->Fill(pt,SPDntr);
							fPt_elec_phot1_multV0M->Fill(pt,V0Mmult);
												
							//cout<<" hi 2 "<<endl;
							//getchar();					
							SelectPhotonicElectronR(iTracks, track, fMCmotherindex, pdg,elec_source, V0Mmult, SPDntr);
							//fNonHFE->SetHistMassBack(fInvmassLS_MC);
							//fNonHFE->SetHistMass(fInvmassULS_MC);
						  fNonHFE->FindNonHFE(iTracks,track,fAOD);
						  		
						  Int_t fNULS_MC = fNonHFE->GetNULS();
						  Int_t fNLS_MC = fNonHFE->GetNLS();
						 	
						 	//cout<<" fNULS_MC = "<<fNULS_MC<<" fNULS = "<<fNULS<<" fNLS_MC = "<< fNLS_MC <<" fNLS = "<<fNLS<<endl;
						 	//getchar();
						  if(fNonHFE->IsULS())
						  {
								fPte_ULS_MC->Fill(track->Pt(),fNULS_MC);
								fPte_ULS_MC_multV0M->Fill(track->Pt(),V0Mmult,fNULS_MC);
								fPte_ULS_MC_multSPD->Fill(track->Pt(),SPDntr,fNULS_MC);
							}

						  if(fNonHFE->IsLS())
						  {
								fPte_LS_MC->Fill(track->Pt(),fNLS_MC);
								fPte_LS_MC_multV0M->Fill(track->Pt(),V0Mmult,fNLS_MC);
								fPte_LS_MC_multSPD->Fill(track->Pt(),SPDntr,fNLS_MC);
							} 
							
							if((elec_source==kPi0NoFeedDown)){
								fPt_elec_pi0_MB->Fill(pt);  
								if(fNonHFE->IsULS()) fPt_elec_pi0_MB_ULS->Fill(pt,fNULS_MC);
								if(fNonHFE->IsLS()) fPt_elec_pi0_MB_LS->Fill(pt,fNLS_MC);
							}
							if((elec_source==kEtaNoFeedDown)){
								fPt_elec_eta_MB->Fill(pt);  
								if(fNonHFE->IsULS()) fPt_elec_eta_MB_ULS->Fill(pt,fNULS_MC);
								if(fNonHFE->IsLS()) fPt_elec_eta_MB_LS->Fill(pt,fNLS_MC);
							}
							if((elec_source==kGPi0NoFeedDown)||(elec_source==kGEtaNoFeedDown)){
								fPt_elec_gamma_MB->Fill(pt);  
								if(fNonHFE->IsULS()) fPt_elec_gamma_MB_ULS->Fill(pt,fNULS_MC);
								if(fNonHFE->IsLS()) fPt_elec_gamma_MB_LS->Fill(pt,fNLS_MC);
							}	
						}			
					}	
					else fpt_e_nonphot_MC->Fill(track->Pt());	
							
					Int_t IsElecHf=GetHFE(fMCparticle,fMCarray);
					if((IsElecHf==kBeauty) || (IsElecHf==kCharm))
					{
						fPtHFEMC_reco->Fill(pt);
						fPtHFEMC_reco_true->Fill(ptMC);
						fPtHFEMC_reco_SPD->Fill(pt,SPDntr);
						fPtHFEMC_reco_V0M->Fill(pt,V0Mmult);
					}		
				} //------------fIsMC Loop------------------
			}	//------------TPC n sigma Cut (-1,3)--------
		}//---TMath::Abs(TOFnsigma) < 3 cut------------
	}//---------------END of track loop----------------------------
	////////////////////////////////////////////////////////////////
	
	fHistMult->Fill(num);
	delete fListOfmotherkink;
}
//--------------------END of UserExec----------------------------------


//_________________________________________
/*
void AliAnalysisTaskQAHFE::SelectPhotonicElectronR(Int_t itrack, AliAODTrack *track, Int_t motherindex, Int_t pdg1, Int_t source,Double_t V0Mmult1 , Double_t SPDntr1)
{
 
    // load MC array
    AliMCEvent* mcEvent;
    AliMCEventHandler* eventHandler;
    AliAODMCHeader *mcHeader;
    TClonesArray* mcArray;
    
    mcArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!mcArray){
            AliError("Array of MC particles not found");
            return;
    }
        
    mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!mcHeader) {
         AliError("Could not find MC Header in AOD");
         return;
    }
   
	for(Int_t jTracks = 0; jTracks<fAOD->GetNumberOfTracks(); jTracks++)
  {
		AliVParticle* VtrackAsso = fAOD->GetTrack(jTracks);
		if (!VtrackAsso)
		{
            printf("ERROR: Could not receive track %d\n", jTracks);
            continue;
    }
        
    AliAODTrack *trackAsso = dynamic_cast<AliAODTrack*>(VtrackAsso);
        
   //track cuts applied
    Int_t pdgass = 0;

       
		AliAODTrack *atrackAsso = dynamic_cast<AliAODTrack*>(VtrackAsso);
		if(!atrackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
		if(atrackAsso->GetTPCNcls() < fAssoTPCCluster) continue;
		if(fAssoITSRefit)
		{
			if(!(atrackAsso->GetStatus()&AliESDtrack::kITSrefit)) continue;
		}
		if(!(atrackAsso->GetStatus()&AliESDtrack::kTPCrefit)) continue;
		AliAODMCParticle *MCass = (AliAODMCParticle*)mcArray->At(TMath::Abs(atrackAsso->GetLabel()));
		if(!(TMath::Abs(MCass->GetPdgCode())==11)) continue;
		pdgass=MCass->GetPdgCode();
		Int_t indexass = MCass->GetMother();
		if(TMath::Abs(indexass-motherindex)>0.8) continue;
 
		if(jTracks==itrack) continue;
        
		Double_t  ptAsso=-999., nsigma=-999.0;
		Double_t mass=-999., width = -999;
		Bool_t fFlagULS=kFALSE;
		//Double_t fAssopTMin=0.1;
        
		nsigma = fpidResponse->NumberOfSigmasTPC(trackAsso, AliPID::kElectron); // fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC
		ptAsso = trackAsso->Pt();
		Int_t chargeAsso = trackAsso->Charge();
		Int_t charge = track->Charge();
        
        //---------------pt and track cuts-----------------
        if(ptAsso < fAssopTMin) continue;
		if(TMath::Abs(trackAsso->Eta())>fAssoEtarange) continue;
		if(TMath::Abs(nsigma) > fAssoTPCnsig ) continue;

        //-------------------AliKFParticle-------------------
        Int_t PDGe1 = 11; Int_t PDGe2 = 11;
        if(charge>0) PDGe1 = -11;
        if(chargeAsso>0) PDGe2 = -11;
        
        if((pdgass*pdg1)<0.) fFlagULS = kTRUE;
        
        AliKFParticle ge1(*track, PDGe1);
        AliKFParticle ge2(*trackAsso, PDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        recg.GetMass(mass,width);
        // Mother associated track
        
        Double_t invms[4];
    	  invms[0]=track->Pt();
		    invms[1]=source;
        invms[2]=mass;
        invms[3]=SPDntr1;
        
        if(fFlagULS){
        	fInvMULS->Fill(track->Pt(),mass);
        	fInvMULS_multV0M->Fill(track->Pt(),mass,V0Mmult1);
        	fInvMULS_multSPD->Fill(track->Pt(),mass,SPDntr1);
        	
        	fInvMULSnSp->Fill(invms);
			if(source==kPi0NoFeedDown)fInvMULSpi0->Fill(track->Pt(),mass);
			if(source==kEtaNoFeedDown)fInvMULSeta->Fill(track->Pt(),mass);
			if((source==kGPi0NoFeedDown)||(source==kGEtaNoFeedDown))fInvMULSgamma->Fill(track->Pt(),mass);
		}
	}
}*/

void AliAnalysisTaskQAHFE::SelectPhotonicElectronR(Int_t itrack, AliAODTrack *track, Int_t motherindex, Int_t pdg1, Int_t source,Double_t V0Mmult1 , Double_t SPDntr1)
{
 
    // load MC array
    AliMCEvent* mcEvent;
    AliMCEventHandler* eventHandler;
    AliAODMCHeader *mcHeader;
    TClonesArray* mcArray;
    //cout<<" hi 3"<<endl;
     // getchar(); 
    mcArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!mcArray){
            AliError("Array of MC particles not found");
            return;
    }
        
    mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!mcHeader) {
         AliError("Could not find MC Header in AOD");
         return;
    }
   
	for(Int_t jTracks = 0; jTracks<fAOD->GetNumberOfTracks(); jTracks++)
  {
		AliVParticle* VtrackAsso = fAOD->GetTrack(jTracks);
		if (!VtrackAsso)
		{
            printf("ERROR: Could not receive track %d\n", jTracks);
            continue;
    }
        
    AliAODTrack *trackAsso = dynamic_cast<AliAODTrack*>(VtrackAsso);
        
   //track cuts applied
    Int_t pdgass = 0;
	//cout<<" hi "<<endl;
     // getchar(); 
		AliAODTrack *atrackAsso = dynamic_cast<AliAODTrack*>(VtrackAsso);
		if(!atrackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
		if(atrackAsso->GetTPCNcls() < fAssoTPCCluster) continue;
		if(fAssoITSRefit)
		{
			if(!(atrackAsso->GetStatus()&AliESDtrack::kITSrefit)) continue;
		}
		if(!(atrackAsso->GetStatus()&AliESDtrack::kTPCrefit)) continue;
		AliAODMCParticle *MCass = (AliAODMCParticle*)mcArray->At(TMath::Abs(atrackAsso->GetLabel()));
		if(!(TMath::Abs(MCass->GetPdgCode())==11)) continue;
		pdgass=MCass->GetPdgCode();
		
        
		Double_t  ptAsso=-999., nsigma=-999.0;
		Double_t mass=-999., width = -999;
		Bool_t fFlagULS=kFALSE;
		//Double_t fAssopTMin=0.1;
        
		nsigma = fpidResponse->NumberOfSigmasTPC(trackAsso, AliPID::kElectron); // fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC
		ptAsso = trackAsso->Pt();
		Int_t chargeAsso = trackAsso->Charge();
		Int_t charge = track->Charge();
        
        //---------------pt and track cuts-----------------
        if(ptAsso < fAssopTMin) continue;
		if(TMath::Abs(trackAsso->Eta())>fAssoEtarange) continue;
		if(TMath::Abs(nsigma) > fAssoTPCnsig ) continue;

        //-------------------AliKFParticle-------------------
        Int_t PDGe1 = 11; Int_t PDGe2 = 11;
        if(charge>0) PDGe1 = -11;
        if(chargeAsso>0) PDGe2 = -11;
        
        if((pdgass*pdg1)<0.) fFlagULS = kTRUE;
        
        AliKFParticle ge1(*track, PDGe1);
        AliKFParticle ge2(*trackAsso, PDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        recg.GetMass(mass,width);
        // Mother associated track
        
        //cout<<" mass = "<<mass<<endl;
        //getchar();
        //COMBINATIONS--------------
        
        if(mass<fInvmassCut)
            {
                if(charge*chargeAsso<0)
                {
                    fPte_ULS_KF->Fill(track->Pt());
										//fPte_ULS_multV0M_KF->Fill(track->Pt(),V0Mmult1);
										//fPte_ULS_multSPD_KF->Fill(track->Pt(),SPDntr1);
                    
                }
                if(charge*chargeAsso>0)
                {
                    
                    fPte_LS_KF->Fill(track->Pt());
										//fPte_LS_multV0M_KF->Fill(track->Pt(),V0Mmult1);
										//fPte_LS_multSPD_KF->Fill(track->Pt(),SPDntr1);
                }
            }
        
     //TRUE PAIRS-----------------   
        Int_t indexass = MCass->GetMother();
				if(TMath::Abs(indexass-motherindex)<0.8)
				{
						if(jTracks==itrack) continue;
				    Double_t invms[4];
					  invms[0]=track->Pt();
						invms[1]=source;
				    invms[2]=mass;
				    invms[3]=SPDntr1;
				    
				    if(fFlagULS){
				    	fInvMULS->Fill(track->Pt(),mass);
				    	//fInvMULS_multV0M->Fill(track->Pt(),mass,V0Mmult1);
				    	//fInvMULS_multSPD->Fill(track->Pt(),mass,SPDntr1);
				    	
				    	fInvMULSnSp->Fill(invms);
					if(source==kPi0NoFeedDown)fInvMULSpi0->Fill(track->Pt(),mass);
					if(source==kEtaNoFeedDown)fInvMULSeta->Fill(track->Pt(),mass);
					if((source==kGPi0NoFeedDown)||(source==kGEtaNoFeedDown))fInvMULSgamma->Fill(track->Pt(),mass);
			}
			
		}
	}
}
//=================================================================================================================================
Int_t AliAnalysisTaskQAHFE::ClassifyTrack(AliAODTrack* track,const AliVVertex* pVtx)
{  
  
	Double_t pt = track->Pt();
	Double_t eta = track->Eta();
	Double_t phi = track->Phi();
	Float_t dx,dy,dxy, dz;

	//if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return 0; //fitler bit 
  if(pt< 0.3 ) return 0;  
	//if (TMath::Abs(eta)>=fEtarange) return 0; 
	
	Double_t nclus = track->GetTPCNcls();  // TPC cluster information
	Double_t nclusF = track->GetTPCNclsF();
 	Double_t nclusN = track->GetTPCsignalN();  // TPC cluster information findable
 	Double_t RatioTPCclusters=nclusN/nclusF;
 	Double_t nclusS = track->GetTPCnclsS();
 	
	if(track->GetTPCNcls() < fTPCNclus) return 0; //TPC N clusters
	if(track->GetITSNcls() < fITSNclus) return 0; // ITS N clusters
	if(nclusN< fTPCNclusPID) return 0 ;
	if(RatioTPCclusters<0.6) return 0;
	
	if((!(track->GetStatus()&AliESDtrack::kITSrefit)|| (!(track->GetStatus()&AliESDtrack::kTPCrefit)))) return 0; // ITS and TPC refit
	
	if(fSPDBoth){ if(!(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1))) return 0;} //Hit on first and second SPD layer
	else if(fSPDAny){ if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) return 0;} //Hit on any layer
	else if(fSPDFirst){ if(!(track->HasPointOnITSLayer(0))) return 0;} //Hit on first and second SPD layer

	//==DCA Cut ======
	Double_t d0z0[2]={-999,-999}, cov[3];
	if(track->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov))
	fDCAxy= d0z0[0];
	fDCAz= d0z0[1];
	if(TMath::Abs(d0z0[0]) > fDCAxyCut || TMath::Abs(d0z0[1]) > fDCAzCut) return 0; 
  
	Double_t chi2ndf = track->Chi2perNDF();
	if(chi2ndf>4.0) return 0;
			
	return 1;
}
//====================================================================================================================================
Int_t AliAnalysisTaskQAHFE::GetElecSourceType(AliAODMCParticle *electron, TClonesArray *mcArray,Double_t &ptm)
{
    //
    // Return what type of gammax it is
    //
    
    // Mother
    Int_t motherlabel = electron->GetMother();
    if(motherlabel<0) return kNoMotherE;
    else {
        
        AliAODMCParticle *mother = (AliAODMCParticle*)mcArray->At(motherlabel);
        Int_t motherpdg = TMath::Abs(mother->GetPdgCode());
        ptm=mother->Pt();
        if(motherpdg == 111) {
            Int_t typepi0eta = GetPi0EtaType(mother,mcArray);
            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) 	return kPi0NoFeedDown;
        }
        if(motherpdg == 221) {
            Int_t typepi0eta = GetPi0EtaType(mother,mcArray);
            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kEtaNoFeedDown;
        }
        if(motherpdg == 22) {
            Int_t gmotherlabel = mother->GetMother();
            if(gmotherlabel<0) return kDirectGamma;
            else {
                AliAODMCParticle *gmother = (AliAODMCParticle*)mcArray->At(gmotherlabel);
                ptm=gmother->Pt();
                Int_t gmotherpdg = TMath::Abs(gmother->GetPdgCode());
                if(gmotherpdg == 111) {
                    Int_t typepi0eta = GetPi0EtaType(gmother,mcArray);
                    if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGPi0NoFeedDown;
                }
                if(gmotherpdg == 221) {
                    Int_t typepi0eta = GetPi0EtaType(gmother,mcArray);
                    if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGEtaNoFeedDown;
                }
                if(gmotherpdg == 22) {
                    Int_t ggmotherlabel = gmother->GetMother();
                    if(ggmotherlabel<0) return kDirectGamma;
                    else {
                        AliAODMCParticle *ggmother = (AliAODMCParticle*)mcArray->At(ggmotherlabel);
                        ptm=ggmother->Pt();
                        Int_t ggmotherpdg = TMath::Abs(ggmother->GetPdgCode());
                        if(ggmotherpdg == 111) {
                            Int_t typepi0eta = GetPi0EtaType(ggmother,mcArray);
                            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGPi0NoFeedDown;
                        }
                        if(ggmotherpdg == 221) {
                            Int_t typepi0eta = GetPi0EtaType(ggmother,mcArray);
                            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGEtaNoFeedDown;
                        }
                    }
                }
            }
        }
    }
    
    return kOthersE;

}
Int_t AliAnalysisTaskQAHFE::GetPi0EtaType(AliAODMCParticle *pi0eta, TClonesArray *mcArray)
{
   
    // Return what type of pi0, eta it is
    
    // IsPrimary
    Bool_t primMC = pi0eta->IsPrimary();
    if(!primMC) return kNoIsPrimary;
    
    // Mother
    Int_t motherlabel = pi0eta->GetMother();
    if(motherlabel<0) return kNoMother;
    else {
        
        AliAODMCParticle *mother = (AliAODMCParticle*)mcArray->At(motherlabel);
        Int_t motherpdg = TMath::Abs(mother->GetPdgCode());
        //    if(motherpdg == 111 || motherpdg == 221 || motherpdg == 223 || motherpdg == 333 || motherpdg == 331 || motherpdg == 113) return kLightMesons;
        if(motherpdg == 111 || motherpdg == 221 || motherpdg == 223 || motherpdg == 333 || motherpdg == 331 || motherpdg == 113 || motherpdg == 213 || motherpdg == 313 || motherpdg == 323) return kLightMesons;
        
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty;
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm;
        return kNoFeedDown;
        
    }
}

Int_t AliAnalysisTaskQAHFE::GetHFE(AliAODMCParticle *electron, TClonesArray *mcArray)
{
	Int_t motherindex=electron->GetMother(); //Getting Electron Mother
	if(motherindex<0) kNoMother;
	AliAODMCParticle *mother = (AliAODMCParticle*)mcArray->At(motherindex);					
	Int_t motherpdg = mother->GetPdgCode();		  
	if ( (int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty;
	if ( (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm;
	return kOthersE;
  
}

//====================================================================================================================================
void AliAnalysisTaskQAHFE::SelectPhotonicElectronMC(Int_t itrack,AliAODTrack *track, AliPIDResponse *pidResponse,Int_t V0Mult1,Int_t nAcceta1)
{
	//Identify non-heavy flavour electrons using Invariant mass method using DCA method
  	fNonHFE = new AliSelectNonHFE();
  	fNonHFE->SetAODanalysis(kTRUE);
  	fNonHFE->SetInvariantMassCut(fInvmassCut);
  	fNonHFE->SetAlgorithm("DCA"); //KF
  	fNonHFE->SetPIDresponse(pidResponse);
  	fNonHFE->SetTrackCuts(-3.5,3.5); //TPCnsigma cuts

	//fNonHFE->SetChi2OverNDFCut(4.0);
	//fNonHFE->SetDCACut(fDCAcut);
	//fNonHFE-> SetDCAPartnerCuts(1,2);   //SetDCAPartnerCuts(Double_t xy, Double_t z)
	fNonHFE->SetAdditionalCuts(0.1,60);  //SetAdditionalCuts(fPtMinAsso,fTpcNclsAsso);

	//fNonHFE->SetHistMassBack(fInvmassLS1);
  	//fNonHFE->SetHistMass(fInvmassULS1);
  	fNonHFE->FindNonHFE(itrack,track,fAOD);
  		
  	Int_t fNULS = fNonHFE->GetNULS();
  	Int_t fNLS = fNonHFE->GetNLS();
 	
 	if(fIsMC)
 	{
  	if(fNonHFE->IsULS())
  	{
		for(Int_t k =0; k< fNonHFE->GetNULS(); k++){
			fPte_ULS_MC->Fill(track->Pt());
			fPte_ULS_MC_multV0M->Fill(track->Pt(),V0Mult1);
			fPte_ULS_MC_multSPD->Fill(track->Pt(),nAcceta1);
		}
	}

  	if(fNonHFE->IsLS())
  	{
		for(Int_t k=0; k < fNonHFE->GetNLS(); k++){
			fPte_LS_MC->Fill(track->Pt());
			fPte_LS_MC_multV0M->Fill(track->Pt(),V0Mult1);
			fPte_LS_MC_multSPD->Fill(track->Pt(),nAcceta1);
		}
	} 
	} 
}


Int_t AliAnalysisTaskQAHFE::GetNcharged(){
    //counts all tracks in eta<1 with charge!=0
	
    Int_t Nch = 0;

		
    if(!fIsMC) return Nch; // if no MC info return 0

    // loop over all tracks 
    for (Int_t igen = 0; igen < fMCarray->GetEntriesFast(); igen++){
        AliAODMCParticle *mctrack=(AliAODMCParticle*)fMCarray->UncheckedAt(igen);
        Int_t charge = mctrack->Charge();
        Double_t eta = mctrack->Eta();
        Bool_t isPhysPrim = mctrack->IsPhysicalPrimary();
        if(charge!=0){
            if(eta > -1.0 && eta < 1.0){
                if(isPhysPrim){
                    Nch++;
                }
            }
        }
    }
    return Nch;
}
//=======================================================================

//===================================================================================
Double_t AliAnalysisTaskQAHFE::Beta(AliAODTrack *track)
{
  	Double_t stoptime=track->GetTOFsignal();
  
  	Double_t c=TMath::C()*1.E-9;// m/ns
  	Float_t startTime= fpidResponse->GetTOFResponse().GetStartTime(((AliAODTrack*)track)->P());//in ps
  	Double_t length= fpidResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kElectron)*1E-3*c;
  	stoptime -= startTime;
  	Double_t scaleStopTime= stoptime*1E-3;
  	scaleStopTime=scaleStopTime*c;
  	return length/scaleStopTime;
}


//===================================================================================
void AliAnalysisTaskQAHFE::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  
	fOutputList = dynamic_cast<TList*> (GetOutputData(1));
	if (!fOutputList)return;

}

