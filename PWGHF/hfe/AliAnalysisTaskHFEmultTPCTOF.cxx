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
#include "AliAnalysisTaskHFEmultTPCTOF.h"
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
//#include  "AliHelperPID.h"

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
#include "TRandom.h"
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskHFEmultTPCTOF)

//________________________________________________________________________//
AliAnalysisTaskHFEmultTPCTOF::AliAnalysisTaskHFEmultTPCTOF(): 
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
fnSigmaVsPt_TPC_cut(0),
fnSigmaVsP_TOF(0), 
fnBetaVsP_TOF(0),
fnSigmaVsEta_TPC(0), 
fnSigmaVsEta_TPC_cut(0), 

fdEdxVsPt_TPC(0),
fnSigmaVsPt_TPC(0),
fnSigmaVsPt_TOF(0), 
fnBetaVsPt_TOF(0), 

fLandau(0),fErr(0),fHadCot_Landau(0),fHadCot_Err(0),fPt_incl_e_Landau(0),fPt_incl_e_Err(0),func_MCnchCorr(0),

//Multiplicity dep

fPteV0M_Lan(0),fPteV0M_Err(0), fPteSPD_Lan(0),fPteSPD_Err(0),
fPteV0M(0),fPteSPD(0),
fHadV0MLan(0),fHadSPDLan(0),fHadV0MErr(0),fHadSPDErr(0),fnSigmaVsP_TPC_cut_mult(0),
//fnSigmaVsP_TPC_cut_multV0M(0),fnSigmaVsP_TPC_cut_multSPD(0),
//fnSigmaVsPt_TPC_cut_multV0M(0),fnSigmaVsPt_TPC_cut_multSPD(0),

//------------Ntracklets Correction and Multiplicity------------------

fMultiPercentile(0),fMultiPercentileSPD(0),
fV0Mult(0),fV0MultCorr(0),fV0AMult(0),fV0CMult(0),fSPD_tracklet_NoEtaCut(0),fSPD_tracklet(0), 
fCent(0),fCentSPD(0), 
fMultV0M_vs_alimult(0),fMultSPD_vs_alimult(0),
fileEstimator(0),
fSPDCorrMultDist_min(0),fSPDCorrMultDist_max(0),fSPDWeightedCorrMultDist_max(0),
fHistNtrVsZvtx(0),fHistNtrCorrVsZvtx_min(0),fHistNtrCorrVsZvtx_max(0),Profile_Mean(0),Profile_MeanCorr(0),
fRefMult(13.079),
fV0MultVsZvtx(0),fV0MultCorrVsZvtx(0),
fSPDCorrMultDist_min_vs_AliSPDMult(0),
fNchVsZvtx(0),SPDNtrCorrVsNch(0),SPDNtrCorrVsNchweights(0),V0MCorrVsNch(0),fSPDNtrCorrVsV0MCorrVsNch(0),
 
    fPtHFEMC_aftertrackcut(0),
    fPtHFEMC_aftertofcut(0),
    fPtHFEMC_aftertpcnsigcut(0),
    fPtHFEMC_aftertoftpcnsigcut(0),
    fPtHFEMC_aftertrackcut_SPD(0),
    fPtHFEMC_aftertofcut_SPD(0),
    fPtHFEMC_aftertpcnsigcut_SPD(0),
    fPtHFEMC_aftertoftpcnsigcut_SPD(0),
  
//------------------------Non-Hfe
fNonHFE(new AliSelectNonHFE()),
fInvmassLS1(0),fInvmassULS1(0),
fPte_ULS(0),fPte_LS(0),
fPte_ULS_multV0M(0),fPte_LS_multV0M(0),
fPte_ULS_multSPD(0),fPte_LS_multSPD(0),

//------------------------------MC 
fIsMC(kFALSE),
fMCarray(0),fMCparticle(0),fMCmother(0),

fPtHFEMCtot(0),fPtHFEMCtot_SPD(0),fPtHFEMCtot_V0M(0),
fPtHFEMC(0),
fPtHFEMC_SPD(0),
fPtHFEMC_V0M(0),

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
fInvMULSnSpwg(0),
fInvMULSpi0(0),
fInvMULSeta(0),
fInvMULSgamma(0),

fPtHFEMC_reco(0),
fPtHFEMC_reco_SPD(0),
fPtHFEMC_reco_V0M(0),

fPi0EtaSpectra(0),
fPi0EtaSpectraSp(0),
fPt_elec_phot2(0),
fPtElec_ULS_MC2(0),
fPtElec_LS_MC2(0),
fPt_elec_phot2_multSPD(0),
fPtElec_ULS_MC2_multSPD(0),
fPtElec_LS_MC2_multSPD(0),
fPt_elec_from_pi02(0),
fPtElec_ULS_MC_from_pi02(0),
fPtElec_LS_MC_from_pi02(0),
fPt_elec_from_eta2(0),
fPtElec_ULS_MC_from_eta2(0),
fPtElec_LS_MC_from_eta2(0),
fPt_elec_from_gamma2(0),
fPtElec_ULS_MC_from_gamma2(0),
fPtElec_LS_MC_from_gamma2(0),
f0a(0),f0b(0),f1(0),f1a(0),f1b(0)

{
  fPID = new AliHFEpid("hfePid");
  for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=0;
}


AliAnalysisTaskHFEmultTPCTOF::AliAnalysisTaskHFEmultTPCTOF(const char *name) 
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
fnSigmaVsPt_TPC_cut(0),
fnSigmaVsP_TOF(0), 
fnBetaVsP_TOF(0),
fnSigmaVsEta_TPC(0), 
fnSigmaVsEta_TPC_cut(0),  

fdEdxVsPt_TPC(0),
fnSigmaVsPt_TPC(0),
fnSigmaVsPt_TOF(0), 
fnBetaVsPt_TOF(0),

fLandau(0),fErr(0),fHadCot_Landau(0),fHadCot_Err(0),fPt_incl_e_Landau(0),fPt_incl_e_Err(0),func_MCnchCorr(0),

//Multiplicity dep

fPteV0M_Lan(0),fPteV0M_Err(0), fPteSPD_Lan(0),fPteSPD_Err(0),
fPteV0M(0),fPteSPD(0),
fHadV0MLan(0),fHadSPDLan(0),fHadV0MErr(0),fHadSPDErr(0),fnSigmaVsP_TPC_cut_mult(0),
//fnSigmaVsP_TPC_cut_multV0M(0),fnSigmaVsP_TPC_cut_multSPD(0),
//fnSigmaVsPt_TPC_cut_multV0M(0),fnSigmaVsPt_TPC_cut_multSPD(0),

//------------Ntracklets Correction and Multiplicity------------------

fMultiPercentile(0),fMultiPercentileSPD(0),
fV0Mult(0),fV0MultCorr(0),fV0AMult(0),fV0CMult(0),fSPD_tracklet_NoEtaCut(0),fSPD_tracklet(0), 
fCent(0),fCentSPD(0), 
fMultV0M_vs_alimult(0),fMultSPD_vs_alimult(0),
fileEstimator(0),
fSPDCorrMultDist_min(0),fSPDCorrMultDist_max(0),fSPDWeightedCorrMultDist_max(0),
fHistNtrVsZvtx(0),fHistNtrCorrVsZvtx_min(0),fHistNtrCorrVsZvtx_max(0),Profile_Mean(0),Profile_MeanCorr(0),
fRefMult(13.079),
fV0MultVsZvtx(0),fV0MultCorrVsZvtx(0),
fSPDCorrMultDist_min_vs_AliSPDMult(0),
fNchVsZvtx(0),SPDNtrCorrVsNch(0),SPDNtrCorrVsNchweights(0),V0MCorrVsNch(0),fSPDNtrCorrVsV0MCorrVsNch(0),

    fPtHFEMC_aftertrackcut(0),
    fPtHFEMC_aftertofcut(0),
    fPtHFEMC_aftertpcnsigcut(0),
    fPtHFEMC_aftertoftpcnsigcut(0),
    fPtHFEMC_aftertrackcut_SPD(0),
    fPtHFEMC_aftertofcut_SPD(0),
    fPtHFEMC_aftertpcnsigcut_SPD(0),
    fPtHFEMC_aftertoftpcnsigcut_SPD(0),
  
//------------------------Non-Hfe
fNonHFE(new AliSelectNonHFE()),
fInvmassLS1(0),fInvmassULS1(0),
fPte_ULS(0),fPte_LS(0),
fPte_ULS_multV0M(0),fPte_LS_multV0M(0),
fPte_ULS_multSPD(0),fPte_LS_multSPD(0),

//------------------------------MC 
fIsMC(kFALSE),
fMCarray(0),fMCparticle(0),fMCmother(0),

fPtHFEMCtot(0),fPtHFEMCtot_SPD(0),fPtHFEMCtot_V0M(0),
fPtHFEMC(0),
fPtHFEMC_SPD(0),
fPtHFEMC_V0M(0),
  
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
fInvMULSnSpwg(0),
fInvMULSpi0(0),
fInvMULSeta(0),
fInvMULSgamma(0),

fPtHFEMC_reco(0),
fPtHFEMC_reco_SPD(0),
fPtHFEMC_reco_V0M(0),

fPi0EtaSpectra(0),
fPi0EtaSpectraSp(0),

fPt_elec_phot2(0),
fPtElec_ULS_MC2(0),
fPtElec_LS_MC2(0),
fPt_elec_phot2_multSPD(0),
fPtElec_ULS_MC2_multSPD(0),
fPtElec_LS_MC2_multSPD(0),
fPt_elec_from_pi02(0),
fPtElec_ULS_MC_from_pi02(0),
fPtElec_LS_MC_from_pi02(0),
fPt_elec_from_eta2(0),
fPtElec_ULS_MC_from_eta2(0),
fPtElec_LS_MC_from_eta2(0),
fPt_elec_from_gamma2(0),
fPtElec_ULS_MC_from_gamma2(0),
fPtElec_LS_MC_from_gamma2(0),
f0a(0),f0b(0),f1(0),f1a(0),f1b(0)



{
  // Constructor
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container

  fPID = new AliHFEpid("hfePid");
  DefineInput(0, TChain::Class());  
  DefineOutput(1, TList::Class());
  DefineOutput(2, TH1F::Class());
  for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=0;
}

//_________________________Destructer_____________________________________
AliAnalysisTaskHFEmultTPCTOF::~AliAnalysisTaskHFEmultTPCTOF()
{
  delete fPID;
  if (fOutputList) { delete fOutputList; fOutputList = 0;}
  if (fNentries){ delete fNentries; fNentries = 0;}
  if (fNentries2){ delete fNentries2; fNentries2 = 0;}
  
     for(Int_t i=0; i<4; i++) {
      if (fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i];
  }
}

//_______________________________________________________________________
void AliAnalysisTaskHFEmultTPCTOF::Init()
{
  // Initialization	
  if(fDebug > 1) printf("AliAnalysisTaskHFEmultTPCTOF::Init() \n");
  return;
}




//________________________________________________________________________
void AliAnalysisTaskHFEmultTPCTOF::UserCreateOutputObjects()
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

  f0a = new TF1("f0a","[0] / (TMath::Power(TMath::Exp( - [1] * x - [2] * x * x ) + x / [3], [4]))",0.5,2.5);
  f0b = new TF1("f0b","[0] / (TMath::Power(TMath::Exp( - [1] * x - [2] * x * x ) + x / [3], [4]))",2.5,10);
  
  f0a->FixParameter(0,5.80666e+01);
  f0a->FixParameter(1,-1.64680e+01);
  f0a->FixParameter(2,1.25223e+01);
  f0a->FixParameter(3,1.25515e-03);
  f0a->FixParameter(4,4.71063e-01 );
  
  f0b->FixParameter(0,1.54389e+00);
  f0b->FixParameter(1,-1.43139e+03);
  f0b->FixParameter(2,4.58994e+02);
  f0b->FixParameter(3,2.04094e+03);
  f0b->FixParameter(4,-9.34144e-05);
  
  
  f1 = new TF1("f1","[0] / (TMath::Power(TMath::Exp( - [1] * x - [2] * x * x ) + x / [3], [4]))",0.5,2.3);
  f1a = new TF1("f1a","[0] / (TMath::Power(TMath::Exp( - [1] * x - [2] * x * x ) + x / [3], [4]))",2.3,4.1);
  f1b = new TF1("f1b","[0] / (TMath::Power(TMath::Exp( - [1] * x - [2] * x * x ) + x / [3], [4]))",4.0,10.);
 
  f1->FixParameter(0,2.14780e+00);
  f1->FixParameter(1,-2.10239e+00);
  f1->FixParameter(2,4.14568e+00);
  f1->FixParameter(3,1.87072e+00);
  f1->FixParameter(4,3.06386e-01);
  
  f1a->FixParameter(0,9.59262e+00 );
  f1a->FixParameter(1,-1.43061e+01);
  f1a->FixParameter(2,2.60984e+01 );
  f1a->FixParameter(3,1.74946e-04);
  f1a->FixParameter(4,1.63876e-01);
  
  f1b->FixParameter(0,1.58238e+00);
  f1b->FixParameter(1,-6.25753e+01);
  f1b->FixParameter(2,7.27244e+00);
  f1b->FixParameter(3,1.83220e-14);
  f1b->FixParameter(4,-1.11985e-03);
	
  func_MCnchCorr= new TF1("func_MCnchCorr","expo(0)+expo(2)", 40.,120.);
  func_MCnchCorr->SetParameters(-1.88813e+00, 3.43870e-02,1.18702e+00,-5.74526e-02);
  
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
	
  Double_t pi=TMath::Pi();

  fcount = new TH1F("fcount", "fcount", 10, 0.0, 10.0);
  fcount->Sumw2();  
  
  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 150, 0.1, 15.1);
  //fOutputList->Add(fHistPt);
  fHistMult = new TH1F("fHistMult", "Multiplicity distribution", 500, 0, 500);
  //fOutputList->Add(fHistMult);
  fTPCSignal = new TH1F("fTPCSignal", "fTPCSignal distribution", 100, 0, 100);
  //fOutputList->Add(fTPCSignal);
	
	//--------------------------Vertex Dist------------------------------------
	fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
  fOutputList->Add(fVtxZ);

	fVtxZ_corr=new TH1F("fVtxZ_corr","Z vertex position;Vtx_{z};counts",1000,-50,50);
  fOutputList->Add(fVtxZ_corr);

  fVtxY = new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",1000,-50,50);
  //fOutputList->Add(fVtxY);

  fVtxX = new TH1F("fVtxX","X vertex position;Vtx_{x};counts",1000,-50,50);
  //fOutputList->Add(fVtxX);
  
 
  //--------------------------DCA Dist------------------------------------
	fDCAZvsPt = new TH2F("fDCAZvsPt","DCAZ ;p_{T};DCA_{z}",300,0,100,100,-5,5);
  //fOutputList->Add(fDCAZvsPt);

  fDCAXYvsPt = new TH2F("fDCAXYvsPt","DCAXY ;p_{T};DCA_{xy}",300,0,100,100,-5,5);
  //fOutputList->Add(fDCAXYvsPt);

	
	//--------------------------TPC dedx nsig vs p ----------------------------------------------
  fdEdxVsP_TPC = new TH2F("fdEdxVsP_TPC", "fdEdxVsP_TPC distribution",300,0,15,750,10,160);
  fdEdxVsP_TPC->Sumw2();
  //fOutputList->Add(fdEdxVsP_TPC);
  
  fnSigmaVsP_TPC= new TH2F("fnSigmaVsP_TPC", "fnSigmaVsP_TPC distribution",300,0.,15.,750,-15.,15.);
  fnSigmaVsP_TPC->Sumw2(); 
  fOutputList->Add(fnSigmaVsP_TPC);
  
  fdEdxVsP_TPC_cut = new TH2F("fdEdxVsP_TPC_cut", "fdEdxVsP_TPC_cut distribution",300,0,15,750,10,160);
  fdEdxVsP_TPC_cut->Sumw2();
  //fOutputList->Add(fdEdxVsP_TPC_cut);
  
 	fnSigmaVsP_TPC_cut = new TH2F("fnSigmaVsP_TPC_cut", "fnSigmaVsP_TPC_cut distribution",300,0,15,250,-15,10); 
  fnSigmaVsP_TPC_cut->Sumw2(); 
  fOutputList->Add(fnSigmaVsP_TPC_cut);
  
  	fnSigmaVsPt_TPC_cut = new TH2F("fnSigmaVsPt_TPC_cut", "fnSigmaVsPt_TPC_cut distribution",300,0,15,250,-15,10); 
  fnSigmaVsPt_TPC_cut->Sumw2(); 
  fOutputList->Add(fnSigmaVsPt_TPC_cut);

  
  fnSigmaVsP_TOF= new TH2F("fnSigmaVsP_TOF", "fnSigmaVsP_TOF distribution",300,0.,15.,1000,-10.,20.);
  fnSigmaVsP_TOF->Sumw2(); 
  fOutputList->Add(fnSigmaVsP_TOF);
  
  fnBetaVsP_TOF= new TH2F("fnBetaVsP_TOF", "fnBetaVsP_TOF distribution",300,0.,15.,1000,-0.5,1.5);
	fnBetaVsP_TOF->Sumw2(); 
	//fOutputList->Add(fnBetaVsP_TOF);


 //--------------------------TPC dedx nsig vs pt ----------------------------------------------
  fdEdxVsPt_TPC = new TH2F("fdEdxVsPt_TPC", "fdEdxVsPt_TPC distribution",300,0,15,750,10,160);
  fdEdxVsPt_TPC->Sumw2();
  //fOutputList->Add(fdEdxVsPt_TPC);
  
  fnSigmaVsPt_TPC= new TH2F("fnSigmaVsPt_TPC", "fnSigmaVsPt_TPC distribution",300,0.,15.,750,-15.,15.);
  fnSigmaVsPt_TPC->Sumw2(); 
  //fOutputList->Add(fnSigmaVsPt_TPC);
  
  fnSigmaVsPt_TOF= new TH2F("fnSigmaVsPt_TOF", "fnSigmaVsPt_TOF distribution",300,0.,15.,1000,-10.,20.);
  fnSigmaVsPt_TOF->Sumw2(); 
  //fOutputList->Add(fnSigmaVsPt_TOF);
  
  fnBetaVsPt_TOF= new TH2F("fnBetaVsPt_TOF", "fnBetaVsPt_TOF distribution",300,0.,15.,1000,-0.5,1.5);
	fnBetaVsPt_TOF->Sumw2(); 
	//fOutputList->Add(fnBetaVsPt_TOF);


 
 //--------------------------TPC  nsig vs eta ----------------------------------------------
  fnSigmaVsEta_TPC = new TH2F("fnSigmaVsEta_TPC", "fnSigmaVsEta_TPC distribution",1600,-0.8,0.8,750,-15.,15.);
  fnSigmaVsEta_TPC->Sumw2();
  //fOutputList->Add(fnSigmaVsEta_TPC);
  
  fnSigmaVsEta_TPC_cut = new TH2F("fnSigmaVsEta_TPC_cut", "fnSigmaVsEta_TPC_cut distribution",1600,-0.8,0.8,750,-15.,15.);
  fnSigmaVsEta_TPC_cut->Sumw2(); 
  //fOutputList->Add(fnSigmaVsEta_TPC_cut);
  
 
  //-----------------Incl e Spectrum and Hadron Contamination---------------------------------------------
  //fHadCot_Landau = new TH1F("fHadCot_Landau", "fHadCot_Landau Distribution",300,0,15);
  //fHadCot_Landau->Sumw2();
  //fOutputList->Add(fHadCot_Landau); 
  
  fHadCot_Err = new TH1F("fHadCot_Err", "fHadCot_Err Distribution",300,0,15);
  fHadCot_Err->Sumw2();
  fOutputList->Add(fHadCot_Err);
  
  fPt_incl_e_Landau = new TH1F("fPt_incl_e_Landau", "fPt_incl_e_Landau Distribution",300,0,15);
  fPt_incl_e_Landau->Sumw2();
  //fOutputList->Add(fPt_incl_e_Landau); 
  
  fPt_incl_e_Err = new TH1F("fPt_incl_e_Err", "fPt_incl_e_Err Distribution",300,0,15);
  fPt_incl_e_Err->Sumw2();
  fOutputList->Add(fPt_incl_e_Err);

  //=================================================================================

  	fPtHFEMC_aftertrackcut  = new TH1F("fPtHFEMC_aftertrackcut", "fPtHFEMC_aftertrackcut Distribution",300,0,15);
  	fPtHFEMC_aftertrackcut->Sumw2();
  	if(fIsMC) fOutputList->Add(fPtHFEMC_aftertrackcut);

	fPtHFEMC_aftertofcut  = new TH1F("ffPtHFEMC_aftertofcut", "fPtHFEMC_aftertofcut",300,0,15);
	fPtHFEMC_aftertofcut->Sumw2();
  	if(fIsMC) fOutputList->Add(fPtHFEMC_aftertofcut);
	
	fPtHFEMC_aftertpcnsigcut  = new TH1F("fPtHFEMC_aftertpcnsigcut", "fPtHFEMC_aftertpcnsigcut",300,0,15);
	fPtHFEMC_aftertpcnsigcut->Sumw2();
  	if(fIsMC) fOutputList->Add(fPtHFEMC_aftertpcnsigcut);

	fPtHFEMC_aftertoftpcnsigcut  = new TH1F("fPtHFEMC_aftertoftpcnsigcut", "fPtHFEMC_aftertoftpcnsigcut Distribution",300,0,15);
	fPtHFEMC_aftertoftpcnsigcut->Sumw2();
  	if(fIsMC) fOutputList->Add(fPtHFEMC_aftertoftpcnsigcut);

  	fPtHFEMC_aftertrackcut_SPD = new TH2F("fPtHFEMC_aftertrackcut_SPD", "fPtHFEMC_aftertrackcut_SPD Distribution",300,0,15,-0.5,299.5);
  	fPtHFEMC_aftertrackcut_SPD->Sumw2();
  	if(fIsMC) fOutputList->Add(fPtHFEMC_aftertrackcut_SPD);

	fPtHFEMC_aftertofcut_SPD  = new TH2F("ffPtHFEMC_aftertofcut_SPD", "fPtHFEMC_aftertofcut_SPD",300,0,15,-0.5,299.5);
	fPtHFEMC_aftertofcut_SPD->Sumw2();
  	if(fIsMC) fOutputList->Add(fPtHFEMC_aftertofcut_SPD);
	
	fPtHFEMC_aftertpcnsigcut_SPD  = new TH2F("fPtHFEMC_aftertpcnsigcut_SPD", "fPtHFEMC_aftertpcnsigcut_SPD",300,0,15,-0.5,299.5);
	fPtHFEMC_aftertpcnsigcut_SPD->Sumw2();
  	if(fIsMC) fOutputList->Add(fPtHFEMC_aftertpcnsigcut_SPD);

	fPtHFEMC_aftertoftpcnsigcut_SPD  = new TH2F("fPtHFEMC_aftertoftpcnsigcut_SPD", "fPtHFEMC_aftertoftpcnsigcut_SPD Distribution",300,0,15,-0.5,299.5);
	fPtHFEMC_aftertoftpcnsigcut_SPD->Sumw2();
  	if(fIsMC) fOutputList->Add(fPtHFEMC_aftertoftpcnsigcut_SPD);

//------------------Centrality and Multiplicity------------------------------------------------- 
 
 
  
  

  //................................
  
  
  //fPteSPD_Lan = new TH2F("fPteSPD_Lan","Track multiplicity vs pt",300,0,15.,200,0,200.);
  //fPteSPD_Lan->Sumw2();
  //fOutputList->Add(fPteSPD_Lan);
  
  fPteSPD_Err = new TH2F("fPteSPD_Err","Track multiplicity vs pt",300,0,15.,300,-0.5,299.5);
  fPteSPD_Err->Sumw2();
  fOutputList->Add(fPteSPD_Err);
 
 //................................
 
  
  fPteSPD= new TH2F("fPteSPD","Track multiplicity vs pt",300,0,15.,300,-0.5,299.5);
	fPteSPD->Sumw2();
  fOutputList->Add(fPteSPD);
  	
  //................................
 
  
  //fHadSPDLan= new TH2F("fHadSPDLan","Track multiplicity vs pt",300,0,15.,200,0,200.);
	//fHadSPDLan->Sumw2();
  //fOutputList->Add(fHadSPDLan);
  
  fHadSPDErr= new TH2F("fHadSPDErr","Track multiplicity vs pt",300,0,15.,300,-0.5,299.5);
	fHadSPDErr->Sumw2();
  fOutputList->Add(fHadSPDErr);
  

//-------------------------------------------------------------------------------------------------------------
  
  Int_t nbinsTPCnsigmaVsPcut[4] = {100, 250, 300,1000};
  Double_t binlowTPCnsigmaVsPcut[4] = {0., -15, -0.5, 0.};
  Double_t binupTPCnsigmaVsPcut[4] = {10., 10, 299.5, 1000};

  fnSigmaVsP_TPC_cut_mult = new THnSparseF("fnSigmaVsP_TPC_cut_mult", "fnSigmaVsP_TPC_cut_mult;p;nsig;SPD;V0M", 4, nbinsTPCnsigmaVsPcut, binlowTPCnsigmaVsPcut, binupTPCnsigmaVsPcut);
  fnSigmaVsP_TPC_cut_mult->Sumw2();
  //fOutputList->Add(fnSigmaVsP_TPC_cut_mult);
  
  //-----------------NonHFE-----------------------------------------------------------------------------------------
  fPte_ULS = new TH1F("fPte_ULS", "ULS electron pt",300,0,15);
  fPte_ULS->Sumw2();
  fOutputList->Add(fPte_ULS);
  
  fPte_LS = new TH1F("fPte_LS", "LS electron pt",300,0,15);
  fPte_LS->Sumw2();
  fOutputList->Add(fPte_LS);
  
  fInvmassLS1 = new TH1F("fInvmassLS1","Inv mass of LS (e,e) for pt^{e}; mass(GeV/c^2); counts;",1000,0,1.0);
	fOutputList->Add(fInvmassLS1);
	
	fInvmassULS1 = new TH1F("fInvmassULS1","Inv mass of ULS (e,e) for pt^{e}; mass(GeV/c^2); counts;",1000,0,1.0);
	fOutputList->Add(fInvmassULS1);
	
	//....................................................................................
	
  
 
  fPte_ULS_multSPD = new TH2F("fPte_ULS_multSPD", "ULS electron pt in percentile",300,0,15,300,-0.5,299.5);
  fPte_ULS_multSPD->Sumw2();
  fOutputList->Add(fPte_ULS_multSPD);
  
  fPte_LS_multSPD = new TH2F("fPte_LS_multSPD", "LS electron pt in percentile",300,0,15,300,-0.5,299.5);
  fPte_LS_multSPD->Sumw2();
  fOutputList->Add(fPte_LS_multSPD);
  

  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  
 
  fPtHFEMC = new TH1F("fPtHFEMC",";p_{t} (GeV/c)",300,0,15);
  fPtHFEMC->Sumw2();
  if(fIsMC) fOutputList->Add(fPtHFEMC);
  
  fPtHFEMC_SPD = new TH2F("fPtHFEMC_SPD",";p_{t} (GeV/c)",300,0,15,300,-0.5,299.5);
  fPtHFEMC_SPD->Sumw2();
  if(fIsMC) fOutputList->Add(fPtHFEMC_SPD);
 

  
   
  pi0MC= new TH1F("pi0MC",";p_{t} (GeV/c)",300,0,15.);
  pi0MC->Sumw2();
  //fOutputList->Add(pi0MC);
  
  etaMC= new TH1F("etaMC",";p_{t} (GeV/c)",300,0,15.);
  etaMC->Sumw2();
  //fOutputList->Add(etaMC);
  
  gammaMC= new TH1F("gammaMC",";p_{t} (GeV/c)",300,0,15.);
  gammaMC->Sumw2();
  //fOutputList->Add(gammaMC);
  
  
  pi0MC_1= new TH1F("pi0MC_1",";p_{t} (GeV/c)",300,0,15.);
  pi0MC_1->Sumw2();
  //fOutputList->Add(pi0MC_1);
  
  etaMC_1= new TH1F("etaMC_1",";p_{t} (GeV/c)",300,0,15.);
  etaMC_1->Sumw2();
  //fOutputList->Add(etaMC_1);
  
  gammaMC_1= new TH1F("gammaMC_1",";p_{t} (GeV/c)",300,0,15.);
  gammaMC_1->Sumw2();
  //fOutputList->Add(gammaMC_1);
  
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
  
  fPt_elec_phot1= new TH1F("fPt_elec_phot1", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_phot1->Sumw2();
  if(fIsMC) fOutputList->Add(fPt_elec_phot1);
  
  fpt_e_nonphot_MC= new TH1F("fpt_e_nonphot_MC", "; p_{T}(GeV/c); counts;",300,0,15.);
  fpt_e_nonphot_MC->Sumw2();
  //fOutputList->Add(fpt_e_nonphot_MC);
  
  fElecNos = new TH1F("fElecNos","no of electrons",10,-0.5,9.5);
  fElecNos->GetXaxis()->SetBinLabel(1,"Total no. of Particles");
  fElecNos->GetXaxis()->SetBinLabel(2,"Inclusive Electrons");
  fElecNos->GetXaxis()->SetBinLabel(3,"Non-hfe Electrons");
  //fOutputList->Add(fElecNos);
 
  fPT_elec_MCtrue = new TH1F("fPT_elec_MCtrue","Pt of e from MC_True",300,0.,15.);
  fPT_elec_MCtrue->Sumw2();
  //fOutputList->Add(fPT_elec_MCtrue);
  
  fpt_incl_elec_MC= new TH1F("fpt_incl_elec_MC","Pt of inclusive e from MC",300,0.,15.);
  fpt_incl_elec_MC->Sumw2();
  //fOutputList->Add(fpt_incl_elec_MC);
  
  fPi0_Pt_MC = new TH1F("fPi0_Pt_MC","Pt of e mother Pi0 from MC",300,0.,15.);
  fPi0_Pt_MC->Sumw2();
  //fOutputList->Add(fPi0_Pt_MC);
  
  fEta_Pt_MC = new TH1F("fEta_Pt_MC","Pt of e mother Eta from MC",300,0.,15.);
  fEta_Pt_MC->Sumw2();
  //fOutputList->Add(fEta_Pt_MC);
  
  fGamma_Pt_MC = new TH1F("fGamma_Pt_MC","Pt of e mother Gamma from MC",300,0.,15.);
  fGamma_Pt_MC->Sumw2();
  //fOutputList->Add(fGamma_Pt_MC);
  
  fPte_ULS_MC = new TH1F("fPte_ULS_MC","Pt of ULS from MC",300,0.,15.);
  fPte_ULS_MC->Sumw2(); 
  if(fIsMC)fOutputList->Add(fPte_ULS_MC);
  
  fPte_LS_MC = new TH1F("fPte_LS_MC","Pt of LS from MC",300,0.,15.);
  fPte_LS_MC->Sumw2();
  if(fIsMC)fOutputList->Add(fPte_LS_MC);
 
  fPte_ULS_MC_multSPD = new TH2F("fPte_ULS_MC_multSPD","Pt of ULS from MC",300,0.,15.,300,-0.5,299.5);
  fPte_ULS_MC_multSPD->Sumw2(); 
  if(fIsMC)fOutputList->Add(fPte_ULS_MC_multSPD);
  
  fPte_LS_MC_multSPD = new TH2F("fPte_LS_MC_multSPD","Pt of LS from MC",300,0.,15.,300,-0.5,299.5);
  fPte_LS_MC_multSPD->Sumw2();
  if(fIsMC)fOutputList->Add(fPte_LS_MC_multSPD);
  
  
  fInvMULS = new TH2F("fInvMULS","Pt and InvMass from MC",300,0.,15.,1000,0,1.0);
  fInvMULS->Sumw2();
  if(fIsMC)fOutputList->Add(fInvMULS);
  

 
  fPt_elec_phot1_multSPD= new TH2F("fPt_elec_phot1_multSPD", "; p_{T}(GeV/c); counts;",300,0,15.,300,-0.5,299.5);
  fPt_elec_phot1_multSPD->Sumw2();
  if(fIsMC)fOutputList->Add(fPt_elec_phot1_multSPD);
  
  //----pi0---------------
  fPt_elec_pi0_MB= new TH1F("fPt_elec_pi0_MB", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_pi0_MB->Sumw2();
  //fOutputList->Add(fPt_elec_pi0_MB);
  
  fPt_elec_pi0_MB_LS= new TH1F("fPt_elec_pi0_MB_LS", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_pi0_MB_LS->Sumw2();
  //fOutputList->Add(fPt_elec_pi0_MB_LS);
  
  fPt_elec_pi0_MB_ULS= new TH1F("fPt_elec_pi0_MB_ULS", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_pi0_MB_ULS->Sumw2();
  //fOutputList->Add(fPt_elec_pi0_MB_ULS);
  
  //----eta---------------
  fPt_elec_eta_MB= new TH1F("fPt_elec_eta_MB", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_eta_MB->Sumw2();
  //fOutputList->Add(fPt_elec_eta_MB);
  
  fPt_elec_eta_MB_LS= new TH1F("fPt_elec_eta_MB_LS", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_eta_MB_LS->Sumw2();
  //fOutputList->Add(fPt_elec_eta_MB_LS);
  
  fPt_elec_eta_MB_ULS= new TH1F("fPt_elec_eta_MB_ULS", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_eta_MB_ULS->Sumw2();
  //fOutputList->Add(fPt_elec_eta_MB_ULS);
  
  
  //----gamma---------------
  fPt_elec_gamma_MB= new TH1F("fPt_elec_gamma_MB", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_gamma_MB->Sumw2();
  //fOutputList->Add(fPt_elec_gamma_MB);
  
  fPt_elec_gamma_MB_LS= new TH1F("fPt_elec_gamma_MB_LS", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_gamma_MB_LS->Sumw2();
  //fOutputList->Add(fPt_elec_gamma_MB_LS);
  
  fPt_elec_gamma_MB_ULS= new TH1F("fPt_elec_gamma_MB_ULS", "; p_{T}(GeV/c); counts;",300,0,15.);
  fPt_elec_gamma_MB_ULS->Sumw2();
  //fOutputList->Add(fPt_elec_gamma_MB_ULS); 
  
  //Weights============================================
  fPt_elec_phot2 = new TH1F("fPt_elec_phot2","",300,0,15.);
  fPt_elec_phot2->Sumw2();
	if(fIsMC)fOutputList->Add(fPt_elec_phot2);
	
   fPtElec_ULS_MC2 = new TH1F("fPtElec_ULS_MC2","p_{T} [GeV/c]",300,0,15.);
   fPtElec_ULS_MC2->Sumw2();
	if(fIsMC)fOutputList->Add(fPtElec_ULS_MC2);
		
	fPtElec_LS_MC2 = new TH1F("fPtElec_LS_MC2","p_{T} [GeV/c]",300,0,15.);
	fPtElec_LS_MC2->Sumw2();
   if(fIsMC)fOutputList->Add(fPtElec_LS_MC2);
	
	
	fPt_elec_phot2_multSPD= new TH2F("fPt_elec_phot2_multSPD","",300,0,15.,300,-0.5,299.5);
	fPt_elec_phot2_multSPD->Sumw2();
	if(fIsMC)fOutputList->Add(fPt_elec_phot2_multSPD);
	
	
   fPtElec_ULS_MC2_multSPD = new TH2F("fPtElec_ULS_MC2_multSPD","p_{T} [GeV/c]",300,0,15.,300,-0.5,299.5);
   fPtElec_ULS_MC2_multSPD->Sumw2();
	if(fIsMC)fOutputList->Add(fPtElec_ULS_MC2_multSPD);
	
   fPtElec_LS_MC2_multSPD = new TH2F("fPtElec_LS_MC2_multSPD","p_{T} [GeV/c]",300,0,15.,300,-0.5,299.5);
   fPtElec_LS_MC2_multSPD->Sumw2();
    if(fIsMC)fOutputList->Add(fPtElec_LS_MC2_multSPD);
	
	
	/*
	fPt_elec_from_pi02 = new TH1F("fPt_elec_from_pi02","",300,0,15.);
	fOutputList->Add(fPt_elec_from_pi02);
	//fPt_elec_from_pi02->Sumw2();
	
	fPtElec_ULS_MC_from_pi02 = new TH1F("fPtElec_ULS_MC_from_pi02","p_{T} [GeV/c]",300,0,15.);
	fOutputList->Add(fPtElec_ULS_MC_from_pi02);
	//fPtElec_ULS_MC_from_pi02->Sumw2();
	
	fPtElec_LS_MC_from_pi02 = new TH1F("fPtElec_LS_MC_from_pi02","p_{T} [GeV/c]",300,0,15.);
    fOutputList->Add(fPtElec_LS_MC_from_pi02);
	//fPtElec_LS_MC_from_pi02->Sumw2();
	
	fPt_elec_from_eta2 = new TH1F("fPt_elec_from_eta2","",300,0,15.);
	fOutputList->Add(fPt_elec_from_eta2);
	fPt_elec_from_eta2->Sumw2();
	
	fPtElec_ULS_MC_from_eta2 = new TH1F("fPtElec_ULS_MC_from_eta2","p_{T} [GeV/c]",300,0,15.);
	fOutputList->Add(fPtElec_ULS_MC_from_eta2);
	fPtElec_ULS_MC_from_eta2->Sumw2();
	
	fPtElec_LS_MC_from_eta2 = new TH1F("fPtElec_LS_MC_from_eta2","p_{T} [GeV/c]",300,0,15.);
    fOutputList->Add(fPtElec_LS_MC_from_eta2);
	fPtElec_LS_MC_from_eta2->Sumw2();
	
	fPt_elec_from_gamma2 = new TH1F("fPt_elec_from_gamma2","",300,0,15.);
	fOutputList->Add(fPt_elec_from_gamma2);
	fPt_elec_from_gamma2->Sumw2();
	
	fPtElec_ULS_MC_from_gamma2 = new TH1F("fPtElec_ULS_MC_from_gamma2","p_{T} [GeV/c]",300,0,15.);
	fOutputList->Add(fPtElec_ULS_MC_from_gamma2);
	fPtElec_ULS_MC_from_gamma2->Sumw2();
	
	fPtElec_LS_MC_from_gamma2 = new TH1F("fPtElec_LS_MC_from_gamma2","p_{T} [GeV/c]",300,0,15.);
    fOutputList->Add(fPtElec_LS_MC_from_gamma2);
	fPtElec_LS_MC_from_gamma2->Sumw2();*/
//===========================================================================================================
  
  //-------------------------------------------------------------------------------------------------------------
  
  Int_t nbinsInvMULS[4] = {100, 6, 100, 300};
  Double_t binlowInvMULS[4] = {0., 0, 0., -0.5};
  Double_t binupInvMULS[4] = {10, 6., 1., 299.5};

	fInvMULSnSp = new THnSparseF("fInvMULSnSp", "fInvMULSnSp;pt;source;mass;multdepSPD", 4, nbinsInvMULS, binlowInvMULS, binupInvMULS);
	fInvMULSnSp->Sumw2();
	if(fIsMC)fOutputList->Add(fInvMULSnSp);
	
	fInvMULSnSpwg = new THnSparseF("fInvMULSnSpwg", "fInvMULSnSpwg;pt;source;mass;multdepSPD", 4, nbinsInvMULS, binlowInvMULS, binupInvMULS);
	fInvMULSnSpwg->Sumw2();
	if(fIsMC)fOutputList->Add(fInvMULSnSpwg);
	
	/*
  fInvMULSpi0 = new TH2F("fInvMULSpi0","Pt and InvMass from MC",300,0.,15.,100,0,1.0);
  fInvMULSpi0->Sumw2();
  fOutputList->Add(fInvMULSpi0);
  
  fInvMULSeta = new TH2F("fInvMULSeta","Pt and InvMass from MC",300,0.,15.,100,0,1.0);
  fInvMULSeta->Sumw2();
  fOutputList->Add(fInvMULSeta);
  
  fInvMULSgamma = new TH2F("fInvMULSgamma","Pt and InvMass from MC",300,0.,15.,100,0,1.0);
  fInvMULSgamma->Sumw2();
  fOutputList->Add(fInvMULSgamma);
  
  */

  //--------------------------------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------------------------------------
  fPtHFEMC_reco = new TH1F("fPtHFEMC_reco",";p_{t} (GeV/c)",300,0.,15.0);
  fPtHFEMC_reco->Sumw2();
  if(fIsMC)fOutputList->Add(fPtHFEMC_reco);
  
  fPtHFEMC_reco_SPD = new TH2F("fPtHFEMC_reco_SPD",";p_{t} (GeV/c)",300,0.,15.0,300,-0.5,299.5);
  fPtHFEMC_reco_SPD->Sumw2();
  if(fIsMC)fOutputList->Add(fPtHFEMC_reco_SPD);
  

 //--------------------------------------------------------------------------------------------------
 //--------SPD correction---------------------------------
  
  const char *estimatorName="tr"; //SPD  
  
  fSPD_tracklet = new TH1F("fSPD_tracklet", "fSPD_tracklet Distribution",300,-0.5,299.5);
  fSPD_tracklet->Sumw2(); 
  fOutputList->Add(fSPD_tracklet); 
  
  fSPDCorrMultDist_min = new TH1F("fSPDCorrMultDist_min",Form("Corrected Mult Dist;  N_{%s};Counts;",estimatorName),300,-0.5,299.5); //
  fSPDCorrMultDist_min->Sumw2();
  fOutputList->Add(fSPDCorrMultDist_min);
  
  fSPDCorrMultDist_max = new TH1F("fSPDCorrMultDist_max",Form("Corrected Mult Dist; N_{%s};Counts;",estimatorName),300,-0.5,299.5); //
  fSPDCorrMultDist_max->Sumw2();
  fOutputList->Add(fSPDCorrMultDist_max);
  
  fSPDWeightedCorrMultDist_max= new TH1F("fSPDWeightedCorrMultDist_max",Form("Corrected Mult Dist; N_{%s};Counts;",estimatorName),300,-0.5,299.5); //
  fSPDWeightedCorrMultDist_max->Sumw2();
  if(fIsMC)fOutputList->Add(fSPDWeightedCorrMultDist_max);
 
  fHistNtrVsZvtx = new TH2F("hNtrVsZvtx",Form("N%s vs VtxZ; Z_{vtx};N_{%s};",estimatorName,estimatorName),300,-15.0,15.0,300,-0.5,299.5); //
  fHistNtrVsZvtx->Sumw2();
  fOutputList->Add(fHistNtrVsZvtx);
  
  fHistNtrCorrVsZvtx_min = new TH2F("hNtrCorrVsZvtx_min",Form("N%s vs VtxZ; Z_{vtx};N_{%s};",estimatorName,estimatorName),300,-15.0,15.0,300,-0.5,299.5); //
  fHistNtrCorrVsZvtx_min->Sumw2();
  fOutputList->Add(fHistNtrCorrVsZvtx_min);
  					
  fHistNtrCorrVsZvtx_max = new TH2F("hNtrCorrVsZvtx_max",Form("N%s vs VtxZ; Z_{vtx};N_{%s};",estimatorName,estimatorName),300,-15.0,15.0,300,-0.5,299.5); //
  fHistNtrCorrVsZvtx_max->Sumw2();
  fOutputList->Add(fHistNtrCorrVsZvtx_max);

  Profile_Mean = new TProfile("Profile_Mean","Profile_Mean",300,-15.0,15.0);
  Profile_Mean->Sumw2();
  fOutputList->Add(Profile_Mean); 
  
  Profile_MeanCorr = new TProfile("Profile_MeanCorr","Profile_MeanCorr",300,-15.0,15.0);
  Profile_MeanCorr->Sumw2();
  fOutputList->Add(Profile_MeanCorr);
 
  fNchVsZvtx= new TH2F("fNchVsZvtx","Corrected Mult Dist;Nch;Counts;",300,-15.0,15.0,400,-0.5,399.5); //
  fNchVsZvtx->Sumw2();
  if(fIsMC)fOutputList->Add(fNchVsZvtx);
  
  SPDNtrCorrVsNch = new TH2F("SPDNtrCorrVsNch","",300,-0.5,299.5,400,-0.5,399.5); //
  SPDNtrCorrVsNch->Sumw2();
  if(fIsMC)fOutputList->Add(SPDNtrCorrVsNch);	
  
  SPDNtrCorrVsNchweights= new TH2F("SPDNtrCorrVsNchweights","",300,-0.5,299.5,400,-0.5,399.5); //
  SPDNtrCorrVsNchweights->Sumw2();
  if(fIsMC)fOutputList->Add(SPDNtrCorrVsNchweights);
   

  Int_t nbinsSPDV0MNch[3] = {200, 1000,400};
  Double_t binlowSPDV0MNch[3] = {0., 0, 0.};
  Double_t binupSPDV0MNch[3] = {200., 1000.,400.};

  fSPDNtrCorrVsV0MCorrVsNch = new THnSparseF("fSPDNtrCorrVsV0MCorrVsNch", "fSPDNtrCorrVsV0MCorrVsNch;SPD;V0M;Nch", 3, nbinsSPDV0MNch,binlowSPDV0MNch,binupSPDV0MNch);
  fSPDNtrCorrVsV0MCorrVsNch->Sumw2();
  //fOutputList->Add(fSPDNtrCorrVsV0MCorrVsNch);

  
  	
	// ----- weights for tagging efficiency -----
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
2.28928,2.58223,2.91267,3.2854,3.70582,4.18004,4.71494,5.3183,5.99886,6.76651,7.6324,8.60909,9.71076,10.9534,12.3551,13.9361,15.7195,17.731,20};
    
    const Int_t nDima=4;
    Int_t nBina[nDima] = {44,nBinspdg,nBinsg,nBinstype};
    fPi0EtaSpectra = new THnSparseF("fPi0EtaSpectra","fPi0EtaSpectra",nDima,nBina);
    fPi0EtaSpectra->SetBinEdges(0,binLimpt);
    fPi0EtaSpectra->SetBinEdges(1,binLimpdg);
    fPi0EtaSpectra->SetBinEdges(2,binLimg);
    fPi0EtaSpectra->SetBinEdges(3,binLimtype);
    fPi0EtaSpectra->Sumw2();
    //fOutputList->Add(fPi0EtaSpectra);
  
  
  
  
  Int_t nbinspt[3] = {1000, 3, 7};
  Double_t binlow[3] = {0., 0, -1.};
  Double_t binup[3] = {10, 3, 6.};
	fPi0EtaSpectraSp = new THnSparseF("fPi0EtaSpectraSp", "fPi0EtaSpectraSp;pt;source;type", 3, nbinspt, binlow,binup);
	fPi0EtaSpectraSp->Sumw2();
   //fOutputList->Add(fPi0EtaSpectraSp);
   
 //===================================================================================================//
 
 
  fNentries2=new TH1F("CutSet", "", 19,-0.5,18.5);
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
   fOutputList->Add(fNentries2);
  
  fNentries2->SetBinContent(1,ftrigger);
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
   
  //==================================================================================================//
  //==================================================================================================//
  					
  const char* nameoutput=GetOutputSlot(2)->GetContainer()->GetName();  // gets the output slot from mgr-> in AddTask
  fNentries=new TH1F(nameoutput, "", 5,-0.5,4.5);
  fNentries->GetXaxis()->SetBinLabel(1,"No. of Events Analyzed");
  fNentries->GetXaxis()->SetBinLabel(2,"No. of Events Accepted");
  fNentries->GetXaxis()->SetBinLabel(3,"After Pile Up cuts");
  fNentries->GetXaxis()->SetBinLabel(4,"After Z vtx cut");
  fNentries->GetXaxis()->SetBinLabel(5,"After N contributors");
  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);
  fNentries->Sumw2();
  fNentries->SetMinimum(0);
  
  fNentries2->GetXaxis()->SetNdivisions(1,kFALSE);
  fNentries2->Sumw2();
  fNentries2->SetMinimum(0);
 
  PostData(1,fOutputList);
  PostData(2,fNentries);

}

//////////////=================================================///////////////////////

void AliAnalysisTaskHFEmultTPCTOF::UserExec(Option_t *) 
{
	fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
	if (!fAOD) { return;}
  
	fNentries->Fill(0);                //No. of Events Analyzed 
	Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & ftrigger);
  
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

/*	//For Pb-Pb ---------------
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
  
  
   /*//====DhananJaya's Code
   const AliVVertex* spdVtx = fAOD->GetPrimaryVertexSPD();
   fVtxZ->Fill(spdVtx->GetZ());
    if (!spdVtx || spdVtx->GetNContributors()<=0) return;
    Float_t zvtx = spdVtx->GetZ();
    Double_t Zvertex= spdVtx->GetZ();
    if (spdVtx->GetNContributors()<=0) return;
    Double_t cov[6]={0};
    spdVtx->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);

    if (zRes>0.25) return;
    //if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;
    //if (TMath::Abs(zvtx) > 10) return;
    if (TMath::Abs(spdVtx->GetZ()) > 10) return;
*/

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

	
	
	fpidResponse = fInputHandler->GetPIDResponse();
	if(!fpidResponse)
	{
		AliDebug(1, "Using default PID Response");
		fpidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
		fPID->SetPIDResponse(fpidResponse);
	}
	fPID->SetPIDResponse(fpidResponse);
	/*
	//===============V0 Multiplicity===========================
	AliAODVZERO *vzeroAOD = dynamic_cast<AliAODVZERO *>( dynamic_cast<AliAODEvent *>(fAOD)->GetVZEROData());
	Int_t V0AMult = static_cast<Int_t>(vzeroAOD->GetMTotV0A());
	Int_t V0CMult = static_cast<Int_t>(vzeroAOD->GetMTotV0C());
	fV0AMult->Fill(V0AMult);
	fV0CMult->Fill(V0CMult);
	Int_t V0Mult=V0AMult+V0CMult;


	//V0M Correction
  Int_t V0AMultCorr=V0AMult, V0CMultCorr=V0CMult, V0MultCorr=V0Mult;
  V0AMultCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(V0AMult,Zvertex));
  V0CMultCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0C(V0CMult,Zvertex));
 	 V0MultCorr = V0AMultCorr + V0CMultCorr;

	Double_t V0Mmult = V0MultCorr;
 */
 Double_t V0Mmult=1;
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

	fSPD_tracklet->Fill(nAcceta);
 	

  	
	Double_t counteta1Corr_min=nAcceta;
	Double_t N_corr_tr_eta1=nAcceta;

  Double_t fRefMult1_min=9.016;

  
  TProfile* estimatorAvg = GetEstimatorHistogram(fAOD);
  if(estimatorAvg){
       counteta1Corr_min=static_cast<Int_t>(GetCorrectedNtracklets(estimatorAvg,nAcceta,Zvertex,fRefMult1_min));
	     N_corr_tr_eta1=static_cast<Int_t>(GetCorrectedNtracklets(estimatorAvg,nAcceta,Zvertex,fRefMult)); //48.75 fRefMult
  } 
  
	
	
  Double_t SPDntr = N_corr_tr_eta1;
  
 
  

  fSPDCorrMultDist_min->Fill(counteta1Corr_min);
  fSPDCorrMultDist_max->Fill(SPDntr);
  	
  fHistNtrVsZvtx->Fill(Zvertex,nAcceta);
  Profile_Mean->Fill(Zvertex,nAcceta);
  fHistNtrCorrVsZvtx_min->Fill(Zvertex,counteta1Corr_min);
  fHistNtrCorrVsZvtx_max->Fill(Zvertex,SPDntr);
  Profile_MeanCorr->Fill(Zvertex,SPDntr);
  
   
		
	
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

			///For the reconstruction efficiency of HFE:-----------
			if(TMath::Abs(fMCparticle->Eta()) <=fEtarange  )
			{					
				if(TMath::Abs(pdg) == 11)
				{		
					Int_t IsElecHf=GetHFE(fMCparticle,fMCarray);
					if((IsElecHf==kBeauty) || (IsElecHf==kCharm)){
						fPtHFEMC->Fill(fMCparticle->Pt());
						fPtHFEMC_SPD->Fill(fMCparticle->Pt(),SPDntr);
					
					}
				}
			} //eta <0.8 condition
			
			//pt spectra for pi0 and eta
			//if(TMath::Abs(fMCparticle->Y())>0.5 ) continue;
			if(fMCparticle->Y()<-0.8 || fMCparticle->Y()>0.8) continue;
				
			
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
					
				//- Fill
				//if(qaweights[1]>0.) fPi0EtaSpectra->Fill(qaweights);
				//if(pi0etaweights[1]>0.) fPi0EtaSpectraSp->Fill(pi0etaweights);
				
			
				///-----------------------------------------------------
		
		       			
					
					
		////====================================================
		} //For loop
	}


   if(fIsMC) 
	{
		Double_t weigtsMC=1.;
		Double_t SPDntrnoweights=SPDntr;
		if(SPDntr>40 && SPDntr<=115) weigtsMC=func_MCnchCorr->Eval(SPDntr);
		else if (SPDntr<=40) weigtsMC=WeightMCncCorr(SPDntr);
	 
	  fSPDWeightedCorrMultDist_max->Fill(SPDntr,weigtsMC);
	  Int_t nch = GetNcharged();
		/*Double_t SPDV0MNch[3];
		SPDV0MNch[0]=SPDntr;
		SPDV0MNch[1]=V0Mmult;
		SPDV0MNch[2]=nch;
      */
		SPDNtrCorrVsNch->Fill(SPDntr,nch);
		SPDNtrCorrVsNchweights->Fill(SPDntr,nch,weigtsMC);
		//fSPDNtrCorrVsV0MCorrVsNch->Fill(SPDV0MNch);
		fNchVsZvtx->Fill(Zvertex,nch);
	}
	
	
	
	Double_t weight_eta=0., weight_pi0=0.;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////======================================TRACK LOOP==================================================/////
//////==================================================================================================/////

	for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) 
  {
		AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
   //cout<<" iTracks "<<iTracks<<endl;
    if(!track) AliFatal("Not a standard AOD");
  		
    Int_t tracktypeTrig=ClassifyTrack(track,pVtx);  //track classify
    
    fDCAZvsPt->Fill(track->Pt(),fDCAz);
	 fDCAXYvsPt->Fill(track->Pt(),fDCAxy);
    	
    if(tracktypeTrig!=1) continue;  //if(tracktype==0) continue; if(tracktype==1)  	
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
		
	
		num++;
    	
		Double_t pt=0., p=0., eta=-999, dEdx_TPC=-999., Signal_TOF=-999., TOFbeta=-999., fTPCnSigma=-999.0, fTOFnSigma=-999.0;    
    
		pt = track->Pt();
		p = track->P();  
		eta=track->Eta();
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

		if(fIsMC && track->GetLabel()>=0)
     	{		
			fMCparticle=(AliAODMCParticle*)fMCarray->At(track->GetLabel());
			pdg = fMCparticle->GetPdgCode();
			if(TMath::Abs(pdg) == 11)
			{													
     			Int_t IsElecHf=GetHFE(fMCparticle,fMCarray);
				if((IsElecHf==kBeauty) || (IsElecHf==kCharm))
				{
					fPtHFEMC_aftertrackcut->Fill(pt);
					fPtHFEMC_aftertrackcut_SPD->Fill(pt,SPDntr);
					
			
					if(TMath::Abs(fTOFnSigma) < fTOFnsig){
						fPtHFEMC_aftertofcut->Fill(pt);
						fPtHFEMC_aftertofcut_SPD->Fill(pt,SPDntr);
					}
					if(fTPCnSigma>fTPCnsigmin && fTPCnSigma<fTPCnsigmax) {
						fPtHFEMC_aftertpcnsigcut->Fill(pt);
						fPtHFEMC_aftertpcnsigcut_SPD->Fill(pt,SPDntr);
					}
					if(TMath::Abs(fTOFnSigma) < fTOFnsig){
						if(fTPCnSigma>fTPCnsigmin && fTPCnSigma<fTPCnsigmax) {
							fPtHFEMC_aftertoftpcnsigcut->Fill(pt);
							fPtHFEMC_aftertoftpcnsigcut_SPD->Fill(pt,SPDntr);
						}
					}
				}
			}	
		}


		if(TMath::Abs(fTOFnSigma) < fTOFnsig) 
		{
			fdEdxVsP_TPC_cut->Fill(p,dEdx_TPC);
			fnSigmaVsP_TPC_cut->Fill(p,fTPCnSigma);
			fnSigmaVsPt_TPC_cut->Fill(pt,fTPCnSigma);
			fnSigmaVsEta_TPC_cut->Fill(eta,fTPCnSigma);
       			
			Double_t TPCnsigmult[4];
    	  		TPCnsigmult[0]=p;
		    	TPCnsigmult[1]=fTPCnSigma;
        		TPCnsigmult[2]=SPDntr;
        		//TPCnsigmult[3]=V0Mmult;
			//fnSigmaVsP_TPC_cut_mult->Fill(TPCnsigmult);

	       	
	
			weight_lan=fLandau->Eval(p);
			weight_err=fErr->Eval(p);
			weight_lan_inv=1-weight_lan;
			weight_err_inv=1-weight_err;
       		       		
			if(fTPCnSigma>fTPCnsigmin && fTPCnSigma<fTPCnsigmax) 
			{ 		
				//fHadCot_Landau->Fill(pt,weight_lan);
				fHadCot_Err->Fill(pt,weight_err);
				//fPt_incl_e_Landau->Fill(pt,weight_lan_inv);
				fPt_incl_e_Err->Fill(pt,weight_err_inv);
  						
				//---------Multiplicity----------------------
			
				//fPteSPD_Lan->Fill(pt,SPDntr,weight_lan_inv);
				fPteSPD_Err->Fill(pt,SPDntr,weight_err_inv);
				fPteSPD->Fill(pt,SPDntr);
				//fHadSPDLan->Fill(pt,SPDntr,weight_lan);
 				fHadSPDErr->Fill(pt,SPDntr,weight_err);	
  				
			
  				
				fNonHFE = new AliSelectNonHFE();
				fNonHFE->SetAODanalysis(kTRUE);
				fNonHFE->SetInvariantMassCut(fInvmassCut);
				fNonHFE->SetAlgorithm("DCA"); //KF,DCA
				fNonHFE->SetPIDresponse(fpidResponse);
				fNonHFE->SetTrackCuts(-1*fAssoTPCnsig,fAssoTPCnsig); //TPCnsigma cuts
				//fNonHFE->SetChi2OverNDFCut(4.0);
				//fNonHFE->SetDCACut(fDCAcut);
				//fNonHFE-> SetDCAPartnerCuts(1,2);   //SetDCAPartnerCuts(Double_t xy, Double_t z)
				fNonHFE->SetAdditionalCuts(fAssopTMin,fAssoTPCCluster);  //
				/*fNonHFE->SetHistMassBack(fInvmassLS1);
				fNonHFE->SetHistMass(fInvmassULS1);
					if(pt>=0.5 && pt<4.5){
				fNonHFE->SetHistMassBack(fInvmassLS1);
				fNonHFE->SetHistMass(fInvmassULS1);
				}*/
				if(pt>=0.5 && pt<1.5){
				  
				fNonHFE->SetHistMassBack(fInvmassLS1);
				fNonHFE->SetHistMass(fInvmassULS1);
				}
				/*
				else if(pt>=1.5 && pt<2.5){
			
				fNonHFE->SetHistMassBack(fInvmassLS1_2);
				fNonHFE->SetHistMass(fInvmassULS1_2);
				}
				else if(pt>=2.5 && pt<3.5){
				fNonHFE->SetHistMassBack(fInvmassLS1_3);
				fNonHFE->SetHistMass(fInvmassULS1_3);
				}
				else if(pt>=3.5 && pt<4.5){
				fNonHFE->SetHistMassBack(fInvmassLS1_4);
				fNonHFE->SetHistMass(fInvmassULS1_4);
				}*/
				fNonHFE->FindNonHFE(iTracks,track,fAOD);
			
				Int_t fNULS = fNonHFE->GetNULS();
				Int_t fNLS = fNonHFE->GetNLS();
			 
				if(fNonHFE->IsULS())
				{
					fPte_ULS->Fill(track->Pt(),fNULS);
					fPte_ULS_multSPD->Fill(track->Pt(),SPDntr,fNULS);
				}
				if(fNonHFE->IsLS())
				{
					fPte_LS->Fill(track->Pt(),fNLS);
					fPte_LS_multSPD->Fill(track->Pt(),SPDntr,fNLS);
				}
  
				if(fIsMC && track->GetLabel()>=0)
   				{
				
					fElecNos->Fill(0);   						
					fMCparticle=(AliAODMCParticle*)fMCarray->At(track->GetLabel());
					pdg = fMCparticle->GetPdgCode();
					Float_t ptMC= fMCparticle->Pt();
					//fPT_elec_MCtrue->Fill(ptMC);
					
					if(TMath::Abs(pdg)!=11) continue ;      //Is electron:
					
					//fpt_incl_elec_MC->Fill(pt);
					fElecNos->Fill(1);
								
					Int_t fMCmotherindex=fMCparticle->GetMother(); //Getting Electron Mother
					if(fMCmotherindex<0) continue ;
					
					fMCmother= (AliAODMCParticle*)fMCarray->At(fMCparticle->GetMother());
					Int_t pdg_mother = fMCmother->GetPdgCode();
								
					//----If mother is Pi0 , eta or gamma --------------------------------------------------------------
					if(TMath::Abs(pdg_mother) == 111 || TMath::Abs(pdg_mother) == 22 || TMath::Abs(pdg_mother) == 221)
					{
						fElecNos->Fill(2);
							
						if(TMath::Abs(pdg_mother) == 111) fPi0_Pt_MC->Fill(fMCmother->Pt());
						if(TMath::Abs(pdg_mother) == 221) fEta_Pt_MC->Fill(fMCmother->Pt());
						if(TMath::Abs(pdg_mother) == 22) fGamma_Pt_MC->Fill(fMCmother->Pt());
							
						Double_t ptmotherw = -1.;	
						Int_t elec_source = GetElecSourceType(fMCparticle,fMCarray,ptmotherw);
										
						if((elec_source==kPi0NoFeedDown) || (elec_source==kGPi0NoFeedDown) || (elec_source==kEtaNoFeedDown) || (elec_source==kGEtaNoFeedDown))
						{						
							fPt_elec_phot1->Fill(pt);
							fPt_elec_phot1_multSPD->Fill(pt,SPDntr);
												
							SelectPhotonicElectronR(iTracks, track, fMCmotherindex, pdg,elec_source, V0Mmult, SPDntr,ptmotherw);
		
						  fNonHFE->FindNonHFE(iTracks,track,fAOD);
						  		
						  Int_t fNULS_MC = fNonHFE->GetNULS();
						  Int_t fNLS_MC = fNonHFE->GetNLS();
						 	
						 	
						  if(fNonHFE->IsULS())
						  {
								fPte_ULS_MC->Fill(track->Pt(),fNULS_MC);
										fPte_ULS_MC_multSPD->Fill(track->Pt(),SPDntr,fNULS_MC);
							}

						  if(fNonHFE->IsLS())
						  {
								fPte_LS_MC->Fill(track->Pt(),fNLS_MC);
											fPte_LS_MC_multSPD->Fill(track->Pt(),SPDntr,fNLS_MC);
							} 
							
						  if((elec_source==kPi0NoFeedDown)){
								/*fPt_elec_pi0_MB->Fill(pt);  
								if(fNonHFE->IsULS()) fPt_elec_pi0_MB_ULS->Fill(pt,fNULS_MC);
								if(fNonHFE->IsLS()) fPt_elec_pi0_MB_LS->Fill(pt,fNLS_MC);
								*/
								
								weight_pi0=0;
								if(ptmotherw >= 0.5 && ptmotherw <= 2.5) weight_pi0 = f0a->Eval(ptmotherw);
								if(ptmotherw >2.5 && ptmotherw < 15) weight_pi0 = f0b->Eval(ptmotherw);
								
								/*fPt_elec_from_pi02->Fill(pt,weight_pi0);
								if(fNonHFE->IsULS()) fPtElec_ULS_MC_from_pi02->Fill(pt,fNonHFE->GetNULS()*weight_pi0);
								if(fNonHFE->IsLS()) fPtElec_LS_MC_from_pi02->Fill(pt,fNonHFE->GetNLS()*weight_pi0);
								*/
								fPt_elec_phot2->Fill(pt,weight_pi0);
								fPt_elec_phot2_multSPD->Fill(pt,SPDntr,weight_pi0);
						
							  
								if(fNonHFE->IsULS()){ fPtElec_ULS_MC2->Fill(pt,fNonHFE->GetNULS()*weight_pi0);
								   fPtElec_ULS_MC2_multSPD->Fill(pt,SPDntr,fNonHFE->GetNULS()*weight_pi0);
								}
								if(fNonHFE->IsLS()){ fPtElec_LS_MC2->Fill(pt,fNonHFE->GetNLS()*weight_pi0);
								  fPtElec_LS_MC2_multSPD->Fill(pt,SPDntr,fNonHFE->GetNLS()*weight_pi0);
								}
								
							}
							if((elec_source==kEtaNoFeedDown)){
								/*fPt_elec_eta_MB->Fill(pt);  
								if(fNonHFE->IsULS()) fPt_elec_eta_MB_ULS->Fill(pt,fNULS_MC);
								if(fNonHFE->IsLS()) fPt_elec_eta_MB_LS->Fill(pt,fNLS_MC);
								*/
								weight_eta=0;
								if(ptmotherw >= 0.5 && ptmotherw < 2.3) weight_eta = f1->Eval(ptmotherw);
								if(ptmotherw >= 2.3 && ptmotherw < 4.1) weight_eta = f1a->Eval(ptmotherw);
								if(ptmotherw >= 4.1 && ptmotherw < 15) weight_eta = f1b->Eval(ptmotherw);
								
								/*fPt_elec_from_eta2->Fill(pt,weight_eta);
								if(fNonHFE->IsULS()) fPtElec_ULS_MC_from_eta2->Fill(pt,fNonHFE->GetNULS()*weight_eta);
								if(fNonHFE->IsLS()) fPtElec_LS_MC_from_eta2->Fill(pt,fNonHFE->GetNLS()*weight_eta);
								*/
								fPt_elec_phot2->Fill(pt,weight_eta);
								fPt_elec_phot2_multSPD->Fill(pt,SPDntr,weight_eta);
							
							    
								if(fNonHFE->IsULS()){ fPtElec_ULS_MC2->Fill(pt,fNonHFE->GetNULS()*weight_eta);
								   fPtElec_ULS_MC2_multSPD->Fill(pt,SPDntr,fNonHFE->GetNULS()*weight_eta);
								}
								if(fNonHFE->IsLS()){ fPtElec_LS_MC2->Fill(pt,fNonHFE->GetNLS()*weight_eta);
								  fPtElec_LS_MC2_multSPD->Fill(pt,SPDntr,fNonHFE->GetNLS()*weight_eta);
								}
								
							}
							if((elec_source==kGPi0NoFeedDown)){
								/*fPt_elec_gamma_MB->Fill(pt);  
								if(fNonHFE->IsULS()) fPt_elec_gamma_MB_ULS->Fill(pt,fNULS_MC);
								if(fNonHFE->IsLS()) fPt_elec_gamma_MB_LS->Fill(pt,fNLS_MC);
								*/
								weight_pi0=0;
								if(ptmotherw >= 0.5 && ptmotherw <= 2.5) weight_pi0 = f0a->Eval(ptmotherw);
								if(ptmotherw >2.5 && ptmotherw < 15) weight_pi0 = f0b->Eval(ptmotherw);
								
								/*fPt_elec_from_gamma2->Fill(pt,weight_pi0);
								if(fNonHFE->IsULS()) fPtElec_ULS_MC_from_gamma2->Fill(pt,fNonHFE->GetNULS()*weight_pi0);
								if(fNonHFE->IsLS()) fPtElec_LS_MC_from_gamma2->Fill(pt,fNonHFE->GetNLS()*weight_pi0);
								*/
								fPt_elec_phot2->Fill(pt,weight_pi0); 
								fPt_elec_phot2_multSPD->Fill(pt,SPDntr,weight_pi0);
							
								if(fNonHFE->IsULS()){ fPtElec_ULS_MC2->Fill(pt,fNonHFE->GetNULS()*weight_pi0);
								   fPtElec_ULS_MC2_multSPD->Fill(pt,SPDntr,fNonHFE->GetNULS()*weight_pi0);
								}
								if(fNonHFE->IsLS()){ fPtElec_LS_MC2->Fill(pt,fNonHFE->GetNLS()*weight_pi0);
								  fPtElec_LS_MC2_multSPD->Fill(pt,SPDntr,fNonHFE->GetNLS()*weight_pi0);
								}
							}
							if((elec_source==kGEtaNoFeedDown)){
								/*fPt_elec_gamma_MB->Fill(pt);  
								if(fNonHFE->IsULS()) fPt_elec_gamma_MB_ULS->Fill(pt,fNULS_MC);
								if(fNonHFE->IsLS()) fPt_elec_gamma_MB_LS->Fill(pt,fNLS_MC);
								*/
								weight_eta=0;
								if(ptmotherw >= 0.5 && ptmotherw < 2.3) weight_eta = f1->Eval(ptmotherw);
								if(ptmotherw >= 2.3 && ptmotherw < 4.1) weight_eta = f1a->Eval(ptmotherw);
								if(ptmotherw >= 4.1 && ptmotherw < 15) weight_eta = f1b->Eval(ptmotherw);
								/*fPt_elec_from_gamma2->Fill(pt,weight_eta);
								if(fNonHFE->IsULS()) fPtElec_ULS_MC_from_gamma2->Fill(pt,fNonHFE->GetNULS()*weight_eta);
								if(fNonHFE->IsLS()) fPtElec_LS_MC_from_gamma2->Fill(pt,fNonHFE->GetNLS()*weight_eta);
								*/
								fPt_elec_phot2->Fill(pt,weight_eta);
								fPt_elec_phot2_multSPD->Fill(pt,SPDntr,weight_eta);
						
							    
								if(fNonHFE->IsULS()){ fPtElec_ULS_MC2->Fill(pt,fNonHFE->GetNULS()*weight_eta);
								   fPtElec_ULS_MC2_multSPD->Fill(pt,SPDntr,fNonHFE->GetNULS()*weight_eta);
								}
								if(fNonHFE->IsLS()){ fPtElec_LS_MC2->Fill(pt,fNonHFE->GetNLS()*weight_eta);
								  fPtElec_LS_MC2_multSPD->Fill(pt,SPDntr,fNonHFE->GetNLS()*weight_eta);
								}  
							}	
						}			
					}	
				//	else fpt_e_nonphot_MC->Fill(track->Pt());	
							
					Int_t IsElecHf=GetHFE(fMCparticle,fMCarray);
					if((IsElecHf==kBeauty) || (IsElecHf==kCharm))
					{
						fPtHFEMC_reco->Fill(pt);
						fPtHFEMC_reco_SPD->Fill(pt,SPDntr);
						
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

void AliAnalysisTaskHFEmultTPCTOF::SelectPhotonicElectronR(Int_t itrack, AliAODTrack *track, Int_t motherindex, Int_t pdg1, Int_t source,Double_t V0Mmult1 , Double_t SPDntr1,Double_t ptmotherwg)
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
        
        Double_t invms[5];
    	  invms[0]=track->Pt();
	     invms[1]=source;
        invms[2]=mass;
        invms[3]=SPDntr1;
	    // invms[4]=V0Mmult1;
        
							
        if(fFlagULS){
        	fInvMULS->Fill(track->Pt(),mass);
        	//fInvMULS_multV0M->Fill(track->Pt(),mass,V0Mmult1);
        	//fInvMULS_multSPD->Fill(track->Pt(),mass,SPDntr1);
        	
        	fInvMULSnSp->Fill(invms);
			/*if(source==kPi0NoFeedDown)fInvMULSpi0->Fill(track->Pt(),mass);
			if(source==kEtaNoFeedDown)fInvMULSeta->Fill(track->Pt(),mass);
			if((source==kGPi0NoFeedDown)||(source==kGEtaNoFeedDown))fInvMULSgamma->Fill(track->Pt(),mass);*/
		}
		
		  Double_t weight_pe=0.;//,weight_eta=0.;
		  if(  (source==kPi0NoFeedDown)  || (source==kGPi0NoFeedDown) )
		  {
		  		if(ptmotherwg >= 0.5 && ptmotherwg <= 2.5) weight_pe = f0a->Eval(ptmotherwg);
				if(ptmotherwg >2.5 && ptmotherwg < 15) weight_pe = f0b->Eval(ptmotherwg);
		  }					
		  if((source==kEtaNoFeedDown)  || (source==kGEtaNoFeedDown)  )
		  {      
				  if(ptmotherwg >= 0.5 && ptmotherwg < 2.3) weight_pe = f1->Eval(ptmotherwg);
				  if(ptmotherwg >= 2.3 && ptmotherwg < 4.1) weight_pe= f1a->Eval(ptmotherwg);
				  if(ptmotherwg >= 4.1 && ptmotherwg < 15) weight_pe = f1b->Eval(ptmotherwg);
		 }
		   if(fFlagULS){
        	fInvMULSnSpwg->Fill(invms,weight_pe);
        	}
	}
}

Double_t AliAnalysisTaskHFEmultTPCTOF::WeightMCncCorr(Int_t multSPDtr)
{
Double_t  MCnchweight=1.;
if(multSPDtr==0) MCnchweight= 1.81234;
else if(multSPDtr==1) MCnchweight= 1.13862;
else if(multSPDtr==2) MCnchweight= 0.979267;
else if(multSPDtr==3) MCnchweight= 0.901707;
else if(multSPDtr==4) MCnchweight= 0.878856;
else if(multSPDtr==5) MCnchweight= 0.889499;
else if(multSPDtr==6) MCnchweight= 0.918609;
else if(multSPDtr==7) MCnchweight= 0.953161;
else if(multSPDtr==8) MCnchweight= 0.981919;
else if(multSPDtr==9) MCnchweight= 1.008;
else if(multSPDtr==10) MCnchweight= 1.02534;
else if(multSPDtr==11) MCnchweight= 1.04052;
else if(multSPDtr==12) MCnchweight= 1.05301;
else if(multSPDtr==13) MCnchweight= 1.05925;
else if(multSPDtr==14) MCnchweight= 1.06661;
else if(multSPDtr==15) MCnchweight= 1.06843;
else if(multSPDtr==16) MCnchweight= 1.06544;
else if(multSPDtr==17) MCnchweight= 1.06622;
else if(multSPDtr==18) MCnchweight= 1.06464;
else if(multSPDtr==19) MCnchweight= 1.05136;
else if(multSPDtr==20) MCnchweight= 1.04541;
else if(multSPDtr==21) MCnchweight= 1.03675;
else if(multSPDtr==22) MCnchweight= 1.02735;
else if(multSPDtr==23) MCnchweight= 1.01762;
else if(multSPDtr==24) MCnchweight= 1.00628;
else if(multSPDtr==25) MCnchweight= 0.998685;
else if(multSPDtr==26) MCnchweight= 0.989993;
else if(multSPDtr==27) MCnchweight= 0.973123;
else if(multSPDtr==28) MCnchweight= 0.973184;
else if(multSPDtr==29) MCnchweight= 0.96248;
else if(multSPDtr==30) MCnchweight= 0.950655;
else if(multSPDtr==31) MCnchweight= 0.944755;
else if(multSPDtr==32) MCnchweight= 0.932536;
else if(multSPDtr==33) MCnchweight= 0.928509;
else if(multSPDtr==34) MCnchweight= 0.925921;
else if(multSPDtr==35) MCnchweight= 0.923119;
else if(multSPDtr==36) MCnchweight= 0.921073;
else if(multSPDtr==37) MCnchweight= 0.916071;
else if(multSPDtr==38) MCnchweight= 0.928057;
else if(multSPDtr==39) MCnchweight= 0.916297;
else if(multSPDtr==40) MCnchweight= 0.921829;
else if(multSPDtr==41) MCnchweight= 0.940199;

return MCnchweight;



}
//=================================================================================================================================
Int_t AliAnalysisTaskHFEmultTPCTOF::ClassifyTrack(AliAODTrack* track,const AliVVertex* pVtx)
{  
  
	Double_t pt = track->Pt();
	Double_t eta = track->Eta();
	Double_t phi = track->Phi();
	Float_t dx,dy,dxy, dz;

	if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return 0; //fitler bit 
  if(pt< 0.3 ) return 0;  
	if (TMath::Abs(eta)>=fEtarange) return 0; 
	
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
Int_t AliAnalysisTaskHFEmultTPCTOF::GetElecSourceType(AliAODMCParticle *electron, TClonesArray *mcArray,Double_t &ptm)
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
Int_t AliAnalysisTaskHFEmultTPCTOF::GetPi0EtaType(AliAODMCParticle *pi0eta, TClonesArray *mcArray)
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

Int_t AliAnalysisTaskHFEmultTPCTOF::GetHFE(AliAODMCParticle *electron, TClonesArray *mcArray)
{
	Int_t motherindex=electron->GetMother(); //Getting Electron Mother
	if(motherindex<0) return kNoMother;
	AliAODMCParticle *mother = (AliAODMCParticle*)mcArray->At(motherindex);					
	Int_t motherpdg = mother->GetPdgCode();		  
	if ( (int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty;
	if ( (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm;
	return kOthersE;
  
}

//====================================================================================================================================


Int_t AliAnalysisTaskHFEmultTPCTOF::GetNcharged(){
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

TProfile* AliAnalysisTaskHFEmultTPCTOF::GetEstimatorHistogram(const AliAODEvent* fAOD)
{
    
  Int_t runNo  = fAOD->GetRunNumber();
 //cout<<"run number"<<runNo<<endl;
  Int_t period = -1; 
   
        
  if (runNo>=258962 && runNo<=259888){
   period = 0;
   if(fIsMC) period = 1;
  } 
  if (runNo>=256504 && runNo<=258537){
   period = 2;
   if(fIsMC) period = 3;
  } 
  if (period < 0 || period > 3) return 0;
   
    
  return fMultEstimatorAvg[period];
}
Double_t AliAnalysisTaskHFEmultTPCTOF::GetCorrectedNtracklets(TProfile* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult) {
  //
  // Correct the number of accepted tracklets based on a TProfile Hist
  //
  //

  if(TMath::Abs(vtxZ)>10.0){
    //    printf("ERROR: Z vertex out of range for correction of multiplicity\n");
    return uncorrectedNacc;
  }

  if(!estimatorAvg){
    printf("ERROR: Missing TProfile for correction of multiplicity\n");
    return uncorrectedNacc;
  }

  Double_t localAvg = estimatorAvg->GetBinContent(estimatorAvg->FindBin(vtxZ));

  Double_t deltaM = uncorrectedNacc*(refMult/localAvg - 1);

  Double_t correctedNacc = uncorrectedNacc + (deltaM>0 ? 1 : -1) * gRandom->Poisson(TMath::Abs(deltaM));

  return correctedNacc;
}

//===================================================================================
Double_t AliAnalysisTaskHFEmultTPCTOF::Beta(AliAODTrack *track)
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
void AliAnalysisTaskHFEmultTPCTOF::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  
	fOutputList = dynamic_cast<TList*> (GetOutputData(1));
	if (!fOutputList)return;

}

