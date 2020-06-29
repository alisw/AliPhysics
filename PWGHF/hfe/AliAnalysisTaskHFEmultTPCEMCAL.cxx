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
#include "AliAnalysisTaskHFEmultTPCEMCAL.h"
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

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
//#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskHFEmultTPCEMCAL)

//________________________________________________________________________//
AliAnalysisTaskHFEmultTPCEMCAL::AliAnalysisTaskHFEmultTPCEMCAL(): 
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
	fEtarange(0.7),
	fTPCnsigmin(-1),
	fTPCnsigmax(3),
	IsM20(kFALSE),
	fCutM20Min(0.02),
	fCutM20Max1(0.9),
	fCutM20Max2(0.7),
	fCutM20Max3(0.5),
	fCutEopEMin(0.8),
	fCutEopEMax(1.2),


	fInvmassCut(0.14),	
	fAssoTPCCluster(60),
	fAssoITSRefit(kTRUE),
	fAssopTMin(0.1),
	fAssoEtarange(0.9),
	fAssoTPCnsig(3.5),
	fdeltaeta(0.01),
   fdeltaphi(0.01),

	//-------------Analysis
	fHistPt(0), fHistMult(0), fTPCSignal(0), feta(0), fNentries(0), fNentries2(0),
	fVtxZ(0),fVtxZ_corr(0),
	fClusEtaPhi(0),
   fClusT(0),
   fNCells(0), 
   fClusE(0),
   fClusEvsnTracklets(0),
	fHistPtMatch(0),
	fEMCTrkMatch(0),
	fEMCTrkPt(0),
	fEMCTrketa(0),
	fEMCTrkphi(0),
	fEMCTPCnsig(0),
	fClsEAftMatch(0),
	fClsEAftMatch_SPD(0),
	fClsEopAftMatch(0),
	fClsEtaPhiAftMatch(0),
	fSparseElectron(0),
	fvalueElectron(0),
	// flag for emcal dcal
    fFlagClsTypeEMC(kTRUE),
    fFlagClsTypeDCAL(kTRUE),
    // trigger events selection
  	fEMCEG1(kFALSE),
  	fEMCEG2(kFALSE),
  	fDCalDG1(kFALSE),
  	fDCalDG2(kFALSE),
  	 	// emcal correction
  	fUseTender(kTRUE),
  	fTenderClusterName("caloClusters"),
  	fTenderTrackName("tracks"),
  	fTracks_tender(0),
  	fCaloClusters_tender(0),
  		//--SPD
  	  fSPD_tracklet(0),
  fSPDCorrMultDist_max(0),
  fSPDWeightedCorrMultDist_max(0),
  Profile_Mean(0),
  Profile_MeanCorr(0),
   fPeriod(0),
   SPDntr(0),
  //------------------------Non-Hfe
fNonHFE(new AliSelectNonHFE()),
fInvmassLS1(0),fInvmassULS1(0),
fPte_ULS(0),fPte_LS(0),
fPte_ULS_multSPD(0),fPte_LS_multSPD(0),

//------------------------------MC 
fIsMC(kFALSE),
fMCArray(0),fMCHeader(0),
fMCparticle(0),fMCmother(0),
pdg(0),
fPtHFEMC(0),
fPtHFEMC_SPD(0),
fPtHFEMC_reco(0),
fPtHFEMC_reco_SPD(0),
  fPtHFEMC_trackcutreco(0),
  fPtHFEMC_trackcutreco_SPD(0),
  fPtHFEMC_trackmatchreco(0),
  fPtHFEMC_trackmatchreco_SPD(0),
  fPtHFEMC_TPCEMCreco(0),
  fPtHFEMC_TPCEMCreco_SPD(0),
   fPtHFEMC_SScutreco(0),
 fPtHFEMC_SScutreco_SPD(0),
    fPtHFEMC_TPCreco(0),
  fPtHFEMC_TPCreco_SPD(0),
fPT_elec_MC(0),
fPT_PIDCaloCut_MC(0),fPT_PIDCaloCut_realistic(0),fPT_elec_realistic(0),
fPt_elec_phot1(0),
fPt_elec_phot1_multSPD(0),
fInvMULS(0),
fInvMULSnSp(0),
fPi0EtaSpectraSp(0),
pi0MC(0),
etaMC(0),
gammaMC(0),
fNTotMCpart(0),fNpureMC(0),  fNembMCpi0(0), fNembMCeta(0),
fCalculateWeight(kFALSE), fSprsPi0EtaWeightCal(0),
fIsFrmEmbPi0(kFALSE),
fIsFrmEmbEta(kFALSE),
fIsFrmPi0(kFALSE),
fIsFrmEta(kFALSE),
ftype(-1),
fWeight(1),
fWeightPi0(1),
fWeightEta(1),
fPi0Weight(0),
fEtaWeight(0),
fRealInclsElecPt(0),
fNonHFeTrkPt(0),
fNonHFeEmbTrkPt(0),
fNonHFeEmbWeightTrkPt(0),
fNonHFeEmbTrkPt_SPD(0),
fNonHFeEmbWeightTrkPt_SPD(0),
fPi0eEmbWeightTrkPt(0),
fEtaeEmbWeightTrkPt(0),
fRecoNonHFeTrkPt(0),
fRecoNonHFeEmbTrkPt(0),
fRecoNonHFeEmbWeightTrkPt(0),
fRecoNonHFeTrkPt_SPD(0),
fRecoNonHFeEmbTrkPt_SPD(0),
fRecoNonHFeEmbWeightTrkPt_SPD(0),
fRecoPi0eEmbWeightTrkPt(0),
fRecoEtaeEmbWeightTrkPt(0),
fInvmassULSPt(0),
fInvmassLSPt(0),
fULSElecPt(0),
fLSElecPt(0),

fNonHFePairInvmassLS(0),
fNonHFePairInvmassULS(0),
fNonHFeEmbInvmassLS(0),
fNonHFeEmbInvmassULS(0),
fNonHFeEmbWeightInvmassLS(0),
fNonHFeEmbWeightInvmassULS(0),
fPi0EmbInvmassLS(0),
fPi0EmbInvmassULS(0),
fPi0EmbWeightInvmassLS(0),
fPi0EmbWeightInvmassULS(0),
fEtaEmbInvmassLS(0),
fEtaEmbInvmassULS(0),
fEtaEmbWeightInvmassLS(0),
fEtaEmbWeightInvmassULS(0),
fRecoLSeEmbTrkPt(0),
fRecoLSeEmbWeightTrkPt(0),
fRecoPi0LSeEmbWeightTrkPt(0),
fRecoEtaLSeEmbWeightTrkPt(0),
fRecoULSeEmbTrkPt(0),
fRecoULSeEmbWeightTrkPt(0),
fRecoPi0ULSeEmbWeightTrkPt(0),
fRecoEtaULSeEmbWeightTrkPt(0)
	{
	  fvalueElectron = new Double_t[9];
	  fPID = new AliHFEpid("hfePid");
	  for(Int_t i=0; i<5; i++) fMultEstimatorAvg[i]=0;
	}


AliAnalysisTaskHFEmultTPCEMCAL::AliAnalysisTaskHFEmultTPCEMCAL(const char *name) 
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
	fEtarange(0.7),
	fTPCnsigmin(-1),
	fTPCnsigmax(3),
	IsM20(kFALSE),
	fCutM20Min(0.02),
	fCutM20Max1(0.9),
	fCutM20Max2(0.7),
	fCutM20Max3(0.5),
	fCutEopEMin(0.8),
	fCutEopEMax(1.2),


	fInvmassCut(0.14),
	fAssoTPCCluster(60),
	fAssoITSRefit(kTRUE),
	fAssopTMin(0.1),
	fAssoEtarange(0.9),
	fAssoTPCnsig(3.5),
	fdeltaeta(0.01),
   fdeltaphi(0.01),
   
	//-------------Analysis
	fHistPt(0), fHistMult(0), fTPCSignal(0), feta(0),fNentries(0), fNentries2(0),
	fVtxZ(0),fVtxZ_corr(0),
	   fClusEtaPhi(0),
   fClusT(0),
   fNCells(0), 
   fClusE(0),
   fClusEvsnTracklets(0),
	fHistPtMatch(0),
	fEMCTrkMatch(0),
	fEMCTrkPt(0),
	fEMCTrketa(0),
	fEMCTrkphi(0),
	fEMCTPCnsig(0),
	fClsEAftMatch(0),
	fClsEAftMatch_SPD(0),
	fClsEopAftMatch(0),
	fClsEtaPhiAftMatch(0),
	fSparseElectron(0),
	fvalueElectron(0),
    // flag for emcal dcal
    fFlagClsTypeEMC(kTRUE),
	fFlagClsTypeDCAL(kTRUE),
	// trigger events selection
  	fEMCEG1(kFALSE),
  	fEMCEG2(kFALSE),
  	fDCalDG1(kFALSE),
  	fDCalDG2(kFALSE),
  	// emcal correction
  	fUseTender(kTRUE),
  	fTenderClusterName("caloClusters"),
  	fTenderTrackName("tracks"),
  	fTracks_tender(0),
  	fCaloClusters_tender(0),
  	//--SPD
  	  fSPD_tracklet(0),
  fSPDCorrMultDist_max(0),
  fSPDWeightedCorrMultDist_max(0),
  Profile_Mean(0),
  Profile_MeanCorr(0),
   fPeriod(0),
  SPDntr(0),
  
  //------------------------Non-Hfe
fNonHFE(new AliSelectNonHFE()),
fInvmassLS1(0),fInvmassULS1(0),
fPte_ULS(0),fPte_LS(0),
fPte_ULS_multSPD(0),fPte_LS_multSPD(0),

//------------------------------MC 
fIsMC(kFALSE),
fMCArray(0),fMCHeader(0),
fMCparticle(0),fMCmother(0),
pdg(0),
fPtHFEMC(0),
fPtHFEMC_SPD(0),
fPtHFEMC_reco(0),
fPtHFEMC_reco_SPD(0),
  fPtHFEMC_trackcutreco(0),
  fPtHFEMC_trackcutreco_SPD(0),
  fPtHFEMC_trackmatchreco(0),
  fPtHFEMC_trackmatchreco_SPD(0),
  fPtHFEMC_TPCEMCreco(0),
  fPtHFEMC_TPCEMCreco_SPD(0),
  fPtHFEMC_SScutreco(0),
 fPtHFEMC_SScutreco_SPD(0),
    fPtHFEMC_TPCreco(0),
  fPtHFEMC_TPCreco_SPD(0),
fPT_elec_MC(0),
fPT_PIDCaloCut_MC(0),fPT_PIDCaloCut_realistic(0),fPT_elec_realistic(0),
fPt_elec_phot1(0),
fPt_elec_phot1_multSPD(0),
fInvMULS(0),
fInvMULSnSp(0),
fPi0EtaSpectraSp(0),
pi0MC(0),
etaMC(0),
gammaMC(0),
fNTotMCpart(0),fNpureMC(0),  fNembMCpi0(0), fNembMCeta(0),
fCalculateWeight(kFALSE),  fSprsPi0EtaWeightCal(0),
fIsFrmEmbPi0(kFALSE),
fIsFrmEmbEta(kFALSE),
fIsFrmPi0(kFALSE),
fIsFrmEta(kFALSE),
ftype(-1),
fWeight(1),
fWeightPi0(1),
fWeightEta(1),
fPi0Weight(0),
fEtaWeight(0),
fRealInclsElecPt(0),
fNonHFeTrkPt(0),
fNonHFeEmbTrkPt(0),
fNonHFeEmbWeightTrkPt(0),
fNonHFeEmbTrkPt_SPD(0),
fNonHFeEmbWeightTrkPt_SPD(0),
fPi0eEmbWeightTrkPt(0),
fEtaeEmbWeightTrkPt(0),
fRecoNonHFeTrkPt(0),
fRecoNonHFeEmbTrkPt(0),
fRecoNonHFeEmbWeightTrkPt(0),
fRecoNonHFeTrkPt_SPD(0),
fRecoNonHFeEmbTrkPt_SPD(0),
fRecoNonHFeEmbWeightTrkPt_SPD(0),
fRecoPi0eEmbWeightTrkPt(0),
fRecoEtaeEmbWeightTrkPt(0),
fInvmassULSPt(0),
fInvmassLSPt(0),
fULSElecPt(0),
fLSElecPt(0),


fNonHFePairInvmassLS(0),
fNonHFePairInvmassULS(0),
fNonHFeEmbInvmassLS(0),
fNonHFeEmbInvmassULS(0),
fNonHFeEmbWeightInvmassLS(0),
fNonHFeEmbWeightInvmassULS(0),
fPi0EmbInvmassLS(0),
fPi0EmbInvmassULS(0),
fPi0EmbWeightInvmassLS(0),
fPi0EmbWeightInvmassULS(0),
fEtaEmbInvmassLS(0),
fEtaEmbInvmassULS(0),
fEtaEmbWeightInvmassLS(0),
fEtaEmbWeightInvmassULS(0),
fRecoLSeEmbTrkPt(0),
fRecoLSeEmbWeightTrkPt(0),
fRecoPi0LSeEmbWeightTrkPt(0),
fRecoEtaLSeEmbWeightTrkPt(0),
fRecoULSeEmbTrkPt(0),
fRecoULSeEmbWeightTrkPt(0),
fRecoPi0ULSeEmbWeightTrkPt(0),
fRecoEtaULSeEmbWeightTrkPt(0)

{
  // Constructor
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  fvalueElectron = new Double_t[9];
  fPID = new AliHFEpid("hfePid");
  DefineInput(0, TChain::Class());  
  DefineOutput(1, TList::Class());
  DefineOutput(2, TH1F::Class());
  for(Int_t i=0; i<5; i++) fMultEstimatorAvg[i]=0;
}

//_________________________Destructer_____________________________________
AliAnalysisTaskHFEmultTPCEMCAL::~AliAnalysisTaskHFEmultTPCEMCAL()
{
  delete fPID;
  delete fSparseElectron;
  delete []fvalueElectron;	
    delete fSprsPi0EtaWeightCal;
  if (fOutputList) { delete fOutputList; fOutputList = 0;}
  if (fNentries){ delete fNentries; fNentries = 0;}
  if (fNentries2){ delete fNentries2; fNentries2 = 0;}
  
     for(Int_t i=0; i<5; i++) {
      if (fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i];
  }
}

//_______________________________________________________________________
void AliAnalysisTaskHFEmultTPCEMCAL::Init()
{
  // Initialization	
  if(fDebug > 1) printf("AliAnalysisTaskHFEmultTPCEMCAL::Init() \n");
  return;
}




//________________________________________________________________________
void AliAnalysisTaskHFEmultTPCEMCAL::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
	//Initialize PID

	 
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
	
  Double_t pi=TMath::Pi();

	
    fPi0Weight = new TF1("fPi0Weight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    fEtaWeight = new TF1("fEtaWeight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    fPi0Weight->SetParameters(2.91546e+03,-8.25772e-02,5.50975e-04,9.57099e-01,4.67815e+00);
	//2.60217e+03,-4.73807e-02, -1.73834e-04, 9.94091e-01,4.77806e+00);
    fEtaWeight->SetParameters(4.20560e+02,-4.71568e-02,5.38450e-03,2.06379e+00,5.88947e+00);
    //3.43821e+02,-6.29809e-02,6.65830e-03,2.83334e+00,6.73711e+00);
  
  
  fcount = new TH1D("fcount", "fcount", 10, 0.0, 10.0);
  fcount->Sumw2();  
  
  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 150, 0.1, 15.1);
  //fOutputList->Add(fHistPt);
  fHistMult = new TH1F("fHistMult", "Multiplicity distribution", 500, 0, 500);
  //fOutputList->Add(fHistMult);
  fTPCSignal = new TH1F("fTPCSignal", "fTPCSignal distribution", 100, 0, 100);
  //fOutputList->Add(fTPCSignal);
	
	 feta = new TH1F("feta","#eta distribution of tracks ;#eta;counts",100,-1.5,1.5);
  fOutputList->Add(feta);
	//--------------------------Vertex Dist------------------------------------
	fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
  fOutputList->Add(fVtxZ);

	fVtxZ_corr=new TH1F("fVtxZ_corr","Z vertex position;Vtx_{z};counts",1000,-50,50);
  fOutputList->Add(fVtxZ_corr);
  
 	
	//--------------------------TPC dedx nsig vs p ----------------------------------------------
  fdEdxVsP_TPC = new TH2F("fdEdxVsP_TPC", "fdEdxVsP_TPC distribution",300,0,15,750,10,160);
  fdEdxVsP_TPC->Sumw2();
  //fOutputList->Add(fdEdxVsP_TPC);
  
  fnSigmaVsP_TPC= new TH2F("fnSigmaVsP_TPC", "fnSigmaVsP_TPC distribution",300,0.,15.,750,-15.,15.);
  fnSigmaVsP_TPC->Sumw2(); 
  //fOutputList->Add(fnSigmaVsP_TPC);
  
  
 //--------------------------TPC dedx nsig vs pt ----------------------------------------------
   
  fnSigmaVsPt_TPC= new TH2F("fnSigmaVsPt_TPC", "fnSigmaVsPt_TPC distribution",300,0.,15.,750,-15.,15.);
  fnSigmaVsPt_TPC->Sumw2(); 
  //fOutputList->Add(fnSigmaVsPt_TPC);  
  
  //=================================EMCalCluster Properties=============================
   fClusEtaPhi		= new TH2F( "fClusEtaPhi","Cluster Eta Phi distribution; #eta ; #phi",50,-2,2,100,0.,7); 
    fOutputList->Add(fClusEtaPhi);
   fClusT     		= new TH1F( "fClusT","Cluster time distribution ; Time(ns) ; counts",500,-1000,1000);
    fOutputList->Add(fClusT);
   fNCells   		= new TH1F("fNCells","ncells distribution ; cell counts ; cluster counts", 50,-10,40);
    fOutputList->Add(fNCells);
   fClusE   		= new TH1F("fClusE","Cluster Energy ; Energy(GeV); counts",200,0.,100.);
    fOutputList->Add(fClusE);
   fClusEvsnTracklets= new TH2F("fClusEvsnTracklets","Cluster Energy vs Multiplicity ; Energy(GeV); SPD Tracklets",200,0,100,300,-0.5,299.5);
    fOutputList->Add(fClusEvsnTracklets);

 //===================================Track Matched to Cluster================================================================//
 
  fHistPtMatch = new TH1F("fHistPtMatch", "p_{T} distribution of tracks matched to EMCAL;p_{T} (GeV/c);counts",1000, 0.0, 100.0);
  fOutputList->Add(fHistPtMatch);

  fEMCTrkMatch = new TH2F("fEMCTrkMatch","Distance of EMCAL cluster to its closest track;#phi;z",100,-0.3,0.3,100,-0.3,0.3);
  fOutputList->Add(fEMCTrkMatch);

  fEMCTrkPt = new TH1F("fEMCTrkPt","p_{T} distribution of tracks with EMCAL cluster;p_{T} (GeV/c);counts",1000,0,100.);
  fOutputList->Add(fEMCTrkPt);

  fEMCTrketa = new TH1F("fEMCTrketa","#eta distribution of tracks matched to EMCAL;#eta;counts",100,-1.5,1.5);
  fOutputList->Add(fEMCTrketa);

  fEMCTrkphi = new TH1F("fEMCTrkphi","#phi distribution of tracks matched to EMCAL;#phi;counts",100,0,2*pi);
  fOutputList->Add(fEMCTrkphi);


  fEMCTPCnsig = new TH2F("fEMCTPCnsig","TPC Nsigma distribution of tracks matched to EMCAL;p (GeV/c);#sigma_{TPC-dE/dx}",250,0,50,200,-10,10);
  fOutputList->Add(fEMCTPCnsig);

  fClsEAftMatch = new TH1F("fClsEAftMatch", "EMCAL cluster energy distribution after track matching; Cluster E;counts", 100, 0.0, 50.0);
  fOutputList->Add(fClsEAftMatch);
  
   fClsEAftMatch_SPD = new TH2F("fClsEAftMatch_SPD", "EMCAL cluster energy distribution after track matching; Cluster E;counts", 100, 0.0, 50.0,300,-0.5,299.5);
  fOutputList->Add(fClsEAftMatch_SPD);
  
  fClsEopAftMatch = new TH1F("fClsEopAftMatch", "EMCAL cluster energy over p distribution after track matching; Cluster E/p;counts", 100, 0.0, 2.0);
  fOutputList->Add(fClsEopAftMatch);

  fClsEtaPhiAftMatch = new TH2F("fClsEtaPhiAftMatch","EMCAL cluster #eta and #phi distribution after track matching;#eta;#phi",100,-0.9,0.9,200,0,6.3);
  fOutputList->Add(fClsEtaPhiAftMatch);

  
  /*Int_t bins[9]		=   {380,   180, 100, 100, 100, 480, 1000,   200,  1000};
  Double_t xmin[9]	=	{  2,  -10,   0,   0,   0,   2,   0,      0 ,   -50};
  Double_t xmax[9]	=	{  40,   8,   2,   2,  2,   50,  1000,   100,    50};

  fSparseElectron 	= new THnSparseD ("Electron","Electron;pT;nSigma;E/P;m02;m20;p;SPDTracklets;ClusterEnergy;zvtx;",9 ,bins,xmin,xmax);
   fSparseElectron->Sumw2();
   fOutputList->Add(fSparseElectron);
		 

	
	Int_t binselec[6]		= {380,   180,  100,  100,  100,   300};
  Double_t xminelec[6]	=	{  2,   -10,    0,    0,    0,     0};
  Double_t xmaxelec[6]	=	{ 40,     8,    2,    2,    2,   300};

  fSparseElectron 	= new THnSparseD ("Electron","Electron;pT;nSigma;E/P;m02;m20;SPDTracklets;",6,binselec,xminelec,xmaxelec);
  fSparseElectron->GetAxis(0)->SetName("pT");
  fSparseElectron->GetAxis(1)->SetName("nSigma");
  fSparseElectron->GetAxis(2)->SetName("E/P");
  fSparseElectron->GetAxis(3)->SetName("m02");
  fSparseElectron->GetAxis(4)->SetName("m20");
  fSparseElectron->GetAxis(5)->SetName("SPDTracklets");
  fSparseElectron->Sumw2();
  fOutputList->Add(fSparseElectron);
 */ 

  Int_t binselec[5]		=      {380,   180,  200,  100,    300};
  Double_t xminelec[5]	=	{  2,   -10,    0,    0,       0};
  Double_t xmaxelec[5]	=	{ 40,     8,    2,    2,     300};

  fSparseElectron 	= new THnSparseD ("Electron","Electron;pT;nSigma;E/P;m02;SPDTracklets;",5,binselec,xminelec,xmaxelec);
  fSparseElectron->GetAxis(0)->SetName("pT");
  fSparseElectron->GetAxis(1)->SetName("nSigma");
  fSparseElectron->GetAxis(2)->SetName("E/P");
  fSparseElectron->GetAxis(3)->SetName("m02");
  fSparseElectron->GetAxis(4)->SetName("SPDTracklets");
  fSparseElectron->Sumw2();
  fOutputList->Add(fSparseElectron);
/*
	Int_t binsls[3]	=      	{280, 1000, 1000};
  Double_t xminls[3]	=	{  2, 0, 0};
  Double_t xmaxls[3]	=	{  30, 2000, 1000};

  fSparseLSElectron 	= new THnSparseD ("LSElectron","LSElectron;pT;V0M;SPDTracklets;",3 ,binsls,xminls,xmaxls);
  fSparseULSElectron 	= new THnSparseD ("ULSElectron","ULSElectron;pT;V0M;SPDTracklets;",3 ,binsls,xminls,xmaxls);
*/
 //-----------------NonHFE-----------------------------------------------------------------------------------------
  
  fInvmassLS1 = new TH1F("fInvmassLS1","Inv mass of LS (e,e) for pt^{e}; mass(GeV/c^2); counts;",1000,0,1.0);
	fOutputList->Add(fInvmassLS1);
	
	fInvmassULS1 = new TH1F("fInvmassULS1","Inv mass of ULS (e,e) for pt^{e}; mass(GeV/c^2); counts;",1000,0,1.0);
	fOutputList->Add(fInvmassULS1);
	
	//....................................................................................
	
   fPte_ULS = new TH1F("fPte_ULS", "ULS electron pt",380,2.,40.);
  fPte_ULS->Sumw2();
  fOutputList->Add(fPte_ULS);
  
  fPte_LS = new TH1F("fPte_LS", "LS electron pt",380,2.,40.);
  fPte_LS->Sumw2();
  fOutputList->Add(fPte_LS);
 
 
  fPte_ULS_multSPD = new TH2F("fPte_ULS_multSPD", "ULS electron pt in percentile",380,2.,40.,300,-0.5,299.5);
  fPte_ULS_multSPD->Sumw2();
  fOutputList->Add(fPte_ULS_multSPD);
  
  fPte_LS_multSPD = new TH2F("fPte_LS_multSPD", "LS electron pt in percentile",380,2.,40.,300,-0.5,299.5);
  fPte_LS_multSPD->Sumw2();
  fOutputList->Add(fPte_LS_multSPD);
  
 //+++++++++++++++++++++++++++++++++++++++++MC++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  
 
  fPtHFEMC = new TH1F("fPtHFEMC",";p_{t} (GeV/c)",380,2.,40.);
  fPtHFEMC->Sumw2();
  if(fIsMC) fOutputList->Add(fPtHFEMC);
  
  fPtHFEMC_SPD = new TH2F("fPtHFEMC_SPD",";p_{t} (GeV/c)",380,2.,40.,300,-0.5,299.5);
  fPtHFEMC_SPD->Sumw2();
  if(fIsMC) fOutputList->Add(fPtHFEMC_SPD);
  
   fPtHFEMC_reco = new TH1F("fPtHFEMC_reco",";p_{t} (GeV/c)",380,2.,40.);
  fPtHFEMC_reco->Sumw2();
  if(fIsMC)fOutputList->Add(fPtHFEMC_reco);
  
  fPtHFEMC_reco_SPD = new TH2F("fPtHFEMC_reco_SPD",";p_{t} (GeV/c)",380,2.,40.,300,-0.5,299.5);
  fPtHFEMC_reco_SPD->Sumw2();
  if(fIsMC)fOutputList->Add(fPtHFEMC_reco_SPD);
  
   fPtHFEMC_trackcutreco = new TH1F("fPtHFEMC_trackcutreco",";p_{t} (GeV/c)",380,2.,40.);
  fPtHFEMC_trackcutreco->Sumw2();
  if(fIsMC)fOutputList->Add(fPtHFEMC_trackcutreco);
  
  fPtHFEMC_trackcutreco_SPD = new TH2F("fPtHFEMC_trackcutreco_SPD",";p_{t} (GeV/c)",380,2.,40.,300,-0.5,299.5);
  fPtHFEMC_trackcutreco_SPD->Sumw2();
  if(fIsMC)fOutputList->Add(fPtHFEMC_trackcutreco_SPD);
  
   fPtHFEMC_trackmatchreco = new TH1F("fPtHFEMC_trackmatchreco",";p_{t} (GeV/c)",380,2.,40.);
  fPtHFEMC_trackmatchreco->Sumw2();
  if(fIsMC)fOutputList->Add(fPtHFEMC_trackmatchreco);
  
  fPtHFEMC_trackmatchreco_SPD = new TH2F("fPtHFEMC_trackmatchreco_SPD",";p_{t} (GeV/c)",380,2.,40.,300,-0.5,299.5);
  fPtHFEMC_trackmatchreco_SPD->Sumw2();
  if(fIsMC)fOutputList->Add(fPtHFEMC_trackmatchreco_SPD);
  
	
  fPtHFEMC_TPCEMCreco = new TH1F("fPtHFEMC_TPCEMCreco",";p_{t} (GeV/c)",380,2.,40.);
  fPtHFEMC_TPCEMCreco->Sumw2();
  if(fIsMC)fOutputList->Add(fPtHFEMC_TPCEMCreco);

		
  fPtHFEMC_TPCEMCreco_SPD = new TH2F("fPtHFEMC_TPCEMCreco_SPD",";p_{t} (GeV/c)",380,2.,40.,300,-0.5,299.5);
  fPtHFEMC_TPCEMCreco_SPD->Sumw2();
  if(fIsMC)fOutputList->Add(fPtHFEMC_TPCEMCreco_SPD);
  
   fPtHFEMC_SScutreco = new TH1F("fPtHFEMC_SScutreco",";p_{t} (GeV/c)",380,2.,40.);
  fPtHFEMC_SScutreco->Sumw2();
  if(fIsMC)fOutputList->Add(fPtHFEMC_SScutreco);

		
  fPtHFEMC_SScutreco_SPD = new TH2F("fPtHFEMC_SScutreco_SPD",";p_{t} (GeV/c)",380,2.,40.,300,-0.5,299.5);
  fPtHFEMC_SScutreco_SPD->Sumw2();
  if(fIsMC)fOutputList->Add(fPtHFEMC_SScutreco_SPD);
  
   fPtHFEMC_TPCreco = new TH1F("fPtHFEMC_TPCreco",";p_{t} (GeV/c)",380,2.,40.);
  fPtHFEMC_TPCreco->Sumw2();
  if(fIsMC)fOutputList->Add(fPtHFEMC_TPCreco);


		
  fPtHFEMC_TPCreco_SPD = new TH2F("fPtHFEMC_TPCreco_SPD",";p_{t} (GeV/c)",380,2.,40.,300,-0.5,299.5);
  fPtHFEMC_TPCreco_SPD->Sumw2();
  if(fIsMC)fOutputList->Add(fPtHFEMC_TPCreco_SPD);

  fPT_elec_MC = new TH1F("fPT_elec_MC","Pt from MC",380,2.,40.);
  fPT_elec_MC->Sumw2();
  if(fIsMC)fOutputList->Add(fPT_elec_MC);
  
   fPT_elec_realistic= new TH1F("fPT_elec_realistic","Pt from MC",380,2.,40.);
  fPT_elec_realistic->Sumw2();
  if(fIsMC)fOutputList->Add(fPT_elec_realistic);
  
  fPT_PIDCaloCut_MC= new TH1F("fPT_PIDCaloCut_MC","Pt from MC",380,2.,40.);
  fPT_PIDCaloCut_MC->Sumw2();
  if(fIsMC)fOutputList->Add(fPT_PIDCaloCut_MC);
  
  fPT_PIDCaloCut_realistic= new TH1F("fPT_PIDCaloCut_realistic","Pt from MC",380,2.,40.);
  fPT_PIDCaloCut_realistic->Sumw2();
  if(fIsMC)fOutputList->Add(fPT_PIDCaloCut_realistic);
  
  fPt_elec_phot1= new TH1F("fPt_elec_phot1", "; p_{T}(GeV/c); counts;",380,2.,40.);
  fPt_elec_phot1->Sumw2();
  if(fIsMC) fOutputList->Add(fPt_elec_phot1);
  
  fPt_elec_phot1_multSPD= new TH2F("fPt_elec_phot1_multSPD", "; p_{T}(GeV/c); counts;",380,2.,40.,300,-0.5,299.5);
  fPt_elec_phot1_multSPD->Sumw2();
  if(fIsMC)fOutputList->Add(fPt_elec_phot1_multSPD);

	 fInvMULS = new TH2F("fInvMULS","Pt and InvMass from MC",380,2.,40.,1000,0,1.0);
  fInvMULS->Sumw2();
  if(fIsMC)fOutputList->Add(fInvMULS);
 

	Int_t nbinsInvMULS[4] = {380, 6, 200, 300};
  Double_t binlowInvMULS[4] = {2., 0, 0., -0.5};
  Double_t binupInvMULS[4] = {40, 6., 2., 299.5};

	fInvMULSnSp = new THnSparseF("fInvMULSnSp", "fInvMULSnSp;pt;source;mass;multdepSPD;", 4, nbinsInvMULS, binlowInvMULS, binupInvMULS);
	fInvMULSnSp->Sumw2();
	if(fIsMC)fOutputList->Add(fInvMULSnSp);

	

  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	Int_t nbinspt[3] = {400, 3, 7};
  Double_t binlow[3] = {0., 0, -1.};
  Double_t binup[3] = {40, 3, 6.};
	fPi0EtaSpectraSp = new THnSparseF("fPi0EtaSpectraSp", "fPi0EtaSpectraSp;pt;source;type;", 3, nbinspt, binlow,binup);
	fPi0EtaSpectraSp->Sumw2();
   if(fIsMC)fOutputList->Add(fPi0EtaSpectraSp);

	pi0MC= new TH1F("pi0MC",";p_{t} (GeV/c)",400,0.,40.);
  pi0MC->Sumw2();
  if(fIsMC)fOutputList->Add(pi0MC);
  
  etaMC= new TH1F("etaMC",";p_{t} (GeV/c)",400,0,40.);
  etaMC->Sumw2();
  if(fIsMC)fOutputList->Add(etaMC);
  
  gammaMC= new TH1F("gammaMC",";p_{t} (GeV/c)",400,0,40.);
  gammaMC->Sumw2();
  if(fIsMC)fOutputList->Add(gammaMC);

   if(fIsMC){
  //if(fCalculateWeight){
      /*  Int_t bin[4] = {500,3,2,7}; //pT, PDG, EnhancedSigOrNot, pi0etaType
        Double_t xmin[4] = {0,0,0,-1};
        Double_t xmax[4] = {50,3,2,6};
    
        fSprsPi0EtaWeightCal = new THnSparseD("fSprsPi0EtaWeightCal","Sparse to calculate #pi^{0} and #eta weight;p_{T};PDG ID;EnhanceSigOrNot;pi0etaType;",4,bin,xmin,xmax);
        fOutputList->Add(fSprsPi0EtaWeightCal);
       */
       Int_t bin[5] = {500,3,2,7,300}; //pT, PDG, EnhancedSigOrNot, pi0etaType, SPDntr
        Double_t xmin[5] = {0,0,0,-1,-0.5};
        Double_t xmax[5] = {50,3,2,6,299.5};
    
        fSprsPi0EtaWeightCal = new THnSparseD("fSprsPi0EtaWeightCal","Sparse to calculate #pi^{0} and #eta weight;p_{T};PDG ID;EnhanceSigOrNot;pi0etaType;SPDntrCorr;",5,bin,xmin,xmax);
        fOutputList->Add(fSprsPi0EtaWeightCal);
    //}
    
    
    
       fRealInclsElecPt = new TH1F("fRealInclsElecPt","p_{T} distribution of MC tagged inclusive electrons;p_{T} (GeV/c);counts",250,0,50);
        fOutputList->Add(fRealInclsElecPt);
     
       fNonHFeTrkPt = new TH1F("fNonHFeTrkPt","Non-HF electrons from all generators;p_{T} (GeV/c);counts",250,0,50);
        fNonHFeTrkPt->Sumw2();
        fOutputList->Add(fNonHFeTrkPt);
        
			fNonHFeEmbTrkPt = new TH1F("fNonHFeEmbTrkPt","Non-HF electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,50);
			fNonHFeEmbTrkPt->Sumw2();
         fOutputList->Add(fNonHFeEmbTrkPt);
            
          fNonHFeEmbWeightTrkPt = new TH1F("fNonHFeEmbWeightTrkPt","Non-HF electrons from embedded #pi^{0} and #eta + No mom with weight + No mom;p_{T} (GeV/c);counts",250,0,50);
          fNonHFeEmbWeightTrkPt->Sumw2();
          fOutputList->Add(fNonHFeEmbWeightTrkPt);
           
           fNonHFeEmbTrkPt_SPD= new TH2F("fNonHFeEmbTrkPt_SPD","Non-HF electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,50,300,-0.5,299.5);
			fNonHFeEmbTrkPt_SPD->Sumw2();
         fOutputList->Add(fNonHFeEmbTrkPt_SPD);
         
			fNonHFeEmbWeightTrkPt_SPD= new TH2F("fNonHFeEmbWeightTrkPt_SPD","Non-HF electrons from embedded #pi^{0} and #eta + No mom with weight + No mom;p_{T} (GeV/c);counts",250,0,50,300,-0.5,299.5);
          fNonHFeEmbWeightTrkPt_SPD->Sumw2();
          fOutputList->Add(fNonHFeEmbWeightTrkPt_SPD);
          
			fPi0eEmbWeightTrkPt = new TH1F("fPi0eEmbWeightTrkPt","Non-HF electrons from embedded #pi^{0} + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
         fPi0eEmbWeightTrkPt->Sumw2();
         fOutputList->Add(fPi0eEmbWeightTrkPt);
        
         fEtaeEmbWeightTrkPt = new TH1F("fEtaeEmbWeightTrkPt","Non-HF electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
         fEtaeEmbWeightTrkPt->Sumw2();
         fOutputList->Add(fEtaeEmbWeightTrkPt);
         
          fRecoNonHFeTrkPt = new TH1F("fRecoNonHFeTrkPt"," Reco Non-HF electrons from all generators;p_{T} (GeV/c);counts",250,0,50);
        fRecoNonHFeTrkPt->Sumw2();
        fOutputList->Add(fRecoNonHFeTrkPt);
        
          fRecoNonHFeEmbTrkPt = new TH1F("fRecoNonHFeEmbTrkPt","Reco Non-HF electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,50);
            fRecoNonHFeEmbTrkPt->Sumw2();
            fOutputList->Add(fRecoNonHFeEmbTrkPt);
        
            fRecoNonHFeEmbWeightTrkPt = new TH1F("fRecoNonHFeEmbWeightTrkPt","Reco Non-HF electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
            fRecoNonHFeEmbWeightTrkPt->Sumw2();
            fOutputList->Add(fRecoNonHFeEmbWeightTrkPt);
        
        fRecoNonHFeTrkPt_SPD= new TH2F("fRecoNonHFeTrkPt_SPD"," Reco Non-HF electrons from all generators;p_{T} (GeV/c);counts",250,0,50,300,-0.5,299.5);
        fRecoNonHFeTrkPt_SPD->Sumw2();
        fOutputList->Add(fRecoNonHFeTrkPt_SPD);
        
			fRecoNonHFeEmbTrkPt_SPD= new TH2F("fRecoNonHFeEmbTrkPt_SPD","Reco Non-HF electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,50,300,-0.5,299.5);
            fRecoNonHFeEmbTrkPt_SPD->Sumw2();
            fOutputList->Add(fRecoNonHFeEmbTrkPt_SPD);
            
			fRecoNonHFeEmbWeightTrkPt_SPD= new TH2F("fRecoNonHFeEmbWeightTrkPt_SPD","Reco Non-HF electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50,300,-0.5,299.5);
            fRecoNonHFeEmbWeightTrkPt_SPD->Sumw2();
            fOutputList->Add(fRecoNonHFeEmbWeightTrkPt_SPD);
        
        
            fRecoPi0eEmbWeightTrkPt = new TH1F("fRecoPi0eEmbWeightTrkPt","Reco Non-HF electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
            fRecoPi0eEmbWeightTrkPt->Sumw2();
            fOutputList->Add(fRecoPi0eEmbWeightTrkPt);
        
            fRecoEtaeEmbWeightTrkPt = new TH1F("fRecoEtaeEmbWeightTrkPt","Reco Non-HF electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
            fRecoEtaeEmbWeightTrkPt->Sumw2();
            fOutputList->Add(fRecoEtaeEmbWeightTrkPt);
        
        
         fInvmassLSPt = new TH2F("fInvmassLSPt", "Invmass of LS (e,e) for pt^{e}>1; p_{T}(GeV/c); mass(GeV/c^2); counts;", 500,0,50,500,0,1.0);
    fOutputList->Add(fInvmassLSPt);
    
    fInvmassULSPt = new TH2F("fInvmassULSPt", "Invmass of ULS (e,e) for pt^{e}>1; p_{T}(GeV/c); mass(GeV/c^2); counts;", 500,0,50,500,0,1.0);
    fOutputList->Add(fInvmassULSPt);
    
     fULSElecPt  = new TH1F("fULSElecPt","p_{T} distribution of ULS electrons;p_{T} (GeV/c);counts",500,0,50);
    fOutputList->Add(fULSElecPt);
    
    fLSElecPt= new TH1F("fLSElecPt","p_{T} distribution of LS electrons;p_{T} (GeV/c);counts",500,0,50);
    fOutputList->Add(fLSElecPt);





	 fNonHFePairInvmassLS = new TH1F("fNonHFePairInvmassLS", "Inv mass of LS (e,e) if both e- are Non-HFE; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fNonHFePairInvmassLS);
        
        fNonHFePairInvmassULS = new TH1F("fNonHFePairInvmassULS", "Inv mass of ULS (e,e) if both e- are Non-HFE; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fNonHFePairInvmassULS);
    
	//====================if(fEffiFromEnhMC)==================

       fNonHFeEmbInvmassLS = new TH1F("fNonHFeEmbInvmassLS", "Inv mass of LS (e,e) for Non-HFE from embedded #pi^{0} and #eta; mass(GeV/c^2); counts;",  50,0,0.5);
            fNonHFeEmbInvmassLS->Sumw2();
            fOutputList->Add(fNonHFeEmbInvmassLS);
        
            fNonHFeEmbInvmassULS = new TH1F("fNonHFeEmbInvmassULS", "Inv mass of ULS (e,e) for Non-HFE from embedded #pi^{0} and #eta; mass(GeV/c^2); counts;",  50,0,0.5);
            fNonHFeEmbInvmassULS->Sumw2();
            fOutputList->Add(fNonHFeEmbInvmassULS);
        
            fNonHFeEmbWeightInvmassLS = new TH1F("fNonHFeEmbWeightInvmassLS", "Inv mass of LS (e,e) for Non-HFE from embedded #pi^{0} and #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
            fNonHFeEmbWeightInvmassLS->Sumw2();
            fOutputList->Add(fNonHFeEmbWeightInvmassLS);
        
            fNonHFeEmbWeightInvmassULS = new TH1F("fNonHFeEmbWeightInvmassULS", "Inv mass of ULS (e,e) for Non-HFE from embedded #pi^{0} and #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
            fNonHFeEmbWeightInvmassULS->Sumw2();
            fOutputList->Add(fNonHFeEmbWeightInvmassULS);
        
            fPi0EmbInvmassLS = new TH1F("fPi0EmbInvmassLS", "Inv mass of LS (e,e) for ele from embedded #pi^{0}; mass(GeV/c^2); counts;",  50,0,0.5);
            fPi0EmbInvmassLS->Sumw2();
            fOutputList->Add(fPi0EmbInvmassLS);
        
            fPi0EmbInvmassULS  = new TH1F("fPi0EmbInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #pi^{0}; mass(GeV/c^2); counts;",  50,0,0.5);
            fPi0EmbInvmassULS->Sumw2();
            fOutputList->Add(fPi0EmbInvmassULS);
        
            fPi0EmbWeightInvmassLS = new TH1F("fPi0EmbWeightInvmassLS", "Inv mass of LS (e,e) for ele from embedded #pi^{0} with weight; mass(GeV/c^2); counts;",  50,0,0.5);
            fPi0EmbWeightInvmassLS->Sumw2();
            fOutputList->Add(fPi0EmbWeightInvmassLS);
        
            fPi0EmbWeightInvmassULS  = new TH1F("fPi0EmbWeightInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #pi^{0} with weight; mass(GeV/c^2); counts;",  50,0,0.5);
            fPi0EmbWeightInvmassULS->Sumw2();
            fOutputList->Add(fPi0EmbWeightInvmassULS);
        
            fEtaEmbInvmassLS = new TH1F("fEtaEmbInvmassLS", "Inv mass of LS (e,e) for ele from embedded #eta; mass(GeV/c^2); counts;",  50,0,0.5);
            fEtaEmbInvmassLS->Sumw2();
            fOutputList->Add(fEtaEmbInvmassLS);
        
            fEtaEmbInvmassULS = new TH1F("fEtaEmbInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #eta; mass(GeV/c^2); counts;",  50,0,0.5);
            fEtaEmbInvmassULS->Sumw2();
            fOutputList->Add(fEtaEmbInvmassULS);
        
            fEtaEmbWeightInvmassLS = new TH1F("fEtaEmbWeightInvmassLS", "Inv mass of LS (e,e) for ele from embedded #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
            fEtaEmbWeightInvmassLS->Sumw2();
            fOutputList->Add(fEtaEmbWeightInvmassLS);
        
            fEtaEmbWeightInvmassULS  = new TH1F("fEtaEmbWeightInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
            fEtaEmbWeightInvmassULS->Sumw2();
            fOutputList->Add(fEtaEmbWeightInvmassULS);
        
            fRecoLSeEmbTrkPt  = new TH1F("fRecoLSeEmbTrkPt","Reco LS electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,50);
            fRecoLSeEmbTrkPt->Sumw2();
            fOutputList->Add(fRecoLSeEmbTrkPt);
        
            fRecoLSeEmbWeightTrkPt = new TH1F("fRecoLSeEmbWeightTrkPt","Reco LS electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
            fRecoLSeEmbWeightTrkPt->Sumw2();
            fOutputList->Add(fRecoLSeEmbWeightTrkPt);
        
            fRecoPi0LSeEmbWeightTrkPt = new TH1F("fRecoPi0LSeEmbWeightTrkPt","Reco LS electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
            fRecoPi0LSeEmbWeightTrkPt->Sumw2();
            fOutputList->Add(fRecoPi0LSeEmbWeightTrkPt);
        
            fRecoEtaLSeEmbWeightTrkPt  = new TH1F("fRecoEtaLSeEmbWeightTrkPt","Reco LS electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
            fRecoEtaLSeEmbWeightTrkPt->Sumw2();
            fOutputList->Add(fRecoEtaLSeEmbWeightTrkPt);
        
            fRecoULSeEmbTrkPt = new TH1F("fRecoULSeEmbTrkPt","Reco ULS electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,50);
            fRecoULSeEmbTrkPt->Sumw2();
            fOutputList->Add(fRecoULSeEmbTrkPt);
        
            fRecoULSeEmbWeightTrkPt = new TH1F("fRecoULSeEmbWeightTrkPt","Reco ULS electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
            fRecoULSeEmbWeightTrkPt->Sumw2();
            fOutputList->Add(fRecoULSeEmbWeightTrkPt);
        
            fRecoPi0ULSeEmbWeightTrkPt = new TH1F("fRecoPi0ULSeEmbWeightTrkPt","Reco ULS electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
            fRecoPi0ULSeEmbWeightTrkPt->Sumw2();
            fOutputList->Add(fRecoPi0ULSeEmbWeightTrkPt);
        
            fRecoEtaULSeEmbWeightTrkPt = new TH1F("fRecoEtaULSeEmbWeightTrkPt","Reco ULS electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
            fRecoEtaULSeEmbWeightTrkPt->Sumw2();
            fOutputList->Add(fRecoEtaULSeEmbWeightTrkPt);
       
 
     
     } 
   //--------------------------------------------------------------------------------------------------
 //--------SPD correction---------------------------------
  
  const char *estimatorName="tr"; //SPD  
  
  fSPD_tracklet = new TH1F("fSPD_tracklet", "fSPD_tracklet Distribution",300,-0.5,299.5);
  fSPD_tracklet->Sumw2(); 
  fOutputList->Add(fSPD_tracklet);
  
  fSPDCorrMultDist_max = new TH1F("fSPDCorrMultDist_max",Form("Corrected Mult Dist; N_{%s};Counts;",estimatorName),300,-0.5,299.5); //
  fSPDCorrMultDist_max->Sumw2();
  fOutputList->Add(fSPDCorrMultDist_max);
  
  fSPDWeightedCorrMultDist_max= new TH1F("fSPDWeightedCorrMultDist_max",Form("Corrected Mult Dist; N_{%s};Counts;",estimatorName),300,-0.5,299.5); //
  fSPDWeightedCorrMultDist_max->Sumw2();
  if(fIsMC)fOutputList->Add(fSPDWeightedCorrMultDist_max);
  
   Profile_Mean = new TProfile("Profile_Mean","Profile_Mean",300,-15.0,15.0);
  Profile_Mean->Sumw2();
  fOutputList->Add(Profile_Mean); 
  
  Profile_MeanCorr = new TProfile("Profile_MeanCorr","Profile_MeanCorr",300,-15.0,15.0);
  Profile_MeanCorr->Sumw2();
  fOutputList->Add(Profile_MeanCorr);
   
  //-----------------------------------------------------------------
  fNentries2=new TH1F("CutSet", "", 32,-0.5,31.5);
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
  fNentries2->GetXaxis()->SetBinLabel(13,"InvmassCut");
  fNentries2->GetXaxis()->SetBinLabel(14,"AssoTPCCluster");  
  fNentries2->GetXaxis()->SetBinLabel(15,"AssoITSRefit");
  fNentries2->GetXaxis()->SetBinLabel(16,"AssopTMin");
  fNentries2->GetXaxis()->SetBinLabel(17,"AssoEtarange");
  fNentries2->GetXaxis()->SetBinLabel(18,"AssoTPCnsig");
  fNentries2->GetXaxis()->SetBinLabel(19,"fEMCEG1");
  fNentries2->GetXaxis()->SetBinLabel(20,"fDCalDG1");
  fNentries2->GetXaxis()->SetBinLabel(21,"fEMCEG2");
  fNentries2->GetXaxis()->SetBinLabel(22,"fDCalDG2");
  fNentries2->GetXaxis()->SetBinLabel(23,"fUseTender");
  fNentries2->GetXaxis()->SetBinLabel(24,"IsM20");
  fNentries2->GetXaxis()->SetBinLabel(25,"fCutM20Min");
  fNentries2->GetXaxis()->SetBinLabel(26,"fCutM20Max1");
  fNentries2->GetXaxis()->SetBinLabel(27,"fCutM20Max2");
  fNentries2->GetXaxis()->SetBinLabel(28,"fCutM20Max3");
  fNentries2->GetXaxis()->SetBinLabel(29,"fCutEopEMin");
  fNentries2->GetXaxis()->SetBinLabel(30,"fCutEopEMax");
  fNentries2->GetXaxis()->SetBinLabel(31,"deltaeta");
  fNentries2->GetXaxis()->SetBinLabel(32,"deltaphi");
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
  fNentries2->SetBinContent(13,fInvmassCut);
  fNentries2->SetBinContent(14,fAssoTPCCluster);  
  fNentries2->SetBinContent(15,fAssoITSRefit);
  fNentries2->SetBinContent(16,fAssopTMin);
  fNentries2->SetBinContent(17,fAssoEtarange);
  fNentries2->SetBinContent(18,fAssoTPCnsig);
  fNentries2->SetBinContent(19,fEMCEG1);
  fNentries2->SetBinContent(20,fDCalDG1);
  fNentries2->SetBinContent(21,fEMCEG2);
  fNentries2->SetBinContent(22,fDCalDG2);
  fNentries2->SetBinContent(23,fUseTender);
  fNentries2->SetBinContent(24,IsM20);
  fNentries2->SetBinContent(25,fCutM20Min);
  fNentries2->SetBinContent(26,fCutM20Max1);
  fNentries2->SetBinContent(27,fCutM20Max2);
  fNentries2->SetBinContent(28,fCutM20Max3);
  fNentries2->SetBinContent(29,fCutEopEMin);
  fNentries2->SetBinContent(30,fCutEopEMax);
  fNentries2->SetBinContent(31,fdeltaeta);
  fNentries2->SetBinContent(32,fdeltaphi);
 

   
  //==================================================================================================//
  //==================================================================================================//
  					
  const char* nameoutput=GetOutputSlot(2)->GetContainer()->GetName();  // gets the output slot from mgr-> in AddTask
  fNentries=new TH1F(nameoutput, "", 6,-0.5,5.5);
  fNentries->GetXaxis()->SetBinLabel(1,"No. of Events Analyzed");
  fNentries->GetXaxis()->SetBinLabel(2,"No. of Events Accepted");
  fNentries->GetXaxis()->SetBinLabel(3,"trigger check");
  fNentries->GetXaxis()->SetBinLabel(4,"After Pile Up cuts");
  fNentries->GetXaxis()->SetBinLabel(5,"After Z vtx cut");
  fNentries->GetXaxis()->SetBinLabel(6,"After N contributors");
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

void AliAnalysisTaskHFEmultTPCEMCAL::UserExec(Option_t *) 
{
   
	fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
	if (!fAOD) { return;}
  
	fNentries->Fill(0);                //No. of Events Analyzed 
	Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & ftrigger);
	if(!isSelected) return;  
	fNentries->Fill(1); 
	if(fUseTender){
    	fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTenderTrackName));
    	fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTenderClusterName)); //emcal correction
  	}
	  /////////////////
  	//trigger check//
  	/////////////////
  	TString firedTrigger;
  	TString TriggerEG1 = "EG1", TriggerEG2 = "EG2", TriggerDG1 = "DG1", TriggerDG2 = "DG2";
  	if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();

  	if(fEMCEG2 && fDCalDG2) if(!firedTrigger.Contains(TriggerEG2) && !firedTrigger.Contains(TriggerDG2)) return;
  	if(fEMCEG1 && fDCalDG1) if(!firedTrigger.Contains(TriggerEG1) && !firedTrigger.Contains(TriggerDG1)) return;

	//cout<<" ftrigger "<<ftrigger<<" fEMCEG1 "<<fEMCEG1<<"   fDCalDG1 "<<fDCalDG1<<"   fEMCEG2 "<<fEMCEG2<<"     fDCalDG2 "<<fDCalDG2<<endl;
	//getchar();   
	
	fNentries->Fill(2);                 //No. of Events Accepted
		
	const AliAODVertex *pVtx= fAOD->GetPrimaryVertex();
	if(!pVtx){return;}
	TString vtxTtl = pVtx->GetTitle();
	//if(!vtxTtl.Contains("VertexerTracks")) {printf("No PrimaryVertexTracks\n");}
	
	Double_t Zvertex=-999,Xvertex=-999,Yvertex=-999;
  
	Zvertex = pVtx->GetZ();
	Yvertex = pVtx->GetY();
	Xvertex = pVtx->GetX();
	fVtxZ->Fill(Zvertex);
	  	
	if(pVtx->GetNContributors()<2) return;
	
	fNentries->Fill(3);
	
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
	fNentries->Fill(4);
		
	if(TMath::Abs(Zvertex)>10){return;} 
	fNentries->Fill(5);
  
	fVtxZ_corr->Fill(Zvertex);
	
			
	
	fpidResponse = fInputHandler->GetPIDResponse();
	if(!fpidResponse)
	{
		AliDebug(1, "Using default PID Response");
		fpidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
		fPID->SetPIDResponse(fpidResponse);
	}
	fPID->SetPIDResponse(fpidResponse);
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////==================================SPD Tracklets==============================================/////
//////==================================================================================================/////
	Int_t nAcceta=0.;
	AliAODTracklets *tracklets = ((AliAODEvent*)fAOD)->GetTracklets();
	Int_t nTracklets = tracklets->GetNumberOfTracklets();
	for (Int_t nn = 0; nn < nTracklets; nn++)
	{
		Double_t theta = tracklets->GetTheta(nn);
		Double_t eta = -TMath::Log(TMath::Tan(theta/2.0));
		if (TMath::Abs(eta) < 1.) nAcceta++;
	}


 	Double_t N_corr_tr_eta1=nAcceta;  
  TProfile* estimatorAvg = GetEstimatorHistogram(fAOD);
  if(estimatorAvg){
	     N_corr_tr_eta1=static_cast<Int_t>(GetCorrectedNtracklets(estimatorAvg,nAcceta,Zvertex,fRefMult)); 
  } 
  SPDntr = N_corr_tr_eta1;
  
  fSPD_tracklet->Fill(nAcceta);
  fSPDCorrMultDist_max->Fill(SPDntr);
  	
  Profile_Mean->Fill(Zvertex,nAcceta);
  Profile_MeanCorr->Fill(Zvertex,SPDntr);
  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////========================================MC LOOP===================================================/////
//////==================================================================================================/////  
	if(fIsMC) 
	{    
		Double_t qaweights[5];
		Double_t pi0etaweights[3];
		
		fMCArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
		fMCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
		 if(fMCHeader){
        ////////////////////////////////
        //Get number of Gen particles //
        ////////////////////////////////
        GetNMCPartProduced();
        
        /////////////////////////////////
        //Calculate Pi0 and Eta weight //
        /////////////////////////////////
        fCalculateWeight=kFALSE;
        //if(fCalculateWeight) 
        GetPi0EtaWeight(fSprsPi0EtaWeightCal);
    
    }
		if(!fMCArray)
		{
			AliError("Array of MC particles not found");
		 	return;
		}
		for(Int_t iMC = 0; iMC < fMCArray->GetEntries(); iMC++)
		{

			fMCparticle = (AliAODMCParticle*) fMCArray->At(iMC);
			pdg = fMCparticle->GetPdgCode();

			///For the reconstruction efficiency of HFE:-----------
			if(TMath::Abs(fMCparticle->Eta()) <=fEtarange  )
			{			
				if(TMath::Abs(pdg) == 11)
				{		
					Int_t IsElecHf=GetHFE(fMCparticle,fMCArray);
					if((IsElecHf==kBeauty) || (IsElecHf==kCharm)){
						fPtHFEMC->Fill(fMCparticle->Pt());
						fPtHFEMC_SPD->Fill(fMCparticle->Pt(),SPDntr);
					
					}
				}
			} //eta <0.7 condition
			
			//pt spectra for pi0 and eta
			//if(TMath::Abs(fMCparticle->Y())>0.5 ) continue;
			if(fMCparticle->Y()<-0.8 || fMCparticle->Y()>0.8) continue;
			
 			if(fMCparticle->IsPrimary()) ///Does not include particles from weak decays or created in an interaction with the material 
			{ 
	
			
					if (TMath::Abs(fMCparticle->GetPdgCode())==111) pi0MC->Fill(fMCparticle->Pt());  // pdg=111 pi0
					if (TMath::Abs(fMCparticle->GetPdgCode())==221) etaMC->Fill(fMCparticle->Pt());  // pdg=221 eta	
					if (TMath::Abs(fMCparticle->GetPdgCode())==22) gammaMC->Fill(fMCparticle->Pt());   // pdg=22  gamma
			
					Int_t type = GetPi0EtaType(fMCparticle);
				
			
			
					///Using thnSparse--------------------------------------
					//Pt
				
					pi0etaweights[0] = fMCparticle->Pt();
						
					// What pdg
					pi0etaweights[1]=-1.;
					if (TMath::Abs(fMCparticle->GetPdgCode())==111) pi0etaweights[1]=0.2;  // pdg=111 pi0
					if (TMath::Abs(fMCparticle->GetPdgCode())==221) pi0etaweights[1]=1.2;  // pdg=221 eta
					if (TMath::Abs(fMCparticle->GetPdgCode())==22) pi0etaweights[1]=2.2;   // pdg=22  gamma
				
								
					// What type
					//Int_t type = GetPi0EtaType(fMCparticle);
				
					pi0etaweights[2]=type;
					
				
					if(pi0etaweights[1]>0.) fPi0EtaSpectraSp->Fill(pi0etaweights);
				
			
				///-----------------------------------------------------
		
		      } 	//IsPrimary()loop	
					
					
		////====================================================
		} //For loop
	}

	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////==================================EMCAL Cluster LOOP==============================================/////
//////==================================================================================================/////

	

  Double_t cluphi = -999.0;
  Double_t clueta =-999.0 ;
  Int_t ncells = -999.0;
  Float_t energy = -999.0;
  Float_t clut = -999.0;
  Double_t  energycell = -999.0;
  Double_t CellId =0;
  Int_t Nclust = -999; 
	if(!fUseTender) Nclust = fAOD->GetNumberOfCaloClusters(); 
  if(fUseTender) Nclust = fCaloClusters_tender->GetEntries();

  Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;
  for ( Int_t index = 0; index < Nclust ; index++ ) 
  {
				AliAODCaloCluster * clu =0x0;
				if(!fUseTender) clu  = (AliAODCaloCluster*)fAOD->GetCaloCluster(index) ; 
				if(fUseTender) clu = dynamic_cast<AliAODCaloCluster*>(fCaloClusters_tender->At(index));
				if(!clu) continue;
			
				fClsTypeEMC = kFALSE; 
				fClsTypeDCAL = kFALSE;
				if (clu->IsEMCAL())
				{
	
				  AliAODCaloCells &cells = *(fAOD->GetEMCALCells());
	
				  Float_t  x[3]; // cluster pos
				  clu->GetPosition(x);
				  TVector3 clustposi(x[0],x[1],x[2]);
			 
				  cluphi = clustposi.Phi();
				  clueta = clustposi.Eta();
				  if(cluphi < 0) cluphi = cluphi+(2*TMath::Pi());
				  if(cluphi > 1.39 && cluphi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
				  if(cluphi > 4.53 && cluphi < 5.708) fClsTypeDCAL = kTRUE;//DCAL  : 260 < phi < 327

				  if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
					if(!fClsTypeEMC) continue; //selecting only EMCAL clusters

				  if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
					if(!fClsTypeDCAL) continue; //selecting only DCAL clusters

				  clut = clu->GetTOF()*1e9 ;
				  energy = clu->E();
				  ncells= clu->GetNCells();
				  
				  //fClusPhi->Fill(cluphi);
				  //fClusEta->Fill(clueta);
				  fClusEtaPhi->Fill(clueta,cluphi); 
				  fNCells->Fill(ncells);
				  fClusE->Fill(energy);
				  fClusT->Fill(clut);
				  fClusEvsnTracklets->Fill(energy,SPDntr);//For Rejection Factor
					/*
				  fvalueCluE[0] = energy;
					fvalueCluE[1] = vzeroMultCorr; //V0M, Multiplicity information
					fvalueCluE[2] = correctednAcc1; //SPD Tracklets

					fSparseClusE->Fill(fvalueCluE); //For Rejection Factor
					*/
				}  

	}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////======================================TRACK LOOP==================================================/////
//////==================================================================================================/////

	Int_t ntracks = -999;
  if(!fUseTender)ntracks = fAOD->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries(); 
  
  Int_t runNumber = fAOD->GetRunNumber();
	Int_t num=0;
	for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) 
	{
			AliAODTrack* track = 0x0;
		 	if(!fUseTender) track = (AliAODTrack*)fAOD->GetTrack(iTracks);
    	   if(fUseTender) track =  dynamic_cast<AliAODTrack*>(fTracks_tender->At(iTracks));
			if(!track) AliFatal("Not a standard AOD"); 
			
			
  			Int_t tracktypeTrig=0;
    		tracktypeTrig=ClassifyTrack(track,pVtx);  //track classify
   		if(tracktypeTrig!=1) continue;  //if(tracktype==0) continue; if(tracktype==1)  	
     		
			num++;
    	
			Double_t pt=0., p=0., eta=-999, dEdx_TPC=-999., Signal_TOF=-999., TOFbeta=-999., fTPCnSigma=-999.0, fTOFnSigma=-999.0;    
		  
			pt = track->Pt();
			p = track->P();  
			eta=track->Eta();
			
			
																		//reconstructed level & have proper TPC PID response(?)
     		if(fIsMC && track->GetLabel()>=0)
     		{		
					//cout<<" track->GetLabel() = "<<track->GetLabel()<<endl;
					//getchar();
					fMCparticle=(AliAODMCParticle*)fMCArray->At(track->GetLabel());
					pdg = fMCparticle->GetPdgCode();
					if(TMath::Abs(pdg) == 11)
					{													
     				Int_t IsElecHf=GetHFE(fMCparticle,fMCArray);
					if((IsElecHf==kBeauty) || (IsElecHf==kCharm))
					{
						fPtHFEMC_trackcutreco->Fill(pt);
						fPtHFEMC_trackcutreco_SPD->Fill(pt,SPDntr);
					}	
					}	
			} 	
			
			feta->Fill(eta);
			dEdx_TPC = track->GetTPCsignal();
			Double_t phi = track->Phi();
		   	
			//TOFbeta=Beta(track);
			Double_t weight_lan=0.,weight_lan_inv=0.;
			Double_t weight_err=0.,weight_err_inv=0.;
	
			fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
		
			fdEdxVsP_TPC->Fill(p,dEdx_TPC);
			fnSigmaVsP_TPC->Fill(p,fTPCnSigma);
			fnSigmaVsPt_TPC->Fill(pt,fTPCnSigma);
			
			///////////////////////////
			//Track matching to EMCAL//
			//////////////////////////
		  Int_t EMCalIndex = -1;
			EMCalIndex = track->GetEMCALcluster();
			if(EMCalIndex < 0) continue;
			fHistPtMatch->Fill(pt);
		
			
			AliVCluster *clustMatch=0x0;
		  if(!fUseTender) if(EMCalIndex >= 0) clustMatch = (AliAODCaloCluster*)fAOD->GetCaloCluster(EMCalIndex) ; 
    	if(fUseTender) if(EMCalIndex >= 0)clustMatch = dynamic_cast<AliAODCaloCluster*>(fCaloClusters_tender->At(EMCalIndex));
    
    	Double_t emcphi = -999, emceta=-999;
    	Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;
    	if(clustMatch && clustMatch->IsEMCAL())
    	{
    	 
    	 	
		  Double_t fPhiDiff = -999, fEtaDiff = -999;
		  GetTrkClsEtaPhiDiff(track, clustMatch, fPhiDiff, fEtaDiff);
		  fEMCTrkMatch->Fill(fPhiDiff,fEtaDiff);

		  if(TMath::Abs(fPhiDiff) > fdeltaeta || TMath::Abs(fEtaDiff)> fdeltaeta) continue;

		  /////////////////////////////////
		  //Select EMCAL or DCAL clusters//
		  /////////////////////////////////
		  Float_t  emcx[3]; // cluster pos
		  clustMatch->GetPosition(emcx);
		  
		  TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
		  emcphi = clustpos.Phi();
		  emceta = clustpos.Eta();
		  if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
		  if(emcphi > 1.39 && emcphi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
		  if(emcphi > 4.53 && emcphi < 5.708) fClsTypeDCAL = kTRUE;//DCAL  : 260 < phi < 327

		  //----selects EMCAL+DCAL clusters when fFlagClsTypeEMC and fFlagClsTypeDCAL is kTRUE
      if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
		  if(!fClsTypeEMC) continue; //selecting only EMCAL clusters

      if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
      if(!fClsTypeDCAL) continue; //selecting only DCAL clusters

      //Double_t clustTime = clustMatch->GetTOF()*1e+9; // ns;

		  //if(fEMCClsTimeCut)
		  //	if(TMath::Abs(clustTime) > 50) continue;

		  /////////////////////////////////////////////
		  //Properties of tracks matched to the EMCAL//
		  /////////////////////////////////////////////
		  
		  	
			if(fIsMC && track->GetLabel()>=0)
     		{			
					pdg = fMCparticle->GetPdgCode();
					if(TMath::Abs(pdg) == 11)
					{										
     				Int_t IsElecHf=GetHFE(fMCparticle,fMCArray);
					if((IsElecHf==kBeauty) || (IsElecHf==kCharm))
					{
						fPtHFEMC_trackmatchreco->Fill(pt);
						fPtHFEMC_trackmatchreco_SPD->Fill(pt,SPDntr);
						
					}	
					}
			} 
	
		  fEMCTrkPt->Fill(pt);
		  fEMCTrketa->Fill(eta);
		  fEMCTrkphi->Fill(phi);
		  fEMCTPCnsig->Fill(p,fTPCnSigma);
		  
		  
		  Double_t clustMatchE = -999.0, Eoptrk = -999.0 , M02trkmatch = -999.0, M20trkmatch = -999.0;
  
		  clustMatchE = clustMatch->E();
		  if(p>0)Eoptrk = clustMatchE/p;
		  
		  fClsEAftMatch->Fill(clustMatchE); // For RF
		  fClsEAftMatch_SPD->Fill(clustMatchE,SPDntr); // For RF
		  fClsEopAftMatch->Fill(Eoptrk);
		  fClsEtaPhiAftMatch->Fill(emceta,emcphi);
	
		  M02trkmatch = clustMatch->GetM02();
		  M20trkmatch = clustMatch->GetM20();
 
		
       if(pt < 2) continue;   
		  fvalueElectron[0] = pt; //matched tracks pt
		  fvalueElectron[1] = fTPCnSigma; // tpc n sigma
		  fvalueElectron[2] = Eoptrk; //E/P
		  fvalueElectron[3] = M02trkmatch; // shower shape cut
		  //fvalueElectron[4] = M20trkmatch;
		  fvalueElectron[4] = SPDntr; //
		  /*fvalueElectron[6] = 0; //SPD Tracklets
		  fvalueElectron[7] = clustMatchE;  //cluster energy after matching
		  fvalueElectron[8] = Zvertex;*/
						
			fSparseElectron->Fill(fvalueElectron);   //Electron information sparse         
	 	  //////////////////
		  //Apply EID cuts//
		  //////////////////	
	
			Bool_t fElectTrack = kFALSE;
			
			 if(Eoptrk < fCutEopEMin || Eoptrk > fCutEopEMax) continue;			
			if(fIsMC && track->GetLabel()>=0)
     		{			
					pdg = fMCparticle->GetPdgCode();
					if(TMath::Abs(pdg) == 11)
					{											
     				Int_t IsElecHf=GetHFE(fMCparticle,fMCArray);
					
					if((IsElecHf==kBeauty) || (IsElecHf==kCharm))
					{
						
						fPtHFEMC_TPCEMCreco->Fill(pt);
						fPtHFEMC_TPCEMCreco_SPD->Fill(pt,SPDntr);
					}	
					}	
			}



			Double_t nsigma_ele=-999;
  			nsigma_ele = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
  			if(nsigma_ele < fTPCnsigmin || nsigma_ele > fTPCnsigmax) continue;
 
 			if(fIsMC && track->GetLabel()>=0)
     		{			
					pdg = fMCparticle->GetPdgCode();
					if(TMath::Abs(pdg) == 11)
					{											
     				Int_t IsElecHf=GetHFE(fMCparticle,fMCArray);
					
					if((IsElecHf==kBeauty) || (IsElecHf==kCharm))
					{
						
						fPtHFEMC_TPCreco->Fill(pt);
						fPtHFEMC_TPCreco_SPD->Fill(pt,SPDntr);
					}	
					}	
			} 
			
	  
  
		  
		  
		 

		 //if(M02trkmatch < fCutM20Min || M02trkmatch > fCutM20Max) return kFALSE;
		
		  if(pt<12){if(M02trkmatch < fCutM20Min || M02trkmatch > fCutM20Max1) continue;}
		  if(pt>=12 && pt<20){if(M02trkmatch < fCutM20Min || M02trkmatch > fCutM20Max2) continue;}
		  if(pt>=20){if(M02trkmatch < fCutM20Min || M02trkmatch > fCutM20Max3) continue;}
		  
		  
		  if(fIsMC && track->GetLabel()>=0)
     		{			
					pdg = fMCparticle->GetPdgCode();
					if(TMath::Abs(pdg) == 11)
					{											
     				Int_t IsElecHf=GetHFE(fMCparticle,fMCArray);
					
					if((IsElecHf==kBeauty) || (IsElecHf==kCharm))
					{
						
						fPtHFEMC_SScutreco->Fill(pt);
						fPtHFEMC_SScutreco_SPD->Fill(pt,SPDntr);
					}	
					}	
			}  
      	
				//fElectTrack = PassEIDCuts(track, clustMatch);     	
      	//if(!fElectTrack) continue;
      
           //////////////////
		  //=======Non-HFE================= //
		  //////////////////      
		 
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
				/*if(pt>=0.5 && pt<1.5){
				  
				fNonHFE->SetHistMassBack(fInvmassLS1);
				fNonHFE->SetHistMass(fInvmassULS1);
				}
				
				else if(pt>=1.5 && pt<2.5){
			
				fNonHFE->SetHistMassBack(fInvmassLS1_2);
				fNonHFE->SetHistMass(fInvmassULS1_2);
				}
				else if(pt>=2.5 && pt<3.5){
				fNonHFE->SetHistMassBack(fInvmassLS1_3);
				fNonHFE->SetHistMass(fInvmassULS1_3);
				}
				*/
				if(pt>=3.5 && pt<4.5){
				fNonHFE->SetHistMassBack(fInvmassLS1);
				fNonHFE->SetHistMass(fInvmassULS1);
				}
				fNonHFE->FindNonHFE(iTracks,track,fAOD,fTracks_tender,fUseTender);
				//fNonHFE->FindNonHFE(iTracks,track,fAOD);
				
			
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
  	
  				
  	//--------------------MC---------------------------------------
  				if(fIsMC && track->GetLabel()>=0)
   			{
   				
   				fMCparticle=(AliAODMCParticle*)fMCArray->At(track->GetLabel());
					pdg = fMCparticle->GetPdgCode();
					Float_t ptMC= fMCparticle->Pt();
					
					fPT_PIDCaloCut_MC->Fill(ptMC);
					fPT_PIDCaloCut_realistic->Fill(pt);
					
					if(TMath::Abs(pdg)!=11) continue ;      //Is electron:
					fPT_elec_MC->Fill(ptMC);
					fPT_elec_realistic->Fill(pt);					
					
					//cout<<" ptMC "<<ptMC<<"   pt "<<pt<<endl; 
					//Example of outputs=======
					//Usually they are close  but sometimes there is a massive difference.
					// ptMC 2.97147   pt 2.70653
					//ptMC 3.76586   pt 3.80316
					//ptMC 17.4227   pt 5.58553
					//getchar();
					
					Int_t fMCmotherindex=fMCparticle->GetMother(); //Getting Electron Mother
					if(fMCmotherindex<0) continue ;
					fMCmother= (AliAODMCParticle*)fMCArray->At(fMCparticle->GetMother());
					Int_t pdg_mother = fMCmother->GetPdgCode();
					/*
					//----If mother is Pi0 , eta or gamma --------------------------------------------------------------
					if(TMath::Abs(pdg_mother) == 111 || TMath::Abs(pdg_mother) == 22 || TMath::Abs(pdg_mother) == 221)
					{
						Double_t ptmotherw = -1.;	
						Int_t elec_source = GetElecSourceType(fMCparticle,ptmotherw);
						if((elec_source==kPi0NoFeedDown) || (elec_source==kGPi0NoFeedDown) || (elec_source==kEtaNoFeedDown) || (elec_source==kGEtaNoFeedDown))
						{						
							fPt_elec_phot1->Fill(pt);
							fPt_elec_phot1_multSPD->Fill(pt,SPDntr);
							SelectPhotonicElectronR(iTracks, track, fMCmotherindex, pdg,elec_source, SPDntr,ptmotherw);
						}
					}
  					*/
  					Int_t IsElecHf=GetHFE(fMCparticle,fMCArray);
					
					if((IsElecHf==kBeauty) || (IsElecHf==kCharm))
					{
						
						fPtHFEMC_reco->Fill(pt);
						fPtHFEMC_reco_SPD->Fill(pt,SPDntr);
						
					}	
				
				
				
            //////////////////////////////////
            //Non-HFE efficiency calculation//
            //////////////////////////////////
            Bool_t EffiDenom = kFALSE;
            Bool_t EffiNumTag = kFALSE;
            //if(fMCHeader && fCalculateNonHFEEffi){
            if(fMCHeader ){
                //if(fEffiFromEnhMC) 
                EffiDenom = GetNonHFEEffiDenom(track);
                //if(!fEffiFromEnhMC) EffiDenom = GetNonHFEEffiDenomGenPurMC(track);
            }
            
            ////////////////////
            //NonHFE selection//
            ////////////////////
            Bool_t fFlagNonHFE=kFALSE;
            Int_t pidM = -1;
            
            SelectPhotonicElectron(iTracks,track,fFlagNonHFE,pidM);
             
           
            //////////////////////////////////
            //Non-HFE efficiency calculation//
            //////////////////////////////////
            //if(fMCHeader && fCalculateNonHFEEffi){
            if(fMCHeader){
                if(fFlagNonHFE){
                    //if(fEffiFromEnhMC) 
                    EffiNumTag = GetNonHFEEffiRecoTag(track);
                    //if(!fEffiFromEnhMC) EffiNumTag = GetNonHFEEffiRecoTagGenPurMC(track);
                }
            }
          } //------------fIsMC Loop------------------  
            
		}//EMCAL track match
    

	}//---------------END of track loop----------------------------
	////////////////////////////////////////////////////////////////
	
	//fHistMult->Fill(num);
	//delete fListOfmotherkink;
}
void AliAnalysisTaskHFEmultTPCEMCAL::GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff)
{
	// Calculate phi and eta difference between a track and a cluster. The position of the track is obtained on the EMCAL surface

	phidiff = 999;
	etadiff = 999;

	if (!t||!v) return;

	Double_t veta = t->GetTrackEtaOnEMCal();
	Double_t vphi = t->GetTrackPhiOnEMCal();

	Float_t pos[3] = {0};
	v->GetPosition(pos);
	TVector3 cpos(pos);
	Double_t ceta     = cpos.Eta();
	Double_t cphi     = cpos.Phi();
	etadiff=veta-ceta;
	phidiff=TVector2::Phi_mpi_pi(vphi-cphi);
}

Bool_t AliAnalysisTaskHFEmultTPCEMCAL::PassEIDCuts(AliAODTrack *track, AliVCluster *clust)
{
  //apply electron identification cuts

  Double_t eop = -1.0;
  Double_t m02 = -999,m20 = -999;
  Double_t clustE = clust->E();
  Double_t TrkPt = track->Pt();
  //Double_t nsigma_ele=-999;
  //nsigma_ele = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
  if(track->P()>0)eop = clustE/track->P();
 
  //if(IsM20)m20 =clust->GetM20();
  //else if(!IsM20)m20 =clust->GetM02();
  
  m02 =clust->GetM02();
  
  //if(nsigma_ele < fTPCnsigmin || nsigma_ele > fTPCnsigmax) return kFALSE;
  
  
  //if(m02 < fCutM20Min || m02 > fCutM20Max) return kFALSE;
  if(TrkPt<12){if(m02 < fCutM20Min || m02 > 0.9) return kFALSE;}
  if(TrkPt>12){if(m02 < fCutM20Min || m02 > 0.7) return kFALSE;}
  if(eop < fCutEopEMin || eop > fCutEopEMax) return kFALSE;

  return kTRUE;
}



//--------------------END of UserExec----------------------------------
Int_t AliAnalysisTaskHFEmultTPCEMCAL::ClassifyTrack(AliAODTrack* track,const AliVVertex* pVtx)
{  
  
	Double_t pt = track->Pt();
	Double_t eta = track->Eta();
	Double_t phi = track->Phi();
	Float_t dx,dy,dxy, dz;
 	
 	
 	//====kink daughters
  Int_t numberofvertices = 100;
  numberofvertices = fAOD->GetNumberOfVertices();
  Double_t listofmotherkink[numberofvertices];
  Int_t numberofmotherkink = 0;
  for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
    AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
    if(!aodvertex) continue;
    if(aodvertex->GetType()==AliAODVertex::kKink) {
      AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
      if(!mother) continue;
      Int_t idmother = mother->GetID();
      listofmotherkink[numberofmotherkink] = idmother;
      numberofmotherkink++;
    }
  }
 
  //reject kink
  Bool_t kinkmotherpass = kTRUE;
  for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
    if(track->GetID() == listofmotherkink[kinkmother]) {
      kinkmotherpass = kFALSE;
      continue;
    }
  }
  if(!kinkmotherpass) return kFALSE;


 	//=====fitler bit
	if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return 0;
	//=====Pt cut  
  if(pt< 0.5) return 0; 
  //=====Eta cut   
	if (TMath::Abs(eta)>fEtarange) return 0; 
	//if (TMath::Abs(eta)>=0.7) return 0; 
	
	Double_t nclus = track->GetTPCNcls();  // TPC cluster information
	Double_t nclusF = track->GetTPCNclsF();
 	Double_t nclusN = track->GetTPCsignalN();  // TPC cluster information findable
 	Double_t RatioTPCclusters=nclusN/nclusF;

 	
 	//=====TPC Cluster, TPC PID cut, ITS clsuter, RatioTPCcluster============= 
	if(track->GetTPCNcls() < fTPCNclus) return 0; //TPC N clusters
	if(track->GetITSNcls() < fITSNclus) return 0; // ITS N clusters
	if(nclusN< fTPCNclusPID) return 0 ;
	if(RatioTPCclusters<0.6) return 0;
	
	//=========ITS TPC Refit=============
	if((!(track->GetStatus()&AliESDtrack::kITSrefit)|| (!(track->GetStatus()&AliESDtrack::kTPCrefit)))) return 0; // ITS and TPC refit
	
	
	//=====Hits on SPD layers=============
	
       if(fSPDBoth){ if(!(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1))) return 0;} //Hit on first and second SPD layer
	else if(fSPDAny){ if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) return 0;} //Hit on any layer
	else if(fSPDFirst){ if(!(track->HasPointOnITSLayer(0))) return 0;} //Hit on first and second SPD layer
  
   //if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) return 0;
    
	//=========DCA Cut ==================
	Double_t d0z0[2]={-999,-999}, cov[3];
	if(track->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov)){
		fDCAxy= d0z0[0];
		fDCAz= d0z0[1];
	}
		if(TMath::Abs(d0z0[0]) > fDCAxyCut || TMath::Abs(d0z0[1]) > fDCAzCut) return 0; 
	
  //==========chi2Ndf=================
	Double_t chi2ndf = track->Chi2perNDF();
	if(chi2ndf>4.0) return 0;
			
	return 1;
}

//_________________________________________

void AliAnalysisTaskHFEmultTPCEMCAL::SelectPhotonicElectronR(Int_t itrack, AliAODTrack *track, Int_t motherindex, Int_t pdg1, Int_t source , Double_t SPDntr1,Double_t ptmotherwg)
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
	
	Int_t ntracks = -999;
  if(!fUseTender)ntracks = fAOD->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries(); 

   
	for(Int_t jTracks = 0; jTracks< ntracks; jTracks++)
  	{
		AliVParticle* VtrackAsso = 0x0;
		if(!fUseTender) VtrackAsso = (AliVParticle*)fAOD->GetTrack(jTracks);
    	if(fUseTender) VtrackAsso  =  dynamic_cast<AliVParticle*>(fTracks_tender->At(jTracks));
	
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
        
        Double_t invms[4];
    	  invms[0]=track->Pt();
	     invms[1]=source;
        invms[2]=mass;
        invms[3]=SPDntr1;
	    // invms[4]=V0Mmult1;
        
							
        if(fFlagULS){
        	fInvMULS->Fill(track->Pt(),mass);
           	fInvMULSnSp->Fill(invms);
		}
		/*   if(fFlagULS){
        	fInvMULSnSpwg->Fill(invms);
        	}
		*/
	}
}


void AliAnalysisTaskHFEmultTPCEMCAL::SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC)
{
    ///////////////////////////////////////////
    //////Non-HFE - Invariant mass method//////
    ///////////////////////////////////////////
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    //const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    Double_t d0z0[2]={-999,-999}, cov[3];
    Double_t DCAxyCut = 0.25, DCAzCut = 1;
    
    Bool_t flagPhotonicElec = kFALSE, flagLSElec = kFALSE;
    Double_t ptAsso=-999., nsigma=-999.0;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
    
    Int_t ntracks = -999;
    //if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
    if(!fUseTender)ntracks = fAOD->GetNumberOfTracks();
    if(fUseTender) ntracks = fTracks_tender->GetEntries();
    
    for (Int_t jtrack = 0; jtrack < ntracks; jtrack++) {
        AliVParticle* VAssotrack = 0x0;
        //if(!fUseTender) VAssotrack  = fVevent->GetTrack(jtrack);
		if(!fUseTender) VAssotrack  = fAOD->GetTrack(jtrack);
        if(fUseTender) VAssotrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jtrack)); //take tracks from Tender list
        
        if (!VAssotrack) {
            printf("ERROR: Could not receive track %d\n", jtrack);
            continue;
        }
        
        AliVTrack *Assotrack = dynamic_cast<AliVTrack*>(VAssotrack);
        AliESDtrack *eAssotrack = dynamic_cast<AliESDtrack*>(VAssotrack);
        AliAODTrack *aAssotrack = dynamic_cast<AliAODTrack*>(VAssotrack);
        
        //------reject same track
        if(jtrack==itrack) continue;
        
        Double_t mass=-999., width = -999;
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
        
        nsigma = fpidResponse->NumberOfSigmasTPC(Assotrack, AliPID::kElectron);
        ptAsso = Assotrack->Pt();
        
        //------track cuts applied
        //if(fAOD) {
            if(!aAssotrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            if(aAssotrack->GetTPCNcls() < fAssoTPCCluster) continue;
            if((!(aAssotrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(aAssotrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
            
            //if(aAssotrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov))
			//if(aAssotrack->PropagateToDCA(pVtx,fAOD->GetMagneticField(), 20., d0z0, cov))
         //       if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;
        //}
        
        //-------loose cut on partner electron
        if(ptAsso < fAssopTMin) continue;
		if(TMath::Abs(aAssotrack->Eta())>fAssoEtarange) continue;
		if(TMath::Abs(nsigma) > fAssoTPCnsig ) continue;

        
        
        Int_t chargeAsso = Assotrack->Charge();
        Int_t charge = track->Charge();
        if(charge>0) fPDGe1 = -11;
        if(chargeAsso>0) fPDGe2 = -11;
        
        fFlagLS=kFALSE; fFlagULS=kFALSE;
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        //-------define KFParticle to get mass
        //AliKFParticle::SetField(fVevent->GetMagneticField());
		AliKFParticle::SetField(fAOD->GetMagneticField());
        AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
        AliKFParticle ge2 = AliKFParticle(*Assotrack, fPDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        Int_t MassCorrect;
        MassCorrect = recg.GetMass(mass,width);
        
        if(fFlagLS && track->Pt()>1) fInvmassLSPt->Fill(track->Pt(),mass);
        if(fFlagULS && track->Pt()>1) fInvmassULSPt->Fill(track->Pt(),mass);
        
        //////////////////////////////////
        //Non-HFE efficiency calculation//
        //////////////////////////////////
        Bool_t EffiNumULSLS = kFALSE;
        //if(fMCHeader && fCalculateNonHFEEffi)
        if(fMCHeader)
        {
           //if(fEffiFromEnhMC) 
           EffiNumULSLS = GetNonHFEEffiULSLS(track, Assotrack, fFlagLS, fFlagULS, mass);
            //if(!fEffiFromEnhMC) EffiNumULSLS = GetNonHFEEffiULSLSGenPurMC(track, Assotrack, fFlagLS, fFlagULS, mass);
        }

        Double_t TrkPt = track->Pt();
        if(mass < fInvmassCut){
            if(fFlagLS)
                fLSElecPt->Fill(TrkPt);

            if(fFlagULS)
                fULSElecPt->Fill(TrkPt);
        }
        
        if(mass < fInvmassCut && fFlagULS && !flagPhotonicElec)
            flagPhotonicElec = kTRUE; //Tag Non-HFE (random mass cut, not optimised)
    }
    fFlagPhotonicElec = flagPhotonicElec;
}
/*
Double_t AliAnalysisTaskHFEmultTPCEMCAL::WeightMCncCorr(Int_t multSPDtr)
{
Double_t  MCnchweight=1.;
if(multSPDtr==0) MCnchweight= 1.40268;
else if(multSPDtr==1) MCnchweight= 0.971848;
else if(multSPDtr==2) MCnchweight= 0.895842;
else if(multSPDtr==3) MCnchweight= 0.878258;
else if(multSPDtr==4) MCnchweight= 0.889699;
else if(multSPDtr==5) MCnchweight= 0.923275;
else if(multSPDtr==6) MCnchweight= 0.968181;
else if(multSPDtr==7) MCnchweight= 1.00688;
else if(multSPDtr==8) MCnchweight= 1.0429;
else if(multSPDtr==9) MCnchweight= 1.06309;
else if(multSPDtr==10) MCnchweight= 1.08015;
else if(multSPDtr==11) MCnchweight= 1.08931;
else if(multSPDtr==12) MCnchweight= 1.10245;
else if(multSPDtr==13) MCnchweight= 1.10636;
else if(multSPDtr==14) MCnchweight= 1.10981;
else if(multSPDtr==15) MCnchweight= 1.10993;
else if(multSPDtr==16) MCnchweight= 1.10245;
else if(multSPDtr==17) MCnchweight= 1.09729;
else if(multSPDtr==18) MCnchweight= 1.09078;
else if(multSPDtr==19) MCnchweight= 1.07516;
else if(multSPDtr==20) MCnchweight= 1.06905;
else if(multSPDtr==21) MCnchweight= 1.05377;
else if(multSPDtr==22) MCnchweight= 1.04431;
else if(multSPDtr==23) MCnchweight= 1.03231;
else if(multSPDtr==24) MCnchweight= 1.02218;
else if(multSPDtr==25) MCnchweight= 1.01132;
else if(multSPDtr==26) MCnchweight= 0.99894;
else if(multSPDtr==27) MCnchweight= 0.994848;
else if(multSPDtr==28) MCnchweight= 0.98253;
else if(multSPDtr==29) MCnchweight= 0.973244;
else if(multSPDtr==30) MCnchweight= 0.968428;
else if(multSPDtr==31) MCnchweight= 0.960728;
else if(multSPDtr==32) MCnchweight= 0.953542;
else if(multSPDtr==33) MCnchweight= 0.95198;
else if(multSPDtr==34) MCnchweight= 0.955185;
else if(multSPDtr==35) MCnchweight= 0.952798;
else if(multSPDtr==36) MCnchweight= 0.956194;
else if(multSPDtr==37) MCnchweight= 0.966226;
else if(multSPDtr==38) MCnchweight= 0.967056;
else if(multSPDtr==39) MCnchweight= 0.971078;
else if(multSPDtr==40) MCnchweight= 0.982089;
return MCnchweight;



}
//=================================================================================================================================
*/
//====================================================================================================================================
Int_t AliAnalysisTaskHFEmultTPCEMCAL::GetElecSourceType(AliAODMCParticle *electron,Double_t &ptm)
{
    //
    // Return what type of gammax it is
    //
    
    // Mother
    Int_t motherlabel = electron->GetMother();
    if(motherlabel<0) return kNoMotherE;
    else {
        
        AliAODMCParticle *mother = (AliAODMCParticle*)fMCArray->At(motherlabel);
        Int_t motherpdg = TMath::Abs(mother->GetPdgCode());
        ptm=mother->Pt();
        if(motherpdg == 111) {
            Int_t typepi0eta = GetPi0EtaType(mother);
            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) 	return kPi0NoFeedDown;
        }
        if(motherpdg == 221) {
            Int_t typepi0eta = GetPi0EtaType(mother);
            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kEtaNoFeedDown;
        }
        if(motherpdg == 22) {
            Int_t gmotherlabel = mother->GetMother();
            if(gmotherlabel<0) return kDirectGamma;
            else {
                AliAODMCParticle *gmother = (AliAODMCParticle*)fMCArray->At(gmotherlabel);
                ptm=gmother->Pt();
                Int_t gmotherpdg = TMath::Abs(gmother->GetPdgCode());
                if(gmotherpdg == 111) {
                    Int_t typepi0eta = GetPi0EtaType(mother);
                    if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGPi0NoFeedDown;
                }
                if(gmotherpdg == 221) {
                    Int_t typepi0eta = GetPi0EtaType(mother);
                    if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGEtaNoFeedDown;
                }
                if(gmotherpdg == 22) {
                    Int_t ggmotherlabel = gmother->GetMother();
                    if(ggmotherlabel<0) return kDirectGamma;
                    else {
                        AliAODMCParticle *ggmother = (AliAODMCParticle*)fMCArray->At(ggmotherlabel);
                        ptm=ggmother->Pt();
                        Int_t ggmotherpdg = TMath::Abs(ggmother->GetPdgCode());
                        if(ggmotherpdg == 111) {
                            Int_t typepi0eta = GetPi0EtaType(mother);
                            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGPi0NoFeedDown;
                        }
                        if(ggmotherpdg == 221) {
                            Int_t typepi0eta = GetPi0EtaType(mother);
                            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGEtaNoFeedDown;
                        }
                    }
                }
            }
        }
    }
    
    return kOthersE;

}

Int_t AliAnalysisTaskHFEmultTPCEMCAL::GetPi0EtaType(AliAODMCParticle *pi0eta)
{
   
    // Return what type of pi0, eta it is
    
    // IsPrimary
    Bool_t primMC = pi0eta->IsPrimary();
    if(!primMC) return kNoIsPrimary;
    
    // Mother
    Int_t motherlabel = pi0eta->GetMother();
    if(motherlabel<0) return kNoMother;
    else {
        
        AliAODMCParticle *mother = (AliAODMCParticle*)fMCArray->At(motherlabel);
        Int_t motherpdg = TMath::Abs(mother->GetPdgCode());
        //    if(motherpdg == 111 || motherpdg == 221 || motherpdg == 223 || motherpdg == 333 || motherpdg == 331 || motherpdg == 113) return kLightMesons;
        if(motherpdg == 111 || motherpdg == 221 || motherpdg == 223 || motherpdg == 333 || motherpdg == 331 || motherpdg == 113 || motherpdg == 213 || motherpdg == 313 || motherpdg == 323) return kLightMesons;
        
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty;
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm;
        return kNoFeedDown;
        
    }
}


Int_t AliAnalysisTaskHFEmultTPCEMCAL::GetHFE(AliAODMCParticle *electron, TClonesArray *mcArray)
{
	Int_t motherindex=electron->GetMother(); //Getting Electron Mother
	if(motherindex<0) return kNoMother;
	AliAODMCParticle *mother = (AliAODMCParticle*)mcArray->At(motherindex);					
	Int_t motherpdg = mother->GetPdgCode();	
	
	if ( (int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty;
	if ( (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm;
	return kOthersE;
  
}
/*
//====================================================================================================================================


Int_t AliAnalysisTaskHFEmultTPCEMCAL::GetNcharged(){
    //counts all tracks in eta<1 with charge!=0
	
    Int_t Nch = 0;

		
    if(!fIsMC) return Nch; // if no MC info return 0

    // loop over all tracks 
    for (Int_t igen = 0; igen < fMCArray->GetEntriesFast(); igen++){
        AliAODMCParticle *mctrack=(AliAODMCParticle*)fMCArray->UncheckedAt(igen);
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
*/
//=======================================================================


Bool_t AliAnalysisTaskHFEmultTPCEMCAL::GetNMCPartProduced()
{
    //Get number of MC particles produced by generators.
    
    TList *lh = fMCHeader->GetCocktailHeaders();
    fNTotMCpart = 0;
    fNembMCpi0 = 0;
    fNembMCeta = 0;
    fNpureMC = 0;
    TString MCgen;
    TString embpi0("pi");
    TString embeta("eta");
    
    if(!lh){
        AliError("no MC header");
        return (0);
    }
    
    for(int igene=0; igene<lh->GetEntries(); igene++)
    {
        AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igene);
        if(!gh) continue;
        
        MCgen =  gh->GetName();
        //cout << "Gen name, N produced = " << gh->GetName() << ", " << gh->NProduced() << endl;
        if(igene==0) fNpureMC = gh->NProduced();  // generated by MB
        
        //   if(MCgen.Contains(embpi0))cout << MCgen << endl;
        //   if(MCgen.Contains(embeta))cout << MCgen << endl;
        
        if(MCgen.Contains(embpi0))fNembMCpi0 = fNTotMCpart;
        if(MCgen.Contains(embeta))fNembMCeta = fNTotMCpart;
        fNTotMCpart += gh->NProduced();
    }
    //cout << "fNpureMC, fNembMCpi0, fNembMCeta, fNTotMCpart : " <<fNpureMC << ", " << fNembMCpi0 << ", " << fNembMCeta << ", " << fNTotMCpart << endl;
    //getchar();
    return kTRUE;
}
void AliAnalysisTaskHFEmultTPCEMCAL::GetPi0EtaWeight(THnSparse *SparseWeight)
{
    //Get pi0 and eta information for weight calculation
    
    Double_t fvalue[5] = {-999,-999,-999,-999,-999};
    
    for(int imc=0; imc< fNTotMCpart; imc++)
    {
        AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCArray->At(imc);
        if(TMath::Abs(AODMCtrack->Eta()) > 0.9) continue;
        
        //-------Get PDG
        Int_t TrackPDG = TMath::Abs(AODMCtrack->GetPdgCode());
        if((TrackPDG != 111) && (TrackPDG != 221) && (TrackPDG != 22)) continue;
        
        Double_t fPartPDGid = -999;
        if (TrackPDG == 111) fPartPDGid = 0.2;
        if (TrackPDG == 221) fPartPDGid = 1.2;
        if (TrackPDG == 22) fPartPDGid = 2.2;
        
        Double_t fTrkPt = AODMCtrack->Pt();
        
        //-------Check if the particle is from Enhanced signal or not
        Bool_t fFromEnhance = kMB;
        if(imc >= fNpureMC)fFromEnhance = kEnhance;
        
        //------Get type of the particle
        Int_t fType = GetPi0EtaType(AODMCtrack);
        
        fvalue[0] = fTrkPt;
        fvalue[1] = fPartPDGid;
        fvalue[2] = fFromEnhance;
        fvalue[3] = fType;
        fvalue[4] =SPDntr;
        
        SparseWeight->Fill(fvalue);
    }
}




//_________________________________________
Bool_t AliAnalysisTaskHFEmultTPCEMCAL::GetNonHFEEffiDenom(AliVTrack *track)
//Bool_t AliAnalysisTaskHFEmultTPCEMCAL::GetNonHFEEffiDenom(AliAODTrack *track)
{
    //Calculate Non-HFE efficiency demoninator
    
    fIsFrmEmbPi0 = kFALSE, fIsFrmEmbEta = kFALSE;
    ftype = -1, fWeightPi0 = 1.0, fWeightEta = 1.0, fWeight=1.0;
    Bool_t fFromMB = kTRUE;
    
    Int_t MomPDG = -999, GMomPDG=-999, GGMomPDG=-999, GGGMomPDG=-999;
    Int_t iMCmom = -999, iMCgmom = -999, iMCggmom = -999, iMCgggmom = -999;
    Double_t MomPt =-999.0;
    
    AliAODMCParticle *MCPart = 0;
    AliAODMCParticle *MCPartMom = 0;
    AliAODMCParticle *MCPartGMom = 0;
    AliAODMCParticle *MCPartGGMom = 0;
    AliAODMCParticle *MCPartGGGMom = 0;
    
    Double_t TrkPt = track->Pt();
    Int_t iTrklabel = TMath::Abs(track->GetLabel());
    if(iTrklabel == 0) return kFALSE;
    
    MCPart = (AliAODMCParticle*)fMCArray->At(iTrklabel);
    if(TMath::Abs(MCPart->GetPdgCode())!=11) return kFALSE;
    fRealInclsElecPt->Fill(TrkPt);
    
    Bool_t fNonHFE = IsNonHFE(MCPart, fFromMB, ftype, iMCmom, MomPDG, MomPt);
    if(!fNonHFE) return kFALSE;
    fNonHFeTrkPt->Fill(TrkPt);
    
    MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
    iMCgmom = MCPartMom->GetMother();
    if(iMCgmom > 0){
        MCPartGMom = (AliAODMCParticle*)fMCArray->At(iMCgmom);
        GMomPDG = TMath::Abs(MCPartGMom->GetPdgCode());
        
        iMCggmom = MCPartGMom->GetMother();
        if(iMCggmom > 0){
            MCPartGGMom = (AliAODMCParticle*)fMCArray->At(iMCggmom);
            GGMomPDG = TMath::Abs(MCPartGGMom->GetPdgCode());
            
            iMCgggmom = MCPartGGMom->GetMother();
            if(iMCgggmom > 0){
                MCPartGGGMom = (AliAODMCParticle*)fMCArray->At(iMCgggmom);
                GGGMomPDG = TMath::Abs(MCPartGGGMom->GetPdgCode());
            }
        }
    }
    
    //cases to consider: eta->e, eta->pi0->e, eta->gamma->e, eta->pi0->gamma->e, pi0->e, pi0->gamma->e
    if(MomPDG == 221){
        if(iMCmom >= fNembMCeta && iMCmom < fNTotMCpart) { //from eta event
            fIsFrmEmbEta = kTRUE; //eta->e
            fWeightEta = fEtaWeight->Eval(MCPartMom->Pt());
        }
    }
    
    if(MomPDG == 111) {
        if(iMCmom >= fNembMCpi0 && iMCmom < fNembMCeta){ //from pi0 event
            fIsFrmEmbPi0 = kTRUE; //pi0 -> e
            fWeightPi0 = fPi0Weight->Eval(MCPartMom->Pt());
        }
        
        if(GMomPDG == 221){
            if(iMCgmom >= fNembMCeta && iMCgmom < fNTotMCpart) { //from eta event
                fIsFrmEmbEta = kTRUE; //eta->pi0-> e
                fWeightEta = fEtaWeight->Eval(MCPartGMom->Pt());
            }
        }
    }
    
    if(MomPDG == 22){
        if(GMomPDG == 221){
            if(iMCgmom >= fNembMCeta && iMCgmom < fNTotMCpart) { //from eta event
                fIsFrmEmbEta = kTRUE; //eta->gamma-> e
                fWeightEta = fEtaWeight->Eval(MCPartGMom->Pt());
            }
        }
        
        if(GMomPDG == 111){
            if(iMCgmom >= fNembMCpi0 && iMCgmom < fNembMCeta) { //from pi0 event
                fIsFrmEmbPi0 = kTRUE; //pi0-> gamma-> e
                fWeightPi0 = fPi0Weight->Eval(MCPartGMom->Pt());
            }
            
            if(GGMomPDG == 221){
                if(iMCggmom >= fNembMCeta && iMCggmom < fNTotMCpart) { //from eta event
                    fIsFrmEmbEta = kTRUE; //eta->pi0->gamma-> e
                    fWeightEta = fEtaWeight->Eval(MCPartGGMom->Pt());
                }
            }
        }
    }
    
    //   cout << "PDG of M, GM, GGM, GGGM of ele: "<< MomPDG << ", " << GMomPDG << ", " << GGMomPDG << ", " << GGGMomPDG << endl;
    //   cout << "==============" <<endl;
    
    if(fIsFrmEmbPi0 || fIsFrmEmbEta){
    		//cout<<" Denom "<<SPDntr<<endl;
    		//getchar();
        fNonHFeEmbTrkPt->Fill(TrkPt);
        fNonHFeEmbTrkPt_SPD->Fill(TrkPt,SPDntr);
       
        if(fIsFrmEmbPi0) {
            fWeight = fWeightPi0;
            fPi0eEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
            fNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
            fNonHFeEmbWeightTrkPt_SPD->Fill(TrkPt,SPDntr,fWeightPi0);
        }
        if(fIsFrmEmbEta){
            fWeight = fWeightEta;
            fEtaeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
            fNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
            fNonHFeEmbWeightTrkPt_SPD->Fill(TrkPt,SPDntr,fWeightEta);
        }
    }
    
    return kTRUE;
}

//=======================================================================
Bool_t  AliAnalysisTaskHFEmultTPCEMCAL::IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromMB, Int_t &type, Int_t &iMCmom, Int_t &MomPDG, Double_t &MomPt)
{
    //Is electron from pi0, eta and gamma
    
    iMCmom = MCPart->GetMother();
    AliAODMCParticle *MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
    MomPDG = TMath::Abs(MCPartMom->GetPdgCode());
    MomPt = MCPartMom->Pt();
    
    if((MomPDG == 111) || (MomPDG == 221) || (MomPDG == 22)){
        if(iMCmom >= fNpureMC)fFromMB = kFALSE;
        type = GetPi0EtaType(MCPartMom);
        //cout<<"IsNonHFE--> iMCmom "<<iMCmom<<"  fNpureMC   "<<fNpureMC<<endl;
        return kTRUE;
    }
    else return kFALSE;
}

Bool_t AliAnalysisTaskHFEmultTPCEMCAL::GetNonHFEEffiRecoTag(AliVTrack *track)
{
    //Tagging method

    Double_t TrkPt = track->Pt();
    
    fRecoNonHFeTrkPt->Fill(TrkPt);
    fRecoNonHFeTrkPt_SPD->Fill(TrkPt,SPDntr);
    
    if(fIsFrmEmbPi0 || fIsFrmEmbEta){
        fRecoNonHFeEmbTrkPt->Fill(TrkPt);
        fRecoNonHFeEmbTrkPt_SPD->Fill(TrkPt,SPDntr);
        
        if(fIsFrmEmbPi0) {
            fRecoPi0eEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
            fRecoNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
            fRecoNonHFeEmbWeightTrkPt_SPD->Fill(TrkPt,SPDntr,fWeightPi0);
        }
        if(fIsFrmEmbEta){
            fRecoEtaeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
            fRecoNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
            fRecoNonHFeEmbWeightTrkPt_SPD->Fill(TrkPt,SPDntr,fWeightEta);
        }
    }
        
    return kTRUE;
}

//_________________________________________
Bool_t AliAnalysisTaskHFEmultTPCEMCAL::GetNonHFEEffiULSLS(AliVTrack *track, AliVTrack *Assotrack, Bool_t fFlagLS, Bool_t fFlagULS, Double_t mass)
{
    //ULS-LS method

    Double_t TrkPt = track->Pt();
    
    //Track information
    Int_t iTrklabel = TMath::Abs(track->GetLabel());
    if(iTrklabel == 0) return kFALSE;
    AliAODMCParticle *MCPart = (AliAODMCParticle*)fMCArray->At(iTrklabel);
    
    if(TMath::Abs(MCPart->GetPdgCode())!=11) return kFALSE;
    Bool_t fFromMB = kTRUE;
    Int_t iMCmom=-999, MomPDG = -999, type=-1;
    Double_t MomPt =-999;
    Bool_t fNonHFE = IsNonHFE(MCPart, fFromMB, type, iMCmom, MomPDG, MomPt);
    
    //Associated partner information
    Int_t iTrkAssolabel = TMath::Abs(Assotrack->GetLabel());
    if(iTrkAssolabel == 0) return kFALSE;
    AliAODMCParticle *MCPartAsso = (AliAODMCParticle*)fMCArray->At(iTrkAssolabel);
    
    if(TMath::Abs(MCPartAsso->GetPdgCode())!=11) return kFALSE; // check origin of asso elec
    
    Bool_t fAssoFromMB = kTRUE;
    Int_t iMCAssomom=-999, AssoMomPDG = -999, fAssotype=-1;
    Double_t AssoMomPt =-999;
    Bool_t fAssoNonHFE = IsNonHFE(MCPartAsso, fAssoFromMB, fAssotype, iMCAssomom, AssoMomPDG, AssoMomPt);
    
    //cout << "Asso ele mom : " << iMCAssomom << ", " << AssoMomPDG << ", " << iMCmom << ", " << MomPDG << ", " << fIsFrmEmbPi0 << ", " << fIsFrmEmbEta << ", " << type << endl;
    
    if(!fAssoNonHFE) return kFALSE;
    if(iMCmom != iMCAssomom) return kFALSE; //ensure electron and partner comes from same mother
    
    if(fFlagLS) fNonHFePairInvmassLS->Fill(mass);
    if(fFlagULS) fNonHFePairInvmassULS->Fill(mass);
    
    if((fIsFrmEmbPi0 || fIsFrmEmbEta) && ftype==kNoMother){ //If parent e from embedded pi0/eta + NoMom
        if(fFlagLS) fNonHFeEmbInvmassLS->Fill(mass);
        if(fFlagULS) fNonHFeEmbInvmassULS->Fill(mass);
        if(fFlagLS) fNonHFeEmbWeightInvmassLS->Fill(mass, fWeight);
        if(fFlagULS) fNonHFeEmbWeightInvmassULS->Fill(mass, fWeight);
        
        if(fIsFrmEmbPi0){ //if from pi0
            if(fFlagLS) fPi0EmbInvmassLS->Fill(mass);
            if(fFlagULS) fPi0EmbInvmassULS->Fill(mass);
            if(fFlagLS) fPi0EmbWeightInvmassLS->Fill(mass, fWeight);
            if(fFlagULS) fPi0EmbWeightInvmassULS->Fill(mass, fWeight);
        }
        if(fIsFrmEmbEta){ //if from eta
            if(fFlagLS) fEtaEmbInvmassLS->Fill(mass);
            if(fFlagULS) fEtaEmbInvmassULS->Fill(mass);
            if(fFlagLS) fEtaEmbWeightInvmassLS->Fill(mass, fWeight);
            if(fFlagULS) fEtaEmbWeightInvmassULS->Fill(mass, fWeight);
        }
    }
    
    if(mass < fInvmassCut){
        if(fFlagLS){
            //new method
            if(fIsFrmEmbPi0 || fIsFrmEmbEta) {
                fRecoLSeEmbTrkPt->Fill(TrkPt);
                
                if(fIsFrmEmbPi0) {
                    fRecoLSeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
                    fRecoPi0LSeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
                }
                if(fIsFrmEmbEta){
                    fRecoLSeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
                    fRecoEtaLSeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
                }
            }
            
        }
        
        if(fFlagULS){
            //new method
            if(fIsFrmEmbPi0 || fIsFrmEmbEta) {
                fRecoULSeEmbTrkPt->Fill(TrkPt);
                
                if(fIsFrmEmbPi0) {
                    fRecoULSeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
                    fRecoPi0ULSeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
                }
                
                if(fIsFrmEmbEta){
                    fRecoULSeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
                    fRecoEtaULSeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
                }
            }
            
        }
    }
    
    return kTRUE;
}




//===============================================================================
TProfile* AliAnalysisTaskHFEmultTPCEMCAL::GetEstimatorHistogram(const AliAODEvent* fAOD)
{

  if (fPeriod < 0 || fPeriod > 5) return 0;   
  return fMultEstimatorAvg[fPeriod];
}
/*
TProfile* AliAnalysisTaskHFEmultTPCEMCAL::GetEstimatorHistogram(const AliAODEvent* fAOD)
{
    
  Int_t runNo  = fAOD->GetRunNumber();
 //cout<<"run number"<<runNo<<endl;
  Int_t period = -1; 
   
        
  if (runNo>=258919 && runNo<=259888){
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
*/
Double_t AliAnalysisTaskHFEmultTPCEMCAL::GetCorrectedNtracklets(TProfile* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult) {
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
/*
//===================================================================================
Double_t AliAnalysisTaskHFEmultTPCEMCAL::Beta(AliAODTrack *track)
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

*/
//===================================================================================
void AliAnalysisTaskHFEmultTPCEMCAL::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  
	fOutputList = dynamic_cast<TList*> (GetOutputData(1));
	if (!fOutputList)return;

}


