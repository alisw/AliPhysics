
// Task to estimate the number of gamma-hadron
// statistic available in the Pb+Pb run.
//
// Authors: E. Epple, M. Oliver, based on code by  B. Sahlmueller and C. Loizides

#include <Riostream.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include "TCustomBinning.h"

#include "AliAnalysisTaskGammaHadron.h"

#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliVEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliExternalTrackParam.h"
#include "AliTrackerBase.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"
#include "AliPicoTrack.h"
#include "AliVVZERO.h"
#include "AliESDUtils.h"
#include "AliEventPoolManager.h"
#include "AliMCEvent.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskGammaHadron)
//
//  This class inherits from AliAnalysisTaskEmcal() -->fcent,fVertex,fNVertCont,fVertexSPD,fNVertSPDCont,fTriggers is defined in there
//  And AliAnalysisTaskEmcal() inherits from AliAnalysisTaskSE()
//________________________________________________________________________
AliAnalysisTaskGammaHadron::AliAnalysisTaskGammaHadron()
  : AliAnalysisTaskEmcal("AliAnalysisTaskGammaHadron", kTRUE),
  fEventCuts(0),fFiducialCuts(0x0),fFiducialCellCut(0x0),fFlowQnVectorMgr(0x0),
  fGammaOrPi0(0),fSEvMEv(0),fSaveTriggerPool(0),fDownScaleMT(1.0),fSidebandChoice(0),
  fDebug(0),fSavePool(0),fPlotQA(0),
  fUseManualEventCuts(0),fCorrectEff(0),fEventWeightChoice(0),
  fRtoD(0),fSubDetector(0),
  fTriggerPtCut(5.),fMaxPi0Pt(23.),fClShapeMin(0),fClShapeMax(10),fClEnergyMin(2),fOpeningAngleCut(0.017),fMaxNLM(10),
  fRmvMTrack(0),fClusEnergyType(0),fHadCorr(0),fHadCorrConstant(0.236),fTrackMatchEta(0),fTrackMatchPhi(0),fTrackMatchEOverPLow(0.6),fTrackMatchEOverPHigh(1.4),
  fMixBCent(0),fMixBZvtx(0),fMixBEMCalMult(0),fMixBClusZvtx(0),
  fPoolMgr(0x0),fTrackDepth(0),fTargetFraction(0.1),fClusterDepth(0),fPoolSize(0),fEventPoolOutputList(0),
  fTriggerType(AliVEvent::kINT7),fPi0MassSelection(3), fMixingEventType(AliVEvent::kINT7),fCurrentEventTrigger(0),fVetoTrigger(AliVEvent::kEMCEGA),
  fApplyPatchCandCut(0),
  fQnCorrEventPlaneAngle(0.0),fQnCorrEventPlane3Angle(0.0),fQnCorrEventPlane4Angle(0.0),
  fParticleLevel(kFALSE),fIsMC(0),fMCEmbedReweightMode(0),fUseMCReactionPlane(0),fMCHeader(0),fMCParticles(0),fMCPi0List(0),fMCReactionPlaneAngle(0),
  fEventCutList(0),fOutputListQA(0),
  fEPAngleV0M(0),fEPAngleTPCA(0),fEPAngleTPCC(0),
  fEP3AngleV0M(0),fEP3AngleTPCA(0),fEP3AngleTPCC(0),
  fEP4AngleV0M(0),fEP4AngleTPCA(0),fEP4AngleTPCC(0),
  fEPR_CosD1(0),fEPR_CosD2(0),fEPR_CosD3(0),
  fEP3R_CosD1(0),fEP3R_CosD2(0),fEP3R_CosD3(0),
  fEP4R_CosD1(0),fEP4R_CosD2(0),fEP4R_CosD3(0),
  fHistMCPi0_PtEtaMult(0),fHistMCPi0_PtEtaEP(0),fEtaPhiMCPion(0),
  fHistClusPairInvarMasspT(0),fHistPi0(0),fMAngle(0),fPtAngle(0),fMassPionRej(0),
  fPtEPAnglePionAcc(0),fPtEPAnglePionAccCent(0),fPtEPAngleMCPion(0),fPtEPAngleTrueRecMCPion(0),
  fPtEP3AnglePionAcc(0),fPtEP3AnglePionAccCent(0),fPtEP3AngleMCPion(0),fPtEP3AngleTrueRecMCPion(0),
  fPtEP4AnglePionAcc(0),fPtEP4AnglePionAccCent(0),fPtEP4AngleMCPion(0),fPtEP4AngleTrueRecMCPion(0),
  fHistTrackPsiEPPtCent(0),fHistTrackPsiEP3PtCent(0),fHistTrackPsiEP4PtCent(0),fMCReactionPlane(0),fPtRPAnglePionAcc(0),fPtRPAngleMCPion(0),fPtRPAngleTrueRecMCPion(0),fHistTrackPsiRPPtCent(0),
  fEtaPhiPionAcc(0),fMassPtPionAcc(0),fMassPtPionRej(0),fMassPtCentPionAcc(0),fMassPtCentPionRej(0),
  fMatchDeltaEtaTrackPt(0),fMatchDeltaPhiTrackPt(0),fMatchCondDeltaEtaTrackPt(0),fMatchCondDeltaPhiTrackPt(0),fClusterEnergyMatchedTracks(0),fHistEOverPvE(0),fHistPOverEvE(0),
  fHistPSDistU(0),fHistPSDistV(0),
  fRand(0),
  fClusEnergy(0),fAccClusEtaPhi(0),fAccClusEtaPhiZvtx(0),bEnableClusPairRot(0),fDoRotBkg(0),fDoClusMixing(0),fDoPosSwapMixing(0),fNRotBkgSamples(1),fPi0Cands(0),fHistEventHash(0),
  bEnablePosSwapHists(false),bLogPSMod(true),fPSMassPtMap(0),fESMassPtMap(0),fUScaleMatrix(0),fVScaleMatrix(0),
  fEMCalMultvZvtx(0),
  fHistClusMCDE(0),fHistClusMCDPhiDEta(0),fHistPi0MCDPt(0),fHistEtaMCDPt(0),fHistPi0MCDPhiDEta(0),fHistEtaMCDPhiDEta(0),
  fUseParamMassSigma(0),fPi0NSigma(2.),fPi0AsymCut(1.0),
  fEffCorrectionCheck(0),
  fHistEvsPt(0),fHistBinCheckPt(0),fHistBinCheckZt(0),fHistBinCheckXi(0), fHistBinCheckEvtPl(0), fHistBinCheckEvtPl2(0),
  fHistDEtaDPhiGammaQA(0),fHistDEtaDPhiTrackQA(0), fHistClusterTime(0),
  fCorrVsManyThings(0),fTriggerHist(0),fClusterProp(0x0),
  fHPoolReady(0x0)
{
	//..Initialize by defult for
	//..AliAnalysisTaskGammaHadron(0,0);
	InitArrays();
}
//________________________________________________________________________
AliAnalysisTaskGammaHadron::AliAnalysisTaskGammaHadron(Int_t InputGammaOrPi0,Int_t InputSeMe,Bool_t InputMCorData)
  : AliAnalysisTaskEmcal("AliAnalysisTaskGammaHadron", kTRUE),
  fEventCuts(0),fFiducialCuts(0x0),fFiducialCellCut(0x0),fFlowQnVectorMgr(0x0),
  fGammaOrPi0(0),fSEvMEv(0),fSaveTriggerPool(0),fDownScaleMT(1.0),fSidebandChoice(0),
  fDebug(0),fSavePool(0),fPlotQA(0),
  fUseManualEventCuts(0),fCorrectEff(0),fEventWeightChoice(0),
  fRtoD(0),fSubDetector(0),
  fTriggerPtCut(5.),fMaxPi0Pt(23.),fClShapeMin(0),fClShapeMax(10),fClEnergyMin(2),fOpeningAngleCut(0.017),fMaxNLM(10),
  fRmvMTrack(0),fClusEnergyType(0),fHadCorr(0),fHadCorrConstant(0.236),fTrackMatchEta(0),fTrackMatchPhi(0),fTrackMatchEOverPLow(0.6),fTrackMatchEOverPHigh(1.4),
  fMixBCent(0),fMixBZvtx(0),fMixBEMCalMult(0),fMixBClusZvtx(0),
  fPoolMgr(0x0),fTrackDepth(0),fTargetFraction(0.1),fClusterDepth(0),fPoolSize(0),fEventPoolOutputList(0),
  fTriggerType(AliVEvent::kINT7),fPi0MassSelection(3), fMixingEventType(AliVEvent::kINT7),fCurrentEventTrigger(0),fVetoTrigger(AliVEvent::kEMCEGA),
  fApplyPatchCandCut(0),
  fQnCorrEventPlaneAngle(0.0),fQnCorrEventPlane3Angle(0.0),fQnCorrEventPlane4Angle(0.0),
  fParticleLevel(kFALSE),fIsMC(InputMCorData),fMCEmbedReweightMode(0),fUseMCReactionPlane(0),fMCHeader(0),fMCParticles(0),fMCPi0List(0),fMCReactionPlaneAngle(0),
  fEventCutList(0),fOutputListQA(0),
  fEPAngleV0M(0),fEPAngleTPCA(0),fEPAngleTPCC(0),
  fEP3AngleV0M(0),fEP3AngleTPCA(0),fEP3AngleTPCC(0),
  fEP4AngleV0M(0),fEP4AngleTPCA(0),fEP4AngleTPCC(0),
  fEPR_CosD1(0),fEPR_CosD2(0),fEPR_CosD3(0),
  fEP3R_CosD1(0),fEP3R_CosD2(0),fEP3R_CosD3(0),
  fEP4R_CosD1(0),fEP4R_CosD2(0),fEP4R_CosD3(0),
  fHistMCPi0_PtEtaMult(0),fHistMCPi0_PtEtaEP(0),fEtaPhiMCPion(0),
  fHistClusPairInvarMasspT(0),fHistPi0(0),fMAngle(0),fPtAngle(0),fMassPionRej(0),
  fPtEPAnglePionAcc(0),fPtEPAnglePionAccCent(0),
  fPtEPAngleMCPion(0),fPtEPAngleTrueRecMCPion(0),
  fPtEP3AnglePionAcc(0),fPtEP3AnglePionAccCent(0),
  fPtEP3AngleMCPion(0),fPtEP3AngleTrueRecMCPion(0),
  fPtEP4AnglePionAcc(0),fPtEP4AnglePionAccCent(0),
  fPtEP4AngleMCPion(0),fPtEP4AngleTrueRecMCPion(0),
  fHistTrackPsiEPPtCent(0),fHistTrackPsiEP3PtCent(0),fHistTrackPsiEP4PtCent(0),fMCReactionPlane(0),fPtRPAnglePionAcc(0),fPtRPAngleMCPion(0),fPtRPAngleTrueRecMCPion(0),fHistTrackPsiRPPtCent(0),
  fEtaPhiPionAcc(0),fMassPtPionAcc(0),fMassPtPionRej(0),fMassPtCentPionAcc(0),fMassPtCentPionRej(0),
  fMatchDeltaEtaTrackPt(0),fMatchDeltaPhiTrackPt(0),fMatchCondDeltaEtaTrackPt(0),fMatchCondDeltaPhiTrackPt(0),fClusterEnergyMatchedTracks(0),fHistEOverPvE(0),fHistPOverEvE(0),
  fHistPSDistU(0),fHistPSDistV(0),
  fRand(0),
  fClusEnergy(0),fAccClusEtaPhi(0),fAccClusEtaPhiZvtx(0),bEnableClusPairRot(0),fDoRotBkg(0),fDoClusMixing(0),fDoPosSwapMixing(0),fNRotBkgSamples(1),fPi0Cands(0),fHistEventHash(0),
  bEnablePosSwapHists(false),bLogPSMod(true),fPSMassPtMap(0),fESMassPtMap(0),fUScaleMatrix(0),fVScaleMatrix(0),
  fEMCalMultvZvtx(0),
  fHistClusMCDE(0),fHistClusMCDPhiDEta(0),fHistPi0MCDPt(0),fHistEtaMCDPt(0),fHistPi0MCDPhiDEta(0),fHistEtaMCDPhiDEta(0),
  fUseParamMassSigma(0),fPi0NSigma(2.),fPi0AsymCut(1.0),
  fEffCorrectionCheck(0),
  fHistEvsPt(0),fHistBinCheckPt(0),fHistBinCheckZt(0),fHistBinCheckXi(0), fHistBinCheckEvtPl(0), fHistBinCheckEvtPl2(0),
  fHistDEtaDPhiGammaQA(0),fHistDEtaDPhiTrackQA(0), fHistClusterTime(0),
  fCorrVsManyThings(0),fTriggerHist(0),fClusterProp(0x0),
  fHPoolReady(0x0)
{
	InitArrays();
	//..set input variables
	fGammaOrPi0        =InputGammaOrPi0;
	fSEvMEv            =InputSeMe;
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::InitArrays()
{
	//..Initialize by defult for
	//..AliAnalysisTaskGammaHadron(0,1);

	//..set input variables
	fGammaOrPi0        =0; //= 0 ( Gamma analysis ), 1 (pi0 analysis)
	fSEvMEv            =0; //= 0 (do only same event analyis with correct triggers), =1 (do event mixing)

	fPlotQA            =0; //= 0 do not plot additional QA histograms, >0 plot them
	fDebug             =0; //set only 1 for debugging
	fSavePool          =0; //= 0 do not save the pool by default. Use the set function to do this.
	fUseManualEventCuts=0; //= 0 use automatic setting from AliEventCuts. =1 load manual cuts
	fCorrectEff        =1; //= 0 no efficiency correction
	//..These two items are set in AliAnalysisTaskEmcal::RetrieveEventObjects()
	//fCent, zVertex

	//..Initialize the arrays to 0
	/*for(Int_t i=0; i<kNIdentifier;i++)
	{
		fHistptAssHadronG[i] =0;
		fHistptAssHadronZt[i]=0;
		fHistptAssHadronXi[i]=0;
		fHistptTriggG[i] =0;
		fHistptTriggZt[i]=0;
		fHistptTriggXi[i]=0;


		for(Int_t j=0; j<kNIdentifier;j++)
		{
			if(j<kNoGammaBins+1)fHistDEtaDPhiG[i][j]  = 0;
			if(j<kNoZtBins+1)   fHistDEtaDPhiZT[i][j] = 0;
			if(j<kNoXiBins+1)   fHistDEtaDPhiXI[i][j] = 0;
		}
	}*/

	fRtoD=180.0/TMath::Pi();

	fFiducialCuts    = new AliFiducialCut();
	fFiducialCellCut = new AliEMCALRecoUtils();
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Define bins in which the 2D histograms are plotted
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//Double_t fZtStep =1.0/(7-1.0);  // Bin width for the zT histograms
	Double_t fXiStep =2.5/(8-1.0);  // Bin width for the Xi histograms

	Double_t fArray_G_BinsValue[kNoGammaBins+1] ={5,7,9,11,14,17,20,23,30,60};
	Double_t fArray_ZT_BinsValue[kNoZtBins+1]   ={0.03,0.08,0.16,0.29,0.5,0.84,1.39,2.};
//	Double_t fArray_ZT_BinsValue[kNoZtBins+1]   ={0,fZtStep,2*fZtStep,3*fZtStep,4*fZtStep,5*fZtStep,6*fZtStep,20};
	Double_t fArray_XI_BinsValue[kNoXiBins+1]   ={-10,0,fXiStep,2*fXiStep,3*fXiStep,4*fXiStep,5*fXiStep,6*fXiStep,10};
	Double_t fArray_HPT_BinsValue[kNoHPtBins+1]   ={0.2,0.4,0.8,1.5,2.5,4,7,11,17};

	memcpy (fArray_G_Bins,  fArray_G_BinsValue,  sizeof (fArray_G_Bins));
	memcpy (fArray_ZT_Bins, fArray_ZT_BinsValue, sizeof (fArray_ZT_Bins));
	memcpy (fArray_XI_Bins, fArray_XI_BinsValue, sizeof (fArray_XI_Bins));
	memcpy (fArray_HPT_Bins, fArray_HPT_BinsValue, sizeof (fArray_HPT_Bins));
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Define vertex and centrality bins for the ME background
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//..if desired one can add a set function to set these values in the add task function
	//Double_t centmix[kNcentBins+1] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0};
	Double_t centmix[kNcentBins+1] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 90.0};
	fMixBCent = new TAxis(kNcentBins,centmix);

	//static const Int_t NvertBins=8;
	//Double_t zvtxmix[kNvertBins+1] = {-10,-6,-4,-2,0,2,4,6,10};
	//Double_t zvtxmix[kNvertBins+1] = {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10};
// Good:
	Double_t zvtxmix[kNvertBins+1] = {-10,-8,-6,-4,-2,0,2,4,6,8,10};
	memcpy (fArrayNVertBins, zvtxmix, sizeof (fArrayNVertBins));
	fMixBZvtx = new TAxis(kNvertBins,zvtxmix);

	// Additional event plane and pt bins for mixed events (only should be used for mixed trigger mode)


	//..Raymond/Megan gives more mixed event yield - don't know about the quality though
	//fTrackDepth     = 100;      //Hanseul sets it to 100! Q:: is this good? Maximum number of tracks??
	fTrackDepth     = 50000;    //Raymonds/Megans value

	fClusterDepth   = 10000;

	//..!!
	//.. fPoolSize is an input that is ignored in the PoolManager Anyway
	//fPoolSize       = 1;     //1000 - Raymond/Megan value, says it is ignored anyway
	fPoolSize       = -1; // fPoolSize is no longer ignored ?  Must be -1 or the max number of events to mix in each pool

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Define Vertex and EMCal multiplicity bins for the Mixed Cluster Pion background
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Double_t emcalMix[kNEMCalMultBins+1] = {0.0, 100, 200, 300, 1000};
	fMixBEMCalMult = new TAxis(kNEMCalMultBins,emcalMix);

	Double_t zClusvtxmix[kNClusVertBins+1] = {-10,-6.,-3.,-1.,1.,3.,6.,10.};
///	memcpy (fArrayNVertBins, zvtxmix, sizeof (fArrayNVertBins));
	fMixBClusZvtx = new TAxis(kNClusVertBins,zClusvtxmix);

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//..Efficiency correction function From March 2018
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	for(Int_t i=0;i<4;i++)
	{
		funcpT_low[i]    = new TF1(Form("pT_low_%i",i),"([0]+[1]*(-1)/x)+[2]*(TMath::Gaus(x,[3],[4],0))",0.15,3.6); //3=mean, 4=width
		funcpT_high[i]   = new TF1(Form("pT_high_%i",i), "([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x)",3.4, 30.); //0.5,10);
		funcpEta_left[i] = new TF1(Form("funcpEtaLeft_%i",i),"([0]*exp(-pow([1]/TMath::Abs(x+0.91),[2]))+[3]*x+[4]*(TMath::Gaus(x,-0.04,[5],0)))",-0.9,0);
		funcpEta_right[i]= new TF1(Form("funcpEtaRight_%i",i),"([0]*exp(-pow([1]/TMath::Abs(-x+0.91),[2]))+[3]*x+[4]*(TMath::Gaus(x,-0.04,[5],0)))",-0.06,0.9);
	}

	Double_t paramSetPt_1[10]  = {0.8350, 0.0621, 0.0986, 0.2000, 1.0124, 0.7568, 0.0277, -0.0034, 0.1506*0.001, -0.0023*0.001 };
	Double_t paramSetPt_2[10]  = {0.8213, 0.0527, 0.0867, 0.1970, 1.1518, 0.7469, 0.0300, -0.0038, 0.1704*0.001, -0.0026*0.001 };
	Double_t paramSetPt_3[10]  = {0.8381, 0.0648, 0.1052, 0.1478, 1.0320, 0.7628, 0.0263, -0.0032, 0.1443*0.001, -0.0023*0.001 };
	Double_t paramSetPt_4[10]  = {0.8437, 0.0668, 0.1083, 0.2000, 0.9741, 0.7677, 0.0255, -0.0030, 0.1260*0.001, -0.0019*0.001 };
	Double_t paramSetEta_1[12] = {1.0086, 0.0074, 0.2404, -0.1230, -0.0107, 0.0427, 0.8579, 0.0088, 0.4697, 0.0772, -0.0352, 0.0645 };
	Double_t paramSetEta_2[12] = {0.9726, 0.0066, 0.2543, -0.1167, -0.0113, 0.0400, 0.8729, 0.0122, 0.4537, 0.0965, -0.0328, 0.0623 };
	Double_t paramSetEta_3[12] = {0.9076, 0.0065, 0.3216, -0.1130, -0.0107, 0.0456, 0.8521, 0.0073, 0.4764, 0.0668, -0.0363, 0.0668 };
	Double_t paramSetEta_4[12] = {1.1259, 0.0105, 0.1961, -0.1330, -0.0103, 0.0440, 0.8421, 0.0066, 0.5061, 0.0580, -0.0379, 0.0651 };

	for(Int_t param=0;param<6;param++)
	{
		if(param<5)
		{
			funcpT_low[0]->SetParameter(param,paramSetPt_1[param]);
			funcpT_high[0]->SetParameter(param,paramSetPt_1[param+5]);
			funcpT_low[1]->SetParameter(param,paramSetPt_2[param]);
			funcpT_high[1]->SetParameter(param,paramSetPt_2[param+5]);
			funcpT_low[2]->SetParameter(param,paramSetPt_3[param]);
			funcpT_high[2]->SetParameter(param,paramSetPt_3[param+5]);
			funcpT_low[3]->SetParameter(param,paramSetPt_4[param]);
			funcpT_high[3]->SetParameter(param,paramSetPt_4[param+5]);
		}
		funcpEta_left[0]->SetParameter(param,paramSetEta_1[param]);
		funcpEta_right[0]->SetParameter(param,paramSetEta_1[param+6]);
		funcpEta_left[1]->SetParameter(param,paramSetEta_2[param]);
		funcpEta_right[1]->SetParameter(param,paramSetEta_2[param+6]);
		funcpEta_left[2]->SetParameter(param,paramSetEta_3[param]);
		funcpEta_right[2]->SetParameter(param,paramSetEta_3[param+6]);
		funcpEta_left[3]->SetParameter(param,paramSetEta_4[param]);
		funcpEta_right[3]->SetParameter(param,paramSetEta_4[param+6]);
	}

  SetPi0MassSelection(fPi0MassSelection);

	// Pi0 Mass and Sigma Fit parameters (for mass window)
	Double_t fPi0MassFitParsValue[5] = {10.49,0.13852,-1.17e-4,2.861e-3,0};
	memcpy (fPi0MassFitPars, fPi0MassFitParsValue, sizeof(fPi0MassFitPars));
	Double_t fPi0SigmaFitParsValue[5] = {8.34,9.90e-3,-1.09e-4,6.86e-4,0};
	memcpy (fPi0SigmaFitPars, fPi0SigmaFitParsValue, sizeof(fPi0SigmaFitPars));

	//..member function of AliAnalysisTaskEmcal
	SetMakeGeneralHistograms(kTRUE);
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::SetPi0MassSelection(Int_t input) {

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//..Pi0 Cut values
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// Bins 5-7,7-9,9-11,11-14,14-17,17-20,20-23,23-30,30-60
	//  Double_t fPi0MassFixedValue[kNoGammaBins] = {0.135,0.135,0.135,
	//                                          0.135,0.135,0.135,
	//                                          0.135,0.135,0.135}; //9
	//  Double_t fPi0SigmaFixedValue[kNoGammaBins] = {0.01,0.01,0.01,
	//                                           0.01,0.01,0.01,
	//                                           0.01,0.01,0.01}; //9

	// These cuts correspond to
	// Lambda Range: [0.10 - 0.40]       Bins: 1 3
	// Energy Range: [2.00 - 119.94]     Bins: 5 6
	// Asym Range:   [0.00 - 0.80]       Bins: 1 4
	//  Double_t fPi0MassFixedValue[kNoGammaBins] = {0.136189, 0.132716, 0.137215,
	//                                               0.144112, 0.155093, 0.167641,
	//                                               0.192909, 0.219976, 0.219976};
	//  Double_t fPi0SigmaFixedValue[kNoGammaBins] = {0.013780,0.016556, 0.015154,
	//                                                0.014779, 0.017486, 0.018040,
	//                                                0.021053, 0.029528, 0.029528}; //9
	// These are the values used prior to 20180824
	// These cuts correspond to
	// Lambda Range: [0.10 - 0.50]       Bins: 1 5
	// Energy Range: [1.50 - 119.94]     Bins: 4 6
	// Asym Range:   [0.00 - 1.0]        Bins: 1 5
	// Angle Min Range:   [0.017 - 3.142]    Bins: 6 23
  // Not split up by cent bins (all cent bins identical
	Double_t fPi0MassFixedValue_0[kNMainCentBins][kNoGammaBins] = {
    { 0.138350, 0.130986, 0.137813, 0.145594, 0.156744, 0.175125, 0.220000, 0.220000, 0.220000},
    { 0.138350, 0.130986, 0.137813, 0.145594, 0.156744, 0.175125, 0.220000, 0.220000, 0.220000},
    { 0.138350, 0.130986, 0.137813, 0.145594, 0.156744, 0.175125, 0.220000, 0.220000, 0.220000},
    { 0.138350, 0.130986, 0.137813, 0.145594, 0.156744, 0.175125, 0.220000, 0.220000, 0.220000}
  };
	Double_t fPi0SigmaFixedValue_0[kNMainCentBins][kNoGammaBins] = {
    { 0.012870, 0.021483, 0.015919, 0.016042, 0.017068, 0.021500, 0.031488, 0.031488, 0.031488},
    { 0.012870, 0.021483, 0.015919, 0.016042, 0.017068, 0.021500, 0.031488, 0.031488, 0.031488},
    { 0.012870, 0.021483, 0.015919, 0.016042, 0.017068, 0.021500, 0.031488, 0.031488, 0.031488},
    { 0.012870, 0.021483, 0.015919, 0.016042, 0.017068, 0.021500, 0.031488, 0.031488, 0.031488}
  }; //9
	
	// These are the values add for GA triggered data on 20180824
  // These were not done in different centrality bins, so all cent bins are identical
	// Lambda Range: [0.10 - 0.50]     Bins: 1 5
	// Energy Range: [1.50 - 119.94]     Bins: 4 6 (old bins)
	// Asym Range:   [0.00 - 0.80]     Bins: 1 4
	// OpeningAngle Range:   [0.017 - 3.142]     Bins: 6 29
	// EOverP Cut    [0.60 - 1.40]
	Double_t fPi0MassFixedValue_1[kNMainCentBins][kNoGammaBins] = {
    { 0.140730, 0.135875, 0.139245, 0.146724, 0.158606, 0.173285, 0.190328, 0.200000, 0.200000},
    { 0.140730, 0.135875, 0.139245, 0.146724, 0.158606, 0.173285, 0.190328, 0.200000, 0.200000},
    { 0.140730, 0.135875, 0.139245, 0.146724, 0.158606, 0.173285, 0.190328, 0.200000, 0.200000},
    { 0.140730, 0.135875, 0.139245, 0.146724, 0.158606, 0.173285, 0.190328, 0.200000, 0.200000}
  };
	Double_t fPi0SigmaFixedValue_1[kNMainCentBins][kNoGammaBins] = {
    { 0.016273, 0.016512, 0.015727, 0.016063, 0.017414, 0.017614, 0.018374, 0.011323, 0.011323},
    { 0.016273, 0.016512, 0.015727, 0.016063, 0.017414, 0.017614, 0.018374, 0.011323, 0.011323},
    { 0.016273, 0.016512, 0.015727, 0.016063, 0.017414, 0.017614, 0.018374, 0.011323, 0.011323},
    { 0.016273, 0.016512, 0.015727, 0.016063, 0.017414, 0.017614, 0.018374, 0.011323, 0.011323}
  }; //9
	// These are the values add for MB data on 20180824
	// Lambda Range: [0.10 - 0.50]     Bins: 1 5
	// Energy Range: [1.50 - 119.94]     Bins: 4 6
	// Asym Range:   [0.00 - 0.80]     Bins: 1 4
	// OpeningAngle Range:   [0.017 - 3.142]     Bins: 6 29
  // EOverP Cut    [0.60 - 1.40]
	Double_t fPi0MassFixedValue_2[kNMainCentBins][kNoGammaBins] = {
    { 0.138988, 0.138666, 0.140109, 0.147061, 0.158484, 0.200000, 0.200000, 0.197369, 0.197369},
    { 0.138988, 0.138666, 0.140109, 0.147061, 0.158484, 0.200000, 0.200000, 0.197369, 0.197369},
    { 0.138988, 0.138666, 0.140109, 0.147061, 0.158484, 0.200000, 0.200000, 0.197369, 0.197369},
    { 0.138988, 0.138666, 0.140109, 0.147061, 0.158484, 0.200000, 0.200000, 0.197369, 0.197369}
  };
	Double_t fPi0SigmaFixedValue_2[kNMainCentBins][kNoGammaBins] = {
    { 0.013807, 0.014188, 0.014031, 0.015245, 0.014975, 0.160000, 0.160000, 0.020498, 0.020498},
    { 0.013807, 0.014188, 0.014031, 0.015245, 0.014975, 0.160000, 0.160000, 0.020498, 0.020498},
    { 0.013807, 0.014188, 0.014031, 0.015245, 0.014975, 0.160000, 0.160000, 0.020498, 0.020498},
    { 0.013807, 0.014188, 0.014031, 0.015245, 0.014975, 0.160000, 0.160000, 0.020498, 0.020498}
  };

  // These are the values for GA Data on 20191202
  // Lambda Range: [0.10 - 0.70]     Bins: 1 7
  // Energy Range: [2.00 - 100.00]     Bins: 1 3
  // Asym Range:   [0.00 - 0.70]     Bins: 1 1
  // OpeningAngle Range:   [0.017 - 3.142]     Bins: 3 12
  // Track Cluster Correction: Subtraction

  // note: last three bins are basically fake
	Double_t fPi0MassFixedValue_3[kNMainCentBins][kNoGammaBins] = {
    { 0.129265, 0.130740, 0.138696, 0.147779, 0.160138, 0.177054, 0.177054, 0.177054, 0.177054},
    { 0.131305, 0.131078, 0.137113, 0.143671, 0.153884, 0.171041, 0.171041, 0.171041, 0.171041},
    { 0.130332, 0.129157, 0.134172, 0.140263, 0.151399, 0.162242, 0.162242, 0.162242, 0.162242},
    { 0.127595, 0.128214, 0.132793, 0.137396, 0.148063, 0.167322, 0.167322, 0.167322, 0.167322}
  };
	Double_t fPi0SigmaFixedValue_3[kNMainCentBins][kNoGammaBins] = {
    { 0.019000, 0.019000, 0.017793, 0.017149, 0.017567, 0.011105, 0.011105, 0.011105, 0.011105},
    { 0.019000, 0.019000, 0.013467, 0.012001, 0.014112, 0.013643, 0.011105, 0.011105, 0.011105},
    { 0.019000, 0.016856, 0.010454, 0.010690, 0.011508, 0.011226, 0.011226, 0.011226, 0.011226},
    { 0.019000, 0.003000, 0.009186, 0.011264, 0.014117, 0.012515, 0.012515, 0.012515, 0.012515}
  };

  // These values correspond to MB Data on 20191202
  // Energy Range: [2.00 - 100.00]     Bins: 1 3
  // Asym Range:   [0.00 - 0.70]     Bins: 1 1
  // OpeningAngle Range:   [0.023 - 3.142]     Bins: 6 12
  // Track Cluster Correction: Subtraction
  // Lambda cuts vary by Centrality
  // Cent 0 : [0.10 - 0.50]     Bins: 1 5
  // Cent 1 : [0.10 - 0.50]     Bins: 1 5
  // Cent 2 : [0.10 - 0.60]     Bins: 1 6
  // Cent 3 : [0.10 - 0.40]     Bins: 1 3


	Double_t fPi0MassFixedValue_4[kNMainCentBins][kNoGammaBins] = {
    {  0.137945, 0.138191, 0.141075,  0.147679, 0.148967, 0.183629  ,0.183629,0.183629,0.183629 },
    {  0.136807, 0.136838, 0.138779,  0.143453, 0.150812, 0.137669  ,0.137669,0.137669,0.137669 },
    {  0.135318, 0.134796, 0.136419, 0.139745, 0.141703, 0.111335   ,0.111335,0.111335,0.111335 },
    {  0.135276, 0.135471, 0.136128, 0.141905, 0.101601, 0.126029   ,0.126029,0.126029,0.126029 }
  };
	Double_t fPi0SigmaFixedValue_4[kNMainCentBins][kNoGammaBins] = {
    {  0.012025, 0.012813, 0.012214, 0.011307, 0.019000, 0.019000  ,0.019000,0.019000,0.019000},
    {  0.010155, 0.010165, 0.009131, 0.009708, 0.019000, 0.018995  ,0.018995,0.018995,0.018995},
    {  0.008285, 0.008858, 0.008698, 0.010906, 0.018227, 0.016506  ,0.016506,0.016506,0.016506},
    {  0.006975, 0.006388, 0.005972, 0.007425, 0.005017, 0.010000  ,0.010000,0.010000,0.010000}
  };


  // These values correspond to GA Data on 20200913
  // Energy Range: [2.00 - 100.00]     Bins: 1 3
  // Asym Range:   [0.00 - 0.70]     Bins: 1 1
  // OpeningAngle Range:   [0.017 - 3.142]     Bins: 3 12
  // Track Cluster Correction: Subtraction
  // Lambda cuts vary by Centrality
  // Cent 0 : [0.10 - 0.50]     Bins: 1 5
  // Cent 1 : [0.10 - 0.50]     Bins: 1 5
  // Cent 2 : [0.10 - 0.50]     Bins: 1 5
  // Cent 3 : [0.10 - 0.50]     Bins: 1 5


	Double_t fPi0MassFixedValue_5[kNMainCentBins][kNoGammaBins] = {
    {  0.122980, 0.134701, 0.137247, 0.145853, 0.159203, 0.168348   ,0.15,0.15,0.15 },
    {  0.126280, 0.135862, 0.137113, 0.143644, 0.149588, 0.169486   ,0.15,0.15,0.15 },
    {  0.126438, 0.132486, 0.136093, 0.141735, 0.151793, 0.160520   ,0.15,0.15,0.15 },
    {  0.135224, 0.134146, 0.135899, 0.141734, 0.151186, 0.164889   ,0.15,0.15,0.15 }
  };
	Double_t fPi0SigmaFixedValue_5[kNMainCentBins][kNoGammaBins] = {
    {  0.005000, 0.009081, 0.010044, 0.011585, 0.014136, 0.012628  ,0.01,0.01,0.01},
    {  0.005000, 0.014835, 0.009292, 0.010112, 0.005000, 0.008841  ,0.01,0.01,0.01},
    {  0.009758, 0.007585, 0.008212, 0.008947, 0.009782, 0.007195  ,0.01,0.01,0.01},
    {  0.005971, 0.006350, 0.006648, 0.008418, 0.011478, 0.005000  ,0.01,0.01,0.01}
  };

// MB Corr4 20210204
  // Energy Range: [2.00 - 100.00]     Bins: 1 3
  // Asym Range:   [0.00 - 1.0]     Bins: 1 4
  // OpeningAngle Range:   [0.023 - 3.142]     Bins: 6 12
  // Track Cluster Correction: Subtraction
  // Lambda cuts vary by Centrality
  // Cent 0 : [0.10 - 0.50]     Bins: 1 5
  // Cent 1 : [0.10 - 0.50]     Bins: 1 5
  // Cent 2 : [0.10 - 0.50]     Bins: 1 5
  // Cent 3 : [0.10 - 0.50]     Bins: 1 5

	Double_t fPi0MassFixedValue_6[kNMainCentBins][kNoGammaBins] = {
    { 0.138998, 0.139368, 0.141303, 0.148221, 0.117443, 0.15, 0.15, 0.15, 0.15},
    { 0.137308, 0.137660, 0.139410, 0.144601, 0.107966, 0.15, 0.15, 0.15, 0.15},
    { 0.136108, 0.135959, 0.137340, 0.141563, 0.111171, 0.15, 0.15, 0.15, 0.15},
    { 0.135356, 0.135408, 0.135935, 0.141102, 0.115666, 0.15, 0.15, 0.15, 0.15}
  };
	Double_t fPi0SigmaFixedValue_6[kNMainCentBins][kNoGammaBins] = {
    { 0.011146, 0.011566, 0.011338, 0.011077, 0.013047, 0.01, 0.01, 0.01, 0.01},
    { 0.009240, 0.009207, 0.009374, 0.011105, 0.016010, 0.01, 0.01, 0.01, 0.01},
    { 0.007751, 0.008359, 0.007867, 0.006892, 0.019000, 0.01, 0.01, 0.01, 0.01},
    { 0.006977, 0.006724, 0.007251, 0.007834, 0.019000, 0.01, 0.01, 0.01, 0.01}
  };

// GA Corr4 20210204
  // Energy Range: [2.00 - 100.00]     Bins: 1 3
  // Asym Range:   [0.00 - 1.00]     Bins: 1 4
  // OpeningAngle Range:   [0.017 - 3.142]     Bins: 3 12
  // Track Cluster Correction: Subtraction
  // Lambda cuts vary by Centrality
  // Cent 0 : [0.10 - 0.50]     Bins: 1 5
  // Cent 1 : [0.10 - 0.50]     Bins: 1 5
  // Cent 2 : [0.10 - 0.50]     Bins: 1 5
  // Cent 3 : [0.10 - 0.50]     Bins: 1 5

	Double_t fPi0MassFixedValue_7[kNMainCentBins][kNoGammaBins] = {
    { 0.135693, 0.132496, 0.139298, 0.147529, 0.158139, 0.15, 0.15, 0.15, 0.15},
    { 0.135277, 0.132225, 0.137970, 0.144848, 0.153996, 0.15, 0.15, 0.15, 0.15},
    { 0.133821, 0.132782, 0.136673, 0.142265, 0.151094, 0.15, 0.15, 0.15, 0.15},
    { 0.133668, 0.132907, 0.135557, 0.140655, 0.150016, 0.15, 0.15, 0.15, 0.15}
  };
	Double_t fPi0SigmaFixedValue_7[kNMainCentBins][kNoGammaBins] = {
    { 0.010231, 0.017349, 0.012613, 0.012747, 0.013998, 0.01, 0.01, 0.01, 0.01},
    { 0.007734, 0.015694, 0.010948, 0.010878, 0.012375, 0.01, 0.01, 0.01, 0.01},
    { 0.006950, 0.005000, 0.008406, 0.008908, 0.009553, 0.01, 0.01, 0.01, 0.01},
    { 0.006978, 0.005000, 0.007047, 0.008444, 0.011205, 0.01, 0.01, 0.01, 0.01}
  };


	switch (input) {
    case 7:
			memcpy (fPi0MassFixed , fPi0MassFixedValue_7, sizeof(fPi0MassFixed));
			memcpy (fPi0SigmaFixed, fPi0SigmaFixedValue_7, sizeof(fPi0SigmaFixed));
      break;
    case 6:
			memcpy (fPi0MassFixed , fPi0MassFixedValue_6, sizeof(fPi0MassFixed));
			memcpy (fPi0SigmaFixed, fPi0SigmaFixedValue_6, sizeof(fPi0SigmaFixed));
      break;
    case 5:
			memcpy (fPi0MassFixed , fPi0MassFixedValue_5, sizeof(fPi0MassFixed));
			memcpy (fPi0SigmaFixed, fPi0SigmaFixedValue_5, sizeof(fPi0SigmaFixed));
      break;
    case 4:
			memcpy (fPi0MassFixed , fPi0MassFixedValue_4, sizeof(fPi0MassFixed));
			memcpy (fPi0SigmaFixed, fPi0SigmaFixedValue_4, sizeof(fPi0SigmaFixed));
      break;
    case 3:
			memcpy (fPi0MassFixed , fPi0MassFixedValue_3, sizeof(fPi0MassFixed));
			memcpy (fPi0SigmaFixed, fPi0SigmaFixedValue_3, sizeof(fPi0SigmaFixed));
      break;
		case 2:
			memcpy (fPi0MassFixed , fPi0MassFixedValue_2 , sizeof(fPi0MassFixed));
			memcpy (fPi0SigmaFixed, fPi0SigmaFixedValue_2, sizeof(fPi0SigmaFixed));
			break;
		case 1: 
			memcpy (fPi0MassFixed , fPi0MassFixedValue_1 , sizeof(fPi0MassFixed));
			memcpy (fPi0SigmaFixed, fPi0SigmaFixedValue_1, sizeof(fPi0SigmaFixed));
			break;
		case 0:
		default:
			memcpy (fPi0MassFixed , fPi0MassFixedValue_0 , sizeof(fPi0MassFixed));
			memcpy (fPi0SigmaFixed, fPi0SigmaFixedValue_0, sizeof(fPi0SigmaFixed));
	}
	fPi0MassSelection = input;
}
//________________________________________________________________________
AliAnalysisTaskGammaHadron::~AliAnalysisTaskGammaHadron()
{
	// Destructor

	//Copied from chris yaldo. Ask Salvatore about it!
	// Destructor. Clean-up the output list, but not the histograms that are put inside
	// (the list is owner and will clean-up these histograms). Protect in PROOF case.
	if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())
	{
		//delete fOutputList1;
	}
	//copied from hanseul
	/*if (fPoolMgr)
	{
		delete fPoolMgr;
	}*/
}
//________________________________________________________________________
TF1* AliAnalysisTaskGammaHadron::GetEffFunction(Int_t no,Int_t cent)
{
	TF1* function = 0;
	if(no==0)function=funcpT_low[cent];
	if(no==1)function=funcpT_high[cent];
	if(no==2)function=funcpEta_left[cent];
	if(no==3)function=funcpEta_right[cent];

	return function;
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::UserCreateOutputObjects()
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::UserCreateOutputObjects()"<<endl;

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Create fOutput list and histograms (fHistZVertex, fHistEventRejection, fHistEventRejection, fHistEventCount, fHistCentrality, fHistEventPlane)
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	AliAnalysisTaskEmcal::UserCreateOutputObjects();

  //fHistEventRejection->GetXaxis()->SetBinLabel(2,"Cent");


	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Add output for AliEventCuts
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	fEventCutList = new TList();
	fEventCutList ->SetOwner();
	fEventCutList ->SetName("EventCutOutput");

	fEventCuts.OverrideAutomaticTriggerSelection(fOffTrigger); //..otherwise only kINT7 events are used for the analysis
	if(fUseManualEventCuts==1)
	{
		//..Enable manual mode.
		//..this just means that the automatic cut settings
		//..are not loaded every time the event is checked
		fEventCuts.SetManualMode();
		fEventCuts.fMC = false; //FixMe substitute by a real flag in the task!
		fEventCuts.SetupLHC15o();
		fEventCuts.fUseVariablesCorrelationCuts = true; //..That is specifically for LHC15o!
	}

	fEventCuts.AddQAplotsToList(fEventCutList);
	fOutput->Add(fEventCutList);
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Create mixed event pools
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if(fSEvMEv==1 || fPoolMgr) //do this for either a mixed event analysis or when an external pool is given
	{
		InitEventMixer();
	}
	// Additional case: SE, but with Save Trigger
	if (fSEvMEv==0 && fSaveTriggerPool) {
		InitEventMixer(1);
	}
	if (fSEvMEv==2) {
		if (!fPoolMgr) {
			AliInfo("Mix Trigger Mode Called without Pool loaded.");
		}
		InitEventMixer(1);
	}

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //   Set up event plane objects for QnVectorFramework
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  AliWarning("Attempting to load Flow QnVector Task\n");
  AliAnalysisTaskFlowVectorCorrections * fFlowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections *> (AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
  if (fFlowQnVectorTask != NULL) {
    fFlowQnVectorMgr = fFlowQnVectorTask->GetAliQnCorrectionsManager();
    AliInfo("Successfully loaded QnVector Corrections");
  } else {
    AliError("Flow Qn Vector correction object not found. Will use uncorrected event plane angle from VZEROM.");
  }

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Define sublists/folders for a better organisation of the figures
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	fOutputListQA   = new TList();
	fOutputListQA   ->SetOwner();
	fOutputListQA   ->SetName("QA_histograms");

	Double_t pi = TMath::Pi();

	//common bins for the histograms
	Int_t nbins[6]  = {0};
	Double_t min[6] = {0};
	Double_t max[6] = {0};

	//settings for p_t cluster distributon
	nbins[0] = 500;
	min[0] = 0;
	max[0] = 100;
	//settings for p_t hadron distribution
	nbins[1] = 60;  //do 1/2 GeV bins so that you can see the 0.5 cut to set as a minimum pT to combine hadron and gamma
	min[1] = 0;
	max[1] = 30;
	//settings for delta phi (g-h) distribution
	nbins[2] = 45;
	min[2] = -90;
	max[2] = 270;
	//settings for delta eta (g-h) distribution
	nbins[3] = 80;
	min[3] = -2;
	max[3] = 2;
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//
	//   Create Histograms
	//
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  if (fIsMC) {
    // Cluster Matching
    fHistClusMCDE = new TH2D("fHistClusMCDE","fHistClusMCDE;E_{clus} (GeV);#Delta E (GeV)",100,0.,20.,100,-2.5,2.5);
    fHistClusMCDPhiDEta = new TH2D("fHistClusMCDPhiDEta","fHistClusMCDPhiDEta;#Delta #phi;#Delta #eta",100,-0.25,0.25,100,-0.25,0.25);
  }

  // Event Hash tracking histogram
  fHistEventHash = new TH1F("HistEventHash","Event Hash Value;Hash Value",11,-0.5,10.5);
  fOutput->Add(fHistEventHash);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   THn Sparse for the 2D histograms
	//   Dimensions are eta,phi
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Int_t dimThn = 0;
    TString titleThn[20];
    Int_t nbinsThn[20] = {0};
    Double_t minThn[20] = {0.};
    Double_t maxThn[20] = {0.};
    Double_t *binEdgesThn[20] = {0};

    titleThn[dimThn] = "#Delta #varphi";
    nbinsThn[dimThn] = 54;
    Double_t deltaPhiArray[54+1];
    binEdgesThn[dimThn] = deltaPhiArray;
    GenerateFixedBinArray(54,-pi/2.,3.*pi/2.,deltaPhiArray);
    minThn[dimThn] = -pi/2.;
    maxThn[dimThn] = 3.*pi/2.;
    dimThn++;

    titleThn[dimThn] = "#Delta #eta";
    nbinsThn[dimThn] = 80;
    //binEdgesThn[dim] = 80.;
    Double_t deltaEtaArray[80+1];
    binEdgesThn[dimThn] = deltaEtaArray;
    GenerateFixedBinArray(80,-2.,2.,deltaEtaArray);
    minThn[dimThn] = -2.;
    maxThn[dimThn] = 2.;
    dimThn++;

    titleThn[dimThn] = "E_{T} gamma";
    nbinsThn[dimThn] = kNoGammaBins;
    binEdgesThn[dimThn] = fArray_G_Bins;
    minThn[dimThn] = fArray_G_Bins[0];
    maxThn[dimThn] = fArray_G_Bins[kNoGammaBins];
    dimThn++;

    titleThn[dimThn] = "z_{T}";
    nbinsThn[dimThn] = kNoZtBins;
    binEdgesThn[dimThn] = fArray_ZT_Bins;
    minThn[dimThn] = fArray_ZT_Bins[0];
    maxThn[dimThn] = fArray_ZT_Bins[kNoZtBins];
    dimThn++;

    if (bEnableTrackPtAxis) {
      titleThn[dimThn] = "Track p_{T}";
      nbinsThn[dimThn] = kNoHPtBins;
      binEdgesThn[dimThn] = fArray_HPT_Bins;
      minThn[dimThn] = fArray_HPT_Bins[0];
      maxThn[dimThn] = fArray_HPT_Bins[kNoHPtBins];
      dimThn++;
    } else {
      titleThn[dimThn] = "#Xi";
      nbinsThn[dimThn] = kNoXiBins;
      binEdgesThn[dimThn] = fArray_XI_Bins;
      minThn[dimThn] = fArray_XI_Bins[0];
      maxThn[dimThn] = fArray_XI_Bins[kNoXiBins];
      dimThn++;
    }

    titleThn[dimThn] = "z vertex position";
    nbinsThn[dimThn] = kNvertBins;
    binEdgesThn[dimThn] = fArrayNVertBins;
    minThn[dimThn] = fArrayNVertBins[0];
    maxThn[dimThn] = fArrayNVertBins[kNvertBins];
    dimThn++;
/*
    titleThn[dimThn] = "ME(0) or SE(1)";
    static const Int_t nSeMeBins=2;
 	Double_t SeMeBinArray[nSeMeBins+1]  = {0,1,2};
    nbinsThn[dimThn] = nSeMeBins;
    binEdgesThn[dimThn] = SeMeBinArray;
    minThn[dimThn] = SeMeBinArray[0];
    maxThn[dimThn] = SeMeBinArray[nSeMeBins];
    dimThn++;
*/
    static const Int_t nEvtPlaneBins=3;
    Double_t evtPlaneArray[nEvtPlaneBins+1] = {0,1,2,3};// = {In Plane, MP, Out of Plane};
    static const Int_t nCentHistBins=4;
    Double_t centBinArray[nCentHistBins+1]  = {0.0,10.0,30.0,50.0,90.0};
    if(fForceBeamType != AliAnalysisTaskEmcal::kpp)
    {
     	//..Event plane
     	titleThn[dimThn] = "IP(0), MP(1), OP(2)";
     	nbinsThn[dimThn] = nEvtPlaneBins;
     	binEdgesThn[dimThn] = evtPlaneArray;
        minThn[dimThn] = evtPlaneArray[0];
     	maxThn[dimThn] = evtPlaneArray[nEvtPlaneBins];
     	dimThn++;

     	//..Centrality
     	titleThn[dimThn] = "Centrality %";
     	nbinsThn[dimThn] = nCentHistBins;
     	binEdgesThn[dimThn] = centBinArray;
     	minThn[dimThn] = centBinArray[0];
     	maxThn[dimThn] = centBinArray[nCentHistBins];
     	dimThn++;
    }
    int nCorrMCStatusBins = 3;
    double corrMCStatusArray[3]; // 0 = (background), 1 = (true eta->2gamma), 2 = (true pi0->2gamma)
    if (fIsMC) {
      titleThn[dimThn] = "MC Status";
      nbinsThn[dimThn] = nCorrMCStatusBins;
      GenerateFixedBinArray(3,0,3,corrMCStatusArray);
      binEdgesThn[dimThn] = corrMCStatusArray;
      minThn[dimThn] = corrMCStatusArray[0];
      maxThn[dimThn] = corrMCStatusArray[nCorrMCStatusBins];
      dimThn++;
    }

    if(fPlotQA==0)
    {
     	fCorrVsManyThings   = new THnSparseF("CorrVsManyThings", "CorrVsManyThings", dimThn, nbinsThn, minThn, maxThn);
  
     	for(Int_t i=0;i<dimThn;i++)
     	{
     		fCorrVsManyThings->GetAxis(i)->SetTitle(titleThn[i]);
     		fCorrVsManyThings->SetBinEdges(i, binEdgesThn[i]);
     	}

      // Creating the trigger thnSparse:
      const Int_t dimThnTrig = 5;
      //TString titleThnTrig[20];
      Int_t nbinsThnTrig[dimThnTrig]  = {0};
      Double_t minThnTrig[dimThnTrig] = {0.};
      Double_t maxThnTrig[dimThnTrig] = {0.};

      Int_t trigMap[dimThnTrig] = {2,5,6,7,8};  // Which axes to include in trigger ThnSparse

      Int_t usedDimThnTrig = 5; // variable one, so the MC one can be removed
      if (!fIsMC) { usedDimThnTrig = 4;}

      for(Int_t i = 0; i < usedDimThnTrig; i++)
      {
        nbinsThnTrig[i] = nbinsThn[trigMap[i]];
        minThnTrig[i]   = minThn[trigMap[i]];
        maxThnTrig[i]   = maxThn[trigMap[i]];
      }

      fTriggerHist = new THnSparseF("TriggerHist","TriggerHist",usedDimThnTrig,nbinsThnTrig,minThnTrig,maxThnTrig);
  
      for(Int_t i = 0; i < usedDimThnTrig; i++)
      {
        fTriggerHist->GetAxis(i)->SetTitle(titleThn[trigMap[i]]);
        fTriggerHist->SetBinEdges(i,binEdgesThn[trigMap[i]]);
      }      

     	//fCorrVsManyThings->Sumw2();
     	fOutput->Add(fCorrVsManyThings);
      fOutput->Add(fTriggerHist);
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   THn Sparse for the Pi0 Candidates
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Int_t dimThnPi0 = 0;
    TString titleThnPi0[11];
    Int_t nBinsThnPi0[11] = {0};
    Double_t minThnPi0[11] = {0.};
    Double_t maxThnPi0[11] = {0.};
    Double_t *binEdgesThnPi0[11] = {0};

    titleThnPi0[dimThnPi0] = "p_{T}^{#gamma#gamma}";
    nBinsThnPi0[dimThnPi0] = 60;
    Double_t pTArray[60+1];
    binEdgesThnPi0[dimThnPi0] = pTArray;
    GenerateFixedBinArray(60,0,30,pTArray);
    minThnPi0[dimThnPi0] = 0;
    maxThnPi0[dimThnPi0] = 30;
    dimThnPi0++;

    titleThnPi0[dimThnPi0] = "M_{#gamma#gamma}";
    nBinsThnPi0[dimThnPi0] = 1000;
    Double_t mGGArray[1000+1];
    binEdgesThnPi0[dimThnPi0] = mGGArray;
    GenerateFixedBinArray(1000,0,1.0,mGGArray);
    minThnPi0[dimThnPi0] = 0;
    maxThnPi0[dimThnPi0] = 1.0;
/*    nBinsThnPi0[dimThnPi0] = 750;
    Double_t mGGArray[750+1];
    binEdgesThnPi0[dimThnPi0] = mGGArray;
    GenerateFixedBinArray(750,0,0.75,mGGArray);
    minThnPi0[dimThnPi0] = 0;
    maxThnPi0[dimThnPi0] = 0.75;*/
    dimThnPi0++;

    titleThnPi0[dimThnPi0] = "Opening Angle [rad]";
    nBinsThnPi0[dimThnPi0] = 12;
    // Small Scan
    //Double_t openAngleArray[6+1] = {0,0.015,0.017,0.019,0.021,0.23,pi};
    // Medium Scan
    Double_t openAngleArray[12+1] = {0,0.015,0.017,0.019,0.021,0.023,0.025,0.027,0.029,0.031,0.033,0.035,pi};
    // Large Scan
//    Double_t openAngleArray[29+1] = {0,0.009,0.011,0.013,0.015,0.017,0.019,0.021,0.023,0.025,0.027,0.029,0.031,0.033,0.035,0.037,0.04,0.044,0.048,0.055,0.07,0.08,0.09,0.1,0.125,0.15,pi/16.,pi/8.,pi/4.,pi};
    binEdgesThnPi0[dimThnPi0] = openAngleArray;
  //  GenerateFixedBinArray(20,0,TMath::Pi(),openAngleArray);
    minThnPi0[dimThnPi0] = 0;
    maxThnPi0[dimThnPi0] = pi;
    dimThnPi0++;

    titleThnPi0[dimThnPi0] = "Max Lambda_{0}^{2}";
    nBinsThnPi0[dimThnPi0] = 8;
    Double_t MaxM02Array[8+1] = {0.1, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 1.0};
//    Double_t MaxM02Array[6+1] = {0.1,0.3,0.35,0.4,0.45,0.5,10};
    binEdgesThnPi0[dimThnPi0] = MaxM02Array;
   // GenerateFixedBinArray(5,0,);
    minThnPi0[dimThnPi0] = 0.1;
    maxThnPi0[dimThnPi0] = 1.0;
    dimThnPi0++;

    titleThnPi0[dimThnPi0] = "Min Cluster Energy";
    nBinsThnPi0[dimThnPi0] = 3;
    Double_t MinClusEnergyArray[3+1] = {2,2.5,3,100};
    binEdgesThnPi0[dimThnPi0] = MinClusEnergyArray;
    minThnPi0[dimThnPi0] = 2;
    maxThnPi0[dimThnPi0] = 100;
/*    nBinsThnPi0[dimThnPi0] = 5;
    Double_t MinClusEnergyArray[5+1] = {0.3,0.5,1,1.5,2,100};
    binEdgesThnPi0[dimThnPi0] = MinClusEnergyArray;
    minThnPi0[dimThnPi0] = 0.3;
    maxThnPi0[dimThnPi0] = 100;*/
    dimThnPi0++;

    titleThnPi0[dimThnPi0] = "Asymmetry";
    nBinsThnPi0[dimThnPi0] = 4;
    Double_t asymmetryArray[4+1] = {0.,0.7,0.8,0.9,1.0};
    binEdgesThnPi0[dimThnPi0] = asymmetryArray;
    // More useful range
    // Old Range
    //GenerateFixedBinArray(5,0,1,asymmetryArray);
    minThnPi0[dimThnPi0] = 0;
    maxThnPi0[dimThnPi0] = 1;
    dimThnPi0++;

    //..ID array -  0 - real, 1 - rotated, 2 - mixed event, 3 - pos-swapped
    titleThnPi0[dimThnPi0] = "Type";
    nBinsThnPi0[dimThnPi0] = 4;
    Double_t mRotArray[4+1];
    binEdgesThnPi0[dimThnPi0] = mRotArray;
    GenerateFixedBinArray(4,0,4,mRotArray);
    minThnPi0[dimThnPi0] = 0;
    maxThnPi0[dimThnPi0] = 4;
    dimThnPi0++;

    // FIXME add number of matched tracks?
 /*   titleThnPi0[dimThnPi0] = "NMatched";
    nBinsThnPi0[dimThnPi0] = 10;
    Double_t mNMatchedArray[10+1] = {0,1,2,3,4,5,6,8,10,15,20};
    binEdgesThnPi0[dimThnPi0] = mNMatchedArray;
    minThnPi0[dimThnPi0] = 0;
    maxThnPi0[dimThnPi0] = 10;
    dimThnPi0++;*/

    // Event Plane Angle
    titleThnPi0[dimThnPi0] = "EPAngle";
    nBinsThnPi0[dimThnPi0] = 3;
    Double_t epArray[3+1] = {0,1,2,3};
    binEdgesThnPi0[dimThnPi0] = epArray;
    minThnPi0[dimThnPi0] = 0;
    maxThnPi0[dimThnPi0] = 3;
    dimThnPi0++;

    Double_t mcStatusArray[14+1];
    if (fIsMC) {
      //..MC Status array:
			// 0 - no Match, 1 - Single Particle to two clusters, 2 - pi0 to two gamma
			// 3 - pi0 Dalitz, 4 - Eta to two gamma, 5 - Eta to three pi0
			// 6 - eta to pi0pi+pi0 or pi0twogamma, 7 - eta to gammapi+pi- or gammaf+f-
			// 8 - gamma to e+e-, 9 - other shared ancestor (not in previous categories
      // 10 PosSwapped pi0->2gamma (pair energies conserved)
      // 11 PosSwapped pi0->2gamma (pair position conserved)
      // 12 PosSwapped eta->2gamma (pair energies conserved)
      // 13 PosSwapped eta->2gamma (pair position conserved)
      titleThnPi0[dimThnPi0] = "MC Match Status";
      nBinsThnPi0[dimThnPi0] = 14;
      binEdgesThnPi0[dimThnPi0] = mcStatusArray;
      GenerateFixedBinArray(14,0,14,mcStatusArray);
      minThnPi0[dimThnPi0] = 0;
      maxThnPi0[dimThnPi0] = 14;
      dimThnPi0++;
    }

    Double_t fPatchStatusArray[2+1];
    if (fApplyPatchCandCut) {
      // 1 if cluster pair is candidate for GA trigger patch (0 if not)
      titleThnPi0[dimThnPi0] = "Patch Status";
      nBinsThnPi0[dimThnPi0] = 2;
      binEdgesThnPi0[dimThnPi0] = fPatchStatusArray;
      GenerateFixedBinArray(2,0,2,fPatchStatusArray);
      minThnPi0[dimThnPi0]=0;
      maxThnPi0[dimThnPi0]=2;
      dimThnPi0++;
    }

    // Information for modification histograms
    // Axes for modification arrays
    Int_t dimThnMod = 0;
    TString titleThnMod[8];
    Int_t nBinsThnMod[8] = {0};
    Double_t minThnMod[8] = {0.};
    Double_t maxThnMod[8] = {0.};
    Double_t *binEdgesThnMod[8] = {0};
/*
    titleThnMod[dimThnMod] = "Modification"; // could replace later
    nBinsThnMod[dimThnMod] = 1000;
    Double_t ModArray[1000+1];
    binEdgesThnMod[dimThnMod] = ModArray;
    GenerateFixedBinArray(1000,0,10.,ModArray);
    minThnMod[dimThnMod] = 0.0;
    maxThnMod[dimThnMod] = 10.;
    dimThnMod++;
*/
    titleThnMod[dimThnMod] = "Initial Mass";
    nBinsThnMod[dimThnMod] = 500;
    Double_t InitialMassArray[500+1];
    binEdgesThnMod[dimThnMod] = InitialMassArray;
    GenerateFixedBinArray(500,0,1.,InitialMassArray);
    minThnMod[dimThnMod] = 0.0;
    maxThnMod[dimThnMod] = 1.0;
    dimThnMod++;

    titleThnMod[dimThnMod] = "Initial pT";
    nBinsThnMod[dimThnMod] = 60;
    Double_t InitialPtArray[60+1];
    binEdgesThnMod[dimThnMod] = InitialPtArray;
    GenerateFixedBinArray(60,0,30,InitialPtArray);
    minThnMod[dimThnMod] =  0.0;
    maxThnMod[dimThnMod] = 30.0;
    dimThnMod++;

    titleThnMod[dimThnMod] = "Final Mass";
    nBinsThnMod[dimThnMod] = 500;
    Double_t FinalMassArray[500+1];
    binEdgesThnMod[dimThnMod] = FinalMassArray;
    GenerateFixedBinArray(500,0,1.,FinalMassArray);
    minThnMod[dimThnMod] = 0.0;
    maxThnMod[dimThnMod] = 1.0;
    dimThnMod++;

    titleThnMod[dimThnMod] = "Final pT";
    nBinsThnMod[dimThnMod] = 60;
    Double_t FinalPtArray[60+1];
    binEdgesThnMod[dimThnMod] = FinalPtArray;
    GenerateFixedBinArray(60,0,30,FinalPtArray);
    minThnMod[dimThnMod] =  0.0;
    maxThnMod[dimThnMod] = 30.0;
    dimThnMod++;


    titleThnMod[dimThnMod] = "Max Lambda_{0}^{2}";
    nBinsThnMod[dimThnMod] = 8;
    Double_t ModMaxM02Array[8+1] = {0.1, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 1.0};
    binEdgesThnMod[dimThnMod] = ModMaxM02Array;
    minThnMod[dimThnMod] = 0.1;
    maxThnMod[dimThnMod] = 1.0;
    dimThnMod++;

    titleThnMod[dimThnMod] = "Min Cluster Energy";
    nBinsThnMod[dimThnMod] = 3;
    binEdgesThnMod[dimThnMod] = MinClusEnergyArray;
    minThnMod[dimThnMod] = 2;
    maxThnMod[dimThnMod] = 100;
    dimThnMod++;

    // Event Plane Angle
    titleThnMod[dimThnMod] = "EPAngle";
    nBinsThnMod[dimThnMod] = 3;
    Double_t ModEPArray[3+1] = {0,1,2,3};
    binEdgesThnMod[dimThnMod] = ModEPArray;
    minThnMod[dimThnMod] = 0;
    maxThnMod[dimThnMod] = 3;
    dimThnMod++;

    // Axes for 2D Scaling Modification Arrays
    Int_t dimThnModMatrix = 0;
    TString titleThnModMatrix[7];
    Int_t nBinsThnModMatrix[7] = {0};
    Double_t minThnModMatrix[7] = {0.};
    Double_t maxThnModMatrix[7] = {0.};
    Double_t *binEdgesThnModMatrix[7] = {0};

    const Int_t nBins2DMod = 500;
    nBinsThnModMatrix[dimThnModMatrix] = nBins2DMod;
    Double_t ModMatrixArray[nBins2DMod+1];
    binEdgesThnModMatrix[dimThnModMatrix] = ModMatrixArray;
    if (bLogPSMod) {
      // Log Version
      titleThnModMatrix[dimThnModMatrix] = "Log Mass Scaling";
      GenerateFixedBinArray(nBins2DMod,-5.,5.,ModMatrixArray);
      minThnModMatrix[dimThnModMatrix] = -5.;
      maxThnModMatrix[dimThnModMatrix] = 5.;
    } else {
      // Linear Version
      titleThnModMatrix[dimThnModMatrix] = "Mass Scaling";
      GenerateFixedBinArray(nBins2DMod,0,10.,ModMatrixArray);
      minThnModMatrix[dimThnModMatrix] = 0.0;
      maxThnModMatrix[dimThnModMatrix] = 10.;
    }
    dimThnModMatrix++;

    nBinsThnModMatrix[dimThnModMatrix] = nBins2DMod;
    Double_t ModMatrixArray2[nBins2DMod+1];
    binEdgesThnModMatrix[dimThnModMatrix] = ModMatrixArray2;
    if (bLogPSMod) {
      // Log Version
      titleThnModMatrix[dimThnModMatrix] = "Log pT Scaling";
      GenerateFixedBinArray(nBins2DMod,-5.,5.,ModMatrixArray2);
      minThnModMatrix[dimThnModMatrix] = -5.;
      maxThnModMatrix[dimThnModMatrix] = 5.;
    } else {
      // Linear Version
      titleThnModMatrix[dimThnModMatrix] = "pT Scaling";
      GenerateFixedBinArray(nBins2DMod,0,10.,ModMatrixArray2);
      minThnModMatrix[dimThnModMatrix] = 0.0;
      maxThnModMatrix[dimThnModMatrix] = 10.;
    }
    dimThnModMatrix++;

    titleThnModMatrix[dimThnModMatrix] = "Max Lambda_{0}^{2}";
    nBinsThnModMatrix[dimThnModMatrix] = 8;
    Double_t ModMatrixMaxM02Array[8+1] = {0.1, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 1.0};
    binEdgesThnModMatrix[dimThnModMatrix] = ModMatrixMaxM02Array;
    minThnModMatrix[dimThnModMatrix] = 0.1;
    maxThnModMatrix[dimThnModMatrix] = 1.0;
    dimThnModMatrix++;

    titleThnModMatrix[dimThnModMatrix] = "Min Cluster Energy";
    nBinsThnModMatrix[dimThnModMatrix] = 3;
    binEdgesThnModMatrix[dimThnModMatrix] = MinClusEnergyArray;
    minThnModMatrix[dimThnModMatrix] = 2;
    maxThnModMatrix[dimThnModMatrix] = 100;
    dimThnModMatrix++;

    if ( fGammaOrPi0>0  && fPlotQA==1)
    {
      fPi0Cands= new THnSparseF("Pi0Cands", "Pi0Cands", dimThnPi0, nBinsThnPi0, minThnPi0, maxThnPi0);
      for(Int_t i=0;i<dimThnPi0;i++)
      {
        fPi0Cands->GetAxis(i)->SetTitle(titleThnPi0[i]);
        fPi0Cands->SetBinEdges(i, binEdgesThnPi0[i]);
      }
      fOutput->Add(fPi0Cands);

      // Initializing the PS Modification histograms
      if (bEnablePosSwapHists) {
        fPSMassPtMap = new THnSparseF("fPSMassPtMap","fPSMassPtMap",dimThnMod,nBinsThnMod,minThnMod,maxThnMod);
        fESMassPtMap = new THnSparseF("fESMassPtMap","fESMassPtMap",dimThnMod,nBinsThnMod,minThnMod,maxThnMod);
        for(Int_t i=0;i<dimThnMod;i++)
        {
          fPSMassPtMap->GetAxis(i)->SetTitle(titleThnMod[i]);
          fESMassPtMap->GetAxis(i)->SetTitle(titleThnMod[i]);

          fPSMassPtMap->SetBinEdges(i, binEdgesThnMod[i]);
          fESMassPtMap->SetBinEdges(i, binEdgesThnMod[i]);
        }
        fOutput->Add(fPSMassPtMap);
        fOutput->Add(fESMassPtMap);

        fUScaleMatrix = new THnSparseF("UScaleMatrix","UScaleMatrix",dimThnModMatrix,nBinsThnModMatrix,minThnModMatrix,maxThnModMatrix);
        fVScaleMatrix = new THnSparseF("VScaleMatrix","VScaleMatrix",dimThnModMatrix,nBinsThnModMatrix,minThnModMatrix,maxThnModMatrix);
        for(Int_t i=0;i<dimThnModMatrix;i++)
        {
          fUScaleMatrix->GetAxis(i)->SetTitle(titleThnModMatrix[i]);
          fVScaleMatrix->GetAxis(i)->SetTitle(titleThnModMatrix[i]);

          fUScaleMatrix->SetBinEdges(i, binEdgesThnModMatrix[i]);
          fVScaleMatrix->SetBinEdges(i, binEdgesThnModMatrix[i]);
        }
        fOutput->Add(fUScaleMatrix);
        fOutput->Add(fVScaleMatrix);
      }
//			Double_t fBinsMixedClusZvtx[nBinsMixedClusZvtx+1] = {-10., -5.,-3.,-1.,1.,3.,5.,10.};
//			Double_t fBinsEMCalMult[nBinsEMCalMult + 1] = {0.,50.,100.,150.,200.,250.,300,500,700,900,1200};
//			fEMCalMultvZvtx = new TH2D("EMCalMultvZvtx","fEMCalMultvZvtx",nBinsMixedClusZvtx,fBinsMixedClusZvtx,nBinsEMCalMult,fBinsEMCalMult);

			if (fIsMC) {
				// Cluster Matching
				//fHistClusMCDE = new TH2D("fHistClusMCDE","fHistClusMCDE;E_{clus} (GeV);#Delta E (GeV)",100,0.,20.,100,-2.5,2.5);
				//fHistClusMCDPhiDEta = new TH2D("fHistClusMCDPhiDEta","fHistClusMCDPhiDEta;#Delta #phi;#Delta #eta",100,-0.25,0.25,100,-0.25,0.25);

				// Cluster Pair Matching
				fHistPi0MCDPt = new TH2D("fHistMatchPi0MCDPt","fHistMatchPi0MCDPt;p_{T}^{#gamma#gamma} (GeV/c);#Delta p_{T} (GeV/c)",100,0.,20.,100,-2.5,2.5);
				fHistEtaMCDPt = new TH2D("fHistMatchEtaMCDPt","fHistMatchEtaMCDPt;p_{T}^{#gamma#gamma} (GeV/c);#Delta p_{T} (GeV/c)",100,0.,20.,100,-2.5,2.5);
				fHistPi0MCDPhiDEta = new TH2D("fHistPi0MCDPhiDEta","fHistPi0MCDPhiDEta;#Delta #phi;#Delta #eta",100,-0.25,0.25,100,-0.25,0.25);
				fHistEtaMCDPhiDEta = new TH2D("fHistEtaMCDPhiDEta","fHistEtaMCDPhiDEta;#Delta #phi;#Delta #eta",100,-0.25,0.25,100,-0.25,0.25);

				fOutput->Add(fHistClusMCDE);
				fOutput->Add(fHistClusMCDPhiDEta);

				fOutput->Add(fHistPi0MCDPt);
				fOutput->Add(fHistEtaMCDPt);
				fOutput->Add(fHistPi0MCDPhiDEta);
				fOutput->Add(fHistEtaMCDPhiDEta);
			}

			if (fMixBClusZvtx && fMixBEMCalMult) {
			fEMCalMultvZvtx = new TH2D("EMCalMultvZvtx","fEMCalMultvZvtx",kNClusVertBins,fMixBClusZvtx->GetXbins()->GetArray(),kNEMCalMultBins,fMixBEMCalMult->GetXbins()->GetArray());
			fOutput->Add(fEMCalMultvZvtx);
			}
			// Initialize EventPoolManager for Cluster Mixing
			if (fDoClusMixing) {
				InitClusMixer();
			}
    }


	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //   THn Sparse for the Cluster properties
    //
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Int_t dimThnQA = 0;
    TString titleThnQA[11];
    Int_t nbinsThnQA[11] = {0};
    Double_t minThnQA[11] = {0.};
    Double_t maxThnQA[11] = {0.};
    Double_t *binEdgesThnQA[11] = {0};

    //..E axis
    TCustomBinning xBinningE;
    xBinningE.SetMinimum(0);
    xBinningE.AddStep(16,0.4); //40 bins
    xBinningE.AddStep(25,0.6); //15 bins
    TArrayD xbinsArrayE;
    xBinningE.CreateBinEdges(xbinsArrayE);

    titleThnQA[dimThnQA] = "E_{#gamma}";
    nbinsThnQA[dimThnQA] = 55;
    Double_t EgArray[55+1];
	for(Int_t i=0;i<56;i++)
    {
		EgArray[i]=xbinsArrayE.At(i);
    }
    binEdgesThnQA[dimThnQA] = EgArray;
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 25;
    dimThnQA++;

    //..L0 Axis
    //..Create the L0 axis with increasing binwidth
    TCustomBinning xBinning;
    xBinning.SetMinimum(0);
    xBinning.AddStep(0.4,0.01);  //40 (100)..first entries of the array are the set ranges and bins
    xBinning.AddStep(1,0.025);   //24..expand the previously defined range but increase the bin width
    xBinning.AddStep(4,0.1);     //30..expand the previously defined range but increase the bin width
    TArrayD xbinsArray;
    xBinning.CreateBinEdges(xbinsArray);

    titleThnQA[dimThnQA] = "#lambda_{0}";
    nbinsThnQA[dimThnQA] = 94;
    Double_t ShapeArray[94+1];
	for(Int_t i=0;i<95;i++)
    {
    		ShapeArray[i]=xbinsArray.At(i);
    }
    binEdgesThnQA[dimThnQA] = ShapeArray;
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 4;
    dimThnQA++;

    Double_t distanceArray[5+1];
    Double_t distanceArray2[5+1];
    if(fPlotQA==1)
    {
    	titleThnQA[dimThnQA] = "cell distance to bad channel";
    	nbinsThnQA[dimThnQA] = 5;
    	binEdgesThnQA[dimThnQA] = distanceArray;
    	GenerateFixedBinArray(5,0,5,distanceArray);
    	minThnQA[dimThnQA] = 0;
    	maxThnQA[dimThnQA] = 5;
    	dimThnQA++;

    	titleThnQA[dimThnQA] = "cell distance to SM edge";
    	nbinsThnQA[dimThnQA] = 5;
    	binEdgesThnQA[dimThnQA] = distanceArray2;
    	GenerateFixedBinArray(5,0,5,distanceArray2);
    	minThnQA[dimThnQA] = 0;
    	maxThnQA[dimThnQA] = 5;
    	dimThnQA++;
    }

    TCustomBinning xBinning2;
    xBinning2.SetMinimum(0);
    xBinning2.AddStep(0.04,0.002);     //..first entries of the array are the set ranges and bins
    xBinning2.AddStep(0.1,0.004);      //..expand the previously defined range but increase the bin width
    Double_t etaArrayDistMatched[35+1];
    Double_t phiArrayDistMatched[35+1];
    TArrayD xbinsArray2;
    xBinning2.CreateBinEdges(xbinsArray2);

    for(Int_t i=0;i<36;i++)
    {
       	etaArrayDistMatched[i]=xbinsArray2.At(i);
       	phiArrayDistMatched[i]=xbinsArray2.At(i);
    }

    if(fPlotQA==2)
    {
    	titleThnQA[dimThnQA] = "#Delta #eta^{match. track-cluster}";
    	nbinsThnQA[dimThnQA] = 35;
    	binEdgesThnQA[dimThnQA] = etaArrayDistMatched;
    	minThnQA[dimThnQA] = 0;
    	maxThnQA[dimThnQA] = 35;
    	dimThnQA++;

    	titleThnQA[dimThnQA] = "#Delta #varphi^{match. track-cluster}";
    	nbinsThnQA[dimThnQA] = 35;
    	binEdgesThnQA[dimThnQA] = phiArrayDistMatched;
    	minThnQA[dimThnQA] = 0;
    	maxThnQA[dimThnQA] = 35;
    	dimThnQA++;
    }

    Double_t etaArray[142+1];
    Double_t phiArray[311+1];
    if(fPlotQA==1)
    {
    	titleThnQA[dimThnQA] = "#eta^{cluster}";
    	nbinsThnQA[dimThnQA] = 142;
    	binEdgesThnQA[dimThnQA] = etaArray;
    	GenerateFixedBinArray(142,-0.71,0.71,etaArray);
    	minThnQA[dimThnQA] = 0;
    	maxThnQA[dimThnQA] = 142;
    	dimThnQA++;

    	titleThnQA[dimThnQA] = "#varphi^{cluster}";
    	nbinsThnQA[dimThnQA] = 311;
    	binEdgesThnQA[dimThnQA] = phiArray;
    	GenerateFixedBinArray(311,75.9,331.1,phiArray);
    	minThnQA[dimThnQA] = 0;
    	maxThnQA[dimThnQA] = 311;
    	dimThnQA++;
    }

    Double_t EPArray[100+1];
//    static const Int_t nCentHistBins=4; // These are define already now
//    Double_t centBinArray[nCentHistBins+1]  = {0.0,10.0,30.0,60.0,100.0};

    if(fPlotQA==2)
    {
    		//..E/p for electron identification
        titleThnQA[dimThnQA] = "E/p";
        nbinsThnQA[dimThnQA] = 100;
        binEdgesThnQA[dimThnQA] = EPArray;
        GenerateFixedBinArray(100,0,2,EPArray);
        minThnQA[dimThnQA] = 0;
        maxThnQA[dimThnQA] = 2;
        dimThnQA++;

        //..Centrality
    		titleThnQA[dimThnQA] = "Centrality %";
    		nbinsThnQA[dimThnQA] = nCentHistBins;
    		binEdgesThnQA[dimThnQA] = centBinArray;
    		minThnQA[dimThnQA] = centBinArray[0];
    		maxThnQA[dimThnQA] = centBinArray[nCentHistBins];
    		dimThnQA++;
    }
    /*//..ID array -  0 - leading, 1- track matched, 2 - leading & track matched
    titleThnQA[dimThnQA] = "ID code of photon";
    nbinsThnQA[dimThnQA] = 5;
    Double_t IdArray[5+1];
    binEdgesThnQA[dimThnQA] = IdArray;
    GenerateFixedBinArray(5,0,5,IdArray);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 5;
    dimThnQA++;
     */
    /*
        titleThnQA[dimThnQA] = "NLM";
        nbinsThnQA[dimThnQA] = 4;
        Double_t NLMArray[4+1];
        binEdgesThnQA[dimThnQA] = NLMArray;
        GenerateFixedBinArray(4,0,4,NLMArray);
        minThnQA[dimThnQA] = 0;
        maxThnQA[dimThnQA] = 4;
        dimThnQA++;
    */
    /*
        titleThnQA[dimThnQA] = "#Cells";
        nbinsThnQA[dimThnQA] = 10;
        Double_t nCellsArray[10+1];
        binEdgesThnQA[dimThnQA] = nCellsArray;
        GenerateFixedBinArray(10,0,10,nCellsArray);
        minThnQA[dimThnQA] = 0;
        maxThnQA[dimThnQA] = 10;
        dimThnQA++;
     */
     /*
        titleThnQA[dimThnQA] = "M_{#gamma#gamma x}";
        nbinsThnQA[dimThnQA] = 1000;
        Double_t MArray[1000+1];
        binEdgesThnQA[dimThnQA] = MArray;
        GenerateFixedBinArray(1000,0,10,MArray);
        minThnQA[dimThnQA] = 0;
        maxThnQA[dimThnQA] = 10;
        dimThnQA++;
    */
    //..additional things to put inside: time
    if(fPlotQA>0 && fGammaOrPi0==0)
    {
    		fClusterProp= new THnSparseF("ClusterProp", "ClusterProp", dimThnQA, nbinsThnQA, minThnQA, maxThnQA);
    		for(Int_t i=0;i<dimThnQA;i++)
    		{
    			fClusterProp->GetAxis(i)->SetTitle(titleThnQA[i]);
    			fClusterProp->SetBinEdges(i, binEdgesThnQA[i]);
    		}
    		fOutput->Add(fClusterProp);
    }

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//    Histograms for common use
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//..Initialize
	fHistBinCheckPt       = new TH1*[kNIdentifier];
	fHistBinCheckZt       = new TH1*[kNIdentifier];
	fHistBinCheckXi       = new TH1*[kNIdentifier];
	fHistBinCheckEvtPl    = new TH1*[kNIdentifier];
	fHistBinCheckEvtPl2   = new TH1*[kNIdentifier];
	fHistDEtaDPhiGammaQA  = new TH2*[kNIdentifier+3]; //Why +2?? -> because I want more than 3 QA versions
	fHistDEtaDPhiTrackQA  = new TH2*[kNIdentifier+3];
	fHistClusterTime      = new TH2*[kNIdentifier+3];

	//..by the identifier different histograms can be filled under different cut conditions
	//..Can eg. later be modified to contain certain delta phi or centrality bins
	for(Int_t identifier=0;identifier<kNIdentifier;identifier++)
	{
		fHistBinCheckPt[identifier] = new TH1D(Form("fHistBinCheckPt_%0d",identifier),Form("fHistBinCheckPt_%0d",identifier), nbins[0], min[0], max[0]);
		fHistBinCheckPt[identifier]->GetXaxis()->SetTitle("p_{T}^{#gamma}");
		fHistBinCheckPt[identifier]->GetYaxis()->SetTitle("Entries");
		fOutput->Add(fHistBinCheckPt[identifier]);

		fHistBinCheckZt[identifier] = new TH1D(Form("fHistBinCheckZt_%0d",identifier),Form("fHistBinCheckZt_%0d",identifier), 1500, 0, 60);
		fHistBinCheckZt[identifier]->GetXaxis()->SetTitle("z_{T}^{#gamma-h}");
		fHistBinCheckZt[identifier]->GetYaxis()->SetTitle("Entries");
		fOutput->Add(fHistBinCheckZt[identifier]);

		fHistBinCheckXi[identifier] = new TH1D(Form("fHistBinCheckXi_%0d",identifier),Form("fHistBinCheckXi_%0d",identifier), 500, -20, 20);
		fHistBinCheckXi[identifier]->GetXaxis()->SetTitle("#xi^{#gamma-h}");
		fHistBinCheckXi[identifier]->GetYaxis()->SetTitle("Entries");
		fOutput->Add(fHistBinCheckXi[identifier]);

		fHistBinCheckEvtPl[identifier] = new TH1D(Form("fHistBinCheckEvtPl_%0d",identifier),Form("fHistBinCheckEvtPl_%0d",identifier), 182, -182, 182);
		fHistBinCheckEvtPl[identifier]->GetXaxis()->SetTitle("#Delta#varphi^{#gamma-EvtPl}");
		fHistBinCheckEvtPl[identifier]->GetYaxis()->SetTitle("Entries");
		fOutput->Add(fHistBinCheckEvtPl[identifier]);

		fHistBinCheckEvtPl2[identifier]= new TH1D(Form("fHistBinCheckEvtPl2_%0d",identifier),Form("fHistBinCheckEvtPl2_%0d",identifier), 47, -2, 92);
		fHistBinCheckEvtPl2[identifier]->GetXaxis()->SetTitle("#Delta#varphi^{#gamma-EvtPl}");
		fHistBinCheckEvtPl2[identifier]->GetYaxis()->SetTitle("Entries");
		fOutput->Add(fHistBinCheckEvtPl2[identifier]);
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Special QA histograms (also to get more info what is going on in mixed event for trigger data)
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	for(Int_t identifier=0;identifier<kNIdentifier+3;identifier++)
	{
		//..geometrical hit distribution of clusters
		fHistDEtaDPhiGammaQA[identifier] = new TH2F(Form("fHistDEtaDPhiGammaQA%d_Id%d",0,identifier),Form("fHistDEtaDPhiGammaQA%d_Id%d",0,identifier),142,-0.71,0.71,311,75.9,331.1);
		fHistDEtaDPhiGammaQA[identifier]->GetXaxis()->SetTitle("#eta^{#gamma}");
		fHistDEtaDPhiGammaQA[identifier]->GetYaxis()->SetTitle("#varphi^{#gamma}");
		fOutputListQA->Add(fHistDEtaDPhiGammaQA[identifier]);

		//..geometrical hit distribution of tracks
		fHistDEtaDPhiTrackQA[identifier] = new TH2F(Form("fHistDEtaDPhiTrackQA%d_Id%d",0,identifier),Form("fHistDEtaDPhiTrackQA%d_Id%d",0,identifier),182,-0.91,0.91,360,0,360);
		fHistDEtaDPhiTrackQA[identifier]->GetXaxis()->SetTitle("#eta^{hadron}");
		fHistDEtaDPhiTrackQA[identifier]->GetYaxis()->SetTitle("#varphi^{hadron}");
		fOutputListQA->Add(fHistDEtaDPhiTrackQA[identifier]);

		//..Time information
		fHistClusterTime[identifier] = new TH2F(Form("fHistClusterTime%d_Id%d",0,identifier),Form("fHistClusterTime%d_Id%d",0,identifier),2000,-100,100,200,0,40);
		fHistClusterTime[identifier]->GetXaxis()->SetTitle("time [ns]");
		fHistClusterTime[identifier]->GetYaxis()->SetTitle("pT");
		fOutputListQA->Add(fHistClusterTime[identifier]);
	}
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Create histograms with the Efficiency model used in this analysis
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Int_t fNPtHistBins = 52;
	Double_t* fPtHistBins = new Double_t[fNPtHistBins+1];
	GenerateFixedBinArray(6, 0, 0.3, fPtHistBins);
	GenerateFixedBinArray(7, 0.3, 1, fPtHistBins+6);
	GenerateFixedBinArray(10, 1, 3, fPtHistBins+13);
	GenerateFixedBinArray(14, 3, 10, fPtHistBins+23);
	GenerateFixedBinArray(15, 10, 30, fPtHistBins+37);
	fEffCorrectionCheck      = new TH2*[4];

	for(Int_t cent=0;cent<4;cent++)
	{
		fEffCorrectionCheck[cent]= new TH2F(Form("fEfficiencyModel_%i",cent),Form("fEfficiencyModel_%i",cent),fNPtHistBins,fPtHistBins,100,-1,1);
//		fEffCorrectionCheck[cent]= new TH2F(Form("fEfficiencyModel_%i",cent),Form("fEfficiencyModel_%i",cent),30,0,30,100,-1,1);
		fEffCorrectionCheck[cent]->GetXaxis()->SetTitle("p_T");
		fEffCorrectionCheck[cent]->GetYaxis()->SetTitle("#eta");
		fOutputListQA->Add(fEffCorrectionCheck[cent]);

		//..Fill them imediatley with the set functions
		//..loop over eta-pt combination and evaluate the efficiency at this point
		//..retrieve the maximum of the eta function
		for(Int_t i=0;i<4;i++)
		{
			if(funcpEta_left[cent]->GetParameter(0)!=0 && funcpEta_right[cent]->GetParameter(0)!=0)
			{
				fscaleEta[i]  = funcpEta_left[i]->GetMaximum();
				if(funcpEta_right[i]->GetMaximum()>fscaleEta[i])fscaleEta[i]  = funcpEta_right[i]->GetMaximum();
			}
			else fscaleEta[i]=1;
		}

		Double_t efficiencyPt=1;
		Double_t DetectionEff;
		Double_t pT,eta;
		Double_t nYBins = fEffCorrectionCheck[0]->GetYaxis()->GetNbins();
		for(Int_t j=1;j<fNPtHistBins+1;j++)
		{
			pT = fEffCorrectionCheck[0]->GetXaxis()->GetBinCenter(j);
			efficiencyPt=1;
			if(pT<=3.5) efficiencyPt=funcpT_low[cent]->Eval(pT);
			else        efficiencyPt=funcpT_high[cent]->Eval(pT);

			for(Int_t k=1;k<nYBins+1;k++)
			{
				eta = fEffCorrectionCheck[0]->GetYaxis()->GetBinCenter(k);
				DetectionEff=efficiencyPt;
				//..only if there is an eta dependency given
				if(funcpEta_left[cent]->GetParameter(0)!=0 && funcpEta_right[cent]->GetParameter(0)!=0)
				{
					if(TMath::Abs(eta) < 0.9) { //functions have singularity at .91
						if(eta<=-0.04)DetectionEff*=funcpEta_left[cent]->Eval(eta,0,0);
						else          DetectionEff*=funcpEta_right[cent]->Eval(eta,0,0);
					}
				}
				if(fCorrectEff==1)fEffCorrectionCheck[cent]->SetBinContent(j,k,DetectionEff/fscaleEta[cent]);
				if(fCorrectEff==0)fEffCorrectionCheck[cent]->SetBinContent(j,k,1);
			}
		}
	}
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Michael's Special Histograms
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Int_t nEventPlaneBins = 128; // For histograms spanning delta psi_{EP,RP}
  Double_t fEventPlaneMin = 0;
  Double_t fEventPlaneMax = TMath::Pi()/2.;
  Double_t fEventPlaneBinArray[nEventPlaneBins+1];
  GenerateFixedBinArray(nEventPlaneBins,fEventPlaneMin,fEventPlaneMax,fEventPlaneBinArray);

  Double_t fEventPlane3Max = TMath::Pi()/3.;
  Double_t fEventPlane3BinArray[nEventPlaneBins+1];
  GenerateFixedBinArray(nEventPlaneBins,fEventPlaneMin,fEventPlane3Max,fEventPlane3BinArray);

  Double_t fEventPlane4Max = TMath::Pi()/4.;
  Double_t fEventPlane4BinArray[nEventPlaneBins+1];
  GenerateFixedBinArray(nEventPlaneBins,fEventPlaneMin,fEventPlane4Max,fEventPlane4BinArray);



  // pt bins for tracks
  Int_t nTrackPtBins = 200-1; // Bin Size = 150 MeV/c
  Double_t fTrackPtMin = 0.15; // Min Pt Cut
  Double_t fTrackPtMax = 30;
  Double_t fTrackPtArray[nTrackPtBins+1];
  GenerateFixedBinArray(nTrackPtBins,fTrackPtMin,fTrackPtMax,fTrackPtArray);

  fHistTrackPsiEPPtCent = new TH3F("fHistTrackPsiEPPtCent","Track #Delta#Psi_{EP};#Delta#Psi_{EP};p_{T} (GeV/c);Cent (%)",nEventPlaneBins,fEventPlaneBinArray,nTrackPtBins,fTrackPtArray,nCentHistBins,centBinArray);
  fOutput->Add(fHistTrackPsiEPPtCent);

  fHistTrackPsiEP3PtCent = new TH3F("fHistTrackPsiEP3PtCent","Track #Delta#Psi_{EP,3};#Delta#Psi_{EP,3};p_{T} (GeV/c);Cent (%)",nEventPlaneBins,fEventPlane3BinArray,nTrackPtBins,fTrackPtArray,nCentHistBins,centBinArray);
  fOutput->Add(fHistTrackPsiEP3PtCent);

  fHistTrackPsiEP4PtCent = new TH3F("fHistTrackPsiEP4PtCent","Track #Delta#Psi_{EP,4};#Delta#Psi_{EP,4};p_{T} (GeV/c);Cent (%)",nEventPlaneBins,fEventPlane4BinArray,nTrackPtBins,fTrackPtArray,nCentHistBins,centBinArray);
  fOutput->Add(fHistTrackPsiEP4PtCent);

  fHistTrackPsiRPPtCent = new TH3F("fHistTrackPsiRPPtCent","Track #Delta#Psi_{RP};#Delta#Psi_{RP};p_{T} (GeV/c);Cent (%)",nEventPlaneBins,fEventPlaneBinArray,nTrackPtBins,fTrackPtArray,nCentHistBins,centBinArray);
  fOutput->Add(fHistTrackPsiRPPtCent);

  int nClusEtaBins = 2;
  int nClusPhiBins = 18;
  fAccClusEtaPhi  = new TH2F("fAccClusEtaPhi","Accepted Cluster;#eta;#phi",nClusEtaBins,-0.7,0.7,nClusPhiBins,0.,2*TMath::Pi());
  fOutput->Add(fAccClusEtaPhi);
  fAccClusEtaPhiZvtx  = new TH3F("fAccClusEtaPhiZvtx","Accepted Cluster;#eta;#phi;z_{vtx}",nClusEtaBins,-0.7,0.7,nClusPhiBins,0.,2*TMath::Pi(),10,-10.,10.);
  fOutput->Add(fAccClusEtaPhiZvtx);

	if ( fGammaOrPi0>0 ) {  // Don't necessarily need this for Gamma analysis
		fClusEnergy = new TH1F("ClusEnergy","Cluster Energy",1000,0,50);
		fClusEnergy->GetXaxis()->SetTitle("E (GeV)");
		fOutput->Add(fClusEnergy);
		
		Int_t nMassBinsAccRej = 3000;
		Double_t binEdgesMassAccRej[nMassBinsAccRej+1];
		GenerateFixedBinArray(nMassBinsAccRej,0.,0.75,binEdgesMassAccRej);

    Int_t nBinsPtForEP = 30;
    Double_t binsPtForEP[nBinsPtForEP+1];
    GenerateFixedBinArray(nBinsPtForEP,0.,30.,binsPtForEP);

		//fMassPionAcc = new TH1F("fMassPtPionAcc","Accepted Pi0 Candidates;M_{#gamma#gamma} (GeV/c^2);p_{T} (GeV/c)",3000,0,0.75,250,0,50);
		//fOutput->Add(fMassPionRej);
		fMassPionRej = new TH1F("fMassPionRej","Rejected Pi0 Candidates;M_{#gamma#gamma} (GeV/c^2)",3000,0,0.75);
		fOutput->Add(fMassPionRej);

    fPtEPAnglePionAcc = new TH2F("PtEPAnglePionAcc","PtEPAnglePionAcc;#Delta#Psi_{EP}",nEventPlaneBins,fEventPlaneMin,fEventPlaneMax,60,0,30);
    fOutput->Add(fPtEPAnglePionAcc);

    fPtEPAnglePionAccCent = new TH3F("PtEPAnglePionAccCent","PtEPAnglePionAccCent;#Delta#Psi_{EP};p_{T} (GeV/c);Cent (%)",nEventPlaneBins,fEventPlaneBinArray,nBinsPtForEP,binsPtForEP,nCentHistBins,centBinArray);
    fOutput->Add(fPtEPAnglePionAccCent);

    fPtEPAngleMCPion = new TH2F("PtEPAngleMCPion","PtEPAngleMCPion;#Delta#Psi_{EP}",nEventPlaneBins,fEventPlaneMin,fEventPlaneMax,60,0,30);
    fOutput->Add(fPtEPAngleMCPion);
    fPtEPAngleTrueRecMCPion = new TH2F("PtEPAngleTrueRecMCPion","PtEPAngleTrueRecMCPion;#Delta#Psi_{EP}",nEventPlaneBins,fEventPlaneMin,fEventPlaneMax,60,0,30);
    fOutput->Add(fPtEPAngleTrueRecMCPion);


    fPtEP3AnglePionAcc = new TH2F("PtEP3AnglePionAcc","PtEP3AnglePionAcc;#Delta#Psi_{EP,3}",nEventPlaneBins,fEventPlaneMin,fEventPlane3Max,60,0,30);
    fOutput->Add(fPtEP3AnglePionAcc);

    fPtEP3AnglePionAccCent = new TH3F("PtEP3AnglePionAccCent","PtEP3AnglePionAccCent;#Delta#Psi_{EP,3};p_{T} (GeV/c);Cent (%)",nEventPlaneBins,fEventPlane3BinArray,nBinsPtForEP,binsPtForEP,nCentHistBins,centBinArray);
    fOutput->Add(fPtEP3AnglePionAccCent);

    fPtEP3AngleMCPion = new TH2F("PtEP3AngleMCPion","PtEP3AngleMCPion;#Delta#Psi_{EP,3}",nEventPlaneBins,fEventPlaneMin,fEventPlane3Max,60,0,30);
    fOutput->Add(fPtEP3AngleMCPion);
    fPtEP3AngleTrueRecMCPion = new TH2F("PtEP3AngleTrueRecMCPion","PtEP3AngleTrueRecMCPion;#Delta#Psi_{EP,3}",nEventPlaneBins,fEventPlaneMin,fEventPlane3Max,60,0,30);
    fOutput->Add(fPtEP3AngleTrueRecMCPion);


    fPtEP4AnglePionAcc = new TH2F("PtEP4AnglePionAcc","PtEP4AnglePionAcc;#Delta#Psi_{EP,4}",nEventPlaneBins,fEventPlaneMin,fEventPlane4Max,60,0,30);
    fOutput->Add(fPtEP4AnglePionAcc);

    fPtEP4AnglePionAccCent = new TH3F("PtEP4AnglePionAccCent","PtEP4AnglePionAccCent;#Delta#Psi_{EP,4};p_{T} (GeV/c);Cent (%)",nEventPlaneBins,fEventPlane4BinArray,nBinsPtForEP,binsPtForEP,nCentHistBins,centBinArray);
    fOutput->Add(fPtEP4AnglePionAccCent);

    fPtEP4AngleMCPion = new TH2F("PtEP4AngleMCPion","PtEP4AngleMCPion;#Delta#Psi_{EP,4}",nEventPlaneBins,fEventPlaneMin,fEventPlane4Max,60,0,30);
    fOutput->Add(fPtEP4AngleMCPion);
    fPtEP4AngleTrueRecMCPion = new TH2F("PtEP4AngleTrueRecMCPion","PtEP4AngleTrueRecMCPion;#Delta#Psi_{EP,4}",nEventPlaneBins,fEventPlaneMin,fEventPlane4Max,60,0,30);
    fOutput->Add(fPtEP4AngleTrueRecMCPion);



    fMCReactionPlane = new TH1F("MCReactionPlane","Reaction Plane (MC Truth)",256,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fMCReactionPlane);
    fPtRPAnglePionAcc = new TH2F("PtRPAnglePionAcc","Accepted Pion vs Reaction Plane;#Delta#Psi_{RP}",nEventPlaneBins,fEventPlaneMin,fEventPlaneMax,60,0,30);
    fOutput->Add(fPtRPAnglePionAcc);
    fPtRPAngleMCPion = new TH2F("PtRPAngleMCPion","PtEPAngleMCPion;#Delta#Psi_{RP}",nEventPlaneBins,fEventPlaneMin,fEventPlaneMax,60,0,30);
    fOutput->Add(fPtRPAngleMCPion);
    fPtRPAngleTrueRecMCPion = new TH2F("PtRPAngleTrueRecMCPion","PtRPAngleTrueRecMCPion;#Delta#Psi_{RP}",nEventPlaneBins,fEventPlaneMin,fEventPlaneMax,60,0,30);
    fOutput->Add(fPtRPAngleTrueRecMCPion);

		fEtaPhiMCPion = new TH2F("fEtaPhiMCPion","MC Pi0 Coordinates (within acceptance);#eta;#phi",200,-0.9,0.9,200,0,2*TMath::Pi());
		fOutput->Add(fEtaPhiMCPion);
		fEtaPhiPionAcc = new TH2F("fEtaPhiPionAcc","Accepted Pi0 Coordinates;#eta;#phi",200,-0.9,0.9,200,0,2*TMath::Pi());
		fOutput->Add(fEtaPhiPionAcc);

		fMassPtPionAcc = new TH2F("fMassPtPionAcc","Accepted Pi0 Candidates;M_{#gamma#gamma} (GeV/c^2);p_{T} (GeV/c)",3000,0,0.75,250,0,50);
		fOutput->Add(fMassPtPionAcc);
		fMassPtPionRej = new TH2F("fMassPtPionRej","Rejected Pi0 Candidates;M_{#gamma#gamma} (GeV/c^2);p_{T} (GeV/c)",3000,0,0.75,250,0,50);
		fOutput->Add(fMassPtPionRej);

		fMassPtCentPionAcc = new TH3F("fMassPtCentPionAcc","Accepted Pi0 Candidates;M_{#gamma#gamma} (GeV/c^2);p_{T} (GeV/c); Cent (%)",nMassBinsAccRej,binEdgesMassAccRej,kNoGammaBins,fArray_G_Bins,nCentHistBins,centBinArray);
		fOutput->Add(fMassPtCentPionAcc);
		fMassPtCentPionRej = new TH3F("fMassPtCentPionRej","Rejected Pi0 Candidates;M_{#gamma#gamma} (GeV/c^2);p_{T} (GeV/c); Cent (%)",nMassBinsAccRej,binEdgesMassAccRej,kNoGammaBins,fArray_G_Bins,nCentHistBins,centBinArray);
		fOutput->Add(fMassPtCentPionRej);

	}

  fHistEOverPvE = new TH2F("fHistEOverPvE","fHistEOverPvE;E_{#gamma} (GeV);E_{#gamma}/p_{Track}", 400, 0, 40, 400, 0, 40); 
  fOutput->Add(fHistEOverPvE);
  fHistPOverEvE = new TH2F("fHistPOverEvE","fHistPOverEvE;E_{#gamma} (GeV);p_{Track}/E_{#gamma}", 400, 0, 40, 200, 0, 20); 
  fOutput->Add(fHistPOverEvE);
	fMatchDeltaEtaTrackPt = new TH2F("fMatchDeltaEtaTrackPt","fMatchDeltaEtaTrackPt;p_{T}^{track} (GeV/c);#Delta#eta",200,0,20,200,-0.1,0.1);
	fOutput->Add(fMatchDeltaEtaTrackPt);
	fMatchDeltaPhiTrackPt = new TH2F("fMatchDeltaPhiTrackPt","fMatchDeltaPhiTrackPt;p_{T}^{track} (GeV/c);#Delta#phi",200,0,20,200,-0.1,0.1);
	fOutput->Add(fMatchDeltaPhiTrackPt);
	fMatchCondDeltaEtaTrackPt = new TH2F("fMatchCondDeltaEtaTrackPt","fMatchCondDeltaEtaTrackPt;p_{T}^{track} (GeV/c);#Delta#eta",200,0,20,200,-0.1,0.1);
	fOutput->Add(fMatchCondDeltaEtaTrackPt);
	fMatchCondDeltaPhiTrackPt = new TH2F("fMatchCondDeltaPhiTrackPt","fMatchCondDeltaPhiTrackPt;p_{T}^{track} (GeV/c);#Delta#phi",200,0,20,200,-0.1,0.1);
	fOutput->Add(fMatchCondDeltaPhiTrackPt);

  fClusterEnergyMatchedTracks = new TH2F("ClusterEnergyMatchedTracks","ClusterEnergyMatchedTracks;E_{clus} (GeV);# Matched Tracks",400,0,40,20,0,20);
  fOutput->Add(fClusterEnergyMatchedTracks);

  Double_t fMCPi0Bins[20+1];
  GenerateFixedBinArray(20,3,23,fMCPi0Bins);
  Double_t fMCPi0EtaBins[7+1];
  GenerateFixedBinArray(7,-0.7,0.7,fMCPi0EtaBins);
  if (fIsMC) {
    fHistMCPi0_PtEtaMult = new TH3F("fHistMCPi0_PtEtaMult","fHistMCPi0_PtEtaMult;p_{T}^{#pi^{0}};#eta;EMCAL Cluster Mult.",20,fMCPi0Bins,7,fMCPi0EtaBins,kNEMCalMultBins,fMixBEMCalMult->GetXbins()->GetArray());
    fOutput->Add(fHistMCPi0_PtEtaMult);
  }

  // Profiles for calculating Event Plane Resolution
  // 2nd Order Event Plane
  fEPAngleV0M = new TH1F("EPAngleV0M","EPAngleV0M;#psi_{2}^{V0M}",270,-TMath::Pi(),2*TMath::Pi());
  fEPAngleTPCA = new TH1F("EPAngleTPCA","EPAngleTPCA;#psi_{2}^{TPCA}",270,-TMath::Pi(),2*TMath::Pi());
  fEPAngleTPCC = new TH1F("EPAngleTPCC","EPAngleTPCC;#psi_{2}^{TPCC}",270,-TMath::Pi(),2*TMath::Pi());
  fOutput->Add(fEPAngleV0M);
  fOutput->Add(fEPAngleTPCA);
  fOutput->Add(fEPAngleTPCC);
  // 3rd Order Event Plane
  fEP3AngleV0M = new TH1F("EP3AngleV0M","EP3AngleV0M;#psi_{3}^{V0M}",270,-TMath::Pi(),2*TMath::Pi());
  fEP3AngleTPCA = new TH1F("EP3AngleTPCA","EP3AngleTPCA;#psi_{3}^{TPCA}",270,-TMath::Pi(),2*TMath::Pi());
  fEP3AngleTPCC = new TH1F("EP3AngleTPCC","EP3AngleTPCC;#psi_{3}^{TPCC}",270,-TMath::Pi(),2*TMath::Pi());
  fOutput->Add(fEP3AngleV0M);
  fOutput->Add(fEP3AngleTPCA);
  fOutput->Add(fEP3AngleTPCC);
  // 4th Order Event Plane
  fEP4AngleV0M = new TH1F("EP4AngleV0M","EP4AngleV0M;#psi_{4}^{V0M}",270,-TMath::Pi(),2*TMath::Pi());
  fEP4AngleTPCA = new TH1F("EP4AngleTPCA","EP4AngleTPCA;#psi_{4}^{TPCA}",270,-TMath::Pi(),2*TMath::Pi());
  fEP4AngleTPCC = new TH1F("EP4AngleTPCC","EP4AngleTPCC;#psi_{4}^{TPCC}",270,-TMath::Pi(),2*TMath::Pi());
  fOutput->Add(fEP4AngleV0M);
  fOutput->Add(fEP4AngleTPCA);
  fOutput->Add(fEP4AngleTPCC);


  fEPR_CosD1 = new TProfile2D*[kNumEPROrders];
  fEPR_CosD2 = new TProfile2D*[kNumEPROrders];
  fEPR_CosD3 = new TProfile2D*[kNumEPROrders];
  TString sEPRName = "EPR_CosD%d_N%d";
  TString sEPRTitle = "<Cos(%d[#Delta#Psi_{%d,2}])>;z_{vtx} (cm);Cent";

  for (Int_t iOrder = 0; iOrder < kNumEPROrders; iOrder++) {
    fEPR_CosD1[iOrder] = new TProfile2D(Form(sEPRName.Data(),1,iOrder+1),Form(sEPRTitle.Data(),iOrder+1,1),kNvertBins,fArrayNVertBins,4,centBinArray);
    fEPR_CosD2[iOrder] = new TProfile2D(Form(sEPRName.Data(),2,iOrder+1),Form(sEPRTitle.Data(),iOrder+1,2),kNvertBins,fArrayNVertBins,4,centBinArray);
    fEPR_CosD3[iOrder] = new TProfile2D(Form(sEPRName.Data(),3,iOrder+1),Form(sEPRTitle.Data(),iOrder+1,3),kNvertBins,fArrayNVertBins,4,centBinArray);

    fOutput->Add(fEPR_CosD1[iOrder]);
    fOutput->Add(fEPR_CosD2[iOrder]);
    fOutput->Add(fEPR_CosD3[iOrder]);
  }

  fEP3R_CosD1 = new TProfile2D*[kNumEPROrders];
  fEP3R_CosD2 = new TProfile2D*[kNumEPROrders];
  fEP3R_CosD3 = new TProfile2D*[kNumEPROrders];
  sEPRName = "EP3R_CosD%d_N%d";
  sEPRTitle = "<Cos(%d[#Delta#Psi_{%d,3}])>;z_{vtx} (cm);Cent";

  for (Int_t iOrder = 0; iOrder < kNumEPROrders; iOrder++) {
    fEP3R_CosD1[iOrder] = new TProfile2D(Form(sEPRName.Data(),1,iOrder+1),Form(sEPRTitle.Data(),iOrder+1,1),kNvertBins,fArrayNVertBins,4,centBinArray);
    fEP3R_CosD2[iOrder] = new TProfile2D(Form(sEPRName.Data(),2,iOrder+1),Form(sEPRTitle.Data(),iOrder+1,2),kNvertBins,fArrayNVertBins,4,centBinArray);
    fEP3R_CosD3[iOrder] = new TProfile2D(Form(sEPRName.Data(),3,iOrder+1),Form(sEPRTitle.Data(),iOrder+1,3),kNvertBins,fArrayNVertBins,4,centBinArray);

    fOutput->Add(fEP3R_CosD1[iOrder]);
    fOutput->Add(fEP3R_CosD2[iOrder]);
    fOutput->Add(fEP3R_CosD3[iOrder]);
  }

  fEP4R_CosD1 = new TProfile2D*[kNumEPROrders];
  fEP4R_CosD2 = new TProfile2D*[kNumEPROrders];
  fEP4R_CosD3 = new TProfile2D*[kNumEPROrders];
  sEPRName = "EP4R_CosD%d_N%d";
  sEPRTitle = "<Cos(%d[#Delta#Psi_{%d,4}])>;z_{vtx} (cm);Cent";

  for (Int_t iOrder = 0; iOrder < kNumEPROrders; iOrder++) {
    fEP4R_CosD1[iOrder] = new TProfile2D(Form(sEPRName.Data(),1,iOrder+1),Form(sEPRTitle.Data(),iOrder+1,1),kNvertBins,fArrayNVertBins,4,centBinArray);
    fEP4R_CosD2[iOrder] = new TProfile2D(Form(sEPRName.Data(),2,iOrder+1),Form(sEPRTitle.Data(),iOrder+1,2),kNvertBins,fArrayNVertBins,4,centBinArray);
    fEP4R_CosD3[iOrder] = new TProfile2D(Form(sEPRName.Data(),3,iOrder+1),Form(sEPRTitle.Data(),iOrder+1,3),kNvertBins,fArrayNVertBins,4,centBinArray);

    fOutput->Add(fEP4R_CosD1[iOrder]);
    fOutput->Add(fEP4R_CosD2[iOrder]);
    fOutput->Add(fEP4R_CosD3[iOrder]);
  }



	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //   Tyler's Special Histograms
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fHistClusPairInvarMasspT= new TH2F("fHistClusPairInvarMasspT","fHistClusPairInvarMasspT", 3000, 0, 20.0, 250, 0, 50);
    fHistClusPairInvarMasspT->GetXaxis()->SetTitle("M_{#gamma#gamma}");
	fHistClusPairInvarMasspT->GetYaxis()->SetTitle("p_{T}_{#pi^{0}}");
	fOutput->Add(fHistClusPairInvarMasspT);

	fMAngle= new TH2F("fMAngle","fMAngle", 2000, 0, 3, 500, 0, 3.14159);
	fMAngle->GetXaxis()->SetTitle("M_{#gamma#gamma}");
	fMAngle->GetYaxis()->SetTitle("Opening Angle [rad]");
	fOutput->Add(fMAngle);

	fPtAngle= new TH2F("fPtAngle","fPtAngle", 250, 0, 50, 500, 0, 3.14159);
	fPtAngle->GetXaxis()->SetTitle("p_{T}_{#pi^{0}}");
	fPtAngle->GetYaxis()->SetTitle("Opening Angle [rad]");
	fOutput->Add(fPtAngle);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//  Eliane's Special Histograms
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	fHistEvsPt= new TH2F("fHistEvsPt","fHistEvsPt", nbins[0], min[0], max[0], 250, 0, 50);
	fHistEvsPt->GetXaxis()->SetTitle("p_{T}_{#gamma}");
	fHistEvsPt->GetYaxis()->SetTitle("E_{#gamma}");
	fOutput->Add(fHistEvsPt);

	//test!!
	fHistPi0 = new TH1F(Form("fHistPi0_%0d",1),Form("fHistPi0_%0d",1), 500, 0, 0.5);
	fHistPi0->GetXaxis()->SetTitle("M_{#gamma#gamma}");
	fHistPi0->GetYaxis()->SetTitle("Entries");
	fOutput->Add(fHistPi0);

	//..The END
	fOutput->Add(fOutputListQA);

	PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::InitEventMixer(Int_t MixMode)
{
	if(fDebug==1){cout<<"Inside of: AliAnalysisTaskGammaHadron::InitEventMixer()"<<endl;
		printf("Event pool parameters: fTrackDepth %d fTargetFraction: %f\n",fTrackDepth,fTargetFraction);
	}
	//--The effective pool size in events is set by trackDepth, so more
	//--low-mult events are required to maintain the threshold than
	//--high-mult events. Centrality pools are indep. of data histogram
	//--binning, no need to match.

	//..Centrality pools
	Int_t nCentBins=fMixBCent->GetNbins();
	Double_t centBins[nCentBins+1];
	centBins[0] = fMixBCent->GetBinLowEdge(1);
	for(Int_t i=1; i<=nCentBins; i++)
	{
		centBins[i] = fMixBCent->GetBinUpEdge(i);
	}

	//..Z-vertex pools
	Int_t nZvtxBins=fMixBZvtx->GetNbins();
	Double_t zvtxbin[nZvtxBins+1];
	zvtxbin[0] = fMixBZvtx->GetBinLowEdge(1);
	for(Int_t i=1; i<=nZvtxBins; i++)
	{
		zvtxbin[i] = fMixBZvtx->GetBinUpEdge(i);
	}

	//..Event plane Pools
	// using evtPlaneArray [nEvtPlaneBins]
//	const Int_t nUsedEvtPlaneBins = 3;
//	Double_t fEventPlaneArray[nUsedEvtPlaneBins+1] = {0,1,2,3}; // In Plane, MP, Out of Plane}
	// Fake event plane array:
	//const Int_t nUsedEvtPlaneBins = 1;

	Double_t fEventPlaneArray[3] = {-0.5,0.5,1.5}; //0 and 1 used for even/odd events
  Int_t nUsedEvtPlaneBins = 1;
  if (bEnableEventHashMixing) nUsedEvtPlaneBins = 2;

	//..Pt Pools
	// using fArray_G_BinsValue [5+1]
//	Int_t nUsedTriggerPtBins = 5;

	//..in case no external pool is provided create one here
	if(!fPoolMgr)
	{
		if (MixMode == 0) {
			fPoolMgr = new AliEventPoolManager(fPoolSize, fTrackDepth, nCentBins, centBins, nZvtxBins, zvtxbin);
			AliInfo("....  Pool Manager Created for Mixed Tracks ....");
		} else { //MixMode == 1
			fPoolMgr = new AliEventPoolManager(fPoolSize, fTrackDepth, nCentBins, centBins, nZvtxBins, zvtxbin, nUsedEvtPlaneBins, fEventPlaneArray, kUsedPi0TriggerPtBins, fArray_G_Bins);
			AliInfo("....  Pool Manager Created for Mixed Triggers ....");
		}
		fPoolMgr->SetTargetValues(fTrackDepth, fTargetFraction, 5);  //pool is ready at 0.1*fTrackDepth = 5000 or events =5
		//save this pool by default
	}
	else
	{
		//..lock all pools
		//..clears empty pools and sets them locked
		//..(is only possible because all save flags are ture in my case  - NASTY NASTY)
		if (fSEvMEv != 2 ) {
			fPoolMgr->ClearPools();
		}
		AliInfo("....  Pool Manager Provided From File ....");
	}

	//..Check binning of pool manager (basic dimensional check for the time being) to see whether external pool fits the here desired one??
	if( (fPoolMgr->GetNumberOfMultBins() != nCentBins) || (fPoolMgr->GetNumberOfZVtxBins() != nZvtxBins) )
	{
		AliFatal("Binning of given pool manager not compatible with binning of correlation task!");
	}
	if (MixMode == 1) {
		if( (fPoolMgr->GetNumberOfPsiBins() != nUsedEvtPlaneBins) || (fPoolMgr->GetNumberOfPtBins() != kUsedPi0TriggerPtBins) )
		{
			AliFatal("Binning of given pool manager for mixed triggers not compatible with binning of correlation task!");
		}
	}
	//if you want to save the pool:
	// If some bins of the pool should be saved, fEventPoolOutputList must be given
	// using AddEventPoolToOutput() (to increase size of fEventPoolOutputList)
	// Note that this is in principle also possible, if an external poolmanager was given
//	if(fEventPoolOutputList.size())

	if(fSavePool==1)
	{
		//? is this an option to save only specific pools instead of the full pool manager?
		//for(Int_t i = 0; i < fEventPoolOutputList.size(); i++)
		{
			/*Double_t minCent = fEventPoolOutputList[i][0];
			Double_t maxCent = fEventPoolOutputList[i][1];
			Double_t minZvtx = fEventPoolOutputList[i][2];
			Double_t maxZvtx = fEventPoolOutputList[i][3];
			Double_t minPt   = fEventPoolOutputList[i][4];
			Double_t maxPt   = fEventPoolOutputList[i][5];
            */
		    //If the pool fulfills the given criteria the saveflag is set to true
			//the flag is used in the ClearPools function to not delete the pool content
			fPoolMgr->SetSaveFlag(-1, 10000, -10000, 100000, 0, 0, -1, 10000000);

			/*
			In case you don't want to store all pools but only the ones specified above
			you have to rund at the very end of filling these lines.

			// Clear unnecessary pools before saving and locks the pool
			fPoolMgr->ClearPools();
			*/
		}
		fOutput->Add(fPoolMgr);
	}

	//..Basic checks and printing of pool properties

	fPoolMgr->Validate();
}

//________________________________________________________________________
void AliAnalysisTaskGammaHadron::InitClusMixer()
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::InitClusMixer()"<<endl;

	if (fPoolMgr) {  // If already exists, likely because event mixing for track mixing was also suggested.
		if (fSEvMEv == 1) AliFatal("Error: Event Pool Manager already exists. Do not run MixedEvent mode at the same time as (Pi0Cands + Cluster Mixing (don't use PlotQA = 1, GammaOrPi0 = 1 and Mixing = 1 at the same time)).");
    else AliFatal("Error: Event Pool Manager already exists.  Do not try to load event pools for cluster mixing.  Cluster mixing from loaded pools not implemented yet.");
	}

	//..EMCal Multiplicity bins for cluster mixing
	Int_t nEMCalMultBins=fMixBEMCalMult->GetNbins();
	Double_t emcalMultBins[nEMCalMultBins+1];
	emcalMultBins[0] = fMixBEMCalMult->GetBinLowEdge(1);
	for(Int_t i=1; i<=nEMCalMultBins; i++)
	{
		emcalMultBins[i] = fMixBEMCalMult->GetBinUpEdge(i);
	}

	//..Z-vertex bins for cluster mixing
	Int_t nClusZvtxBins=fMixBClusZvtx->GetNbins();
	Double_t zClusvtxbin[nClusZvtxBins+1];
	zClusvtxbin[0] = fMixBClusZvtx->GetBinLowEdge(1);
	for(Int_t i=1; i<=nClusZvtxBins; i++)
	{
		zClusvtxbin[i] = fMixBClusZvtx->GetBinUpEdge(i);
	}

	//Using same trackdepth, etc., as mixed event mode.
	AliInfo("....  Pool Manager Created for cluster mixing ....");
	//fPoolMgr = new AliEventPoolManager(fPoolSize,fTrackDepth,nCentBins,centBins,nZvtxBins,zvtxbin);
	fPoolMgr = new AliEventPoolManager(fPoolSize,fClusterDepth,nEMCalMultBins,emcalMultBins,nClusZvtxBins,zClusvtxbin);
	fPoolMgr->SetTargetValues(fClusterDepth,0.05,5); //pool is ready at 0.05*fClusterDepth = 500 or events =5

	// Can still add option to save event pools out.

	//..Basic checks and printing of pool properties
	fPoolMgr->Validate();
}

///
/// Saves event pool to be used later
/// This might be a possibility to increase statistic
/// in the mixed event pool
//____________________________________________________________________
void AliAnalysisTaskGammaHadron::AddEventPoolsToOutput(Double_t minCent, Double_t maxCent,  Double_t minZvtx, Double_t maxZvtx, Double_t minPt, Double_t maxPt)
{
	//..This allows you to add only specific pools and not the full pool manager to the output file
	std::vector<Double_t> binVec;
	binVec.push_back(minCent);
	binVec.push_back(maxCent);
	binVec.push_back(minZvtx);
	binVec.push_back(maxZvtx);
	binVec.push_back(minPt);
	binVec.push_back(maxPt);
	fEventPoolOutputList.push_back(binVec);
}

//________________________________________________________________________
void AliAnalysisTaskGammaHadron::ExecOnce()
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::ExecOnce()"<<endl;

	//..This function does...
	AliAnalysisTaskEmcal::ExecOnce();

}
///Overwrites the AliAnalysisTaskEmcal::IsEventSelected()
///function
///
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::IsEventSelected()
{
	//..checks: rejects DAQ incompletes, magnetic field selection, offline trigger
	//..vertex selection, SPD pile-up (if enabled), centrality cuts
	if (!fEventCuts.AcceptEvent(InputEvent()))
	{
		PostData(1, fOutput);
		return kFALSE;
	}

	//.. .. .. .. .. .. .. .. .. .. ..
	//..Start of copy part from AliAnalysisTaskEmcal
	if (!fTrigClass.IsNull())
	{
		TString fired;

		const AliAODEvent *aev = dynamic_cast<const AliAODEvent*>(InputEvent());
		if (aev)
		{
			fired = aev->GetFiredTriggerClasses();
		}
		else cout<<"Error analysis only for AODs"<<endl;

		if (!fired.Contains("-B-"))
		{
			if (fGeneralHistograms) fHistEventRejection->Fill("trigger",1);
			return kFALSE;
		}
     /* // For some wired reason gives an error in alibild when doing pull request
        // need to check that later again because it's an exact copy from "AliAnalysisTaskEmcal"
		std::unique_ptr<TObjArray> arr(fTrigClass.Tokenize("|"));
		if (!arr)
		{
			if (fGeneralHistograms) fHistEventRejection->Fill("trigger",1);
			return kFALSE;
		}
		Bool_t match = 0;
		for (Int_t i=0;i<arr->GetEntriesFast();++i)
		{
			TObject *obj = arr->At(i);
			if (!obj)
				continue;

			//Check if requested trigger was fired
			TString objStr = obj->GetName();
			if(fEMCalTriggerMode == kOverlapWithLowThreshold &&
					(objStr.Contains("J1") || objStr.Contains("J2") || objStr.Contains("G1") || objStr.Contains("G2"))) {
				// This is relevant for EMCal triggers with 2 thresholds
				// If the kOverlapWithLowThreshold was requested than the overlap between the two triggers goes with the lower threshold trigger
				TString trigType1 = "J1";
				TString trigType2 = "J2";
				if(objStr.Contains("G"))
				{
					trigType1 = "G1";
					trigType2 = "G2";
				}
				if(objStr.Contains(trigType2) && fired.Contains(trigType2.Data()))
				{ //requesting low threshold + overlap
					match = 1;
					break;
				}
				else if(objStr.Contains(trigType1) && fired.Contains(trigType1.Data()) && !fired.Contains(trigType2.Data())) { //high threshold only
					match = 1;
					break;
				}
			}
			else
			{
				// If this is not an EMCal trigger, or no particular treatment of EMCal triggers was requested,
				// simply check that the trigger was fired
				if (fired.Contains(obj->GetName()))
				{
					match = 1;
					break;
				}
			}
		}
		if (!match)
		{
			if (fGeneralHistograms) fHistEventRejection->Fill("trigger",1);
			return kFALSE;
		}*/
	}

	if (fTriggerTypeSel != kND)
	{
		if (!HasTriggerType(fTriggerTypeSel))
		{
			if (fGeneralHistograms) fHistEventRejection->Fill("trigTypeSel",1);
			return kFALSE;
		}
	}

  if ((fMinCent != -999) && (fMaxCent != -999)) {
    if (fCent<fMinCent || fCent>fMaxCent) {
//      if (fGeneralHistograms) fHistEventRejection->Fill("Cent",1); // Disabling cent bin
      return kFALSE;
    }
  }

    /*
	//.. Maybe these two as well .. .. ..
	if (fSelectPtHardBin != -999 && fSelectPtHardBin != fPtHardBin)  {
		if (fGeneralHistograms) fHistEventRejection->Fill("SelPtHardBin",1);
		return kFALSE;
	}

	// Reject filter for MC data
	if (!CheckMCOutliers()) return kFALSE;
    */
	//..End of copy part from AliAnalysisTaskEmcal
	//.. .. .. .. .. .. .. ..
	return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::Run()
{
	//..This function is called in AliAnalysisTaskEmcal::UserExec.
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::Run()"<<endl;
	//..Determine the trigger for the current event
	fCurrentEventTrigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

	if (fIsMC) {
    fMCHeader = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));
    if (!fMCHeader) AliWarning("Missing MC Header");
		fMCParticles = GetMCParticleContainer("mcparticles");
		if (!fMCParticles) AliWarning("Missing \"mcparticles\" container"); 
    else {
      fMCParticles->SetParticleEtaLimits(-0.8, 0.8); // Loose cut for all MC particles
    }



	}

	//..Get the ClusterContainer for the event
	//AliClusterContainer* clusters  = GetClusterContainer(0);  //how do I know which cells are selected
	if (!fCaloClusters)
	{
		fCaloClusters = (TClonesArray*)GetClusterContainer(0);
		//cout<<"load calo clusters"<<endl;
	}
	//..Get the emcal cells
	if (!fCaloCells)
	{
		if (fCaloCellsName.IsNull())
		{
			fCaloCells = InputEvent()->GetEMCALCells();
		}
		else
		{
			fCaloCells =  dynamic_cast<AliVCaloCells*>(InputEvent()->FindListObject(fCaloCellsName));
			if (!fCaloCells) AliError(Form("%s: Could not retrieve cells %s!", GetName(), fCaloCellsName.Data()));
		}
		//cout<<"load calo cells"<<endl;
	}

	if (!fGeom)
	{
		AliWarning(Form("%s - AliAnalysisTaskGammaHadron::Run - Geometry is not available!", GetName()));
		return kFALSE;
	}
	//..check here some basic properties of the event
	//is centrality and zvertex in range? - see if this is actually done in IsSelected in the EMCal Task

	/*if(fCurrentEventTrigger & fTriggerType)
	{
		if(fCurrentEventTrigger & fMixingEventType)cout<<"*********************contains both triggers!!"<<endl;
	}*/

	//..for same event only analyse events when there is a cluster inside
	//..and when the event has the correct trigger
	if (fSEvMEv==0 && !fCaloClusters)                         return kFALSE;
	if (fSEvMEv==0 && !(fCurrentEventTrigger & fTriggerType)) return kFALSE;

	return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::FillHistograms()
{
	//..This function is called in AliAnalysisTaskEmcal::UserExec.
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::FillHistograms()"<<endl;

  Int_t iEventHash = CalculateEventHash();
  fHistEventHash->Fill((float) iEventHash);

  // Getting corrected event plane, saving information
  LoadQnCorrectedEventPlane();
  if (fIsMC && fMCHeader) {
    fMCReactionPlaneAngle = fMCHeader->GetReactionPlaneAngle();
    fMCReactionPlane->Fill(fMCReactionPlaneAngle);
  }

  AliClusterContainer *clusters = GetClusterContainer(0);
	if (!clusters) return 0;

  clusters->SetDefaultClusterEnergy(fClusEnergyType);

	// 1. First get an event pool corresponding in mult (cent) and
	//    zvertex to the current event. Once initialized, the pool
	//    should contain nMix (reduced) events. This routine does not
	//    pre-scan the chain. The first several events of every chain
	//    will be skipped until the needed pools are filled to the
	//    specified depth. If the pool categories are not too rare, this
	//    should not be a problem. If they are rare, you could lose
	//    statistics.

	// 2. Collect the whole pool's content of tracks into one TObjArray
	//    (bgTracks), which is effectively a single background super-event.

	// 3. The reduced and bgTracks arrays must both be passed into
	//    FillCorrelations(). Also nMix should be passed in, so a weight
	//    of 1./nMix can be applied.

	//..Get pool containing tracks from other events like this one
	Double_t zVertex = fVertex[2];
	AliParticleContainer* tracks =0x0;
	tracks   = GetParticleContainer(0);

  // Make sure Cluster Acceptance array is zeroed at start of event
  for (int i = 0; i < kCLUS_BUF_SIZE; i++) fClusterAcceptanceStatus[i] = 0;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//    Mixed event section
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if(fSEvMEv==1)
	{
		AliEventPool* pool = 0x0;
		pool = fPoolMgr->GetEventPool(fCent, zVertex);
		if (!pool)
		{
			AliWarning(Form("No pool found. Centrality %f, ZVertex %f",fCent, zVertex));
			return kFALSE;
		}

        //. . . . . . . . . . . . .
		//..Start combining triggers (from fTriggerType events)
		//..with a pool filled with tracks (from fMixingEventType)
		if(pool->IsReady() && (fCurrentEventTrigger & fTriggerType))
		{
			//..get number of current events in pool
			Int_t nMix = pool->GetCurrentNEvents();

//			cout<<"number of events in pool: "<<nMix<<endl;
			for(Int_t jMix=0; jMix<nMix; jMix++)
			{
				TObjArray* bgTracks=0x0;
				bgTracks = pool->GetEvent(jMix);

				if(!bgTracks)
				{
					AliError("could not retrieve TObjArray from EventPool!");
				}
				//..Loop over clusters and fill histograms
				if(fGammaOrPi0==0) CorrelateClusterAndTrack(0,bgTracks,0,1.0/nMix);//correlate with mixed event
				else               CorrelatePi0AndTrack(0,bgTracks,0,1.0/nMix);    //correlate with mixed event
			}
		}
		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		//    Update the pool
		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		//..update pool only with tracks from event type fMixingEventType,
		//..and do NOT add tracks from GA triggered events (BIT(15))
//		if((fCurrentEventTrigger & fMixingEventType) && ((fCurrentEventTrigger & AliVEvent::kEMCEGA)==0))
		if((fCurrentEventTrigger & fMixingEventType) && ((fCurrentEventTrigger & fVetoTrigger)==0))
		{
			TObjArray* tracksClone=0x0;
			tracksClone = CloneToCreateTObjArray(tracks);

			//..if there is no track object or the pool is locked do not update
			if(tracksClone && !pool->GetLockFlag())
			{
				pool->UpdatePool(tracksClone);
			}
		}
	}
  Double_t fEventWeight = 1;
  if (fEventWeightChoice != 0) {
    fEventWeight = 1;
  }
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//    Same event section
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//..Loop over clusters and fill histograms
	//..Do this only for events that are of fTriggerType
	if(fSEvMEv==0 && fCurrentEventTrigger & fTriggerType)
	{
		if(fGammaOrPi0==0) CorrelateClusterAndTrack(tracks,0,1,fEventWeight);//correlate with same event
		else               CorrelatePi0AndTrack(tracks,0,1,fEventWeight);    //correlate with same event
    FillTrackHistograms(tracks);
	}
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//    Mixed Trigger section
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//..Do this only for events that are of fTriggerType
	if(fSEvMEv==2 && fCurrentEventTrigger & fTriggerType)
	{
		if(fGammaOrPi0==0) CorrelateClusterAndTrack(tracks,0,1,fEventWeight);//correlate with same event
		else               CorrelatePi0AndTrack(tracks,0,1,fEventWeight);    //correlate with same event
	}


	return kTRUE;
}
///
/// Clone the tracks to create an object reduced
/// in size to be filled in the event pool
///
//________________________________________________________________________
TObjArray* AliAnalysisTaskGammaHadron::CloneToCreateTObjArray(AliParticleContainer* tracks)
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::CloneToCreateTObjArray()"<<endl;
	//..clones a track list
	if(!tracks)                            return 0;
	if(tracks->GetNAcceptedParticles()==0) return 0;

	TObjArray* tracksClone = new TObjArray;
	tracksClone->SetOwner(kTRUE);

	Int_t NoOfTracksInEvent =tracks->GetNParticles();
	AliVParticle* track=0;

	for(Int_t NoTrack = 0; NoTrack < NoOfTracksInEvent; NoTrack++)
	{
		track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
		if(!track)continue; //check if the track is a good track
		//tracksClone->Add((AliVParticle*)track);  //only add accepted tracks
		tracksClone->Add(new AliPicoTrack(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0, 0, 0, 0));
	}
	if(tracksClone->GetEntries()!=tracks->GetNAcceptedParticles())cout<<"!!!!!!! Major error!!!! "<<"Accepted tracks in event: "<<tracks->GetNAcceptedParticles()<<", Tracks in TObjArray: "<<tracksClone->GetEntries()<<endl;

	return tracksClone;
}
void AliAnalysisTaskGammaHadron::FillTrackHistograms(AliParticleContainer* tracks) {
  double pi = TMath::Pi();

  // fHistTrackPsiRPPtCent fHistTrackPsiEPPtCent
	Int_t NoOfTracksInEvent =tracks->GetNParticles();
	AliVParticle* track=0;
	for(Int_t NoTrack = 0; NoTrack < NoOfTracksInEvent; NoTrack++)
	{
		track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
		if(!track)continue; //check if the track is a good track

    TLorentzVector fLocalVector;
    tracks->GetMomentumFromParticle(fLocalVector,track);

    Double_t fLocalPhi = track->Phi();
    if (fLocalPhi < 0) fLocalPhi += TMath::TwoPi(); // LocalPhi in [0,2pi]

    // Reconstructed event plane angle
    Double_t fDeltaPsiEP = abs(DeltaPhi(fLocalVector,fQnCorrEventPlaneAngle)); // DeltaPhi in [-pi,pi]
    if ((TMath::Pi() - fDeltaPsiEP) < fDeltaPsiEP) fDeltaPsiEP = TMath::Pi() - fDeltaPsiEP; //DeltaPsi in [0,pi/2]
    fHistTrackPsiEPPtCent->Fill(fDeltaPsiEP,track->Pt(),fCent);

    // Angles w.r.t. 3rd order EP
    Double_t fDeltaPsiEP3 = abs(DeltaPhi(fLocalVector,fQnCorrEventPlane3Angle));
    fDeltaPsiEP3 = fmod(fDeltaPsiEP3,2.*pi/3.);
    if (pi/3. <= fDeltaPsiEP3) fDeltaPsiEP3 = 2.*pi/3. - fDeltaPsiEP3;
    fHistTrackPsiEP3PtCent->Fill(fDeltaPsiEP3,track->Pt(),fCent);

    // Angles w.r.t. 4th order EP
    Double_t fDeltaPsiEP4 = abs(DeltaPhi(fLocalVector,fQnCorrEventPlane4Angle));
    fDeltaPsiEP4 = fmod(fDeltaPsiEP4,pi/2.);
    if (pi/4. <= fDeltaPsiEP4) fDeltaPsiEP4 = pi/2. - fDeltaPsiEP4;
    fHistTrackPsiEP4PtCent->Fill(fDeltaPsiEP4,track->Pt(),fCent);

    // MC True ReactionPlaneAngle
    Double_t fDeltaPsiRP = abs(DeltaPhi(fLocalVector,fMCReactionPlaneAngle));
    if ((TMath::Pi() - fDeltaPsiRP) < fDeltaPsiRP) fDeltaPsiRP = TMath::Pi() - fDeltaPsiRP;
    fHistTrackPsiRPPtCent->Fill(fDeltaPsiRP,track->Pt(),fCent);

  }
  return;
}
///
/// Select tracks and clusters to correlate with each other
///
//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::CorrelateClusterAndTrack(AliParticleContainer* tracks,TObjArray* bgTracksArray,Bool_t SameMix, Double_t InputWeight)
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::CorrelateClusterAndTrack()"<<endl;

	//...........................................
	//..Do cluster loop.
	AliClusterContainer* clusters  = GetClusterContainer(0);  //how do I know which cells are selected
	if (!clusters) return 0;
	Int_t NoOfClustersInEvent = clusters->GetNClusters();
	//	Double_t EffWeight_Gamma;
	Double_t EffWeight_Hadron=1.0;
	Double_t Weight=1;    //weight to normalize mixed and same event distributions individually

	AliVCluster* cluster = 0;
	AliVCluster* cluster2= 0;
	AliVParticle* trackNULL=0;

	//...........................................
	//..for mixed events normalize per events in pool
	if(SameMix==0)
	{
		Weight=InputWeight;
	}
	//..for same events normalize later by counting the entries in spec. norm histograms
	if(SameMix==1)
	{
		Weight=1;
	}
	//...........................................
	//..run the loop for filling the histograms
	Int_t GammaCounter=0;
	for(Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
	{
		cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
		if(!cluster)continue; //check if the cluster is a good cluster

		TLorentzVector CaloClusterVec;
		clusters->GetMomentum(CaloClusterVec, cluster);
		AliTLorentzVector aliCaloClusterVec = AliTLorentzVector(CaloClusterVec); //..can acess phi from
		//	EffWeight_Gamma=GetEff(aliCaloClusterVec); //commenting out until used

		//------------------------------------------------
		//..This section is for the moment to test
		//..cluster distributions without cuts
		if(SameMix==1)FillQAHistograms(0,clusters,cluster,trackNULL,Weight);

		fFiducialCellCut->SetNumberOfCellsFromEMCALBorder(0);
		if(!AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster
		if(SameMix==1)FillQAHistograms(1,clusters,cluster,trackNULL,Weight);
		//------------------------------------------------

		if(SameMix==1)
		{
			fHistEvsPt->Fill(CaloClusterVec.Pt(),CaloClusterVec.E()); //the .pt only works for gammas (E=M) for other particle this is wrong
		}

		//........................
		//..Applying Trigger Pt Cut
		if ( CaloClusterVec.Pt() >= fTriggerPtCut )
		{
			//...........................................
			//..combine gammas with same event tracks
			GammaCounter++;
			if(SameMix==1)
			{
				if(!tracks)  return 0;
				Int_t NoOfTracksInEvent =tracks->GetNParticles();
				AliVParticle* track=0;

				if(NoOfTracksInEvent!=0) FillTriggerHist(aliCaloClusterVec,0,Weight);
				Int_t trackCounter=0;
				for(Int_t NoTrack = 0; NoTrack < NoOfTracksInEvent; NoTrack++)
				{
					track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
					if(!track)continue; //check if the track is a good track
					trackCounter++;

					EffWeight_Hadron=GetTrackEff(track->Pt(),track->Eta());
          // Could use the MC status of the cluster
					FillGhHistograms(0,aliCaloClusterVec,track,0,Weight/EffWeight_Hadron);
					if(GammaCounter==1)FillQAHistograms(2,clusters,cluster,track,Weight/EffWeight_Hadron); //fill only once per track (first gamma) - good for each track
					if(trackCounter==1)FillQAHistograms(3,clusters,cluster,track,Weight/EffWeight_Hadron); //fill only once per gamma (first track) - good for gamma distr.
				}
				//...........................................
				//..double cluster loop for testing an anti pi0 cut
				for( Int_t NoCluster2 = 0; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
				{
					if(NoCluster1!=NoCluster2 && NoCluster1<NoCluster2) //..don't combine same clusters and don't combine them twice
					{
						cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);
						if(!cluster2 || !AccClusterForAna(clusters,cluster2))continue; //check if the cluster is a good cluster

						TLorentzVector CaloClusterVec2;
						TLorentzVector CaloClusterVecpi0;
						clusters->GetMomentum(CaloClusterVec2, cluster2);
						if(cluster2->GetUserDefEnergy(fClusEnergyType)>fClEnergyMin && cluster->GetUserDefEnergy(fClusEnergyType)>fClEnergyMin)
						{
							CaloClusterVecpi0=CaloClusterVec+CaloClusterVec2;
							fHistPi0->Fill(CaloClusterVecpi0.M());
							fHistClusPairInvarMasspT->Fill(CaloClusterVecpi0.M(),CaloClusterVecpi0.Pt());
						}
					}
				}
			}
			//...........................................
			//..combine gammas with mixed event tracks
			if(SameMix==0)
			{
				Int_t Nbgtrks = bgTracksArray->GetEntries();
				if(Nbgtrks!=0) FillTriggerHist(aliCaloClusterVec,0,Weight);
				for(Int_t ibg=0; ibg<Nbgtrks; ibg++)
				{
					AliPicoTrack* track = static_cast<AliPicoTrack*>(bgTracksArray->At(ibg));
					if(!track) continue;

					EffWeight_Hadron=GetTrackEff(track->Pt(),track->Eta());
					FillGhHistograms(0,aliCaloClusterVec,track,0,Weight/EffWeight_Hadron);
				}
			}
		}
		//...........................................
		//..Additional histograms to test fiducial cell cuts
		fFiducialCellCut->SetNumberOfCellsFromEMCALBorder(1);
		if(!AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster
		if(SameMix==1)FillQAHistograms(4,clusters,cluster,trackNULL,Weight);

		fFiducialCellCut->SetNumberOfCellsFromEMCALBorder(2);
		if(!AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster
		if(SameMix==1)FillQAHistograms(5,clusters,cluster,trackNULL,Weight);
	}
	return GammaCounter;
}
///
/// Select tracks and pi0/decay-gammas to correlate with each other
///
//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::CorrelatePi0AndTrack(AliParticleContainer* tracks,TObjArray* bgTracksArray,Bool_t SameMix, Double_t InputWeight)
{
  if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::CorrelatePi0AndTrack()"<<endl;

	Double_t pi = TMath::Pi();
	//...........................................
	//--Do cluster loop.
	AliClusterContainer* clusters  = GetClusterContainer(0);
	if (!clusters) return 0;
	Int_t NoOfClustersInEvent =clusters->GetNClusters();
	Int_t nAccClusters = 0;
	Int_t nAccPi0Clusters = 0;
//	Double_t Pi0Mass = 0.13487; // Hi Michael -> this center value should also be made flexible
//	Double_t Pi0Window = 0.02;  //0.03 // Hi Michael -> the width will vary with pT
//	Double_t EffWeight_Gamma;
	Double_t EffWeight_Hadron=1.0;
	Double_t Weight;    //weight to normalize mixed and same event distributions individually

	// Double_t ClusterEnergyCut = 1; // On top of cuts in GetAcceptCluster,AccClusterForAna

	AliVCluster* cluster = 0;
	AliVCluster* cluster2= 0;
	AliVParticle* trackNULL=0;

	Weight = InputWeight; // Good enough for now.

	Double_t zVertex = fVertex[2];

	// Info for Mixed Clusters
	if ( fGammaOrPi0>0  && fPlotQA==1) {
		fEMCalMultvZvtx->Fill(zVertex,NoOfClustersInEvent); //  #clusters with E > 0.3 GeV vs Z-vertez
	}

  // Info for MC Pi0 Information
	if ( fGammaOrPi0>0  && fPlotQA==1 && fIsMC) {
    // Build array of accepted MC pi0s
    // have function Int_t                       CheckAcceptanceStatus(AliTLorentzVector Vec);
    fMCPi0List = {};
  //  for (Int_t i = 0; i < fMCParticles->Get
    AliAODMCParticle   * fMCParticle = 0;

  // FIXME Get container cuts to work
   // while ((fMCParticle = fMCParticles->GetNextAcceptMCParticle())) {
    for (Int_t i = 0; i < fMCParticles->GetNParticles(); i++) {
      fMCParticle = fMCParticles->GetMCParticle(i);
      //AliInfo(Form("Investigating MC particle (code = %d) with eta phi = (%.3f,%.1f)\n",fMCParticle->GetPdgCode(),fMCParticle->Eta(),fMCParticle->Phi()));
      if (fMCParticle->GetPdgCode() != 111) continue;
      Double_t fMCEtaAbs = fabs(fMCParticle->Eta()); // assuming all cuts sym. w.r.t sign of eta
      Double_t fMCPhi = fMCParticle->Phi();
      if (fMCPhi < 0) fMCPhi += TMath::TwoPi();
      fMCPhi = fMCPhi * fRtoD; // Degrees [0,360)
      if (fMCEtaAbs > clusters->GetMaxEta()) continue; //0.7 or so cut
      // 80 < $\varphi$ < 187 degrees
      if ((fMCPhi > 80.) && (fMCPhi < 187.)) { // EMCal Main
        fMCPi0List.push_back(fMCParticle);
        continue;
      }
      //0.22< |$\eta$|<0.7, 260 < $\varphi$ < 320 degrees
      //|$\eta$|<0.7, 320 < $\varphi$ < 327 degrees
      if (fMCPhi > 260.) {
        if (fMCPhi < 320.) { // DCal 2/3
          if (fMCEtaAbs > 0.22) {
            fMCPi0List.push_back(fMCParticle);
            continue;
          }
        } else if (fMCPhi < 327.) { // DCal 1/3 SM
          fMCPi0List.push_back(fMCParticle);
          continue;
        }
      }
    }
    FillMCPi0Hists(NoOfClustersInEvent);
  }

  //...........................................
  //..Combine pi0s with tracks in the same event
	if(SameMix==1)
	{
		AliEventPool* pool = 0x0;
		//..Load Mixed Cluser pool for mixed clusters in pi0s
		if (fPlotQA && fDoClusMixing) {
			Double_t fEMCalMultiplicity = NoOfClustersInEvent;

			pool = fPoolMgr->GetEventPool(fEMCalMultiplicity, zVertex);
			if (!pool)
			{
				AliWarning(Form("No pool found. EMCal Multiplicity %f, ZVertex %f",fCent, zVertex));
//				return kFALSE;
			}
		}

    // Also building Pi0Cand array, if fPlotQA == 1
		for(Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
		{
			cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
			if(!cluster || !AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster

			fClusEnergy->Fill(cluster->GetUserDefEnergy(fClusEnergyType),Weight);
			Int_t iMCIndexClus1 = -1;
			if (fIsMC && fPlotQA) {
				iMCIndexClus1 = FindMCPartForClus(cluster);
			}

			TLorentzVector CaloClusterVec;
			clusters->GetMomentum(CaloClusterVec, cluster);
			//acc if pi0 candidate
			nAccClusters++;

			if (fPlotQA && fDoClusMixing && pool && pool->IsReady()) {
				//.. Get current number of events in pool
				Int_t nMix = pool->GetCurrentNEvents();
				//cout<<"number of events in pools: "<<nMix<<endl;

				for (Int_t jMix = 0; jMix < nMix; jMix++)
				{
					TObjArray * mixedClusters = 0;
					mixedClusters = pool->GetEvent(jMix);

					if (!mixedClusters) {
						cout<<"Could not retrieve TObjArray from EventPool!"<<endl;
						continue;
					}
					Int_t nMixedClusters = mixedClusters->GetEntries();

					for(Int_t NoCluster2 = 0; NoCluster2 < nMixedClusters; NoCluster2++)
					{
						AliVCluster * cluster2 = (AliVCluster *) mixedClusters->At(NoCluster2);
						if (!cluster2) continue;  // No need to check acc cluster (already checked when building pool)

						TLorentzVector CaloClusterVec2;
						TLorentzVector CaloClusterVecpi0;
						clusters->GetMomentum(CaloClusterVec2, cluster2);
						Double_t fMaxClusM02 = TMath::Max(cluster->GetM02(),cluster2->GetM02());

						CaloClusterVecpi0=CaloClusterVec+CaloClusterVec2;
            // Fill Pi0 Cands for Mixed Events
						FillPi0CandsHist(CaloClusterVec,CaloClusterVec2,CaloClusterVecpi0,fMaxClusM02,Weight,1);
					}
				}
			}

			//for(Int_t NoCluster2 = 0; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
			for(Int_t NoCluster2 = NoCluster1 + 1; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
			{
				//				if(NoCluster1!=NoCluster2)
				//				{
				cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);
				if(!cluster2 || !AccClusterForAna(clusters,cluster2))continue; //check if the cluster is a good cluster


				Int_t iMCIndexClus2 = -1;
				if (fIsMC && fPlotQA) {
					iMCIndexClus2 = FindMCPartForClus(cluster2);
				}

				TLorentzVector CaloClusterVec2;
				TLorentzVector CaloClusterVecpi0;

				Double_t fMaxClusM02 = TMath::Max(cluster->GetM02(),cluster2->GetM02());

				//old framework				cluster2->GetMomentum(CaloClusterVec2, fVertex);
				clusters->GetMomentum(CaloClusterVec2, cluster2);
				//	if(cluster2->E()>2 && cluster->E()>2)
				// if(cluster2->GetNonLinCorrEnergy()>ClusterEnergyCut && cluster->GetNonLinCorrEnergy()>ClusterEnergyCut)
				//	{
				CaloClusterVecpi0=CaloClusterVec+CaloClusterVec2;
				fHistPi0->Fill(CaloClusterVecpi0.M());
				if (fPlotQA) {
					FillPi0CandsHist(CaloClusterVec,CaloClusterVec2,CaloClusterVecpi0,fMaxClusM02,Weight,0,iMCIndexClus1,iMCIndexClus2);

          // Filling Pi0 Cands with position swapped clusters
          if (fDoPosSwapMixing > 0 && NoOfClustersInEvent > 2) {

            // Filling the position swap modification functions
            Double_t fMinClusEnergy = TMath::Min(cluster->GetUserDefEnergy(fClusEnergyType),cluster2->GetUserDefEnergy(fClusEnergyType));
            // May move to separate function
            Double_t ModArray[8];
            /*
            ModArray[1] = fMaxClusM02;
            ModArray[2] = fMinClusEnergy;
            ModArray[0] = TMath::Sqrt(cluster->GetNonLinCorrEnergy()/cluster2->GetNonLinCorrEnergy());
            fVDist->Fill(ModArray,Weight);
            // Redo for switch
            ModArray[0] = TMath::Sqrt(cluster2->GetNonLinCorrEnergy() / cluster->GetNonLinCorrEnergy());
            fVDist->Fill(ModArray,Weight);
            */
            Int_t iSelectClusNo = -1;
            Int_t nTimes = 1; // 1 if no looping

            if (!fRand) fRand = new TRandom3(0);
            if (fDoPosSwapMixing == 1) {
              iSelectClusNo = fRand->Integer(NoOfClustersInEvent);
            } else {
              nTimes = NoOfClustersInEvent; // number of available clusters to mix
            }
            Int_t iRandomSamples = 0;

            //AliInfo(Form("Starting PosSwap Cycle for pair %d %d",NoCluster1,NoCluster2));

            for (Int_t j = 0; j < nTimes; j++) {
              if (nTimes > 1) { // If not random. This still works if NoOfClusters = 3
                iSelectClusNo = j;
                if ((iSelectClusNo == NoCluster1) || (iSelectClusNo == NoCluster2)) continue;
              }
              if (nTimes == 1) {
                if ((iSelectClusNo == NoCluster1) || (iSelectClusNo == NoCluster2)) {
                  iSelectClusNo = fRand->Integer(NoOfClustersInEvent);
                  j--; // redo the loop
                  iRandomSamples++; // count to avoid getting stuck in inf loop
                  if (iRandomSamples >= 5) break;
                  else continue;
                }
              }

              AliVCluster * cluster3 = 0;
              cluster3 = clusters->GetAcceptCluster(iSelectClusNo);
              //AliVCluster * cluster3 = clusters->GetAcceptCluster(iSelectClusNo);
              Int_t iMCIndexClus3 = -1;

              if (cluster3 && AccClusterForAna(clusters,cluster3)) {
                if (fIsMC) {
                  iMCIndexClus3 = FindMCPartForClus(cluster3);
                }
                // Do energy swap  on cluster first, before getting momentum
                TLorentzVector CaloClusterVecSwap;
                TLorentzVector CaloClusterVecPi0Swap;
                TLorentzVector CaloClusterVec3; // The third cluster
                clusters->GetMomentum(CaloClusterVec3, cluster3);
//                clusters->GetMomentum(CaloClusterVecPi0Swap, cluster3); //recycling pi0swap
                // Fill the mass modification histograms

                Double_t f3MaxClusM02 = TMath::Max(fMaxClusM02,cluster3->GetM02());
                Double_t f3MinClusEnergy = TMath::Min(fMinClusEnergy,cluster3->GetUserDefEnergy(fClusEnergyType));

                // Theta31 vs Theta32?
                Double_t Theta31 = CaloClusterVec3.Angle(CaloClusterVec.Vect());
                Double_t Theta32 = CaloClusterVec3.Angle(CaloClusterVec2.Vect());


                if ((Theta31 < 1e-4) || (Theta32 < 1e-4)) continue; // avoiding risk of FP error.

                // 2D Scalin distributions
                if (bEnablePosSwapHists) {
                  ModArray[2] = f3MaxClusM02;
                  ModArray[3] = f3MinClusEnergy;
                  //UMatrix
  //                ModArray[0] = TMath::Sqrt((1-TMath::Cos(Theta31))/(1-TMath::Cos(Theta32)));
                  ModArray[1] = TMath::Cos(Theta31/2.) / TMath::Cos(Theta32/2.);
                  if (bLogPSMod) {
                    // Could simplify computation further with log
                    ModArray[0] = 0.5*(TMath::Log(1-TMath::Cos(Theta31)) - TMath::Log(1-TMath::Cos(Theta32)));
                    //ModArray[0] = TMath::Log(ModArray[0]);
                    ModArray[1] = TMath::Log(ModArray[1]);
                  } else {
                    ModArray[0] = TMath::Sqrt((1-TMath::Cos(Theta31))/(1-TMath::Cos(Theta32)));
                  }
                  fUScaleMatrix->Fill(ModArray,Weight);
                  if (bLogPSMod) {
                    ModArray[0] = -ModArray[0];
                    ModArray[1] = -ModArray[1];
                  } else {
                    ModArray[0] = 1./ModArray[0];
                    ModArray[1] = 1./ModArray[1];
                  }
                  fUScaleMatrix->Fill(ModArray,Weight);

                  //VMatrix
  //                ModArray[0] = TMath::Sqrt(cluster->GetNonLinCorrEnergy()/cluster2->GetNonLinCorrEnergy());
                  ModArray[1] = (cluster3->GetUserDefEnergy(fClusEnergyType) + cluster->GetUserDefEnergy(fClusEnergyType)) /
                  (cluster3->GetUserDefEnergy(fClusEnergyType) + cluster2->GetUserDefEnergy(fClusEnergyType));
                  if (bLogPSMod) {
                    ModArray[0] = 0.5 * (TMath::Log(cluster->GetUserDefEnergy(fClusEnergyType)) - TMath::Log(cluster2->GetUserDefEnergy(fClusEnergyType)));
                    //ModArray[0] = TMath::Log(ModArray[0]);
                    ModArray[1] = TMath::Log(ModArray[1]);
                  } else {
                    ModArray[0] = TMath::Sqrt(cluster->GetUserDefEnergy(fClusEnergyType)/cluster2->GetUserDefEnergy(fClusEnergyType));
                  }
                  fVScaleMatrix->Fill(ModArray,Weight);
                  if (bLogPSMod) {
                    ModArray[0] = -ModArray[0];
                    ModArray[1] = -ModArray[1];
                  } else {
                    ModArray[0] = 1./ModArray[0];
                    ModArray[1] = 1./ModArray[1];
                  }
                  fVScaleMatrix->Fill(ModArray,Weight);
                }


                // ======================================================================
                // Set swap cluster to have energy of cluster2 (swapping pos 3 with pos 2)
                // ======================================================================
                // If cluster pair rotation is enabled the pos swap is replaced with that.
                // Do the cluster pair rotation around the axis of the combined 4-vector
                // ======================================================================
                Double_t evtPlaneAngle = 0;

               // Note: have CaloClusterVecpi0=CaloClusterVec+CaloClusterVec2;
                // have v.Rotate(TMath::Pi()/4., v1); // rotation around v1
                TLorentzVector fRotatedClus1,fRotatedClus2;
                fRotatedClus1 = CaloClusterVec;
                fRotatedClus2 = CaloClusterVec2;

                // Set to false in rotated pair mode
                Bool_t fAcceptableCluster = true;

                if (bEnableClusPairRot) {
                  // randomly decide which way
                  Double_t fRotationAngle = TMath::Pi() / 2. + TMath::Pi() * fRand->Integer(2);
                  fRotatedClus1.Rotate(fRotationAngle,CaloClusterVecpi0.Vect());
                  fRotatedClus2.Rotate(fRotationAngle,CaloClusterVecpi0.Vect());
                  //AliInfo(Form("DEBUG RotPair: angle = %.1f, combined = (%.2f,%.2f) clus1 (%.2f,%.2f) -> (%.2f,%.2f), clus2 (%.2f,%.2f) -> (%.2f,%.2f)",fRotationAngle,CaloClusterVecpi0.Eta(),CaloClusterVecpi0.Phi(),CaloClusterVec.Eta(),CaloClusterVec.Phi(),fRotatedClus1.Eta(),fRotatedClus1.Phi(),CaloClusterVec2.Eta(),CaloClusterVec2.Phi(),fRotatedClus2.Eta(),fRotatedClus2.Phi()));
                  CaloClusterVecPi0Swap = fRotatedClus1 + CaloClusterVec3;
                  fAcceptableCluster = QuickCheckAccClus(fRotatedClus1);
                } else {
                  CaloClusterVecSwap = CaloClusterVec2;
                  CaloClusterVecSwap.SetPhi(CaloClusterVec3.Phi());
                  CaloClusterVecSwap.SetTheta(CaloClusterVec3.Theta());
  //                CaloClusterVecSwap.SetPhi(CaloClusterVecPi0Swap.Phi());
  //                CaloClusterVecSwap.SetTheta(CaloClusterVecPi0Swap.Theta());

                  CaloClusterVecPi0Swap = CaloClusterVec + CaloClusterVecSwap;
                }

                evtPlaneAngle=DeltaPhi(CaloClusterVecPi0Swap,fQnCorrEventPlaneAngle); //fEPV0);




                Int_t evtPlaneCategory=-1;
                Double_t angleFromAxis;
                // Calculate EP category for the PosSwapped cluster pair
                angleFromAxis=fabs(evtPlaneAngle);
                if((pi-angleFromAxis)<angleFromAxis)angleFromAxis = pi-angleFromAxis;
                if(angleFromAxis>=0 && angleFromAxis<pi/6.)           evtPlaneCategory=0;
                else if (angleFromAxis>=pi/6. && angleFromAxis<pi/3.) evtPlaneCategory=1;
                else if (angleFromAxis>=pi/3. && angleFromAxis<=pi/2.)evtPlaneCategory=2;

                // Saving the Pos Swap Mapping information
                if (bEnablePosSwapHists) {
                  ModArray[4] = f3MaxClusM02;
                  ModArray[5] = f3MinClusEnergy;
                  ModArray[6] = evtPlaneCategory;
                  // UMap
                  ModArray[0] = CaloClusterVecpi0.M();// Initial Mass
                  ModArray[1] = CaloClusterVecpi0.Pt();// Initial Pt
                  ModArray[2] = CaloClusterVecPi0Swap.M();// Final Mass
                  ModArray[3] = CaloClusterVecPi0Swap.Pt();// Final Pt
                  fPSMassPtMap->Fill(ModArray,Weight);
                }
                if (fAcceptableCluster) {
                  // for MC, fill once with energy match MC index and once with position match index (0.5 weight each time)

                  if (bEnableClusPairRot) {
                    if (fIsMC) {
                      FillPi0CandsHist(fRotatedClus1,CaloClusterVec3,CaloClusterVecPi0Swap,fMaxClusM02,Weight,2,iMCIndexClus1,iMCIndexClus3,1);
                    } else {
                      FillPi0CandsHist(fRotatedClus1,CaloClusterVec3,CaloClusterVecPi0Swap,fMaxClusM02,Weight,2,iMCIndexClus1,iMCIndexClus3);
                    }
                  } else {

                    if (fIsMC) {
                      // Keeping the MC indices of the original pair (A,B) (Energy pair conserved for MC Info)
                      // Final argument 1 means MC indices for for E_A,E_B (A location is changed)
                      FillPi0CandsHist(CaloClusterVec,CaloClusterVecSwap,CaloClusterVecPi0Swap,fMaxClusM02,0.5*Weight,2,iMCIndexClus1,iMCIndexClus2,1);
                      // Alternate, use index of 3rd cluster
                      // MC id for (A,C) (Position pair conserved)
                      // Final argument 2 means MC indices for for x_A,x_C (An energy is changed)
                      FillPi0CandsHist(CaloClusterVec,CaloClusterVecSwap,CaloClusterVecPi0Swap,fMaxClusM02,0.5*Weight,2,iMCIndexClus1,iMCIndexClus3,2);
                    } else {
                      FillPi0CandsHist(CaloClusterVec,CaloClusterVecSwap,CaloClusterVecPi0Swap,fMaxClusM02,Weight,2,iMCIndexClus1,iMCIndexClus2);
                    }
                  }
                }
                // Energy Swap Map
                if (bEnablePosSwapHists) {
                  CaloClusterVecSwap = CaloClusterVec3;
                  CaloClusterVecSwap.SetPhi(CaloClusterVec2.Phi());
                  CaloClusterVecSwap.SetTheta(CaloClusterVec2.Theta());
                  CaloClusterVecPi0Swap = CaloClusterVec + CaloClusterVecSwap;

                  ModArray[2] = CaloClusterVecPi0Swap.M();// Final Mass
                  ModArray[3] = CaloClusterVecPi0Swap.Pt();// Final Pt
                  fESMassPtMap->Fill(ModArray,Weight);
                }


                // ======================================================================
                // Now, do it again using the energy of cluster 1 (swapping pos 3 with pos 1)
                // ======================================================================
                // If doing clus pair rotation, use the 2nd cluster of the rotated pair
                // ======================================================================
                if (bEnableClusPairRot) {
                  CaloClusterVecPi0Swap = fRotatedClus2 + CaloClusterVec3;
                  if (!QuickCheckAccClus(fRotatedClus2)) continue; // save time
                } else {
                  CaloClusterVecSwap = CaloClusterVec;
                  CaloClusterVecSwap.SetPhi(CaloClusterVec3.Phi());
                  CaloClusterVecSwap.SetTheta(CaloClusterVec3.Theta());

                  CaloClusterVecPi0Swap = CaloClusterVecSwap + CaloClusterVec2;
                }



                evtPlaneAngle=DeltaPhi(CaloClusterVecPi0Swap,fQnCorrEventPlaneAngle); //fEPV0);
                evtPlaneCategory=-1;

                // Calculate EP category for the PosSwapped cluster pair
                angleFromAxis=fabs(evtPlaneAngle);
                if((pi-angleFromAxis)<angleFromAxis)angleFromAxis = pi-angleFromAxis;
                if(angleFromAxis>=0 && angleFromAxis<pi/6.)           evtPlaneCategory=0;
                else if (angleFromAxis>=pi/6. && angleFromAxis<pi/3.) evtPlaneCategory=1;
                else if (angleFromAxis>=pi/3. && angleFromAxis<=pi/2.)evtPlaneCategory=2;

                if (bEnablePosSwapHists) {
                  ModArray[6] = evtPlaneCategory;
                  // Saving the Pos Swap Mapping information
                  // VMap
                  ModArray[0] = CaloClusterVecpi0.M();// Initial Mass
                  ModArray[1] = CaloClusterVecpi0.Pt();// Initial Pt
                  ModArray[2] = CaloClusterVecPi0Swap.M();// Final Mass
                  ModArray[3] = CaloClusterVecPi0Swap.Pt();// Final Pt
                  fPSMassPtMap->Fill(ModArray,Weight);
                }
                if (bEnableClusPairRot) {
                  if (fIsMC) {
                    FillPi0CandsHist(fRotatedClus2,CaloClusterVec3,CaloClusterVecPi0Swap,fMaxClusM02,Weight,2,iMCIndexClus2,iMCIndexClus3,1);
                  } else {
                    FillPi0CandsHist(fRotatedClus2,CaloClusterVec3,CaloClusterVecPi0Swap,fMaxClusM02,Weight,2,iMCIndexClus2,iMCIndexClus3);
                  }
                } else {
                  if (fIsMC) {
                    // MC info for energy pair (A,B)
                    FillPi0CandsHist(CaloClusterVecSwap,CaloClusterVec2,CaloClusterVecPi0Swap,fMaxClusM02,0.5*Weight,2,iMCIndexClus2,iMCIndexClus1,1);
                    // Now use MC id (C,B) for Position Pair
                    FillPi0CandsHist(CaloClusterVecSwap,CaloClusterVec2,CaloClusterVecPi0Swap,fMaxClusM02,0.5*Weight,2,iMCIndexClus3,iMCIndexClus2,2);
                  } else {
                    FillPi0CandsHist(CaloClusterVecSwap,CaloClusterVec2,CaloClusterVecPi0Swap,fMaxClusM02,Weight,2,iMCIndexClus2,iMCIndexClus1);
                  }
                }

                // Energy Swap Map
                CaloClusterVecSwap = CaloClusterVec3;
                CaloClusterVecSwap.SetPhi(CaloClusterVec.Phi());
                CaloClusterVecSwap.SetTheta(CaloClusterVec.Theta());

                CaloClusterVecPi0Swap = CaloClusterVecSwap + CaloClusterVec2;

                if (bEnablePosSwapHists) {
                  ModArray[2] = CaloClusterVecPi0Swap.M();// Final Mass
                  ModArray[3] = CaloClusterVecPi0Swap.Pt();// Final Pt
                  fESMassPtMap->Fill(ModArray,Weight);
                }
              }
            }
          }
				}
				fHistClusPairInvarMasspT->Fill(CaloClusterVecpi0.M(),CaloClusterVecpi0.Pt());
				fMAngle->Fill(CaloClusterVecpi0.M(), CaloClusterVec.Angle(CaloClusterVec2.Vect()),0.5);
				fPtAngle->Fill(CaloClusterVecpi0.Pt(), CaloClusterVec.Angle(CaloClusterVec2.Vect()),0.5);
				//        if(AccClusPairForAna(cluster,cluster2,CaloClusterVecpi0)) {
				//				if((CaloClusterVecpi0.M()>=Pi0Mass-Pi0Window) && (CaloClusterVecpi0.M()<=Pi0Mass+Pi0Window)){
				//					nAccPi0Clusters++;
				//          fMassPtPionAcc->Fill(CaloClusterVecpi0.M(),CaloClusterVecpi0.Pt());
				//				}
				//        else {
				//          fMassPtPionRej->Fill(CaloClusterVecpi0.M(),CaloClusterVecpi0.Pt());
				//        }
				//	}
			}
      //AliInfo("  Done analyzing with this Event");
		}
		if (fPlotQA && fDoClusMixing)
		{
      //AliInfo("  Saving clusters to pool ...");
			TObjArray * accClusterArr = new TObjArray();
			for(Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
			{
				cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
				if(!cluster || !AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster
				//N.B.: only stored acceptable clusters in the mixed cluster pools

				AliVCluster * accClus = (AliVCluster *) cluster->Clone();
				accClusterArr->Add(accClus);
			}
			if (!pool->GetLockFlag()) {
				pool->UpdatePool(accClusterArr);
			}
     // AliInfo("       ... Done!");
		}
	}
	//...........................................
	//..for mixed events normalize per events in pool
	if(SameMix==0)
	{
		Weight=InputWeight;
	}
	//..for same events normalize later by counting the entries in spec. norm histograms
	if(SameMix==1)
	{
		Weight=1;
	}
//	vector<TLorentzVector> AcceptTriggerArray = {};
	//...........................................
	//run the loop for filling the histograms
	if (fSEvMEv == 2) { // Use Mixed Triggers from Pool
    if (fDownScaleMT != 1.0) {
      //Check if random gen. ready
      if (!fRand) fRand = new TRandom3(0);
    }
		if (!tracks) return 0;
		Int_t NoOfTracksInEvent = tracks->GetNParticles();
		AliVParticle * track = 0x0;
		// Select From Available pool.  Have fCent, zVertex.
		AliEventPool * pool = 0x0;
		// maybe iterate over pt bins?
		Double_t EventPlaneAngle = 0.;

    // Use event hash to split data set, avoid autocorrelation
    Int_t iEventHash = CalculateEventHash();
    if (bEnableEventHashMixing) {
      EventPlaneAngle = (Double_t) !(iEventHash); // use pool 0 for 1, pool 1 for 0
    }

		for (Int_t PtIndex = 0; PtIndex < kUsedPi0TriggerPtBins; PtIndex++) {
			pool = fPoolMgr->GetEventPool(fCent, zVertex,EventPlaneAngle,PtIndex);
			if (!pool) {
				AliWarning(Form("No pool found. Centrality %f, ZVertex %f, PtIndex %d",fCent, zVertex,PtIndex));
			}
			Int_t nMix = pool->GetCurrentNEvents();
			for (Int_t jMix = 0; jMix < nMix; jMix++) {
				TObjArray * mixedTriggers = 0;
				// Potential to introduce random sampling for downscaling here
        if (fDownScaleMT != 1.0) {
          if (fRand->Rndm() > fDownScaleMT) continue;
        }
				mixedTriggers = pool->GetEvent(jMix);
				if (!mixedTriggers) {
					AliWarning("Could not retrieve TObjArray from EventPool!");
					continue;
				}
				Int_t nTriggers = mixedTriggers->GetEntries();
        for (Int_t iTrack = 0; iTrack < NoOfTracksInEvent; iTrack++) {
          track = (AliVParticle *) tracks->GetAcceptParticle(iTrack);
          if (!track) continue;
          EffWeight_Hadron=GetTrackEff(track->Pt(),track->Eta());
          for (Int_t iTrigger = 0; iTrigger < nTriggers; iTrigger++) {
            AliTLorentzVector * Pi0CandVector = (AliTLorentzVector *) mixedTriggers->At(iTrigger);
						FillGhHistograms(0,*Pi0CandVector,track,0,Weight/EffWeight_Hadron);
					}
				}
			}
		}

	} else {
	for( Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
	{
		cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
		if(!cluster || !AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster

		TLorentzVector CaloClusterVec;
		clusters->GetMomentum(CaloClusterVec,cluster);
		AliTLorentzVector aliCaloClusterVec = AliTLorentzVector(CaloClusterVec); //..can acess phi from

		FillQAHistograms(0,clusters,cluster,trackNULL,Weight);

		for( Int_t NoCluster2 = NoCluster1 + 1; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
		{
			cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);
			if(!cluster2 || !AccClusterForAna(clusters,cluster2))continue; //check if the cluster is a good cluster

			TLorentzVector CaloClusterVec2;
			TLorentzVector aliCaloClusterVecpi0;
			//old framework			cluster2->GetMomentum(CaloClusterVec2, fVertex);
			clusters->GetMomentum(CaloClusterVec2,cluster2); /// Vec+=2 2.1.17
			AliTLorentzVector aliCaloClusterVec2 = AliTLorentzVector(CaloClusterVec2); //..can acess phi from

			aliCaloClusterVecpi0=aliCaloClusterVec+aliCaloClusterVec2;

			//........................
			//..Applying Trigger Pt Cut
			if ( aliCaloClusterVecpi0.Pt() < fTriggerPtCut) continue;

			//........................
			//..Applying Opening Angle Cut
			if ( aliCaloClusterVec.Angle(aliCaloClusterVec2.Vect()) < fOpeningAngleCut) continue;

			if(AccClusPairForAna(cluster,cluster2,aliCaloClusterVecpi0))
			{
        Double_t fLocalPhi = aliCaloClusterVecpi0.Phi();
        // Angle relative to 2nd order event plane
        if (fLocalPhi < 0) fLocalPhi += 2*pi;
        Double_t evtPlaneAngle=abs(DeltaPhi(aliCaloClusterVecpi0,fQnCorrEventPlaneAngle));
        if ((pi - evtPlaneAngle) < evtPlaneAngle) evtPlaneAngle = pi - evtPlaneAngle;
        fPtEPAnglePionAcc->Fill(evtPlaneAngle,aliCaloClusterVecpi0.Pt());
        fPtEPAnglePionAccCent->Fill(evtPlaneAngle,aliCaloClusterVecpi0.Pt(),fCent);

        // Angle relative to 3rd order event plane

        Double_t evtPlane3Angle=abs(DeltaPhi(aliCaloClusterVecpi0,fQnCorrEventPlane3Angle));
        evtPlane3Angle  = fmod(evtPlane3Angle,2.*pi/3); // now in range [0,2pi/3]
        if (pi/3. <= evtPlane3Angle) evtPlane3Angle = 2.*pi/3. - evtPlane3Angle;
        fPtEP3AnglePionAcc->Fill(evtPlane3Angle,aliCaloClusterVecpi0.Pt());
        fPtEP3AnglePionAccCent->Fill(evtPlane3Angle,aliCaloClusterVecpi0.Pt(),fCent);

        // Angle relative to 4th order event plane

        Double_t evtPlane4Angle=abs(DeltaPhi(aliCaloClusterVecpi0,fQnCorrEventPlane4Angle));
        evtPlane4Angle = fmod(evtPlane4Angle,pi/2.); // now in range [0,pi/2.]
        if (pi/4. <= evtPlane4Angle) evtPlane4Angle = pi/2. - evtPlane4Angle;
        fPtEP4AnglePionAcc->Fill(evtPlane4Angle,aliCaloClusterVecpi0.Pt());
        fPtEP4AnglePionAccCent->Fill(evtPlane4Angle,aliCaloClusterVecpi0.Pt(),fCent);


        Double_t reactionPlaneAngle=abs(DeltaPhi(aliCaloClusterVecpi0,fMCReactionPlaneAngle));
        if ((pi - reactionPlaneAngle) < reactionPlaneAngle) reactionPlaneAngle = pi - reactionPlaneAngle;
        fPtRPAnglePionAcc->Fill(reactionPlaneAngle,aliCaloClusterVecpi0.Pt());

				fEtaPhiPionAcc->Fill(aliCaloClusterVecpi0.Eta(),fLocalPhi);
				fMassPtPionAcc->Fill(aliCaloClusterVecpi0.M(),aliCaloClusterVecpi0.Pt());
				fMassPtCentPionAcc->Fill(aliCaloClusterVecpi0.M(),aliCaloClusterVecpi0.Pt(),fCent);
				nAccPi0Clusters++;

				if (fSaveTriggerPool) { //update an array of accepted triggers
					// TObjArray of TLorentzVectors
//					TObjArray * AcceptTriggerArray = 0x0;
					TObjArray * AcceptTriggerArray = new TObjArray(1);
					AcceptTriggerArray->SetOwner(kTRUE);
					TLorentzVector * Pi0CandVector = (TLorentzVector *) aliCaloClusterVecpi0.Clone();
					AcceptTriggerArray->Add(Pi0CandVector);

//					vector<TLorentzVector> AcceptTriggerArray = {};
//					AcceptTriggerArray.push_back(aliCaloClusterVecpi0);

//					Double_t evtPlaneAngle= DeltaPhi(aliCaloClusterVecpi0,fEPV0);
					Double_t evtPlaneCategory=0;

          Int_t iEventHash = CalculateEventHash();
          if (bEnableEventHashMixing) {
            evtPlaneCategory = (Double_t) iEventHash; // are you a 1 or a 0?
          }
/*
					Double_t angleFromAxis;
					//..fold around 0 axis
					angleFromAxis=fabs(evtPlaneAngle);
					//..fold once more arounad pi/2 axis
					if((pi-angleFromAxis)<angleFromAxis)angleFromAxis = pi-angleFromAxis;
					//.. --> now we have only one quadrant left 0 to  <pi/2
					if(angleFromAxis>=0 && angleFromAxis<pi/6.)           evtPlaneCategory=0;
					else if (angleFromAxis>=pi/6. && angleFromAxis<pi/3.) evtPlaneCategory=1;
					else if (angleFromAxis>=pi/3. && angleFromAxis<=pi/2.)evtPlaneCategory=2;*/
					Double_t Pi0Pt = aliCaloClusterVecpi0.Pt();
					Int_t Pi0PtBin = -1;
					for (Int_t k = 0; k < kUsedPi0TriggerPtBins; k++)
					{
						if ((Pi0Pt >= fArray_G_Bins[k]) && (Pi0Pt < fArray_G_Bins[k+1]))
						{
							Pi0PtBin = k;
							break;
						}
					}
					if (Pi0PtBin >= 0) {
						// centrality, z-vertex, event plane angle, pt
						AliEventPool * pool = 0x0;
						pool = fPoolMgr->GetEventPool(fCent,zVertex,evtPlaneCategory,Pi0PtBin);
						if (pool) {
							pool->UpdatePool(AcceptTriggerArray);
						}
					}
				}


        // FIXME Check the MC information Here
        Int_t iMCIndexClus1 = -1;
        Int_t iMCIndexClus2 = -1;
        Int_t iCorrMCStatus = -1;
				if (fIsMC) {
					iMCIndexClus1 = FindMCPartForClus(cluster);
					iMCIndexClus2 = FindMCPartForClus(cluster2);

          Int_t aMCTreeHeight1 = 0;
          Int_t aMCRootPartClus1 = FindMCRootPart(iMCIndexClus1,&aMCTreeHeight1);
          Int_t aMCTreeHeight2 = 0;
          Int_t aMCRootPartClus2 = FindMCRootPart(iMCIndexClus2,&aMCTreeHeight2);

          if (aMCRootPartClus1 == aMCRootPartClus2) {
            Int_t iLCA = FindMCLowComAnc(iMCIndexClus1,iMCIndexClus2); // Lowest Common Ancestor of clusters
            if (iLCA == -1) { // impossible case
              iCorrMCStatus = 0;
            } else {
              AliAODMCParticle * pLCA = fMCParticles->GetMCParticle(iLCA);
              Int_t iLCAPdg = pLCA->GetPdgCode();

              if (iLCAPdg == 111) iCorrMCStatus = 2; // pi0 -> 2 gamma
              else if (iLCAPdg == 221) { //eta
                // Check that it is eta -> 2gamma. Otherwise, consider it background.
                iCorrMCStatus = 0;
                Int_t nLCADaughters = pLCA->GetNDaughters();
                if (nLCADaughters == 2) { // 2 Gammas
                  AliAODMCParticle * pDaughter1 = fMCParticles->GetMCParticle(pLCA->GetDaughterLabel(0));
                  AliAODMCParticle * pDaughter2 = fMCParticles->GetMCParticle(pLCA->GetDaughterLabel(0)+1);
                  if (pDaughter1 && pDaughter2 && pDaughter1->GetPdgCode() == 22 && pDaughter2->GetPdgCode() == 22) {
                    iCorrMCStatus = 1; // eta -> 2gamma
                  }
                }
              }
            }
          } else {
            // True Background pair
            iCorrMCStatus = 0;
          }
				}


				//...........................................
				//..combine gammas with same event tracks
				if(SameMix==1)
				{
					//cout<<"SameMix==1"<<endl;
					if(!tracks)  return 0;
					Int_t NoOfTracksInEvent =tracks->GetNParticles();
					AliVParticle* track=0;
					FillTriggerHist(aliCaloClusterVecpi0,iCorrMCStatus,Weight);

					for(Int_t NoTrack = 0; NoTrack < NoOfTracksInEvent; NoTrack++)
					{
						track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
						if(!track)continue; //check if the track is a good track

						//..fill here eventually a pi0 four-vector instead of CaloClusterVec
						EffWeight_Hadron=GetTrackEff(track->Pt(),track->Eta());
						FillGhHistograms(0,aliCaloClusterVecpi0,track,iCorrMCStatus,Weight/EffWeight_Hadron);
					}
				}
				//...........................................
				//..combine gammas with mixed event tracks
				if(SameMix==0)
				{
					Int_t Nbgtrks = bgTracksArray->GetEntries();
					if(Nbgtrks!=0) FillTriggerHist(aliCaloClusterVecpi0,0,Weight);
					for(Int_t ibg=0; ibg<Nbgtrks; ibg++)
					{
						AliPicoTrack* track = static_cast<AliPicoTrack*>(bgTracksArray->At(ibg));
						if(!track) continue;

						//**fill here eventually a pi0 four-vector instead of CaloClusterVec
						EffWeight_Hadron=GetTrackEff(track->Pt(),track->Eta());
						FillGhHistograms(0,aliCaloClusterVecpi0,track,0,Weight/EffWeight_Hadron);
  // FIXME need to inser corrmc status
					}
				}
			}
		}
	}
	}
//	if (fSEvME == 0 && fSaveTriggerPool) { // Saving trigger information
//		pool->UpdatePool(accClusterArr);
//	}
	return nAccPi0Clusters/2;
}
///
/// Fill Pi0 Cand THnSparse with relevant info
/// To Do: add in rotation method
///
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::FillPi0CandsHist(AliTLorentzVector CaloClusterVec, AliTLorentzVector CaloClusterVec2, AliTLorentzVector CaloClusterVecPi0, Double_t fMaxClusM02, Double_t Weight, Int_t isMixed, Int_t mcIndex1, Int_t mcIndex2, Int_t PosSwapStatus)
{
	Double_t pi = TMath::Pi();
  Double_t FillWeight = Weight; // May want to modify weight

	Double_t valueArray[11];
	valueArray[0]=CaloClusterVecPi0.Pt();
	valueArray[1]=CaloClusterVecPi0.M();
	valueArray[2]=CaloClusterVec.Angle(CaloClusterVec2.Vect());

	Double_t fE1 = CaloClusterVec.E();
	Double_t fE2 = CaloClusterVec2.E();
	Double_t fAsym = (fE1+fE2 > 0.000001) ? TMath::Abs(fE2-fE1)/(fE1+fE2) : 0; //Don't divide by zero

	valueArray[3]=fMaxClusM02;
	valueArray[4]=TMath::Min(fE1,fE2);
	valueArray[5]=fAsym;

	if (isMixed == 1) valueArray[6] = 2;
  else if (isMixed == 2) valueArray[6] = 3; // cluster swapping
	else valueArray[6]=0;

  // Event Plane will be 7
  Double_t evtPlaneAngle=DeltaPhi(CaloClusterVecPi0,fQnCorrEventPlaneAngle); //fEPV0);
	Int_t evtPlaneCategory=-1;
	Double_t angleFromAxis;
	angleFromAxis=fabs(evtPlaneAngle);
	if((pi-angleFromAxis)<angleFromAxis)angleFromAxis = pi-angleFromAxis;
	if(angleFromAxis>=0 && angleFromAxis<pi/6.)           evtPlaneCategory=0;
	else if (angleFromAxis>=pi/6. && angleFromAxis<pi/3.) evtPlaneCategory=1;
	else if (angleFromAxis>=pi/3. && angleFromAxis<=pi/2.)evtPlaneCategory=2;

  valueArray[7]=evtPlaneCategory;


	// MC Status determination
	Int_t MCMatchStatus = 0; // 0 for no matches
  // Could change code here to allow MC id for Position Swap
//	if (!isMixed && fIsMC && fMCParticles) {
	if ((isMixed != 1) && fIsMC && fMCParticles) {
		if (mcIndex1 < 0 || mcIndex2 < 0) {
		}
		else if (mcIndex1 == mcIndex2) {
			 MCMatchStatus = 1; // 2 clusters from 1 MC Part
		} else {
			Int_t aMCTreeHeight1 = 0;
			Int_t aMCRootPartClus1 = FindMCRootPart(mcIndex1,&aMCTreeHeight1);
			Int_t aMCTreeHeight2 = 0;
			Int_t aMCRootPartClus2 = FindMCRootPart(mcIndex2,&aMCTreeHeight2);

			if (aMCRootPartClus1 == aMCRootPartClus2) { //The MC Parts still have a common ancestor

				Int_t iLCA = FindMCLowComAnc(mcIndex1,mcIndex2); // Lowest Common Ancestor of clusters
				if (iLCA != -1) {
					AliAODMCParticle * pLCA = fMCParticles->GetMCParticle(iLCA);
          TLorentzVector fLCAVector;
          fMCParticles->GetMomentumFromParticle(fLCAVector,pLCA);
					Int_t iLCAPdg = pLCA->GetPdgCode();

          Int_t iLCAStatus = pLCA->MCStatusCode(); // 1 for embedded (0 otherwise)

          if (fMCEmbedReweightMode > 0) {
            if (iLCAStatus == 1) { // is embedded particle
            }
            if (fMCEmbedReweightMode == 2 && iLCAPdg == 221) return; // Vetoing on embedded eta
          }


					if (iLCAPdg == 111) {
						Int_t nLCADaughters = pLCA->GetNDaughters();
						if (nLCADaughters == 2) { // 2 Gammas
              AliAODMCParticle * pDaughter1 = fMCParticles->GetMCParticle(pLCA->GetDaughterLabel(0));
              AliAODMCParticle * pDaughter2 = fMCParticles->GetMCParticle(pLCA->GetDaughterLabel(0)+1);
							if (pDaughter1 && pDaughter2 && pDaughter1->GetPdgCode() == 22 && pDaughter2->GetPdgCode() == 22) {
								MCMatchStatus = 2;

                // Reaction Plane Angle information
                Double_t fTrueAngleFromEventPlane = abs(DeltaPhi(fLCAVector,fQnCorrEventPlaneAngle));
                if ((TMath::Pi() - fTrueAngleFromEventPlane) < fTrueAngleFromEventPlane) fTrueAngleFromEventPlane = TMath::Pi() - fTrueAngleFromEventPlane;
                Double_t fTrueAngleFromReactionPlane = abs(DeltaPhi(fLCAVector,fMCReactionPlaneAngle));
                if ((TMath::Pi() - fTrueAngleFromReactionPlane) < fTrueAngleFromReactionPlane) fTrueAngleFromReactionPlane = TMath::Pi() - fTrueAngleFromReactionPlane;

                fPtEPAngleTrueRecMCPion->Fill(fTrueAngleFromEventPlane,pLCA->Pt());
                fPtRPAngleTrueRecMCPion->Fill(fTrueAngleFromReactionPlane,pLCA->Pt());

                if (PosSwapStatus == 1) {
                  MCMatchStatus = 10;
                } else if (PosSwapStatus == 2) {
                  MCMatchStatus = 11;
                }
								// Check Reconstructed Pi0 DeltaPt and DeltaPhiDeltaEta
								Double_t fPi0_Pt = CaloClusterVecPi0.Pt();
								Double_t fDeltaPt = fPi0_Pt - pLCA->Pt();

								Double_t fDeltaPhi = DeltaPhi(CaloClusterVecPi0,pLCA->Phi());
								Double_t fDeltaEta = CaloClusterVecPi0.Eta() - pLCA->Eta();

								fHistPi0MCDPt->Fill(fPi0_Pt,fDeltaPt);
								fHistPi0MCDPhiDEta->Fill(fDeltaPhi,fDeltaEta);
							}
						} else { // other pi0, most likely gamma,e+,e-
							MCMatchStatus = 3;
						}
					} else if (iLCAPdg == 221) {
						// Categorize Etas here
						// Check daughters  2 gamma, 3 pi0, or 1pi0,1pi+,ipi-
						Int_t nLCADaughters = pLCA->GetNDaughters();
						if (nLCADaughters == 2) { // 2 Gammas
							AliAODMCParticle * pDaughter1 = fMCParticles->GetMCParticle(pLCA->GetDaughterLabel(0));
							AliAODMCParticle * pDaughter2 = fMCParticles->GetMCParticle(pLCA->GetDaughterLabel(0)+1);
							if (pDaughter1 && pDaughter2 && pDaughter1->GetPdgCode() == 22 && pDaughter2->GetPdgCode() == 22) {

								MCMatchStatus = 4;
                if (PosSwapStatus == 1) {
                  MCMatchStatus = 12;
                } else if (PosSwapStatus == 2) {
                  MCMatchStatus = 13;
                }
								// Check Reconstructed Eta DeltaPt and DeltaPhiDeltaEta
								Double_t fPi0_Pt = CaloClusterVecPi0.Pt();
								Double_t fDeltaPt = fPi0_Pt - pLCA->Pt();

								Double_t fDeltaPhi = DeltaPhi(CaloClusterVecPi0,pLCA->Phi());
								Double_t fDeltaEta = CaloClusterVecPi0.Eta() - pLCA->Eta();

								fHistEtaMCDPt->Fill(fPi0_Pt,fDeltaPt);
								fHistEtaMCDPhiDEta->Fill(fDeltaPhi,fDeltaEta);
							}
						} else if (nLCADaughters == 3) { // 3 pi0 or pi0,pi+,pi-
							AliAODMCParticle * pDaughter1 = fMCParticles->GetMCParticle(pLCA->GetDaughterLabel(0));
							AliAODMCParticle * pDaughter2 = fMCParticles->GetMCParticle(pLCA->GetDaughterLabel(0)+1);
							AliAODMCParticle * pDaughter3 = fMCParticles->GetMCParticle(pLCA->GetDaughterLabel(0)+2);
							Int_t iPdg1 = pDaughter1->GetPdgCode();
							Int_t iPdg2 = pDaughter2->GetPdgCode();
							Int_t iPdg3 = pDaughter3->GetPdgCode();
							if (iPdg1 == 111 && iPdg2 == 111 && iPdg3 == 111) MCMatchStatus = 5;
								// eta -> 3pi0
							else if (iPdg1 + iPdg2 + iPdg3 == 111 && iPdg1*iPdg2*iPdg3 == 111*211*(-211)) {
								// eta -> pi0,pi+,pi-
								MCMatchStatus = 6;
							} else if (iPdg1 + iPdg2 + iPdg3 == 22) {
								// eta -> gamma,pi+,pi-
								// eta -> gamma,e+,e-
								// eta -> gamma,mu+,mu-
								MCMatchStatus = 7;
							} else if (iPdg1 + iPdg2 + iPdg3 == 111+22+22 && iPdg1*iPdg2*iPdg3 == 111*22*22) {
								// eta-> pi0,gamma,gamma
								MCMatchStatus = 6; // similar to pi0,pi+,pi-
							} else {
								MCMatchStatus = 9; // group this with other shared ancestor (probably never happen)
							}
						}
					}
					else if (iLCAPdg == 22) MCMatchStatus = 8;
					else MCMatchStatus = 9; // FIXME check if it is a parton? Then jet.

				} else {
					MCMatchStatus = 0;
				}
			}
		}
	}
	if (fIsMC) {
		valueArray[8] = MCMatchStatus;
		if (fApplyPatchCandCut) valueArray[9] = DetermineGAPatchCand(CaloClusterVec,CaloClusterVec2);
	} else { // if no MC, then axis 7 is patch status
		if (fApplyPatchCandCut) valueArray[8] = DetermineGAPatchCand(CaloClusterVec,CaloClusterVec2);
	}


	fPi0Cands->Fill(valueArray,FillWeight);

	if (!fDoRotBkg) return;
	if (isMixed) return; // don't do rotational background if off or if this is looking at mixed cluster pairs
	// Rotational Background
	//  const Double_t fOpeningAngleCut = 0.017;

  // Can replace this with Joshua Koenig's rotated background here
  // add in a reference to the cluster container, so the second loop can be done here?
  // or it might be better to do that in the position swap code

	if (!fRand) fRand = new TRandom3(0);

	Double_t fParA = 1.3;

	for (int i = 0; i < fNRotBkgSamples; i++) {
		Double_t fEta,fPhi;
		Double_t fOpeningAngle = 0;
		Int_t nLoopTrials = 15; // avoid too many trials
		for (int j = 0; j < nLoopTrials; j++) {
			fEta = fRand->Uniform(-0.7,0.7);  // change to eta cut maybe
			// GetClusterContainer("caloClusters")->SetMinEta()
			// GetClusterContainer("caloClusters")->SetMaxEta()
			fPhi = fRand->Uniform(80,254); // pretend DCAL next to EMCAL
			if(fPhi > 187) {
				// Check PHOS hole
				if (TMath::Abs(fEta) < .22 && fPhi < 247) continue;
				fPhi+= 73;  // shift DCAL points
			}
			fPhi = fPhi * 3.141592653589793 / 180.;

			// Opening Angle Cut
			CaloClusterVec2.SetPhi(fPhi);
			CaloClusterVec2.SetTheta(2.*TMath::ATan(TMath::Exp(-fEta)));
			fOpeningAngle = CaloClusterVec.Angle(CaloClusterVec2.Vect());

			// Weighting towards lower angles (higher pT)
			Double_t fR = fRand->Rndm();
			if (fR > 0.1 + 0.9 * TMath::Exp(-TMath::Power(fOpeningAngle/fParA,3))) continue;

			break;
		}

		CaloClusterVec2.SetPhi(fPhi);
		CaloClusterVec2.SetTheta(2.*TMath::ATan(TMath::Exp(-fEta)));

		CaloClusterVecPi0=CaloClusterVec+CaloClusterVec2;

		valueArray[0]=CaloClusterVecPi0.Pt();
		valueArray[1]=CaloClusterVecPi0.M();
		valueArray[2]=fOpeningAngle;
		valueArray[3]=fMaxClusM02;
		valueArray[4]=TMath::Min(fE1,fE2);
		valueArray[5]=fAsym;
		valueArray[6]=1;

    // Event Plane Angle For Rot Bkg
    evtPlaneAngle=DeltaPhi(CaloClusterVecPi0,fQnCorrEventPlaneAngle); //fEPV0);
    evtPlaneCategory=-1;
    angleFromAxis=fabs(evtPlaneAngle);
    if((pi-angleFromAxis)<angleFromAxis)angleFromAxis = pi-angleFromAxis;
    if(angleFromAxis>=0 && angleFromAxis<pi/6.)           evtPlaneCategory=0;
    else if (angleFromAxis>=pi/6. && angleFromAxis<pi/3.) evtPlaneCategory=1;
    else if (angleFromAxis>=pi/3. && angleFromAxis<=pi/2.)evtPlaneCategory=2;

    valueArray[7]=evtPlaneCategory;

		if (fIsMC) {
			valueArray[7] = 0; // MC info is meaningless for rot bkg
			if (fApplyPatchCandCut) valueArray[8] = DetermineGAPatchCand(CaloClusterVec,CaloClusterVec2);
		} else { // if no MC, then axis 7 is patch status
			if (fApplyPatchCandCut) valueArray[7] = DetermineGAPatchCand(CaloClusterVec,CaloClusterVec2);
		}

		fPi0Cands->Fill(valueArray,Weight);
	}
}
///
/// Fill histogram with Trigger information
///
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::FillTriggerHist(AliTLorentzVector ClusterVec, Int_t CorrMCStatus, Double_t Weight)
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::FillTriggerHist()"<<endl;

	Double_t pi = TMath::Pi();

	Double_t G_PT_Value = ClusterVec.Pt();

	Double_t zVertex = fVertex[2];
	//..all from EMcal base class : fEPV0,fEPV0A,fEPV0C
	//fEPV0  = aliEP->GetEventplane("V0" ,InputEvent());
    //fEPV0A = aliEP->GetEventplane("V0A",InputEvent());
    //fEPV0C = aliEP->GetEventplane("V0C",InputEvent());

	Double_t evtPlaneAngle= DeltaPhi(ClusterVec,fQnCorrEventPlaneAngle); //fEPV0);
	Int_t evtPlaneCategory=-1;

  Double_t angleFromAxis;
  //..fold around 0 axis
  angleFromAxis=fabs(evtPlaneAngle);
  //..fold once more arounad pi/2 axis
  if((pi-angleFromAxis)<angleFromAxis)angleFromAxis = pi-angleFromAxis;
  //.. --> now we have only one quadrant left 0 to  <pi/2
  if(angleFromAxis>=0 && angleFromAxis<pi/6.)           evtPlaneCategory=0;
  else if (angleFromAxis>=pi/6. && angleFromAxis<pi/3.) evtPlaneCategory=1;
  else if (angleFromAxis>=pi/3. && angleFromAxis<=pi/2.)evtPlaneCategory=2;

  double fCorrMCStatus = CorrMCStatus;

	Double_t valueArray[5];
	valueArray[0]=G_PT_Value;
	valueArray[1]=zVertex;
	valueArray[2]=evtPlaneCategory;
	valueArray[3]=fCent;
  valueArray[4]=fCorrMCStatus;

	if(fPlotQA==0)
	{
		fTriggerHist->Fill(valueArray,Weight);
	}
}
///
/// Fill histograms with cluster and track information
///
//________________________________________________________________________
//void AliAnalysisTaskGammaHadron::FillGhHistograms(Int_t identifier,AliTLorentzVector ClusterVec,AliVParticle* TrackVec, Double_t Weight)
void AliAnalysisTaskGammaHadron::FillGhHistograms(Int_t identifier,AliTLorentzVector ClusterVec,AliVParticle* TrackVec, Int_t CorrMCStatus, Double_t Weight)
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::FillGhHistograms()"<<endl;

	Double_t pi = TMath::Pi();

	//..This function fills several histograms under different cut conditions.
	//..it is run within a cluster{ track{}} loop to get all combinations.

	//..A word to the weight - for mixed events it devides by the number of events in the current pool 1/nEvents
	//..                     - for same events you devide by the number of triggers (but not now ->later in the post analysis)
	//..                     - for both you have to take into account the efficiency of your correlated pair
	Double_t G_PT_Value = ClusterVec.Pt();

	Double_t deltaEta   = ClusterVec.Eta()-TrackVec->Eta();
	Double_t deltaPhi   = DeltaPhi(ClusterVec,TrackVec);
	//Double_t ZT_Value   = TMath::Cos(deltaPhi)*TrackVec->P()/ClusterVec.P(); //   TrackVec->Pt()/G_PT_Value;
	Double_t ZT_Value   = TrackVec->Pt()/G_PT_Value; //   TrackVec->Pt()/G_PT_Value;
	//..Careful here: usually this is done for an opening angle (hadron-jet axis) of less than 90¡. Due to
	//..resolution momentum smearing (our guess - check that!) there are particles appearing at angles greater than 90¡
	Double_t XI_Value=-50;
	if(ZT_Value>0)
	{
		XI_Value   = TMath::Log(1.0/ZT_Value);
	}
	Double_t zVertex = fVertex[2];
	//..all from EMCal base class : fEPV0,fEPV0A,fEPV0C
	//fEPV0  = aliEP->GetEventplane("V0" ,InputEvent());
	//fEPV0A = aliEP->GetEventplane("V0A",InputEvent());
	//fEPV0C = aliEP->GetEventplane("V0C",InputEvent());
	Double_t evtPlaneAngle=DeltaPhi(ClusterVec,fQnCorrEventPlaneAngle); //fEPV0);
	Int_t evtPlaneCategory=-1;
	Double_t angleFromAxis;
	//..fold around 0 axis
	angleFromAxis=fabs(evtPlaneAngle);
	//..fold once more arounad pi/2 axis
	if((pi-angleFromAxis)<angleFromAxis)angleFromAxis = pi-angleFromAxis;
	//.. --> now we have only one quadrant left 0 to <pi/2
	if(angleFromAxis>=0 && angleFromAxis<pi/6.)           evtPlaneCategory=0;
	else if (angleFromAxis>=pi/6. && angleFromAxis<pi/3.) evtPlaneCategory=1;
	else if (angleFromAxis>=pi/3. && angleFromAxis<=pi/2.)evtPlaneCategory=2;

  double fCorrMCStatus = CorrMCStatus;
 // if (fIsMC) {
 // }

//	Double_t valueArray[8];
	Double_t valueArray[9];
	valueArray[0]=deltaPhi;
	valueArray[1]=deltaEta;
	valueArray[2]=G_PT_Value;
	valueArray[3]=ZT_Value;
  if (bEnableTrackPtAxis) valueArray[4]=TrackVec->Pt();
	else valueArray[4]=XI_Value;
	valueArray[5]=zVertex;
	valueArray[6]=evtPlaneCategory;
	valueArray[7]=fCent;
  valueArray[8]=fCorrMCStatus;

	if(identifier==0 && fPlotQA==0)fCorrVsManyThings->Fill(valueArray,Weight);

	//..Histograms to test the binning
	//fHistBinCheckEvtPl[identifier] ->Fill(evtPlaneAngle*fRtoD,Weight);
	fHistBinCheckEvtPl[identifier] ->Fill(fQnCorrEventPlaneAngle*fRtoD,Weight);
	fHistBinCheckEvtPl2[identifier]->Fill(angleFromAxis*fRtoD,Weight);
	fHistBinCheckPt[identifier] ->Fill(G_PT_Value,Weight);
	fHistBinCheckZt[identifier] ->Fill(ZT_Value,Weight);
	fHistBinCheckXi[identifier] ->Fill(XI_Value,Weight);
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::FillQAHistograms(Int_t identifier,AliClusterContainer* clusters,AliVCluster* caloCluster,AliVParticle* TrackVec,Double_t weight)
{
	TLorentzVector caloClusterVec;
	clusters->GetMomentum(caloClusterVec,caloCluster);
	AliTLorentzVector aliCaloClusterVec = AliTLorentzVector(caloClusterVec); //..can acess phi from
	Double_t energy= caloCluster->GetUserDefEnergy(fClusEnergyType);

	if(identifier==0 && fPlotQA>0 && fGammaOrPi0==0 && energy >=1)
	{
		//..Get leading cluster
		AliVCluster* leadingClus = GetLeadingCluster("Pt",clusters);

		Double_t etaDistMatched=-1;
		Double_t phiDistMatched=-1;

		//..Distance to SM border
		Int_t etaCellDist;
		Int_t phiCellDist;
		GetDistanceToSMBorder(caloCluster,etaCellDist,phiCellDist);
		if(etaCellDist<0)cout<<"etaCellDist: "<<etaCellDist<<endl;
		if(phiCellDist<0)cout<<"phiCellDist: "<<phiCellDist<<endl;
		Int_t minCellDistance=etaCellDist;
		if(etaCellDist>phiCellDist)minCellDistance=phiCellDist;

		//..ID array -  1 - leading, 2- track matched, 3 - leading & track matched
		Int_t gammaInfo=0;
		if(caloCluster==leadingClus)gammaInfo=1;
		//cout<<"cluster ID: "<<caloCluster->GetID()<<", lead cluster ID: "<<leadingClus->GetID()<<endl;
		if(DetermineMatchedTrack(caloCluster,etaDistMatched,phiDistMatched))gammaInfo=2;
		if(gammaInfo==2 && caloCluster==leadingClus)gammaInfo=3;
		//cout<<"eta distance matched: "<<etaDistMatched<<", phi dist matched: "<<phiDistMatched<<endl;

		//Eg, lambda0,NLM, ncells, distance to bad ,e/p, Mgg
		if(fPlotQA==1)
		{
			Double_t valueArray[6];
			valueArray[0] = energy;
			valueArray[1] = caloCluster->GetM02();
			valueArray[2] = caloCluster->GetDistanceToBadChannel()-1; //..shift to -1 since it starts at 1 and not at 0
			valueArray[3] = minCellDistance;                          //..closest distance to SM border
			valueArray[4] = caloClusterVec.Eta();
			valueArray[5] = aliCaloClusterVec.Phi_0_2pi()*fRtoD;
			fClusterProp->Fill(valueArray); //..all clusters - no cuts
		}

		if(fPlotQA==2)
		{
			Int_t Ntrks = caloCluster->GetNTracksMatched();
			//Double_t etadiff=100;
			//Double_t phidiff=100;
			Double_t etadiffTrack;
			Double_t phidiffTrack;
			//Double_t mom=0;
			Double_t momTrack;
			Double_t valueArray[6];

			//..loop over matched tracks
			for (Int_t i = 0; i < Ntrks; ++i)
			{
				AliVTrack* track = static_cast<AliVTrack*>(caloCluster->GetTrackMatched(i));

				Double_t veta = track->GetTrackEtaOnEMCal();
				Double_t vphi = track->GetTrackPhiOnEMCal();
				momTrack      = track->P();

				Float_t pos[3] = {0};
				caloCluster->GetPosition(pos);
				TVector3 cpos(pos);
				Double_t ceta     = cpos.Eta();
				Double_t cphi     = cpos.Phi();

				etadiffTrack  = veta-ceta;
				phidiffTrack  = TVector2::Phi_mpi_pi(vphi-cphi);
				/* only closest track
			 	if(sqrt(pow(etadiffTrack,2)+pow(phidiffTrack,2))<sqrt(pow(etadiff,2)+pow(phidiff,2)))
				{
					etadiff = etadiffTrack;
					phidiff = phidiffTrack;
					mom     = momTrack;
				}
			    }*/
				//..Fill for every track
				valueArray[0] = energy;
				valueArray[1] = caloCluster->GetM02();
				valueArray[2] = fabs(etadiffTrack);
				valueArray[3] = fabs(phidiffTrack);
				valueArray[4] = energy/momTrack;
				valueArray[5] = fCent;
				fClusterProp->Fill(valueArray); //..all clusters - no cuts
			}
			if(Ntrks==0)
			{
				valueArray[0] = energy;
				valueArray[1] = caloCluster->GetM02();
				valueArray[2] = 1;
				valueArray[3] = 1;
				valueArray[4] = 0;
				valueArray[5] = fCent;
				fClusterProp->Fill(valueArray); //..all clusters - no cuts
			}
		}
	}
	/*do similar test here?*/fHistDEtaDPhiGammaQA[identifier] ->Fill(caloClusterVec.Eta(),aliCaloClusterVec.Phi_0_2pi()*fRtoD,weight);
	if(TrackVec)             fHistDEtaDPhiTrackQA[identifier] ->Fill(TrackVec->Eta(),TrackVec->Phi()*fRtoD,weight);
	fHistClusterTime[identifier]  ->Fill(caloCluster->GetTOF()*1000000000,caloCluster->GetUserDefEnergy(fClusEnergyType),weight);
}
//
// Fills histograms for MC Pi0 information
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::FillMCPi0Hists(Int_t fMultiplicity) {
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::FillMCPi0Hists()"<<endl;
  Int_t nMCPi0s = fMCPi0List.size();
  for (Int_t i = 0; i < nMCPi0s; i++) {
    AliAODMCParticle * fMCPi0_Part = fMCPi0List[i];

    AliTLorentzVector vMCPi0;
    fMCPi0_Part->Momentum(vMCPi0);

    Double_t fLocalPhi = fMCPi0_Part->Phi();
    if (fLocalPhi < 0) fLocalPhi += TMath::TwoPi();

    fEtaPhiMCPion->Fill(fMCPi0_Part->Eta(),fLocalPhi);

    // Reconstructed angle
    Double_t evtPlaneAngle=abs(DeltaPhi(vMCPi0,fQnCorrEventPlaneAngle)); // DeltaPhi is [-pi,pi]
    //if (evtPlaneAngle < 0) evtPlaneAngle += TMath::TwoPi();
    if ((TMath::Pi() - evtPlaneAngle) < evtPlaneAngle) evtPlaneAngle = TMath::Pi() - evtPlaneAngle;
    fPtEPAngleMCPion->Fill(evtPlaneAngle,fMCPi0_Part->Pt());

    // MC Truth angle
    Double_t reactionPlaneAngle=abs(DeltaPhi(vMCPi0,fMCReactionPlaneAngle));
    // Map to [0,-pi/2]
    //if (reactionPlaneAngle < 0) reactionPlaneAngle += TMath::TwoPi();
    if ((TMath::Pi() - reactionPlaneAngle) < reactionPlaneAngle) reactionPlaneAngle = TMath::Pi() - reactionPlaneAngle;
    fPtRPAngleMCPion->Fill(reactionPlaneAngle,fMCPi0_Part->Pt());

    //fMCPi0s = new TH3F("fHistMCPi0s","fHistMCPi0s;p_{T}^{#pi^{0}};",20,fMCPi0Bins,7,fMCPi0EtaBins,kNEMCalMultBins,fMixBEMCalMult->GetXbins()->GetArray());
    fHistMCPi0_PtEtaMult->Fill(fMCPi0_Part->Pt(),fMCPi0_Part->Eta(),fMultiplicity);
  }
}
//
// Accept cluster for analysis. More cuts besides in ApplyClusterCuts and ApplyKinematicCuts
// Apply a hadronic corrections (subtracing MIP energy) if enabled via SetHadronicCorrection
//
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::AccClusterForAna(AliClusterContainer* clusters, AliVCluster* caloCluster)
{
  int iClusterID = caloCluster->GetID();
  // This assumes each cluster has a unique ID
  // Check if already checked:
  if(fClusterAcceptanceStatus[iClusterID] != 0) {
    if (fClusterAcceptanceStatus[iClusterID] == -1) return 0;
    else return 1;
  }

	TLorentzVector caloClusterVec;
	clusters->GetMomentum(caloClusterVec,caloCluster);
	//!!!! eventually transform to AliTLorentzvector

  fClusterAcceptanceStatus[iClusterID] = -1; // This assumes the cluster will be rejected. changed at end

	//..Accepts clusters if certain conditions are fulfilled
	Bool_t Accepted=1; //..By default accepted

	//!!double check these cuts carefully with the experts!!

	//-----------------------------
	//..Check if the cluster energy is above set threshold
  //..This is check is done first for efficiency
  if (caloCluster->GetUserDefEnergy(fClusEnergyType) < fClEnergyMin) {
    return 0;
  }

	//-----------------------------
	//..at least 2 cells in cluster
	if(caloCluster->GetNCells()<2)
	{
		//..Reject the cluster as a good candidate for your analysis
		return 0;
	}
	//-----------------------------
	//..number of local maxima should be 1 or 0 (for cluster splitting this has to be changed)
	if(caloCluster->GetNExMax()>fMaxNLM)
	{
		//..Reject the cluster as a good candidate for your analysis
		return 0;
	}
	//-----------------------------
	//..cut on the cluster shape
	if(fClShapeMin>0 && fClShapeMax>0
	   && (caloCluster->GetM02()<fClShapeMin || caloCluster->GetM02()>fClShapeMax))
	{
		//..Reject the cluster as a good candidate for your analysis
		return 0;
	}
	//-----------------------------
	//..remove clusters with a matched track
	Double_t etaDiff = 0;
	Double_t phiDiff = 0;
	if(fRmvMTrack==1 && DetermineMatchedTrack(caloCluster,etaDiff,phiDiff,0) > 0)
	{
		return 0;
	}
  if (fHadCorr > 0) {
    Double_t fClusterEnergyBeforeCorrection = caloCluster->GetNonLinCorrEnergy();
    Int_t nMatchedTracks = DetermineMatchedTrack(caloCluster,etaDiff,phiDiff,1);    // 1 applies had corr
    fClusterEnergyMatchedTracks->Fill(fClusterEnergyBeforeCorrection,nMatchedTracks);
  }
  //-----------------------------
  //if(fHadCorr>0 && DetermineMatchedTrack(caloCluster,etaDiff,phiDiff))
  if(fHadCorr>0) // Had Corr may lower cluster E below threshold
  {
    //..Check if the cluster energy is now below the threshold
    if (caloCluster->GetUserDefEnergy(fClusEnergyType) < fClEnergyMin) {
      return 0;
    }
  }
	//-----------------------------
	//..Fiducial volume cut. If it is located neither in EMCal nor in DCal reject
	//..Additionally, do we want all EMCal Clusters, ECal only, or DCal only?
	//.. If fSubdetector = 0, accept from ECal and DCal.
	//.. If fSubdetector = 1, accept from only ECal.
	//.. If fSubdetector = 2, accept from only DCal.
	switch (fSubDetector) {
		case 2:
			if(!fFiducialCuts->IsInFiducialCut(caloClusterVec.Eta(),caloClusterVec.Phi(),AliFiducialCut::kDCAL)) return 0;
			break;
		case 1:
			if(!fFiducialCuts->IsInFiducialCut(caloClusterVec.Eta(),caloClusterVec.Phi(),AliFiducialCut::kEMCAL)) return 0;
			break;
		default:
		case 0:
			if(!fFiducialCuts->IsInFiducialCut(caloClusterVec.Eta(),caloClusterVec.Phi(),AliFiducialCut::kEMCAL) &&
			!fFiducialCuts->IsInFiducialCut(caloClusterVec.Eta(),caloClusterVec.Phi(),AliFiducialCut::kDCAL)) return 0;
	}
  //-----------------------------
  //..Fiducial volume cut II. Cuts on the distance to the EMCal border last+first row and last+first collumn
  //fFiducialCellCut->SetNumberOfCellsFromEMCALBorder(1); //ELI this could be momentum dependent and also different for merged clusters!!!
  if(!fFiducialCellCut->CheckCellFiducialRegion(fGeom,caloCluster,fCaloCells))
  {
   return 0;
  }
  //-----------------------------
  //..Exotic cells

  double ClusterPhi2Pi = caloClusterVec.Phi();
  if (ClusterPhi2Pi < 0) ClusterPhi2Pi += 2*TMath::Pi();

  fClusterAcceptanceStatus[iClusterID] = 1; // This marks the cluster as accepted
  fAccClusEtaPhi->Fill(caloClusterVec.Eta(),ClusterPhi2Pi);
  fAccClusEtaPhiZvtx->Fill(caloClusterVec.Eta(),ClusterPhi2Pi,fVertex[2]);
  return Accepted;
}
//
// Accept cluster pair cut for Pi0 analysis. More cuts besides in ApplyClusterCuts and ApplyKinematicCuts
//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::AccClusPairForAna(AliVCluster* cluster1, AliVCluster * cluster2, TLorentzVector vecPi0)
{
	//..Accepts clusters if certain conditions are fulfilled
	Bool_t Accepted=1; //..By default accepted

	Double_t fE1 = cluster1->GetUserDefEnergy(fClusEnergyType);
	Double_t fE2 = cluster2->GetUserDefEnergy(fClusEnergyType);

  // This check should be redundant now ...
  if (fE1<fClEnergyMin || fE2<fClEnergyMin) return 0;  // Check cluster energy min.

	Double_t fAsym = (fE1+fE2 > 0.000001) ? TMath::Abs(fE2-fE1)/(fE1+fE2) : 0; //Don't divide by zero

  if (fAsym > fPi0AsymCut) return 0; // Apply asymmetry cut

  // For reference
  //fPi0MassFitPars  = {d,e,m1,m2,0}
  //fPi0SigmaFitPars = {d,e,m1,m2,0}

  Double_t Pi0Pt    = vecPi0.Pt();
  Double_t Pi0Mass  = 0;
  Double_t Pi0Sigma = 0;
  Double_t SBsplit  = 0;
  // Defining the total sideband region
  Double_t SBLowerRange = 0.16; // To be overidden by pi0mass + 5 sigma
  Double_t SBUpperRange = 0.45;

  if (Pi0Pt > fMaxPi0Pt) return 0;

  if (fUseParamMassSigma)
  {
    // Estimate Mass Peak
    if (Pi0Pt < fPi0MassFitPars[0])
    {
      Pi0Mass = fPi0MassFitPars[2]*Pi0Pt + fPi0MassFitPars[1] - fPi0MassFitPars[2]*fPi0MassFitPars[0];
    } else {
      Pi0Mass = fPi0MassFitPars[3]*Pi0Pt + fPi0MassFitPars[1] - fPi0MassFitPars[3]*fPi0MassFitPars[0];
    }
    // Estimate Mass Sigma
    if (Pi0Pt < fPi0SigmaFitPars[0])
    {
      Pi0Sigma = fPi0SigmaFitPars[2]*Pi0Pt + fPi0SigmaFitPars[1] - fPi0SigmaFitPars[2]*fPi0SigmaFitPars[0];
    } else {
      Pi0Sigma = fPi0SigmaFitPars[3]*Pi0Pt + fPi0SigmaFitPars[1] - fPi0SigmaFitPars[3]*fPi0SigmaFitPars[0];
    }
  }
  else
  { // Using fixed mass windows
    // Finding pT Bin:
    Int_t ptBin = 0; // Default is to use lowest bin

    for (Int_t k = 0; k < kNoGammaBins-1; k++) {
      if ((Pi0Pt >= fArray_G_Bins[k]) && (Pi0Pt < fArray_G_Bins[k+1]))
      {
        ptBin = k;
        break;
      }
    }

    Double_t centBinArray[kNMainCentBins+1]  = {0.0,10.0,30.0,50.0,90.0};
    Int_t iCentBin = 0;

    for (Int_t k = 0; k < kNMainCentBins; k++) {
      if ((fCent >= centBinArray[k]) && (fCent < centBinArray[k+1]))
      {
        iCentBin = k;
        break;
      }
    }
    Pi0Mass  = fPi0MassFixed[iCentBin][ptBin];
    Pi0Sigma = fPi0SigmaFixed[iCentBin][ptBin];
    SBLowerRange = Pi0Mass + Pi0Sigma*5; //..Lower range  = Pi0Peak + 5 sigma
    SBsplit  = (SBUpperRange-SBLowerRange)*0.5;
  }
  //..if you select the pi0 peak region
  if(fGammaOrPi0==1 && TMath::Abs(vecPi0.M() - Pi0Mass) > fPi0NSigma * Pi0Sigma)
  {
    fMassPionRej->Fill(vecPi0.M());
    fMassPtPionRej->Fill(vecPi0.M(),vecPi0.Pt());
    fMassPtCentPionRej->Fill(vecPi0.M(),vecPi0.Pt(),fCent);
    return 0;
  }
	// Sidebands
	if(fGammaOrPi0 >= 2) {  // >= used currently for backwards compatibility.
    Double_t fSBMin = 0;
    Double_t fSBMax = 0;
    switch (fSidebandChoice) {
      case 0: // Largest sideband
      default:
        fSBMin = SBLowerRange;
        fSBMax = SBUpperRange;
        break;
      case 1: // Sideband 1 (1/2 SB 1)
        fSBMin = SBLowerRange;
        fSBMax = SBUpperRange - SBsplit;
        break;
      case 2: // Sideband 2 (1/2 SB 2)
        fSBMin = SBLowerRange + SBsplit;
        fSBMax = SBUpperRange;
        break;
      case 3: // Sideband 3 (1/4 SB 1)
        fSBMin = SBLowerRange;
        fSBMax = (SBUpperRange - 1.5*SBsplit);
        break;
      case 4: // Sideband 4 (1/4 SB 2)
        fSBMin = SBLowerRange + 0.5*SBsplit;
        fSBMax = (SBUpperRange - SBsplit);
        break;
      case 5: // Sideband 5 (1/4 SB 3)
        fSBMin = SBLowerRange + SBsplit;
        fSBMax = (SBUpperRange - 0.5*SBsplit);
        break;
      case 6: // Sideband 6 (1/4 SB 4)
        fSBMin = SBLowerRange + 1.5*SBsplit;
        fSBMax = SBUpperRange;
        break;
      case 7: // Sideband 7 (1/3 SB 1)
        fSBMin = SBLowerRange;
        fSBMax = SBUpperRange - (4./3.)*SBsplit;
        break;
      case 8: // Sideband 8 (1/3 SB 2)
        fSBMin = SBLowerRange + (2./3.)*SBsplit;
        fSBMax = SBUpperRange - (2./3.)*SBsplit;
        break;
      case 9: // Sideband 9 (1/3 SB 3)
        fSBMin = SBLowerRange + (4./3.)*SBsplit;
        fSBMax = SBUpperRange;
        break;
    }
    if (vecPi0.M() < fSBMin || vecPi0.M() >= fSBMax) {
      fMassPionRej->Fill(vecPi0.M());
      fMassPtPionRej->Fill(vecPi0.M(),vecPi0.Pt());
      fMassPtCentPionRej->Fill(vecPi0.M(),vecPi0.Pt(),fCent);
      return 0;
    }
  }
  return Accepted;
}

Bool_t AliAnalysisTaskGammaHadron::QuickCheckAccClus(TLorentzVector ClusterVec) {
  if (!fAccClusEtaPhi) return true; // default to allowing all clusters if histogram missing

// FIXME
  return true;
  double fClusPhi = ClusterVec.Phi();
  if (fClusPhi < 0) fClusPhi += 2 * TMath::Pi();
  // Todo: add option for fAccClusEtaPhiZvtx
  float fCellContent = fAccClusEtaPhi->GetBinContent(ClusterVec.Eta(),fClusPhi);
  if (fCellContent > 0) AliInfo(Form("DEBUG: Quick accept found nonzero bin content at (%f,%f)= %f\n",ClusterVec.Eta(),fClusPhi,fCellContent));
  AliInfo(Form("DEBUG: Quick accept found bin content at (%f,%f)= %f\n",ClusterVec.Eta(),fClusPhi,fCellContent));
  return fCellContent > 0;
}

//________________________________________________________________________
Double_t AliAnalysisTaskGammaHadron::DeltaPhi(AliTLorentzVector ClusterVec,AliVParticle* TrackVec)
{
	Double_t Phi_g = ClusterVec.Phi_0_2pi();
	Double_t Phi_h = TrackVec->Phi();

	Double_t dPhi = -999;
	Double_t pi = TMath::Pi();

	dPhi = Phi_g-Phi_h;
	//--cut the away side peak on the left of the NS peak
	//--and insert it to the very right: \-^-/  ---> -^-/\.
	//--to create a correlation histogram that starts at -pi/2 and ends at 3/2pi
	if (dPhi <= -pi/2)    dPhi += 2*pi;
	if (dPhi > 3.0*pi/2.0)dPhi -= 2*pi;

	return dPhi;
}
//________________________________________________________________________
Double_t AliAnalysisTaskGammaHadron::DeltaPhi(AliTLorentzVector ClusterVec,Double_t phi_EVP)
{
	Double_t phi_g = ClusterVec.Phi_0_2pi();

	Double_t dPhi = -999;
	Double_t pi = TMath::Pi();

	//..delta phi between gamma and evt plane
	//..phi_g is [80,328] and phi_EVP [-90,90]
	dPhi = phi_g-phi_EVP;
	//..to create an angle that starts at 0 and ends at 2pi
	if (dPhi >= 2*pi)dPhi -= 2*pi;
	if (dPhi < 0)    dPhi += 2*pi;
	//..to create an angle that starts at -pi and ends at pi
	if (dPhi >= pi)dPhi -= 2*pi;

	return dPhi;
}
//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::FindMCPartForClus(AliVCluster * caloCluster) {
	if (!fMCParticles) {
		return -1;
	}
	TLorentzVector caloClusterVec;
	caloCluster->GetMomentum(caloClusterVec,0);

	Double_t fClus_E   = caloCluster->GetUserDefEnergy(fClusEnergyType);
	Double_t fClus_Eta = caloClusterVec.Eta();

	Int_t iLabel = caloCluster->GetLabel();
	AliAODMCParticle * mcPart = 0;
	if (iLabel > -1) {
		mcPart = fMCParticles->GetMCParticle(iLabel);
		Double_t fDeltaPhi = DeltaPhi(caloClusterVec,mcPart->Phi());
		Double_t fDeltaEta = fClus_Eta - mcPart->Eta();
		Double_t fDeltaE   = fClus_E  -  mcPart->E();

		fHistClusMCDE->Fill(fClus_E,fDeltaE);
		fHistClusMCDPhiDEta->Fill(fDeltaPhi,fDeltaEta);
	}

	return caloCluster->GetLabel();
}
// Return the id of the root particle for the given (by index) MC Particle
// If the particle is its own root, return the particle's ID
//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::FindMCRootPart(Int_t iMCIndex, Int_t * iMCTreeHeight) {
	if (iMCIndex == -1) return -1;
	if (!fMCParticles) {
		return -1;
	}
	AliAODMCParticle * mcPart = fMCParticles->GetMCParticle(iMCIndex);
	if (!mcPart) {return -2;}

	*iMCTreeHeight = 1;

	AliAODMCParticle * pMotherParticle = 0;

	Int_t iCurrentAnswer = iMCIndex;
	Int_t iMother = mcPart->GetMother();

	for (Int_t i = 0; i < 22; i++) {
		if (iMother == -1) return iCurrentAnswer;
		iCurrentAnswer = iMother;
		(*iMCTreeHeight)++;
		pMotherParticle = fMCParticles->GetMCParticle(iMother);
		if (!pMotherParticle) return iCurrentAnswer;
		iMother = pMotherParticle->GetMother();
	}

	return iCurrentAnswer;
}

// Return the id of the Lowest Common Ancestor (LCA) of the two given root particles.
// Return -1 if they do not have an LCA.
//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::FindMCLowComAnc(Int_t iMCIndex1, Int_t iMCIndex2) {
	if (iMCIndex1 == -1 || iMCIndex2 == -1) return -1;
	if (!fMCParticles) {
		return -1;
	}
	AliAODMCParticle * mcPart1 = (AliAODMCParticle *) fMCParticles->GetMCParticle(iMCIndex1);
	AliAODMCParticle * mcPart2 = (AliAODMCParticle *) fMCParticles->GetMCParticle(iMCIndex2);
	if (!mcPart1 || !mcPart2) return -1;

	// Lists of ancestors;
	std::vector<Int_t> lAncPart_1 = {iMCIndex1};
	std::vector<Int_t> lAncPart_2 = {iMCIndex2};

	// Avoid infinite loops if MC particles aren't set up right
	Int_t iCutOff = 50;
	Int_t iMinHeight = 1;

	// Building Ancestor Lists
	AliAODMCParticle * fTempMC1 = mcPart1;
	AliAODMCParticle * fTempMC2 = mcPart2;

	for (Int_t i = 0; i < iCutOff; i++) {
		Int_t iMother1 = -1;
		Int_t iMother2 = -1;
		if (fTempMC1) iMother1 = fTempMC1->GetMother();
		if (fTempMC2) iMother2 = fTempMC2->GetMother();

		if (iMother1 == -1) {
			fTempMC1 = 0;
    } else {
			lAncPart_1.push_back(iMother1);
			fTempMC1 = fMCParticles->GetMCParticle(iMother1);
		}
		if (iMother2 == -1) {
			fTempMC2 = 0;
		} else {
			fTempMC2 = fMCParticles->GetMCParticle(iMother2);
			lAncPart_2.push_back(iMother2);
		}
		if ((iMother1 == -1) && (iMother2 == -1)) break;
	}

	Int_t iSize1 = lAncPart_1.size();
	Int_t iSize2 = lAncPart_2.size();

	Int_t iLastCommonAncestor = -1;

	iMinHeight = std::min(lAncPart_1.size(),lAncPart_2.size());
	for (Int_t i = 0; i < iMinHeight; i++) {
		if (lAncPart_1[iSize1-i] == lAncPart_2[iSize2-i]) {
			iLastCommonAncestor = lAncPart_1[iSize1-i];
		}
	}

	return iLastCommonAncestor;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::DetermineGAPatchCand(AliTLorentzVector CaloClusterVec, AliTLorentzVector CaloClusterVec2) {

	Double_t fE1 = CaloClusterVec.E();
	Double_t fE2 = CaloClusterVec2.E();

	// From OCDB, could get these dynamically
	if (fE1 + fE2 < 10) return 0;// 10 is the GA thr in 15o
	if ((fE1 >= 10) || (fE2 >= 10)) return 1; // one cluster is a patch cand on its own

	Double_t fDeltaPhi = DeltaPhi(CaloClusterVec2,CaloClusterVec.Phi());
	Double_t fDeltaEta = CaloClusterVec2.Eta() - CaloClusterVec.Eta();

	if (abs(fDeltaPhi) > 4*0.0143) return 0;
	if (abs(fDeltaEta) > 4*0.0143) return 0;

	return 1;
}
//
//  Determine if this cluster has a track matched to it
//  If bApplyHadCorr == 1, then hadronic correction is applied during this step, and the
//  return value will be the number of matched tracks
//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::DetermineMatchedTrack(AliVCluster* caloCluster,Double_t &etadiff,Double_t & phidiff, Bool_t bApplyHadCorr)
{
	Int_t foundTrackMatched=0;
	Int_t Ntrks = caloCluster->GetNTracksMatched();
	if(Ntrks==0) return foundTrackMatched; //..if no matched track removal is wanted set it to 0.

	//..loop over matched tracks
	for (Int_t i = 0; i < Ntrks; ++i)
	{
		AliVTrack* track = static_cast<AliVTrack*>(caloCluster->GetTrackMatched(i));

		Double_t veta = track->GetTrackEtaOnEMCal();
		Double_t vphi = track->GetTrackPhiOnEMCal();

		Float_t pos[3] = {0};
		caloCluster->GetPosition(pos);
		TVector3 cpos(pos);
		Double_t ceta     = cpos.Eta();
		Double_t cphi     = cpos.Phi();
		Double_t fEtaDiff=veta-ceta;
		Double_t fPhiDiff=TVector2::Phi_mpi_pi(vphi-cphi);

		//?  // check if track also points to cluster
		//?   Int_t cid = track->GetEMCALcluster();
		//?   if (fDoTrackClus && (cid != icluster)) return energyclus;

		//..check if the track was matched within a stricter criteria
		Double_t etaCut = fTrackMatchEta;
		Double_t phiCut = fTrackMatchPhi;
		Double_t trackPt = track->Pt();

		fMatchDeltaPhiTrackPt->Fill(trackPt,fPhiDiff);
		fMatchDeltaEtaTrackPt->Fill(trackPt,fEtaDiff);
		//
		// For input -1 or 0, use the parametrized track matching cuts given in https://alice-notes.web.cern.ch/node/813
    // FIXME make these (changeable) parameters
		if (etaCut <= 0) etaCut = 0.010 + TMath::Power((trackPt + 4.07), -2.5);
		if (phiCut <= 0) phiCut = 0.015 + TMath::Power((trackPt + 3.65), -2.);

		if(TMath::Abs(fEtaDiff)<etaCut) fMatchCondDeltaPhiTrackPt->Fill(trackPt,fPhiDiff);
		if(TMath::Abs(fPhiDiff)<phiCut) fMatchCondDeltaEtaTrackPt->Fill(trackPt,fEtaDiff);

		if(TMath::Abs(fEtaDiff)<etaCut && TMath::Abs(fPhiDiff)<phiCut)
		{
			if (!fIsMC || (fMinMCLabel <= 0 || TMath::Abs(track->GetLabel()) > fMinMCLabel)) // label check copied from AliClusterContainer
			{
        Double_t fP = track->P();
        Double_t fE = caloCluster->GetUserDefEnergy(AliVCluster::kNonLinCorr); // Use the non linearity corrected energy here
        if (fP > 1e-6) fHistEOverPvE->Fill(fE,fE/fP);
        if (fE > 1e-6) fHistPOverEvE->Fill(fE,fP/fE);
//					if (track->P() > 1e-6) fHistEOverPvE->Fill(caloCluster->GetNonLinCorrEnergy()/track->P(),caloCluster->GetNonLinCorrEnergy()); 
//					if (caloCluster->GetNonLinCorrEnergy() > 1e-6) fHistPOverEvE->Fill(track->P()/caloCluster->GetNonLinCorrEnergy(),caloCluster->GetNonLinCorrEnergy()); 					
//					if ((fP > 1e-6) && (fE/fP > fTrackMatchEOverPLow) && (fE/fP < fTrackMatchEOverPHigh)) {
        if (fP > 1e-6 && fTrackMatchEOverPHigh > 0.) {
//          if ((fE/fP < fTrackMatchEOverPLow) || (fE/fP > fTrackMatchEOverPHigh && fTrackMatchEOverPHigh > 0.)) {
          if ((fE/fP < fTrackMatchEOverPLow) || (fE/fP > fTrackMatchEOverPHigh)) {
            continue; // Don't accept match
          }
        }
        foundTrackMatched++;
				etadiff = fEtaDiff;
				phidiff = fPhiDiff;
        if (!bApplyHadCorr) break; // If we are not applying the hadronic correction, break at first match
        caloCluster->SetHadCorrEnergy(fE - fHadCorrConstant);
        // may get negative energy clusters
			}
		}
	}
	return foundTrackMatched;
}
//
// Calculate the distance to a supermodule border
// calculated is the smaller distance in either direction of the max. energy cell of the cluster
//
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::GetDistanceToSMBorder(AliVCluster* caloCluster,Int_t &etaCellDist,Int_t &phiCellDist)
{
	//..celldist 0 means it is at the border of the EMCal
	//..celldist 1 means there is one row/collum between the cell and the border
	etaCellDist=-1;
	phiCellDist=-1;
	//..partially copied from RecoUtils
	Int_t absIdMax  = -1, iSM =-1, ieta = -1, iphi = -1;
	Bool_t shared = kFALSE;
	fFiducialCellCut->GetMaxEnergyCell(fGeom, fCaloCells, caloCluster, absIdMax,  iSM, ieta, iphi, shared);
	if (absIdMax==-1) cout<<"In GetDistanceToSMBorder  why??"<<endl;

	//cout<<"max cell eta collumn: "<<ieta<<", max cell phi row: "<<iphi<<endl;
	//..Define Borders
	Int_t lastRow,firstRow;
	Int_t lastCollumn,firstCollumn;
	//..phi
	firstRow=0;
	lastRow=24;
	if      (fGeom->GetSMType(iSM) == AliEMCALGeometry::kEMCAL_Half) lastRow /= 2;
	else if (fGeom->GetSMType(iSM) == AliEMCALGeometry::kEMCAL_3rd ) lastRow /= 3;// 1/3 sm case
	else if (fGeom->GetSMType(iSM) == AliEMCALGeometry::kDCAL_Ext  ) lastRow /= 3;// 1/3 sm case
	//..eta
	firstCollumn=0;
	lastCollumn=48;
	if(fGeom->GetSMType(iSM) == AliEMCALGeometry::kDCAL_Standard )  lastCollumn = lastCollumn*2/3;

	lastRow    =lastRow-1;      //..range starts from 0 and only goes up to 23
	lastCollumn=lastCollumn-1;  //..range starts from 0

	//..Calculate smallest distance
	//..phi
	phiCellDist=iphi; //eli in case it starts from 1
    if((iphi-firstRow)>=(lastRow-iphi))phiCellDist=lastRow-iphi;

    //..eta
    if(iSM%2==0)
    {
    	  //..only outer border (not at eta=0) except DCal modules
      etaCellDist=ieta;
      if(lastCollumn<47 && (ieta-firstCollumn)>=(lastCollumn-ieta))etaCellDist=lastCollumn-ieta;
    }
    else
    {
    	  //..only outer border (not at eta=0) except DCal modules
    	  etaCellDist=lastCollumn-ieta;
      if(lastCollumn<47 && (ieta-firstCollumn)<=(lastCollumn-ieta))etaCellDist=ieta-firstCollumn;
    }
    //cout<<"supermodule: "<<iSM<<", last collumn: "<<lastCollumn<<endl;
    //cout<<"eta distance: "<<etaCellDist<<" phi distance: "<<phiCellDist<<endl;
}
//
// Get leading cluster of clusters that pass
// the standard cuts for the analysis
//________________________________________________________________________
AliVCluster* AliAnalysisTaskGammaHadron::GetLeadingCluster(const char* opt,AliClusterContainer* clusters)
{
	//..2 options are leading by pT and leading by E
	TString option(opt);

	AliVCluster *clusterMaxE = 0;
	AliVCluster *clusterMaxPt = 0;
	AliVCluster *cluster = 0;
	Double_t et=0;
	Double_t etmax=0;

	Bool_t first=1;
	for(Int_t NoCluster1 = 0; NoCluster1 < clusters->GetNClusters(); NoCluster1++ )
	{
		cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
		if(!cluster)continue;                            //check if the cluster is a good cluster
		if(!AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster

		if(first==1)
		{
			clusterMaxE=cluster;
			clusterMaxPt=cluster;
			first=0;
		}

		if (cluster->E() > clusterMaxE->E()) clusterMaxE = cluster;

		TLorentzVector mom;
		clusters->GetMomentum(mom,cluster);
		et = mom.Et();
		if (et > etmax)
		{
			clusterMaxPt = cluster;
			etmax = et;
		}
	}
	//if (option.Contains("Pt"))    cout<<"selected pt option"<<endl;
	//else if (option.Contains("E"))cout<<"selected E option"<<endl;
	//else  cout<<"no option selected"<<endl;

	if (option.Contains("Pt"))    return clusterMaxPt;
	else if (option.Contains("E"))return clusterMaxE;
	else return 0;
}
//________________________________________________________________________
Double_t AliAnalysisTaskGammaHadron::GetTrackEff(Double_t pT, Double_t eta)
{
	Double_t DetectionEff=1;
	if(fCorrectEff==0)return DetectionEff;

	//..Check which centrality
	Int_t centBin;
	if(fCent<=0.1)centBin=0;
	if(fCent>0.1 && fCent<=0.3)centBin=1;
	if(fCent>0.3 && fCent<=0.5)centBin=2;
	if(fCent>0.5 && fCent<=0.9)centBin=3;

	if(pT<=3.5)
	{
		DetectionEff=funcpT_low[centBin]->Eval(pT);
	}
	else
	{
		DetectionEff=funcpT_high[centBin]->Eval(pT);
	}
	//..eta part
	if(eta<=-0.04)
	{
		DetectionEff*=funcpEta_left[centBin]->Eval(eta,0,0)/fscaleEta[centBin];
	}
	else
	{
		DetectionEff*=funcpEta_right[centBin]->Eval(eta,0,0)/fscaleEta[centBin];
	}

	return DetectionEff;
}
//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::CalculateEventHash() {

  //UInt_t fTimeStamp = InputEvent()->GetTimeStamp();
  Int_t nTracks = InputEvent()->GetNumberOfTracks();
  Int_t nClusters = InputEvent()->GetNumberOfCaloClusters();
  //Double_t fT0TOF = (Double_t) InputEvent()->GetT0TOF()[0];

  UInt_t iOrbitNumber = InputEvent()->GetOrbitNumber();
  //Int_t iEventNumberInFile = InputEvent()->GetEventNumberInFile();

  if (!fIsMC) { // Data
    return (iOrbitNumber) % 2;
  } else { // MC
    // orbit number and bunch cross not defined in MC,
    return (nTracks-nClusters) % 2;
  }
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::LoadQnCorrectedEventPlane() {
  //..This function is called at the beginning of FillHistograms
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::LoadQnCorrectedEventPlane()"<<endl;

  Double_t fZVertex = fVertex[2];

  fQnCorrEventPlaneAngle = fEPV0; // Default to uncorrected V0M if QnVectors are not set up correctly

  if (fIsMC && fUseMCReactionPlane && fMCHeader) {
    fQnCorrEventPlaneAngle = fMCHeader->GetReactionPlaneAngle();
    return;
  }

  if (fFlowQnVectorMgr == 0) return;

  // Want to set fQnCorrEventPlaneAngle

  Int_t iHarmonic = 2;

  const AliQnCorrectionsQnVector * fV0MQnVector;
  const AliQnCorrectionsQnVector * fTPCAQnVector;
  const AliQnCorrectionsQnVector * fTPCCQnVector;
  Double_t fV0MQnEP = 0.0;
  Double_t fTPCAQnEP = 0.0;
  Double_t fTPCCQnEP = 0.0;

  fV0MQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROQoverM");
  fTPCAQnVector = fFlowQnVectorMgr->GetDetectorQnVector("TPCPosEtaQoverM");
  fTPCCQnVector = fFlowQnVectorMgr->GetDetectorQnVector("TPCNegEtaQoverM");

  if (fV0MQnVector != NULL) fV0MQnEP = fV0MQnVector->EventPlane(iHarmonic); else return;
  if (fTPCAQnVector != NULL) fTPCAQnEP = fTPCAQnVector->EventPlane(iHarmonic); else return;
  if (fTPCCQnVector != NULL) fTPCCQnEP = fTPCCQnVector->EventPlane(iHarmonic); else return;

  // Second Order Event Plane

  Double_t fDPsi1 = fV0MQnEP - fTPCAQnEP;
  Double_t fDPsi2 = fV0MQnEP - fTPCCQnEP;
  Double_t fDPsi3 = fTPCAQnEP - fTPCCQnEP; // A-side is probably eta > 0

  fEPAngleV0M->Fill(fV0MQnEP);
  fEPAngleTPCA->Fill(fTPCAQnEP);
  fEPAngleTPCC->Fill(fTPCCQnEP);

  for (Int_t iOrder = 0; iOrder < kNumEPROrders; iOrder++) {
    fEPR_CosD1[iOrder]->Fill(fZVertex,fCent,TMath::Cos((iOrder+1)*fDPsi1));
    fEPR_CosD2[iOrder]->Fill(fZVertex,fCent,TMath::Cos((iOrder+1)*fDPsi2));
    fEPR_CosD3[iOrder]->Fill(fZVertex,fCent,TMath::Cos((iOrder+1)*fDPsi3));
  }

  fQnCorrEventPlaneAngle = fV0MQnEP; // We use V0 Combination

  // Third Order Event Plane
  iHarmonic = 3;
  fV0MQnEP = fV0MQnVector->EventPlane(iHarmonic);
  fTPCAQnEP = fTPCAQnVector->EventPlane(iHarmonic);
  fTPCCQnEP = fTPCCQnVector->EventPlane(iHarmonic);

  fEP3AngleV0M->Fill(fV0MQnEP);
  fEP3AngleTPCA->Fill(fTPCAQnEP);
  fEP3AngleTPCC->Fill(fTPCCQnEP);

  fDPsi1 = fV0MQnEP - fTPCAQnEP;
  fDPsi2 = fV0MQnEP - fTPCCQnEP;
  fDPsi3 = fTPCAQnEP - fTPCCQnEP;

  for (Int_t iOrder = 0; iOrder < kNumEPROrders; iOrder++) {
    fEP3R_CosD1[iOrder]->Fill(fZVertex,fCent,TMath::Cos((iOrder+1)*fDPsi1));
    fEP3R_CosD2[iOrder]->Fill(fZVertex,fCent,TMath::Cos((iOrder+1)*fDPsi2));
    fEP3R_CosD3[iOrder]->Fill(fZVertex,fCent,TMath::Cos((iOrder+1)*fDPsi3));
  }

  fQnCorrEventPlane3Angle = fV0MQnEP; // We use V0 Combination

  // Fourth Order Event Plane
  iHarmonic = 4;
  fV0MQnEP = fV0MQnVector->EventPlane(iHarmonic);
  fTPCAQnEP = fTPCAQnVector->EventPlane(iHarmonic);
  fTPCCQnEP = fTPCCQnVector->EventPlane(iHarmonic);

  fEP4AngleV0M->Fill(fV0MQnEP);
  fEP4AngleTPCA->Fill(fTPCAQnEP);
  fEP4AngleTPCC->Fill(fTPCCQnEP);

  fDPsi1 = fV0MQnEP - fTPCAQnEP;
  fDPsi2 = fV0MQnEP - fTPCCQnEP;
  fDPsi3 = fTPCAQnEP - fTPCCQnEP;

  for (Int_t iOrder = 0; iOrder < kNumEPROrders; iOrder++) {
    fEP4R_CosD1[iOrder]->Fill(fZVertex,fCent,TMath::Cos((iOrder+1)*fDPsi1));
    fEP4R_CosD2[iOrder]->Fill(fZVertex,fCent,TMath::Cos((iOrder+1)*fDPsi2));
    fEP4R_CosD3[iOrder]->Fill(fZVertex,fCent,TMath::Cos((iOrder+1)*fDPsi3));
  }

  fQnCorrEventPlane4Angle = fV0MQnEP; // We use V0 Combination


}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
AliAnalysisTaskGammaHadron* AliAnalysisTaskGammaHadron::AddTaskGammaHadron(
		Int_t       InputGammaOrPi0,//..gamma analysis=0, pi0 analyis=1, pi0 SB1=2, pi0 SB2=3,
		Int_t       InputSeMe,      //..same event=0 mixed event =1 mixed trigger = 2
		Bool_t      InputMCorData,  // 0->Data, 1->MC 
		UInt_t      evtTriggerType, //..use this type of events to combine gammas(trigger) with hadrons
		UInt_t      evtMixingType,  //..use only this type of events to fill your mixed event pool with tracks
		Bool_t      isRun2,         //..changes some settigs and cuts depending on 2013 or 2015/2016 data
		Double_t    trackptcut,     //..
		Double_t    clusEcut,      //..
		Bool_t      SavePool,       //..saves a mixed event pool to the output event
		const char *trackName,
		const char *clusName,
		const char *taskname,
		const char *suffix)
{
	//cout<<"in AddTaskGammaHadron.C(...)"<<endl;
	// Get the pointer to the existing analysis manager via the static access method.
	//==============================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr)
	{
		::Error("AddTaskGammaHadron", "No analysis manager to connect to.");
		return 0;
	}

	// Check the analysis type using the event handlers connected to the analysis manager.
	//==============================================================================
	AliVEventHandler* handler = mgr->GetInputEventHandler();
	if (!handler)
	{
		::Error("AddTaskGammaHadron", "This task requires an input event handler");
		return 0;
	}
	if (handler->InheritsFrom("AliESDInputHandler"))
	{
		::Error("AddTaskGammaHadron", "We have never taken care if this works for ESDs");
		return 0;
	}

	//..in case of AOD the default names are:
	//if(trackName=="usedefault")trackName = "tracks";
	//if(clusName =="usedefault")clusName  = "caloClusters";
	if(strncmp(trackName,"usedefault",10)==0)trackName = "tracks";
	if(strncmp(clusName,"usedefault",10)==0) clusName  = "caloClusters";

	//-------------------------------------------------------
	// Built the name of the Task together
	//-------------------------------------------------------
	TString GammaPi0Name;
	if(InputGammaOrPi0 == 0)
	{
		GammaPi0Name += "GH";
	}
	else if(InputGammaOrPi0 == 1)
	{
		GammaPi0Name += "Pi0H";
	}
	else if(InputGammaOrPi0 == 2)
	{
		GammaPi0Name += "Pi0H_SB1";
	}
	else if(InputGammaOrPi0 == 3)
	{
		GammaPi0Name += "Pi0H_SB2";
	}
	else if(InputGammaOrPi0 == 4)
	{
		GammaPi0Name += "Pi0H_SB";
	}
	else if(InputGammaOrPi0 == 5)
	{
		GammaPi0Name += "Pi0H_SB3";
	}
	else if(InputGammaOrPi0 == 6)
	{
		GammaPi0Name += "Pi0H_SB4";
	}
	else if(InputGammaOrPi0 == 7)
	{
		GammaPi0Name += "Pi0H_SB5";
	}
	else if(InputGammaOrPi0 == 8)
	{
		GammaPi0Name += "Pi0H_SB6";
	}
	TString SameMixName;
	if(InputSeMe == 0)
	{
		SameMixName += "SE";
	}
	else
	{
		SameMixName += "ME";
	}

	TString combinedName;
	combinedName.Form("%s_%s_%s_%s_%s",taskname,(const char*)GammaPi0Name,(const char*)SameMixName,trackName,clusName);
	if(strncmp(suffix,"",1)!=0)
	//if(suffix!="")
	{
		combinedName += "_";
		combinedName += suffix;
	}
	cout<<"combinedName: "<<combinedName<<endl;
	TString contName(combinedName);
	contName += "_histos";

	//-------------------------------------------------------
	// Init the task and do settings
	//-------------------------------------------------------
	AliAnalysisTaskGammaHadron* AnalysisTask = new AliAnalysisTaskGammaHadron(InputGammaOrPi0,InputSeMe,InputMCorData);

	//..Add the containers and set the names
	AnalysisTask->AddClusterContainer(clusName);
	if (InputMCorData) {
		AnalysisTask->AddMCParticleContainer("mcparticles");
	}
	if (strncmp(trackName,"mcparticles",11)==0)
	{
		AliMCParticleContainer* mcpartCont = AnalysisTask->AddMCParticleContainer(trackName);
		mcpartCont->SelectPhysicalPrimaries(kTRUE);
	}
	else if (strncmp(trackName,"tracks",6))
	{
		AliTrackContainer* trackCont = AnalysisTask->AddTrackContainer(trackName);
		trackCont->SetFilterHybridTracks(kTRUE); //gives me Hyprid tracks
	}
	else  //implemented for testing correction framework
	{
		AliTrackContainer* trackCont = AnalysisTask->AddTrackContainer(trackName);
		trackCont->SetFilterHybridTracks(kTRUE); //gives me Hyprid tracks
	}
	//..check that condition!! maybe for mixed events its different!!!!!!
	if(!AnalysisTask->GetTrackContainer(trackName) || !AnalysisTask->GetClusterContainer(clusName))
	{
		cout<<"Task can not run like this!"<<endl;
		return 0;
	}

	//-------------------------------------------------------
	// Add some selection criteria
	//-------------------------------------------------------
	//..set the beamtype and the run2 flag
	Double_t    trackEta   = 0.9;    //..+- eta range for track acceptance
	Double_t    clusterEta = 0.7;    //..+- eta range for cluster acceptance

	AnalysisTask->SetOffTrigger(evtTriggerType|evtMixingType); //..select only evets of type evtTriggerType and evtMixingType
	AnalysisTask->SetNeedEmcalGeom(kTRUE);
	//..for Run1 pPb
	if(isRun2==0)
	{
		AnalysisTask->SetUseManualEvtCuts(kTRUE);
		AnalysisTask->SetUseAliAnaUtils(kTRUE);
		AnalysisTask->SetVzRange(-10,10);
		AnalysisTask->SetCentRange(0.0,100.0);
		// AnalysisTask->SetCentralityEstimator("ZNA");
	}
	//..new task for run2
	if(isRun2==1)
	{
		AnalysisTask->SetNCentBins(5);
		AnalysisTask->SetUseNewCentralityEstimation(kTRUE);
	}

	if(AnalysisTask->GetTrackContainer(trackName))
	{
		AnalysisTask->GetTrackContainer(trackName)->SetParticleEtaLimits(-trackEta,trackEta); //..Eta limits (-0.8,0.8 as in Pi0-h publication)
		AnalysisTask->GetTrackContainer(trackName)->SetPtLimits(trackptcut,30.0);             //..pT limits for accepted tracks
	}
	if(AnalysisTask->GetClusterContainer(clusName))
	{
		AnalysisTask->GetClusterContainer(clusName)->SetClusECut(0);                 //by default set to 0
		AnalysisTask->GetClusterContainer(clusName)->SetClusPtCut(0);                //by default set to 0.15
		AnalysisTask->GetClusterContainer(clusName)->SetClusUserDefEnergyCut(AliVCluster::kNonLinCorr,clusEcut);
		AnalysisTask->GetClusterContainer(clusName)->SetDefaultClusterEnergy(AliVCluster::kNonLinCorr);
		AnalysisTask->GetClusterContainer(clusName)->SetEtaLimits(-clusterEta,clusterEta);
		AnalysisTask->GetClusterContainer(clusName)->SetClusTimeCut(-50e-9,50e9);
	}

	//..some additional input for the analysis
	AnalysisTask->SetSavePool(SavePool);
	AnalysisTask->SetEvtTriggerType(evtTriggerType);   //..Trigger to be used for filling same event histograms
	AnalysisTask->SetEvtMixType(evtMixingType);        //..Trigger to be used to fill tracks into the pool (no GA trigger!!)
	AnalysisTask->SetNLM(1);                           //..Maximum of number of local maxima
	if(InputGammaOrPi0==0)
	{
		AnalysisTask->SetM02(0.1,0.4);                 //..Ranges of allowed cluster shapes in the analysis
	}

	//for later AnalysisTask->SetEffHistGamma(THnF *h);
	//for later AnalysisTask->SetEffHistHadron(THnF *h);

	//-------------------------------------------------------
	// Final settings, pass to manager and set the containers
	//-------------------------------------------------------
	mgr->AddTask(AnalysisTask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(),TList::Class(),
			AliAnalysisManager::kOutputContainer,
			Form("%s", AliAnalysisManager::GetCommonFileName()));
	mgr->ConnectInput  (AnalysisTask, 0,  cinput1 );
	mgr->ConnectOutput (AnalysisTask, 1, coutput1 );

	return AnalysisTask;
}
//the end
