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

/* AliAnaysisTaskFlowTPCEMCalRun2
 *
 */
#include <TGrid.h>
#include <TSpline.h>
#include <TString.h>
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TList.h"
#include "THnSparse.h"
#include "AliMultSelection.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliAODPid.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliKFParticle.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"
//#include "AliAnalysisTaskQnVectorAnalysis.h"
//#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsManager.h"
#include "AliVertexingHFUtils.h"
//#include "AliQnCorrectionsQnVector.h"
//#include "AliAnalysisTaskQnVectorAnalysis.h"
#include "AliAnalysisTaskSEHFTenderQnVectors.h"
#include "AliHFQnVectorHandler.h"
//#include "AliAODHandler.h"
//#include "AliAODVertex.h"

#include "AliAnalysisTaskFlowTPCEMCalRun2.h"

//using std::cout;
//using std::endl;

//class AliAnalysisTaskFlowTPCEMCalRun2;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskFlowTPCEMCalRun2) // classimp: necessary for root

//_____________________________________________________________________________
AliAnalysisTaskFlowTPCEMCalRun2::AliAnalysisTaskFlowTPCEMCalRun2(const char *name) : AliAnalysisTaskSE(name),
	fAOD(0),
	fOutputList(0),
	fVevent(0),
	fpidResponse(0),
	fHistPt(0),
	fNevents(0),
	fCent(0),
	fVtxZ(0),
	fVtxX(0),
	fVtxY(0),
	fTrkPt(0),
        fTrkPtbef(0),
	fTrketa(0),
	fTrkphi(0),
	fTrketa2(0), //1GeVでのカットを加えた
	fdEdx(0),
	fTrkP(0),
	fHistClustE(0),
	fEMCClsEtaPhi(0),
	fHistNoCells(0),
	fHistCalCell(0),
	fHistNCls(0),
	fHistNClsE1(0),
	fHistNClsE2(0),
	fHistNClsE3(0),
	fEMCTrkPt(0),
	fEMCTrketa(0),
	fEMCTrkphi(0),
	fClsEtaPhiAftMatch(0),
	fTPCnsig(0),
        fTOFnsig(0),
        fITSnsig(0),
        fTPCnsig_TOFnsig(0),
        fHistele_TOFcuts(0),
        fHisthad_TOFcuts(0),
	fHisteop(0),
	fM20(0),
	fM02(0),
	fHistNsigEop(0),
	fEMCTrkMatchPhi(0),
	fEMCTrkMatchEta(0),
	fHistelectron(0),
	fHisthadron(0),
	fInvmassLS(0),
	fInvmassULS(0),
	fInvmassLS_2D(0),
	fInvmassULS_2D(0),
	//fUseTender(kTRUE),
	fTracks_tender(0),
        fCaloClusters_tender(0),
	fMCparticle(0),
	fMCcheckMother(0),
	fMCarray(0),
	fMCTrackpart(0),
	fMCheader(0),
	fCheckEtaMC(0),
	NembMCeta(0),
	NembMCpi0(0),
	Nch(0),
	NpureMCproc(0),
	NpureMC(0),
	fHistMCorgPi0(0),
	fHistMCorgEta(0),
	fHistMCorgB(0),
	fHistMCorgD(0),
	fPt_Btoe(0),
	fNDB(0),
	fHist_eff_HFE(0),
	fHist_eff_TPC(0),
	fHistPhoReco0(0),
	fHistPhoReco1(0),
	fHistPhoReco2(0),
        fPi010_0(0),
        fPi010_1(0),
        fEta010(0),
	fPi3050_0(0),
	fPi3050_1(0),
	fEta3050(0),
	fHistPhoPi0(0),
	fHistPhoEta(0),
	fHistPhoPi01(0),
	fHistPhoEta1(0),
	flPercentile(0),
	//fFlowQnVectorMgr(0),
	//flowQnVectorTask(0),
	//fmyEventPlane(0),
	//fmyEventPlane2(0),
	//fmyEventPlane3(0),
	fEPcorV0AC(0),
	fEPcorV0ATPC(0),
	fEPcorV0CTPC(0),
	//fTrkPhiEPFullTPC(0),
	//fTrkPhiEPFullV0A(0),
	//fTrkPhiEPFullTPC_Pt(0),
	fTrkPhiEPV0A_Pt(0),
	fTrkPhiEPV0A_Pt_ele(0),
        fTrkPhiEPV0A_Pt_ele_lowpt(0),
	//fTrkPhiEP2(0),
	//fTrkPhiEP_Pt(0),
	//fTrkPhiEP2_Pt(0),
	fTrkPhicos2(0),
        fTrkPhisin2(0),
        fTrkPhicos2_elelow(0),
        fTrkPhisin2_elelow(0),
        fTrkPhicos2_elehigh(0),
        fTrkPhisin2_elehigh(0),
        fTrkPhicos2_hfehigh(0),
        fTrkPhisin2_hfehigh(0),
        fTrkPhicos2_hadhigh(0),
        fTrkPhisin2_hadhigh(0),
        fTrkPhicos2_phoLShigh(0),
        fTrkPhisin2_phoLShigh(0),
        fTrkPhicos2_phoULShigh(0),
        fTrkPhisin2_phoULShigh(0),
	//fInplane(0),
	//fOutplane(0),
        fInplane_ele(0),
        fOutplane_ele(0),
        fInplane_hfe(0),
        fOutplane_hfe(0),
        fInplane_LSpho(0),
        fOutplane_LSpho(0),
        fInplane_ULSpho(0),
        fOutplane_ULSpho(0),
	DCAxy(3.0),
	DCAz(3.0),
	fDCAxy_Pt_ele(0),
	fDCAxy_Pt_had(0),
	fDCAxy_Pt_hfe(0),
	fDCAxy_Pt_Inplane_ele(0),
	fDCAxy_Pt_Outplane_ele(0),
	fDCAxy_Pt_Inplane_hfe(0),
	fDCAxy_Pt_Outplane_hfe(0),
	fDCAxy_Pt_Inplane_had(0),
	fDCAxy_Pt_Outplane_had(0),
	fHistPt_HFE_MC_D(0),
	fHistPt_HFE_MC_B(0),
	fDCAxy_Pt_D(0),
	fDCAxy_Pt_D_WeightNew(0),
	fDCAxy_Pt_D_WeightVar1(0),
	fDCAxy_Pt_D_WeightVar2(0),
	fDCAxy_Pt_Dpm(0),
	fDCAxy_Pt_Dpm_WeightNew(0),
	fDCAxy_Pt_Dpm_WeightVar1(0),
	fDCAxy_Pt_Dpm_WeightVar2(0),
	fDCAxy_Pt_D0(0),
	fDCAxy_Pt_D0_WeightNew(0),
	//fDCAxy_Pt_D0_WeightVar1(0),
	//fDCAxy_Pt_D0_WeightVar2(0),
	fDCAxy_Pt_Ds(0),
	fDCAxy_Pt_Ds_WeightNew(0),
	fDCAxy_Pt_Ds_WeightVar1(0),
	fDCAxy_Pt_Ds_WeightVar2(0),
	fDCAxy_Pt_lambda(0),
	fDCAxy_Pt_lambda_WeightNew(0),
	fDCAxy_Pt_lambda_WeightVar1(0),
	fDCAxy_Pt_lambda_WeightVar2(0),
	fDCAxy_Pt_Bmeson(0),
	fDCAxy_Pt_Bmeson_WeightNew(0),
	fDCAxy_Pt_Bmeson_WeightVar1(0),
	fDCAxy_Pt_Bmeson_WeightVar2(0),
	fDCAxy_Pt_Bbaryon(0),
	fDCAxy_Pt_Bbaryon_WeightNew(0),
	fDCAxy_Pt_Bbaryon_WeightVar1(0),
	fDCAxy_Pt_Bbaryon_WeightVar2(0),
        ftpcnsig(-1.0),
        femceop(0.9),
        femcss_mim(0.01),
        femcss_max(0.35),
        finvmass(0.1),
        finvmass_pt(0.15),
	massMin(0.1),
	fDCAxy_Pt_LS(0),
	fDCAxy_Pt_ULS(0),
	fDCAxy_Pt_Inplane_LS(0),
	fDCAxy_Pt_Outplane_LS(0),
	fDCAxy_Pt_Inplane_ULS(0),
	fDCAxy_Pt_Outplane_ULS(0),
	fsubV0ACcos2(0),
	fsubV0ATPCcos2(0),
	fsubV0CTPCcos2(0),
	fcorTrkPhicent_charge(0),
	fcorTrkPhicent_ele(0),
	fcorcentcos2_charge(0),
	fcorcentInplane(0),
	fcorcentOutplane(0),
	//fQx1(0),
	//fQy1(0)
	//fQx2(0),
	//fQy2(0),
	//fQx3(0),
	//fQy3(0),
	fPercentileqn(false),
	fLoadedSplines(false),
	fHarmonic(2),
	//fCalibType(AliHFQnVectorHandler::kQnFrameworkCalib),
	fCalibType(AliHFQnVectorHandler::kQnCalib),
	fNormMethod(AliHFQnVectorHandler::kQoverM),
	//fOADBFileName("alien:////alice/cern.ch/user/f/fgrosa/QnVectorCalibrations/calibV0TrklTPCNoEtaCutRun218rVtx14MRP2New.root"),
	fOADBFileName(""),
	fFlowMethod(kEvShapeEP),
	fEvPlaneDet(kFullV0),
	fSubEvDetA(kPosTPC),
	fSubEvDetB(kNegTPC),
	//fRDCuts(nullptr),
	fqnMeth(kq2TPC),
	fTenderTaskName("HFTenderQnVectorsV2"),
	//fqnSplineFileName("alien:///alice/cern.ch/user/f/fgrosa/q2Splines/Splines_q2_3050_1centbin_LHC18r.root"),
	fqnSplineFileName(""),
	fEtaGapInTPCHalves(0),
	fScalProdLimit(0.4),
	fMinCentr(30.),
	fMaxCentr(50.),
        iCentral(kFALSE),
        iSemiCentral(kTRUE),
        iBevt(kFALSE),
	fDWeightNew(0),//weight pT at D mesons and B mesons
	fDWeightVar1(0),
	fDWeightVar2(0),
	fDPlusWeightVar1(0),
	fDsWeightVar1(0),
	fLcWeightVar1(0),
	fLcWeightVar2(0),
	//fBWeight(0),
	fBWeightNew(0),
	fBWeightVar1(0),
	fBWeightVar2(0),
	fSparseElectron(0),
	fvalueElectron(0),
        iTree(kFALSE)


{
	// Default constructor
for(int iHisto=0; iHisto<3; iHisto++){

	fHistEPResolVsCentrVsqn[iHisto] = nullptr;

}


	for(int iDet=0; iDet<6; iDet++){

		fqnSplinesList[iDet] = nullptr;
		fHistqnVsCentrPercCalib[iDet] = nullptr;
		fHistEvPlaneQncorr[iDet] = nullptr;
		fHistEvPlane[iDet] = nullptr;
	}


        fvalueElectron = new Double_t[6];

	DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
	DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms

}


//_____________________________________________________________________________
AliAnalysisTaskFlowTPCEMCalRun2::AliAnalysisTaskFlowTPCEMCalRun2() : AliAnalysisTaskSE("HFEflowRun2"),
	fAOD(0),
	fOutputList(0),
	fVevent(0),
	fpidResponse(0),
	fHistPt(0),
	fNevents(0),
	fCent(0),
	fVtxZ(0),
	fVtxX(0),
	fVtxY(0),
	fTrkPt(0),
        fTrkPtbef(0),
	fTrketa(0),
	fTrkphi(0),
	fTrketa2(0), //1GeVでのカットを加えた
	fdEdx(0),
	fTrkP(0),
	fHistClustE(0),
	fEMCClsEtaPhi(0),
	fHistNoCells(0),
	fHistCalCell(0),
	fHistNCls(0),
	fHistNClsE1(0),
	fHistNClsE2(0),
	fHistNClsE3(0),
	fEMCTrkPt(0),
	fEMCTrketa(0),
	fEMCTrkphi(0),
	fClsEtaPhiAftMatch(0),
	fTPCnsig(0),
        fTOFnsig(0),
        fITSnsig(0),
        fTPCnsig_TOFnsig(0),
        //fTrkPt_TPCnsig_TOFnsig(0),
        fHistele_TOFcuts(0),
        fHisthad_TOFcuts(0),
	fHisteop(0),
	fM20(0),
	fM02(0),
	fHistNsigEop(0),
	fEMCTrkMatchPhi(0),
	fEMCTrkMatchEta(0),
	fHistelectron(0),
	fHisthadron(0),
	fInvmassLS(0),
	fInvmassULS(0),
	fInvmassLS_2D(0),
	fInvmassULS_2D(0),
	//fUseTender(kTRUE),
	fTracks_tender(0),
        fCaloClusters_tender(0),
	fMCparticle(0),
	fMCcheckMother(0),
	fMCarray(0),
	fMCTrackpart(0),
	fMCheader(0),
	fCheckEtaMC(0),
	NembMCeta(0),
	NembMCpi0(0),
	Nch(0),
	NpureMCproc(0),
	NpureMC(0),
	fHistMCorgPi0(0),
	fHistMCorgEta(0),
	fHistMCorgD(0),
	fHistMCorgB(0),
	fPt_Btoe(0),
	fNDB(0),
	fHist_eff_HFE(0),
	fHist_eff_TPC(0),
	fHistPhoReco0(0),
	fHistPhoReco1(0),
	fHistPhoReco2(0),
        fPi010_0(0),
        fPi010_1(0),
        fEta010(0),
	fPi3050_0(0),
	fPi3050_1(0),
	fEta3050(0),
	fHistPhoPi0(0),
	fHistPhoEta(0),
	fHistPhoPi01(0),
	fHistPhoEta1(0),
	flPercentile(0),
	//fFlowQnVectorMgr(0),
	//flowQnVectorTask(0),
	//fmyEventPlane(0),
	//fmyEventPlane2(0),
	//fmyEventPlane3(0),
	fEPcorV0AC(0),
	fEPcorV0ATPC(0),
	fEPcorV0CTPC(0),
	//fTrkPhiEPFullTPC(0),
	//fTrkPhiEPFullV0A(0),
	//fTrkPhiEPFullTPC_Pt(0),
	fTrkPhiEPV0A_Pt(0),
	fTrkPhiEPV0A_Pt_ele(0),
        fTrkPhiEPV0A_Pt_ele_lowpt(0),
	//fTrkPhiEP2(0),
	//fTrkPhiEP_Pt(0),
	//fTrkPhiEP2_Pt(0),
	fTrkPhicos2(0),
        fTrkPhisin2(0),
        fTrkPhicos2_elelow(0),
        fTrkPhisin2_elelow(0),
        fTrkPhicos2_elehigh(0),
        fTrkPhisin2_elehigh(0),
        fTrkPhicos2_hfehigh(0),
        fTrkPhisin2_hfehigh(0),
        fTrkPhicos2_hadhigh(0),
        fTrkPhisin2_hadhigh(0),
        fTrkPhicos2_phoLShigh(0),
        fTrkPhisin2_phoLShigh(0),
        fTrkPhicos2_phoULShigh(0),
        fTrkPhisin2_phoULShigh(0),
	//fInplane(0),
	//fOutplane(0),
        fInplane_ele(0),
        fOutplane_ele(0),
        fInplane_hfe(0),
        fOutplane_hfe(0),
        fInplane_LSpho(0),
        fOutplane_LSpho(0),
        fInplane_ULSpho(0),
        fOutplane_ULSpho(0),
	DCAxy(3.0),
	DCAz(3.0),
	fDCAxy_Pt_ele(0),
	fDCAxy_Pt_had(0),
	fDCAxy_Pt_hfe(0),
	fDCAxy_Pt_Inplane_ele(0),
	fDCAxy_Pt_Outplane_ele(0),
	fDCAxy_Pt_Inplane_hfe(0),
	fDCAxy_Pt_Outplane_hfe(0),
	fDCAxy_Pt_Inplane_had(0),
	fDCAxy_Pt_Outplane_had(0),
	fHistPt_HFE_MC_D(0),
	fHistPt_HFE_MC_B(0),
	fDCAxy_Pt_D(0),
	fDCAxy_Pt_D_WeightNew(0),
	fDCAxy_Pt_D_WeightVar1(0),
	fDCAxy_Pt_D_WeightVar2(0),
	fDCAxy_Pt_Dpm(0),
	fDCAxy_Pt_Dpm_WeightNew(0),
	fDCAxy_Pt_Dpm_WeightVar1(0),
	fDCAxy_Pt_Dpm_WeightVar2(0),
	fDCAxy_Pt_D0(0),
	fDCAxy_Pt_D0_WeightNew(0),
	//fDCAxy_Pt_D0_WeightVar1(0),
	//fDCAxy_Pt_D0_WeightVar2(0),
	fDCAxy_Pt_Ds(0),
	fDCAxy_Pt_Ds_WeightNew(0),
	fDCAxy_Pt_Ds_WeightVar1(0),
	fDCAxy_Pt_Ds_WeightVar2(0),
	fDCAxy_Pt_lambda(0),
	fDCAxy_Pt_lambda_WeightNew(0),
	fDCAxy_Pt_lambda_WeightVar1(0),
	fDCAxy_Pt_lambda_WeightVar2(0),
	fDCAxy_Pt_Bmeson(0),
	fDCAxy_Pt_Bmeson_WeightNew(0),
	fDCAxy_Pt_Bmeson_WeightVar1(0),
	fDCAxy_Pt_Bmeson_WeightVar2(0),
	fDCAxy_Pt_Bbaryon(0),
	fDCAxy_Pt_Bbaryon_WeightNew(0),
	fDCAxy_Pt_Bbaryon_WeightVar1(0),
	fDCAxy_Pt_Bbaryon_WeightVar2(0),
        ftpcnsig(-1.0),
        femceop(0.9),
        femcss_mim(0.01),
        femcss_max(0.35),
        finvmass(0.1),
        finvmass_pt(0.15),
	massMin(0.1),
	fDCAxy_Pt_LS(0),
	fDCAxy_Pt_ULS(0),
	fDCAxy_Pt_Inplane_LS(0),
	fDCAxy_Pt_Outplane_LS(0),
	fDCAxy_Pt_Inplane_ULS(0),
	fDCAxy_Pt_Outplane_ULS(0),
	fsubV0ACcos2(0),
	fsubV0ATPCcos2(0),
	fsubV0CTPCcos2(0),
	fcorTrkPhicent_charge(0),
	fcorTrkPhicent_ele(0),
	fcorcentcos2_charge(0),
	fcorcentInplane(0),
	fcorcentOutplane(0),
	//fQx1(0),
	//fQy1(0)
	//fQx2(0),
	//fQy2(0),
	//fQx3(0),
	//fQy3(0),
	fPercentileqn(false),
	fLoadedSplines(false),
	fHarmonic(2),
	//fCalibType(AliHFQnVectorHandler::kQnFrameworkCalib),
	fCalibType(AliHFQnVectorHandler::kQnCalib),
	fNormMethod(AliHFQnVectorHandler::kQoverM),
	//fOADBFileName("alien:////alice/cern.ch/user/f/fgrosa/QnVectorCalibrations/calibV0TrklTPCNoEtaCutRun218rVtx14MRP2New.root"),
	fOADBFileName(""),
	fFlowMethod(kEvShapeEP),
	fEvPlaneDet(kFullV0),
	fSubEvDetA(kPosTPC),
	fSubEvDetB(kNegTPC),
	//fRDCuts(rdCuts),
	fqnMeth(kq2TPC),
	fTenderTaskName("HFTenderQnVectorsV2"),
	//fTenderTaskName(""),
	//fqnSplineFileName("alien:///alice/cern.ch/user/f/fgrosa/q2Splines/Splines_q2_3050_1centbin_LHC18r.root"),
	fqnSplineFileName(""),
	fEtaGapInTPCHalves(0),
	fScalProdLimit(0.4),
	fMinCentr(30.),
	fMaxCentr(50.),
        iCentral(kFALSE),
        iSemiCentral(kTRUE),
        iBevt(kFALSE),
	fDWeightNew(0),//weight pT at D mesons and B mesons
	fDWeightVar1(0),
	fDWeightVar2(0),
	fDPlusWeightVar1(0),
	fDsWeightVar1(0),
	fLcWeightVar1(0),
	fLcWeightVar2(0),
	//fBWeight(0),
	fBWeightNew(0),
	fBWeightVar1(0),
	fBWeightVar2(0),
	fSparseElectron(0),
	fvalueElectron(0),
        iTree(kFALSE)

	// Standard constructor
{

	for(int iHisto=0; iHisto<3; iHisto++){

		fHistEPResolVsCentrVsqn[iHisto] = nullptr;

	}

	for(int iDet=0; iDet<6; iDet++){

		fqnSplinesList[iDet] = nullptr;
		fHistqnVsCentrPercCalib[iDet] = nullptr;
		fHistEvPlaneQncorr[iDet] = nullptr;
		fHistEvPlane[iDet] = nullptr;

	}


        fvalueElectron = new Double_t[6];

	// default constructor, don't allocate memory here!
	// this is used by root for IO purposes, it needs to remain empty
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
AliAnalysisTaskFlowTPCEMCalRun2::~AliAnalysisTaskFlowTPCEMCalRun2()
{
	// destructor
  //   if(fOutputList) {
  // 	    delete fOutputList;
  // 	//    delete fTracks_tender;
  // // at the end of your task, it is deleted from memory by calling this function
  //   }

	for(int iDet=0; iDet<6; iDet++){

		delete fHistEvPlaneQncorr[iDet];
		delete fHistEvPlane[iDet];
		delete fHistqnVsCentrPercCalib[iDet];
                delete []fvalueElectron;

	}

	for(int iHisto=0; iHisto<3; iHisto++){

		delete fHistEPResolVsCentrVsqn[iHisto];

	}

	//delete fRDCuts;

	for(int iDet=0; iDet<6;iDet++){
		if(fqnSplinesList[iDet] && fLoadedSplines) delete fqnSplinesList[iDet];
	}


}
//____________________________________________________________________________
//void AliAnalysisTaskFlowTPCEMCalRun2::LocalInit()
//{

//Initialization
//fRDCuts->SetMinCentrality(fMinCentr);
//fRDCuts->SetMaxCentrality(fMaxCentr);

//}

//_____________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalRun2::UserCreateOutputObjects()
{
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

    // example of a histogram
    fHistPt = new TH1F("fHistPt", "fHistPt", 100, 0, 10);       // create your histogra
    fOutputList->Add(fHistPt);          // don't forget to add it to the list! the list will be written to file, so if you want
                                        // your histogram in the output file, add it to the list!

//No of events
fNevents = new TH1F("fNevents","No of events",3,-0.5,2.5);
fOutputList->Add(fNevents);
fNevents->GetYaxis()->SetTitle("counts");
//fNevents->GetXaxis()->SetBinlabel(1,"All");
//fNevents->GetXaxis()->SetBinlabel(2,"With >2 Trks");
//fNevents->GetXaxis()->SetBinlabel(3,"Vtx_{z}<10cm");

fCent = new TH1F("fCent","Centrality",100,0,100);
fOutputList->Add(fCent);

fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",500,-25,25);
fOutputList->Add(fVtxZ);

fVtxX = new TH1F("fVtxX","X vertex position;Vtx_{x};counts",500,-25,25);
fOutputList->Add(fVtxX);

fVtxY = new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",500,-25,25);
fOutputList->Add(fVtxY);

fTrkPt = new TH1F("fTrkPt","p_{T} distribution of all tracks;p_{T} (GeV/c);counts",500,0,100);
fOutputList->Add(fTrkPt);

fTrkPtbef = new TH1F("fTrkPtbef","p_{T} distribution of all tracks before DCA cuts;p_{T} (GeV/c);counts",500,0,100);
fOutputList->Add(fTrkPtbef);

fTrketa = new TH1F("fTrketa","All Track #eta distribution;#eta;counts",100,-1.5,1.5);
fOutputList->Add(fTrketa);

fTrkphi = new TH1F("fTrkphi","All Track #phi distribution;#phi;counts",100,0,6.3);
fOutputList->Add(fTrkphi);

fTrketa2 = new TH1F("fTrketa2","All Track #eta 2 distribution;#eta 2;counts",100,-1.5,1.5);
fOutputList->Add(fTrketa2);

fTrkP = new TH1F("fTrkP","All Track P distribution;P;counts",500,0,50);
fOutputList->Add(fTrkP);

fdEdx = new TH2F("fdEdx","All Track dE/dx distribution;p (GeV/c);de/dx",150,0,15,500,0,160);
fOutputList->Add(fdEdx);

fHistClustE = new TH1F("fHistClustE","EMCAL cluster energy distribution; Cluster E;counts",500,0.0,50.0);
fOutputList->Add(fHistClustE);

fEMCClsEtaPhi = new TH2F("fEMCClsEtaPhi","EMCAL cluster #eta and #phi distribution;#eta;#phi",100,-0.9,0.9,200,0,6.3);
fOutputList->Add(fEMCClsEtaPhi);

fHistNoCells = new TH2F("fHistNoCells","Number of EMCAL cells in a cluster;Cluster E;N^{EMC}_{cells}",500,0,50,30,0,30);
fOutputList->Add(fHistNoCells);

fHistCalCell = new TH2F("fHistCalCell","Energy of EMCAL cells;cell ID;E (GeV)",15000,-0.5,14999.5,150,0,30);
fOutputList->Add(fHistCalCell);

fHistNCls = new TH1F("fHistNCls","Number of EMCAL cluster in the event;N^{EMC}_{cls};counts",100,0,100);
fOutputList->Add(fHistNCls);

fHistNClsE1 = new TH1F("fHistNClsE1","Number of EMCAL cluster in the event (E>0.1 GeV);N^{EMC}_{cls};counts",100,0,100);
fOutputList->Add(fHistNClsE1);

fHistNClsE2 = new TH1F("fHistNClsE2","Number of EMCAL cluster in the event (E>0.2 GeV);N^{EMC}_{cls};counts",100,0,100);
fOutputList->Add(fHistNClsE2);

fHistNClsE3 = new TH1F("fHistNClsE3","Number of EMCAL cluster in the event (E>0.5 GeV);N^{EMC}_{cls};counts",100,0,100);
fOutputList->Add(fHistNClsE3);

fEMCTrkPt = new TH1F("fEMCTrkPt","p_{T} distribution of tracks with EMCAL cluster;p_{T};counts",500,0,50);
fOutputList->Add(fEMCTrkPt);

fEMCTrketa = new TH1F("fEMCTrketa","#eta distribution of tracks matched EMCAL cluster;#eta;counts",60,-1.5,1.5);
fOutputList->Add(fEMCTrketa);

fEMCTrkphi = new TH1F("fEMCTrkphi","#phi distribution of tracks matched EMCAL cluster;#phi;counts",100,0,6.3);
fOutputList->Add(fEMCTrkphi);

fClsEtaPhiAftMatch = new TH2F("fClsEtaPhiAftMatch","EMCAL cluster #eta and #phi distribution after track matching;#eta;#phi",100,-0.9,0.9,200,0,6.3);
fOutputList->Add(fClsEtaPhiAftMatch);

fTPCnsig = new TH2F("fTPCnsig","All Track TPC Nsigma distribution;p (GeV/c);#sigma_{TPC-dE/dx}",500,0,50,200,-10,10);
fOutputList->Add(fTPCnsig);

fTOFnsig = new TH2F("fTOFnsig","All Track TOF Nsigma distribution;p (GeV/c);#sigma_{TOF-dE/dx}",500,0,50,200,-10,10);
fOutputList->Add(fTOFnsig);

fITSnsig = new TH2F("fITSnsig","All Track ITS Nsigma distribution;p (GeV/c);#sigma_{ITS-dE/dx}",500,0,50,200,-10,10);
fOutputList->Add(fITSnsig);

fTPCnsig_TOFnsig = new TH2F("fTPCnsig_TOFnsig","All Track TOF vs TPC Nsigma distribution;#sigma {TOF-dE/dx};#sigma_{TPC-dE/dx}",200,-10,10,200,-10,10);
fOutputList->Add(fTPCnsig_TOFnsig);

fHistele_TOFcuts = new TH2F("fHistele_TOFcuts","TPC Nsigma vs #it{p}_{T} with -3 < TOF Nsigma < 3 ;#sigma {TPC-dE/dx};#it{p}_{T}",200,-10,10,500,0,50);
fOutputList->Add(fHistele_TOFcuts);

fHisthad_TOFcuts = new TH2F("fHisthad_TOFcuts","TPC Nsigma vs #it{p}_{T} with TOF Nsigma < -3.5;#sigma {TPC-dE/dx};#it{p}_{T}",200,-10,10,500,0,50);
fOutputList->Add(fHisthad_TOFcuts);

fHisteop = new TH1F("fHisteop","E/p with -1<n#sigma<3;E/p;counts",40,0,2);
fOutputList->Add(fHisteop);

fM20 = new TH2F("fM20","M20 vs Pt distribution",500,0,50,200,0,2);
fOutputList->Add(fM20);

fM02 = new TH2F("fM02","M02 vs Pt distribution",500,0,50,200,0,2);
fOutputList->Add(fM02);

fHistNsigEop = new TH2F("fHistNsigEop","TPCnsigma vs E/p, p_{T}>2 GeV/c;E/p;#sigma_{TPC-dE/dx}",40,0.0,2.0,200,-10,10);
fOutputList->Add(fHistNsigEop);

fEMCTrkMatchPhi = new TH1F("fEMCTrkMatchPhi","Distance of EMCAL cluster to its closest track;#delta#phi;counts",100,-0.3,0.3);
fOutputList->Add(fEMCTrkMatchPhi);

fEMCTrkMatchEta = new TH1F("fEMCTrkMatchEta","Distance of EMCAL cluster to its closest track;#delta#eta;counts",100,-0.3,0.3);
fOutputList->Add(fEMCTrkMatchEta);

fHistelectron = new TH2F("fHistelectron","P_{T} vs E/P with -1<n#sigma<3;E/P;P_{T}",40,0,2,100,0,20);
fOutputList->Add(fHistelectron);

fHisthadron = new TH2F("fHisthadron","P_{T} vs E/P with n#sigma<-3;E/P;P_{T}",40,0,2,100,0,20);
fOutputList->Add(fHisthadron);

fInvmassLS = new TH1F("fInvmassLS","Invmass of LS (e,e) for pt^{e}>1;mass(GeV/c^2);counts",500,0,1.0);
fOutputList->Add(fInvmassLS);

fInvmassULS = new TH1F("fInvmassULS","Invmass of ULS (e,e) for py^{e}>1;mass(GeV/c^2);counts",500,0,1.0);
fOutputList->Add(fInvmassULS);

fInvmassLS_2D = new TH2F("fInvmassLS_2D","Invmass of LS vs P_{T};mass(GeV/c^2);P_{T} (GeV/c)",500,0,1.0,400,0.0,20.0);
fOutputList->Add(fInvmassLS_2D);

fInvmassULS_2D = new TH2F("fInvmassULS_2D","Invmass of ULS vs P_{T};mass(GeV/c^2);P_{T} (GeV/c)",500,0,1.0,400,0.0,20.0);
fOutputList->Add(fInvmassULS_2D);

fCheckEtaMC = new TH1F("fCheckEtaMC","check Eta range cut in MC",160,-0.8,0.8);
fOutputList->Add(fCheckEtaMC);

fHistMCorgPi0 = new TH2F("fHistMCorgPi0","MC org Pi0",2,-0.5,1.5,100,0,50);
fOutputList->Add(fHistMCorgPi0);

fHistMCorgEta = new TH2F("fHistMCorgEta","MC org Eta",2,-0.5,1.5,100,0,50);
fOutputList->Add(fHistMCorgEta);

fHistMCorgB = new TH1F("fHistMCorgB","MC org B",600,0,60);
fOutputList->Add(fHistMCorgB);

fHistMCorgD = new TH1F("fHistMCorgD","MC org D",600,0,60);
fOutputList->Add(fHistMCorgD);

Double_t eop_range[8] = {2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};
fPt_Btoe = new TH2F("fPt_Btoe","B meson vs electron:electron p_{t} (GeV/c)",20,eop_range,600,0,60);
fOutputList->Add(fPt_Btoe);

fNDB = new TH1F("fNDB","No of events",2,-0.5,1.5);
fOutputList->Add(fNDB);

fHist_eff_HFE = new TH1F("fHist_eff_HFE","efficiency::HFE",600,0,60);
fOutputList->Add(fHist_eff_HFE);

fHist_eff_TPC = new TH1F("fHist_eff_TPC","efficiency::TPC",600,0,60);
fOutputList->Add(fHist_eff_TPC);

fMCcheckMother = new TH1F("fMCcheckMother","Mother MC PDG",1000,-0.5,999.5);
fOutputList->Add(fMCcheckMother);

fHistPhoReco0 = new TH1F("fHistPhoReco0","P_{T} (allPho)",500,0,100);
fOutputList->Add(fHistPhoReco0);

fHistPhoReco1 = new TH1F("fHistPhoReco1","P_{T} (NonHFE)",500,0,100);
fOutputList->Add(fHistPhoReco1);

fHistPhoReco2 = new TH1F("fHistPhoReco2","P_{T} (HFE)",500,0,100);
fOutputList->Add(fHistPhoReco2);

//Pi0 Weight

fPi010_0 = new TF1("fPi010_0","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
fPi010_0->SetParameters(6.77248e+02, 8.59098e-01, 1.77145e-01, 1.45314e+00, 4.75048e+00);
fPi010_1 = new TF1("fPi010_1","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
fPi010_1->SetParameters(1.81210e+05, -3.41909e-01, 2.32972e-02 , 5.46069e-01, 4.66485e+00);

fPi3050_0 = new TF1("fPi3050_0","[0]*x/pow([1]+x/[2],[3])");
fPi3050_0->SetParameters(0.937028,0.674846,9.02659,10.);
fPi3050_1 = new TF1("fPi3050_1","[0]*x/pow([1]+x/[2],[3])");
fPi3050_1->SetParameters(2.7883,0.,2.5684,5.63827);


//Eta Weight

fEta010 = new TF1("fEta010","[0]*x/pow([1]+x/[2]+x*x/[3],[4])");
fEta010->SetParameters(2.50883e-02,-1.63341e+00,6.58911e+00,8.07446e-01,3.12257e+00);

fEta3050 = new TF1("fEta3050","[0]*x/pow([1]+x/[2]+x*x/[3],[4])",0,100);
fEta3050 -> SetParameters(5.87918e+01,2.14009e-01,4.03579e+00,2.38693e+00,2.52382e+00);
//fOutputList -> Add(fEta3050);

fHistPhoPi0 = new TH1F("fHistPhoPi0","total pi0 in sample;p_{T}(GeV/c)",500,0,100);
fHistPhoPi0->Sumw2();
fOutputList->Add(fHistPhoPi0);

fHistPhoEta = new TH1F("fHistPhoEta","total Eta in sample;p_{T}(GeV/c)",500,0,100);
fHistPhoEta->Sumw2();
fOutputList->Add(fHistPhoEta);

fHistPhoPi01 = new TH1F("fHistPhoPi01","reco Pi0 in sample;p_{T}(GeV/c)",500,0,100);
fHistPhoPi01->Sumw2();
fOutputList->Add(fHistPhoPi01);

fHistPhoEta1 = new TH1F("fHistPhoEta1","reco Eta in sample;p_{T}(GeV/c)",500,0,100);
fHistPhoEta1->Sumw2();
fOutputList->Add(fHistPhoEta1);

flPercentile = new TH1F("flPercentile","Centrality",600,0,60);
fOutputList->Add(flPercentile);

//fmyEventPlane = new TH1F("fmyEventPlane","EventPlane",320,-1.6,1.6);
//fOutputList->Add(fmyEventPlane);
//
//fmyEventPlane2 = new TH1F("fmyEventPlane2","EventPlane2",320,-1.6,1.6);
//fOutputList->Add(fmyEventPlane2);
//
//fmyEventPlane3 = new TH1F("fmyEventPlane3","EventPlane3",320,-1.6,1.6);
//fOutputList->Add(fmyEventPlane3);
//
fEPcorV0AC = new TH2F("fEPcorV0AC","EPcorV0AC",320,-1.6,1.6,320,-1.6,1.6);
fOutputList->Add(fEPcorV0AC);
//
fEPcorV0ATPC = new TH2F("fEPcorV0ATPC","EPcorV0ATPC",320,-1.6,1.6,320,-1.6,1.6);
fOutputList->Add(fEPcorV0ATPC);
//
fEPcorV0CTPC = new TH2F("fEPcorV0CTPC","EPcorV0CTPC",320,-1.6,1.6,320,-1.6,1.6);
fOutputList->Add(fEPcorV0CTPC);
//
//fTrkPhiEPFullTPC = new TH1F("fTrkPhiEPFullTPC","TrkPhi-EventPlane(FullTPC)",640,-1.6,1.6);
//fOutputList->Add(fTrkPhiEPFullTPC);

//fTrkPhiEPFullV0A = new TH1F("fTrkPhiEPFullV0A","TrkPhi-EventPlane(V0A)",640,-1.6,1.6);
//fOutputList->Add(fTrkPhiEPFullV0A);

//fTrkPhiEPFullTPC_Pt = new TH2F("fTrkPhiEPFullTPC_Pt","TrkPhiFullTPC_eachPt",500,0,100,640,-1.6,1.6);
//fOutputList->Add(fTrkPhiEPFullTPC_Pt);

fTrkPhiEPV0A_Pt = new TH2F("fTrkPhiEPV0A_Pt","TrkPhiV0A_eachPt",500,0,100,640,-1.6,1.6);
fOutputList->Add(fTrkPhiEPV0A_Pt);

fTrkPhiEPV0A_Pt_ele = new TH2F("fTrkPhiEPV0A_Pt_ele","TrkPhiV0A_eachelePt",500,0,100,640,-1.6,1.6);
fOutputList->Add(fTrkPhiEPV0A_Pt_ele);

fTrkPhiEPV0A_Pt_ele_lowpt = new TH2F("fTrkPhiEPV0A_Pt_ele_lowpt","TrkPhiV0A_eachelePt at lowpt",500,0,100,640,-1.6,1.6);
fOutputList->Add(fTrkPhiEPV0A_Pt_ele_lowpt);

//fTrkPhiEP2 = new TH1F("fTrkPhiEP2","TrkPhi-EventPlane",640,-1.6,1.6);
//fOutputList->Add(fTrkPhiEP2);
//
//fTrkPhiEP_Pt = new TH2F("fTrkPhiEP_Pt","TrkPhi_eachPt",500,0,100,640,-1.6,1.6);
//fOutputList->Add(fTrkPhiEP_Pt);
//
//fTrkPhiEP2_Pt = new TH2F("fTrkPhiEP2_Pt","TrkPhi_eachPt2",500,0,100,640,-1.6,1.6);
//fOutputList->Add(fTrkPhiEP2_Pt);
//
fTrkPhicos2 = new TH2F("fTrkPhicos2","fTrkPhicos2",200,0,20,200,-1,1);
fOutputList->Add(fTrkPhicos2);

fTrkPhisin2 = new TH2F("fTrkPhisin2","fTrkPhisin2",200,0,20,200,-1,1);
fOutputList->Add(fTrkPhisin2);

fTrkPhicos2_elelow = new TH2F("fTrkPhicos2_elelow","fTrkPhicos2_elelow",200,0,20,200,-1,1);
fOutputList->Add(fTrkPhicos2_elelow);
//
fTrkPhisin2_elelow = new TH2F("fTrkPhisin2_elelow","fTrkPhisin2_elelow",200,0,20,200,-1,1);
fOutputList->Add(fTrkPhisin2_elelow);

fTrkPhicos2_elehigh = new TH2F("fTrkPhicos2_elehigh","fTrkPhicos2_elehigh",200,0,20,200,-1,1);
fOutputList->Add(fTrkPhicos2_elehigh);
//
fTrkPhisin2_elehigh = new TH2F("fTrkPhisin2_elehigh","fTrkPhisin2_elehigh",200,0,20,200,-1,1);
fOutputList->Add(fTrkPhisin2_elehigh);

fTrkPhicos2_hfehigh = new TH2F("fTrkPhicos2_hfehigh","fTrkPhicos2_hfehigh",200,0,20,200,-1,1);
fOutputList->Add(fTrkPhicos2_hfehigh);
//
fTrkPhisin2_hfehigh = new TH2F("fTrkPhisin2_hfehigh","fTrkPhisin2_hfehigh",200,0,20,200,-1,1);
fOutputList->Add(fTrkPhisin2_hfehigh);

fTrkPhicos2_hadhigh = new TH2F("fTrkPhicos2_hadhigh","fTrkPhicos2_hadhigh",200,0,20,200,-1,1);
fOutputList->Add(fTrkPhicos2_hadhigh);
//
fTrkPhisin2_hadhigh = new TH2F("fTrkPhisin2_hadhigh","fTrkPhisin2_hadhigh",200,0,20,200,-1,1);
fOutputList->Add(fTrkPhisin2_hadhigh);

fTrkPhicos2_phoLShigh = new TH2F("fTrkPhicos2_phoLShigh","fTrkPhicos2_phoLShigh",200,0,20,200,-1,1);
fOutputList->Add(fTrkPhicos2_phoLShigh);
//
fTrkPhisin2_phoLShigh = new TH2F("fTrkPhisin2_phoLShigh","fTrkPhisin2_phoLShigh",200,0,20,200,-1,1);
fOutputList->Add(fTrkPhisin2_phoLShigh);

fTrkPhicos2_phoULShigh = new TH2F("fTrkPhicos2_phoULShigh","fTrkPhicos2_phoULShigh",200,0,20,200,-1,1);
fOutputList->Add(fTrkPhicos2_phoULShigh);

fTrkPhisin2_phoULShigh = new TH2F("fTrkPhisin2_phoULShigh","fTrkPhisin2_phoULShigh",200,0,20,200,-1,1);
fOutputList->Add(fTrkPhisin2_phoULShigh);

//
//fInplane = new TH1F("fInplane","p_{T} distribution of inplane",200,0,20);
//fOutputList->Add(fInplane);
//
//fOutplane = new TH1F("fOutplane","p_{T} distribution of outplane",200,0,20);
//fOutputList->Add(fOutplane);

fInplane_ele = new TH1F("fInplane_ele","p_{T} distribution of inplane_ele",200,0,20);
fOutputList->Add(fInplane_ele);

fOutplane_ele = new TH1F("fOutplane_ele","p_{T} distribution of outplane_ele",200,0,20);
fOutputList->Add(fOutplane_ele);

fInplane_hfe = new TH1F("fInplane_hfe","p_{T} distribution of inplane_hfe",200,0,20);
fOutputList->Add(fInplane_hfe);

fOutplane_hfe = new TH1F("fOutplane_hfe","p_{T} distribution of outplane_hfe",200,0,20);
fOutputList->Add(fOutplane_hfe);

fInplane_LSpho = new TH1F("fInplane_LSpho","p_{T} distribution of inplane_LSpho",200,0,20);
fOutputList->Add(fInplane_LSpho);

fOutplane_LSpho = new TH1F("fOutplane_LSpho","p_{T} distribution of outplane_LSpho",200,0,20);
fOutputList->Add(fOutplane_LSpho);

fInplane_ULSpho = new TH1F("fInplane_ULSpho","p_{T} distribution of inplane_ULSpho",200,0,20);
fOutputList->Add(fInplane_ULSpho);

fOutplane_ULSpho = new TH1F("fOutplane_ULSpho","p_{T} distribution of outplane_ULSpho",200,0,20);
fOutputList->Add(fOutplane_ULSpho);

Double_t ptbin[14] = {0.0,0.5,1.0,1.5,2,2.5,3,4,6,8,10,12,14,20};

fDCAxy_Pt_ele = new TH2F("fDCAxy_Pt_ele","DCA_{xy} vs Pt (electron);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_ele);

fDCAxy_Pt_had = new TH2F("fDCAxy_Pt_had","DCA_{xy} vs Pt (hadron);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_had);

fDCAxy_Pt_hfe = new TH2F("fDCAxy_Pt_hfe","DCA_{xy} vs Pt (hfe);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_hfe);

fDCAxy_Pt_LS = new TH2F("fDCAxy_Pt_LS","DCA_{xy} vs Pt LS pairs;p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_LS);

fDCAxy_Pt_ULS = new TH2F("fDCAxy_Pt_ULS","DCA_{xy} vs Pt ULS pairs;p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_ULS);

fDCAxy_Pt_Inplane_LS = new TH2F("fDCAxy_Pt_Inplane_LS","DCA_{xy} vs Pt In-plane LS pairs;p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Inplane_LS);

fDCAxy_Pt_Outplane_LS = new TH2F("fDCAxy_Pt_Outplane_LS","DCA_{xy} vs Pt Out-plane LS pairs;p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Outplane_LS);

fDCAxy_Pt_Inplane_ULS = new TH2F("fDCAxy_Pt_Inplane_ULS","DCA_{xy} vs Pt In-plane ULS pairs;p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Inplane_ULS);

fDCAxy_Pt_Outplane_ULS = new TH2F("fDCAxy_Pt_Outplane_ULS","DCA_{xy} vs Pt Out-plane ULS pairs;p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Outplane_ULS);

fDCAxy_Pt_Inplane_ele = new TH2F("fDCAxy_Pt_Inplane_ele","DCA_{xy} vs Pt Inplane electron;p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Inplane_ele);

fDCAxy_Pt_Outplane_ele = new TH2F("fDCAxy_Pt_Outplane_ele","DCA_{xy} vs Pt Outplane electron;p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Outplane_ele);

fDCAxy_Pt_Inplane_hfe = new TH2F("fDCAxy_Pt_Inplane_hfe","DCA_{xy} vs Pt Inplane hfe;p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Inplane_hfe);

fDCAxy_Pt_Outplane_hfe = new TH2F("fDCAxy_Pt_Outplane_hfe","DCA_{xy} vs Pt Outplane hfe;p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Outplane_hfe);

fDCAxy_Pt_Inplane_had = new TH2F("fDCAxy_Pt_Inplane_had","DCA_{xy} vs Pt Inplane hadron;p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Inplane_had);

fDCAxy_Pt_Outplane_had = new TH2F("fDCAxy_Pt_Outplane_had","DCA_{xy} vs Pt Outplane hadron;p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Outplane_had);

fHistPt_HFE_MC_D  = new TH1F("fHistPt_HFE_MC_D","HFE from D MC",13,ptbin);
fOutputList->Add(fHistPt_HFE_MC_D);

fHistPt_HFE_MC_B  = new TH1F("fHistPt_HFE_MC_B","HFE fron B MC",13,ptbin);
fOutputList->Add(fHistPt_HFE_MC_B);

fDCAxy_Pt_D = new TH2F("fDCAxy_Pt_D","DCA_{xy} vs Pt D(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_D);

fDCAxy_Pt_D_WeightNew = new TH2F("fDCAxy_Pt_D_WeightNew","DCA_{xy} vs Pt D(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_D_WeightNew->Sumw2();
fOutputList->Add(fDCAxy_Pt_D_WeightNew);

fDCAxy_Pt_D_WeightVar1 = new TH2F("fDCAxy_Pt_D_WeightVar1","DCA_{xy} vs Pt D(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_D_WeightVar1->Sumw2();
fOutputList->Add(fDCAxy_Pt_D_WeightVar1);

fDCAxy_Pt_D_WeightVar2 = new TH2F("fDCAxy_Pt_D_WeightVar2","DCA_{xy} vs Pt D(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_D_WeightVar2->Sumw2();
fOutputList->Add(fDCAxy_Pt_D_WeightVar2);

fDCAxy_Pt_Dpm = new TH2F("fDCAxy_Pt_Dpm","DCA_{xy} vs Pt D+-(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Dpm);

fDCAxy_Pt_Dpm_WeightNew = new TH2F("fDCAxy_Pt_Dpm_WeightNew","DCA_{xy} vs Pt D+-(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_Dpm_WeightNew->Sumw2();
fOutputList->Add(fDCAxy_Pt_Dpm_WeightNew);

fDCAxy_Pt_Dpm_WeightVar1 = new TH2F("fDCAxy_Pt_Dpm_WeightVar1","DCA_{xy} vs Pt D+-(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_Dpm_WeightVar1->Sumw2();
fOutputList->Add(fDCAxy_Pt_Dpm_WeightVar1);

fDCAxy_Pt_Dpm_WeightVar2 = new TH2F("fDCAxy_Pt_Dpm_WeightVar2","DCA_{xy} vs Pt D+-(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_Dpm_WeightVar2->Sumw2();
fOutputList->Add(fDCAxy_Pt_Dpm_WeightVar2);

fDCAxy_Pt_D0 = new TH2F("fDCAxy_Pt_D0","DCA_{xy} vs Pt D0(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_D0);

fDCAxy_Pt_D0_WeightNew = new TH2F("fDCAxy_Pt_D0_WeightNew","DCA_{xy} vs Pt D0(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_D0_WeightNew->Sumw2();
fOutputList->Add(fDCAxy_Pt_D0_WeightNew);

//fDCAxy_Pt_D0_WeightVar1 = new TH2F("fDCAxy_Pt_D0_WeightVar1","DCA_{xy} vs Pt D0(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
//fDCAxy_Pt_D0_WeightVar1->Sumw2();
//fOutputList->Add(fDCAxy_Pt_D0_WeightVar1);
//
//fDCAxy_Pt_D0_WeightVar2  = new TH2F("fDCAxy_Pt_D0_WeightVar2","DCA_{xy} vs Pt D0(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
//fDCAxy_Pt_D0_WeightVar2->Sumw2();
//fOutputList->Add(fDCAxy_Pt_D0_WeightVar2);

fDCAxy_Pt_Ds = new TH2F("fDCAxy_Pt_Ds","DCA_{xy} vs Pt Ds(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Ds);

fDCAxy_Pt_Ds_WeightNew = new TH2F("fDCAxy_Pt_Ds_WeightNew","DCA_{xy} vs Pt Ds(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_Ds_WeightNew->Sumw2();
fOutputList->Add(fDCAxy_Pt_Ds_WeightNew);

fDCAxy_Pt_Ds_WeightVar1 = new TH2F("fDCAxy_Pt_Ds_WeightVar1","DCA_{xy} vs Pt Ds(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_Ds_WeightVar1->Sumw2();
fOutputList->Add(fDCAxy_Pt_Ds_WeightVar1);

fDCAxy_Pt_Ds_WeightVar2 = new TH2F("fDCAxy_Pt_Ds_WeightVar2","DCA_{xy} vs Pt Ds(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_Ds_WeightVar2->Sumw2();
fOutputList->Add(fDCAxy_Pt_Ds_WeightVar2);

fDCAxy_Pt_lambda = new TH2F("fDCAxy_Pt_lambda","DCA_{xy} vs Pt lambda(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_lambda);

fDCAxy_Pt_lambda_WeightNew = new TH2F("fDCAxy_Pt_lambda_WeightNew","DCA_{xy} vs Pt lambda(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_lambda_WeightNew->Sumw2();
fOutputList->Add(fDCAxy_Pt_lambda_WeightNew);

fDCAxy_Pt_lambda_WeightVar1 = new TH2F("fDCAxy_Pt_lambda_WeightVar1","DCA_{xy} vs Pt lambda(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_lambda_WeightVar1->Sumw2();
fOutputList->Add(fDCAxy_Pt_lambda_WeightVar1);

fDCAxy_Pt_lambda_WeightVar2 = new TH2F("fDCAxy_Pt_lambda_WeightVar2","DCA_{xy} vs Pt lambda(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_lambda_WeightVar2->Sumw2();
fOutputList->Add(fDCAxy_Pt_lambda_WeightVar2);

fDCAxy_Pt_Bmeson = new TH2F("fDCAxy_Pt_Bmeson","DCA_{xy} vs Pt all B meson(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Bmeson);

fDCAxy_Pt_Bmeson_WeightNew = new TH2F("fDCAxy_Pt_Bmeson_WeightNew","DCA_{xy} vs Pt all B meson(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_Bmeson_WeightNew->Sumw2();
fOutputList->Add(fDCAxy_Pt_Bmeson_WeightNew);

fDCAxy_Pt_Bmeson_WeightVar1 = new TH2F("fDCAxy_Pt_Bmeson_WeightVar1","DCA_{xy} vs Pt all B meson(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_Bmeson_WeightVar1->Sumw2();
fOutputList->Add(fDCAxy_Pt_Bmeson_WeightVar1);

fDCAxy_Pt_Bmeson_WeightVar2 = new TH2F("fDCAxy_Pt_Bmeson_WeightVar2","DCA_{xy} vs Pt all B meson(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_Bmeson_WeightVar2->Sumw2();
fOutputList->Add(fDCAxy_Pt_Bmeson_WeightVar2);

fDCAxy_Pt_Bbaryon = new TH2F("fDCAxy_Pt_Bbaryon","DCA_{xy} vs Pt all B baryon(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Bbaryon);

fDCAxy_Pt_Bbaryon_WeightNew = new TH2F("fDCAxy_Pt_Bbaryon_WeightNew","DCA_{xy} vs Pt all B baryon(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_Bbaryon_WeightNew->Sumw2();
fOutputList->Add(fDCAxy_Pt_Bbaryon_WeightNew);

fDCAxy_Pt_Bbaryon_WeightVar1 = new TH2F("fDCAxy_Pt_Bbaryon_WeightVar1","DCA_{xy} vs Pt all B baryon(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_Bbaryon_WeightVar1->Sumw2();
fOutputList->Add(fDCAxy_Pt_Bbaryon_WeightVar1);

fDCAxy_Pt_Bbaryon_WeightVar2 = new TH2F("fDCAxy_Pt_Bbaryon_WeightVar2","DCA_{xy} vs Pt all B baryon(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fDCAxy_Pt_Bbaryon_WeightVar2->Sumw2();
fOutputList->Add(fDCAxy_Pt_Bbaryon_WeightVar2);

fsubV0ACcos2 = new TH2F("fsubV0ACcos2","fsubV0ACcos2 vs cetrality",40,0,80,200,-1,1);
fOutputList->Add(fsubV0ACcos2);

fsubV0ATPCcos2 = new TH2F("fsubV0ATPCcos2","fsubV0ATPCcos2 vs cetrality",40,0,80,200,-1,1);
fOutputList->Add(fsubV0ATPCcos2);

fsubV0CTPCcos2 = new TH2F("fsubV0CTPCcos2","fsubV0CTPCcos2 vs cetrality",40,0,80,200,-1,1);
fOutputList->Add(fsubV0CTPCcos2);

fcorTrkPhicent_charge = new TH2F("fcorTrkPhicent_charge","TrkPhi vs centrality",320,-1.6,1.6,40,0,80);
fOutputList->Add(fcorTrkPhicent_charge);

fcorTrkPhicent_ele = new TH2F("fcorTrkPhicent_ele","TrkPhi vs centrality",320,-1.6,1.6,40,0,80);
fOutputList->Add(fcorTrkPhicent_ele);

fcorcentcos2_charge = new TH2F("fcorcentcos2_charge","centrality vs cos2#delta#phi",40,0,80,200,-1,1);
fOutputList->Add(fcorcentcos2_charge);

fcorcentInplane = new TH2F("fcorcentInplane","centrality vs p_{T} of Inplane",40,0,80,200,0,20);
fOutputList->Add(fcorcentInplane);

fcorcentOutplane = new TH2F("fcorcentOutplane","centrality vs p_{T} of Outplane",40,0,80,200,0,20);
fOutputList->Add(fcorcentOutplane);

//fQx1 = new TH2F("fQx1","centrality vs Qx1",40,0,80,320,-1.6,1.6);
//fOutputList->Add(fQx1);

//fQy1 = new TH2F("fQy1","centrality vs Qy1",40,0,80,320,-1.6,1.6);
//fOutputList->Add(fQy1);
//
//fQx2 = new TH2F("fQx2","centrality vs Qx2",40,0,80,320,-1.6,1.6);
//fOutputList->Add(fQx2);
//
//fQy2 = new TH2F("fQy2","centrality vs Qy2",40,0,80,320,-1.6,1.6);
//fOutputList->Add(fQy2);
//
//fQx3 = new TH2F("fQx3","centrality vs Qx3",40,0,80,320,-1.6,1.6);
//fOutputList->Add(fQx3);
//
//fQy3 = new TH2F("fQy3","centrality vs Qy3",40,0,80,320,-1.6,1.6);
//fOutputList->Add(fQy3);


//D Meson pt weighting
    Int_t nbins = 13;
    Double_t xbins[14] = {1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,16.,24.,36.,50.};
    //Double_t err[13] = {};
    //fDWeight = new TH1F("fDWeight","D^{0}_data/AllD_MC;p_{T} (GeV/c);Weight;",nbins,xbins);
    fDWeightNew = new TH1F("fDWeightNew","D^{0}_data/AllD_MCNew;p_{T} (GeV/c);Weight;",nbins,xbins);
    fDWeightVar1 = new TH1F("fDWeightVar1","D^{0}_data/AllD_MC;p_{T} (GeV/c);Weight Var1;",nbins,xbins);
    fDWeightVar2 = new TH1F("fDWeightVar2","D^{0}_data/AllD_MC;p_{T} (GeV/c);Weight Var2;",nbins,xbins);
    
    fDPlusWeightVar1 = new TH1F("fDPlusWeightVar1","(D^{+}/D^{0})_data*(D^{0}_data/D^{+}_{MC});p_{T} (GeV/c);Weight Var1;",nbins,xbins);
    fDsWeightVar1 = new TH1F("fDsWeightVar1","(D^{s}/D^{0})_data*(D^{0}_data/D^{s}_{MC});p_{T} (GeV/c);Weight Var1;",nbins,xbins);
    
    fLcWeightVar1 = new TH1F("fLcWeightVar1","Lc weight, Lc/D0 from model;p_{T} (GeV/c);Weight Var1;",nbins,xbins);
    fLcWeightVar2 = new TH1F("fLcWeightVar2","Lc weight, Lc/D0 from data;p_{T} (GeV/c);Weight Var2;",nbins,xbins);
    //Double_t ratio[13] = {2.03552,1.0201,0.45925,0.211574,0.11987,0.0898116,0.0631282,0.0546798,0.0477205,0.0410021,0.0307936,0.0398483,0.0175335};
    //Double_t err[13] = {0.541651,0.146443,0.0498454,0.024907,0.01438,0.0107908,0.00848616,0.0061723,0.00587082,0.00566712,0.00597994,0.00811015,0.00693105};
    
    if (fMinCentr==0 && fMaxCentr==10 ) {
        Double_t ratio[13] = {0.106888,0.0650239,0.0343858,0.0172579,0.00957876,0.00640323,0.00399907,0.00269269,0.00163078,0.000942387,0.000441093,0.000353811,0.000143011};
        Double_t err[13] = {0.0284416,0.0093333,0.00373075,0.00203067,0.00114824,0.000768388,0.00053676,0.000303334,0.000199878,0.000129785,8.53822e-05,7.13313e-05,5.61316e-05};
        //Double_t ratioNew[13] = {0.197449,0.118714,0.0627949,0.0321233,0.0182153,0.0124903,0.00801369,0.00553768,0.00340667,0.00193131,0.00089526,0.000678224,0.00026223};
        Double_t ratioNew[13] = {0.382299,0.223983,0.116971,0.0573891,0.0292338,0.0192381,0.011778,0.00863768,0.00534462,0.00301279,0.00159646,0.00131757,0.000523965};
        Double_t ratioVar1[13] = {0.249988,0.132914,0.0673368,0.0340132,0.0189431,0.01274,0.00801369,0.00543372,0.00326751,0.00179833,0.000779735,0.000564287,0.000159311};
        Double_t ratioVar2[13] = {0.144911,0.104514,0.058253,0.0302335,0.0174875,0.0122405,0.00801369,0.00564164,0.00354583,0.00206429,0.00101078,0.00079216,0.000365149};
        Double_t wLcVar1[13] = {1.57532,1.46238,0.948202,0.497431,0.250927,0.151636,0.0827256,0.0487993,0.0218917,0.0076259,0.00180605,0.00055039,8.79344e-05};
        Double_t wLcVar2[13] = {3.47398,1.72259,0.839403,0.514328,0.261815,0.145702,0.0895206,0.0380316,0.0234082,0.0051148,0.00266179,0.00049151,0};
        for (int idata=1; idata<14; idata++) {
            //fDWeight->SetBinContent(idata,ratio[idata-1]);
            //fDWeight->SetBinError(idata,err[idata-1]);
            
            fDWeightNew->SetBinContent(idata,ratioNew[idata-1]);
            fDWeightVar1->SetBinContent(idata,ratioVar1[idata-1]);
            fDWeightVar2->SetBinContent(idata,ratioVar2[idata-1]);
            fLcWeightVar1->SetBinContent(idata,wLcVar1[idata-1]);
            fLcWeightVar2->SetBinContent(idata,wLcVar2[idata-1]);
            fDPlusWeightVar1->SetBinContent(idata,ratioVar1[idata-1]);
            fDsWeightVar1->SetBinContent(idata,ratioVar1[idata-1]);
        }
    }else if (fMinCentr==30 && fMaxCentr==50 ) {
        Double_t ratio[13] = {0.079428,0.0402934,0.0258836,0.0165168,0.0117076,0.00807683,0.00545914,0.00413535,0.00218055,0.00147282,0.000578039,0.000286482,0.000286482};
        Double_t err[13] = {0.0112672,0.00212923,0.000929744,0.000665396,0.000505178,0.00041285,0.000325958,0.000210603,0.000152176,0.000110271,6.60032e-05,5.6821e-05,0};
        //Double_t ratioNew[13] = {0.0611667,0.0299809,0.0186145,0.0114545,0.00789462,0.0053161,0.00355211,0.002672,0.00141544,0.000970815,0.000391017,0.000201985,0.000107469};
        Double_t ratioNew[13] = {0.171312,0.0872016,0.0553106,0.034332,0.0236477,0.0159165,0.0106375,0.00797653,0.00422176,0.00292325,0.00118981,0.000624289,0.0003415285};
        Double_t ratioVar1[13] = {0.0705071,0.031752,0.019312,0.0118197,0.00807804,0.00538886,0.00355211,0.0026429,0.00137442,0.000931529,0.000360442,0.000168603,0.0000726835};
        Double_t ratioVar2[13] = {0.0518263,0.0282098,0.0179171,0.0110893,0.00771121,0.00524334,0.00355211,0.00270111,0.00145645,0.0010101,0.000421592,0.000235367,0.0001422545};
        Double_t wLcVar1[13] = {2.02089,0.889089,0.490217,0.266978,0.158602,0.0931538,0.0536931,0.0323656,0.0131682,0.00593832,0.00123755,0.000273619,0.};
        Double_t wLcVar2[13] = {0.621365,0.207712,0.0824139,0.0377918,0.0179345,0.0140968,0.00861207,0.00516433,0.00249769,0.0017025,0.000517206,6.8216e-05,0.};
        Double_t wDPlusVar1[13] = {0.125876,0.0642961,0.0414463,0.026278,0.01833,0.0124752,0.00830149,0.00628302,0.00333208,0.00227467,0.000910711,0.000469747,0.000028783};
        Double_t wDsVar1[13] = {0.317507,0.147254,0.0890879,0.0533869,0.0365636,0.0245759,0.0163403,0.0121869,0.00646731,0.00444424,0.00180416,0.000929189,0.000054218};
        /*Double_t ratio[13] = {0.566977,0.233989,0.0909109,0.0346338,0.0155742,0.00734675,0.00362088,0.00189595,0.000670163,0.000284606,5.13181e-05,8.84883e-06,8.84883e-06};
        Double_t err[13] = {0.0804276,0.0123637,0.00326455,0.0013947,0.000671652,0.000375317,0.000216074,9.65006e-05,4.67481e-05,2.13019e-05,5.85888e-06,1.75491e-06,0};
        Double_t ratioNew[13] = {0.566977,0.233989,0.0909109,0.0346338,0.0155742,0.00734675,0.00362088,0.00189595,0.000670163,0.000284606,5.13181e-05,8.84883e-06,8.84883e-06};
        Double_t ratioVar1[13] = {0.566977,0.233989,0.0909109,0.0346338,0.0155742,0.00734675,0.00362088,0.00189595,0.000670163,0.000284606,5.13181e-05,8.84883e-06,8.84883e-06};
        Double_t ratioVar2[13] = {0.566977,0.233989,0.0909109,0.0346338,0.0155742,0.00734675,0.00362088,0.00189595,0.000670163,0.000284606,5.13181e-05,8.84883e-06,8.84883e-06};*/
        for (int idata=1; idata<14; idata++) {
            //fDWeight->SetBinContent(idata,ratio[idata-1]);
            //fDWeight->SetBinError(idata,err[idata-1]);
            fDWeightNew->SetBinContent(idata,ratioNew[idata-1]);
            fDWeightVar1->SetBinContent(idata,ratioVar1[idata-1]);
            fDWeightVar2->SetBinContent(idata,ratioVar2[idata-1]);
            fLcWeightVar1->SetBinContent(idata,wLcVar1[idata-1]);
            fLcWeightVar2->SetBinContent(idata,wLcVar2[idata-1]);
            fDPlusWeightVar1->SetBinContent(idata,wDPlusVar1[idata-1]);
            fDsWeightVar1->SetBinContent(idata,wDsVar1[idata-1]);
        }
    }else{
        Double_t ratio[13] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
        Double_t err[13];
        Double_t ratioNew[13] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
        Double_t ratioVar1[13] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
        Double_t ratioVar2[13] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
        Double_t wLcVar1[13] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
        Double_t wLcVar2[13] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
        Double_t wDPlusVar1[13] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
        Double_t wDsVar1[13] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
        
        for (int idata=1; idata<14; idata++) {
            //fDWeight->SetBinContent(idata,ratio[idata-1]);
            //fDWeight->SetBinError(idata,err[idata-1]);
            fDWeightNew->SetBinContent(idata,ratioNew[idata-1]);
            fDWeightVar1->SetBinContent(idata,ratioVar1[idata-1]);
            fDWeightVar2->SetBinContent(idata,ratioVar2[idata-1]);
            fLcWeightVar1->SetBinContent(idata,wLcVar1[idata-1]);
            fLcWeightVar2->SetBinContent(idata,wLcVar2[idata-1]);
            fDPlusWeightVar1->SetBinContent(idata,wDPlusVar1[idata-1]);
            fDsWeightVar1->SetBinContent(idata,wDsVar1[idata-1]);
        }
    }
    
    //fDWeight->Sumw2();
    //fOutputList->Add(fDWeight);
    fOutputList->Add(fDWeightNew);
    fOutputList->Add(fDWeightVar1);
    fOutputList->Add(fDWeightVar2);
    fOutputList->Add(fLcWeightVar1);
    fOutputList->Add(fLcWeightVar2);
    fOutputList->Add(fDPlusWeightVar1);
    fOutputList->Add(fDsWeightVar1);
    
    
    //B Meson pt weighting
    Int_t nbinsB = 250;
    Double_t xbinsB[251] = {0.,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7,7.2,7.4,7.6,7.8,8,8.2,8.4,8.6,8.8,9,9.2,9.4,9.6,9.8,10,10.2,10.4,10.6,10.8,11,11.2,11.4,11.6,11.8,12,12.2,12.4,12.6,12.8,13,13.2,13.4,13.6,13.8,14,14.2,14.4,14.6,14.8,15,15.2,15.4,15.6,15.8,16,16.2,16.4,16.6,16.8,17,17.2,17.4,17.6,17.8,18,18.2,18.4,18.6,18.8,19,19.2,19.4,19.6,19.8,20,20.2,20.4,20.6,20.8,21,21.2,21.4,21.6,21.8,22,22.2,22.4,22.6,22.8,23,23.2,23.4,23.6,23.8,24,24.2,24.4,24.6,24.8,25,25.2,25.4,25.6,25.8,26,26.2,26.4,26.6,26.8,27,27.2,27.4,27.6,27.8,28,28.2,28.4,28.6,28.8,29,29.2,29.4,29.6,29.8,30,30.2,30.4,30.6,30.8,31,31.2,31.4,31.6,31.8,32,32.2,32.4,32.6,32.8,33,33.2,33.4,33.6,33.8,34,34.2,34.4,34.6,34.8,35,35.2,35.4,35.6,35.8,36,36.2,36.4,36.6,36.8,37,37.2,37.4,37.6,37.8,38,38.2,38.4,38.6,38.8,39,39.2,39.4,39.6,39.8,40,40.2,40.4,40.6,40.8,41,41.2,41.4,41.6,41.8,42,42.2,42.4,42.6,42.8,43,43.2,43.4,43.6,43.8,44,44.2,44.4,44.6,44.8,45,45.2,45.4,45.6,45.8,46,46.2,46.4,46.6,46.8,47,47.2,47.4,47.6,47.8,48,48.2,48.4,48.6,48.8,49,49.2,49.4,49.6,49.8,50.};
    //fBWeight = new TH1F("fBWeight","TAMU RAA x FONLL/MC;p_{T} (GeV/c);Weight;",nbinsB,xbinsB);
    fBWeightNew = new TH1F("fBWeightNew","TAMU RAA x FONLL(New)/MC;p_{T} (GeV/c);Weight;",nbinsB,xbinsB);
    fBWeightVar1 = new TH1F("fBWeightVar1","TAMU RAA(Max) x FONLL(New)/MC;p_{T} (GeV/c);Weight;",nbinsB,xbinsB);
    fBWeightVar2 = new TH1F("fBWeightVar2","TAMU RAA(Min) x FONLL(New)/MC;p_{T} (GeV/c);Weight;",nbinsB,xbinsB);
    
    if (fMinCentr==0 && fMaxCentr==10 ) {
        Double_t ratioB[250] = {0.498558,0.62782,0.672398,0.718151,0.731677,0.734297,0.770069,0.776389,0.809054,0.825417,0.859501,0.874636,0.898427,0.923358,0.92293,0.926125,0.926343,0.916683,0.917184,0.911809,0.886661,0.868603,0.852465,0.837679,0.834915,0.813126,0.776888,0.762151,0.745226,0.726064,0.689477,0.675392,0.657326,0.639663,0.610833,0.599204,0.574208,0.551114,0.534405,0.508084,0.503657,0.468409,0.459998,0.447575,0.429903,0.413056,0.396865,0.385802,0.36568,0.358106,0.351127,0.342732,0.324592,0.325566,0.317839,0.306017,0.292589,0.292114,0.28069,0.274267,0.268707,0.260392,0.261063,0.248518,0.245039,0.237828,0.226948,0.221301,0.216985,0.214584,0.210782,0.197669,0.195025,0.190053,0.186186,0.179784,0.171765,0.169136,0.159454,0.161479,0.157118,0.153996,0.14812,0.141091,0.138405,0.137488,0.134021,0.127508,0.126735,0.119841,0.117923,0.11417,0.109935,0.110609,0.106261,0.102496,0.0988583,0.0988566,0.0962846,0.0971962,0.0912081,0.0906153,0.0909837,0.089602,0.0851854,0.0829322,0.0854181,0.0793186,0.0828916,0.079228,0.0765657,0.0755108,0.0746448,0.0749831,0.0726684,0.0741421,0.0725614,0.0725931,0.0708997,0.0691765,0.0675163,0.0638974,0.0636498,0.0696325,0.0621972,0.0654069,0.0613243,0.0604615,0.0588438,0.0602627,0.0602133,0.0562468,0.058719,0.0581855,0.0582384,0.0562759,0.053968,0.053529,0.0549016,0.0524293,0.0523908,0.0529672,0.0515757,0.0493422,0.0481627,0.0482736,0.0459697,0.0458283,0.0453141,0.0433415,0.045522,0.043535,0.0411514,0.0456961,0.0441323,0.0430352,0.0428453,0.0426922,0.0436595,0.040605,0.0381453,0.03843,0.0410627,0.0374946,0.0391381,0.0379844,0.0367569,0.036893,0.0372399,0.0355532,0.0336221,0.0344674,0.0330855,0.0332478,0.0308431,0.0318583,0.0313297,0.03175,0.0314691,0.0320693,0.0312073,0.0297707,0.0285189,0.0289505,0.0291229,0.0291655,0.0305036,0.0281665,0.0276302,0.0291551,0.026976,0.027417,0.027251,0.0245767,0.0246722,0.024308,0.0247152,0.0273125,0.0242685,0.0232144,0.0220058,0.023307,0.0224514,0.0219948,0.0208518,0.0216804,0.0215576,0.0195631,0.0195621,0.0211566,0.0185032,0.0196948,0.0190276,0.0193687,0.0191364,0.0195232,0.0190942,0.0177003,0.0186524,0.0180672,0.017813,0.0173329,0.0163168,0.0163818,0.0158271,0.0169163,0.0164239,0.0157851,0.0169332,0.0180851,0.0152819,0.0149112,0.0141721,0.015732,0.0163534,0.014337,0.0143059,0.014786,0.0140978,0.0145358,0.0145749,0.0134334,0.0130966,0.0138116,0.013262,0.0143179,0.0132258,0.0135375,0.0132521,0.0132662};
        Double_t ratioBNew[250] = {0.527854,0.676468,0.697991,0.711774,0.722505,0.744184,0.771669,0.791643,0.815676,0.842682,0.869755,0.886323,0.900165,0.921397,0.927047,0.93348,0.936192,0.926124,0.917184,0.909077,0.889903,0.876369,0.858462,0.838199,0.821877,0.805409,0.785763,0.762809,0.741829,0.72182,0.704218,0.683978,0.660144,0.639746,0.617747,0.596334,0.573373,0.553575,0.533483,0.516636,0.49255,0.475004,0.460341,0.441742,0.4269,0.412142,0.40238,0.387698,0.37121,0.360204,0.352943,0.341443,0.332603,0.323606,0.315956,0.305826,0.295788,0.291485,0.28239,0.275661,0.268983,0.260012,0.25416,0.24767,0.24152,0.234847,0.229889,0.2255,0.216294,0.21002,0.204699,0.200694,0.194504,0.189242,0.185419,0.179002,0.173614,0.169286,0.162982,0.158871,0.15605,0.149979,0.146925,0.143822,0.138876,0.135658,0.131925,0.128432,0.125932,0.121191,0.11921,0.116162,0.111859,0.11006,0.106977,0.105162,0.103074,0.100297,0.0993633,0.0962979,0.0940254,0.0921904,0.0906611,0.0891579,0.0880537,0.0863541,0.0855137,0.0817595,0.0811291,0.0799165,0.0793575,0.0781609,0.0766887,0.0750305,0.0738189,0.073269,0.072447,0.0706597,0.0707376,0.0681811,0.0680168,0.0670508,0.065782,0.0663052,0.0645379,0.0638748,0.0625307,0.0614706,0.0613937,0.0600465,0.0599109,0.0578851,0.0589849,0.0577687,0.0577089,0.0564293,0.0548461,0.0541342,0.0537008,0.0534945,0.0522089,0.0527365,0.0517209,0.0505692,0.0501077,0.0491331,0.0493484,0.0485754,0.0477572,0.0471697,0.0461767,0.0455698,0.0463311,0.0454191,0.0441808,0.0445879,0.0432137,0.043078,0.0429534,0.041554,0.0404588,0.040246,0.0399505,0.0400038,0.0393302,0.0378221,0.0377077,0.0375886,0.0364395,0.0361215,0.0359106,0.0354513,0.035101,0.0343332,0.0340229,0.033936,0.0329639,0.0324529,0.0320335,0.0311993,0.0318253,0.0314773,0.0295871,0.0297416,0.0294255,0.0293431,0.0287639,0.028382,0.0275853,0.0278151,0.0270327,0.026975,0.0267576,0.0258747,0.0254334,0.0253068,0.0250752,0.0246926,0.0236031,0.0242901,0.0232739,0.0232222,0.0229169,0.0232431,0.0216698,0.0220481,0.0212689,0.0212692,0.0212703,0.0207071,0.0207659,0.0199056,0.019779,0.0194439,0.0190936,0.0189077,0.0188995,0.017954,0.0178651,0.0179176,0.0177795,0.0169843,0.0170956,0.0171893,0.0168782,0.0165971,0.017014,0.0161936,0.0161965,0.0158903,0.0157885,0.0152633,0.0151995,0.0154153,0.0155522,0.0153605,0.0147381,0.0147705,0.0144276,0.0145398,0.0142471,0.0144064,0.0138917,0.0140977,0.014126,0.0138707,0.0134346,0.0137999,0.013464,0.0130542};
        Double_t ratioBVar1[250] = {0.54281,0.695606,0.717742,0.731959,0.743084,0.765527,0.794016,0.814866,0.839998,0.86832,0.896867,0.914756,0.930023,0.953155,0.96042,0.968762,0.97354,0.965327,0.958592,0.953066,0.936273,0.925755,0.910985,0.894068,0.881733,0.869652,0.854532,0.83615,0.820235,0.805694,0.794136,0.779845,0.761552,0.747232,0.730978,0.715231,0.697306,0.682818,0.667487,0.655676,0.633963,0.619851,0.608774,0.591691,0.578793,0.565203,0.557723,0.542693,0.524328,0.512981,0.506383,0.493153,0.483236,0.472629,0.463576,0.450508,0.437227,0.432141,0.419709,0.410571,0.401326,0.388495,0.380185,0.370806,0.36184,0.352007,0.344679,0.338147,0.324347,0.314906,0.306865,0.300774,0.291389,0.283384,0.277522,0.267772,0.259559,0.252931,0.243351,0.237051,0.232676,0.22346,0.218746,0.213963,0.206445,0.201502,0.195801,0.190463,0.186603,0.179431,0.176352,0.171701,0.165204,0.162412,0.15773,0.154925,0.151721,0.147509,0.146013,0.14139,0.137937,0.135132,0.132778,0.130467,0.128743,0.126153,0.12482,0.119241,0.118222,0.116358,0.115448,0.113612,0.11138,0.108881,0.107034,0.106149,0.104871,0.1022,0.102228,0.0984524,0.0981344,0.0966614,0.0947546,0.0954302,0.0928108,0.0917825,0.0897781,0.0881844,0.0880027,0.0860019,0.0857383,0.0827723,0.0842771,0.0824731,0.0823216,0.0804317,0.0781125,0.0770371,0.0763594,0.0760055,0.07412,0.0748097,0.0733108,0.0716217,0.070912,0.0694779,0.0697275,0.0685814,0.0673733,0.0664924,0.0650417,0.0641368,0.0651575,0.0638252,0.0620368,0.0625599,0.060585,0.060348,0.0601271,0.0581234,0.0565479,0.0562073,0.0557518,0.0557834,0.0548022,0.0526606,0.0524614,0.0522559,0.05062,0.0501403,0.0498098,0.0491355,0.0486133,0.0475141,0.0470493,0.0468939,0.0455165,0.0447775,0.0441656,0.0429835,0.0438132,0.043302,0.0406715,0.0408535,0.0403895,0.0402465,0.0394231,0.038871,0.037752,0.0380386,0.0369416,0.0368358,0.0365122,0.0352817,0.0346547,0.034457,0.034117,0.0335721,0.0320675,0.0329771,0.0315746,0.0314818,0.0310455,0.0314648,0.029314,0.0298044,0.0287306,0.0287104,0.0286915,0.0279118,0.0279712,0.0267935,0.0266041,0.026135,0.025646,0.0253783,0.0253495,0.0240644,0.0239284,0.023982,0.0237805,0.0227009,0.0228338,0.022943,0.0225121,0.0221218,0.0226618,0.0215541,0.0215431,0.0211211,0.0209714,0.0202598,0.0201613,0.0204336,0.0206009,0.0203331,0.0194958,0.0195254,0.0190591,0.0191943,0.0187951,0.0189924,0.0183014,0.0185603,0.018585,0.0182369,0.0176516,0.0181193,0.0176665,0.0171173};
        Double_t ratioBVar2[250] = {0.512898,0.65733,0.67824,0.69159,0.701926,0.722841,0.749322,0.768421,0.791354,0.817044,0.842643,0.85789,0.870306,0.889639,0.893675,0.898198,0.898843,0.886921,0.875776,0.865088,0.843533,0.826984,0.80594,0.782329,0.762022,0.741167,0.716994,0.689468,0.663422,0.637946,0.6143,0.588111,0.558735,0.53226,0.504516,0.477438,0.449439,0.424332,0.399479,0.377596,0.351136,0.330158,0.311908,0.291794,0.275008,0.259082,0.247036,0.232703,0.218092,0.207428,0.199503,0.189734,0.181971,0.174584,0.168336,0.161143,0.154349,0.150828,0.145071,0.14075,0.13664,0.131529,0.128135,0.124534,0.121201,0.117687,0.1151,0.112853,0.108241,0.105134,0.102532,0.100614,0.0976182,0.0951008,0.0933164,0.0902324,0.0876688,0.0856413,0.0826121,0.080691,0.0794234,0.0764967,0.0751029,0.0736805,0.0713075,0.0698141,0.0680495,0.0664012,0.0652605,0.0629505,0.0620669,0.0606222,0.0585142,0.0577086,0.0562239,0.0554,0.0544271,0.0530847,0.0527133,0.051206,0.0501136,0.0492492,0.0485438,0.0478486,0.0473641,0.0465557,0.0462072,0.0442785,0.0440359,0.0434749,0.0432671,0.0427093,0.0419975,0.0411799,0.0406037,0.0403892,0.0400228,0.0391198,0.0392472,0.0379098,0.0378991,0.0374402,0.0368093,0.0371801,0.0362649,0.0359672,0.0352834,0.0347569,0.0347848,0.0340911,0.0340835,0.0329978,0.0336927,0.0330643,0.0330962,0.0324268,0.0315796,0.0312313,0.0310422,0.0309835,0.0302978,0.0306634,0.030131,0.0295167,0.0293034,0.0287882,0.0289693,0.0285695,0.0281411,0.027847,0.0273117,0.0270028,0.0275047,0.027013,0.0263247,0.0266158,0.0258425,0.025808,0.0257798,0.0249847,0.0243697,0.0242847,0.0241492,0.0242241,0.0238581,0.0229835,0.022954,0.0229212,0.022259,0.0221028,0.0220115,0.021767,0.0215887,0.0211522,0.0209965,0.0209781,0.0204113,0.0201284,0.0199013,0.0194151,0.0198373,0.0196527,0.0185027,0.0186296,0.0184616,0.0184396,0.0181047,0.017893,0.0174185,0.0175915,0.0171238,0.0171143,0.017003,0.0164678,0.0162122,0.0161565,0.0160335,0.0158132,0.0151387,0.0156032,0.0149732,0.0149626,0.0147882,0.0150214,0.0140256,0.0142919,0.0138073,0.013828,0.0138492,0.0135023,0.0135606,0.0130178,0.0129538,0.0127529,0.0125413,0.012437,0.0124495,0.0118436,0.0118018,0.0118533,0.0117786,0.0112676,0.0113574,0.0114356,0.0112443,0.0110724,0.0113663,0.0108331,0.01085,0.0106594,0.0106056,0.0102667,0.0102377,0.0103971,0.0105035,0.010388,0.00998035,0.0100156,0.00979607,0.0098853,0.00969908,0.0098204,0.00948191,0.00963511,0.00966696,0.00950457,0.0092176,0.0094804,0.00926154,0.00899106};
        for (int idata=1; idata<251; idata++) {
            //fBWeight->SetBinContent(idata,ratioB[idata-1]);
            fBWeightNew->SetBinContent(idata,ratioBNew[idata-1]);
            fBWeightVar1->SetBinContent(idata,ratioBVar1[idata-1]);
            fBWeightVar2->SetBinContent(idata,ratioBVar2[idata-1]);
        }
    }else if (fMinCentr==30 && fMaxCentr==50 ) {
        Double_t ratioB[250] = {0.613763,0.784834,0.810916,0.836203,0.830989,0.868295,0.883505,0.913221,0.934821,0.961457,0.985556,1.01031,1.03038,1.04399,1.04894,1.05124,1.04879,1.0374,1.02127,1.00945,0.982878,0.956372,0.934072,0.906477,0.887155,0.864324,0.842351,0.811834,0.784108,0.757401,0.730328,0.709669,0.681074,0.661868,0.631406,0.60966,0.590193,0.567588,0.54466,0.531031,0.516403,0.489309,0.472915,0.463391,0.447918,0.434322,0.422943,0.406763,0.397189,0.386448,0.376006,0.367809,0.360888,0.348902,0.337959,0.332709,0.326533,0.309368,0.305947,0.29407,0.29256,0.281315,0.277646,0.268647,0.25997,0.252176,0.245948,0.237849,0.230473,0.226727,0.219848,0.215782,0.207203,0.201134,0.196046,0.190369,0.188318,0.181247,0.17419,0.168473,0.162579,0.157984,0.15497,0.150132,0.144187,0.143221,0.139225,0.134864,0.129645,0.129193,0.123249,0.12106,0.119186,0.115747,0.113364,0.108896,0.107027,0.103195,0.101673,0.0990433,0.09838,0.0941805,0.0918415,0.0901931,0.0891644,0.0870442,0.0855439,0.0830778,0.0831012,0.0818048,0.0817035,0.0784107,0.0771535,0.0767316,0.0750625,0.0731386,0.0727318,0.0727746,0.0723,0.0707661,0.0694893,0.067761,0.0659144,0.065337,0.063405,0.0611122,0.062151,0.0613304,0.0606145,0.0604442,0.0577981,0.057869,0.0572931,0.0565873,0.0549912,0.0541402,0.0534318,0.0533835,0.0517918,0.0516942,0.051115,0.0508591,0.0495408,0.0494513,0.0499725,0.0483911,0.0479715,0.0464031,0.0463092,0.0458376,0.0445392,0.0433594,0.0431873,0.0424705,0.0415192,0.04182,0.0414201,0.0406436,0.041191,0.0404933,0.0395873,0.0392351,0.0372448,0.0381429,0.035662,0.0358194,0.0363743,0.0352793,0.0348819,0.0344932,0.0342752,0.0337997,0.0320732,0.0319788,0.031834,0.0309246,0.0303152,0.0300479,0.0296263,0.0287396,0.0286671,0.0280178,0.0277545,0.0288423,0.027348,0.0273515,0.0254216,0.0258547,0.0251372,0.0254501,0.0239668,0.0245532,0.0241397,0.0232339,0.0232454,0.0222495,0.0223577,0.021668,0.021833,0.0208466,0.0208539,0.0211989,0.0208092,0.0197581,0.0198031,0.0191019,0.0196992,0.0190488,0.019162,0.0189245,0.0187646,0.0176105,0.0178578,0.0170691,0.0170069,0.0161859,0.0169536,0.0156688,0.0155009,0.0164949,0.0150437,0.0150612,0.0148613,0.0148058,0.0142766,0.0148055,0.0145336,0.0150874,0.0144156,0.0138823,0.0139659,0.0130704,0.0133007,0.0138011,0.0131655,0.0131122,0.0133701,0.0131427,0.0125634,0.012817,0.0123759,0.0123669,0.0120754,0.0116675,0.0115908,0.0114222,0.0114497,0.0113251,0.012022,0.011139};
        Double_t ratioBNew[250] = {0.629879,0.788083,0.813433,0.837649,0.83473,0.871466,0.88658,0.917419,0.935438,0.961492,0.987841,1.01214,1.03055,1.04675,1.05152,1.05398,1.05146,1.03976,1.02127,1.01217,0.985053,0.956285,0.936065,0.908407,0.889645,0.863525,0.844254,0.814193,0.788496,0.759985,0.73241,0.710624,0.68301,0.663172,0.632443,0.610222,0.59217,0.569207,0.545748,0.531126,0.516459,0.489473,0.472985,0.465141,0.448287,0.435402,0.42446,0.408086,0.397874,0.386006,0.376618,0.368262,0.36165,0.349324,0.339087,0.332359,0.328113,0.309907,0.305498,0.294595,0.294602,0.282511,0.27876,0.270066,0.260252,0.253447,0.247503,0.238828,0.230165,0.226219,0.219926,0.217022,0.208488,0.201821,0.196422,0.190846,0.1887,0.181658,0.174466,0.169161,0.162862,0.158847,0.154782,0.150732,0.1443,0.143006,0.139508,0.135351,0.129897,0.129112,0.123742,0.121913,0.11981,0.116241,0.113353,0.108844,0.107325,0.10318,0.102095,0.099406,0.0986112,0.0947324,0.0920699,0.0910207,0.0894683,0.087799,0.0859247,0.0831379,0.083139,0.0824289,0.0816908,0.07866,0.0770417,0.0769627,0.0754568,0.0732459,0.0730796,0.0729844,0.0725179,0.0702804,0.0694963,0.0676671,0.0658698,0.0653475,0.0636846,0.0609512,0.0621307,0.061395,0.0607265,0.0608395,0.0579561,0.0583977,0.0574282,0.0567415,0.0553473,0.0545031,0.0534561,0.0533514,0.0518442,0.0519239,0.0513135,0.0509472,0.0498151,0.049784,0.0502312,0.0484379,0.0479072,0.0464586,0.046161,0.0456644,0.0446926,0.0436307,0.0437529,0.0425312,0.0417728,0.0419363,0.0415331,0.0405329,0.0409199,0.0405974,0.0398935,0.0395114,0.0375411,0.0380874,0.0360252,0.0357682,0.0366069,0.0356826,0.0350621,0.0344909,0.0342995,0.0336487,0.0321502,0.0318939,0.0315562,0.0312153,0.0307212,0.0299212,0.0297284,0.0287068,0.0285927,0.0281519,0.0276306,0.0287455,0.027473,0.0273274,0.0257916,0.0256825,0.0252547,0.0255688,0.0242197,0.0247519,0.0242898,0.0233756,0.0233464,0.0223576,0.0222939,0.0216614,0.0219721,0.0209188,0.0210538,0.0211145,0.0208338,0.0198707,0.0201319,0.0192302,0.0196697,0.0191208,0.0193517,0.0191076,0.0189315,0.0179213,0.0181558,0.0172182,0.0169685,0.0161772,0.0170894,0.0156943,0.0156453,0.0167149,0.015151,0.0151612,0.0150915,0.0149835,0.0143214,0.0148453,0.0144279,0.01508,0.0145104,0.0138127,0.0140673,0.0131155,0.0132893,0.0139538,0.0131358,0.0130515,0.0133015,0.0131268,0.0125072,0.0127295,0.0124395,0.0124556,0.0118749,0.0117609,0.0116434,0.0114951,0.0116084,0.011298,0.0121,0.0111589};
        Double_t ratioBVar1[250] = {0.573077,0.717766,0.741751,0.764901,0.763463,0.798537,0.814112,0.84448,0.863463,0.890334,0.918046,0.944505,0.966183,0.98657,0.997,1.00608,1.01129,1.00854,1,1.00155,0.986125,0.969671,0.962596,0.948575,0.944539,0.933351,0.930149,0.915439,0.905731,0.892737,0.880535,0.874948,0.861608,0.857325,0.837869,0.828285,0.82316,0.809796,0.793977,0.789405,0.783337,0.756718,0.744393,0.744272,0.728337,0.717369,0.708315,0.688908,0.67871,0.664664,0.653959,0.644245,0.636892,0.618813,0.603806,0.594538,0.589303,0.558567,0.552316,0.534034,0.535292,0.514361,0.508419,0.493305,0.475994,0.464062,0.453608,0.438059,0.422455,0.415447,0.404079,0.398897,0.383331,0.371165,0.361304,0.351096,0.347183,0.334247,0.321023,0.311261,0.299663,0.292261,0.284764,0.277289,0.265432,0.263023,0.256558,0.248881,0.238819,0.237341,0.227436,0.224038,0.220139,0.213545,0.208205,0.199887,0.197062,0.189418,0.187392,0.182424,0.180933,0.173784,0.168868,0.166912,0.164035,0.160944,0.157478,0.152342,0.152315,0.150985,0.149604,0.144026,0.141036,0.140865,0.138082,0.13401,0.13368,0.13348,0.132602,0.128485,0.127027,0.12366,0.120352,0.119375,0.116315,0.111301,0.113433,0.112068,0.110827,0.111011,0.10573,0.106515,0.104726,0.103454,0.100893,0.0993345,0.0974076,0.097198,0.094434,0.0945608,0.0934312,0.0927463,0.0906679,0.0905938,0.0913901,0.0881104,0.0871283,0.0844774,0.0839202,0.0830013,0.0812194,0.0792744,0.079481,0.0772469,0.0758548,0.0761371,0.0753906,0.0735609,0.0742491,0.0736497,0.0723589,0.071652,0.0680658,0.069043,0.0652924,0.0648141,0.0663211,0.0646342,0.063498,0.0624515,0.0620931,0.0609033,0.0581799,0.057705,0.0570831,0.0564556,0.0555514,0.0540944,0.0537356,0.051879,0.051663,0.0508568,0.0499054,0.0519092,0.0496018,0.0493295,0.0465484,0.0463426,0.045562,0.0461198,0.043678,0.0446292,0.0437878,0.0421317,0.0420709,0.0402814,0.0401591,0.0390122,0.0395643,0.0376604,0.0378963,0.0379983,0.0374861,0.0357463,0.0362092,0.0345808,0.0353644,0.034371,0.0347796,0.0343342,0.0340113,0.0321904,0.0326054,0.0309157,0.0304616,0.0290355,0.0306669,0.0281581,0.0280648,0.0299779,0.0271678,0.0271809,0.0270508,0.0268521,0.0256608,0.0265943,0.0258418,0.0270047,0.0259797,0.0247258,0.0251768,0.0234689,0.0237754,0.0249596,0.0234919,0.0233367,0.0237791,0.0234624,0.0223507,0.0227437,0.0222214,0.0222459,0.0212047,0.0209973,0.0207835,0.0205149,0.0207132,0.0201556,0.0215822,0.0198998};
        Double_t ratioBVar2[250] = {0.360924,0.463659,0.490858,0.517781,0.527756,0.562621,0.583404,0.614127,0.635702,0.661931,0.687458,0.71048,0.728122,0.742842,0.748039,0.750202,0.747529,0.737199,0.721136,0.711003,0.687762,0.663216,0.644631,0.621141,0.604108,0.582581,0.566292,0.543469,0.524326,0.504077,0.485197,0.470848,0.453267,0.441402,0.42274,0.41011,0.400569,0.387892,0.374939,0.368077,0.361179,0.345513,0.337034,0.33457,0.325438,0.318942,0.313639,0.304058,0.298802,0.292062,0.286966,0.282445,0.279074,0.271095,0.264535,0.260545,0.258367,0.245036,0.242464,0.234625,0.235382,0.226387,0.223988,0.217547,0.210128,0.205075,0.200668,0.193996,0.187287,0.18438,0.179529,0.17742,0.170682,0.165445,0.161225,0.156841,0.155261,0.149639,0.143875,0.139651,0.134594,0.131412,0.128179,0.12495,0.119736,0.118777,0.115983,0.112633,0.108196,0.107642,0.103261,0.101827,0.100162,0.097266,0.094935,0.0912398,0.0900464,0.0866456,0.0858103,0.0836242,0.0830286,0.0798326,0.0776566,0.0768387,0.0755938,0.0742479,0.0727258,0.0704279,0.0704896,0.0699478,0.069381,0.0668643,0.0655448,0.0655336,0.0643063,0.0624754,0.0623867,0.0623584,0.0620125,0.0601501,0.0595294,0.0580116,0.0565186,0.0561177,0.0547358,0.0524305,0.0534901,0.0529011,0.052369,0.0525104,0.0500636,0.0504872,0.0496905,0.0491372,0.0479698,0.0472773,0.0464077,0.0463552,0.045083,0.0451896,0.0446953,0.0444128,0.0434617,0.0434703,0.043897,0.0423645,0.0419348,0.0407001,0.0404724,0.0400697,0.0392491,0.0383478,0.0384864,0.0374422,0.0368044,0.0369785,0.0366527,0.0357989,0.03617,0.0359139,0.0353197,0.0350095,0.0332905,0.0338021,0.0319976,0.0317948,0.0325664,0.0317695,0.031242,0.0307575,0.0306113,0.0300544,0.0287388,0.0285323,0.0282526,0.0279695,0.0275486,0.0268524,0.0267005,0.0258032,0.025721,0.0253444,0.0248946,0.0259194,0.0247914,0.0246794,0.0233107,0.0232302,0.0228611,0.0231634,0.0219583,0.0224583,0.0220562,0.0212425,0.0212324,0.0203489,0.0203066,0.0197457,0.0200444,0.0190981,0.0192363,0.0193066,0.0190645,0.0181971,0.0184505,0.0176375,0.0180544,0.017564,0.0177897,0.0175787,0.0174299,0.0165124,0.0167412,0.0158887,0.0156701,0.0149506,0.0158056,0.0145263,0.0144918,0.0154943,0.0140551,0.0140751,0.0140209,0.013931,0.0133255,0.0138232,0.0134447,0.0140628,0.0135417,0.0129002,0.0131477,0.0122673,0.0124391,0.0130708,0.0123136,0.0122436,0.0124873,0.0123325,0.011759,0.0119768,0.0117126,0.0117363,0.0111974,0.0110981,0.0109952,0.0108631,0.0109782,0.0106925,0.0114598,0.0105762};
        /*Double_t ratioB[250] = {0.528468,0.676528,0.69985,0.722598,0.719079,0.75247,0.766871,0.794031,0.81433,0.839234,0.862171,0.885957,0.905921,0.920518,0.92776,0.932943,0.93418,0.927708,0.917184,0.910725,0.891085,0.871534,0.855825,0.835218,0.822144,0.805686,0.789807,0.765574,0.74352,0.721928,0.69942,0.682462,0.657241,0.640435,0.612102,0.591613,0.572795,0.55046,0.527428,0.513092,0.497552,0.469891,0.452484,0.441642,0.425178,0.410609,0.398269,0.381579,0.371263,0.360024,0.349237,0.340697,0.333486,0.32174,0.311097,0.305812,0.299775,0.283748,0.280412,0.269394,0.267931,0.257601,0.254249,0.246051,0.238174,0.231127,0.225531,0.218232,0.211603,0.208313,0.202149,0.198573,0.190843,0.18542,0.180898,0.175828,0.174104,0.167734,0.161368,0.156232,0.150923,0.146812,0.144165,0.139813,0.134422,0.133665,0.130076,0.126139,0.12139,0.121099,0.115653,0.113723,0.112086,0.108971,0.106844,0.102746,0.101093,0.0975802,0.0962457,0.0938589,0.0933322,0.0894457,0.0873193,0.0858455,0.0849587,0.0830286,0.081686,0.0794172,0.0795254,0.0783695,0.078357,0.0752801,0.0741529,0.0738267,0.0722983,0.0705208,0.0702036,0.0703201,0.069936,0.0685253,0.0673606,0.065755,0.064031,0.0635373,0.0617238,0.0595547,0.060631,0.0598935,0.0592567,0.0591523,0.0566221,0.056751,0.0562451,0.0556102,0.054098,0.0533164,0.0526735,0.0526806,0.0511629,0.0511194,0.0505989,0.0503977,0.0491419,0.0491037,0.0496723,0.0481498,0.0477813,0.0462665,0.0462201,0.0457961,0.0445443,0.0434086,0.0432802,0.0426051,0.041693,0.0420377,0.0416778,0.0409378,0.041531,0.0408687,0.0399945,0.0396785,0.0377035,0.0386513,0.0361735,0.0363695,0.0369697,0.0358926,0.0355235,0.0351626,0.0349751,0.034524,0.0327929,0.0327288,0.0326127,0.0317123,0.031118,0.0308738,0.0304706,0.0295875,0.0295417,0.0289009,0.0286572,0.0298095,0.0282926,0.0283237,0.0263508,0.0268257,0.0261065,0.0264569,0.0249391,0.0255739,0.0251674,0.0242464,0.0242816,0.0232636,0.0233992,0.0226991,0.0228937,0.0218802,0.0219088,0.0222925,0.0219034,0.0208167,0.020884,0.0201636,0.0208137,0.0201455,0.0202843,0.0200517,0.019901,0.0186946,0.0189748,0.0181538,0.0181046,0.0172467,0.0180816,0.0167268,0.016563,0.0176415,0.0161043,0.016138,0.0159385,0.0158937,0.0153398,0.0159227,0.0156447,0.0162558,0.0155462,0.0149848,0.0150888,0.0141342,0.0143964,0.0149517,0.0142761,0.0142313,0.0145244,0.0142902,0.0136728,0.0139614,0.0134931,0.0134955,0.0131893,0.0127552,0.0126828,0.0125095,0.0125509,0.0124254,0.0132018,0.0122431};
        Double_t ratioBNew[250] = {0.528468,0.676528,0.69985,0.722598,0.719079,0.75247,0.766871,0.794031,0.81433,0.839234,0.862171,0.885957,0.905921,0.920518,0.92776,0.932943,0.93418,0.927708,0.917184,0.910725,0.891085,0.871534,0.855825,0.835218,0.822144,0.805686,0.789807,0.765574,0.74352,0.721928,0.69942,0.682462,0.657241,0.640435,0.612102,0.591613,0.572795,0.55046,0.527428,0.513092,0.497552,0.469891,0.452484,0.441642,0.425178,0.410609,0.398269,0.381579,0.371263,0.360024,0.349237,0.340697,0.333486,0.32174,0.311097,0.305812,0.299775,0.283748,0.280412,0.269394,0.267931,0.257601,0.254249,0.246051,0.238174,0.231127,0.225531,0.218232,0.211603,0.208313,0.202149,0.198573,0.190843,0.18542,0.180898,0.175828,0.174104,0.167734,0.161368,0.156232,0.150923,0.146812,0.144165,0.139813,0.134422,0.133665,0.130076,0.126139,0.12139,0.121099,0.115653,0.113723,0.112086,0.108971,0.106844,0.102746,0.101093,0.0975802,0.0962457,0.0938589,0.0933322,0.0894457,0.0873193,0.0858455,0.0849587,0.0830286,0.081686,0.0794172,0.0795254,0.0783695,0.078357,0.0752801,0.0741529,0.0738267,0.0722983,0.0705208,0.0702036,0.0703201,0.069936,0.0685253,0.0673606,0.065755,0.064031,0.0635373,0.0617238,0.0595547,0.060631,0.0598935,0.0592567,0.0591523,0.0566221,0.056751,0.0562451,0.0556102,0.054098,0.0533164,0.0526735,0.0526806,0.0511629,0.0511194,0.0505989,0.0503977,0.0491419,0.0491037,0.0496723,0.0481498,0.0477813,0.0462665,0.0462201,0.0457961,0.0445443,0.0434086,0.0432802,0.0426051,0.041693,0.0420377,0.0416778,0.0409378,0.041531,0.0408687,0.0399945,0.0396785,0.0377035,0.0386513,0.0361735,0.0363695,0.0369697,0.0358926,0.0355235,0.0351626,0.0349751,0.034524,0.0327929,0.0327288,0.0326127,0.0317123,0.031118,0.0308738,0.0304706,0.0295875,0.0295417,0.0289009,0.0286572,0.0298095,0.0282926,0.0283237,0.0263508,0.0268257,0.0261065,0.0264569,0.0249391,0.0255739,0.0251674,0.0242464,0.0242816,0.0232636,0.0233992,0.0226991,0.0228937,0.0218802,0.0219088,0.0222925,0.0219034,0.0208167,0.020884,0.0201636,0.0208137,0.0201455,0.0202843,0.0200517,0.019901,0.0186946,0.0189748,0.0181538,0.0181046,0.0172467,0.0180816,0.0167268,0.016563,0.0176415,0.0161043,0.016138,0.0159385,0.0158937,0.0153398,0.0159227,0.0156447,0.0162558,0.0155462,0.0149848,0.0150888,0.0141342,0.0143964,0.0149517,0.0142761,0.0142313,0.0145244,0.0142902,0.0136728,0.0139614,0.0134931,0.0134955,0.0131893,0.0127552,0.0126828,0.0125095,0.0125509,0.0124254,0.0132018,0.0122431};
        Double_t ratioBVar1[250] = {0.543441,0.695667,0.719653,0.743089,0.73956,0.774051,0.78908,0.817324,0.838612,0.864768,0.889046,0.914378,0.93597,0.952245,0.961158,0.968204,0.971448,0.966978,0.958592,0.954794,0.937517,0.920647,0.908186,0.890889,0.882019,0.86995,0.858929,0.839181,0.822106,0.805814,0.788725,0.778117,0.758203,0.748037,0.724299,0.709568,0.696604,0.678976,0.659911,0.651178,0.640401,0.613179,0.598384,0.591556,0.576458,0.563099,0.552026,0.534127,0.524402,0.512725,0.501066,0.492075,0.484518,0.469903,0.456447,0.450488,0.443121,0.420671,0.416769,0.401237,0.399756,0.384892,0.380318,0.368382,0.356827,0.346432,0.338145,0.327248,0.317312,0.312347,0.303043,0.297595,0.285905,0.27766,0.270755,0.263024,0.260292,0.250612,0.240941,0.233114,0.225032,0.218743,0.214638,0.208,0.199823,0.198542,0.193057,0.187062,0.179873,0.179295,0.17109,0.168097,0.165539,0.160804,0.157534,0.151365,0.148804,0.143514,0.141432,0.137809,0.13692,0.131108,0.127884,0.12562,0.124218,0.121294,0.119233,0.115824,0.115885,0.114106,0.113992,0.109425,0.107697,0.107134,0.104829,0.102167,0.101624,0.101708,0.10107,0.0989493,0.0971877,0.0947934,0.0922324,0.0914466,0.088764,0.0855749,0.0870506,0.0859218,0.0849395,0.0847212,0.0810318,0.0811507,0.0803625,0.0793915,0.0771707,0.0759947,0.0750183,0.0749686,0.0727507,0.072631,0.0718343,0.0714919,0.0696553,0.069546,0.0702958,0.0680875,0.0675132,0.0653215,0.0652048,0.0645562,0.0627425,0.061095,0.0608669,0.0598709,0.0585437,0.0589818,0.0584316,0.0573498,0.0581359,0.0571648,0.055899,0.0554147,0.0526161,0.0538975,0.0504038,0.0506381,0.0514347,0.0498981,0.0493476,0.0488092,0.0485122,0.0478504,0.0454168,0.0452938,0.0450992,0.0438211,0.0429676,0.0425987,0.0420108,0.0407629,0.0406695,0.0397577,0.0393932,0.0409468,0.0388344,0.0388484,0.0361157,0.0367395,0.0357283,0.0361813,0.0340805,0.0349224,0.0343423,0.0330614,0.0330853,0.0316751,0.0318366,0.0308616,0.0311037,0.0297053,0.0297226,0.0302213,0.0296726,0.0281802,0.0282509,0.0272569,0.0281156,0.0271936,0.0273614,0.0270285,0.0268063,0.0251633,0.0255225,0.0244009,0.0243176,0.0231489,0.0242524,0.0224195,0.0221844,0.0236123,0.0215398,0.0215699,0.0212883,0.0212137,0.0204601,0.0212229,0.0208379,0.0216368,0.0206781,0.0199176,0.020042,0.0187612,0.0190961,0.019819,0.0189105,0.0188383,0.0192131,0.0188905,0.018062,0.0184308,0.0178004,0.0177915,0.0173761,0.0167928,0.0166862,0.0164472,0.0164905,0.0163147,0.0173225,0.0160538};
        Double_t ratioBVar2[250] = {0.513494,0.657388,0.680046,0.702107,0.698598,0.73089,0.744663,0.770738,0.790048,0.813701,0.835295,0.857536,0.875872,0.888791,0.894362,0.897681,0.896912,0.888437,0.875776,0.866656,0.844654,0.822421,0.803464,0.779547,0.762269,0.741421,0.720684,0.691967,0.664935,0.638041,0.610115,0.586807,0.556278,0.532833,0.499906,0.473658,0.448987,0.421944,0.394945,0.375005,0.354702,0.326604,0.306584,0.291727,0.273898,0.258118,0.244512,0.22903,0.218123,0.207324,0.197408,0.189319,0.182453,0.173577,0.165747,0.161136,0.156429,0.146825,0.144055,0.13755,0.136105,0.130309,0.12818,0.12372,0.119522,0.115823,0.112918,0.109216,0.105894,0.104279,0.101255,0.099551,0.095781,0.0931801,0.0910409,0.0886326,0.0879164,0.0848564,0.0817939,0.0793508,0.0768141,0.0748818,0.0736922,0.071627,0.0690201,0.0687883,0.0670956,0.0652155,0.0629067,0.0629027,0.060215,0.0593496,0.0586328,0.0571375,0.0561542,0.054127,0.0533809,0.0516468,0.0510594,0.0499091,0.0497441,0.0477829,0.0467545,0.0460709,0.0456993,0.0447628,0.0441389,0.0430099,0.0431655,0.0426333,0.0427216,0.0411352,0.0406088,0.0405192,0.0397673,0.0388742,0.0387834,0.0389317,0.0388025,0.0381012,0.0375335,0.0367166,0.0358295,0.0356281,0.0346837,0.0335346,0.0342115,0.0338651,0.033574,0.0335835,0.0322125,0.0323513,0.0321277,0.0318289,0.0310254,0.030638,0.0303287,0.0303927,0.0295751,0.0296079,0.0293635,0.0293035,0.0286286,0.0286613,0.0290488,0.0282121,0.0280493,0.0272114,0.0272353,0.0270361,0.0263462,0.0257221,0.0256935,0.0253394,0.0248424,0.0250935,0.024924,0.0245258,0.0249261,0.0245726,0.0240901,0.0239423,0.0227909,0.0234051,0.0219433,0.0221008,0.0225047,0.021887,0.0216995,0.021516,0.021438,0.0211977,0.0201691,0.0201638,0.0201262,0.0196034,0.0192683,0.019149,0.0189303,0.0184121,0.018414,0.0180441,0.0179212,0.0186722,0.0177508,0.017799,0.0165858,0.0169118,0.0164847,0.0167326,0.0157976,0.0162253,0.0159925,0.0154314,0.015478,0.0148521,0.0149618,0.0145365,0.0146837,0.0140551,0.0140949,0.0143636,0.0141343,0.0134533,0.013517,0.0130702,0.0135118,0.0130974,0.0132072,0.013075,0.0129958,0.0122258,0.0124272,0.0119067,0.0118916,0.0113445,0.0119107,0.0110341,0.0109416,0.0116706,0.0106688,0.0107062,0.0105887,0.0105737,0.0102194,0.0106225,0.0104515,0.0108747,0.0104143,0.010052,0.0101356,0.00950728,0.00969676,0.0100843,0.00964163,0.00962429,0.00983562,0.00968996,0.00928355,0.00949209,0.00918577,0.00919945,0.0090025,0.00871755,0.00867932,0.00857183,0.00861128,0.00853617,0.00908116,0.00843244};*/
        for (int idata=1; idata<251; idata++) {
            //fBWeight->SetBinContent(idata,ratioB[idata-1]);
            fBWeightNew->SetBinContent(idata,ratioBNew[idata-1]);
            fBWeightVar1->SetBinContent(idata,ratioBVar1[idata-1]);
            fBWeightVar2->SetBinContent(idata,ratioBVar2[idata-1]);
        }
    }else{
        Double_t ratioB[250];
        Double_t ratioBNew[250];
        Double_t ratioBVar1[250];
        Double_t ratioBVar2[250];
        for (int idata=1; idata<251; idata++) {
            //fBWeight->SetBinContent(idata,ratioB[idata-1]);
            fBWeightNew->SetBinContent(idata,ratioBNew[idata-1]);
            fBWeightVar1->SetBinContent(idata,ratioBVar1[idata-1]);
            fBWeightVar2->SetBinContent(idata,ratioBVar2[idata-1]);
        }
    }

    //fOutputList->Add(fBWeight);
    fOutputList->Add(fBWeightNew);
    fOutputList->Add(fBWeightVar1);
    fOutputList->Add(fBWeightVar2);

//add by sudo
Int_t Sparsebins[6]={  100, 100, 100,  50,  50, 200}; // trigger;pT;nSigma;eop;m20;m02;sqrtm02m20;eID;iSM;cent
Double_t Sparsexmin[6]={ 0,   0, -10,  -5,  -5,   0};
Double_t Sparsexmax[6]={10,  10,  10,   5,   5,   2};
fSparseElectron = new THnSparseD ("fSparseElectron","correlation;Pt;P;TPCnsigma;ITSnsigma;TOFnsigma;E/p;",6,Sparsebins,Sparsexmin,Sparsexmax);
if(iTree)fOutputList -> Add(fSparseElectron);


const int ncentbins = static_cast<int>(fMaxCentr-fMinCentr);

TString detConfName[6] = {"TPC","TPCPosEta","TPCNegEta","V0","V0A","V0C"};

TString qnaxisname=Form("#it{q}_{%d}",fHarmonic);
//cout << "fqnMeth = " << fqnMeth << endl;
switch(fqnMeth) {
	case kq2TPC:
		qnaxisname=Form("#it{q}_{%d}^{TPC}",fHarmonic);
		break;
	case kq2PosTPC:
		qnaxisname=Form("#it{q}_{%d}^{TPCPosEta}",fHarmonic);
		break;
	case kq2NegTPC:
		qnaxisname=Form("#it{q}_{%d}^{TPCNegEta}",fHarmonic);
		break;
	case kq2VZERO:
		qnaxisname=Form("#it{q}_{%d}^{V0}",fHarmonic);
		break;
	case kq2VZEROA:
		qnaxisname=Form("#it{q}_{%d}^{V0A}",fHarmonic);
		break;
	case kq2VZEROC:
		qnaxisname=Form("#it{q}_{%d}^{V0C}",fHarmonic);
		break;
}

//cout << "fPercentileqn = " << fPercentileqn << endl;

TString qnpercaxisname = qnaxisname + " (%)";
TString qnaxisnamefill = qnaxisname;
if(fPercentileqn)
	qnaxisnamefill = qnpercaxisname;

	int nqnbins=1; //single bin if unbiased analysis
	double qnmin = 0.;
	double qnmax = 15.;
	if(fFlowMethod==kEvShapeEP || fFlowMethod==kEvShapeSP || fFlowMethod==kEvShapeEPVsMass) {
		if(!fPercentileqn) {
			nqnbins=300;
		}
		else {
			nqnbins=100;
			qnmin = 0.;
			qnmax = 100.;
		}
	}



for(int iDet = 0; iDet < 6; iDet++){

	//fHistEvPlane[iDet] = new TH1F(Form("fHistEvPlane%s",detConfName[iDet].Data()),Form("hEvPlane%s",detConfName[iDet].Data()),200,0.,TMath::Pi());
	fHistEvPlane[iDet] = new TH1F(Form("fHistEvPlane%s",detConfName[iDet].Data()),Form("hEvPlane%s",detConfName[iDet].Data()),200,-TMath::Pi()/2.,TMath::Pi()/2.);
	fOutputList->Add(fHistEvPlane[iDet]);

	fHistEvPlaneQncorr[iDet] = new TH2F(Form("fHistEvPlaneQncorr%sVsCent",detConfName[iDet].Data()),Form("hEvPlaneQncorr%s;centrality(%%);#psi_{%d}",detConfName[iDet].Data(),fHarmonic),ncentbins,fMinCentr,fMaxCentr,200,0.,TMath::Pi());

	//fHistEvPlaneQncorr[iDet] = new TH3F(Form("fHistEvPlaneQncorr%sVsqnVsCent",detConfName[iDet].Data()),Form("hEvPlaneQncorr%s;centrality(%%);%s;#psi_{%d}",detConfName[iDet].Data(),qnaxisnamefill.Data(),fHarmonic),ncentbins,fMinCentr,fMaxCentr,nqnbins,qnmin,qnmax,200,0.,TMath::Pi());
	fOutputList->Add(fHistEvPlaneQncorr[iDet]);

	//histos for qn vs. centrality with fine binning (for qn percentiles calibration)
		if(fFlowMethod == kEvShapeEP || fFlowMethod == kEvShapeSP || fFlowMethod == kEvShapeEPVsMass){

			fHistqnVsCentrPercCalib[iDet] = new TH2F(Form("fHistqnVsCentr%s",detConfName[iDet].Data()),Form("#it{q}_{%d}^{%s} vs. centrality;centrality(%%);#it{q}_{%d}^{%s}",fHarmonic,detConfName[iDet].Data(),fHarmonic,detConfName[iDet].Data()),ncentbins,fMinCentr,fMaxCentr,15000,0.,15.);
		}

	}

//EP / Qn resolutions
TString detLabels[3][2] = {{"A","B"},{"A","C"},{"B","C"}};

for(int iResoHisto=0; iResoHisto<3; iResoHisto++){

	if(fFlowMethod == kEP || fFlowMethod == kEvShapeEP || fFlowMethod == kEPVsMass || fFlowMethod == kEvShapeEPVsMass){
		fHistEPResolVsCentrVsqn[iResoHisto] = new TH2F(Form("fHistEvPlaneReso%d",iResoHisto+1),Form("Event plane angle Resolution;centrality (%%);%s;cos2(#psi_{%s}-#psi_{%s})",qnaxisnamefill.Data(),detLabels[iResoHisto][0].Data(),detLabels[iResoHisto][1].Data()),ncentbins,fMinCentr,fMaxCentr,220,-1.1,1.1);
		//fHistEPResolVsCentrVsqn[iResoHisto] = new TH3F(Form("fHistEvPlaneReso%d",iResoHisto+1),Form("Event plane angle Resolution;centrality (%%);%s;cos2(#psi_{%s}-#psi_{%s})",qnaxisnamefill.Data(),detLabels[iResoHisto][0].Data(),detLabels[iResoHisto][1].Data()),ncentbins,fMinCentr,fMaxCentr,nqnbins,qnmin,qnmax,220,-1.1,1.1);
	}

	else if(fFlowMethod == kSP || fFlowMethod == kEvShapeSP){
		fHistEPResolVsCentrVsqn[iResoHisto] = new TH2F(Form("hScalProdQnVectors%d",iResoHisto+1),Form("Scalar product between Q-vectors;centrality (%%);%s;(Q{%s}Q{%s})",qnaxisnamefill.Data(),detLabels[iResoHisto][0].Data(),detLabels[iResoHisto][1].Data()),ncentbins,fMinCentr,fMaxCentr,200,-fScalProdLimit*fScalProdLimit,fScalProdLimit*fScalProdLimit);
		//fHistEPResolVsCentrVsqn[iResoHisto] = new TH3F(Form("hScalProdQnVectors%d",iResoHisto+1),Form("Scalar product between Q-vectors;centrality (%%);%s;(Q{%s}Q{%s})",qnaxisnamefill.Data(),detLabels[iResoHisto][0].Data(),detLabels[iResoHisto][1].Data()),ncentbins,fMinCentr,fMaxCentr,nqnbins,qnmin,qnmax,200,-fScalProdLimit*fScalProdLimit,fScalProdLimit*fScalProdLimit);

	}

	fOutputList->Add(fHistEPResolVsCentrVsqn[iResoHisto]);

}



PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the
// fOutputList object. the manager will in the end take care of writing your output to file
// so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalRun2::UserExec(Option_t *)
{

//============= Systematic Parameter ================
//Track cut
Double_t CutDCAxy = DCAxy;
Double_t CutDCAz = DCAz;

//PID cut
Double_t CutEopHad = -3.5;

    cout<<"fMinCentr ="<<fMinCentr<<endl;
    cout<<"fMaxCentr ="<<fMaxCentr<<endl;
cout << "cut selections ---------------------" << endl;
cout << "tpcnsig = " << ftpcnsig << endl;
cout << "emceop = " << femceop << endl;
cout << "emcss_mim = " << femcss_mim << endl;
cout << "emcss_max = " << femcss_max << endl;
cout << "invmass = " << finvmass << endl;
cout << "invmass_pt = " << finvmass_pt << endl;
cout << "iCentralt = " << iCentral << endl;
cout << "iSemiCentralt = " << iSemiCentral << endl;
cout << "-------------------------------------" << endl;

//===================================================


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
    fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));


    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fVevent) {
	    printf("ERROR: fVEvent not available\n");
	    return;
    }


  //TString firedTrigger;

////////////if Tender////////////////
//if(fUseTender){
////new branches with calibrted tracks and clusters
//if(IsAODanalysis()) fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("tracks"));
//
//fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("caloClusters"));
//}

    fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("tracks"));
    fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("caloClusters"));

////////////PID initialised///////
    fpidResponse = fInputHandler->GetPIDResponse();

    /////////Centrality////////
     Float_t lPercentile = 300;
     AliMultSelection *MultSelection = 0x0;
      MultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
      if( !MultSelection) {
      	    //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
      	    AliWarning("AliMultSelection object not found!");
      }else{
      	    lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
      }

    flPercentile -> Fill(lPercentile);

    //if(TMath::Abs(lPercentile)<30 || TMath::Abs(lPercentile)>50)return;
    //cout <<"lPercentile = "<< lPercentile << endl;
    if(TMath::Abs(lPercentile)<fMinCentr || TMath::Abs(lPercentile)>fMaxCentr)return;

    ///////Event Plane for 2015//////
//   flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections *> (AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
//
//   if (flowQnVectorTask != NULL) {
//           /* AliQnCorrectionsManager *fFlowQnVectorMgr; shall be a member of the user's analysis task */
//           /* and store the framework manager */
//           fFlowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
//   }
//   else {
//           AliFatal("Flow Qn vector corrections framework needed but it is not present. ABORTING!!!");
//   }
//
//  //select the desired harmonic
//   Int_t myHarmonic = 2;
//   const AliQnCorrectionsQnVector *myQnVector;//VOA
//   const AliQnCorrectionsQnVector *myQnVector2;//VOC
//   const AliQnCorrectionsQnVector *myQnVector3;//TPC
//
//   Double_t myEventPlane = 0.0;
//   Double_t myEventPlane2 = 0.0;
//   Double_t myEventPlane3 = 0.0;
//
//   Float_t Qx1 = 0.0;
//   Float_t Qy1 = 0.0;
//   //Float_t Qx2 = 0.0;
//   //Float_t Qy2 = 0.0;
//   //Float_t Qx3 = 0.0;
//   //Float_t Qy3 = 0.0;
//
//   //get the fully corrected Qn vector from VZEROA sub-detector
//   myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA");
//   myQnVector2 = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC");
//   myQnVector3 = fFlowQnVectorMgr->GetDetectorQnVector("TPC");
//
//   //cout << myQnVector << endl;
//
//   //always check for a valid returned Qn vector
//   if (myQnVector != NULL){
//	   myEventPlane = myQnVector->EventPlane(myHarmonic);
//	   Qx1 = myQnVector->Qx(myHarmonic);
//	   Qy1 = myQnVector->Qy(myHarmonic);
//   }
//
//   // cout << myEventPlane << endl;
//   cout << "Qx1 = " << Qx1 << endl;
//   cout << "Qy1 = " << Qy1 << endl;
//
//
//   if (myQnVector2 !=NULL)
//	   myEventPlane2 = myQnVector2->EventPlane(myHarmonic);
//   //Qx2 = myQnVector2->Qx(myHarmonic);
//   //Qy2 = myQnVector2->Qy(myHarmonic);
//
//   if (myQnVector3 !=NULL)
//	   myEventPlane3 = myQnVector3->EventPlane(myHarmonic);
//   //Qx3 = myQnVector3->Qx(myHarmonic);
//   //Qy3 = myQnVector3->Qy(myHarmonic);
//
//   //fmyEventPlane->Fill(myEventPlane);
//
//   TList *myQnFrameworkList = fFlowQnVectorMgr -> GetQnVectorList();
   //const TList *fmdcQnList = fFlowQnVectorMgr -> GetDetectorQnVectorList("FMDC");

   /* friendly printout of the Qn vectors lists */
//  fFlowQnVectorMgr->GetQnVectorList()->Print("",-1);

   /* friendly printout of the FMDC Qn vectors list */
   //fFlowQnVectorMgr->GetDetectorQnVectorList("FMDC")->Print("",-1);

   // Double_t epV0A = 0, epV0C = 0, epV0 = 0;
   // Double_t qxV0A = 0, qyV0A = 0, qxV0C = 0, qyV0C = 0, qxV0 = 0, qyV0 = 0;


   // if(AOD){

   //          epV0 = TVector2::Phi_0_2pi(fAOD->GetEventplane()->CalculateVZEROEventPlane(fAOD,10,2,qxV0,qyV0));

   //	    epV0A = TVector2::Phi_0_2pi(fAOD->GetEventplane()->CalculateVZEROEventPlane(fAOD,8,2,qxV0A,qyV0A));

   //	    epV0C = TVector2::Phi_0_2pi(fAOD->GetEventplane()->CalculateVZEROEventPlane(fAOD,9,2,qxV0C,qyV0C));

   // }

   //if(epV0 > TMath::Pi()) epV0 = epV0 - TMath::Pi();
   //if(epV0A > TMath::Pi()) epV0A = epV0A - TMath::Pi();
   //if(epV0C > TMath::Pi()) epV0C = epV0C - TMath::Pi();


   //////////////////////////////////////Event vertex////////////////////////////////////
    fNevents -> Fill(0);//all events
    Double_t Zvertex = -100, Xvertex = -100, Yvertex = -100;
    const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    Double_t NcontV = pVtx->GetNContributors();
    if(NcontV<2)return;
    fNevents -> Fill(1);//events with 2 tracks

    Zvertex = pVtx->GetZ();
    Xvertex = pVtx->GetX();
    Yvertex = pVtx->GetY();
    fVtxZ->Fill(Zvertex);
    fVtxX->Fill(Xvertex);
    fVtxY->Fill(Yvertex);

    //Event selection
    if(TMath::Abs(Zvertex)>10.0)return;
    //if(myEventPlane == 0.0 && myEventPlane2 == 0.0)return;//to reject NULL events
    fNevents->Fill(2);//event after z vtx cut
    //fCent->Fill(centrality);//centrality dist.

    //Get sign of B field
    int Bsign = 0;
    if(fAOD->GetMagneticField() < 0)Bsign = -1;
    if(fAOD->GetMagneticField() > 0)Bsign = 1;

    //cout << "Bsign = " << Bsign << endl;

    //fmyEventPlane->Fill(myEventPlane);     //EP of V0A
    //fmyEventPlane2->Fill(myEventPlane2);   //EP of V0C
    //fmyEventPlane3->Fill(myEventPlane3);   //EP of TPC

    //fQx1->Fill(lPercentile,Qx1);
    //fQy1->Fill(lPercentile,Qy1);
    //fQx2->Fill(lPercentile,Qx2);
    //fQy2->Fill(lPercentile,Qy2);
    //fQx3->Fill(lPercentile,Qx3);
    //fQy3->Fill(lPercentile,Qy3);

    //fEPcorV0AC->Fill(myEventPlane,myEventPlane2);    //correlation of V0A and V0C
    //fEPcorV0ATPC->Fill(myEventPlane,myEventPlane3);  //correlation of V0A and TPC
    //fEPcorV0CTPC->Fill(myEventPlane2,myEventPlane3); //correlation of V0C and TPC

    //Double_t subV0AC = -999.;
    //Double_t subV0ATPC = -999.;
    //Double_t subV0CTPC = -999.;

    //subV0AC = myEventPlane - myEventPlane2;
    //subV0ATPC = myEventPlane - myEventPlane3;
    //subV0CTPC = myEventPlane2 - myEventPlane3;

	//Double_t const PI = TMath::Pi();
	//Double_t const PI_2 = TMath::Pi()/2.;

        ////V0AC

	//if(subV0AC > PI_2){
	//	subV0AC -= PI;
	//}

	//if(subV0AC < -PI_2){
	//	subV0AC += PI;
	//}
	//Double_t subV0ACcos2 = -999.;

	//subV0ACcos2 = cos(2*subV0AC);
	//fsubV0ACcos2 -> Fill(lPercentile,subV0ACcos2);

	////V0ATPC

	//if(subV0ATPC > PI_2){
	//	subV0ATPC -= PI;
	//}

	//if(subV0ATPC < -PI_2){
	//	subV0ATPC += PI;
	//}
	//Double_t subV0ATPCcos2 = -999.;

	//subV0ATPCcos2 = cos(2*subV0ATPC);
	//fsubV0ATPCcos2 -> Fill(lPercentile,subV0ATPCcos2);

	////V0CTPC

	//if(subV0CTPC > PI_2){
	//	subV0CTPC -= PI;
	//}

	//if(subV0CTPC < -PI_2){
	//	subV0CTPC += PI;
	//}
	//Double_t subV0CTPCcos2 = -999.;

	//subV0CTPCcos2 = cos(2*subV0CTPC);
	//fsubV0CTPCcos2 -> Fill(lPercentile,subV0CTPCcos2);

	///////Event Plane for 2018//////
    double evCentr = lPercentile;
    //double evCentr = fRDCuts->GetCentrality(fAOD);

	//Get Qn-vectors from tender task
    AliHFQnVectorHandler *HFQnVectorHandler = nullptr;
    bool isHandlerFound = false;

    cout << "<----------- fTenderTaskName = " << fTenderTaskName << endl;
    AliAnalysisTaskSEHFTenderQnVectors *HFQnVectorTask = dynamic_cast<AliAnalysisTaskSEHFTenderQnVectors*>(AliAnalysisManager::GetAnalysisManager()->GetTask(fTenderTaskName.Data()));

    cout << "<----------- HFQnVectorTask = " << HFQnVectorTask << endl;

    if(HFQnVectorTask){

	    //cout << "HFQnVectorTask = " << HFQnVectorTask << endl;

	    HFQnVectorHandler = HFQnVectorTask->GetQnVectorHandler();

	    if(fPercentileqn && !fqnSplinesList[0]){
		    for(int iDet=0; iDet<6; iDet++){
			    fqnSplinesList[iDet] = dynamic_cast<TList*>(HFQnVectorTask->GetSplineForqnPercentileList(iDet));

		    }
	    }
    }

    if(HFQnVectorHandler){

	    isHandlerFound = true;


	    if(HFQnVectorHandler->GetHarmonic()!=fHarmonic){
		    AliWarning("Harmonic of task and Qn-vector handler not consistent!");
		    return;
	    }

	    if(HFQnVectorHandler->GetCalibrationType() != fCalibType){
		    AliWarning("Calibration strategy of task and Qn-vector handler not consistent");
		    return;

	    }

	    if(HFQnVectorHandler->GetNormalisationMethod() != fNormMethod){
		    AliWarning("Normalisation method of task and Qn-vector handler not consistent");
		    return;
	    }

	    if(fCalibType==AliHFQnVectorHandler::kQnCalib && HFQnVectorHandler->GetCalibrationsOADBFileName() != fOADBFileName){
		    AliWarning("OADB file name for calibrations of task and Qn-vector handler not consistent");
		    AliWarning(HFQnVectorHandler->GetCalibrationsOADBFileName());
		    return;
	    }


    }

    else{//create a new handler if not found in tender task

	    //cout << "isHandlerFound =" << isHandlerFound << endl;

	    AliWarning("Qn-vector tender task not found! Create a new one");
	    HFQnVectorHandler = new AliHFQnVectorHandler(fCalibType,fNormMethod,fHarmonic,fOADBFileName);
	    HFQnVectorHandler->SetAODEvent(fAOD);
	    HFQnVectorHandler->ComputeCalibratedQnVectorV0();
	    HFQnVectorHandler->ComputeCalibratedQnVectorTPC();
    }

    if(fPercentileqn && !fqnSplinesList[0]){
	    fLoadedSplines=LoadSplinesForqnPercentile();
    }

    double QnFullTPC[2], QnPosTPC[2],QnNegTPC[2];
    double QnFullV0[2],QnV0A[2],QnV0C[2];
    double MultQnFullTPC = -1.,MultQnPosTPC = -1.,MultQnNegTPC = -1.;
    double MultQnFullV0 = -1.,MultQnV0A = -1.,MultQnV0C = -1.;
    double PsinFullTPC = -1.,PsinPosTPC = -1.,PsinNegTPC = -1.;
    double PsinFullV0 = -1.,PsinV0A = -1.,PsinV0C = -1.;

    //cout << "HFQnVectorHandler" << HFQnVectorHandler << endl;

    //get the unnormalised Qn-vectors --> normalization can be done in the task
    HFQnVectorHandler->GetUnNormQnVecTPC(QnFullTPC,QnPosTPC,QnNegTPC);
    HFQnVectorHandler->GetUnNormQnVecV0(QnFullV0,QnV0A,QnV0C);

    HFQnVectorHandler->GetMultQnVecTPC(MultQnFullTPC,MultQnPosTPC,MultQnNegTPC);
    HFQnVectorHandler->GetMultQnVecV0(MultQnFullV0,MultQnV0A,MultQnV0C);

    HFQnVectorHandler->GetEventPlaneAngleTPC(PsinFullTPC,PsinPosTPC,PsinNegTPC);
    HFQnVectorHandler->GetEventPlaneAngleV0(PsinFullV0,PsinV0A,PsinV0C);

    //cout << "PsinFullTPC = " << PsinFullTPC << endl;
    //cout << "PsinPosTPC = " << PsinPosTPC << endl;
    //cout << "PsinNegTPC = " << PsinNegTPC << endl;
    //cout << "PsinFullV0 = " << PsinFullV0 << endl;
    //cout << "PsinV0A = " << PsinV0A << endl;
    //cout << "PsinV0C = " << PsinV0C << endl;

    double mainQn[2], SubAQn[2],SubBQn[2],SubCQn[2];
    double mainPsin = -1.,SubAPsin = -1.,SubBPsin = -1.,SubCPsin = -1.;
    double mainMultQn = -1.,SubAMultQn = -1.,SubBMultQn = -1.,SubCMultQn = -1.;
    GetMainQnVectorInfo(mainPsin,mainMultQn,mainQn,SubAPsin,SubAMultQn,SubAQn,SubBPsin,SubBMultQn,SubBQn,HFQnVectorHandler);

    SubCPsin = mainPsin;
    SubCMultQn = mainMultQn;
    SubCQn[0] = mainQn[0];
    SubCQn[1] = mainQn[1];

    int nsubevents = 3;
    if(fEvPlaneDet==fSubEvDetA || fEvPlaneDet==fSubEvDetB)
	    nsubevents = 2;

    double qnFullTPC = -1.,qnPosTPC = -1.,qnNegTPC = -1.;
    double qnFullV0 = -1.,qnV0A = -1.,qnV0C = -1.;
    double mainqn = -1.,mainpercqn = -1.;
    TSpline3* qnspline = nullptr;

    if(fFlowMethod==kEvShapeEP || fFlowMethod==kEvShapeSP || fFlowMethod==kEvShapeEPVsMass) {
    //cout << "qnFullTPC0 = " << qnFullTPC << endl;
	    HFQnVectorHandler->GetqnTPC(qnFullTPC,qnPosTPC,qnNegTPC);
    //cout << "qnFullTPC1 = " << qnFullTPC << endl;
	    HFQnVectorHandler->GetqnV0(qnFullV0,qnV0A,qnV0C);
    }

    //cout << "qnFullTPC = " << qnFullTPC << endl;
    //cout << "qnPosTPC = " << qnPosTPC << endl;
    //cout << "qnNegTPC = " << qnNegTPC << endl;
    //cout << "qnFullV0 = " << qnFullV0 << endl;
    //cout << "qnV0A = " << qnV0A << endl;
    //cout << "qnV0C = " << qnV0C << endl;

    TAxis* centraxis = fHistqnVsCentrPercCalib[0]->GetXaxis();
    int centrbin = centraxis->FindBin(evCentr);
    double centrbinmin = centraxis->GetBinLowEdge(centrbin);
    double centrbinmax = centrbinmin+centraxis->GetBinWidth(centrbin);

    //cout << "fqnMeth =" << fqnMeth << endl;
    //cout << "kq2TPC = " << kq2TPC << endl;
    //cout << "qnFullTPC = " << qnFullTPC << endl;


    switch(fqnMeth){

	    case kq2TPC:
	    mainqn = qnFullTPC;
            //cout << "case, fqnMeth =" << fqnMeth << endl;
            //cout << "case, mainqn =" << mainqn << endl;
	    if(fPercentileqn) qnspline = static_cast<TSpline3*>(fqnSplinesList[0]->FindObject(Form("sq2Int_centr_%0.f_%0.f",centrbinmin,centrbinmax)));
	    break;

	    case kq2PosTPC:
	    mainqn = qnPosTPC;
	    if(fPercentileqn) qnspline = static_cast<TSpline3*>(fqnSplinesList[1]->FindObject(Form("sq2Int_centr_%0.f_%0.f",centrbinmin,centrbinmax)));
	    break;

	    case kq2NegTPC:
	    mainqn = qnNegTPC;
	    if(fPercentileqn) qnspline = static_cast<TSpline3*>(fqnSplinesList[2]->FindObject(Form("sq2Int_centr_%0.f_%0.f",centrbinmin,centrbinmax)));
	    break;

	    case kq2VZERO:
	    mainqn = qnFullV0;
	    if(fPercentileqn) qnspline = static_cast<TSpline3*>(fqnSplinesList[3]->FindObject(Form("sq2Int_centr_%0.f_%0.f",centrbinmin,centrbinmax)));
	    break;

	    case kq2VZEROA:
	    mainqn = qnV0A;
	    if(fPercentileqn) qnspline = static_cast<TSpline3*>(fqnSplinesList[4]->FindObject(Form("sq2Int_centr_%0.f_%0.f",centrbinmin,centrbinmax)));
	    break;

	    case kq2VZEROC:
	    mainqn = qnV0C;
	    if(fPercentileqn) qnspline = static_cast<TSpline3*>(fqnSplinesList[5]->FindObject(Form("sq2Int_centr_%0.f_%0.f",centrbinmin,centrbinmax)));
	    break;
    }

    if(fPercentileqn){
	    if(qnspline)
		    mainpercqn = qnspline->Eval(mainqn);
	    else
		    AliWarning("Centrality binning and centrality intervals of qn splines do not match!");
    }
    else{
	    mainpercqn = mainqn;
    }

    //cout << "mainqn = " << mainqn << endl;

    //phi axis of EP is changed from 0 ~ Pi to -Pi/2 ~ Pi/2
    if( TMath::Pi()/2. < PsinFullTPC && PsinFullTPC  < TMath::Pi() ){
	    PsinFullTPC = PsinFullTPC - TMath::Pi();
    }

    if( TMath::Pi()/2. < PsinFullV0 && PsinFullV0  < TMath::Pi() ){
	    PsinFullV0 = PsinFullV0 - TMath::Pi();
    }

    if( TMath::Pi()/2. < PsinPosTPC && PsinPosTPC  < TMath::Pi() ){
	    PsinPosTPC = PsinPosTPC - TMath::Pi();
    }

    if( TMath::Pi()/2. < PsinNegTPC && PsinNegTPC  < TMath::Pi() ){
	    PsinNegTPC = PsinNegTPC - TMath::Pi();
    }

    if( TMath::Pi()/2. < PsinV0A && PsinV0A  < TMath::Pi() ){
	    PsinV0A = PsinV0A - TMath::Pi();
    }

    if( TMath::Pi()/2. < PsinV0C && PsinV0C  < TMath::Pi() ){
	    PsinV0C = PsinV0C - TMath::Pi();
    }
    //cout << "PsinFullTPC = " <<PsinFullTPC << endl;

    //Fill event-based histograms
//EP
    fHistEvPlane[0]->Fill(PsinFullTPC);
    fHistEvPlane[1]->Fill(PsinPosTPC);
    fHistEvPlane[2]->Fill(PsinNegTPC);
    fHistEvPlane[3]->Fill(PsinFullV0);
    fHistEvPlane[4]->Fill(PsinV0A);
    fHistEvPlane[5]->Fill(PsinV0C);


    //centrality vs EP
    fHistEvPlaneQncorr[0]->Fill(evCentr,PsinFullTPC);
    fHistEvPlaneQncorr[1]->Fill(evCentr,PsinPosTPC);
    fHistEvPlaneQncorr[2]->Fill(evCentr,PsinNegTPC);
    fHistEvPlaneQncorr[3]->Fill(evCentr,PsinFullV0);
    fHistEvPlaneQncorr[4]->Fill(evCentr,PsinV0A);
    fHistEvPlaneQncorr[5]->Fill(evCentr,PsinV0C);

    fEPcorV0AC->Fill(PsinV0A,PsinV0C);    //correlation of V0A and V0C
    fEPcorV0ATPC->Fill(PsinV0A,PsinFullTPC);  //correlation of V0A and TPC
    fEPcorV0CTPC->Fill(PsinV0C,PsinFullTPC); //correlation of V0C and TPC

    Double_t subV0AC = -999.;
    Double_t subV0ATPC = -999.;
    Double_t subV0CTPC = -999.;

    subV0AC = PsinV0A - PsinV0C;
    subV0ATPC = PsinV0A - PsinFullTPC;
    subV0CTPC = PsinV0C - PsinFullTPC;

    Double_t subV0ACcos2 = -999.;
    Double_t subV0ATPCcos2 = -999.;
    Double_t subV0CTPCcos2 = -999.;

    subV0ACcos2 = cos(2*subV0AC);
    fsubV0ACcos2 -> Fill(evCentr,subV0ACcos2);

    subV0ATPCcos2 = cos(2*subV0ATPC);
    fsubV0ATPCcos2 -> Fill(evCentr,subV0ATPCcos2);

    subV0CTPCcos2 = cos(2*subV0CTPC);
    fsubV0CTPCcos2 -> Fill(evCentr,subV0CTPCcos2);


    //centrality vs mainpercqn vs EP
    //fHistEvPlaneQncorr[0]->Fill(evCentr,mainpercqn,PsinFullTPC);
    //fHistEvPlaneQncorr[1]->Fill(evCentr,mainpercqn,PsinPosTPC);
    //fHistEvPlaneQncorr[2]->Fill(evCentr,mainpercqn,PsinNegTPC);
    //fHistEvPlaneQncorr[3]->Fill(evCentr,mainpercqn,PsinFullV0);
    //fHistEvPlaneQncorr[4]->Fill(evCentr,mainpercqn,PsinV0A);
    //fHistEvPlaneQncorr[5]->Fill(evCentr,mainpercqn,PsinV0C);

    //cout << "PsinFullV0 =" << PsinFullV0 << endl;
    //cout << "PsinV0A =" << PsinV0A << endl;
    //cout << "evCentr = " << evCentr << endl;
    //cout << "mainpercqn = " << mainpercqn << endl;

    //EP / Qn resolution histograms
    if(fFlowMethod==kEP || fFlowMethod==kEvShapeEP || fFlowMethod==kEPVsMass || fFlowMethod==kEvShapeEPVsMass){
	    fHistEPResolVsCentrVsqn[0]->Fill(evCentr,TMath::Cos(fHarmonic*GetDeltaPsiSubInRange(SubAPsin,SubBPsin)));
	    //fHistEPResolVsCentrVsqn[0]->Fill(evCentr,mainpercqn,TMath::Cos(fHarmonic*GetDeltaPsiSubInRange(SubAPsin,SubBPsin)));
	    if(nsubevents==3){
		    fHistEPResolVsCentrVsqn[1]->Fill(evCentr,TMath::Cos(fHarmonic*GetDeltaPsiSubInRange(SubAPsin,SubCPsin)));
		    //fHistEPResolVsCentrVsqn[1]->Fill(evCentr,mainpercqn,TMath::Cos(fHarmonic*GetDeltaPsiSubInRange(SubAPsin,SubCPsin)));
		    fHistEPResolVsCentrVsqn[2]->Fill(evCentr,TMath::Cos(fHarmonic*GetDeltaPsiSubInRange(SubBPsin,SubCPsin)));
		    //fHistEPResolVsCentrVsqn[2]->Fill(evCentr,mainpercqn,TMath::Cos(fHarmonic*GetDeltaPsiSubInRange(SubBPsin,SubCPsin)));

	    }
    }

    else if(fFlowMethod==kSP || fFlowMethod==kEvShapeSP){

	    fHistEPResolVsCentrVsqn[0]->Fill(evCentr,(SubAQn[0]*SubBQn[0]+SubAQn[1]*SubBQn[1])/(SubAMultQn * SubBMultQn));
	    //fHistEPResolVsCentrVsqn[0]->Fill(evCentr,mainpercqn,(SubAQn[0]*SubBQn[0]+SubAQn[1]*SubBQn[1])/(SubAMultQn * SubBMultQn));

	    if(nsubevents==3){

		    fHistEPResolVsCentrVsqn[1]->Fill(evCentr,(SubAQn[0]*SubCQn[0]+SubAQn[1]*SubCQn[1])/(SubAMultQn * SubCMultQn));
		    //fHistEPResolVsCentrVsqn[1]->Fill(evCentr,mainpercqn,(SubAQn[0]*SubCQn[0]+SubAQn[1]*SubCQn[1])/(SubAMultQn * SubCMultQn));
		    fHistEPResolVsCentrVsqn[2]->Fill(evCentr,(SubBQn[0]*SubCQn[0]+SubBQn[1]*SubCQn[1])/(SubBMultQn * SubCMultQn));
		    //fHistEPResolVsCentrVsqn[2]->Fill(evCentr,mainpercqn,(SubBQn[0]*SubCQn[0]+SubBQn[1]*SubCQn[1])/(SubBMultQn * SubCMultQn));
	    }
    }

    iBevt = kFALSE;
    if(fMCarray)CheckMCgen(fMCheader,0.6);
    //cout << "Bevt 0 = " << iBevt << endl;

    //cout << "--------- NembMCpi0 ; NembMCeta " << endl;
    //cout << NembMCpi0 << " ;  " << NembMCeta << endl;

////////////////////////////// EMCAL cluster loop////////////////////////////////////

    Int_t Nclust = fVevent->GetNumberOfCaloClusters();

    int NclustAll = 0;
    int NclustE1 = 0;//Number of cluster E>0.1
    int NclustE2 = 0;//Number of cluster E>0.2
    int NclustE3 = 0;//Number of cluster E>0.5

    for(Int_t icl=0; icl<Nclust; icl++)
    {
	    AliVCluster *clust = 0x0;
	    //clust = (AliVCluster*)fVevent->GetCaloCluster(icl); // address cluster matched to track
	    clust = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(icl));

	    if(clust && clust->IsEMCAL())
	    {
		    Double_t clustE = clust->E();

		    //Select EMCAL or DCAL clusters
		    Float_t emcx[3];//cluster pos
		    clust->GetPosition(emcx);
		    TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
		    Double_t emcphi = clustpos.Phi();
		    Double_t emceta = clustpos.Eta();

		    fEMCClsEtaPhi->Fill(emceta,emcphi);

		    fHistClustE->Fill(clustE);

		    fHistNoCells->Fill(clustE,clust->GetNCells());

		    // printf("%f \n",clustE);


		    NclustAll++;
		    if(clustE>0.1)NclustE1++;
		    if(clustE>0.2)NclustE2++;
		    if(clustE>0.5)NclustE3++;

	    }

    }

    fHistNCls->Fill(NclustAll);
    fHistNClsE1->Fill(NclustE1);
    fHistNClsE2->Fill(NclustE2);
    fHistNClsE3->Fill(NclustE3);

    //cell information
    AliVCaloCells *fCaloCells = fVevent->GetEMCALCells();

    Short_t cellAddr, nSACell;
    Int_t mclabel;
    Short_t iSACell;
    Double_t cellAmp=-1., cellTimeT=-1., clusterTime=-1., efrac=-1.;

    nSACell = fCaloCells->GetNumberOfCells();
    for(iSACell = 0; iSACell < nSACell; iSACell++){
	    Bool_t haveCell = fCaloCells->GetCell(iSACell,cellAddr,cellAmp,cellTimeT,mclabel,efrac);//cellの様々な情報を取ってくる
	    if(haveCell)fHistCalCell->Fill(cellAddr,cellAmp);
    }


    /////////////////////////////////Track loop///////////////////////////////


    //Int_t iTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event
    Int_t iTracks(fTracks_tender->GetEntries());           // see how many tracks there are in the event
    for(Int_t i = 0; i < iTracks; i++) {                 // loop over all these tracks

	    //AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
	    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fTracks_tender->At(i));         // get a track (type AliAODTrack) from the event

	    if(!track) continue;                            // if we failed, skip this track
	    fHistPt->Fill(track->Pt());                     // plot the pt value of the track in a histogram

	    //////////////Apply track cuts//////////
	    if(track->GetTPCNcls() < 80)continue;
	    if(track->GetITSNcls() < 3)continue;
	    if((!(track->GetStatus()&AliESDtrack::kITSrefit) || (!(track->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
	    if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)))continue;

	    //////////////Track properties//////////
	    Bool_t fFlagNonHFE=kFALSE;

	    Double_t dEdx=-999, fTPCnSigma=-999, fTOFnSigma=-999, fITSnSigma=-999;
	    Double_t TrkPhi=-999, TrkPt=-999,TrkEta=-999,TrkP=-999;

	    TrkPhi = track->Phi();
	    TrkPt = track->Pt();
	    TrkEta = track->Eta();
	    if(TrkEta < -0.6 || TrkEta > 0.6)continue;
	    TrkP = track->P();
	    dEdx = track->GetTPCsignal();
	    fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
	    fTOFnSigma = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
	    fITSnSigma = fpidResponse->NumberOfSigmasITS(track, AliPID::kElectron);

	    fTrkPt -> Fill(TrkPt);
	    fTrkphi -> Fill(TrkPhi);
	    fTrketa -> Fill(TrkEta);
	    fTrkP -> Fill(TrkP);
	    fdEdx -> Fill(TrkP,dEdx);
	    fTPCnsig -> Fill(TrkP,fTPCnSigma);
	    fTOFnsig -> Fill(TrkP,fTOFnSigma);
	    fITSnsig -> Fill(TrkP,fITSnSigma);
	    fTPCnsig_TOFnsig -> Fill(fTOFnSigma,fTPCnSigma);

	    fvalueElectron[0] = TrkPt;
	    fvalueElectron[1] = TrkP;
	    fvalueElectron[2] = fTPCnSigma;
	    fvalueElectron[3] = fITSnSigma;
	    fvalueElectron[4] = fTOFnSigma;
	    fvalueElectron[5] = -1.0;

	    ////Charged Particle v2////

	    Double_t TrkPhiPI = -999.;

	    TrkPhiPI = TVector2::Phi_mpi_pi(TrkPhi);

	    //Double_t TrkPhiEPTPC = -999.;
	    Double_t TrkPhiEPV0A = -999.;
	    Double_t TrkPhiEPV0A_ele = -999.;
	    Double_t TrkPhiEPV0A_had = -999.;
	    Double_t TrkPhiEPV0A_ele_lowpt = -999.;

	    //Double_t TrkPhiEP2 = -999;
	    //Double_t const PI = TMath::Pi();
	    //Double_t const PI_2 = TMath::Pi()/2.;
	    //Double_t const PI_4 = TMath::Pi()/4.;


	    Double_t TrkPhicos2 = -999;
	    Double_t TrkPhisin2 = -999;

	    ////if(track->Pt() > 1.5){

	    //TrkPhiEPTPC = TrkPhiPI - PsinFullTPC;
	    TrkPhiEPV0A = TrkPhiPI - PsinV0A;

	    //if(TrkPhiEPTPC > TMath::Pi()/2.){
	    //	TrkPhiEPTPC -= TMath::Pi();
	    //}

	    //if(TrkPhiEPTPC < -TMath::Pi()/2.){
	    //	TrkPhiEPTPC += TMath::Pi();
	    //}

	    if(TrkPhiEPV0A > TMath::Pi()/2.){
		    TrkPhiEPV0A -= TMath::Pi();
	    }

	    if(TrkPhiEPV0A < -TMath::Pi()/2.){
		    TrkPhiEPV0A += TMath::Pi();
	    }

	    //cout << "TrkPhiEPTPC =" << TrkPhiEPTPC << endl;
	    //cout << "TrkPhiEPV0A =" << TrkPhiEPV0A << endl;

	    //fTrkPhiEPFullTPC -> Fill(TrkPhiEPTPC);
	    //fTrkPhiEPFullV0A -> Fill(TrkPhiEPV0A);

	    //fTrkPhiEPFullTPC_Pt -> Fill(track->Pt(),TrkPhiEPTPC);
	    fTrkPhiEPV0A_Pt -> Fill(track->Pt(),TrkPhiEPV0A);

	    //fcorTrkPhicent_charge -> Fill(TrkPhiEP,lPercentile);

	    TrkPhicos2 = cos(2*TrkPhiEPV0A);
	    TrkPhisin2 = sin(2*TrkPhiEPV0A);

	    fTrkPhicos2 -> Fill(track->Pt(),TrkPhicos2);
	    fTrkPhisin2 -> Fill(track->Pt(),TrkPhisin2);

	    //fcorcentcos2_charge -> Fill(lPercentile,TrkPhicos2);


	    //if(TrkPhiEP > -PI_4 && TrkPhiEP < PI_4){

	    //	fInplane -> Fill(track->Pt());
	    //	fcorcentInplane -> Fill(lPercentile,track->Pt());

	    //}

	    //if(TrkPhiEP > PI_4 || TrkPhiEP < -PI_4){

	    //	fOutplane -> Fill(track->Pt());
	    //	fcorcentOutplane -> Fill(lPercentile,track->Pt());

	    //}
	    //////Charged Particle v2 each Pt////

	    //fTrkPhiEP_Pt -> Fill(track->Pt(),TrkPhiEP);
	    //}

	    //Event cut for eta distribution
	    //if(TMath::Abs(TrkPt)>1.0)return;
	    //fTrketa2 ->Fill(TrkEta);

	    Int_t EMCalIndex = -1;
	    EMCalIndex = track->GetEMCALcluster();  // get index of EMCal cluster which matched to track
	    //cout << "EMCalIndex = " << EMCalIndex << endl;

	    //cout << "EMCal Index = " << EMCalIndex << endl;

	    /////track cut/////

	    //===DCA cut===
	    Double_t DCA[2] = {-999.,-999.}, covar[3];
	    if(track -> PropagateToDCA(pVtx,fVevent -> GetMagneticField(),20.,DCA,covar))
	    {

		    //cout << "DCA[0] = "<< DCA[0] << endl;
		    //cout << "DCA[1] = "<< DCA[1] << endl;

		    //if(TMath::Abs(DCA[0]) > CutDCAxy || TMath::Abs(DCA[1]) > CutDCAz)continue;
		    if(TMath::Abs(DCA[0]) > 2.4 || TMath::Abs(DCA[1]) > 3.2)continue;

	    }

	    //cout << "EMCal Index = " << EMCalIndex << endl;

	    //////Get MC information///////
	    Int_t ilabel = TMath::Abs(track->GetLabel());
	    Int_t pdg = -999;
	    Double_t pid_ele = 0.0;
	    Double_t pTmom = -1.0;
	    Int_t pidM = -1;
	    Int_t ilabelM = -1;

	    Bool_t pid_eleD = kFALSE;
	    Bool_t pid_eleB = kFALSE;
	    Bool_t pid_eleP = kFALSE;

	    Int_t ilabelGM = -1;
	    Int_t pidGM = -1;
	    Double_t pTGMom = -1;

	    Bool_t iEmbPi0 = kFALSE;
	    Bool_t iEmbEta = kFALSE;


	    if(ilabel>0 && fMCarray){

		    fMCTrackpart = (AliAODMCParticle*) fMCarray->At(ilabel);

		    pdg = fMCTrackpart->GetPdgCode();

		    if(TMath::Abs(pdg)==11)pid_ele = 1.0;

		    if(pid_ele==1.0)FindMother(fMCTrackpart, ilabelM, pidM, pTmom);

		    pid_eleB = IsBdecay(pidM);

		    pid_eleP = IsPdecay(pidM);

		    pid_eleD = IsDdecay(pidM);
                    //if(pidM==4122)cout << "pidM = " << pidM << " ; pid_eleD = " << pid_eleD << endl;

		    if(pid_eleD || pid_eleB)fNDB->Fill(0);

		    if(pid_eleD){

			    AliAODMCParticle* fMCTrackpartMom = (AliAODMCParticle*) fMCarray->At(ilabelM);

			    FindMother(fMCTrackpartMom,ilabelGM,pidGM,pTGMom);

			    if(IsBdecay(pidGM)){
				    //cout << pid_eleD << " ; "<< pid_eleB << "; before pidM = " << pidM << ";   pTmom = " << pTmom <<endl;
				    pid_eleB = IsBdecay(pidGM);

				    pid_eleD = kFALSE;

				    pidM = pidGM;

				    pTmom = pTGMom;

				    //cout << pid_eleD << " ; "<< pid_eleB << "; afetr pidM = " << pidM << ";   pTmom = " << pTmom <<endl;
			    }

		    } //pid_eleD

                    if(iBevt)pid_eleD = kFALSE;
                    //cout << "pidM = " << pidM << " ; " << pid_eleD << endl;

		    if(pid_eleD || pid_eleB)fNDB->Fill(1);

		    if(pidM==111){
			    if(ilabelM>=NembMCpi0 && ilabelM<NembMCeta)iEmbPi0 = kTRUE;
			    if(ilabelM>=NembMCeta && ilabelM<NpureMCproc)iEmbEta = kTRUE;

		    }
		    if(pidM==221){

			    if(ilabelM>=NembMCeta && ilabelM<NpureMCproc)iEmbEta = kTRUE;

		    }

		    if(pidM==22){

			    AliAODMCParticle* fMCparticleM = (AliAODMCParticle*) fMCarray->At(ilabelM);
			    FindMother(fMCparticleM, ilabelM, pidM, pTmom);

			    if(pidM==111){
				    if(ilabelM>=NembMCpi0 && ilabelM<NembMCeta)iEmbPi0 = kTRUE;
				    if(ilabelM>=NembMCeta && ilabelM<NpureMCproc)iEmbEta = kTRUE;

			    }

			    if(pidM==221){

				    if(ilabelM>=NembMCeta && ilabelM<NpureMCproc)iEmbEta = kTRUE;

			    }

		    } //pidM = 22

		    fMCcheckMother->Fill(abs(pidM));
	    }


	    //if(pidM==443)continue;
	    //if(pidM==-99)continue;

	    if(pid_eleB || pid_eleD) {

		    fHist_eff_HFE->Fill(TrkPt);

		    if(fTPCnSigma>-1 && fTPCnSigma<3) fHist_eff_TPC->Fill(TrkPt);

	    }


	    ///// calucurate weight of photon for MC /////
	    Double_t WeightPho = -999.0;

	    if(iEmbPi0){

		    if(iCentral)
		    {
			    if(pTmom<6.75)
			    {
				    WeightPho = fPi010_0->Eval(pTmom);
			    }
			    else
			    {
				    WeightPho = fPi010_1->Eval(pTmom);
			    }
		    }

		    if(iSemiCentral)
		    {
			    if(pTmom<4.0)
			    {
				    WeightPho = fPi3050_0 -> Eval(pTmom);
			    }
			    else
			    {
				    WeightPho = fPi3050_1 -> Eval(pTmom);
			    }
		    }
	    }

	    if(iEmbEta){

		    if(iCentral)
		    {
			    WeightPho = fEta010->Eval(pTmom);
		    }
		    if(iSemiCentral)
		    {
			    WeightPho = fEta3050 -> Eval(pTmom);
		    }
	    }

	    //////// calculation of electron v2 at low pt using TPC and TOF info. //////

	    Double_t TrkPhicos2_elelow = -999;
	    Double_t TrkPhisin2_elelow = -999;

	    if(fTOFnSigma > -3 && fTOFnSigma < 3){ //TOF nsigma cut

		    fHistele_TOFcuts->Fill(fTPCnSigma,track->Pt());

		    if(fTPCnSigma > 0 && fTPCnSigma <3){ //TPC nsigma cut

			    TrkPhiEPV0A_ele_lowpt = TrkPhiPI - PsinV0A;

			    if(TrkPhiEPV0A_ele_lowpt > TMath::Pi()/2.){
				    TrkPhiEPV0A_ele_lowpt -= TMath::Pi();
			    }

			    if(TrkPhiEPV0A_ele_lowpt < -TMath::Pi()/2.){
				    TrkPhiEPV0A_ele_lowpt += TMath::Pi();
			    }

			    fTrkPhiEPV0A_Pt_ele_lowpt -> Fill(track->Pt(),TrkPhiEPV0A_ele_lowpt);

			    TrkPhicos2_elelow = cos(2*TrkPhiEPV0A_ele_lowpt);
			    TrkPhisin2_elelow = sin(2*TrkPhiEPV0A_ele_lowpt);

			    fTrkPhicos2_elelow -> Fill(track->Pt(),TrkPhicos2_elelow);
			    fTrkPhisin2_elelow -> Fill(track->Pt(),TrkPhisin2_elelow);
		    }

	    }

	    if(fTOFnSigma < -3.5){ //TOF nsigma cut

		    fHisthad_TOFcuts->Fill(fTPCnSigma,track->Pt());

	    }

	    //////Track matching to EMCAL//////

	    AliVCluster *clustMatch=0x0;
	    //cout << "EMCalIndex = " << EMCalIndex << endl;
	    //if(EMCalIndex>=0)clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex); // address cluster matched to track
	    if(EMCalIndex>=0) clustMatch = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(EMCalIndex));

	    //cout << "Charge = " << track -> Charge() << endl;

	    if(clustMatch && clustMatch->IsEMCAL())
	    {
		    Double_t clustMatchE = clustMatch->E();
		    if(track->P()>0) fvalueElectron[5] = clustMatchE/track->P();
	    }
	    if(iTree)fSparseElectron->Fill(fvalueElectron);



	    Double_t emcphi = -999, emceta = -999;

	    if(clustMatch && clustMatch->IsEMCAL())
	    {

		    //shower shape cut
		    Double_t fPhiDiff = -999, fEtaDiff = -999;
		    GetTrkClsEtaPhiDiff(track,clustMatch,fPhiDiff,fEtaDiff);
		    fEMCTrkMatchPhi->Fill(fPhiDiff);
		    fEMCTrkMatchEta->Fill(fEtaDiff);

		    //cout << "fPhiDiff = "<< fPhiDiff << endl;

		    if(TMath::Abs(fPhiDiff)>0.05 || TMath::Abs(fEtaDiff)>0.05)continue;


		    //Select EMCAL or DCAL clusters
		    Float_t emcx[3];//cluster pos
		    clustMatch->GetPosition(emcx);
		    TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
		    emcphi = clustpos.Phi();
		    emceta = clustpos.Eta();

		    fClsEtaPhiAftMatch->Fill(emceta,emcphi);


		    //////Properties of tracks matched to the EMCAL/////
		    fEMCTrkPt->Fill(TrkPt);

		    if(TrkPt>1.0){
			    fEMCTrketa->Fill(TrkEta);
			    fEMCTrkphi->Fill(TrkPhi);
		    }

		    Double_t clustMatchE = clustMatch->E();


		    /////EMCAL EID info///////
		    Double_t eop = -1.0;
		    Double_t m02 = -99999, m20 = -99999;

		    if(track->P()>0)eop = clustMatchE/track->P();
		    m02 = clustMatch->GetM02();
		    m20 = clustMatch->GetM20();

		    //add by sudo

		    if(track->Pt()>2.0){
			    fHistNsigEop->Fill(eop,fTPCnSigma);
		    }

		    fM20->Fill(track->Pt(),clustMatch->GetM20());
		    fM02->Fill(track->Pt(),clustMatch->GetM02());

		    //cout << "eop" << eop << endl;

		    Double_t TrkPhicos2_elehigh = -999;
		    Double_t TrkPhisin2_elehigh = -999;

		    //if((fTPCnSigma > -1 && fTPCnSigma <3) && (m20 > 0.01 && m20 < 0.3)){ //TPC nsigma & shower shape cut
		    if(fTPCnSigma > ftpcnsig && fTPCnSigma <3){ //TPC nsigma & shower shape cut

			    if(track->Pt()<3.0){
				    //if( (fTOFnSigma<-1 || fTOFnSigma>1) || (fITSnSigma<-3 || fITSnSigma>1) || (fTPCnSigma<0) )continue;
				    if( (fTOFnSigma<-1 || fTOFnSigma>1) || (fITSnSigma<-3 || fITSnSigma>1) )continue;
			    }


			    //if(eop>0.9 && eop<1.3){ //eop cut
			    if(eop>femceop && eop<1.3 && (m20 > femcss_mim && m20 < femcss_max)){ //eop cut

				    //SelectPhotonicElectron(iTracks,track,fFlagNonHFE,TrkPt,DCAxy,Bsign,TrkPhiPI,PsinV0A);
				    SelectPhotonicElectron(iTracks,track,fFlagNonHFE,TrkPt,DCA[0],Bsign,TrkPhiPI,PsinV0A);
				    ////electron v2////

				    //if(track->Pt() > 1.5){

				    TrkPhiEPV0A_ele = TrkPhiPI - PsinV0A;

				    if(TrkPhiEPV0A_ele > TMath::Pi()/2.){
					    TrkPhiEPV0A_ele -= TMath::Pi();
				    }

				    if(TrkPhiEPV0A_ele < -TMath::Pi()/2.){
					    TrkPhiEPV0A_ele += TMath::Pi();
				    }
				    //if(TrkPhiEP2 > PI_2){
				    //	TrkPhiEP2 -= PI;
				    //}

				    //if(TrkPhiEP2 < -PI_2){
				    //	TrkPhiEP2 += PI;
				    //}
				    //fTrkPhiEP2 -> Fill(TrkPhiEP2);
				    //fcorTrkPhicent_ele -> Fill(TrkPhiEP2,lPercentile);

				    ////electron v2 each Pt/////

				    fTrkPhiEPV0A_Pt_ele -> Fill(track->Pt(),TrkPhiEPV0A_ele);

				    TrkPhicos2_elehigh = cos(2*TrkPhiEPV0A_ele);
				    TrkPhisin2_elehigh = sin(2*TrkPhiEPV0A_ele);

				    fTrkPhicos2_elehigh -> Fill(track->Pt(),TrkPhicos2_elehigh);
				    fTrkPhisin2_elehigh -> Fill(track->Pt(),TrkPhisin2_elehigh);
				    //}

				    if(!fFlagNonHFE){

					    fTrkPhicos2_hfehigh -> Fill(track->Pt(),TrkPhicos2_elehigh);
					    fTrkPhisin2_hfehigh -> Fill(track->Pt(),TrkPhisin2_elehigh);
					    fDCAxy_Pt_hfe -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());

				    }


				    if(TrkPhiEPV0A_ele > -TMath::Pi()/4. && TrkPhiEPV0A_ele < TMath::Pi()/4.){

					    fInplane_ele -> Fill(track->Pt());
					    fDCAxy_Pt_Inplane_ele -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());

					    if(!fFlagNonHFE){

						    fInplane_hfe -> Fill(track->Pt());//using this at pt < 3 GeV
						    fDCAxy_Pt_Inplane_hfe -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());

					    }


				    }



				    if(TrkPhiEPV0A_ele > TMath::Pi()/4. || TrkPhiEPV0A_ele < -TMath::Pi()/4.){

					    fOutplane_ele -> Fill(track->Pt());
					    fDCAxy_Pt_Outplane_ele -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());

					    if(!fFlagNonHFE){

						    fOutplane_hfe -> Fill(track->Pt());//using this at pt < 3 GeV
						    fDCAxy_Pt_Outplane_hfe -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());

					    }

				    }


				    /////DCA distribution////
				    //fDCAxy_Pt_ele->Fill(TrkPt,DCA[0]*Bsign*track->Charge());

				    /////Identify Non-HFE/////
				    //SelectPhotonicElectron(iTracks,track,fFlagNonHFE,TrkPt,DCAxy,Bsign,TrkPhiPI,PsinV0A);

				    if(pid_eleP){

					    fHistPhoReco0 -> Fill(track->Pt());

					    //cout << "iSemiCentral = " << iSemiCentral << " ; iEmbPi0 = " << iEmbPi0 << " ; WeightPho1 = " << WeightPho << endl;

					    if(iEmbPi0)fHistPhoPi0->Fill(track->Pt(),WeightPho);
					    if(iEmbEta)fHistPhoEta->Fill(track->Pt(),WeightPho);


					    if(fFlagNonHFE){
						    fHistPhoReco1 -> Fill(track->Pt());
						    if(iEmbPi0)fHistPhoPi01->Fill(track->Pt(),WeightPho);
						    if(iEmbEta)fHistPhoEta1->Fill(track->Pt(),WeightPho);

					    }

					    else{
						    fHistPhoReco2 -> Fill(track->Pt());
					    }
				    }

			    }

			    fHisteop->Fill(eop);
			    fHistelectron->Fill(eop,track->Pt());
		    }

		    Double_t dWeight = -99;
		    Double_t bWeight = -99;
		    if(fTPCnSigma > ftpcnsig && fTPCnSigma <3 && m20 > 0.02 && m20 < 0.3){ // TPC nsigma & shower shape cut

			    if(eop>femceop && eop<1.3){ // E/p cut
				    if(pid_eleB) fHistPt_HFE_MC_B -> Fill(track->Pt());
				    if(pid_eleD) fHistPt_HFE_MC_D -> Fill(track->Pt());

				    fDCAxy_Pt_ele->Fill(TrkPt,DCA[0]*Bsign*track->Charge());
				    //D meson : pidM>400, <499, =421, =413, >430, <436 D+ : =411 Ds : =431 Lc : 4122 from AliAnalysisTaskTPCCalBeauty.cxx
				    
                                    if(pid_eleD)
				    {
					    if(TMath::Abs(pidM)> 400 || TMath::Abs(pidM)< 499 || TMath::Abs(pidM)== 421 || TMath::Abs(pidM)== 413 || TMath::Abs(pidM)> 430 || TMath::Abs(pidM)< 436){//if from D meson
						    if (pTmom>1 && pTmom<50.) { //in proper pt range

							    dWeight = fDWeightNew->GetBinContent(fDWeightNew->FindBin(pTmom));
							    fDCAxy_Pt_D_WeightNew->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),dWeight);

							    dWeight = fDWeightVar1->GetBinContent(fDWeightVar1->FindBin(pTmom));
							    fDCAxy_Pt_D_WeightVar1->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),dWeight);

							    dWeight = fDWeightVar2->GetBinContent(fDWeightVar2->FindBin(pTmom));
							    fDCAxy_Pt_D_WeightVar2->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),dWeight);

							    fDCAxy_Pt_D->Fill(TrkPt,DCA[0]*Bsign*track->Charge());
						    }
					    }if (TMath::Abs(pidM)== 411) { //if from D+ meson
						    if (pTmom>1 && pTmom<50.) { //in proper pt range

							    dWeight = fDWeightNew->GetBinContent(fDWeightNew->FindBin(pTmom));
							    fDCAxy_Pt_Dpm_WeightNew->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),dWeight);

							    dWeight = fDPlusWeightVar1->GetBinContent(fDPlusWeightVar1->FindBin(pTmom));
							    fDCAxy_Pt_Dpm_WeightVar1->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),dWeight);

							    dWeight = fDWeightVar2->GetBinContent(fDWeightVar2->FindBin(pTmom));
							    fDCAxy_Pt_Dpm_WeightVar2->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),dWeight);

							    fDCAxy_Pt_Dpm->Fill(TrkPt,DCA[0]*Bsign*track->Charge());
						    }

					    }if (TMath::Abs(pidM)== 421) { //if from D0 meson
						    if (pTmom>1 && pTmom<50.) { //in proper pt range

							    dWeight = fDWeightNew->GetBinContent(fDWeightNew->FindBin(pTmom));
							    fDCAxy_Pt_D0_WeightNew->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),dWeight);
							    fDCAxy_Pt_D0->Fill(TrkPt,DCA[0]*Bsign*track->Charge());
						    }

					    }if (TMath::Abs(pidM)== 431) { //if from Ds meson
						    if (pTmom>1 && pTmom<50.) { //in proper pt range

							    dWeight = fDWeightNew->GetBinContent(fDWeightNew->FindBin(pTmom));
							    fDCAxy_Pt_Ds_WeightNew->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),dWeight);

							    dWeight = fDsWeightVar1->GetBinContent(fDsWeightVar1->FindBin(pTmom));
							    fDCAxy_Pt_Ds_WeightVar1->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),dWeight);

							    dWeight = fDWeightVar2->GetBinContent(fDWeightVar2->FindBin(pTmom));
							    fDCAxy_Pt_Ds_WeightVar2->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),dWeight);

							    fDCAxy_Pt_Ds->Fill(TrkPt,DCA[0]*Bsign*track->Charge());
						    }

					    }if (TMath::Abs(pidM)==4122) { //if from Lc
						    if (pTmom>1 && pTmom<50.) { //in proper pt range
							    dWeight = fDWeightNew->GetBinContent(fDWeightNew->FindBin(pTmom));
							    fDCAxy_Pt_lambda_WeightNew->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),dWeight);

							    dWeight = fLcWeightVar1->GetBinContent(fLcWeightVar1->FindBin(pTmom));
							    fDCAxy_Pt_lambda_WeightVar1->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),dWeight);

							    dWeight = fLcWeightVar2->GetBinContent(fLcWeightVar2->FindBin(pTmom));
							    fDCAxy_Pt_lambda_WeightVar2->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),dWeight);

							    fDCAxy_Pt_lambda->Fill(TrkPt,DCA[0]*Bsign*track->Charge());
						    }

                                     } // end if pid_eleD


				    }if (pid_eleB && (TMath::Abs(pidM)> 500 || TMath::Abs(pidM)< 599) ) {//if from B meson
					    //cout<<"TESTING5"<<endl;
					    if (pTmom>0. && pTmom<50.) { //in proper pt range

						    bWeight = fBWeightNew->GetBinContent(fBWeightNew->FindBin(pTmom));
						    fDCAxy_Pt_Bmeson_WeightNew->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),bWeight);

						    bWeight = fBWeightVar1->GetBinContent(fBWeightVar1->FindBin(pTmom));
						    fDCAxy_Pt_Bmeson_WeightVar1->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),bWeight);

						    bWeight = fBWeightVar2->GetBinContent(fBWeightVar2->FindBin(pTmom));
						    fDCAxy_Pt_Bmeson_WeightVar2->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),bWeight);

						    fDCAxy_Pt_Bmeson->Fill(TrkPt,DCA[0]*Bsign*track->Charge());
					    }

				    }if (pid_eleB && (TMath::Abs(pidM)> 5000 || TMath::Abs(pidM)< 5999)) {//if from B baryon
					    if (pTmom>0. && pTmom<50.) { //in proper pt range

						    bWeight = fBWeightNew->GetBinContent(fBWeightNew->FindBin(pTmom));
						    fDCAxy_Pt_Bbaryon_WeightNew->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),bWeight);

						    bWeight = fBWeightVar1->GetBinContent(fBWeightVar1->FindBin(pTmom));
						    fDCAxy_Pt_Bbaryon_WeightVar1->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),bWeight);

						    bWeight = fBWeightVar2->GetBinContent(fBWeightVar2->FindBin(pTmom));
						    fDCAxy_Pt_Bbaryon_WeightVar2->Fill(TrkPt,DCA[0]*Bsign*track->Charge(),bWeight);

						    fDCAxy_Pt_Bbaryon->Fill(TrkPt,DCA[0]*Bsign*track->Charge());
					    }

				    }


				    // 411 : D+, 421 :  D0, 413 : D*+, 423 : D*0, 431 : D_s+, 433 : D_s*+
				    //if(pid_eleD){
				    //				if(TMath::Abs(pidM)==411 || TMath::Abs(pidM)== 413){fDCAxy_Pt_Dpm->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
				    //				if(TMath::Abs(pidM)==421 || TMath::Abs(pidM)== 423){fDCAxy_Pt_D0->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
				    //				if(TMath::Abs(pidM)==431 || TMath::Abs(pidM)== 433){fDCAxy_Pt_Ds->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
				    //}
				    //if(TMath::Abs(pidM)==4122){fDCAxy_Pt_lambda->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
				    //if(pid_eleB){fDCAxy_Pt_B->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
			    }
		    }

		    Double_t TrkPhicos2_hadhigh = -999;
		    Double_t TrkPhisin2_hadhigh = -999;

		    if((fTPCnSigma < -3) && (m20 > 0.01 && m20 < 3.0)){

			    if(track->Pt()<3.0){
				    if( (fTOFnSigma<-1 || fTOFnSigma>1) || (fITSnSigma<-3 || fITSnSigma>1) )continue;
			    }

			    TrkPhiEPV0A_had = TrkPhiPI - PsinV0A;

			    if(TrkPhiEPV0A_had > TMath::Pi()/2.){
				    TrkPhiEPV0A_had -= TMath::Pi();
			    }

			    if(TrkPhiEPV0A_had < -TMath::Pi()/2.){
				    TrkPhiEPV0A_had += TMath::Pi();
			    }

			    TrkPhicos2_hadhigh = cos(2*TrkPhiEPV0A_had);
			    TrkPhisin2_hadhigh = sin(2*TrkPhiEPV0A_had);

			    fTrkPhicos2_hadhigh -> Fill(track->Pt(),TrkPhicos2_hadhigh);
			    fTrkPhisin2_hadhigh -> Fill(track->Pt(),TrkPhisin2_hadhigh);

			    fHisthadron -> Fill(eop,track->Pt());
			    fDCAxy_Pt_had -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());
			    //cout<<"DCA= "<<DCA[i]<<endl;
			    //cout<<"Bsign= "<<Bsign<<endl;
			    //cout<<"charge= "<<track->Charge()<<endl;
			    if(TrkPhiEPV0A_had > -TMath::Pi()/4. && TrkPhiEPV0A_had < TMath::Pi()/4.){
				    fDCAxy_Pt_Inplane_had -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());
			    }
			    if(TrkPhiEPV0A_had > TMath::Pi()/4. || TrkPhiEPV0A_had < -TMath::Pi()/4.){
				    fDCAxy_Pt_Outplane_had -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());
			    }


		    }


		    }

	    }                                                   // continue until all the tracks are processed
	    PostData(1, fOutputList);                           // stream the results the analysis of this event to
	    // it to a file
	    }

//_____________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalRun2::Terminate(Option_t *)
{
	// terminate
	// called at the END of the analysis (when all events are processed)
}
//____________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalRun2::FindMother(AliAODMCParticle* part, int &label, int &pid, double &ptmom){

	if(part->GetMother()>-1){
		label = part->GetMother();
		AliAODMCParticle *partM = (AliAODMCParticle*)fMCarray->At(label);
		pid = partM->GetPdgCode();
		ptmom = partM->Pt();
	}
	else
	{
		pid = -99;
	}
}
//____________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalRun2::GetTrkClsEtaPhiDiff(AliVTrack *t,AliVCluster *v,Double_t &phidiff, Double_t &etadiff)
{
	phidiff = 999;
	etadiff = 999;

	if(!t||!v)return;

	Double_t veta = t->GetTrackEtaOnEMCal();
	Double_t vphi = t->GetTrackPhiOnEMCal();

	Float_t pos[3] = {0};
	v->GetPosition(pos);
	TVector3 cpos(pos);
	Double_t ceta = cpos.Eta();
	Double_t cphi = cpos.Phi();
	etadiff=veta-ceta;
	phidiff=TVector2::Phi_mpi_pi(vphi-cphi);
}
//_____________________________________________________________________________
//void AliAnalysisTaskFlowTPCEMCalRun2::SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Double_t TrkPt, Double_t DCAxy, Int_t Bsign)
void AliAnalysisTaskFlowTPCEMCalRun2::SelectPhotonicElectron(Int_t itrack, AliAODTrack *track, Bool_t &fFlagPhotonicElec, Double_t TrkPt, Double_t DCAxy, Int_t Bsign, Double_t TrkPhiPI, double PsinV0A)
{
	///////////////////////////////////////////
	//////Non-HFE - Invariant mass method//////
	///////////////////////////////////////////


//##################### Set cone radius ######################//

Double_t CutmassMin = massMin;

//############################################################//



	Bool_t flagPhotonicElec = kFALSE;

	Int_t ntracks = -999;
	//ntracks = fVevent->GetNumberOfTracks();
	ntracks = fTracks_tender->GetEntries();

	for (Int_t jtrack = 0; jtrack < ntracks; jtrack++) {
		AliVParticle* VAssotrack = 0x0;
	        // VAssotrack  = fVevent->GetTrack(jtrack);
                 VAssotrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jtrack)); //take tracks from Tender list

        if (!VAssotrack) {
            printf("ERROR: Could not receive track %d\n", jtrack);
            continue;
        }

        AliVTrack *Assotrack = dynamic_cast<AliVTrack*>(VAssotrack);
        AliAODTrack *aAssotrack = dynamic_cast<AliAODTrack*>(VAssotrack);

        //------reject same track
	if(jtrack==itrack) continue;
	if(aAssotrack->Px()==track->Px() && aAssotrack->Py()==track->Py() && aAssotrack->Pz()==track->Pz())continue;

        Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
        Double_t ptAsso=-999., nsigma=-999.0, width = -999., mass=-999.;
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;//electron

        nsigma = fpidResponse->NumberOfSigmasTPC(aAssotrack, AliPID::kElectron);
        ptAsso = aAssotrack->Pt();
        Int_t chargeAsso = aAssotrack->Charge();
        Int_t charge = track->Charge();
        if(charge>0) fPDGe1 = -11;//positron
        if(chargeAsso>0) fPDGe2 = -11;//positron
	if(charge == chargeAsso) fFlagLS = kTRUE;
	if(charge != chargeAsso) fFlagULS = kTRUE;
	//------track cuts applied
	if(fAOD) {
		if(!aAssotrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            if(aAssotrack->GetTPCNcls() < 70) continue;
            if((!(aAssotrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(aAssotrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
        }
        //else{
        //}

        //-------loose cut on partner electron
        //if(ptAsso <0.2) continue;
        if(ptAsso <finvmass_pt) continue;
        if(aAssotrack->Eta()<-0.6 || aAssotrack->Eta()>0.6) continue;
        if(nsigma < -3 || nsigma > 3) continue;

        //-------define KFParticle to get mass
        AliKFParticle::SetField(fVevent->GetMagneticField());
        AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
        AliKFParticle ge2 = AliKFParticle(*aAssotrack, fPDGe2);
        AliKFParticle recg(ge1, ge2);

        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;

	//-------Get mass
	Int_t MassCorrect;
	MassCorrect = recg.GetMass(mass,width);
	Double_t TrkPhiEPV0A_phoLS = -999.;
	Double_t TrkPhiEPV0A_phoULS = -999.;
	Double_t TrkPhicos2_phoLShigh = -999., TrkPhisin2_phoLShigh = -999.;
	Double_t TrkPhicos2_phoULShigh = -999., TrkPhisin2_phoULShigh = -999.;

	if(fFlagLS && track->Pt()>1){

		fInvmassLS->Fill(mass);
		fInvmassLS_2D->Fill(mass,track->Pt());

		//if(mass<CutmassMin){
		if(mass<finvmass){

			fDCAxy_Pt_LS -> Fill(TrkPt,DCAxy*charge*Bsign);
			TrkPhiEPV0A_phoLS = TrkPhiPI - PsinV0A;

			if(TrkPhiEPV0A_phoLS > TMath::Pi()/2.){
				TrkPhiEPV0A_phoLS -= TMath::Pi();
			}

			if(TrkPhiEPV0A_phoLS < -TMath::Pi()/2.){
				TrkPhiEPV0A_phoLS += TMath::Pi();
			}

			TrkPhicos2_phoLShigh = cos(2*TrkPhiEPV0A_phoLS);
			TrkPhisin2_phoLShigh = sin(2*TrkPhiEPV0A_phoLS);

			fTrkPhicos2_phoLShigh -> Fill(track->Pt(),TrkPhicos2_phoLShigh);
			fTrkPhisin2_phoLShigh -> Fill(track->Pt(),TrkPhisin2_phoLShigh);

			if(TrkPhiEPV0A_phoLS > -TMath::Pi()/4. && TrkPhiEPV0A_phoLS < TMath::Pi()/4.){

				fInplane_LSpho -> Fill(track->Pt());
                fDCAxy_Pt_Inplane_LS -> Fill(TrkPt,DCAxy*charge*Bsign);

			}

			if(TrkPhiEPV0A_phoLS > TMath::Pi()/4. || TrkPhiEPV0A_phoLS < -TMath::Pi()/4.){

				fOutplane_LSpho -> Fill(track->Pt());
                fDCAxy_Pt_Outplane_LS -> Fill(TrkPt,DCAxy*charge*Bsign);

			}

		}
	}

	if(fFlagULS && track->Pt()>1){

		fInvmassULS->Fill(mass);
		fInvmassULS_2D->Fill(mass,track->Pt());

		//cout<<"DCAULS= "<<DCAxy<<endl;
		//cout<<"chargeULS= "<<charge<<endl;
		//cout<<"BsignULS= "<<Bsign<<endl;
		

		if(mass<CutmassMin){

			fDCAxy_Pt_ULS -> Fill(TrkPt,DCAxy*charge*Bsign);

			TrkPhiEPV0A_phoULS = TrkPhiPI - PsinV0A;

			//cout << "TrkPhiEPV0A_phoULS = " << TrkPhiEPV0A_phoULS << endl;

			if(TrkPhiEPV0A_phoULS > TMath::Pi()/2.){
				TrkPhiEPV0A_phoULS -= TMath::Pi();
			}

			if(TrkPhiEPV0A_phoULS < -TMath::Pi()/2.){
				TrkPhiEPV0A_phoULS += TMath::Pi();
			}

			//cout << "TrkPhiEPV0A_phoULS =" << TrkPhiEPV0A_phoULS << endl;

			TrkPhicos2_phoULShigh = cos(2*TrkPhiEPV0A_phoULS);
			TrkPhisin2_phoULShigh = sin(2*TrkPhiEPV0A_phoULS);

			//cout << "TrkPhicos2_phoULShigh =" << TrkPhicos2_phoULShigh << endl;

			fTrkPhicos2_phoULShigh -> Fill(track->Pt(),TrkPhicos2_phoULShigh);
			fTrkPhisin2_phoULShigh -> Fill(track->Pt(),TrkPhisin2_phoULShigh);


			if(TrkPhiEPV0A_phoULS > -TMath::Pi()/4. && TrkPhiEPV0A_phoULS < TMath::Pi()/4.){

				fInplane_ULSpho -> Fill(track->Pt());
                fDCAxy_Pt_Inplane_ULS -> Fill(TrkPt,DCAxy*charge*Bsign);

			}

			if(TrkPhiEPV0A_phoULS > TMath::Pi()/4. || TrkPhiEPV0A_phoULS < -TMath::Pi()/4.){

				fOutplane_ULSpho -> Fill(track->Pt());
                fDCAxy_Pt_Outplane_ULS -> Fill(TrkPt,DCAxy*charge*Bsign);

			}


		}
	}

	if(mass<0.1 && fFlagULS && !flagPhotonicElec){
		flagPhotonicElec = kTRUE; //Tag Non-HFE (random mass cut, not optimised)
	}
	fFlagPhotonicElec = flagPhotonicElec;
	}
}

//____________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalRun2::CheckMCgen(AliAODMCHeader* fMCheader,Double_t CutEta)
{
	TList *lh = fMCheader -> GetCocktailHeaders();
	NpureMC = 0;
	NpureMCproc = 0;
	NembMCpi0 = 0;
	NembMCeta = 0;
	Nch = 0;
	TString MCgen;
	TString embpi0("pi");
	TString embeta("eta");
	TString embbeauty("bele");

	if(lh)
{
		//for(int igene=0; igene<lh->GetEntries(); igene++)
		for(int igene=0; igene<lh->GetEntries()-1; igene++)
		{
			AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igene);
			if(gh)
			{
				MCgen = gh->GetName();
                                //cout << "MCgen = " << MCgen << endl;
				if(igene==0)NpureMC = gh->NProduced(); //generated by PYTHIA or HIJING

				if(MCgen.Contains(embpi0))NembMCpi0 = NpureMCproc;
				if(MCgen.Contains(embeta))NembMCeta = NpureMCproc;
				if(MCgen.Contains(embbeauty))
                                   {
                                    iBevt = kTRUE;
                                    }  

				NpureMCproc += gh->NProduced(); //generated by PYTHIA or HIJING
                                //cout << "NpureMCproc = " << NpureMCproc << endl;
			}
		}
}

for(int imc=0; imc<NpureMCproc; imc++)
{
	Bool_t iEnhance = kTRUE;
	//if(imc>=NpureMC)iEnhance = kFALSE;
	Int_t iHijing = 1; //select particles from Hijing or PYTHIA

	fMCparticle = (AliAODMCParticle*) fMCarray->At(imc);
	Int_t pdgGen = TMath::Abs(fMCparticle->GetPdgCode());
	Double_t pdgEta = fMCparticle->Eta();
	Double_t pTtrue = fMCparticle->Pt();
	Int_t chargetrue = fMCparticle->Charge();
	Bool_t isPhysPrim = fMCparticle->IsPhysicalPrimary();

	//------------ Get N charged ------------------
	if(chargetrue!=0){
		if(TMath::Abs(pdgEta)<1.0){
			if(isPhysPrim){

				Nch++;

			}
		}
	}

	if(TMath::Abs(pdgEta)>CutEta)continue;

	fCheckEtaMC->Fill(pdgEta);

	Int_t pdgMom = -99;
	Int_t labelMom = -1;
	Double_t pTmom = -1.0;

	FindMother(fMCparticle,labelMom,pdgMom,pTmom);
	if(pdgMom==-99 && iEnhance)iHijing = 0; //particle from enhance
	if(pdgMom>0 && iEnhance)iHijing = -1; //particles from enhance but feeddown

	if(iHijing>-1)
	{
		if(pdgGen==111)fHistMCorgPi0->Fill(iHijing,pTtrue);
		if(pdgGen==221)fHistMCorgEta->Fill(iHijing,pTtrue);
	}

	if(TMath::Abs(pdgGen)!=11)continue;
	if(pTtrue<2.0)continue;

	Int_t pdgGM = -99;
	Int_t labelGM = -1;
	Double_t pTGMom = -1.0;

	if(pdgMom!=0)
	{
		AliAODMCParticle* fMCparticleMom = (AliAODMCParticle*) fMCarray->At(labelMom);

		if(IsDdecay(pdgMom)){

			fHistMCorgD->Fill(fMCparticle->Pt());

			FindMother(fMCparticleMom,labelGM,pdgGM,pTGMom);

			if(IsBdecay(pdgGM)){

				fPt_Btoe->Fill(fMCparticle->Pt(),pTGMom);

			}
		}

		if(IsBdecay(pdgMom)){

			fHistMCorgB->Fill(fMCparticle->Pt());

			fPt_Btoe->Fill(fMCparticle->Pt(),pTmom);
		}
	}


}
return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalRun2::IsPdecay(int mpid){

int abmpid = TMath::Abs(mpid);
if(abmpid==22 || abmpid==111 || abmpid==221){

	return kTRUE;

}
else
{
	return kFALSE;

}
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalRun2::IsDdecay(int mpid){

	int abmpid = TMath::Abs(mpid);
	if(abmpid==411 || abmpid==421 || abmpid==413 || abmpid==423 || abmpid==431 || abmpid==433 || abmpid==4122){

		return kTRUE;

	}

	else
	{
		return kFALSE;
		}

}

//411:D+, 421:D0, 413:D*+, 423:D*0, 431:D_s+, 433:D_s*

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowTPCEMCalRun2::IsBdecay(int mpid){

int abmpid = TMath::Abs(mpid);
if(abmpid==511 || abmpid==521 || abmpid==513 || abmpid==523 || abmpid==531 || abmpid==533){

	return kTRUE;

}

else
{
return kFALSE;
}

}

//511:B0, 521:B+, 513:B*0, 523:B*+, 531:B_s0, 533:B_s*

//_____________________________________________________________________________
bool AliAnalysisTaskFlowTPCEMCalRun2::LoadSplinesForqnPercentile()
{

	//load splines from file

	TString listname[6] = {"SplineListq2TPC","SplineListq2TPCPosEta","SplineListq2TPCNegEta","SplineListq2V0","SplineListq2V0A","SplineListq2V0C"};

	if(!gGrid){

		TGrid::Connect("alien//");

	}

        cout << "fqnSplineFileName = "<< fqnSplineFileName << endl;
        cout << "fqnSplineFileName = "<< fqnSplineFileName.Data() << endl;
        cout << "fOADBFileName = "<< fOADBFileName << endl;

	TFile* splinesfile = TFile::Open(fqnSplineFileName.Data());
	if(!splinesfile){
		AliFatal("File with splines for qn percentiles not found!");
		return false;
	}

	for(int iDet=0; iDet<6; iDet++){
		fqnSplinesList[iDet] = (TList*)splinesfile->Get(listname[iDet].Data());

		if(!fqnSplinesList[iDet]){
			AliFatal("TList with splines for qn percentiles not found in the spline file!");
			return false;
		}

		fqnSplinesList[iDet]->SetOwner();
	}

	splinesfile->Close();

	return true;

}
//________________________________________________________________________________
double AliAnalysisTaskFlowTPCEMCalRun2::GetDeltaPsiSubInRange(double psi1,double psi2)
{

	double delta = psi1 - psi2;
	if(TMath::Abs(delta) > TMath::Pi() / fHarmonic){
		if(delta>0.)delta -= 2.*TMath::Pi() / fHarmonic;
		else delta += 2.*TMath::Pi() / fHarmonic;
	}

	return delta;

}

//________________________________________________________________________________
void AliAnalysisTaskFlowTPCEMCalRun2::GetMainQnVectorInfo(double &mainPsin,double &mainMultQn,double mainQn[2],double &SubAPsin,double &SubAMultQn,double SubAQn[2],double &SubBPsin,double &SubBMultQn,double SubBQn[2],AliHFQnVectorHandler*HFQnVectorHandler)
{

	double QnFullTPC[2], QnPosTPC[2], QnNegTPC[2];
	double QnFullV0[2], QnV0A[2], QnV0C[2];
	double MultQnFullTPC = -1., MultQnPosTPC = -1., MultQnNegTPC = -1.;
	double MultQnFullV0 = -1., MultQnV0A = -1., MultQnV0C = -1.;
	double PsinFullTPC = -1., PsinPosTPC = -1., PsinNegTPC = -1.;
	double PsinFullV0 = -1., PsinV0A = -1., PsinV0C = -1.;

	//get the unnormalised Qn-vectors --> normalisation can be done in the task
	HFQnVectorHandler->GetUnNormQnVecTPC(QnFullTPC,QnPosTPC,QnNegTPC);
	HFQnVectorHandler->GetUnNormQnVecV0(QnFullV0,QnV0A,QnV0C);

	HFQnVectorHandler->GetMultQnVecTPC(MultQnFullTPC,MultQnPosTPC,MultQnNegTPC);
	HFQnVectorHandler->GetMultQnVecV0(MultQnFullV0,MultQnV0A,MultQnV0C);

	HFQnVectorHandler->GetEventPlaneAngleTPC(PsinFullTPC,PsinPosTPC,PsinNegTPC);
	HFQnVectorHandler->GetEventPlaneAngleV0(PsinFullV0,PsinV0A,PsinV0C);

	if(fEtaGapInTPCHalves>0.) {
		vector<AliAODTrack*> trackstoremove;
		for(int iTrack=0; iTrack<fAOD->GetNumberOfTracks(); iTrack++) {
			AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
			if(TMath::Abs(track->Eta())<fEtaGapInTPCHalves/2)
				trackstoremove.push_back(track);
		}
		HFQnVectorHandler->RemoveTracksFromQnTPC(trackstoremove, QnFullTPC, QnPosTPC, QnNegTPC, MultQnFullTPC, MultQnPosTPC, MultQnNegTPC, true);
		PsinFullTPC = (TMath::Pi()+TMath::ATan2(-QnFullTPC[1],-QnFullTPC[0]))/2;
		PsinPosTPC = (TMath::Pi()+TMath::ATan2(-QnPosTPC[1],-QnPosTPC[0]))/2;
		PsinNegTPC = (TMath::Pi()+TMath::ATan2(-QnNegTPC[1],-QnNegTPC[0]))/2;
	}

	switch(fEvPlaneDet) {
		case kFullTPC:
			mainPsin   = PsinFullTPC;
			mainMultQn = MultQnFullTPC;
			mainQn[0]  = QnFullTPC[0];
			mainQn[1]  = QnFullTPC[1];
			break;
		case kPosTPC:
			mainPsin   = PsinPosTPC;
			mainMultQn = MultQnPosTPC;
			mainQn[0]  = QnPosTPC[0];
			mainQn[1]  = QnPosTPC[1];
			break;
		case kNegTPC:
			mainPsin   = PsinNegTPC;
			mainMultQn = MultQnNegTPC;
			mainQn[0]  = QnNegTPC[0];
			mainQn[1]  = QnNegTPC[1];
			break;
		case kFullV0:
			mainPsin   = PsinFullV0;
			mainMultQn = MultQnFullV0;
			mainQn[0]  = QnFullV0[0];
			mainQn[1]  = QnFullV0[1];
			break;
		case kV0A:
			mainPsin   = PsinV0A;
			mainMultQn = MultQnV0A;
			mainQn[0]  = QnV0A[0];
			mainQn[1]  = QnV0A[1];
			break;
		case kV0C:
			mainPsin   = PsinV0C;
			mainMultQn = MultQnV0C;
			mainQn[0]  = QnV0C[0];
			mainQn[1]  = QnV0C[1];
			break;
	}

	switch(fSubEvDetA) {
		case kFullTPC:
			SubAPsin   = PsinFullTPC;
			SubAMultQn = MultQnFullTPC;
			SubAQn[0]  = QnFullTPC[0];
			SubAQn[1]  = QnFullTPC[1];
			break;
		case kPosTPC:
			SubAPsin   = PsinPosTPC;
			SubAMultQn = MultQnPosTPC;
			SubAQn[0]  = QnPosTPC[0];
			SubAQn[1]  = QnPosTPC[1];
			break;
		case kNegTPC:
			SubAPsin   = PsinNegTPC;
			SubAMultQn = MultQnNegTPC;
			SubAQn[0]  = QnNegTPC[0];
			SubAQn[1]  = QnNegTPC[1];
			break;
		case kFullV0:
			SubAPsin   = PsinFullV0;
			SubAMultQn = MultQnFullV0;
			SubAQn[0]  = QnFullV0[0];
			SubAQn[1]  = QnFullV0[1];
			break;
		case kV0A:
			SubAPsin   = PsinV0A;
			SubAMultQn = MultQnV0A;
			SubAQn[0]  = QnV0A[0];
			SubAQn[1]  = QnV0A[1];
			break;
		case kV0C:
			SubAPsin   = PsinV0C;
			SubAMultQn = MultQnV0C;
			SubAQn[0]  = QnV0C[0];
			SubAQn[1]  = QnV0C[1];
			break;
	}

	switch(fSubEvDetB) {
		case kFullTPC:
			SubBPsin   = PsinFullTPC;
			SubBMultQn = MultQnFullTPC;
			SubBQn[0]  = QnFullTPC[0];
			SubBQn[1]  = QnFullTPC[1];
			break;
		case kPosTPC:
			SubBPsin   = PsinPosTPC;
			SubBMultQn = MultQnPosTPC;
			SubBQn[0]  = QnPosTPC[0];
			SubBQn[1]  = QnPosTPC[1];
			break;
		case kNegTPC:
			SubBPsin   = PsinNegTPC;
			SubBMultQn = MultQnNegTPC;
			SubBQn[0]  = QnNegTPC[0];
			SubBQn[1]  = QnNegTPC[1];
			break;
		case kFullV0:
			SubBPsin   = PsinFullV0;
			SubBMultQn = MultQnFullV0;
			SubBQn[0]  = QnFullV0[0];
			SubBQn[1]  = QnFullV0[1];
			break;
		case kV0A:
			SubBPsin   = PsinV0A;
			SubBMultQn = MultQnV0A;
			SubBQn[0]  = QnV0A[0];
			SubBQn[1]  = QnV0A[1];
			break;
		case kV0C:
			SubBPsin   = PsinV0C;
			SubBMultQn = MultQnV0C;
			SubBQn[0]  = QnV0C[0];
			SubBQn[1]  = QnV0C[1];
			break;
	}


}
