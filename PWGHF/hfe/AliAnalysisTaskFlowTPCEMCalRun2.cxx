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
	fDCAxy_Pt_Inplane_ele(0),
	fDCAxy_Pt_Outplane_ele(0),
	fDCAxy_Pt_Inplane_hfe(0),
	fDCAxy_Pt_Outplane_hfe(0),
	fDCAxy_Inplane_ele(0),
	fDCAxy_Outplane_ele(0),
	fDCAxy_Inplane_hfe(0),
	fDCAxy_Outplane_hfe(0),
	fHistPt_HFE_MC_D(0),
	fHistPt_HFE_MC_B(0),
	fDCAxy_Pt_Dpm(0),
	fDCAxy_Pt_D0(0),
	fDCAxy_Pt_Ds(0),
	fDCAxy_Pt_lambda(0),
	fDCAxy_Pt_B(0),
        ftpcnsig(-1.0),
        femceop(0.9),
        femcss_mim(0.01),
        femcss_max(0.35),
        finvmass(0.1),
        finvmass_pt(0.15),
	massMin(0.1),
	fDCAxy_Pt_LS(0),
	fDCAxy_Pt_ULS(0),
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
	fDCAxy_Pt_Inplane_ele(0),
	fDCAxy_Pt_Outplane_ele(0),
	fDCAxy_Pt_Inplane_hfe(0),
	fDCAxy_Pt_Outplane_hfe(0),
	fDCAxy_Inplane_ele(0),
	fDCAxy_Outplane_ele(0),
	fDCAxy_Inplane_hfe(0),
	fDCAxy_Outplane_hfe(0),
	fHistPt_HFE_MC_D(0),
	fHistPt_HFE_MC_B(0),
	fDCAxy_Pt_Dpm(0),
	fDCAxy_Pt_D0(0),
	fDCAxy_Pt_Ds(0),
	fDCAxy_Pt_lambda(0),
	fDCAxy_Pt_B(0),
        ftpcnsig(-1.0),
        femceop(0.9),
        femcss_mim(0.01),
        femcss_max(0.35),
        finvmass(0.1),
        finvmass_pt(0.15),
	massMin(0.1),
	fDCAxy_Pt_LS(0),
	fDCAxy_Pt_ULS(0),
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

fDCAxy_Pt_ele = new TH2F("fDCAxy_Pt_ele","DCA_{xy} vs Pt (electron);p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_ele);

fDCAxy_Pt_had = new TH2F("fDCAxy_Pt_had","DCA_{xy} vs Pt (hadron);p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_had);

fDCAxy_Pt_LS = new TH2F("fDCAxy_Pt_LS","DCA_{xy} vs Pt LS pairs;p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_LS);

fDCAxy_Pt_ULS = new TH2F("fDCAxy_Pt_ULS","DCA_{xy} vs Pt ULS pairs;p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_ULS);

fDCAxy_Pt_Inplane_ele = new TH2F("fDCAxy_Pt_Inplane_ele","DCA_{xy} vs Pt Inplane electron;p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Inplane_ele);

fDCAxy_Pt_Outplane_ele = new TH2F("fDCAxy_Pt_Outplane_ele","DCA_{xy} vs Pt Outplane electron;p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Outplane_ele);

fDCAxy_Pt_Inplane_hfe = new TH2F("fDCAxy_Pt_Inplane_hfe","DCA_{xy} vs Pt Inplane hfe;p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Inplane_hfe);

fDCAxy_Pt_Outplane_hfe = new TH2F("fDCAxy_Pt_Outplane_hfe","DCA_{xy} vs Pt Outplane hfe;p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Outplane_hfe);

Double_t ptbin[14] = {0.0,0.5,1.0,1.5,2,2.5,3,4,6,8,10,12,14,20};

fDCAxy_Inplane_ele = new TH2F("fDCAxy_Inplane_ele","DCA_{xy} vs Pt Inplane electron;DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Inplane_ele);

fDCAxy_Outplane_ele = new TH2F("fDCAxy_Outplane_ele","DCA_{xy} vs Pt Outplane electron;DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Outplane_ele);

fDCAxy_Inplane_hfe = new TH2F("fDCAxy_Inplane_hfe","DCA_{xy} vs Pt Inplane hfe;DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Inplane_hfe);

fDCAxy_Outplane_hfe = new TH2F("fDCAxy_Outplane_hfe","DCA_{xy} vs Pt Outplane hfe;DCAxy*charge*Bsign",13,ptbin,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Outplane_hfe);

fHistPt_HFE_MC_D  = new TH1F("fHistPt_HFE_MC_D","HFE from D MC",600,0,60);
fOutputList->Add(fHistPt_HFE_MC_D);

fHistPt_HFE_MC_B  = new TH1F("fHistPt_HFE_MC_B","HFE fron B MC",600,0,60);
fOutputList->Add(fHistPt_HFE_MC_B);

fDCAxy_Pt_Dpm = new TH2F("fDCAxy_Pt_Dpm","DCA_{xy} vs Pt D+-(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Dpm);

fDCAxy_Pt_D0= new TH2F("fDCAxy_Pt_D0","DCA_{xy} vs Pt D0(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_D0);

fDCAxy_Pt_Ds= new TH2F("fDCAxy_Pt_Ds","DCA_{xy} vs Pt Ds(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_Ds);

fDCAxy_Pt_lambda = new TH2F("fDCAxy_Pt_lambda","DCA_{xy} vs Pt lambda(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_lambda);

fDCAxy_Pt_B= new TH2F("fDCAxy_Pt_B","DCA_{xy} vs Pt all B meson(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
fOutputList->Add(fDCAxy_Pt_B);

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


    if(fMCarray)CheckMCgen(fMCheader,0.6);

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

		if(pid_eleD || pid_eleB)fNDB->Fill(0);

		if(pid_eleD){

			AliAODMCParticle* fMCTrackpartMom = (AliAODMCParticle*) fMCarray->At(ilabelM);

			FindMother(fMCTrackpartMom,ilabelGM,pidGM,pTGMom);

			if(IsBdecay(pidGM)){

				pid_eleB = IsBdecay(pidGM);

				pid_eleD = kFALSE;

			}

		} //pid_eleD

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

                                }


				if(TrkPhiEPV0A_ele > -TMath::Pi()/4. && TrkPhiEPV0A_ele < TMath::Pi()/4.){

					fInplane_ele -> Fill(track->Pt());
					fDCAxy_Inplane_ele -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());
					fDCAxy_Pt_Inplane_ele -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());

                                         if(!fFlagNonHFE){

                                                 fInplane_hfe -> Fill(track->Pt());//using this at pt < 3 GeV
												 fDCAxy_Inplane_hfe -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());
												 fDCAxy_Pt_Inplane_hfe -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());

                                         }


				}



				if(TrkPhiEPV0A_ele > TMath::Pi()/4. || TrkPhiEPV0A_ele < -TMath::Pi()/4.){

					fOutplane_ele -> Fill(track->Pt());
					fDCAxy_Outplane_ele -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());
					fDCAxy_Pt_Outplane_ele -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());

                                        if(!fFlagNonHFE){

                                                fOutplane_hfe -> Fill(track->Pt());//using this at pt < 3 GeV
												fDCAxy_Outplane_hfe -> Fill(TrkPt,DCA[0]*Bsign*track->Charge());
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

		if(fTPCnSigma > ftpcnsig && fTPCnSigma <3 && m20 > 0.02 && m20 < 0.3){ // TPC nsigma & shower shape cut
						
						if(eop>femceop && eop<1.3){ // E/p cut
								if(pid_eleB) fHistPt_HFE_MC_B -> Fill(track->Pt());
								if(pid_eleD) fHistPt_HFE_MC_D -> Fill(track->Pt());
								
								fDCAxy_Pt_ele->Fill(TrkPt,DCA[0]*Bsign*track->Charge());
										// 411 : D+, 421 :  D0, 413 : D*+, 423 : D*0, 431 : D_s+, 433 : D_s*+
										if(pid_eleD){
														if(TMath::Abs(pidM)==411 || TMath::Abs(pidM)== 413){fDCAxy_Pt_Dpm->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
														if(TMath::Abs(pidM)==421 || TMath::Abs(pidM)== 423){fDCAxy_Pt_D0->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
														if(TMath::Abs(pidM)==431 || TMath::Abs(pidM)== 433){fDCAxy_Pt_Ds->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
										}
										if(TMath::Abs(pidM)==4122){fDCAxy_Pt_lambda->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
										if(pid_eleB){fDCAxy_Pt_B->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
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

			}

			if(TrkPhiEPV0A_phoLS > TMath::Pi()/4. || TrkPhiEPV0A_phoLS < -TMath::Pi()/4.){

				fOutplane_LSpho -> Fill(track->Pt());

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

			}

			if(TrkPhiEPV0A_phoULS > TMath::Pi()/4. || TrkPhiEPV0A_phoULS < -TMath::Pi()/4.){

				fOutplane_ULSpho -> Fill(track->Pt());

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
	if(abmpid==411 || abmpid==421 || abmpid==413 || abmpid==423 || abmpid==431 || abmpid==433){

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
