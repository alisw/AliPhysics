/**********************************************************************************
 * Date   :: 08/02/2018								  *
 * Author :: Pranjal Sarma & B. Bhattacharjee, Gauhati University, India	  *
 * MC Analysis Task for Pi,Ka & Pr vs Multiplicity in pp @ 13 TeV with TOF	  *
 **********************************************************************************/

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "AliPIDResponse.h"
#include <AliPID.h>
#include <TSpline.h>
#include "AliESDtrack.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliVTrack.h"
#include "AliAnalysisTaskTOFMC.h"
#include "AliVEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliVEvent.h"
#include "AliMCParticle.h"

#include "AliCentrality.h"
#include "TProfile.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"

#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliMultSelectionTask.h"


#include "TParticle.h"
#include "TClonesArray.h"
#include "AliTrackReference.h"
#include "iostream"
// Authors: Pranjal Sarma(Date modified=29/08/16)

Double_t v0mpc,v0mpc_inel,dtPion,dtKaon,dtProton;

ClassImp(AliAnalysisTaskTOFMC)

//________________________________________________________________________
AliAnalysisTaskTOFMC::AliAnalysisTaskTOFMC()
 : AliAnalysisTaskSE(),fESD(0), fOutputList(0),fPIDResponse(0),fesdTrackCuts(0x0),fesdTrackCuts_no_dca(0x0),
fTrigSel(AliVEvent::kINT7),fMultSelection(0x0),fMultSelection_INEL(0x0),
fdEdxP(0),fdEdxPt(0),fdEdxPq(0),fdEdxPtq(0),fbetaAllPt(0),fbetaAllP(0),fbetaAllPtq(0),fbetaAllPq(0),
fPtVsTPion(0),fPtVsTKaon(0),fPtVsTProton(0),
fEventCounter(0),fEventPS(0),fEventVtx(0),fEventVtx10(0),fZVertex(0),fZVertexRec(0),fZVertexEff(0),fZVertex10(0),
fCorrRefMultVsNch(0),fCorrRefMultVsV0M(0),fV0MPC(0),

fPtV0MGenPion(0),fPtV0MTPCRecPion(0),fPtV0MTOFRecPion(0),
fPtV0MGenPionP(0),fPtV0MTPCRecPionP(0),fPtV0MTOFRecPionP(0),
fPtV0MGenPionM(0),fPtV0MTPCRecPionM(0),fPtV0MTOFRecPionM(0),

fPtV0MGenKaon(0),fPtV0MTPCRecKaon(0),fPtV0MTOFRecKaon(0),
fPtV0MGenKaonP(0),fPtV0MTPCRecKaonP(0),fPtV0MTOFRecKaonP(0),
fPtV0MGenKaonM(0),fPtV0MTPCRecKaonM(0),fPtV0MTOFRecKaonM(0),

fPtV0MGenProton(0),fPtV0MTPCRecProton(0),fPtV0MTOFRecProton(0),
fPtV0MGenProtonP(0),fPtV0MTPCRecProtonP(0),fPtV0MTOFRecProtonP(0),
fPtV0MGenProtonM(0),fPtV0MTPCRecProtonM(0),fPtV0MTOFRecProtonM(0),

fPtV0MDCAxyTOFPriPion(0),fPtV0MDCAxyTOFPriPionP(0),fPtV0MDCAxyTOFPriPionM(0),
fPtV0MDCAxyTOFWeakPion(0),fPtV0MDCAxyTOFWeakPionP(0),fPtV0MDCAxyTOFWeakPionM(0),
fPtV0MDCAxyTOFMatPion(0),fPtV0MDCAxyTOFMatPionP(0),fPtV0MDCAxyTOFMatPionM(0),

fPtV0MDCAxyTOFPriProton(0),fPtV0MDCAxyTOFPriProtonP(0),fPtV0MDCAxyTOFPriProtonM(0),
fPtV0MDCAxyTOFWeakProton(0),fPtV0MDCAxyTOFWeakProtonP(0),fPtV0MDCAxyTOFWeakProtonM(0),
fPtV0MDCAxyTOFMatProton(0),fPtV0MDCAxyTOFMatProtonP(0),fPtV0MDCAxyTOFMatProtonM(0),


fPPVsMultUtils(new AliPPVsMultUtils()),
fPtV0MGenPion_kINT7(0),fPtV0MGenPion_inel(0),fPtV0MGenPion_signal_loss(0),
fPtV0MGenKaon_kINT7(0),fPtV0MGenKaon_inel(0),fPtV0MGenKaon_signal_loss(0),
fPtV0MGenProton_kINT7(0),fPtV0MGenProton_inel(0),fPtV0MGenProton_signal_loss(0),


fPtV0MTOFRecPion_nSigma(0),fPtV0MTOFRecPionP_nSigma(0),fPtV0MTOFRecPionM_nSigma(0),
fPtV0MTOFRecKaon_nSigma(0),fPtV0MTOFRecKaonP_nSigma(0),fPtV0MTOFRecKaonM_nSigma(0),
fPtV0MTOFRecProton_nSigma(0),fPtV0MTOFRecProtonP_nSigma(0),fPtV0MTOFRecProtonM_nSigma(0),

fPtV0MTOFRecPion_nSigma_excl(0),fPtV0MTOFRecPionP_nSigma_excl(0),fPtV0MTOFRecPionM_nSigma_excl(0),
fPtV0MTOFRecKaon_nSigma_excl(0),fPtV0MTOFRecKaonP_nSigma_excl(0),fPtV0MTOFRecKaonM_nSigma_excl(0),
fPtV0MTOFRecProton_nSigma_excl(0),fPtV0MTOFRecProtonP_nSigma_excl(0),fPtV0MTOFRecProtonM_nSigma_excl(0),

fTOFTimeV0MPtPi(0),fTOFTimeV0MPtK(0),fTOFTimeV0MPtP(0),
fTOFTimeV0MPtMismatchDecayPi(0),fTOFTimeV0MPtMismatchDecayK(0),fTOFTimeV0MPtMismatchDecayP(0), 

fEventV0MPS(0),fEventV0MVtx(0),fEventV0M(0),
fTOFLabel(),

fV0MPC_vertexcut(0),
fPtTPC_AllP(0),fPtTOF_AllP(0), fPtTPC_AllN(0),fPtTOF_AllN(0),

fTPC_CR(0), fChi2TPCcluster(0), fDCAZ(0),fDCAxy(0),

fMinTPCcr(0),fMaxChi2PerTPC(0),fMaxDCAz(0),fMaxDCAxy(0)

{}

//________________________________________________________________________
//AliAnalysisTaskTOFMC::AliAnalysisTaskTOFMC(const char *name)
AliAnalysisTaskTOFMC::AliAnalysisTaskTOFMC(const char *PeriodName, Int_t nTPC_CR, Int_t Chi2_TPCcluser, Int_t DCAz, Int_t DCAxy)
     : AliAnalysisTaskSE("name"),fESD(0), fOutputList(0),fPIDResponse(0),fesdTrackCuts(0x0),fesdTrackCuts_no_dca(0x0),
fTrigSel(AliVEvent::kINT7),fMultSelection(0x0),fMultSelection_INEL(0x0),
fdEdxPt(0),fdEdxP(0),fdEdxPtq(0),fdEdxPq(0),fbetaAllPt(0),fbetaAllP(0),fbetaAllPtq(0),fbetaAllPq(0),
fPtVsTPion(0),fPtVsTKaon(0),fPtVsTProton(0),
fEventCounter(0),fEventPS(0),fEventVtx(0),fEventVtx10(0),fZVertex(0),fZVertexRec(0),fZVertexEff(0),fZVertex10(0),


fCorrRefMultVsNch(0),fCorrRefMultVsV0M(0),fV0MPC(0),
fPtV0MGenPion(0),fPtV0MTPCRecPion(0),fPtV0MTOFRecPion(0),
fPtV0MGenPionP(0),fPtV0MTPCRecPionP(0),fPtV0MTOFRecPionP(0),
fPtV0MGenPionM(0),fPtV0MTPCRecPionM(0),fPtV0MTOFRecPionM(0),

fPtV0MGenKaon(0),fPtV0MTPCRecKaon(0),fPtV0MTOFRecKaon(0),
fPtV0MGenKaonP(0),fPtV0MTPCRecKaonP(0),fPtV0MTOFRecKaonP(0),
fPtV0MGenKaonM(0),fPtV0MTPCRecKaonM(0),fPtV0MTOFRecKaonM(0),

fPtV0MGenProton(0),fPtV0MTPCRecProton(0),fPtV0MTOFRecProton(0),
fPtV0MGenProtonP(0),fPtV0MTPCRecProtonP(0),fPtV0MTOFRecProtonP(0),
fPtV0MGenProtonM(0),fPtV0MTPCRecProtonM(0),fPtV0MTOFRecProtonM(0),

fPtV0MDCAxyTOFPriPion(0),fPtV0MDCAxyTOFPriPionP(0),fPtV0MDCAxyTOFPriPionM(0),
fPtV0MDCAxyTOFWeakPion(0),fPtV0MDCAxyTOFWeakPionP(0),fPtV0MDCAxyTOFWeakPionM(0),
fPtV0MDCAxyTOFMatPion(0),fPtV0MDCAxyTOFMatPionP(0),fPtV0MDCAxyTOFMatPionM(0),

fPtV0MDCAxyTOFPriProton(0),fPtV0MDCAxyTOFPriProtonP(0),fPtV0MDCAxyTOFPriProtonM(0),
fPtV0MDCAxyTOFWeakProton(0),fPtV0MDCAxyTOFWeakProtonP(0),fPtV0MDCAxyTOFWeakProtonM(0),
fPtV0MDCAxyTOFMatProton(0),fPtV0MDCAxyTOFMatProtonP(0),fPtV0MDCAxyTOFMatProtonM(0),


fPPVsMultUtils(new AliPPVsMultUtils()),
fPtV0MGenPion_kINT7(0),fPtV0MGenPion_inel(0),fPtV0MGenPion_signal_loss(0),
fPtV0MGenKaon_kINT7(0),fPtV0MGenKaon_inel(0),fPtV0MGenKaon_signal_loss(0),
fPtV0MGenProton_kINT7(0),fPtV0MGenProton_inel(0),fPtV0MGenProton_signal_loss(0),

fPtV0MTOFRecPion_nSigma(0),fPtV0MTOFRecPionP_nSigma(0),fPtV0MTOFRecPionM_nSigma(0),
fPtV0MTOFRecKaon_nSigma(0),fPtV0MTOFRecKaonP_nSigma(0),fPtV0MTOFRecKaonM_nSigma(0),
fPtV0MTOFRecProton_nSigma(0),fPtV0MTOFRecProtonP_nSigma(0),fPtV0MTOFRecProtonM_nSigma(0),

fPtV0MTOFRecPion_nSigma_excl(0),fPtV0MTOFRecPionP_nSigma_excl(0),fPtV0MTOFRecPionM_nSigma_excl(0),
fPtV0MTOFRecKaon_nSigma_excl(0),fPtV0MTOFRecKaonP_nSigma_excl(0),fPtV0MTOFRecKaonM_nSigma_excl(0),
fPtV0MTOFRecProton_nSigma_excl(0),fPtV0MTOFRecProtonP_nSigma_excl(0),fPtV0MTOFRecProtonM_nSigma_excl(0),

fTOFTimeV0MPtPi(0),fTOFTimeV0MPtK(0),fTOFTimeV0MPtP(0),
fTOFTimeV0MPtMismatchDecayPi(0),fTOFTimeV0MPtMismatchDecayK(0),fTOFTimeV0MPtMismatchDecayP(0),

fEventV0MPS(0),fEventV0MVtx(0),fEventV0M(0),
fTOFLabel(),

fV0MPC_vertexcut(0),
fPtTPC_AllP(0),fPtTOF_AllP(0), fPtTPC_AllN(0),fPtTOF_AllN(0),

fTPC_CR(0), fChi2TPCcluster(0), fDCAZ(0),fDCAxy(0),

fMinTPCcr(0),fMaxChi2PerTPC(0),fMaxDCAz(0),fMaxDCAxy(0)



{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 Recoiner
  DefineOutput(1, TList::Class());

	for (Int_t i = 0; i < 3; i++) fTOFLabel[i] = 999;
}

//________________________________________________________________________
void AliAnalysisTaskTOFMC::UserCreateOutputObjects()
{


  //pid response object
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

inputHandler->SetNeedField();

	 Double_t nPtbinsAll=39;
        Double_t PtbinsAll[] = {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5};

	Double_t nPtbins=59;
	const Double_t Ptbins[] = {0.01, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0};//nicolo binning


        const Int_t nRefMbins = 120;
        const Double_t RefMbins[121] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120};

	 Double_t nV0Mbins=11;
        Double_t V0Mbins[12] ={0.0, 0.1, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};//for next run i.e 10aug
        //Double_t V0Mbins[11] ={0.0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};//9aug


	double DCAxybins[] = {-3, -2.95, -2.9, -2.85, -2.8, -2.75, -2.7, -2.65, -2.6, -2.55, -2.5, -2.45, -2.4, -2.35, -2.3, -2.25, -2.2, -2.15, -2.1, -2.05, -2, -1.95, -1.9, -1.85, -1.8, -1.75, -1.7, -1.65, -1.6, -1.55, -1.5, -1.45, -1.4, -1.35, -1.3, -1.25, -1.2, -1.15, -1.1, -1.05, -1, -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0 , 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3};



	 Int_t nTimebins=500;
        Double_t Timebins[501];

  // Create histograms
  // Called once
  fOutputList = new TList();

                                                                                
	fEventCounter = new TH1F( "fEventCounter", ";Evt. Sel. Step;Count",9,0,9);
	fEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fEventCounter->GetXaxis()->SetBinLabel(2, "Pass Physics selection and trigger");
        fEventCounter->GetXaxis()->SetBinLabel(3, "INEL>0");
        fEventCounter->GetXaxis()->SetBinLabel(4, "Incomplete DAQ");
        fEventCounter->GetXaxis()->SetBinLabel(5, "SPD cluster vs tracklets cut");
        fEventCounter->GetXaxis()->SetBinLabel(6, "PileupSPDInMultBins");
        fEventCounter->GetXaxis()->SetBinLabel(7, "Reconstructed Vertex");
        fEventCounter->GetXaxis()->SetBinLabel(8, "Vertexcut<10");
        fEventCounter->GetXaxis()->SetBinLabel(9, "Selected by Analysis");
	fOutputList->Add(fEventCounter);

	fEventPS = new TH1F("fEventPS","Event after Physics selection",3,0,3);
	fOutputList->Add(fEventPS);

        fEventVtx = new TH1F("fEventVtx","Event after Phy sel & w/o vertex cut",3,0,3);
        fOutputList->Add(fEventVtx);

        fEventVtx10 = new TH1F("fEventVtx10","Event after Phy sel & w vertex cut",3,0,3);
        fOutputList->Add(fEventVtx10);

        fZVertex = new TH1F("fZVertex","Z vertex dist;Vtx_{z};Counts",40,-20,20);
        fOutputList->Add(fZVertex);

        fZVertexRec = new TH1F("fZVertexRec","Z vertex dist;Vtx_{z};Counts",40,-20,20);
        fOutputList->Add(fZVertexRec);

        fZVertexEff = new TH1F("fZVertexEff","Z vertex dist;Vtx_{z};Counts",40,-20,20);
        fOutputList->Add(fZVertexEff);

        fZVertex10 = new TH1F("fZVertex10","Z vertex dist;Vtx_{z};Counts",40,-20,20);
        fOutputList->Add(fZVertex10);

	fV0MPC = new TH1F("fV0MPC","V0M PC;V0M PC;count",120,0,120);
	//fV0MPC = new TH1F("fV0MPC","V0M PC;V0M PC;count",nV0Mbins,V0Mbins);
	fOutputList->Add(fV0MPC);

fCorrRefMultVsNch = new TH2F("fCorrRefMultVsNch","Ref Mult vs N^{Gen}_{ch};N^{Gen}_{ch};Ref Mult_{|#eta|<0.8}",120,0,120,nRefMbins,RefMbins);
        fOutputList->Add(fCorrRefMultVsNch);
fCorrRefMultVsV0M = new TH2F("fCorrRefMultVsV0M","Ref Mult vs V0M PC;V0M PC;Ref Mult_{|#eta|<0.8}",nV0Mbins,V0Mbins,nRefMbins,RefMbins);
        fOutputList->Add(fCorrRefMultVsV0M);

	fPtV0MGenPion = new TH2F("fPtV0MGenPion","p_{T} vs V0M of generated;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenPion);

	fPtV0MTPCRecPion = new TH2F("fPtV0MTPCRecPion","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTPCRecPion);

	fPtV0MTOFRecPion = new TH2F("fPtV0MTOFRecPion","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecPion);

	fPtV0MGenPionP = new TH2F("fPtV0MGenPionP","p_{T} vs V0M  of generated;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenPionP);

	fPtV0MTPCRecPionP = new TH2F("fPtV0MTPCRecPionP","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTPCRecPionP);

	fPtV0MTOFRecPionP = new TH2F("fPtV0MTOFRecPionP","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecPionP);

	fPtV0MGenPionM = new TH2F("fPtV0MGenPionM","p_{T} vs V0M  of generated;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenPionM);

	fPtV0MTPCRecPionM = new TH2F("fPtV0MTPCRecPionM","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
	fOutputList->Add(fPtV0MTPCRecPionM);

	fPtV0MTOFRecPionM = new TH2F("fPtV0MTOFRecPionM","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecPionM);
        


	fPtV0MGenKaon = new TH2F("fPtV0MGenKaon","p_{T} vs V0M  of generated;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
	fOutputList->Add(fPtV0MGenKaon);
        fPtV0MTPCRecKaon = new TH2F("fPtV0MTPCRecKaon","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTPCRecKaon);
        fPtV0MTOFRecKaon = new TH2F("fPtV0MTOFRecKaon","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecKaon);
        fPtV0MGenKaonP = new TH2F("fPtV0MGenKaonP","p_{T} vs V0M  of generated;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenKaonP);
        fPtV0MTPCRecKaonP = new TH2F("fPtV0MTPCRecKaonP","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTPCRecKaonP);
        fPtV0MTOFRecKaonP = new TH2F("fPtV0MTOFRecKaonP","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecKaonP);
        fPtV0MGenKaonM = new TH2F("fPtV0MGenKaonM","p_{T} vs V0M  of generated;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenKaonM);
        fPtV0MTPCRecKaonM = new TH2F("fPtV0MTPCRecKaonM","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTPCRecKaonM);
        fPtV0MTOFRecKaonM = new TH2F("fPtV0MTOFRecKaonM","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecKaonM);



	fPtV0MGenProton = new TH2F("fPtV0MGenProton","p_{T} vs V0M  of generated;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
	fOutputList->Add(fPtV0MGenProton);
        fPtV0MTPCRecProton = new TH2F("fPtV0MTPCRecProton","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTPCRecProton);
        fPtV0MTOFRecProton = new TH2F("fPtV0MTOFRecProton","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecProton);
        fPtV0MGenProtonP = new TH2F("fPtV0MGenProtonP","p_{T} vs V0M  of generated;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenProtonP);
        fPtV0MTPCRecProtonP = new TH2F("fPtV0MTPCRecProtonP","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTPCRecProtonP);
        fPtV0MTOFRecProtonP = new TH2F("fPtV0MTOFRecProtonP","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecProtonP);
        fPtV0MGenProtonM = new TH2F("fPtV0MGenProtonM","p_{T} vs V0M of generated;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenProtonM);
        fPtV0MTPCRecProtonM = new TH2F("fPtV0MTPCRecProtonM","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTPCRecProtonM);
        fPtV0MTOFRecProtonM = new TH2F("fPtV0MTOFRecProtonM","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecProtonM);


	Double_t DCAbin[1201];
	DCAbin[0]=-3.00;
        for (Int_t i=1;i<1201;i++) DCAbin[i]=DCAbin[i-1]+0.005;

	fPtV0MDCAxyTOFPriPion = new TH3F("fPtV0MDCAxyTOFPriPion","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
       fOutputList->Add(fPtV0MDCAxyTOFPriPion);
	fPtV0MDCAxyTOFWeakPion = new TH3F("fPtV0MDCAxyTOFWeakPion","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFWeakPion);
	fPtV0MDCAxyTOFMatPion = new TH3F("fPtV0MDCAxyTOFMatPion","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFMatPion);
	fPtV0MDCAxyTOFPriPionP = new TH3F("fPtV0MDCAxyTOFPriPionP","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFPriPionP);
	fPtV0MDCAxyTOFWeakPionP = new TH3F("fPtV0MDCAxyTOFWeakPionP","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFWeakPionP);
	fPtV0MDCAxyTOFMatPionP = new TH3F("fPtV0MDCAxyTOFMatPionP","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFMatPionP);
	fPtV0MDCAxyTOFPriPionM = new TH3F("fPtV0MDCAxyTOFPriPionM","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFPriPionM);
	fPtV0MDCAxyTOFWeakPionM = new TH3F("fPtV0MDCAxyTOFWeakPionM","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFWeakPionM);
	fPtV0MDCAxyTOFMatPionM = new TH3F("fPtV0MDCAxyTOFMatPionM","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFMatPionM);


	fPtV0MDCAxyTOFPriProton = new TH3F("fPtV0MDCAxyTOFPriProton","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFPriProton);

	fPtV0MDCAxyTOFWeakProton = new TH3F("fPtV0MDCAxyTOFWeakProton","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFWeakProton);
	fPtV0MDCAxyTOFMatProton = new TH3F("fPtV0MDCAxyTOFMatProton","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFMatProton);
	fPtV0MDCAxyTOFPriProtonP = new TH3F("fPtV0MDCAxyTOFPriProtonP","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFPriProtonP);
	fPtV0MDCAxyTOFWeakProtonP = new TH3F("fPtV0MDCAxyTOFWeakProtonP","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFWeakProtonP);
	fPtV0MDCAxyTOFMatProtonP = new TH3F("fPtV0MDCAxyTOFMatProtonP","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFMatProtonP);
	fPtV0MDCAxyTOFPriProtonM = new TH3F("fPtV0MDCAxyTOFPriProtonM","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFPriProtonM);
	fPtV0MDCAxyTOFWeakProtonM = new TH3F("fPtV0MDCAxyTOFWeakProtonM","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFWeakProtonM);
	fPtV0MDCAxyTOFMatProtonM = new TH3F("fPtV0MDCAxyTOFMatProtonM","Pt vs V0M vs DCAxy;p_{T} (GeV/c);V0M PC;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fPtV0MDCAxyTOFMatProtonM);


	fdEdxP = new TH2F("fdEdxP","TPC dE/dx vs p_{TPC};p_{TPC};dE/dx", 500,0,5,500,0,500);
	fOutputList->Add(fdEdxP);
        fdEdxPt = new TH2F("fdEdxPt","TPC dE/dx vs p_{T};p_{T};dE/dx", 500,0,5,500,0,500);
        fOutputList->Add(fdEdxPt);
        fdEdxPq = new TH2F("fdEdxPq","TPC dE/dx vs p_{TPC}/q;p_{TPC}/q;dE/dx", 1000,-5,5,500,0,500);
        fOutputList->Add(fdEdxPq);
        fdEdxPtq = new TH2F("fdEdxPtq","TPC dE/dx vs p_{T}/q;p_{T}/q;dE/dx", 1000,-5,5,500,0,500);
        fOutputList->Add(fdEdxPtq);
	fbetaAllPt = new TH2F("fbetaAllPt","#beta distribution of all particles;p_{T};#beta",500,0,5.,1100,0.,1.1);
	fOutputList->Add(fbetaAllPt);
        fbetaAllP = new TH2F("fbetaAllP","#beta distribution of all particles;p;#beta",500,0,5.,1100,0.,1.1);
        fOutputList->Add(fbetaAllP);
	fbetaAllPtq = new TH2F("fbetaAllPtq","#beta distribution of all particles;p_{T}./q;#beta",1000,-5,5.,1100,0.,1.1);
	fOutputList->Add(fbetaAllPtq);
        fbetaAllPq = new TH2F("fbetaAllPq","#beta distribution of all particles;p/q;#beta",1000,-5,5.,1100,0.,1.1);
        fOutputList->Add(fbetaAllPq);
	fPtVsTPion = new TH2F("fPtVsTPion",";#Delta t_{pi};p_{T}",nPtbins,Ptbins,250,-1000,6000);
	fOutputList->Add(fPtVsTPion);
        fPtVsTKaon = new TH2F("fPtVsTKaon",";#Delta t_{K};p_{T}",nPtbins,Ptbins,500,-10000,10000);
        fOutputList->Add(fPtVsTKaon);
        fPtVsTProton = new TH2F("fPtVsTProton",";#Delta t_{p};p_{T}",nPtbins,Ptbins,500,-18000,2000);
        fOutputList->Add(fPtVsTProton);


	fPtV0MGenPion_kINT7 = new TH2F("fPtV0MGenPion_kINT7","Pt vs vom",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenPion_kINT7);
	fPtV0MGenPion_inel = new TH2F("fPtV0MGenPion_inel","Pt vs vom",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenPion_inel);
	fPtV0MGenPion_signal_loss = new TH2F("fPtV0MGenPion_signal_loss","Pt vs vom",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenPion_signal_loss);

	fPtV0MGenKaon_kINT7 = new TH2F("fPtV0MGenKaon_kINT7","Pt vs vom",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenKaon_kINT7);
	fPtV0MGenKaon_inel = new TH2F("fPtV0MGenKaon_inel","Pt vs vom",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenKaon_inel);
	fPtV0MGenKaon_signal_loss = new TH2F("fPtV0MGenKaon_signal_loss","Pt vs vom",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenKaon_signal_loss);

	fPtV0MGenProton_kINT7 = new TH2F("fPtV0MGenProton_kINT7","Pt vs vom",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenProton_kINT7);
	fPtV0MGenProton_inel = new TH2F("fPtV0MGenProton_inel","Pt vs vom",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenProton_inel);
	fPtV0MGenProton_signal_loss = new TH2F("fPtV0MGenProton_signal_loss","Pt vs v0m",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MGenProton_signal_loss);

	
fPtV0MTOFRecPion_nSigma = new TH2F("fPtV0MTOFRecPion_nSigma","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecPion_nSigma);
fPtV0MTOFRecPionP_nSigma = new TH2F("fPtV0MTOFRecPionP_nSigma","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecPionP_nSigma);
fPtV0MTOFRecPionM_nSigma = new TH2F("fPtV0MTOFRecPionM_nSigma","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecPionM_nSigma);
fPtV0MTOFRecKaon_nSigma = new TH2F("fPtV0MTOFRecKaon_nSigma","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecKaon_nSigma);
fPtV0MTOFRecKaonP_nSigma = new TH2F("fPtV0MTOFRecKaonP_nSigma","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecKaonP_nSigma);
fPtV0MTOFRecKaonM_nSigma = new TH2F("fPtV0MTOFRecKaonM_nSigma","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecKaonM_nSigma);
fPtV0MTOFRecProton_nSigma = new TH2F("fPtV0MTOFRecProton_nSigma","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecProton_nSigma);
fPtV0MTOFRecProtonP_nSigma = new TH2F("fPtV0MTOFRecProtonP_nSigma","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecProtonP_nSigma);
fPtV0MTOFRecProtonM_nSigma = new TH2F("fPtV0MTOFRecProtonM_nSigma","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecProtonM_nSigma);


	fPtV0MTOFRecPion_nSigma_excl = new TH2F("fPtV0MTOFRecPion_nSigma_excl","p_{T} vs V0M  of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecPion_nSigma_excl);
	fPtV0MTOFRecPionP_nSigma_excl = new TH2F("fPtV0MTOFRecPionP_nSigma_excl","p_{T} vs V0M  of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecPionP_nSigma_excl);
	fPtV0MTOFRecPionM_nSigma_excl = new TH2F("fPtV0MTOFRecPionM_nSigma_excl","p_{T} vs V0M  of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecPionM_nSigma_excl);
	fPtV0MTOFRecKaon_nSigma_excl = new TH2F("fPtV0MTOFRecKaon_nSigma_excl","p_{T} vs V0M  of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecKaon_nSigma_excl);
	fPtV0MTOFRecKaonP_nSigma_excl = new TH2F("fPtV0MTOFRecKaonP_nSigma_excl","p_{T} vs V0M  of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecKaonP_nSigma_excl);
	fPtV0MTOFRecKaonM_nSigma_excl = new TH2F("fPtV0MTOFRecKaonM_nSigma_excl","p_{T} vs V0M  of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecKaonM_nSigma_excl);
	fPtV0MTOFRecProton_nSigma_excl = new TH2F("fPtV0MTOFRecProton_nSigma_excl","p_{T} vs V0M  of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecProton_nSigma_excl);
	fPtV0MTOFRecProtonP_nSigma_excl = new TH2F("fPtV0MTOFRecProtonP_nSigma_excl","p_{T} vs V0M  of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecProtonP_nSigma_excl);
	fPtV0MTOFRecProtonM_nSigma_excl = new TH2F("fPtV0MTOFRecProtonM_nSigma_excl","p_{T} vs V0M of reco;p_{T} (GeV/c);V0M PC",nPtbins,Ptbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fPtV0MTOFRecProtonM_nSigma_excl);


	 Timebins[0]=-10000.0;
        for (Int_t i=1;i<501;i++){
        Timebins[i]=Timebins[i-1]+40.0;
}	

	fTOFTimeV0MPtPi=new TH3F("fTOFTimeV0MPtPi","TOF Time vs pT #pi;p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFTimeV0MPtPi);
	fTOFTimeV0MPtK=new TH3F("fTOFTimeV0MPtK","TOF Time vs pT K;p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
        fOutputList->Add(fTOFTimeV0MPtK);
	fTOFTimeV0MPtP=new TH3F("fTOFTimeV0MPtP","TOF Time vs pT P;p_{T} (GeV/c);T-T_{0}-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
        fOutputList->Add(fTOFTimeV0MPtP);


	//mismatch
	fTOFTimeV0MPtMismatchDecayPi=new TH3F("fTOFTimeV0MPtMismatchDecayPi"," Mis decay TOF Time vs pT #pi;p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFTimeV0MPtMismatchDecayPi);
	fTOFTimeV0MPtMismatchDecayK=new TH3F("fTOFTimeV0MPtMismatchDecayK","Mis decay TOF Time vs pT K;p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
        fOutputList->Add(fTOFTimeV0MPtMismatchDecayK);
	fTOFTimeV0MPtMismatchDecayP=new TH3F("fTOFTimeV0MPtMismatchDecayP","Mis decay TOF Time vs pT P;p_{T} (GeV/c);T-T_{0}-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
        fOutputList->Add(fTOFTimeV0MPtMismatchDecayP);

	Double_t eventbins[]={0,1,2,3};
        fEventV0MPS = new TH2F("fEventV0MPS","Event vs V0M after PS ;Events;V0M percentile",3,eventbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fEventV0MPS);
        fEventV0MVtx = new TH2F("fEventV0MVtx","Event vs V0M after reco vertex;Events;V0M percentile",3,eventbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fEventV0MVtx);
        fEventV0M = new TH2F("fEventV0M","Event vs V0M after reco vertex |Z|<10;Events;V0M percentile",3,eventbins,nV0Mbins,V0Mbins);
        fOutputList->Add(fEventV0M);


        fV0MPC_vertexcut=new TH1F("fV0MPC_vertexcut","V0M percentile;V0M PC",120,0,120);
        fOutputList->Add(fV0MPC_vertexcut);

        fPtTPC_AllP=new TH1F("fPtTPC_AllP","TPC pt distribution;p_{T} (GeV/c)",nPtbins,Ptbins);
        fOutputList->Add(fPtTPC_AllP);
        fPtTOF_AllP=new TH1F("fPtTOF_AllP","TOF pt distribution;p_{T} (GeV/c)",nPtbins,Ptbins);
        fOutputList->Add(fPtTOF_AllP);
        fPtTPC_AllN=new TH1F("fPtTPC_AllN","TPC pt distribution;p_{T} (GeV/c)",nPtbins,Ptbins);
        fOutputList->Add(fPtTPC_AllN);
        fPtTOF_AllN=new TH1F("fPtTOF_AllN","TOF pt distribution;p_{T} (GeV/c)",nPtbins,Ptbins);
        fOutputList->Add(fPtTOF_AllN);

	fTPC_CR=new TH1F("fTPC_CR","TPC cr distribution;# of CR",200,0,200);
        fOutputList->Add(fTPC_CR);
        fChi2TPCcluster=new TH1F("fChi2TPCcluster","Chi2 /TPC cluster distribution;# of CR",100,0,10);
        fOutputList->Add(fChi2TPCcluster);
        fDCAZ=new TH1F("fDCAZ","DCA z distribution;# of CR",100,0,10);
        fOutputList->Add(fDCAZ);
        fDCAxy=new TH1F("fDCAxy","DCA xy distribution;# of CR",1200,-3,3);
        fOutputList->Add(fDCAxy);

	 //ESD track cut
        fesdTrackCuts =  AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
	fesdTrackCuts->SetMinNCrossedRowsTPC(fMinTPCcr);
	fesdTrackCuts->SetMaxChi2PerClusterTPC(fMaxChi2PerTPC);
	fesdTrackCuts->SetMaxDCAToVertexZ(fMaxDCAz);
	fesdTrackCuts->SetMaxDCAToVertexXYPtDep(fMaxDCAxy);

	//no DCAxy cut for secondaries
        fesdTrackCuts_no_dca =  AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);// no DCA xy cut
	fesdTrackCuts_no_dca->SetMinNCrossedRowsTPC(fMinTPCcr);
        fesdTrackCuts_no_dca->SetMaxChi2PerClusterTPC(fMaxChi2PerTPC);
        fesdTrackCuts_no_dca->SetMaxDCAToVertexZ(fMaxDCAz);



fPtV0MGenPion_signal_loss->Sumw2();
fPtV0MGenKaon_signal_loss->Sumw2();
fPtV0MGenProton_signal_loss->Sumw2();







	PostData(1, fOutputList);
}
//_________________________________________________________________________
bool IsPionNSigmaTOF(float pt, float nSigmaTOFPion, float TOFtime)
{
        if(pt>0.5){
        if (TMath::Abs(nSigmaTOFPion)<3)  return true;
}
  return false;
}

bool IsPionNSigmaTPC(float pt, float nSigmaTPCPion)
{
        if (TMath::Abs(nSigmaTPCPion)<3)  return true;

  return false;
}

bool IsKaonNSigmaTOF(float pt, float nSigmaTOFKaon, float TOFtime)
{
        if(pt>0.5){
     //rejection of unwanted contamination
      if(pt>1 && TOFtime<-400)  return false;

      if (TMath::Abs(nSigmaTOFKaon)<3)  return true;
}
  return false;
}

bool IsKaonNSigmaTPC(float pt, float nSigmaTPCKaon)
{
        if (TMath::Abs(nSigmaTPCKaon)<3) return true;

  return false;
}

bool IsProtonNSigmaTOF(float pt, float nSigmaTOFProton, float TOFtime)
{

        if(pt>0.5){
      if(pt>1.8 && TOFtime<-300)  return false;
        if (TMath::Abs(nSigmaTOFProton)<3)  return true;
}
        return false;
}

bool IsProtonNSigmaTPC(float pt, float nSigmaTPCProton)
{
        if (TMath::Abs(nSigmaTPCProton)<3) return true;
  return false;
}

//________________________________________________________________________
void AliAnalysisTaskTOFMC::UserExec(Option_t *)
{
  // Main loop
  // Called for each event
        AliESDEvent *fESD = 0x0;
          fESD = dynamic_cast<AliESDEvent*>(InputEvent());
          if (!fESD) {
            printf("ERROR: fESD not available\n");
            return;
}

	 AliESDVZERO *esdV0=fESD->GetVZEROData();

        fEventCounter->Fill(0.5);

	AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
	AliMCEvent* fMC = eventHandler->MCEvent();
        AliStack* fStack = fMC->Stack();

	 //pileup (type-1)
        Bool_t ispileup = fESD->IsPileupFromSPDInMultBins();
        if(!ispileup){
        //incomplete DAQ
        if (!fESD->IsIncompleteDAQ()){
        //tracklet vs cluster cut
        AliAnalysisUtils *AnalysisUtils = new AliAnalysisUtils();
        Double_t IsCluVstrk = AnalysisUtils->IsSPDClusterVsTrackletBG(fESD);
        if(!IsCluVstrk){
        
	//z vertex cut<10
        Bool_t IsVertex = selectVertex2015pp(fESD,kTRUE,kFALSE,kTRUE);
        if(IsVertex){

	const AliESDVertex * fVertex = fESD->GetPrimaryVertex();
        if (TMath::Abs(fVertex->GetZ())<10){

	Double_t INELgt0=(AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1.0) >= 1);
        if(INELgt0){

	fMultSelection_INEL = (AliMultSelection*) fESD->FindListObject("MultSelection"); // Esto es para 13 TeV
        if (!fMultSelection_INEL)
           //cout<<"------- No AliMultSelection Object Found --------"<<fMultSelection<<endl;
           AliWarning("No AliMultSelection Object Found --------");
        else
          v0mpc_inel = fMultSelection_INEL->GetMultiplicityPercentile("V0M",kFALSE);

//	if(AliPPVsMultUtils::IsINELgtZERO(fESD)){
	//if(fPPVsMultUtils->IsINELgtZERO(fESD)){
	Double_t yGenPi2=0;
        Double_t yGenK2=0;
        Double_t yGenP2=0;
	Int_t nTrack=0;

	for(int iPart = 1; iPart < (fMC->GetNumberOfTracks()); iPart++) {
        AliMCParticle *mcPart  = (AliMCParticle*)fMC->GetTrack(iPart);

        Double_t eta=mcPart->Eta();
        Double_t pt=mcPart->Pt();
        Int_t pdgcode=mcPart->PdgCode();
        Int_t label = mcPart->GetLabel();
	Double_t pz=mcPart->Pz();

	Double_t p=TMath::Sqrt(mcPart->Px()*mcPart->Px()+mcPart->Py()*mcPart->Py()+mcPart->Pz()*mcPart->Pz());
        Double_t eK = TMath::Sqrt(p*p + AliPID::ParticleMass(AliPID::kKaon)*AliPID::ParticleMass(AliPID::kKaon));
        Double_t ePi = TMath::Sqrt(p*p + AliPID::ParticleMass(AliPID::kPion)*AliPID::ParticleMass(AliPID::kPion));
        Double_t eP = TMath::Sqrt(p*p + AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton));
        if (eK -pz!=0.) yGenK2 =0.5*TMath::Log((eK + pz)/(eK -pz));
        if (ePi -pz!=0.) yGenPi2 =0.5*TMath::Log((ePi + pz)/(ePi -pz));
        if (eP -pz!=0.) yGenP2 =0.5*TMath::Log((eP + pz)/(eP -pz));

        if(fStack->IsPhysicalPrimary(label)){
	if (TMath::Abs(eta)<0.8){
	if (TMath::Abs(yGenPi2)<0.5){
	if(TMath::Abs(pdgcode)==211){
	fPtV0MGenPion_inel->Fill(pt,v0mpc_inel);
}//pdg
}//y
	if (TMath::Abs(yGenK2)<0.5){
        if(TMath::Abs(pdgcode)==321){
	fPtV0MGenKaon_inel->Fill(pt,v0mpc_inel);
}//pdg
}//y
	if (TMath::Abs(yGenP2)<0.5){
        if(TMath::Abs(pdgcode)==2212){
	fPtV0MGenProton_inel->Fill(pt,v0mpc_inel);
}//pdg
}//y
}//eta cut
}//primary
}//track loop

}//inel
}//vertex cut
}//vertex 10cut
}//BG cut
}//DAQ cut
}//Pile up
//#####################################################physics selection###################################################################
        UInt_t maskPhysSel = ((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
        maskPhysSel &= fTrigSel;
        TString firedTriggerClasses = fESD->GetFiredTriggerClasses();
        if (maskPhysSel != fTrigSel)
          return ;
        fEventCounter->Fill(1.5);

	Double_t INELgt0=(AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1.0) >= 1);
        if(!INELgt0)
           return;
        fEventCounter->Fill(2.5);


        //incomplete DAQ
        if (fESD->IsIncompleteDAQ())
          return;
        fEventCounter->Fill(3.5);

        //tracklet vs cluster cut
        AliAnalysisUtils *AnalysisUtils = new AliAnalysisUtils();
        Double_t IsCluVstrk = AnalysisUtils->IsSPDClusterVsTrackletBG(fESD);
        if(IsCluVstrk)
          return;
        fEventCounter->Fill(4.5);
	
	//pileup (type-1)
        //Bool_t ispileup = fESD->IsPileupFromSPDInMultBins();
        if(ispileup)
          return ;
        fEventCounter->Fill(5.5);

        fMultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection"); // Esto es para 13 TeV
        if (!fMultSelection)
           //cout<<"------- No AliMultSelection Object Found --------"<<fMultSelection<<endl;
           AliWarning("No AliMultSelection Object Found --------");
        else
          v0mpc = fMultSelection->GetMultiplicityPercentile("V0M",kFALSE);

        const AliESDVertex * fVertex = fESD->GetPrimaryVertex();

        fZVertex->Fill(fVertex->GetZ());
        fEventV0MPS->Fill(1,v0mpc);
        fEventPS->Fill(1);

	fV0MPC->Fill(v0mpc);

        Bool_t IsVertex = selectVertex2015pp(fESD,kTRUE,kFALSE,kTRUE);
        if(!IsVertex)
          return ;

        fEventCounter->Fill(6.5);
        fEventVtx->Fill(1);
        fEventV0MVtx->Fill(1,v0mpc);
        fZVertexRec->Fill(fVertex->GetZ());
  	
	fV0MPC_vertexcut->Fill(v0mpc);


//	const AliESDVertex * fVertex = fESD->GetPrimaryVertex();
        if (TMath::Abs(fVertex->GetZ())>10) return;

        //z vertex cut<10
        fEventCounter->Fill(7.5);
        fEventVtx10->Fill(1);
        fEventV0M->Fill(1,v0mpc);

        fZVertex10->Fill(fVertex->GetZ());
        fEventCounter->Fill(8.5);


        Double_t refmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);

        fCorrRefMultVsV0M->Fill(v0mpc,refmult08);

//=======================================================Rec============================================================

	
	for (Int_t i = 0; i < 3; i++) fTOFLabel[i];//=track->fTOFLabel[i];

	fPIDResponse->SetTOFResponse(fESD,AliPIDResponse::kBest_T0);

	Double_t yRecPi=0.;
        Double_t yRecK=0.;
        Double_t yRecP=0.;
	for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
        AliESDtrack *track= fESD->GetTrack(iTracks);
                if (!track) {
                 printf("ERROR: Could not receive track %d\n", iTracks);
         continue;
}

	Int_t label=TMath::Abs(track->GetLabel());
        TParticle *trackMC = fStack->Particle(TMath::Abs(track->GetLabel()));
	if(!trackMC) return;

        Int_t pdg = trackMC->GetPdgCode();

	Double_t pt=track->Pt();
	Double_t eta=track->Eta();


	Float_t dcaxy[2];
        Float_t dcaz[3];
        track->GetImpactParameters(dcaxy,dcaz);

	yRecPi=Rapidity(track,AliPID::ParticleMass(AliPID::kPion));
	yRecK=Rapidity(track,AliPID::ParticleMass(AliPID::kKaon));
	yRecP=Rapidity(track,AliPID::ParticleMass(AliPID::kProton));

	Double_t nSigmaTOFPion=fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);
        Double_t nSigmaTOFKaon=fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon);
        Double_t nSigmaTOFProton=fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton);

	Double_t nSigmaTPCPion=fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion);
        Double_t nSigmaTPCKaon=fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon);
        Double_t nSigmaTPCProton=fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton);
        Double_t nSigmaTPCElectron=fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);

	
	Bool_t TOFPIDStatus = TOFPID(track);
	Bool_t TPCPIDStatus = TPCPID(track);
	
        if(fesdTrackCuts_no_dca->AcceptTrack(track)){//for contamination
	if(TOFPIDStatus){
	if(TMath::Abs(eta)<0.8){
        if(TMath::Abs(yRecPi)<0.5){
        if(pt>=0.5 && pt<10.){
	if(TMath::Abs(pdg)==211){

	if(fStack->IsPhysicalPrimary(label)) fPtV0MDCAxyTOFPriPion->Fill(pt,v0mpc,dcaxy[0]);
	if(fStack->IsSecondaryFromWeakDecay(label)) fPtV0MDCAxyTOFWeakPion->Fill(pt,v0mpc,dcaxy[0]);
	if(fStack->IsSecondaryFromMaterial(label)) fPtV0MDCAxyTOFMatPion->Fill(pt,v0mpc,dcaxy[0]);

	if(track->Charge()>0){
	if(fStack->IsPhysicalPrimary(label)) fPtV0MDCAxyTOFPriPionP->Fill(pt,v0mpc,dcaxy[0]);
	if(fStack->IsSecondaryFromWeakDecay(label)) fPtV0MDCAxyTOFWeakPionP->Fill(pt,v0mpc,dcaxy[0]);
        if(fStack->IsSecondaryFromMaterial(label)) fPtV0MDCAxyTOFMatPionP->Fill(pt,v0mpc,dcaxy[0]);
}//positive
	if(track->Charge()<0){
	if(fStack->IsPhysicalPrimary(label)) fPtV0MDCAxyTOFPriPionM->Fill(pt,v0mpc,dcaxy[0]);
        if(fStack->IsSecondaryFromWeakDecay(label)) fPtV0MDCAxyTOFWeakPionM->Fill(pt,v0mpc,dcaxy[0]);
        if(fStack->IsSecondaryFromMaterial(label)) fPtV0MDCAxyTOFMatPionM->Fill(pt,v0mpc,dcaxy[0]);
}//negative
}//pdg
}//pt
}//yPion

	if(TMath::Abs(yRecP)<0.5){
        if(pt>=0.5 && pt<10.){
	if(TMath::Abs(pdg)==2212){

        if(fStack->IsPhysicalPrimary(label)) fPtV0MDCAxyTOFPriProton->Fill(pt,v0mpc,dcaxy[0]);
        if(fStack->IsSecondaryFromWeakDecay(label)) fPtV0MDCAxyTOFWeakProton->Fill(pt,v0mpc,dcaxy[0]);
        if(fStack->IsSecondaryFromMaterial(label)) fPtV0MDCAxyTOFMatProton->Fill(pt,v0mpc,dcaxy[0]);

        if(track->Charge()>0){
        if(fStack->IsPhysicalPrimary(label)) fPtV0MDCAxyTOFPriProtonP->Fill(pt,v0mpc,dcaxy[0]);
        if(fStack->IsSecondaryFromWeakDecay(label)) fPtV0MDCAxyTOFWeakProtonP->Fill(pt,v0mpc,dcaxy[0]);
        if(fStack->IsSecondaryFromMaterial(label)) fPtV0MDCAxyTOFMatProtonP->Fill(pt,v0mpc,dcaxy[0]);
}//positive

        if(track->Charge()<0){
        if(fStack->IsPhysicalPrimary(label)) fPtV0MDCAxyTOFPriProtonM->Fill(pt,v0mpc,dcaxy[0]);
        if(fStack->IsSecondaryFromWeakDecay(label)) fPtV0MDCAxyTOFWeakProtonM->Fill(pt,v0mpc,dcaxy[0]);
        if(fStack->IsSecondaryFromMaterial(label)) fPtV0MDCAxyTOFMatProtonM->Fill(pt,v0mpc,dcaxy[0]);
}//negative
}//pdg
}//pt
}//y
}//eta
}//tofpid
}//no dca xy cut for contamination estimation


        if(fesdTrackCuts->AcceptTrack(track)){	
	
	fTPC_CR->Fill(track->GetTPCCrossedRows());
	fChi2TPCcluster->Fill(track->GetTPCchi2()/track->GetTPCNcls());
	fDCAZ->Fill(dcaxy[1]);
	fDCAxy->Fill(dcaxy[0]);	

	if(TMath::Abs(track->Eta())<0.8){
        if(TPCPIDStatus){
        if(track->Charge()>0.) fPtTPC_AllP->Fill(track->Pt());
        if(track->Charge()<0.) fPtTPC_AllN->Fill(track->Pt());
}
        if(TOFPIDStatus){
        if(track->Charge()>0.) fPtTOF_AllP->Fill(track->Pt());
        if(track->Charge()<0.) fPtTOF_AllN->Fill(track->Pt());
}
}

	double pidTime[5]; track->GetIntegratedTimes(pidTime);
        Double_t t=track->GetTOFsignal();
        bool isPionNsigmaTOF = 0;
        bool isKaonNsigmaTOF = 0;
        bool isProtonNsigmaTOF  = 0;
        bool isPionNsigmaTOF_excl = 0;
        bool isKaonNsigmaTOF_excl = 0;
        bool isProtonNsigmaTOF_excl  = 0;


        isPionNsigmaTOF = (IsPionNSigmaTOF(track->Pt(), nSigmaTOFPion, t-pidTime[2]));
        isKaonNsigmaTOF = (IsKaonNSigmaTOF(track->Pt(),nSigmaTOFKaon, t-pidTime[3]));
        isProtonNsigmaTOF = (IsProtonNSigmaTOF(track->Pt(),nSigmaTOFProton, t-pidTime[4]));

         isPionNsigmaTOF_excl = (IsPionNSigmaTOF(track->Pt(), nSigmaTOFPion, t-pidTime[2]) && !IsKaonNSigmaTOF(track->Pt(),nSigmaTOFKaon, t-pidTime[3]) && !IsProtonNSigmaTOF(track->Pt(),nSigmaTOFProton, t-pidTime[4]));
        isKaonNsigmaTOF_excl = (!IsPionNSigmaTOF(track->Pt(),nSigmaTOFPion, t-pidTime[2]) && IsKaonNSigmaTOF(track->Pt(),nSigmaTOFKaon, t-pidTime[3]) && !IsProtonNSigmaTOF(track->Pt(), nSigmaTOFProton, t-pidTime[4]));
        isProtonNsigmaTOF_excl = (!IsPionNSigmaTOF(track->Pt(),nSigmaTOFPion, t-pidTime[2]) && !IsKaonNSigmaTOF(track->Pt(),nSigmaTOFKaon, t-pidTime[3]) && IsProtonNSigmaTOF(track->Pt(),nSigmaTOFProton, t-pidTime[4]));



	if(TPCPIDStatus){

	Double_t pt=track->Pt();
        Double_t pTPC=track->GetTPCmomentum();
        Double_t q=track->Charge();

        fdEdxP->Fill(pTPC,track->GetTPCsignal());
        fdEdxPt->Fill(pt,track->GetTPCsignal());
        fdEdxPq->Fill(pTPC/q,track->GetTPCsignal());
        fdEdxPtq->Fill(pt/q,track->GetTPCsignal());
}	

	if(TOFPIDStatus){
	 Double_t fmom=track->GetP();
        Double_t texpPion=fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kPion);
        Double_t texpKaon=fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kKaon);
        Double_t texpProton=fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kProton);
        Double_t t=track->GetTOFsignal();
        //Double_t t0 = fPIDResponse->SetTOFResponse(fESD,AliPIDResponse::kBest_T0);
        Double_t t0 = fPIDResponse->GetTOFResponse().GetStartTime(fmom); // T0best time

        Double_t dtPion=(t-texpPion-t0);
        Double_t dtKaon=(t-texpKaon-t0);
        Double_t dtProton=(t-texpProton-t0);

	fPtVsTPion->Fill(pt,dtPion);
	fPtVsTKaon->Fill(pt,dtKaon);
	fPtVsTProton->Fill(pt,dtProton);


	Double_t tofexptime=fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kElectron);
	Double_t tof= t*1E-3; // ns, average T0 fill subtracted, no info from T0detector     
	Double_t c=TMath::C()*1.E-9;// m/ns
	Double_t  length=tofexptime*1E-3*c;//m
	tof=tof*c;
	Double_t beta=length/tof;
        Double_t q=track->Charge();

        fbetaAllPt->Fill(track->Pt(),beta);
        fbetaAllP->Fill(track->P(),beta);
        fbetaAllPtq->Fill(track->Pt()/q,beta);
        fbetaAllPq->Fill(track->P()/q,beta);
}//tof pid



	if(fStack->IsPhysicalPrimary(label)){
	if(TMath::Abs(eta)<0.8){
        if(TMath::Abs(yRecPi)<0.5){
        if(pt>=0.2 && pt<10.){
        if(TPCPIDStatus){
	if( TMath::Abs(pdg)==211){
	
	fPtV0MTPCRecPion->Fill(pt,v0mpc);
	if(track->Charge()>0) fPtV0MTPCRecPionP->Fill(pt,v0mpc);
	if(track->Charge()<0) fPtV0MTPCRecPionM->Fill(pt,v0mpc);
}//pid
}//tpc

       	if(TOFPIDStatus){
	if( TMath::Abs(pdg)==211){

        fPtV0MTOFRecPion->Fill(pt,v0mpc);
        if(track->Charge()>0) fPtV0MTOFRecPionP->Fill(pt,v0mpc);
        if(track->Charge()<0) fPtV0MTOFRecPionM->Fill(pt,v0mpc);
  
	if (isPionNsigmaTOF){      
        fPtV0MTOFRecPion_nSigma->Fill(pt,v0mpc);
        if(track->Charge()>0) fPtV0MTOFRecPionP_nSigma->Fill(pt,v0mpc);
        if(track->Charge()<0) fPtV0MTOFRecPionM_nSigma->Fill(pt,v0mpc);
}
	if (isPionNsigmaTOF_excl){
        fPtV0MTOFRecPion_nSigma_excl->Fill(pt,v0mpc);
        if(track->Charge()>0) fPtV0MTOFRecPionP_nSigma_excl->Fill(pt,v0mpc);
        if(track->Charge()<0) fPtV0MTOFRecPionM_nSigma_excl->Fill(pt,v0mpc);
}
}//pid
}//tof
}//pt
}//rap

	if(TMath::Abs(yRecK)<0.5){
        if(pt>=0.2 && pt<10.){
        if(TPCPIDStatus){
	if( TMath::Abs(pdg)==321){

        fPtV0MTPCRecKaon->Fill(pt,v0mpc);
        if(track->Charge()>0) fPtV0MTPCRecKaonP->Fill(pt,v0mpc);
        if(track->Charge()<0) fPtV0MTPCRecKaonM->Fill(pt,v0mpc);
}//pid
}//tpc

        if(TOFPIDStatus){
	if( TMath::Abs(pdg)==321){

        fPtV0MTOFRecKaon->Fill(pt,v0mpc);
        if(track->Charge()>0) fPtV0MTOFRecKaonP->Fill(pt,v0mpc);
        if(track->Charge()<0) fPtV0MTOFRecKaonM->Fill(pt,v0mpc);
        
	if (isKaonNsigmaTOF){
        fPtV0MTOFRecKaon_nSigma->Fill(pt,v0mpc);
        if(track->Charge()>0) fPtV0MTOFRecKaonP_nSigma->Fill(pt,v0mpc);
        if(track->Charge()<0) fPtV0MTOFRecKaonM_nSigma->Fill(pt,v0mpc);
}
	
	if (isKaonNsigmaTOF_excl){	
        fPtV0MTOFRecKaon_nSigma_excl->Fill(pt,v0mpc);
        if(track->Charge()>0) fPtV0MTOFRecKaonP_nSigma_excl->Fill(pt,v0mpc);
        if(track->Charge()<0) fPtV0MTOFRecKaonM_nSigma_excl->Fill(pt,v0mpc);
}
}//pid
}//tof
}//pt
}//y
	if(TMath::Abs(yRecP)<0.5){
        if(pt>=0.2 && pt<10.){
        if(TPCPIDStatus){
	if( TMath::Abs(pdg)==2212){

        fPtV0MTPCRecProton->Fill(pt,v0mpc);
        if(track->Charge()>0) fPtV0MTPCRecProtonP->Fill(pt,v0mpc);
        if(track->Charge()<0) fPtV0MTPCRecProtonM->Fill(pt,v0mpc);
}//pid
}//tpcpid

        if(TOFPIDStatus){
        if( TMath::Abs(pdg)==2212){

        fPtV0MTOFRecProton->Fill(pt,v0mpc);
        if(track->Charge()>0) fPtV0MTOFRecProtonP->Fill(pt,v0mpc);
        if(track->Charge()<0) fPtV0MTOFRecProtonM->Fill(pt,v0mpc);
        
	if (isProtonNsigmaTOF){
        fPtV0MTOFRecProton_nSigma->Fill(pt,v0mpc);
        if(track->Charge()>0) fPtV0MTOFRecProtonP_nSigma->Fill(pt,v0mpc);
        if(track->Charge()<0) fPtV0MTOFRecProtonM_nSigma->Fill(pt,v0mpc);
}
	if (isProtonNsigmaTOF_excl){
        fPtV0MTOFRecProton_nSigma_excl->Fill(pt,v0mpc);
        if(track->Charge()>0) fPtV0MTOFRecProtonP_nSigma_excl->Fill(pt,v0mpc);
        if(track->Charge()<0) fPtV0MTOFRecProtonM_nSigma_excl->Fill(pt,v0mpc);
}
}//pid
}//tof
}//pt
}//y
}//eta
}//primary
		

	if(TMath::Abs(eta)<0.8){
        if(TOFPIDStatus){
	//PID with the TOF
        Double_t fTOFTime = track->GetTOFsignal(); //Gets the TOF signal
        Double_t fT0TrkTime = fPIDResponse->GetTOFResponse().GetStartTime(track->P()); // T0best time

        Double_t pidTime[6]; track->GetIntegratedTimes(pidTime,6);
        //Double_t TExpTimeEl=fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kElectron); or
        Double_t fTOFExpTimePi = pidTime[AliPID::kPion];
        Double_t fTOFExpTimeK = pidTime[AliPID::kKaon];
        Double_t fTOFExpTimeP = pidTime[AliPID::kProton];

	const Double_t TdiffK = fTOFTime-fT0TrkTime-fTOFExpTimeK;
	const Double_t TdiffP = fTOFTime-fT0TrkTime-fTOFExpTimeP;


	//newly added on 12may
	Double_t fMCTOFTime=0;

	if (fMC && fTOFLabel[0] > 0) {
	track->GetTOFLabel(fTOFLabel);

	TParticle *particle;
	TClonesArray *arrayTR;
	AliTrackReference *trackRef;

	
	fMC->GetParticleAndTR(fTOFLabel[0], particle, arrayTR);
        //if (fMC->GetParticleAndTR(fTOFLabel[0], particle, arrayTR)<0) return
        if (fMC->GetParticleAndTR(fTOFLabel[0], particle, arrayTR)>0){
    	for (Int_t itr = 0; itr < arrayTR->GetEntries(); itr++) {
	trackRef = (AliTrackReference *)arrayTR->At(itr);
	if (!trackRef || trackRef->DetectorId() != AliTrackReference::kTOF) continue;

	fMCTOFTime = trackRef->GetTime() * 1.e12; /* s -> ps */
//      fMCTOFLength = trackRef->GetLength();
      /* break as soon as we get it */
      break;
}
    }
}
	

	if(TMath::Abs(yRecPi)<0.5){
	const Double_t TdiffPi = fMCTOFTime-fT0TrkTime-fTOFExpTimePi;
        if(( TMath::Abs(pdg)==211) && fStack->IsPhysicalPrimary(label)) fTOFTimeV0MPtPi->Fill(track->Pt(),TdiffPi,v0mpc);
        if(( TMath::Abs(pdg)==211) && !fStack->IsPhysicalPrimary(label)) fTOFTimeV0MPtMismatchDecayPi->Fill(track->Pt(),TdiffPi,v0mpc);
}//rap
	
	if(TMath::Abs(yRecK)<0.5){
        const Double_t TdiffK = fMCTOFTime-fT0TrkTime-fTOFExpTimeK;
        if(( TMath::Abs(pdg)==321) && fStack->IsPhysicalPrimary(label)) fTOFTimeV0MPtK->Fill(track->Pt(),TdiffK,v0mpc);
        if(( TMath::Abs(pdg)==321) && !fStack->IsPhysicalPrimary(label)) fTOFTimeV0MPtMismatchDecayK->Fill(track->Pt(),TdiffK,v0mpc);
}//rap
	
	if(TMath::Abs(yRecP)<0.5){
        const Double_t TdiffP = fMCTOFTime-fT0TrkTime-fTOFExpTimeP;
        if(( TMath::Abs(pdg)==2212) && fStack->IsPhysicalPrimary(label)) fTOFTimeV0MPtP->Fill(track->Pt(),TdiffP,v0mpc);
        if(( TMath::Abs(pdg)==2212) && !fStack->IsPhysicalPrimary(label)) fTOFTimeV0MPtMismatchDecayP->Fill(track->Pt(),TdiffP,v0mpc);
}//rap

}//tof
}//eta
}//track cut with dcaxy
}//track loop
//=======================================================Gen============================================================

	Double_t yGenPi=0;
        Double_t yGenK=0;
        Double_t yGenP=0;
	Double_t nTrack=0;

	for(int iPart = 1; iPart < (fMC->GetNumberOfTracks()); iPart++) {
        AliMCParticle *mcPart  = (AliMCParticle*)fMC->GetTrack(iPart);

        Double_t eta=mcPart->Eta();
        Double_t pt=mcPart->Pt();
        Int_t pdgcode=mcPart->PdgCode();
        Int_t label = mcPart->GetLabel();
	Double_t pz=mcPart->Pz();

	Double_t p=TMath::Sqrt(mcPart->Px()*mcPart->Px()+mcPart->Py()*mcPart->Py()+mcPart->Pz()*mcPart->Pz());

        Double_t eK = TMath::Sqrt(p*p + AliPID::ParticleMass(AliPID::kKaon)*AliPID::ParticleMass(AliPID::kKaon));
        Double_t ePi = TMath::Sqrt(p*p + AliPID::ParticleMass(AliPID::kPion)*AliPID::ParticleMass(AliPID::kPion));
        Double_t eP = TMath::Sqrt(p*p + AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton));
        if (eK -pz!=0.) yGenK =0.5*TMath::Log((eK + pz)/(eK -pz));
        if (ePi -pz!=0.) yGenPi =0.5*TMath::Log((ePi + pz)/(ePi -pz));
        if (eP -pz!=0.) yGenP =0.5*TMath::Log((eP + pz)/(eP -pz));

	if (TMath::Abs(eta)<0.8){
        if(fStack->IsPhysicalPrimary(label)){
	if(pt>0)	nTrack++;

}
}//eta


        if(fStack->IsPhysicalPrimary(label)){
	if (TMath::Abs(eta)<0.8){
	if (TMath::Abs(yGenPi)<0.5){
	if(TMath::Abs(pdgcode)==211){
	
	fPtV0MGenPion_kINT7->Fill(pt,v0mpc);

        if (pt>=0.2 && pt<10. ) {
	
	fPtV0MGenPion->Fill(pt,v0mpc);
	if(mcPart->Charge()<0) fPtV0MGenPionM->Fill(pt,v0mpc);
	if(mcPart->Charge()>0)	fPtV0MGenPionP->Fill(pt,v0mpc);
}//pt
}//pdg
}//y
	if (TMath::Abs(yGenK)<0.5){
        if(TMath::Abs(pdgcode)==321){

	fPtV0MGenKaon_kINT7->Fill(pt,v0mpc);

        if (pt>=0.2 && pt<10. ) {
        fPtV0MGenKaon->Fill(pt,v0mpc);
        if(mcPart->Charge()>0) fPtV0MGenKaonP->Fill(pt,v0mpc);
        if(mcPart->Charge()<0) fPtV0MGenKaonM->Fill(pt,v0mpc);
}//pdg
}//pt
}//y
	if (TMath::Abs(yGenP)<0.5){
        if(TMath::Abs(pdgcode)==2212){

	fPtV0MGenProton_kINT7->Fill(pt,v0mpc);
        if (pt>=0.2 && pt<10. ) {
        fPtV0MGenProton->Fill(pt,v0mpc);
        if(mcPart->Charge()>0) fPtV0MGenProtonP->Fill(pt,v0mpc);
        if(mcPart->Charge()<0) fPtV0MGenProtonM->Fill(pt,v0mpc);
}//pt
}//pdg
}//y

}//eta cut
}//primary
}//track loop

if(nTrack>0 && refmult08>0)	fCorrRefMultVsNch->Fill(nTrack,refmult08);

  PostData(1, fOutputList);
}
//_________________________________________________________________________
Float_t AliAnalysisTaskTOFMC::GetVertex(AliESDEvent* esd) const
{
  Float_t zvtx = -999;
  //const AliAODVertex* vtxAOD = aod->GetPrimaryVertex();
 const AliESDVertex* vtxESD = esd->GetPrimaryVertex();
  if (!vtxESD)
    return zvtx;
  if(vtxESD->GetNContributors()>0)
    zvtx = vtxESD->GetZ();
  return zvtx;
}
//________________________________________________________________________
void AliAnalysisTaskTOFMC::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
	fEventCounter = dynamic_cast<TH1F*> (fOutputList->At(0));
	fEventPS = dynamic_cast<TH1F*> (fOutputList->At(1));
        fEventVtx = dynamic_cast<TH1F*> (fOutputList->At(2));
        fEventVtx10 = dynamic_cast<TH1F*> (fOutputList->At(3));
        fZVertex = dynamic_cast<TH1F*> (fOutputList->At(4));
        fZVertexRec = dynamic_cast<TH1F*> (fOutputList->At(5));
        fZVertexEff = dynamic_cast<TH1F*> (fOutputList->At(6));
        fZVertex10 = dynamic_cast<TH1F*> (fOutputList->At(7));
        fV0MPC = dynamic_cast<TH1F*> (fOutputList->At(8));

	fCorrRefMultVsNch = dynamic_cast<TH2F*> (fOutputList->At(9));
	fCorrRefMultVsV0M = dynamic_cast<TH2F*> (fOutputList->At(10));
	fPtV0MGenPion = dynamic_cast<TH2F*> (fOutputList->At(11));
	fPtV0MTPCRecPion = dynamic_cast<TH2F*> (fOutputList->At(12));
	fPtV0MTOFRecPion = dynamic_cast<TH2F*> (fOutputList->At(13));
	fPtV0MGenPionP = dynamic_cast<TH2F*> (fOutputList->At(14));
	fPtV0MTPCRecPionP = dynamic_cast<TH2F*> (fOutputList->At(15));
	fPtV0MTOFRecPionP = dynamic_cast<TH2F*> (fOutputList->At(16));
	fPtV0MGenPionM = dynamic_cast<TH2F*> (fOutputList->At(17));
	fPtV0MTPCRecPionM = dynamic_cast<TH2F*> (fOutputList->At(18));
	fPtV0MTOFRecPionM = dynamic_cast<TH2F*> (fOutputList->At(19));

	fPtV0MGenKaon = dynamic_cast<TH2F*> (fOutputList->At(20));
	fPtV0MTPCRecKaon = dynamic_cast<TH2F*> (fOutputList->At(21));
	fPtV0MTOFRecKaon = dynamic_cast<TH2F*> (fOutputList->At(22));
	fPtV0MGenKaonP = dynamic_cast<TH2F*> (fOutputList->At(23));
	fPtV0MTPCRecKaonP = dynamic_cast<TH2F*> (fOutputList->At(24));
	fPtV0MTOFRecKaonP = dynamic_cast<TH2F*> (fOutputList->At(25));
	fPtV0MGenKaonM = dynamic_cast<TH2F*> (fOutputList->At(26));
	fPtV0MTPCRecKaonM = dynamic_cast<TH2F*> (fOutputList->At(27));
	fPtV0MTOFRecKaonM = dynamic_cast<TH2F*> (fOutputList->At(28));

	fPtV0MGenProton = dynamic_cast<TH2F*> (fOutputList->At(29));
	fPtV0MTPCRecProton = dynamic_cast<TH2F*> (fOutputList->At(30));
	fPtV0MTOFRecProton = dynamic_cast<TH2F*> (fOutputList->At(31));
	fPtV0MGenProtonP = dynamic_cast<TH2F*> (fOutputList->At(32));
	fPtV0MTPCRecProtonP = dynamic_cast<TH2F*> (fOutputList->At(33));
	fPtV0MTOFRecProtonP = dynamic_cast<TH2F*> (fOutputList->At(34));
	fPtV0MGenProtonM = dynamic_cast<TH2F*> (fOutputList->At(35));
	fPtV0MTPCRecProtonM = dynamic_cast<TH2F*> (fOutputList->At(36));
	fPtV0MTOFRecProtonM = dynamic_cast<TH2F*> (fOutputList->At(37));

	fPtV0MDCAxyTOFPriPion = dynamic_cast<TH3F*> (fOutputList->At(38));
	fPtV0MDCAxyTOFWeakPion = dynamic_cast<TH3F*> (fOutputList->At(39));
	fPtV0MDCAxyTOFMatPion = dynamic_cast<TH3F*> (fOutputList->At(40));
	fPtV0MDCAxyTOFPriPionP = dynamic_cast<TH3F*> (fOutputList->At(41));
	fPtV0MDCAxyTOFWeakPionP = dynamic_cast<TH3F*> (fOutputList->At(42));
	fPtV0MDCAxyTOFMatPionP = dynamic_cast<TH3F*> (fOutputList->At(43));
	fPtV0MDCAxyTOFPriPionM = dynamic_cast<TH3F*> (fOutputList->At(44));
	fPtV0MDCAxyTOFWeakPionM = dynamic_cast<TH3F*> (fOutputList->At(45));
	fPtV0MDCAxyTOFMatPionM = dynamic_cast<TH3F*> (fOutputList->At(46));

	fPtV0MDCAxyTOFPriProton = dynamic_cast<TH3F*> (fOutputList->At(47));
	fPtV0MDCAxyTOFWeakProton = dynamic_cast<TH3F*> (fOutputList->At(48));
	fPtV0MDCAxyTOFMatProton = dynamic_cast<TH3F*> (fOutputList->At(49));
	fPtV0MDCAxyTOFPriProtonP = dynamic_cast<TH3F*> (fOutputList->At(50));
	fPtV0MDCAxyTOFWeakProtonP = dynamic_cast<TH3F*> (fOutputList->At(51));
	fPtV0MDCAxyTOFMatProtonP = dynamic_cast<TH3F*> (fOutputList->At(52));
	fPtV0MDCAxyTOFPriProtonM = dynamic_cast<TH3F*> (fOutputList->At(53));
	fPtV0MDCAxyTOFWeakProtonM = dynamic_cast<TH3F*> (fOutputList->At(54));
	fPtV0MDCAxyTOFMatProtonM = dynamic_cast<TH3F*> (fOutputList->At(55));


	fdEdxP = dynamic_cast<TH2F*> (fOutputList->At(56));
        fdEdxPt = dynamic_cast<TH2F*> (fOutputList->At(57));
        fdEdxPq = dynamic_cast<TH2F*> (fOutputList->At(58));
        fdEdxPtq = dynamic_cast<TH2F*> (fOutputList->At(59));
	fbetaAllPt = dynamic_cast<TH2F*> (fOutputList->At(60));
        fbetaAllP = dynamic_cast<TH2F*> (fOutputList->At(61));
	fbetaAllPtq = dynamic_cast<TH2F*> (fOutputList->At(62));
        fbetaAllPq = dynamic_cast<TH2F*> (fOutputList->At(63));

	fPtVsTPion = dynamic_cast<TH2F*> (fOutputList->At(64));
        fPtVsTKaon = dynamic_cast<TH2F*> (fOutputList->At(65));
        fPtVsTProton = dynamic_cast<TH2F*> (fOutputList->At(66));

        fPtV0MGenPion_kINT7 = dynamic_cast<TH2F*> (fOutputList->At(67));
        fPtV0MGenPion_inel = dynamic_cast<TH2F*> (fOutputList->At(68));
        fPtV0MGenPion_signal_loss = dynamic_cast<TH2F*> (fOutputList->At(69));
        fPtV0MGenKaon_kINT7 = dynamic_cast<TH2F*> (fOutputList->At(70));
        fPtV0MGenKaon_inel = dynamic_cast<TH2F*> (fOutputList->At(71));
        fPtV0MGenKaon_signal_loss = dynamic_cast<TH2F*> (fOutputList->At(72));
        fPtV0MGenProton_kINT7 = dynamic_cast<TH2F*> (fOutputList->At(73));
        fPtV0MGenProton_inel = dynamic_cast<TH2F*> (fOutputList->At(74));
        fPtV0MGenProton_signal_loss = dynamic_cast<TH2F*> (fOutputList->At(75));
	


	fPtV0MTOFRecPion_nSigma = dynamic_cast<TH2F*> (fOutputList->At(76));
	fPtV0MTOFRecPionP_nSigma = dynamic_cast<TH2F*> (fOutputList->At(77));
	fPtV0MTOFRecPionM_nSigma = dynamic_cast<TH2F*> (fOutputList->At(78));
	fPtV0MTOFRecKaon_nSigma = dynamic_cast<TH2F*> (fOutputList->At(79));
	fPtV0MTOFRecKaonP_nSigma = dynamic_cast<TH2F*> (fOutputList->At(80));
	fPtV0MTOFRecKaonM_nSigma = dynamic_cast<TH2F*> (fOutputList->At(81));
	fPtV0MTOFRecProton_nSigma = dynamic_cast<TH2F*> (fOutputList->At(82));
	fPtV0MTOFRecProtonP_nSigma = dynamic_cast<TH2F*> (fOutputList->At(83));
	fPtV0MTOFRecProtonM_nSigma = dynamic_cast<TH2F*> (fOutputList->At(84));

	fPtV0MTOFRecPion_nSigma_excl = dynamic_cast<TH2F*> (fOutputList->At(85));
	fPtV0MTOFRecPionP_nSigma_excl = dynamic_cast<TH2F*> (fOutputList->At(86));
	fPtV0MTOFRecPionM_nSigma_excl = dynamic_cast<TH2F*> (fOutputList->At(87));
	fPtV0MTOFRecKaon_nSigma_excl = dynamic_cast<TH2F*> (fOutputList->At(88));
	fPtV0MTOFRecKaonP_nSigma_excl = dynamic_cast<TH2F*> (fOutputList->At(89));
	fPtV0MTOFRecKaonM_nSigma_excl = dynamic_cast<TH2F*> (fOutputList->At(90));
	fPtV0MTOFRecProton_nSigma_excl = dynamic_cast<TH2F*> (fOutputList->At(91));
	fPtV0MTOFRecProtonP_nSigma_excl = dynamic_cast<TH2F*> (fOutputList->At(92));
	fPtV0MTOFRecProtonM_nSigma_excl = dynamic_cast<TH2F*> (fOutputList->At(93));


	fTOFTimeV0MPtPi = dynamic_cast<TH3F*> (fOutputList->At(94));
	fTOFTimeV0MPtK = dynamic_cast<TH3F*> (fOutputList->At(95));
	fTOFTimeV0MPtP = dynamic_cast<TH3F*> (fOutputList->At(96));
	fTOFTimeV0MPtMismatchDecayPi = dynamic_cast<TH3F*> (fOutputList->At(97));
	fTOFTimeV0MPtMismatchDecayK = dynamic_cast<TH3F*> (fOutputList->At(98));
	fTOFTimeV0MPtMismatchDecayP = dynamic_cast<TH3F*> (fOutputList->At(99));

	fEventV0MPS= dynamic_cast<TH2F*> (fOutputList->At(100));
        fEventV0MVtx= dynamic_cast<TH2F*> (fOutputList->At(101));
        fEventV0M= dynamic_cast<TH2F*> (fOutputList->At(102));

        fV0MPC_vertexcut= dynamic_cast<TH1F*> (fOutputList->At(103));

        fPtTPC_AllP= dynamic_cast<TH1F*> (fOutputList->At(104));
        fPtTPC_AllN= dynamic_cast<TH1F*> (fOutputList->At(105));
        fPtTOF_AllP= dynamic_cast<TH1F*> (fOutputList->At(106));
        fPtTOF_AllN= dynamic_cast<TH1F*> (fOutputList->At(107));
        
	fTPC_CR= dynamic_cast<TH1F*> (fOutputList->At(108));
	fChi2TPCcluster= dynamic_cast<TH1F*> (fOutputList->At(109));
	fDCAZ= dynamic_cast<TH1F*> (fOutputList->At(110));
	fDCAxy= dynamic_cast<TH1F*> (fOutputList->At(111));



fPtV0MGenPion_signal_loss->Divide(fPtV0MGenPion_inel,fPtV0MGenPion_kINT7);
fPtV0MGenKaon_signal_loss->Divide(fPtV0MGenKaon_inel,fPtV0MGenKaon_kINT7);
fPtV0MGenProton_signal_loss->Divide(fPtV0MGenProton_inel,fPtV0MGenProton_kINT7);


	TCanvas *c1 = new TCanvas();
	c1->Divide(3,2);
  	c1->cd(1);
	fEventCounter->Draw();
	c1->cd(2);
        fEventPS->Draw();
        c1->cd(3);
        fEventVtx->Draw();
        c1->cd(4);
        fEventVtx10->Draw();


	TCanvas *c4=new TCanvas();
        c4->Divide(2,2);
        c4->cd(1);
        fCorrRefMultVsNch->Draw("colz");
        c4->cd(2);
        fCorrRefMultVsV0M->Draw("colz");
        c4->cd(3);
        fV0MPC->Draw();

	TCanvas *c5=new TCanvas();
        c5->Divide(3,2);
        c5->cd(1);
        fPtVsTPion->Draw();
        c5->cd(2);
        fPtVsTKaon->Draw();
        c5->cd(3);
        fPtVsTProton->Draw();
        c5->cd(4);
        fPtV0MGenPion_signal_loss->Draw();
        c5->cd(5);
        fPtV0MGenKaon_signal_loss->Draw();
        c5->cd(6);
        fPtV0MGenProton_signal_loss->Draw();

	TCanvas *cR5=new TCanvas();
        cR5->Divide(3,2);
        cR5->cd(1);
        fTOFTimeV0MPtPi->Draw();
        cR5->cd(2);
        fTOFTimeV0MPtMismatchDecayPi->Draw();
        cR5->cd(3);
        fTOFTimeV0MPtMismatchDecayK->Draw();
        cR5->cd(4);
        fTOFTimeV0MPtMismatchDecayP->Draw();


/*
TFile *f=new TFile("result/12dec/MC_final_TOF_output.root","recreate");//v0m bin different
f->cd();

fTPC_CR->Write();
fChi2TPCcluster->Write();
fDCAZ->Write();
fDCAxy->Write();

fEventCounter->Write();
fEventPS->Write();
fEventVtx->Write();
fEventVtx10->Write();
fZVertex->Write();
fZVertexRec->Write();
fZVertexEff->Write();
fZVertex10->Write();
fV0MPC->Write();
fV0MPC_vertexcut->Write();
fEventV0MPS->Write();
fEventV0MVtx->Write();
fEventV0M->Write();


fCorrRefMultVsNch->Write();
fCorrRefMultVsV0M->Write();
fPtVsTPion->Write();
fPtVsTKaon->Write();
fPtVsTProton->Write();
fdEdxP->Write();
fdEdxPt->Write();
fdEdxPq->Write();
fdEdxPtq->Write();
fbetaAllPt->Write();
fbetaAllP->Write();
fbetaAllPtq->Write();
fbetaAllPq->Write();

fPtGenPion_signal_loss->Write();
fPtGenKaon_signal_loss->Write();
fPtGenProton_signal_loss->Write();



fPtV0MGenPion->Write();
fPtV0MTPCRecPion->Write();
fPtV0MTOFRecPion->Write();
fPtV0MGenPionP->Write();
fPtV0MTPCRecPionP->Write();
fPtV0MTOFRecPionP->Write();
fPtV0MGenPionM->Write();
fPtV0MTPCRecPionM->Write();
fPtV0MTOFRecPionM->Write();
fPtV0MGenKaon->Write();
fPtV0MTPCRecKaon->Write();
fPtV0MTOFRecKaon->Write();
fPtV0MGenKaonP->Write();
fPtV0MTPCRecKaonP->Write();
fPtV0MTOFRecKaonP->Write();
fPtV0MGenKaonM->Write();
fPtV0MTPCRecKaonM->Write();
fPtV0MTOFRecKaonM->Write();
fPtV0MGenProton->Write();
fPtV0MTPCRecProton->Write();
fPtV0MTOFRecProton->Write();
fPtV0MGenProtonP->Write();
fPtV0MTPCRecProtonP->Write();
fPtV0MTOFRecProtonP->Write();
fPtV0MGenProtonM->Write();
fPtV0MTPCRecProtonM->Write();
fPtV0MTOFRecProtonM->Write();

fPtV0MDCAxyTOFPriPion->Write();
fPtV0MDCAxyTOFWeakPion->Write();
fPtV0MDCAxyTOFMatPion->Write();
fPtV0MDCAxyTOFPriPionP->Write();
fPtV0MDCAxyTOFWeakPionP->Write();
fPtV0MDCAxyTOFMatPionP->Write();
fPtV0MDCAxyTOFPriPionM->Write();
fPtV0MDCAxyTOFWeakPionM->Write();
fPtV0MDCAxyTOFMatPionM->Write();
fPtV0MDCAxyTOFPriProton->Write();
fPtV0MDCAxyTOFWeakProton->Write();
fPtV0MDCAxyTOFMatProton->Write();
fPtV0MDCAxyTOFPriProtonP->Write();
fPtV0MDCAxyTOFWeakProtonP->Write();
fPtV0MDCAxyTOFMatProtonP->Write();
fPtV0MDCAxyTOFPriProtonM->Write();
fPtV0MDCAxyTOFWeakProtonM->Write();
fPtV0MDCAxyTOFMatProtonM->Write();

fPtV0MTOFRecPion_nSigma->Write();
fPtV0MTOFRecPionP_nSigma->Write();
fPtV0MTOFRecPionM_nSigma->Write();
fPtV0MTOFRecKaon_nSigma->Write();
fPtV0MTOFRecKaonP_nSigma->Write();
fPtV0MTOFRecKaonM_nSigma->Write();
fPtV0MTOFRecProton_nSigma->Write();
fPtV0MTOFRecProtonP_nSigma->Write();
fPtV0MTOFRecProtonM_nSigma->Write();

fPtV0MTOFRecPion_nSigma_excl->Write();
fPtV0MTOFRecPionP_nSigma_excl->Write();
fPtV0MTOFRecPionM_nSigma_excl->Write();
fPtV0MTOFRecKaon_nSigma_excl->Write();
fPtV0MTOFRecKaonP_nSigma_excl->Write();
fPtV0MTOFRecKaonM_nSigma_excl->Write();
fPtV0MTOFRecProton_nSigma_excl->Write();
fPtV0MTOFRecProtonP_nSigma_excl->Write();
fPtV0MTOFRecProtonM_nSigma_excl->Write();

fTOFTimeV0MPtPi->Write();
fTOFTimeV0MPtK->Write();
fTOFTimeV0MPtP->Write();
fTOFTimeV0MPtMismatchDecayPi->Write();
fTOFTimeV0MPtMismatchDecayK->Write();
fTOFTimeV0MPtMismatchDecayP->Write();

fPtTPC_AllP->Write();
fPtTPC_AllN->Write();
fPtTOF_AllP->Write();
fPtTOF_AllN->Write();	
*/

}
//-----------------------------------------------------------------
Bool_t AliAnalysisTaskTOFMC::selectVertex2015pp(AliESDEvent *esd,
                          Bool_t checkSPDres, //enable check on vtx resolution
                          Bool_t requireSPDandTrk, //ask for both trk and SPD vertex 
                          Bool_t checkProximity) //apply cut on relative position of spd and trk verteces 
{

  if (!esd) return kFALSE;

  const AliESDVertex * trkVertex = esd->GetPrimaryVertexTracks();
  const AliESDVertex * spdVertex = esd->GetPrimaryVertexSPD();
  Bool_t hasSPD = spdVertex->GetStatus();
  Bool_t hasTrk = trkVertex->GetStatus();

  //Note that AliVertex::GetStatus checks that N_contributors is > 0
  //reject events if both are explicitly requested and none is available
  if (requireSPDandTrk && !(hasSPD && hasTrk)) return kFALSE;

  //reject events if none between the SPD or track verteces are available
  //if no trk vertex, try to fall back to SPD vertex;
  if (!hasTrk) {
    if (!hasSPD) return kFALSE;
    //on demand check the spd vertex resolution and reject if not satisfied
    if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
  } else {
    if (hasSPD) {
      //if enabled check the spd vertex resolution and reject if not satisfied
      //if enabled, check the proximity between the spd vertex and trak vertex, and reject if not satisfied
      if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
      if ((checkProximity && TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())>0.5)) return kFALSE;
    }
  }

/*  //Cut on the vertex z position
  const AliESDVertex *vertex = esd->GetPrimaryVertex();
  if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
  return kTRUE;
*/
}
//_________________________________________________________________________________________________
Bool_t AliAnalysisTaskTOFMC::IsGoodSPDvertexRes(const AliESDVertex * spdVertex)
{
  if (!spdVertex) return kFALSE;
  if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25)) return kFALSE;
  return kTRUE;
}

//_________________________________________________________________________________________________

Bool_t AliAnalysisTaskTOFMC::TPCPID(AliESDtrack *track)
{
    if ((track->GetStatus() & AliESDtrack::kTPCin   ) == 0) return kFALSE;
    if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) return kFALSE;
    if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) return kFALSE;
    return kTRUE;
}

// =============== For TOF PID CHECK =================

Bool_t AliAnalysisTaskTOFMC:: TOFPID(AliESDtrack * track)
{
    Double_t TOFPtRange =0.5;         // Controling Agent 1

// *************** Check if the particle has TOF Matching ***************

    UInt_t status;
    status=track->GetStatus();

    if ((status & AliESDtrack::kITSrefit) == 0) return kFALSE;
    if ((status & AliESDtrack::kTPCrefit) == 0) return kFALSE;

    if((status&AliVTrack::kTOFout)==0 || (status&AliVTrack::kTIME)==0)
        return kFALSE;
 // TPC TOF mismatch is be implemented
    Float_t length = track->GetIntegratedLength();
    if (length < 350.)        
        return kFALSE;

// -------  in addition to KTOFout and kTIME we look at the pt  ------

//    if(track->Pt()<TOFPtRange) return kFALSE;
    return kTRUE;
}

//================================================================
Double_t AliAnalysisTaskTOFMC::Rapidity(AliESDtrack *track , Double_t mass)
{
    Double_t E,rap,pz,pt;
    pt=track->Pt();
    pz = track->Pz();
    E = TMath::Sqrt(pt*pt+pz*pz+mass*mass);
    rap = 0.5 * TMath::Log ((E+pz)/(E-pz));
    return rap;
}

