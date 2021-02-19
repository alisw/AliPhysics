#include "AliAnalysisTaskNonlinearFlow.h"
#include "AliGFWWeights.h"
#include "CorrelationCalculator.h"

#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TComplex.h>
#include <TBits.h>
// AliRoot includes
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"

// ROOT includes
#include <TList.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TMatrixDSym.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TComplex.h>

// AliRoot includes
#include "AliEventCuts.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisUtils.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliAODITSsaTrackCuts.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliCentrality.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliVVertex.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "AliMultSelection.h"
#include "AliInputEventHandler.h"

// STL includes
#include <iostream>
using std::cout;
using std::endl;
ClassImp(AliAnalysisTaskNonlinearFlow)
//___________________________________________________________________________
AliAnalysisTaskNonlinearFlow::AliAnalysisTaskNonlinearFlow():
AliAnalysisTaskSE(),
fAOD(0),
fitssatrackcuts(0),
fFilterbit(96),
fFilterbitDefault(96),
fEtaCut(0.8),
fVtxCut(10.0),
fVtxCutDefault(10.0),
fMinPt(0.2),
fMaxPt(3.0),
fTPCclusters(70),
fChi2PerTPCcluster(10000),
fMinITSClus(0),
fMaxChi2(0),
fUseDCAzCut(0),
fDCAz(0),
fUseDCAxyCut(0),
fDCAxy(0),
fSample(1),
fCentFlag(0),
fTrigger(0),
fLS(false),
fNUE(0),
fNUA(0),
fNtrksName("Mult"),
//....
fPeriod("LHC15o"),

fListOfObjects(0),

fMultTOFLowCut(0),
fMultTOFHighCut(0),
fMultCentLowCut(0),

fTrackEfficiency(0),
hTrackEfficiency(0),
hTrackEfficiencyRun(0),

fFlowRunByRunWeights(false),
fFlowPeriodWeights(false),
fFlowUse3Dweights(false),
fFlowWeightsList(nullptr),

fPhiWeight(0),
fPhiWeightPlus(0),
fPhiWeightMinus(0),
hPhiWeight(0),
hPhiWeight1D(0),
hPhiWeight_LHC15i_part1(0),
hPhiWeight_LHC15i_part2(0),
hPhiWeight_LHC15j_part1(0),
hPhiWeight_LHC15j_part2(0),
hPhiWeight_LHC15l_part1(0),
hPhiWeight_LHC15l_part2(0),
hPhiWeight_LHC15l_part3(0),
hPhiWeightPlus_LHC15i_part1(0),
hPhiWeightPlus_LHC15i_part2(0),
hPhiWeightPlus_LHC15j_part1(0),
hPhiWeightPlus_LHC15j_part2(0),
hPhiWeightPlus_LHC15l_part1(0),
hPhiWeightPlus_LHC15l_part2(0),
hPhiWeightPlus_LHC15l_part3(0),
hPhiWeightMinus_LHC15i_part1(0),
hPhiWeightMinus_LHC15i_part2(0),
hPhiWeightMinus_LHC15j_part1(0),
hPhiWeightMinus_LHC15j_part2(0),
hPhiWeightMinus_LHC15l_part1(0),
hPhiWeightMinus_LHC15l_part2(0),
hPhiWeightMinus_LHC15l_part3(0),

hEventCount(0),
hMult(0),
fVtxAfterCuts(0),
fCentralityDis(0),
fV0CentralityDis(0),
hMultV0vsNtrksAfterCuts(0),
hMultSPDvsNtrksAfterCuts(0),
hNtrksVSmultPercentile(0),
fCentralityV0MCL1(0),
fCentralityV0MCL0(0),
fCentralityCL0CL1(0),
fMultvsCentr(0),
fMult128vsCentr(0),
fMultTPCvsTOF(0),
fMultTPCvsESD(0),

hSPDClsVsTrk(0),
hV0C012vsTkl(0),
hV0C012vsV0C3(0),
hV0MOnVsOf(0),
hSPDOnVsOf(0),

fPhiDis1D(0),
fPhiDis(0),
fEtaDis(0),
fEtaBefore(0),
fPtDis(0),
fPtBefore(0),
hDCAxyBefore(0),
hDCAzBefore(0),
hITSclustersBefore(0),
hChi2Before(0),
hDCAxy(0),
hDCAz(0),
hITSclusters(0),
hChi2(0),
rand(32213) 
 {

}
//______________________________________________________________________________
AliAnalysisTaskNonlinearFlow::AliAnalysisTaskNonlinearFlow(const char *name):
	AliAnalysisTaskSE(name),
	fAOD(0),
	fitssatrackcuts(0),
	fFilterbit(96),
	fEtaCut(0.8),
	fVtxCut(10.0),
	fMinPt(0.2),
	fMaxPt(3.0),
	fTPCclusters(70),
	fMinITSClus(0),
	fMaxChi2(0),
	fUseDCAzCut(0),
	fDCAz(0),
	fDCAzDefault(0),
	fUseDCAxyCut(0),
	fDCAxy(0),
	fDCAxyDefault(0),
	fSample(1),
	fCentFlag(0),
	fTrigger(0),
	fAliTrigger(0),
	fLS(false),
	fNUE(0),
	fNUA(0),
        fNtrksName("Mult"),
	//....
	fPeriod("LHC15o"),

	fListOfObjects(0),

	fMultTOFLowCut(0),
	fMultTOFHighCut(0),
	fMultCentLowCut(0),

	fTrackEfficiency(0),
	hTrackEfficiency(0),
	hTrackEfficiencyRun(0),

        fFlowRunByRunWeights(false),
        fFlowPeriodWeights(false),
        fFlowUse3Dweights(false),
        fFlowWeightsList(nullptr),

	fPhiWeight(0),
	fPhiWeightPlus(0),
	fPhiWeightMinus(0),
	hPhiWeight(0),
	hPhiWeight1D(0),
	hPhiWeight_LHC15i_part1(0),
	hPhiWeight_LHC15i_part2(0),
	hPhiWeight_LHC15j_part1(0),
	hPhiWeight_LHC15j_part2(0),
	hPhiWeight_LHC15l_part1(0),
	hPhiWeight_LHC15l_part2(0),
	hPhiWeight_LHC15l_part3(0),
	hPhiWeightPlus_LHC15i_part1(0),
	hPhiWeightPlus_LHC15i_part2(0),
	hPhiWeightPlus_LHC15j_part1(0),
	hPhiWeightPlus_LHC15j_part2(0),
	hPhiWeightPlus_LHC15l_part1(0),
	hPhiWeightPlus_LHC15l_part2(0),
	hPhiWeightPlus_LHC15l_part3(0),
	hPhiWeightMinus_LHC15i_part1(0),
	hPhiWeightMinus_LHC15i_part2(0),
	hPhiWeightMinus_LHC15j_part1(0),
	hPhiWeightMinus_LHC15j_part2(0),
	hPhiWeightMinus_LHC15l_part1(0),
	hPhiWeightMinus_LHC15l_part2(0),
	hPhiWeightMinus_LHC15l_part3(0),

	hEventCount(0),
	hMult(0),
	fVtxAfterCuts(0),
	fCentralityDis(0),
	fV0CentralityDis(0),
	hMultV0vsNtrksAfterCuts(0),
	hMultSPDvsNtrksAfterCuts(0),
	hNtrksVSmultPercentile(0),
	fCentralityV0MCL1(0),
	fCentralityV0MCL0(0),
	fCentralityCL0CL1(0),
	fMultvsCentr(0),
	fMult128vsCentr(0),
	fMultTPCvsTOF(0),
	fMultTPCvsESD(0),

	hSPDClsVsTrk(0),
	hV0C012vsTkl(0),
	hV0C012vsV0C3(0),
	hV0MOnVsOf(0),
	hSPDOnVsOf(0),

	fPhiDis1D(0),
	fPhiDis(0),
	fEtaDis(0),
	fEtaBefore(0),
	fPtDis(0),
	fPtBefore(0),
	hDCAxyBefore(0),
	hDCAzBefore(0),
	hITSclustersBefore(0),
	hChi2Before(0),
	hDCAxy(0),
	hDCAz(0),
	hITSclusters(0),
	hChi2(0),
	rand(32213) {

		// Output slot #1 writes into a TList
		DefineOutput(1, TList::Class());
		// DefineOutput(2, TList::Class());
		DefineInput(1, TFile::Class());
		DefineInput(2, TFile::Class());
	}

//_____________________________________________________________________________
AliAnalysisTaskNonlinearFlow::~AliAnalysisTaskNonlinearFlow()
{
	// Destructor
	// histograms are in the output list and deleted when the output
	// list is deleted by the TSelector dtor
	if (fListOfObjects)
		delete fListOfObjects;

}

//______________________________________________________________________________
void AliAnalysisTaskNonlinearFlow::UserCreateOutputObjects()
{

	//OpenFile(1);
	fListOfObjects = new TList();
	fListOfObjects->SetOwner();

	//..Settings for AliEventCuts:
	//..This adds QA plots to the output
	fEventCuts.AddQAplotsToList(fListOfObjects);
	//..kINT7 is set in the class as default, if I want to have kHigHMultV0 in pp, I have to switch to manual mode

	fEventCuts.SetManualMode();
	fEventCuts.fRequireTrackVertex = false; // !!
	fEventCuts.fMinVtz = -10.f;
	fEventCuts.fMaxVtz = 10.f;
	fEventCuts.fMaxResolutionSPDvertex = 0.25f;
	// Distance between track and SPD vertex < 0.2 cm
	fEventCuts.fPileUpCutMV = true;

	// range on Xaxis:
	// double xbins_tmp[] = {50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300,
        // 		1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000 };
       

        if (fNtrksName == "Mult") {
	    nn = 3000;
            for (int i = 0; i <= 3000; i++) {
                xbins[i] = i;
            }
        } else {
            nn = 10;
            for (int i = 0; i <= 10; i++) {
                xbins[i] = i * 10;
            }
        }


	hEventCount = new TH1D("hEventCount", "; centrality;;", 1, 0, 1);
	fListOfObjects->Add(hEventCount);

	hMult = new TH1F("hMult", ";number of tracks; entries", nn, xbins);
	hMult->Sumw2();
	fListOfObjects->Add(hMult);

	fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (after cuts); Vtx z [cm]; Counts", 120, -30, 30);
	fVtxAfterCuts->Sumw2();
	fListOfObjects->Add(fVtxAfterCuts);

	fCentralityDis = new TH1F("fCentralityDis", "centrality distribution; centrality; Counts", 100, 0, 100);
	fListOfObjects->Add(fCentralityDis);

	fV0CentralityDis = new TH1F("fV0CentralityDis", "centrality V0/<V0> distribution; centrality; Counts", 100, 0, 10);
	fListOfObjects->Add(fV0CentralityDis);

	hMultV0vsNtrksAfterCuts = new TH2F("hMultV0vsNtrksAfterCuts","V0 mult vs. number of tracks; V0 mult; number of tracks", 100, 0, 10, nn, xbins);
	fListOfObjects->Add(hMultV0vsNtrksAfterCuts);

	hMultSPDvsNtrksAfterCuts = new TH2F("hMultSPDvsNtrksAfterCuts","SPD mult vs. number of tracks; SPD mult; number of tracks", 100, 0, 10, nn, xbins);
	fListOfObjects->Add(hMultSPDvsNtrksAfterCuts);

	hNtrksVSmultPercentile = new TH2F("hNtrksVSmultPercentile", ";Multiplicity percentile;ITSsa tracks", 100, 0, 100, 1000, 0, 2000);
	fListOfObjects->Add(hNtrksVSmultPercentile);

	Int_t inSlotCounter=1;
	if(fNUA) {
                fFlowWeightsList = (TList*) GetInputData(inSlotCounter);
		inSlotCounter++;
	};
	if(fNUE) {
		fTrackEfficiency = (TFile*)GetInputData(inSlotCounter);
		inSlotCounter++;
	};

	// Physics profiles
	//	NL response
	InitProfile(multProfile, "");
        for (int i = 0; i < 10; i++) InitProfile(multProfile_bin[i], Form("_%d", i));

	// Post output data.
	PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskNonlinearFlow::UserExec(Option_t *)
{
	bootstrap_value = rand.Integer(10);

	//..apply physics selection
	UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
	Bool_t isTrigselected = false;
        if (fTrigger == 0) {
	    isTrigselected = fSelectMask&AliVEvent::kINT7;
	    fAliTrigger = AliVEvent::kINT7;
        } else if (fTrigger == 1) {
	    isTrigselected = fSelectMask&AliVEvent::kHighMultV0;
	    fAliTrigger = AliVEvent::kHighMultV0;
        }
	if(isTrigselected == false) return;

	//..check if I have AOD
	fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
	if(!fAOD){
		Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
		return;
	}

	if (fPeriod.EqualTo("LHC15o")) {
		if(!fEventCuts.AcceptEvent(fAOD)) { // automatic event selection for Run2
			PostData(1,fListOfObjects);
			return;
		}
	} else {
		if(fAOD->IsPileupFromSPDInMultBins() ) { return; }

		AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
		if (!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return; }

		if(!multSelection->GetThisEventIsNotPileup() || !multSelection->GetThisEventIsNotPileupInMultBins() || !multSelection->GetThisEventHasNoInconsistentVertices() || !multSelection->GetThisEventPassesTrackletVsCluster()) { return; }

		Int_t nTracksPrim = fAOD->GetPrimaryVertex()->GetNContributors();
		if(nTracksPrim < 0.5) { return; }
	}
	hEventCount->Fill("after fEventCuts", 1.);

	//..filling Vz distribution
	AliVVertex *vtx = fAOD->GetPrimaryVertex();
	float fVtxZ = vtx->GetZ();
	if(TMath::Abs(fVtxZ) > fVtxCut) return;
        NTracksCalculation(fInputEvent);
	if(TMath::Abs(fVtxZ) > fVtxCutDefault) return;
	fVtxAfterCuts->Fill(fVtxZ);

	//..standard event plots (cent. percentiles, mult-vs-percentile)

	const auto pms(static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection")));
	const auto dCentrality(pms->GetMultiplicityPercentile("V0M"));
	float fMultV0Meq = 0;
	float fMultMeanV0M = 0;
	float fMultSPD = 0;
	float fMultMeanSPD = 0;
	float centrV0 = 0;
	float cent = dCentrality;
	float centSPD = 0;
	float v0Centr = 0;
	float cl1Centr = 0;
	float cl0Centr = 0;

	fCentralityDis->Fill(centrV0);
	fV0CentralityDis->Fill(cent);


        // checking the run number for aplying weights & loading TList with weights
        //
        // if (fCurrSystFlag == 0) {
        if ( !fPeriod.EqualTo("LHC15o") ) {
		if (fNUA && !LoadWeights()) { AliFatal("Weights not loaded! at pp"); return; }
        } else {
		if (fNUA && !LoadWeightsSystematics()) { AliFatal("Weights not loaded! for LHC15o"); return; }
        }



	//..all charged particles
	if(fLS == true){
		AnalyzeAOD(fInputEvent, centrV0, cent, centSPD, fVtxZ, false);
		AnalyzeAOD(fInputEvent, centrV0, cent, centSPD, fVtxZ, true);
	} else AnalyzeAOD(fInputEvent, centrV0, cent, centSPD, fVtxZ, false);

	// Post output data.
	PostData(1, fListOfObjects);
}

void AliAnalysisTaskNonlinearFlow::NTracksCalculation(AliVEvent* aod) {
	const int nAODTracks = aod->GetNumberOfTracks();
	NtrksCounter = 0;

	//..for DCA
	double pos[3], vz, vx, vy;
	vz = aod->GetPrimaryVertex()->GetZ();
	vx = aod->GetPrimaryVertex()->GetX();
	vy = aod->GetPrimaryVertex()->GetY();

	//..LOOP OVER TRACKS........
	//........................................
	for(Int_t nt = 0; nt < nAODTracks; nt++)
	{
		AliAODTrack *aodTrk = (AliAODTrack*) fInputEvent->GetTrack(nt);

		if (!aodTrk){
			delete aodTrk;
			continue;
		}

		double pos[3];
		aodTrk->GetXYZ(pos);
		double dcaZ = 100;
		double dcaX = 100;
		double dcaY = 100;
		double dcaXY = 100;
		dcaZ = pos[2] - vz;
		dcaX = pos[0] - vx;
		dcaY = pos[1] - vy;
		dcaXY = TMath::Sqrt(dcaX*dcaX + dcaY*dcaY);

		int nClustersITS = 0;
		nClustersITS = aodTrk->GetITSNcls();
		float chi2PerClusterITS = -1;
		if(nClustersITS != 0) chi2PerClusterITS = aodTrk->GetITSchi2()/float(nClustersITS);

		if (!(aodTrk->TestFilterBit(fFilterbitDefault))) { continue; }
		if (fFilterbitDefault == 96) {
			if (TMath::Abs(dcaZ) > fDCAzDefault) continue;
		}

		if(aodTrk->Pt() < fMinPt) continue;
		if(aodTrk->Pt() > fMaxPt) continue;

		if(TMath::Abs(aodTrk->Eta()) > fEtaCut) continue;

		NtrksCounter += 1;
	} // end loop of all track
}

//________________________________________________________________________
void AliAnalysisTaskNonlinearFlow::AnalyzeAOD(AliVEvent* aod, float centrV0, float cent, float centSPD, float fVtxZ, bool fPlus)
{

	const int nAODTracks = aod->GetNumberOfTracks();

	// Init the number of tracks
	double NtrksBefore = 0;
        NtrksAfter = 0;
	NtrksAfterGap10M = 0;
	NtrksAfterGap10P = 0;
	NtrksAfterGap14M = 0;
	NtrksAfterGap14P = 0;
	NtrksAfter3subL = 0;
	NtrksAfter3subM = 0;
	NtrksAfter3subR = 0;


	//..for DCA
	double pos[3], vz, vx, vy;
	vz = aod->GetPrimaryVertex()->GetZ();
	vx = aod->GetPrimaryVertex()->GetX();
	vy = aod->GetPrimaryVertex()->GetY();

	double Qcos[20][20] = {0};
	double Qsin[20][20] = {0};
	double QcosGap10M[20][20] = {0};
	double QsinGap10M[20][20] = {0};
	double QcosGap10P[20][20] = {0};
	double QsinGap10P[20][20] = {0};
	double QcosGap14M[20][20] = {0};
	double QsinGap14M[20][20] = {0};
	double QcosGap14P[20][20] = {0};
	double QsinGap14P[20][20] = {0};
	double QcosSubLeft[20][20] = {0};
	double QsinSubLeft[20][20] = {0};
	double QcosSubMiddle[20][20] = {0};
	double QsinSubMiddle[20][20] = {0};
	double QcosSubRight[20][20] = {0};
	double QsinSubRight[20][20] = {0};



	int run = GetRunPart(fInputEvent->GetRunNumber());
	double runNumber = fInputEvent->GetRunNumber();

	//..LOOP OVER TRACKS........
	//........................................
	for(Int_t nt = 0; nt < nAODTracks; nt++)
	{

		AliAODTrack *aodTrk = (AliAODTrack*) fInputEvent->GetTrack(nt);

		if (!aodTrk){
			delete aodTrk;
			continue;
		}

		aodTrk->GetXYZ(pos);
		double dcaZ = 100;
		double dcaX = 100;
		double dcaY = 100;
		double dcaXY = 100;
		dcaZ = pos[2] - vz;
		dcaX = pos[0] - vx;
		dcaY = pos[1] - vy;
		dcaXY = TMath::Sqrt(dcaX*dcaX + dcaY*dcaY);

                // nClustersITS cut
		int nClustersITS = 0;
		nClustersITS = aodTrk->GetITSNcls();
		float chi2PerClusterITS = -1;
		if(nClustersITS != 0) chi2PerClusterITS = aodTrk->GetITSchi2()/float(nClustersITS);

                // nClustersTPC cut
                int nClustersTPC = 0;
		nClustersTPC = aodTrk->GetTPCNcls();
		if (nClustersTPC < fTPCclusters) { continue; }
		float chi2PerClusterTPC = -1;
		if (nClustersTPC != 0) chi2PerClusterTPC = aodTrk->GetTPCchi2()/float(nClustersTPC);
		if (chi2PerClusterTPC > fChi2PerTPCcluster) { continue; }

		if (!(aodTrk->TestFilterBit(fFilterbit))) { continue; }
		if (fFilterbit == 96) {
			if (TMath::Abs(dcaZ) > fDCAz) continue;
		}
                if(fUseDCAzCut && dcaZ > fDCAz) { continue; }
                if(fUseDCAxyCut && dcaXY > (0.0105+0.0350/pow(aodTrk->Pt(),1.1))*fDCAxy) { continue; }

		if(aodTrk->Pt() < fMinPt) continue;
		if(aodTrk->Pt() > fMaxPt) continue;

		if(TMath::Abs(aodTrk->Eta()) > fEtaCut) continue;

		NtrksAfter += 1;

		//..get phi-weight for NUA correction
		double weight = 1;
		if(fNUA == 1) {
                        // if (fCurrSystFlag == 0) {
			if ( !fPeriod.EqualTo("LHC15o") ) {
				weight = GetFlowWeight(aodTrk, fVtxZ, kRefs);
			} else {
				weight = GetFlowWeightSystematics(aodTrk, fVtxZ, kRefs);
			}
		}
		double weightPt = 1;
		if(fNUE == 1) weightPt = GetPtWeight(aodTrk->Pt(), aodTrk->Eta(), fVtxZ, runNumber);
		NtrksBefore += weightPt;

		//..calculate Q-vectors
		//..no eta gap
		for(int iharm=0; iharm<20; iharm++)
		{
			for(int ipow=0; ipow<20; ipow++)
			{
				Qcos[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
				Qsin[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
			}
		}

		//..Gap > 1.0
		if(aodTrk->Eta() < -0.5)
		{
			NtrksAfterGap10M++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap10M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
					QsinGap10M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
				}
			}
		}
		if(aodTrk->Eta() > 0.5)
		{
			NtrksAfterGap10P++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap10P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
					QsinGap10P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
				}
			}
		}

		//..Gap > 1.4
		if(aodTrk->Eta() < -0.7)
		{
			NtrksAfterGap14M++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap14M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
					QsinGap14M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
				}
			}
		}
		if(aodTrk->Eta() > 0.7)
		{
			NtrksAfterGap14P++;
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					QcosGap14P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
					QsinGap14P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
				}
			}
		}
		//..3-subevent method
		if(aodTrk->Eta() < -0.4)
		{//..left part
			NtrksAfter3subL += 1;
			for(int iharm=0; iharm<10; iharm++)
			{
				for(int ipow=0; ipow<10; ipow++)
				{
					QcosSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
					QsinSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
				}
			}
		}
		if(aodTrk->Eta() >= -0.4 && aodTrk->Eta() <= 0.4)
		{//..middle part
			NtrksAfter3subM += 1;
			for(int iharm=0; iharm<10; iharm++)
			{
				for(int ipow=0; ipow<10; ipow++)
				{
					QcosSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
					QsinSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
				}
			}
		}
		if(aodTrk->Eta() > 0.4)
		{//..right part
			NtrksAfter3subR += 1;
			for(int iharm=0; iharm<10; iharm++)
			{
				for(int ipow=0; ipow<10; ipow++)
				{
					QcosSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
					QsinSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
				}
			}
		}


	} // end loop of all track

	//............................
	//..GENERIC FRAMEWORK RP
	//............................

	//..calculate Q-vector for each harmonics n and power p
        correlator.FillQVector(correlator.Qvector, Qcos, Qsin); 
        correlator.FillQVector(correlator.Qvector10M, QcosGap10M, QsinGap10M); 
        correlator.FillQVector(correlator.Qvector10P, QcosGap10P, QsinGap10P); 
        correlator.FillQVector(correlator.Qvector14M, QcosGap14M, QsinGap14M); 
        correlator.FillQVector(correlator.Qvector14P, QcosGap14P, QsinGap14P); 
        correlator.FillQVector(correlator.QvectorSubLeft, QcosSubLeft, QsinSubLeft); 
        correlator.FillQVector(correlator.QvectorSubRight, QcosSubRight, QsinSubRight); 
        correlator.FillQVector(correlator.QvectorSubMiddle, QcosSubMiddle, QsinSubMiddle); 

	hMult->Fill(NtrksCounter);

	// CalculateProfile(centProfile, cent);
	if (fNtrksName == "Mult") {
	    CalculateProfile(multProfile, NtrksCounter);
	    CalculateProfile(multProfile_bin[bootstrap_value], NtrksCounter);
        } else {
	    CalculateProfile(multProfile, cent);
	    CalculateProfile(multProfile_bin[bootstrap_value], cent);
        }
	// CalculateProfile(centProfile_bin[bootstrap_value], cent);

}

//____________________________________________________________________
//	END OF MAIN PROGRAM
//____________________________________________________________________
Bool_t AliAnalysisTaskNonlinearFlow::IsGoodPSEvent(AliVEvent* event)
{

	IsSPDClusterVsTrackletBG(event, true);
	IsV0C012vsTklBG(event, true);
	IsV0Casym(event, true);
	IsV0MOnVsOfPileup(event, true);
	IsSPDOnVsOfPileup(event, true);

	bool is = true;

	if(IsSPDClusterVsTrackletBG(event, false)) is = false;
	if(IsV0C012vsTklBG(event, false)) is = false;
	if(IsV0Casym(event, false)) is = false;
	if(IsV0MOnVsOfPileup(event, false)) is = false;
	if(IsSPDOnVsOfPileup(event, false)) is = false;
	if(IsV0PFPileup(event)) is = false;

	return is;

}
//-----------------------------------------------------------------------------
Bool_t AliAnalysisTaskNonlinearFlow::IsSPDClusterVsTrackletBG(const AliVEvent* event, bool fillHist){
	// rejects BG based on SPD tracklets vs. clusters correlation
	// returns true if the event is BG
	const AliVMultiplicity* mult = event->GetMultiplicity();

	Int_t nTkl = mult->GetNumberOfTracklets();
	Int_t nCls = event->GetNumberOfITSClusters(0) + event->GetNumberOfITSClusters(1);

	if(fillHist == true) hSPDClsVsTrk->Fill(nTkl, nCls);

	return nCls > 65 + nTkl*4;
}
//-------------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskNonlinearFlow::IsV0C012vsTklBG(const AliVEvent* event, bool fillHist){
	// rejects BG based on V0C012 vs tracklet correlation
	// returns true if the event is BG
	const AliVMultiplicity* mult = event->GetMultiplicity();
	AliVVZERO* vzero = event->GetVZEROData();

	Float_t nTkl       = mult->GetNumberOfTracklets();
	Float_t multV0C012 = vzero->GetMTotV0C() - vzero->GetMRingV0C(3);

	if(fillHist == true) hV0C012vsTkl->Fill(nTkl, multV0C012);

	return nTkl < 6 && multV0C012 > 150 + nTkl*20;
}
//-------------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskNonlinearFlow::IsV0Casym(const AliVEvent* event, bool fillHist){
	// rehect BG based on V0C012 vs. V0C3 mult. correlation
	// returns true if the event is BG
	AliVVZERO* vzero = event->GetVZEROData();

	Float_t multV0C012 = vzero->GetMRingV0C(0)+vzero->GetMRingV0C(1)+vzero->GetMRingV0C(2);
	Float_t multV0C3   = vzero->GetMRingV0C(3);

	if(fillHist == true) hV0C012vsV0C3->Fill(multV0C012, multV0C3);

	return (multV0C3 < -25 + 0.15*multV0C012);
}
//-------------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskNonlinearFlow::IsV0MOnVsOfPileup(const AliVEvent* event, bool fillHist){
	// rejects pileup based on V0M online vs offline correlation
	// return true if the event is pileup
	AliVVZERO* vzero = event->GetVZEROData();

	// V0A0 excluded from online V0A charge sum => excluding also from offline sum for consistency
	Float_t on = vzero->GetTriggerChargeA() + vzero->GetTriggerChargeC();
	Float_t of = vzero->GetMTotV0A() - vzero->GetMRingV0A(0) + vzero->GetMTotV0C();

	if(fillHist == true) hV0MOnVsOf->Fill(of, on);

	return (on < -145 + 7.2*of);
}
//-------------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskNonlinearFlow::IsSPDOnVsOfPileup(const AliVEvent* event, bool fillHist){
	// rejects pileup based on SPD online vs. offline correlation
	// returns true if the event is pileup
	AliVMultiplicity* mult = event->GetMultiplicity();
	TBits onMap = mult->GetFastOrFiredChips();
	TBits ofMap = mult->GetFiredChipMap();

	Int_t on = onMap.CountBits(0);
	Int_t of = ofMap.CountBits(0);

	if(fillHist == true) hSPDOnVsOf->Fill(of, on);

	return (on < -4.16 + 0.84*of);
}
//-------------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskNonlinearFlow::IsV0PFPileup(const AliVEvent* event){
	// return true if the event is pileup

	int fVIRBBAflags = 10;
	int fVIRBBCflags = 10;
	int fVIRBGAflags = 33;
	int fVIRBGCflags = 33;

	AliVVZERO* vzero = event->GetVZEROData();

	Bool_t vir[21] = {0};
	UChar_t bcMod4 = event->GetBunchCrossNumber()%4;

	for (Int_t bc=0;bc<=20;bc++) {
		UChar_t nBBA=0;
		UChar_t nBBC=0;
		UChar_t nBGA=0;
		UChar_t nBGC=0;
		if (fVIRBBAflags<33) for (Int_t i=0;i<32;i++) nBBA+=vzero->GetPFBBFlag(i+32,bc);
		if (fVIRBBCflags<33) for (Int_t i=0;i<32;i++) nBBC+=vzero->GetPFBBFlag(i   ,bc);
		if (fVIRBGAflags<33) for (Int_t i=0;i<32;i++) nBGA+=vzero->GetPFBGFlag(i+32,bc);
		if (fVIRBGCflags<33) for (Int_t i=0;i<32;i++) nBGC+=vzero->GetPFBGFlag(i   ,bc);
		vir[bc] |= nBBA>=fVIRBBAflags;
		vir[bc] |= nBBC>=fVIRBBCflags;
		vir[bc] |= nBGA>=fVIRBGAflags;
		vir[bc] |= nBGC>=fVIRBGCflags;
	}

	// clock index is counting from future to past
	Int_t bcMin = 10 - 7 + bcMod4;
	Int_t bcMax = 10 + 4 + bcMod4;
	for (Int_t bc=bcMin;bc<=bcMax;bc++) {
		if (bc==10) continue; // skip current bc
		if (bc < 0) continue;
		if (bc >20) continue;
		if (vir[bc]) return kTRUE;
	}

	return kFALSE;
}
//____________________________________________________________________
int AliAnalysisTaskNonlinearFlow::GetRunPart(int run)
{

	int fRun = 0;

	//..LHC15i, part 1
	if(run == 236137 || run == 236138 || run == 236150 || run == 236151 || run == 236153
			|| run == 236158 || run == 236159 || run == 236163 || run == 236164 || run == 236203
			|| run == 236204 || run == 236222 || run == 236227 || run == 236234 || run == 236238
			|| run == 236240 || run == 236242 || run == 236244 || run == 236246 || run == 236248)
		fRun = 1;
	//..LHC15i, part2
	if(run == 236281 || run == 236284 || run == 236285 || run == 236331 || run == 236334
			|| run == 236337 || run == 236348 || run == 236349 || run == 236352 || run == 236353
			|| run == 236354 || run == 236356 || run == 236357 || run == 236359 || run == 236360
			|| run == 236386 || run == 236389 || run == 236393 || run == 236395 || run == 236397
			|| run == 236441 || run == 236443 || run == 236444 || run == 236446 || run == 236453
			|| run == 236459 || run == 236462 || run == 236541 || run == 236554 || run == 236556
			|| run == 236558 || run == 236562 || run == 236563 || run == 236564 || run == 236565
			|| run == 236569)
		fRun = 2;

	//..LHC15j, part1
	if(run == 238091 || run == 238097 || run == 238129 || run == 238131 || run == 238132
			|| run == 238133 || run == 238136 || run == 238139 || run == 238140 || run == 238142
			|| run == 238144 || run == 238145 || run == 238147 || run == 238148 || run == 238159
			|| run == 238160 || run == 238164 || run == 238170 || run == 238570)
		fRun = 3;
	//..LHC15j, part2
	if(run == 237029 || run == 237406 || run == 237408 || run == 237409 || run == 237507
			|| run == 237512 || run == 237515 || run == 237645 || run == 237670 || run == 237671
			|| run == 237675 || run == 237676 || run == 237678 || run == 237681 || run == 237684
			|| run == 237691 || run == 237698 || run == 237699 || run == 237705 || run == 237706
			|| run == 237707 || run == 237708 || run == 237710 || run == 237711 || run == 237713
			|| run == 237765 || run == 237768 || run == 237777 || run == 237779 || run == 237780
			|| run == 237782 || run == 237787 || run == 237789 || run == 237790 || run == 237791
			|| run == 237793 || run == 237795 || run == 237796 || run == 237806 || run == 237842
			|| run == 237844 || run == 237845 || run == 237847 || run == 237945 || run == 237948
			|| run == 237969 || run == 237978 || run == 237982 || run == 237983 || run == 238073
			|| run == 238176 || run == 238179 || run == 238184 || run == 238185 || run == 238187
			|| run == 238395 || run == 238451 || run == 238454 || run == 238455 || run == 238456
			|| run == 238457 || run == 238458 || run == 238459 || run == 238460 || run == 238472
			|| run == 238474 || run == 238604 || run == 238606 || run == 238607 || run == 238610
			|| run == 238614 || run == 238621)
		fRun = 4;

	//..LHC15l, part1
	if(run == 241257 || run == 241261 || run == 241263 || run == 241267 || run == 241268
			|| run == 241269 || run == 241281 || run == 241288 || run == 241295 || run == 241296)
		fRun = 5;
	//..LHC15l, part2
	if(run == 240069) fRun = 6;
	//..LHC15l, part3
	if(run == 239319 || run == 239324 || run == 239518 || run == 239519 || run == 240183
			|| run == 240194 || run == 240196 || run == 240201 || run == 240204 || run == 240212
			|| run == 240220 || run == 240241 || run == 240250 || run == 240256 || run == 240262
			|| run == 240263 || run == 240265 || run == 240271 || run == 240274 || run == 240293
			|| run == 240303 || run == 240312 || run == 240376 || run == 240380 || run == 240381
			|| run == 240382 || run == 240385 || run == 240392 || run == 240394 || run == 240404
			|| run == 240411 || run == 240443 || run == 240444 || run == 240447 || run == 240450
			|| run == 240452 || run == 240610 || run == 240612 || run == 240845 || run == 240854
			|| run == 240860 || run == 240864 || run == 240872 || run == 240874 || run == 240875
			|| run == 240880 || run == 241001 || run == 241010 || run == 241014 || run == 241021
			|| run == 241032 || run == 241043 || run == 241047 || run == 241050 || run == 241054
			|| run == 241055 || run == 241056 || run == 241057 || run == 241062 || run == 241069
			|| run == 241075 || run == 241141 || run == 241144 || run == 241354 || run == 241360
			|| run == 241361 || run == 241393 || run == 241396 || run == 241407 || run == 241412)
		fRun = 7;


	return fRun;

}
//____________________________________________________________________
double AliAnalysisTaskNonlinearFlow::GetPtWeight(double pt, double eta, float vz, double runNumber)
{
        return 1;
	hTrackEfficiencyRun = (TH3F*)fTrackEfficiency->Get(Form("eff_LHC15o_HIJING_%.0lf", runNumber));
	double binPt = hTrackEfficiencyRun->GetXaxis()->FindBin(pt);
	double binEta = hTrackEfficiencyRun->GetYaxis()->FindBin(eta);
	double binVz = hTrackEfficiencyRun->GetZaxis()->FindBin(vz);
	//..take into account error on efficiency: randomly get number from gaussian distribution of eff. where width = error
	double eff = hTrackEfficiencyRun->GetBinContent(binPt, binEta, binVz);
	double error = hTrackEfficiencyRun->GetBinError(binPt, binEta, binVz);

	double weight = 1;
	if((eff < 0.03) || ((error/eff) > 0.1)) weight = 1;
	else{
		TRandom3 r(0);
		double efficiency = 0;
		efficiency = r.Gaus(eff, error);
		weight = 1./efficiency; //..taking into account errors
		//weight = 1./eff;
	}

	return weight;

}
//____________________________________________________________________
double AliAnalysisTaskNonlinearFlow::GetWeight(double phi, double eta, double pt, int fRun, bool fPlus, double vz, double runNumber)
{
	TList* weights_list = dynamic_cast<TList*>(fPhiWeight);
        // cout << "weights_list" << weights_list << endl;
        // weights_list->ls();
        
	TList* averaged_list = dynamic_cast<TList*>(weights_list->FindObject("averaged"));
        // cout << "averaged_list" << averaged_list << endl;
	TH2D* hPhiWeightRun = dynamic_cast<TH2D*>(averaged_list->FindObject("Charged"));
        // cout << "hist_list" << hPhiWeightRun << endl;

	double weight = hPhiWeightRun->GetBinContent(hPhiWeightRun->GetXaxis()->FindBin(phi),
			hPhiWeightRun->GetYaxis()->FindBin(eta));
			// , hPhiWeightRun->GetZaxis()->FindBin(vz));
	return weight;
}

const char* AliAnalysisTaskNonlinearFlow::GetSpeciesName(const PartSpecies species) const
{
  const char* name;

  switch(species) {
    case kRefs: name = "Refs"; break;
    case kCharged: name = "Charged"; break;
    case kPion: name = "Pion"; break;
    case kKaon: name = "Kaon"; break;
    case kProton: name = "Proton"; break;
    case kCharUnidentified: name = "UnidentifiedCharged"; break;
    case kK0s: name = "K0s"; break;
    case kLambda: name = "Lambda"; break;
    case kPhi: name = "Phi"; break;
    default: name = "Unknown";
  }

  return name;
}

Bool_t AliAnalysisTaskNonlinearFlow::LoadWeightsSystematics()
{
	if(fCurrSystFlag == 0) fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i",fAOD->GetRunNumber()));
        else fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i_SystFlag%i_",fAOD->GetRunNumber(), fCurrSystFlag));
        if(!fWeightsSystematics)
        {
            printf("Weights could not be found in list!\n");
            return kFALSE;
        }
        fWeightsSystematics->CreateNUA();
        return kTRUE;
}

Double_t AliAnalysisTaskNonlinearFlow::GetFlowWeightSystematics(const AliVParticle* track, double fVtxZ, const PartSpecies species)
{

    double dPhi = track->Phi();
    double dEta = track->Eta();
    double dVz = fVtxZ;
    double dWeight = 1.0;
    dWeight = fWeightsSystematics->GetNUA(dPhi, dEta, dVz);
    return dWeight;
}

Bool_t AliAnalysisTaskNonlinearFlow::LoadWeights()
{
  // (Re-) Loading of flow vector weights
  // ***************************************************************************
  if(!fFlowWeightsList) { AliError("Flow weights list not found! Terminating!"); return kFALSE; }

  TList* listFlowWeights = nullptr;
  
  TString fFlowWeightsTag = "";
  if(!fFlowWeightsTag.IsNull()) {
      // using weights Tag if provided (systematics)
      listFlowWeights = (TList*) fFlowWeightsList->FindObject(fFlowWeightsTag.Data());
      if(!listFlowWeights) { AliError(Form("TList with tag '%s' not found!",fFlowWeightsTag.Data())); fFlowWeightsList->ls(); return kFALSE; }
  } else {
      if(!fFlowRunByRunWeights && !fFlowPeriodWeights) {
          // loading run-averaged weights
          listFlowWeights = (TList*) fFlowWeightsList->FindObject("averaged");
          if(!listFlowWeights) { AliError("TList with flow run-averaged weights not found."); fFlowWeightsList->ls(); return kFALSE; }
      } else if(fFlowPeriodWeights){
        // loading period-specific weights
        listFlowWeights = (TList*) fFlowWeightsList->FindObject(ReturnPPperiod(fAOD->GetRunNumber()));
        if(!listFlowWeights) { AliError("Loading period weights failed!"); fFlowWeightsList->ls(); return kFALSE; }
      }
      else {
          // loading run-specific weights
          listFlowWeights = (TList*) fFlowWeightsList->FindObject(Form("%d",fAOD->GetRunNumber()));

          if(!listFlowWeights) {
              // run-specific weights not found for this run; loading run-averaged instead
              AliWarning(Form("TList with flow weights (run %d) not found. Using run-averaged weights instead (as a back-up)", fAOD->GetRunNumber()));
              listFlowWeights = (TList*) fFlowWeightsList->FindObject("averaged");
              if(!listFlowWeights) { AliError("Loading run-averaged weights failed!"); fFlowWeightsList->ls(); return kFALSE; }
          }
      }
  }


  for(Int_t iSpec(0); iSpec <= kRefs; ++iSpec) {
    if(fFlowUse3Dweights) {
      fh3Weights[iSpec] = (TH3D*) listFlowWeights->FindObject(Form("%s3D",GetSpeciesName(PartSpecies(iSpec))));
      if(!fh3Weights[iSpec]) { AliError(Form("Weight 3D (%s) not found",GetSpeciesName(PartSpecies(iSpec)))); return kFALSE; }
    } else {
      fh2Weights[iSpec] = (TH2D*) listFlowWeights->FindObject(GetSpeciesName(PartSpecies(iSpec)));
      if(!fh2Weights[iSpec]) { AliError(Form("Weight 2D (%s) not found",GetSpeciesName(PartSpecies(iSpec)))); return kFALSE; }
    }
  }

  return kTRUE;
}

Double_t AliAnalysisTaskNonlinearFlow::GetFlowWeight(const AliVParticle* track, double fVtxZ, const PartSpecies species)
{
  // if not applying for reconstructed
  // if(!fFlowWeightsApplyForReco && HasMass(species)) { return 1.0; }

  Double_t dWeight = 1.0;
  if(fFlowUse3Dweights) {
    Int_t iBin = fh3Weights[species]->FindFixBin(track->Phi(),track->Eta(),fVtxZ);
    dWeight = fh3Weights[species]->GetBinContent(iBin);
  } else {
    Int_t iBin = fh2Weights[species]->FindFixBin(track->Phi(),track->Eta());
    dWeight = fh2Weights[species]->GetBinContent(iBin);
  }

  if(dWeight <= 0.0) { dWeight = 1.0; }
  return dWeight;
}

void AliAnalysisTaskNonlinearFlow::InitProfile(PhysicsProfile& multProfile, TString label) {

	for(int h=0; h<6; h++)
	{
		multProfile.fChcn2[h] = new TProfile(Form("fChc%d{2}%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
		multProfile.fChcn2[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2[h]);

		multProfile.fChcn2_Gap10[h] = new TProfile(Form("fChc%d{2}_Gap10%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
		multProfile.fChcn2_Gap10[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_Gap10[h]);

		multProfile.fChcn2_Gap14[h] = new TProfile(Form("fChc%d{2}_Gap14%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
		multProfile.fChcn2_Gap14[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_Gap14[h]);

		multProfile.fChcn2_3subLM[h] = new TProfile(Form("fChc%d{2}_3subLM%s", h+2, label.Data()), "<<2>> Re for 3-subevent method, left+middle; # of tracks", nn, xbins);
		multProfile.fChcn2_3subLM[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_3subLM[h]);

		multProfile.fChcn2_3subRM[h] = new TProfile(Form("fChc%d{2}_3subRM%s", h+2, label.Data()), "<<2>> Re for 3-subevent method, right+middle; # of tracks", nn, xbins);
		multProfile.fChcn2_3subRM[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_3subRM[h]);

		multProfile.fChcn2_3subLR[h] = new TProfile(Form("fChc%d{2}_3subLR%s", h+2, label.Data()), "<<2>> Re for 3-subevent method, left+right; # of tracks", nn, xbins);
		multProfile.fChcn2_3subLR[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn2_3subLR[h]);

		multProfile.fChcn4[h] = new TProfile(Form("fChc%d{4}%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
		multProfile.fChcn4[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4[h]);

		multProfile.fChcn4_Gap10[h] = new TProfile(Form("fChc%d{4}_Gap10%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
		multProfile.fChcn4_Gap10[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_Gap10[h]);

		multProfile.fChcn4_3subLLMR[h] = new TProfile(Form("fChc%d{4}_3subLLMR%s", h+2, label.Data()), "<<4>> 3-subevent method; # of tracks", nn, xbins);
		multProfile.fChcn4_3subLLMR[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_3subLLMR[h]);

		multProfile.fChcn4_3subRRML[h] = new TProfile(Form("fChc%d{4}_3subRRML%s", h+2, label.Data()), "<<4>> 3-subevent method; # of tracks", nn, xbins);
		multProfile.fChcn4_3subRRML[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_3subRRML[h]);

		multProfile.fChcn4_3subMMLR[h] = new TProfile(Form("fChc%d{4}_3subMMLR%s", h+2, label.Data()), "<<4>> 3-subevent method; # of tracks", nn, xbins);
		multProfile.fChcn4_3subMMLR[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_3subMMLR[h]);

		multProfile.fChcn4_3subGap2[h] = new TProfile(Form("fChc%d{4}_3subGap2%s", h+2, label.Data()), "<<4>> 3-subevent method; # of tracks", nn, xbins);
		multProfile.fChcn4_3subGap2[h]->Sumw2();
		fListOfObjects->Add(multProfile.fChcn4_3subGap2[h]);

	} // harmonics


	multProfile.fChc422 = new TProfile(Form("fChc422%s", label.Data()), "", nn, xbins);
	multProfile.fChc422->Sumw2();
	fListOfObjects->Add(multProfile.fChc422);

	multProfile.fChc532 = new TProfile(Form("fChc532%s", label.Data()), "", nn, xbins);
	multProfile.fChc532->Sumw2();
	fListOfObjects->Add(multProfile.fChc532);

	multProfile.fChc422_Gap10A = new TProfile(Form("fChc422_Gap10A%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap10A->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap10A);

	multProfile.fChc422_Gap10B = new TProfile(Form("fChc422_Gap10B%s", label.Data()), "", nn, xbins);
	multProfile.fChc422_Gap10B->Sumw2();
	fListOfObjects->Add(multProfile.fChc422_Gap10B);

	multProfile.fChc532_Gap10A = new TProfile(Form("fChc532_Gap10A%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap10A->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap10A);

	multProfile.fChc532_Gap10B = new TProfile(Form("fChc532_Gap10B%s", label.Data()), "", nn, xbins);
	multProfile.fChc532_Gap10B->Sumw2();
	fListOfObjects->Add(multProfile.fChc532_Gap10B);

	// SC(n,m): SC(3,2)
	multProfile.fChsc3232 = new TProfile(Form("fChsc3232%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232);

	multProfile.fChsc3232_Gap10 = new TProfile(Form("fChsc3232_Gap10%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_Gap10->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_Gap10);

	multProfile.fChsc3232_3subMMLRA = new TProfile(Form("fChsc3232_3subMMLRA%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_3subMMLRA->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_3subMMLRA);

	multProfile.fChsc3232_3subMMLRB = new TProfile(Form("fChsc3232_3subMMLRB%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_3subMMLRB->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_3subMMLRB);

	multProfile.fChsc3232_3subLLMRA = new TProfile(Form("fChsc3232_3subLLMRA%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_3subLLMRA->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_3subLLMRA);

	multProfile.fChsc3232_3subLLMRB = new TProfile(Form("fChsc3232_3subLLMRB%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_3subLLMRB->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_3subLLMRB);

	multProfile.fChsc3232_3subRRMLA = new TProfile(Form("fChsc3232_3subRRMLA%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_3subRRMLA->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_3subRRMLA);

	multProfile.fChsc3232_3subRRMLB = new TProfile(Form("fChsc3232_3subRRMLB%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc3232_3subRRMLB->Sumw2();
	fListOfObjects->Add(multProfile.fChsc3232_3subRRMLB);

	// SC(n,m): SC(4,2)
	multProfile.fChsc4242 = new TProfile(Form("fChsc4242%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242);

	multProfile.fChsc4242_Gap10 = new TProfile(Form("fChsc4242_Gap10%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_Gap10->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_Gap10);

	multProfile.fChsc4242_3subMMLRA = new TProfile(Form("fChsc4242_3subMMLRA%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_3subMMLRA->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_3subMMLRA);

	multProfile.fChsc4242_3subMMLRB = new TProfile(Form("fChsc4242_3subMMLRB%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_3subMMLRB->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_3subMMLRB);

	multProfile.fChsc4242_3subLLMRA = new TProfile(Form("fChsc4242_3subLLMRA%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_3subLLMRA->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_3subLLMRA);

	multProfile.fChsc4242_3subLLMRB = new TProfile(Form("fChsc4242_3subLLMRB%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_3subLLMRB->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_3subLLMRB);

	multProfile.fChsc4242_3subRRMLA = new TProfile(Form("fChsc4242_3subRRMLA%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_3subRRMLA->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_3subRRMLA);

	multProfile.fChsc4242_3subRRMLB = new TProfile(Form("fChsc4242_3subRRMLB%s", label.Data()), "# of tracks", nn, xbins);
	multProfile.fChsc4242_3subRRMLB->Sumw2();
	fListOfObjects->Add(multProfile.fChsc4242_3subRRMLB);

}

void AliAnalysisTaskNonlinearFlow::CalculateProfile(PhysicsProfile& profile, double Ntrks) {
	//..calculate 2-particle correlations
	//..................................
	double Dn2 = correlator.Two(0, 0).Re();
	double Dn2Gap10 = correlator.TwoGap10(0, 0).Re();
	double Dn2Gap14 = correlator.TwoGap14(0, 0).Re();
	double Dn2_3subLM = correlator.Two_3SubLM(0, 0).Re();
	double Dn2_3subRM = correlator.Two_3SubRM(0, 0).Re();
	double Dn2_3subLR = correlator.Two_3SubLR(0, 0).Re();

	if(NtrksAfter > 1 && Dn2 != 0) {
		//..v2{2} = <cos2(phi1 - phi2)>
		TComplex v22 = correlator.Two(2, -2);
		double v22Re = v22.Re()/Dn2;
		profile.fChcn2[0]->Fill(Ntrks, v22Re, Dn2);

		//..v3{2} = <cos3(phi1 - phi2)>
		TComplex v32 = correlator.Two(3, -3);
		double v32Re = v32.Re()/Dn2;
		profile.fChcn2[1]->Fill(Ntrks, v32Re, Dn2);

		//..v4{2} = <cos4(phi1 - phi2)>
		TComplex v42 = correlator.Two(4, -4);
		double v42Re = v42.Re()/Dn2;
		profile.fChcn2[2]->Fill(Ntrks, v42Re, Dn2);

	}

	if(NtrksAfterGap10M > 0 && NtrksAfterGap10P > 0 && Dn2Gap10 != 0)
	{
		//..v2{2} with eta Gap > 1.0
		TComplex v22Gap10 = correlator.TwoGap10(2, -2);
		double v22ReGap10 = v22Gap10.Re()/Dn2Gap10;
		profile.fChcn2_Gap10[0]->Fill(Ntrks, v22ReGap10, Dn2Gap10);

		//..v3{2} with eta Gap > 1.0
		TComplex v32Gap10 = correlator.TwoGap10(3, -3);
		double v32ReGap10 = v32Gap10.Re()/Dn2Gap10;
		profile.fChcn2_Gap10[1]->Fill(Ntrks, v32ReGap10, Dn2Gap10);

		//..v4{2} with eta Gap > 1.0
		TComplex v42Gap10 = correlator.TwoGap10(4, -4);
		double v42ReGap10 = v42Gap10.Re()/Dn2Gap10;
		profile.fChcn2_Gap10[2]->Fill(Ntrks, v42ReGap10, Dn2Gap10);
	}

	if(NtrksAfterGap14M > 0 && NtrksAfterGap14P > 0 && Dn2Gap14 != 0)
	{
		//..v2{2} with eta Gap > 1.4
		TComplex v22Gap14 = correlator.TwoGap14(2, -2);
		double v22ReGap14 = v22Gap14.Re()/Dn2Gap14;
		profile.fChcn2_Gap14[0]->Fill(Ntrks, v22ReGap14, Dn2Gap14);

		//..v3{2} with eta Gap > 1.4
		TComplex v32Gap14 = correlator.TwoGap14(3, -3);
		double v32ReGap14 = v32Gap14.Re()/Dn2Gap14;
		profile.fChcn2_Gap14[1]->Fill(Ntrks, v32ReGap14, Dn2Gap14);

		//..v4{2} with eta Gap > 1.4
		TComplex v42Gap14 = correlator.TwoGap14(4, -4);
		double v42ReGap14 = v42Gap14.Re()/Dn2Gap14;
		profile.fChcn2_Gap14[2]->Fill(Ntrks, v42ReGap14, Dn2Gap14);
	}

	//..for 3-subevent method, Gap0
	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 0 && Dn2_3subLM != 0)
	{//..left+middle
		TComplex v22_3subLM = correlator.Two_3SubLM(2, -2);
		double v22Re_3subLM = v22_3subLM.Re()/Dn2_3subLM;
		profile.fChcn2_3subLM[0]->Fill(Ntrks, v22Re_3subLM, Dn2_3subLM);

		TComplex v32_3subLM = correlator.Two_3SubLM(3, -3);
		double v32Re_3subLM = v32_3subLM.Re()/Dn2_3subLM;
		profile.fChcn2_3subLM[1]->Fill(Ntrks, v32Re_3subLM, Dn2_3subLM);

		TComplex v42_3subLM = correlator.Two_3SubLM(4, -4);
		double v42Re_3subLM = v42_3subLM.Re()/Dn2_3subLM;
		profile.fChcn2_3subLM[2]->Fill(Ntrks, v42Re_3subLM, Dn2_3subLM);
	}

	if(NtrksAfter3subM > 0 && NtrksAfter3subR > 0 && Dn2_3subRM != 0)
	{//..right+middle
		TComplex v22_3subRM = correlator.Two_3SubRM(2, -2);
		double v22Re_3subRM = v22_3subRM.Re()/Dn2_3subRM;
		profile.fChcn2_3subRM[0]->Fill(Ntrks, v22Re_3subRM, Dn2_3subRM);

		TComplex v32_3subRM = correlator.Two_3SubRM(3, -3);
		double v32Re_3subRM = v32_3subRM.Re()/Dn2_3subRM;
		profile.fChcn2_3subRM[1]->Fill(Ntrks, v32Re_3subRM, Dn2_3subRM);

		TComplex v42_3subRM = correlator.Two_3SubRM(4, -4);
		double v42Re_3subRM = v42_3subRM.Re()/Dn2_3subRM;
		profile.fChcn2_3subRM[2]->Fill(Ntrks, v42Re_3subRM, Dn2_3subRM);
	}

	if(NtrksAfter3subL > 0 && NtrksAfter3subR > 0 && Dn2_3subLR != 0)
	{//..right+middle
		TComplex v22_3subLR = correlator.Two_3SubLR(2, -2);
		double v22Re_3subLR = v22_3subLR.Re()/Dn2_3subLR;
		profile.fChcn2_3subLR[0]->Fill(Ntrks, v22Re_3subLR, Dn2_3subLR);

		TComplex v32_3subLR = correlator.Two_3SubLR(3, -3);
		double v32Re_3subLR = v32_3subLR.Re()/Dn2_3subLR;
		profile.fChcn2_3subLR[1]->Fill(Ntrks, v32Re_3subLR, Dn2_3subLR);

		TComplex v42_3subLR = correlator.Two_3SubLR(4, -4);
		double v42Re_3subLR = v42_3subLR.Re()/Dn2_3subLR;
		profile.fChcn2_3subLR[2]->Fill(Ntrks, v42Re_3subLR, Dn2_3subLR);
	}

	//..calculate 3-particle correlations
	//................................
	double Dn3 = correlator.Three(0, 0, 0).Re();
	double Dn3Gap10A = correlator.ThreeGap10A(0, 0, 0).Re();
	double Dn3Gap10B = correlator.ThreeGap10B(0, 0, 0).Re();

	if(NtrksAfter > 2 && Dn3 != 0)
	{
		//..v4{psi2}
		TComplex v422 = correlator.Three(4, -2, -2);
		double v422Re = v422.Re()/Dn3;
		profile.fChc422->Fill(Ntrks, v422Re, Dn3);

		//..v5{psi32}
		TComplex v532 = correlator.Three(5, -3, -2);
		double v532Re = v532.Re()/Dn3;
		profile.fChc532->Fill(Ntrks, v532Re, Dn3);
	}

        // A-type
	if(NtrksAfterGap10M > 0 && NtrksAfterGap10P > 1 && Dn3Gap10A != 0)
	{

		TComplex v422Gap10A = correlator.ThreeGap10A(4, -2, -2);
		double v422Gap10ARe = v422Gap10A.Re()/Dn3Gap10A;
		profile.fChc422_Gap10A->Fill(Ntrks, v422Gap10ARe, Dn3Gap10A);

		TComplex v532Gap10A = correlator.ThreeGap10A(5, -3, -2);
		double v532Gap10ARe = v532Gap10A.Re()/Dn3Gap10A;
		profile.fChc532_Gap10A->Fill(Ntrks, v532Gap10ARe, Dn3Gap10A);
	}

        // B-type
	if(NtrksAfterGap10P > 0 && NtrksAfterGap10M > 1 && Dn3Gap10B != 0)
	{

		TComplex v422Gap10B = correlator.ThreeGap10B(4, -2, -2);
		double v422Gap10BRe = v422Gap10B.Re()/Dn3Gap10B;
		profile.fChc422_Gap10B->Fill(Ntrks, v422Gap10BRe, Dn3Gap10B);

		TComplex v532Gap10B = correlator.ThreeGap10B(5, -3, -2);
		double v532Gap10BRe = v532Gap10B.Re()/Dn3Gap10B;
		profile.fChc532_Gap10B->Fill(Ntrks, v532Gap10BRe, Dn3Gap10B);
	}

	//..calculate 4-particle correlations
	//................................
	double Dn4 = correlator.Four(0, 0, 0, 0).Re();
	double Dn4Gap10 = correlator.FourGap10(0, 0, 0, 0).Re();
	double Dn4_3subMMLR = correlator.Four_3SubMMLR(0, 0, 0, 0).Re();
	double Dn4_3subLLMR = correlator.Four_3SubLLMR(0, 0, 0, 0).Re();
	double Dn4_3subRRML = correlator.Four_3SubRRML(0, 0, 0, 0).Re();

	if(NtrksAfter > 3 && Dn4 != 0)
	{

		TComplex v24 = correlator.Four(2, 2, -2, -2);
		double v24Re = v24.Re()/Dn4;
		profile.fChcn4[0]->Fill(Ntrks, v24Re, Dn4);
		// fcn4Ntrks1bin[0][fBin]->Fill(NtrksAfter, v24Re, Dn4);

		TComplex v34 = correlator.Four(3, 3, -3, -3);
		double v34Re = v34.Re()/Dn4;
		profile.fChcn4[1]->Fill(Ntrks, v34Re, Dn4);
		// fcn4Ntrks1bin[1][fBin]->Fill(NtrksAfter, v34Re, Dn4);

		TComplex v44 = correlator.Four(4, 4, -4, -4);
		double v44Re = v44.Re()/Dn4;
		profile.fChcn4[2]->Fill(Ntrks, v44Re, Dn4);
		// fcn4Ntrks1bin[2][fBin]->Fill(NtrksAfter, v44Re, Dn4);

		//..SC(3,2,-3,-2)
		TComplex sc3232 = correlator.Four(3, 2, -3, -2);
		double sc3232Re = sc3232.Re()/Dn4;
		profile.fChsc3232->Fill(Ntrks, sc3232Re, Dn4);

		//..SC(4,2,-4,-2)
		TComplex sc4242 = correlator.Four(4, 2, -4, -2);
		double sc4242Re = sc4242.Re()/Dn4;
		profile.fChsc4242->Fill(Ntrks, sc4242Re, Dn4);

	}

	if(NtrksAfterGap10M > 1 && NtrksAfterGap10P > 1 && Dn4Gap10 !=0)
	{
		TComplex v24Gap10 = correlator.FourGap10(2, 2, -2, -2);
		double v24Gap10Re = v24Gap10.Re()/Dn4Gap10;
		profile.fChcn4_Gap10[0]->Fill(Ntrks, v24Gap10Re, Dn4Gap10);

		TComplex v34Gap10 = correlator.FourGap10(3, 3, -3, -3);
		double v34Gap10Re = v34Gap10.Re()/Dn4Gap10;
		profile.fChcn4_Gap10[1]->Fill(Ntrks, v34Gap10Re, Dn4Gap10);

		TComplex v44Gap10 = correlator.FourGap10(4, 4, -4, -4);
		double v44Gap10Re = v44Gap10.Re()/Dn4Gap10;
		profile.fChcn4_Gap10[2]->Fill(Ntrks, v44Gap10Re, Dn4Gap10);

		TComplex sc3232Gap10 = correlator.FourGap10(3, 2, -3, -2);
		double sc3232Gap10Re = sc3232Gap10.Re()/Dn4Gap10;
		profile.fChsc3232_Gap10->Fill(Ntrks, sc3232Gap10Re, Dn4Gap10);

		TComplex sc4242Gap10 = correlator.FourGap10(4, 2, -4, -2);
		double sc4242Gap10Re = sc4242Gap10.Re()/Dn4Gap10;
		profile.fChsc4242_Gap10->Fill(Ntrks, sc4242Gap10Re, Dn4Gap10);
	}

	//..3-subevent method
	if(NtrksAfter3subL > 0 && NtrksAfter3subR > 0 && NtrksAfter3subM > 1 && Dn4_3subMMLR != 0)
	{
		TComplex v24_3sub = correlator.Four_3SubMMLR(2, 2, -2, -2);
		double v24_3subRe = v24_3sub.Re()/Dn4_3subMMLR;
		profile.fChcn4_3subMMLR[0]->Fill(Ntrks, v24_3subRe, Dn4_3subMMLR);

		TComplex v34_3sub = correlator.Four_3SubMMLR(3, 3, -3, -3);
		double v34_3subRe = v34_3sub.Re()/Dn4_3subMMLR;
		profile.fChcn4_3subMMLR[1]->Fill(Ntrks, v34_3subRe, Dn4_3subMMLR);

		TComplex v44_3sub = correlator.Four_3SubMMLR(4, 4, -4, -4);
		double v44_3subRe = v44_3sub.Re()/Dn4_3subMMLR;
		profile.fChcn4_3subMMLR[2]->Fill(Ntrks, v44_3subRe, Dn4_3subMMLR);

		TComplex sc3232_3subA = correlator.Four_3SubMMLR(3, 2, -3, -2);
		double sc3232_3subARe = sc3232_3subA.Re()/Dn4_3subMMLR;
		profile.fChsc3232_3subMMLRA->Fill(Ntrks, sc3232_3subARe, Dn4_3subMMLR);

		TComplex sc3232_3subB = correlator.Four_3SubMMLR(3, 2, -2, -3);
		double sc3232_3subBRe = sc3232_3subB.Re()/Dn4_3subMMLR;
		profile.fChsc3232_3subMMLRB->Fill(Ntrks, sc3232_3subBRe, Dn4_3subMMLR);

		TComplex sc4242_3subA = correlator.Four_3SubMMLR(4, 2, -4, -2);
		double sc4242_3subARe = sc4242_3subA.Re()/Dn4_3subMMLR;
		profile.fChsc4242_3subMMLRA->Fill(Ntrks, sc4242_3subARe, Dn4_3subMMLR);

		TComplex sc4242_3subB = correlator.Four_3SubMMLR(4, 2, -2, -4);
		double sc4242_3subBRe = sc4242_3subB.Re()/Dn4_3subMMLR;
		profile.fChsc4242_3subMMLRB->Fill(Ntrks, sc4242_3subBRe, Dn4_3subMMLR);
	}

	//..3-subevent method
	if(NtrksAfter3subL > 1 && NtrksAfter3subR > 0 && NtrksAfter3subM > 0 && Dn4_3subLLMR != 0)
	{
		TComplex v24_3sub = correlator.Four_3SubLLMR(2, 2, -2, -2);
		double v24_3subRe = v24_3sub.Re()/Dn4_3subLLMR;
		profile.fChcn4_3subLLMR[0]->Fill(Ntrks, v24_3subRe, Dn4_3subLLMR);

		TComplex v34_3sub = correlator.Four_3SubLLMR(3, 3, -3, -3);
		double v34_3subRe = v34_3sub.Re()/Dn4_3subLLMR;
		profile.fChcn4_3subLLMR[1]->Fill(Ntrks, v34_3subRe, Dn4_3subLLMR);

		TComplex v44_3sub = correlator.Four_3SubLLMR(4, 4, -4, -4);
		double v44_3subRe = v44_3sub.Re()/Dn4_3subLLMR;
		profile.fChcn4_3subLLMR[2]->Fill(Ntrks, v44_3subRe, Dn4_3subLLMR);

		TComplex sc3232_3subA = correlator.Four_3SubLLMR(3, 2, -3, -2);
		double sc3232_3subARe = sc3232_3subA.Re()/Dn4_3subLLMR;
		profile.fChsc3232_3subLLMRA->Fill(Ntrks, sc3232_3subARe, Dn4_3subLLMR);

		TComplex sc3232_3subB = correlator.Four_3SubLLMR(3, 2, -2, -3);
		double sc3232_3subBRe = sc3232_3subB.Re()/Dn4_3subLLMR;
		profile.fChsc3232_3subLLMRB->Fill(Ntrks, sc3232_3subBRe, Dn4_3subLLMR);

		TComplex sc4242_3subA = correlator.Four_3SubLLMR(4, 2, -4, -2);
		double sc4242_3subARe = sc4242_3subA.Re()/Dn4_3subLLMR;
		profile.fChsc4242_3subLLMRA->Fill(Ntrks, sc4242_3subARe, Dn4_3subLLMR);

		TComplex sc4242_3subB = correlator.Four_3SubLLMR(4, 2, -2, -4);
		double sc4242_3subBRe = sc4242_3subB.Re()/Dn4_3subLLMR;
		profile.fChsc4242_3subLLMRB->Fill(Ntrks, sc4242_3subBRe, Dn4_3subLLMR);
	}

	//..3-subevent method
	if(NtrksAfter3subL > 0 && NtrksAfter3subR > 1 && NtrksAfter3subM > 0 && Dn4_3subRRML != 0)
	{
		TComplex v24_3sub = correlator.Four_3SubRRML(2, 2, -2, -2);
		double v24_3subRe = v24_3sub.Re()/Dn4_3subRRML;
		profile.fChcn4_3subRRML[0]->Fill(Ntrks, v24_3subRe, Dn4_3subRRML);

		TComplex v34_3sub = correlator.Four_3SubRRML(3, 3, -3, -3);
		double v34_3subRe = v34_3sub.Re()/Dn4_3subRRML;
		profile.fChcn4_3subRRML[1]->Fill(Ntrks, v34_3subRe, Dn4_3subRRML);

		TComplex v44_3sub = correlator.Four_3SubRRML(4, 4, -4, -4);
		double v44_3subRe = v44_3sub.Re()/Dn4_3subRRML;
		profile.fChcn4_3subRRML[2]->Fill(Ntrks, v44_3subRe, Dn4_3subRRML);

		TComplex sc3232_3subA = correlator.Four_3SubRRML(3, 2, -3, -2);
		double sc3232_3subARe = sc3232_3subA.Re()/Dn4_3subRRML;
		profile.fChsc3232_3subRRMLA->Fill(Ntrks, sc3232_3subARe, Dn4_3subRRML);

		TComplex sc3232_3subB = correlator.Four_3SubRRML(3, 2, -2, -3);
		double sc3232_3subBRe = sc3232_3subB.Re()/Dn4_3subRRML;
		profile.fChsc3232_3subRRMLB->Fill(Ntrks, sc3232_3subBRe, Dn4_3subRRML);

		TComplex sc4242_3subA = correlator.Four_3SubRRML(4, 2, -4, -2);
		double sc4242_3subARe = sc4242_3subA.Re()/Dn4_3subRRML;
		profile.fChsc4242_3subRRMLA->Fill(Ntrks, sc4242_3subARe, Dn4_3subRRML);

		TComplex sc4242_3subB = correlator.Four_3SubRRML(4, 2, -2, -4);
		double sc4242_3subBRe = sc4242_3subB.Re()/Dn4_3subRRML;
		profile.fChsc4242_3subRRMLB->Fill(Ntrks, sc4242_3subBRe, Dn4_3subRRML);
	}

}

const char* AliAnalysisTaskNonlinearFlow::ReturnPPperiod(const Int_t runNumber) const
{
  Bool_t isHM = kFALSE;
  if(fAliTrigger == AliVEvent::kHighMultV0) isHM = kTRUE;

  if(runNumber >= 252235 && runNumber <= 264347){ // LHC16
    if(!isHM && runNumber >= 252235 && runNumber <= 252375) return "LHC16de"; //d
    if(!isHM && runNumber >= 253437 && runNumber <= 253591) return "LHC16de"; //e
    if(runNumber >= 254128 && runNumber <= 254332) return "LHC16ghi"; //g
    if(runNumber >= 254604 && runNumber <= 255467) return "LHC16ghi"; //h
    if(runNumber >= 255539 && runNumber <= 255618) return "LHC16ghi"; //i
    if(runNumber >= 256219 && runNumber <= 256418) return "LHC16j";
    if(runNumber >= 256941 && runNumber <= 258537) return "LHC16k";
    if(runNumber >= 258962 && runNumber <= 259888) return "LHC16l";
    if(runNumber >= 262424 && runNumber <= 264035) return "LHC16o";
    if(runNumber >= 264076 && runNumber <= 264347) return "LHC16p";
  }

  if(runNumber >= 270581 && runNumber <= 282704){ // LHC17
    if(!isHM && runNumber >= 270581 && runNumber <= 270667) return "LHC17ce";
    if(runNumber >= 270822 && runNumber <= 270830){
      if(isHM) return "averaged";
      else return "LHC17ce";
    }
    if(runNumber >= 270854 && runNumber <= 270865){
      if(isHM) return "averaged";
      else return "LHC17f";
    }
    if(runNumber >= 271870 && runNumber <= 273103) return "LHC17h";
    if(runNumber >= 273591 && runNumber <= 274442) return "LHC17i";
    if(!isHM && runNumber >= 274593 && runNumber <= 274671) return "LHC17j";
    if(runNumber >= 274690 && runNumber <= 276508) return "LHC17k";
    if(runNumber >= 276551 && runNumber <= 278216) return "LHC17l";
    if(runNumber >= 278914 && runNumber <= 280140) return "LHC17m";
    if(runNumber >= 280282 && runNumber <= 281961) return "LHC17o";
    if(runNumber >= 282528 && runNumber <= 282704) return "LHC17r";
  }

  if(runNumber >= 285009 && runNumber <= 294925){ // LHC18
    if(runNumber >= 285009 && runNumber <= 285396){
      if(isHM) return "LHC18bd";
      else return "LHC18b";
    }
    if(runNumber >= 285978 && runNumber <= 286350){
      if(isHM) return "LHC18bd";
      else return "LHC18d";
    }
    if(runNumber >= 286380 && runNumber <= 286937) return "LHC18e";
    if(runNumber >= 287000 && runNumber <= 287658) return "LHC18f";
    if(runNumber >= 288804 && runNumber <= 288806){
      if(isHM) return "LHC18hjk";
      else return "LHC18ghijk";
    }
    if(runNumber == 288943){
      if(isHM) return "LHC18hjk";
      else return "LHC18ghijk";
    }
    if(runNumber >= 289165 && runNumber <= 289201){
      if(isHM) return "LHC18hjk";
      else return "LHC18ghijk";
    }
    if(!isHM && runNumber >= 288619 && runNumber <= 288750) return "LHC18ghijk"; //g, no HM event, only MB
    if(!isHM && runNumber >= 288861 && runNumber <= 288909) return "LHC18ghijk"; //i, no HM event, only MB
    if(runNumber >= 289240 && runNumber <= 289971) return "LHC18l";
    if(runNumber >= 290323 && runNumber <= 292839){
      if(isHM) return "LHC18m";
      else return "LHC18mn";
    }
    if(!isHM && runNumber >= 293357 && runNumber <= 293359) return "LHC18mn"; //n, no HM event, only MB
    if(runNumber >= 293475 && runNumber <= 293898) return "LHC18o";
    if(runNumber >= 294009 && runNumber <= 294925) return "LHC18p";
  }

  AliWarning("Unknown period! Returning averaged weights");
  return "averaged";
}


//_____________________________________________________________________________
void AliAnalysisTaskNonlinearFlow::Terminate(Option_t *)
{
	// Terminate loop
	Printf("Terminate()");
}

ClassImp(PhysicsProfile);
PhysicsProfile::PhysicsProfile() :
		fChsc4242(nullptr),
		fChsc4242_Gap0(nullptr),
		fChsc4242_Gap2(nullptr),
		fChsc4242_Gap4(nullptr),
		fChsc4242_Gap6(nullptr),
		fChsc4242_Gap8(nullptr),      
		fChsc4242_Gap10(nullptr),    
		fChsc4242_3sub(nullptr),	
		fChsc4242_3subMMLRA(nullptr),
		fChsc4242_3subMMLRB(nullptr),							
		fChsc4242_3subLLMRA(nullptr),							
		fChsc4242_3subLLMRB(nullptr),							
		fChsc4242_3subRRMLA(nullptr),							
		fChsc4242_3subRRMLB(nullptr),							
		fChsc4224_3sub(nullptr),							
		fChsc4242_3subGap2(nullptr),					
		fChsc4224_3subGap2(nullptr),					
		fChsc3232(nullptr),									
		fChsc3232_Gap0(nullptr),							
		fChsc3232_Gap2(nullptr),							
		fChsc3232_Gap4(nullptr),                            
		fChsc3232_Gap6(nullptr),                            
		fChsc3232_Gap8(nullptr),                            
		fChsc3232_Gap10(nullptr),                            
		fChsc3232_3sub(nullptr),							
		fChsc3232_3subMMLRA(nullptr),							
		fChsc3232_3subMMLRB(nullptr),							
		fChsc3232_3subLLMRA(nullptr),							
		fChsc3232_3subLLMRB(nullptr),							
		fChsc3232_3subRRMLA(nullptr),							
		fChsc3232_3subRRMLB(nullptr),							
		fChsc3223_3sub(nullptr),							
		fChsc3232_3subGap2(nullptr),					
		fChsc3223_3subGap2(nullptr),					
		fChc422(nullptr), 
		fChc532(nullptr),
		fChc422_Gap0A(nullptr),   
		fChc422_Gap0B(nullptr),   
		fChc532_Gap0A(nullptr),   
		fChc532_Gap0B(nullptr),   
		fChc422_Gap2A(nullptr),   
		fChc422_Gap2B(nullptr),   
		fChc532_Gap2A(nullptr),   
		fChc532_Gap2B(nullptr),   
		fChc422_Gap4A(nullptr),   
		fChc422_Gap4B(nullptr),   
		fChc532_Gap4A(nullptr),   
		fChc532_Gap4B(nullptr),   
		fChc422_Gap6A(nullptr),   
		fChc422_Gap6B(nullptr),   
		fChc532_Gap6A(nullptr),   
		fChc532_Gap6B(nullptr),   
		fChc422_Gap8A(nullptr),   
		fChc422_Gap8B(nullptr),   
		fChc532_Gap8A(nullptr),   
		fChc532_Gap8B(nullptr),   
		fChc422_Gap10A(nullptr), 
		fChc422_Gap10B(nullptr), 
		fChc532_Gap10A(nullptr), 
		fChc532_Gap10B(nullptr)
{
		memset(fChcn2, 0, sizeof(fChcn2));
		memset(fChcn2_Gap0, 0, sizeof(fChcn2_Gap0));
		memset(fChcn2_Gap2, 0, sizeof(fChcn2_Gap2));
		memset(fChcn2_Gap4, 0, sizeof(fChcn2_Gap4));
		memset(fChcn2_Gap6, 0, sizeof(fChcn2_Gap6));
		memset(fChcn2_Gap8, 0, sizeof(fChcn2_Gap8));
		memset(fChcn2_Gap10, 0, sizeof(fChcn2_Gap10));
		memset(fChcn2_Gap14, 0, sizeof(fChcn2_Gap14));
		memset(fChcn2_Gap16, 0, sizeof(fChcn2_Gap16));
		memset(fChcn2_Gap18, 0, sizeof(fChcn2_Gap18));

		memset(fChcn2_3subLM, 0, sizeof(fChcn2_3subLM));
		memset(fChcn2_3subRM, 0, sizeof(fChcn2_3subRM));
		memset(fChcn2_3subLR, 0, sizeof(fChcn2_3subLR));
		memset(fChcn2_3subGap2LM, 0, sizeof(fChcn2_3subGap2LM));
		memset(fChcn2_3subGap2RM, 0, sizeof(fChcn2_3subGap2RM));

		memset(fChcn4, 0, sizeof(fChcn4));
		memset(fChcn4_Gap0, 0, sizeof(fChcn4_Gap0));
		memset(fChcn4_Gap2, 0, sizeof(fChcn4_Gap2));
		memset(fChcn4_Gap4, 0, sizeof(fChcn4_Gap4));
		memset(fChcn4_Gap6, 0, sizeof(fChcn4_Gap6));
		memset(fChcn4_Gap8, 0, sizeof(fChcn4_Gap8));
		memset(fChcn4_Gap10, 0, sizeof(fChcn4_Gap10));
		memset(fChcn4_3subMMLR, 0, sizeof(fChcn4_3subMMLR));
		memset(fChcn4_3subLLMR, 0, sizeof(fChcn4_3subLLMR));
		memset(fChcn4_3subRRML, 0, sizeof(fChcn4_3subRRML));
		memset(fChcn4_3subGap2, 0, sizeof(fChcn4_3subGap2));
}
