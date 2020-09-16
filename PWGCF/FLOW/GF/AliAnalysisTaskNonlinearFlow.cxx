#include "AliAnalysisTaskNonlinearFlow.h"

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
		fFilterbit(768),
		fEtaCut(0.8),
		fVtxCut(10.0),
		fMinPt(0.2),
		fMaxPt(3.0),
		fTPCclusters(70),
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
		//....
		fPeriod("LHC15o"),

		fListOfObjects(0),

		fMultTOFLowCut(0),
		fMultTOFHighCut(0),
		fMultCentLowCut(0),

		fTrackEfficiency(0),
		hTrackEfficiency(0),
		hTrackEfficiencyRun(0),

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

		hNtrksPt0530Pt0230(0),
		hNtrksPt0730Pt0230(0),
		hNtrksEta09Eta10(0),
		hNtrksEta08Eta10(0),
		hNtrksAllNtrksLS(0),
		hNtrksNoGapGap0(0),
		hNtrksNoGapGap2(0),
		hNtrksNoGapGap4(0),
		hNtrksNoGapGap6(0),
		hNtrksNoGapGap8(0),
		hNtrksNoGapGap(0),
		hNtrksNoGapGap14(0),
		hNtrksNoGapGap16(0),
		hNtrksNoGapGap18(0),
		hNtrksNoGap3sub(0),
		hNtrksNoGap3subGap(0),

		hReco(0),
		hRecoPion(0),
		hRecoKaon(0),
		hRecoProton(0),
		hRecoElectron(0),
		hRecoMuon(0),
		hRecoLSplus(0),
		hRecoLSminus(0),
		hPtRecoNtrks(0),
		hEtaRecoNtrks(0),
		hVzRecoNtrks(0),
		hPtRecoNtrksReco(0),
		hEtaRecoNtrksReco(0),
		hVzRecoNtrksReco(0),
		hTruth(0),
		hTruthPion(0),
		hTruthKaon(0),
		hTruthProton(0),
		hTruthElectron(0),
		hTruthMuon(0),
		hTruthLSplus(0),
		hTruthLSminus(0),
		hPtTruthNtrks(0),
		hEtaTruthNtrks(0),
		hVzTruthNtrks(0),
		hPtTruthNtrksReco(0),
		hEtaTruthNtrksReco(0),
		hVzTruthNtrksReco(0),
		hNtrksRecoNtrksTruth(0),
		hNtrksRecoCorrNtrksTruth(0),

		hPrimary(0),
		hPions(0),
		hDCAptMC(0),
		hDCAptMC_material(0),
		hDCAptMC_weak(0) {

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
	fUseDCAxyCut(0),
	fDCAxy(0),
	fSample(1),
	fCentFlag(0),
	fTrigger(0),
	fLS(false),
	fNUE(0),
	fNUA(0),
	//....
	fPeriod("LHC15o"),

	fListOfObjects(0),

	fMultTOFLowCut(0),
	fMultTOFHighCut(0),
	fMultCentLowCut(0),

	fTrackEfficiency(0),
	hTrackEfficiency(0),
	hTrackEfficiencyRun(0),

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

	hNtrksPt0530Pt0230(0),
	hNtrksPt0730Pt0230(0),
	hNtrksEta09Eta10(0),
	hNtrksEta08Eta10(0),
	hNtrksAllNtrksLS(0),
	hNtrksNoGapGap0(0),
	hNtrksNoGapGap2(0),
	hNtrksNoGapGap4(0),
	hNtrksNoGapGap6(0),
	hNtrksNoGapGap8(0),
	hNtrksNoGapGap(0),
	hNtrksNoGapGap14(0),
	hNtrksNoGapGap16(0),
	hNtrksNoGapGap18(0),
	hNtrksNoGap3sub(0),
	hNtrksNoGap3subGap(0),

	//..MC
	hMultMC(0),
	fPhiDisTruth(0),
	fEtaDisTruth(0),
	fPtDisTruth(0),

	hReco(0),
	hRecoPion(0),
	hRecoKaon(0),
	hRecoProton(0),
	hRecoElectron(0),
	hRecoMuon(0),
	hRecoLSplus(0),
	hRecoLSminus(0),
	hPtRecoNtrks(0),
	hEtaRecoNtrks(0),
	hVzRecoNtrks(0),
	hPtRecoNtrksReco(0),
	hEtaRecoNtrksReco(0),
	hVzRecoNtrksReco(0),
	hTruth(0),
	hTruthPion(0),
	hTruthKaon(0),
	hTruthProton(0),
	hTruthElectron(0),
	hTruthMuon(0),
	hTruthLSplus(0),
	hTruthLSminus(0),
	hPtTruthNtrks(0),
	hEtaTruthNtrks(0),
	hVzTruthNtrks(0),
	hPtTruthNtrksReco(0),
	hEtaTruthNtrksReco(0),
	hVzTruthNtrksReco(0),
	hNtrksRecoNtrksTruth(0),
	hNtrksRecoCorrNtrksTruth(0),

	hPrimary(0),
	hPions(0),
	hDCAptMC(0),
	hDCAptMC_material(0),
	hDCAptMC_weak(0), rand(32213) {

		// Output slot #1 writes into a TList
		DefineOutput(1, TList::Class());
		DefineOutput(2, TList::Class());

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

	nn = 32;
	double xbins_tmp[] = {50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300,
		1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000 };
        for (int i = 0; i <= nn; i++) {
                xbins[i] = xbins_tmp[i];
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

	//..phi weight: it is done run-by-run, the histograms are obtained in GetWeight() function
	//    default
	fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/PhiWeight_LHC15o_HIR.root");
	//    for systematics
	if(fFilterbit == 768) fPhiWeight = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/PhiWeight_LHC15o_HIR_FB768.root");

	//    default
	fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_LHC15o_HIR.root");
	// for systematics
	if(fFilterbit == 768) fTrackEfficiency = TFile::Open("alien:///alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2015/TrackingEfficiency_PbPb5TeV_LHC15o_HIR_FB768.root");

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
	isTrigselected = fSelectMask&AliVEvent::kINT7;
	if(isTrigselected == false) return;

	//..check if I have AOD
	fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
	if(!fAOD){
		Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
		return;
	}
	if(!fEventCuts.AcceptEvent(fAOD)) { // automatic event selection for Run2
		PostData(1,fListOfObjects);
		return;
	}
	hEventCount->Fill("after fEventCuts", 1.);

	//..filling Vz distribution
	AliVVertex *vtx = fAOD->GetPrimaryVertex();
	float fVtxZ = vtx->GetZ();
	if(TMath::Abs(fVtxZ) > fVtxCut) return;
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

	//..all charged particles
	if(fLS == true){
		AnalyzeAOD(fInputEvent, centrV0, cent, centSPD, fVtxZ, false);
		AnalyzeAOD(fInputEvent, centrV0, cent, centSPD, fVtxZ, true);
	} else AnalyzeAOD(fInputEvent, centrV0, cent, centSPD, fVtxZ, false);

	// Post output data.
	PostData(1, fListOfObjects);
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

		int nClustersITS = 0;
		nClustersITS = aodTrk->GetITSNcls();
		float chi2PerClusterITS = -1;
		if(nClustersITS != 0) chi2PerClusterITS = aodTrk->GetITSchi2()/float(nClustersITS);

		if (!(aodTrk->TestFilterBit(fFilterbit))) { continue; }
		if (fFilterbit == 96) {
			if (TMath::Abs(dcaZ) > fDCAz) continue;
		}

		if(aodTrk->Pt() < fMinPt) continue;
		if(aodTrk->Pt() > fMaxPt) continue;

		if(TMath::Abs(aodTrk->Eta()) > fEtaCut) continue;

		NtrksAfter += 1;

		//..get phi-weight for NUA correction
		double weight = 1;
		if(fNUA == 1) {
			weight = GetWeight(aodTrk->Phi(), aodTrk->Eta(), aodTrk->Pt(), run, fPlus, fVtxZ, runNumber);
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
	for(int iharm=0; iharm<20; iharm++)
	{
		for(int ipow=0; ipow<20; ipow++)
		{
			Qvector[iharm][ipow] = TComplex(Qcos[iharm][ipow], Qsin[iharm][ipow]);
			Qvector10M[iharm][ipow] = TComplex(QcosGap10M[iharm][ipow], QsinGap10M[iharm][ipow]);
			Qvector10P[iharm][ipow] = TComplex(QcosGap10P[iharm][ipow], QsinGap10P[iharm][ipow]);
			Qvector14M[iharm][ipow] = TComplex(QcosGap14M[iharm][ipow], QsinGap14M[iharm][ipow]);
			Qvector14P[iharm][ipow] = TComplex(QcosGap14P[iharm][ipow], QsinGap14P[iharm][ipow]);
			QvectorSubLeft[iharm][ipow] = TComplex(QcosSubLeft[iharm][ipow], QsinSubLeft[iharm][ipow]);
			QvectorSubRight[iharm][ipow] = TComplex(QcosSubRight[iharm][ipow], QsinSubRight[iharm][ipow]);
			QvectorSubMiddle[iharm][ipow] = TComplex(QcosSubMiddle[iharm][ipow], QsinSubMiddle[iharm][ipow]);
		}
	}

	hMult->Fill(NtrksAfter);

	// CalculateProfile(centProfile, cent);
	CalculateProfile(multProfile, NtrksAfter);
	// CalculateProfile(centProfile_bin[bootstrap_value], cent);
	CalculateProfile(multProfile_bin[bootstrap_value], NtrksAfter);

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

	double weight = 1;
	hTrackEfficiencyRun = (TH3F*)fTrackEfficiency->Get(Form("eff_LHC15o_HIJING_%.0lf", runNumber));
	double binPt = hTrackEfficiencyRun->GetXaxis()->FindBin(pt);
	double binEta = hTrackEfficiencyRun->GetYaxis()->FindBin(eta);
	double binVz = hTrackEfficiencyRun->GetZaxis()->FindBin(vz);
	//..take into account error on efficiency: randomly get number from gaussian distribution of eff. where width = error
	double eff = hTrackEfficiencyRun->GetBinContent(binPt, binEta, binVz);
	double error = hTrackEfficiencyRun->GetBinError(binPt, binEta, binVz);

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
	double weight = 1;
	hPhiWeightRun = (TH3F*)fPhiWeight->Get(Form("fPhiWeight_%0.lf", runNumber));
	weight = hPhiWeightRun->GetBinContent(hPhiWeightRun->GetXaxis()->FindBin(phi),
			hPhiWeightRun->GetYaxis()->FindBin(eta),
			hPhiWeightRun->GetZaxis()->FindBin(vz));
	return weight;
}
//_____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Q(int n, int p)
{

	if(n>=0) return Qvector[n][p];
	else return TComplex::Conjugate(Qvector[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap10M(int n, int p)
{

	if(n>=0) return Qvector10M[n][p];
	else return TComplex::Conjugate(Qvector10M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap10P(int n, int p)
{

	if(n>=0) return Qvector10P[n][p];
	else return TComplex::Conjugate(Qvector10P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap0M(int n, int p)
{

	if(n>=0) return Qvector0M[n][p];
	else return TComplex::Conjugate(Qvector0M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap0P(int n, int p)
{

	if(n>=0) return Qvector0P[n][p];
	else return TComplex::Conjugate(Qvector0P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap2M(int n, int p)
{

	if(n>=0) return Qvector2M[n][p];
	else return TComplex::Conjugate(Qvector2M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap2P(int n, int p)
{

	if(n>=0) return Qvector2P[n][p];
	else return TComplex::Conjugate(Qvector2P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap4M(int n, int p)
{

	if(n>=0) return Qvector4M[n][p];
	else return TComplex::Conjugate(Qvector4M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap4P(int n, int p)
{

	if(n>=0) return Qvector4P[n][p];
	else return TComplex::Conjugate(Qvector4P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap6M(int n, int p)
{

	if(n>=0) return Qvector6M[n][p];
	else return TComplex::Conjugate(Qvector6M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap6P(int n, int p)
{

	if(n>=0) return Qvector6P[n][p];
	else return TComplex::Conjugate(Qvector6P[-n][p]);

}

//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap8M(int n, int p)
{

	if(n>=0) return Qvector8M[n][p];
	else return TComplex::Conjugate(Qvector8M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap8P(int n, int p)
{

	if(n>=0) return Qvector8P[n][p];
	else return TComplex::Conjugate(Qvector8P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::pGap10M(int n, int p)
{

	if(n>=0) return pvectorM[n][p];
	else return TComplex::Conjugate(pvectorM[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::pGap10P(int n, int p)
{

	if(n>=0) return pvectorP[n][p];
	else return TComplex::Conjugate(pvectorP[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap14M(int n, int p)
{

	if(n>=0) return Qvector14M[n][p];
	else return TComplex::Conjugate(Qvector14M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap14P(int n, int p)
{

	if(n>=0) return Qvector14P[n][p];
	else return TComplex::Conjugate(Qvector14P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap16M(int n, int p)
{

	if(n>=0) return Qvector16M[n][p];
	else return TComplex::Conjugate(Qvector16M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap16P(int n, int p)
{

	if(n>=0) return Qvector16P[n][p];
	else return TComplex::Conjugate(Qvector16P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap18M(int n, int p)
{

	if(n>=0) return Qvector18M[n][p];
	else return TComplex::Conjugate(Qvector18M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QGap18P(int n, int p)
{

	if(n>=0) return Qvector18P[n][p];
	else return TComplex::Conjugate(Qvector18P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QsubLeft(int n, int p)
{

	if(n>=0) return QvectorSubLeft[n][p];
	else return TComplex::Conjugate(QvectorSubLeft[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QsubRight(int n, int p)
{

	if(n>=0) return QvectorSubRight[n][p];
	else return TComplex::Conjugate(QvectorSubRight[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QsubMiddle(int n, int p)
{

	if(n>=0) return QvectorSubMiddle[n][p];
	else return TComplex::Conjugate(QvectorSubMiddle[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QsubGap2Left(int n, int p)
{

	if(n>=0) return QvectorSubGap2Left[n][p];
	else return TComplex::Conjugate(QvectorSubGap2Left[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QsubGap2Right(int n, int p)
{

	if(n>=0) return QvectorSubGap2Right[n][p];
	else return TComplex::Conjugate(QvectorSubGap2Right[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::QsubGap2Middle(int n, int p)
{

	if(n>=0) return QvectorSubGap2Middle[n][p];
	else return TComplex::Conjugate(QvectorSubGap2Middle[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::p(int n, int p)
{

	if(n>=0) return pvector[n][p];
	else return TComplex::Conjugate(pvector[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::q(int n, int p)
{

	if(n>=0) return qvector[n][p];
	else return TComplex::Conjugate(qvector[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::pGap0M(int n, int p)
{

	if(n>=0) return pvector0M[n][p];
	else return TComplex::Conjugate(pvector0M[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::pGap0P(int n, int p)
{

	if(n>=0) return pvector0P[n][p];
	else return TComplex::Conjugate(pvector0P[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::pGap4M(int n, int p)
{

	if(n>=0) return pvector4M[n][p];
	else return TComplex::Conjugate(pvector4M[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::pGap4P(int n, int p)
{

	if(n>=0) return pvector4P[n][p];
	else return TComplex::Conjugate(pvector4P[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::pGap8M(int n, int p)
{

	if(n>=0) return pvector8M[n][p];
	else return TComplex::Conjugate(pvector8M[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::pGap8P(int n, int p)
{

	if(n>=0) return pvector8P[n][p];
	else return TComplex::Conjugate(pvector8P[n][p]);

}
//____________________________________________________________________
void AliAnalysisTaskNonlinearFlow::ResetQ(const int nMaxHarm, const int nMaxPow)
{

	for(int i=0; i<nMaxHarm; i++)
	{
		for(int j=0; j<nMaxPow; j++)
		{
			Qvector[i][j] = TComplex(0.,0.);
		}
	}

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Two(int n1, int n2)
{

	TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap0(int n1, int n2)
{

	TComplex formula = QGap0M(n1,1)*QGap0P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap2(int n1, int n2)
{

	TComplex formula = QGap2M(n1,1)*QGap2P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap4(int n1, int n2)
{

	TComplex formula = QGap4M(n1,1)*QGap4P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap6(int n1, int n2)
{

	TComplex formula = QGap6M(n1,1)*QGap6P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap8(int n1, int n2)
{

	TComplex formula = QGap8M(n1,1)*QGap8P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap10(int n1, int n2)
{

	TComplex formula = QGap10M(n1,1)*QGap10P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap14(int n1, int n2)
{

	TComplex formula = QGap14M(n1,1)*QGap14P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap16(int n1, int n2)
{

	TComplex formula = QGap16M(n1,1)*QGap16P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap18(int n1, int n2)
{

	TComplex formula = QGap18M(n1,1)*QGap18P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Two_3SubLM(int n1, int n2)
{

	TComplex formula = QsubLeft(n1,1)*QsubMiddle(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Two_3SubRM(int n1, int n2)
{

	TComplex formula = QsubMiddle(n1,1)*QsubRight(n2,1);
	return formula;

}
//
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Two_3SubLR(int n1, int n2)
{

	TComplex formula = QsubLeft(n1,1)*QsubRight(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Two_3SubGap2LM(int n1, int n2)
{

	TComplex formula = QsubGap2Left(n1,1)*QsubGap2Middle(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Two_3SubGap2RM(int n1, int n2)
{

	TComplex formula = QsubGap2Middle(n1,1)*QsubGap2Right(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap0M(int n1, int n2)
{

	TComplex formula = QGap0M(n1,1)*QGap0M(n2,1) - QGap0M(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap2M(int n1, int n2)
{

	TComplex formula = QGap2M(n1,1)*QGap2M(n2,1) - QGap2M(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap4M(int n1, int n2)
{

	TComplex formula = QGap4M(n1,1)*QGap4M(n2,1) - QGap4M(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap6M(int n1, int n2)
{

	TComplex formula = QGap6M(n1,1)*QGap6M(n2,1) - QGap6M(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap8M(int n1, int n2)
{

	TComplex formula = QGap8M(n1,1)*QGap8M(n2,1) - QGap8M(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap10M(int n1, int n2)
{

	TComplex formula = QGap10M(n1,1)*QGap10M(n2,1) - QGap10M(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap0P(int n1, int n2)
{

	TComplex formula = QGap0P(n1,1)*QGap0P(n2,1) - QGap0P(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap2P(int n1, int n2)
{

	TComplex formula = QGap2P(n1,1)*QGap2P(n2,1) - QGap2P(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap4P(int n1, int n2)
{

	TComplex formula = QGap4P(n1,1)*QGap4P(n2,1) - QGap4P(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap6P(int n1, int n2)
{

	TComplex formula = QGap6P(n1,1)*QGap6P(n2,1) - QGap6P(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap8P(int n1, int n2)
{

	TComplex formula = QGap8P(n1,1)*QGap8P(n2,1) - QGap8P(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::TwoGap10P(int n1, int n2)
{

	TComplex formula = QGap10P(n1,1)*QGap10P(n2,1) - QGap10P(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Three(int n1, int n2, int n3)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
		- Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap0A(int n1, int n2, int n3)
{

	TComplex formula = QGap0M(n1,1)*QGap0P(n2,1)*QGap0P(n3,1)-QGap0M(n1,1)*QGap0P(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap0B(int n1, int n2, int n3)
{

	TComplex formula = QGap0P(n1,1)*QGap0M(n2,1)*QGap0M(n3,1)-QGap0P(n1,1)*QGap0M(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap2A(int n1, int n2, int n3)
{

	TComplex formula = QGap2M(n1,1)*QGap2P(n2,1)*QGap2P(n3,1)-QGap2M(n1,1)*QGap2P(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap2B(int n1, int n2, int n3)
{

	TComplex formula = QGap2P(n1,1)*QGap2M(n2,1)*QGap2M(n3,1)-QGap2P(n1,1)*QGap2M(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap4A(int n1, int n2, int n3)
{

	TComplex formula = QGap4M(n1,1)*QGap4P(n2,1)*QGap4P(n3,1)-QGap4M(n1,1)*QGap4P(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap4B(int n1, int n2, int n3)
{

	TComplex formula = QGap4P(n1,1)*QGap4M(n2,1)*QGap4M(n3,1)-QGap4P(n1,1)*QGap4M(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap6A(int n1, int n2, int n3)
{

	TComplex formula = QGap6M(n1,1)*QGap6P(n2,1)*QGap6P(n3,1)-QGap6M(n1,1)*QGap6P(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap6B(int n1, int n2, int n3)
{

	TComplex formula = QGap6P(n1,1)*QGap6M(n2,1)*QGap6M(n3,1)-QGap6P(n1,1)*QGap6M(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap8A(int n1, int n2, int n3)
{

	TComplex formula = QGap8M(n1,1)*QGap8P(n2,1)*QGap8P(n3,1)-QGap8M(n1,1)*QGap8P(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap8B(int n1, int n2, int n3)
{

	TComplex formula = QGap8P(n1,1)*QGap8M(n2,1)*QGap8M(n3,1)-QGap8P(n1,1)*QGap8M(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap10A(int n1, int n2, int n3)
{

	TComplex formula = QGap10M(n1,1)*QGap10P(n2,1)*QGap10P(n3,1)-QGap10M(n1,1)*QGap10P(n2+n3,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap10B(int n1, int n2, int n3)
{

	TComplex formula = QGap10P(n1,1)*QGap10M(n2,1)*QGap10M(n3,1)-QGap10P(n1,1)*QGap10M(n2+n3,2);
	return formula;

}


//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap0_subM(int n1, int n2, int n3)
{

	TComplex formula = QGap0M(n1,1)*QGap0M(n2,1)*QGap0M(n3,1)-QGap0M(n1+n2,2)*QGap0M(n3,1)-QGap0M(n2,1)*QGap0M(n1+n3,2)
		- QGap0M(n1,1)*QGap0M(n2+n3,2)+2.*QGap0M(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap0_subP(int n1, int n2, int n3)
{

	TComplex formula = QGap0P(n1,1)*QGap0P(n2,1)*QGap0P(n3,1)-QGap0P(n1+n2,2)*QGap0P(n3,1)-QGap0P(n2,1)*QGap0P(n1+n3,2)
		- QGap0P(n1,1)*QGap0P(n2+n3,2)+2.*QGap0P(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap2_subM(int n1, int n2, int n3)
{

	TComplex formula = QGap2M(n1,1)*QGap2M(n2,1)*QGap2M(n3,1)-QGap2M(n1+n2,2)*QGap2M(n3,1)-QGap2M(n2,1)*QGap2M(n1+n3,2)
		- QGap2M(n1,1)*QGap2M(n2+n3,2)+2.*QGap2M(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap2_subP(int n1, int n2, int n3)
{

	TComplex formula = QGap2P(n1,1)*QGap2P(n2,1)*QGap2P(n3,1)-QGap2P(n1+n2,2)*QGap2P(n3,1)-QGap2P(n2,1)*QGap2P(n1+n3,2)
		- QGap2P(n1,1)*QGap2P(n2+n3,2)+2.*QGap2P(n1+n2+n3,3);
	return formula;

}

//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap4_subM(int n1, int n2, int n3)
{

	TComplex formula = QGap4M(n1,1)*QGap4M(n2,1)*QGap4M(n3,1)-QGap4M(n1+n2,2)*QGap4M(n3,1)-QGap4M(n2,1)*QGap4M(n1+n3,2)
		- QGap4M(n1,1)*QGap4M(n2+n3,2)+2.*QGap4M(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap4_subP(int n1, int n2, int n3)
{

	TComplex formula = QGap4P(n1,1)*QGap4P(n2,1)*QGap4P(n3,1)-QGap4P(n1+n2,2)*QGap4P(n3,1)-QGap4P(n2,1)*QGap4P(n1+n3,2)
		- QGap4P(n1,1)*QGap4P(n2+n3,2)+2.*QGap4P(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap6_subM(int n1, int n2, int n3)
{

	TComplex formula = QGap6M(n1,1)*QGap6M(n2,1)*QGap6M(n3,1)-QGap6M(n1+n2,2)*QGap6M(n3,1)-QGap6M(n2,1)*QGap6M(n1+n3,2)
		- QGap6M(n1,1)*QGap6M(n2+n3,2)+2.*QGap6M(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap6_subP(int n1, int n2, int n3)
{

	TComplex formula = QGap6P(n1,1)*QGap6P(n2,1)*QGap6P(n3,1)-QGap6P(n1+n2,2)*QGap6P(n3,1)-QGap6P(n2,1)*QGap6P(n1+n3,2)
		- QGap6P(n1,1)*QGap6P(n2+n3,2)+2.*QGap6P(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap8_subM(int n1, int n2, int n3)
{

	TComplex formula = QGap8M(n1,1)*QGap8M(n2,1)*QGap8M(n3,1)-QGap8M(n1+n2,2)*QGap8M(n3,1)-QGap8M(n2,1)*QGap8M(n1+n3,2)
		- QGap8M(n1,1)*QGap8M(n2+n3,2)+2.*QGap8M(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap8_subP(int n1, int n2, int n3)
{

	TComplex formula = QGap8P(n1,1)*QGap8P(n2,1)*QGap8P(n3,1)-QGap8P(n1+n2,2)*QGap8P(n3,1)-QGap8P(n2,1)*QGap8P(n1+n3,2)
		- QGap8P(n1,1)*QGap8P(n2+n3,2)+2.*QGap8P(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap10_subM(int n1, int n2, int n3)
{

	TComplex formula = QGap10M(n1,1)*QGap10M(n2,1)*QGap10M(n3,1)-QGap10M(n1+n2,2)*QGap10M(n3,1)-QGap10M(n2,1)*QGap10M(n1+n3,2)
		- QGap10M(n1,1)*QGap10M(n2+n3,2)+2.*QGap10M(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::ThreeGap10_subP(int n1, int n2, int n3)
{

	TComplex formula = QGap10P(n1,1)*QGap10P(n2,1)*QGap10P(n3,1)-QGap10P(n1+n2,2)*QGap10P(n3,1)-QGap10P(n2,1)*QGap10P(n1+n3,2)
		- QGap10P(n1,1)*QGap10P(n2+n3,2)+2.*QGap10P(n1+n2+n3,3);
	return formula;

}

//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Four(int n1, int n2, int n3, int n4)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
		- Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
		+ Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
		+ 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FourGap0(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap0P(n1,1)*QGap0P(n2,1)*QGap0M(n3,1)*QGap0M(n4,1)-QGap0P(n1+n2,2)*QGap0M(n3,1)*QGap0M(n4,1)
		-QGap0P(n1,1)*QGap0P(n2,1)*QGap0M(n3+n4,2)+QGap0P(n1+n2,2)*QGap0M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FourGap0M(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap0M(n1,1)*QGap0M(n2,1)*QGap0M(n3,1)*QGap0M(n4,1)-QGap0M(n1+n2,2)*QGap0M(n3,1)*QGap0M(n4,1)-QGap0M(n2,1)*QGap0M(n1+n3,2)*QGap0M(n4,1)
		- QGap0M(n1,1)*QGap0M(n2+n3,2)*QGap0M(n4,1)+2.*QGap0M(n1+n2+n3,3)*QGap0M(n4,1)-QGap0M(n2,1)*QGap0M(n3,1)*QGap0M(n1+n4,2)
		+ QGap0M(n2+n3,2)*QGap0M(n1+n4,2)-QGap0M(n1,1)*QGap0M(n3,1)*QGap0M(n2+n4,2)+QGap0M(n1+n3,2)*QGap0M(n2+n4,2)
		+ 2.*QGap0M(n3,1)*QGap0M(n1+n2+n4,3)-QGap0M(n1,1)*QGap0M(n2,1)*QGap0M(n3+n4,2)+QGap0M(n1+n2,2)*QGap0M(n3+n4,2)
		+ 2.*QGap0M(n2,1)*QGap0M(n1+n3+n4,3)+2.*QGap0M(n1,1)*QGap0M(n2+n3+n4,3)-6.*QGap0M(n1+n2+n3+n4,4);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FourGap0P(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap0P(n1,1)*QGap0P(n2,1)*QGap0P(n3,1)*QGap0P(n4,1)-QGap0P(n1+n2,2)*QGap0P(n3,1)*QGap0P(n4,1)-QGap0P(n2,1)*QGap0P(n1+n3,2)*QGap0P(n4,1)
		- QGap0P(n1,1)*QGap0P(n2+n3,2)*QGap0P(n4,1)+2.*QGap0P(n1+n2+n3,3)*QGap0P(n4,1)-QGap0P(n2,1)*QGap0P(n3,1)*QGap0P(n1+n4,2)
		+ QGap0P(n2+n3,2)*QGap0P(n1+n4,2)-QGap0P(n1,1)*QGap0P(n3,1)*QGap0P(n2+n4,2)+QGap0P(n1+n3,2)*QGap0P(n2+n4,2)
		+ 2.*QGap0P(n3,1)*QGap0P(n1+n2+n4,3)-QGap0P(n1,1)*QGap0P(n2,1)*QGap0P(n3+n4,2)+QGap0P(n1+n2,2)*QGap0P(n3+n4,2)
		+ 2.*QGap0P(n2,1)*QGap0P(n1+n3+n4,3)+2.*QGap0P(n1,1)*QGap0P(n2+n3+n4,3)-6.*QGap0P(n1+n2+n3+n4,4);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FourGap2(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap2P(n1,1)*QGap2P(n2,1)*QGap2M(n3,1)*QGap2M(n4,1)-QGap2P(n1+n2,2)*QGap2M(n3,1)*QGap2M(n4,1)
		-QGap2P(n1,1)*QGap2P(n2,1)*QGap2M(n3+n4,2)+QGap2P(n1+n2,2)*QGap2M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FourGap4(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap4P(n1,1)*QGap4P(n2,1)*QGap4M(n3,1)*QGap4M(n4,1)-QGap4P(n1+n2,2)*QGap4M(n3,1)*QGap4M(n4,1)
		-QGap4P(n1,1)*QGap4P(n2,1)*QGap4M(n3+n4,2)+QGap4P(n1+n2,2)*QGap4M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FourGap6(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap6P(n1,1)*QGap6P(n2,1)*QGap6M(n3,1)*QGap6M(n4,1)-QGap6P(n1+n2,2)*QGap6M(n3,1)*QGap6M(n4,1)
		-QGap6P(n1,1)*QGap6P(n2,1)*QGap6M(n3+n4,2)+QGap6P(n1+n2,2)*QGap6M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FourGap8(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap8P(n1,1)*QGap8P(n2,1)*QGap8M(n3,1)*QGap8M(n4,1)-QGap8P(n1+n2,2)*QGap8M(n3,1)*QGap8M(n4,1)
		-QGap8P(n1,1)*QGap8P(n2,1)*QGap8M(n3+n4,2)+QGap8P(n1+n2,2)*QGap8M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FourGap10(int n1, int n2, int n3, int n4)
{

	TComplex formula = QGap10P(n1,1)*QGap10P(n2,1)*QGap10M(n3,1)*QGap10M(n4,1)-QGap10P(n1+n2,2)*QGap10M(n3,1)*QGap10M(n4,1)
		-QGap10P(n1,1)*QGap10P(n2,1)*QGap10M(n3+n4,2)+QGap10P(n1+n2,2)*QGap10M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Four_3SubMMLR(int n1, int n2, int n3, int n4)
{
	TComplex formula = QsubMiddle(n1,1)*QsubMiddle(n2,1)*QsubLeft(n3,1)*QsubRight(n4,1)-QsubMiddle(n1+n2,2)*QsubLeft(n3,1)*QsubRight(n4,1);
	return formula;
}

//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Four_3SubLLMR(int n1, int n2, int n3, int n4)
{
	TComplex formula = QsubLeft(n1,1)*QsubLeft(n2,1)*QsubMiddle(n3,1)*QsubRight(n4,1)-QsubLeft(n1+n2,2)*QsubMiddle(n3,1)*QsubRight(n4,1);
	return formula;
}

//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Four_3SubRRML(int n1, int n2, int n3, int n4)
{
	TComplex formula = QsubRight(n1,1)*QsubRight(n2,1)*QsubMiddle(n3,1)*QsubLeft(n4,1)-QsubRight(n1+n2,2)*QsubMiddle(n3,1)*QsubLeft(n4,1);
	return formula;
}

//____________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Four_3SubGap2Evts(int n1, int n2, int n3, int n4)
{
	TComplex formula = QsubGap2Middle(n1,1)*QsubGap2Middle(n2,1)*QsubGap2Left(n3,1)*QsubGap2Right(n4,1)
		-QsubGap2Middle(n1+n2,2)*QsubGap2Left(n3,1)*QsubGap2Right(n4,1);
	return formula;
}
//___________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Five(int n1, int n2, int n3, int n4, int n5)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)
		- Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)
		+ 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)
		+ Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)
		+ Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)
		- Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)
		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)
		- 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)
		+ Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)
		+ Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)
		- Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)
		+ Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)
		- 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)
		- 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)
		+ Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)
		+ Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)
		+ 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)
		+ 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)
		- 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)
		+ Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)
		+ Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)
		+ 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)
		+ 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)
		- 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)
		- 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)
		- 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)+24.*Q(n1+n2+n3+n4+n5,5);
	return formula;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FiveGap0A(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	TComplex formula = TwoGap0M(n1, n2) * ThreeGap0_subP(n3, n4, n5);
	//TComplex formula = (QGap0M(n1,1)*QGap0M(n2,1) - QGap0M(n1+n2,2)) * (QGap0P(n3,1)*QGap0P(n4,1)*QGap0P(n5,1)-QGap0P(n3+n4,2)*QGap0P(n5,1)-QGap0P(n4,1)*QGap0P(n3+n5,2) - QGap0P(n3,1)*QGap0P(n4+n5,2)+2.*QGap0P(n3+n4+n5,3));
	return formula;

}

//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FiveGap0A_2(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	//TComplex formula = TwoGap0M(n1, n2) * ThreeGap0_subP(n3, n4, n5);
	TComplex formula = (QGap0M(n1,1)*QGap0M(n2,1) - QGap0M(n1+n2,2)) * (QGap0P(n3,1)*QGap0P(n4,1)*QGap0P(n5,1)-QGap0P(n3+n4,2)*QGap0P(n5,1)-QGap0P(n4,1)*QGap0P(n3+n5,2) - QGap0P(n3,1)*QGap0P(n4+n5,2)+2.*QGap0P(n3+n4+n5,3));
	return formula;

}

//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FiveGap0B(int n1, int n2, int n3, int n4, int n5) // (- - -, + +)
{

	TComplex formula = TwoGap0P(n1, n2) * ThreeGap0_subM(n3, n4, n5);
	return formula;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FiveGap2A(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	TComplex formula = TwoGap2M(n1, n2) * ThreeGap2_subP(n3, n4, n5);
	return formula;

}

//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FiveGap2B(int n1, int n2, int n3, int n4, int n5) // (- - -, + +)
{

	TComplex formula = TwoGap2P(n1, n2) * ThreeGap2_subM(n3, n4, n5);
	return formula;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FiveGap4A(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	TComplex formula = TwoGap4M(n1, n2) * ThreeGap4_subP(n3, n4, n5);
	return formula;

}

//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FiveGap4B(int n1, int n2, int n3, int n4, int n5) // (- - -, + +)
{

	TComplex formula = TwoGap4P(n1, n2) * ThreeGap4_subM(n3, n4, n5);
	return formula;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FiveGap6A(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	TComplex formula = TwoGap6M(n1, n2) * ThreeGap6_subP(n3, n4, n5);
	return formula;

}

//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FiveGap6B(int n1, int n2, int n3, int n4, int n5) // (- - -, + +)
{

	TComplex formula = TwoGap6P(n1, n2) * ThreeGap6_subM(n3, n4, n5);
	return formula;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FiveGap8A(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	TComplex formula = TwoGap8M(n1, n2) * ThreeGap8_subP(n3, n4, n5);
	return formula;

}

//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FiveGap8B(int n1, int n2, int n3, int n4, int n5) // (- - -, + +)
{

	TComplex formula = TwoGap8P(n1, n2) * ThreeGap8_subM(n3, n4, n5);
	return formula;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FiveGap10A(int n1, int n2, int n3, int n4, int n5) // ( + +,- - - )
{

	TComplex formula = TwoGap10M(n1, n2) * ThreeGap10_subP(n3, n4, n5);
	return formula;

}

//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::FiveGap10B(int n1, int n2, int n3, int n4, int n5) // (- - -, + +)
{

	TComplex formula = TwoGap10P(n1, n2) * ThreeGap10_subM(n3, n4, n5);
	return formula;

}
//___________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Six(int n1, int n2, int n3, int n4, int n5, int n6)
{


	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
		- Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
		+ 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)
		+ Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
		+ Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
		- Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
		- 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)
		+ Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)
		+ Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)
		- Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
		+ Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
		- 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)
		- 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
		+ Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
		+ Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
		+ 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)
		+ 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
		- 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
		+ Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
		+ Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
		+ 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)
		+ 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
		- 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
		- 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)
		- 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)
		- Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)
		+ Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)
		- 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)
		- Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)
		- Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)
		+ Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)
		- 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)
		+ 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
		+ Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
		+ Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
		+ Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)
		+ Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
		- 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
		- Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)
		- 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)
		+ 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)
		- 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)
		+ 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
		+ Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
		+ Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
		+ Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)
		+ Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
		- 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
		- Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)
		- 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)
		+ 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)
		- 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)
		+ 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
		- 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)
		- 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)
		- 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)
		- Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
		+ Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
		- 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)
		- Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
		- Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)
		+ Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
		- 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
		+ 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)
		- 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)
		- 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)
		+ 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
		- 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
		+ 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)
		+ 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
		- 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)
		- 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)
		- 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)
		- 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)
		+ 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
		+ Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
		+ Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
		+ Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)
		+ Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
		- 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
		- Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)
		- 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)
		+ 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)
		- 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)
		+ 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
		- 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)
		- 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)
		- 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)
		+ 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
		- 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
		+ 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)
		+ 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
		+ 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)
		+ 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
		- 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
		+ 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)
		+ 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
		+ 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)
		- 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)
		+ 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)
		- 120.*Q(n1+n2+n3+n4+n5+n6,6);
	return formula;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::SixGap0(int n1, int n2, int n3, int n4, int n5, int n6)
{

	TComplex formula = ThreeGap0_subM(n1, n2, n3)*ThreeGap0_subP(n4, n5, n6);
	return formula;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::SixGap2(int n1, int n2, int n3, int n4, int n5, int n6)
{

	TComplex formula = ThreeGap2_subM(n1, n2, n3)*ThreeGap2_subP(n4, n5, n6);
	return formula;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::SixGap4(int n1, int n2, int n3, int n4, int n5, int n6)
{

	TComplex formula = ThreeGap4_subM(n1, n2, n3)*ThreeGap4_subP(n4, n5, n6);
	return formula;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::SixGap6(int n1, int n2, int n3, int n4, int n5, int n6)
{

	TComplex formula = ThreeGap6_subM(n1, n2, n3)*ThreeGap6_subP(n4, n5, n6);
	return formula;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::SixGap8(int n1, int n2, int n3, int n4, int n5, int n6)
{

	TComplex formula = ThreeGap8_subM(n1, n2, n3)*ThreeGap8_subP(n4, n5, n6);
	return formula;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::SixGap10(int n1, int n2, int n3, int n4, int n5, int n6)
{

	TComplex formula = ThreeGap10_subM(n1, n2, n3)*ThreeGap10_subP(n4, n5, n6);
	return formula;

}
//_________________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7)
{

	TComplex Correlation = {0, 0};
	int Narray[] = {n1, n2, n3, n4, n5, n6};

	for(int k=7; k-->0; )
	{// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

		int array[6] = {0,1,2,3,4,5};
		int iPerm = 0;
		int argument = 0;
		int count = 0;

		// k==6: there is just one combination, we can add it manually
		if(k==6){
			Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
				Six(n1, n2, n3, n4, n5, n6)*Q(n7, 7-k);
		}// k==6

		else if(k==5){
			do{
				iPerm += 1;
				if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
					count += 1;
					Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
						Five(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
								Narray[int(array[3])], Narray[int(array[4])])*
						Q(Narray[int(array[5])]+n7, 7-k);
				}
			}while(std::next_permutation(array, array+6));
		}// k==5

		else if(k==4){
			do{
				iPerm += 1;
				if(iPerm%2 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
							Four(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
									Narray[int(array[3])])*
							Q(Narray[int(array[4])]+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==4

		else if(k==3){
			do{
				iPerm += 1;
				if(iPerm%6 == 1){
					if(array[0] < array[1] && array[1] < array[2]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
							Three(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])])*
							Q(Narray[int(array[3])]+Narray[int(array[4])]+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==3

		else if(k==2){
			do{
				iPerm += 1;
				if(iPerm%24 == 1){
					if(array[0] < array[1]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
							Two(Narray[int(array[0])], Narray[int(array[1])])*
							Q(Narray[int(array[2])]+Narray[int(array[3])]+Narray[int(array[4])]
									+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==2

		else if(k == 1){
			Correlation = Correlation
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7, 7-k)
				+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7, 7-k);
		}// k==1

		else if(k == 0){
			Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1+n2+n3+n4+n5+n6+n7, 7-k);
		}// k==0

		else{
			cout<<"invalid range of k"<<endl;
			return {0,0};
		}

	}// loop over k

	return Correlation;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{

	TComplex Correlation = {0, 0};
	int Narray[] = {n1, n2, n3, n4, n5, n6, n7};

	for(int k=8; k-->0; )
	{// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

		int array[7] = {0,1,2,3,4,5,6};
		int iPerm = 0;
		int argument = 0;
		int count = 0;

		// k==7: there is just one combination, we can add it manually
		if(k==7){
			Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
				Seven(n1, n2, n3, n4, n5, n6, n7)*Q(n8, 8-k);
		}// k==7

		else if(k==6){
			do{
				iPerm += 1;
				if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4] && array[4] < array[5]){
					count += 1;
					Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
						Six(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
								Narray[int(array[3])], Narray[int(array[4])], Narray[int(array[5])])*
						Q(Narray[int(array[6])]+n8, 8-k);
				}
			}while(std::next_permutation(array, array+7));
		}// k==6

		else if(k==5){
			do{
				iPerm += 1;
				if(iPerm%2 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Five(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
									Narray[int(array[3])], Narray[int(array[4])])*
							Q(Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==5

		else if(k==4){
			do{
				iPerm += 1;
				if(iPerm%6 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Four(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])], Narray[int(array[3])])*
							Q(Narray[int(array[4])]+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==4

		else if(k==3){
			do{
				iPerm += 1;
				if(iPerm%24 == 1){
					if(array[0] < array[1] && array[1] < array[2]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Three(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])])*
							Q(Narray[int(array[3])]+Narray[int(array[4])]+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==3

		else if(k==2){
			do{
				iPerm += 1;
				if(iPerm%120 == 1){
					if(array[0] < array[1]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
							Two(Narray[int(array[0])], Narray[int(array[1])])*
							Q(Narray[int(array[2])]+Narray[int(array[3])]+Narray[int(array[4])]
									+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==2

		else if(k == 1){
			Correlation = Correlation
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7+n8, 8-k)
				+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n7, 1)*Q(n1+n2+n3+n4+n5+n6+n8, 8-k);
		}// k==1

		else if(k == 0){
			Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1+n2+n3+n4+n5+n6+n7+n8, 8-k);
		}// k==0

		else{
			cout<<"invalid range of k"<<endl;
			return {0,0};
		}

	}// loop over k

	return Correlation;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskNonlinearFlow::EightGap0(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{

	TComplex formula = FourGap0M(n1, n2, n3, n4)*FourGap0P(n5, n6, n7, n8);
	return formula;

}
//_____________________________________________________________________________
Short_t AliAnalysisTaskNonlinearFlow::GetCentrCode(AliVEvent* ev)
{
	Short_t centrCode = -1;
	Float_t lPercentile = 0;
	Float_t V0M_Cent = 0, SPD_Cent = 0;

	//if (fAnalysisType == "AOD"){
	//AliAODEvent* aod = (AliAODEvent*)ev;
	AliMultSelection *MultSelection = 0;
	MultSelection = (AliMultSelection * ) ev->FindListObject("MultSelection");
	if(!MultSelection){
		lPercentile = -100;
	}
	else{
		if (fCentFlag == 0)
			lPercentile = MultSelection->GetMultiplicityPercentile("V0M");

		if (fCentFlag == 1)
			lPercentile = MultSelection->GetMultiplicityPercentile("CL0");

		if (fCentFlag == 2)
			lPercentile = MultSelection->GetMultiplicityPercentile("CL1");

		V0M_Cent = MultSelection->GetMultiplicityPercentile("V0M");
		SPD_Cent = MultSelection->GetMultiplicityPercentile("CL1");

	}

	if ((lPercentile > 0) && (lPercentile <= 5.0))
		centrCode = 0;
	else if ((lPercentile > 5.0) && (lPercentile <= 10.0))
		centrCode = 1;
	else if ((lPercentile > 10.0) && (lPercentile <= 20.0))
		centrCode = 2;
	else if ((lPercentile > 20.0) && (lPercentile <= 30.0))
		centrCode = 3;
	else if ((lPercentile > 30.0) && (lPercentile <= 40.0))
		centrCode = 4;
	else if ((lPercentile > 40.0) && (lPercentile <= 50.0))
		centrCode = 5;
	else if ((lPercentile > 50.0) && (lPercentile <= 60.0))
		centrCode = 6;
	else if ((lPercentile > 60.0) && (lPercentile <= 70.0))
		centrCode = 7;
	else if ((lPercentile > 70.0) && (lPercentile <= 80.0))
		centrCode = 8;
	else if ((lPercentile > 80.0) && (lPercentile <= 90.0))
		centrCode = 9;


	return centrCode;

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
	double Dn2 = Two(0, 0).Re();
	double Dn2Gap10 = TwoGap10(0, 0).Re();
	double Dn2Gap14 = TwoGap14(0, 0).Re();
	double Dn2_3subLM = Two_3SubLM(0, 0).Re();
	double Dn2_3subRM = Two_3SubRM(0, 0).Re();
	double Dn2_3subLR = Two_3SubRM(0, 0).Re();

	if(NtrksAfter > 1 && Dn2 != 0) {
		//..v2{2} = <cos2(phi1 - phi2)>
		TComplex v22 = Two(2, -2);
		double v22Re = v22.Re()/Dn2;
		profile.fChcn2[0]->Fill(Ntrks, v22Re, Dn2);

		//..v3{2} = <cos3(phi1 - phi2)>
		TComplex v32 = Two(3, -3);
		double v32Re = v32.Re()/Dn2;
		profile.fChcn2[1]->Fill(Ntrks, v32Re, Dn2);

		//..v4{2} = <cos4(phi1 - phi2)>
		TComplex v42 = Two(4, -4);
		double v42Re = v42.Re()/Dn2;
		profile.fChcn2[2]->Fill(Ntrks, v42Re, Dn2);

	}

	if(NtrksAfterGap10M > 0 && NtrksAfterGap10P > 0 && Dn2Gap10 != 0)
	{
		//..v2{2} with eta Gap > 1.0
		TComplex v22Gap10 = TwoGap10(2, -2);
		double v22ReGap10 = v22Gap10.Re()/Dn2Gap10;
		profile.fChcn2_Gap10[0]->Fill(Ntrks, v22ReGap10, Dn2Gap10);

		//..v3{2} with eta Gap > 1.0
		TComplex v32Gap10 = TwoGap10(3, -3);
		double v32ReGap10 = v32Gap10.Re()/Dn2Gap10;
		profile.fChcn2_Gap10[1]->Fill(Ntrks, v32ReGap10, Dn2Gap10);

		//..v4{2} with eta Gap > 1.0
		TComplex v42Gap10 = TwoGap10(4, -4);
		double v42ReGap10 = v42Gap10.Re()/Dn2Gap10;
		profile.fChcn2_Gap10[2]->Fill(Ntrks, v42ReGap10, Dn2Gap10);
	}

	if(NtrksAfterGap14M > 0 && NtrksAfterGap14P > 0 && Dn2Gap14 != 0)
	{
		//..v2{2} with eta Gap > 1.4
		TComplex v22Gap14 = TwoGap14(2, -2);
		double v22ReGap14 = v22Gap14.Re()/Dn2Gap14;
		profile.fChcn2_Gap14[0]->Fill(Ntrks, v22ReGap14, Dn2Gap14);

		//..v3{2} with eta Gap > 1.4
		TComplex v32Gap14 = TwoGap14(3, -3);
		double v32ReGap14 = v32Gap14.Re()/Dn2Gap14;
		profile.fChcn2_Gap14[1]->Fill(Ntrks, v32ReGap14, Dn2Gap14);

		//..v4{2} with eta Gap > 1.4
		TComplex v42Gap14 = TwoGap14(4, -4);
		double v42ReGap14 = v42Gap14.Re()/Dn2Gap14;
		profile.fChcn2_Gap14[2]->Fill(Ntrks, v42ReGap14, Dn2Gap14);
	}

	//..for 3-subevent method, Gap0
	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 0 && Dn2_3subLM != 0)
	{//..left+middle
		TComplex v22_3subLM = Two_3SubLM(2, -2);
		double v22Re_3subLM = v22_3subLM.Re()/Dn2_3subLM;
		profile.fChcn2_3subLM[0]->Fill(Ntrks, v22Re_3subLM, Dn2_3subLM);

		TComplex v32_3subLM = Two_3SubLM(3, -3);
		double v32Re_3subLM = v32_3subLM.Re()/Dn2_3subLM;
		profile.fChcn2_3subLM[1]->Fill(Ntrks, v32Re_3subLM, Dn2_3subLM);

		TComplex v42_3subLM = Two_3SubLM(4, -4);
		double v42Re_3subLM = v42_3subLM.Re()/Dn2_3subLM;
		profile.fChcn2_3subLM[2]->Fill(Ntrks, v42Re_3subLM, Dn2_3subLM);
	}

	if(NtrksAfter3subM > 0 && NtrksAfter3subR > 0 && Dn2_3subRM != 0)
	{//..right+middle
		TComplex v22_3subRM = Two_3SubRM(2, -2);
		double v22Re_3subRM = v22_3subRM.Re()/Dn2_3subRM;
		profile.fChcn2_3subRM[0]->Fill(Ntrks, v22Re_3subRM, Dn2_3subRM);

		TComplex v32_3subRM = Two_3SubRM(3, -3);
		double v32Re_3subRM = v32_3subRM.Re()/Dn2_3subRM;
		profile.fChcn2_3subRM[1]->Fill(Ntrks, v32Re_3subRM, Dn2_3subRM);

		TComplex v42_3subRM = Two_3SubRM(4, -4);
		double v42Re_3subRM = v42_3subRM.Re()/Dn2_3subRM;
		profile.fChcn2_3subRM[2]->Fill(Ntrks, v42Re_3subRM, Dn2_3subRM);
	}

	if(NtrksAfter3subL > 0 && NtrksAfter3subR > 0 && Dn2_3subLR != 0)
	{//..right+middle
		TComplex v22_3subLR = Two_3SubLR(2, -2);
		double v22Re_3subLR = v22_3subLR.Re()/Dn2_3subLR;
		profile.fChcn2_3subLR[0]->Fill(Ntrks, v22Re_3subLR, Dn2_3subLR);

		TComplex v32_3subLR = Two_3SubLR(3, -3);
		double v32Re_3subLR = v32_3subLR.Re()/Dn2_3subLR;
		profile.fChcn2_3subLR[1]->Fill(Ntrks, v32Re_3subLR, Dn2_3subLR);

		TComplex v42_3subLR = Two_3SubLR(4, -4);
		double v42Re_3subLR = v42_3subLR.Re()/Dn2_3subLR;
		profile.fChcn2_3subLR[2]->Fill(Ntrks, v42Re_3subLR, Dn2_3subLR);
	}

	//..calculate 3-particle correlations
	//................................
	double Dn3 = Three(0, 0, 0).Re();
	double Dn3Gap10A = ThreeGap10A(0, 0, 0).Re();
	double Dn3Gap10B = ThreeGap10B(0, 0, 0).Re();

	if(NtrksAfter > 2 && Dn3 != 0)
	{
		//..v4{psi2}
		TComplex v422 = Three(4, -2, -2);
		double v422Re = v422.Re()/Dn3;
		profile.fChc422->Fill(Ntrks, v422Re, Dn3);

		//..v5{psi32}
		TComplex v532 = Three(5, -3, -2);
		double v532Re = v532.Re()/Dn3;
		profile.fChc532->Fill(Ntrks, v532Re, Dn3);
	}

        // A-type
	if(NtrksAfterGap10M > 0 && NtrksAfterGap10P > 1 && Dn3Gap10A != 0)
	{

		TComplex v422Gap10A = ThreeGap10A(4, -2, -2);
		double v422Gap10ARe = v422Gap10A.Re()/Dn3Gap10A;
		profile.fChc422_Gap10A->Fill(Ntrks, v422Gap10ARe, Dn3Gap10A);

		TComplex v532Gap10A = ThreeGap10A(5, -3, -2);
		double v532Gap10ARe = v532Gap10A.Re()/Dn3Gap10A;
		profile.fChc532_Gap10A->Fill(Ntrks, v532Gap10ARe, Dn3Gap10A);
	}

        // B-type
	if(NtrksAfterGap10P > 0 && NtrksAfterGap10M > 1 && Dn3Gap10B != 0)
	{

		TComplex v422Gap10B = ThreeGap10B(4, -2, -2);
		double v422Gap10BRe = v422Gap10B.Re()/Dn3Gap10B;
		profile.fChc422_Gap10B->Fill(Ntrks, v422Gap10BRe, Dn3Gap10B);

		TComplex v532Gap10B = ThreeGap10B(5, -3, -2);
		double v532Gap10BRe = v532Gap10B.Re()/Dn3Gap10B;
		profile.fChc532_Gap10B->Fill(Ntrks, v532Gap10BRe, Dn3Gap10B);
	}

	//..calculate 4-particle correlations
	//................................
	double Dn4 = Four(0, 0, 0, 0).Re();
	double Dn4Gap10 = FourGap10(0, 0, 0, 0).Re();
	double Dn4_3subMMLR = Four_3SubMMLR(0, 0, 0, 0).Re();
	double Dn4_3subLLMR = Four_3SubLLMR(0, 0, 0, 0).Re();
	double Dn4_3subRRML = Four_3SubRRML(0, 0, 0, 0).Re();

	if(NtrksAfter > 3 && Dn4 != 0)
	{

		TComplex v24 = Four(2, 2, -2, -2);
		double v24Re = v24.Re()/Dn4;
		profile.fChcn4[0]->Fill(Ntrks, v24Re, Dn4);
		// fcn4Ntrks1bin[0][fBin]->Fill(NtrksAfter, v24Re, Dn4);

		TComplex v34 = Four(3, 3, -3, -3);
		double v34Re = v34.Re()/Dn4;
		profile.fChcn4[1]->Fill(Ntrks, v34Re, Dn4);
		// fcn4Ntrks1bin[1][fBin]->Fill(NtrksAfter, v34Re, Dn4);

		TComplex v44 = Four(4, 4, -4, -4);
		double v44Re = v44.Re()/Dn4;
		profile.fChcn4[2]->Fill(Ntrks, v44Re, Dn4);
		// fcn4Ntrks1bin[2][fBin]->Fill(NtrksAfter, v44Re, Dn4);

		//..SC(3,2,-3,-2)
		TComplex sc3232 = Four(3, 2, -3, -2);
		double sc3232Re = sc3232.Re()/Dn4;
		profile.fChsc3232->Fill(Ntrks, sc3232Re, Dn4);

		//..SC(4,2,-4,-2)
		TComplex sc4242 = Four(4, 2, -4, -2);
		double sc4242Re = sc4242.Re()/Dn4;
		profile.fChsc4242->Fill(Ntrks, sc4242Re, Dn4);

	}

	if(NtrksAfterGap10M > 1 && NtrksAfterGap10P > 1 && Dn4Gap10 !=0)
	{
		TComplex v24Gap10 = FourGap10(2, 2, -2, -2);
		double v24Gap10Re = v24Gap10.Re()/Dn4Gap10;
		profile.fChcn4_Gap10[0]->Fill(Ntrks, v24Gap10Re, Dn4Gap10);

		TComplex v34Gap10 = FourGap10(3, 3, -3, -3);
		double v34Gap10Re = v34Gap10.Re()/Dn4Gap10;
		profile.fChcn4_Gap10[1]->Fill(Ntrks, v34Gap10Re, Dn4Gap10);

		TComplex v44Gap10 = FourGap10(4, 4, -4, -4);
		double v44Gap10Re = v44Gap10.Re()/Dn4Gap10;
		profile.fChcn4_Gap10[2]->Fill(Ntrks, v44Gap10Re, Dn4Gap10);

		TComplex sc3232Gap10 = FourGap10(3, 2, -3, -2);
		double sc3232Gap10Re = sc3232Gap10.Re()/Dn4Gap10;
		profile.fChsc3232_Gap10->Fill(Ntrks, sc3232Gap10Re, Dn4Gap10);

		TComplex sc4242Gap10 = FourGap10(4, 2, -4, -2);
		double sc4242Gap10Re = sc4242Gap10.Re()/Dn4Gap10;
		profile.fChsc4242_Gap10->Fill(Ntrks, sc4242Gap10Re, Dn4Gap10);
	}

	//..3-subevent method
	if(NtrksAfter3subL > 0 && NtrksAfter3subR > 0 && NtrksAfter3subM > 1 && Dn4_3subMMLR != 0)
	{
		TComplex v24_3sub = Four_3SubMMLR(2, 2, -2, -2);
		double v24_3subRe = v24_3sub.Re()/Dn4_3subMMLR;
		profile.fChcn4_3subMMLR[0]->Fill(Ntrks, v24_3subRe, Dn4_3subMMLR);

		TComplex v34_3sub = Four_3SubMMLR(3, 3, -3, -3);
		double v34_3subRe = v34_3sub.Re()/Dn4_3subMMLR;
		profile.fChcn4_3subMMLR[1]->Fill(Ntrks, v34_3subRe, Dn4_3subMMLR);

		TComplex v44_3sub = Four_3SubMMLR(4, 4, -4, -4);
		double v44_3subRe = v44_3sub.Re()/Dn4_3subMMLR;
		profile.fChcn4_3subMMLR[2]->Fill(Ntrks, v44_3subRe, Dn4_3subMMLR);

		TComplex sc3232_3subA = Four_3SubMMLR(3, 2, -3, -2);
		double sc3232_3subARe = sc3232_3subA.Re()/Dn4_3subMMLR;
		profile.fChsc3232_3subMMLRA->Fill(Ntrks, sc3232_3subARe, Dn4_3subMMLR);

		TComplex sc3232_3subB = Four_3SubMMLR(3, 2, -2, -3);
		double sc3232_3subBRe = sc3232_3subB.Re()/Dn4_3subMMLR;
		profile.fChsc3232_3subMMLRB->Fill(Ntrks, sc3232_3subBRe, Dn4_3subMMLR);

		TComplex sc4242_3subA = Four_3SubMMLR(4, 2, -4, -2);
		double sc4242_3subARe = sc4242_3subA.Re()/Dn4_3subMMLR;
		profile.fChsc4242_3subMMLRA->Fill(Ntrks, sc4242_3subARe, Dn4_3subMMLR);

		TComplex sc4242_3subB = Four_3SubMMLR(4, 2, -2, -4);
		double sc4242_3subBRe = sc4242_3subB.Re()/Dn4_3subMMLR;
		profile.fChsc4242_3subMMLRB->Fill(Ntrks, sc4242_3subBRe, Dn4_3subMMLR);
	}

	//..3-subevent method
	if(NtrksAfter3subL > 1 && NtrksAfter3subR > 0 && NtrksAfter3subM > 0 && Dn4_3subLLMR != 0)
	{
		TComplex v24_3sub = Four_3SubLLMR(2, 2, -2, -2);
		double v24_3subRe = v24_3sub.Re()/Dn4_3subLLMR;
		profile.fChcn4_3subLLMR[0]->Fill(Ntrks, v24_3subRe, Dn4_3subLLMR);

		TComplex v34_3sub = Four_3SubLLMR(3, 3, -3, -3);
		double v34_3subRe = v34_3sub.Re()/Dn4_3subLLMR;
		profile.fChcn4_3subLLMR[1]->Fill(Ntrks, v34_3subRe, Dn4_3subLLMR);

		TComplex v44_3sub = Four_3SubLLMR(4, 4, -4, -4);
		double v44_3subRe = v44_3sub.Re()/Dn4_3subLLMR;
		profile.fChcn4_3subLLMR[2]->Fill(Ntrks, v44_3subRe, Dn4_3subLLMR);

		TComplex sc3232_3subA = Four_3SubLLMR(3, 2, -3, -2);
		double sc3232_3subARe = sc3232_3subA.Re()/Dn4_3subLLMR;
		profile.fChsc3232_3subLLMRA->Fill(Ntrks, sc3232_3subARe, Dn4_3subLLMR);

		TComplex sc3232_3subB = Four_3SubLLMR(3, 2, -2, -3);
		double sc3232_3subBRe = sc3232_3subB.Re()/Dn4_3subLLMR;
		profile.fChsc3232_3subLLMRB->Fill(Ntrks, sc3232_3subBRe, Dn4_3subLLMR);

		TComplex sc4242_3subA = Four_3SubLLMR(4, 2, -4, -2);
		double sc4242_3subARe = sc4242_3subA.Re()/Dn4_3subLLMR;
		profile.fChsc4242_3subLLMRA->Fill(Ntrks, sc4242_3subARe, Dn4_3subLLMR);

		TComplex sc4242_3subB = Four_3SubLLMR(4, 2, -2, -4);
		double sc4242_3subBRe = sc4242_3subB.Re()/Dn4_3subLLMR;
		profile.fChsc4242_3subLLMRB->Fill(Ntrks, sc4242_3subBRe, Dn4_3subLLMR);
	}

	//..3-subevent method
	if(NtrksAfter3subL > 0 && NtrksAfter3subR > 1 && NtrksAfter3subM > 0 && Dn4_3subRRML != 0)
	{
		TComplex v24_3sub = Four_3SubRRML(2, 2, -2, -2);
		double v24_3subRe = v24_3sub.Re()/Dn4_3subRRML;
		profile.fChcn4_3subRRML[0]->Fill(Ntrks, v24_3subRe, Dn4_3subRRML);

		TComplex v34_3sub = Four_3SubRRML(3, 3, -3, -3);
		double v34_3subRe = v34_3sub.Re()/Dn4_3subRRML;
		profile.fChcn4_3subRRML[1]->Fill(Ntrks, v34_3subRe, Dn4_3subRRML);

		TComplex v44_3sub = Four_3SubRRML(4, 4, -4, -4);
		double v44_3subRe = v44_3sub.Re()/Dn4_3subRRML;
		profile.fChcn4_3subRRML[2]->Fill(Ntrks, v44_3subRe, Dn4_3subRRML);

		TComplex sc3232_3subA = Four_3SubRRML(3, 2, -3, -2);
		double sc3232_3subARe = sc3232_3subA.Re()/Dn4_3subRRML;
		profile.fChsc3232_3subRRMLA->Fill(Ntrks, sc3232_3subARe, Dn4_3subRRML);

		TComplex sc3232_3subB = Four_3SubRRML(3, 2, -2, -3);
		double sc3232_3subBRe = sc3232_3subB.Re()/Dn4_3subRRML;
		profile.fChsc3232_3subRRMLB->Fill(Ntrks, sc3232_3subBRe, Dn4_3subRRML);

		TComplex sc4242_3subA = Four_3SubRRML(4, 2, -4, -2);
		double sc4242_3subARe = sc4242_3subA.Re()/Dn4_3subRRML;
		profile.fChsc4242_3subRRMLA->Fill(Ntrks, sc4242_3subARe, Dn4_3subRRML);

		TComplex sc4242_3subB = Four_3SubRRML(4, 2, -2, -4);
		double sc4242_3subBRe = sc4242_3subB.Re()/Dn4_3subRRML;
		profile.fChsc4242_3subRRMLB->Fill(Ntrks, sc4242_3subBRe, Dn4_3subRRML);
	}

}

//_____________________________________________________________________________
void AliAnalysisTaskNonlinearFlow::Terminate(Option_t *)
{
	// Terminate loop
	Printf("Terminate()");
}
