/**************************************************************************
 * Copyright(c) 1998-2022, ALICE Experiment at CERN, All rights reserved. *
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

/*
   Xic0 -> eXi analysis (* newly written, Aug. 2022)

   Chong Kim
   Pusan National University
   kimc@cern.ch
*/

#include "AliAnalysisTaskSEXic0SL.h"
ClassImp(AliAnalysisTaskSEXic0SL);

#include "AliAnalysisManager.h"
#include "AliAODcascade.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"

#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliNormalizationCounter.h"
#include "AliPIDResponse.h"
#include "AliPPVsMultUtils.h"
#include "AliRDHFCutsXictoeleXifromAODtracks.h" //For fEvtCuts
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "THistManager.h" //For fHisto

#include "TBranch.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TLeaf.h"
#include "TString.h"
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

AliAnalysisTaskSEXic0SL::AliAnalysisTaskSEXic0SL(const char* name, const char* option):
	AliAnalysisTaskSE(name), fTaskOpt(option)
{
	ControlAnaObjects(0);
	ControlOutputContainers(0);
	ResetTreeVariables();
}//Constructor

AliAnalysisTaskSEXic0SL::~AliAnalysisTaskSEXic0SL()
{
	ControlAnaObjects(1);
}//Destructor

void AliAnalysisTaskSEXic0SL::Terminate(Option_t *)
{
	cout <<"\nDone!\n";
	return;
}//Terminate

//=======================================================================================
void AliAnalysisTaskSEXic0SL::UserCreateOutputObjects()
{
	if (IsCutsByFile)
	{
		cout <<"\nSet selection cuts by using external file...\n";

	}//IsCutsByFile
	else //Apply hard-corded cuts, includes legacy mode
	{
		//Eventwise cut, for MB (default)
		fEvtCutsMB = new AliRDHFCutsXictoeleXifromAODtracks();
		fEvtCutsMB->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent); //Multi vertexer pileup rejection
		if (IsMC==false) //data
		{
			//fEvtCutsMB->SetTriggerMask(AliVEvent::kINT7); //This is included in SetUseIny7TriggerPP2012()
			fEvtCutsMB->SetUseInt7TriggerPP2012();
		}
		else if (IsLegacy==false) fEvtCutsMB->SetTriggerMask(AliVEvent::kAny); //MC only
		fEvtCutsMB->SetUsePhysicsSelection(true);

		//Eventwise cut, for HMV0 trigger
		fEvtCutsHMV0 = new AliRDHFCutsXictoeleXifromAODtracks();
		fEvtCutsHMV0->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
		if (IsMC==false)
		{
			fEvtCutsHMV0->SetTriggerMask(AliVEvent::kHighMultV0);
			fEvtCutsHMV0->SetTriggerClass(""); //Causes problem on NormalizationCounter if deleted, Jan. 9
		}
		fEvtCutsHMV0->SetUsePhysicsSelection(true);

		//AliESDTrackCut for track filtering
		fTrkCuts = new AliESDtrackCuts();
		fTrkCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kBoth);
		fTrkCuts->SetDCAToVertex2D(true);
		fTrkCuts->SetMaxChi2PerClusterITS(cut_maxChi2PerClusterITS);
		fTrkCuts->SetRequireTPCRefit(true); //TPC refit
		fTrkCuts->SetRequireITSRefit(true); //ITS refit

		if (IsLegacy)
		{
			cout <<"\nWARNING! Legacy conditions are being used!\n";

			fEvtCutsMB->SetTriggerClass("");
			fTrkCuts->SetMaxDCAToVertexXY(cut_maxDCAToVertexXY);
			fTrkCuts->SetMaxDCAToVertexZ(cut_maxDCAToVertexZ);
			fTrkCuts->SetMinNClustersITS(cut_minNClustersITS); //ITS cluster
			//fTrkCuts->SetMinNClustersTPC(fNClustersTPCMin = 70); //TPC min # of clusters - why masked?
		}//IsLegacy
	}//Use Hard-corded cuts

	//Trigger setup
	TrigStore.clear();
	if (TrigOnMB) TrigStore.push_back(AliVEvent::kINT7);
	if (TrigOnHMV0) TrigStore.push_back(AliVEvent::kHighMultV0);
	if (TrigStore.size() == 0) AliFatal("ERROR: no trigger is enabled!");
	for (unsigned int a=0; a<TrigStore.size(); a++)
	{
		if (TrigStore[a] == AliVEvent::kINT7) AliInfo(Form("Adding trigger: kINT7 (bit %i)", TrigStore[a]));
		if (TrigStore[a] == AliVEvent::kHighMultV0) AliInfo(Form("Adding trigger: kHMV0 (bit %i)\n", TrigStore[a]));
	}

	//-----------------------------------------------------

	//Histograms
	fHisto = new THistManager("Histo");
	fHisto->CreateTH1("RunNo", ";run", cut_runNoUp - cut_runNoLo, cut_runNoLo, cut_runNoUp, "s");

	vector<const char*> evtCuts = {"Bfield", "Mult", "PID", "Trig>0", "PileupX>0","VtxZ10", "INEL>0"};
	TH1* H1_cutEff = fHisto->CreateTH1("EvtCut", ";cut", evtCuts.size(), 0, (float)evtCuts.size(), "s");
	for (unsigned int a=0; a<evtCuts.size(); a++) H1_cutEff->GetXaxis()->SetBinLabel(a+1, evtCuts[a]);

	vector<const char*> evtTrigList = {"Any", "MB", "HMV0"};
	TH1* H1_trig = fHisto->CreateTH1("EvtTrig", ";trig", evtTrigList.size(), 0, (float)evtTrigList.size());
	for (unsigned int a=0; a<evtTrigList.size(); a++) H1_trig->GetXaxis()->SetBinLabel(a+1, evtTrigList[a]);

	//xChceks for raw e/Xi yields (after filtering), Nov. 2022
	fHisto->CreateTH1("e_minMassUS", "", 300, 0, 3, "s");
	fHisto->CreateTH1("c_massXi", "", 400, 1.1, 1.5, "s");
	fHisto->CreateTH1("c_massXi_MB", "", 400, 1.1, 1.5, "s");
	fHisto->CreateTH1("c_massXi_HMV0", "", 400, 1.1, 1.5, "s");

	//Tree
	fTree = new TTree("T", Form("Xic0SemiLeptonic%s%s%s", IsPA?"_pA":"", IsMC?"_MC":"", IsLegacy?"_LEGACY":""));
	ControlOutputTree(fTree, IsMC);

	//AliNormalizationCounter
	fANC_MB_0to100   = new AliNormalizationCounter("ANC_MB_0to100");
	fANC_MB_30to100  = new AliNormalizationCounter("ANC_MB_30to100");
	fANC_MB_0p1to30  = new AliNormalizationCounter("ANC_MB_0p1to30");
	fANC_HMV0_0to0p1 = new AliNormalizationCounter("ANC_HMV0_0to0p1");
	fANC_MB_0to100  ->SetStudyMultiplicity(true, 1.); fANC_MB_0to100  ->Init();
	fANC_MB_30to100 ->SetStudyMultiplicity(true, 1.); fANC_MB_30to100 ->Init();
	fANC_MB_0p1to30 ->SetStudyMultiplicity(true, 1.); fANC_MB_0p1to30 ->Init();
	fANC_HMV0_0to0p1->SetStudyMultiplicity(true, 1.); fANC_HMV0_0to0p1->Init();

	fANC_INEL0_MB_0to100   = new AliNormalizationCounter("ANC_INEL0_MB_0to100");
	fANC_INEL0_MB_30to100  = new AliNormalizationCounter("ANC_INEL0_MB_30to100");
	fANC_INEL0_MB_0p1to30  = new AliNormalizationCounter("ANC_INEL0_MB_0p1to30");
	fANC_INEL0_HMV0_0to0p1 = new AliNormalizationCounter("ANC_INEL0_HMV0_0to0p1");
	fANC_INEL0_MB_0to100  ->SetStudyMultiplicity(true, 1.); fANC_INEL0_MB_0to100  ->Init();
	fANC_INEL0_MB_30to100 ->SetStudyMultiplicity(true, 1.); fANC_INEL0_MB_30to100 ->Init();
	fANC_INEL0_MB_0p1to30 ->SetStudyMultiplicity(true, 1.); fANC_INEL0_MB_0p1to30 ->Init();
	fANC_INEL0_HMV0_0to0p1->SetStudyMultiplicity(true, 1.); fANC_INEL0_HMV0_0to0p1->Init();

	ControlOutputContainers(1);
	return;
}//UserCreateOutputObjects

//=======================================================================================
void AliAnalysisTaskSEXic0SL::UserExec(Option_t *)
{
	fEvtID++;
	ResetTreeVariables();

	//Eventwise selection starts
	//-----------------------------------------------------

	//Check AOD event object
	AliVEvent *event = InputEvent();
	if (event->IsA() == AliAODEvent::Class()) fEvt = dynamic_cast<AliAODEvent*>(event);
	else AliFatal("ERROR: input event type is NOT an AOD"); //Return 1

	//Check bfield
	const Double_t t_bfield = (Double_t)fEvt->GetMagneticField();
	if (fabs(t_bfield) < cut_bfield) AliFatal("ERROR: B-field?"); //Return 2
	fHisto->FillTH1("EvtCut", "Bfield", 1); //@

	//Check multiplicity selection is available
	fMultSel = (AliMultSelection*)fEvt->FindListObject("MultSelection");
	if (fMultSel) fEvtMult = fMultSel->GetMultiplicityPercentile(IsPA?"V0A":"V0M");
	else AliFatal("ERROR: MultSelection"); //Return 3
	fHisto->FillTH1("EvtCut", "Mult", 1); //@

	//Get input handler, Check event validity (* Dec. 21, outdated: does nothing, but keep it as legacy)
	fInputHandler = (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
	if (IsLegacy && fInputHandler->IsEventSelected()==false) return; //Return 4

	//Check pID
	fPID = (AliPIDResponse*)fInputHandler->GetPIDResponse();
	if (!fPID) AliFatal("ERROR: PID"); //Return 5
	fHisto->FillTH1("EvtCut", "PID", 1); //@

	//Check run number
	fEvtRunNo = fEvt->GetRunNumber();
	fHisto->FillTH1("RunNo", fEvtRunNo, 1); //@

	//Separate data/MC
	if (IsMC == false) //data
	{
		//Check vtx validity (NOT vtxZ)
		fEvtVtx = fEvt->GetPrimaryVertex();
		if (!fEvtVtx || fEvtVtx->GetNContributors()<cut_vtxNContributors) return; //Return 6d

		//Check trigger
		bool t_trigFired = false;
		fEvtTrig = fInputHandler->IsEventSelected();
		for (unsigned int a=0; a<TrigStore.size(); a++)
		{
			if (fEvtTrig & TrigStore[a]) t_trigFired = true;

			//Outdated but keep it as a legacy (Dec., 2022)
			//LHC16k and LHC16l are CD dedicated runs! CD-online-trigger is used
			if ( IsLegacy && (TrigStore[a] == AliVEvent::kINT7) &&
				 (fTaskOpt.Contains("LHC16k") || fTaskOpt.Contains("LHC16l")) &&
				 (fEvt->GetFiredTriggerClasses().Contains("CINT7-B-NOPF-CENT")) ) t_trigFired = true;
		}//a, TrigStore
		if (!t_trigFired) return; //Return 7d
		else fHisto->FillTH1("EvtCut", "Trig>0", 1); //@

		const bool t_trigFiredMB   = (fEvtTrig & AliVEvent::kINT7)?true:false;
		const bool t_trigFiredHMV0 = (fEvtTrig & AliVEvent::kHighMultV0)?true:false;
		fHisto->FillTH1("EvtTrig", "Any", 1); //@
		if (t_trigFiredMB)   fHisto->FillTH1("EvtTrig", "MB",   1); //@
		if (t_trigFiredHMV0) fHisto->FillTH1("EvtTrig", "HMV0", 1); //@

		//Activate RDHF cut by invoking it
		fEvtGoodMB   = fEvtCutsMB  ->IsEventSelected(fEvt);
		fEvtGoodHMV0 = fEvtCutsHMV0->IsEventSelected(fEvt);

		const bool t_pileup_MB   = fEvtCutsMB  ->IsEventRejectedDueToPileup(); //true = pile up
		const bool t_pileup_HMV0 = fEvtCutsHMV0->IsEventRejectedDueToPileup();
		if (t_pileup_MB != t_pileup_HMV0) AliFatal("Pileup status is different between MB and HMV0");
		else if (t_pileup_MB && t_pileup_HMV0) return;
		else fHisto->FillTH1("EvtCut", "PileupX", 1); //@

		//Check INEL>0
		if (AliPPVsMultUtils::IsINELgtZERO(fEvt)==true)	fEvtINELLgt0 = true;
		if (fEvtINELLgt0) fHisto->FillTH1("EvtCut", "INEL>0", 1); //@

		//Store AliNormalizationCounter: DO NOT CHANGE THIS POSITION! (after vtx found, before |vtxZ| < 10)
		if (t_trigFiredMB)
		{
			if (fEvtMult >= 0.0 && fEvtMult <= 100.0)
			{
				fANC_MB_0to100->StoreEvent(fEvt, fEvtCutsMB, IsMC);
				if (fEvtINELLgt0) fANC_INEL0_MB_0to100->StoreEvent(fEvt, fEvtCutsMB, IsMC);
			}
			if (fEvtMult >= 30.0 && fEvtMult <= 100.0)
			{
				fANC_MB_30to100->StoreEvent(fEvt, fEvtCutsMB, IsMC);
				if (fEvtINELLgt0) fANC_INEL0_MB_30to100->StoreEvent(fEvt, fEvtCutsMB, IsMC);
			}
			if (fEvtMult >= 0.1 && fEvtMult <= 30.0)
			{
				fANC_MB_0p1to30->StoreEvent(fEvt, fEvtCutsMB, IsMC);
				if (fEvtINELLgt0) fANC_INEL0_MB_0p1to30->StoreEvent(fEvt, fEvtCutsMB, IsMC);
			}
		}
		if (t_trigFiredHMV0 && (fEvtMult >= 0.0 && fEvtMult <= 0.1) )
		{
			fANC_HMV0_0to0p1->StoreEvent(fEvt, fEvtCutsHMV0, IsMC);
			if (fEvtINELLgt0) fANC_INEL0_HMV0_0to0p1->StoreEvent(fEvt, fEvtCutsHMV0, IsMC);
		}

		//Check vtxZ
		fEvtVtxZ = fEvtVtx->GetZ();
		if (fabs(fEvtVtxZ) > cut_vtxZ) return; //Return 9d
		else fHisto->FillTH1("EvtCut", "VtxZ10", 1); //@

		//Dec. 21, newly added
		if (IsLegacy==false && fEvtGoodMB==false && fEvtGoodHMV0==false) return; //Return 10d
	}//data
	else //MC
	{
		fMCEvt = fInputHandler->MCEvent();

		//Check vtxZ
		fEvtVtx = fMCEvt->GetPrimaryVertex();
		fEvtVtxZ = fEvtVtx->GetZ();
		if (fabs(fEvtVtxZ) > cut_vtxZ) return; //Return 7mc

		//RDHF event cut, Dec. 27, 2022
		if (IsLegacy==false && fEvtCutsMB->IsEventSelected(fEvt)==false) return; //Return 8mc

		//Check generated MC particles w/ eXi pairs
		int nTruth = 0;
		const int nTracksMC = fMCEvt->GetNumberOfTracks();
		for (int a=0; a<nTracksMC; a++)
		{
			//Possible 4-body decay contaminations:
			//a. Xic+ -> (e+) + (nu) + (Xi*0) -> (e+) + (nu) + (Xi- + pi+)
			//b. Xic0 -> (e-) + (nu_bar) + (Xi*+) -> (e-) + (nu_bar) + (Xi+ + pi0)
			fMCPart = (AliAODMCParticle*)fMCEvt->GetTrack(a);
			if (!fMCPart) continue; //Continue

			//Check if this is a desired particle
			bool t_proceed = false;
			if (IsLegacy)
			{
				if ( (abs(fMCPart->GetPdgCode()) == PDGCode_Xic0) ) t_proceed = true;
			}
			else //For later 4-body decay study
			{
				if ( (abs(fMCPart->GetPdgCode()) == PDGCode_Xic0) ||
					 (abs(fMCPart->GetPdgCode()) == PDGCode_Xicp) ) t_proceed = true;
			}
			if (t_proceed == false) continue; //Continue

			//Search e and Xi among daughters
			std::vector<int> idxEle;
			std::vector<int> idxXi;
			for (int b=fMCPart->GetDaughterFirst(); b<=fMCPart->GetDaughterLast(); b++)
			{
				if (b<0) { cout <<"ERROR: negative index found!\n"; break; }
				const int t_PDG = fMCEvt->GetTrack(b)->PdgCode();
				if (abs(t_PDG) == PDGCode_e) idxEle.push_back(b);
				if (abs(t_PDG) == PDGCode_Xi) idxXi.push_back(b);
			}//b
			if (idxEle.size()==0 || idxXi.size()==0) continue; //Continue

			//#
			fMCOrig [nTruth] = CheckOrigin(fMCEvt, fMCPart);
			fMCLabel[nTruth] = fMCPart->GetLabel();
			fMCPDG  [nTruth] = fMCPart->GetPdgCode();
			fMCPt   [nTruth] = fMCPart->Pt();
			fMCY    [nTruth] = fMCPart->Y();

			if (idxEle.size()>1 || idxXi.size()>1) AliWarning("More than one electron/Xi found in MC truth!");
			else
			{
				fMCElePt[nTruth] = fMCEvt->GetTrack(idxEle[0])->Pt();
				fMCEleY [nTruth] = fMCEvt->GetTrack(idxEle[0])->Y();
				fMCXiPt [nTruth] = fMCEvt->GetTrack(idxXi[0])->Pt();
				fMCXiY  [nTruth] = fMCEvt->GetTrack(idxXi[0])->Y();
			}
			nTruth++;

			if (nTruth > MaxNTruth) cout <<Form("WARNING: nTruth exceeds max in run%i evt%i\n", fEvtRunNo, fEvtID);
		}//a, tracks (MC)
		fMCNum = nTruth;

		if (ValidEvtOnly && nTruth==0) return; //!!
	}//MC

	//Trackwise selection starts
	//-----------------------------------------------------

	//Cascade (charged Xi)
	int nXi = 0;
	const int nCasc = fEvt->GetNumberOfCascades();
	for (int a=0; a<nCasc; a++)
	{
		fCasc = ((AliAODEvent*)fEvt)->GetCascade(a);
		if (!fCasc->GetSecondaryVtx() || !fCasc->GetDecayVertexXi()) continue;
		if (FilterCascade(fCasc, fEvtVtx, fPID) == false) continue;

		//++++++++++++++++++++++++++++++++++++++++++++++++++++

		//#
		fCascChgXi  [nXi] = fCasc->ChargeXi();
		fCascCosPAXi[nXi] = fCasc->CosPointingAngleXi(fEvtVtx->GetX(), fEvtVtx->GetY(), fEvtVtx->GetZ());
		fCascCosPAV0[nXi] = fCasc->CosPointingAngle(fCasc->GetDecayVertexXi());

		fCascDcaBachToPV[nXi] = fCasc->DcaBachToPrimVertex();
		fCascDcaV0ToPV  [nXi] = fCasc->DcaV0ToPrimVertex();
		fCascDcaXiDau   [nXi] = fCasc->DcaXiDaughters();
		fCascDcaV0Dau   [nXi] = fCasc->DcaV0Daughters();
		fCascDcaPosToPV [nXi] = fCasc->DcaPosToPrimVertex();
		fCascDcaNegToPV [nXi] = fCasc->DcaNegToPrimVertex();

		fCascDecayLenXi   [nXi] = fCasc->DecayLengthXi(fEvtVtx->GetX(), fEvtVtx->GetY(), fEvtVtx->GetZ());
		fCascDecayLenV0   [nXi] = fCasc->DecayLengthV0();
		fCascDecayLenXiOld[nXi] = sqrt( pow(fCasc->DecayVertexXiX(), 2) + pow(fCasc->DecayVertexXiY(), 2) );
		fCascDecayLenV0Old[nXi] = sqrt( pow(fCasc->DecayVertexV0X(), 2) + pow(fCasc->DecayVertexV0Y(), 2) );

		fCascMassLmb    [nXi] = fCasc->MassLambda();
		fCascMassLmbAnti[nXi] = fCasc->MassAntiLambda();
		fCascMassOmega  [nXi] = fCasc->MassOmega();
		fCascMassXi     [nXi] = fCasc->MassXi();

		fCascPtXi[nXi] = sqrt(fCasc->Pt2Xi());
		fCascPxXi[nXi] = fCasc->MomXiX();
		fCascPyXi[nXi] = fCasc->MomXiY();
		fCascPzXi[nXi] = fCasc->MomXiZ();

		//++++++++++++++++++++++++++++++++++++++++++++++++++++

		AliAODTrack* fBach_pi = (AliAODTrack*)fCasc->GetDecayVertexXi()->GetDaughter(0);
		AliAODTrack* fV0d_pos = (AliAODTrack*)fCasc->GetDaughter(0);
		AliAODTrack* fV0d_neg = (AliAODTrack*)fCasc->GetDaughter(1);

		fCascPt_BachPi      [nXi] = fBach_pi->Pt();
		fCascPt_V0dPos      [nXi] = fV0d_pos->Pt();
		fCascPt_V0dNeg      [nXi] = fV0d_neg->Pt();
		fCascTPCNxedR_BachPi[nXi] = fBach_pi->GetTPCNCrossedRows();
		fCascTPCNxedR_V0dPos[nXi] = fV0d_pos->GetTPCNCrossedRows();
		fCascTPCNxedR_V0dNeg[nXi] = fV0d_neg->GetTPCNCrossedRows();
		fCascTPCNclsF_BachPi[nXi] = fBach_pi->GetTPCNclsF();
		fCascTPCNclsF_V0dPos[nXi] = fV0d_pos->GetTPCNclsF();
		fCascTPCNclsF_V0dNeg[nXi] = fV0d_neg->GetTPCNclsF();
		if (IsMC)
		{
			fCascMomLabel[nXi] = GetCascMotherID(fMCEvt, fCasc); //For truth mathcing later
			if (fCascMomLabel[nXi] >= 0) fCascMomPDG[nXi] = fMCEvt->GetTrack(fCascMomLabel[nXi])->PdgCode();
		}

		//xCheck, Dec. 2022
		if ( (fCascDecayLenXiOld[nXi] > cutCasc_minDecayLenXi) && (fCascDecayLenV0Old[nXi] > cutCasc_minDecayLenV0) )
		{
			fHisto->FillTH1("c_massXi", fCasc->MassXi());

			if ( (fEvtTrig & AliVEvent::kINT7) && (fEvtMult >= 0.0 && fEvtMult <= 100.0) ) 
			{
				fHisto->FillTH1("c_massXi_MB", fCasc->MassXi());
			}
			if ( (fEvtTrig & AliVEvent::kHighMultV0) && (fEvtMult >= 0.0 && fEvtMult <= 0.1) )
			{
				fHisto->FillTH1("c_massXi_HMV0", fCasc->MassXi());
			}
		}

		nXi++;
		if (nXi > MaxNCasc) cout <<Form("WARNING: nCasc exceeds max in run%i evt%i\n", fEvtRunNo, fEvtID);
	}//a, nCasc
	fCascNum = nXi;
	if (ValidEvtOnly && nXi==0) return; //!!

	//Electron tracks
	int nEle = 0;
	vector<int> Blacklist;
	const int nTracks = fEvt->GetNumberOfTracks();
	for (int a=0; a<nTracks; a++)
	{
		fTrk = (AliAODTrack*)fEvt->GetTrack(a);

		if (FilterTrack(fTrk, fEvtVtx) == false) continue; //Track filter 1, common
		if (FilterTrackElectron(fTrk, fPID) == false) continue; //Track filter 2, electron

		//!! Reject electrons on the blacklist (i.e., pairs from photon conversion)
		if ( ValidEvtOnly &&
			 (std::find(Blacklist.begin(), Blacklist.end(), fTrk->GetID())!=Blacklist.end()) ) continue;

		//Check if this electron is originated from conversion
		//++++++++++++++++++++++++++++++++++++++++++++++++++++

		const Short_t  e1_charge = fTrk->Charge();
		const Int_t    e1_trkID  = fTrk->GetID();
		const Double_t e1_energy = fTrk->E(MassEle);
		const Double_t e1_px     = fTrk->Px();
		const Double_t e1_py     = fTrk->Py();
		const Double_t e1_pz     = fTrk->Pz();

		Double_t minMass_ls = 999.; //Like-Sign
		Double_t minMass_us = 999.; //Unlike-Sign
		for (int b=0; b<nTracks; b++)
		{
			AliAODTrack* tempTrk = (AliAODTrack*)fEvt->GetTrack(b);

			//Apply loose track cut
			if ( !tempTrk || (tempTrk->GetID() == e1_trkID) ||
				 (tempTrk->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA) == false) ||
				 (fabs(fPID->NumberOfSigmasTPC(tempTrk, AliPID::kElectron)) > cutEle_nSigmaTPCAbsConv) ) continue;

			const Short_t  e2_charge = tempTrk->Charge();
			const Double_t e2_energy = tempTrk->E(MassEle);
			const Double_t e2_px     = tempTrk->Px();
			const Double_t e2_py     = tempTrk->Py();
			const Double_t e2_pz     = tempTrk->Pz();
			const Double_t e1e2_mass = sqrt( pow(e1_energy+e2_energy, 2) -
					( pow(e1_px+e2_px, 2) + pow(e1_py+e2_py, 2) + pow(e1_pz+e2_pz, 2) ) );

			if ( (e1_charge*e2_charge > 0) && (e1e2_mass < minMass_ls) ) minMass_ls = e1e2_mass;
			if ( (e1_charge*e2_charge < 0) && (e1e2_mass < minMass_us) ) minMass_us = e1e2_mass;
		}//b, 2nd track loop
		//cout <<Form("%5i %3i %3i %6.4f \n", fEvtID, fTrk, fTrk->GetID(), fTrk->Pt());
		fHisto->FillTH1("e_minMassUS", minMass_us); //xCheck, Nov. 8, 2022

		//!! Track filter 3, electron from photon conversion
		if ( ValidEvtOnly && (minMass_us < cutEle_massConv) ) { Blacklist.push_back(fTrk->GetID()); continue; }

		//++++++++++++++++++++++++++++++++++++++++++++++++++++

		//#
		fEleMinMassLS[nEle] = minMass_ls;
		fEleMinMassUS[nEle] = minMass_us;
		fEleNSigmaTOF[nEle] = fPID->NumberOfSigmasTOF(fTrk, AliPID::kElectron);
		fEleNSigmaTPC[nEle] = fPID->NumberOfSigmasTPC(fTrk, AliPID::kElectron);

		fEleChg[nEle] = fTrk->Charge();
		fEleEta[nEle] = fTrk->Eta();
		fElePhi[nEle] = fTrk->Phi();
		fElePt [nEle] = fTrk->Pt();
		fElePx [nEle] = fTrk->Px();
		fElePy [nEle] = fTrk->Py();
		fElePz [nEle] = fTrk->Pz();
		fEleY  [nEle] = fTrk->Y();

		fEleITSNcls [nEle] = fTrk->GetITSNcls();
		fEleTPCNsig [nEle] = fTrk->GetTPCsignalN();
		fEleTPCNxedR[nEle] = fTrk->GetTPCNCrossedRows();
		fEleTPCNclsF[nEle] = fTrk->GetTPCNclsF();

		if (IsMC)
		{
			const int l_ele = fTrk->GetLabel();

			fEleLabel   [nEle] = l_ele;
			fEleMomLabel[nEle] = fMCEvt->GetTrack( abs(l_ele) )->GetMother();
			fEleMomPDG  [nEle] = fMCEvt->GetTrack( fEleMomLabel[nEle] )->PdgCode();
		}

		nEle++;
		if (nEle > MaxNEle) cout <<Form("WARNING: nEle exceeds max in run%i evt%i\n", fEvtRunNo, fEvtID);
	}//a, nTracks
	fEleNum = nEle;
	if (ValidEvtOnly && nEle==0) return; //!!

	if (nXi>0 && nEle>0) fTree->Fill(); //#
	ControlOutputContainers(1);
	return;
}//UserExec

//=======================================================================================

int AliAnalysisTaskSEXic0SL::CheckOrigin(AliMCEvent* MCEvt, AliAODMCParticle *MCPart)
{
	//Original code from AliVertexingHFUtils::CheckOrigin()
	bool isQuarkFound = false;
	bool isFromB = false;
	//bool isFromC = false;

	int step = 0;
	int momTrk = MCPart->GetMother();
	while (momTrk >= 0)
	{
		AliAODMCParticle* momPart = (AliAODMCParticle*)MCEvt->GetTrack(momTrk);
		if (!momPart) {	cout <<"ERROR: cannot find the mother particle object!\n"; break; }
		else
		{
			//Charm quark: 4
			//Charmed/ccbar mesons: 411 (D+) ~ 445 (chi_x2(1P))
			//Charmed baryons: 4122 (Lambda_c+) ~ 4444 (Omega_ccc++)
			//Bottom quark: 5
			//Bottom/bbar mesons: 511 (B0) ~ 557 (Upsilon_3(1D))
			//Bottom baryons: 5122 (Lambda_b0) ~ 5554 (Omega_bbb-)
			const Int_t momPdg = abs(momPart->GetPdgCode()); //Take absolute number

			//if ((momPdg==4) || (momPdg>400 && momPdg<500) || (momPdg>4000 && momPdg<5000)) isFromC = true;
			if ((momPdg==5) || (momPdg>500 && momPdg<600) || (momPdg>5000 && momPdg<6000)) isFromB = true;
			if (momPdg==4 || momPdg==5)	isQuarkFound = true;

			//cout <<Form("%i %i %4i / C%i B%i\n", fEvtID, step, momPdg, isFromC, isFromB);
			momTrk = momPart->GetMother(); //Climb up one ancetry level
			step++;
		}
	}

	//Last updated Jan. 3, 2023
	if (isQuarkFound == false) return -1;
	else if (isFromB == true) return 5; //Quark found and is bottom
	else return 4; //Quark found and it's NOT bottom, thus charm (only Xic0 or Xic+ supposed to be provided)
}//CheckOrigin

int AliAnalysisTaskSEXic0SL::GetCascMotherID(AliMCEvent* MCEvt, AliAODcascade* Casc)
{
	//CAVEAT: reco level, these indice can be negative if quality is poor, but SHOULD BE USED
	const Int_t l_bachPi = ((AliAODTrack*)Casc->GetDecayVertexXi()->GetDaughter(0))->GetLabel();
	const Int_t l_v0dPos = ((AliAODTrack*)Casc->GetDaughter(0))->GetLabel();
	const Int_t l_v0dNeg = ((AliAODTrack*)Casc->GetDaughter(1))->GetLabel();

	//CAVEAT: truth level hereafter, provide argument wrapped with absolute (otherwise segfault happens, anyway)
	const Int_t m_bachPi = MCEvt->GetTrack( abs(l_bachPi) )->GetMother(); //= Xi (strange, charged)
	const Int_t m_v0dPos = MCEvt->GetTrack( abs(l_v0dPos) )->GetMother(); //= lambda0
	const Int_t m_v0dNeg = MCEvt->GetTrack( abs(l_v0dNeg) )->GetMother(); //= lambda0
	if (m_v0dPos != m_v0dNeg) return -999;

	const Int_t n_v0dPos = MCEvt->GetTrack( m_v0dPos )->GetMother(); //= Xi
	const Int_t n_v0dNeg = MCEvt->GetTrack( m_v0dNeg )->GetMother(); //= Xi
	if (n_v0dPos != n_v0dNeg) return -998;

	if ( (m_bachPi == n_v0dPos) && (m_bachPi == n_v0dNeg) && (n_v0dPos == n_v0dNeg)	)
	{
		const Int_t l_Xic0 = MCEvt->GetTrack( m_bachPi )->GetMother();
		return l_Xic0;
	}
	else return -997;
}//GetCascMotherID

bool AliAnalysisTaskSEXic0SL::FilterTrack(AliAODTrack* Trk, const AliVVertex* Vtx)
{
	if ( (Trk->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA) == false) || //Filterbit 4
		 (Trk->GetTPCsignalN() < cut_TPCsignalN) || //fSetProdTrackTPCNclsPID
		 (fabs(Trk->Eta()) > cut_trkEta) ||
		 (Trk->Pt() < cut_trkPt) ) return false; //Cut 1,2,3,4

	//Track quality
	AliAODTrack* tempTrkAOD = dynamic_cast<AliAODTrack*>(Trk);
	AliESDtrack tempTrkESD(tempTrkAOD);
	tempTrkESD.SetTPCClusterMap(tempTrkAOD->GetTPCClusterMap());
	tempTrkESD.SetTPCSharedMap(tempTrkAOD->GetTPCSharedMap());
	tempTrkESD.SetTPCPointsF(tempTrkAOD->GetTPCNclsF());
	Double_t pos[3]; Vtx->GetXYZ(pos);
	Double_t cvm[6]; Vtx->GetCovarianceMatrix(cvm);

	//AliESDVertex(const Double_t position[3], const Double_t covmatrix[6], Double_t chi2, Int_t nContributors...)
	//RelateToVertex (const AliESDVertex *vtx, Double_t b, Double_t maxd...)
	const AliESDVertex vtxESD(pos, cvm, 100., 100);
	tempTrkESD.RelateToVertex(&vtxESD, 0., 3.);
	if (fTrkCuts->IsSelected(&tempTrkESD) == false) return false; //Cut 5

	return true;
}//FilterTrack

bool AliAnalysisTaskSEXic0SL::FilterTrackElectron(AliAODTrack* Trk, AliPIDResponse* PID)
{
	//Regard all tracks being tested are electron candidates
	const Float_t e_nSigmaTOF = PID->NumberOfSigmasTOF(Trk, AliPID::kElectron); //Default -999
	const Float_t e_nSigmaTPC = PID->NumberOfSigmasTPC(Trk, AliPID::kElectron); //Default ?
	const Double_t e_pT = Trk->Pt();
	const Double_t e_pT_cutVal = (e_pT > 5.0) ? 5.0:e_pT;
	const Double_t cutEle_nSigmaTPCMin = -4.3 + (1.17 * e_pT_cutVal) - 0.094 * pow(e_pT_cutVal, 2);

	if ( (fabs(e_nSigmaTOF) > cutEle_nSigmaTOFAbs) ||
		 (e_nSigmaTPC > cutEle_nSigmaTPCMax) ||
		 (e_nSigmaTPC < cutEle_nSigmaTPCMin) ) return false; //Cut 1, 2

	return true;
}//FilterTrackElectron

bool AliAnalysisTaskSEXic0SL::FilterCascade(AliAODcascade* Casc, const AliVVertex* Vtx, AliPIDResponse* PID)
{
	//Xic0 (m: 2470.91 +- 0.25 MeV):
	//Xic0 -> "e+ Xi-" -> "e+ (pi- lambda0 nu)"     -> "e+ (pi- (p     pi-) nu)" or -> XiNeg
	//Xic0 -> "e- Xi+" -> "e- (pi+ lambda0_bar nu)" -> "e- (pi+ (p_bar pi+) nu)"    -> XiPos
	//i.e., at final level, cascade has "bach_pi+(-) & v0_p & v0_pi-"
	const bool XiNeg = (Casc->ChargeXi()<0)?true:false;

	AliAODTrack* bach_pi = (AliAODTrack*)Casc->GetDecayVertexXi()->GetDaughter(0); //AliAODcascade::GetBachID()
	AliAODTrack* v0d_pos = (AliAODTrack*)Casc->GetDaughter(0); //Positive (checked: data and mc)
	AliAODTrack* v0d_neg = (AliAODTrack*)Casc->GetDaughter(1); //Negative
	if (!bach_pi || !v0d_pos || !v0d_neg) return false; //Return 1

	//nSigmaTPC
	if (XiNeg) //Xic0 => e+ Xi- => e+ (pi- (p pi-))
	{
		if ( (fabs(PID->NumberOfSigmasTPC(bach_pi, AliPID::kPion))   > cutCasc_nSigmaTPCAbs) ||
			 (fabs(PID->NumberOfSigmasTPC(v0d_pos, AliPID::kProton)) > cutCasc_nSigmaTPCAbs) ||
			 (fabs(PID->NumberOfSigmasTPC(v0d_neg, AliPID::kPion))   > cutCasc_nSigmaTPCAbs) ) return false; //2a
	}
	else //XiPos: Xic0 => e- Xi+ => e- (pi+ (p_bar pi+)
	{
		if ( (fabs(PID->NumberOfSigmasTPC(bach_pi, AliPID::kPion))   > cutCasc_nSigmaTPCAbs) ||
			 (fabs(PID->NumberOfSigmasTPC(v0d_pos, AliPID::kPion))   > cutCasc_nSigmaTPCAbs) ||
			 (fabs(PID->NumberOfSigmasTPC(v0d_neg, AliPID::kProton)) > cutCasc_nSigmaTPCAbs) ) return false; //2b
	}

	//!! Old decay length (actually transverse radius) cut: DO NOT USE THIS FOR NORMAL ANALYSIS
	//const Double_t DecayLenXi = sqrt( pow(fCasc->DecayVertexXiX(), 2) + pow(fCasc->DecayVertexXiY(), 2) );
	//const Double_t DecayLenV0 = sqrt( pow(fCasc->DecayVertexV0X(), 2) + pow(fCasc->DecayVertexV0Y(), 2) );
	//if ( (DecayLenXi < cutCasc_minDecayLenXi) || (DecayLenV0 < cutCasc_minDecayLenV0) ) return false;

	//!! Decay length cut: turn off for now, for xCheck (Oct. 5, 2022)
	//const Double_t DecayLenXi = Casc->DecayLengthXi(Vtx->GetX(), Vtx->GetY(), Vtx->GetZ()); //PV to Cascade (Xi)
	//const Double_t DecayLenV0 = Casc->DecayLengthV0(); //Cascade to V0
	//if ( (DecayLenXi < cutCasc_minDecayLenXi) || (DecayLenV0 < cutCasc_minDecayLenV0) ) return false; //3

	//DCA of Bach/V0 to PV
	const Double_t DcaBachToPV = Casc->DcaBachToPrimVertex();
	const Double_t DcaV0ToPV   = Casc->DcaV0ToPrimVertex();
	if ( (DcaBachToPV < cutCasc_minDcaBachToPV) || (DcaV0ToPV < cutCasc_minDcaV0ToPV) ) return false; //4

	//DCA of V0/Xi to its daughters
	const Double_t DcaXiDau = Casc->DcaXiDaughters();
	const Double_t DcaV0Dau = Casc->DcaV0Daughters();
	if ( (DcaXiDau > cutCasc_maxDcaXiDau) || (DcaV0Dau > cutCasc_maxDcaV0Dau) ) return false; //5

	//DCA of V0 daughters to PV (regardless of Xi sign since cut values are the same)
	//* AliAODv0::DcaPosToPrimVertex gets "fd0[0] from AliAODRecoDecay.h: rphi impact params WRT Primary Vtx [cm]"
	const Double_t DcaV0dPosToPV = Casc->DcaPosToPrimVertex();
	const Double_t DcaV0dNegToPV = Casc->DcaNegToPrimVertex();
	if ( (DcaV0dPosToPV < cutCasc_minDcaV0PosToPV) || (DcaV0dNegToPV < cutCasc_minDcaV0NegToPV) ) return false; //6

	//Pointing angles
	//* AliAODRecoDecay::CosPointingAngle(AliAODVertex *vtx1)
	const Double_t CosPAngleXi = Casc->CosPointingAngleXi(Vtx->GetX(), Vtx->GetY(), Vtx->GetZ()); //Xi to PV
	const Double_t CosPAngleV0 = Casc->CosPointingAngle(Casc->GetDecayVertexXi()); //V0 to Xi
	if ( (CosPAngleXi < cutCasc_minCosPAngleXi) || (CosPAngleV0 < cutCasc_minCosPAngleV0) ) return false; //7

	//Close to Lambda0 mass: 1115.683 +- 0.006 (MeV), cut tolerance: 0.008
	const Double_t massLmbPDG  = TDatabasePDG::Instance()->GetParticle(PDGCode_Lambda)->Mass();
	const Double_t massLmbTmp0 = Casc->MassLambda();
	const Double_t massLmbTmp1 = Casc->MassAntiLambda();
	if ( (fabs(massLmbPDG - massLmbTmp0) > cutCasc_massTolLambda) &&
		 (fabs(massLmbPDG - massLmbTmp1) > cutCasc_massTolLambda) ) return false; //8

	//Close to Xi mass: 1321.71 +- 0.07 (MeV) -> 1.32171 +- 0.00007 (GeV), cut tolerance: 0.01
	const Double_t massXiPDG = TDatabasePDG::Instance()->GetParticle(PDGCode_Xi)->Mass();
	const Double_t massXiTmp = Casc->MassXi();
	if ( fabs(massXiPDG - massXiTmp) > cutCasc_massTolXi ) return false; //9

	/*
	//NOT SURE ABOUT FOLLOWING CUTS...
	//Check validity of mother (Xi)
	if (sqrt(Casc->Pt2Xi()) > 999.) return false; //Original code calculates this manually

	//Rule out Omega by mass: Omega- (1672.45 +- 0.29 MeV) -> Lambda_0 + K- (67.8 %)
	const Double_t massOmegaPDG = TDatabasePDG::Instance()->GetParticle(PDGCode_Omega)->Mass();
	const Double_t massOmegaTmp = Casc->MassOmega();
	if ( fabs(massOmegaPDG - massOmegaTmp) < cutCasc_massTolOmega ) return false; //Original cut was "0.0"
	*/

	return true;
}//FilterCascade

//=======================================================================================

void AliAnalysisTaskSEXic0SL::ControlOutputContainers(int option)
{
	if (option == 0)
	{
		DefineOutput(0, TChain::Class());
		DefineOutput(1, TDirectory::Class());
		DefineOutput(2, TTree::Class());
		DefineOutput(3, AliNormalizationCounter::Class());
		DefineOutput(4, AliNormalizationCounter::Class());
		DefineOutput(5, AliNormalizationCounter::Class());
		DefineOutput(6, AliNormalizationCounter::Class());
		DefineOutput(7, AliNormalizationCounter::Class());
		DefineOutput(8, AliNormalizationCounter::Class());
		DefineOutput(9, AliNormalizationCounter::Class());
		DefineOutput(10, AliNormalizationCounter::Class());
	}
	else if (option == 1)
	{
		PostData(1, fHisto->GetListOfHistograms());
		PostData(2, fTree);
		PostData(3, fANC_MB_0to100);
		PostData(4, fANC_MB_30to100);
		PostData(5, fANC_MB_0p1to30);
		PostData(6, fANC_HMV0_0to0p1);
		PostData(7, fANC_INEL0_MB_0to100);
		PostData(8, fANC_INEL0_MB_30to100);
		PostData(9, fANC_INEL0_MB_0p1to30);
		PostData(10, fANC_INEL0_HMV0_0to0p1);
	}
	else AliFatal("ERROR: improper use of ControlOutputContainers!");
	return;
}//ControlOutputContainers

void AliAnalysisTaskSEXic0SL::ControlAnaObjects(int option)
{
	if (option == 0)
	{
		fCasc = 0;
		fMCPart = 0;
		fTrk = 0;
		fTrkCuts = 0;
		fInputHandler = 0;
		fMCEvt = 0;
		fMultSel = 0;
		fANC_MB_0to100 = 0;
		fANC_MB_30to100 = 0;
		fANC_MB_0p1to30 = 0;
		fANC_HMV0_0to0p1 = 0;
		fANC_INEL0_MB_0to100 = 0;
		fANC_INEL0_MB_30to100 = 0;
		fANC_INEL0_MB_0p1to30 = 0;
		fANC_INEL0_HMV0_0to0p1 = 0;
		fPID = 0;
		fEvtCutsMB = 0;
		fEvtCutsHMV0 = 0;
		fEvt = 0;
		fEvtVtx = 0;
		fHisto = 0;
		fTree = 0;
	}
	else if (option == 1)
	{
		delete fCasc;
		delete fMCPart;
		delete fTrk;
		delete fTrkCuts;
		delete fInputHandler;
		delete fMCEvt;
		delete fMultSel;
		delete fANC_MB_0to100;
		delete fANC_MB_30to100;
		delete fANC_MB_0p1to30;
		delete fANC_HMV0_0to0p1;
		delete fANC_INEL0_MB_0to100;
		delete fANC_INEL0_MB_30to100;
		delete fANC_INEL0_MB_0p1to30;
		delete fANC_INEL0_HMV0_0to0p1;
		delete fPID;
		delete fEvtCutsMB;
		delete fEvtCutsHMV0;
		delete fEvt;
		delete fEvtVtx;
		delete fHisto;
		delete fTree;
	}
	else AliFatal("ERROR: improper use of ControlAnaObjects!");

	return;
}//ControlAnaObjects

void AliAnalysisTaskSEXic0SL::ControlOutputTree(TTree* T, bool isMC, bool readOnly)
{
	int iLeaf = 0;

	if (isMC == false) //data
	{
		iLeaf = 0;
		if (!readOnly)
		{
			TString strEvt = "evtID/i:trig/i:runNo/I:mult/F:vtxZ/D:goodMB/O:goodHMV0/O:INEL0/O";
			T->Branch("Event", 0, strEvt);
		}
		((TLeaf*)T->GetBranch("Event")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEvtID);       iLeaf++; //0
		((TLeaf*)T->GetBranch("Event")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEvtTrig);     iLeaf++;
		((TLeaf*)T->GetBranch("Event")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEvtRunNo);    iLeaf++;
		((TLeaf*)T->GetBranch("Event")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEvtMult);     iLeaf++;
		((TLeaf*)T->GetBranch("Event")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEvtVtxZ);     iLeaf++;
		((TLeaf*)T->GetBranch("Event")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEvtGoodMB);   iLeaf++; //5
		((TLeaf*)T->GetBranch("Event")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEvtGoodHMV0); iLeaf++;
		((TLeaf*)T->GetBranch("Event")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEvtINELLgt0); iLeaf++;
	}
	else //MC
	{
		iLeaf = 0;
		if (!readOnly) T->Branch("Event", 0, "runNo/I:vtxZ/D");
		((TLeaf*)T->GetBranch("Event")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEvtRunNo); iLeaf++;
		((TLeaf*)T->GetBranch("Event")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEvtVtxZ);  iLeaf++;

		iLeaf = 0;
		if (!readOnly)
		{
			TString strMCTruth = "mcN/I:mc_label[mcN]/I:mc_orig[mcN]/I:mc_PDG[mcN]/I:mc_Pt[mcN]/D:mc_Y[mcN]/D";
			strMCTruth += ":mc_ElePt[mcN]/D:mc_EleY[mcN]/D";
			strMCTruth += ":mc_XiPt[mcN]/D:mc_XiY[mcN]/D";
			T->Branch("MCTruth", 0, strMCTruth.Data());
		}
		((TLeaf*)T->GetBranch("MCTruth")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fMCNum);   iLeaf++; //0
		((TLeaf*)T->GetBranch("MCTruth")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fMCLabel); iLeaf++;
		((TLeaf*)T->GetBranch("MCTruth")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fMCOrig);  iLeaf++;
		((TLeaf*)T->GetBranch("MCTruth")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fMCPDG);   iLeaf++;
		((TLeaf*)T->GetBranch("MCTruth")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fMCPt);    iLeaf++;
		((TLeaf*)T->GetBranch("MCTruth")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fMCY);     iLeaf++;
		((TLeaf*)T->GetBranch("MCTruth")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fMCElePt); iLeaf++;
		((TLeaf*)T->GetBranch("MCTruth")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fMCEleY);  iLeaf++;
		((TLeaf*)T->GetBranch("MCTruth")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fMCXiPt);  iLeaf++;
		((TLeaf*)T->GetBranch("MCTruth")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fMCXiY);   iLeaf++;
	}

	//+++++++++++++++++++++++++++++++++++++++++++

	iLeaf = 0;
	if (!readOnly)
	{
		TString strEle = "eN/I:e_charge[eN]/I:e_itsNcls[eN]/I";
		strEle += ":e_minMassLS[eN]/F:e_minMassUS[eN]/F:e_nSigmaTOF[eN]/F:e_nSigmaTPC[eN]/F";
		strEle += ":e_eta[eN]/D:e_phi[eN]/D:e_pT[eN]/D:e_px[eN]/D:e_py[eN]/D:e_pz[eN]/D:e_Y[eN]/D";
		strEle += ":e_tpcNsig[eN]/s:e_tpcNxedR[eN]/s:e_tpcNclsF[eN]/s";
		//
		if (IsMC) strEle += ":e_label[eN]/I:e_momLabel[eN]/I:e_momPDG[eN]/I";
		T->Branch("Ele", 0, strEle.Data());
	}
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleNum);       iLeaf++; //0
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleChg);       iLeaf++;
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleITSNcls);   iLeaf++;
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleMinMassLS); iLeaf++;
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleMinMassUS); iLeaf++;
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleNSigmaTOF); iLeaf++; //5
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleNSigmaTPC); iLeaf++;
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleEta);       iLeaf++;
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fElePhi);       iLeaf++;
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fElePt);        iLeaf++;
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fElePx);        iLeaf++; //10
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fElePy);        iLeaf++;
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fElePz);        iLeaf++;
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleY);         iLeaf++;
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleTPCNsig);   iLeaf++;
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleTPCNxedR);  iLeaf++; //15
	((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleTPCNclsF);  iLeaf++;
	//
	if (IsMC)
	{
		((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleLabel);    iLeaf++;
		((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleMomLabel); iLeaf++;
		((TLeaf*)T->GetBranch("Ele")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fEleMomPDG);   iLeaf++;
	}

	//+++++++++++++++++++++++++++++++++++++++++++

	iLeaf = 0;
	if (!readOnly)
	{
		TString strCasc = "cN/I:c_chgXi[cN]/I";
		strCasc += ":c_cosPAXi[cN]/D:c_cosPAV0[cN]/D:c_dcaBachToPV[cN]/D:c_dcaV0ToPV[cN]/D";
		strCasc += ":c_dcaXiDau[cN]/D:c_dcaV0Dau[cN]/D:c_dcaPosToPV[cN]/D:c_dcaNegToPV[cN]/D";
		//
		strCasc += ":c_dLenXi[cN]/D:c_dLenXiOld[cN]/D:c_dLenV0[cN]/D:c_dLenV0Old[cN]/D";
		strCasc += ":c_massLmb[cN]/D:c_massLmbAnti[cN]/D:c_massOmega[cN]/D:c_massXi[cN]/D";
		strCasc += ":c_pTXi[cN]/D:c_pxXi[cN]/D:c_pyXi[cN]/D:c_pzXi[cN]/D";
		//
		strCasc += ":c_pT_BachPi[cN]/D:c_pT_V0dPos[cN]/D:c_pT_V0dNeg[cN]/D";
		strCasc += ":c_tpcNxedR_BachPi[cN]/s:c_tpcNxedR_V0dPos[cN]/s:c_tpcNxedR_V0dNeg[cN]/s";
		strCasc += ":c_tpcNclsF_BachPi[cN]/s:c_tpcNclsF_V0dPos[cN]/s:c_tpcNclsF_V0dNeg[cN]/s";
		//
		if (IsMC) strCasc += ":c_momLabel[cN]/I:c_momPDG[cN]/I";
		T->Branch("Casc", 0, strCasc.Data());
	}
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascNum);           iLeaf++;
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascChgXi);         iLeaf++;
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascCosPAXi);       iLeaf++; //CosPA
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascCosPAV0);       iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascDcaBachToPV);   iLeaf++; //DCA
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascDcaV0ToPV);     iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascDcaXiDau);      iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascDcaV0Dau);      iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascDcaPosToPV);    iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascDcaNegToPV);    iLeaf++; //|
	//
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascDecayLenXi);    iLeaf++; //DLen
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascDecayLenXiOld); iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascDecayLenV0);    iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascDecayLenV0Old); iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascMassLmb);       iLeaf++; //Mass
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascMassLmbAnti);   iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascMassOmega);     iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascMassXi);        iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascPtXi);          iLeaf++; //Xi_pT
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascPxXi);          iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascPyXi);          iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascPzXi);          iLeaf++; //|
	//
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascPt_BachPi);       iLeaf++; //pT
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascPt_V0dPos);       iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascPt_V0dNeg);       iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascTPCNxedR_BachPi); iLeaf++; //TPCx
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascTPCNxedR_V0dPos); iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascTPCNxedR_V0dNeg); iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascTPCNclsF_BachPi); iLeaf++; //TPCf
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascTPCNclsF_V0dPos); iLeaf++; //|
	((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascTPCNclsF_V0dNeg); iLeaf++; //|
	//
	if (IsMC)
	{
		((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascMomLabel); iLeaf++;
		((TLeaf*)T->GetBranch("Casc")->GetListOfLeaves()->At(iLeaf))->SetAddress(&fCascMomPDG);   iLeaf++;
	}

	return;
}//ControlOutputTree

void AliAnalysisTaskSEXic0SL::ResetTreeVariables(void)
{
	fEvtTrig     = -999;
	fEvtRunNo    = -999;
	fEvtMult     = -999.;
	fEvtVtxZ     = -999.;
	fEvtGoodMB   = false;
	fEvtGoodHMV0 = false;
	fEvtINELLgt0 = false;

	fMCNum = -999;
	for (int a=0; a<MaxNTruth; a++)
	{
		fMCLabel[a] = -999;
		fMCOrig [a] = -999;
		fMCPDG  [a] = -999;
		fMCPt   [a] = -999.;
		fMCY    [a] = -999.;
		fMCElePt[a] = -999.;
		fMCEleY [a] = -999.;
		fMCXiPt [a] = -999.;
		fMCXiY  [a] = -999.;
	}//fMC

	fEleNum = -999;
	for (int a=0; a<MaxNEle; a++)
	{
		fEleChg      [a] = -999;
		fEleITSNcls  [a] = -999;
		fEleLabel    [a] = -999;
		fEleMomLabel [a] = -999;
		fEleMomPDG   [a] = -999;
		fEleMinMassLS[a] = -999.;
		fEleMinMassUS[a] = -999.;
		fEleNSigmaTOF[a] = -999.;
		fEleNSigmaTPC[a] = -999.;
		fEleEta      [a] = -999.;
		fElePhi      [a] = -999.;
		fElePt       [a] = -999.;
		fElePx       [a] = -999.;
		fElePy       [a] = -999.;
		fElePz       [a] = -999.;
		fEleY        [a] = -999.;
		fEleTPCNsig  [a] = 999; ///Unsigned short
		fEleTPCNxedR [a] = 999;
		fEleTPCNclsF [a] = 999;
	}//fEle

	fCascNum = -999;
	for (int a=0; a<MaxNCasc; a++)
	{
		fCascChgXi        [a] = -999;
		fCascMomLabel     [a] = -999;
		fCascMomPDG       [a] = -999;
		fCascCosPAXi      [a] = -999.;
		fCascCosPAV0      [a] = -999.;
		fCascDcaBachToPV  [a] = -999.;
		fCascDcaV0ToPV    [a] = -999.;
		fCascDcaXiDau     [a] = -999.;
		fCascDcaV0Dau     [a] = -999.;
		fCascDcaPosToPV   [a] = -999.;
		fCascDcaNegToPV   [a] = -999.;
		//
		fCascDecayLenXi   [a] = -999.;
		fCascDecayLenXiOld[a] = -999.;
		fCascDecayLenV0   [a] = -999.;
		fCascDecayLenV0Old[a] = -999.;
		fCascMassLmb      [a] = -999.;
		fCascMassLmbAnti  [a] = -999.;
		fCascMassOmega    [a] = -999.;
		fCascMassXi       [a] = -999.;
		fCascPtXi         [a] = -999.;
		fCascPxXi         [a] = -999.;
		fCascPyXi         [a] = -999.;
		fCascPzXi         [a] = -999.;
		//
		fCascPt_BachPi      [a] = -999.;
		fCascPt_V0dPos      [a] = -999.;
		fCascPt_V0dNeg      [a] = -999.;
		fCascTPCNxedR_BachPi[a] = 999;
		fCascTPCNxedR_V0dPos[a] = 999;
		fCascTPCNxedR_V0dNeg[a] = 999;
		fCascTPCNclsF_BachPi[a] = 999;
		fCascTPCNclsF_V0dPos[a] = 999;
		fCascTPCNclsF_V0dNeg[a] = 999;
	}//fCasc

	return;
}//ResetDummyIndices
