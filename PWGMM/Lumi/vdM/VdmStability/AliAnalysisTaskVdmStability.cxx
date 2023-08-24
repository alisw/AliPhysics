/*************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
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

// authors: I.Lofnes, ingrid.mckibben.lofnes@cern.ch
//          H.Degenhardt, hermann.franz.degenhardt@cern.ch
// Analysis task for determining the vdm-scan stability
//

#include <iostream>
using namespace std;

#include <TChain.h>
#include <TH1D.h>

#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliESDVZERO.h>

#include "AliAnalysisTaskVdmStability.h"
ClassImp(AliAnalysisTaskVdmStability)

//_________________________________________________________________________________
AliAnalysisTaskVdmStability::AliAnalysisTaskVdmStability():
AliAnalysisTaskSE(),
fOutputList(),
fRunNumber(0),
fVtxZ(-999.),
fnVtxCont(0),
fIsGoodZ(kFALSE),
fSelectPhysics(kFALSE),
fIsV0ANDfired(kFALSE),
fIsT0fired(kFALSE),
fIsEMCALfired(kFALSE),
fIsMUfired(kFALSE),
fIsDIMUfired(kFALSE),
fPileupEvent(kFALSE),
fTV0A(-10240.0f),
fTV0C(-10240.0f),
fGoodTime(kFALSE),
fzCut30(kFALSE),
fzCut10(kFALSE),
fzCut30nCont0(kFALSE),
fzCut10nCont0(kFALSE),
fzCut30nCont1(kFALSE),
fzCut10nCont1(kFALSE),
fV0TimeDiff_min(5.5),
fV0TimeDiff_max(11.5),
fV0TimeSum_min(11.5),
fV0TimeSum_max(17.5),
fEventTag(0),
fEvent(0x0),
fEventStatV0(0x0),
fEventStatT0(0x0),
fEventTree(0x0),
fNRuns(1000),
fNSelectionCases(25),
fFillTTree(false)
{
    // ROOT IO constructor, don't allocate memory here!
}

//_________________________________________________________________________________
AliAnalysisTaskVdmStability::AliAnalysisTaskVdmStability(const char* taskname):
AliAnalysisTaskSE(taskname),
fOutputList(),
fRunNumber(0),
fVtxZ(-999.),
fnVtxCont(0),
fIsGoodZ(kFALSE),
fSelectPhysics(kFALSE),
fIsV0ANDfired(kFALSE),
fIsT0fired(kFALSE),
fIsEMCALfired(kFALSE),
fIsMUfired(kFALSE),
fIsDIMUfired(kFALSE),
fPileupEvent(kFALSE),
fTV0A(-10240.0f),
fTV0C(-10240.0f),
fGoodTime(kFALSE),
fzCut30(kFALSE),
fzCut10(kFALSE),
fzCut30nCont0(kFALSE),
fzCut10nCont0(kFALSE),
fzCut30nCont1(kFALSE),
fzCut10nCont1(kFALSE),
fV0TimeDiff_min(5.5),
fV0TimeDiff_max(11.5),
fV0TimeSum_min(11.5),
fV0TimeSum_max(17.5),
fEventTag(0),
fEvent(0x0),
fEventStatV0(0x0),
fEventStatT0(0x0),
fEventTree(0x0),
fNRuns(1000),
fNSelectionCases(25),
fFillTTree(false)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());    //TList of event statistics
    DefineOutput(2, TTree::Class());    //event tree information
}

//_________________________________________________________________________________
AliAnalysisTaskVdmStability::~AliAnalysisTaskVdmStability()
{
    if(fEventStatV0)       { delete fEventStatV0;       fEventStatV0=0; }
    if(fEventStatT0)       { delete fEventStatT0;       fEventStatT0=0; }
    if (fEventTree)         { delete fEventTree;        fEventTree=0; }
    for (Int_t i = 0; i < 25; ++i) {
      if (v0_H[i])         		{ delete v0_H[i]; }
      if (t0_H[i])         		{ delete t0_H[i]; }
      if (emcal_H[i])         	{ delete emcal_H[i]; }
      if (muon_H[i])         	{ delete muon_H[i]; }
      if (dimuon_H[i])         	{ delete dimuon_H[i]; }
      if (v0_Timing[i])         { delete v0_Timing[i]; }
      if (t0_Timing[i])         { delete t0_Timing[i]; }
      if (emcal_Timing[i])      { delete emcal_Timing[i]; }
      if (muon_Timing[i])       { delete muon_Timing[i]; }
      if (dimuon_Timing[i])     { delete dimuon_Timing[i]; }
    }
}
//defíne output
//_________________________________________________________________________________
void AliAnalysisTaskVdmStability::UserCreateOutputObjects()
{	
    // set ownership of list
    fOutputList.SetOwner();
    // ---| skip list initialisation if already done |----------------------------
    if (!fOutputList.IsEmpty()) return;
    
    //Initialize output tree
    //fEventTree = std::unique_ptr<TTree>(new TTree("event", "event"));
    fEventTree = new TTree("event", "event");
    AddEventTreeVariables(fEventTree);
    
    // create our histo and add it to the list
    Int_t nbins = kNbinsEvent;
    if (!fEventStatV0){
        fEventStatV0=new TH1D("hEventStatV0","Event statistics V0",nbins,0,nbins);
        fEventStatV0->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
        fEventStatV0->GetXaxis()->SetBinLabel(2,"After Phys. Sel.");
        fEventStatV0->GetXaxis()->SetBinLabel(3,"V0 timing");
        fEventStatV0->GetXaxis()->SetBinLabel(4,"PU");
        fEventStatV0->GetXaxis()->SetBinLabel(5,"good z");
        fEventStatV0->GetXaxis()->SetBinLabel(6,"V0 timing + PU");
        fEventStatV0->GetXaxis()->SetBinLabel(7,"V0 timing + good z");
        fEventStatV0->GetXaxis()->SetBinLabel(8,"V0 timing + PU + good z");
        
    }
    fOutputList.Add(fEventStatV0);
    
    if (!fEventStatT0){
        fEventStatT0=new TH1D("hEventStatT0","Event statistics T0",nbins,0,nbins);
        fEventStatT0->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
        fEventStatT0->GetXaxis()->SetBinLabel(2,"After Phys. Sel.");
        fEventStatT0->GetXaxis()->SetBinLabel(3,"V0 timing");
        fEventStatT0->GetXaxis()->SetBinLabel(4,"PU");
        fEventStatT0->GetXaxis()->SetBinLabel(5,"good z");
        fEventStatT0->GetXaxis()->SetBinLabel(6,"V0 timing + PU");
        fEventStatT0->GetXaxis()->SetBinLabel(7,"V0 timing + good z");
        fEventStatT0->GetXaxis()->SetBinLabel(8,"V0 timing + PU + good z");
        
    }
    fOutputList.Add(fEventStatT0);
    
	TString selectionCases[25] = {
		"no_selection",
		"physics_selected",
		"V0_timing_cut",
		"pileup_rejection",
		"z-vertex_cut",
		"z-vertex_30",
		"z-vertex_10",
		"z-vertex_30_nContCut1",
		"z-vertex_10_nContCut1",
		"z-vertex_30_nContCut2",
		"z-vertex_10_nContCut2",
		"V0_timing_cut_pileup_rejection",
		"V0_timing_cut_z-vertex_cut",
		"V0_timing_cut_pileup_rejection_z-vertex_cut",
		"V0_timing_cut_pileup_rejection_z-vertex_30",
		"V0_timing_cut_pileup_rejection_z-vertex_10",
		"V0_timing_cut_pileup_rejection_z-vertex_30nContCut1",
		"V0_timing_cut_pileup_rejection_z-vertex_10nContCut1",
		"V0_timing_cut_pileup_rejection_z-vertex_30nContCut2",
		"V0_timing_cut_pileup_rejection_z-vertex_10nContCut2",
		"V0_timing_cut_z-vertex_30nContCut1",
		"TIME",
		"TIME_PS",
		"TIME_PU",
		"TIME_PS_PU"
	};
	
    for (Int_t iCase = 0; iCase < fNSelectionCases; iCase++){
        v0_H[iCase] = new TH1D(Form("v0_H_%s",selectionCases[iCase].Data()),Form("V0 events, %s",selectionCases[iCase].Data()),fNRuns,-0.5,fNRuns-0.5);
        t0_H[iCase] = new TH1D(Form("t0_H_%s",selectionCases[iCase].Data()),Form("T0 events, %s",selectionCases[iCase].Data()),fNRuns,-0.5,fNRuns-0.5);
        emcal_H[iCase] = new TH1D(Form("emcal_H_%s",selectionCases[iCase].Data()),Form("EMCAL events, %s",selectionCases[iCase].Data()),fNRuns,-0.5,fNRuns-0.5);
        muon_H[iCase] = new TH1D(Form("muon_H_%s",selectionCases[iCase].Data()),Form("MUON events, %s",selectionCases[iCase].Data()),fNRuns,-0.5,fNRuns-0.5);
        dimuon_H[iCase] = new TH1D(Form("dimuon_H_%s",selectionCases[iCase].Data()),Form("DIMUON events, %s",selectionCases[iCase].Data()),fNRuns,-0.5,fNRuns-0.5);
        v0_Timing[iCase] = new TH2D(Form("v0_Timing_%s",selectionCases[iCase].Data()),"sum vs diff",500,-30,40,500,-20,50);
        t0_Timing[iCase] = new TH2D(Form("t0_Timing_%s",selectionCases[iCase].Data()),"sum vs diff",500,-30,40,500,-20,50);
        emcal_Timing[iCase] = new TH2D(Form("emcal_Timing_%s",selectionCases[iCase].Data()),"sum vs diff",500,-30,40,500,-20,50);
        muon_Timing[iCase] = new TH2D(Form("muon_Timing_%s",selectionCases[iCase].Data()),"sum vs diff",500,-30,40,500,-20,50);
        dimuon_Timing[iCase] = new TH2D(Form("dimuon_Timing_%s",selectionCases[iCase].Data()),"sum vs diff",500,-30,40,500,-20,50);
		fOutputList.Add(v0_H[iCase]);
		fOutputList.Add(t0_H[iCase]);
		fOutputList.Add(emcal_H[iCase]);
		fOutputList.Add(muon_H[iCase]);
		fOutputList.Add(dimuon_H[iCase]);
		fOutputList.Add(v0_Timing[iCase]);
		fOutputList.Add(t0_Timing[iCase]);
		fOutputList.Add(emcal_Timing[iCase]);
		fOutputList.Add(muon_Timing[iCase]);
		fOutputList.Add(dimuon_Timing[iCase]);
    }
    
    // add the list to our output file
    PostData(1,&fOutputList);
    if (fFillTTree) PostData(2,fEventTree);
}

//Event loop
//_________________________________________________________________________________
void AliAnalysisTaskVdmStability::UserExec(Option_t *)
{
    //Get analysis manager
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
    Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
    
    if (isAOD) {
        cout << "Error: AOD not inititalized yet " << endl;
        return;
    }
    
    //Get input handler
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    if (!inputHandler) return;
    
    //get event
    fEvent = static_cast<AliESDEvent*>(InputEvent());
    
    //Check if V0 event
    TString ftc = fEvent->GetFiredTriggerClasses();
    Bool_t isCINT7 = ftc.Contains("CINT7-B");
    Bool_t isCENT = ftc.Contains("CENT");
    if (isCINT7 && isCENT) fIsV0ANDfired = kTRUE;
    else fIsV0ANDfired = kFALSE;
    if (!fIsV0ANDfired) return;
    
    //check if which triggers were fired
    fIsT0fired = fEvent->GetHeader()->IsTriggerInputFired("0TVX");
	fIsEMCALfired = fEvent->GetHeader()->IsTriggerInputFired("0EMC");
	fIsMUfired = fEvent->GetHeader()->IsTriggerInputFired("0MSL");
	fIsDIMUfired = fEvent->GetHeader()->IsTriggerInputFired("0MUL");
    
    //Check physics selection
    fSelectPhysics = kTRUE;
    if(inputHandler){
        if((isESD && inputHandler->GetEventSelection())){
            fSelectPhysics = inputHandler->IsEventSelected();
        }
    }
    
    //set run number
    fRunNumber = static_cast<UInt_t >(InputEvent()->GetRunNumber());
    
    //Set z-vertex information
    const AliESDVertex* trackVtx = fEvent->GetPrimaryVertexTracks();
    fVtxZ = trackVtx->GetZ();
    fnVtxCont = trackVtx->GetNContributors();
    fIsGoodZ = CheckZVtx(fVtxZ, fnVtxCont, 10., kTRUE, 1);
    
    //Set V0 timing
    AliESDVZERO* esdV0 = fEvent->GetVZEROData();
    fTV0A = esdV0->GetV0ATime();
    fTV0C = esdV0->GetV0CTime();
    fGoodTime = CheckTime(fTV0A,fTV0C);
    
    //set pileup
    fPileupEvent = fEvent->IsPileupFromSPD(5,0.8,3.,2.,5.);
    
    //Fill event statistics histograms
    //-----------------------------------------------------
    //Trigger information
    fEventStatV0->Fill(kAllEvents);
    if (fIsT0fired) fEventStatT0->Fill(kAllEvents);
        
    //Physics selection
    if (fSelectPhysics) {
        fEventStatV0->Fill(kSelectedEvents);
        if (fIsT0fired) fEventStatT0->Fill(kSelectedEvents);
    }
    //V0 timing
    if ((fSelectPhysics) && (fGoodTime)) {
        fEventStatV0->Fill(kV0TimeEvents);
        if (fIsT0fired) fEventStatT0->Fill(kV0TimeEvents);
    }
    //pileup
    if ((fSelectPhysics) && !(fPileupEvent)) {
        fEventStatV0->Fill(kPileupEvents);
        if (fIsT0fired) fEventStatT0->Fill(kPileupEvents);
    }
    
    fzCut30 		= CheckZVtx(fVtxZ, fnVtxCont, 30., kFALSE, -999);
    fzCut10 		= CheckZVtx(fVtxZ, fnVtxCont, 10., kFALSE, -999);
    fzCut30nCont0 = CheckZVtx(fVtxZ, fnVtxCont, 30.,  kTRUE, 0);
    fzCut10nCont0 = CheckZVtx(fVtxZ, fnVtxCont, 10.,  kTRUE, 0);
    fzCut30nCont1 = CheckZVtx(fVtxZ, fnVtxCont, 30.,  kTRUE, 1);
    fzCut10nCont1 = CheckZVtx(fVtxZ, fnVtxCont, 10.,  kTRUE, 1);
    
    Bool_t badTimeRun = BadTimingRun(fRunNumber);
    Float_t tV0sum  = fTV0A + fTV0C;
    Float_t tV0diff = fTV0A - fTV0C;
    TString binLabel = Form("%d",fRunNumber);
            
    for (Int_t iCase = 0; iCase < fNSelectionCases; iCase++){
		v0_H[iCase]->GetXaxis()->FindBin(binLabel.Data());
		t0_H[iCase]->GetXaxis()->FindBin(binLabel.Data());
		
		if (IsSelected(iCase)){
			v0_H[iCase]->Fill(binLabel.Data(),1);
			if (!badTimeRun) v0_Timing[iCase]->Fill(tV0diff,tV0sum);
			if (fIsT0fired){
				t0_H[iCase]->Fill(binLabel.Data(),1);
				if (!badTimeRun) t0_Timing[iCase]->Fill(tV0diff,tV0sum);
			} 
			if (fIsEMCALfired){
				emcal_H[iCase]->Fill(binLabel.Data(),1);
				if (!badTimeRun) emcal_Timing[iCase]->Fill(tV0diff,tV0sum);
			}
			if (fIsMUfired){
				muon_H[iCase]->Fill(binLabel.Data(),1);
				if (!badTimeRun) muon_Timing[iCase]->Fill(tV0diff,tV0sum);
			}
			if (fIsDIMUfired){
				dimuon_H[iCase]->Fill(binLabel.Data(),1);
				if (!badTimeRun) dimuon_Timing[iCase]->Fill(tV0diff,tV0sum);
			}
		}
	}
    
    //Fill tree with event information
    if (fFillTTree) fEventTree->Fill();
    
    PostData(1,&fOutputList);
    if (fFillTTree) PostData(2,fEventTree);
    
}

// Set branch names for event tree variables
//_________________________________________________________________________________
void AliAnalysisTaskVdmStability::AddEventTreeVariables(TTree* &tree)
{
    tree->Branch("RunNumber", &fRunNumber);
    //tree->Branch("VtxZ", &fVtxZ);
    //tree->Branch("nVtxCont", &fnVtxCont);
    tree->Branch("PhysSelected",&fSelectPhysics);
    tree->Branch("V0Trigger",&fIsV0ANDfired);
    tree->Branch("T0Trigger",&fIsT0fired);
    //tree->Branch("EMCALTrigger",&fIsEMCALfired);
    //tree->Branch("MUONTrigger",&fIsMUfired);
    //tree->Branch("DIMUONTrigger",&fIsDIMUfired);
    tree->Branch("PileupEvent",&fPileupEvent);
    tree->Branch("timeV0A",&fTV0A);
    tree->Branch("timeV0C",&fTV0C);
}

//Terminate
//_________________________________________________________________________________
void AliAnalysisTaskVdmStability::Terminate(Option_t *)
{
    
}

//Check if good V0 timing
Bool_t AliAnalysisTaskVdmStability::CheckTime(Float_t timeA, Float_t timeC){
    Float_t tV0sum = timeA + timeC;
    Float_t tV0diff = timeA - timeC;
    //return (((5.5 < tV0diff) && (tV0diff < 11.5)) && ((11.5 < tV0sum) && (tV0sum < 17.5)));
    return (((fV0TimeDiff_min < tV0diff) && (tV0diff < fV0TimeDiff_max)) && ((fV0TimeSum_min < tV0sum) && (tV0sum < fV0TimeSum_max)));
}

//Timing Run
Bool_t AliAnalysisTaskVdmStability::BadTimingRun(Int_t run){
	return (run == 253978 || run == 253961 || run == 253958 || run == 253957 || run==253956 || run==253951);
}

//Check zVtx
Bool_t AliAnalysisTaskVdmStability::CheckZVtx(Double_t zVtx, Int_t nCont, Double_t zCut, Bool_t contCut, Int_t nContCut){
    if (TMath::Abs(zVtx) > zCut) return kFALSE; //check zCut, return false if bad z
    if (contCut) return nCont > nContCut;					//if nContCut, return the check
    return kTRUE;
}

Bool_t AliAnalysisTaskVdmStability::IsSelected(Int_t nSelection){
	// NO selection
	if (nSelection == 0) return kTRUE;
	
    //physics selected
	if (nSelection == 1) return (fSelectPhysics);
	
	//V0 timing cut
	if (nSelection == 2) return (fSelectPhysics && fGoodTime);
	
	//pileup rejection
	if (nSelection == 3) return (fSelectPhysics && !fPileupEvent);
	
	//z-vertex cut (standart: +/-10 nCont > 1)
	if (nSelection == 4) return (fSelectPhysics && fIsGoodZ);
	
	//z-vertex cut (+/-30)
	if (nSelection == 5) return (fSelectPhysics && fzCut30);
	
	//z-vertex cut (+/-10)
	if (nSelection == 6) return (fSelectPhysics && fzCut10);
	
	//z-vertex cut (+/-30) + nCont > 0
	if (nSelection == 7) return (fSelectPhysics && fzCut30nCont0);
	
	//z-vertex cut (+/-10) + nCont > 0
	if (nSelection == 8) return (fSelectPhysics && fzCut10nCont0);
	
	//z-vertex cut (+/-30) + nCont > 1
	if (nSelection == 9) return (fSelectPhysics && fzCut30nCont1);
	
	//z-vertex cut (+/-10) + nCont > 1
	if (nSelection == 10) return (fSelectPhysics && fzCut10nCont1);
	
	//V0 timing cut + pileup rejection
	if (nSelection == 11) return (fSelectPhysics && fGoodTime && !fPileupEvent);
                
    //V0 timing cut + z-veretex cut
	if (nSelection == 12) return (fSelectPhysics && fGoodTime && fIsGoodZ);
	
	//V0 timing cut + pileup rejection + z-vertex cut
	if (nSelection == 13) return (fSelectPhysics && fGoodTime && !fPileupEvent && fIsGoodZ);
	
	//V0 timing cut + pileup rejection + z-vertex cut (+/-30)
	if (nSelection == 14) return (fSelectPhysics && fGoodTime && !fPileupEvent && fzCut30);
	
	//V0 timing cut + pileup rejection + z-vertex cut (+/-10)
	if (nSelection == 15) return (fSelectPhysics && fGoodTime && !fPileupEvent && fzCut10);
	
	//V0 timing cut + pileup rejection + z-vertex cut (+/-30) + nCont > 0
	if (nSelection == 16) return (fSelectPhysics && fGoodTime && !fPileupEvent && fzCut30nCont0);
	
	//V0 timing cut + pileup rejection + z-vertex cut (+/-10) + nCont > 0
	if (nSelection == 17) return (fSelectPhysics && fGoodTime && !fPileupEvent && fzCut10nCont0);
	
	//V0 timing cut + pileup rejection + z-vertex cut (+/-30) + nCont > 1
	if (nSelection == 18) return (fSelectPhysics && fGoodTime && !fPileupEvent && fzCut30nCont1);
	
	//V0 timing cut + pileup rejection + z-vertex cut (+/-10) + nCont > 1
	if (nSelection == 19) return (fSelectPhysics && fGoodTime && !fPileupEvent && fzCut10nCont1);
	
	//V0 timing cut + z-vertex cut (+/-30) + nCont > 0
	if (nSelection == 20) return (fSelectPhysics && fGoodTime && fzCut30nCont0);
	
	//V0 timing cut
	if (nSelection == 21) return (fGoodTime);
	
	//V0 timing cut + PS
	if (nSelection == 22) return (fSelectPhysics && fGoodTime);
	
	//V0 timing cut + PU
	if (nSelection == 23) return (fGoodTime && !fPileupEvent);
	
	//V0 timing cut + PU + PS
	if (nSelection == 24) return (fSelectPhysics && fGoodTime && !fPileupEvent);
	
	printf("\n\n!!! WARNING - NO SELECTION CUTS FOR SELECTION CASE %i!!!\n\n",nSelection);
	return kFALSE;
}
