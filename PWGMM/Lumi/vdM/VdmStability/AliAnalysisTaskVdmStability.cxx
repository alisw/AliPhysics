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

// author: I.Lofnes, ingrid.mckibben.lofnes@cern.ch
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
#include <AliTriggerAnalysis.h>
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
fPileupEvent(kFALSE),
fTV0A(-10240.0f),
fTV0C(-10240.0f),
fGoodTime(kFALSE),
fEventTag(0),
fEvent(0x0),
fEventStatV0(0x0),
fEventStatT0(0x0),
fEventTree(0x0),
fNRuns(1000),
fNSelectionCases(20)
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
fPileupEvent(kFALSE),
fTV0A(-10240.0f),
fTV0C(-10240.0f),
fGoodTime(kFALSE),
fEventTag(0),
fEvent(0x0),
fEventStatV0(0x0),
fEventStatT0(0x0),
fEventTree(0x0),
fNRuns(1000),
fNSelectionCases(20)
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
    if (v0_H)         		{ delete v0_H; }
    if (t0_H)         		{ delete t0_H; }
    if (v0_Timing)         { delete v0_Timing; }
    if (t0_Timing)         { delete t0_Timing; }
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
    
	TString selectionCases[fNSelectionCases] = {
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
		"V0_timing_cut_pileup_rejection_z-vertex_10nContCut2"
	};
	
    for (Int_t iCase = 0; iCase < fNSelectionCases; iCase++){
        v0_H[iCase] = new TH1D(Form("v0_H_%s",selectionCases[iCase].Data()),Form("V0 events, %s",selectionCases[iCase].Data()),fNRuns,-0.5,fNRuns-0.5);
        t0_H[iCase] = new TH1D(Form("t0_H_%s",selectionCases[iCase].Data()),Form("T0 events, %s",selectionCases[iCase].Data()),fNRuns,-0.5,fNRuns-0.5);
        v0_Timing[iCase] = new TH2D(Form("v0_Timing_%s",selectionCases[iCase].Data()),"sum vs diff",500,-30,40,500,-20,50);
        t0_Timing[iCase] = new TH2D(Form("t0_Timing_%s",selectionCases[iCase].Data()),"sum vs diff",500,-30,40,500,-20,50);
		fOutputList.Add(v0_H[iCase]);
		fOutputList.Add(t0_H[iCase]);
		fOutputList.Add(v0_Timing[iCase]);
		fOutputList.Add(t0_Timing[iCase]);
    }
    
    // add the list to our output file
    PostData(1,&fOutputList);
    PostData(2,fEventTree);
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
    if (ftc.Contains("CINT7-B"))
        fIsV0ANDfired = kTRUE;
    else fIsV0ANDfired = kFALSE;
    if (!fIsV0ANDfired) return;
    
    //check if T0 trigger was fired
    fIsT0fired = fEvent->GetHeader()->IsTriggerInputFired("0TVX");
    
    //Check physics selection
    fSelectPhysics = kTRUE;
    if( inputHandler){
        if((isESD && inputHandler->GetEventSelection())){
            fSelectPhysics = inputHandler->IsEventSelected();
        }
    }
    
    //set run number
    fRunNumber = static_cast<UInt_t >(InputEvent()->GetRunNumber());
    
    //Set z-vertex information
    const AliESDVertex* trackVtx    = fEvent->GetPrimaryVertexTracks();
    fVtxZ = trackVtx->GetZ();
    fnVtxCont = trackVtx->GetNContributors();
    fIsGoodZ = CheckZVtx(fVtxZ, fnVtxCont, 10., kTRUE, 1);
    
    //Set V0 timing
    AliESDVZERO * esdV0 = fEvent->GetVZEROData();
    fTV0A = esdV0->GetV0ATime();
    fTV0C = esdV0->GetV0CTime();
    fGoodTime = CheckTime(fTV0A,fTV0C);
    
    //set pileup
    fPileupEvent = fEvent->IsPileupFromSPD(3,0.8,3.,2.,5.);
    
    //Fill event statistics histograms
    //-----------------------------------------------------
    //Trigger information
    fEventStatV0->Fill(kAllEvents);
    if (fIsT0fired)
        fEventStatT0->Fill(kAllEvents);
        
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
    //Z-vertex
    if ((fSelectPhysics) && (fIsGoodZ)) {
        fEventStatV0->Fill(kGoodZEvents);
        if (fIsT0fired) fEventStatT0->Fill(kGoodZEvents);
    }
    //Timing + PU
    if ((fSelectPhysics) && (fGoodTime) && !(fPileupEvent) ) {
        fEventStatV0->Fill(kV0andPUEvents);
        if (fIsT0fired) fEventStatT0->Fill(kV0andPUEvents);
    }
    //Timing + z
    if ((fSelectPhysics) && (fGoodTime) && (fIsGoodZ) ) {
        fEventStatV0->Fill(kV0andZEvents);
        if (fIsT0fired) fEventStatT0->Fill(kV0andZEvents);
    }
    //Timing + PU + z
    if ((fSelectPhysics) && (fGoodTime) && !(fPileupEvent) && (fIsGoodZ) ) {
        fEventStatV0->Fill(kV0andPUandZEvents);
        if (fIsT0fired) fEventStatT0->Fill(kV0andPUandZEvents);
    }
    
    Bool_t zCut30 		= CheckZVtx(fVtxZ, fnVtxCont, 30., kFALSE, -999);
    Bool_t zCut10 		= CheckZVtx(fVtxZ, fnVtxCont, 10., kFALSE, -999);
    Bool_t zCut30nCont0 = CheckZVtx(fVtxZ, fnVtxCont, 30.,  kTRUE, 0);
    Bool_t zCut10nCont0 = CheckZVtx(fVtxZ, fnVtxCont, 10.,  kTRUE, 0);
    Bool_t zCut30nCont1 = CheckZVtx(fVtxZ, fnVtxCont, 30.,  kTRUE, 1);
    Bool_t zCut10nCont1 = CheckZVtx(fVtxZ, fnVtxCont, 10.,  kTRUE, 1);
    Bool_t badTimeRun = BadTimingRun(fRunNumber);
    Float_t tV0sum  = fTV0A + fTV0C;
    Float_t tV0diff = fTV0A - fTV0C;
    TString binLabel = Form("%d",fRunNumber);
            
    for (Int_t iCase; iCase < fNSelectionCases; iCase++){
		v0_H[iCase]->GetXaxis()->FindBin(binLabel.Data());
		t0_H[iCase]->GetXaxis()->FindBin(binLabel.Data());
	}
	
	// NO selection
	v0_H[0]->Fill(binLabel.Data(),1);
	if (!badTimeRun) v0_Timing[0]->Fill(tV0diff,tV0sum);
	if (fIsT0fired){
		t0_H[0]->Fill(binLabel.Data(),1);
		if (!badTimeRun) t0_Timing[0]->Fill(tV0diff,tV0sum);
	}  
    //physics selected
    if (fSelectPhysics){
		v0_H[1]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[1]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[1]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[1]->Fill(tV0diff,tV0sum);
		}
	}
	//V0 timing cut
	if (fSelectPhysics && fGoodTime){
		v0_H[2]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[2]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[2]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[2]->Fill(tV0diff,tV0sum);
		}
	}
	//pileup rejection
	if (fSelectPhysics && !fPileupEvent){
		v0_H[3]->Fill(binLabel.Data(),1);
		v0_Timing[3]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[3]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[3]->Fill(tV0diff,tV0sum);
		}
	}
	//z-vertex cut (standart: +/-10 nCont > 1)
	if (fSelectPhysics && fIsGoodZ){
		v0_H[4]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[4]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[4]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[4]->Fill(tV0diff,tV0sum);
		}
	}
	//z-vertex cut (+/-30)
	if (fSelectPhysics && zCut30){
		v0_H[5]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[5]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[5]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[5]->Fill(tV0diff,tV0sum);
		}
	}
	//z-vertex cut (+/-10)
	if (fSelectPhysics && zCut10){
		v0_H[6]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[6]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[6]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[6]->Fill(tV0diff,tV0sum);
		}
	}
	//z-vertex cut (+/-30) + nCont > 0
	if (fSelectPhysics && zCut30nCont0){
		v0_H[7]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[7]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[7]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[7]->Fill(tV0diff,tV0sum);
		}
	}
	//z-vertex cut (+/-10) + nCont > 0
	if (fSelectPhysics && zCut10nCont0){
		v0_H[8]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[8]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[8]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[8]->Fill(tV0diff,tV0sum);
		}
	}
	//z-vertex cut (+/-30) + nCont > 1
	if (fSelectPhysics && zCut30nCont1){
		v0_H[9]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[9]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[9]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[9]->Fill(tV0diff,tV0sum);
		}
	}
	//z-vertex cut (+/-10) + nCont > 1
	if (fSelectPhysics && zCut10nCont1){
		v0_H[10]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[10]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[10]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[10]->Fill(tV0diff,tV0sum);
		}
	}
	//V0 timing cut + pileup rejection
	if (fSelectPhysics && fGoodTime && !fPileupEvent){
		v0_H[11]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[11]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[11]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[11]->Fill(tV0diff,tV0sum);
		}
	}
                
    //V0 timing cut + z-veretex cut
    if (fSelectPhysics && fGoodTime && fIsGoodZ){
		v0_H[12]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[12]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[12]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[12]->Fill(tV0diff,tV0sum);
		}
	}
	
	//V0 timing cut + pileup rejection + z-vertex cut
	if (fSelectPhysics && fGoodTime && !fPileupEvent && fIsGoodZ){
		v0_H[13]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[13]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[13]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[13]->Fill(tV0diff,tV0sum);
		}
	}
	
	//V0 timing cut + pileup rejection + z-vertex cut (+/-30)
	if (fSelectPhysics && fGoodTime && !fPileupEvent && zCut30){
		v0_H[14]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[14]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[14]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[14]->Fill(tV0diff,tV0sum);
		}
	}
	
	//V0 timing cut + pileup rejection + z-vertex cut (+/-10)
	if (fSelectPhysics && fGoodTime && !fPileupEvent && zCut10){
		v0_H[15]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[15]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[15]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[15]->Fill(tV0diff,tV0sum);
		}
	}
	
	//V0 timing cut + pileup rejection + z-vertex cut (+/-30) + nCont > 0
	if (fSelectPhysics && fGoodTime && !fPileupEvent && zCut30nCont0){
		v0_H[16]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[16]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[16]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[16]->Fill(tV0diff,tV0sum);
		}
	}
	
	//V0 timing cut + pileup rejection + z-vertex cut (+/-10) + nCont > 0
	if (fSelectPhysics && fGoodTime && !fPileupEvent && zCut10nCont0){
		v0_H[17]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[17]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[17]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[17]->Fill(tV0diff,tV0sum);
		}
	}
	
	//V0 timing cut + pileup rejection + z-vertex cut (+/-30) + nCont > 1
	if (fSelectPhysics && fGoodTime && !fPileupEvent && zCut30nCont1){
		v0_H[18]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[18]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[18]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[18]->Fill(tV0diff,tV0sum);
		}
	}
	
	//V0 timing cut + pileup rejection + z-vertex cut (+/-10) + nCont > 1
	if (fSelectPhysics && fGoodTime && !fPileupEvent && zCut10nCont1){
		v0_H[19]->Fill(binLabel.Data(),1);
		if (!badTimeRun) v0_Timing[19]->Fill(tV0diff,tV0sum);
		if (fIsT0fired){
			t0_H[19]->Fill(binLabel.Data(),1);
			if (!badTimeRun) t0_Timing[19]->Fill(tV0diff,tV0sum);
		}
	}
    
    //Fill tree with event information
    fEventTree->Fill();
    
    PostData(1,&fOutputList);
    PostData(2,fEventTree);
    
}

// Set branch names for event tree variables
//_________________________________________________________________________________
void AliAnalysisTaskVdmStability::AddEventTreeVariables(TTree* &tree)
{
    tree->Branch("RunNumber", &fRunNumber);
    tree->Branch("VtxZ", &fVtxZ);
    tree->Branch("nVtxCont", &fnVtxCont);
    tree->Branch("PhysSelected",&fSelectPhysics);
    tree->Branch("V0Trigger",&fIsV0ANDfired);
    tree->Branch("T0Trigger",&fIsT0fired);
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
    return (((4 < tV0diff) && (tV0diff < 12)) && ((10 < tV0sum) && (tV0sum < 18)));
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
