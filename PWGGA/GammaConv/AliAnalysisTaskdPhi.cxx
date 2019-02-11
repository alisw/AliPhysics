/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *																          *
 * Authors: Svein Lindal                                                  *
 * Version 1.0			                                                  *
 *																          *
 *																          *
 * Permission to use, copy, modify and distribute this software and its	  *
 * documentation strictly for non-commercial purposes is hereby granted	  *
 * without fee, provided that the above copyright notice appears in all	  *
 * copies and that both the copyright notice and this permission notice	  *
 * appear in the supporting documentation. The authors make no claims	  *
 * about the suitability of this software for any purpose. It is		  *
 * provided "as is" without express or implied warranty.				  *
 **************************************************************************/

////////////////////////////////////////////////
//--------------------------------------------- 
// Class doing conversion gamma dPhi correlations
// Gamma Conversion analysis
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskdPhi.h"

#include <TH2I.h>
#include <TList.h>
#include <TChain.h>
#include <TFile.h>
#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <TGeoGlobalMagField.h>

#include "AliConversionTrackCuts.h"
#include "AliConversionCuts.h"
#include "AliConversionMesonCuts.h"
#include "AliAODConversionPhoton.h"
#include "AliAODConversionMother.h"
#include "TGrid.h"
// #include "AliAnaConvCorrPhoton.h"
// #include "AliAnaConvCorrPion.h"
// #include "AliAnaConvIsolation.h"
#include "AliV0ReaderV1.h"
// Author Svein Lindal <slindal@fys.uio.no>
using namespace std;

ClassImp(AliAnalysisTaskdPhi)


//________________________________________________________________________
AliAnalysisTaskdPhi::AliAnalysisTaskdPhi(const char *name) : AliAnalysisTaskSE(name),
	fHistograms(NULL),
	fCorrSparse(NULL),
	fTrigSparse(NULL),
	fTrackSparse(NULL),
	fMassSparse(NULL),
	fV0Reader(NULL),
	fSaveReaderHists(kFALSE),
	fV0FilterEvent(NULL),
	fV0FilterPhoton(NULL),
	fV0Filters(),
	fEventFilters(),
	fEventFilter(NULL),
	fPhotonFilter(NULL),
	fMesonFilter(NULL),
	fMesonFilters(),
	fTrackFilter(NULL),
	fTrackFilters(),
	fGammas(),
	fTracks(),
	hMEvents(NULL),
	hTrackCent(NULL),
	hTrigPt(NULL),
	hTrackPt(NULL), 
	hTrigPhi(NULL),
	fDeltaAODBranchName("AliAODGammaConversion_gamma"), 
	fAxistPt(),
	fAxiscPt(),
	fAxisdEta(),
	fAxisTrigEta(),
	fAxisAssEta(),
	fAxisdPhi(),
	fAxisCent(),
	fAxisZ(), 
	fAxisPiM(), 
	fAxisTrackFilters(),							     
	fAxisV0Filters(),							     
	fAxisMesonFilters(),							     
	fkTrackAxis(kFALSE),
	fkV0Axis(kFALSE),
	fkPionAxis(kFALSE),
	fAxesList(), 
	fTrigAxesList(), 
	fTrackAxesList(), 
	fMassAxesList(),
	fDoPhoton(kFALSE), 
	fCorrectionMap(NULL)
{
	//constructor
	SetUpBins();

	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());

	fGammas.SetOwner(kTRUE);
	fTracks.SetOwner(kTRUE);
}



//________________________________________________________________________
AliAnalysisTaskdPhi::~AliAnalysisTaskdPhi(){
	//destructor

	if(fV0FilterEvent)
		delete fV0FilterEvent;
	fV0FilterEvent = NULL;

	if(fV0FilterPhoton)
		delete fV0FilterPhoton;
	fV0FilterPhoton = NULL;

	if(fMesonFilter)
		delete fMesonFilter;
	fMesonFilter = NULL;

	if(fEventFilter)
		delete fEventFilter;
	fEventFilter = NULL;

	if(fPhotonFilter)
		delete fPhotonFilter;
	fPhotonFilter = NULL;

	if(fHistograms)
		delete fHistograms;
	fHistograms = NULL;

	if(fTrackFilter)
		delete fTrackFilter;
	fTrackFilter = NULL;

	fGammas.Delete();
	fTracks.Delete();

}

///________________________________________________________________________
void AliAnalysisTaskdPhi::SetUpBins() {

	fAxisTrigEta.SetNameTitle("tEta", "Eta");
	fAxisTrigEta.Set(320, -0.8, 0.8);

	fAxisAssEta.SetNameTitle("aEta", "Eta");
	fAxisAssEta.Set(360, -0.9, 0.9);

	fAxisdEta.SetNameTitle("dEta", "Eta");
	fAxisdEta.Set(34, -1.7, 1.7);

	fAxisdPhi.SetNameTitle("dPhi", "delta Phi");
	fAxisdPhi.Set(32, -TMath::PiOver2(), 3*TMath::PiOver2());

	Double_t tptbins[10] = {0.1, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15, 50, 100};                                    fAxistPt.SetNameTitle("tPt", "tPt");
	fAxistPt.Set(9, tptbins);

	Double_t cptbins[13] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 25, 50, 100};                     fAxiscPt.SetNameTitle("cPt", "cPt");
	fAxiscPt.Set(12, cptbins);

	fAxisZ.SetNameTitle("vertexz", "Z");
	fAxisZ.Set(4, -10, 10);

	Double_t centbins[5] = {0, 10, 30, 60, 100.1};
	fAxisCent.SetNameTitle("centrality", "Cent");
	fAxisCent.Set(4, centbins);

	fAxisPiM.SetNameTitle("InvMassPi0", "Invariant mass");
	Double_t mbins[7] = {0.1, 0.11, 0.12, 0.15, 0.16, 0.18, 0.2};
	fAxisPiM.Set(6, mbins);
}


///________________________________________________________________________
// void AliAnalysisTaskdPhi::SetUpCorrObjects() {
   // //Set up corr objects
   // AliDebug(AliLog::kDebug + 5, "Set Up corr objects");

   // if(fDoPhoton) {
   //   fPhotonCorr = new AliAnaConvCorrPhoton("PhotonCorr","photon %s");
   //   SetUpCorrAxes(fPhotonCorr);
   //   fPhotonCorr->CreateHistograms();
   //   fHistograms->Add(fPhotonCorr->GetHistograms());
   // }

   // fPionCorr = new AliAnaConvCorrPion("PionCorr", "pion");
   // SetUpCorrAxes(fPionCorr);
   // fPionCorr->GetAxisM().Set(fAxisPiM.GetNbins(), fAxisPiM.GetXbins()->GetArray());
   // fPionCorr->CreateHistograms();
   // fHistograms->Add(fPionCorr->GetHistograms());
   //}

 ///________________________________________________________________________    
// void AliAnalysisTaskdPhi::SetUpCorrAxes(AliAnaConvCorrBase * corr) {
   ///Set up axes in corr object
   // corr->GetAxisCent().Set(fAxisCent.GetNbins(), fAxisCent.GetXbins()->GetArray());
   // const Double_t * zbins = fAxisZ.GetXbins()->GetArray();
   // if(zbins) {
   //   corr->GetAxisZ().Set(fAxisZ.GetNbins(), fAxisZ.GetXbins()->GetArray());
   // } else {
   //   corr->GetAxisZ().Set(fAxisZ.GetNbins(), fAxisZ.GetBinLowEdge(1), fAxisZ.GetBinUpEdge(fAxisZ.GetNbins()));
   // }

   // corr->GetAxistPt().Set(fAxistPt.GetNbins(), fAxistPt.GetXbins()->GetArray());
   // corr->GetAxiscPt().Set(fAxiscPt.GetNbins(), fAxiscPt.GetXbins()->GetArray());
   //   corr->GetAxisdEta().Set(fAxisdEta.GetNbins(), fAxisdEta.GetBinLowEdge(1), fAxisdEta.GetBinUpEdge(fAxisdEta.GetNbins()));
   // corr->GetAxisTrigEta().Set(fAxisTrigEta.GetNbins(), fAxisTrigEta.GetBinLowEdge(1), fAxisTrigEta.GetBinUpEdge(fAxisTrigEta.GetNbins()));
   // corr->GetAxisAssEta().Set(fAxisAssEta.GetNbins(), fAxisAssEta.GetBinLowEdge(1), fAxisAssEta.GetBinUpEdge(fAxisAssEta.GetNbins()));
// }    


//________________________________________________________________________
void AliAnalysisTaskdPhi::UserCreateOutputObjects() {
	// Create histograms
	// TGrid::Connect("alien://",0,0,"t");
	// if(!gGrid) AliWarning("no GGrid");
	// TFile *tfile = TFile::Open("alien:///alice/cern.ch/user/s/slindal/trackMap.root", "READ");
	// if(tfile) {
	//   THnF * corrmap = dynamic_cast<THnF*>(tfile->Get("hTrackCorr"));
	//   if (corrmap) {
	//     fCorrectionMap = dynamic_cast<THnF*>(THn::CreateHn("corr", "corr", corrmap));
	//     for(Int_t i = 0; i < fCorrectionMap->GetNdimensions(); i++) {
	// 	TAxis * axis = fCorrectionMap->GetAxis(i);
	// 	axis->SetRange(1, axis->GetNbins());
	//     }
		
	//     cout << "yessssssssssssssssssssssssssssssssssssssssssssssssss"<<endl;
	//   } else {
	//     cout << "xxxxxxxxxxxxxxxxx xxxxxxxxxxxxxxx xxxxxxxxxxxxxxxx"<<endl;
	//   }
	//   tfile->Close();
	// } else {
	//   cout << "no tfile shit shit shit "<<endl;
	//   AliFatal("file not ther!!!");
	// }

	
	fHistograms = new TList();
	fHistograms->SetName("dPhi_histograms");
	fHistograms->SetOwner(kTRUE);
	
	
	if(!fV0Reader){
		fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
	}
	
	if(!fV0Reader){
		printf("Error: No V0 Reader");
		return;
	} // GetV0Reader
	
	
	if(fSaveReaderHists) {
		AliConversionPhotonCuts * v0cuts = fV0Reader->GetConversionCuts();
		if(v0cuts) {
			TList * histograms = v0cuts->GetCutHistograms();
			if(!histograms) {
				AliWarning("initializing v0 reader hists");
				v0cuts->InitCutHistograms("V0Reader", kTRUE);
			}
			histograms = v0cuts->GetCutHistograms();
			if(histograms) {
				fHistograms->Add(histograms);
			}
		}
	}
	
	for(Int_t igf = 0; igf < fV0Filters[0].GetEntriesFast(); igf ++){
		AliConversionPhotonCuts * f = dynamic_cast<AliConversionPhotonCuts*>(fV0Filters[0].At(igf));
		if(f) {
			TList * histograms = f->GetCutHistograms();
			if(histograms) fHistograms->Add(f->GetCutHistograms());
		}
	}

	for(Int_t igf = 0; igf < fV0Filters[1].GetEntriesFast(); igf ++){
		AliConversionPhotonCuts * f = dynamic_cast<AliConversionPhotonCuts*>(fV0Filters[1].At(igf));
		if(f) {
			TList * histograms = f->GetCutHistograms();
			if(histograms) fHistograms->Add(f->GetCutHistograms());
		}
	}

	for(Int_t igf = 0; igf < fMesonFilters[0].GetEntriesFast(); igf ++){
		AliConversionMesonCuts * f = dynamic_cast<AliConversionMesonCuts*>(fMesonFilters[0].At(igf));
		if(f) {
			TList * histograms = f->GetCutHistograms();
			if(histograms) fHistograms->Add(f->GetCutHistograms());
		}
	}

	for(Int_t igf = 0; igf < fMesonFilters[1].GetEntriesFast(); igf ++){
		AliConversionMesonCuts * f = dynamic_cast<AliConversionMesonCuts*>(fMesonFilters[1].At(igf));
		if(f) {
			TList * histograms = f->GetCutHistograms();
			if(histograms) fHistograms->Add(f->GetCutHistograms());
		}
	}

	if(fV0FilterEvent) {
		fV0FilterEvent->InitCutHistograms("V0FilterEvent", kFALSE);
		fHistograms->Add(fV0FilterEvent->GetCutHistograms());
	}
	if(fV0FilterPhoton) {
		fV0FilterPhoton->InitCutHistograms("V0FilterPhoton", kFALSE);
		fHistograms->Add(fV0FilterPhoton->GetCutHistograms());
	}
	if(fMesonFilter) {
		fMesonFilter->InitCutHistograms("PionFilter", kFALSE);
		fHistograms->Add(fMesonFilter->GetCutHistograms());
	}
	if(fPhotonFilter) {
		fPhotonFilter->InitCutHistograms("PhotonFilter", kFALSE);
		fHistograms->Add(fPhotonFilter->GetCutHistograms());
	}


	AliConversionTrackCuts * tc = dynamic_cast<AliConversionTrackCuts*>(fTrackFilter);
	if(tc) {
		fHistograms->Add(tc->CreateHistograms());
	}


	// for(Int_t i = 0; i < fTrackFilters.GetEntriesFast(); i++) {
	//   AliConversionTrackCuts * tc = dynamic_cast<AliConversionTrackCuts*>(fTrackFilters.At(i));
	//   if(tc) fHistograms->Add(tc->CreateHistograms());
	// }

	//SetUpCorrObjects();

	///Set up ME histograms
	TList * MEHistograms = new TList();
	MEHistograms->SetName("MEHistograms");
	MEHistograms->SetOwner(kTRUE);
	fHistograms->Add(MEHistograms);
	
	hMEvents = new TH2I("hMEvents", "Nevents vs centrality vertexz",
				fAxisZ.GetNbins(), fAxisZ.GetXbins()->GetArray(),
				fAxisCent.GetNbins(), fAxisCent.GetXbins()->GetArray());
	MEHistograms->Add(hMEvents);

	hTrackCent = new TH2I("hTrackCent", "N accepted tracks vs centrality",
				fAxisCent.GetNbins() > 1 ? (int) (10*(fAxisCent.GetXmax() - fAxisCent.GetXmin()))  : 1, 
				fAxisCent.GetXmin(), fAxisCent.GetXmax(),
				fAxisCent.GetNbins() > 1 ? 900 : 50, 
				0,
				fAxisCent.GetNbins() > 1 ? 1800 : 50);
	MEHistograms->Add(hTrackCent);

	hTrigPt = new TH3F("hTrigPt", "trigger pt", 100, 0., 10., 
				10, 0., 50., 
				5,  0.05, 0.2);
	MEHistograms->Add(hTrigPt);
	hTrackPt = new TH2F("hTrackPt", "track pt", 100, 0, 10, 10, 0, 50);//fAxisCent.GetNbins(), fAxisCent.GetXbins()->GetArray()); 
	MEHistograms->Add(hTrackPt);
	hTrigPhi = new TH1F("hTrigPhi", "trigger pt", 32, 0, 2*TMath::Pi());
	MEHistograms->Add(hTrigPhi);



	Int_t ntrackfilters[2] = {fTrackFilters[0].GetEntriesFast(), fTrackFilters[1].GetEntriesFast()};
	fkTrackAxis = kTRUE;
	fAxisTrackFilters.SetNameTitle("trackCuts", "trackCuts");
	fAxisTrackFilters.Set(ntrackfilters[0] + ntrackfilters[1] + 1, -ntrackfilters[0] -0.5, ntrackfilters[1]  + 0.5);

	Int_t nV0filters[2] = {fV0Filters[0].GetEntriesFast(), fV0Filters[1].GetEntriesFast()};
	fkV0Axis = kTRUE;
	fAxisV0Filters.SetNameTitle("V0Cuts", "V0Cuts");
	fAxisV0Filters.Set(nV0filters[0] + nV0filters[1] + 1, -nV0filters[0] -0.5, nV0filters[1]  + 0.5);
	
	Int_t nmesonfilters[2] = {fMesonFilters[0].GetEntriesFast(), fMesonFilters[1].GetEntriesFast()};
	fkPionAxis = kTRUE;
	fAxisMesonFilters.SetNameTitle("mesonCuts", "mesonCuts");
	fAxisMesonFilters.Set(nmesonfilters[0] + nmesonfilters[1] + 1, -nmesonfilters[0] -0.5, nmesonfilters[1]  + 0.5);
	
	fAxesList.AddAt(&fAxisdEta, 0);
	fAxesList.AddAt(&fAxisdPhi, 1);
	fAxesList.AddAt(&fAxistPt, 2);
	fAxesList.AddAt(&fAxiscPt, 3);
	fAxesList.AddAt(&fAxisCent, 4);
	fAxesList.AddAt(&fAxisZ, 5);
	fAxesList.AddAt(&fAxisPiM, 6);
	fAxesList.AddAt(&fAxisTrackFilters, 7);
	fAxesList.AddAt(&fAxisV0Filters, 8);
	fAxesList.AddAt(&fAxisMesonFilters, 9);
	
	fTrackAxesList.AddAt(&fAxisAssEta, 0);
	fTrackAxesList.AddAt(&fAxistPt, 1);
	fTrackAxesList.AddAt(&fAxiscPt, 2);
	fTrackAxesList.AddAt(&fAxisCent, 3);
	fTrackAxesList.AddAt(&fAxisZ, 4);
	//fTrackAxesList.AddAt(&fAxisPiM, 5);
	fTrackAxesList.AddAt(&fAxisTrackFilters, 5);
	fTrackAxesList.AddAt(&fAxisV0Filters, 6);
	fTrackAxesList.AddAt(&fAxisMesonFilters, 7);
	
	fTrigAxesList.AddAt(&fAxisTrigEta, 0);
	fTrigAxesList.AddAt(&fAxistPt, 1);
	fTrigAxesList.AddAt(&fAxisCent, 2);
	fTrigAxesList.AddAt(&fAxisZ, 3);
	fTrigAxesList.AddAt(&fAxisPiM, 4);
	fTrigAxesList.AddAt(&fAxisV0Filters, 5);
	fTrigAxesList.AddAt(&fAxisMesonFilters, 6);
	
	TList masslist;
	TAxis massax;
	massax.SetNameTitle("mass", "mass");
	massax.Set(360, 0.04, 0.4); //hardcoded! change also in filling!
	
	masslist.AddAt(&massax, 0);
	masslist.AddAt(&fAxistPt, 1);
	masslist.AddAt(&fAxisCent, 2);
	masslist.AddAt(&fAxisV0Filters, 3);
	masslist.AddAt(&fAxisMesonFilters, 4);
	
	fCorrSparse = CreateSparse(TString("pionSparse"), TString("pionSparse"), &fAxesList);
	fTrackSparse = CreateSparse(TString("trackSparse"), TString("trackSparse"), &fTrackAxesList);
	fTrigSparse = CreateSparse(TString("trigSparse"), TString("trigSparse"), &fTrigAxesList);
	fMassSparse = CreateSparse("massSparse", "massSparse", &masslist);
	
	fHistograms->Add(fCorrSparse);
	fHistograms->Add(fTrackSparse);
	fHistograms->Add(fTrigSparse);
	fHistograms->Add(fMassSparse);


	///Add gamma and track containers:
	for(Int_t i = 0; i < fV0Filters[1].GetEntriesFast() + 1; i++) {
		fGammas.Add(new TObjArray());
	}

	for(Int_t i = 0; i < fTrackFilters[1].GetEntriesFast() + 1; i++) {
		fTracks.Add(new TObjArray());
	}
	PostData(1, fHistograms);
}

///________________________________________________________________________
THnSparseF * AliAnalysisTaskdPhi::CreateSparse(TString nameString, TString titleString, TList * axesList) {
	//Create sparse
	const Int_t dim = axesList->GetSize();
	
	TAxis * axes[dim];
	Int_t   bins[dim];
	Double_t min[dim];
	Double_t max[dim];

	for(Int_t i = 0; i<dim; i++) {
		TAxis * axis = dynamic_cast<TAxis*>(axesList->At(i));
		if(axis) axes[i] = axis;
		else {
		cout << "AliAnalysisTaskdPhi::CreateSparse: Error error, all the axes are not present in axis list" << endl;
		return NULL;
		}
	}
	
	for(Int_t i = 0; i<dim; i++) {
		//cout << axes[i]->GetTitle() << endl;
		bins[i] = axes[i]->GetNbins(); 
		min[i] = axes[i]->GetBinLowEdge(1);
		max[i] = axes[i]->GetBinUpEdge(axes[i]->GetNbins());
	}

	THnSparseF * sparse = new THnSparseF(Form("%s", nameString.Data()), 
						Form("%s", titleString.Data()), 
						dim, bins, min, max);
	
	for(Int_t i = 0; i<dim; i++) {
		sparse->GetAxis(i)->SetNameTitle(axes[i]->GetName(), axes[i]->GetTitle() );
		if(axes[i]->GetXbins()->GetSize() > 0) {
		sparse->SetBinEdges(i, axes[i]->GetXbins()->GetArray() );
		}
	}
	return sparse;
}

//________________________________________________________________________
void AliAnalysisTaskdPhi::UserExec(Option_t *) {
	///User exec. 
	///This is a very ugly function, cut the complexity of the logic demands it.
	

	AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
	if(esdEvent) {
		if (!TGeoGlobalMagField::Instance()->GetField()) esdEvent->InitMagneticField(); 
	}


	//if(! fV0Filter->EventIsSelected(fInputEvent)) return;
	if(!fV0Reader){
		AliError("Error: No V0 Reader");
		return;
	} // GetV0Reader
	
	if(!fV0Reader->IsEventSelected()) {
		return;
	}

	AliDebug(5, "Processing event");
	
	
	for(Int_t i = 0; i < fGammas.GetEntriesFast(); i++) {
		static_cast<TObjArray*>(fGammas.At(i))->Clear();
	}

	for(Int_t i = 0; i < fTracks.GetEntriesFast(); i++) {
		static_cast<TObjArray*>(fTracks.At(i))->Clear();
	}


	
	AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
	Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
	
	AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
	if (!inputHandler) {
		cout << "cout no input event handler"<<endl;
		return;
	}

	for(Int_t igf = 0; igf < fV0Filters[0].GetEntriesFast(); igf++) {
		AliConversionPhotonCuts * f = dynamic_cast<AliConversionPhotonCuts*>(fV0Filters[0].At(igf));
		if ( f && !f->GetPIDResponse() ) {
			if ( inputHandler->GetPIDResponse() ){
				f->SetPIDResponse( inputHandler->GetPIDResponse() );
			} 
		} else {
			break;
		}
	}

	for(Int_t igf = 0; igf < fV0Filters[1].GetEntriesFast(); igf++) {
		AliConversionPhotonCuts * f = dynamic_cast<AliConversionPhotonCuts*>(fV0Filters[1].At(igf));
		if ( f && !f->GetPIDResponse() ) {
			if ( inputHandler->GetPIDResponse() ){
				f->SetPIDResponse( inputHandler->GetPIDResponse() );
			} 
		} else {
			break;
		}
	}

	
	if ( fV0FilterPhoton && !fV0FilterPhoton->GetPIDResponse() ) {
		if ( inputHandler->GetPIDResponse() ){
			fV0FilterPhoton->SetPIDResponse( inputHandler->GetPIDResponse() );
		}
	}


	///Initialize track cuts. Delete tracks that have been constrained to vertex (copies)
	AliConversionTrackCuts * tc = dynamic_cast<AliConversionTrackCuts*>(fTrackFilter);
	if(tc) {
		tc->SetEvent(fInputEvent);
	}
	
	for(Int_t i = 0; i < fTrackFilters[0].GetEntriesFast(); i++){
		AliConversionTrackCuts * tct = dynamic_cast<AliConversionTrackCuts*>(fTrackFilters[0].At(i));
		if(tct) {
			tct->SetEvent(fInputEvent);
		}
	}
	for(Int_t i = 0; i < fTrackFilters[1].GetEntriesFast(); i++){
		AliConversionTrackCuts * tct = dynamic_cast<AliConversionTrackCuts*>(fTrackFilters[1].At(i));
		if(tct) {
			tct->SetEvent(fInputEvent);
		}
	}
	
	Double_t centrality = 0.0;
	Double_t vertexz = fInputEvent->GetPrimaryVertex()->GetZ();
	if(isAOD) {
		AliAODHeader * header = static_cast<AliAODHeader*>(fInputEvent->GetHeader());
		centrality = header->GetCentrality();
	} else {
		centrality = static_cast<AliESDEvent*>(fInputEvent)->GetCentrality()->GetCentralityPercentile("V0M");
	}
	
	const Int_t centBin = GetBin(fAxisCent, centrality);
	const Int_t vertexBin = GetBin(fAxisZ, vertexz);
	
	if(DebugLevel () > 4) {
		cout << "centrality: " << centrality <<  " " << GetBin(fAxisCent, centrality) << endl;
		cout << "vertexz: " << vertexz <<  " " << GetBin(fAxisZ, vertexz) << endl;
	}
	
	if(centBin < 0 || vertexBin < 0) {
		//	AliError("bin out of range");
		return;
	}
	
	
	//TClonesArray * aodGammas = GetConversionGammas(isAOD);
	TClonesArray * aodGammas = fV0Reader->GetReconstructedGammas();
	if(!aodGammas) {
		AliError("no aod gammas found!");
		return;
	}
	
	TObjArray * ggammas = static_cast<TObjArray*>(fGammas.At(0));

	///Fill arrays of accepted gammas
	if(DebugLevel() > 1) printf("Number of conversion gammas %d \n", aodGammas->GetEntriesFast());
	for(Int_t ig = 0; ig < aodGammas->GetEntriesFast(); ig++) {
		AliAODConversionPhoton * photon = dynamic_cast<AliAODConversionPhoton*>(aodGammas->At(ig));
		if(!photon) continue;
		if(!fV0FilterPhoton || fV0FilterPhoton->PhotonIsSelected(photon, fInputEvent)) {
			ggammas->Add(photon);
		} else {
			for(Int_t igf = 0; igf < fV0Filters[1].GetEntriesFast(); igf++) {
				AliConversionPhotonCuts * gfilter = dynamic_cast<AliConversionPhotonCuts*>(fV0Filters[1].At(igf));
				if(gfilter && gfilter->PhotonIsSelected(photon, fInputEvent)) {
					static_cast<TObjArray*>(fGammas.At(igf+1))->Add(photon);
				}
			}
		}
	}

	Bool_t lowgmap[fV0Filters[0].GetEntriesFast()][ggammas->GetEntriesFast()];

	for(Int_t ig = 0; ig < ggammas->GetEntriesFast(); ig++ ) {
		AliAODConversionPhoton * gamma = static_cast<AliAODConversionPhoton*>(ggammas->At(ig));
		for(Int_t igf = 0; igf < fV0Filters[0].GetEntriesFast(); igf++) {
			
			AliConversionPhotonCuts * v0cuts = static_cast<AliConversionPhotonCuts*>(fV0Filters[0].At(igf));
			if(!(v0cuts->PhotonIsSelected(gamma, fInputEvent))) {
				lowgmap[igf][ig] = kTRUE;
			} else {
				lowgmap[igf][ig] = kFALSE;
			}
		}
	}
	

	if(DebugLevel() > 4) printf("Number of accepted gammas %d \n", ggammas->GetEntriesFast());
	hMEvents->Fill(vertexz, centrality);
	
	///create track array
	const Int_t ntrackfilters[2] = { fTrackFilters[0].GetEntriesFast(), fTrackFilters[1].GetEntriesFast()};
	
	TObjArray * ttracks = static_cast<TObjArray*>(fTracks.At(0));
	const Double_t aetalim[2] = { fAxisAssEta.GetXmin(), fAxisAssEta.GetXmax()};
	const Double_t aptlim[2] = { fAxiscPt.GetXmin(), fAxiscPt.GetXmax()};
	for(Int_t iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); iTrack++) {
		AliVTrack * track = static_cast<AliVTrack*>(fInputEvent->GetTrack(iTrack));
		if(track->Pt() < aptlim[0] || track->Pt() > aptlim[1]) continue;
		if(track->Eta() < aetalim[0] || track->Eta() > aetalim[1]) continue;
		if(fTrackFilter->IsSelected(track)) {
			hTrackPt->Fill(track->Pt(), centrality);
			ttracks->Add(track);
		} else {
			///upside cuts
			for(Int_t itf = 1; itf < ntrackfilters[1] + 1; itf++) {
				AliAnalysisCuts * trackCuts = static_cast<AliAnalysisCuts*>(fTrackFilters[1].At(itf -1));
				if(trackCuts->IsSelected(track)) {
					static_cast<TObjArray*>(fTracks.At(itf))->Add(track);
				}
			}
		}
	}


	Bool_t lowtrackmap[ntrackfilters[0]][ttracks->GetEntriesFast()];
	///Check lowside cuts
	for(Int_t iTrack = 0; iTrack < ttracks->GetEntriesFast(); iTrack++ ) {
		AliVTrack * track = static_cast<AliVTrack*>(ttracks->At(iTrack));
		for(Int_t itf = 0; itf < ntrackfilters[0]; itf++) {
			AliAnalysisCuts * trackCuts = static_cast<AliAnalysisCuts*>(fTrackFilters[0].At(itf));
			if(!trackCuts->IsSelected(track)) {
				lowtrackmap[itf][iTrack] = kTRUE;
			} else {
				lowtrackmap[itf][iTrack] = kFALSE;
			}
		}
	}

	hTrackCent->Fill(centrality, ttracks->GetEntriesFast());
	
	const Double_t etalim[2] = { fAxisTrigEta.GetXmin(), fAxisTrigEta.GetXmax()};
	if(DebugLevel() > 4) printf("Number of accepted gammas, tracks %d  %d \n", ggammas->GetEntriesFast(), ttracks->GetEntriesFast());
	
	//AliAnaConvCorrBase * gCorr = fPhotonCorr; //GetCorrObject(vertexBin, centBin, fPhotonCorr);
	//  AliAnaConvCorrPion * piCorr = fPionCorr; //static_cast<AliAnaConvCorrPion*>(GetCorrObject(vertexBin, centBin, fPionCorr));
	// if(!piCorr) {
	//   AliError("corr object missing");
	//   return;
	// }
	
	TObjArray pions;

	///corr values
	Double_t dphivalues[fAxesList.GetSize()];
	dphivalues[4] = centrality;
	dphivalues[5] = vertexz;

	///Trigger me counters
	Double_t trigValues[7];
	trigValues[2] = centrality;
	trigValues[3] = vertexz;

	///Mass histogram
	Double_t massval[5];
	massval[2] = centrality;
	
	///Set up track me counters and initialize
	const Int_t nbins = fAxistPt.GetNbins();
	Bool_t tmap[fAxistPt.GetNbins()];
	Bool_t lv0tmap[fAxistPt.GetNbins()][fV0Filters[0].GetEntriesFast()];
	Bool_t uv0tmap[fAxistPt.GetNbins()][fV0Filters[1].GetEntriesFast()];

	Bool_t lpitmap[fAxistPt.GetNbins()][fMesonFilters[0].GetEntriesFast()];
	Bool_t upitmap[fAxistPt.GetNbins()][fMesonFilters[1].GetEntriesFast()];
	
	for(Int_t igf = 0; igf < fV0Filters[0].GetEntriesFast(); igf++) {
		for(Int_t ptbin = 0; ptbin < fAxistPt.GetNbins(); ptbin++) {
			lv0tmap[ptbin][igf] = kFALSE;
		}
	}

	for(Int_t igf = 0; igf < fV0Filters[1].GetEntriesFast(); igf++) {
		for(Int_t ptbin = 0; ptbin < fAxistPt.GetNbins(); ptbin++) {
			uv0tmap[ptbin][igf] = kFALSE;
		}
	}

	for(Int_t igf = 0; igf < fMesonFilters[0].GetEntriesFast(); igf++) {
		for(Int_t ptbin = 0; ptbin < fAxistPt.GetNbins(); ptbin++) {
			lpitmap[ptbin][igf] = kFALSE;
		}
	}

	for(Int_t igf = 0; igf < fMesonFilters[1].GetEntriesFast(); igf++) {
		for(Int_t ptbin = 0; ptbin < fAxistPt.GetNbins(); ptbin++) {
			upitmap[ptbin][igf] = kFALSE;
		}
	}

	for(Int_t ptbin = 0; ptbin < nbins; ptbin++){
		tmap[ptbin] = kFALSE;
	}

	/////////////////////////////////////////


	for(Int_t igf1 = 0; igf1 < fV0Filters[1].GetEntriesFast() + 1; igf1++) {
		TObjArray * gamm1 = static_cast<TObjArray*>(fGammas.At(igf1));
		for(Int_t i1 = 0; i1 < gamm1->GetEntriesFast(); i1++) {
			AliAODConversionPhoton * ph1 = static_cast<AliAODConversionPhoton*>(gamm1->UncheckedAt(i1));
			Int_t tIDs[4] = {ph1->GetLabel(0), ph1->GetLabel(1), -1, -1};
			
			///Combine gamma into pions
			Int_t igmax = 0;
			for(Int_t igf2 = 0; igf2 <= igf1; igf2++) {
			TObjArray * gamm2 = NULL;
			if(igf2 == igf1) {
				gamm2 = gamm1;
				igmax = i1;
			} else {
				gamm2 = static_cast<TObjArray*>(fGammas.At(igf2));
				igmax = gamm2->GetEntriesFast();
			}

			for(Int_t i2 = 0; i2 < igmax; i2++) {
				AliAODConversionPhoton * ph2 = static_cast<AliAODConversionPhoton*>(gamm2->UncheckedAt(i2));
				
				if( ph2->GetTrackLabelPositive()==ph1->GetTrackLabelPositive() 
					|| ph2->GetTrackLabelNegative()==ph1->GetTrackLabelNegative()
					|| ph2->GetTrackLabelNegative()==ph1->GetTrackLabelPositive()
					|| ph2->GetTrackLabelPositive()==ph1->GetTrackLabelNegative()) {
					continue;
				}
				
				AliAODConversionMother * pion = new AliAODConversionMother(ph1, ph2);
				if(pion->Eta() < etalim[0] || pion->Eta() > etalim[1]) continue;
				pion->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());	  


				tIDs[2] = ph2->GetLabel(0);
				tIDs[3] = ph2->GetLabel(1);

				massval[0] = pion->M();
				massval[1] = pion->Pt();
				massval[3] = igf1;
				massval[4] = 0;

				dphivalues[2] = pion->Pt();
				dphivalues[6] = pion->M();
				dphivalues[8] = igf1;
				dphivalues[9] = 0;
				
				trigValues[0] = pion->Eta();
				trigValues[1] = pion->Pt();
				trigValues[4] = pion->M();
				trigValues[5] = igf1; 
				trigValues[6] = 0;
				

				///Check that particle is in histo phase space
				if(pion->Pt() > fAxistPt.GetBinLowEdge(1) && pion->Pt() < fAxistPt.GetXmax() &&
					pion->Eta() > etalim[0] && pion->Eta() < etalim[1] &&
					pion->M() > 0.04 && pion->M() < 0.4
					) {
						
						
						if(fMesonFilter->MesonIsSelected(pion, kTRUE)) {

							fMassSparse->Fill(massval);

							Bool_t lpimap[fMesonFilters[0].GetEntriesFast()];
							///See if it passes lowside cuts
							if(igf1 == 0 && igf2 == 0) {

								///Low side pion
								massval[3] = 0;
								for(Int_t ilpf = 0; ilpf < fMesonFilters[0].GetEntriesFast(); ilpf++) {
									if(!(static_cast<AliConversionMesonCuts*>(fMesonFilters[0].At(ilpf))->MesonIsSelected(pion, kTRUE))) {
										lpimap[ilpf] = kTRUE;
										massval[4] = -(ilpf + 1);
										fMassSparse->Fill(massval);
									} else {
										lpimap[ilpf] = kFALSE;
									}
									massval[4] = 0;
								}
								///Lowside v0
								for(Int_t iglf = 0; iglf < fV0Filters[0].GetEntriesFast(); iglf++) {
									if(lowgmap[iglf][i1] || lowgmap[iglf][i2]){
										massval[3] = -(iglf+1);
										fMassSparse->Fill(massval);
									}
								}
								massval[3] = 0;
							}  /// End lowside mass histo fillers
									
							
							///Check that particle is in histo phase space
							if(pion->Pt() > fAxistPt.GetBinLowEdge(1) && pion->Pt() < fAxistPt.GetXmax() &&
								pion->M() > fAxisPiM.GetXmin() && pion->M() < fAxisPiM.GetXmax() &&
								pion->Eta() > etalim[0] && pion->Eta() < etalim[1]) {

								const Int_t tbin = fAxistPt.FindFixBin(pion->Pt());
								///Fill standard triggers including upside v0 filters
								fTrigSparse->Fill(trigValues);
								if(igf1 == 0 && igf2 == 0) {
									// piCorr->FillTriggerCounters(pion);
									// piCorr->CorrelateWithTracks(pion, &tracks[0], tIDs, centrality, vertexz);								
									////Only mix events with pion in signal region
									hTrigPt->Fill(pion->Pt(), centrality, pion->M());
									if(pion->M() > 0.1 && pion->M() < 0.15) {
										hTrigPhi->Fill(pion->Phi());
										///Check trigger bin
										if (tbin > 0 && tbin < (nbins + 1)) {
											tmap[tbin-1] = kTRUE;
										}
										///Check if trigger also in low side (both gamma present in low side!)
										for(Int_t ilgf = 0; ilgf < fV0Filters[0].GetEntriesFast(); ilgf++) {
											if(!lowgmap[ilgf][i1] || !lowgmap[ilgf][i2]) {
												lv0tmap[tbin-1][ilgf] = kTRUE;
											}
										}
										///See if the lowside pion filter also passes this, if not
										for(Int_t ilpf = 0; ilpf < fMesonFilters[0].GetEntriesFast(); ilpf++) {
											if(!lpimap[ilpf]) {
												lpitmap[tbin-1][ilpf] = kTRUE;
											}
										}
									}
								} else {
									if(pion->M() > 0.1 && pion->M() < 0.15) {
										uv0tmap[tbin-1][igf1 - 1] = kTRUE;
									}
								}
										
									
								///Fill the triggers not selected in lowside filters only if passsing standard v0 filter
								if(igf1 == 0 && igf2 == 0) {	
									///Lowside v0 filters
									for(Int_t iglf = 0; iglf < fV0Filters[0].GetEntriesFast(); iglf++) {
										if(lowgmap[iglf][i1] || lowgmap[iglf][i2]){
										trigValues[5] = -(iglf+1);
										fTrigSparse->Fill(trigValues);
										}
									}
									////Low side pion filters
									trigValues[5] = 0;
									for(Int_t iplf = 0; iplf < fMesonFilters[0].GetEntriesFast(); iplf ++) {
										if(lpimap[iplf]) {
										trigValues[6] = -(iplf + 1);
										fTrigSparse->Fill(trigValues);
										}
									}
								} // ifg1 == 0 
							
								trigValues[5] = igf1; 
								trigValues[6] = 0;

								///////////////////////////////////////////////
								/// Correlate with tracks
								///////////////////////////////////////////////

								Int_t ntf = 1;
								if(igf1 == 0 && igf2 == 0) { 
									ntf = fTrackFilters[1].GetEntriesFast() + 1;
								}

								for(Int_t itf = 0; itf < ntf; itf++) {
									TObjArray * tracks = static_cast<TObjArray*>(fTracks.At(itf));
									for(int ij = 0; ij < tracks->GetEntriesFast(); ij++) {
										AliVTrack * track = static_cast<AliVTrack*>(tracks->At(ij));
										Int_t tid = track->GetID();

										if(tid == tIDs[0] || tid == tIDs[1] || tid == tIDs[2] || tid == tIDs[3])  {
											continue;
										}
										
										if(tid < 0) {
											if(-tid-1 == tIDs[0]+1 || -tid-1 == tIDs[1]+1 || -tid-1 == tIDs[2]+1 || -tid-1 == tIDs[3]+1)  {
												continue;
											}
										}
									
										dphivalues[0] = pion->Eta() - track->Eta();
										dphivalues[1] = GetDPhi(pion->Phi() - track->Phi());
										dphivalues[3] = track->Pt();
										dphivalues[7] = itf;
										dphivalues[8] = igf1;
										dphivalues[9] = 0;
										fCorrSparse->Fill(dphivalues, GetTrackCorrection(vertexz, track));
								
										if(itf == 0 && igf1 == 0 && igf2 == 0) {
											///Fill the low side track filters
											for(Int_t itlf = 0; itlf < fTrackFilters[0].GetEntriesFast(); itlf++) {
												if(lowtrackmap[itlf][ij]){
													dphivalues[7] = -(itlf+1);
													fCorrSparse->Fill(dphivalues, GetTrackCorrection(vertexz, track));
												}
											}
											///Fill the low side v0 filters
											dphivalues[7] = 0;
											for(Int_t iglf = 0; iglf < fV0Filters[0].GetEntriesFast(); iglf++) {
												if(lowgmap[iglf][i1] || lowgmap[iglf][i2]){
													dphivalues[8] = -(iglf+1);
													fCorrSparse->Fill(dphivalues, GetTrackCorrection(vertexz, track));
												}
											}

											///Fill the low side pi filter
											dphivalues[7] = 0;
											dphivalues[8] = 0;
											for(Int_t iplf = 0; iplf < fMesonFilters[0].GetEntriesFast(); iplf ++) {
												if(lpimap[iplf]) {
													dphivalues[9] = -(iplf + 1);
													fCorrSparse->Fill(dphivalues, GetTrackCorrection(vertexz, track));
												}
											}
										}  /// end non standard filters track corr
									
									} // end for tracks
								} // end trackfilters loop
							} //end check pion in histogram range to prevent overflow
						} else {
							/////////////////////////////
							//// Not passing standard meson cuts, check upside filters	
							////////////////////////////
							
							///Only check the pions from standard v0 filter 
							if(igf1 == 0 && igf2 == 0) {
								for(Int_t ipuf = 0; ipuf < fMesonFilters[1].GetEntriesFast(); ipuf++) {
									if(static_cast<AliConversionMesonCuts*>(fMesonFilters[1].At(ipuf))->MesonIsSelected(pion, kTRUE)) {
										///Fill invariant mass hist
										massval[4] = (ipuf + 1);
										fMassSparse->Fill(massval);
										///Check that particle is in histo phase space --- redundant!
										
										

										if(pion->Pt() > fAxistPt.GetBinLowEdge(1) && pion->Pt() < fAxistPt.GetXmax() &&
										pion->M() > fAxisPiM.GetXmin() && pion->M() < fAxisPiM.GetXmax() &&
										pion->Eta() > etalim[0] && pion->Eta() < etalim[1]) {
											
											
											////Only mix events with pion in signal region
											if(pion->M() > 0.1 && pion->M() < 0.15) {
												const Int_t tbin = fAxistPt.FindFixBin(pion->Pt());
												upitmap[tbin-1][ipuf] = kTRUE;
											}
											
											///Fill trigger counters
											trigValues[6] = (ipuf + 1);
											fTrigSparse->Fill(trigValues);
											
											///Correlate with standard tracks
											for(int ij = 0; ij < ttracks->GetEntriesFast(); ij++) {
												AliVTrack * track = static_cast<AliVTrack*>(ttracks->At(ij));
												Int_t tid = track->GetID();
												
												if(tid == tIDs[0] || tid == tIDs[1] || tid == tIDs[2] || tid == tIDs[3] ) {
													continue;
												}

												if(tid < 0) {
													if(-tid-1 == tIDs[0]+1 || -tid-1 == tIDs[1]+1 || -tid-1 == tIDs[2]+1 || -tid-1 == tIDs[3]+1)  {
														continue;
													}
												}
										
												dphivalues[0] = pion->Eta() - track->Eta();
												dphivalues[1] = GetDPhi(pion->Phi() - track->Phi());
												dphivalues[3] = track->Pt();
												dphivalues[7] = 0; // track filter
												dphivalues[8] = 0; // v0 filter
												dphivalues[9] = ipuf + 1; // pion filter
												fCorrSparse->Fill(dphivalues, GetTrackCorrection(vertexz, track));
											} /// end track corr
										}
									} // MesonIsSelected
								}
							}
						} /// end else ..  end upside meson filters
							/////////////////////////////////////////////
					} ///Etalim && pt limits
				} // i2 second gamma
			}
		} // i1 first gamma
	}
	//FillCounters(&pions, tracks, ntrackfilters, centrality, vertexz);
	
	///////Fill track counters after entire event has been passed through
	////
	
	Double_t trackValues[fTrackAxesList.GetSize()];
	trackValues[3] = centrality;
	trackValues[4] = vertexz;
	//trackValues[5] = particle->M(); remove !!!
	
	for(Int_t tbin = 0; tbin < fAxistPt.GetNbins(); tbin++) {
		trackValues[1] = fAxistPt.GetBinCenter(tbin+1);
		
		if(tmap[tbin]) {
			
			for(Int_t itf = 0; itf < fTrackFilters[1].GetEntriesFast() + 1; itf++) {
				TObjArray * tracks = static_cast<TObjArray*>(fTracks.At(itf));
				for(Int_t iTrack = 0; iTrack < tracks->GetEntriesFast(); iTrack++) {
					AliVTrack * track = static_cast<AliVTrack*>(tracks->At(iTrack));
					trackValues[0] = track->Eta();
					trackValues[2] = track->Pt();
					trackValues[5] = itf;
					trackValues[6] = 0;  ///v0 filter
					trackValues[7] = 0; ////Pi filter
					fTrackSparse->Fill(trackValues, GetTrackCorrection(vertexz, track));

					if(itf == 0) {
						for(Int_t itlf = 0; itlf < fTrackFilters[0].GetEntriesFast(); itlf++) {
							if(lowtrackmap[itlf][iTrack]) {
								trackValues[5] = -(itlf + 1);
								fTrackSparse->Fill(trackValues, GetTrackCorrection(vertexz, track) );
							}
						}
						trackValues[5] = 0;
						
						///Check lowside gamma
						for(Int_t iglf = 0; iglf < fV0Filters[0].GetEntriesFast(); iglf++) {
							if(!lv0tmap[tbin][iglf]) {
								trackValues[6] = -(iglf + 1);
								fTrackSparse->Fill(trackValues, GetTrackCorrection(vertexz, track));
							}
						}
						trackValues[6] = 0;

						////Check lowside pion
						for(Int_t iplf = 0; iplf < fMesonFilters[0].GetEntriesFast(); iplf++) {
							if(!lpitmap[tbin][iplf]) {
								trackValues[7] = -(iplf + 1);
								fTrackSparse->Fill(trackValues, GetTrackCorrection(vertexz, track));
							}
						}
					} // itf == 00
				}
			}
		} else {
			///If not in main, see if in upside filters
			///Do upside v0 filters
			for(Int_t iguf = 0; iguf < fV0Filters[1].GetEntriesFast(); iguf++) {
				if (uv0tmap[tbin][iguf] ) {

					//cout << "c vtx " <<  centrality << vertexz << endl;

					for(Int_t iTrack = 0; iTrack < ttracks->GetEntriesFast(); iTrack++) {
						AliVTrack * track = static_cast<AliVTrack*>(ttracks->At(iTrack));
						trackValues[0] = track->Eta();
						trackValues[2] = track->Pt();
						trackValues[5] = 0;
						trackValues[6] = iguf+1;  ///v0 filter
						trackValues[7] = 0; ////Pi filter
						fTrackSparse->Fill(trackValues, GetTrackCorrection(vertexz, track));
					}
				}
			}
				///Do upside pi filter
			for(Int_t ipuf = 0; ipuf < fMesonFilters[1].GetEntriesFast(); ipuf++) {
				if (upitmap[tbin][ipuf] ) {
					//cout << "c2 vtx2 " <<  centrality << vertexz << endl;
					for(Int_t iTrack = 0; iTrack < ttracks->GetEntriesFast(); iTrack++) {
						AliVTrack * track = static_cast<AliVTrack*>(ttracks->At(iTrack));
						trackValues[0] = track->Eta();
						trackValues[2] = track->Pt();
						trackValues[5] = 0;
						trackValues[6] = 0;  ///v0 filter
						trackValues[7] = ipuf+1; ////Pi filter
						fTrackSparse->Fill(trackValues, GetTrackCorrection(vertexz, track));
					}
				}
			}
		}
	}
	//////
	//////

	//piCorr->FillCounters(&pions, tracks, centrality, vertexz);
	
	PostData(1, fHistograms);

}

//_______________________________________________________________________________
// void AliAnalysisTaskdPhi::FillCounters(TObjArray * particles, TObjArray tracks[], Int_t ntrackfilters, Float_t cent, Float_t vtxz) {
  
//   return;
  
  
//  //Fill ME Counters
//   const Int_t nbins = fAxistPt.GetNbins();
//   Bool_t tmap[nbins];
//   for(Int_t ptbin = 0; ptbin < nbins; ptbin++){
//       tmap[ptbin] = kFALSE;
//     }
//   }


//   Double_t trackValues[fTrackAxesList.GetSize()];
//   trackValues[3] = cent;
//   trackValues[4] = vtxz;

//   for(Int_t ip = 0; ip < particles->GetEntriesFast(); ip++){
//     AliAODConversionParticle * particle = static_cast<AliAODConversionParticle*>(particles->At(ip));

//     Int_t tbin = fAxistPt.FindFixBin(particle->Pt());
//     if (tbin > 0 && tbin < nbins + 1) {
//       if(tmap[tbin - 1] == kTRUE) {
//  	continue;
//       } else {
//  	tmap[tbin -1 ] = kTRUE;

//  	trackValues[5] = particle->M();
//  	trackValues[1] = particle->Pt();
	
//  	for(Int_t itf = 0; itf < ntrackfilters; itf++) {
//  	  if(fkTrackAxis) trackValues[6] = itf;
//  	  for(int ij = 0; ij < tracks[itf].GetEntriesFast(); ij++) {
//  	    AliVTrack * track = static_cast<AliVTrack*>(tracks->UncheckedAt(ij));
//  	    trackValues[0] = track->Eta();
//  	    trackValues[2] = track->Pt();
//  	    fTrackSparse->Fill(trackValues);	
//  	  }
//  	}
//       }
//     }
//   }
// }




// //________________________________________________________________
// void AliAnalysisTaskdPhi::CorrelateWithTracks(AliAODConversionParticle * particle, TObjArray tracks[], Int_t ntrackfilters, Bool_t ** lowtrackmap, Int_t nltf, Int_t const tIDs[4], Double_t dphivalues[]) {
//   //Correlate particle with tracks
//   ///Correlate with tracks

// }

//________________________________________________________________________
void AliAnalysisTaskdPhi::Terminate(Option_t *) {
	
	// Draw result to the screen
	// Called once at the end of the query
}

//________________________________________________________________________
TClonesArray * AliAnalysisTaskdPhi::GetConversionGammas(Bool_t isAOD) {
	if(isAOD) {

		TClonesArray * gammas = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(fDeltaAODBranchName.Data()));
		if(gammas) {
			return gammas;
		}
		
		FindDeltaAODBranchName(fInputEvent);
		gammas = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(fDeltaAODBranchName.Data()));
		return gammas;
		
	} else {
		TClonesArray * gammas = dynamic_cast<TClonesArray*>(GetInputData(1));
		return gammas;
	}
  
}

//________________________________________________________________________
void AliAnalysisTaskdPhi::FindDeltaAODBranchName(AliVEvent * event){
	///Find aod branch
	TList *list=event->GetList();
	for(Int_t ii=0;ii<list->GetEntries();ii++){
		TString name((list->At(ii))->GetName());
		if(name.BeginsWith("GammaConv")&&name.EndsWith("gamma")){
			fDeltaAODBranchName=name;
			AliDebug(AliLog::kDebug + 5, Form("Set DeltaAOD BranchName to: %s",fDeltaAODBranchName.Data()));
			return;
		}
	}
}
  

//________________________________________________________________________
Double_t AliAnalysisTaskdPhi::GetTrackCorrection(Double_t vtxz, AliVTrack * track) {
	////Get track correction from map
	Int_t coord[4] = {-1, -1, -1, -1};
	if(fCorrectionMap) {
		Double_t values[4] = { vtxz, track->Pt(), track->Eta(), track->Phi() };
		Double_t correction = fCorrectionMap->GetBinContent(fCorrectionMap->GetBin(values, kFALSE), coord);
		if (fCorrectionMap->IsInRange(coord)) {
			return correction;
		} 
	}
	return 1.0;
}


