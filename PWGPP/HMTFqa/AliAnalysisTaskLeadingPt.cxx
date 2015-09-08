/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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
 * provided "as is" without express or implied warranty.                  * **************************************************************************/

//Authors: Antonio Ortiz Velasquez, antonio.ortiz@nucleares.unam.mx
//modified by: Sergio Iga , sergio.arturo.iga.buitron@cern.ch
//ptleading distributions, vs nch analysis

#include "AliAnalysisTaskLeadingPt.h"

// ROOT includes
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TProfile.h>
#include <TParticle.h>
#include <TFile.h>

// AliRoot includes
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>

#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>

#include <TTreeStream.h>

#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>

#include "AliCentrality.h" 
#include <AliESDv0.h>
#include <AliKFVertex.h>
#include <AliAODVertex.h>

#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliAODMCHeader.h> 
#include <AliAODHeader.h>
#include <AliESDUtils.h>

#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"

// STL includes
#include <iostream>
using namespace std;



ClassImp(AliAnalysisTaskLeadingPt)
	//_____________________________________________________________________________
	//AliAnalysisTaskLeadingPt::AliAnalysisTaskQAHighPtDeDx(const char *name):
AliAnalysisTaskLeadingPt::AliAnalysisTaskLeadingPt():
	AliAnalysisTaskSE(),
	fESD(0x0),
	fAOD(0x0),
	fTrackFilterGolden(0x0),
	fTrackFilterTPC(0x0),
	fAnalysisType("ESD"),
	ftrigBit(0x0),
	fRandom(0x0),
	fPileUpRej(kFALSE),
	fPileUpRejMV(kFALSE),
	fAliAnalysisUtils(0),
	fnContributors(5),
	fVtxCut(10.0),  
	fEtaCut(0.9),  
	fTriggeredEventMB(-999),
	fVtxStatus(-999),
	fZvtx(-999),
	fRun(-999),
	fEventId(-999),
	fPtCutLeading(0),
	fListOfObjects(0),
	fEvents(0x0),
	fVtxBeforeCuts(0x0), 
	fVtxAfterCuts(0x0),
	fn1(0x0),
	ftrackVsClusters(0),
	ftrackVsClusters_SPD(0),
	ftrackVsClusters_MV(0),
	ftrackVsClusters_SPDMV(0)


{
	//default constructor, seven multiplicity classes
	for(Int_t i=0;i<7;++i){
		//W/o PileUp Rejection
		fnchL[i]=0;
		fptLL[i]=0;
		fptL[i]=0;
		fnchH[i]=0;
		fptLH[i]=0;
		fptH[i]=0;

		//SPD PileUp Rejection
		fnchL_SPD[i]=0;
		fptLL_SPD[i]=0;
		fptL_SPD[i]=0;
		fnchH_SPD[i]=0;
		fptLH_SPD[i]=0;
		fptH_SPD[i]=0;

		//MV PileUp Rejection
		fnchL_MV[i]=0;
		fptLL_MV[i]=0;
		fptL_MV[i]=0;
		fnchH_MV[i]=0;
		fptLH_MV[i]=0;
		fptH_MV[i]=0;

		//SPD && MV PileUp Rejection
		fnchL_SPDMV[i]=0;
		fptLL_SPDMV[i]=0;
		fptL_SPDMV[i]=0;
		fnchH_SPDMV[i]=0;
		fptLH_SPDMV[i]=0;
		fptH_SPDMV[i]=0;

	}



}


AliAnalysisTaskLeadingPt::AliAnalysisTaskLeadingPt(const char *name):
	AliAnalysisTaskSE(name),
	fESD(0x0),
	fAOD(0x0),
	fTrackFilterGolden(0x0),
	fTrackFilterTPC(0x0),
	fAnalysisType("ESD"),
	ftrigBit(0x0),
	fRandom(0x0),
	fPileUpRej(kFALSE),
	fPileUpRejMV(kFALSE),
	fAliAnalysisUtils(0),
	fnContributors(5),
	fVtxCut(10.0),  
	fEtaCut(0.9),  
	fTriggeredEventMB(-999),
	fVtxStatus(-999),
	fZvtx(-999),
	fRun(-999),
	fEventId(-999),
	fPtCutLeading(0),
	fListOfObjects(0),
	fEvents(0x0), 
	fVtxBeforeCuts(0x0), 
	fVtxAfterCuts(0x0),
	fn1(0x0),
	ftrackVsClusters(0),
	ftrackVsClusters_SPD(0),
	ftrackVsClusters_MV(0),
	ftrackVsClusters_SPDMV(0)

{
	// Default constructor (should not be used)
	for(Int_t i=0;i<7;++i){
		//W/o PileUp Rejection
		fnchL[i]=0;
		fptLL[i]=0;
		fptL[i]=0;
		fnchH[i]=0;
		fptLH[i]=0;
		fptH[i]=0;

		//SPD PileUp Rejection
		fnchL_SPD[i]=0;
		fptLL_SPD[i]=0;
		fptL_SPD[i]=0;
		fnchH_SPD[i]=0;
		fptLH_SPD[i]=0;
		fptH_SPD[i]=0;

		//MV PileUp Rejection
		fnchL_MV[i]=0;
		fptLL_MV[i]=0;
		fptL_MV[i]=0;
		fnchH_MV[i]=0;
		fptLH_MV[i]=0;
		fptH_MV[i]=0;

		//SPD && MV PileUp Rejection
		fnchL_SPDMV[i]=0;
		fptLL_SPDMV[i]=0;
		fptL_SPDMV[i]=0;
		fnchH_SPDMV[i]=0;
		fptLH_SPDMV[i]=0;
		fptH_SPDMV[i]=0;

	}

	DefineOutput(1, TList::Class());//esto es nuevo

}




AliAnalysisTaskLeadingPt::~AliAnalysisTaskLeadingPt() {
	//
	// Destructor
	//

	if (fListOfObjects) {
		delete fListOfObjects;
		fListOfObjects = 0x0;
	}

}

//______________________________________________________________________________
void AliAnalysisTaskLeadingPt::UserCreateOutputObjects()
{ 
	// This method is called once per worker node
	// Here we define the output: histograms and debug tree if requested 
	// We also create the random generator here so it might get different seeds...

	fRandom = new TRandom(0); // 0 means random seed

	//OpenFile(1);
	fListOfObjects = new TList();
	fListOfObjects->SetOwner();

	//
	// Histograms
	//  
	fEvents = new TH1I("fEvents","Number of analyzed events; Events; Counts", 3, 0, 3);
	fListOfObjects->Add(fEvents);

	fVtxBeforeCuts = new TH1F("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
	fListOfObjects->Add(fVtxBeforeCuts);

	fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
	fListOfObjects->Add(fVtxAfterCuts);

	fn1=new TH1F("fn1","fn1",11,-1,10);
	fListOfObjects->Add(fn1);

	const Int_t nPtBins = 68;
	Double_t xBins[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
		0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
		1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
		2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
		4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
		11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
		26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0 };

	const Char_t * legendsNch[7] = {"all #it{z}", "0 < #it{z} #leq 1", "1 < #it{z} #leq 2", 
		"2 < #it{z} #leq 3", "3 < #it{z} #leq 4", "4 < #it{z} #leq 5", "#it{z} > 5" };
	
	//W/o PileUp Rejection

	for( Int_t i_z = 0; i_z < 7; ++i_z ){

		//soft below the ptleading cut
		fnchL[i_z] = new TH1F( Form( "fnchL%d",i_z ),
				Form("%s multiplicity distribution #it{p}_{T}^{#it{l}} < %1.1f; #it{N}_{ch} (combined |#eta|<0.8); entries", legendsNch[i_z], fPtCutLeading ), 
				126, -6, 120 );
		fListOfObjects->Add( fnchL[i_z] );

		fptLL[i_z] = new TH1F( Form( "fptLL%d",i_z ),
				Form("%s #it{p}_{T}^{leading} distribution; #it{p}_{T}^{#it{l}} (GeV/#it{c}); entries", legendsNch[i_z] ), 
				nPtBins, xBins );
		fListOfObjects->Add( fptLL[i_z] );

		fptL[i_z] = new TH1F( Form( "fptL%d",i_z ),
				Form("%s #it{p}_{T} distribution #it{p}_{T}^{#it{l}} < %1.1f; #it{p}_{T} (GeV/#it{c}); entries", legendsNch[i_z], fPtCutLeading ),
				nPtBins, xBins );
		fListOfObjects->Add( fptL[i_z] );


		//hard above the ptleading cut
		fnchH[i_z] = new TH1F( Form( "fnchH%d",i_z ),
				Form("%s multiplicity distribution #it{p}_{T}^{#it{l}} > %1.1f; #it{N}_{ch} (combined |#eta|<0.8); entries", legendsNch[i_z], fPtCutLeading ),
				126, -6, 120 );
		fListOfObjects->Add( fnchH[i_z] );

		fptLH[i_z] = new TH1F( Form( "fptLH%d",i_z ),
				Form("%s #it{p}_{T}^{leading} distribution; #it{p}_{T}^{#it{l}} (GeV/#it{c}); entries", legendsNch[i_z] ),
				nPtBins, xBins );
		fListOfObjects->Add( fptLH[i_z] );

		fptH[i_z] = new TH1F( Form( "fptH%d",i_z ),
				Form("%s #it{p}_{T} distribution #it{p}_{T}^{#it{l}} > %1.1f; #it{p}_{T} (GeV/#it{c}); entries", legendsNch[i_z], fPtCutLeading ),
				nPtBins, xBins );
		fListOfObjects->Add( fptH[i_z] );


	}

	//SPD PileUp Rejection

	for( Int_t i_z = 0; i_z < 7; ++i_z ){

		//soft below the ptleading cut
		fnchL_SPD[i_z] = new TH1F( Form( "fnchL_SPD%d",i_z ),
				Form("%s multiplicity distribution #it{p}_{T}^{#it{l}} < %1.1f; #it{N}_{ch} (combined |#eta|<0.8); entries", legendsNch[i_z], fPtCutLeading ), 
				126, -6, 120 );
		fListOfObjects->Add( fnchL_SPD[i_z] );

		fptLL_SPD[i_z] = new TH1F( Form( "fptLL_SPD%d",i_z ),
				Form("%s #it{p}_{T}^{leading} distribution; #it{p}_{T}^{#it{l}} (GeV/#it{c}); entries", legendsNch[i_z] ), 
				nPtBins, xBins );
		fListOfObjects->Add( fptLL_SPD[i_z] );

		fptL_SPD[i_z] = new TH1F( Form( "fptL_SPD%d",i_z ),
				Form("%s #it{p}_{T} distribution #it{p}_{T}^{#it{l}} < %1.1f; #it{p}_{T} (GeV/#it{c}); entries", legendsNch[i_z], fPtCutLeading ),
				nPtBins, xBins );
		fListOfObjects->Add( fptL_SPD[i_z] );


		//hard above the ptleading cut
		fnchH_SPD[i_z] = new TH1F( Form( "fnchH_SPD%d",i_z ),
				Form("%s multiplicity distribution #it{p}_{T}^{#it{l}} > %1.1f; #it{N}_{ch} (combined |#eta|<0.8); entries", legendsNch[i_z], fPtCutLeading ),
				126, -6, 120 );
		fListOfObjects->Add( fnchH_SPD[i_z] );

		fptLH_SPD[i_z] = new TH1F( Form( "fptLH_SPD%d",i_z ),
				Form("%s #it{p}_{T}^{leading} distribution; #it{p}_{T}^{#it{l}} (GeV/#it{c}); entries", legendsNch[i_z] ),
				nPtBins, xBins );
		fListOfObjects->Add( fptLH_SPD[i_z] );

		fptH_SPD[i_z] = new TH1F( Form( "fptH_SPD%d",i_z ),
				Form("%s #it{p}_{T} distribution #it{p}_{T}^{#it{l}} > %1.1f; #it{p}_{T} (GeV/#it{c}); entries", legendsNch[i_z], fPtCutLeading ),
				nPtBins, xBins );
		fListOfObjects->Add( fptH_SPD[i_z] );


	}

	// MV PileUpRejection

	for( Int_t i_z = 0; i_z < 7; ++i_z ){

		//soft below the ptleading cut
		fnchL_MV[i_z] = new TH1F( Form( "fnchL_MV%d",i_z ),
				Form("%s multiplicity distribution #it{p}_{T}^{#it{l}} < %1.1f; #it{N}_{ch} (combined |#eta|<0.8); entries", legendsNch[i_z], fPtCutLeading ), 
				126, -6, 120 );
		fListOfObjects->Add( fnchL_MV[i_z] );

		fptLL_MV[i_z] = new TH1F( Form( "fptLL_MV%d",i_z ),
				Form("%s #it{p}_{T}^{leading} distribution; #it{p}_{T}^{#it{l}} (GeV/#it{c}); entries", legendsNch[i_z] ), 
				nPtBins, xBins );
		fListOfObjects->Add( fptLL_MV[i_z] );

		fptL_MV[i_z] = new TH1F( Form( "fptL_MV%d",i_z ),
				Form("%s #it{p}_{T} distribution #it{p}_{T}^{#it{l}} < %1.1f; #it{p}_{T} (GeV/#it{c}); entries", legendsNch[i_z], fPtCutLeading ),
				nPtBins, xBins );
		fListOfObjects->Add( fptL_MV[i_z] );


		//hard above the ptleading cut
		fnchH_MV[i_z] = new TH1F( Form( "fnchH_MV%d",i_z ),
				Form("%s multiplicity distribution #it{p}_{T}^{#it{l}} > %1.1f; #it{N}_{ch} (combined |#eta|<0.8); entries", legendsNch[i_z], fPtCutLeading ),
				126, -6, 120 );
		fListOfObjects->Add( fnchH_MV[i_z] );

		fptLH_MV[i_z] = new TH1F( Form( "fptLH_MV%d",i_z ),
				Form("%s #it{p}_{T}^{leading} distribution; #it{p}_{T}^{#it{l}} (GeV/#it{c}); entries", legendsNch[i_z] ),
				nPtBins, xBins );
		fListOfObjects->Add( fptLH_MV[i_z] );

		fptH_MV[i_z] = new TH1F( Form( "fptH_MV%d",i_z ),
				Form("%s #it{p}_{T} distribution #it{p}_{T}^{#it{l}} > %1.1f; #it{p}_{T} (GeV/#it{c}); entries", legendsNch[i_z], fPtCutLeading ),
				nPtBins, xBins );
		fListOfObjects->Add( fptH_MV[i_z] );


	}

	// SPD & MV PileUpRejection

	for( Int_t i_z = 0; i_z < 7; ++i_z ){

		//soft below the ptleading cut
		fnchL_SPDMV[i_z] = new TH1F( Form( "fnchL_SPDMV%d",i_z ),
				Form("%s multiplicity distribution #it{p}_{T}^{#it{l}} < %1.1f; #it{N}_{ch} (combined |#eta|<0.8); entries", legendsNch[i_z], fPtCutLeading ), 
				126, -6, 120 );
		fListOfObjects->Add( fnchL_SPDMV[i_z] );

		fptLL_SPDMV[i_z] = new TH1F( Form( "fptLL_SPDMV%d",i_z ),
				Form("%s #it{p}_{T}^{leading} distribution; #it{p}_{T}^{#it{l}} (GeV/#it{c}); entries", legendsNch[i_z] ), 
				nPtBins, xBins );
		fListOfObjects->Add( fptLL_SPDMV[i_z] );

		fptL_SPDMV[i_z] = new TH1F( Form( "fptL_SPDMV%d",i_z ),
				Form("%s #it{p}_{T} distribution #it{p}_{T}^{#it{l}} < %1.1f; #it{p}_{T} (GeV/#it{c}); entries", legendsNch[i_z], fPtCutLeading ),
				nPtBins, xBins );
		fListOfObjects->Add( fptL_SPDMV[i_z] );


		//hard above the ptleading cut
		fnchH_SPDMV[i_z] = new TH1F( Form( "fnchH_SPDMV%d",i_z ),
				Form("%s multiplicity distribution #it{p}_{T}^{#it{l}} > %1.1f; #it{N}_{ch} (combined |#eta|<0.8); entries", legendsNch[i_z], fPtCutLeading ),
				126, -6, 120 );
		fListOfObjects->Add( fnchH_SPDMV[i_z] );

		fptLH_SPDMV[i_z] = new TH1F( Form( "fptLH_SPDMV%d",i_z ),
				Form("%s #it{p}_{T}^{leading} distribution; #it{p}_{T}^{#it{l}} (GeV/#it{c}); entries", legendsNch[i_z] ),
				nPtBins, xBins );
		fListOfObjects->Add( fptLH_SPDMV[i_z] );

		fptH_SPDMV[i_z] = new TH1F( Form( "fptH_SPDMV%d",i_z ),
				Form("%s #it{p}_{T} distribution #it{p}_{T}^{#it{l}} > %1.1f; #it{p}_{T} (GeV/#it{c}); entries", legendsNch[i_z], fPtCutLeading ),
				nPtBins, xBins );
		fListOfObjects->Add( fptH_SPDMV[i_z] );


	}

	ftrackVsClusters = new TH2F("ftrackVsClusters","nTracklets vs SPD clusters",180,0,180,700,0,700);
	fListOfObjects->Add( ftrackVsClusters );
	ftrackVsClusters_SPD = new TH2F("ftrackVsClusters_SPD","nTracklets vs SPD clusters",180,0,180,700,0,700);
	fListOfObjects->Add(ftrackVsClusters_SPD );
	ftrackVsClusters_MV = new TH2F("ftrackVsClusters_MV","nTracklets vs SPD clusters",180,0,180,700,0,700);
	fListOfObjects->Add( ftrackVsClusters_MV );
	ftrackVsClusters_SPDMV = new TH2F("ftrackVsClusters_SPDMV","nTracklets vs SPD clusters",180,0,180,700,0,700);
	fListOfObjects->Add( ftrackVsClusters_SPDMV );

	// Post output data.
	PostData(1, fListOfObjects);

	if(!fAliAnalysisUtils)
	{
		fAliAnalysisUtils = new AliAnalysisUtils();
		SetPileUpMvSettings(fAliAnalysisUtils);
	}

}

//______________________________________________________________________________
void AliAnalysisTaskLeadingPt::UserExec(Option_t *) 
{
	// Main loop

	//
	// First we make sure that we have valid input(s)!
	//



	AliVEvent *event = InputEvent();
	if (!event) {
		Error("UserExec", "Could not retrieve event");
		return;
	}



	if (fAnalysisType == "ESD"){
		fESD = dynamic_cast<AliESDEvent*>(event);
		if(!fESD){
			Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}    
	} else {
		fAOD = dynamic_cast<AliAODEvent*>(event);
		if(!fAOD){
			Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}    
	}





	// Get trigger decision
	fTriggeredEventMB = 0; //init

	fn1->Fill(0);

	((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler() ))->SetNeedField();

	if( ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
			->IsEventSelected() & ftrigBit ){
		fTriggeredEventMB = 1;  //event triggered as minimum bias
	}


	// real data that are not triggered we skip
	if(!fTriggeredEventMB)
		return; 

	//AliESDUtils::RefitESDVertexTracks(fESD); // para LHC11c_p1 solamente

	fn1->Fill(1);



	if (fAnalysisType == "ESD"){

		const AliESDVertex *vtxESD = fESD->GetPrimaryVertexTracks();
		if(vtxESD->GetNContributors()<1) {
			// SPD vertex
			vtxESD = fESD->GetPrimaryVertexSPD();
			/* quality checks on SPD-vertex */
			TString vertexType = vtxESD->GetTitle();
			if (vertexType.Contains("vertexer: Z") && (vtxESD->GetDispersion() > 0.04 || vtxESD->GetZRes() > 0.25))  
				fZvtx  = -1599; //vertex = 0x0; //
			else if (vtxESD->GetNContributors()<1) 
				fZvtx  = -999; //vertex = 0x0; //
			else
				fZvtx = vtxESD->GetZ();
		}  
		else
			fZvtx = vtxESD->GetZ();

	}
	else // AOD
		fZvtx = GetVertex(fAOD);

	fVtxBeforeCuts->Fill(fZvtx);

	//cut on the z position of vertex
	if (TMath::Abs(fZvtx) > fVtxCut) {	
		return;
	}
	fn1->Fill(2);






	Float_t centrality = -10;

	if (fAnalysisType == "ESD"){
		AnalyzeESD(fESD);
	} else {
		AnalyzeAOD(fAOD);
	}
	

	fVtxAfterCuts->Fill(fZvtx);


	// Post output data.
	PostData(1, fListOfObjects);
}

//________________________________________________________________________
void AliAnalysisTaskLeadingPt::AnalyzeESD(AliESDEvent* esdEvent)
{
	fRun  = esdEvent->GetRunNumber();
	fEventId = 0;
	if(esdEvent->GetHeader())
		fEventId = GetEventIdAsLong(esdEvent->GetHeader());

	Int_t trackmult08=0;
	trackmult08=AliESDtrackCuts::GetReferenceMultiplicity(esdEvent, AliESDtrackCuts::kTrackletsITSTPC, 0.8);

    const AliVMultiplicity* mult = esdEvent->GetMultiplicity();
    if (!mult) { cout<<"No multiplicity object"<<endl; return; }
    int ntracklet   = mult->GetNumberOfTracklets();
    int spdClusters = esdEvent->GetNumberOfITSClusters(0) + esdEvent->GetNumberOfITSClusters(1);
	
	fn1->Fill(3);

	// accepted event
	fEvents->Fill(0);

	const Int_t nESDTracks = esdEvent->GetNumberOfTracks();

	//First loop to get ptleading
	Double_t ptL = 0;
	Int_t iTL = -1;
	for(Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = esdEvent->GetTrack(iT);

		//only golden track cuts
		UInt_t selectDebug = 0;
		if (fTrackFilterGolden) {
			selectDebug = fTrackFilterGolden->IsSelected(esdTrack);
			if (!selectDebug) {
				//cout<<"this is not a golden track"<<endl;
				continue;
			}
		}

		Double_t eta  = esdTrack->Eta();
		Double_t pt   = esdTrack->Pt();

		if( TMath::Abs(eta) > fEtaCut )
			continue;
		//ptmin cut, only > 150 Mev/c
		if( pt < 0.15 )
			continue;
		//extracting the ptleading
		if( pt > ptL ){
			ptL = pt;
			iTL = iT;
		}

	}//end of track loop

	Bool_t isHard = kFALSE;
	if( ptL >=  fPtCutLeading )
		isHard = kTRUE;

	Float_t z_08 = trackmult08 / 10.19;
	//extracting the multiplicity bin
	Int_t bin_z_08 = 0;
	if( z_08 > 0.0 && z_08 <= 1.0 )
		bin_z_08 = 1;
	if( z_08 > 1.0 && z_08 <= 2.0 )
		bin_z_08 = 2;
	if( z_08 > 2.0 && z_08 <= 3.0 )
		bin_z_08 = 3;
	if( z_08 > 3.0 && z_08 <= 4.0 )
		bin_z_08 = 4;
	if( z_08 > 4.0 && z_08 <= 5.0 )
		bin_z_08 = 5;
	if( z_08 > 5.0 )
		bin_z_08 = 6;

	//Filling multiplicity distribution
	isHard ? fnchH[0]->Fill(trackmult08) : fnchL[0]->Fill(trackmult08);
	isHard ? fptLH[0]->Fill(ptL) : fptLL[0]->Fill(ptL);
	if( bin_z_08 > 0 ){
		isHard  ? fnchH[bin_z_08]->Fill(trackmult08) : fnchL[bin_z_08]->Fill(trackmult08);
		isHard  ? fptLH[bin_z_08]->Fill(ptL) : fptLL[bin_z_08]->Fill(ptL);
	}


	for(Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = esdEvent->GetTrack(iT);

		//only golden track cuts
		UInt_t selectDebug = 0;
		if (fTrackFilterGolden) {
			selectDebug = fTrackFilterGolden->IsSelected(esdTrack);
			if (!selectDebug) {
				//cout<<"this is not a golden track"<<endl;
				continue;
			}
		}

		Double_t eta  = esdTrack->Eta();
		Double_t pt   = esdTrack->Pt();

		if( TMath::Abs(eta) > fEtaCut )
			continue;
		//ptmin cut, only > 150 Mev/c
		if( pt < 0.15 )
			continue;

		isHard ? fptH[0]->Fill(pt) : fptL[0]->Fill(pt);		
		if( bin_z_08 > 0 ){
			isHard ? fptH[bin_z_08]->Fill(pt) : fptL[bin_z_08]->Fill(pt);
		}

	}//end of track loop

	if(bin_z_08 == 6)
		ftrackVsClusters->Fill(ntracklet,spdClusters);


	Bool_t isPileup = kTRUE;
	if(fPileUpRej)
		isPileup = esdEvent->IsPileupFromSPD();

	if(!isPileup){ //if SPD pileup rejected
		fn1->Fill(4);

		// accepted event
		fEvents->Fill(0);

		const Int_t nESDTracks = esdEvent->GetNumberOfTracks();

		//First loop to get ptleading
		Double_t ptL = 0;
		Int_t iTL = -1;
		for(Int_t iT = 0; iT < nESDTracks; iT++) {

			AliESDtrack* esdTrack = esdEvent->GetTrack(iT);

			//only golden track cuts
			UInt_t selectDebug = 0;
			if (fTrackFilterGolden) {
				selectDebug = fTrackFilterGolden->IsSelected(esdTrack);
				if (!selectDebug) {
					//cout<<"this is not a golden track"<<endl;
					continue;
				}
			}

			Double_t eta  = esdTrack->Eta();
			Double_t pt   = esdTrack->Pt();

			if( TMath::Abs(eta) > fEtaCut )
				continue;
			//ptmin cut, only > 150 Mev/c
			if( pt < 0.15 )
				continue;
			//extracting the ptleading
			if( pt > ptL ){
				ptL = pt;
				iTL = iT;
			}

		}//end of track loop

		Bool_t isHard = kFALSE;
		if( ptL >=  fPtCutLeading )
			isHard = kTRUE;

		Float_t z_08 = trackmult08 / 10.19;
		//extracting the multiplicity bin
		Int_t bin_z_08 = 0;
		if( z_08 > 0.0 && z_08 <= 1.0 )
			bin_z_08 = 1;
		if( z_08 > 1.0 && z_08 <= 2.0 )
			bin_z_08 = 2;
		if( z_08 > 2.0 && z_08 <= 3.0 )
			bin_z_08 = 3;
		if( z_08 > 3.0 && z_08 <= 4.0 )
			bin_z_08 = 4;
		if( z_08 > 4.0 && z_08 <= 5.0 )
			bin_z_08 = 5;
		if( z_08 > 5.0 )
			bin_z_08 = 6;

		//Filling multiplicity distribution
		isHard ? fnchH_SPD[0]->Fill(trackmult08) : fnchL_SPD[0]->Fill(trackmult08);
		isHard ? fptLH_SPD[0]->Fill(ptL)         : fptLL_SPD[0]->Fill(ptL);
		if( bin_z_08 > 0 ){
			isHard ? fnchH_SPD[bin_z_08]->Fill(trackmult08) : fnchL_SPD[bin_z_08]->Fill(trackmult08);
			isHard ? fptLH_SPD[bin_z_08]->Fill(ptL)         : fptLL_SPD[bin_z_08]->Fill(ptL);
		}

		for(Int_t iT = 0; iT < nESDTracks; iT++) {

			AliESDtrack* esdTrack = esdEvent->GetTrack(iT);

			//only golden track cuts
			UInt_t selectDebug = 0;
			if (fTrackFilterGolden) {
				selectDebug = fTrackFilterGolden->IsSelected(esdTrack);
				if (!selectDebug) {
					//cout<<"this is not a golden track"<<endl;
					continue;
				}
			}

			Double_t eta  = esdTrack->Eta();
			Double_t pt   = esdTrack->Pt();

			if( TMath::Abs(eta) > fEtaCut )
				continue;
			//ptmin cut, only > 150 Mev/c
			if( pt < 0.15 )
				continue;

			isHard ? fptH_SPD[0]->Fill(pt) : fptL_SPD[0]->Fill(pt);		
			if( bin_z_08 > 0 ){
				isHard ? fptH_SPD[bin_z_08]->Fill(pt) : fptL_SPD[bin_z_08]->Fill(pt);
			}

		}//end of track loop
		if(bin_z_08 == 6)
			ftrackVsClusters_SPD->Fill(ntracklet,spdClusters);

	}
	
	Bool_t isPileUpMV = kTRUE;
	if(fPileUpRejMV)	
		isPileUpMV = fAliAnalysisUtils->IsPileUpMV(esdEvent);


	if(!isPileUpMV){
		fn1->Fill(5);

		// accepted event
		fEvents->Fill(0);

		const Int_t nESDTracks = esdEvent->GetNumberOfTracks();

		//First loop to get ptleading
		Double_t ptL = 0;
		Int_t iTL = -1;
		for(Int_t iT = 0; iT < nESDTracks; iT++) {

			AliESDtrack* esdTrack = esdEvent->GetTrack(iT);

			//only golden track cuts
			UInt_t selectDebug = 0;
			if (fTrackFilterGolden) {
				selectDebug = fTrackFilterGolden->IsSelected(esdTrack);
				if (!selectDebug) {
					//cout<<"this is not a golden track"<<endl;
					continue;
				}
			}

			Double_t eta  = esdTrack->Eta();
			Double_t pt   = esdTrack->Pt();

			if( TMath::Abs(eta) > fEtaCut )
				continue;
			//ptmin cut, only > 150 Mev/c
			if( pt < 0.15 )
				continue;
			//extracting the ptleading
			if( pt > ptL ){
				ptL = pt;
				iTL = iT;
			}

		}//end of track loop

		Bool_t isHard = kFALSE;
		if( ptL >=  fPtCutLeading )
			isHard = kTRUE;

		Float_t z_08 = trackmult08 / 10.19;
		//extracting the multiplicity bin
		Int_t bin_z_08 = 0;
		if( z_08 > 0.0 && z_08 <= 1.0 )
			bin_z_08 = 1;
		if( z_08 > 1.0 && z_08 <= 2.0 )
			bin_z_08 = 2;
		if( z_08 > 2.0 && z_08 <= 3.0 )
			bin_z_08 = 3;
		if( z_08 > 3.0 && z_08 <= 4.0 )
			bin_z_08 = 4;
		if( z_08 > 4.0 && z_08 <= 5.0 )
			bin_z_08 = 5;
		if( z_08 > 5.0 )
			bin_z_08 = 6;

		//Filling multiplicity distribution
		isHard ? fnchH_MV[0]->Fill(trackmult08) : fnchL_MV[0]->Fill(trackmult08);
		isHard ? fptLH_MV[0]->Fill(ptL)         : fptLL_MV[0]->Fill(ptL);
		if( bin_z_08 > 0 ){
			isHard ? fnchH_MV[bin_z_08]->Fill(trackmult08) : fnchL_MV[bin_z_08]->Fill(trackmult08);
			isHard ? fptLH_MV[bin_z_08]->Fill(ptL)         : fptLL_MV[bin_z_08]->Fill(ptL);
		}

		for(Int_t iT = 0; iT < nESDTracks; iT++) {

			AliESDtrack* esdTrack = esdEvent->GetTrack(iT);

			//only golden track cuts
			UInt_t selectDebug = 0;
			if (fTrackFilterGolden) {
				selectDebug = fTrackFilterGolden->IsSelected(esdTrack);
				if (!selectDebug) {
					//cout<<"this is not a golden track"<<endl;
					continue;
				}
			}

			Double_t eta  = esdTrack->Eta();
			Double_t pt   = esdTrack->Pt();

			if( TMath::Abs(eta) > fEtaCut )
				continue;
			//ptmin cut, only > 150 Mev/c
			if( pt < 0.15 )
				continue;

			isHard ? fptH_MV[0]->Fill(pt) : fptL_MV[0]->Fill(pt);		
			if( bin_z_08 > 0 ){
				isHard ? fptH_MV[bin_z_08]->Fill(pt) : fptL_MV[bin_z_08]->Fill(pt);
			}

		}//end of track loop
		if(bin_z_08 == 6)
			ftrackVsClusters_MV->Fill(ntracklet,spdClusters);

	}

	if(!isPileUpMV && !isPileup){
		fn1->Fill(6);	

		// accepted event
		fEvents->Fill(0);

		const Int_t nESDTracks = esdEvent->GetNumberOfTracks();

		//First loop to get ptleading
		Double_t ptL = 0;
		Int_t iTL = -1;
		for(Int_t iT = 0; iT < nESDTracks; iT++) {

			AliESDtrack* esdTrack = esdEvent->GetTrack(iT);

			//only golden track cuts
			UInt_t selectDebug = 0;
			if (fTrackFilterGolden) {
				selectDebug = fTrackFilterGolden->IsSelected(esdTrack);
				if (!selectDebug) {
					//cout<<"this is not a golden track"<<endl;
					continue;
				}
			}

			Double_t eta  = esdTrack->Eta();
			Double_t pt   = esdTrack->Pt();

			if( TMath::Abs(eta) > fEtaCut )
				continue;
			//ptmin cut, only > 150 Mev/c
			if( pt < 0.15 )
				continue;
			//extracting the ptleading
			if( pt > ptL ){
				ptL = pt;
				iTL = iT;
			}

		}//end of track loop

		Bool_t isHard = kFALSE;
		if( ptL >=  fPtCutLeading )
			isHard = kTRUE;

		Float_t z_08 = trackmult08 / 10.19;
		//extracting the multiplicity bin
		Int_t bin_z_08 = 0;
		if( z_08 > 0.0 && z_08 <= 1.0 )
			bin_z_08 = 1;
		if( z_08 > 1.0 && z_08 <= 2.0 )
			bin_z_08 = 2;
		if( z_08 > 2.0 && z_08 <= 3.0 )
			bin_z_08 = 3;
		if( z_08 > 3.0 && z_08 <= 4.0 )
			bin_z_08 = 4;
		if( z_08 > 4.0 && z_08 <= 5.0 )
			bin_z_08 = 5;
		if( z_08 > 5.0 )
			bin_z_08 = 6;

		//Filling multiplicity distribution
		isHard ? fnchH_SPDMV[0]->Fill(trackmult08) : fnchL_SPDMV[0]->Fill(trackmult08);
		isHard ? fptLH_SPDMV[0]->Fill(ptL)         : fptLL_SPDMV[0]->Fill(ptL);
		if( bin_z_08 > 0 ){
			isHard ? fnchH_SPDMV[bin_z_08]->Fill(trackmult08) : fnchL_SPDMV[bin_z_08]->Fill(trackmult08);
			isHard ? fptLH_SPDMV[bin_z_08]->Fill(ptL)         : fptLL_SPDMV[bin_z_08]->Fill(ptL);
		}

		for(Int_t iT = 0; iT < nESDTracks; iT++) {

			AliESDtrack* esdTrack = esdEvent->GetTrack(iT);

			//only golden track cuts
			UInt_t selectDebug = 0;
			if (fTrackFilterGolden) {
				selectDebug = fTrackFilterGolden->IsSelected(esdTrack);
				if (!selectDebug) {
					//cout<<"this is not a golden track"<<endl;
					continue;
				}
			}

			Double_t eta  = esdTrack->Eta();
			Double_t pt   = esdTrack->Pt();

			if( TMath::Abs(eta) > fEtaCut )
				continue;
			//ptmin cut, only > 150 Mev/c
			if( pt < 0.15 )
				continue;

			isHard ? fptH_SPDMV[0]->Fill(pt) : fptL_SPDMV[0]->Fill(pt);		
			if( bin_z_08 > 0 ){
				isHard ? fptH_SPDMV[bin_z_08]->Fill(pt) : fptL_SPDMV[bin_z_08]->Fill(pt);
			}

		}//end of track loop

		if(bin_z_08 == 6)
			ftrackVsClusters_SPDMV->Fill(ntracklet,spdClusters);

	}

}

//________________________________________________________________________
void AliAnalysisTaskLeadingPt::AnalyzeAOD(AliAODEvent* aodEvent)
{
	fRun  = aodEvent->GetRunNumber();
	fEventId = 0;
	if(aodEvent->GetHeader())
		fEventId = GetEventIdAsLong(aodEvent->GetHeader());

	Int_t     trackmult08 = 0; // no pt cuts
	trackmult08 =  (dynamic_cast<AliAODHeader*>(aodEvent->GetHeader()))->GetRefMultiplicityComb08();

	//UInt_t    time      = 0; // Missing AOD info? aodEvent->GetTimeStamp();

	//Int_t     trackmult = 0; // no pt cuts
	//Int_t     nadded    = 0;

	Bool_t isPileup = aodEvent->IsPileupFromSPD();
	if(fPileUpRej)
		if(isPileup)
			return;
	fn1->Fill(4);



	if(fTriggeredEventMB) {// Only MC case can we have not triggered events

		// accepted event
		fEvents->Fill(0);

		Int_t nAODTracks = aodEvent->GetNumberOfTracks();
		//First loop to get ptleading
		Double_t ptL = 0;

		for(Int_t iT = 0; iT < nAODTracks; iT++) {

			//AliAODTrack* aodTrack = AODevent->GetTrack(iT);
			AliVTrack   *trackTmp = (AliVTrack *)aodEvent->GetTrack(iT);
			AliAODTrack * aodTrack  = dynamic_cast<AliAODTrack *>(trackTmp);

			if (fTrackFilterGolden) {     
				// "Global track RAA analysis QM2011 + Chi2ITS<36"; bit 1024
				if(!aodTrack->TestFilterBit(1024))
					continue;
			}

			Double_t eta  = aodTrack->Eta();
			Double_t pt   = aodTrack->Pt();

			if( TMath::Abs(eta) > fEtaCut )
				continue;
			//ptmin cut, only > 150 Mev/c
			if( pt < 0.15 )
				continue;
			//extracting the ptleading
			if( pt > ptL )
				ptL = pt;

		}//end of track loop

		Bool_t isHard = kFALSE;
		if( ptL >=  fPtCutLeading )
			isHard = kTRUE;

		Float_t z_08 = trackmult08 / 10.19;
		//extracting the multiplicity bin
		Int_t bin_z_08 = 0;
		if( z_08 > 0.0 && z_08 <= 1.0 )
			bin_z_08 = 1;
		if( z_08 > 1.0 && z_08 <= 2.0 )
			bin_z_08 = 2;
		if( z_08 > 2.0 && z_08 <= 3.0 )
			bin_z_08 = 3;
		if( z_08 > 3.0 && z_08 <= 4.0 )
			bin_z_08 = 4;
		if( z_08 > 4.0 && z_08 <= 5.0 )
			bin_z_08 = 5;
		if( z_08 > 5.0 )
			bin_z_08 = 6;

		//Filling multiplicity distribution
		isHard ? fnchH[0]->Fill(trackmult08) : fnchL[0]->Fill(trackmult08);
		isHard ? fptLH[0]->Fill(ptL) : fptLL[0]->Fill(ptL);
		if( bin_z_08 > 0 ){
			isHard ? fnchH[bin_z_08]->Fill(trackmult08) : fnchL[bin_z_08]->Fill(trackmult08);
			isHard ? fptLH[bin_z_08]->Fill(ptL) : fptLL[bin_z_08]->Fill(ptL);
		}

		for(Int_t iT = 0; iT < nAODTracks; iT++) {

			//AliAODTrack* aodTrack = AODevent->GetTrack(iT);
			AliVTrack   *trackTmp = (AliVTrack *)aodEvent->GetTrack(iT);
			AliAODTrack * aodTrack  = dynamic_cast<AliAODTrack *>(trackTmp);

			if (fTrackFilterGolden) {     
				// "Global track RAA analysis QM2011 + Chi2ITS<36"; bit 1024
				if(!aodTrack->TestFilterBit(1024))
					continue;
			}

			Double_t eta  = aodTrack->Eta();
			Double_t pt   = aodTrack->Pt();

			if( TMath::Abs(eta) > fEtaCut )
				continue;
			//ptmin cut, only > 150 Mev/c
			if( pt < 0.15 )
				continue;

			isHard ? fptH[0]->Fill(pt) : fptL[0]->Fill(pt);
			if( bin_z_08 > 0 ){
				isHard ? fptH[bin_z_08]->Fill(pt) : fptL[bin_z_08]->Fill(pt);
			}

		}//end of track loop


	} // end if triggered

}

//_____________________________________________________________________________
Float_t AliAnalysisTaskLeadingPt::GetVertex(const AliVEvent* event) const
{
	Float_t zvtx = -999;

	const AliVVertex* primaryVertex = event->GetPrimaryVertex(); 

	if(primaryVertex->GetNContributors()>0)
		zvtx = primaryVertex->GetZ();

	return zvtx;
}
//_____________________________________________________________________________
ULong64_t AliAnalysisTaskLeadingPt::GetEventIdAsLong(AliVHeader* header) const
{
	// To have a unique id for each event in a run!
	// Modified from AliRawReader.h
	return ((ULong64_t)header->GetBunchCrossNumber()+
			(ULong64_t)header->GetOrbitNumber()*3564+
			(ULong64_t)header->GetPeriodNumber()*16777215*3564);
}

void AliAnalysisTaskLeadingPt::SetPileUpMvSettings(AliAnalysisUtils* fAliAnalysisUtils)
{
	fAliAnalysisUtils->SetMinPlpContribMV(fnContributors);
	fAliAnalysisUtils->SetMaxPlpChi2MV(5.0);
	fAliAnalysisUtils->SetMinWDistMV(15.0);
	fAliAnalysisUtils->SetCheckPlpFromDifferentBCMV(kFALSE);
}

