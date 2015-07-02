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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// This is a simple task designed to test AliPPVsMultUtils functionality.
// --- david.dobrigkeit.chinellato@cern.ch
// Also has been added analysis for dN/dEta and Multiplicity + Multiplicity correlations with V0 estimators+ V0 sectors
// ---- hector.bello.martinez@cern.ch
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


#include <Riostream.h>
#include "TList.h"
#include <TChain.h>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
//#define LOG_NO_INFO
//#define LOG_NO_DEBUG
//#define LOG_NO_WARNING
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include <AliHeader.h>
#include "AliCentrality.h"

#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"
#include "AliAnalysisMultCorrTaskQA.h"

#include <AliESDVertex.h>
#include <AliMultiplicity.h>

#include <TTree.h>
#include <TDirectory.h>
#include <TBits.h>

#include "AliESDInputHandler.h"
#include <AliAnalysisFilter.h>

using std::cout;
using std::endl;

ClassImp(AliAnalysisMultCorrTaskQA)

AliAnalysisMultCorrTaskQA::AliAnalysisMultCorrTaskQA()
    : AliAnalysisTaskSE(), 
    fListHist(0), 
    fCuts(0),
    fTrackFilterSKE(0x0),
    fHistEventCounter(0), 
    fHistRefMult08(0), 
    fHistRefMult05(0), 
    fHistV0M(0), 
    fHistV0A(0), 
    fHistV0C(0), 
    fdNdeta(0), 
    fcorrRef05Ref08(0), 
    fcorrV0ARef08(0), 
    fcorrV0CRef08(0), 
    fcorrV0MRef08(0), 
    fHistV0Aamp(0), 
    fHistV0Camp(0), 
    fHistV0Mamp(0), 
    fcorrV0AampRef08(0), 
    fcorrV0CampRef08(0), 
    fcorrV0MampRef08(0),
    fModulesV0(0), 
    fPPVsMultUtils(0)
{

}

AliAnalysisMultCorrTaskQA::AliAnalysisMultCorrTaskQA(const char *name)
    : AliAnalysisTaskSE(name), 
    fListHist(0), 
    fCuts(0), 
    fTrackFilterSKE(0x0),
    fHistEventCounter(0), 
    fHistRefMult08(0), 
    fHistRefMult05(0), 
    fHistV0M(0), 
    fHistV0A(0), 
    fHistV0C(0), 
    fdNdeta(0), 
    fcorrRef05Ref08(0), 
    fcorrV0ARef08(0), 
    fcorrV0CRef08(0),
    fcorrV0MRef08(0), 
    fHistV0Aamp(0), 
    fHistV0Camp(0), 
    fHistV0Mamp(0), 
    fcorrV0AampRef08(0),
    fcorrV0CampRef08(0), 
    fcorrV0MampRef08(0),
    fModulesV0(0), 
    fPPVsMultUtils(0)
{
    DefineOutput(1, TList::Class()); // List of Histograms
}


AliAnalysisMultCorrTaskQA::~AliAnalysisMultCorrTaskQA()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------

    if (fListHist) {
        delete fListHist;
        fListHist = 0x0;
    }
    if (fPPVsMultUtils) {
        delete fPPVsMultUtils;
        fPPVsMultUtils = 0x0;
    }
}

//________________________________________________________________________
void AliAnalysisMultCorrTaskQA::UserCreateOutputObjects()
{
    //Helper
    if(! fPPVsMultUtils ) {
        fPPVsMultUtils = new AliPPVsMultUtils();
    }
    //... cuts .....
    fCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    
    inputHandler->SetNeedField();

    // Create histograms
    fListHist = new TList();
    fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner
    if(! fHistEventCounter ){
        fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",7,0,7);
        fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fHistEventCounter->GetXaxis()->SetBinLabel(2, "Selected by Analysis");
        fHistEventCounter->GetXaxis()->SetBinLabel(3, "PileupSPDInMultBins");
        fHistEventCounter->GetXaxis()->SetBinLabel(4, "NotINEL>0");
        fHistEventCounter->GetXaxis()->SetBinLabel(5, "Not inVertexcut");
        fHistEventCounter->GetXaxis()->SetBinLabel(6, "Inconsistent SPD&TrackVertx");
        fHistEventCounter->GetXaxis()->SetBinLabel(7, "Not MinimumBias");
        fListHist->Add(fHistEventCounter);
    }
    //Histogram Output: Event-by-Event
    if(! fHistRefMult08 ) {
        fHistRefMult08 = new TH1D( "fHistRefMult08", "Multiplicity |#eta| < 0.8 ;Ref. Mult. |#eta| < 0.8;Count",200,0,200);
        fListHist->Add(fHistRefMult08);
    }
    if(! fHistRefMult05 ) {
        fHistRefMult05 = new TH1D( "fHistRefMult05", "Multiplicity |#eta| < 0.5 ;Ref. Mult. |#eta| < 0.5;Count",200,0,200);
        fListHist->Add(fHistRefMult05);
    }
    if(! fHistV0M ) {
        fHistV0M = new TH1D( "fHistV0M", "Multiplicity V0M;V0M Percentile;Count",100,0,100);
        fListHist->Add(fHistV0M);
    }
    if(! fHistV0A ) {
        fHistV0A = new TH1D( "fHistV0A", "Multiplicity V0A;V0A Percentile;Count",100,0,100);
        fListHist->Add(fHistV0A);
    }
    if(! fHistV0C ) {
        fHistV0C = new TH1D( "fHistV0C", "Multiplicity V0C;V0C Percentile;Count",100,0,100);
        fListHist->Add(fHistV0C);
    }
    if(! fHistV0Aamp ) {
        fHistV0Aamp = new TH1D( "fHistV0Aamp", "Multiplicity V0A;V0A Amplitud;Count",700,0,700);
        fListHist->Add(fHistV0Aamp);
    }
    if(! fHistV0Camp ) {
        fHistV0Camp = new TH1D( "fHistV0Camp", "Multiplicity V0C;V0C Amplitud;Count",700,0,700);
        fListHist->Add(fHistV0Camp);
    }
    if(! fHistV0Mamp ) {
        fHistV0Mamp = new TH1D( "fHistV0Mamp", "Multiplicity V0M;V0M Amplitud;Count",700,0,700);
        fListHist->Add(fHistV0Mamp);
    }
    if(!fdNdeta){
        fdNdeta= new TH1D( "fdNdeta","dN/d#eta (|#eta| < 0.8);#eta (rads);Count",100,-1,1);
        fListHist->Add(fdNdeta);
    }

    if(!fcorrRef05Ref08){
        fcorrRef05Ref08= new TH2D( "fcorrRef05Ref08"," Multiplicity Correlation (|#eta| < 0.5 and |#eta| < 0.8); Nch (|#eta| < 0.5); Nch (|#eta| < 0.8)",200,0,200,200,0,200);
        fListHist->Add(fcorrRef05Ref08);
    }

    if(!fcorrV0ARef08){
	    fcorrV0ARef08= new TH2D( "fcorrV0ARef08"," Multiplicity Correlation (V0A percentile and |#eta| < 0.8); V0A Percentile; Nch (|#eta| < 0.8)",100,0,100,200,0,200);
	    fListHist->Add(fcorrV0ARef08);
    }
    if(!fcorrV0CRef08){
	    fcorrV0CRef08= new TH2D( "fcorrV0CRef08"," Multiplicity Correlation (V0C percentile and |#eta| < 0.8); V0C Percentile; Nch (|#eta| < 0.8)",100,0,100,200,0,200);
	    fListHist->Add(fcorrV0CRef08);
    }
    if(!fcorrV0AampRef08){
	    fcorrV0AampRef08= new TH2D( "fcorrV0AampRef08"," Multiplicity Correlation (V0A and |#eta| < 0.8); V0A Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
	    fListHist->Add(fcorrV0AampRef08);
    }
    if(!fcorrV0CampRef08){
	    fcorrV0CampRef08= new TH2D( "fcorrV0CampRef08"," Multiplicity Correlation (V0C and |#eta| < 0.8); V0C Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
	    fListHist->Add(fcorrV0CampRef08);
    }
    if(!fcorrV0MRef08){
	    fcorrV0MRef08= new TH2D( "fcorrV0MRef08"," Multiplicity Correlation (V0M and |#eta| < 0.8); V0M Percentile; Nch (|#eta| < 0.8)",100,0,100,200,0,200);
	    fListHist->Add(fcorrV0MRef08);
    }
    if(!fcorrV0MampRef08){
	    fcorrV0MampRef08= new TH2D( "fcorrV0MampRef08"," Multiplicity Correlation (V0M and |#eta| < 0.8); V0M Amplitud; Nch (|#eta| < 0.8)",700,0,700,200,0,200);
	    fListHist->Add(fcorrV0MampRef08);
    }
    if(!fModulesV0){
	    fModulesV0= new TH2D( "fModulesV0"," Multiplicity vs cell; V0 Sector; counts ",64,0,64,1000,0,1000);
	    fListHist->Add(fModulesV0);
    }


    PostData(1, fListHist);
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisMultCorrTaskQA::UserExec(Option_t *)
{
	// Main loop
	// Called for each event
	AliESDEvent *lESDevent = 0x0;

	lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
	if (!lESDevent) {
		AliWarning("ERROR: lESDevent not available \n");
		return;
	}
	//------------------------------------------------
	// Selection Investigation with AliPPVsMultUtils
	//------------------------------------------------
	fHistEventCounter->Fill(0.5);

	if( AliPPVsMultUtils::IsNotPileupSPDInMultBins(lESDevent) == kFALSE){fHistEventCounter->Fill(2.5);}
	if( AliPPVsMultUtils::IsINELgtZERO( lESDevent ) == kFALSE){fHistEventCounter->Fill(3.5);}
	if( AliPPVsMultUtils::IsAcceptedVertexPosition( lESDevent ) == kFALSE){fHistEventCounter->Fill(4.5);}
	if( AliPPVsMultUtils::HasNoInconsistentSPDandTrackVertices( lESDevent ) == kFALSE){fHistEventCounter->Fill(5.5);}
	if( AliPPVsMultUtils::IsMinimumBias( lESDevent ) == kFALSE){fHistEventCounter->Fill(6.5);}


	if( !AliPPVsMultUtils::IsEventSelected(lESDevent) ) {
		PostData(1, fListHist);// Event isn't selected, post output data, done here
		return;
	}


	fHistEventCounter->Fill(1.5);
	//This should NOT have any underflow: 
	// ---> all error codes should have been removed by IsEventSelected already 

	//Reference Multiplicity
	fHistRefMult08->Fill(AliPPVsMultUtils::GetStandardReferenceMultiplicity( lESDevent ) );
	fHistRefMult05->Fill(AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.5) );
	//V0 Percentiles
	fHistV0M->Fill( fPPVsMultUtils->GetMultiplicityPercentile( lESDevent, "V0M" ) );
	fHistV0A->Fill( fPPVsMultUtils->GetMultiplicityPercentile( lESDevent, "V0A" ) );
	fHistV0C->Fill( fPPVsMultUtils->GetMultiplicityPercentile( lESDevent, "V0C" ) );
	LoopESD(lESDevent);
	AliVVZERO* esdV0 = lESDevent->GetVZEROData();

	Float_t multV0A= esdV0->GetMTotV0A();
	Float_t multV0C= esdV0->GetMTotV0C();
	Float_t multV0M= multV0A+multV0C;
	
        fHistV0Aamp->Fill( esdV0->GetMTotV0A() );
	fHistV0Camp->Fill( esdV0->GetMTotV0C() );


	fHistV0Mamp->Fill( multV0M );

	//Correlations
	fcorrRef05Ref08->Fill(AliESDtrackCuts::GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 0.5), AliPPVsMultUtils::GetStandardReferenceMultiplicity( lESDevent ) );
	fcorrV0ARef08->Fill( fPPVsMultUtils->GetMultiplicityPercentile( lESDevent, "V0A" ), AliPPVsMultUtils::GetStandardReferenceMultiplicity( lESDevent ));
	fcorrV0CRef08->Fill( fPPVsMultUtils->GetMultiplicityPercentile( lESDevent, "V0C" ), AliPPVsMultUtils::GetStandardReferenceMultiplicity( lESDevent ));
	fcorrV0AampRef08->Fill( esdV0->GetMTotV0A(), AliPPVsMultUtils::GetStandardReferenceMultiplicity( lESDevent ));
	fcorrV0CampRef08->Fill( esdV0->GetMTotV0C(), AliPPVsMultUtils::GetStandardReferenceMultiplicity( lESDevent ));
	fcorrV0MampRef08->Fill( esdV0->GetMTotV0A()+esdV0->GetMTotV0C(), AliPPVsMultUtils::GetStandardReferenceMultiplicity( lESDevent ));
	fcorrV0MRef08->Fill( fPPVsMultUtils->GetMultiplicityPercentile( lESDevent, "V0M" ), AliPPVsMultUtils::GetStandardReferenceMultiplicity( lESDevent ));



	for(Int_t ich =0; ich < 64; ich++){
		fModulesV0->Fill(ich,esdV0->GetMultiplicity(ich));
	}
	// Post 
	PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisMultCorrTaskQA::Terminate(Option_t *)
{
	// Draw result to the screen
	/*
	   TList *cRetrievedList = 0x0;
	   cRetrievedList = (TList*)GetOutputData(1);
	   if(!cRetrievedList) {
	   Printf("ERROR - AliAnalysisMultCorrTaskQA : ouput data container list not available\n");
	   return;
	   }

	   fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
	   if (!fHistEventCounter) {
	   Printf("ERROR - AliAnalysisMultCorrTaskQA : fHistEventCounter not available");
	   return;
	   }

	   TCanvas *canCheck = new TCanvas("AliAnalysisMultCorrTaskQA","Control Histo",10,10,510,510);
	   canCheck->cd(1)->SetLogy();

	   fHistEventCounter->SetMarkerStyle(22);
	   fHistEventCounter->DrawCopy("E");*/
}

//----------------------------------------------
void AliAnalysisMultCorrTaskQA::LoopESD( AliESDEvent *lESDevent)
{

	Int_t ntracks = lESDevent->GetNumberOfTracks();
	if(ntracks==0){
		return;
	}
	Int_t nesdtracks=0;
	TObjArray* acceptedtracks = new TObjArray();
	for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
		AliESDtrack *track = (AliESDtrack *)lESDevent->GetTrack(iTracks);

		if (!track) {
			Error("UserExec", "Could not receive track %d", iTracks);
			continue;
		}

		if(!fCuts->AcceptTrack(track))continue;

		nesdtracks++;
		fdNdeta->Fill(track->Eta());
		acceptedtracks->Add(track);

	} 
	delete acceptedtracks;

}


// - david.dobrigkeit.chinellato@cern.ch
// --hector.bello.martinez@cern.ch
