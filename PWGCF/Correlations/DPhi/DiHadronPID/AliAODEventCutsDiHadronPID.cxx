// -------------------------------------------------------------------------
// Object managing event cuts, and holding QA histograms.
// -------------------------------------------------------------------------
// Author: Misha Veldhoen (misha.veldhoen@cern.ch)

#include <iostream>
#include "TList.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNamed.h"
#include "TIterator.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODEventCutsDiHadronPID.h"

using namespace std;

ClassImp(AliAODEventCutsDiHadronPID);

// -------------------------------------------------------------------------
AliAODEventCutsDiHadronPID::AliAODEventCutsDiHadronPID():
	TNamed(),
	fIsPbPb(kTRUE),
	fIsMC(kFALSE),
	fTrigger(AliVEvent::kMB),
	fMinCentrality(5.),
	fMaxCentrality(0.),
	fCentralityEstimator("V0M"),
	fMaxVertexZ(10.),
	fMinRefMult(0),
	fTestTrigger(kFALSE),
	fTestCentrality(kFALSE),
	fTestVertexZ(kFALSE),
	fTestMinRefMult(kFALSE),
	fSelectedEventQAHistos(0x0),
	fAllEventQAHistos(0x0),
	fHistTrigger(0x0),
	fHistRefMultiplicity(0x0),
	fHistCentrality(0x0),
	fHistCentralityQuality(0x0),
	fHistVertexZ(0x0),
	fDebug(0)

{

	// 
	// Default Constructor
	//

	cout<<"AliAODEventCutsDiHadronPID Default Constructor Called."<<endl;
	if (fDebug > 1) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

}

// -------------------------------------------------------------------------
AliAODEventCutsDiHadronPID::AliAODEventCutsDiHadronPID(const char* name):
	TNamed(name,"AOD Event Cuts"),
	fIsPbPb(kTRUE),
	fIsMC(kFALSE),	
	fTrigger(AliVEvent::kMB),
	fMinCentrality(5.),
	fMaxCentrality(0.),
	fCentralityEstimator("V0M"),
	fMaxVertexZ(10.),
	fMinRefMult(0),	
	fTestTrigger(kFALSE),
	fTestCentrality(kFALSE),
	fTestVertexZ(kFALSE),
	fTestMinRefMult(kFALSE),	
	fSelectedEventQAHistos(0x0),
	fAllEventQAHistos(0x0),
	fHistTrigger(0x0),
	fHistRefMultiplicity(0x0),
	fHistCentrality(0x0),
	fHistCentralityQuality(0x0),
	fHistVertexZ(0x0),
	fDebug(0)
	
{

	//
	// Named Constructor
	//

	cout<<"AliAODEventCutsDiHadronPID Named Constructor Called."<<endl;
	if (fDebug > 1) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

}

// -------------------------------------------------------------------------
AliAODEventCutsDiHadronPID::~AliAODEventCutsDiHadronPID() {

	//
	// Destructor
	//

	cout<<"AliAODEventCutsDiHadronPID Destructor Called."<<endl;
	if (fDebug > 1) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	if (fSelectedEventQAHistos) delete fSelectedEventQAHistos;
	fSelectedEventQAHistos = 0x0;
	if (fAllEventQAHistos) delete fAllEventQAHistos;
	fAllEventQAHistos = 0x0;

}

// -------------------------------------------------------------------------
Long64_t AliAODEventCutsDiHadronPID::Merge(TCollection* list) {

	//
	// Merger. 
	// 

	cout<<"AliAODEventCutsDiHadronPID Merger Called."<<endl;
	if (fDebug > 1) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	if (!list) return 0;
	if (list->IsEmpty()) return 1;

	if (!fSelectedEventQAHistos||!fAllEventQAHistos) {
		cout<<"AliAODEventCutsDiHadronPID::Merge() - Warning, current object's histograms are missing... Generating."<<endl;
		CreateHistos();
	}

	TIterator* iter = list->MakeIterator();
	TObject* obj;

	// List of collections
	TList collection_fSelectedEventQAHistos;
	TList collection_fAllEventQAHistos;

	Int_t count = 0;

  	while ((obj = iter->Next())) {
    	AliAODEventCutsDiHadronPID* entry = dynamic_cast<AliAODEventCutsDiHadronPID*> (obj);
    	if (entry == 0) continue;

    	// Check if the object to be merged really has the same name! (FIXME!)

    	// Getting the lists from obj.
    	TList* list_fSelectedEventQAHistos = entry->GetListOfSelectedEventQAHistos();
    	TList* list_fAllEventQAHistos = entry->GetListOfAllEventQAHistos();

    	// Adding the retrieved lists to the collection.
    	if (list_fSelectedEventQAHistos) collection_fSelectedEventQAHistos.Add(list_fSelectedEventQAHistos);
    	if (list_fAllEventQAHistos) collection_fAllEventQAHistos.Add(list_fAllEventQAHistos);

    	count++;
    }

    // Merging. Note that we require the original list to exist.
    //  * Assume that if the collection happens to be empty, then nothing will happen.
    //  * All other variables are taken from the original object.
    if (fSelectedEventQAHistos) fSelectedEventQAHistos->Merge(&collection_fSelectedEventQAHistos);
    if (fAllEventQAHistos) fAllEventQAHistos->Merge(&collection_fAllEventQAHistos);

    delete iter;

	return count+1;

}

// -------------------------------------------------------------------------
void AliAODEventCutsDiHadronPID::CreateHistos() {

	cout<<"AliAODEventCutsDiHadronPID - Creating histograms"<<endl;
	if (fDebug > 1) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Create list of Event related QA histograms (selected events).
	fSelectedEventQAHistos = new TList();
	fSelectedEventQAHistos->SetName("SelectedEventQAHistos");
	fSelectedEventQAHistos->SetOwner(kTRUE);

	// The same, but for all events.
	fAllEventQAHistos = new TList();
	fAllEventQAHistos->SetName("AllEventQAHistos");
	fAllEventQAHistos->SetOwner(kTRUE);

	// Creating arrays of pointers to the QA histos.
	fHistTrigger = new TH1F*[2];
	fHistRefMultiplicity = new TH1F*[2];
	fHistCentrality = new TH1F*[2];
	fHistCentralityQuality = new TH1F*[2];
	fHistVertexZ = new TH1F*[2];

	const char* HistType[2] = {"Selected","All"};

	for (Int_t iHistType = 0; iHistType < 2; iHistType++) {

		// Trigger Histogram.
		fHistTrigger[iHistType] = new TH1F(Form("fHistTrigger%s",HistType[iHistType]),"Trigger;;Count",2,-0.5,1.5);
		fHistTrigger[iHistType]->GetXaxis()->SetBinLabel(1,"kMB");
		fHistTrigger[iHistType]->GetXaxis()->SetBinLabel(2,"Other");
		if (iHistType == 0) fSelectedEventQAHistos->Add(fHistTrigger[iHistType]);
		else fAllEventQAHistos->Add(fHistTrigger[iHistType]);

		// Ref Multiplicity Histogram.
		fHistRefMultiplicity[iHistType] = new TH1F(Form("fHistRefMultiplicity%s",HistType[iHistType]),"Reference Multiplicity;N_{tracks};Count",100,0.,10000.);
		if (iHistType == 0) fSelectedEventQAHistos->Add(fHistRefMultiplicity[iHistType]);
		else fAllEventQAHistos->Add(fHistRefMultiplicity[iHistType]);

		// Centrality Histogram.
		fHistCentrality[iHistType] = new TH1F(Form("fHistCentrality%s",HistType[iHistType]),"Centrality;Centrality;Count",20,0,100);
		if (iHistType == 0) fSelectedEventQAHistos->Add(fHistCentrality[iHistType]);
		else fAllEventQAHistos->Add(fHistCentrality[iHistType]);

		// Centrality Quality.
		fHistCentralityQuality[iHistType] = new TH1F(Form("fHistCentralityQuality%s",HistType[iHistType]),"Centrality Quality;Quality;Count",2,-0.5,1.5);
		fHistCentralityQuality[iHistType]->GetXaxis()->SetBinLabel(1,"0");
		fHistCentralityQuality[iHistType]->GetXaxis()->SetBinLabel(2,"Other");
		if (iHistType == 0) fSelectedEventQAHistos->Add(fHistCentralityQuality[iHistType]);
		else fAllEventQAHistos->Add(fHistCentralityQuality[iHistType]);

		// VertexZ Histogram.
		fHistVertexZ[iHistType] = new TH1F(Form("fHistVertexZ%s",HistType[iHistType]),"VertexZ;z (cm);Count",60,-15.,15.);
		if (iHistType == 0) fSelectedEventQAHistos->Add(fHistVertexZ[iHistType]);
		else fAllEventQAHistos->Add(fHistVertexZ[iHistType]);

	}

}

// -------------------------------------------------------------------------
Bool_t AliAODEventCutsDiHadronPID::IsSelected(AliAODEvent* event) {
	
	if (fDebug > 1) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	if (!event) return kFALSE;

	if (!fAllEventQAHistos||!fSelectedEventQAHistos) {cout<<"AliAODEventCutsDiHadronPID - Histograms were not created, you should have called CreateHistos()..."<<endl;}

	// Input the event handler.
	AliInputEventHandler* InputHandler = (AliInputEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
	if (!InputHandler) return kFALSE;

	Bool_t select = kTRUE;

	// Test Trigger.
	UInt_t trigger = InputHandler->IsEventSelected();
	Int_t triggerselect = 0; // 0 = selected.
    if (fTestTrigger) {
		if (!(trigger & fTrigger)) {
			select = kFALSE;
			triggerselect = 1;	// 1 = not selected.
		}
	}	

	AliCentrality* CurrentCentrality = 0x0;
	Int_t CurrentCentralityQuality = -999;
	Float_t percentile = -999.;

	if (fIsPbPb) {

		// Get the centrality object.
	    CurrentCentrality = event->GetCentrality();
	    if (!CurrentCentrality) select = kFALSE;

	    // Check the quality of the centrality estimation.
	    // If 0 then quality is OK, c.f. TOF/PbPb276/macros/TOFmatchEff.C
	    CurrentCentralityQuality = CurrentCentrality->GetQuality();
	    //cout<<"Centrality: "<<CurrentCentrality->GetCentralityPercentile(fCentralityEstimator.Data())<<" Quality: "<<CurrentCentrality->GetQuality()<<endl;
	    if (CurrentCentralityQuality) select = kFALSE;

		// Test Centrality.
	    percentile = CurrentCentrality->GetCentralityPercentile(fCentralityEstimator.Data());
		if (fTestCentrality) {
	    	if ((percentile < fMaxCentrality)||(percentile > fMinCentrality)) select = kFALSE;
		}

	}

	// Get the primary vertex.
	AliAODVertex* CurrentPrimaryVertex = event->GetPrimaryVertex();
    if (!CurrentPrimaryVertex) select = kFALSE;

	// Test Vertex Z.
    Double_t vtxz = CurrentPrimaryVertex->GetZ();
	if (fTestVertexZ) {
    	if (TMath::Abs(vtxz) > fMaxVertexZ) select = kFALSE; 
	}

	// Get the event header.
	AliAODHeader* CurrentHeader = event->GetHeader();

	// Test minimum reference multiplicity.
	Int_t CurrentRefMultiplicity = CurrentHeader->GetRefMultiplicity();
	if (fTestMinRefMult) {
		if (CurrentRefMultiplicity < fMinRefMult) select = kFALSE;
	}

	// Fill the histograms for selected events.
	if (select) {
		fHistTrigger[0]->Fill(triggerselect);
		fHistRefMultiplicity[0]->Fill(CurrentHeader->GetRefMultiplicity());
		if (fIsPbPb) fHistCentrality[0]->Fill(percentile);
		if (fIsPbPb) fHistCentralityQuality[0]->Fill(CurrentCentralityQuality);
		fHistVertexZ[0]->Fill(vtxz);
	}

	// Fill the histograms for all events.
	fHistTrigger[1]->Fill(triggerselect);
	fHistRefMultiplicity[1]->Fill(CurrentHeader->GetRefMultiplicity());
	if (fIsPbPb) fHistCentrality[1]->Fill(percentile);
	if (fIsPbPb) fHistCentralityQuality[1]->Fill(CurrentCentralityQuality);
	fHistVertexZ[1]->Fill(vtxz);

	cout<<"Event Selected: "<<select<<endl;

	return select;

}

// -------------------------------------------------------------------------
void AliAODEventCutsDiHadronPID::PrintCuts() {

	// NOT IMPLEMENTED.
	if (fDebug > 1) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	return;

}
