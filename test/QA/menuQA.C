/*
 *  menuQA.C
 *  
 *
 *  Created by schutz on 01/08/08.
 *  Copyright 2008 CERN. All rights reserved.
 *fQA
 */
gROOT->Reset("a") ; 

TControlBar * fQA   = NULL ; 
TControlBar * fDet  = NULL ; 
TControlBar * fHist = NULL ; 
const Int_t fRun = 48868 ; 
AliQA * fQAResult = NULL ; 
TCanvas * fCa = NULL ; 

void menuQA()
{
	cout << AliQA::GetQAResultFileName()  << endl ; 
	TFile * qaResultFile = TFile::Open("QA.root") ; //AliQA::GetQAResultFileName() ) ; 
	if ( ! qaResultFile ) {
		printf("File %s not found in current directory\n", AliQA::GetQAResultFileName() ) ; 
		return ; 
	}
	fQAResult = (AliQA *)qaResultFile->Get("QA") ; 
	if ( ! fQAResult ) {
		printf("QA object not found in %s\n", AliQA::GetQAResultFileName() ) ; 
		return ; 
	}

	TFile * qaDataFile = TFile::Open(Form("Merged.QA.%d.root", fRun)) ; 
	if ( ! qaDataFile ) {
		printf("File Merged.QA.%d.root not found in current directory\n", fRun) ; 
		return ; 
	}
	if (fQA)
		delete fQA ; 
	TCanvas * cc = new TCanvas("test", "", 400, 800) ; 
	fQA = new TControlBar("vertical", Form("Active detectors in Run %d", fRun), 100, 20);
	fQA->SetButtonWidth(24) ;
	fQA->AddButton("Clean Screen", "Cls()", "Clean the screen");	
	// search active detectors
	TList * listOfDetectors = qaDataFile->GetListOfKeys() ; 
	for (Int_t det = 0 ; det < listOfDetectors->GetEntries() ; det++) {
		char * detName = listOfDetectors->At(det)->GetName() ; 
		if (fQAResult->IsSetAny(AliQA::GetDetIndex(detName)))
			char * buttonName = Form("QA SIGNALLED !! : %s", detName) ; 
		else 
			char * buttonName = Form("QA OK : %s", detName) ; 
	  fQA->AddButton(buttonName, Form("MakeDetMenu(\"%s\")", detName), Form("Display the QA histograms for %s", detName));	
	}
	fQA->Show();
}

void MakeDetMenu(char * detName)
{
	if (fDet) 
		delete fDet ;  
	if (fHist) 
		delete fHist ; 
	if (fCa) 
		delete fCa ; 
	fDet = new TControlBar("vertical", detName, 6000, 8000);
	fDet->SetButtonWidth(100) ; 
	// serach all the QA tasks
	TDirectory * save = gDirectory ; 
	gDirectory->cd(detName) ; 
	TList * listOfTasks = gDirectory->GetListOfKeys() ; 
	for (Int_t task = 0 ; task < listOfTasks->GetEntries() ; task++) {
		char * taskName = listOfTasks->At(task)->GetName() ; 
		if (fQAResult->IsSetAny(AliQA::GetDetIndex(detName), AliQA::GetTaskIndex(taskName)))
			char * buttonName = Form("QA SIGNALLED !! : %s", taskName) ; 
		else 
			char * buttonName = Form("QA OK : %s", taskName) ; 
			fDet->AddButton(buttonName, Form("MakeTaskMenu(\"%s\", \"%s\")", detName, taskName), Form("Display the QA histograms for %s", taskName));	
	}
	fDet->Show() ; 
	gDirectory = save ; 
}

void MakeTaskMenu(char * detName, char * taskName )
{
	if (fHist) 
		delete fHist ; 
	if (fCa)
		delete fCa ; 
	fHist = new TControlBar("vertical", Form("Found QA histo for %s/%s", detName, taskName));
	fHist->SetButtonWidth(100) ; 
	fHist->AddButton("ALL", Form("DisplayAll(\"%s\", \"%s\")", detName, taskName), Form("Display the QA histograms for %s", detName));
	TDirectory * save = gDirectory ; 
	gDirectory->cd(Form("%s/%s", detName, taskName)) ; 
	TList * listOfHistos = gDirectory->GetListOfKeys() ; 
	for (Int_t h = 0 ; h < listOfHistos->GetEntries() ; h++) {
		char * hName = listOfHistos->At(h)->GetName() ; 
		fHist->AddButton(hName, Form("Display(\"%s\")", hName), Form("Display the QA histograms %s", hName));
	}
	fHist->Show() ; 
	gDirectory = save ; 
}


void Display(char * hName)
{
	if (fCa)
		delete fCa ; 
	cout << Display << " " << hName << endl ; 
	fCa = new TCanvas(hName, "test", 800, 600) ;
	fCa->SetLogy() ; 
	TH1 * hh = dynamic_cast<TH1*>(gDirectory->FindObjectAny(hName));
	if (hh) {
		hh->Draw() ; 
		fCa->Modified() ; 
		fCa->Update() ;		
	}
}

void DisplayAll(char * detName, char * taskName)
{
	cout << Display << " " << detName << endl ; 
	TDirectory * save = gDirectory ; 
	gDirectory->cd(Form("%s/%s", detName, taskName)) ; 
	TList * listOfHistos = gDirectory->GetListOfKeys() ; 
	Int_t nHisto = listOfHistos->GetEntries() ; 
	Int_t ny = TMath::Sqrt(nHisto) ; 
	Int_t nx = nHisto / ny + 1 ;
	if (fCa) 
		delete fCa ; 
	fCa = new TCanvas(Form("QA %s in %s", taskName, detName), Form("QA %s in %s", taskName, detName), nx*300, ny*300) ;
	fCa->Divide(nx, ny) ;
	for (Int_t h = 0 ; h < listOfHistos->GetEntries() ; h++) {
		char * hName = listOfHistos->At(h)->GetName() ; 
		TH1 * hh = dynamic_cast<TH1*>(gDirectory->FindObjectAny(hName));
		TPad * pad = fCa->cd(h+1) ;
		pad->SetLogy() ; 
		if (hh) {
			hh->Draw() ; 
			fCa->Modified() ; 
			fCa->Update() ;			
		}
	}
	fCa->Modified() ; 
	fCa->Update() ;
	gDirectory = save ; 
}

void Cls()
{
	if (fDet) 
		delete fDet ;  
	if (fHist) 
		delete fHist ; 
	if (fCa) 
		delete fCa ; 
}
