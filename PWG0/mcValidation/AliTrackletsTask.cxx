/* $Id:$ */

// -----------------------------------------------
// Task to extract distributions 
// for traclets paramters
// for a quick comparison 
// between MC and data
// eta, phi, deltaPhi, deltaEta, vtxX, vtxY, vtxZ, SPD GFO
// -----------------------------------------------


#include "AliTrackletsTask.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH1I.h>

#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliAnalysisManager.h>
#include <AliESDInputHandler.h>
#include <AliESDHeader.h>
#include <AliESDVertex.h>
#include <AliMultiplicity.h>

ClassImp(AliTrackletsTask)

AliTrackletsTask::AliTrackletsTask() :
	AliAnalysisTask("AliTrackletsTask", ""),
	fESD(0),
	fOutput(0),
	fNtracks(0x0),
	fPhi(0x0),
	fEtaPhi(0x0),
	fDeltaPhi(0x0),
	fDeltaTheta(0x0),
	fVtxX(0x0),
	fVtxY(0x0),
	fVtxZ(0x0),
	fVtx(0x0),
	fVtxContributors(0x0)
{
	//
	// Constructor. Initialization of pointers
	//
	
	// Define input and output slots here
	DefineInput(0, TChain::Class());
	DefineOutput(0, TList::Class());
}

//------------------------------------------------------------------

AliTrackletsTask::~AliTrackletsTask()
{
	//
	// Destructor
	//
	
	// histograms are in the output list and deleted when the output
	// list is deleted by the TSelector dtor
	
	if (fOutput) {
		delete fOutput;
		fOutput = 0;
	}
}

//--------------------------------------------------------------------

void AliTrackletsTask::ConnectInputData(Option_t *)
{
	//
	// Connect ESD
	// 
	
	Printf("AliTrackletsTask::ConnectInputData called");
	
	AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	
	if (!esdH) {
		Printf("ERROR: Could not get ESDInputHandler");
	} else {
		fESD = esdH->GetEvent();
		
		TString branches("AliESDHeader Tracks AliMultiplicity");
		
		// Enable only the needed branches
		esdH->SetActiveBranches(branches);
	}
}

//---------------------------------------------------------------------

void AliTrackletsTask::CreateOutputObjects()
{
	//
	// create result objects and add to output list
	//
	
	Printf("AliTrackletsTask::CreateOutputObjects");
	
	fOutput = new TList;
	fOutput->SetOwner();
	
	fNtracks = new TH1I("fNtracks", "n. of ESD tracks", 1000, 0, 1000);
	fOutput->Add(fNtracks);

	fPhi = new TH1D("fPhi", "Phi (rad)", 720,0,2*TMath::Pi());
	fOutput->Add(fPhi);

	fEtaPhi = new TH2D("fEtaPhi", "Phi vs Eta;#eta;#phi [rad];count", 80, -4, 4, 18*5,0,2*TMath::Pi());
	fOutput->Add(fEtaPhi);

	fVtxX = new TH1D("fVtxX", "x SPD primary vertex; x [cm]; count", 100, -1, 1);
	fOutput->Add(fVtxX);

	fVtxY = new TH1D("fVtxY", "y SPD primary vertex; y [cm]; count", 100, -1, 1);
	fOutput->Add(fVtxY);

	fVtxZ = new TH1D("fVtxZ", "z SPD primary vertex; z [cm]; count", 100, -30, 30);
	fOutput->Add(fVtxZ);

	fVtx = new TH3D("fVtx", "SPD primary vertex; x [cm]; y [cm]; z [cm]; count", 100, -1, 1, 100,-1,1,100,-30,30);
	fOutput->Add(fVtx);
	
	fVtxContributors = new TH3D("fVtxContributors", "SPD primary vertex with N contributors > 0; x [cm]; y [cm]; z [cm]; count", 100, -1, 1, 100,-1,1,100,-30,30);
	fOutput->Add(fVtxContributors);
	
	fDeltaPhi = new TH1D("fDeltaPhi", "fDeltaPhi;#Delta #phi;Entries", 500, -1, 1);
	fOutput->Add(fDeltaPhi);

	fDeltaTheta = new TH1D("fDeltaTheta", "fDeltaTheta;#Delta #theta;Entries", 500, -0.2, 0.2);
	fOutput->Add(fDeltaTheta);

}

//----------------------------------------------------------------------

void AliTrackletsTask::Exec(Option_t*)
{
	//
	// process the event
	//

	PostData(0, fOutput);

	if (!fESD){
		AliError("ESD branch not available");
		return;
	}
  
	AliESDHeader* esdHeader = fESD->GetHeader();
	if (!esdHeader){
		Printf("ERROR: esdHeader could not be retrieved");
		return;
	}

	Bool_t isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
	if (isSelected == kFALSE) {
		AliInfo("Event not selected by Physics analysis");
		return;
	}

	Int_t ntrk = fESD->GetNumberOfTracks() ;
	fNtracks->Fill(ntrk);

	const AliMultiplicity* mult = fESD->GetMultiplicity();
        if (!mult){
		AliError("AliMultiplicity not found, returning");
		return;
	}

	Int_t nTracklets = mult->GetNumberOfTracklets();

	for (Int_t iTracklet = 0; iTracklet < nTracklets; iTracklet++){
		Float_t eta = mult->GetEta(iTracklet);
                Float_t phi = mult->GetPhi(iTracklet);
                Float_t deltaPhi = mult->GetDeltaPhi(iTracklet);
                Float_t deltaTheta = mult->GetDeltaTheta(iTracklet);
		if (phi < 0) phi += TMath::Pi() * 2;
		fPhi->Fill(phi);
		fEtaPhi->Fill(eta, phi); 
		fDeltaPhi->Fill(deltaPhi); 
		fDeltaTheta->Fill(deltaTheta); 
	}

	// vertex SPD
	const AliESDVertex* vertexSPD = 0;
	vertexSPD = fESD->GetPrimaryVertexSPD();
	if (vertexSPD && (!(vertexSPD->IsFromVertexerZ() && vertexSPD->GetDispersion() > 0.02))) {
		Double_t vtxX = vertexSPD->GetX();
		Double_t vtxY = vertexSPD->GetY();
		Double_t vtxZ = vertexSPD->GetZ();
		fVtxX->Fill(vtxX);
		fVtxY->Fill(vtxY);
		fVtxZ->Fill(vtxZ);
		fVtx->Fill(vtxX,vtxY,vtxZ);
		if (vertexSPD->GetNContributors() > 0){
			fVtxContributors->Fill(vtxX,vtxY,vtxZ);
		}
	}


}

//----------------------------------------------------------------------

void AliTrackletsTask::Terminate(Option_t *)
{
	//
	// Plotting distributions, and saving them in a file
	//

	fOutput = dynamic_cast<TList*> (GetOutputData(0));
	if (!fOutput)
		Printf("ERROR: fOutput not available");
    
	fOutput->Print();

	if (fOutput){
		fNtracks = dynamic_cast<TH1I*>(fOutput->FindObject("fNtracks"));
		fPhi = dynamic_cast<TH1D*>(fOutput->FindObject("fPhi"));
		fEtaPhi = dynamic_cast<TH2D*>(fOutput->FindObject("fEtaPhi"));
		fDeltaPhi = dynamic_cast<TH1D*>(fOutput->FindObject("fDeltaPhi"));
		fDeltaTheta = dynamic_cast<TH1D*>(fOutput->FindObject("fDeltaTheta"));
		fVtxX = dynamic_cast<TH1D*>(fOutput->FindObject("fVtxX"));
		fVtxY = dynamic_cast<TH1D*>(fOutput->FindObject("fVtxY"));
		fVtxZ = dynamic_cast<TH1D*>(fOutput->FindObject("fVtxZ"));
		fVtx = dynamic_cast<TH3D*>(fOutput->FindObject("fVtx"));
		fVtxContributors = dynamic_cast<TH3D*>(fOutput->FindObject("fVtxContributors"));

	        new TCanvas("Ntracks", " Ntracks ",50, 50, 550, 550) ;
		fNtracks->Draw();
	        new TCanvas("phi", " phi ",50, 50, 550, 550) ;
		fPhi->Draw();
	        new TCanvas("etaphi", " etaphi ",50, 50, 550, 550) ;
		fEtaPhi->Draw();
	        new TCanvas("deltaPhi", " deltaPhi ",50, 50, 550, 550) ;
		fDeltaPhi->Draw();
	        new TCanvas("deltaTheta", " deltaTheta ",50, 50, 550, 550) ;
		fDeltaTheta->Draw();
	        new TCanvas("vtxX", " vtxX ",50, 50, 550, 550) ;
		fVtxX->Draw();
	        new TCanvas("vtxY", " vtxY ",50, 50, 550, 550) ;
		fVtxY->Draw();
	        new TCanvas("vtxZ", " vtxZ ",50, 50, 550, 550) ;
		fVtxZ->Draw();
		
		TFile* outputFile = new TFile("histograms.root", "RECREATE");
		fNtracks->Write();
		fPhi->Write();
		fEtaPhi->Write();
		fDeltaPhi->Write();
		fDeltaTheta->Write();
		fVtxX->Write();
		fVtxY->Write();
		fVtxZ->Write();
		fVtx->Write();
		fVtxContributors->Write();
		outputFile->Write();
		outputFile->Close();

	}
}





