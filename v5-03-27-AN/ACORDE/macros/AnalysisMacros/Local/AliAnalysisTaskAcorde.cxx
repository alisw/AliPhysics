/////////////////////////////////////////////////////////////////////////////
//
// AliAnalysisTaskAcorde class
//
// Description:
//
//	Reads the information of ACORDE-ESD
//	Also it is implemented a simple matching between tracks
//	to look for the extrapolation of them to ACORDE modules
//
//	Create some fHistograms and one tree with the information
//	of the matching tracks
//
//  Author: Mario Rodríguez Cahuantzi
//		<mario.rocah@gmail.com>
//		<mrodrigu@mail.cern.ch>
//
//  Created: June 30th. 2010 @ FCFM - BUAP, Puebla, MX
//  Last update: created
//
/////////////////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliAnalysisTaskAcorde.h"

#include "TMath.h"
#include "TArrayI.h"
#include "TDatabasePDG.h"
#include "TVectorD.h"

#include "AliESDtrack.h"
#include "AliTracker.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
ClassImp(AliAnalysisTaskAcorde)
//________________________________________________________________________
AliAnalysisTaskAcorde::AliAnalysisTaskAcorde(const char *name) 
  : AliAnalysisTask(name, ""), 
	fESD(0),
	fCosmicTree(0),
	fnTracks(0),
	fNMatchTracks(0),
	fkMinTPCclusters(30),
	fkMinTrackDist(1000.),
	fkMinCutDir(0.95),
  	fXAco(0),
  	fYAco(0),
  	fZAco(0),
	fTrigger(0),
	fActiveTriggerDetector(0),
	fNSingleTracks(0),
	fNMatchedTracks(0),
	fHisto(0),
 	fNTracks(0),
	fAcordeHitsAll(0),
	fAcordeMultiAll(0),
	fAcordeHitsTPC(0),
	fAcordeMultiTPC(0),
	fNTracksMatched(0),
	fTracksToAcorde(0)

{
	// Constructor of class

  	DefineInput(0, TChain::Class());
    	DefineOutput(0,TTree::Class());
	DefineOutput(1,TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskAcorde::~AliAnalysisTaskAcorde()
{
	// destructor --> avoid watchdog mails?
	
	delete fESD;
	delete fHisto;
	delete fCosmicTree;
}
//________________________________________________________________________
void AliAnalysisTaskAcorde::ConnectInputData(Option_t *) 
{
	// Connect ESD or AOD here
	// Called once

	TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  	if (!tree) 
	{
    		Printf("ERROR: Could not read chain from input slot 0");
  	}else 
	{

    		AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    		if (!esdH) 
		{
      			Printf("ERROR: Could not get ESDInputHandler");
    		}else fESD = esdH->GetEvent();
  	}
}

//________________________________________________________________________
void AliAnalysisTaskAcorde::CreateOutputObjects()
{

	// creates one TList with some fHistograms

	fHisto = new TList();

	fNTracksMatched = new TH1F("fNTracksMatched","fNTracksMatched from matching algorithm implemented",300,1,300);
	fNTracks = new TH1F("fNTracks","fnTracks distribution",500,1,500);
	fAcordeHitsAll = new TH1F("fAcordeHitsAll","Hits of ACORDE from ESD",60,-0.5,59.5);
	fAcordeHitsTPC = new TH1F("fAcordeHitsAllTPC","Hits of ACORDE from ESD (if fNTracks>0)",61,0,60);
	fAcordeMultiAll = new TH1F("fAcordeMultiAll","Multiplicity of ACORDE modules from ESD",60,-0.5,59.5);
	fAcordeMultiTPC = new TH1F("fAcordeMultiAllTPC","Multiplicity of ACORDE modules from ESD (id fNTracks>0)",61,0,60);
	fTracksToAcorde = new TH2F("fTracksToAcorde","Extrapolated tracks to ACORDE x:z",1200,-600,600,1200,-600,600);

	fHisto->Add(fNTracksMatched);
	fHisto->Add(fNTracks);
	fHisto->Add(fAcordeHitsAll);
	fHisto->Add(fAcordeHitsTPC);
	fHisto->Add(fAcordeMultiAll);
	fHisto->Add(fAcordeMultiTPC);
	fHisto->Add(fTracksToAcorde);

	// Create Tree branches
  	// Called just once

  	fCosmicTree = new TTree("fCosmicTree","fCosmicTreeMRC");
	fCosmicTree->Branch("fnTracks",&fnTracks,"fnTracks/I");
	fCosmicTree->Branch("fNMatchTracks",&fNMatchTracks,"fNMatchTracks/I");
	fCosmicTree->Branch("fXAco",&fXAco,"fXAco/F");
	fCosmicTree->Branch("fYAco",&fYAco,"fYAco/F");
	fCosmicTree->Branch("fZAco",&fZAco,"fZAco/F");
	fCosmicTree->Branch("fTrigger",&fTrigger,"fTrigger/I");

}


void AliAnalysisTaskAcorde::Exec(Option_t *) 
{
	// Main loop
  	// Called for each event

	Int_t counterOfMuons = 0;



  	if (!fESD) 
	{
    		Printf("ERROR: fESD not available");
    		return;
  	}



	// Pointer to the information of ACORDE detector

	AliESDACORDE *acordeESD = fESD->GetACORDEData();

	Int_t contMulti = 0;
	Int_t contMultiTPC = 0;

	for(Int_t imod=0;imod<60;imod++)
	{
		if (acordeESD->GetHitChannel(imod)) 
		{
			fAcordeHitsAll->Fill(imod);
			contMulti++;
		}

	}fAcordeMultiAll->Fill(contMulti);

	for(Int_t imod=0;imod<60;imod++)
	{
		if (acordeESD->GetHitChannel(imod))
		{
			fAcordeHitsTPC->Fill(imod);
			contMultiTPC++;
		}
	}fAcordeMultiTPC->Fill(contMultiTPC);


	// Assingment of the type of Trigger to fCosmicTree

	fActiveTriggerDetector = fESD->GetFiredTriggerClasses();

	if (fActiveTriggerDetector.Contains("C0SCO")) fTrigger=0;
	if (fActiveTriggerDetector.Contains("C0SH2")) fTrigger=1;
  	if (fActiveTriggerDetector.Contains("C0AMU")) fTrigger=2;
  	if (fActiveTriggerDetector.Contains("C0LSR")) fTrigger=3;
  	if (fActiveTriggerDetector.Contains("C0SCO1")) fTrigger=4;
  	if (fActiveTriggerDetector.Contains("C0OBE")) fTrigger=5;
  	if (fActiveTriggerDetector.Contains("C0PCO")) fTrigger=6;
  	if (fActiveTriggerDetector.Contains("C0OB3")) fTrigger=7;
  	if (fActiveTriggerDetector.Contains("C0OCP")) fTrigger=8;
    	if (fActiveTriggerDetector.Contains("C0ASL")) fTrigger=9;
	

	// Begins the matching analysis between tracks

	fnTracks = fESD->GetNumberOfTracks();
	if (fnTracks<=0) return;
	fNTracks->Fill(fnTracks);

	fPair = new TArrayI(fnTracks);
	fNSingleTracks=0;
	fNMatchedTracks=0;
	Int_t muons = 0;
	TFile *curfile = AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile();
	for (Int_t iTrack=0; iTrack<fnTracks; iTrack++)
	{

		// Special case: fnTracks == 1
		if (fnTracks ==1)
		{
			AliESDtrack *singleTrack = fESD->GetTrack(iTrack);
			if (!singleTrack || !singleTrack->GetOuterParam() || !singleTrack->GetInnerParam()) continue;
			if (singleTrack->GetTPCNcls()<fkMinTPCclusters) continue;
			Float_t vertX = singleTrack->Xv();
			Float_t vertZ = singleTrack->Zv();
			if (TMath::Abs(vertX)>165 || TMath::Abs(vertZ)>245) continue; // far away from the VERTEX
			if (singleTrack->E()< 1) continue; // less than 1 GeV 
			if (singleTrack->P()< 1) continue; // less than 1 GeV 
			fNMatchTracks=0;
			fNSingleTracks = fnTracks;
			muons = 1;

		} // fnTracks ==1
		else if (fnTracks>1) // fnTracks > 1
		{

			AliESDtrack *track0 = fESD->GetTrack(iTrack);
			(*fPair)[iTrack]=-1;
			if (!track0 || !track0->GetOuterParam() || !track0->GetInnerParam()) continue;

			Double_t trackdir0[3];
			track0->GetDirection(trackdir0);
			Float_t minDist = fkMinTrackDist;	

			const AliExternalTrackParam * trackIn = track0->GetInnerParam();
			const AliExternalTrackParam * trackOut = track0->GetOuterParam();

			if (!trackIn || !trackOut) continue;
			if (fnTracks>4 && TMath::Abs(trackIn->GetTgl())<0.0015) continue; // Filter TPC-Laser

			if (track0->GetTPCNcls()<fkMinTPCclusters) continue;

			for (Int_t jTrack=iTrack; jTrack<fnTracks; ++jTrack)
			{
				AliESDtrack *track1 = fESD->GetTrack(jTrack);
				if (!track1 || !track1->GetOuterParam() || !track1->GetInnerParam()) continue;
				Double_t trackdir1[3];
				track1->GetDirection(trackdir1);
				if (track1->GetTPCNcls()<fkMinTPCclusters) continue;
				Float_t direction = trackdir0[0]*trackdir1[0] + trackdir0[1]*trackdir1[1] + trackdir0[2]*trackdir1[2];
				if (TMath::Abs(direction)<fkMinCutDir) continue; // direction vector product
				Float_t dvertex0[3];
				Float_t dvertex1[3];

				dvertex0[0] = track0->Xv();
				dvertex0[1] = track0->Yv();
				dvertex0[2] = track0->Zv();

				dvertex1[0] = track1->Xv();
				dvertex1[1] = track1->Yv();
				dvertex1[2] = track1->Zv();

				if (TMath::Abs(dvertex0[0])>165 || TMath::Abs(dvertex0[2])>245 || TMath::Abs(dvertex1[0])>165 || TMath::Abs(dvertex0[2])>245) continue;
				// calculate the distance between track0 and track1
				
				Float_t dy = track0->GetY()+track1->GetY();
				Float_t sy2 = track0->GetSigmaY2()+track1->GetSigmaY2();
				Float_t dphi = track0->GetAlpha()-track1->GetAlpha()-TMath::Pi();
				Float_t sphi2 = track0->GetSigmaSnp2()+track1->GetSigmaSnp2();
				Float_t dtheta = track0->GetTgl()+track1->GetTgl();
				Float_t stheta2 = track0->GetSigmaTgl2()+track1->GetSigmaTgl2();
				Float_t normDist = TMath::Sqrt(dy*dy/sy2+dphi+dphi/sphi2+dtheta*dtheta/stheta2);
				if (normDist>minDist) continue;

				Float_t distTracks = TMath::Sqrt((dvertex0[0]-dvertex1[0])*(dvertex0[0]-dvertex1[0]) + (dvertex0[1]-dvertex1[1])*(dvertex0[1]-dvertex1[1]) + (dvertex0[2]-dvertex1[2])*(dvertex0[2]-dvertex1[2]) );
				if (distTracks>2.5) continue;


				minDist = normDist;

				// after all the cuts we should have only tracks that can be matched

				(*fPair)[iTrack] = jTrack;


			}	
		} // fnTracks > 1
	

	} // Loop over all the tracks




	// Matching tracks
	
	for (Int_t itrack=0;itrack<fnTracks;itrack++)
	{
		AliESDtrack *upTrack = 0;
		AliESDtrack *downTrack = 0;

		AliESDtrack *track0 = fESD->GetTrack(itrack);
		if (!track0 || !track0->GetOuterParam() || !track0->GetInnerParam() ) continue;

		Double_t mxyz[3];
		track0->GetOuterParam()->GetXYZ(mxyz);
		
		upTrack = track0;
	

		if ((*fPair)[itrack]>0)
		{
			AliESDtrack *track1 = fESD->GetTrack((*fPair)[itrack]);
			if (!track1 || !track1->GetOuterParam() || !track1->GetInnerParam()) continue;
			Float_t dvertex0[3];
			Float_t dvertex1[3];

			dvertex0[0] = track0->Xv();
			dvertex0[1] = track0->Yv();
			dvertex0[2] = track0->Zv();

			dvertex1[0] = track1->Xv();
			dvertex1[1] = track1->Yv();
			dvertex1[2] = track1->Zv();

			Float_t distTracks = TMath::Sqrt((dvertex0[0]-dvertex1[0])*(dvertex0[0]-dvertex1[0]) + (dvertex0[1]-dvertex1[1])*(dvertex0[1]-dvertex1[1]) + (dvertex0[2]-dvertex1[2])*(dvertex0[2]-dvertex1[2]) );
			if (distTracks>2.5) continue;

			if (track1->GetOuterParam())
			{
				Double_t nxyz[3];
				track1->GetOuterParam()->GetXYZ(nxyz);
				if (nxyz[1]>mxyz[1])
				{
					upTrack = track1;
				}else downTrack = track1;
			}



			// Here we propagate the Up-Tracks to ACORDE


			const Double_t kRL3 = 510; // radious of L3 magnet
			const Double_t kfXAcorde = 850.; // radios of "ACORDE" above the magnet

			Double_t agxyz[3];
			

			AliExternalTrackParam *upper = (AliExternalTrackParam *)(upTrack->GetOuterParam()->Clone());

			Bool_t isOk = upper->PropagateTo(kRL3,fESD->GetMagneticField());
			upper->GetXYZ(agxyz);
			if (agxyz[1]<0) continue;

			upper->GetXYZ(agxyz);
			Double_t galpha = TMath::ATan2(agxyz[1],agxyz[0]);
			Double_t alpha = galpha;
			galpha*=180/TMath::Pi();

			if (galpha<45.) alpha = TMath::Pi()/8;
			if (galpha>135.) alpha = TMath::Pi()*(7/8);
			if (galpha>45 && galpha<135.) alpha = TMath::Pi()/2;

			if (isOk) upper->Rotate(alpha);

			if (isOk) isOk = upper->PropagateTo(kfXAcorde,0);

			TVectorD rgxyz(3);
			AliExternalTrackParam *param = (AliExternalTrackParam *)upper;
			param->GetXYZ(rgxyz.GetMatrixArray());
			fXAco = rgxyz[0];
			fYAco = rgxyz[1];
			fZAco = rgxyz[2];
			fTracksToAcorde->Fill(fXAco,fZAco);
			// Count how many tracks have been matched
			fNMatchedTracks++;
		} else fNSingleTracks++;
	}

	if ( (fNMatchedTracks+fNSingleTracks) <= fnTracks ) muons=fNMatchedTracks+fNSingleTracks;
	fNTracksMatched->Fill(muons);

	fNMatchTracks = fNSingleTracks + fNMatchedTracks;
	
	// Fills the tree

	fCosmicTree->Fill();
	
	// Post output data.

	PostData(0,fCosmicTree);
	PostData(1,fHisto);
}      

//________________________________________________________________________
void AliAnalysisTaskAcorde::Terminate(Option_t *) 
{

//	AliDebug(1,"Do nothig in Terminate");

	//nPart->Draw();
        TCanvas *c1 = new TCanvas("c1","TPC-resolution for P reconstruction .... MRC");

        c1->Divide(2,2);
	c1->cd(1);
	fAcordeHitsAll->Draw();
	c1->cd(2)->SetLogy();
	fAcordeMultiAll->Draw();
	c1->cd(3)->SetLogy();
	//fNTracksMatched->Draw();
	fNTracks->Draw();
	c1->cd(4);
	fTracksToAcorde->Draw();
	
	


}
