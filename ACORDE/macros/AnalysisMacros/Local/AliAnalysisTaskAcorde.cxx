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
//	Create some histograms and one tree with the information
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
	cosmicTree(0),
	nTracks(0),
	nMatchTracks(0),
	minTPCclusters(30),
	minTrackDist(1000.),
	minCutDir(0.95),
  	xAco(0),
  	yAco(0),
  	zAco(0),
	trigger(0),
	ActiveTriggerDetector(0),
	nSingleTracks(0),
	nMatchedTracks(0),
	histo(0),
 	ntracks(0),
	acordeHitsAll(0),
	acordeMultiAll(0),
	acordeHitsTPC(0),
	acordeMultiTPC(0),
	nTracksMatched(0),
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
	delete histo;
	delete cosmicTree;
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

	// creates one TList with some histograms

	histo = new TList();

	nTracksMatched = new TH1F("nTracksMatched","nTracksMatched from matching algorithm implemented",300,1,300);
	ntracks = new TH1F("ntracks","nTracks distribution",500,1,500);
	acordeHitsAll = new TH1F("acordeHitsAll","Hits of ACORDE from ESD",60,-0.5,59.5);
	acordeHitsTPC = new TH1F("acordeHitsAllTPC","Hits of ACORDE from ESD (if ntracks>0)",61,0,60);
	acordeMultiAll = new TH1F("acordeMultiAll","Multiplicity of ACORDE modules from ESD",60,-0.5,59.5);
	acordeMultiTPC = new TH1F("acordeMultiAllTPC","Multiplicity of ACORDE modules from ESD (id ntracks>0)",61,0,60);
	fTracksToAcorde = new TH2F("fTracksToAcorde","Extrapolated tracks to ACORDE x:z",1200,-600,600,1200,-600,600);

	histo->Add(nTracksMatched);
	histo->Add(ntracks);
	histo->Add(acordeHitsAll);
	histo->Add(acordeHitsTPC);
	histo->Add(acordeMultiAll);
	histo->Add(acordeMultiTPC);
	histo->Add(fTracksToAcorde);

	// Create Tree branches
  	// Called just once

  	cosmicTree = new TTree("cosmicTree","cosmicTreeMRC");
	cosmicTree->Branch("nTracks",&nTracks,"nTracks/I");
	cosmicTree->Branch("nMatchTracks",&nMatchTracks,"nMatchTracks/I");
	cosmicTree->Branch("xAco",&xAco,"xAco/F");
	cosmicTree->Branch("yAco",&yAco,"yAco/F");
	cosmicTree->Branch("zAco",&zAco,"zAco/F");
	cosmicTree->Branch("trigger",&trigger,"trigger/I");

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
			acordeHitsAll->Fill(imod);
			contMulti++;
		}

	}acordeMultiAll->Fill(contMulti);

	for(Int_t imod=0;imod<60;imod++)
	{
		if (acordeESD->GetHitChannel(imod))
		{
			acordeHitsTPC->Fill(imod);
			contMultiTPC++;
		}
	}acordeMultiTPC->Fill(contMultiTPC);


	// Assingment of the type of Trigger to cosmicTree

	ActiveTriggerDetector = fESD->GetFiredTriggerClasses();

	if (ActiveTriggerDetector.Contains("C0SCO")) trigger=0;
	if (ActiveTriggerDetector.Contains("C0SH2")) trigger=1;
  	if (ActiveTriggerDetector.Contains("C0AMU")) trigger=2;
  	if (ActiveTriggerDetector.Contains("C0LSR")) trigger=3;
  	if (ActiveTriggerDetector.Contains("C0SCO1")) trigger=4;
  	if (ActiveTriggerDetector.Contains("C0OBE")) trigger=5;
  	if (ActiveTriggerDetector.Contains("C0PCO")) trigger=6;
  	if (ActiveTriggerDetector.Contains("C0OB3")) trigger=7;
  	if (ActiveTriggerDetector.Contains("C0OCP")) trigger=8;
    	if (ActiveTriggerDetector.Contains("C0ASL")) trigger=9;
	

	// Begins the matching analysis between tracks

	nTracks = fESD->GetNumberOfTracks();
	if (nTracks<=0) return;
	ntracks->Fill(nTracks);

	fPair = new TArrayI(nTracks);
	nSingleTracks=0;
	nMatchedTracks=0;
	Int_t muons = 0;
	TFile *curfile = AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile();
	for (Int_t iTrack=0; iTrack<nTracks; iTrack++)
	{

		// Special case: nTracks == 1
		if (nTracks ==1)
		{
			AliESDtrack *singleTrack = fESD->GetTrack(iTrack);
			if (!singleTrack || !singleTrack->GetOuterParam() || !singleTrack->GetInnerParam()) continue;
			if (singleTrack->GetTPCNcls()<minTPCclusters) continue;
			Float_t vertX = singleTrack->Xv();
			Float_t vertZ = singleTrack->Zv();
			if (TMath::Abs(vertX)>165 || TMath::Abs(vertZ)>245) continue; // far away from the VERTEX
			if (singleTrack->E()< 1) continue; // less than 1 GeV 
			if (singleTrack->P()< 1) continue; // less than 1 GeV 
			nMatchTracks=0;
			nSingleTracks = nTracks;
			muons = 1;

		} // nTracks ==1
		else if (nTracks>1) // nTracks > 1
		{

			AliESDtrack *track0 = fESD->GetTrack(iTrack);
			(*fPair)[iTrack]=-1;
			if (!track0 || !track0->GetOuterParam() || !track0->GetInnerParam()) continue;

			Double_t track0_dir[3];
			track0->GetDirection(track0_dir);
			Float_t minDist = minTrackDist;	

			const AliExternalTrackParam * trackIn = track0->GetInnerParam();
			const AliExternalTrackParam * trackOut = track0->GetOuterParam();

			if (!trackIn || !trackOut) continue;
			if (nTracks>4 && TMath::Abs(trackIn->GetTgl())<0.0015) continue; // Filter TPC-Laser

			if (track0->GetTPCNcls()<minTPCclusters) continue;

			for (Int_t jTrack=iTrack; jTrack<nTracks; ++jTrack)
			{
				AliESDtrack *track1 = fESD->GetTrack(jTrack);
				if (!track1 || !track1->GetOuterParam() || !track1->GetInnerParam()) continue;
				Double_t track1_dir[3];
				track1->GetDirection(track1_dir);
				if (track1->GetTPCNcls()<minTPCclusters) continue;
				Float_t direction = track0_dir[0]*track1_dir[0] + track0_dir[1]*track1_dir[1] + track0_dir[2]*track1_dir[2];
				if (TMath::Abs(direction)<minCutDir) continue; // direction vector product
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
		} // nTracks > 1
	

	} // Loop over all the tracks




	// Matching tracks
	
	for (Int_t itrack=0;itrack<nTracks;itrack++)
	{
		AliESDtrack *UpTrack = 0;
		AliESDtrack *DownTrack = 0;

		AliESDtrack *Track0 = fESD->GetTrack(itrack);
		if (!Track0 || !Track0->GetOuterParam() || !Track0->GetInnerParam() ) continue;

		Double_t mxyz[3];
		Track0->GetOuterParam()->GetXYZ(mxyz);
		
		UpTrack = Track0;
	

		if ((*fPair)[itrack]>0)
		{
			AliESDtrack *Track1 = fESD->GetTrack((*fPair)[itrack]);
			if (!Track1 || !Track1->GetOuterParam() || !Track1->GetInnerParam()) continue;
			Float_t dvertex0[3];
			Float_t dvertex1[3];

			dvertex0[0] = Track0->Xv();
			dvertex0[1] = Track0->Yv();
			dvertex0[2] = Track0->Zv();

			dvertex1[0] = Track1->Xv();
			dvertex1[1] = Track1->Yv();
			dvertex1[2] = Track1->Zv();

			Float_t distTracks = TMath::Sqrt((dvertex0[0]-dvertex1[0])*(dvertex0[0]-dvertex1[0]) + (dvertex0[1]-dvertex1[1])*(dvertex0[1]-dvertex1[1]) + (dvertex0[2]-dvertex1[2])*(dvertex0[2]-dvertex1[2]) );
			if (distTracks>2.5) continue;

			if (Track1->GetOuterParam())
			{
				Double_t nxyz[3];
				Track1->GetOuterParam()->GetXYZ(nxyz);
				if (nxyz[1]>mxyz[1])
				{
					UpTrack = Track1;
				}else DownTrack = Track1;
			}



			// Here we propagate the Up-Tracks to ACORDE


			const Double_t kRL3 = 510; // radious of L3 magnet
			const Double_t kxAcorde = 850.; // radios of "ACORDE" above the magnet

			Double_t agxyz[3];
			

			AliExternalTrackParam *upper = (AliExternalTrackParam *)(UpTrack->GetOuterParam()->Clone());

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

			if (isOk) isOk = upper->PropagateTo(kxAcorde,0);

			TVectorD rgxyz(3);
			AliExternalTrackParam *param = (AliExternalTrackParam *)upper;
			param->GetXYZ(rgxyz.GetMatrixArray());
			xAco = rgxyz[0];
			yAco = rgxyz[1];
			zAco = rgxyz[2];
			fTracksToAcorde->Fill(xAco,zAco);
			// Count how many tracks have been matched
			nMatchedTracks++;
		} else nSingleTracks++;
	}

	if ( (nMatchedTracks+nSingleTracks) <= nTracks ) muons=nMatchedTracks+nSingleTracks;
	nTracksMatched->Fill(muons);

	nMatchTracks = nSingleTracks + nMatchedTracks;
	
	// Fills the tree

	cosmicTree->Fill();
	
	// Post output data.

	PostData(0,cosmicTree);
	PostData(1,histo);
}      

//________________________________________________________________________
void AliAnalysisTaskAcorde::Terminate(Option_t *) 
{

//	AliDebug(1,"Do nothig in Terminate");

	//nPart->Draw();
        TCanvas *c1 = new TCanvas("c1","TPC-resolution for P reconstruction .... MRC");

        c1->Divide(2,2);
	c1->cd(1);
	acordeHitsAll->Draw();
	c1->cd(2)->SetLogy();
	acordeMultiAll->Draw();
	c1->cd(3)->SetLogy();
	//nTracksMatched->Draw();
	ntracks->Draw();
	c1->cd(4);
	fTracksToAcorde->Draw();
	
	


}
