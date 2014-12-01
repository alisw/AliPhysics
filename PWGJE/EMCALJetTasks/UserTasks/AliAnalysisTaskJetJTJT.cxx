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
//
//
// Jet fragmentation transverse momentum (j_T) analysis task
//
// Author: T.Snellman

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TVector.h>
#include <TGraphErrors.h>
#include <TGrid.h>
#include <TSystem.h>
#include <TFile.h>

#include "AliCentrality.h"



#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"

#include "AliAnalysisTaskJetJTJT.h"


ClassImp(AliAnalysisTaskJetJTJT)

	//________________________________________________________________________
	AliAnalysisTaskJetJTJT::AliAnalysisTaskJetJTJT() : 
		AliAnalysisTaskEmcalJet("AliAnalysisTaskJetJTJT", kTRUE),
		fHistTracksPt(0),
		fHistTracksJt(0),
		fHistClustersPt(0),
		fHistLeadingJetPt(0),
		fHistJetsPt(0),
		fHistBackgroundDone(0),
		fHistJTPta(0),
		fHistLogJTPta(0),
		fHistJTPta_all(0),
		fHistJTBg(0),
		fHistLogJTBg(0),
		fHistBgMulti(0),
		fHistBgPt(0),
		fHistJetEta(0),
		fHistJetMulti(0),
		fHistJetTracksPt(0),
		fhTrackingEfficiency(0),
		fNpttBins(1),
		fNptaBins(1),
		fJetsCont(0),
		//fJetsConts(0),
		//nJetsConts(0),
		fTracksCont(0),
		fCaloClustersCont(0),
		fTracks(0),
		fTrackArrayName("nonejk"),
		runPeriod(""),
		fEfficiency(0),
		debug(0)


{
	// Default constructor.

	fHistTracksPt       = new TH1*[fNcentBins];
	fHistTracksJt       = new TH1*[fNcentBins];
	fHistClustersPt     = new TH1*[fNcentBins];
	fHistLeadingJetPt   = new TH1*[fNcentBins];
	fHistJetsPt         = new TH1**[fNcentBins];
	fHistBackgroundDone = new TH1**[fNcentBins];
	fHistJTPta          = new TH1***[fNcentBins];
	fHistLogJTPta       = new TH1***[fNcentBins];
	fHistJTPta_all      = new TH1***[fNcentBins];
	fHistJTBg           = new TH1***[fNcentBins];
	fHistLogJTBg        = new TH1***[fNcentBins];
	fHistBgMulti        = new TH1**[fNcentBins];
	fHistBgPt           = new TH1**[fNcentBins];	 
	fHistJetEta	    = new TH1**[fNcentBins];
	fHistJetMulti       = new TH1**[fNcentBins];
	fHistJetTracksPt     = new TH1**[fNcentBins];
	fhTrackingEfficiency = new TProfile*[fNcentBins];
	//CentBinBorders      = new Double_t[10];


	for (Int_t i = 0; i < fNcentBins; i++) {
		fHistJTPta[i] = 0;
		fHistLogJTPta[i] = 0;
		fHistJTPta_all[i] = 0;
		fHistJTBg[i] = 0;
		fHistLogJTBg[i] = 0;
		fHistBackgroundDone[i] = 0;
		fHistTracksPt[i] = 0;
		fHistTracksJt[i] = 0;
		fHistClustersPt[i] = 0;
		fHistLeadingJetPt[i] = 0;
		fHistJetsPt[i] = 0;
		fHistBgMulti[i] = 0;
		fHistBgPt[i] = 0;  
		fHistJetEta[i] = 0; 
		fHistJetMulti[i] = 0;
		fHistJetTracksPt[i] = 0;
		fhTrackingEfficiency[i] = 0;
	}

	/*for(Int_t i = 0; i < nJetsConts; i++){
	  fJetsConts[i] = 0;
	  }*/
	SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskJetJTJT::AliAnalysisTaskJetJTJT(const char *name) : 
	AliAnalysisTaskEmcalJet(name, kTRUE),
	fHistTracksPt(0),
	fHistTracksJt(0),
	fHistClustersPt(0),
	fHistLeadingJetPt(0),
	fHistJetsPt(0),
	fHistBackgroundDone(0),
	fHistJTPta(0),
	fHistLogJTPta(0),
	fHistJTPta_all(0),
	fHistJTBg(0),
	fHistLogJTBg(0),
	fHistBgMulti(0),
	fHistBgPt(0),
	fHistJetEta(0),
	fHistJetMulti(0),
	fHistJetTracksPt(0),
	fhTrackingEfficiency(0),
	fNpttBins(1),
	fNptaBins(1),
	fJetsCont(0),
	//fJetsConts(0),
	//nJetsConts(0),
	fTracksCont(0),
	fCaloClustersCont(0),
	fTracks(0),
	fTrackArrayName("nonejk"),
	runPeriod(""),
	fEfficiency(0),
	debug(0)
{
	// Standard constructor.
	fHistTracksPt       = new TH1*[fNcentBins];
	fHistTracksJt       = new TH1*[fNcentBins];
	fHistClustersPt     = new TH1*[fNcentBins];
	fHistLeadingJetPt   = new TH1*[fNcentBins];
	fHistJetsPt         = new TH1**[fNcentBins];
	fHistBackgroundDone = new TH1**[fNcentBins];	
	fHistJTPta	    = new TH1***[fNcentBins];	
	fHistLogJTPta	    = new TH1***[fNcentBins];	
	fHistJTPta_all	    = new TH1***[fNcentBins];	
	fHistJTBg	    = new TH1***[fNcentBins];	
	fHistLogJTBg	    = new TH1***[fNcentBins];	
	fHistBgMulti        = new TH1**[fNcentBins];
	fHistBgPt           = new TH1**[fNcentBins];
	fHistJetEta         = new TH1**[fNcentBins];
	fHistJetMulti       = new TH1**[fNcentBins];
	fHistJetTracksPt     = new TH1**[fNcentBins];
	fhTrackingEfficiency = new TProfile*[fNcentBins];
	//CentBinBorders	    = new Double_t[10];


	for (Int_t i = 0; i < fNcentBins; i++) {
		fHistJTPta[i] = 0;
		fHistLogJTPta[i] = 0;
		fHistJTPta_all[i] = 0;
		fHistJTBg[i] = 0;
		fHistLogJTBg[i] = 0;
		fHistBackgroundDone[i] = 0;
		fHistTracksPt[i] = 0;
		fHistTracksJt[i] = 0;
		fHistClustersPt[i] = 0;
		fHistLeadingJetPt[i] = 0;
		fHistJetsPt[i] = 0;
		fHistBgMulti[i] = 0;        
		fHistBgPt[i] = 0;           
		fHistJetEta[i] = 0;         
		fHistJetMulti[i] = 0; 
		fHistJetTracksPt[i] = 0; 
		fhTrackingEfficiency[i] = 0;
	}

	/*for(Int_t i = 0; i < nJetsConts; i++){
	  fJetsConts[i] = 0;
	  }*/
	SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskJetJTJT::~AliAnalysisTaskJetJTJT()
{
	// Destructor.
}


void AliAnalysisTaskJetJTJT::setCentBinBorders( int n, Double_t *c){
	fNcentBins=n;  
	if(debug > 0){
		cout << "AliAnalysisTaskJetJTJT::setCentBinBorders: " << endl;
	}
	for(int i= 0 ; i < fNcentBins; i++){
		CentBinBorders[i]= c[i];
		if(debug > 0)
			cout << CentBinBorders[i] << endl;
	}	
}

void AliAnalysisTaskJetJTJT::setTriggPtBorders( int n, Double_t *c){
	fNpttBins=n;  
	if(debug > 0)
		cout << "AliAnalysisTaskJetJTJT::setTriggPtBorders: " << endl;
	for(int i= 0 ; i < fNpttBins; i++){
		TriggPtBorders[i]= c[i];
		if(debug > 0)
			cout << TriggPtBorders[i] << endl;
	}
}

void AliAnalysisTaskJetJTJT::setAssocPtBorders( int n, Double_t *c){
	fNptaBins=n;  
	if(debug > 0)
		cout << "AliAnalysisTaskJetJTJT::setAssocPtBorders: " << endl;
	for(int i= 0 ; i < fNptaBins; i++){
		AssocPtBorders[i]= c[i];
		if(debug > 0)
			cout << AssocPtBorders[i] << endl;
	}
}


//________________________________________________________________________
void AliAnalysisTaskJetJTJT::UserCreateOutputObjects()
{
	// Create user output.

	AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
	if(debug > 0)
		cout << "Creating Histograms" << endl;

	fJetsCont           = GetJetContainer(0);
	/*for(int i = 0; i < nJetsConts ; i++){
	  fJetsConts[0]	    = GetJetContainer(i);
	  }*/
	/*if(fJetsCont) { //get particles and clusters connected to jets
	  fTracksCont       = fJetsCont->GetParticleContainer();
	  fCaloClustersCont = fJetsCont->GetClusterContainer();
	  } else {        //no jets, just analysis tracks and clusters
	  fTracksCont       = GetParticleContainer(0);
	  fCaloClustersCont = GetClusterContainer(0);
	  }*/
	fTracksCont       = GetParticleContainer(0);
	fCaloClustersCont = GetClusterContainer(0);
	fTracksCont->SetClassName("AliVTrack");
	fCaloClustersCont->SetClassName("AliAODCaloCluster");

	TString histname;
	//Int_t fMinBinJt = 0;
	//Int_t fMaxBinJt = 5;

	Int_t NBINSJt=150;
	double LogBinsJt[NBINSJt+1], LimLJt=0.01, LimHJt=10;
	double logBWJt = (TMath::Log(LimHJt)-TMath::Log(LimLJt))/(NBINSJt-1);
	LogBinsJt[0] = 0;
	for(int ij=1;ij<=NBINSJt;ij++) LogBinsJt[ij]=LimLJt*exp(ij*logBWJt);

	int NBINSJtW=150;
	double LimLJtW=TMath::Log(0.01), LimHJtW=TMath::Log(10);

	//==== Efficiency ====
	if(debug > 0)
		cout << "AliAnalysisTaskJetJTJT::UserCreateOutputObjects: Creating efficiency" << endl;
	fEfficiency = new JTJTEfficiency;
	fEfficiency->SetMode( 1 ); // 0:NoEff, 1:Period 2:RunNum 3:Auto
	fEfficiency->SetDataPath("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data"); // Efficiency root file location local or alien



	for (Int_t ic = 0; ic < fNcentBins; ic++) {
		if (fParticleCollArray.GetEntriesFast()>0) {
			histname = "fHistTracksPt_";
			histname += ic;
			fHistTracksPt[ic] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt / 4);
			fHistTracksPt[ic]->GetXaxis()->SetTitle("p_{T,track} (GeV/c)");
			fHistTracksPt[ic]->GetYaxis()->SetTitle("counts");
			fOutput->Add(fHistTracksPt[ic]);
		}

		if (fParticleCollArray.GetEntriesFast()>0) {
			histname = "fHistTracksJt_";
			histname += ic;
			fHistTracksJt[ic] = new TH1F(histname.Data(), histname.Data(), NBINSJt, LogBinsJt);
			fHistTracksJt[ic]->GetXaxis()->SetTitle("J_{T,track} (GeV/c)");
			fHistTracksJt[ic]->GetYaxis()->SetTitle("counts");
			fOutput->Add(fHistTracksJt[ic]);
		}

		if (fParticleCollArray.GetEntriesFast()>0) {
			histname = "fhTrackingEfficiency_";
			histname += ic;
			fhTrackingEfficiency[ic] = new TProfile(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2);
			fhTrackingEfficiency[ic]->GetXaxis()->SetTitle("p_{T,track} (GeV/c)");
			fhTrackingEfficiency[ic]->GetYaxis()->SetTitle("counts");
			fOutput->Add(fhTrackingEfficiency[ic]);
		}
		fHistJTPta[ic] = new TH1**[fNpttBins];
		fHistLogJTPta[ic] = new TH1**[fNpttBins];
		fHistJTPta_all[ic] = new TH1**[fNpttBins];
		fHistJTBg[ic] = new TH1**[fNpttBins];
		fHistLogJTBg[ic] = new TH1**[fNpttBins];
		fHistJetsPt[ic] = new TH1*[fNpttBins];
		fHistBackgroundDone[ic] = new TH1*[fNpttBins];
		fHistBgMulti[ic] = new TH1*[fNpttBins];
		fHistBgPt[ic] = new TH1*[fNpttBins];
		fHistJetEta[ic] = new TH1*[fNpttBins];
		fHistJetMulti[ic] = new TH1*[fNpttBins];
		fHistJetTracksPt[ic] = new TH1*[fNpttBins];
		for(Int_t j=0; j < fNpttBins; j++){
			fHistJTPta[ic][j] = new TH1*[fNptaBins];
			fHistLogJTPta[ic][j] = new TH1*[fNptaBins];
			fHistJTPta_all[ic][j] = new TH1*[fNptaBins];
			fHistJTBg[ic][j] = new TH1*[fNptaBins];
			fHistLogJTBg[ic][j] = new TH1*[fNptaBins];
			for(Int_t k=0; k < fNptaBins; k++){
				fHistJTPta[ic][j][k] = 0;
				fHistLogJTPta[ic][j][k] = 0;
				fHistJTPta_all[ic][j][k] = 0;
				fHistJTBg[ic][j][k] = 0;
				fHistLogJTBg[ic][j][k] = 0;
			}
			fHistJetsPt[ic][j] = 0;
			fHistBackgroundDone[ic][j] = 0;
			fHistBgMulti[ic][j] = 0;
			fHistBgPt[ic][j] = 0;
			fHistJetEta[ic][j] = 0;
			fHistJetMulti[ic][j] =0;
			fHistJetTracksPt[ic][j] = 0;
		}


		if (fParticleCollArray.GetEntriesFast()>0) {
			for(Int_t iptt = 0 ; iptt <  fNpttBins; iptt++){
				for(Int_t ipta = 0 ; ipta < fNptaBins; ipta++){
					histname = "hJTPtaD00C";
					//histname += ic;
					histname += Form("C%02dT%02dA%02d", ic, iptt, ipta);
					if(debug > 1)
						cout << histname << endl;
					fHistJTPta[ic][iptt][ipta] = new TH1F(histname.Data(), histname.Data(),NBINSJt, LogBinsJt);
					fHistJTPta[ic][iptt][ipta]->GetXaxis()->SetTitle("J_{T,track} (GeV/c)");
					fHistJTPta[ic][iptt][ipta]->GetYaxis()->SetTitle("counts");
					fOutput->Add(fHistJTPta[ic][iptt][ipta]);

					histname = "hLogJTPtaD00C";
					//histname += ic;
					histname += Form("C%02dT%02dA%02d", ic, iptt, ipta);
					if(debug > 1)
						cout << histname << endl;
					fHistLogJTPta[ic][iptt][ipta] = new TH1F(histname.Data(), histname.Data(), NBINSJtW, LimLJtW, LimHJtW);
					fHistLogJTPta[ic][iptt][ipta]->GetXaxis()->SetTitle("ln(J_{T,track}/ (GeV/c))");
					fHistLogJTPta[ic][iptt][ipta]->GetYaxis()->SetTitle("counts");
					fOutput->Add(fHistLogJTPta[ic][iptt][ipta]);

					histname = "hJTPta_allD00C";
					//histname += ic;
					histname += Form("C%02dT%02dA%02d", ic, iptt, ipta);
					if(debug > 1)
						cout << histname << endl;
					fHistJTPta_all[ic][iptt][ipta] = new TH1F(histname.Data(), histname.Data(), NBINSJt, LogBinsJt);
					fHistJTPta_all[ic][iptt][ipta]->GetXaxis()->SetTitle("J_{T,track} (GeV/c)");
					fHistJTPta_all[ic][iptt][ipta]->GetYaxis()->SetTitle("counts");
					fOutput->Add(fHistJTPta_all[ic][iptt][ipta]);

					histname = "hJTBg";
					//histname += ic;
					histname += Form("C%02dT%02dA%02d", ic, iptt, ipta);
					if(debug > 1)
						cout << histname << endl;
					fHistJTBg[ic][iptt][ipta] = new TH1F(histname.Data(), histname.Data(), NBINSJt, LogBinsJt);
					fHistJTBg[ic][iptt][ipta]->GetXaxis()->SetTitle("J_{T,track} (GeV/c)");
					fHistJTBg[ic][iptt][ipta]->GetYaxis()->SetTitle("counts");
					fOutput->Add(fHistJTBg[ic][iptt][ipta]);

					histname = "hLogJTBg";
					//histname += ic;
					histname += Form("C%02dT%02dA%02d", ic, iptt, ipta);
					if(debug > 1)
						cout << histname << endl;
					fHistLogJTBg[ic][iptt][ipta] = new TH1F(histname.Data(), histname.Data(), NBINSJtW, LimLJtW, LimHJtW);
					fHistLogJTBg[ic][iptt][ipta]->GetXaxis()->SetTitle("ln(J_{T,track}/ (GeV/c))");
					fHistLogJTBg[ic][iptt][ipta]->GetYaxis()->SetTitle("counts");
					fOutput->Add(fHistLogJTBg[ic][iptt][ipta]);


				}
			}
		}

		if (fClusterCollArray.GetEntriesFast()>0) {
			histname = "fHistClustersPt_";
			histname += ic;
			fHistClustersPt[ic] = new TH1F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2);
			fHistClustersPt[ic]->GetXaxis()->SetTitle("p_{T,clus} (GeV/c)");
			fHistClustersPt[ic]->GetYaxis()->SetTitle("counts");
			fOutput->Add(fHistClustersPt[ic]);
		}

		if (fJetCollArray.GetEntriesFast()>0) {
			histname = "fHistLeadingJetPt_";
			histname += ic;
			fHistLeadingJetPt[ic] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
			fHistLeadingJetPt[ic]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
			fHistLeadingJetPt[ic]->GetYaxis()->SetTitle("counts");
			fOutput->Add(fHistLeadingJetPt[ic]);

			for(Int_t iptt = 0 ; iptt <  fNpttBins; iptt++){

				histname = "fHistJetsPt_";
				histname += Form("C%02dT%02d", ic, iptt);
				fHistJetsPt[ic][iptt] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
				fHistJetsPt[ic][iptt]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
				fHistJetsPt[ic][iptt]->GetYaxis()->SetTitle("counts");
				fOutput->Add(fHistJetsPt[ic][iptt]);

				histname = "fHistBackgroundDone_";
				histname += Form("C%02dT%02d", ic, iptt);;
				fHistBackgroundDone[ic][iptt] = new TH1F(histname.Data(), histname.Data(), 2, -1, 2);
				fHistBackgroundDone[ic][iptt]->GetXaxis()->SetTitle("Number of jets");
				fHistBackgroundDone[ic][iptt]->GetYaxis()->SetTitle("0 = not done, 1 = done");
				fOutput->Add(fHistBackgroundDone[ic][iptt]);

				histname = "fHistJetEta_";
				histname += Form("C%02dT%02d", ic, iptt);
				fHistJetEta[ic][iptt] = new TH1F(histname.Data(), histname.Data(), fNbins, -5, 5);
				fHistJetEta[ic][iptt]->GetXaxis()->SetTitle("#eta");
				fHistJetEta[ic][iptt]->GetYaxis()->SetTitle("jets");
				fOutput->Add(fHistJetEta[ic][iptt]);

				histname = "fHistJetMulti_";
				histname += Form("C%02dT%02d", ic, iptt);
				fHistJetMulti[ic][iptt] = new TH1F(histname.Data(), histname.Data(), 200, 0, 200);
				fHistJetMulti[ic][iptt]->GetXaxis()->SetTitle("Multiplicity");
				fHistJetMulti[ic][iptt]->GetYaxis()->SetTitle("jets");
				fOutput->Add(fHistJetMulti[ic][iptt]);

				histname = "fHistBgMulti_";
				histname += Form("C%02dT%02d", ic, iptt);
				fHistBgMulti[ic][iptt] = new TH1F(histname.Data(), histname.Data(), 200, 0, 200);
				fHistBgMulti[ic][iptt]->GetXaxis()->SetTitle("Multiplicity");
				fHistBgMulti[ic][iptt]->GetYaxis()->SetTitle("Events");
				fOutput->Add(fHistBgMulti[ic][iptt]);

				histname = "fHistBgPt_";
				histname += Form("C%02dT%02d", ic, iptt);
				fHistBgPt[ic][iptt] = new TH1F(histname.Data(), histname.Data(), fNbins,fMinBinPt, fMaxBinPt/4);
				fHistBgPt[ic][iptt]->GetXaxis()->SetTitle("Multiplicity");
				fHistBgPt[ic][iptt]->GetYaxis()->SetTitle("jets");
				fOutput->Add(fHistBgPt[ic][iptt]);

				histname = "fHistJetTracksPt_";
				//histname += ic;
				histname += Form("C%02dT%02d", ic, iptt);
				if(debug > 1)
					cout << histname << endl;
				fHistJetTracksPt[ic][iptt] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt/4);
				fHistJetTracksPt[ic][iptt]->GetXaxis()->SetTitle("p_{T,track} (GeV/c)");
				fHistJetTracksPt[ic][iptt]->GetYaxis()->SetTitle("counts");
				fOutput->Add(fHistJetTracksPt[ic][iptt]);
			}

			/*
			   if (!(GetJetContainer()->GetRhoName().IsNull())) {
			   histname = "fHistJetsCorrPtArea_";
			   histname += i;
			   fHistJetsCorrPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins*2, -fMaxBinPt, fMaxBinPt, 30, 0, 3);
			   fHistJetsCorrPtArea[i]->GetXaxis()->SetTitle("p_{T}^{corr} [GeV/c]");
			   fHistJetsCorrPtArea[i]->GetYaxis()->SetTitle("area");
			   fOutput->Add(fHistJetsCorrPtArea[i]);
			   }
			   */
		}
	}

	PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetJTJT::FillHistograms()
{
	// Fill histograms.

	AliCentrality *aliCent = InputEvent()->GetCentrality();
	fCentBin = 0;
	if (aliCent) {
		//fCent = aliCent->GetCentralityPercentile(fCentEst.Data());
		fCent = aliCent->GetCentralityPercentile("V0M");
		/*if(debug > 0){
		  cout << "Centrality " << fCent << endl;
		  }*/
		for(int ic = 0; ic < fNcentBins; ic++){
			/*if(debug > 0){
			  cout << "ic: " << ic << " / " << fNcentBins << endl;
			  cout << "Centrality bin " << fCentBin << endl;
			  cout << "Border: " << CentBinBorders[ic] << endl;
			  } */
			if(fCent > CentBinBorders[ic]){
				fCentBin = ic;
			}
		}
		//cout << "Centrality bin: " << fCentBin << endl;
	} else {
		AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
		fCentBin = 3;
	}
	int fHadronSelectionCut = 5; //5=Hybrid cut
	if (fTracksCont) {
		AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle(0)); 
		while(track) {
			double ptt = track->Pt();

			//<<<<<<<<<<<< Efficiency >>>>>>>>>>>
			//AliJBaseTrack *triggTr = (AliJBaseTrack*)fchargedHadronList->At(i);

			double effCorr = 1./fEfficiency->GetCorrection(ptt, fHadronSelectionCut, fCent);  // here you generate warning if ptt>30
			//double effCorr = 1.;
			fhTrackingEfficiency[fCentBin]->Fill( ptt, 1./effCorr );
			//triggTr->SetTrackEff( 1./effCorr );
			//<<<<<<<<<<<< Efficiency >>>>>>>>>>>

			if(ptt > 0 && 1.0/ptt > 0){
				fHistTracksPt[fCentBin]->Fill(ptt,1./ptt*effCorr); 
			}


			track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
		}
	}

	if (fCaloClustersCont) {
		AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(0); 
		while(cluster) {
			TLorentzVector nPart;
			cluster->GetMomentum(nPart, fVertex);
			fHistClustersPt[fCentBin]->Fill(nPart.Pt());

			cluster = fCaloClustersCont->GetNextAcceptCluster();
		}
	}

	Int_t fPttBin, fPtaBin;
	fPtaBin = 0;


	if (fJetsCont) {
		//Int_t Njets = fJetsCont->GetNJets();
		Int_t Njets = 0;
		if(debug > 1){
			cout << "Number of Jets: " << Njets << endl;
		}

		//Make a array to hold jets to be tested in background jT
		Float_t jetPhis[200] = {};
		Float_t jetEtas[200] = {};
		AliEmcalJet *jet = fJetsCont->GetNextAcceptJet(0);
		Int_t ij = 0;
		while(jet) {
			//cout << "Jet found " << ij << " pt: " << jet->Pt() << endl;
			if(jet->Pt() > 5){    //Only consider jets with pT > 5 GeV
				jetPhis[ij] = jet->Phi();
				jetEtas[ij] = jet->Eta();
				ij++;
				Njets++;
				if(debug > 1)
					cout << "i: " << ij << " jetPhi: " << jetPhis[ij] << " jetEta: " << jetEtas[ij] << endl;
			}else{
				//jetPhis[ij] = 100;
				//jetEtas[ij] = 100;
				//Njets--;
				if(debug > 1)
					cout << "jetPt: " << jet->Pt() << " jetPhi: " << jet->Phi() << " jetEta: " << jet->Eta() << endl;
			}
			//i++; 
			jet = fJetsCont->GetNextAcceptJet();
		}


		//fTracks =dynamic_cast <TClonesArray*>(InputEvent()->FindListObject( fTrackArrayName.Data() ));
		jet = fJetsCont->GetNextAcceptJet(0); 
		while(jet) {
			if(jet->Pt() > 5){
				if(jet->Eta() < -0.4 || jet->Eta() > 0.4){
					if(debug > 0)
						cout << "Jet outside eta range, Eta: " << jet->Eta() << endl;
					jet = fJetsCont->GetNextAcceptJet();
					continue;
				}
				//cout << "Jet found " << ij << " pt: " << jet->Pt() << endl;
				//Get the trigger pT bin
				fPttBin = 0;
				for(int iptt = 0 ; iptt < fNpttBins; iptt++){
					if(jet->Pt() > TriggPtBorders[iptt]){
						fPttBin = iptt;
					}

				}
				fHistJetEta[fCentBin][fPttBin]->Fill(jet->Eta());
				if(jet->Pt() > 0 && 1.0/jet->Pt() > 0){
					fHistJetsPt[fCentBin][fPttBin]->Fill(jet->Pt(),1.0/jet->Pt());  //Fill jet dN/(pT dpT)
				}

				/*if (fHistJetsCorrPtArea[fCentBin]) {
				  Float_t corrPt = jet->Pt() - fJetsCont->GetRhoVal() * jet->Area();
				  fHistJetsCorrPtArea[fCentBin]->Fill(corrPt, jet->Area());
				  }*/
				//Float_t jetp = sqrt(jet->Px()*jet->Px()+jet->Py()*jet->Py()+jet->Pz()*jet->Pz()); //Jet pT norm


				Int_t nTrack = jet->GetNumberOfTracks();
				if (debug > 0)			
					cout << "Number of tracks " << nTrack << " Jet Pt: " << jet->Pt() << endl;
				fHistJetMulti[fCentBin][fPttBin]->Fill(nTrack);
				for(Int_t it = 0; it < nTrack; it++ ){
					AliVParticle *track = (AliVParticle*)jet->TrackAt( it, fTracks );
					if( !track ){
						cout << "No Track found" << endl;
						continue;
					}
					fPtaBin = 0; //Get the associated pT bin
					for(int ipta = 0 ; ipta < fNptaBins; ipta++){
						if(track->Pt() > AssocPtBorders[ipta]){
							fPtaBin = ipta;
						}
					}
					fHistJetTracksPt[fCentBin][fPttBin]->Fill(track->Pt());
					if(debug > 2)
						cout << "Filling fHistJetTracksPt C" << fCentBin << " T" << fPttBin << endl;
					/*Float_t dotproduct = track->Px()*jet->Px()+track->Py()*jet->Py()+track->Pz()*jet->Pz(); // p_track dot p_jet
					Float_t constp = sqrt(track->Px()*track->Px()+track->Py()*track->Py()+track->Pz()*track->Pz()); // track pT norm
					Float_t normproduct = constp*jetp; // jet pT norm times track pT norm
					Float_t costheta2 = dotproduct/normproduct; 
					//Float_t sintheta = sqrt(1-costheta2*costheta2);
					Float_t jt = constp*sqrt(1-costheta2*costheta2);*/
					Float_t jt = getJt(track,jet,0);
					double effCorr = 1./fEfficiency->GetCorrection(track->Pt(), fHadronSelectionCut, fCent);  // here you generate warning if ptt>30
					//double effCorr = 1.;
					if(jt > 0 && 1.0/jt > 0){
						fHistTracksJt[fCentBin]->Fill(jt,1.0/jt*effCorr); //Fill dN/(djT jT)
						fHistJTPta[fCentBin][fPttBin][fPtaBin]->Fill(jt,1.0/jt*effCorr); //Fill dN/(djT jT)
						fHistLogJTPta[fCentBin][fPttBin][fPtaBin]->Fill(TMath::Log(jt),1.0/jt*effCorr); //Fill logarithmic dN/(djT jT)
					}
					if(debug > 1)
						cout << "Filling JT C" << fCentBin << "T" <<  fPttBin << "A" << fPtaBin << " jt:" << jt << " with " << 1.0/jt*effCorr << endl;
				}

				//Get Jet azimuth and rapidity of jet
				Float_t jetAngle = jet->Phi();
				Float_t jetRap = jet->Eta();

				//Rotate jet angle for background cone
				Float_t rotatedAngle = jetAngle+TMath::Pi()/2;
				if(rotatedAngle > TMath::Pi()*2){
					rotatedAngle = rotatedAngle- TMath::Pi()*2;
				}
				Float_t jetArea = jet->Area();
				Float_t testRadius = TMath::Sqrt(jetArea/TMath::Pi());

				Bool_t doBg = 1;

				//Test if there are jets in the background test cone
				for(int i_j = 0; i_j < Njets; i_j++){
					Float_t diffR = TMath::Sqrt(TMath::Power(jetPhis[i_j]-rotatedAngle,2)+TMath::Power(jetEtas[i_j]-jetRap,2));
					if(debug > 1){
						cout << "i_j: " << i_j << " JetPhi: " << jetPhis[i_j] << " jetEta: " << jetEtas[i_j] << endl;
						cout << "DiffR: " << diffR << " doBG: " << doBg <<endl;
					}
					if(diffR < testRadius *2){ //Jets muts be at least 2*cone radius away from the background cone axis
						doBg =0;
						break;
					}

				}

				// Do jT for all particles in respect to jet axis
				if (fTracksCont) {
					int counter = 0;
					AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle(0)); 
					while(track) {
						/*Float_t dotproduct = track->Px()*jet->Px()+track->Py()*jet->Py()+track->Pz()*jet->Pz();
						  Float_t constp = sqrt(track->Px()*track->Px()+track->Py()*track->Py()+track->Pz()*track->Pz());
						  Float_t normproduct = constp*jetp;
						  Float_t costheta2 = dotproduct/normproduct;
						//Float_t sintheta = sqrt(1-costheta2*costheta2);
						Float_t jt = constp*sqrt(1-costheta2*costheta2);*/
						Double_t jt = getJt(track,jet,0);
						double effCorr = 1./fEfficiency->GetCorrection(track->Pt(), fHadronSelectionCut, fCent);  // here you generate warning if ptt>30
						if(jt > 0 && 1.0/jt > 0){
							fHistJTPta_all[fCentBin][fPttBin][fPtaBin]->Fill(jt,1.0/jt*effCorr);
						}
						for(int ipta = 0 ; ipta < fNptaBins; ipta++){
							if(track->Pt() > AssocPtBorders[ipta]){
								fPtaBin = ipta;
							}
						}

						//If background is to be filled
						if(doBg){
							Float_t diffR = TMath::Sqrt(TMath::Power(track->Phi()-rotatedAngle,2)+TMath::Power(track->Eta()-jetRap,2));
							//Particles in the rotated cone
							if(diffR < testRadius){
								counter++;
								fHistBgPt[fCentBin][fPttBin]->Fill(track->Pt());
								/*dotproduct = -track->Px()*jet->Py()+track->Py()*jet->Px()+track->Pz()*jet->Pz();
								  constp = sqrt(track->Px()*track->Px()+track->Py()*track->Py()+track->Pz()*track->Pz());
								  normproduct = constp*jetp;
								  costheta2 = dotproduct/normproduct;
								//sintheta = sqrt(1-costheta2*costheta2);
								jt = constp*sqrt(1-costheta2*costheta2);*/
								jt = getJt(track,jet,1);
								if(jt > 0 && 1.0/jt > 0){
									fHistJTBg[fCentBin][fPttBin][fPtaBin]->Fill(jt,1.0/jt*effCorr);
									fHistLogJTBg[fCentBin][fPttBin][fPtaBin]->Fill(TMath::Log(jt),1.0/jt*effCorr);
								}
								if(debug > 1)
									cout << "Filling Background C" << fCentBin << "T" <<  fPttBin << "A" << fPtaBin << " jt:" << jt << " with " << 1.0/jt*effCorr << endl;
								//Fill background jT

							}
						}
						track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
					}
					if(doBg){
						fHistBgMulti[fCentBin][fPttBin]->Fill(counter);
					}
				}
				if(doBg){
					fHistBackgroundDone[fCentBin][fPttBin]->Fill(1);
				}else{
					fHistBackgroundDone[fCentBin][fPttBin]->Fill(0);
				}

			}
			jet = fJetsCont->GetNextAcceptJet(); 

		}
		jet = fJetsCont->GetLeadingJet();
		if(jet){
			if(jet->Pt() > 0 && 1.0/jet->Pt() > 0){
				fHistLeadingJetPt[fCentBin]->Fill(jet->Pt(),1.0/jet->Pt());
			}
		}
	}

	CheckClusTrackMatching();

	return kTRUE;
}


//-----------------------------------------------------------------------
Double_t AliAnalysisTaskJetJTJT::getJt(AliVTrack *track, AliEmcalJet *jet,int reverse){
	Float_t dotproduct = 0;
	Float_t jetp = sqrt(jet->Px()*jet->Px()+jet->Py()*jet->Py()+jet->Pz()*jet->Pz()); //Jet pT norm
	if(reverse){
		dotproduct = -track->Px()*jet->Py()+track->Py()*jet->Px()+track->Pz()*jet->Pz();
	} else{
		dotproduct = track->Px()*jet->Px()+track->Py()*jet->Py()+track->Pz()*jet->Pz();
	}
	Float_t constp = sqrt(track->Px()*track->Px()+track->Py()*track->Py()+track->Pz()*track->Pz());
	Float_t normproduct = constp*jetp;
	Float_t costheta2 = dotproduct/normproduct;
	//Float_t sintheta = sqrt(1-costheta2*costheta2);
	Float_t jt = constp*sqrt(1-costheta2*costheta2);
	return jt;
}

Double_t AliAnalysisTaskJetJTJT::getJt(AliVParticle *track, AliEmcalJet *jet,int reverse){
	Float_t dotproduct = 0;
	Float_t jetp = sqrt(jet->Px()*jet->Px()+jet->Py()*jet->Py()+jet->Pz()*jet->Pz()); //Jet pT norm
	if(reverse){
		dotproduct = -track->Px()*jet->Py()+track->Py()*jet->Px()+track->Pz()*jet->Pz();
	} else{
		dotproduct = track->Px()*jet->Px()+track->Py()*jet->Py()+track->Pz()*jet->Pz();
	}
	Float_t constp = sqrt(track->Px()*track->Px()+track->Py()*track->Py()+track->Pz()*track->Pz());
	Float_t normproduct = constp*jetp;
	Float_t costheta2 = dotproduct/normproduct;
	//Float_t sintheta = sqrt(1-costheta2*costheta2);
	Float_t jt = constp*sqrt(1-costheta2*costheta2);
	return jt;
}

//________________________________________________________________________
void AliAnalysisTaskJetJTJT::CheckClusTrackMatching()
{

	if(!fTracksCont || !fCaloClustersCont)
		return;

	Double_t deta = 999;
	Double_t dphi = 999;

	//Get closest cluster to track
	AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle(0)); 
	while(track) {
		//Get matched cluster
		Int_t emc1 = track->GetEMCALcluster();
		if(fCaloClustersCont && emc1>=0) {
			AliVCluster *clusMatch = fCaloClustersCont->GetCluster(emc1);
			if(clusMatch) {
				AliPicoTrack::GetEtaPhiDiff(track, clusMatch, dphi, deta);
				//fHistPtDEtaDPhiTrackClus->Fill(track->Pt(),deta,dphi);
			}
		}
		track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
	}

	//Get closest track to cluster
	AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(0); 
	while(cluster) {
		TLorentzVector nPart;
		cluster->GetMomentum(nPart, fVertex);
		fHistClustersPt[fCentBin]->Fill(nPart.Pt());

		//Get matched track
		AliVTrack *mt = NULL;      
		AliAODCaloCluster *acl = dynamic_cast<AliAODCaloCluster*>(cluster);
		if(acl) {
			if(acl->GetNTracksMatched()>1)
				mt = static_cast<AliVTrack*>(acl->GetTrackMatched(0));
		}
		else {
			AliESDCaloCluster *ecl = dynamic_cast<AliESDCaloCluster*>(cluster);
			Int_t im = ecl->GetTrackMatchedIndex();
			if(fTracksCont && im>=0) {
				mt = static_cast<AliVTrack*>(fTracksCont->GetParticle(im));
			}
		}
		if(mt) {
			AliPicoTrack::GetEtaPhiDiff(mt, cluster, dphi, deta);
			//fHistPtDEtaDPhiClusTrack->Fill(nPart.Pt(),deta,dphi);

			//debugging
			/*
			   if(mt->IsEMCAL()) {
			   Int_t emc1 = mt->GetEMCALcluster();
			   Printf("current id: %d  emc1: %d",fCaloClustersCont->GetCurrentID(),emc1);
			   AliVCluster *clm = fCaloClustersCont->GetCluster(emc1);
			   AliPicoTrack::GetEtaPhiDiff(mt, clm, dphi, deta);
			   Printf("deta: %f dphi: %f",deta,dphi);
			   }
			   */    
		}
		cluster = fCaloClustersCont->GetNextAcceptCluster();
	}
}

//________________________________________________________________________
/*AliJetContainer* AliAnalysisTaskJetJTJT::AddJetContainer(const char *n, TString defaultCutType, Float_t jetRadius) {

  AliAnalysisTaskEmcalJet::ExecOnce();
  nJetsConts++;
  AliJetContainer *cont = 0x0;
  cont = AliAnalysisTaskEmcalJet::AddJetContainer(n,defaultCutType,jetRadius);
  return cont;
  }*/


//________________________________________________________________________
void AliAnalysisTaskJetJTJT::ExecOnce() {

	if(debug > 0){
	cout << "AliAnalysisTaskJetJTJT::ExecOnce(): " << endl;
	cout << "Get fTracks from " << fTrackArrayName.Data() << endl;
	}
	fTracks =dynamic_cast <TClonesArray*>(InputEvent()->FindListObject( fTrackArrayName.Data() ));

	AliAnalysisTaskEmcalJet::ExecOnce();
	if(debug > 1)
		cout << "Efficiency: Set Run Period Name " << runPeriod << endl;
	fEfficiency->SetPeriodName(runPeriod);
	if(debug > 1)
		cout << "Efficiency: Set Run number " << InputEvent()->GetRunNumber() << endl;
	fEfficiency->SetRunNumber( InputEvent()->GetRunNumber() ); //TODO Get run Number
	if(debug > 1)
		cout << "Efficiency: Load()" << endl;
	fEfficiency->Load();
	if(debug > 1)
		cout << "fEfficiency loaded" << endl;



	if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
	if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
	if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetJTJT::Run()
{
	// Run analysis code here, if needed. It will be executed before FillHistograms().

	return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliAnalysisTaskJetJTJT::Terminate(Option_t *) 
{
	// Called once at the end of the analysis.
}



//________________________________________________________________________
JTJTEfficiency::JTJTEfficiency():
	fMode(kAuto),
	fPeriod(-1),
	fDataPath(""),
	fName(""),
	fPeriodStr(""),
	fMCPeriodStr(""),
	fRunNumber(0),
	fTag(""),
	fInputRootName(""),
	fInputRoot(NULL),
	fCentBinAxis(0x0)
{

}

JTJTEfficiency::JTJTEfficiency(const JTJTEfficiency& obj) :
  fMode(obj.fMode),
  fPeriod(obj.fPeriod),
  fDataPath(obj.fDataPath),
  fName(obj.fName),
  fPeriodStr(obj.fPeriodStr),
  fMCPeriodStr(obj.fMCPeriodStr),
  fRunNumber(obj.fRunNumber),
  fTag(obj.fTag),
  fInputRootName(obj.fInputRootName),
  fInputRoot(obj.fInputRoot),
  fCentBinAxis(obj.fCentBinAxis)
{
	// copy constructor TODO: handling of pointer members
	JUNUSED(obj);
}

JTJTEfficiency& JTJTEfficiency::operator=(const JTJTEfficiency& obj){
	// equal sign operator TODO: content
	JUNUSED(obj);
	return *this;
}


//________________________________________________________________________
double JTJTEfficiency::GetCorrection( double pt, int icut , double cent ) const {
	if( fMode == kNotUse ) return 1;
	int icent = fCentBinAxis->FindBin( cent ) -1 ;
	if( icent < 0 || icent > fCentBinAxis->GetNbins()-1 ) {
		cout<<"J_WARNING : Centrality "<<cent<<" is out of CentBinBorder"<<endl;
		return 1;
	}
	// TODO error check for icent;
	int ivtx = 0;
	if( ! fCorrection[ivtx][icent][icut] ) {
		cout<<"J_WARNING : No Eff Info "<<pt<<"\t"<<icut<<"\t"<<cent<<"\t"<<icent<<endl;
		return 1;
	}
	TGraphErrors * gr = fCorrection[ivtx][icent][icut];
	//=== TEMPERORY SETTING. IT will be removed soon.
	if( pt > 30 ) pt = 30; // Getting eff of 30GeV for lager pt
	double cor = gr->Eval(pt);
	if ( cor < 0.2 ) cor = 0.2;
	return cor;
}


TString JTJTEfficiency::GetEffName() {
	/*
	   1. kNotUse : no Load, efficiency is 1 always
	   2. has fInputRootName : Load that or crash
	   3. has fName : Load fName [+runnumber] or crash
	   4. has runnumber : Find Good MC period from AliJRunTable, or crash
	   3. has period : Find Good MC period from AliJRunTable, or crash

	   }
	   */
if(fPeriodStr == "LHC10b"){
	fInputRootName = "Eff--LHC10b-LHC10d1-0-.root";
}
if(fPeriodStr == "LHC10c"){
	fInputRootName = "Eff--LHC10c-LHC10d4-0-.root";
}
if(fPeriodStr == "LHC10d"){
	fInputRootName = "Eff--LHC10d-LHC10f6a-0-.root";
}
if(fPeriodStr == "LHC10e"){
	fInputRootName = "Eff--LHC10e-LHC10e20-0-.root";
}
if(fPeriodStr == "LHC10h"){
	fInputRootName = "Eff--LHC10h-LHC11a10a_bis-0-.root";
}
if(fPeriodStr == "LHC11a"){
	fInputRootName = "Eff--LHC11a-LHC11b10a-0-.root";
}
if(fPeriodStr == "LHC13b"){
	fInputRootName = "Eff--LHC13b-LHC13b2-efix_p1-0-.root";
}

if(fPeriodStr == "LHC13c"){
	fInputRootName = "Eff--LHC13c-LHC13b2-efix_p1-0-.root";
}
if(fPeriodStr == "LHC13d"){
	fInputRootName = "Eff--LHC13d-LHC13b2-efix_p1-0-.root";
}
if(fPeriodStr == "LHC13e"){
	fInputRootName = "Eff--LHC13e-LHC13b2-efix_p1-0-.root";
}

return fInputRootName;
}

TString JTJTEfficiency::GetEffFullName() {
	GetEffName();
	fInputRootName = fDataPath + "/" + fInputRootName;
	return fInputRootName;
}


//________________________________________________________________________
bool JTJTEfficiency::Load(){
	// Load Efficiency File based on fMode
	if( fMode == kNotUse ) {
		cout<<"J_WARNING : Eff Mode is \"NOTUSE\". eff is 1 !!!"<<endl;
		return true;
	}
	GetEffFullName();
	if (TString(fInputRootName).BeginsWith("alien:"))  TGrid::Connect("alien:");
	fInputRoot = TFile::Open( fInputRootName);
	//fInputRoot = new TFile( fInputRootName,"READ");
	if( !fInputRoot ) {
		cout<< "J_ERROR : "<<fInputRootName <<" does not exist"<<endl;
		gSystem->Exit(1);
	}

	//fEffDir[0] = (TDirectory*)fInputRoot->Get("EffRE");
	///fEffDir[1] = (TDirectory*)fInputRoot->Get("EffMC");
	fEffDir[2] = (TDirectory*)fInputRoot->Get("Efficiency");
	//iif( fEffDir[0] && fEffDir[1] && fEffDir[2] )
	if( !fEffDir[2] )
	{
		cout<< "J_ERROR : Directory EFF is not exist"<<endl;
		gSystem->Exit(1);
	}

	fCentBinAxis = (TAxis*)fEffDir[2]->Get("CentralityBin");
	if( !fCentBinAxis ){
		cout<< "J_ERROR : No CentralityBin in directory"<<endl;
		gSystem->Exit(1);
	}


	int nVtx = 1;
	int nCentBin = fCentBinAxis->GetNbins();
	for( int ivtx=0;ivtx<nVtx;ivtx++ ){
		for( int icent=0;icent<nCentBin;icent++ ){
			for( int icut=0;icut<kJNTrackCuts;icut++ ){
				fCorrection[ivtx][icent][icut]
					= (TGraphErrors*) fEffDir[2]->Get(Form("gCor%02d%02d%02d", ivtx,icent,icut));
				//cout<<"J_LOG : Eff graph - "<<Form("gCor%02d%02d%02d", ivtx,icent,icut)<<" - "<<g<<endl;
			}
		}
	}
	cout<<"J_LOG : Eff file is "<<fInputRootName<<endl;
	cout<<"J_LOG : Eff Cent Bins are ";
	for( int i=0;i<=nCentBin;i++ ){
		//cout<<fCentBinAxis->GetXbins()->At(i)<<" ";
	}
	//cout<<endl;
	return true;
}



