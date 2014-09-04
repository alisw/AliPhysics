/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *					                           					          *
 * Authors: Friederike Bock, Baldo Sahlmueller 		                      *
 * Version 1.0                        	       							  *
 *                           						 					  *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is    	  *
 * provided "as is" without express or implied warranty.       			  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class handling all kinds of selection cuts for
// Photon from EMCAL clusters
//---------------------------------------------
////////////////////////////////////////////////

#include "AliCaloPhotonCuts.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODHandler.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "AliStack.h"
#include "AliAODConversionMother.h"
#include "AliAODConversionPhoton.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "TList.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliV0ReaderV1.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliPicoTrack.h"
#include "AliEMCALRecoUtils.h"
#include "AliTrackerBase.h"

class iostream;

using namespace std;

ClassImp(AliCaloPhotonCuts)


const char* AliCaloPhotonCuts::fgkCutNames[AliCaloPhotonCuts::kNCuts] = {
	"ClusterType",          //0   0: all,    1: EMCAL,   2: PHOS
	"EtaMin",               //1   0: -10,    1: -0.6687, 2: -0,5, 3: -2
	"EtaMax",               //2   0: 10,     1: 0.66465, 2: 0.5,  3: 2
	"PhiMin",               //3   0: -10000, 1: 1.39626
	"PhiMax",               //4   0: 10000, 1: 3.125
	"DistanceToBadChannel",	//5   0: 0,      1: 5
	"Timing",               //6   0: no cut
	"TrackMatching",        //7   0: 0,      1: 5
	"ExoticCell",           //8   0: no cut
	"MinEnergy",            //9   0: no cut, 1: 0.05,    2: 0.1,  3: 0.15, 4: 0.2, 5: 0.3, 6: 0.5, 7: 0.75, 8: 1, 9: 1.25 (all GeV)
	"MinNCells",            //10  0: no cut, 1: 1,       2: 2,    3: 3,    4: 4,   5: 5,   6: 6
	"MinM02",				//11
	"MaxM02",				//12
	"MinM20",				//13
	"MaxM20",				//14
	"MaximumDispersion",	//15
	"NLM"					//16
};


//________________________________________________________________________
AliCaloPhotonCuts::AliCaloPhotonCuts(const char *name,const char *title) :
	AliAnalysisCuts(name,title),
	fHistograms(NULL),	
	fClusterType(0),
	fMinEtaCut(-10),
	fMaxEtaCut(10),
	fUseEtaCut(0),
	fMinPhiCut(-10000),
	fMaxPhiCut(-10000),
	fUsePhiCut(0),
	fMinDistanceToBadChannel(0),
	fUseDistanceToBadChannel(0),
	fMaxTimeDiff(10e10),
	fUseTimeDiff(0),
	fMinDistTrackToCluster(0),
	fUseDistTrackToCluster(0),
	fExoticCell(0),
	fUseExoticCell(0),
	fMinEnergy(0),
	fUseMinEnergy(0),
	fMinNCells(0),
	fUseNCells(0),
	fMaxM02(1000),
	fMinM02(0),
	fUseM02(0),
	fMaxM20(1000),
	fMinM20(0),
	fUseM20(0),
	fMaxDispersion(1000),
	fUseDispersion(0),
	fMinNLM(0),
	fMaxNLM(1000),
	fUseNLM(0),
	fCutString(NULL),
	fHistCutIndex(NULL),
	fHistAcceptanceCuts(NULL),
	fHistClusterIdentificationCuts(NULL),
	fHistClusterEtavsPhiBeforeAcc(NULL),
	fHistClusterEtavsPhiAfterAcc(NULL),
	fHistClusterEtavsPhiAfterQA(NULL),
	fHistDistanceToBadChannelBeforeAcc(NULL),
	fHistDistanceToBadChannelAfterAcc(NULL),
    fHistClusterRBeforeQA(NULL),
    fHistClusterRAfterQA(NULL),
    fHistClusterdEtadPhiBeforeQA(NULL),
    fHistClusterdEtadPhiAfterQA(NULL),
	fHistClusterTimevsEBeforeQA(NULL),
	fHistClusterTimevsEAfterQA(NULL),
	fHistExoticCellBeforeQA(NULL),
	fHistExoticCellAfterQA(NULL),
	fHistNMatchedTracks(NULL),
	fHistDistanceTrackToClusterBeforeQA(NULL),
	fHistDistanceTrackToClusterAfterQA(NULL),
	fHistEnergyOfClusterBeforeQA(NULL),
	fHistEnergyOfClusterAfterQA(NULL),
	fHistNCellsBeforeQA(NULL),
	fHistNCellsAfterQA(NULL),
	fHistM02BeforeQA(NULL),
	fHistM02AfterQA(NULL),
	fHistM20BeforeQA(NULL),
	fHistM20AfterQA(NULL),
	fHistDispersionBeforeQA(NULL),
	fHistDispersionAfterQA(NULL),
	fHistNLMBeforeQA(NULL),
	fHistNLMAfterQA(NULL)
{
   for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
   fCutString=new TObjString((GetCutNumber()).Data());
}

//________________________________________________________________________
AliCaloPhotonCuts::AliCaloPhotonCuts(const AliCaloPhotonCuts &ref) :
   AliAnalysisCuts(ref),
 	fHistograms(NULL),	
 	fClusterType(ref.fClusterType),
	fMinEtaCut(ref.fMinEtaCut),
	fMaxEtaCut(ref.fMaxEtaCut),
	fUseEtaCut(ref.fUseEtaCut),
	fMinPhiCut(ref.fMinPhiCut),
	fMaxPhiCut(ref.fMaxPhiCut),
	fUsePhiCut(ref.fUsePhiCut),
	fMinDistanceToBadChannel(ref.fMinDistanceToBadChannel),
	fUseDistanceToBadChannel(ref.fUseDistanceToBadChannel),
	fMaxTimeDiff(ref.fMaxTimeDiff),
	fUseTimeDiff(ref.fUseTimeDiff),
	fMinDistTrackToCluster(ref.fMinDistTrackToCluster),
	fUseDistTrackToCluster(ref.fUseDistTrackToCluster),
	fExoticCell(ref.fExoticCell),
	fUseExoticCell(ref.fUseExoticCell),
	fMinEnergy(ref.fMinEnergy),
	fUseMinEnergy(ref.fUseMinEnergy),
	fMinNCells(ref.fMinNCells),
	fUseNCells(ref.fUseNCells),
	fMaxM02(ref.fMaxM02),
	fMinM02(ref.fMinM02),
	fUseM02(ref.fUseM02),
	fMaxM20(ref.fMaxM20),
	fMinM20(ref.fMinM20),
	fUseM20(ref.fUseDispersion),
	fMaxDispersion(ref.fMaxDispersion),
	fUseDispersion(ref.fUseDispersion),
	fMinNLM(ref.fMinNLM),
	fMaxNLM(ref.fMaxNLM),
	fUseNLM(ref.fUseNLM),
	fCutString(NULL),
	fHistCutIndex(NULL),
	fHistAcceptanceCuts(NULL),
	fHistClusterIdentificationCuts(NULL),
	fHistClusterEtavsPhiBeforeAcc(NULL),
	fHistClusterEtavsPhiAfterAcc(NULL),
	fHistClusterEtavsPhiAfterQA(NULL),
	fHistDistanceToBadChannelBeforeAcc(NULL),
	fHistDistanceToBadChannelAfterAcc(NULL),
	fHistClusterRBeforeQA(NULL),
	fHistClusterRAfterQA(NULL),
	fHistClusterdEtadPhiBeforeQA(NULL),
	fHistClusterdEtadPhiAfterQA(NULL),
	fHistClusterTimevsEBeforeQA(NULL),
	fHistClusterTimevsEAfterQA(NULL),
	fHistExoticCellBeforeQA(NULL),
	fHistExoticCellAfterQA(NULL),
	fHistNMatchedTracks(NULL),
	fHistDistanceTrackToClusterBeforeQA(NULL),
	fHistDistanceTrackToClusterAfterQA(NULL),
	fHistEnergyOfClusterBeforeQA(NULL),
	fHistEnergyOfClusterAfterQA(NULL),
	fHistNCellsBeforeQA(NULL),
	fHistNCellsAfterQA(NULL),
	fHistM02BeforeQA(NULL),
	fHistM02AfterQA(NULL),
	fHistM20BeforeQA(NULL),
	fHistM20AfterQA(NULL),
	fHistDispersionBeforeQA(NULL),
	fHistDispersionAfterQA(NULL),
	fHistNLMBeforeQA(NULL),
	fHistNLMAfterQA(NULL)
{
   // Copy Constructor
   for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=ref.fCuts[jj];}
   fCutString=new TObjString((GetCutNumber()).Data());

}


//________________________________________________________________________
AliCaloPhotonCuts::~AliCaloPhotonCuts() {
   // Destructor
   //Deleting fHistograms leads to seg fault it it's added to output collection of a task
   // if(fHistograms)
   //    delete fHistograms;
   // fHistograms = NULL;
   if(fCutString != NULL){
      delete fCutString;
      fCutString = NULL;
   }
}

//________________________________________________________________________
void AliCaloPhotonCuts::InitCutHistograms(TString name){

	// Initialize Cut Histograms for QA (only initialized and filled if function is called)
	TH1::AddDirectory(kFALSE);

	if(fHistograms != NULL){
		delete fHistograms;
		fHistograms=NULL;
	}
	if(fHistograms==NULL){
		fHistograms=new TList();
		fHistograms->SetOwner(kTRUE);
		if(name=="")fHistograms->SetName(Form("CaloCuts_%s",GetCutNumber().Data()));
		else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
	}

	// IsPhotonSelected
	fHistCutIndex=new TH1F(Form("IsPhotonSelected %s",GetCutNumber().Data()),"IsPhotonSelected",5,-0.5,4.5);
	fHistCutIndex->GetXaxis()->SetBinLabel(kPhotonIn+1,"in");
	fHistCutIndex->GetXaxis()->SetBinLabel(kDetector+1,"detector");
	fHistCutIndex->GetXaxis()->SetBinLabel(kAcceptance+1,"acceptance");
	fHistCutIndex->GetXaxis()->SetBinLabel(kClusterQuality+1,"cluster QA");
	fHistCutIndex->GetXaxis()->SetBinLabel(kPhotonOut+1,"out");
	fHistograms->Add(fHistCutIndex);

	// Acceptance Cuts
	fHistAcceptanceCuts=new TH1F(Form("AcceptanceCuts %s",GetCutNumber().Data()),"AcceptanceCuts",5,-0.5,4.5);
	fHistAcceptanceCuts->GetXaxis()->SetBinLabel(1,"in");
	fHistAcceptanceCuts->GetXaxis()->SetBinLabel(2,"eta");
	fHistAcceptanceCuts->GetXaxis()->SetBinLabel(3,"phi");
	fHistAcceptanceCuts->GetXaxis()->SetBinLabel(4,"distance to bad channel");
	fHistAcceptanceCuts->GetXaxis()->SetBinLabel(5,"out");
	fHistograms->Add(fHistAcceptanceCuts);

	// Cluster Cuts
	fHistClusterIdentificationCuts =new TH1F(Form("ClusterQualityCuts %s",GetCutNumber().Data()),"ClusterQualityCuts",11,-0.5,10.5);
	fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(1,"in");
	fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(2,"timing");
	fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(3,"track matching");
	fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(4,"Exotics");
	fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(5,"minimum energy");
	fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(6,"minimum NCells");
	fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(7,"M02");
	fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(8,"M20");
	fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(9,"dispersion");
	fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(10,"NLM");
	fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(11,"out");
	fHistograms->Add(fHistClusterIdentificationCuts);

	// Acceptance related histogramms
	fHistClusterEtavsPhiBeforeAcc=new TH2F(Form("EtaPhi_beforeAcceptance %s",GetCutNumber().Data()),"EtaPhi_beforeAcceptance",462,-TMath::Pi(),TMath::Pi(),110,-0.7,0.7);
	fHistograms->Add(fHistClusterEtavsPhiBeforeAcc);
	fHistClusterEtavsPhiAfterAcc=new TH2F(Form("EtaPhi_afterAcceptance %s",GetCutNumber().Data()),"EtaPhi_afterAcceptance",462,-TMath::Pi(),TMath::Pi(),110,-0.7,0.7);
	fHistograms->Add(fHistClusterEtavsPhiAfterAcc);
	fHistClusterEtavsPhiAfterQA=new TH2F(Form("EtaPhi_afterClusterQA %s",GetCutNumber().Data()),"EtaPhi_afterClusterQA",462,-TMath::Pi(),TMath::Pi(),110,-0.7,0.7);
	fHistograms->Add(fHistClusterEtavsPhiAfterQA);
	fHistDistanceToBadChannelBeforeAcc = new TH1F(Form("DistanceToBadChannel_beforeAcceptance %s",GetCutNumber().Data()),"DistanceToBadChannel_beforeAcceptance",200,0,40);
	fHistograms->Add(fHistDistanceToBadChannelBeforeAcc);
	fHistDistanceToBadChannelAfterAcc = new TH1F(Form("DistanceToBadChannel_afterAcceptance %s",GetCutNumber().Data()),"DistanceToBadChannel_afterAcceptance",200,0,40);
	fHistograms->Add(fHistDistanceToBadChannelAfterAcc);
	
	// Cluster quality related histograms
	fHistClusterRBeforeQA = new TH1F(Form("R_Cluster_beforeClusterQA %s",GetCutNumber().Data()),"R of cluster",200,400,500);
	fHistograms->Add(fHistClusterRBeforeQA);
	fHistClusterRAfterQA = new TH1F(Form("R_Cluster_afterClusterQA %s",GetCutNumber().Data()),"R of cluster_matched",200,400,500);
	fHistograms->Add(fHistClusterRAfterQA);
	fHistClusterdEtadPhiBeforeQA=new TH2F(Form("dEtaVsdPhi_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_beforeClusterQA",240,-0.3,0.3,240,-0.3,0.3);
	fHistograms->Add(fHistClusterdEtadPhiBeforeQA);
	fHistClusterdEtadPhiAfterQA=new TH2F(Form("dEtaVsdPhi_afterClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_afterClusterQA",240,-0.3,0.3,240,-0.3,0.3);
	fHistograms->Add(fHistClusterdEtadPhiAfterQA);
	fHistClusterTimevsEBeforeQA=new TH2F(Form("ClusterTimeVsE_beforeClusterQA %s",GetCutNumber().Data()),"ClusterTimeVsE_beforeClusterQA",400,-10e-6,10e-6,100,0.,40);
	fHistograms->Add(fHistClusterTimevsEBeforeQA);
	fHistClusterTimevsEAfterQA=new TH2F(Form("ClusterTimeVsE_afterClusterQA %s",GetCutNumber().Data()),"ClusterTimeVsE_afterClusterQA",400,-10e-6,10e-6,100,0.,40);
	fHistograms->Add(fHistClusterTimevsEAfterQA);
	fHistExoticCellBeforeQA=new TH2F(Form("ExoticCell_beforeClusterQA %s",GetCutNumber().Data()),"ExoticCell_beforeClusterQA",400,0,40,50,0.75,1);
	fHistograms->Add(fHistExoticCellBeforeQA);
	fHistExoticCellAfterQA=new TH2F(Form("ExoticCell_afterClusterQA %s",GetCutNumber().Data()),"ExoticCell_afterClusterQA",400,0,40,50,0.75,1);
	fHistograms->Add(fHistExoticCellAfterQA);
	fHistNMatchedTracks = new TH1F(Form("NMatchedTracks_%s",GetCutNumber().Data()),"NMatchedTracks",22,-1.5,20.5);
	fHistograms->Add(fHistNMatchedTracks);
	fHistDistanceTrackToClusterBeforeQA = new TH1F(Form("DistanceToTrack_beforeClusterQA %s",GetCutNumber().Data()),"DistanceToTrack_beforeClusterQA",200,0,2);
	fHistograms->Add(fHistDistanceTrackToClusterBeforeQA);
	fHistDistanceTrackToClusterAfterQA = new TH1F(Form("DistanceToTrack_afterClusterQA %s",GetCutNumber().Data()),"DistanceToTrack_afterClusterQA",200,0,2);
	fHistograms->Add(fHistDistanceTrackToClusterAfterQA);
	fHistEnergyOfClusterBeforeQA = new TH1F(Form("EnergyOfCluster_beforeClusterQA %s",GetCutNumber().Data()),"EnergyOfCluster_beforeClusterQA",300,0,30);
	fHistograms->Add(fHistEnergyOfClusterBeforeQA);
	fHistEnergyOfClusterAfterQA = new TH1F(Form("EnergyOfCluster_afterClusterQA %s",GetCutNumber().Data()),"EnergyOfCluster_afterClusterQA",300,0,30);
	fHistograms->Add(fHistEnergyOfClusterAfterQA);
	fHistNCellsBeforeQA = new TH1F(Form("NCellPerCluster_beforeClusterQA %s",GetCutNumber().Data()),"NCellPerCluster_beforeClusterQA",50,0,50);
	fHistograms->Add(fHistNCellsBeforeQA);
	fHistNCellsAfterQA = new TH1F(Form("NCellPerCluster_afterClusterQA %s",GetCutNumber().Data()),"NCellPerCluster_afterClusterQA",50,0,50);
	fHistograms->Add(fHistNCellsAfterQA);
	fHistM02BeforeQA = new TH1F(Form("M02_beforeClusterQA %s",GetCutNumber().Data()),"M02_beforeClusterQA",100,0,5);
	fHistograms->Add(fHistM02BeforeQA);
	fHistM02AfterQA = new TH1F(Form("M02_afterClusterQA %s",GetCutNumber().Data()),"M02_afterClusterQA",100,0,5);
	fHistograms->Add(fHistM02AfterQA);
	fHistM20BeforeQA = new TH1F(Form("M20_beforeClusterQA %s",GetCutNumber().Data()),"M20_beforeClusterQA",100,0,2.5);
	fHistograms->Add(fHistM20BeforeQA);
	fHistM20AfterQA = new TH1F(Form("M20_afterClusterQA %s",GetCutNumber().Data()),"M20_afterClusterQA",100,0,2.5);
	fHistograms->Add(fHistM20AfterQA);
	fHistDispersionBeforeQA = new TH1F(Form("Dispersion_beforeClusterQA %s",GetCutNumber().Data()),"Dispersion_beforeClusterQA",100,0,4);
	fHistograms->Add(fHistDispersionBeforeQA);
	fHistDispersionAfterQA = new TH1F(Form("Dispersion_afterClusterQA %s",GetCutNumber().Data()),"Dispersion_afterClusterQA",100,0,4);
	fHistograms->Add(fHistDispersionAfterQA);
	fHistNLMBeforeQA = new TH1F(Form("NLM_beforeClusterQA %s",GetCutNumber().Data()),"NLM_beforeClusterQA",10,0,10);
	fHistograms->Add(fHistNLMBeforeQA);
	fHistNLMAfterQA = new TH1F(Form("NLM_afterClusterQA %s",GetCutNumber().Data()),"NLM_afterClusterQA",10,0,10);
	fHistograms->Add(fHistNLMAfterQA);
	
	TH1::AddDirectory(kTRUE);
}


///________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelectedMC(TParticle *particle,AliStack *fMCStack){
   // MonteCarlo Photon Selection

	if(!fMCStack)return kFALSE;

	if (particle->GetPdgCode() == 22){

		if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
		if ( particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut ) return kFALSE;
		
		if(particle->GetMother(0) >-1 && fMCStack->Particle(particle->GetMother(0))->GetPdgCode() == 22){
			return kFALSE; // no photon as mothers!
		}
		if(particle->GetMother(0) >= fMCStack->GetNprimary()){
			return kFALSE; // the gamma has a mother, and it is not a primary particle
		}
		return kTRUE;
	}
	return kFALSE;
}
///________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelectedAODMC(AliAODMCParticle *particle,TClonesArray *aodmcArray){
	// MonteCarlo Photon Selection

	if(!aodmcArray)return kFALSE;
	if (particle->GetPdgCode() == 22){
		if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
		if ( particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut ) return kFALSE;
		if(particle->GetMother() > -1){
			if((static_cast<AliAODMCParticle*>(aodmcArray->At(particle->GetMother())))->GetPdgCode() == 22){
				return kFALSE; // no photon as mothers!
			}
			if(!(static_cast<AliAODMCParticle*>(aodmcArray->At(particle->GetMother()))->IsPrimary())){
				return kFALSE; // the gamma has a mother, and it is not a primary particle
			}
		}
		return kTRUE; // return in case of accepted gamma
	}
	return kFALSE;
}



///________________________________________________________________________
// This function selects the clusters based on their quality criteria
///________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterQualityCuts(AliVCluster* cluster, AliVEvent *event, Bool_t isMC)
{   // Specific Photon Cuts

	Int_t cutIndex = 0;
	if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);
	cutIndex++;

// 	Double_t minR = 999.0;
	// get the minimum radius of tracks to cluster
// 	if(fHistDistanceTrackToClusterBeforeQA || fHistDistanceTrackToClusterAfterQA){
// 		Float_t pos[3];
// 		cluster->GetPosition(pos);  // Get cluster position
// 		
// 		TVector3 cp(pos);
// 		int NtrMatched = 0;
// 		NtrMatched = cluster->GetNTracksMatched();
// 		fHistNMatchedTracks->Fill(NtrMatched);
// 		//loop over tracks for Jet QA
// 		TList *l = event->GetList();
// 		TClonesArray *tracks = dynamic_cast<TClonesArray*>(l->FindObject("Tracks"));
// 		for(int itrack = 0; itrack < NtrMatched; itrack++){
// 			AliVTrack *trackcluster = static_cast<AliVTrack*>(tracks->At(itrack));
// 			if (! trackcluster) {
// 				AliError(Form("Couldn't get ESD track %d\n", itrack));
// 				continue;
// 			}
// 			Double_t dphi = -999.0;
// 			Double_t deta = -999.0;
// 			AliPicoTrack::GetEtaPhiDiff(trackcluster, cluster, dphi, deta);
// 			cout << "here" << endl;
// 			Double_t dr = sqrt(dphi*dphi + deta+deta);
// 			if(dr < minR)
// 				minR = dr;
// 		}//loop over tracks
// 	}
	
	// Fill Histos before Cuts
	if(fHistClusterTimevsEBeforeQA) fHistClusterTimevsEBeforeQA->Fill(cluster->GetTOF(), cluster->E());
// 	if(fHistExoticCellBeforeQA) fHistExoticCellBeforeQA->Fill(cluster->E(), );
// 	if(fHistDistanceTrackToClusterBeforeQA) fHistDistanceTrackToClusterBeforeQA->Fill(minR);
	if(fHistEnergyOfClusterBeforeQA) fHistEnergyOfClusterBeforeQA->Fill(cluster->E());
	if(fHistNCellsBeforeQA) fHistNCellsBeforeQA->Fill(cluster->GetNCells());
	if(fHistM02BeforeQA) fHistM02BeforeQA->Fill(cluster->GetM02());
	if(fHistM20BeforeQA) fHistM20BeforeQA->Fill(cluster->GetM20());
	if(fHistDispersionBeforeQA) fHistDispersionBeforeQA->Fill(cluster->GetDispersion());
// 	if(fHistNLMBeforeQA) fHistNLMBeforeQA->Fill(cluster->GetNExMax());
	
	// Check wether timing is ok
	if (fUseTimeDiff){
		if(abs(cluster->GetTOF()) > fMaxTimeDiff && !isMC){
			if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex); //1
			return kFALSE;
		}
	}	
	cutIndex++; //2, next cut

	// Minimum distance to track
// 	if (fUseDistTrackToCluster){
//  		Float_t pos[3];
//  		cluster->GetPosition(pos);  // Get cluster position
//  		TVector3 cp(pos);
//  		int NtrMatched = 0;
//  		NtrMatched = cluster->GetNTracksMatched();
// 		fHistNMatchedTracks->Fill(NtrMatched);
//  		if(NtrMatched > 0){
// 			//loop over tracks for QA
// 			TList *l = event->GetList();
// 			TClonesArray *tracks = dynamic_cast<TClonesArray*>(l->FindObject("Tracks"));
// 			
// 			Double_t dphi = 999.0;
// 			Double_t deta = 999.0;
// 			Double_t dr2 = 999.0;
// 
// 			for(int itrack = 0; itrack < NtrMatched; itrack++){
// 				AliVTrack *trackcluster = NULL;
// 				trackcluster = static_cast<AliVTrack*>(tracks->At(itrack));
// 				if (! trackcluster) {
// 				AliError(Form("Couldn't get ESD track %d\n", itrack));
// 				continue;
// 				}
// 				AliPicoTrack::GetEtaPhiDiff(trackcluster, cluster, dphi, deta);
// 				dr2 = dphi*dphi + deta+deta;
// 				//cout << dr << endl;
// 				if(dr2 < fMinDistTrackToCluster*fMinDistTrackToCluster){
// 		//        if(dphi < fMinDistTrackToCluster || deta < fMinDistTrackToCluster){
// 				if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex); //2
// 				return kFALSE;
// 				}
// 				
// 			}//loop over tracks
// 		}
// //		if(cluster->GetEmcCpvDistance() < fMinDistTrackToCluster){
// //			if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex); //2
// //			return kFALSE;
// //		}
// 	}
	cutIndex++;//3, next cut

	// exotic cell cut --IMPLEMENT LATER---
// 	if(!AcceptanceCuts(photon)){
// 		if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex); //3
// 		return kFALSE;
// 	}
	cutIndex++; //4, next cut
	
	// minimum cell energy cut
	if (fUseMinEnergy){
		if(cluster->E() < fMinEnergy){
			if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex); //4
			return kFALSE;
		}
	}	
	cutIndex++; //5, next cut
	
	// minimum number of cells
	if (fUseNCells){
		if(cluster->GetNCells() < fMinNCells) {
			if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex); //5
			return kFALSE;
		}
	}	
	cutIndex++; //6, next cut
	
	// M02 cut
	if (fUseM02){
		if( cluster->GetM02()< fMinM02 || cluster->GetM02() > fMaxM02 ) {
			if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex); //6
			return kFALSE;
		}
	}	
	cutIndex++; //7, next cut
	
	// M20 cut
	if (fUseM20){
		if( cluster->GetM20()< fMinM20 || cluster->GetM20() > fMaxM20 ) {
			if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex); //7
			return kFALSE;
		}
	}	
	cutIndex++; //8, next cut
	
	// dispersion cut
	if (fUseDispersion){
		if( cluster->GetDispersion()> fMaxDispersion) {
			if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex); //8
			return kFALSE;
		}
	}	
	cutIndex++; //9, next cut
	
	// NLM cut --IMPLEMENT LATER---
// 	if (fUseNLM){
// 		if( cluster->GetDispersion()> fMaxDispersion) {
// 			if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex); //9
// 			return kFALSE;
// 		}
// 	}	
	cutIndex++; //9, next cut
	
	// DONE with selecting photons
	if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex); //10

	// Histos after Cuts
	Double_t vertex[3] = {0};
	event->GetPrimaryVertex()->GetXYZ(vertex);
	// TLorentzvector with cluster
	TLorentzVector clusterVector;
	cluster->GetMomentum(clusterVector,vertex);
	Double_t etaCluster = clusterVector.Eta();
	Double_t phiCluster = clusterVector.Phi();
	
	if(fHistClusterEtavsPhiAfterQA) fHistClusterEtavsPhiAfterQA->Fill(phiCluster,etaCluster);
	if(fHistClusterTimevsEAfterQA) fHistClusterTimevsEAfterQA->Fill(cluster->GetTOF(), cluster->E());
// 	if(fHistExoticCellAfterQA) fHistExoticCellAfterQA->Fill(cluster->E(), );
// 	if(fHistDistanceTrackToClusterAfterQA) fHistDistanceTrackToClusterAfterQA->Fill(minR);
	if(fHistEnergyOfClusterAfterQA) fHistEnergyOfClusterAfterQA->Fill(cluster->E());
	if(fHistNCellsAfterQA) fHistNCellsAfterQA->Fill(cluster->GetNCells());
	if(fHistM02AfterQA) fHistM02AfterQA->Fill(cluster->GetM02());
	if(fHistM20AfterQA) fHistM20AfterQA->Fill(cluster->GetM20());
	if(fHistDispersionAfterQA) fHistDispersionAfterQA->Fill(cluster->GetDispersion());
// 	if(fHistNLMBeforeQA) fHistNLMAfterQA->Fill(cluster->GetNExMax());

	return kTRUE;

}


///________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelected(AliVCluster *cluster, AliVEvent * event, Bool_t isMC)
{
	//Selection of Reconstructed photon clusters with Calorimeters

	FillClusterCutIndex(kPhotonIn);

	Double_t vertex[3] = {0};
	event->GetPrimaryVertex()->GetXYZ(vertex);
	// TLorentzvector with cluster
	TLorentzVector clusterVector;
	cluster->GetMomentum(clusterVector,vertex);
	Double_t etaCluster = clusterVector.Eta();
	Double_t phiCluster = clusterVector.Phi();

	// Histos before cuts
	if(fHistClusterEtavsPhiBeforeAcc) fHistClusterEtavsPhiBeforeAcc->Fill(phiCluster,etaCluster);
	
	// Cluster Selection - 0= accept any calo cluster
	if (fClusterType > 0){
		//Select EMCAL cluster
		if (fClusterType == 1 && !cluster->IsEMCAL()){
			FillClusterCutIndex(kDetector);
			return kFALSE;
		}
		//Select PHOS cluster
		if (fClusterType == 2 && !cluster->IsPHOS()){
			FillClusterCutIndex(kDetector);
			return kFALSE;
		}
	}
	
	// Acceptance Cuts
	if(!AcceptanceCuts(cluster,event)){
		FillClusterCutIndex(kAcceptance);
		return kFALSE;
	}
	// Cluster Quality Cuts
	if(!ClusterQualityCuts(cluster,event,isMC)){
		FillClusterCutIndex(kClusterQuality);
		return kFALSE;
	}

	// Photon passed cuts
	FillClusterCutIndex(kPhotonOut);
	return kTRUE;
}


///________________________________________________________________________
Bool_t AliCaloPhotonCuts::AcceptanceCuts(AliVCluster *cluster, AliVEvent* event) 
{
   // Exclude certain areas for photon reconstruction

	Int_t cutIndex=0;
	if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
	cutIndex++;

	
	Double_t vertex[3] = {0};
	event->GetPrimaryVertex()->GetXYZ(vertex);
	// TLorentzvector with cluster
	TLorentzVector clusterVector;
	cluster->GetMomentum(clusterVector,vertex);
	Double_t etaCluster = clusterVector.Eta();
	Double_t phiCluster = clusterVector.Phi();
	
	// check eta range
	if (fUseEtaCut){
		if (etaCluster < fMinEtaCut || etaCluster > fMaxEtaCut){
			if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
			return kFALSE;
		}
	}
	cutIndex++;
	
	// check phi range
	if (fUsePhiCut ){
		if (phiCluster < fMinPhiCut || phiCluster > fMaxEtaCut){
			if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
			return kFALSE;
		}
	}
	cutIndex++;
	
	// check distance to bad channel
	if (fUseDistanceToBadChannel){
		if (cluster->GetDistanceToBadChannel() < fMinDistanceToBadChannel){
			if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
			return kFALSE;
		}	
	}
	cutIndex++;
	if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);

	// Histos after cuts
	if(fHistClusterEtavsPhiAfterAcc) fHistClusterEtavsPhiAfterAcc->Fill(phiCluster,etaCluster);
	
	return kTRUE;
}

Bool_t AliCaloPhotonCuts::MatchConvPhotonToCluster(AliAODConversionPhoton* convPhoton, AliVCluster* cluster, AliVEvent* event ){

	if (!fUseDistTrackToCluster) return kFALSE;
	
	AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(event);
	AliAODEvent *aodev = 0;
	if (!esdev) {
		aodev = dynamic_cast<AliAODEvent*>(event);
		if (!aodev) {
			AliError("Task needs AOD or ESD event, returning");
			return kFALSE;
		}
	}

    Double_t vertex[3] = {0,0,0};
	event->GetPrimaryVertex()->GetXYZ(vertex);

    if(!cluster->IsEMCAL() && !cluster->IsPHOS()){AliError("Cluster is neither EMCAL nor PHOS, returning"); return kFALSE;}

    Float_t clusterPosition[3] = {0,0,0};
    cluster->GetPosition(clusterPosition);
    Double_t clusterR = TMath::Sqrt( clusterPosition[0]*clusterPosition[0] + clusterPosition[1]*clusterPosition[1] );
    if(fHistClusterRBeforeQA) fHistClusterRBeforeQA->Fill(clusterR);

//cout << "+++++++++ Cluster: x, y, z, R" << clusterPosition[0] << ", " << clusterPosition[1] << ", " << clusterPosition[2] << ", " << clusterR << "+++++++++" << endl;

	Bool_t matched = kFALSE;
	for (Int_t i = 0; i < 2; i++){
		Int_t tracklabel = convPhoton->GetLabel(i);
		AliVTrack *inTrack = 0x0;	
		if(esdev) {
			if(tracklabel > event->GetNumberOfTracks() ) continue;
			inTrack = esdev->GetTrack(tracklabel);
		} else {
			if(((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1"))->AreAODsRelabeled()){
				inTrack = dynamic_cast<AliVTrack*>(event->GetTrack(tracklabel));	
			} else {
				for(Int_t ii=0; ii<event->GetNumberOfTracks(); ii++) {
					inTrack = dynamic_cast<AliVTrack*>(event->GetTrack(ii));
					if(inTrack){
						if(inTrack->GetID() == tracklabel) {
							continue;
						}
					}
				}
			}
		}
	if( inTrack->Pt() < 0.005 ) continue;

		AliESDtrack *esdt = dynamic_cast<AliESDtrack*>(inTrack);
		AliAODTrack *aodt = 0;
		if (!esdt) {
			aodt = dynamic_cast<AliAODTrack*>(inTrack);
			if (!aodt){AliError("Track is neither ESD nor AOD, continue"); continue;}
		}

		AliExternalTrackParam *trackParam = 0;
		if (esdt) {
			const AliExternalTrackParam *in = esdt->GetInnerParam();
			if (!in){AliError("Could not get InnerParam of Track, continue"); continue;}
			trackParam = new AliExternalTrackParam(*in);
		} else {
			Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
			aodt->PxPyPz(pxpypz);
			aodt->XvYvZv(xyz);
			aodt->GetCovarianceXYZPxPyPz(cv);
			trackParam = new AliExternalTrackParam(xyz,pxpypz,cv,aodt->Charge());
		}

		Bool_t propagated = kFALSE;
		AliExternalTrackParam emcParam(*trackParam);
		Float_t dPhi = 0;
		Float_t dEta = 0;

		if(cluster->IsEMCAL()){
			Float_t eta = 0; Float_t phi = 0; Float_t pt = 0;
			propagated = AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 430, 0.000510999, 20, eta, phi, pt);
			if(propagated){
				propagated = AliEMCALRecoUtils::ExtrapolateTrackToCluster(&emcParam, cluster, 0.000510999, 5, dEta, dPhi);
			}
		}
		if(cluster->IsPHOS()){
			propagated = AliTrackerBase::PropagateTrackToBxByBz(&emcParam, clusterR, 0.000510999, 20, kTRUE, 0.8, -1);
			if (propagated){
				Double_t trkPos[3] = {0,0,0};
				emcParam.GetXYZ(trkPos);
				TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);
				TVector3 clsPosVec(clusterPosition);
				dPhi = clsPosVec.DeltaPhi(trkPosVec);
				dEta = clsPosVec.Eta()-trkPosVec.Eta();
			}
		}

		if (propagated){
			Float_t dR2 = dPhi*dPhi + dEta*dEta;
			if (fHistDistanceTrackToClusterBeforeQA)fHistDistanceTrackToClusterBeforeQA->Fill(TMath::Sqrt(dR2));
			if (fHistClusterdEtadPhiBeforeQA) fHistClusterdEtadPhiBeforeQA->Fill(dEta, dPhi);
			if(dR2 < fMinDistTrackToCluster*fMinDistTrackToCluster){
				matched = kTRUE;
				if (fHistDistanceTrackToClusterAfterQA)fHistDistanceTrackToClusterAfterQA->Fill(TMath::Sqrt(dR2));
				if (fHistClusterdEtadPhiAfterQA) fHistClusterdEtadPhiAfterQA->Fill(dEta, dPhi);
				if (fHistClusterRAfterQA) fHistClusterRAfterQA->Fill(clusterR);
			}

		}
		delete trackParam;
	}

	return matched;

}

//____________________________________________________________________________________________


///________________________________________________________________________
Bool_t AliCaloPhotonCuts::UpdateCutString() {
   ///Update the cut string (if it has been created yet)

   if(fCutString && fCutString->GetString().Length() == kNCuts) {
      fCutString->SetString(GetCutNumber());
   } else {
      return kFALSE;
   }
   return kTRUE;
}

///________________________________________________________________________
Bool_t AliCaloPhotonCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
	// Initialize Cuts from a given Cut string
	AliInfo(Form("Set CaloCut Number: %s",analysisCutSelection.Data()));
	if(analysisCutSelection.Length()!=kNCuts) {
		AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
		return kFALSE;
	}
	if(!analysisCutSelection.IsDigit()){
		AliError("Cut selection contains characters");
		return kFALSE;
	}

	const char *cutSelection = analysisCutSelection.Data();
	#define ASSIGNARRAY(i)  fCuts[i] = cutSelection[i] - '0'
	for(Int_t ii=0;ii<kNCuts;ii++){
		ASSIGNARRAY(ii);
	}

	// Set Individual Cuts
	for(Int_t ii=0;ii<kNCuts;ii++){
		if(!SetCut(cutIds(ii),fCuts[ii]))return kFALSE;
	}
	PrintCutsWithValues();
	return kTRUE;
}

///________________________________________________________________________
Bool_t AliCaloPhotonCuts::SetCut(cutIds cutID, const Int_t value) {
	///Set individual cut ID

	switch (cutID) {		
		
		case kClusterType:
			if( SetClusterTypeCut(value)) {
				fCuts[kClusterType] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;
		
		case kEtaMin:
			if( SetMinEtaCut(value)) {
				fCuts[kEtaMin] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;

		case kEtaMax:
			if( SetMaxEtaCut(value)) {
				fCuts[kEtaMax] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;

		case kPhiMin:
			if( SetMinPhiCut(value)) {
				fCuts[kPhiMin] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;

		case kPhiMax:
			if( SetMaxPhiCut(value)) {
				fCuts[kPhiMax] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;

		case kDistanceToBadChannel:
			if( SetDistanceToBadChannelCut(value)) {
				fCuts[kDistanceToBadChannel] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;

		case kTiming:
			if( SetTimingCut(value)) {
				fCuts[kTiming] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;

		case kTrackMatching:
			if( SetTrackMatchingCut(value)) {
				fCuts[kTrackMatching] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;

		case kExoticCell:
			if( SetExoticCellCut(value)) {
				fCuts[kExoticCell] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;

		case kMinEnery:
			if( SetMinEnergyCut(value)) {
				fCuts[kMinEnery] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;

		case kNMinCells:
			if( SetMinNCellsCut(value)) {
				fCuts[kNMinCells] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;
			
		case kMinM02:
			if( SetMinM02(value)) {
				fCuts[kMinM02] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;

		case kMaxM02:
			if( SetMaxM02(value)) {
				fCuts[kMaxM02] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;
		
		case kMinM20:
			if( SetMinM20(value)) {
				fCuts[kMinM20] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;

		case kMaxM20:
			if( SetMaxM20(value)) {
				fCuts[kMaxM20] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;

		case kDispersion:
			if( SetDispersion(value)) {
				fCuts[kDispersion] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;

		case kNLM:
			if( SetNLM(value)) {
				fCuts[kNLM] = value;
				UpdateCutString();
				return kTRUE;
			} else return kFALSE;

		case kNCuts:
			AliError("Cut id out of range");
			return kFALSE;
	}

	AliError("Cut id %d not recognized");
	return kFALSE;


}
///________________________________________________________________________
void AliCaloPhotonCuts::PrintCuts() {
   // Print out current Cut Selection
   for(Int_t ic = 0; ic < kNCuts; ic++) {
      printf("%-30s : %d \n", fgkCutNames[ic], fCuts[ic]);
   }
}

void AliCaloPhotonCuts::PrintCutsWithValues() {
	// Print out current Cut Selection with value
	printf("\nCluster cutnumber \n");
	for(Int_t ic = 0; ic < kNCuts; ic++) {
		printf("%d",fCuts[ic]);
	}
	printf("\n\n");

	printf("Acceptance cuts: \n");
	if (fClusterType == 0) printf("\tall calorimeter clusters are used\n");
	if (fClusterType == 1) printf("\tEMCAL calorimeter clusters are used\n");
	if (fClusterType == 2) printf("\tPHOS calorimeter clusters are used\n");
	if (fUseEtaCut) printf("\t%3.2f < eta_{cluster} < %3.2f\n", fMinEtaCut, fMaxEtaCut );
	if (fUsePhiCut) printf("\t%3.2f < phi_{cluster} < %3.2f\n", fMinPhiCut, fMaxPhiCut );
	if (fUseDistanceToBadChannel) printf("\tcut on exotics applied \n");
	
	printf("Cluster Quality cuts: \n");
	if (fUseTimeDiff) printf("\t time difference < %3.2f\n", fMaxTimeDiff );
	if (fUseDistTrackToCluster) printf("\tmin distance to track > %3.2f\n", fMinDistTrackToCluster );
	if (fUseExoticCell)printf("\t min distance to track > %3.2f\n", fMinDistTrackToCluster );
    if (fUseMinEnergy)printf("\t E_{cluster} > %3.2f\n", fMinEnergy );
	if (fUseNCells) printf("\t number of cells per cluster > %d\n", fMinNCells );
	if (fUseM02) printf("\t %3.2f < M02 < %3.2f\n", fMinM02, fMaxM02 );
	if (fUseM20) printf("\t %3.2f < M20 < %3.2f\n", fMinM20, fMaxM20 );
	if (fUseDispersion) printf("\t dispersion < %3.2f\n", fMaxDispersion );
	if (fUseNLM) printf("\t %d < NLM < %d\n", fMinNLM, fMaxNLM );
	
}

// EMCAL acceptance 2011
// 1.39626, 3.125 (phi)
// -0.66687,,0.66465


///________________________________________________________________________
Bool_t AliCaloPhotonCuts::SetClusterTypeCut(Int_t clusterType)
{   // Set Cut
	switch(clusterType){
	case 0: // all clusters
		fClusterType=0;
		break;
	case 1: // EMCAL clusters
		fClusterType=1;
		break;
	case 2: // PHOS clusters
		fClusterType=2;
		break;
	default:
		AliError(Form("ClusterTypeCut not defined %d",clusterType));
		return kFALSE;
	}
	return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinEtaCut(Int_t minEta)
{
	switch(minEta){
	case 0:
		if (!fUseEtaCut) fUseEtaCut=0;
		fMinEtaCut=-10.;
		break;
	case 1:
		if (!fUseEtaCut) fUseEtaCut=1;
		fMinEtaCut=-0.6687;
		break;
	case 2: 
		if (!fUseEtaCut) fUseEtaCut=1;
		fMinEtaCut=-0.5;
		break;
	case 3: 
		if (!fUseEtaCut) fUseEtaCut=1;
		fMinEtaCut=-2;
		break;
	default:
		AliError(Form("MinEta Cut not defined %d",minEta));
		return kFALSE;
	}
	return kTRUE;
}


//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMaxEtaCut(Int_t maxEta)
{
	switch(maxEta){
	case 0: 
		if (!fUseEtaCut) fUseEtaCut=0;
		fMaxEtaCut=10;
		break;		
	case 1:
		if (!fUseEtaCut) fUseEtaCut=1;
		fMaxEtaCut=0.66465;
		break;
	case 2: 
		if (!fUseEtaCut) fUseEtaCut=1;
		fMaxEtaCut=0.5;
		break;
	case 3: 
		if (!fUseEtaCut) fUseEtaCut=1;
		fMaxEtaCut=2;
		break;
	default:
		AliError(Form("MaxEta Cut not defined %d",maxEta));
		return kFALSE;
	}
	return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinPhiCut(Int_t minPhi)
{
	switch(minPhi){
	case 0: 
		if (!fUsePhiCut) fUsePhiCut=0;
		fMinPhiCut=-10000;
		break;
	case 1: 
		if (!fUsePhiCut) fUsePhiCut=1;
		fMinPhiCut=1.39626;
		break;
	default:
		AliError(Form("MinPhi Cut not defined %d",minPhi));
		return kFALSE;
	}
	return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMaxPhiCut(Int_t maxPhi)
{
	switch(maxPhi){
	case 0: 
		if (!fUsePhiCut) fUsePhiCut=0;
		fMaxPhiCut=10000;
		break;
	case 1: 
		if (!fUsePhiCut) fUsePhiCut=1;
		fMaxPhiCut=3.125;
		break;
	default:
		AliError(Form("Max Phi Cut not defined %d",maxPhi));
		return kFALSE;
	}
	return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetDistanceToBadChannelCut(Int_t distanceToBadChannel)
{
	switch(distanceToBadChannel){
	case 0: 
		fUseDistanceToBadChannel=0;
		fMinDistanceToBadChannel=0;
		break;
	case 1: 
		if (!fUseDistanceToBadChannel) fUseDistanceToBadChannel=1;
		fMinDistanceToBadChannel=5;
		break;
	default:
		AliError(Form("minimum distance to bad channel Cut not defined %d",distanceToBadChannel));
		return kFALSE;
	}
	return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetTimingCut(Int_t timing)
{
	switch(timing){
	case 0: 
		fUseTimeDiff=0;
		fMaxTimeDiff=500;
		break;
	case 1: 
		if (!fUseTimeDiff) fUseTimeDiff=1;
		fMaxTimeDiff=10e-7; //1000ns
		break;
	case 2: 
		if (!fUseTimeDiff) fUseTimeDiff=1;
		fMaxTimeDiff=50e-8; //500ns
		break;
	case 3: 
		if (!fUseTimeDiff) fUseTimeDiff=1;
		fMaxTimeDiff=20e-8; //200ns
		break;
	case 4: 
		if (!fUseTimeDiff) fUseTimeDiff=1;
		fMaxTimeDiff=10e-8; //100ns
		break;
	case 5: 
		if (!fUseTimeDiff) fUseTimeDiff=1;
		fMaxTimeDiff=50e-9; //50ns
		break;

	default:
		AliError(Form("Timing Cut not defined %d",timing));
		return kFALSE;
	}
	return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetTrackMatchingCut(Int_t trackMatching)
{
	switch(trackMatching){
	case 0: 
		fUseDistTrackToCluster=0;
		fMinDistTrackToCluster=0;
		break;
	case 1: 
		if (!fUseDistTrackToCluster) fUseDistTrackToCluster=1;
		fMinDistTrackToCluster= 0.03;	//0.04;  
		break;
	case 2: 
		if (!fUseDistTrackToCluster) fUseDistTrackToCluster=1;
		fMinDistTrackToCluster= 0.035;	//0.05; 
		break;
	case 3: 
		if (!fUseDistTrackToCluster) fUseDistTrackToCluster=1;
		fMinDistTrackToCluster= 0.04; 	//0.1;  
		break;
	case 4: 
		if (!fUseDistTrackToCluster) fUseDistTrackToCluster=1;
		fMinDistTrackToCluster=	0.045;	//0.13; 
		break;
	case 5: 
		if (!fUseDistTrackToCluster) fUseDistTrackToCluster=1;
		fMinDistTrackToCluster= 0.05;	//0.15 
		break;
	case 6: 
		if (!fUseDistTrackToCluster) fUseDistTrackToCluster=1;
		fMinDistTrackToCluster= 0.055;	//0.2; 
		break;
	case 7: 
		if (!fUseDistTrackToCluster) fUseDistTrackToCluster=1;
		fMinDistTrackToCluster= 0.06; 	//0.3; 
		break;
	case 8: 
		if (!fUseDistTrackToCluster) fUseDistTrackToCluster=1;
		fMinDistTrackToCluster= 0.07;	//0.4; 
		break;
	case 9: 
		if (!fUseDistTrackToCluster) fUseDistTrackToCluster=1;
		fMinDistTrackToCluster= 0.1;	//0.5; 
		break;

	default:
		AliError(Form("Track Matching Cut not defined %d",trackMatching));
		return kFALSE;
	}
	return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetExoticCellCut(Int_t exoticCell)
{
	switch(exoticCell){
	case 0: 
		fUseExoticCell=0;
		fExoticCell=0;
		break;
	case 1: 
		if (!fUseExoticCell) fUseExoticCell=1;
		fExoticCell=5; 
		break;
	default:
		AliError(Form("Exotic cell Cut not defined %d",exoticCell));
		return kFALSE;
	}
	return kTRUE;
}
		
//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinEnergyCut(Int_t minEnergy)
{
	switch(minEnergy){
	case 0: 
		if (!fUseMinEnergy) fUseMinEnergy=0;
		fMinEnergy=0;
		break;
	case 1: 
		if (!fUseMinEnergy) fUseMinEnergy=1;
		fMinEnergy=0.05; 
		break;
	case 2: 
		if (!fUseMinEnergy) fUseMinEnergy=1;
		fMinEnergy=0.1; 
		break;
	case 3: 
		if (!fUseMinEnergy) fUseMinEnergy=1;
		fMinEnergy=0.15; 
		break;
	case 4: 
		if (!fUseMinEnergy) fUseMinEnergy=1;
		fMinEnergy=0.2; 
		break;
	case 5: 
		if (!fUseMinEnergy) fUseMinEnergy=1;
		fMinEnergy=0.3; 
		break;
	case 6: 
		if (!fUseMinEnergy) fUseMinEnergy=1;
		fMinEnergy=0.5; 
		break;
	case 7: 
		if (!fUseMinEnergy) fUseMinEnergy=1;
		fMinEnergy=0.75; 
		break;
	case 8: 
		if (!fUseMinEnergy) fUseMinEnergy=1;
		fMinEnergy=1.; 
		break;
	case 9: 
		if (!fUseMinEnergy) fUseMinEnergy=1;
		fMinEnergy=1.25; 
		break;
	default:
		AliError(Form("Minimum Energy Cut not defined %d",minEnergy));
		return kFALSE;
	}
	return kTRUE;
}
		
//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinNCellsCut(Int_t minNCells)
{
	switch(minNCells){
	case 0:
		if (!fUseNCells) fUseNCells=0;
		fMinNCells=0;
		break;
	case 1: 
		if (!fUseNCells) fUseNCells=1;
		fMinNCells=1; 
		break;
	case 2: 
		if (!fUseNCells) fUseNCells=1;
		fMinNCells=2; 
		break;
	case 3: 
		if (!fUseNCells) fUseNCells=1;
		fMinNCells=3; 
		break;
	case 4: 
		if (!fUseNCells) fUseNCells=1;
		fMinNCells=4; 
		break;
	case 5: 
		if (!fUseNCells) fUseNCells=1;
		fMinNCells=5; 
		break;
	case 6: 
		if (!fUseNCells) fUseNCells=1;
		fMinNCells=6; 
		break;

	default:
		AliError(Form("Min N cells Cut not defined %d",minNCells));
		return kFALSE;
	}
	return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMaxM02(Int_t maxM02)
{
	switch(maxM02){
	case 0: 
		if (!fUseM02) fUseM02=0;
		fMaxM02=100;
		break;
	case 1: 
		if (!fUseM02) fUseM02=1;
		fMaxM02=1.; 
		break;
	case 2: 
		if (!fUseM02) fUseM02=1;
		fMaxM02=0.7; 
		break;
	case 3: 
		if (!fUseM02) fUseM02=1;
		fMaxM02=0.5; 
		break;
	case 4: 
		if (!fUseM02) fUseM02=1;
		fMaxM02=0.4; 
		break;
	default:
		AliError(Form("Max M02 Cut not defined %d",maxM02));
		return kFALSE;
	}
	return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinM02(Int_t minM02)
{
	switch(minM02){
	case 0: 
		if (!fUseM02) fUseM02=0;
		fMinM02=0;
		break;
	case 1: 
		if (!fUseM02) fUseM02=1;
		fMinM02=0.002; 
		break;
	default:
		AliError(Form("Min M02 not defined %d",minM02));
		return kFALSE;
	}
	return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMaxM20(Int_t maxM20)
{
	switch(maxM20){
	case 0: 
		if (!fUseM20) fUseM20=0;
		fMaxM20=100;
		break;
	case 1: 
		if (!fUseM20) fUseM20=1;
		fMaxM20=0.5; 
		break;
	default:
		AliError(Form("Max M20 Cut not defined %d",maxM20));
		return kFALSE;
	}
	return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinM20(Int_t minM20)
{
	switch(minM20){
	case 0: 
		if (!fUseM20) fUseM20=0;
		fMinM20=0;
		break;
	case 1: 
		if (!fUseM20) fUseM20=1;
		fMinM20=0.002; 
		break;
	default:
		AliError(Form("Min M20 Cut not defined %d",minM20));
		return kFALSE;
	}
	return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetDispersion(Int_t dispersion)
{
	switch(dispersion){
	case 0: 
		if (!fUseDispersion) fUseDispersion=0;
		fMaxDispersion =100;
		break;
	case 1: 
		if (!fUseDispersion) fUseDispersion=1;
		fMaxDispersion=2.; 
		break;
	default:
		AliError(Form("Maximum Dispersion Cut not defined %d",dispersion));
		return kFALSE;
	}
	return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetNLM(Int_t nlm)
{
	switch(nlm){
	case 0: 
		if (!fUseNLM) fUseNLM=0;
		fMinNLM =0;
		fMaxNLM =100;
		break;
	case 1: 
		if (!fUseNLM) fUseNLM=1;
		fMinNLM =0;
		fMaxNLM =1;
		break;
	default:
		AliError(Form("NLM Cut not defined %d",nlm));
		return kFALSE;
	}
	return kTRUE;
}
	
///________________________________________________________________________
TString AliCaloPhotonCuts::GetCutNumber(){
   // returns TString with current cut number
   TString a(kNCuts);
   for(Int_t ii=0;ii<kNCuts;ii++){
      a.Append(Form("%d",fCuts[ii]));
   }
   return a;
}
	
	
