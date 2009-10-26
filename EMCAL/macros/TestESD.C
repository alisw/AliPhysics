#if !defined(__CINT__) || defined(__MAKECINT__)

//Root include files 
#include <Riostream.h>
#include <TFile.h>
#include <TChain.h>
#include <TParticle.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1F.h>
#include <TVector.h>
#include <TParticle.h>
#include <TRefArray.h>
#include <TArrayS.h>

//AliRoot include files 
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliPID.h"
#include "AliLog.h"

#endif

//Change the bool depending on what information you want to print
// when all FALSE, prints minimum cluster information.
Bool_t kPrintKine = kTRUE; //Do not use for raw data.
Bool_t kPrintCaloCells = kFALSE;
Bool_t kPrintTrackMatches = kFALSE;
Bool_t kPrintClusterCells = kFALSE;

void TestESD() {

  TFile* f = new TFile("AliESDs.root");
  TTree* esdTree = (TTree*)f->Get("esdTree");
  
  AliESDEvent* esd = new AliESDEvent();
  esd->ReadFromTree(esdTree);

  Int_t nEvt = esdTree->GetEntries();
  Float_t pos[3] = {0.,0.,0.};

  for(Int_t iev = 0; iev < nEvt; iev++) {
    cout << "<<<< Event: " << iev+1 << "/" << nEvt << " >>>>"<<endl;
    esdTree->GetEvent(iev);
	  
	//In case you want to play with MC data
	AliStack *stack = 0;
	if(kPrintKine){  
		AliRunLoader *rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),  "read");
		rl->LoadKinematics();  
		rl->GetEvent(iev);
		stack = rl->Stack();
	}  
  	  
	//get reconstructed vertex position
	Double_t vertex_position[3];
	esd->GetVertex()->GetXYZ(vertex_position);
	  
	//GetCellsClusters Array  
    AliESDCaloCells &cells= *(esd->GetEMCALCells());
	  // Uncomment to see the full list of digit amplitudes and times.
	if(kPrintCaloCells){  
		Int_t nTotalCells = cells.GetNumberOfCells() ;  
		Int_t type        = cells.GetType();
		for (Int_t icell=  0; icell <  nTotalCells; icell++) {
			cout<<"Cell   : "<<icell<<"/"<<nTotalCells<<" ID: "<<cells.GetCellNumber(icell)<<" Amplitude: "<<cells.GetAmplitude(icell)<<" Time: "<<cells.GetTime(icell)<<endl;	  
		}// cell loop
	}
	  
	//GetCaloClusters Array
    TRefArray* caloClusters = new TRefArray();
    esd->GetEMCALClusters(caloClusters);

    //loop over clusters
    Int_t nclus = caloClusters->GetEntries();
    for (Int_t icl = 0; icl < nclus; icl++) {

      AliESDCaloCluster* clus = (AliESDCaloCluster*)caloClusters->At(icl);
      Float_t energy = clus->E();
      clus->GetPosition(pos);
      TVector3 vpos(pos[0],pos[1],pos[2]);
      TLorentzVector p;
      clus->GetMomentum(p,vertex_position);
      Double_t cphi = vpos.Phi();
      Double_t ceta = vpos.Eta();

      Int_t nMatched = clus->GetNTracksMatched();
      Int_t trackIndex = clus->GetTrackMatched();
      Int_t nLabels = clus->GetNLabels();
      Int_t labelIndex = clus->GetLabel();

      Int_t nCells = clus->GetNCells();

      //For later: ADD CHECK THAT CLUSTER IS WITHIN SM FIDUCIAL VOLUME

      cout << "Cluster: " << icl+1 << "/" << nclus << " Energy: " << energy << " Phi: " 
	   << cphi << " Eta: " << ceta << " NCells: " << nCells << " #Matches: " << nMatched 
	   << " Index: " << trackIndex << " #Labels: " << nLabels << " Index: " 
	   << labelIndex << endl;

	if(kPrintTrackMatches && trackIndex >= 0) {
		AliESDtrack* track = esd->GetTrack(trackIndex);
		Double_t tphi = track->GetOuterParam()->Phi();
		Double_t teta = track->GetOuterParam()->Eta();
		Double_t tmom = track->GetOuterParam()->P();
		cout << "\t Track Momentum: " << tmom << " phi: " << tphi << " eta: " << teta << endl;

		Double_t deta = teta - ceta;
		Double_t dphi = tphi - cphi;
		if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
		if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
		Double_t dR = sqrt(dphi*dphi + deta*deta);

		Double_t pOverE = tmom/energy;

		if(dR < 0.02 && pOverE < 1.8 && nCells > 1) {
			cout << "\n\t Excellent MATCH! dR = " << dR << " p/E = " << pOverE << " nCells = " << nCells << endl;
		}

	}
		
	//Get CaloCells of cluster and print
	if(kPrintClusterCells){	
		UShort_t * index = clus->GetCellsAbsId() ;
		Double_t * fraction = clus->GetCellsAmplitudeFraction() ;
		for(Int_t i = 0; i < nCells ; i++){
			Int_t absId =   index[i]; // or clus->GetCellNumber(i) ;
			Double_t ampFract =  fraction[i];
			Float_t amp       = cells.GetCellAmplitude(absId) ;
			Float_t time      = cells.GetCellTime(absId);
			cout<<"         Cluster Cell: AbsID : "<< absId << "; Amplitude "<< amp << "; Fraction "<<ampFract<<"; Time " <<time<<endl;
		}
	}
		
	//Print primary info
	if(!stack || !kPrintKine) continue;
	if(labelIndex >= 0 && labelIndex < stack->GetNtrack()){
		TParticle * particle = stack->Particle(labelIndex);
		//Print primary values
		cout<<"         More  contributing primary: "<<particle->GetName()<< "; Energy "<<particle->Energy()<<endl;   
		for(Int_t i = 1; i < nLabels; i++){
			particle = stack->Particle(clus->(GetLabels()->At(i)));
			cout<<"         Other contributing primary: "<<particle->GetName()<< "; Energy "<<particle->Energy()<<endl;
		}
	}
	else if( labelIndex >= stack->GetNtrack()) cout <<"PROBLEM, label is too large : "<<labelIndex<<" >= particles in stack "<< stack->GetNtrack() <<endl;
	else cout<<"Negative label!!!  : "<<labelIndex<<endl;
		
	}


  }



}
