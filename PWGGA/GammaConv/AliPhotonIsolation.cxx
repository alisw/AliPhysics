
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 **
 * Authors: Helena Horstmann, Daniel Mühlheim*
 * Version 1.0*
 **
 * Permission to use, copy, modify and distribute this software and its *
 * documentation strictly for non-commercial purposes is hereby granted *
 * without fee, provided that the above copyright notice appears in all *
 * copies and that both the copyright notice and this permission notice *
 * appear in the supporting documentation. The authors make no claims*
 * about the suitability of this software for any purpose. It is*
 * provided "as is" without express or implied warranty.*
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Basic Isolation Class
//---------------------------------------------
////////////////////////////////////////////////


#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliPhotonIsolation.h"
#include "AliV0ReaderV1.h"
#include "AliVParticle.h"
#include "TChain.h"
#include "TH1F.h"
#include "TAxis.h"

#include <vector>
#include <map>
#include <utility>


class iostream;

using namespace std;


ClassImp(AliPhotonIsolation)

//_______________________________________________________________________
AliPhotonIsolation::AliPhotonIsolation(const char *name, Int_t photonType) : AliAnalysisTaskSE(name),

  fPhotonType(photonType),
  fV0ReaderName(""),
  fCorrTaskSetting(""),
  fMapClustertoPtR1(),
  fMapClustertoPtR2(),
  fMapClustertoPtR3(),
  fMapClustertoPtR4(),
  fListHistos(NULL),
  fHistTest(NULL),
  fHistClusterEnergy(NULL),
  fHistIso(NULL)
{
  // Default constructor
  DefineInput(0, TChain::Class());
}


//________________________________________________________________________
AliPhotonIsolation::~AliPhotonIsolation(){
  // default deconstructor
  fMapClustertoPtR1.clear();
  fMapClustertoPtR2.clear();
  fMapClustertoPtR3.clear();
  fMapClustertoPtR4.clear();

  if(fHistIso) delete fHistIso;
  if(fHistTest) delete fHistTest;
  if(fHistClusterEnergy) delete fHistClusterEnergy;
  if(fListHistos != NULL){
    delete fListHistos;
    fListHistos = NULL; //#################################################################### Why?
  }
}

//________________________________________________________________________
void AliPhotonIsolation::Terminate(Option_t *){
  fMapClustertoPtR1.clear();
  fMapClustertoPtR2.clear();
  fMapClustertoPtR3.clear();
  fMapClustertoPtR4.clear();
}

//________________________________________________________________________
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ OUTPUT
void AliPhotonIsolation::UserCreateOutputObjects(){
  if(fListHistos != NULL){
    delete fListHistos;
    fListHistos = NULL;
  }
  if(fListHistos == NULL){
    fListHistos = new TList();
    fListHistos->SetOwner(kTRUE);
    fListHistos->SetName(Form("PhotonIsolation_%i",fPhotonType));
  }

  //Create user output objects
  fHistIso = new TH1F("fHistM02","M02 Distribution",200,0,50);
  fHistTest = new TH1F("fHistClus","Number of clusters per Event",20,0,20);
  fHistClusterEnergy = new TH1F("fHistEnergy","Pt Spectrum",280,0,140);
  fListHistos->Add(fHistTest);
  fListHistos->Add(fHistIso);
  fListHistos->Add(fHistClusterEnergy);

}

//________________________________________________________________________
void AliPhotonIsolation::Initialize(){
  // Initialize function to be called once before analysis
  fMapClustertoPtR1.clear();
  fMapClustertoPtR2.clear();
  fMapClustertoPtR3.clear();
  fMapClustertoPtR4.clear();
}

//________________________________________________________________________
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ USEREXEC
void AliPhotonIsolation::UserExec(Option_t *){

  // main method of AliPhotonIsolation, first initialize and then process event

  Initialize();

  ProcessEvent(fInputEvent);

  //DebugIsolation();

  return;
}

//________________________________________________________________________
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PROCESS EVENT
void AliPhotonIsolation::ProcessEvent(AliVEvent *event){
  Bool_t debug = 0; //Debugging on/off
  //###########################################################GET NUMBER OF CLUSTERS IN EVENT
  Int_t nClus = 0; 
  TClonesArray * arrClusters = NULL; 
  if(!fCorrTaskSetting.CompareTo("")){
    nClus = event->GetNumberOfCaloClusters();
  } else {
    arrClusters = dynamic_cast<TClonesArray*>(event->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    nClus = arrClusters->GetEntries();
  }

  //#####################################################CREATE CLASS TO GET EVENT INFORMATION
  AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodev = 0;
  if (!esdev) {
    aodev = dynamic_cast<AliAODEvent*>(event);
    if (!aodev) {
      AliError("Task needs AOD or ESD event, returning");
      return;
    }
  }
  
  if(debug){
    cout << "___________________________________________________This is the start of an event___________________________________" << endl;
  }

  fHistTest->Fill(nClus);

  //#########################################################################LOOP OVER CLUSTERS
  for(Int_t iclus=0;iclus < nClus;iclus++){ 

      if(debug){
        cout << "________________________________̣̣This is a new cluster_________________________________" << endl;
      }

    Float_t ptsum1 = 0., ptsum2 = 0., ptsum3 = 0., ptsum4 = 0.;
    //#####################################################INITIALISE AND CREATE CLUSTER OBJECT
    AliVCluster* cluster = NULL; 
    if(arrClusters){
      if(esdev){
	if(arrClusters)
	  cluster = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClusters->At(iclus));
      } else if(aodev){
	if(arrClusters)
	  cluster = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClusters->At(iclus));
      }
    }
    else
      cluster = event->GetCaloCluster(iclus); //create ALiVCluster object
    if (!cluster){ continue;}
    //####################################################DETERMINE ETA, PHI AND E_T FOR THE CLUSTER
    Float_t eta_clus = 0., phi_clus = 0., theta_clus = 0; //define variables for cluster properties
    Float_t clspos[3] = {0.,0.,0.}; //initialise position vector

    cluster->GetPosition(clspos); //get position of cluster
    TVector3 clspos_v3(clspos); //create 3Vector to gain access to geometrical functions
    eta_clus = clspos_v3.Eta(); //determine Pseudo rapidity
    phi_clus = clspos_v3.Phi(); //determine Phi
    theta_clus = clspos_v3.Theta(); //determine Theta in radians

    Float_t ET_clus = (cluster->E())*TMath::Sin(theta_clus); //Determine cluster energy (E_T=E*sin(theta), because m=0 or negligible)

    Float_t m02 = cluster->GetM02();
    //    cout << "M02: " << m02 << endl;

    fHistIso->Fill(m02);

    if(debug){
      cout << "Clusterenergie:   " << ET_clus << endl;
    }

    fHistClusterEnergy->Fill(ET_clus);
      
      //#########################################################################LOOP OVER TRACKS
      for (Int_t itr=0;itr<event->GetNumberOfTracks();itr++){ 
        if(debug){
          cout << "______________________________This is a new track_______________________" << endl;
        }
	
	//#####################################################INITIALISE AND CREATE TRACK OBJECT
	AliVTrack *inTrack = 0x0; //initialise general track variable
	if(esdev){
	  inTrack = esdev->GetTrack(itr); //get track
	  if(!inTrack) continue; //test track
	  //  AliESDtrack *esdt = dynamic_cast<AliESDtrack*>(inTrack); //create class with track information      
	} else if(aodev) { //same for AODs
      inTrack = dynamic_cast<AliVTrack*>(aodev->GetTrack(itr));
      if(!inTrack) continue;
	  //  AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(inTrack);
	}
	 
	//####################################################DETERMINE ETA AND PHI FOR THE TRACK
	Float_t eta_track = 0., phi_track = 0.; //define variables for track properties
	eta_track = inTrack->Eta(); 
	phi_track = inTrack->Phi();
	
	//#####################################CHECK THE DISTANCE/ANGLE BETWEEN TRACK AND CLUSTER
	Float_t R2 = (TMath::Power(eta_clus-eta_track,2)+TMath::Power(phi_clus-phi_track,2)); //calculate the squared radius      
/*
    if(debug){
      cout << "PhiTrack:  " << phi_track <<  endl;
      cout << "EtaTrack:  " << eta_track <<  endl;
      cout << "PhiClus:   " << phi_clus <<  endl;
      cout << "EtaClus:   " << eta_clus <<  endl;
      cout << "R:         " << R2 <<  endl;
    }
*/
/*
    if(debug){
      cout <<  << endl;
    }
*/

	Float_t pT = inTrack->Pt();//get pT of track

	//add up pT values of tracks in cone
	if(R2<0.01){//compare to R=0.1
	  ptsum1 += pT;
	}else if(R2<0.04){//compare to R=0.2
	  ptsum2 += pT;
	}else if(R2<0.09){//compare to R=0.3
	  ptsum3 += pT;
	}else if(R2<0.16){//compare to R=0.4
	  ptsum4 += pT;
	}
	
    if(debug && R2<0.16){
      cout << "R2:     " << R2 << endl;
    }

      } //end of for loop tracks

      ptsum2+=ptsum1; //add pT of innerst cone to second cone
      ptsum3+=ptsum2; //add complete second to third
      ptsum4+=ptsum3; //add complete third to fourth
      
      if(debug){
        cout << "ptsum1:  " << ptsum1 << endl;
        cout << "ptsum2:  " << ptsum2 << endl;
        cout << "ptsum3:  " << ptsum3 << endl;
        cout << "ptsum4:  " << ptsum4 << endl;
      }

      //map Ptsum to clusterID, map all the ptsum that are bigger than 0
      if(ptsum1>0){
	fMapClustertoPtR1.insert(make_pair(cluster->GetID(),ptsum1));
	fMapClustertoPtR2.insert(make_pair(cluster->GetID(),ptsum2));
	fMapClustertoPtR3.insert(make_pair(cluster->GetID(),ptsum3));
	fMapClustertoPtR4.insert(make_pair(cluster->GetID(),ptsum4));
      }else if(ptsum2>0){
	fMapClustertoPtR2.insert(make_pair(cluster->GetID(),ptsum2));
	fMapClustertoPtR3.insert(make_pair(cluster->GetID(),ptsum3));
	fMapClustertoPtR4.insert(make_pair(cluster->GetID(),ptsum4));
      }else if(ptsum3>0){
	fMapClustertoPtR3.insert(make_pair(cluster->GetID(),ptsum3));
	fMapClustertoPtR4.insert(make_pair(cluster->GetID(),ptsum4));	    
      }else if(ptsum4>0){
	fMapClustertoPtR4.insert(make_pair(cluster->GetID(),ptsum4));
      }

      if(debug){
        cout << "Cluster ID: " << iclus << endl;
        cout << "ptsumR4: " << ptsum4 << endl;
      }

  }//end for loop clusters

  if(debug){
    for (map<Int_t,Float_t>::iterator it= fMapClustertoPtR4.begin(); it != fMapClustertoPtR4.end(); ++it) {
      cout << it->first << "\t" << it->second << endl;
    }
  }

  return;
}


//Determine Isolation dependent on ClusterID, Isolation Cone and Isolation Energy

Bool_t AliPhotonIsolation::GetIsolation(Int_t clusterID, Float_t R, Float_t isoPt){

  Bool_t isolated = kFALSE;
  if (R>0.09 && R<0.11){ //check Radius
    if(fMapClustertoPtR1[clusterID]<isoPt) isolated = kTRUE;
    //cout << "MapPt" << fMapClustertoPtR1[clusterID] << endl;
  }else if(R>0.19 && R<0.21){ //check Radius
    if(fMapClustertoPtR2[clusterID]<isoPt) isolated = kTRUE;
  }else if(R>0.29 && R<0.31){ //check Radius
    if(fMapClustertoPtR3[clusterID]<isoPt) isolated = kTRUE;
  }else if(R>0.39 && R<0.41){ //check Radius
    if(fMapClustertoPtR4[clusterID]<isoPt) isolated = kTRUE;
  }else{
    AliFatal("AliPhotonIsolation: GetIsolation - Radius of the Isolation Cone has to be set to R=0.1, R=0.2, R=0.3 or R=0.4");
  }

 // cout << "_______________________________________________________________________Here I am to check out isolated Photons!" << endl;

  return isolated;
}


void AliPhotonIsolation::DebugIsolation(){

  return;
}


//
