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

#include "TObjArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TParticle.h"
#include "TList.h"
#include "THashList.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskEpRatio.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliVCaloCells.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliPIDResponse.h" 

// *********************** New
#include "AliVertex.h"
#include "AliVVertex.h"
#include "AliVertexerTracks.h"
#include "AliExternalTrackParam.h"
#include "AliTrackReference.h"

#include <AliCDBManager.h>
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliPHOSCalibData.h"
#include "AliOADBContainer.h"

// Analysis task to fill histograms with E/p ratios for electrons, positrons, muons
// and hadrons (pions, kaons..), where E is an energy of PHOS cluster
// and p is the momentum of it's matched track.

// Authors: Boris Polishchuk,Tsubasa Okubo
// Date   : 04.07.2013

ClassImp(AliAnalysisTaskEpRatio)
//________________________________________________________________________
AliAnalysisTaskEpRatio::AliAnalysisTaskEpRatio(const char *name) :  AliAnalysisTaskSE(name),
  fRunNumber(-999),
  fOutputContainer(0),
  fPHOSGeo(0),
  fPIDResponse(0) // Yuri
// Output slots #0 write into a TH1 container
{
  DefineOutput(1,TList::Class());
  
}
//________________________________________________________________________
void AliAnalysisTaskEpRatio::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner();
  
  fOutputContainer->Add( new TH1F("h_CellMultEvent","PHOS Cell Multiplicity per Event",200,0,200) );
  
  fOutputContainer->Add( new TH2F("hEp","E/p All charged particles (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_ele","E/p electron & positron (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_ele2","E/p electron & positron (E/p VS. 1/p)",1000,0,2.,500,0,5) );
  fOutputContainer->Add( new TH2F("hEp_other","E/p other (E/p VS. p_{T})",1000,0,2,1000,0,10.) );

  fOutputContainer->Add( new TH2F("hEp_electron","E/p Electron only (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_positron","E/p Positron only (E/p VS. p_{T})",1000,0,2,1000,0,10.) );

  fOutputContainer->Add( new TH2F("hEp_mod1","E/p All charged particles in Module1 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_mod2","E/p All charged particles in Module2 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_mod3","E/p All charged particles in Module3 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_mod4","E/p All charged particles in Module4 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );

  fOutputContainer->Add( new TH2F("hEp_ele_mod1","E/p electron & positron in Module1 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_ele_mod2","E/p electron & positron in Module2 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_ele_mod3","E/p electron & positron in Module3 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_ele_mod4","E/p electron & positron in Module4 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );

  fOutputContainer->Add( new TH2F("hEp_other_mod1","E/p other in Module1 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_other_mod2","E/p other in Module2 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_other_mod3","E/p other in Module3 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_other_mod4","E/p other in Module4 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );

  fOutputContainer->Add( new TH2F("h_dedx_PHOS","dE/dx in PHOS via TPC",1000,0,10,2200,0,220.) );
  fOutputContainer->Add( new TH2F("h_dedx_ele_PHOS","dE/dx electron in PHOS via TPC",1000,0,10,2200,0,220.) );
  fOutputContainer->Add( new TH2F("h_dedx_other_PHOS","dE/dx other in PHOS via TPC",1000,0,10,2200,0,220.) );
  fOutputContainer->Add( new TH2F("h_dedx_pion_PHOS"  ,"dE/dx pion in PHOS via TPC",1000,0,10,2200,0,220.) );
  fOutputContainer->Add( new TH2F("h_dedx_proton_PHOS","dE/dx proton in PHOS via TPC",1000,0,10,2200,0,220.) );
  fOutputContainer->Add( new TH2F("h_dedx_kaon_PHOS"  ,"dE/dx kaon in PHOS via TPC",1000,0,10,2200,0,220.) );
  fOutputContainer->Add( new TH2F("h_dedx_muon_PHOS"  ,"dE/dx muon in PHOS via TPC",1000,0,10,2200,0,220.) );

  fOutputContainer->Add( new TH2F("h_ClusNXZM1","Cluster (X,Z) Module1",   64,0.5,64.5, 56,0.5,56.5) );
  fOutputContainer->Add( new TH2F("h_ClusNXZM2","Cluster (X,Z) Module2",   64,0.5,64.5, 56,0.5,56.5) );
  fOutputContainer->Add( new TH2F("h_ClusNXZM3","Cluster (X,Z) Module3",   64,0.5,64.5, 56,0.5,56.5) );
  fOutputContainer->Add( new TH2F("h_ClusNXZM4","Cluster (X,Z) Module4",   64,0.5,64.5, 56,0.5,56.5) );
  fOutputContainer->Add( new TH2F("h_ClusNXZM1_ele","Electron Cluster (X,Z) Module1", 64,0.5,64.5, 56,0.5,56.5) );
  fOutputContainer->Add( new TH2F("h_ClusNXZM2_ele","Electron Cluster (X,Z) Module2", 64,0.5,64.5, 56,0.5,56.5) );
  fOutputContainer->Add( new TH2F("h_ClusNXZM3_ele","Electron Cluster (X,Z) Module3", 64,0.5,64.5, 56,0.5,56.5) );
  fOutputContainer->Add( new TH2F("h_ClusNXZM4_ele","Electron Cluster (X,Z) Module4", 64,0.5,64.5, 56,0.5,56.5) );

  fOutputContainer->Add( new TH2F("h_EClus_NCell","Energy VS. NCell in PHOS",1000,0,10,30,0,30.) );
  fOutputContainer->Add( new TH1F("h_ECluster","Energy of Cluster in PHOS",1000,0,10.) );
  fOutputContainer->Add( new TH1F("h_NCells","Number of Cells per Cluster in PHOS",30,0,30.) );

  // =====================================================================================================
  fOutputContainer->Add( new TH2F("hEp_had","E/p Hadrons (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_had2","E/p Hadrons (E/p VS. 1/p)",1000,0,2.,500,0,5) );
  fOutputContainer->Add( new TH2F("h_dedx_had_PHOS","dE/dx Hadron in PHOS via TPC",1000,0,10,2200,0,220.) );
  fOutputContainer->Add( new TH2F("h_ClusNXZM1_had","Hadron Cluster (X,Z) Module1", 64,0.5,64.5, 56,0.5,56.5) );
  fOutputContainer->Add( new TH2F("h_ClusNXZM2_had","Hadron Cluster (X,Z) Module2", 64,0.5,64.5, 56,0.5,56.5) );
  fOutputContainer->Add( new TH2F("h_ClusNXZM3_had","Hadron Cluster (X,Z) Module3", 64,0.5,64.5, 56,0.5,56.5) );
  fOutputContainer->Add( new TH2F("h_ClusNXZM4_had","Hadron Cluster (X,Z) Module4", 64,0.5,64.5, 56,0.5,56.5) );
  fOutputContainer->Add( new TH2F("hEp_had_mod1","E/p Hadron in Module1 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_had_mod2","E/p Hadron in Module2 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_had_mod3","E/p Hadron in Module3 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );
  fOutputContainer->Add( new TH2F("hEp_had_mod4","E/p Hadron in Module4 (E/p VS. p_{T})",1000,0,2,1000,0,10.) );

  // =====================================================================================================
  fOutputContainer->Add( new TH1F("h_NClusEvent","Number PHOS Clusters per Event in All Modules",31,-0.5,30.5) );
  fOutputContainer->Add( new TH1F("h_NClusEventmod1","Number PHOS Clusters per Event in Module1",21,-0.5,20.5) );
  fOutputContainer->Add( new TH1F("h_NClusEventmod2","Number PHOS Clusters per Event in Module2",21,-0.5,20.5) );
  fOutputContainer->Add( new TH1F("h_NClusEventmod3","Number PHOS Clusters per Event in Module3",21,-0.5,20.5) );
  fOutputContainer->Add( new TH1F("h_NClusEventmod4","Number PHOS Clusters per Event in Module4",21,-0.5,20.5) );

  fOutputContainer->Add( new TH1F("h_EClusEvent","PHOS Cluster Energy per Event in All Modules",2010,-0.5,20.5) );
  fOutputContainer->Add( new TH1F("h_EClusEventmod1","PHOS Cluster Energy per Event in Module1",2010,-0.5,20.5) );
  fOutputContainer->Add( new TH1F("h_EClusEventmod2","PHOS Cluster Energy per Event in Module2",2010,-0.5,20.5) );
  fOutputContainer->Add( new TH1F("h_EClusEventmod3","PHOS Cluster Energy per Event in Module3",2010,-0.5,20.5) );
  fOutputContainer->Add( new TH1F("h_EClusEventmod4","PHOS Cluster Energy per Event in Module4",2010,-0.5,20.5) );

  fOutputContainer->Add( new TH1F("h_NTrackEvent","Number of TPC Tracks per Event",1001,-0.5,1000.5) );
  fOutputContainer->Add( new TH1F("h_PTrackEvent","TPC Momentum per Event",2010,-0.5,20.5) );
  fOutputContainer->Add( new TH1F("h_NClusTPCEvent","Number of TPC Clusters per Event",301,-0.5,300.5) );
  
  // ======================================================================================================
  
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());

  fPIDResponse = inputHandler->GetPIDResponse();
  if(!fPIDResponse) AliFatal(" !!! FATAL!! No PIDResponse!!");

  PostData(1, fOutputContainer);
  
}

//________________________________________________________________________
void AliAnalysisTaskEpRatio::UserExec(Option_t *) 
{
  //Main function that does all the job.

  AliVEvent* event = InputEvent();
  fRunNumber = event->GetRunNumber(); 

  SetGeometry();
  
  Double_t v[3]={0,0,0}; //vertex;  
  const AliVVertex *primvertex = event->GetPrimaryVertex();

  primvertex->GetXYZ(v);
  TVector3 vtx(v);
  
  if(fabs(v[2])>10.) return;
  
  const Int_t kNumberOfTracks = event->GetNumberOfTracks();
  const Int_t kNumberOfCaloClusters = event->GetNumberOfCaloClusters();
  
  AliVCaloCells *cells = event->GetPHOSCells();
  Int_t multCells = cells->GetNumberOfCells();
  FillHistogram("h_CellMultEvent",multCells);
  
  TLorentzVector pc1; //4-momentum of PHOS cluster 1
  Int_t inPHOS=0;
  
  // ================
  // === TPC dedx ===
  // ================
  FillHistogram("h_NTrackEvent",kNumberOfTracks);

  for(Int_t iTrack=0; iTrack<kNumberOfTracks; iTrack++){

    AliVTrack* esdTrack = dynamic_cast<AliVTrack*> (event->GetTrack(iTrack));
    if(!esdTrack) AliFatal(Form("iTrack %d is Not a AliVTrack!!",iTrack));

    FillHistogram("h_PTrackEvent",esdTrack->P() );
    FillHistogram("h_NClusTPCEvent",esdTrack->GetTPCNcls() );
  }

  // ================
  // === E/p ========
  // ================
  Float_t position[3];
  Int_t mod1, relId[4], cellAbsId, cellX, cellZ;
  Int_t nPHOSCluster = 0, inMod1 = 0, inMod2 = 0, inMod3 = 0,inMod4 = 0 ;
  
  for(Int_t iClu=0; iClu<kNumberOfCaloClusters; iClu++){
    
    AliVCluster *c1 = event->GetCaloCluster(iClu);

    if( !c1->IsPHOS() ) continue;
    if ( c1->GetType() !=AliVCluster::kPHOSNeutral ) continue; // reject CPV clusters
    
    //if( c1->E()<0.3 ) continue;
    FillHistogram("h_EClus_NCell",c1->E(),c1->GetNCells());
    FillHistogram("h_ECluster",c1->E());
    FillHistogram("h_NCells",c1->GetNCells());
    if( c1->GetNCells()<3 ) continue;

    c1->GetPosition(position);
    TVector3 global1(position);
    fPHOSGeo->GlobalPos2RelId(global1,relId);
    mod1  = relId[0];
    cellX = relId[2];
    cellZ = relId[3];

    cellAbsId = c1->GetCellAbsId(0);
    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    mod1 = relId[0];
    
    FillHistogram("h_EClusEvent",c1->E());
    
    if( mod1 == 1 ){
      inMod1++; nPHOSCluster++;
      FillHistogram("h_EClusEventmod1",c1->E());
    }
    else if( mod1 == 2 ){
      inMod2++; nPHOSCluster++;
      FillHistogram("h_EClusEventmod2",c1->E());
    }
    else if( mod1 == 3 ){
      inMod3++; nPHOSCluster++;
      FillHistogram("h_EClusEventmod3",c1->E());
    }
    else if( mod1 == 4 ){
      inMod4++; nPHOSCluster++;
      FillHistogram("h_EClusEventmod4",c1->E());
    }

    Int_t nMatched = c1->GetNTracksMatched();
    
    Double_t Dx = c1->GetTrackDx();
    Double_t Dz = c1->GetTrackDz();
    Double_t r = sqrt(Dx*Dx+Dz*Dz); // unit is [cm]
    
    if(r > 5.) continue;
    if(nMatched<=0) continue;
    
    AliVTrack* esdTrack = dynamic_cast<AliVTrack*> (c1->GetTrackMatched(0));
    // track was not found
    if(!esdTrack) continue;

    if( !(TMath::Abs(esdTrack->Eta())< 0.8) ) continue;

    Short_t charge   = esdTrack->Charge();
    if(fPIDResponse) {

      //if(c1->E()<0.3) continue;  // MIP & noise cut

      // =============================================================================
      FillHistogram("hEp",c1->E()/esdTrack->P(),esdTrack->Pt());
      FillHistogram("h_dedx_PHOS",esdTrack->P(),esdTrack->GetTPCsignal());

      if( mod1 == 1 ){
	FillHistogram("hEp_mod1",c1->E()/esdTrack->P(),esdTrack->Pt());
	FillHistogram("h_ClusNXZM1",cellX,cellZ,1.);
      }
      else if( mod1 == 2 ){
	FillHistogram("hEp_mod2",c1->E()/esdTrack->P(),esdTrack->Pt());
	FillHistogram("h_ClusNXZM2",cellX,cellZ,1.);
      }
      else if( mod1 == 3 ){
        FillHistogram("hEp_mod3",c1->E()/esdTrack->P(),esdTrack->Pt());
	FillHistogram("h_ClusNXZM3",cellX,cellZ,1.);
      }
      else if( mod1 == 4 ){
	FillHistogram("hEp_mod4",c1->E()/esdTrack->P(),esdTrack->Pt());
	FillHistogram("h_ClusNXZM4",cellX,cellZ,1.);
      }


      Bool_t pidPion=kFALSE, pidKaon=kFALSE, pidProton=kFALSE, pidElectron=kFALSE, pidMuon=kFALSE, pidHadron=kFALSE;
      Double_t nsigmaProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack, AliPID::kProton)); // esdTrack << inEvHMain
      Double_t nsigmaKaon     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack, AliPID::kKaon)); 
      Double_t nsigmaPion     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack, AliPID::kPion)); 
      //Double_t nsigmaElectron = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack, AliPID::kElectron));
      Double_t nsigmaElectron = fPIDResponse->NumberOfSigmasTPC(esdTrack, AliPID::kElectron);
      Double_t nsigmaMuon     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack, AliPID::kMuon));
     
      if((nsigmaKaon    <nsigmaPion) && (nsigmaKaon    <nsigmaProton) && (nsigmaKaon    <nsigmaElectron) 
	 && (nsigmaKaon     <nsigmaMuon    ) && (nsigmaKaon    <3)) pidKaon     = kTRUE;
      if((nsigmaPion    <nsigmaKaon) && (nsigmaPion    <nsigmaProton) && (nsigmaPion    <nsigmaElectron)
	 && (nsigmaPion     <nsigmaMuon    ) && (nsigmaPion    <3)) pidPion     = kTRUE;
      if((nsigmaProton  <nsigmaKaon) && (nsigmaProton  <nsigmaPion  ) && (nsigmaProton  <nsigmaElectron)
	 && (nsigmaProton   <nsigmaMuon    ) && (nsigmaProton  <3)) pidProton   = kTRUE;
      if((nsigmaMuon    <nsigmaKaon) && (nsigmaMuon    <nsigmaPion  ) && (nsigmaMuon    <nsigmaProton  )
	 && (nsigmaMuon     <nsigmaElectron) && (nsigmaMuon    <3)) pidMuon     = kTRUE;
      if((nsigmaElectron<nsigmaKaon) && (nsigmaElectron<nsigmaPion  ) && (nsigmaElectron<nsigmaProton  ) 
	 && (nsigmaElectron <nsigmaMuon    ) && (nsigmaElectron<3) && (nsigmaElectron>-2)) pidElectron = kTRUE;
      if( (nsigmaElectron<-3) ) pidHadron = kTRUE;
      
      if (pidPion){
	FillHistogram("h_dedx_pion_PHOS",esdTrack->P(),esdTrack->GetTPCsignal());
      }
      if (pidKaon){
	FillHistogram("h_dedx_kaon_PHOS",esdTrack->P(),esdTrack->GetTPCsignal());
      }
      if (pidProton){
	FillHistogram("h_dedx_proton_PHOS",esdTrack->P(),esdTrack->GetTPCsignal());
      }
      if (pidMuon){
        FillHistogram("h_dedx_muon_PHOS",esdTrack->P(),esdTrack->GetTPCsignal());
      }

      if (pidElectron){
	FillHistogram("hEp_ele",c1->E()/esdTrack->P(),esdTrack->Pt());
        FillHistogram("hEp_ele2",c1->E()/esdTrack->P(),1/esdTrack->P());
	FillHistogram("h_dedx_ele_PHOS",esdTrack->P(),esdTrack->GetTPCsignal());

	// =============================================================================

	if (charge>0) FillHistogram("hEp_positron",c1->E()/esdTrack->P(),esdTrack->Pt());
        if (charge<0) FillHistogram("hEp_electron",c1->E()/esdTrack->P(),esdTrack->Pt());

	if( mod1 == 1 ){
	  FillHistogram("h_ClusNXZM1_ele",cellX,cellZ,1.);
	  FillHistogram("hEp_ele_mod1",c1->E()/esdTrack->P(),esdTrack->Pt());
	}
	else if( mod1 == 2 ){
	  FillHistogram("h_ClusNXZM2_ele",cellX,cellZ,1.);
	  FillHistogram("hEp_ele_mod2",c1->E()/esdTrack->P(),esdTrack->Pt());
	}
	else if( mod1 == 3 ){
	  FillHistogram("h_ClusNXZM3_ele",cellX,cellZ,1.);
	  FillHistogram("hEp_ele_mod3",c1->E()/esdTrack->P(),esdTrack->Pt());
	}
	else if( mod1 == 4 ){
	  FillHistogram("h_ClusNXZM4_ele",cellX,cellZ,1.);
	  FillHistogram("hEp_ele_mod4",c1->E()/esdTrack->P(),esdTrack->Pt());
	}
      }
      else{
	FillHistogram("hEp_other",c1->E()/esdTrack->P(),esdTrack->Pt());
	FillHistogram("h_dedx_other_PHOS",esdTrack->P(),esdTrack->GetTPCsignal());

	if( mod1 == 1 ){
	  FillHistogram("hEp_other_mod1",c1->E()/esdTrack->P(),esdTrack->Pt());
	}
	else if( mod1 == 2 ){
	  FillHistogram("hEp_other_mod2",c1->E()/esdTrack->P(),esdTrack->Pt());
	}
	else if( mod1 == 3 ){
	  FillHistogram("hEp_other_mod3",c1->E()/esdTrack->P(),esdTrack->Pt());
	}
	else if( mod1 == 4 ){
	  FillHistogram("hEp_other_mod4",c1->E()/esdTrack->P(),esdTrack->Pt());
	}

      }

      if (pidHadron){
	FillHistogram("hEp_had",c1->E()/esdTrack->P(),esdTrack->Pt());
	FillHistogram("hEp_had2",c1->E()/esdTrack->P(),1/esdTrack->P());
	FillHistogram("h_dedx_had_PHOS",esdTrack->P(),esdTrack->GetTPCsignal());

	// =============================================================================

	if( mod1 == 1 ){
	  FillHistogram("h_ClusNXZM1_had",cellX,cellZ,1.);
	  FillHistogram("hEp_had_mod1",c1->E()/esdTrack->P(),esdTrack->Pt());
	}
	else if( mod1 == 2 ){
	  FillHistogram("h_ClusNXZM2_had",cellX,cellZ,1.);
	  FillHistogram("hEp_had_mod2",c1->E()/esdTrack->P(),esdTrack->Pt());
	}
	else if( mod1 == 3 ){
	  FillHistogram("h_ClusNXZM3_had",cellX,cellZ,1.);
	  FillHistogram("hEp_had_mod3",c1->E()/esdTrack->P(),esdTrack->Pt());
	}
	else if( mod1 == 4 ){
	  FillHistogram("h_ClusNXZM4_had",cellX,cellZ,1.);
	  FillHistogram("hEp_had_mod4",c1->E()/esdTrack->P(),esdTrack->Pt());
	}
      }
      
    }
    
    // =============================================================================================================================
    inPHOS++;
  }
  
  FillHistogram("h_NClusEvent",nPHOSCluster);
  FillHistogram("h_NClusEventmod1",inMod1);
  FillHistogram("h_NClusEventmod2",inMod2);
  FillHistogram("h_NClusEventmod3",inMod3);
  FillHistogram("h_NClusEventmod4",inMod4);

  nPHOSCluster = 0;
  inMod1 = 0;
  inMod2 = 0;
  inMod3 = 0;
  inMod4 = 0;
 
  //delete caloClustersArr;
  PostData(1,fOutputContainer);
  return ;
}

//_____________________________________________________________________________
void AliAnalysisTaskEpRatio::FillHistogram(const char * key,Double_t x)const{
  //FillHistogram
  TH1I * tmpI = dynamic_cast<TH1I*>(fOutputContainer->FindObject(key)) ;
  if(tmpI){
    tmpI->Fill(x) ;
    return ;
  }
  TH1F * tmpF = dynamic_cast<TH1F*>(fOutputContainer->FindObject(key)) ;
  if(tmpF){
    tmpF->Fill(x) ;
    return ;
  }
  TH1D * tmpD = dynamic_cast<TH1D*>(fOutputContainer->FindObject(key)) ;
  if(tmpD){
    tmpD->Fill(x) ;
    return ;
  }
  AliInfo(Form("can not find histogram <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskEpRatio::FillHistogram(const char * key,Double_t x,Double_t y)const{
  //FillHistogram
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH1F")){
    ((TH1F*)tmp)->Fill(x,y) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y) ;
    return ;
  }
  AliError(Form("Calling FillHistogram with 2 parameters for histo <%s> of type %s",key,tmp->IsA()->GetName())) ;
}

//_____________________________________________________________________________
void AliAnalysisTaskEpRatio::SetGeometry()
{
  //Init geometry.
  
  if(!fPHOSGeo){
    
    AliOADBContainer geomContainer("phosGeo");
    geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
    TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");

    if (fRunNumber < 224994)
      fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP") ;
    else
      fPHOSGeo =  AliPHOSGeometry::GetInstance("Run2") ;
    
    for(Int_t mod=0; mod<5; mod++) {
      if(!matrixes->At(mod)) {
        if( fDebug )
          AliInfo(Form("No PHOS Matrix for mod:%d, geo=%p\n", mod, fPHOSGeo));
        continue;
      }
      else {
        fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
        if( fDebug >1 )
          AliInfo(Form("Adding PHOS Matrix for mod:%d, geo=%p\n", mod, fPHOSGeo));
      }
    }
  } 
  
}

//_____________________________________________________________________________
void AliAnalysisTaskEpRatio::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
  //Fills 1D histograms with key
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y,z) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH3F")){
    ((TH3F*)tmp)->Fill(x,y,z) ;
    return ;
  }
}
