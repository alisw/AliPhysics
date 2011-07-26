
#if !defined(__CINT__) || defined(__MAKECINT__)

//Root include files 
//#include <Riostream.h>
#include <TFile.h>
//#include <TSystem.h>
#include <TH1F.h>
#include <TParticle.h>
#include <TRefArray.h>

//AliRoot include files 
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloTrigger.h"
#include "AliPID.h"
#include "AliEMCALGeometry.h"

#endif

//Change the bool depending on what information you want to print
// when all FALSE, prints minimum cluster information.
Bool_t kPrintKine         = kFALSE; //Do not use for raw data.
Bool_t kPrintCaloCells    = kFALSE;
Bool_t kPrintCaloTrigger  = kFALSE;
Bool_t kPrintTrackMatches = kFALSE;
Bool_t kPrintClusterCells = kFALSE;
Bool_t kPrintClusterPID   = kFALSE;

void TestESD() {
  // Main method to read information stored in AliESDCaloClusters and AliESDCaloCells
	
  // Init some example histograms
  // ESD
  TH1F * hEta  = (TH1F*) new TH1F("hEta","reco eta",1000, -0.71,0.71); 
  TH1F * hPhi  = (TH1F*) new TH1F("hPhi","reco phi",130, 70,200); 
  TH1F * hE    = (TH1F*) new TH1F("hE"  ,"reco e",300, 0,30); 
  TH1F * hTime = (TH1F*) new TH1F("hTime"  ,"reco time",1000, 0,1000); 
  hEta ->SetXTitle("#eta");
  hPhi ->SetXTitle("#phi (deg)");
  hE   ->SetXTitle("E (GeV)");
  hTime->SetXTitle("time (ns)");
  
  // Monte Carlo
  TH1F * hMCEta = (TH1F*) new TH1F("hMCEta","MC eta",1000, -0.71,0.71); 
  TH1F * hMCPhi = (TH1F*) new TH1F("hMCPhi","MC phi",130, 70,200); 
  TH1F * hMCE   = (TH1F*) new TH1F("hMCE"  ,"MC e",300, 0,30); 
  hMCEta->SetXTitle("#eta");
  hMCPhi->SetXTitle("#phi (deg)");
  hMCE  ->SetXTitle("E (GeV)");
  
  //ESD - MonteCarlo
  TH1F * hDEta = (TH1F*) new TH1F("hDEta"," eta cluster - eta MC",500, -0.05,0.05); 
  TH1F * hDPhi = (TH1F*) new TH1F("hDPhi","phi cluster - phi MC",200, -20,20); 
  TH1F * hDE   = (TH1F*) new TH1F("hDE"  ,"e cluster - e MC",200, -10,10); 
  hDEta->SetXTitle("#eta_{reco}-#eta_{MC}");
  hDPhi->SetXTitle("#phi_{reco}-#phi_{MC} (deg)");
  hDE  ->SetXTitle("E_{reco}-E_{MC} (GeV)");
  
  //Open the ESD file, get the tree with events
  TFile* f = new TFile("AliESDs.root");
  TTree* esdTree = (TTree*)f->Get("esdTree");
  AliESDEvent* esd = new AliESDEvent();
  esd->ReadFromTree(esdTree);

  //Init geometry and array that will contain the clusters
  TRefArray* caloClusters = new TRefArray();
  AliEMCALGeometry *geom =  AliEMCALGeometry::GetInstance("EMCAL_FIRSTYEAR") ;  	
  Float_t pos[3];
  
  //Loop of events
  Int_t nEvt = esdTree->GetEntries();
  for(Int_t iev = 0; iev < nEvt; iev++) {
    cout << "<<<< Event: " << iev+1 << "/" << nEvt << " >>>>"<<endl;
    esdTree->GetEvent(iev);
    
    //Pass the geometry transformation matrix from ESDs to geometry
    for(Int_t mod=0; mod < (geom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
      //printf("matrix %d\n",mod);
      if(esd->GetEMCALMatrix(mod)) {
        //printf("EMCAL: mod %d, matrix %p\n",mod, esd->GetEMCALMatrix(mod));
        //(esd->GetEMCALMatrix(mod))->Print("");
        geom->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
      }//matrix
    }//module
    
    //In case you want to play with MC data, get stack
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
    
    //------------------------------------------------------
    //Get Cells Array, all cells in event, print cells info 
    //------------------------------------------------------ 
    AliVCaloCells &cells= *(esd->GetEMCALCells());
    //AliESDCaloCells &cells= *(esd->GetEMCALCells());

    if(kPrintCaloCells){  
      Int_t nTotalCells = cells.GetNumberOfCells() ;  
      //Int_t type        = cells.GetType();
      for (Int_t icell=  0; icell <  nTotalCells; icell++) {
	cout<<"Cell   : "<<icell<<"/"<<nTotalCells<<" ID: "<<cells.GetCellNumber(icell)<<" Amplitude: "<<cells.GetAmplitude(icell)<<" Time: "<<cells.GetTime(icell)*1e9<<endl;	  
      }// cell loop
    }
    
    //------------------------------------------------------
    // Calo Trigger 
    //------------------------------------------------------
	
	  if(kPrintCaloTrigger)
	  {  
		  AliESDCaloTrigger& trg = *(esd->GetCaloTrigger("EMCAL"));
		  
		  trg.Reset();
		  while (trg.Next())
		  {
			  int posX, posY;
			  trg.GetPosition(posX, posY);
			  
			  if (posX > -1 && posY > -1) 
			  {
				  Int_t ts = 0;
				  trg.GetL1TimeSum(ts);
				  
				  cout << "Position: " << posX << " " << posY << " L1 amplitude: " << ts << endl;
			  }
		  }
	  }

    //------------------------------------------------------
    // Calo Clusters 
    //------------------------------------------------------
    
    //Get CaloClusters Array
    caloClusters->Clear();
    esd->GetEMCALClusters(caloClusters);
    
    //loop over clusters
    Int_t nclus = caloClusters->GetEntries();
    for (Int_t icl = 0; icl < nclus; icl++) {
      
      AliVCluster* clus = (AliVCluster*)caloClusters->At(icl);
      //AliESDCluster* clus = (AliESDCluster*)caloClusters->At(icl);
      Float_t energy = clus->E();
      clus->GetPosition(pos);
      TVector3 vpos(pos[0],pos[1],pos[2]);
      
      //We can get a momentum TLorentzVector per cluster, corrected by the vertex position 
      //TLorentzVector p;
      //clus->GetMomentum(p,vertex_position);
      
      Double_t cphi = vpos.Phi();
      Double_t ceta = vpos.Eta();
      
      Int_t nMatched   = clus->GetNTracksMatched();
      Int_t trackIndex = clus->GetTrackMatchedIndex();
      Int_t nLabels    = clus->GetNLabels();
      Int_t labelIndex = clus->GetLabel();
      Int_t nCells     = clus->GetNCells();
		
      //Fill some histograms
      hEta->Fill(ceta);
      hPhi->Fill(cphi*TMath::RadToDeg());
      hE  ->Fill(energy);
      hTime->Fill(clus->GetTOF()*1e9);
      
      //Print basic cluster information
      cout << "Cluster: " << icl+1 << "/" << nclus << " Energy: " << energy << "; Phi: " 
	   << cphi*TMath::RadToDeg() << "; Eta: " << ceta << "; NCells: " << nCells << " ;#Matches: " << nMatched 
	   << "; Index: " << trackIndex << "; #Labels: " << nLabels << " Index: " 
	   << labelIndex << "; Time "<<clus->GetTOF()*1e9<<" ns "<<endl;
      
      //Print primary info
      if(stack && kPrintKine) {
	if(labelIndex >= 0 && labelIndex < stack->GetNtrack()){
	  TParticle * particle = stack->Particle(labelIndex);
	  //Fill histograms with primary info
	  hMCEta->Fill(particle->Eta());
	  hMCPhi->Fill(particle->Phi()*TMath::RadToDeg());
	  hMCE  ->Fill(particle->Energy());
	  hDEta ->Fill(ceta-particle->Eta());
	  hDPhi ->Fill(cphi*TMath::RadToDeg()-particle->Phi()*TMath::RadToDeg());
	  hDE   ->Fill(energy-particle->Energy());
	  //Print primary values
	  cout<<"         More  contributing primary: "<<particle->GetName()<<"; with kinematics: "<<endl;
	  cout<<" \t     Energy: "<<particle->Energy()<<"; Phi: "<<particle->Phi()*TMath::RadToDeg()<<"; Eta: "<<particle->Eta()<<endl;   
	  for(Int_t i = 1; i < nLabels; i++){
	    //particle = stack->Particle((((AliESDCaloCluster*)clus)->GetLabelsArray())->At(i));
	    particle = stack->Particle((clus->GetLabels())[i]);
	    //or Int_t *labels = clus->GetLabels();
	    //particle = stack->Particle(labels[i]);
	    cout<<"         Other contributing primary: "<<particle->GetName()<< "; Energy "<<particle->Energy()<<endl;
	  }
	}
	else if( labelIndex >= stack->GetNtrack()) cout <<"PROBLEM, label is too large : "<<labelIndex<<" >= particles in stack "<< stack->GetNtrack() <<endl;
		  else cout<<"Negative label!!!  : "<<labelIndex<<endl;
      } // play with stack
      
      // Matching results
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
      }// matching
      
      //Get PID weights and print them
      if(kPrintClusterPID){
	const Double_t *pid = clus->GetPID();
	printf("PID weights: ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f, hadrons: pion %0.2f, kaon %0.2f, proton %0.2f , neutron %0.2f, kaon %0.2f \n",
	       pid[AliVCluster::kPhoton],   pid[AliVCluster::kPi0],
	       pid[AliVCluster::kElectron], pid[AliVCluster::kEleCon],
	       pid[AliVCluster::kPion],     pid[AliVCluster::kKaon],   pid[AliVCluster::kProton],
	       pid[AliVCluster::kNeutron],  pid[AliVCluster::kKaon0]);
      }//PID
      
      //Get CaloCells of cluster and print their info, position.
      if(kPrintClusterCells){	
	UShort_t * index    = clus->GetCellsAbsId() ;
	Double_t * fraction = clus->GetCellsAmplitudeFraction() ;
	Int_t sm = -1;
	for(Int_t i = 0; i < nCells ; i++){
	  Int_t absId       =   index[i]; // or clus->GetCellNumber(i) ;
	  Double_t ampFract =  fraction[i];
	  Float_t amp       = cells.GetCellAmplitude(absId) ;
	  Double_t time     = cells.GetCellTime(absId);
	  cout<<"\t Cluster Cell: AbsID : "<< absId << " == "<<clus->GetCellAbsId(i) <<"; Amplitude "<< amp << "; Fraction "<<ampFract<<"; Time " <<time*1e9<<endl;
	  //Geometry methods  
	  Int_t iSupMod =  0 ;
	  Int_t iTower  =  0 ;
	  Int_t iIphi   =  0 ;
	  Int_t iIeta   =  0 ;
	  Int_t iphi    =  0 ;
	  Int_t ieta    =  0 ;
	  if(geom){
	    geom->GetCellIndex(absId,iSupMod,iTower,iIphi,iIeta); 
	    //Gives SuperModule and Tower numbers
	    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
					      iIphi, iIeta,iphi,ieta);
	    //Gives label of cell in eta-phi position per each supermodule
	    Float_t cellPhi = 0;
	    Float_t cellEta = 0;
	    geom->EtaPhiFromIndex(absId,cellEta,cellPhi);
	    cout<< "                SModule "<<iSupMod<<"; Tower "<<iTower <<"; Eta "<<iIeta
		<<"; Phi "<<iIphi<<"; Index: Cell Eta "<<ieta<<"; Cell Phi "<<iphi
		<<"; Global: Cell Eta "<<cellEta<<"; Cell Phi "<<cellPhi*TMath::RadToDeg()<<endl;
	    if(i==0) sm = iSupMod;
	    else{
	      if(sm!=iSupMod) printf("******CLUSTER SHARED BY 2 SuperModules!!!!\n");
	    }	
	  }// geometry on
	}// cluster cell loop
      }// print cell clusters
    } //cluster loop
  } // event loop
  
  
  //Write histograms in a file
  TFile * fhisto = (TFile*) new TFile("histos.root","recreate");
  hEta->Write();
  hPhi->Write();
  hE  ->Write();
  hTime->Write();
  if(kPrintKine){
    hMCEta->Write();
    hMCPhi->Write();
    hMCE  ->Write();
    hDEta->Write();
    hDPhi->Write();
    hDE  ->Write();
  }
  fhisto->Close();
}
