void ReadAOD_MCBranch(char *filein="AliAODs.root"){

//=============================================================================
// Macro to read AOD+MC Branch: 
// in this example the macro reads generated and reconstructed JPsi 
// Plots refer to all MC generated particles, to MC particles which have
// been reconstructed and reconstructed tracks
//
// Each aod event has 3 MC particles (J/Psi, mu1, mu2)
// To access MC muons:   mctrackall->IsPhysicalPrimary()
// To access JPsi muons: mctrackall->IsPrimary() && !mctrackall->IsPhysicalPrimary()

// AliStack::IsPhysicalPrimary() 
// Test if a particle is a physical primary according to the following definition:
// Particles produced in the collision including products of strong and
// electromagnetic decay and excluding feed-down from weak decays of strange
// particles.
//=============================================================================

  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPalette(1); 

// Open AOD file    
  TFile *file = new TFile(filein);
  TTree *aodTree = file->Get("aodTree");
  
// Read AOD event    
  AliAODEvent *event = new AliAODEvent();
  event->ReadFromTree(aodTree);
  
// Read MC info   
  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
  if(!mcarray) {
    printf("No MC info in AOD!\n");
    return;
  }  
  
// Book Histos
  TH2D *hpt_mc_rec = new TH2D("hpt_mc_rec","Pt correlations MC-REC",50,0.,10.,50,0.,10.);
  TH1D *hpt_rec = new TH1D("hpt_rec","Pt REC tracks",20,0.,10.);
  TH1D *hpt_mc = new TH1D("hpt_mc","PT MC particles corresponding to rec tracks",20,0.,10.);
  TH1D *hpt_allmc = new TH1D("hpt_allmc","PT all MC particles",20,0.,10.);
  TH1D *hpt_primarymc = new TH1D("hpt_primarymc","PT MC primaries (J/Psi and muons)",20,0.,10.);
  TH1D *hpt_physprimarymc = new TH1D("hpt_physprimarymc","PT Muons MC",20,0.,10.);
  TH1D *hpt_notphysprimarymc = new TH1D("hpt_notphysprimarymc","Pt JPsi MC",20,0.,10.);
  TH1D *hpt_diff = new TH1D("hpt_diff","Pt: Rec - MC",100,-5.,5.);

  TH1D *hCharge_primarymc = new TH1D("hCharge_primarymc","Charge MC primary particle",10,-5.,5.);
  TH1D *hCharge_physprimmc = new TH1D("hCharge_physprimmc","Charge Muons MC",10,-5.,5.);
  TH1D *hCharge_notphysprimmc = new TH1D("hCharge_notphysprimmc","Charge JPsi MC",10,-5.,5.);
  TH1D *hCharge_rec = new TH1D("hCharge_rec","hCharge_rec",10,-5.,5.);
  TH1D *hCharge_mc = new TH1D("hCharge_mc","hCharge_mc",10,-5.,5.);

  TH1D *hMass_mc = new TH1D("hMass_mc","Mass MC particles",100,0.,5..);
  TH1D *hMass_rec = new TH1D("hMass_rec","Mass rec tracks",100,0.,5..);
  TH1D *hMass_allmc = new TH1D("hMass_allmc","Mass all MC particles",100,0.,5..);
  TH1D *hMass_primarymc = new TH1D("hMass_primarymc","Mass primary particles MC",100,0.,5..);
  TH1D *hMass_physprimarymc = new TH1D("hMass_physprimarymc","Mass muons MC",100,0.,5..);
  TH1D *hMass_notphysprimarymc = new TH1D("hMass_notphysprimarymc","Mass JPsi MC",100,0.,5..);
  TH1D *hMassJPsi_mc = new TH1D("hMassJPsi_mc","Mass JPsi MC",100,0.,5..);
  
// Loop on AOD events
  Int_t nev=0;

  while(aodTree->GetEvent(nev++)){  

// Loop on all MC tracks    
    for(int ii=0;ii<mcarray->GetEntries();ii++){
      AliAODMCParticle *mctrackall = (AliAODMCParticle*) mcarray->At(ii);     
      // all particles
      hpt_allmc->Fill(mctrackall->Pt());
      hMass_allmc->Fill(mctrackall->M());
      
      // Muons
      if(mctrackall->IsPhysicalPrimary()) {
        hCharge_physprimmc->Fill(mctrackall->Charge());
        hpt_physprimarymc->Fill(mctrackall->Pt());
        hMass_physprimarymc->Fill(mctrackall->M());
      }
      	
      if(mctrackall->IsPrimary()) {
        hpt_primarymc->Fill(mctrackall->Pt());
        hMass_primarymc->Fill(mctrackall->M());
	hCharge_primarymc->Fill(mctrackall->Charge());
      }
        
      // JPsi
      if(mctrackall->IsPrimary() && !mctrackall->IsPhysicalPrimary()) {
        hCharge_notphysprimmc->Fill(mctrackall->Charge());
        hMass_notphysprimarymc->Fill(mctrackall->M());
        hpt_notphysprimarymc->Fill(mctrackall->Pt());
      }	
    }
    
// Loop on reconstructed tracks    
    for(int i=0;i<event->GetNumberOfTracks();i++){
      AliAODTrack *track = dynamic_cast<AliAODTrack*>(event->GetTrack(i));
      if(!track) {
        printf("Error: no track! \n");
        continue;
      } 
      hpt_rec->Fill(track->Pt());
      hMass_rec->Fill(track->M());
      hCharge_rec->Fill(track->Charge());
      
// Read MC particle corresponding to reconstructed tracks   
      AliAODMCParticle *mctrack = (AliAODMCParticle*) mcarray->At(track->GetLabel());     
      if(!mctrack){
        printf("Error: no MC track! \n");
	continue;
      }
      //printf("nev=%d i=%d mc=%f rec=%f\n",nev-1,i,mctrack->Pt(),track->Pt());
      hpt_mc->Fill(mctrack->Pt());
      hpt_mc_rec->Fill(track->Pt(),mctrack->Pt());
      hpt_diff->Fill(track->Pt()-mctrack->Pt());
      hMass_mc->Fill(mctrack->M());
      hCharge_mc->Fill(mctrack->Charge());
   }
 }
 
// Plots 
 TCanvas *c = new TCanvas("c","PT spectra",20,20,600,600);
 c->Divide(3,3);
 c->cd(1);
 hpt_primarymc->Draw();
 hpt_primarymc->GetXaxis()->SetTitle("Pt MC primary particles");
 c->cd(2);  
 hpt_mc->GetXaxis()->SetTitle("Pt MC particles");
 hpt_mc->Draw();
 hpt_mc->SetLineColor(4);
 c->cd(3);
 hpt_rec->GetXaxis()->SetTitle("Rec Pt track");
 hpt_rec->Draw();
 hpt_rec->SetLineColor(2);
 c->cd(4);
 hpt_physprimarymc->Draw();
 hpt_physprimarymc->GetXaxis()->SetTitle("Pt MC muon track");
 c->cd(5);
 hpt_notphysprimarymc->Draw();
 hpt_notphysprimarymc->GetXaxis()->SetTitle("Pt MC JPsi track");
 c->cd(7);
 hpt_mc_rec->Draw("colz");
 hpt_mc_rec->GetXaxis()->SetTitle("Pt rec track");
 hpt_mc_rec->GetYaxis()->SetTitle("Pt MC track");
 c->cd(8);
 hpt_diff->Draw();
 hpt_diff->GetXaxis()->SetTitle("Rec-MC Pt");

 TCanvas *c1 = new TCanvas("c1","Mass distributions",40,40,600,600);
 c1->Divide(2,2);
 c1->cd(1);
 hMass_primarymc->Draw();
 hMass_primarymc->GetXaxis()->SetTitle("Mass MC primary particles");
 c1->cd(2);
 hMass_physprimarymc->Draw();
 hMass_physprimarymc->GetXaxis()->SetTitle("Mass MC muon particles");
 c1->cd(3);
 hMass_notphysprimarymc->Draw();
 hMass_notphysprimarymc->GetXaxis()->SetTitle("Mass MC J/Psi");
 c1->cd(4);
 hMass_rec->Draw();
 hMass_rec->GetXaxis()->SetTitle("Mass Rec tracks");

 TCanvas *c2= new TCanvas("c2","Charge distributions",60,60,600,600);
 c2->Divide(2,2);
 c2->cd(1);
 hCharge_primarymc->Draw();
 hCharge_primarymc->GetXaxis()->SetTitle("Charge MC primary particles");
 c2->cd(2);
 hCharge_physprimmc->Draw();
 hCharge_physprimmc->GetXaxis()->SetTitle("Charge muon MC");;
 c2->cd(3);
 hCharge_notphysprimmc->Draw();
 hCharge_notphysprimmc->GetXaxis()->SetTitle("Charge J/Psi MC");;
 c2->cd(4);
 hCharge_rec->Draw();
 hCharge_rec->GetXaxis()->SetTitle("Charge rec tracks");
}
