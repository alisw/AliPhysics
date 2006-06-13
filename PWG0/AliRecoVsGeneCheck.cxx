#include "AliRecoVsGeneCheck.h"

#include <Riostream.h>

//____________________________________________________________________
ClassImp(AliRecoVsGeneCheck);

//____________________________________________________________________
AliRecoVsGeneCheck::AliRecoVsGeneCheck() {

  fhVtzZRecoVsMC     = new TH2F("vtx_z_reco_vs_mc","",100,-20,20,100,-20,20);

  fhVtxZRes          = new TH1F("vtx_z_res","",100,0,0.05);
  fhVtxZResVsZ       = new TH2F("vtx_z_res_vs_z","",100,-20,20,100,0,0.05);
  fhVtxZResVsNPart   = new TH2F("vtx_z_res_vs_npart","",100,-0.5,99.5,100,0,0.05);

  fhVtxDzNorm        = new TH1F("vtx_dz_norm","",100,-10,10);
  fhVtxDzNormVsZ     = new TH2F("vtx_dz_norm_vs_z","",100,-20,20,100,-10,10);
  fhVtxDzNormVsNPart = new TH2F("vtx_dz_norm_vs_npart","",100,-0.5,99.5,100,-10,10);

  fhVtxZMC           = new TH1F("vtx_z_mc","",100,-20,20);
  fhVtxZReco         = new TH1F("vtx_z_reco","",100,-20,20);

  fhNPart            = new TH1F("n_part","",100,-0.5,99.5);

  fhDPtVsPtVsEta     = new TH3F("dpt_vs_pt_vs_eta","", 30,-1.5,1.5,30,0,3,120,-0.3,0.3);
  fhDEtaVsPtVsEta    = new TH3F("deta_vs_pt_vs_eta","",30,-1.5,1.5,30,0,3,120,-0.3,0.3);
  
  fhVtxZMC  ->Sumw2();
  fhVtxZReco->Sumw2();
  fhNPart   ->Sumw2();

  fhVtzZRecoVsMC    ->SetXTitle("vertex z_{mc}"); 
  fhVtzZRecoVsMC    ->SetYTitle("vertex z_{reco}"); 

  fhVtxZRes         ->SetXTitle("vertex #sigma_{z_{reco}}");
  fhVtxZResVsZ      ->SetXTitle("vertex z_{mc}");
  fhVtxZResVsZ      ->SetYTitle("vertex #sigma_{z_{reco}}");  		 
  fhVtxZResVsNPart  ->SetXTitle("n part");
  fhVtxZResVsNPart  ->SetYTitle("vertex #sigma_{z_{reco}}");
    
  fhVtxDzNorm       ->SetXTitle("vertex (z_{mc} - z_{reco}) / #sigma_{z_{reco}}");
  fhVtxDzNormVsZ    ->SetXTitle("vertex z_{mc}");
  fhVtxDzNormVsZ    ->SetYTitle("vertex (z_{mc} - z_{reco}) / #sigma_{z_{reco}}");
  fhVtxDzNormVsNPart->SetXTitle("npart");				     
  fhVtxDzNormVsNPart->SetYTitle("vertex (z_{mc} - z_{reco}) / #sigma_{z_{reco}}"); 

  fhNPart           ->SetXTitle("n part");  		
 
  fhVtxZMC          ->SetXTitle("vertex z_{mc}");
  fhVtxZReco        ->SetXTitle("vertex z_{reco}");

  fhDPtVsPtVsEta ->SetXTitle("#eta");
  fhDPtVsPtVsEta ->SetYTitle("p_{T} [GeV/c]");
  fhDPtVsPtVsEta ->SetZTitle("#deltap_{T} [GeV/c]");
  
  fhDEtaVsPtVsEta ->SetXTitle("#eta");
  fhDEtaVsPtVsEta ->SetYTitle("p_{T} [GeV/c]");
  fhDEtaVsPtVsEta ->SetZTitle("#delta#eta");
  
}

//____________________________________________________________________
void
AliRecoVsGeneCheck::Event(Double_t* vtx, Double_t* vtx_res, Double_t* mcvtx, Int_t n_part) {

  fhVtzZRecoVsMC     ->Fill(mcvtx[2],vtx[2]);     
  		     
  fhVtxZRes          ->Fill(vtx_res[2]);
  fhVtxZResVsZ       ->Fill(mcvtx[2], vtx_res[2]);
  fhVtxZResVsNPart   ->Fill(n_part, vtx_res[2]);
  		     
  if (vtx_res[2]!=0) {
    Float_t dzNorm = (mcvtx[2] - vtx[2])/vtx_res[2];

    fhVtxDzNorm       ->Fill(dzNorm);       
    fhVtxDzNormVsZ    ->Fill(mcvtx[2], dzNorm);         
    fhVtxDzNormVsNPart->Fill(n_part, dzNorm);         

  }  		     
  fhVtxZMC      ->Fill(mcvtx[2]);
  fhVtxZReco    ->Fill(vtx[2]);
  fhNPart       ->Fill(n_part);     
}


//____________________________________________________________________
void
AliRecoVsGeneCheck::Track(AliESDtrack* esdTrack, TParticle* mcParticle) {

  Double_t p[3];
  esdTrack->GetConstrainedPxPyPz(p); // ### TODO or GetInnerPxPyPy / GetOuterPxPyPy
  TVector3 vector(p);
  
  Float_t reco_p      = esdTrack->GetP();
  Float_t reco_theta  = vector.Theta();
  //  Float_t reco_phi    = 
  Float_t reco_eta    = -TMath::Log(TMath::Tan(reco_theta/2.));
  Float_t reco_pt     = reco_p*TMath::Sin(reco_theta);

  Float_t mc_eta = mcParticle->Eta();
  Float_t mc_pt = mcParticle->Pt();

  fhDPtVsPtVsEta->Fill(mc_eta, mc_pt, mc_pt - reco_pt);

  fhDEtaVsPtVsEta->Fill(mc_eta, mc_pt, mc_eta - reco_eta);

}

//____________________________________________________________________
void 
AliRecoVsGeneCheck::SaveHistograms(Char_t* dir) {

  gDirectory->mkdir(dir);
  gDirectory->cd(dir);

  gDirectory->mkdir("event");
  gDirectory->cd("event");

  fhVtzZRecoVsMC      ->Write();
  fhVtxZRes           ->Write();
  fhVtxZResVsZ        ->Write();
  fhVtxZResVsNPart    ->Write();  			  
  fhVtxDzNorm         ->Write();
  fhVtxDzNormVsZ      ->Write();
  fhVtxDzNormVsNPart  ->Write();
  fhVtxZMC            ->Write();
  fhVtxZReco          ->Write();
  fhNPart             ->Write();

  gDirectory->cd("../");


  gDirectory->mkdir("track");
  gDirectory->cd("track");
  
  fhDPtVsPtVsEta      ->Write();
  fhDEtaVsPtVsEta     ->Write();
  
  gDirectory->cd("../");

  gDirectory->cd("../");
}
