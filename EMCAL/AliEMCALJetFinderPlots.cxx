
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


/*
 
$Log$
Revision 1.1.1.1  2003/05/29 18:56:58  horner
Initial import - Mark


*/

//_________________________________________________________________________
//  Class for Filling JetFinder Plots
//
//*-- Author: Mark Horner (LBL/UCT)
//
//


#include "TMath.h"
#include "AliEMCALJetFinderPlots.h"

ClassImp(AliEMCALJetFinderPlots)
	
AliEMCALJetFinderPlots::AliEMCALJetFinderPlots()
{
  fInitialised = kFALSE;
  fNominalEnergy = 0.0;
  fConeRadius = 0.5;
}

void AliEMCALJetFinderPlots::InitPlots()
{
	
    hFragmFcn = new TH1F("hFragmFcn","Fragmentation Function",100,0,1);
    hFragmFcn->Sumw2();
    hPartonFragmFcn = new TH1F("hPartonFragmFcn","Parton Fragmentation Function",100,0,1);
    hPartonFragmFcn->Sumw2();
    hPartonJT = new TH1F("hPartonJT","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
    hPartonJT->Sumw2();
    hPartonPL = new TH1F("hPartonPL","Track Momentum Parallel to Parton Axis ",100,0.,100.);
    hPartonPL->Sumw2();
    hJetJT = new TH1F("hJetJT","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
    hJetJT->Sumw2();
    hJetPL = new TH1F("hJetPL","Track Momentum Parallel to Jet Axis ",100,0.,100.);
    hJetPL->Sumw2();
    hJetEt = new TH1F("hJetEt","E_{T}^{reco}",250,0.,250.);
    hJetEt->Sumw2();
    hJetEta = new	TH1F("hJetEta","#eta_{jet}^{reco}",180,-0.9,0.9);
    hJetEta->Sumw2();
    hJetPhi = new 	TH1F("hJetPhi","#phi_{jet}^{reco}",62,0.,3.1);
    hJetPhi->Sumw2();
    hPartonEta = new	TH1F("hPartonEta","#eta_{Parton}",180,-0.9,0.9);
    hPartonEta->Sumw2();
    hPartonPhi = new 	TH1F("hPartonPhi","#phi_{Parton}",62,0.,3.1);
    hPartonPhi->Sumw2();
    hEtaDiff = new	TH1F("hEtaDiff","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
    hEtaDiff->Sumw2();
    hPhiDiff  = new TH1F("hPhiDiff","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
    hPhiDiff->Sumw2();
    hNJets = new	TH1F("hNJets","N Reconstructed jets",11,-0.5,10.5);
    hNJets->Sumw2();
    hEtaPhiSpread = new TH2F("hEtaPhiSpread","#eta - #phi Distribution of Reconstructed Jets",100,-0.5,0.5,100,-0.5,0.5);
    hEtaPhiSpread->Sumw2();    
  hNJets->SetXTitle("N_{jets}^{reco}/event");
  hNJets->SetYTitle("N_{events}");

  //Jet properties
  hJetEt->SetFillColor(16);
  hJetEt->SetXTitle("E_{T}^{reco}");
  
  hJetEta->SetFillColor(16);
  hJetEta->SetXTitle("#eta_{jet}^{reco}");
  
  hJetPhi->SetFillColor(16);
  hJetPhi->SetXTitle("#phi_{jet}^{reco}");
  
  hPartonEta->SetFillColor(16);
  hPartonEta->SetXTitle("#eta_{parton}");
  
  hPartonPhi->SetFillColor(16);
  hPartonPhi->SetXTitle("#phi_{parton}");

  hPartonPL->SetXTitle("p (GeV/c)");
  hPartonJT->SetXTitle("p (GeV/c)");

  hPartonFragmFcn->SetXTitle("Z = p_{T}^{Chg}/E_{T}^{parton}");

  //Jet component properties

  hJetPL->SetXTitle("p (GeV/c)");
  hJetJT->SetXTitle("p (GeV/c)");
  hFragmFcn->SetXTitle("Z = p_{T}^{Chg}/E_{T}^{reco}");
  hPartonFragmFcn->SetXTitle("Z = p_{T}^{Chg}/E_{T}^{reco}");

  hEtaDiff->SetXTitle("#eta_{jet}^{reco}-#eta_{jet}^{input}");
  hPhiDiff->SetXTitle("#phi_{jet}^{reco}-#phi_{jet}^{input}");
  hEtaPhiSpread->SetXTitle("#eta");
  hEtaPhiSpread->SetYTitle("#phi");
  fInitialised = kTRUE;	
  
}

AliEMCALJetFinderPlots::~AliEMCALJetFinderPlots(){

delete    hFragmFcn;// = new TH1F("hFragmFcn","Fragmentation Function",100,0,1);
delete   hPartonJT;// = new TH1F("hPartonJT","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
delete    hPartonPL;// = new TH1F("hPartonPL","Track Momentum Parallel to Parton Axis ",100,0.,100.);
delete    hJetJT;// = new TH1F("hJetJT","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
delete    hJetPL;// = new TH1F("hJetPL","Track Momentum Parallel to Jet Axis ",100,0.,100.);
delete    hJetEt;// = new TH1F("hJetEt","E_{T}^{reco}",250,0.,250.);
delete    hJetEta;// = new       TH1F("hJetEta","#eta_{jet}^{reco}",180,-0.9,0.9);
delete    hJetPhi;// = new       TH1F("hJetPhi","#phi_{jet}^{reco}",62,0.,3.1);
delete    hPartonEta;// = new    TH1F("hPartonEta","#eta_{Parton}",180,-0.9,0.9);
delete    hPartonPhi;// = new    TH1F("hPartonPhi","#phi_{Parton}",62,0.,3.1);
delete    hEtaDiff;// = new      TH1F("hEtaDiff","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
delete    hPhiDiff;//  = new TH1F("hPhiDiff","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
delete    hNJets;// = new        TH1F("hNJets","N Reconstructed jets",11,-0.5,10.5);
delete 	  hEtaPhiSpread;	    
}	

void AliEMCALJetFinderPlots::FillFromOutput(AliEMCALJetFinderOutput* output)
{
if (!fInitialised) InitPlots();	
  fOutput = output;
	
  Int_t nPartons = fOutput->GetNPartons();
  hNJets->Fill(fOutput->GetNJets());
  if (fOutput->GetNJets()!=1) return;
  AliEMCALParton* parton;
  AliEMCALJet* jet; 
  jet = fOutput->GetJet(0);
  
  hJetEt->Fill(jet->Energy());
  hJetEta->Fill(jet->Eta()  );
  hJetPhi->Fill(jet->Phi() );
  if (nPartons ==0) return;
  parton = fOutput->GetParton(0);
  hPartonEta->Fill( parton->Eta() );
  hPartonPhi->Fill( parton->Phi() );

  //hJetEtDiff->Fill( jet->Energy() - parton->Energy() );
  hEtaDiff->Fill( jet->Eta() - parton->Eta() );
  hPhiDiff->Fill( jet->Phi() - parton->Phi() );
  hEtaPhiSpread->Fill(jet->Eta()-parton->Eta(),jet->Phi() - parton->Phi());
 /* 
  Float_t *pt,*phi,*eta;
  Int_t *pdg;
  pt  = new Float_t[parton->GetNTracks()];
  eta = new Float_t[parton->GetNTracks()];
  phi = new Float_t[parton->GetNTracks()];
  pdg = new Int_t[parton->GetNTracks()];*/



 Float_t pt[2000];
 Float_t eta[2000];
 Float_t phi[2000];
 Int_t pdg[2000];
  
  parton->GetTrackList(pt,eta,phi,pdg);
  for(Int_t iT=0; iT< parton->GetNTracks() ; iT++ ) 
  {
      if ( (eta[iT]-parton->Eta())*(eta[iT]-parton->Eta())+
           (phi[iT]-parton->Phi())*(phi[iT]-parton->Phi()) >fConeRadius * fConeRadius ) continue; 
      Double_t tt = 2.0*atan(exp(-eta[iT])); // These names are short to make the equation manageable
      Double_t rt = 2.0*atan(exp(-parton->Eta()));
      Double_t ctt = cos(tt);
      Double_t crt = cos(rt);
      Double_t stt = sin(tt);
      Double_t srt = sin(rt);
      Double_t ctp = cos(phi[iT]);
      Double_t crp = cos(parton->Phi());
      Double_t stp = sin(phi[iT]);
      Double_t srp = sin(parton->Phi());
      Double_t alpha = acos(crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt);
      Double_t correctp = pt[iT]/stt; 	
      hPartonPL->Fill( correctp*cos(alpha));
      if ( (parton->Eta()-eta[iT])*(parton->Eta()-eta[iT]) +  
      (parton->Phi()-phi[iT])*(parton->Phi()-phi[iT]) < 0.2*0.2 ) 
    	      hPartonJT->Fill( correctp*sin(alpha));
      if (fNominalEnergy == 0.0) {
	      hPartonFragmFcn->Fill(  correctp*sin(tt)/parton->Energy()  );
      }else 
      {
              hPartonFragmFcn->Fill(correctp*sin(tt)/fNominalEnergy);
      }
  }// loop over tracks

/*
  pt  = new Float_t[jet->NTracks()];
  eta = new Float_t[jet->NTracks()];
  phi = new Float_t[jet->NTracks()];
  pdg = new Int_t[jet->NTracks()];*/
  jet->TrackList(pt,eta,phi,pdg);
  for(Int_t iT=0; iT< jet->NTracks() ; iT++ )
  {
      	  Double_t tt = 2.0*atan(exp(-eta[iT])); // These names are short to make the equation manageable
	  Double_t rt = 2.0*atan(exp(-jet->Eta()));
	  Double_t ctt = cos(tt);                   
	  Double_t crt = cos(rt);                         
	  Double_t stt = sin(tt);                               
	  Double_t srt = sin(rt);
	  Double_t ctp = cos(phi[iT]);
	  Double_t crp = cos(jet->Phi());
	  Double_t stp = sin(phi[iT]);
	  Double_t srp = sin(jet->Phi());
	  Double_t alpha = acos(crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt);   
	  Double_t correctp = pt[iT]/stt;
	  hJetPL->Fill( correctp*cos(alpha));      
	  if ( (jet->Eta()-eta[iT])*(jet->Eta()-eta[iT]) +  
			  (jet->Phi()-phi[iT])*(jet->Phi()-phi[iT]) < 0.2*0.2 )   
		  hJetJT->Fill( correctp*sin(alpha));
	  if (fNominalEnergy==0.0){
		  hFragmFcn->Fill(  correctp*sin(tt)/parton->Energy()  );
	  } else
	  {
                  hFragmFcn->Fill(  correctp*sin(tt)/fNominalEnergy );
	  }
  }// loop over tracks
  
}


	
