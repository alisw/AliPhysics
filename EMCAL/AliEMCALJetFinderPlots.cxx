
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


/* $Id$ */


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
  fConeRadius = 0.3;
  fDebug = 0;
  fOutput=0;
    hFragmFcn=0;// = new TH1F("hFragmFcn","Fragmentation Function",100,0,1);
    hPartonFragmFcn=0;// = new TH1F("hPartonFragmFcn","Fragmentation Function",100,0,1);
   hPartonJT=0;// = new TH1F("hPartonJT","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
    hPartonPL=0;// = new TH1F("hPartonPL","Track Momentum Parallel to Parton Axis ",100,0.,100.);
    hJetJT=0;// = new TH1F("hJetJT","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
    hJetPL=0;// = new TH1F("hJetPL","Track Momentum Parallel to Jet Axis ",100,0.,100.);
    hJetEt=0;// = new TH1F("hJetEt","E_{T}^{reco}",250,0.,250.);
    hJetEta=0;// = new       TH1F("hJetEta","#eta_{jet}^{reco}",180,-0.9,0.9);
    hJetPhi=0;// = new       TH1F("hJetPhi","#phi_{jet}^{reco}",62,0.,3.1);
    hPartonEta=0;// = new    TH1F("hPartonEta","#eta_{Parton}",180,-0.9,0.9);
    hPartonPhi=0;// = new    TH1F("hPartonPhi","#phi_{Parton}",62,0.,3.1);
    hEtaDiff=0;// = new      TH1F("hEtaDiff","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
    hPhiDiff=0;//  = new TH1F("hPhiDiff","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
    hNJets=0;// = new        TH1F("hNJets","N Reconstructed jets",11,-0.5,10.5);
  hEtaPhiSpread=0;	    

hFragmFcn2=0;	// ("hFragmFcn2","Fragmentation Function",100,0,1);
hPartonFragmFcn2=0;// ("hFragmFcn2","Parton Fragmentation Function",100,0,1);
hPartonJT2=0;	// ("hPartonJT2","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
hPartonPL2=0;	// ("hPartonPL2","Track Momentum Parallel to Parton Axis ",100,0.,100.);
hJetJT2=0;	// ("hJetJT2","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
hJetPL2=0;	// ("hJetPL2","Track Momentum Parallel to Jet Axis ",100,0.,100.);
hJetEt2=0;	// ("hJetEt2","E_{T}^{reco}",250,0.,250.);
hJetEta2=0;	// ("hJetEta2","#eta_{jet}^{reco}",180,-0.9,0.9);
hJetPhi2=0;	// ("hJetPhi2","#phi_{jet}^{reco}",62,0.,3.1);
hPartonEta2=0;	// ("hPartonEta2","#eta_{Parton}",180,-0.9,0.9);
hPartonPhi2=0;	// ("hPartonPhi2","#phi_{Parton}",62,0.,3.1);
hEtaDiff2=0;	// ("hEtaDiff2","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
hPhiDiff2=0;	// ("hPhiDiff2","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
hEtaPhiSpread2=0;	// ("hEtaPhiSpread2","#eta - #phi Distribution 
						//of Reconstructed Jets",192,-0.7,0.7,288,pi/3,pi);
hNJets2=0;	// ("hNJets2","N Reconstructed jets",11,-0.5,10.5);
hJetEtSecond2=0; //("hJetEtSecond2","E_{T}^{reco}",250,0.,250.); 
hJetEtRatio2=0;  //("hJetEtRatio2","Ratio of Second Highest to Highest",100,0,1);
hEtaPhiDist2=0;  //("hEtaPhiDist2","Angular Distance Between First and Second",100,0,3);

}

void AliEMCALJetFinderPlots::InitPlots()
{
//========================= CASE 1 =======================================	
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

//======================= CASE 2 ======================================
  

hFragmFcn2 		= new TH1F("hFragmFcn2","Fragmentation Function",100,0,1);
hFragmFcn2->Sumw2();
hPartonFragmFcn2 	= new TH1F("hPartonFragmFcn2","Parton Fragmentation Function",100,0,1);
hPartonFragmFcn2->Sumw2();
hPartonJT2 		= new TH1F("hPartonJT2","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
hPartonJT2->Sumw2();
hPartonPL2		= new TH1F("hPartonPL2","Track Momentum Parallel to Parton Axis ",100,0.,100.);
hPartonPL2->Sumw2();
hJetJT2			= new TH1F("hJetJT2","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
hJetJT2->Sumw2();
hJetPL2			= new TH1F("hJetPL2","Track Momentum Parallel to Jet Axis ",100,0.,100.);
hJetPL2->Sumw2();
hJetEt2			= new TH1F("hJetEt2","E_{T}^{reco}",250,0.,250.);
hJetEt2->Sumw2();
hJetEta2		= new TH1F("hJetEta2","#eta_{jet}^{reco}",180,-0.9,0.9);
hJetEta2->Sumw2();
hJetPhi2		= new TH1F("hJetPhi2","#phi_{jet}^{reco}",62,0.,3.1);
hJetPhi2->Sumw2();
hPartonEta2		= new TH1F("hPartonEta2","#eta_{Parton}",180,-0.9,0.9);
hPartonEta2->Sumw2();
hPartonPhi2		= new TH1F("hPartonPhi2","#phi_{Parton}",62,0.,3.1);
hPartonPhi2->Sumw2();
hEtaDiff2		= new TH1F("hEtaDiff2","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
hEtaDiff2->Sumw2();
hPhiDiff2		= new TH1F("hPhiDiff2","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
hPhiDiff2->Sumw2();
hEtaPhiSpread2 = new TH2F("hEtaPhiSpread2","#eta - #phi Distribution of Reconstructed Jets",100,-0.5,0.5,100,-0.5,0.5);
hEtaPhiSpread2->Sumw2();
hNJets2			= new TH1F("hNJets2","N Reconstructed jets",11,-0.5,10.5);
hNJets2->Sumw2();
hJetEtSecond2		= new TH1F("hJetEtSecond2","E_{T}^{reco}",250,0.,250.); 
hJetEtSecond2->Sumw2();
hJetEtRatio2		= new TH1F("hJetEtRatio2","Ratio of Second Highest to Highest",100,0,1);
hJetEtRatio2->Sumw2();
hEtaPhiDist2		= new TH1F("hEtaPhiDist2","Angular Distance Between First and Second",100,0,3);
hEtaPhiDist2->Sumw2();
  
  fInitialised = kTRUE;	
  
}

AliEMCALJetFinderPlots::~AliEMCALJetFinderPlots(){

delete    hFragmFcn;// = new TH1F("hFragmFcn","Fragmentation Function",100,0,1);
delete    hPartonFragmFcn;// = new TH1F("hFragmFcn","Fragmentation Function",100,0,1);
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

	delete				hFragmFcn2;	// ("hFragmFcn2","Fragmentation Function",100,0,1);
	delete				hPartonFragmFcn2;// ("hFragmFcn2","Parton Fragmentation Function",100,0,1);
	delete				hPartonJT2;	// ("hPartonJT2","Track Momentum Perpendicular to Parton Axis",100,0.,10.);
	delete				hPartonPL2;	// ("hPartonPL2","Track Momentum Parallel to Parton Axis ",100,0.,100.);
	delete				hJetJT2;	// ("hJetJT2","Track Momentum Perpendicular to Jet Axis",100,0.,10.);
	delete				hJetPL2;	// ("hJetPL2","Track Momentum Parallel to Jet Axis ",100,0.,100.);
	delete				hJetEt2;	// ("hJetEt2","E_{T}^{reco}",250,0.,250.);
	delete				hJetEta2;	// ("hJetEta2","#eta_{jet}^{reco}",180,-0.9,0.9);
	delete				hJetPhi2;	// ("hJetPhi2","#phi_{jet}^{reco}",62,0.,3.1);
	delete				hPartonEta2;	// ("hPartonEta2","#eta_{Parton}",180,-0.9,0.9);
	delete				hPartonPhi2;	// ("hPartonPhi2","#phi_{Parton}",62,0.,3.1);
	delete				hEtaDiff2;	// ("hEtaDiff2","#eta_{jet}^{reco}-#eta_{jet}^{input}",100,-0.5,0.5);
	delete				hPhiDiff2;	// ("hPhiDiff2","#phi_{jet}^{reco}-#phi_{jet}^{input}",100,-0.5,0.5);
	delete				hEtaPhiSpread2;	// ("hEtaPhiSpread2","#eta - #phi Distribution 
							//of Reconstructed Jets",192,-0.7,0.7,288,pi/3,pi);
	delete				hNJets2;	// ("hNJets2","N Reconstructed jets",11,-0.5,10.5);
	delete				hJetEtSecond2; //("hJetEtSecond2","E_{T}^{reco}",250,0.,250.); 
	delete				hJetEtRatio2;  //("hJetEtRatio2","Ratio of Second Highest to Highest",100,0,1);
	delete 				hEtaPhiDist2;  //("hEtaPhiDist2","Angular Distance Between First and Second",100,0,3);



}	

void AliEMCALJetFinderPlots::FillFromOutput(AliEMCALJetFinderOutput* output)
{
if (!fInitialised) InitPlots();	
  fOutput = output;
  if (!fOutput) return;
  hNJets->Fill(fOutput->GetNJets());
if (fOutput->GetNJets()>1)
{	
//========================= CASE 2 ===========================	
  Int_t nPartons = fOutput->GetNPartons();
  hNJets2->Fill(fOutput->GetNJets());
  AliEMCALParton* parton;
  AliEMCALJet* jethighest=0; 
  AliEMCALJet* jetsecond=0; 
  //  Find Highest and Second Highest Jet
  for (Int_t counter = 0; counter<fOutput->GetNJets();counter++)
  {
	  if (counter==0){ 
		  jethighest = fOutput->GetJet(0);
		  jetsecond  = fOutput->GetJet(1);
	  }
	  if (counter>0)
	  {
		  Float_t energyhighest = jethighest->Energy();
		  Float_t energysecond = jetsecond->Energy();
		  
		  if ((fOutput->GetJet(counter))->Energy()>energyhighest) 
		  {
			  jetsecond=jethighest;
			  jethighest=fOutput->GetJet(counter);
		  }else if ((fOutput->GetJet(counter))->Energy()>energysecond) 
		  {
			  jetsecond=fOutput->GetJet(counter);
		  }
	  }
  }

  // End finding highest and second highest and continue
  hJetEt2->Fill(jethighest->Energy());
  hJetEta2->Fill(jethighest->Eta()  );
  hJetPhi2->Fill(jethighest->Phi() );
  if (nPartons ==0) return;
  parton = fOutput->GetParton(0);
  
  hPartonEta2->Fill( parton->Eta() );
  hPartonPhi2->Fill( parton->Phi() );

  //hJetEtDiff->Fill( jet->Energy() - parton->Energy() );
  hEtaDiff2->Fill( jethighest->Eta() - parton->Eta() );
  hPhiDiff2->Fill( jethighest->Phi() - parton->Phi() );
  hEtaPhiSpread2->Fill(jethighest->Eta()-parton->Eta(),jethighest->Phi() - parton->Phi());
  hJetEtSecond2->Fill(jetsecond->Energy()); 
  hJetEtRatio2->Fill(jetsecond->Energy()/jethighest->Energy()); 
  hEtaPhiDist2->Fill( TMath::Sqrt((jethighest->Eta() - jetsecond->Eta())*(jethighest->Eta() - jetsecond->Eta())
		      + (jethighest->Phi() - jetsecond->Phi())*(jethighest->Phi() - jetsecond->Phi())	  ));  
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
      hPartonPL2->Fill( correctp*cos(alpha));
      if ( (parton->Eta()-eta[iT])*(parton->Eta()-eta[iT]) +  
      (parton->Phi()-phi[iT])*(parton->Phi()-phi[iT]) < 0.2*0.2 ) 
    	      hPartonJT2->Fill( correctp*sin(alpha));
      if (fNominalEnergy == 0.0) {
	      hPartonFragmFcn2->Fill(  correctp*sin(tt)/parton->Energy()  );
      }else 
      {
              hPartonFragmFcn2->Fill(correctp*sin(tt)/fNominalEnergy);
      }
  }// loop over tracks

/*
  pt  = new Float_t[jet->NTracks()];
  eta = new Float_t[jet->NTracks()];
  phi = new Float_t[jet->NTracks()];
  pdg = new Int_t[jet->NTracks()];*/
  jethighest->TrackList(pt,eta,phi,pdg);
  for(Int_t iT=0; iT< jethighest->NTracks() ; iT++ )
  {
      	  Double_t tt = 2.0*atan(exp(-eta[iT])); // These names are short to make the equation manageable
	  Double_t rt = 2.0*atan(exp(-jethighest->Eta()));
	  Double_t ctt = cos(tt);                   
	  Double_t crt = cos(rt);                         
	  Double_t stt = sin(tt);                               
	  Double_t srt = sin(rt);
	  Double_t ctp = cos(phi[iT]);
	  Double_t crp = cos(jethighest->Phi());
	  Double_t stp = sin(phi[iT]);
	  Double_t srp = sin(jethighest->Phi());
	  Double_t alpha = acos(crp*ctp*srt*stt+srp*stp*srt*stt+crt*ctt);   
	  Double_t correctp = pt[iT]/stt;
	  hJetPL2->Fill( correctp*cos(alpha));      
	  if ( (jethighest->Eta()-eta[iT])*(jethighest->Eta()-eta[iT]) +  
			  (jethighest->Phi()-phi[iT])*(jethighest->Phi()-phi[iT]) < 0.2*0.2 )   
		  hJetJT2->Fill( correctp*sin(alpha));
	  if (fNominalEnergy==0.0){
		  hFragmFcn2->Fill(  correctp*sin(tt)/parton->Energy()  );
	  } else
	  {
                  hFragmFcn2->Fill(  correctp*sin(tt)/fNominalEnergy );
	  }
  }// loop over tracks
  }

  if (fOutput->GetNJets()==1)
  {

//========================= CASE 1 ===========================	
  Int_t nPartons = fOutput->GetNPartons();
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
}


