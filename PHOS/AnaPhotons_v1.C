//===============================================================
// In this Macro (which can (must) be compilated), you will find all the 
// analysis functions to build photon spectrum, invariant mass 
// spectrum of photon  pairs and combinatorial background calculations 
// in ALICE electromagnetic calorimeter
// Author: Gines MARTINEZ, Subatech, 15 june 2001
//============================================================== 
#include "TH2.h"
#include "TH1.h"
#include "TFile.h"
#include "TRandom.h"
#include "TObjectTable.h"
#include "AliPHOSIndexToObject.h"
#include "AliPHOSRecParticle.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TF1.h"

TObjectTable * gObjectTable;
TRandom * gRandom;


void AnaMinv(char * filename)
{
  TH2F * h_Minv_lowpT       = new TH2F("h_Minv_lowpT","Minv vs pT low",500,0.0,1.0,40,0.,10.);
  TH2F * h_Minv_highpT      = new TH2F("h_Minv_highpT","Minv vs pT high",500,0.0,1.0,50,0.,100.);
  TH2F * h_Minv_lowpT_back  = new TH2F("h_Minv_lowpT_back","Minv vs pT low back",500,0.0,1.0,40,0.,10.);
  TH2F * h_Minv_highpT_back = new TH2F("h_Minv_highpT_back","Minv vs pT high back",500,0.0,1.0,50,0.,100.);


  TH1F * h_Pseudoeta = new TH1F("h_Pseudoeta","Pseudoeta photons",500,-1.0,1.0);
  TH1F * h_Pt        = new TH1F("h_Pt","Pt photons",400,0.,10.);
  TH2F * h_Peta_Pt   = new TH2F("h_Peta_Pt","Pseudo vs pT",40,0.,10.,50,-1.0,1.0);
  TH1F * h_Phi       = new TH1F("h_Phi","Phi photons",400,-4.,4.);
  TH2F * h_Peta_Phi  = new TH2F("h_Peta_Phi","Pseudo vs Phi",200,-4,4,200,-1.0,1.0);

  TH1F * h_DeltaR   = new TH1F("h_DeltaR","Delta R",400,0,4);
  TH1F * h_Asymmetry= new TH1F("h_Asymmetry","Asymmetry",400, -2., 2.);

  AliPHOSIndexToObject * RecData = AliPHOSIndexToObject::GetInstance(filename)  ;
 
  AliPHOSRecParticle * RecParticle1;
  AliPHOSRecParticle * RecParticle2;
 

  Float_t RelativeRCut = 0.0001 ; 
  Float_t AsymmetryCut = 0.7 ;
  Float_t Asymmetry;

  Int_t iEvent, iRecParticle1, iRecParticle2;
  Int_t nRecParticle;
  Float_t invariant_mass, invariant_mass_mixed;
  TLorentzVector P_photon1, P_photon2, P_photonMixed1, P_photonMixed2 ;
  Float_t average_multiplicity = 0.;

  for(iEvent=0; iEvent<RecData->GetMaxEvent(); iEvent++)
    //  for(iEvent=0; iEvent<1000; iEvent++)
    {
      //     if (iEvent==2) gObjectTable->Print();
      //if (iEvent==15) gObjectTable->Print();
      RecData->GetEvent(iEvent);
      printf(">>> Event %d \n",iEvent);
      nRecParticle=RecData->GimeNRecParticles();
      average_multiplicity += ((Float_t) (nRecParticle) ) / ( (Float_t)RecData->GetMaxEvent() ) ;
      // Construction de la masse invariante des pairs
      if (nRecParticle > 1) 
	{
	  for(iRecParticle1=0; iRecParticle1<nRecParticle; iRecParticle1++)
	    {
	      RecParticle1 = (AliPHOSRecParticle *)  RecData->GimeRecParticle(iRecParticle1);
	      RecParticle1->Momentum(P_photon1);

	      h_Pseudoeta->Fill(P_photon1.PseudoRapidity());
	      h_Pt->Fill(P_photon1.Pt());
	      h_Phi->Fill(P_photon1.Phi());
	      h_Peta_Pt->Fill(P_photon1.Pt(), P_photon1.PseudoRapidity());
	      h_Peta_Phi->Fill(P_photon1.Phi(), P_photon1.PseudoRapidity() );
	 
 	      for(iRecParticle2=iRecParticle1+1; iRecParticle2<nRecParticle; iRecParticle2++)
 		{
 		  RecParticle2 = (AliPHOSRecParticle *)  RecData->GimeRecParticle(iRecParticle2);
 		  RecParticle2->Momentum(P_photon2); 
		  Asymmetry = TMath::Abs((P_photon1.E()-P_photon2.E())/(P_photon1.E()+P_photon2.E()));
  		  if ( (P_photon1 != P_photon2) && 
		       (P_photon1.DeltaR(P_photon2) > RelativeRCut) &&
		       (Asymmetry < AsymmetryCut)                          )
  		    {
		      h_DeltaR->Fill(P_photon1.DeltaR(P_photon2));
		      h_Asymmetry->Fill( Asymmetry );

		      //   printf("A. p1 es %f \n",P_photon1->E());
  		      invariant_mass = (P_photon1 + P_photon2).M();
		      // printf("B. p1 es %f \n",P_photon1->E());
  		      h_Minv_lowpT->Fill(invariant_mass, (P_photon1 + P_photon2).Pt() );
  		      h_Minv_highpT->Fill(invariant_mass,(P_photon1 + P_photon2).Pt() );
 		    }  
 		}
 	    }
	}
    }
  printf(">>> Average Multiplicity is %f \n",average_multiplicity);
  Int_t Background = (Int_t)  (RecData->GetMaxEvent() * average_multiplicity * (average_multiplicity-1.)/2.) ;
  printf(">>> Background is %d \n",Background);

  Double_t Pt_Mixed1, Pt_Mixed2;
  Double_t Y_Mixed1, Y_Mixed2;
  Double_t Phi_Mixed1, Phi_Mixed2;

  for(iEvent=0; iEvent<Background; iEvent++)
    {
      //      printf(">>> Background Event %d \n",iEvent);
      Pt_Mixed1 =  h_Pt->GetRandom(); 
      Pt_Mixed2 =  h_Pt->GetRandom();
      h_Peta_Phi->GetRandom2(Phi_Mixed1, Y_Mixed1);
      h_Peta_Phi->GetRandom2(Phi_Mixed2, Y_Mixed2);
      P_photonMixed1.SetPtEtaPhiM( Pt_Mixed1, Y_Mixed1, Phi_Mixed1, 0.0);
      P_photonMixed2.SetPtEtaPhiM( Pt_Mixed2, Y_Mixed2, Phi_Mixed2, 0.0);
      Asymmetry = TMath::Abs((P_photonMixed1.E()-P_photonMixed2.E())/(P_photonMixed1.E()+P_photonMixed2.E()));
      
      if ( (P_photonMixed1.DeltaR(P_photonMixed2) > RelativeRCut) &&
           (Asymmetry < AsymmetryCut  )                               )
	{
	  invariant_mass_mixed = (P_photonMixed1 + P_photonMixed2).M();
	  h_Minv_lowpT_back->Fill(invariant_mass_mixed, (P_photonMixed1 + P_photonMixed2).Pt() );
	  h_Minv_highpT_back->Fill(invariant_mass_mixed,(P_photonMixed1 + P_photonMixed2).Pt() );
	}  
    }
  

  char outputname[80];
  sprintf(outputname,"%s.Minv",filename);
  TFile output(outputname,"recreate");
  h_Minv_lowpT->Write();
  h_Minv_highpT->Write();
  h_Minv_lowpT_back->Write();  
  h_Minv_highpT_back->Write();
  h_Pseudoeta->Write();
  h_Pt->Write();
  h_Peta_Pt->Write();
  h_Phi->Write();
  h_Peta_Phi->Write();
  h_Asymmetry->Write();
  h_DeltaR->Write();
  
  output.Close();
}


void AnaPtSpectrum(char * filename, Int_t NumberPerPtBin, Option_t * particle, Option_t * opt)
{

  Int_t NumberOfPtBins = NumberPerPtBin;
  Float_t PtCalibration = 0.250;

  TFile * in = new TFile(filename);

  TH2F * h_Minv_pT = 0;  
  TH2F * h_Minv_pT_back = 0; 
  TH2F * frame = 0 ;

  if (strstr(opt,"low"))
    {
      h_Minv_pT      = (TH2F *) in->Get("h_Minv_lowpT"); ;
      h_Minv_pT_back = (TH2F *) in->Get("h_Minv_lowpT_back");
      PtCalibration  = 0.250;   
      frame = new TH2F("PtSpectrumlow","Pt Spectrum low",10, 0.,10.,10,0.1,10000);
    }  
  if (strstr(opt,"high"))
    {
      h_Minv_pT      = (TH2F *) in->Get("h_Minv_highpT"); ;
      h_Minv_pT_back = (TH2F *) in->Get("h_Minv_highpT_back");
      PtCalibration = 2.5;
      frame = new TH2F("PtSpectrumhigh","Pt Spectrum high",10, 0.,100.,10,0.1,10000);
    }

  if ( h_Minv_pT == 0 ) 
    {
      printf(">>> Bad Option! \n");
      return;
    }
  Int_t Norma_1 = 100; Float_t Norma_minv_1 = 0.2;
  Int_t Norma_2 = 200; Float_t Norma_minv_2 = 0.4;

  Int_t Minv_1 = 56;
  Int_t Minv_2 = 76;
  if (strstr(particle,"eta"))
    {
      Minv_1 = 234;
      Minv_2 = 314;
    }

  if (strstr(particle,"norma"))
    {
      Minv_1 = 100;
      Minv_2 = 200;
    }
  
  Int_t NHistos = 40/NumberOfPtBins;
  Int_t iHistos;

  TH1D * signal = 0;
  TH1D * background = 0;
  TH1D * ratio = 0;
  TH1D * difference = 0;

  Float_t Pt[NHistos];
  Float_t PtError[NHistos];
  Float_t Nmesons[NHistos];
  Float_t NmesonsError[NHistos];

  Float_t Ntota, Nback, Norma, NormaError, Renorma;

  for(iHistos=0; iHistos<NHistos; iHistos++)
    {
      signal     = h_Minv_pT->ProjectionX("signal",         NumberOfPtBins*iHistos+1,NumberOfPtBins*(iHistos+1));
      background = h_Minv_pT_back->ProjectionX("background",NumberOfPtBins*iHistos+1,NumberOfPtBins*(iHistos+1));
      //signal->Rebin();
      //background->Rebin();
      ratio = new TH1D(*signal);
      ratio->Sumw2(); 
      ratio->Add(background,-1.0);
      ratio->Divide(background);
      difference = new TH1D(*signal);
      difference->Sumw2();
      ratio->Fit("pol0","","",Norma_minv_1,Norma_minv_2); 
      if (background->Integral(Norma_1,Norma_2) == 0)
	Renorma = 0.;
      else
	Renorma = signal->Integral(Norma_1,Norma_2)/background->Integral(Norma_1,Norma_2);
      difference->Add(background,(-1.)*Renorma); 
       
      //ratio->Draw();
      //      background->Draw("same");
      //     difference->Draw();

      Ntota = signal->Integral(Minv_1,Minv_2);
      Nback = background->Integral(Minv_1,Minv_2);
      Norma = ratio->GetFunction("pol0")->GetParameter(0);
      if (Renorma == 0.)
	NormaError = 0.;
      else
	NormaError =  ratio->GetFunction("pol0")->GetParError(0);
      printf("Ntotal %f Nback %f Norma %f and NormaError %f \n",Ntota, Nback, Norma, NormaError);
      printf("differencia is %f \n",difference->Integral(Minv_1,Minv_2));
      Nmesons[iHistos] = Ntota - Renorma * Nback;
      NmesonsError[iHistos] = TMath::Sqrt( Ntota + Nback*Renorma*Renorma + Nback*Nback*NormaError*NormaError );  
      Pt[iHistos] = (iHistos+0.5)*NumberOfPtBins*PtCalibration;
      PtError[iHistos] = NumberOfPtBins*PtCalibration/2.; 
      //   ratio->Delete("");
      //difference->Delete("");
    }
  // in->Close();


  char filenameout[80];
  sprintf(filenameout,"%s.PtSpectrum_%d_%s_%s",filename, NumberPerPtBin, particle, opt);
  TFile out(filenameout,"recreate");
  TGraphErrors * PtSpectrum = new TGraphErrors(NHistos, Pt, Nmesons, PtError, NmesonsError);
  PtSpectrum->SetName("PtSpectrum");
  PtSpectrum->Write();
  out.Close();

  frame->Draw();
  frame->SetStats(0);
  frame->SetXTitle("Neutral meson pT (GeV/c)");
  frame->SetYTitle("Number of neutral mesons per pT bin");
  PtSpectrum->SetMarkerStyle(27);
  PtSpectrum->Draw("P");

}
