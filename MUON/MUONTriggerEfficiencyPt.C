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


/* Macro to produce trigger single muon efficiency versus pt plots for the 
   3 pt cuts. Results are compared to the reference (red curves).   
   To be used with (at least) 10000 events as follows
   AliGenBox * gener = new AliGenBox(1);
   gener->SetPtRange(0.,10.);
   gener->SetPhiRange(0., 360.);         
   gener->SetThetaRange(171.000,178.001);
   gener->SetPart(13);           // or -13
   gener->SetOrigin(0.,0., 0.);  
   gener->SetSigma(0.0, 0.0, 0.0);     

   Author: Fabien Guerin, LPC Clermont-Ferrand, Jan. 2006
*/


// ROOT includes
#include "TBranch.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TParticle.h"
#include "TTree.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"


// STEER includes
#include "AliRun.h"
#include "AliDetector.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliLoader.h"
#include "AliStack.h"
#include "AliSegmentation.h"
#include "AliMC.h"

// MUON includes
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONHit.h"
#include "AliMUONChamber.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONTrack.h"

Double_t fitArc(Double_t *x,Double_t *par) 
{
  Double_t h1=par[1]*(x[0]+par[2]);
  return par[0]*TMath::ATan(TMath::Power(h1,6))+par[3];
}

Double_t fitArch(Double_t *x,Double_t *par) 
{
  Double_t h0=2*par[0]/(TMath::Pi());
  Double_t h1=par[1]*x[0]+par[2]*x[0]*x[0];
  return h0*TMath::TanH(h1)+par[3];
}

void MUONTriggerEfficiencyPt(char filename[10]="galice.root")
{

// define style
    TStyle *st1 = new TStyle("st1","My STYLE");
    
    st1->SetOptStat(0);
    st1->SetOptFit(111);  
    st1->SetFrameFillColor(10);
    st1->SetCanvasColor(10);
    st1->SetFillColor(10);
    st1->SetTitleBorderSize(0);
    st1->SetTitleOffset(1.2,"XY");
    st1->SetTitleSize(0.05,"XY"); 
    st1->SetLabelSize(0.045,"XY");  
    st1->SetLabelFont(22,"XY");  
    st1->SetTitleFont(22,"XY");     
    st1->SetOptTitle(0);   
    st1->SetStatColor(10);
    st1->SetStatFont(62);     
    st1->SetFitFormat("4.4g");
    st1->SetPadLeftMargin(0.15);
    st1->SetPadRightMargin(0.06); 
    st1->SetPadTopMargin(0.06);
    st1->SetPadBottomMargin(0.15); 
    st1->cd();
    
    gROOT->ForceStyle();
    //TGaxis::SetMaxDigits(3);
    
// beginning of macro    
    TH1F *ptlpt    = new TH1F("ptlpt","",50,0.15,5.); 
    TH1F *pthpt    = new TH1F("pthpt","",50,0.15,5.);
    TH1F *ptcoinc  = new TH1F("ptcoinc","",50,0.15,5.);  
        
    TParticle *particle;
    AliStack* stack; 
    
    Double_t coincmuon=0;
    Double_t lptmuon=0;
    Double_t hptmuon=0;
    Double_t ptmu=0;

// output file
    char digitdat[100];
    sprintf(digitdat,"MUONTriggerEfficiencyPt.out");  
    FILE *fdat=fopen(digitdat,"w");

// Creating Run Loader and openning file containing Hits
    AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
    
    if (RunLoader ==0x0) {
	printf(">>> Error : Error Opening %s file \n",filename);
	return;
    }
    
    Int_t nevents = RunLoader->GetNumberOfEvents();
    
// loading hits
    AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");     
    MUONLoader->LoadHits("READ");
    AliMUONData data_hits(MUONLoader,"MUON","MUON");
    
    TClonesArray * globalTrigger;
    AliMUONGlobalTrigger * gloTrg;

// loading trigger
    MUONLoader->LoadDigits("READ");  
    AliMUONData muondata(MUONLoader,"MUON","MUON");

// Loading kine
    RunLoader->LoadKinematics("READ");
    
    for (Int_t ievent=0; ievent<nevents; ievent++) {  // Event loop
	
	RunLoader->GetEvent(ievent);
	
	if (ievent%500==0) printf("ievent = %d \n",ievent);


// kine
	Int_t iparticle, nparticles;
	stack = RunLoader->Stack();
	nparticles = (Int_t) stack->GetNtrack();        
	for (iparticle=0; iparticle<nparticles; iparticle++) {             
            particle = stack->Particle(iparticle);   
            Float_t pt = particle->Pt();	    
            Int_t pdgcode = particle->GetPdgCode();                        
	    if (pdgcode==-13 || pdgcode==13) ptmu = pt;             
        }

// trigger 
        muondata.SetTreeAddress("D,GLT");
        muondata.GetTriggerD();
    
        globalTrigger = muondata.GlobalTrigger();

        Int_t nglobals = (Int_t) globalTrigger->GetEntriesFast(); // should be 1 
        for (Int_t iglobal=0; iglobal<nglobals; iglobal++) { // Global Trigger
	    gloTrg = static_cast<AliMUONGlobalTrigger*>(globalTrigger->At(iglobal));
	    
	    if (gloTrg->SingleLpt()>=1 ) {
		lptmuon++;
		ptlpt->Fill(ptmu);
	    }       
	    if (gloTrg->SingleHpt()>=1 ) {
		hptmuon++;
		pthpt->Fill(ptmu);
	    }
        } // end of loop on Global Trigger      
        muondata.ResetTrigger();  

// Hits
        RunLoader->GetEvent(ievent);  
        data_hits.SetTreeAddress("H");    
        
        Int_t itrack, ntracks, NbHits[4];
        Int_t SumNbHits;

        for (Int_t j=0; j<4; j++) NbHits[j]=0;
 
        ntracks = (Int_t) data_hits.GetNtracks();      
     
        for (itrack=0; itrack<ntracks; itrack++) { // track loop
         data_hits.GetTrack(itrack); 

         Int_t ihit, nhits;
         nhits = (Int_t) data_hits.Hits()->GetEntriesFast();   
         AliMUONHit* mHit;

    
          for (ihit=0; ihit<nhits; ihit++) {
            mHit = static_cast<AliMUONHit*>(data_hits.Hits()->At(ihit));
            Int_t Nch        = mHit->Chamber(); 
            Int_t hittrack   = mHit->Track();
            Int_t IdPart     = mHit->Particle();	    

            for (Int_t j=0;j<4;j++) {
              Int_t kch=11+j;
	      if (hittrack==0) {
                if (Nch==kch && (IdPart==-13 || IdPart==13) && NbHits[j]==0) NbHits[j]++; 
		}               
            }    
          }
          data_hits.ResetHits();
        } // end track loop
        
// 3/4 coincidence 
        SumNbHits=NbHits[0]+NbHits[1]+NbHits[2]+NbHits[3];

        if (SumNbHits==3 || SumNbHits==4) {
            coincmuon++;
            ptcoinc->Fill(ptmu);
        }            
                
      } // end loop on event
      
      MUONLoader->UnloadHits();
      MUONLoader->UnloadDigits();
      MUONLoader->UnloadRecPoints();
      RunLoader->UnloadKinematics();   

      delete RunLoader;
     
      fprintf(fdat,"\n");    
      fprintf(fdat,"\n");
      fprintf(fdat," Number of events = %d \n",nevents);  
      fprintf(fdat," Number of events with 3/4 coinc = %d \n",(Int_t)coincmuon);
      fprintf(fdat," Nomber of dimuons with 3/4 coinc Lpt cut = %d \n",(Int_t)lptmuon); 
      fprintf(fdat," Number of dimuons with 3/4 coinc Hpt cut = %d \n",(Int_t)hptmuon);  
      fprintf(fdat,"\n"); 
      
      Double_t efficiency,error;

      efficiency=lptmuon/coincmuon;  
      error=efficiency*TMath::Sqrt((lptmuon+coincmuon)/(lptmuon*coincmuon));  
      fprintf(fdat," Efficiency Lpt cut = %4.4f +/- %4.4f\n",efficiency,error);
      
      efficiency=hptmuon/coincmuon; 
      error=efficiency*TMath::Sqrt((hptmuon+coincmuon)/(hptmuon*coincmuon)); 
      fprintf(fdat," Efficiency Hpt cut = %4.4f +/- %4.4f\n",efficiency,error);
      
      fclose(fdat);
      
      Double_t x1,x2,xval,xerr;
      
      for (Int_t i=0;i<50;i++) {
	  x1=ptcoinc->GetBinContent(i+1);
	  
	  x2=ptlpt->GetBinContent(i+1);
	  if (x1!=0 && x2!=0) {
	      xval=x2/x1;
	      xerr=xval*TMath::Sqrt((x1+x2)/x1*x2);   
	      ptlpt->SetBinContent(i+1,xval);
	      ptlpt->SetBinError(i+1,0);
	  } else {
	      ptlpt->SetBinContent(i+1,0.);
	      ptlpt->SetBinError(i+1,0.);     
	  } 
	  
	  x2=pthpt->GetBinContent(i+1);
	  if (x1!=0 && x2!=0) {
	      xval=x2/x1;
	      xerr=xval*TMath::Sqrt((x1+x2)/x1*x2);   
	      pthpt->SetBinContent(i+1,xval);
	      pthpt->SetBinError(i+1,0);
	  } else {
	      pthpt->SetBinContent(i+1,0.);
	      pthpt->SetBinError(i+1,0.);     
	  }              
      }
      
      
      TF1 *fitlpt = new TF1("fitlpt",fitArc,0.,5.,4); 
      TF1 *fithpt = new TF1("fithpt",fitArc,0.,5.,4);      
      
      TCanvas *c1 = new TCanvas("c1","Trigger efficiency vs pt and y for single muon",200,0,900,400); 
      c1->Divide(2,1);      
      
      c1->cd(1);
      ptlpt->SetMinimum(0.);
      ptlpt->SetMaximum(1.05);   
      ptlpt->SetTitle("");
      ptlpt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
      ptlpt->GetYaxis()->SetTitle("Efficiency");     
      //ptlpt->GetXaxis()->SetRange(3);
      ptlpt->Draw(""); 
      fitlpt->SetParameters(0.602,0.774,0.273,0.048);  
      fitlpt->SetLineColor(2);
      fitlpt->SetLineWidth(2);
      fitlpt->Draw("SAME");    
      TLegend * leg = new TLegend(0.5,0.38,0.85,0.53); 
      leg = new TLegend(0.5,0.38,0.85,0.53); 
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->SetTextFont(22);
      leg->SetHeader("Lpt trigger pt cut");
      leg->AddEntry(fitlpt,"Ref.","l");
      leg->AddEntry(ptlpt,"New","l");                 
      leg->Draw("SAME");
      
      c1->cd(2);
      pthpt->SetMinimum(0.);
      pthpt->SetMaximum(1.05);   
      pthpt->SetTitle("");
      pthpt->GetXaxis()->SetTitle("P_{T} (GeV/c)");  
      pthpt->GetYaxis()->SetTitle("Efficiency");     
      //pthpt->GetXaxis()->SetRange(3);
      pthpt->Draw(""); 
      fithpt->SetParameters(0.627,0.393,0.855,0.0081); 
      fithpt->SetLineColor(2);
      fithpt->SetLineWidth(2);
      fithpt->Draw("SAME");    
      leg = new TLegend(0.5,0.38,0.85,0.53); 
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextSize(0.05);
      leg->SetTextFont(22);
      leg->SetHeader("Hpt trigger pt cut");
      leg->AddEntry(fithpt,"Ref.","l");
      leg->AddEntry(pthpt,"New","l");                 
      leg->Draw("SAME");
      
      c1->SaveAs("MUONTriggerEfficiencyPt.gif");
      c1->SaveAs("MUONTriggerEfficiencyPt.eps");
      
}
