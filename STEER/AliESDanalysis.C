//********************************************************************
//     Example (very naive for the moment) of the data analysis 
//                    using the ESD classes
//       It demonstrates the idea of the "combined PID".
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//********************************************************************

#ifndef __CINT__
  #include <Riostream.h>
  #include "TKey.h"
  #include "TFile.h"
  #include "TH2F.h"
  #include "TCanvas.h"
  #include "TStopwatch.h"

  #include "AliESD.h"
#endif

Int_t AliESDanalysis(Int_t nev=1) { 
  TH2F *tpcHist=
     new TH2F("tpcHist","TPC dE/dX vs momentum",100,0.,4.,100,0.,500.);
  tpcHist->SetXTitle("p (GeV/c)"); tpcHist->SetYTitle("dE/dx (Arb. Units)");
     tpcHist->SetMarkerStyle(8); 
     tpcHist->SetMarkerSize(0.3);
 
  TH2F *elHist=new TH2F("elHist","dE/dX vs momentum",100,0.,4.,50,0.,500.);
     elHist->SetMarkerStyle(8); 
     elHist->SetMarkerSize(0.3);

  TH2F *piHist=new TH2F("piHist","dE/dX vs momentum",100,0.,4.,50,0.,500.);
     piHist->SetMarkerColor(2); 
     piHist->SetMarkerStyle(8); 
     piHist->SetMarkerSize(0.3);

  TH2F *kaHist=new TH2F("kaHist","dE/dX vs momentum",100,0.,4.,100,0.,500.);
     kaHist->SetMarkerColor(4); 
     kaHist->SetMarkerStyle(8); 
     kaHist->SetMarkerSize(0.3);

  TH2F *prHist=new 
   TH2F("prHist","Classification into e, pi, K and p",100,0.,4.,100,0.,500.);
   prHist->SetXTitle("p (GeV/c)"); prHist->SetYTitle("dE/dx (Arb. Units)");
     prHist->SetMarkerColor(6); 
     prHist->SetMarkerStyle(8); 
     prHist->SetMarkerSize(0.3);
   
   TFile *ef=TFile::Open("AliESDs.root");
   if (!ef->IsOpen()) {cerr<<"Can't AliESDs.root !\n"; return 1;}

   TStopwatch timer;
   Int_t rc=0,n=0;
   TKey *key=0;
   TIter next(ef->GetListOfKeys());

   //****** Tentative particle type "concentrations"
   Double_t c[5]={0.05, 0., 0.85, 0.10, 0.05};

   //******* The loop over events
   while ((key=(TKey*)next())!=0) {
     cerr<<"Processing event number : "<<n++<<endl;
     AliESD *event=(AliESD*)key->ReadObj();

     Int_t ntrk=event->GetNumberOfTracks();
     cerr<<"Number of ESD tracks : "<<ntrk<<endl; 
     //****** The loop over tracks
     while (ntrk--) {
       AliESDtrack *t=event->GetTrack(ntrk);

       Double_t p=t->GetP();

       if (t->GetStatus()&AliESDtrack::kTPCin) {
	 Double_t dedx=t->GetTPCsignal();
         tpcHist->Fill(p,dedx,1);
       }

       if (t->GetStatus()&AliESDtrack::kESDpid) {
         Double_t dedx=t->GetTPCsignal();
	 Double_t r[10]; t->GetESDpid(r);
         Double_t rc=0.;
         Int_t i;
         for (i=0; i<AliESDtrack::kSPECIES; i++) rc+=(c[i]*r[i]);
         if (rc==0.) continue;

	 //Here we apply Bayes' formula
         Double_t w[10];
         for (i=0; i<AliESDtrack::kSPECIES; i++) w[i]=c[i]*r[i]/rc;

         if (w[4]>w[3] && w[4]>w[2] && w[4]>w[0]) prHist->Fill(p,dedx,1);
         if (w[3]>w[4] && w[3]>w[2] && w[3]>w[0]) kaHist->Fill(p,dedx,1);
         if (w[2]>w[3] && w[2]>w[4] && w[2]>w[0]) piHist->Fill(p,dedx,1);
         if (w[0]>w[3] && w[0]>w[2] && w[0]>w[4]) elHist->Fill(p,dedx,1);
       }

     } 
     delete event;
   }
   timer.Stop(); timer.Print();

   TCanvas *c1=new TCanvas("c1","",0,0,600,1200);
   c1->Divide(1,2);

   c1->cd(1);
   tpcHist->Draw();
   c1->cd(2);
   prHist->Draw();
   kaHist->Draw("same");
   piHist->Draw("same");
   elHist->Draw("same");

   ef->Close();

   return rc;
}
