//********************************************************************
//     Example (very naive for the moment) of using the ESD classes
//        It demonstrates the idea of the "combined PID".
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//********************************************************************

#if !defined( __CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include "TKey.h"
  #include "TFile.h"
  #include "TH1F.h"
  #include "TH2F.h"
  #include "TCanvas.h"
  #include "TStopwatch.h"
  #include "TParticle.h"
  #include "TROOT.h"

  #include "AliRun.h"
  #include "AliStack.h"
  #include "AliRunLoader.h"
  #include "AliLoader.h"

  #include "AliESD.h"
#endif

extern AliRun *gAlice;
extern TROOT *gROOT;

Int_t AliESDcomparison(const Char_t *dir=".") { 
   TH2F *tpcHist=(TH2F*)gROOT->FindObject("tpcHist");
   if (!tpcHist)
     tpcHist=new TH2F("tpcHist","TPC dE/dX vs momentum",100,0.,4.,100,0.,500.);
   tpcHist->SetXTitle("p (GeV/c)"); tpcHist->SetYTitle("dE/dx (Arb. Units)");
   tpcHist->SetMarkerStyle(8); 
   tpcHist->SetMarkerSize(0.3);
 
   const Char_t *hname[]={
    "piG","piR","piF","piGood","piFake",
    "kaG","kaR","kaF","kaGood","kaFake",
    "prG","prR","prF","prGood","prFake"
   };
   Int_t nh=sizeof(hname)/sizeof(const Char_t *);
   TH1F **hprt=new TH1F*[nh]; 

   for (Int_t i=0; i<nh; i++) {
     hprt[i]=(TH1F*)gROOT->FindObject(hname[i]);
     if (hprt[i]==0) {hprt[i]=new TH1F(hname[i],"",20,0.,4.);hprt[i]->Sumw2();}
   }
   TH1F *piG=hprt[0];
   TH1F *piR=hprt[1];
   TH1F *piF=hprt[2];
   TH1F *piGood=hprt[3]; 
        piGood->SetTitle("Combined PID for pions"); 
        piGood->SetLineColor(4); piGood->SetXTitle("p (GeV/c)");
   TH1F *piFake=hprt[4]; 
        piFake->SetLineColor(2);
   TH1F *kaG=hprt[5];
   TH1F *kaR=hprt[6];
   TH1F *kaF=hprt[7];
   TH1F *kaGood=hprt[8];
        kaGood->SetTitle("Combined PID for kaons"); 
        kaGood->SetLineColor(4); kaGood->SetXTitle("p (GeV/c)");
   TH1F *kaFake=hprt[9]; 
        kaFake->SetLineColor(2);
   TH1F *prG=hprt[10];
   TH1F *prR=hprt[11];
   TH1F *prF=hprt[12];
   TH1F *prGood=hprt[13];
        prGood->SetTitle("Combined PID for protons"); 
        prGood->SetLineColor(4); prGood->SetXTitle("p (GeV/c)");
   TH1F *prFake=hprt[14]; 
        prFake->SetLineColor(2);

   delete[] hprt;

   Char_t fname[100];

   if (gAlice) {
      delete gAlice->GetRunLoader();
      delete gAlice;
      gAlice=0;
   }
   sprintf(fname,"%s/galice.root",dir);
   AliRunLoader *rl = AliRunLoader::Open(fname);
   if (rl == 0x0) {
      cerr<<"Can not open session"<<endl;
      return 1;
   }
   if (rl->LoadgAlice()) {
      cerr<<"LoadgAlice returned error"<<endl;
      delete rl;
      return 1;
   }
   if (rl->LoadHeader()) {
      cerr<<"LoadHeader returned error"<<endl;
      delete rl;
      return 1;
   }
   rl->LoadKinematics();
   AliStack *stack = rl->Stack();

   sprintf(fname,"%s/AliESDs.root",dir);
   TFile *ef=TFile::Open(fname);
   if (!ef->IsOpen()) {cerr<<"Can't AliESDs.root !\n"; return 1;}

   TStopwatch timer;
   Int_t rc=0,n=0;
   TKey *key=0;
   TIter next(ef->GetListOfKeys());

   //****** Tentative particle type "concentrations"
   Double_t c[5]={0.01, 0.01, 0.85, 0.10, 0.05};

   //******* The loop over events
   while ((key=(TKey*)next())!=0) {
     rl->GetEvent(n);
     ef->cd();

     cerr<<"Processing event number : "<<n++<<endl;

     AliESD *event=(AliESD*)key->ReadObj();

     Int_t ntrk=event->GetNumberOfTracks();
     cerr<<"Number of ESD tracks : "<<ntrk<<endl; 
     //****** The loop over tracks

     Int_t pisel=0,kasel=0,prsel=0,nsel=0;

     while (ntrk--) {
       AliESDtrack *t=event->GetTrack(ntrk);

       Double_t p=t->GetP();

       if (t->GetStatus()&AliESDtrack::kTPCin) {
	 Double_t dedx=t->GetTPCsignal();
         tpcHist->Fill(p,dedx,1);
       }

       UInt_t status=AliESDtrack::kESDpid;
       status|=AliESDtrack::kTPCpid; 
       status|=AliESDtrack::kTOFpid; 

       if ((t->GetStatus()&status) == status) {
         nsel++;

         Int_t lab=TMath::Abs(t->GetLabel());
         TParticle *part=stack->Particle(lab);
         Int_t code=part->GetPdgCode();

         Double_t r[10]; t->GetESDpid(r);

         Double_t rcc=0.;
         Int_t i;
         for (i=0; i<AliESDtrack::kSPECIES; i++) rcc+=(c[i]*r[i]);
         if (rcc==0.) continue;

	 //Here we apply Bayes' formula
         Double_t w[10];
         for (i=0; i<AliESDtrack::kSPECIES; i++) w[i]=c[i]*r[i]/rcc;

         if (w[4]>w[3] && w[4]>w[2] && w[4]>w[1] && w[4]>w[0]) {//proton
	    prsel++;
	    prG->Fill(p);
            if (TMath::Abs(code)==2212) prR->Fill(p);
            else prF->Fill(p);
         }

         if (w[3]>w[4] && w[3]>w[2] && w[3]>w[1] && w[3]>w[0]) {//kaon
	    kasel++;
	    kaG->Fill(p);
            if (TMath::Abs(code)==321) kaR->Fill(p);
            else kaF->Fill(p);
         }

	 if (w[2]>w[3] && w[2]>w[4] && w[2]>w[1] && w[2]>w[0]) {//pion
	    pisel++;
	    piG->Fill(p);
	    if (TMath::Abs(code)==211) piR->Fill(p);
            else piF->Fill(p);
         }

       }

     } 
     delete event;
     cerr<<"Number of selected ESD tracks : "<<nsel<<endl;
     cerr<<"Number of selected pion ESD tracks : "<<pisel<<endl;
     cerr<<"Number of selected kaon ESD tracks : "<<kasel<<endl;
     cerr<<"Number of selected proton ESD tracks : "<<prsel<<endl;
   }

   timer.Stop(); timer.Print();

   TCanvas *c1=(TCanvas*)gROOT->FindObject("c1");
   if (c1) delete c1; 
   c1=new TCanvas("c1","",0,0,600,1200);
   c1->Divide(1,4);

   c1->cd(1);
   tpcHist->Draw();

   c1->cd(2);
   //piG->Sumw2(); piF->Sumw2(); piR->Sumw2();
   piGood->Divide(piR,piG,1,1,"b");
   piFake->Divide(piF,piG,1,1,"b");
   piGood->Draw("hist");
   piFake->Draw("same");

   c1->cd(3);
   //kaG->Sumw2(); kaF->Sumw2(); kaR->Sumw2();
   kaGood->Divide(kaR,kaG,1,1,"b");
   kaFake->Divide(kaF,kaG,1,1,"b");
   kaGood->Draw("hist");
   kaFake->Draw("same");

   c1->cd(4);
   //prG->Sumw2(); prF->Sumw2(); prR->Sumw2();
   prGood->Divide(prR,prG,1,1,"b");
   prFake->Divide(prF,prG,1,1,"b");
   prGood->Draw("hist");
   prFake->Draw("same");

   ef->Close();

   delete rl;

   return rc;
}
