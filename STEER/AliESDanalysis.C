//********************************************************************
//     Example (very naive for the moment) of the data analysis 
//                    using the ESD classes
//       It demonstrates the idea of the "combined PID".
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
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

  #include "AliRun.h"
  #include "AliStack.h"
  #include "AliRunLoader.h"
  #include "AliLoader.h"

  #include "AliESD.h"
#endif

extern AliRun *gAlice;

Int_t AliESDanalysis(Int_t nev=1) { 
  TH2F *tpcHist=
     new TH2F("tpcHist","TPC dE/dX vs momentum",100,0.,4.,100,0.,500.);
   tpcHist->SetXTitle("p (GeV/c)"); tpcHist->SetYTitle("dE/dx (Arb. Units)");
     tpcHist->SetMarkerStyle(8); 
     tpcHist->SetMarkerSize(0.3);
 
  TH1F *piG=new TH1F("piG","",20,0.,4.);
  TH1F *piR=new TH1F("piR","",20,0.,4.);
  TH1F *piF=new TH1F("piF","",20,0.,4.);
  TH1F *piGood=new TH1F("piGood","Combined PID for pions",20,0.,4.); 
  piGood->SetLineColor(4); piGood->SetXTitle("p (GeV/c)");
  TH1F *piFake=new TH1F("piFake","",20,0.,4.); piFake->SetLineColor(2);

  TH1F *kaG=new TH1F("kaG","",20,0.,4.);
  TH1F *kaR=new TH1F("kaR","",20,0.,4.);
  TH1F *kaF=new TH1F("kaF","",20,0.,4.);
  TH1F *kaGood=new TH1F("kaGood","Combined PID for kaons",20,0.,4.); 
  kaGood->SetLineColor(4); kaGood->SetXTitle("p (GeV/c)");
  TH1F *kaFake=new TH1F("kaFake","",20,0.,4.); kaFake->SetLineColor(2);

  TH1F *prG=new TH1F("prG","",20,0.,4.);
  TH1F *prR=new TH1F("prR","",20,0.,4.);
  TH1F *prF=new TH1F("prF","",20,0.,4.);
  TH1F *prGood=new TH1F("prGood","Combined PID for protons",20,0.,4.); 
  prGood->SetLineColor(4); prGood->SetXTitle("p (GeV/c)");
  TH1F *prFake=new TH1F("prFake","",20,0.,4.); prFake->SetLineColor(2);

   if (gAlice) {
      delete gAlice->GetRunLoader();
      delete gAlice;
      gAlice=0;
   }
   AliRunLoader *rl = AliRunLoader::Open("galice.root");
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
   AliStack* stack = rl->Stack();

   TFile *ef=TFile::Open("AliESDs.root");
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

   TCanvas *c1=new TCanvas("c1","",0,0,600,1200);
   c1->Divide(1,4);

   c1->cd(1);
   tpcHist->Draw();

   c1->cd(2);
   piG->Sumw2(); piF->Sumw2(); piR->Sumw2();
   piGood->Divide(piR,piG,1,1,"b");
   piFake->Divide(piF,piG,1,1,"b");
   piGood->Draw("hist");
   piFake->Draw("same");

   c1->cd(3);
   kaG->Sumw2(); kaF->Sumw2(); kaR->Sumw2();
   kaGood->Divide(kaR,kaG,1,1,"b");
   kaFake->Divide(kaF,kaG,1,1,"b");
   kaGood->Draw("hist");
   kaFake->Draw("same");

   c1->cd(4);
   prG->Sumw2(); prF->Sumw2(); prR->Sumw2();
   prGood->Divide(prR,prG,1,1,"b");
   prFake->Divide(prF,prG,1,1,"b");
   prGood->Draw("hist");
   prFake->Draw("same");

   ef->Close();

   delete rl;

   return rc;
}
