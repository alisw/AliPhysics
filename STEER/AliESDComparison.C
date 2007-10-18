//********************************************************************
//     Example (very naive for the moment) of using the ESD classes
//        It demonstrates the idea of the "combined PID".
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//********************************************************************

#if !defined( __CINT__) || defined(__MAKECINT__)
  #include <TMath.h>
  #include <TError.h>
  #include <Riostream.h>
  #include <TH1F.h>
  #include <TH2F.h>
  #include <TParticle.h>
  #include <TCanvas.h>
  #include <TBenchmark.h>
  #include <TFile.h>
  #include <TTree.h>
  #include <TROOT.h>

  #include <AliStack.h>
  #include <AliRunLoader.h>
  #include <AliRun.h>
  #include <AliESDEvent.h>
#endif

extern AliRun *gAlice;
extern TBenchmark *gBenchmark;
extern TROOT *gROOT;

static Int_t allpisel=0;
static Int_t allkasel=0;
static Int_t allprsel=0;
static Int_t allnsel=0;

Int_t AliESDComparison(const Char_t *dir=".") { 
   gBenchmark->Start("AliESDComparison");

   ::Info("AliESDComparison.C","Doing comparison...");

   Double_t pi=0.2,pa=3;

   TH2F *tpcHist=(TH2F*)gROOT->FindObject("tpcHist");
   if (!tpcHist)
     tpcHist=new TH2F("tpcHist","TPC dE/dX vs momentum",100,pi,pa,100,0.,300.);
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
     if (hprt[i]==0) {hprt[i]=new TH1F(hname[i],"",20,pi,pa);hprt[i]->Sumw2();}
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
      ::Error("AliESDComparison.C","Can not open session !");
      return 1;
   }
   if (rl->LoadgAlice()) {
      ::Error("AliESDComparison.C","LoadgAlice returned error !");
      delete rl;
      return 1;
   }
   if (rl->LoadHeader()) {
      ::Error("AliESDComparison.C","LoadHeader returned error !");
      delete rl;
      return 1;
   }
   rl->LoadKinematics();

   sprintf(fname,"%s/AliESDs.root",dir);
   TFile *ef=TFile::Open(fname);
   if (!ef || !ef->IsOpen()) {
      ::Error("AliESDComparison.C","Can't AliESDs.root !"); 
      delete rl;
      return 1;
   }
   AliESDEvent* event = new AliESDEvent();
   TTree* tree = (TTree*) ef->Get("esdTree");
   if (!tree) {
      ::Error("AliESDComparison.C", "no ESD tree found");
      delete rl;
      return 1;
   }
   event->ReadFromTree(tree);

   //****** Tentative particle type "concentrations"
   Double_t c[5]={0.01, 0.01, 0.85, 0.10, 0.05};
   //Double_t c[5]={0.2, 0.2, 0.2, 0.2, 0.2};
   AliPID::SetPriors(c);


   //******* The loop over events
   Int_t e=0;
   while (tree->GetEvent(e)) {
      cout<<endl<<endl<<"********* Processing event number: "<<e<<"*******\n";

      rl->GetEvent(e);
 
      e++;

      Int_t ntrk=event->GetNumberOfTracks();
      cerr<<"Number of ESD tracks : "<<ntrk<<endl; 

      Int_t pisel=0,kasel=0,prsel=0,nsel=0;

      AliStack *stack = rl->Stack();

      while (ntrk--) {
        AliESDtrack *t=event->GetTrack(ntrk);

	//*** Some track quality cuts ****
        if (t->GetITSclusters(0) < 6 ) continue;
        if (t->GetTPCclusters(0) < 60) continue;
        //if (t->GetTRDclusters(0) < 60) continue;

        if (!t->IsOn(AliESDtrack::kESDpid)) continue;
        if (!t->IsOn(AliESDtrack::kITSpid)) continue;
        if (!t->IsOn(AliESDtrack::kTPCpid)) continue;
        //if (!t->IsOn(AliESDtrack::kTRDpid)) continue;
        if (!t->IsOn(AliESDtrack::kTOFpid)) continue;
        {
           nsel++;

           Double_t p=t->GetP();
	   Double_t dedx=t->GetTPCsignal();
           tpcHist->Fill(p,dedx,1);

           Int_t lab=TMath::Abs(t->GetLabel());
           TParticle *part=stack->Particle(lab);
           Int_t code=part->GetPdgCode();

           Double_t r[10]; t->GetESDpid(r);
           //t->GetITSpid(r);
           //t->GetTPCpid(r);
           //t->GetTRDpid(r);
           //t->GetTOFpid(r);

           AliPID pid(r);

           Double_t w[10];
           w[0]=pid.GetProbability(AliPID::kElectron);
           w[1]=pid.GetProbability(AliPID::kMuon);
           w[2]=pid.GetProbability(AliPID::kPion);
           w[3]=pid.GetProbability(AliPID::kKaon);
           w[4]=pid.GetProbability(AliPID::kProton);


           if (TMath::Abs(code)==2212) prR->Fill(p);
           if (w[4]>w[3] && w[4]>w[2] && w[4]>w[1] && w[4]>w[0]) {//proton
	      prsel++;
	      prG->Fill(p);
              if (TMath::Abs(code)!=2212) prF->Fill(p);
           }

           if (TMath::Abs(code)==321) kaR->Fill(p);
           if (w[3]>w[4] && w[3]>w[2] && w[3]>w[1] && w[3]>w[0]) {//kaon
	      kasel++;
	      kaG->Fill(p);
              if (TMath::Abs(code)!=321) kaF->Fill(p);
           }

	   if (TMath::Abs(code)==211) piR->Fill(p);
	   if (w[2]>w[4] && w[2]>w[3] && w[2]>w[0] && w[2]>w[1]) {//pion
	      pisel++;
	      piG->Fill(p);
	      if (TMath::Abs(code)!=211) piF->Fill(p);
           }
	}
      }
      cout<<"Number of selected ESD tracks : "<<nsel<<endl;
      cout<<"Number of selected pion ESD tracks : "<<pisel<<endl;
      cout<<"Number of selected kaon ESD tracks : "<<kasel<<endl;
      cout<<"Number of selected proton ESD tracks : "<<prsel<<endl;

      allnsel+=nsel; allpisel+=pisel; allkasel+=kasel; allprsel+=prsel;

   } // ***** End of the loop over events

   delete event;
   delete tree;
   ef->Close();

   TCanvas *c1=(TCanvas*)gROOT->FindObject("c1");
   if (!c1) {
      c1=new TCanvas("c1","",0,0,600,1200);
      c1->Divide(1,4);
   }

   c1->cd(1);
   tpcHist->Draw();

   c1->cd(2);
   //piG->Sumw2(); piF->Sumw2(); piR->Sumw2();
   piGood->Add(piG,piF,1,-1);
   piGood->Divide(piGood,piR,1,1,"b");
   piFake->Divide(piF,piG,1,1,"b");
   piGood->Draw("hist");
   piFake->Draw("same");

   c1->cd(3);
   //kaG->Sumw2(); kaF->Sumw2(); kaR->Sumw2();
   kaGood->Add(kaG,kaF,1,-1);
   kaGood->Divide(kaGood,kaR,1,1,"b");
   kaFake->Divide(kaF,kaG,1,1,"b");
   kaGood->Draw("hist");
   kaFake->Draw("same");

   c1->cd(4);
   //prG->Sumw2(); prF->Sumw2(); prR->Sumw2();
   prGood->Add(prG,prF,1,-1);
   prGood->Divide(prGood,prR,1,1,"b");
   prFake->Divide(prF,prG,1,1,"b");
   prGood->Draw("hist");
   prFake->Draw("same");

   c1->Update();

   cout<<endl;
   cout<<endl;
   e=(Int_t)piR->GetEntries();
   Int_t o=(Int_t)piG->GetEntries();
   if (e*o) cout<<"Efficiency (contamination) for pions : "<<
      piGood->GetEntries()/e<<'('<<piF->GetEntries()/o<<')'<<endl; 
   e=(Int_t)kaR->GetEntries();
   o=(Int_t)kaG->GetEntries();
   if (e*o) cout<<"Efficiency (contamination) for kaons : "<<
      kaGood->GetEntries()/e<<'('<<kaF->GetEntries()/o<<')'<<endl; 
   e=(Int_t)prR->GetEntries();
   o=(Int_t)prG->GetEntries();
   if (e*o) cout<<"Efficiency (contamination) for protons : "<<
      prGood->GetEntries()/e<<'('<<prF->GetEntries()/o<<')'<<endl;
   cout<<endl; 

   cout<<"Total number of selected ESD tracks : "<<allnsel<<endl;
   cout<<"Total number of selected pion ESD tracks : "<<allpisel<<endl;
   cout<<"Total number of selected kaon ESD tracks : "<<allkasel<<endl;
   cout<<"Total number of selected proton ESD tracks : "<<allprsel<<endl;

   ef->Close();
   TFile fc("AliESDComparison.root","RECREATE");
   c1->Write();
   fc.Close();

   gBenchmark->Stop("AliESDComparison");
   gBenchmark->Show("AliESDComparison");

   delete rl;
   return 0;
}
