//********************************************************************
//     Example (very naive for the moment) of the data analysis 
//                    using the ESD classes.
//       It demonstrates the idea of the "combined PID" 
//            applied to the Lambda0 reconstruction. 
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//********************************************************************

#if !defined( __CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include <TTree.h>
  #include "TFile.h"
  #include "TH1F.h"
  #include "TH2F.h"
  #include "TCanvas.h"
  #include "TStopwatch.h"
  #include "TParticle.h"

  #include "AliRun.h"

  #include "AliESD.h"

#endif

extern AliRun *gAlice;

Int_t AliESDv0Analysis(Int_t nev=1) { 
   TH1F *hm=new TH1F("hm","Effective Mass",40,1.065,1.165);
   hm->SetXTitle("Mass (GeV/c**2)");

   TFile *ef=TFile::Open("AliESDs.root");
   if (!ef->IsOpen()) {cerr<<"Can't AliESDs.root !\n"; return 1;}
   AliESD* event = new AliESD;
   TTree* tree = (TTree*) ef->Get("esdTree");
   if (!tree) {cerr<<"no ESD tree found\n"; return 1;};
   tree->SetBranchAddress("ESD", &event);

   TStopwatch timer;
   Int_t rc=0,n=0;

   //****** Tentative particle type "concentrations"
   Double_t c[5]={0.0, 0.0, 0.1, 0.1, 0.1};

   //******* The loop over events
   while (tree->GetEvent(n)) {

     cerr<<"Processing event number : "<<n++<<endl;

     Int_t nv0=event->GetNumberOfV0s();
     cerr<<"Number of ESD v0s : "<<nv0<<endl; 

     while (nv0--) {
       AliESDv0 *v0=event->GetV0(nv0);
       Int_t pi=v0->GetPindex();
       AliESDtrack *t=event->GetTrack(pi);
       Int_t isProton=1;
       if ((t->GetStatus()&AliESDtrack::kESDpid)!=0) {
	 Double_t r[10]; t->GetESDpid(r);
         Double_t rcc=0.;
         Int_t i;
         for (i=0; i<AliESDtrack::kSPECIES; i++) rcc+=(c[i]*r[i]);
         if (rcc==0.) continue;
         //Here we apply Bayes' formula
         Double_t w[10];
         for (i=0; i<AliESDtrack::kSPECIES; i++) w[i]=c[i]*r[i]/rcc;

         if (w[4]<w[3]) isProton=0;
         if (w[4]<w[2]) isProton=0;
         if (w[4]<w[1]) isProton=0;
	 if (w[4]<w[0]) isProton=0;
       }
       if (!isProton) continue;
       v0->ChangeMassHypothesis(3122);
       Double_t mass=v0->GetEffMass();
       hm->Fill(mass);
     } 
   }

   delete event;
   ef->Close();

   timer.Stop(); timer.Print();

   hm->Draw();

   ef->Close();


   return rc;
}
