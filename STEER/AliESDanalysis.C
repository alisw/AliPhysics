//********************************************************************
//     Example (very naive for the moment) of the data analysis 
//                    using the ESD classes
//     Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//********************************************************************

#if !defined( __CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include <TTree.h>
  #include "TFile.h"
  #include "TH1F.h"
  #include "TCanvas.h"
  #include "TStyle.h"
  #include "TStopwatch.h"

  #include "AliESD.h"
#endif

extern TStyle *gStyle;

Int_t AliESDanalysis() { 
   TStopwatch timer;

   gStyle->SetOptStat(111110);
   gStyle->SetOptFit(1);

   Double_t V0mass=0.497672, V0width=0.020, V0window=0.05; 
   Double_t mmin=V0mass-V0window, mmax=V0mass+V0window;
   TH1F *hm =new TH1F("hm","K0s",40, mmin, mmax);
   hm->SetXTitle("Mass (GeV/c**2)"); hm->SetLineColor(2);
   TH1F *hp =new TH1F("hp","Momentum of the positive daughter",40, 0, 2);
   hp->SetXTitle("P (GeV/c)"); hp->SetLineColor(4);

//****** File with the ESD
   TFile *ef=TFile::Open("AliESDs.root");
   if (!ef || !ef->IsOpen()) {cerr<<"Can't AliESDs.root !\n"; return 1;}
   AliESD* event = new AliESD;
   TTree* tree = (TTree*) ef->Get("esdTree");
   if (!tree) {cerr<<"no ESD tree found\n"; return 1;};
   tree->SetBranchAddress("ESD", &event);

   Int_t n=0;

//******* The loop over events
   while (tree->GetEvent(n)) {
     cout<<endl<<"Processing event number : "<<n++<<endl;

     Int_t ntrk=event->GetNumberOfTracks();
     cout<<"Number of ESD tracks : "<<ntrk<<endl; 
     Int_t nv0=event->GetNumberOfV0s();
     cout<<"Number of ESD V0s : "<<nv0<<endl; 
     Int_t ncas=event->GetNumberOfCascades();
     cout<<"Number of ESD cascades : "<<ncas<<endl; 

     //****** The loop over tracks
     Int_t nk=0;
     while (ntrk--) {
        AliESDtrack *track=event->GetTrack(ntrk);
        UInt_t status=track->GetStatus();

	//select only tracks with the "combined PID"
        if ((status&AliESDtrack::kESDpid)==0) continue;

        Double_t w[10]; track->GetESDpid(w);
        //count only "Kaon-like" tracks
        if (w[3]>w[4] && w[3]>w[2] && w[3]>w[1] && w[3]>w[0]) nk++;        
     }
     cout<<"Number of \"Kaon-like\" tracks : "<<nk<<endl;

     //****** The loop over V0s
     while (nv0--) {
        AliESDv0 *v0=event->GetV0(nv0);
        v0->ChangeMassHypothesis(310); // K0s
        Double_t mass=v0->GetEffMass();
        hm->Fill(mass);

        Int_t pidx=v0->GetPindex();               // now let's get an access  
        AliESDtrack *track=event->GetTrack(pidx); // to the positive daughter
        Double_t p=track->GetP();
        hp->Fill(p);
     }

     //****** The loop over cascades
     while (ncas--) {
        AliESDcascade *cas=event->GetCascade(ncas);
        Double_t q; //"quality" of the associated Lambda
        cas->ChangeMassHypothesis(q,3312); // Xi-
        // Here you do something with your Xis
        //  ...
        // You can get the access to the daughters
     }

   }

   delete event;
   ef->Close();

   timer.Stop(); timer.Print();

   TCanvas *c1=new TCanvas("c1","",0,0,600,1200);
   c1->Divide(1,2);

   c1->cd(1);
   hm->Fit("gaus","","",V0mass-V0width,V0mass+V0width);

   c1->cd(2);
   hp->Fit("expo","","",0.3,2); 

   return 0;
}
