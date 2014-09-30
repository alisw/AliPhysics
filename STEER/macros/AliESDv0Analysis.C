//********************************************************************
//     Example (very naive for the moment) of the data analysis 
//                    using the ESD classes.
//       It demonstrates the idea of the "combined PID" 
//            applied to the Lambda0 reconstruction. 
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//********************************************************************

#if !defined( __CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include <TROOT.h>
  #include <TTree.h>
  #include <TFile.h>
  #include <TH1F.h>
  #include <TCanvas.h>

  #include "AliESDEvent.h"
  #include "AliESDv0.h"
#endif

extern TROOT *gROOT;

Int_t AliESDv0Analysis(const Char_t *dir=".") { 
   TH1F *hm=(TH1F*)gROOT->FindObject("hm");
   if (!hm) {
      hm=new TH1F("hm","Lambda+LambdaBar Effective Mass",60,1.065,1.165);
      hm->SetXTitle("Mass (GeV/c**2)");
   }
   Char_t fname[100];
   sprintf(fname,"%s/AliESDs.root",dir);
   TFile *ef=TFile::Open(fname);
   if (!ef||!ef->IsOpen()) {cerr<<"Can't AliESDs.root !\n"; return 1;}
   cerr<<"\n****** "<<fname<<" ******\n";

   AliESDEvent* event = new AliESDEvent();


   TTree* tree = (TTree*) ef->Get("esdTree");
   if (!tree) {cerr<<"no ESD tree found\n"; return 1;};
   event->ReadFromTree(tree);

   Int_t rc=0,n=0;

   //****** Tentative particle type "concentrations"
   Double_t c[5]={0.0, 0.0, 1, 0, 1};
   AliPID pid;
   pid.SetPriors(c);

   //******* The loop over events
    while (tree->GetEvent(n)) {

     cerr<<"Processing event number : "<<n++<<endl;

     Int_t nv0=event->GetNumberOfV0s();
     cerr<<"Number of ESD v0s : "<<nv0<<endl; 

     while (nv0--) {
       AliESDv0 *v0=event->GetV0(nv0);
       if (v0->GetOnFlyStatus()) continue;

       Int_t protonIdx=v0->GetPindex();
       Int_t pionIdx  =v0->GetNindex();
      
       v0->ChangeMassHypothesis(3122);
       Double_t mass=v0->GetEffMass();
       if (mass>1.17) {  //check also the LambdaBar hypothesis
          v0->ChangeMassHypothesis(-3122);
          mass=v0->GetEffMass();
          if (mass>1.17) continue;
          Int_t tmp=protonIdx; protonIdx=pionIdx; pionIdx=tmp;
       } 

       AliESDtrack *protonTrk=event->GetTrack(protonIdx);
       AliESDtrack *pionTrk  =event->GetTrack(pionIdx);

       if (protonTrk->GetP()<0.5) continue;

       // Check if the "proton track" is a proton
       if ((protonTrk->GetStatus()&AliESDtrack::kESDpid)!=0) {
	 Double_t r[10]; protonTrk->GetESDpid(r);
         pid.SetProbabilities(r);
         Double_t pp=pid.GetProbability(AliPID::kProton);
         if (pp < pid.GetProbability(AliPID::kElectron)) continue;
         if (pp < pid.GetProbability(AliPID::kMuon)) continue;
         if (pp < pid.GetProbability(AliPID::kPion)) continue;
         if (pp < pid.GetProbability(AliPID::kKaon)) continue;
       }
 
       //Check if the "pion track" is a pion
       if ((pionTrk->GetStatus()&AliESDtrack::kESDpid)!=0) {
	 Double_t r[10]; pionTrk->GetESDpid(r);
         pid.SetProbabilities(r);
         Double_t ppi=pid.GetProbability(AliPID::kPion);
         if (ppi < pid.GetProbability(AliPID::kElectron)) continue;
         if (ppi < pid.GetProbability(AliPID::kMuon)) continue;
         if (ppi < pid.GetProbability(AliPID::kKaon)) continue;
         if (ppi < pid.GetProbability(AliPID::kProton)) continue;
        }

       hm->Fill(mass);
     } 

   }

   delete event;
   delete tree;
   ef->Close();

   TCanvas *c1=(TCanvas*)gROOT->FindObject("c1");
   if (!c1) {
      c1=new TCanvas();
   }

   if (hm->GetEntries()>100) hm->Fit("gaus","","",1.11,1.12);
   else hm->Draw();

   c1->Update();

   return rc;
}
