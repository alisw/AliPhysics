/****************************************************************************
 *           Very important, delicate and rather obscure macro.             *
 *                    ("a la" AliTPCComparison.C)                           *
 *                                                                          *
 *               Creates list of "trackable" tracks,                        *
 *             calculates efficiency, resolutions etc.                      *
 *         (To get the list of the "trackable" tracks one should            *
 *             first run the AliTPCComparison.C macro)                      *
 *                                                                          *
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <TMath.h>
  #include <TError.h>
  #include <Riostream.h>
  #include <TH1F.h>
  #include <TH2F.h>
  #include <TTree.h>
  #include <TParticle.h>
  #include <TCanvas.h>
  #include <TLine.h>
  #include <TText.h>
  #include <TBenchmark.h>
  #include <TStyle.h>
  #include <TFile.h>
  #include <TROOT.h>

  #include "AliHeader.h"
  #include "AliTrackReference.h"
  #include "AliRunLoader.h"
  #include "AliRun.h"
  #include "AliESD.h"
#endif

Int_t GoodTracksTRD(const Char_t *dir=".");

extern AliRun *gAlice;
extern TBenchmark *gBenchmark;
extern TROOT *gROOT;

static Int_t allgood=0;
static Int_t allselected=0;
static Int_t allfound=0;

Int_t AliTRDComparisonV2
(Float_t ptcutl=0.2, Float_t ptcuth=10., const Char_t *dir=".") {
   gBenchmark->Start("AliTRDComparisonV2");

   ::Info("AliTRDComparisonV2.C","Doing comparison...");
   

   TH1F *hp=(TH1F*)gROOT->FindObject("hp");
   if (!hp) {
      hp=new TH1F("hp","PHI resolution",50,-70.,70.); 
      hp->SetFillColor(4);
      hp->SetXTitle("(mrad)"); 
   }
   TH1F *hl=(TH1F*)gROOT->FindObject("hl");
   if (!hl) {
      hl=new TH1F("hl","LAMBDA resolution",50,-70,70);
      hl->SetFillColor(4);
      hl->SetXTitle("(mrad)");
   }
   TH1F *hc=(TH1F*)gROOT->FindObject("hc");
   if (!hc) {
      hc=new TH1F("hc","Number of the assigned clusters",25,110,135);
      hc->SetLineColor(2);
   }
   TH1F *hpt=(TH1F*)gROOT->FindObject("hpt");
   if (!hpt) {
      hpt=new TH1F("hpt","Relative Pt resolution",30,-10.,10.); 
      hpt->SetFillColor(2);
      hpt->SetXTitle("(%)");
   }
   TH1F *hmpt=(TH1F*)gROOT->FindObject("hmpt");
   if (!hmpt) {
      hmpt=new TH1F("hmpt","Y and Z resolution",30,-30,30); 
      hmpt->SetFillColor(6);
      hmpt->SetXTitle("(mm)");
   }
   TH1F *hz=(TH1F*)gROOT->FindObject("hz");
   if (!hz) hz=new TH1F("hz","Z resolution",30,-30,30); 



   TH1F *hgood=(TH1F*)gROOT->FindObject("hgood");
   if (!hgood) hgood=new TH1F("hgood","Good tracks",30,0.2,6.1);
    
   TH1F *hfound=(TH1F*)gROOT->FindObject("hfound");
   if (!hfound) hfound=new TH1F("hfound","Found tracks",30,0.2,6.1);

   TH1F *hfake=(TH1F*)gROOT->FindObject("hfake");
   if (!hfake) hfake=new TH1F("hfake","Fake tracks",30,0.2,6.1);

   TH1F *hg=(TH1F*)gROOT->FindObject("hg");
   if (!hg) hg=new TH1F("hg","Efficiency for good tracks",30,0.2,6.1);
   hg->SetLineColor(4); hg->SetLineWidth(2);

   TH1F *hf=(TH1F*)gROOT->FindObject("hf");
   if (!hf) hf=new TH1F("hf","Efficiency for fake tracks",30,0.2,6.1);
   hf->SetFillColor(1); hf->SetFillStyle(3013); hf->SetLineWidth(2);

   TH1F *he=(TH1F*)gROOT->FindObject("he");
   if (!he) 
      he =new TH1F("he","dE/dX for pions with 0.4<p<0.5 GeV/c",50,0.,1000.);

   TH2F *hep=(TH2F*)gROOT->FindObject("hep");
   if (!hep) hep=new TH2F("hep","dE/dX vs momentum",50,0.,2.,50,0.,2000.);
   hep->SetMarkerStyle(8);
   hep->SetMarkerSize(0.4);


   Char_t fname[100];
   sprintf(fname,"%s/GoodTracksTRD.root",dir);

   TFile *refFile=TFile::Open(fname,"old");
   if (!refFile || !refFile->IsOpen()) {
   ::Info("AliTRDComparisonV2.C","Marking good tracks (will take a while)...");
     if (GoodTracksTRD(dir)) {
        ::Error("AliTRDComparisonV2.C","Can't generate the reference file !");
        return 1;
     }
   }
   refFile=TFile::Open(fname,"old");
   if (!refFile || !refFile->IsOpen()) {
     ::Error("AliTRDComparisonV2.C","Can't open the reference file !");
     return 1;
   }   
  
   TTree *trdTree=(TTree*)refFile->Get("trdTree");
   if (!trdTree) {
     ::Error("AliTRDComparisonV2.C","Can't get the reference tree !");
     return 2;
   }
   TBranch *branch=trdTree->GetBranch("TRD");
   if (!branch) {
     ::Error("AliTRDComparisonV2.C","Can't get the TRD branch !");
     return 3;
   }
   TClonesArray dummy("AliTrackReference",1000), *refs=&dummy;
   branch->SetAddress(&refs);


   sprintf(fname,"%s/AliESDs.root",dir);
   TFile *ef=TFile::Open(fname);
   if ((!ef)||(!ef->IsOpen())) {
      sprintf(fname,"%s/AliESDtrd.root",dir);
      ef=TFile::Open(fname);
      if ((!ef)||(!ef->IsOpen())) {
         ::Error("AliTRDComparisonV2.C","Can't open AliESDtrd.root !");
         return 4;
      }
   }
   AliESD* event = new AliESD;
   TTree* esdTree = (TTree*) ef->Get("esdTree");
   if (!esdTree) {
      ::Error("AliTRDComparisonV2.C", "no ESD tree found");
      return 6;
   }
   esdTree->SetBranchAddress("ESD", &event);


   //******* Loop over events *********
   Int_t e=0;
   while (esdTree->GetEvent(e)) {
     cout<<endl<<endl<<"********* Processing event number: "<<e<<"*******\n";
 
     Int_t nentr=event->GetNumberOfTracks();
     allfound+=nentr;

     if (trdTree->GetEvent(e++)==0) {
        cerr<<"No reconstructable tracks !\n";
        continue;
     }

     Int_t ngood=refs->GetEntriesFast(); 
     allgood+=ngood;

     const Int_t MAX=15000;
     Int_t notf[MAX], nnotf=0;
     Int_t fake[MAX], nfake=0;
     Int_t mult[MAX], numb[MAX], nmult=0;
     Int_t k;
     for (k=0; k<ngood; k++) {
	AliTrackReference *ref=(AliTrackReference*)refs->UncheckedAt(k); 
        Int_t lab=ref->Label(), tlab=-1;
        Float_t ptg=ref->Pt();

        if (ptg<ptcutl) continue;
        if (ptg>ptcuth) continue;

        allselected++;

        hgood->Fill(ptg);

        AliESDtrack *esd=0;
        Int_t cnt=0;
        for (Int_t i=0; i<nentr; i++) {
           AliESDtrack *t=event->GetTrack(i);
	   UInt_t status=t->GetStatus();

           if ((status&AliESDtrack::kTRDout)==0) continue;

           Int_t lbl=t->GetTRDLabel();
           if (lab==TMath::Abs(lbl)) {
	      if (cnt==0) {esd=t; tlab=lbl;}
              cnt++;
           }
        }
        if (cnt==0) {
           notf[nnotf++]=lab;
           continue;
        } else if (cnt>1){
           mult[nmult]=lab;
           numb[nmult]=cnt; nmult++;        
        }

        if (lab==tlab) hfound->Fill(ptg);
        else {
          fake[nfake++]=lab;
          hfake->Fill(ptg); 
        }

        AliExternalTrackParam out(*esd->GetOuterParam());

        Double_t bz=event->GetMagneticField();
        Double_t xg = ref->LocalX();
        out.PropagateTo(xg,bz);

        Float_t phi=TMath::ASin(out.GetSnp()) + out.GetAlpha();
        if (phi < 0) phi+=2*TMath::Pi();
        if (phi >=2*TMath::Pi()) phi-=2*TMath::Pi();
        Float_t phig=ref->Phi();
        hp->Fill((phi - phig)*1000.);

        Float_t lam=TMath::ATan(out.GetTgl()); 
        Float_t lamg=TMath::ATan2(ref->Pz(),ptg);
        hl->Fill((lam - lamg)*1000.);

        hc->Fill(esd->GetTRDclusters(0));

        Float_t pt_1=TMath::Abs(out.Get1Pt());
        hpt->Fill((pt_1 - 1/ptg)/(1/ptg)*100.);

        Float_t y=out.GetY();
        Float_t yg=ref->LocalY();
        hmpt->Fill(10*(y - yg));

        Float_t z=out.GetZ();
        Float_t zg=ref->Z();
        hz->Fill(10*(z - zg));

        Float_t mom=1./(pt_1*TMath::Cos(lam));
        Float_t dedx=esd->GetTRDsignal();
        hep->Fill(mom,dedx,1.);

        Int_t pdg=(Int_t)ref->GetLength();  //this is particle's PDG !

        if (TMath::Abs(pdg)==211) //pions
           if (mom>0.4 && mom<0.5) he->Fill(dedx,1.);

     }

     cout<<"\nList of Not found tracks :\n";
     for (k=0; k<nnotf; k++){
       cout<<notf[k]<<"\t";
       if ((k%9)==8) cout<<"\n";
     }
     cout<<"\n\nList of fake  tracks :\n";
     for (k=0; k<nfake; k++){
       cout<<fake[k]<<"\t";
       if ((k%9)==8) cout<<"\n";
     }
     cout<<"\n\nList of multiple found tracks :\n";
     for (k=0; k<nmult; k++) {
         cout<<"id.   "<<mult[k]
             <<"     found - "<<numb[k]<<"times\n";
     }
     cout<<endl;

     cout<<"Number of found tracks : "<<nentr<<endl;
     cout<<"Number of \"good\" tracks : "<<ngood<<endl;

     refs->Clear();
   } //***** End of the loop over events

   delete event;
   ef->Close();
   
   delete trdTree;
   refFile->Close();

   Stat_t ng=hgood->GetEntries(), nf=hfound->GetEntries();
   if (ng!=0) cout<<"\n\nIntegral efficiency is about "<<nf/ng*100.<<" %\n";
   cout<<"Total number selected of \"good\" tracks ="<<allselected<<endl<<endl;
   cout<<"Total number of found tracks ="<<allfound<<endl;
   cout<<"Total number of \"good\" tracks ="<<allgood<<endl;
   cout<<endl;

   gStyle->SetOptStat(111110);
   gStyle->SetOptFit(1);

   TCanvas *c1=new TCanvas("c1","",0,0,700,850);

   Int_t minc=33; 

   TPad *p1=new TPad("p1","",0,0.3,.5,.6); p1->Draw();
   p1->cd(); p1->SetFillColor(42); p1->SetFrameFillColor(10); 
   if (hp->GetEntries()<minc) hp->Draw(); else hp->Fit("gaus"); c1->cd();

   TPad *p2=new TPad("p2","",0.5,.3,1,.6); p2->Draw(); 
   p2->cd(); p2->SetFillColor(42); p2->SetFrameFillColor(10);
   if (hl->GetEntries()<minc) hl->Draw(); else hl->Fit("gaus"); c1->cd();

   TPad *p3=new TPad("p3","",0,0,0.5,0.3); p3->Draw();
   p3->cd(); p3->SetFillColor(42); p3->SetFrameFillColor(10); 
   if (hpt->GetEntries()<minc) hpt->Draw(); else hpt->Fit("gaus"); c1->cd();

   TPad *p4=new TPad("p4","",0.5,0,1,0.3); p4->Draw();
   p4->cd(); p4->SetFillColor(42); p4->SetFrameFillColor(10);
   if (hmpt->GetEntries()<minc) hmpt->Draw(); else hmpt->Fit("gaus"); 
   hz->Draw("same"); c1->cd();
   

   TPad *p5=new TPad("p5","",0,0.6,1,1); p5->Draw(); p5->cd(); 
   p5->SetFillColor(41); p5->SetFrameFillColor(10);
   hfound->Sumw2(); hgood->Sumw2(); hfake->Sumw2();
   hg->Divide(hfound,hgood,1,1.,"b");
   hf->Divide(hfake,hgood,1,1.,"b");
   hg->SetMaximum(1.4);
   hg->SetYTitle("Tracking efficiency");
   hg->SetXTitle("Pt (GeV/c)");
   hg->Draw();

   TLine *line1 = new TLine(0.2,1.0,6.1,1.0); line1->SetLineStyle(4);
   line1->Draw("same");
   TLine *line2 = new TLine(0.2,0.9,6.1,0.9); line2->SetLineStyle(4);
   line2->Draw("same");

   hf->SetFillColor(1);
   hf->SetFillStyle(3013);
   hf->SetLineColor(2);
   hf->SetLineWidth(2);
   hf->Draw("histsame");
   TText *text = new TText(0.461176,0.248448,"Fake tracks");
   text->SetTextSize(0.05);
   text->Draw();
   text = new TText(0.453919,1.11408,"Good tracks");
   text->SetTextSize(0.05);
   text->Draw();

   TCanvas *c2=new TCanvas("c2","",320,32,530,590);
   TPad *p6=new TPad("p6","",0.,0.,1.,.5); p6->Draw();
   p6->cd(); p6->SetFillColor(42); p6->SetFrameFillColor(10); 
   he->SetFillColor(2); he->SetFillStyle(3005);  
   he->SetXTitle("Arbitrary Units"); 
   if (he->GetEntries()<minc) he->Draw(); else he->Fit("gaus"); c2->cd();

   TPad *p7=new TPad("p7","",0.,0.5,1.,1.); p7->Draw(); 
   p7->cd(); p7->SetFillColor(42); p7->SetFrameFillColor(10);
   hep->SetXTitle("p (Gev/c)"); hep->SetYTitle("dE/dX (Arb. Units)"); 
   hep->Draw(); c1->cd();

   TFile fc("AliTRDComparisonV2.root","RECREATE");
   c1->Write();
   c2->Write();
   fc.Close();

   gBenchmark->Stop("AliTRDComparisonV2");
   gBenchmark->Show("AliTRDComparisonV2");

   return 0;
}



Int_t GoodTracksTRD(const Char_t *dir) {
   if (gAlice) { 
       delete gAlice->GetRunLoader();
       delete gAlice;//if everything was OK here it is already NULL
       gAlice = 0x0;
   }

   Char_t fname[100];
   sprintf(fname,"%s/galice.root",dir);

   AliRunLoader *rl = AliRunLoader::Open(fname,"COMPARISON");
   if (!rl) {
      ::Error("GoodTracksTRD","Can't start session !");
      return 1;
   }

   rl->LoadgAlice();
   rl->LoadHeader();
   rl->LoadTrackRefs();

   Int_t nev=rl->GetNumberOfEvents();
   ::Info("GoodTracksTRD","Number of events : %d\n",nev);  

   sprintf(fname,"%s/GoodTracksTPC.root",dir);
   TFile *tpcFile=TFile::Open(fname);
   if ((!tpcFile)||(!tpcFile->IsOpen())) {
       ::Error("GoodTracksTRD","Can't open the GoodTracksTPC.root !");
       delete rl;
       return 5; 
   }
   TClonesArray dum("AliTrackReference",1000), *tpcRefs=&dum;
   TTree *tpcTree=(TTree*)tpcFile->Get("tpcTree");
   if (!tpcTree) {
       ::Error("GoodTracksTRD","Can't get the TPC reference tree !");
       delete rl;
       return 6;
   }
   TBranch *tpcBranch=tpcTree->GetBranch("TPC");
   if (!tpcBranch) {
      ::Error("GoodTracksTRD","Can't get the TPC reference branch !");
      delete rl;
      return 7;
   }
   tpcBranch->SetAddress(&tpcRefs);

   sprintf(fname,"%s/GoodTracksTRD.root",dir);
   TFile *trdFile=TFile::Open(fname,"recreate");
   TClonesArray dummy2("AliTrackReference",1000), *trdRefs=&dummy2;
   TTree trdTree("trdTree","Info about the reconstructable TRD tracks");
   trdTree.Branch("TRD",&trdRefs);

   //********  Loop over generated events 
   for (Int_t e=0; e<nev; e++) {
     Int_t k;

     rl->GetEvent(e);  trdFile->cd();

     Int_t np = rl->GetHeader()->GetNtrack();
     cout<<"Event "<<e<<" Number of particles: "<<np<<endl;


     TTree *TR=rl->TreeTR();
     TBranch *branch=TR->GetBranch("TRD");
     if (branch==0) {
        ::Error("GoodTracksTRD","No TRD track references !");
        delete rl;
        return 5;
     }
     TClonesArray dummy("AliTrackReference",1000), *refs=&dummy;
     branch->SetAddress(&refs);
     Int_t nr=TR->GetEntries();

     //Preselect the "good" track candidates (TRD info only)
     Int_t nt=0;
     for (Int_t r=0; r<nr; r++) {
         refs->Clear();
         TR->GetEntry(r);
         Int_t n=refs->GetEntriesFast();
 
         if (n<=1) continue;

	 AliTrackReference *ref0=(AliTrackReference *)refs->UncheckedAt(0);
         if (ref0->LocalX() > 300.) continue;

	 AliTrackReference *refn=(AliTrackReference *)refs->UncheckedAt(n-1);
         if (refn->LocalX() < 363.) continue;

         if (TMath::Abs(ref0->Alpha() - refn->Alpha()) > 1e-5) continue;   

         new((*trdRefs)[nt++]) AliTrackReference(*refn);
     }

     //Check if the candidates are "good" from the TPC point of view
     tpcTree->GetEvent(e);
     Int_t ntpc=tpcRefs->GetEntriesFast();
     for (Int_t t=0; t<nt; t++) {
	 AliTrackReference *ref=(AliTrackReference *)trdRefs->UncheckedAt(t);
         Int_t lab=ref->Label();
         for (k=0; k<ntpc; k++) {
             AliTrackReference 
                *tpcRef=(AliTrackReference *)tpcRefs->UncheckedAt(k);
             if (tpcRef->Label()==lab) {
	        ref->SetLength(tpcRef->GetLength());
                break;
             } 
        }
	if (k==ntpc) delete trdRefs->RemoveAt(t); 
     }
     trdRefs->Compress();
     trdTree.Fill();

     trdRefs->Clear();
     tpcRefs->Clear();

   } //*** end of the loop over generated events

   trdTree.Write();
   trdFile->Close();

   delete tpcTree;
   tpcFile->Close();

   delete rl;
   return 0;
}


