/****************************************************************************
 *           Very important, delicate and rather obscure macro.             *
 *                                                                          *
 *               Creates list of "trackable" tracks,                        *
 *             calculates efficiency, resolutions etc.                      *
 *  The ESD tracks must be in an appropriate state (kITSin or kITSrefit)    *
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

  #include "AliStack.h"
  #include "AliHeader.h"
  #include "AliTrackReference.h"
  #include "AliRunLoader.h"
  #include "AliRun.h"
  #include "AliESD.h"

  #include "AliITSRecPoint.h"
  #include "AliITSLoader.h"
#endif

Int_t GoodTracksITS(const Char_t *dir=".");

extern AliRun *gAlice;
extern TBenchmark *gBenchmark;
extern TROOT *gROOT;

static Int_t allgood=0;
static Int_t allselected=0;
static Int_t allfound=0;

Int_t AliITSComparisonV2
(Float_t ptcutl=0.2, Float_t ptcuth=10., const Char_t *dir=".") {
   gBenchmark->Start("AliITSComparisonV2");

   ::Info("AliITSComparisonV2.C","Doing comparison...");
   

   TH1F *hp=(TH1F*)gROOT->FindObject("hp");
   if (!hp) hp=new TH1F("hp","PHI resolution",50,-20.,20.); 
   hp->SetFillColor(4);

   TH1F *hl=(TH1F*)gROOT->FindObject("hl");
   if (!hl) hl=new TH1F("hl","LAMBDA resolution",50,-20,20);
   hl->SetFillColor(4);

   TH1F *hpt=(TH1F*)gROOT->FindObject("hpt");
   if (!hpt) hpt=new TH1F("hpt","Relative Pt resolution",30,-10.,10.); 
   hpt->SetFillColor(2);
 
   TH1F *hmpt=(TH1F*)gROOT->FindObject("hmpt");
   if (!hmpt) 
      hmpt=new TH1F("hmpt","Transverse impact parameter",30,-777,777); 
   hmpt->SetFillColor(6);

   TH1F *hz=(TH1F*)gROOT->FindObject("hz");
   if (!hz) hz=new TH1F("hz","Longitudinal impact parameter",30,-777,777); 



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
      he =new TH1F("he","dE/dX for pions with 0.4<p<0.5 GeV/c",50,0.,100.);

   TH2F *hep=(TH2F*)gROOT->FindObject("hep");
   if (!hep) hep=new TH2F("hep","dE/dX vs momentum",50,0.,2.,50,0.,400.);
   hep->SetMarkerStyle(8);
   hep->SetMarkerSize(0.4);


   Char_t fname[100];
   sprintf(fname,"%s/GoodTracksITS.root",dir);

   TFile *refFile=TFile::Open(fname,"old");
   if (!refFile || !refFile->IsOpen()) {
   ::Info("AliITSComparisonV2.C","Marking good tracks (will take a while)...");
     if (GoodTracksITS(dir)) {
        ::Error("AliITSComparisonV2.C","Can't generate the reference file !");
        return 1;
     }
   }
   refFile=TFile::Open(fname,"old");
   if (!refFile || !refFile->IsOpen()) {
     ::Error("AliITSComparisonV2.C","Can't open the reference file !");
     return 1;
   }   
  
   TTree *itsTree=(TTree*)refFile->Get("itsTree");
   if (!itsTree) {
     ::Error("AliITSComparisonV2.C","Can't get the reference tree !");
     return 2;
   }
   TBranch *branch=itsTree->GetBranch("ITS");
   if (!branch) {
     ::Error("AliITSComparisonV2.C","Can't get the ITS branch !");
     return 3;
   }
   TClonesArray dummy("AliTrackReference",1000), *refs=&dummy;
   branch->SetAddress(&refs);


   sprintf(fname,"%s/AliESDs.root",dir);
   TFile *ef=TFile::Open(fname);
   if ((!ef)||(!ef->IsOpen())) {
      sprintf(fname,"%s/AliESDits.root",dir);
      ef=TFile::Open(fname);
      if ((!ef)||(!ef->IsOpen())) {
         ::Error("AliITSComparisonV2.C","Can't open AliESDits.root !");
         return 4;
      }
   }
   AliESD* event = new AliESD;
   TTree* esdTree = (TTree*) ef->Get("esdTree");
   if (!esdTree) {
      ::Error("AliITSComparison.C", "no ESD tree found");
      return 6;
   }
   esdTree->SetBranchAddress("ESD", &event);


   //******* Loop over events *********
   Int_t e=0;
   while (esdTree->GetEvent(e)) {
     cout<<endl<<endl<<"********* Processing event number: "<<e<<"*******\n";
 
     Int_t nentr=event->GetNumberOfTracks();
     allfound+=nentr;

     if (itsTree->GetEvent(e++)==0) {
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
        Float_t ptg=TMath::Sqrt(ref->Px()*ref->Px() + ref->Py()*ref->Py());

        if (ptg<ptcutl) continue;
        if (ptg>ptcuth) continue;

        allselected++;

        hgood->Fill(ptg);

        AliESDtrack *esd=0;
        Int_t cnt=0;
        for (Int_t i=0; i<nentr; i++) {
           AliESDtrack *t=event->GetTrack(i);
	   UInt_t status=t->GetStatus();

           if ((status&AliESDtrack::kITSrefit)==0) continue;

           Int_t lbl=t->GetLabel();
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

        Double_t alpha=esd->GetAlpha(),xv,par[5]; 
        esd->GetExternalParameters(xv,par);
        Float_t phi=TMath::ASin(par[2]) + alpha;
        if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
        if (phi>=TMath::Pi()) phi-=2*TMath::Pi();
        Float_t lam=TMath::ATan(par[3]); 
        Float_t pt_1=TMath::Abs(par[4]);

        Float_t phig=TMath::ATan2(ref->Py(),ref->Px());
        hp->Fill((phi - phig)*1000.);

        Float_t lamg=TMath::ATan2(ref->Pz(),ptg);
        hl->Fill((lam - lamg)*1000.);

        Float_t d,z; esd->GetImpactParameters(d,z);
        hmpt->Fill(10000*d);
        hz->Fill(10000*z);

        hpt->Fill((pt_1 - 1/ptg)/(1/ptg)*100.);

        Float_t mom=1./(pt_1*TMath::Cos(lam));
        Float_t dedx=esd->GetITSsignal();
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
   
   delete itsTree;
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
   hp->SetFillColor(4);  hp->SetXTitle("(mrad)"); 
   if (hp->GetEntries()<minc) hp->Draw(); else hp->Fit("gaus"); c1->cd();

   TPad *p2=new TPad("p2","",0.5,.3,1,.6); p2->Draw(); 
   p2->cd(); p2->SetFillColor(42); p2->SetFrameFillColor(10);
   hl->SetXTitle("(mrad)");
   if (hl->GetEntries()<minc) hl->Draw(); else hl->Fit("gaus"); c1->cd();

   TPad *p3=new TPad("p3","",0,0,0.5,0.3); p3->Draw();
   p3->cd(); p3->SetFillColor(42); p3->SetFrameFillColor(10); 
   hpt->SetXTitle("(%)");
   if (hpt->GetEntries()<minc) hpt->Draw(); else hpt->Fit("gaus"); c1->cd();

   TPad *p4=new TPad("p4","",0.5,0,1,0.3); p4->Draw();
   p4->cd(); p4->SetFillColor(42); p4->SetFrameFillColor(10);
   hmpt->SetXTitle("(micron)");
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

   TFile fc("AliITSComparisonV2.root","RECREATE");
   c1->Write();
   c2->Write();
   fc.Close();

   gBenchmark->Stop("AliITSComparisonV2");
   gBenchmark->Show("AliITSComparisonV2");

   return 0;
}



Int_t GoodTracksITS(const Char_t *dir) {
   if (gAlice) { 
       delete gAlice->GetRunLoader();
       delete gAlice;//if everything was OK here it is already NULL
       gAlice = 0x0;
   }

   Char_t fname[100];
   sprintf(fname,"%s/galice.root",dir);

   AliRunLoader *rl = AliRunLoader::Open(fname,"COMPARISON");
   if (!rl) {
      ::Error("GoodTracksITS","Can't start session !");
      return 1;
   }

   rl->LoadgAlice();
   rl->LoadHeader();
   rl->LoadKinematics();

   AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0) {
       ::Error("GoodTracksITS","Can not find the ITSLoader");
       delete rl;
       return 4;
   }
   itsl->LoadRecPoints();
  

   Int_t nev=rl->GetNumberOfEvents();
   ::Info("GoodTracksITS","Number of events : %d\n",nev);  

   sprintf(fname,"%s/GoodTracksTPC.root",dir);
   TFile *tpcFile=TFile::Open(fname);
   if ((!tpcFile)||(!tpcFile->IsOpen())) {
       ::Error("GoodTracksITS","Can't open the GoodTracksTPC.root !");
       delete rl;
       return 5; 
   }
   TClonesArray dum("AliTrackReference",1000), *tpcRefs=&dum;
   TTree *tpcTree=(TTree*)tpcFile->Get("tpcTree");
   if (!tpcTree) {
       ::Error("GoodTracksITS","Can't get the TPC reference tree !");
       delete rl;
       return 6;
   }
   TBranch *tpcBranch=tpcTree->GetBranch("TPC");
   if (!tpcBranch) {
      ::Error("GoodTracksITS","Can't get the TPC reference branch !");
      delete rl;
      return 7;
   }
   tpcBranch->SetAddress(&tpcRefs);

   sprintf(fname,"%s/GoodTracksITS.root",dir);
   TFile *itsFile=TFile::Open(fname,"recreate");
   TClonesArray dummy2("AliTrackReference",1000), *itsRefs=&dummy2;
   TTree itsTree("itsTree","Tree with info about the reconstructable ITS tracks");
   itsTree.Branch("ITS",&itsRefs);

   //********  Loop over generated events 
   for (Int_t e=0; e<nev; e++) {
     Int_t k;

     rl->GetEvent(e);  itsFile->cd();

     Int_t np = rl->GetHeader()->GetNtrack();
     cout<<"Event "<<e<<" Number of particles: "<<np<<endl;

     //******** Fill the "good" masks
     Int_t *good=new Int_t[np]; for (k=0; k<np; k++) good[k]=0;

     TTree *cTree=itsl->TreeR();
     if (!cTree) {
        ::Error("GoodTracksITS","Can't get the cluster tree !"); 
        delete rl;
        return 8;
     }
     TBranch *branch=cTree->GetBranch("ITSRecPoints");
     if (!branch) {
        ::Error("GoodTracksITS","Can't get the clusters branch !"); 
        delete rl;
        return 9;
     }
     TClonesArray dummy("AliITSRecPoint",10000), *clusters=&dummy;
     branch->SetAddress(&clusters);

     Int_t entr=(Int_t)cTree->GetEntries();
     for (k=0; k<entr; k++) {
         cTree->GetEvent(k);
         Int_t ncl=clusters->GetEntriesFast(); if (ncl==0) continue;
         while (ncl--) {
            AliITSRecPoint *pnt=(AliITSRecPoint*)clusters->UncheckedAt(ncl);

            Int_t lay=pnt->GetLayer();
            if (lay<0 || lay>5) {
               ::Error("GoodTracksITS","Wrong layer !");
               delete rl;
               return 10;
            }

            Int_t l0=pnt->GetLabel(0);
	       if (l0>=np) {
// 		 cerr<<"Wrong label: "<<l0<<endl;
		 continue;
	       }
            Int_t l1=pnt->GetLabel(1);
	       if (l1>=np) {
// 		 cerr<<"Wrong label: "<<l1<<endl;
		 continue;
	       }
            Int_t l2=pnt->GetLabel(2);
	       if (l2>=np) {
// 		 cerr<<"Wrong label: "<<l2<<endl;
		 continue;
	       }
            Int_t mask=1<<lay;
            if (l0>=0) good[l0]|=mask; 
            if (l1>=0) good[l1]|=mask; 
            if (l2>=0) good[l2]|=mask;
         }
         clusters->Clear();
     }
   

     //****** select tracks which are "good" enough
     AliStack* stack = rl->Stack();

     tpcTree->GetEvent(e);
     Int_t nk=tpcRefs->GetEntriesFast();
     Int_t nt=0;
     for (k=0; k<nk; k++) {
        AliTrackReference *tpcRef=(AliTrackReference *)tpcRefs->UncheckedAt(k);
        Int_t lab=tpcRef->Label();
        if (good[lab] != 0x3F) continue;
        TParticle *p = (TParticle*)stack->Particle(lab);
        if (p == 0x0) {
           cerr<<"Can not get particle "<<lab<<endl;
           continue;
        }

	AliTrackReference *ref=new((*itsRefs)[nt]) AliTrackReference(*tpcRef);
	ref->SetMomentum(p->Px(),p->Py(),p->Pz());
	ref->SetPosition(p->Vx(),p->Vy(),p->Vz());
        nt++;
     }
     tpcRefs->Clear();

     itsTree.Fill();
     itsRefs->Clear();

     delete[] good;

   } //*** end of the loop over generated events

   itsTree.Write();
   itsFile->Close();

   delete tpcTree;
   tpcFile->Close();

   delete rl;
   return 0;
}


