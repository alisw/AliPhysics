#ifndef __CINT__
  #include <iostream.h>
  #include <fstream.h>

  #include "AliRun.h"
  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITStrackerV2.h"
  #include "AliITStrackV2.h"
  #include "AliITSclusterV2.h"

  #include "TFile.h"
  #include "TTree.h"
  #include "TH1.h"
  #include "TObjArray.h"
  #include "TStyle.h"
  #include "TCanvas.h"
  #include "TLine.h"
  #include "TText.h"
  #include "TParticle.h"
#endif

struct GoodTrack {
  Int_t lab;
  Int_t code;
  Float_t px,py,pz;
  Float_t x,y,z;
};
Int_t good_tracks(GoodTrack *gt, Int_t max);

Int_t AliITSComparisonV2() {
   cerr<<"Doing comparison...\n";

   TFile *cf=TFile::Open("AliITSclustersV2.root");
   if (!cf->IsOpen()) {cerr<<"Can't open AliITSclustersV2.root !\n"; return 1;}
   AliITSgeom *geom=(AliITSgeom*)cf->Get("AliITSgeom");
   if (!geom) { cerr<<"Can't get the ITS geometry !\n"; return 2; }
   AliITStrackerV2 tracker(geom);   

// Load tracks
   TFile *tf=TFile::Open("AliITStracksV2.root");
   if (!tf->IsOpen()) {cerr<<"Can't open AliITStracksV2.root !\n"; return 3;}
   TObjArray tarray(2000);
   TTree *tracktree=(TTree*)tf->Get("TreeT");
   TBranch *tbranch=tracktree->GetBranch("tracks");
   Int_t nentr=(Int_t)tracktree->GetEntries(),i;
   for (i=0; i<nentr; i++) {
       AliITStrackV2 *iotrack=new AliITStrackV2;
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);

       Int_t tpcLabel=iotrack->GetLabel();
       tracker.CookLabel(iotrack,0.);
       Int_t itsLabel=iotrack->GetLabel();
       if (itsLabel != tpcLabel) iotrack->SetLabel(-TMath::Abs(itsLabel));
       if (tpcLabel < 0)         iotrack->SetLabel(-TMath::Abs(itsLabel));
       /*
       if (itsLabel==1234) {
         Int_t nc=iotrack->GetNumberOfClusters();
         for (Int_t k=0; k<nc; k++) {
           Int_t index=iotrack->GetClusterIndex(k);
           AliITSclusterV2 *c=tracker.GetCluster(index);
           cout<<c->GetLabel(0)<<' '<<c->GetLabel(1)<<' '<<c->GetLabel(2)<<endl;
         }
       }
       */
       tarray.AddLast(iotrack);
   }   
   tf->Close();
   cf->Close();

/////////////////////////////////////////////////////////////////////////
   GoodTrack gt[15000];
   Int_t ngood=0;
   ifstream in("good_tracks_its");
   if (in) {
      cerr<<"Reading good tracks...\n";
      while (in>>gt[ngood].lab>>gt[ngood].code>>
                 gt[ngood].px>>gt[ngood].py>>gt[ngood].pz>>
                 gt[ngood].x >>gt[ngood].y >>gt[ngood].z) {
         ngood++;
         cerr<<ngood<<'\r';
         if (ngood==15000) {
            cerr<<"Too many good tracks !\n";
            break;
         }
      }
      if (!in.eof()) cerr<<"Read error (good_tracks_its) !\n";
   } else {
      cerr<<"Marking good tracks (this will take a while)...\n";
      ngood=good_tracks(gt,15000);
      ofstream out("good_tracks_its");
      if (out) {
         for (Int_t ngd=0; ngd<ngood; ngd++)            
	    out<<gt[ngd].lab<<' '<<gt[ngd].code<<' '<<
                 gt[ngd].px<<' '<<gt[ngd].py<<' '<<gt[ngd].pz<<' '<<
                 gt[ngd].x <<' '<<gt[ngd].y <<' '<<gt[ngd].z <<endl;
      } else cerr<<"Can not open file (good_tracks_its) !\n";
      out.close();
   }
   cerr<<"Number of good tracks : "<<ngood<<endl;

   TH1F *hp=new TH1F("hp","PHI resolution",50,-20.,20.); hp->SetFillColor(4);
   TH1F *hl=new TH1F("hl","LAMBDA resolution",50,-20,20);hl->SetFillColor(4);
   TH1F *hpt=new TH1F("hpt","Relative Pt resolution",30,-10.,10.); 
   hpt->SetFillColor(2); 
   TH1F *hmpt=new TH1F("hmpt","Transverse impact parameter",30,-300,300); 
   hmpt->SetFillColor(6);
   TH1F *hz=new TH1F("hz","Longitudinal impact parameter",30,-300,300); 
   //hmpt->SetFillColor(6);

   TH1F *hgood=new TH1F("hgood","Good tracks",30,0.1,6.1);    
   TH1F *hfound=new TH1F("hfound","Found tracks",30,0.1,6.1);
   TH1F *hfake=new TH1F("hfake","Fake tracks",30,0.1,6.1);
   TH1F *hg=new TH1F("hg","",30,0.1,6.1); //efficiency for good tracks
   hg->SetLineColor(4); hg->SetLineWidth(2);
   TH1F *hf=new TH1F("hf","Efficiency for fake tracks",30,0.1,6.1);
   hf->SetFillColor(1); hf->SetFillStyle(3013); hf->SetLineWidth(2);

   TH1F *hptw=new TH1F("hptw","Weghted pt",30,0.1,6.1);

   while (ngood--) {
      Int_t lab=gt[ngood].lab, tlab=-1;
      Double_t pxg=gt[ngood].px, pyg=gt[ngood].py, pzg=gt[ngood].pz;
      Double_t ptg=TMath::Sqrt(pxg*pxg+pyg*pyg);

      hgood->Fill(ptg);

      AliITStrackV2 *track=0;
      Int_t j;
      for (j=0; j<nentr; j++) {
          track=(AliITStrackV2*)tarray.UncheckedAt(j);
          tlab=track->GetLabel();
          if (lab==TMath::Abs(tlab)) break;
      }
      if (j==nentr) {
	cerr<<"Track "<<lab<<" was not found !\n";
        continue;
      }
      track->Propagate(track->GetAlpha(),3.,0.1/65.19*1.848,0.1*1.848);
      track->PropagateToVertex();

      if (lab==tlab) hfound->Fill(ptg);
      else { hfake->Fill(ptg); cerr<<lab<<" fake\n";}

      Double_t xv,par[5]; track->GetExternalParameters(xv,par);
      Float_t phi=TMath::ASin(par[2]) + track->GetAlpha();
      if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
      if (phi>=TMath::Pi()) phi-=2*TMath::Pi();
      Float_t lam=TMath::ATan(par[3]); 
      Float_t pt_1=TMath::Abs(par[4]);

      Double_t phig=TMath::ATan2(pyg,pxg);
      hp->Fill((phi - phig)*1000.);

      Double_t lamg=TMath::ATan2(pzg,ptg);
      hl->Fill((lam - lamg)*1000.);

      Double_t d=10000*track->GetD();
      hmpt->Fill(d);

      hptw->Fill(ptg,TMath::Abs(d));

      Double_t z=10000*track->GetZ();
      hz->Fill(z);

//if (TMath::Abs(gt[ngood].code)==11 && ptg>4.)
      hpt->Fill((pt_1 - 1/ptg)/(1/ptg)*100.);

   }

   Stat_t ng=hgood->GetEntries(); cerr<<"Good tracks "<<ng<<endl;
   Stat_t nf=hfound->GetEntries();
   if (ng!=0) 
      cerr<<"Integral efficiency is about "<<nf/ng*100.<<" %\n";

   gStyle->SetOptStat(111110);
   gStyle->SetOptFit(1);

   TCanvas *c1=new TCanvas("c1","",0,0,700,850);

   TPad *p1=new TPad("p1","",0,0.3,.5,.6); p1->Draw();
   p1->cd(); p1->SetFillColor(42); p1->SetFrameFillColor(10); 
   hp->SetFillColor(4);  hp->SetXTitle("(mrad)"); hp->Fit("gaus"); c1->cd();

   TPad *p2=new TPad("p2","",0.5,.3,1,.6); p2->Draw(); 
   p2->cd(); p2->SetFillColor(42); p2->SetFrameFillColor(10);
   hl->SetXTitle("(mrad)"); hl->Fit("gaus"); c1->cd();

   TPad *p3=new TPad("p3","",0,0,0.5,0.3); p3->Draw();
   p3->cd(); p3->SetFillColor(42); p3->SetFrameFillColor(10); 
   hpt->SetXTitle("(%)"); hpt->Fit("gaus"); c1->cd();

   TPad *p4=new TPad("p4","",0.5,0,1,0.3); p4->Draw();
   p4->cd(); p4->SetFillColor(42); p4->SetFrameFillColor(10);
   hmpt->SetXTitle("(micron)"); hmpt->Fit("gaus"); hz->Draw("same"); c1->cd();
   //hfound->Sumw2();
   //hptw->Sumw2(); 
   //hg->SetMaximum(333);
   //hg->SetYTitle("Impact Parameter Resolution (micron)");
   //hg->SetXTitle("Pt (GeV/c)");
   //hg->GetXaxis()->SetRange(0,10);
   //hg->Divide(hptw,hfound,1,1.);
   //hg->DrawCopy(); c1->cd();
   

   TPad *p5=new TPad("p5","",0,0.6,1,1); p5->Draw(); p5->cd(); 
   p5->SetFillColor(41); p5->SetFrameFillColor(10);
   hfound->Sumw2(); hgood->Sumw2(); hfake->Sumw2();
   hg->Divide(hfound,hgood,1,1.,"b");
   hf->Divide(hfake,hgood,1,1.,"b");
   hg->SetMaximum(1.4);
   hg->SetYTitle("Tracking efficiency");
   hg->SetXTitle("Pt (GeV/c)");
   hg->Draw();

   TLine *line1 = new TLine(0,1.0,7,1.0); line1->SetLineStyle(4);
   line1->Draw("same");
   TLine *line2 = new TLine(0,0.9,7,0.9); line2->SetLineStyle(4);
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

   return 0;
}

Int_t good_tracks(GoodTrack *gt, Int_t max) {
   if (gAlice) {delete gAlice; gAlice=0;}

   TFile *file=TFile::Open("galice.root");
   if (!file->IsOpen()) {cerr<<"Can't open galice.root !\n"; exit(4);}
   if (!(gAlice=(AliRun*)file->Get("gAlice"))) {
     cerr<<"gAlice have not been found on galice.root !\n";
     exit(5);
   }

   Int_t np=gAlice->GetEvent(0);

   Int_t *good=new Int_t[np];
   Int_t k;
   for (k=0; k<np; k++) good[k]=0;

   AliITS *ITS=(AliITS*)gAlice->GetDetector("ITS");
   if (!ITS) {
      cerr<<"can't get ITS !\n"; exit(8);
   }
   AliITSgeom *geom=ITS->GetITSgeom();
   if (!geom) {
      cerr<<"cen't get ITS geometry !\n"; exit(9);
   }

   TFile *cf=TFile::Open("AliITSclustersV2.root");
   if (!cf->IsOpen()){
      cerr<<"Can't open AliITSclustersV2.root !\n"; exit(6);
   }
   TTree *cTree=(TTree*)cf->Get("cTree");
   if (!cTree) {
      cerr<<"Can't get cTree !\n"; exit(7);
   }
   TBranch *branch=cTree->GetBranch("Clusters");
   if (!branch) {
      cerr<<"Can't get clusters branch !\n"; exit(8);
   }
   TClonesArray *clusters=new TClonesArray("AliITSclusterV2",10000);
   branch->SetAddress(&clusters);

   Int_t entr=(Int_t)cTree->GetEntries();
   for (k=0; k<entr; k++) {
     if (!cTree->GetEvent(k)) continue;
     Int_t lay,lad,det;  geom->GetModuleId(k-1,lay,lad,det);
     if (lay<1 || lay>6) {
	cerr<<"wrong layer !\n"; exit(10);
     }
     Int_t ncl=clusters->GetEntriesFast();
     while (ncl--) {
        AliITSclusterV2 *pnt=(AliITSclusterV2*)clusters->UncheckedAt(ncl);
        Int_t l0=pnt->GetLabel(0);
        Int_t l1=pnt->GetLabel(1);
        Int_t l2=pnt->GetLabel(2);
        Int_t mask=1<<(lay-1);
        if (l0>=0) good[l0]|=mask; 
        if (l1>=0) good[l1]|=mask; 
        if (l2>=0) good[l2]|=mask;
     }
   }
   clusters->Delete(); delete clusters;
   cf->Close();

   ifstream in("good_tracks_tpc");
   if (!in) {
     cerr<<"can't get good_tracks_tpc !\n"; exit(11);
   }
   Int_t nt=0;
   Double_t px,py,pz,x,y,z;
   Int_t code,lab;
   while (in>>lab>>code>>px>>py>>pz>>x>>y>>z) {
      if (good[lab] != 0x3F) continue;
      TParticle *p = (TParticle*)gAlice->Particle(lab);
      gt[nt].lab=lab;
      gt[nt].code=p->GetPdgCode();
//**** px py pz - in global coordinate system
      gt[nt].px=p->Px(); gt[nt].py=p->Py(); gt[nt].pz=p->Pz();
      gt[nt].x=gt[nt].y=gt[nt].z=0.;
      nt++;
      if (nt==max) {cerr<<"Too many good tracks !\n"; break;}
   }

   delete[] good;

   delete gAlice; gAlice=0;
   file->Close();

   return nt;
}


