//******* Run TPCTestTracking.C(1) for fast simulator data 
//******* Run TPCTestTracking.C(2) for slow simulator data 
void TPCTestTracking(int fast=1) {
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gSystem->Load("libGeant3Dummy.so");      // a dummy version of Geant3
      gSystem->Load("PHOS/libPHOSdummy.so");   // the standard Alice classes 
      gSystem->Load("libgalice.so");           // the standard Alice classes 
   } else {
      delete gAlice;
      gAlice=0;
   }

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!file) file = new TFile("galice.root");

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }

   gAlice->GetEvent(0);

   TClonesArray *particles=gAlice->Particles(); 
   int np=particles->GetEntriesFast();

   AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");

   if (fast==1) {
      cerr<<"Making clusters...\n";
      TPC->Hits2Clusters();
      TClonesArray *clusters=TPC->Clusters();
      if (!clusters) {cerr<<"No clusters found !\n"; return;}
      int n=clusters->GetEntriesFast();
      cerr<<"Number of clusters "<<n<<"                                  \n";

      cerr<<"Marking \"good\" tracks...                                  \n";
      for (int i=0; i<n; i++) {
          AliTPCcluster *c=(AliTPCcluster*)clusters->UncheckedAt(i);
          int lab=c->fTracks[0]; if (lab<0) lab=-lab; lab--;
          if (lab<0) continue; //noise cluster
          int sector=c->fSector, row=c->fPadRow;
          GParticle *p=(GParticle*)particles->UncheckedAt(lab);
          int ks;
          if (row==51) {ks=p->GetKS()|0x100; p->SetKS(ks);}
          if (row==51-8) {ks=p->GetKS()|0x80; p->SetKS(ks);}
          ks=p->GetKS()+1; p->SetKS(ks);
      }

   } else {
      cerr<<"Looking for clusters...\n";
      TPC->Digits2Clusters();
      TClonesArray *clusters=TPC->Clusters();
      if (!clusters) {cerr<<"No clusters found !\n"; return;}
      int n=clusters->GetEntriesFast();
      cerr<<"Number of clusters "<<n<<"                                  \n";

      cerr<<"Marking \"good\" tracks...                                  \n";
      TClonesArray *digits=TPC->Digits();
      TTree *TD=gAlice->TreeD();
      TD->GetBranch("TPC")->SetAddress(&digits);
      int *count = new int[np];
      int i;
      for (i=0; i<np; i++) count[i]=0;
      int sectors_by_rows=(int)TD->GetEntries();
      for (i=0; i<sectors_by_rows; i++) {
          if (!TD->GetEvent(i)) continue;
          int sec, row;
          int ndigits=digits->GetEntriesFast();
          int j;
          for (j=0; j<ndigits; j++) {
              AliTPCdigit *dig = (AliTPCdigit*)digits->UncheckedAt(j);
              int idx=dig->fTracks[0]-1;
              if (idx<0) continue;
              sec=dig->fSector-1; row=dig->fPadRow-1;
              if (dig->fSignal>10) count[idx]+=dig->fSignal;
          }
          for (j=0; j<np; j++) {
              GParticle *p=(GParticle*)particles->UncheckedAt(j);
              if (count[j]>75) {
                 int ks;
                 if (row==51  ) {ks=p->GetKS()|0x100; p->SetKS(ks);}
                 if (row==51-8) {ks=p->GetKS()|0x80;  p->SetKS(ks);}
                 ks=p->GetKS()+1; p->SetKS(ks);
              }
              count[j]=0;
          }
      }
      delete[] count;

   }

   cerr<<"Looking for tracks...\n";
   TPC->Clusters2Tracks();
   TClonesArray *tracks=TPC->Tracks();
   if (!tracks) {cerr<<"No tracks found !\n"; return;}
   int nt=tracks->GetEntriesFast();
   cerr<<"Number of found tracks "<<nt<<endl;

/////////////////////////////////////////////////////////////////////////
   cerr<<"Doing comparison...\n";
   TH1F *hp=new TH1F("hp","PHI resolution",50,-100.,100.); hp->SetFillColor(4); 
   TH1F *hl=new TH1F("hl","LAMBDA resolution",50,-100,100); hl->SetFillColor(4);
   TH1F *hpt=new TH1F("hpt","Relative Pt resolution",50,-10.,10.); 
   hpt->SetFillColor(2); 
   TH1F *hd=new TH1F("hd","Impact parameter distribution ",50,0,25); 
   hd->SetFillColor(6);

   TH1F *hgood=new TH1F("hgood","Good tracks",20,0,2);    
   TH1F *hfound=new TH1F("hfound","Found tracks",20,0,2);
   TH1F *hfake=new TH1F("hfake","Fake tracks",20,0,2);
   TH1F *hg=new TH1F("hg","",20,0,2); //efficiency for good tracks
   hg->SetLineColor(4); hg->SetLineWidth(2);
   TH1F *hf=new TH1F("hf","Efficiency for fake tracks",20,0,2);
   hf->SetFillColor(1); hf->SetFillStyle(3013); hf->SetLineWidth(2);

   for (int i=0; i<np; i++) {
      GParticle *p = (GParticle*)particles->UncheckedAt(i);
      if (p->GetParent()>=0) continue;  //secondary particle
      if (p->GetKS()<0x100+0x80+32) continue; 
      Double_t ptg=p->GetPT(),pxg=p->GetPx(),pyg=p->GetPy(),pzg=p->GetPz();
      if (ptg<0.100) continue;
      if (fabs(pzg/ptg)>0.999) continue;
      hgood->Fill(ptg);
      int found=0;
      for (int j=0; j<nt; j++) {
          AliTPCtrack *track=(AliTPCtrack*)tracks->UncheckedAt(j);
          int lab=track->GetLab();
      if (fabs(lab)!=i) continue;
          found=1;
          Double_t xk=76.;

          track->PropagateTo(xk);
          xk-=0.11;
          track->PropagateTo(xk,42.7,2.27); //C
          xk-=2.6;
          track->PropagateTo(xk,36.2,1.98e-3); //C02
          xk-=0.051;
          track->PropagateTo(xk,42.7,2.27); //C

          xk-=0.4;
          track->PropagateTo(xk,21.82,2.33); //ITS+beam_pipe+etc (approximately)

          track->PropagateToVertex(); //comparison should be done at the vertex

          if (lab==i) hfound->Fill(ptg);
          else hfake->Fill(ptg);
          Double_t px,py,pz,pt=fabs(track->GetPt());track->GetPxPyPz(px,py,pz);
          Double_t phig=TMath::ATan(pyg/pxg);
          Double_t phi =TMath::ATan(py /px );
          hp->Fill((phi - phig)*1000.);
          Double_t lamg=TMath::ATan(pzg/ptg);
          Double_t lam =TMath::ATan(pz /pt );
          hl->Fill((lam - lamg)*1000.);
          hpt->Fill((pt - ptg)/ptg*100.);
          Double_t x,y,z; track->GetXYZ(x,y,z);
          hd->Fill(sqrt(x*x + y*y + z*z)*10.);
          break;
      }
      if (!found) cerr<<"Track number "<<i<<" was not found !\n";
   }
   Stat_t ngood =hgood->GetEntries(); cerr<<"Good tracks "<<ngood<<endl;
   Stat_t nfound=hfound->GetEntries();
   if (ngood!=0) 
      cerr<<"Integral efficiency is about "<<nfound/ngood*100.<<" %\n";

   gStyle->SetOptStat(111110);

   TCanvas *c1=new TCanvas("c1","",0,0,700,700);
   TPad *p1=new TPad("p1","",0,0.5,0.5,1); p1->Draw(); hp->SetXTitle("(mrad)");
   p1->cd(); hp->Draw(); c1->cd();
   TPad *p2=new TPad("p2","",0.5,0.5,1,1); p2->Draw(); hl->SetXTitle("(mrad)");
   p2->cd(); hl->Draw(); c1->cd();
   TPad *p3=new TPad("p3","",0,0,0.5,0.5); p3->Draw(); hpt->SetXTitle("(%)");
   p3->cd(); hpt->Draw(); c1->cd();
   TPad *p4=new TPad("p4","",0.5,0,1,0.5); p4->Draw(); hd->SetXTitle("(mm)");
   p4->cd(); hd->Draw(); c1->cd();

   TCanvas *c2=new TCanvas("c2",""); c2->cd();
   hfound->Sumw2(); hgood->Sumw2(); hfake->Sumw2();
   hg->Divide(hfound,hgood,1,1.,"b");
   hf->Divide(hfake,hgood,1,1.,"b");
   hg->SetMaximum(1.4);
   hg->SetYTitle("Tracking efficiency");
   hg->SetXTitle("Pt (GeV/c)");
   hg->Draw();

   TLine *line1 = new TLine(0,1.0,2,1.0); line1->SetLineStyle(4);
   line1->Draw("same");
   TLine *line2 = new TLine(0,0.9,2,0.9); line2->SetLineStyle(4);
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
}

