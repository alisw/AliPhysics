void AliTPCTestTracking() {
   const char *pname="Param1";
   const char *tname="TreeD0_Param1";

// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
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
   int ver=TPC->IsVersion();
   cerr<<"TPC version "<<ver<<" has been found !\n";

   AliTPCD *digp= (AliTPCD*)file->Get(pname);
   if (digp!=0) TPC->SetDigParam(digp);
   else cerr<<"Warning: default TPC parameters will be used !\n";

   int nrows=TPC->GetDigParam()->GetParam().GetNRowUp()-1;
   switch (ver) {
   case 1:
      cerr<<"Making clusters...\n";
      TPC->Hits2Clusters();
      TClonesArray *clusters=TPC->Clusters();
      if (!clusters) {cerr<<"No clusters found !\n"; return;}
      int n=clusters->GetEntriesFast();
      cerr<<"Number of clusters "<<n<<"                                  \n";
      
      cerr<<"Marking \"good\" tracks...                                  \n";
      for (int i=0; i<n; i++) {
          AliTPCcluster *c=(AliTPCcluster*)clusters->UncheckedAt(i);
          int lab=c->fTracks[0];
          if (lab<0) continue; //noise cluster
          lab=TMath::Abs(lab);
          int sector=c->fSector-1, row=c->fPadRow-1;
          GParticle *p=(GParticle*)particles->UncheckedAt(lab);
          int ks;
          if (row==nrows)   {ks=p->GetKS()|0x1000; p->SetKS(ks);}
          if (row==nrows-8) {ks=p->GetKS()|0x800; p->SetKS(ks);}
          ks=p->GetKS()+1; p->SetKS(ks);
      }
      
      break;
   case 2:
      cerr<<"Looking for clusters...\n";
      TPC->Digits2Clusters();
      TClonesArray *clusters=TPC->Clusters();
      if (!clusters) {cerr<<"No clusters found !\n"; return;}
      int n=clusters->GetEntriesFast();
      cerr<<"Number of clusters "<<n<<"                                  \n";
      
      cerr<<"Marking \"good\" tracks...                                  \n";
      TTree *TD=gDirectory->Get(tname);
      TClonesArray *digits=TPC->Digits();
      TD->GetBranch("Digits")->SetAddress(&digits);

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
              int idx0=dig->fTracks[0];
              int idx1=dig->fTracks[1];
              int idx2=dig->fTracks[2];
              sec=dig->fSector-1; row=dig->fPadRow-1;
              if (idx0>=0 && dig->fSignal>0) count[idx0]+=1;
              if (idx1>=0 && dig->fSignal>0) count[idx1]+=1;
              if (idx2>=0 && dig->fSignal>0) count[idx2]+=1;
          }
          for (j=0; j<np; j++) {
              GParticle *p=(GParticle*)particles->UncheckedAt(j);
              if (count[j]>1) {
                 int ks;
                 if (row==nrows   ) {ks=p->GetKS()|0x1000; p->SetKS(ks);}
                 if (row==nrows-8 ) {ks=p->GetKS()|0x800;  p->SetKS(ks);}
                 ks=p->GetKS()+1; p->SetKS(ks);
              }
              count[j]=0;
          }
      }
      delete[] count;
      
      break;
   default:
      cerr<<"Invalid TPC version !\n";
      return;
   }

   cerr<<"Looking for tracks...\n";
   TPC->Clusters2Tracks();
   int nt=0;
   TClonesArray *tracks=TPC->Tracks();
   if (tracks) nt=tracks->GetEntriesFast();
   cerr<<"Number of found tracks "<<nt<<endl;

/////////////////////////////////////////////////////////////////////////
   cerr<<"Doing comparison...\n";
   TH1F *hp=new TH1F("hp","PHI resolution",50,-100.,100.); hp->SetFillColor(4); 
   TH1F *hl=new TH1F("hl","LAMBDA resolution",50,-100,100); hl->SetFillColor(4);
   TH1F *hpt=new TH1F("hpt","Relative Pt resolution",30,-10.,10.); 
   hpt->SetFillColor(2); 
   TH1F *hd=new TH1F("hd","Impact parameter distribution ",30,0,25); 
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
      if (p->GetKS()<0x1000+0x800+2+30) continue;
      Double_t ptg=p->GetPT(),pxg=p->GetPx(),pyg=p->GetPy(),pzg=p->GetPz();
      if (ptg<0.100) continue;
      if (fabs(pzg/ptg)>0.999) continue;

      //cout<<i<<endl;

      hgood->Fill(ptg);
      int found=0;
      for (int j=0; j<nt; j++) {
          AliTPCtrack *track=(AliTPCtrack*)tracks->UncheckedAt(j);
          int lab=track->GetLab();
	  if (fabs(lab)!=i) continue;
	  //if (lab!=i) continue;
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
          else { hfake->Fill(ptg); cerr<<lab<<" fake\n";}
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
   hd->SetXTitle("(mm)"); hd->Draw(); c1->cd();

   TPad *p5=new TPad("p5","",0,0.6,1,1); p5->Draw(); p5->cd(); 
   /*p5->SetTopMargin(0.25);*/ p5->SetFillColor(41); p5->SetFrameFillColor(10);
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

