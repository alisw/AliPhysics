//#include "alles.h"
struct GoodTrack {
  Int_t lab;
  Int_t code;
  Float_t px,py,pz;
  Float_t x,y,z;
};
Int_t good_tracks(GoodTrack *gt, Int_t max);

Int_t AliTPCComparison() {
  Int_t i;
   gBenchmark->Start("AliTPCComparison");
   cerr<<"Doing comparison...\n";

   TFile *cf=TFile::Open("AliTPCclusters.root");
   if (!cf->IsOpen()) {cerr<<"Can't open AliTPCclusters.root !\n"; return 1;}
   AliTPCParam *digp= (AliTPCParam*)cf->Get("75x40_100x60");
   if (!digp) { cerr<<"TPC parameters have not been found !\n"; return 2; }

// Load clusters
   AliTPCClustersArray *ca=new AliTPCClustersArray;
   ca->Setup(digp);
   ca->SetClusterType("AliTPCcluster");
   ca->ConnectTree("Segment Tree");
   Int_t nentr=Int_t(ca->GetTree()->GetEntries());
   for (i=0; i<nentr; i++) {
     //AliSegmentID *s=;
       ca->LoadEntry(i);
   }

// Load tracks
   TFile *tf=TFile::Open("AliTPCtracks.root");
   if (!tf->IsOpen()) {cerr<<"Can't open AliTPCtracks.root !\n"; return 3;}
   TObjArray tarray(2000);
   TTree *tracktree=(TTree*)tf->Get("TreeT");
   TBranch *tbranch=tracktree->GetBranch("tracks");
   nentr=(Int_t)tracktree->GetEntries();
   
   for (i=0; i<nentr; i++) {
       AliTPCtrack *iotrack=new AliTPCtrack;
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);
       iotrack->CookLabel(ca);
       tarray.AddLast(iotrack);
   }   
   tf->Close();

   delete ca;
   cf->Close();

/////////////////////////////////////////////////////////////////////////
   GoodTrack gt[20000];
   Int_t ngood=0;
   ifstream in("good_tracks");
   if (in) {
      cerr<<"Reading good tracks...\n";
      while (in>>gt[ngood].lab>>gt[ngood].code
	       >>gt[ngood].px >>gt[ngood].py>>gt[ngood].pz
	       >>gt[ngood].x  >>gt[ngood].y >>gt[ngood].z) {
         ngood++;
         cerr<<ngood<<'\r';
         if (ngood==20000) {
            cerr<<"Too many good tracks !\n";
            break;
         }
      }
      if (!in.eof()) cerr<<"Read error (good_tracks) !\n";
   } else {
      cerr<<"Marking good tracks (this will take a while)...\n";
      ngood=good_tracks(gt,20000);
      ofstream out("good_tracks");
      if (out) {
         for (Int_t ngd=0; ngd<ngood; ngd++)            
	    out<<gt[ngd].lab<<' '<<gt[ngd].code<<' '
	       <<gt[ngd].px <<' '<<gt[ngd].py<<' '<<gt[ngd].pz<<' '
	       <<gt[ngd].x  <<' '<<gt[ngd].y <<' '<<gt[ngd].z <<endl;
      } else cerr<<"Can not open file (good_tracks) !\n";
      out.close();
   }
   cerr<<"Number of good tracks : "<<ngood<<endl;

   cerr<<"Doing comparison...\n";
   TH1F *hp=new TH1F("hp","PHI resolution",50,-20.,20.); hp->SetFillColor(4);
   TH1F *hl=new TH1F("hl","LAMBDA resolution",50,-20,20);hl->SetFillColor(4);
   TH1F *hpt=new TH1F("hpt","Relative Pt resolution",30,-10.,10.); 
   hpt->SetFillColor(2); 
   TH1F *hmpt=new TH1F("hmpt","Relative Pt resolution (pt>4GeV/c)",30,-60,60); 
   hmpt->SetFillColor(6);

   TH1F *hgood=new TH1F("hgood","Good tracks",20,0,7);    
   TH1F *hfound=new TH1F("hfound","Found tracks",20,0,7);
   TH1F *hfake=new TH1F("hfake","Fake tracks",20,0,7);
   TH1F *hg=new TH1F("hg","",20,0,7); //efficiency for good tracks
   hg->SetLineColor(4); hg->SetLineWidth(2);
   TH1F *hf=new TH1F("hf","Efficiency for fake tracks",20,0,7);
   hf->SetFillColor(1); hf->SetFillStyle(3013); hf->SetLineWidth(2);

   TH1F *he =new TH1F("he","dE/dX for pions with 0.4<p<0.5 GeV/c",50,0.,100.);
   TH2F *hep=new TH2F("hep","dE/dX vs momentum",50,0.,2.,50,0.,200.);

   Int_t mingood = ngood;  //MI change
   Int_t * track_notfound = new Int_t[ngood];
   Int_t   itrack_notfound =0;
   Int_t * track_fake  = new Int_t[ngood];
   Int_t   itrack_fake = 0;
   Int_t * track_multifound = new Int_t[ngood];
   Int_t * track_multifound_n = new Int_t[ngood];
   Int_t   itrack_multifound =0;

   while (ngood--) {
      Int_t lab=gt[ngood].lab,tlab=-1;
      Float_t ptg=
      TMath::Sqrt(gt[ngood].px*gt[ngood].px + gt[ngood].py*gt[ngood].py);

      hgood->Fill(ptg);

      AliTPCtrack *track=0;
      for (i=0; i<nentr; i++) {
          track=(AliTPCtrack*)tarray.UncheckedAt(i);
          tlab=track->GetLabel();
          if (lab==TMath::Abs(tlab)) break;
      }
      if (i==nentr) {
	//	cerr<<"Track "<<lab<<" was not found !\n";
	track_notfound[itrack_notfound++]=lab;
        continue;
      }
      
      //MI change  - addition
      Int_t micount=0;
      Int_t mi;
      AliTPCtrack * mitrack;
      for (mi=0; mi<nentr; mi++) {
	mitrack=(AliTPCtrack*)tarray.UncheckedAt(mi);          
	if (lab==TMath::Abs(mitrack->GetLabel())) micount++;
      }
      if (micount>1) {
	//cout<<"Track no. "<<lab<<" found "<<micount<<"  times\n";
	track_multifound[itrack_multifound]=lab;
	track_multifound_n[itrack_multifound]=micount;
	itrack_multifound++;
      }
      

      //
      Double_t xk=gt[ngood].x;//digp->GetPadRowRadii(0,0);
      track->PropagateTo(xk);

      if (lab==tlab) hfound->Fill(ptg);
      else { 
	track_fake[itrack_fake++]=lab;
	hfake->Fill(ptg); 
	//cerr<<lab<<" fake\n";
      }

      Double_t px,py,pz,pt=fabs(track->GetPt());track->GetPxPyPz(px,py,pz);

      if (TMath::Abs(gt[ngood].code)==11 && ptg>4.) {
         hmpt->Fill((1/pt - 1/ptg)/(1/ptg)*100.);
      } else {
         Float_t phig=TMath::ATan2(gt[ngood].py,gt[ngood].px);
         Float_t phi =TMath::ATan2(py,px );
         hp->Fill((phi - phig)*1000.);

         Float_t lamg=TMath::ATan2(gt[ngood].pz,ptg);
         Float_t lam =TMath::ATan2(pz ,pt );
         hl->Fill((lam - lamg)*1000.);

         hpt->Fill((1/pt - 1/ptg)/(1/ptg)*100.);
      }

      Float_t mom=track->GetP();
      Float_t dedx=track->GetdEdx();
      hep->Fill(mom,dedx,1.);
      if (TMath::Abs(gt[ngood].code)==211)
	 if (mom>0.4 && mom<0.5) {
            he->Fill(dedx,1.);
         }

   }

   
   Stat_t ng=hgood->GetEntries(); cerr<<"Good tracks "<<ng<<endl;
   Stat_t nf=hfound->GetEntries();
   if (ng!=0) 
      cerr<<"\n\n\nIntegral efficiency is about "<<nf/ng*100.<<" %\n";
   
   //MI change  - addition
   cout<<"Total number of found tracks ="<<nentr<<"\n";
   cout<<"Total number of \"good\" tracks ="<<mingood<<"\n";


   cout<<"\nList of Not found tracks :\n";
   //   Int_t i;
   for ( i = 0; i< itrack_notfound; i++){
     cout<<track_notfound[i]<<"\t";
     if ((i%5)==4) cout<<"\n";
   }
   cout<<"\nList of fake  tracks :\n";
   for ( i = 0; i< itrack_fake; i++){
     cout<<track_fake[i]<<"\t";
     if ((i%5)==4) cout<<"\n";
   }
    cout<<"\nList of multiple found tracks :\n";
   for ( i=0; i<itrack_multifound; i++) {
       cout<<"id.   "<<track_multifound[i]<<"     found - "<<track_multifound_n[i]<<"times\n";
   }
   cout<<"\n\n\n";

   
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
   hmpt->SetXTitle("(%)"); hmpt->Fit("gaus"); c1->cd();

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



   TCanvas *c2=new TCanvas("c2","",320,32,530,590);

   TPad *p6=new TPad("p6","",0.,0.,1.,.5); p6->Draw();
   p6->cd(); p6->SetFillColor(42); p6->SetFrameFillColor(10); 
   he->SetFillColor(2); he->SetFillStyle(3005);  
   he->SetXTitle("Arbitrary Units"); 
   he->Fit("gaus"); c2->cd();

   TPad *p7=new TPad("p7","",0.,0.5,1.,1.); p7->Draw(); 
   p7->cd(); p7->SetFillColor(42); p7->SetFrameFillColor(10);
   hep->SetXTitle("p (Gev/c)"); hep->SetYTitle("dE/dX (Arb. Units)"); 
   hep->Draw(); c1->cd();

   gBenchmark->Stop("AliTPCComparison");
   gBenchmark->Show("AliTPCComparison");


   return 0;
}


Int_t good_tracks(GoodTrack *gt, Int_t max) {
   Int_t nt=0;

   TFile *file=TFile::Open("galice.root");
   if (!file->IsOpen()) {cerr<<"Can't open galice.root !\n"; exit(4);}

   if (!(gAlice=(AliRun*)file->Get("gAlice"))) {
     cerr<<"gAlice have not been found on galice.root !\n";
     exit(5);
   }

   gAlice->GetEvent(0);   

   AliTPC *TPC=(AliTPC*)gAlice->GetDetector("TPC");
   Int_t ver = TPC->IsVersion(); 
   cerr<<"TPC version "<<ver<<" has been found !\n";

   AliTPCParam *digp=(AliTPCParam*)file->Get("75x40_100x60");
   if (!digp) { cerr<<"TPC parameters have not been found !\n"; exit(6); }
   TPC->SetParam(digp);

   Int_t nrow_up=digp->GetNRowUp();
   Int_t nrows=digp->GetNRowLow()+nrow_up;
   Int_t zero=digp->GetZeroSup();
   Int_t gap=Int_t(0.125*nrows);
   Int_t good_number=Int_t(0.4*nrows);

   TClonesArray *particles=gAlice->Particles(); 
   Int_t np=particles->GetEntriesFast();
   Int_t *good=new Int_t[np];
   for (Int_t ii=0; ii<np; ii++) good[ii]=0;


   //MI change to be possible compile macro
   //definition out of the swith statemnet
    Int_t sectors_by_rows=0;
    TTree *TD=0;
    AliSimDigits da, *digits=&da;
    Int_t *count=0;
   switch (ver) {
   case 1:
     {
      TFile *cf=TFile::Open("AliTPCclusters.root");
      if (!cf->IsOpen()){cerr<<"Can't open AliTPCclusters.root !\n";exit(5);}
      AliTPCClustersArray *ca=new AliTPCClustersArray;
      ca->Setup(digp);
      ca->SetClusterType("AliTPCcluster");
      ca->ConnectTree("Segment Tree");
      Int_t nrows=Int_t(ca->GetTree()->GetEntries());
      for (Int_t n=0; n<nrows; n++) {
          AliSegmentID *s=ca->LoadEntry(n);
          Int_t sec,row;
          digp->AdjustSectorRow(s->GetID(),sec,row);
          AliTPCClustersRow &clrow = *ca->GetRow(sec,row);
          Int_t ncl=clrow.GetArray()->GetEntriesFast();
          while (ncl--) {
              AliTPCcluster *c=(AliTPCcluster*)clrow[ncl];
              Int_t lab=c->GetLabel(0);
              if (lab<0) continue; //noise cluster
              lab=TMath::Abs(lab);
              if (sec>=digp->GetNInnerSector())
              if (row==nrow_up-1    ) good[lab]|=0x1000;
              if (sec>=digp->GetNInnerSector())
              if (row==nrow_up-1-gap) good[lab]|=0x800;
              good[lab]++;
          }
          ca->ClearRow(sec,row);
      }
      cf->Close();
     }
      break;
   case 2:
      TD=(TTree*)gDirectory->Get("TreeD_75x40_100x60");
      TD->GetBranch("Segment")->SetAddress(&digits);
      count = new Int_t[np];
      Int_t i;
      for (i=0; i<np; i++) count[i]=0;
      sectors_by_rows=(Int_t)TD->GetEntries();
      for (i=0; i<sectors_by_rows; i++) {
          if (!TD->GetEvent(i)) continue;
          Int_t sec,row;
          digp->AdjustSectorRow(digits->GetID(),sec,row);
          cerr<<sec<<' '<<row<<"                                     \r";
          digits->First();
          while (digits->Next()) {
              Int_t it=digits->CurrentRow(), ip=digits->CurrentColumn();
              Short_t dig = digits->GetDigit(it,ip);
              Int_t idx0=digits->GetTrackID(it,ip,0); 
              Int_t idx1=digits->GetTrackID(it,ip,1);
              Int_t idx2=digits->GetTrackID(it,ip,2);
              if (idx0>=0 && dig>=zero) count[idx0]+=1;
              if (idx1>=0 && dig>=zero) count[idx1]+=1;
              if (idx2>=0 && dig>=zero) count[idx2]+=1;
          }
          for (Int_t j=0; j<np; j++) {
              if (count[j]>1) {
                 if (sec>=digp->GetNInnerSector())
		   if (row==nrow_up-1    ) good[j]|=0x1000;
                 if (sec>=digp->GetNInnerSector())
		   if (row==nrow_up-1-gap) good[j]|=0x800;
                 good[j]++;
              }
              count[j]=0;
          }
      }
      delete[] count;
      break;
   default:
      cerr<<"Invalid TPC version !\n";
      file->Close();
      exit(7);
   }

   TTree *TH=gAlice->TreeH();
   //   TClonesArray *hits=TPC->Hits();
   Int_t npart=(Int_t)TH->GetEntries();

   while (npart--) {
      AliTPChit *hit0=0;

      TPC->ResetHits();
      TH->GetEvent(npart);
      AliTPChit * hit = (AliTPChit*) TPC->FirstHit(-1);
      while (hit){
	if (hit->fQ==0.) break;
	hit =  (AliTPChit*) TPC->NextHit();
      }
      if (hit) {
	hit0 = new AliTPChit(*hit); //Make copy of hit
	hit = hit0;
      }
      else continue;
      AliTPChit *hit1=(AliTPChit*)TPC->NextHit();	
      if (hit1==0) continue;
      /*      Int_t nhits=hits->GetEntriesFast();
      if (nhits==0) continue;
      AliTPChit *hit;
      Int_t j;
      for (j=0; j<nhits-1; j++) {
         hit=(AliTPChit*)hits->UncheckedAt(j);
         if (hit->fQ==0.) break;
      }
      if (j==nhits-1) continue;
      AliTPChit *hit1=(AliTPChit*)hits->UncheckedAt(j+1);
      */
      if (hit1->fQ != 0.) continue;
      Int_t i=hit->Track();
      TParticle *p = (TParticle*)particles->UncheckedAt(i);
      if (p->GetFirstMother()>=0) continue;  //secondary particle
      if (good[i] < 0x1000+0x800+2+good_number) continue;
      if (p->Pt()<0.100) continue;
      if (TMath::Abs(p->Pz()/p->Pt())>0.999) continue;

      gt[nt].lab=i;
      gt[nt].code=p->GetPdgCode();
//**** px py pz - in global coordinate system, x y z - in local !
      gt[nt].px=hit->X(); gt[nt].py=hit->Y(); gt[nt].pz=hit->Z();
      Float_t cs,sn; digp->AdjustCosSin(hit1->fSector,cs,sn);
      gt[nt].x = hit1->X()*cs + hit1->Y()*sn;
      gt[nt].y =-hit1->X()*sn + hit1->Y()*cs;
      gt[nt].z = hit1->Z();
      nt++;	
	if (hit0) delete hit0;
      cerr<<i<<"                \r";
      if (nt==max) {cerr<<"Too many good tracks !\n"; break;}
   }
   delete[] good;

   delete gAlice; gAlice=0;

   file->Close();
   gBenchmark->Stop("AliTPCComparison");
   gBenchmark->Show("AliTPCComparison");   
   return nt;
}



