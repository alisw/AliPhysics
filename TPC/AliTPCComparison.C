/****************************************************************************
 *           Very important, delicate and rather obscure macro.             *
 *                                                                          *
 *               Creates list of "trackable" tracks,                        *
 *             calculates efficiency, resolutions etc.                      *
 *                                                                          *
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 * with several nice improvements by: M.Ivanov, GSI, m.ivanov@gsi.de        *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include "alles.h"
  #include "AliTPCtracker.h"
  #include "AliRunLoader.h"
  #include "AliTPCLoader.h"
  #include "AliMagF.h"
  #include "AliRun.h"
#endif

struct GoodTrackTPC {
  Int_t lab;
  Int_t code;
  Float_t px,py,pz;
  Float_t x,y,z;
};

Int_t 
good_tracks_tpc(GoodTrackTPC *gt, const Int_t max, const char* evfoldname);

extern AliRun *gAlice;

Int_t AliTPCComparison(Bool_t thcut=1.0) {

   if (gAlice) { 
       delete gAlice->GetRunLoader();
       delete gAlice;//if everything was OK here it is already NULL
       gAlice = 0x0;
    }

   cerr<<"Doing comparison...\n";

   const Int_t MAX=20000;
   Int_t good_tracks_tpc(
      GoodTrackTPC *gt, const Int_t max, 
      const char* evfoldname = AliConfig::fgkDefaultEventFolderName
   );//declaration only

   gBenchmark->Start("AliTPCComparison");

   AliRunLoader *rl = AliRunLoader::Open("galice.root","COMPARISON");
   if (!rl) {
       cerr<<"Can't start sesion !\n";
       return 1;
   }

   rl->LoadgAlice();
   if (rl->GetAliRun()) 
   AliKalmanTrack::SetConvConst(
     1000/0.299792458/rl->GetAliRun()->Field()->SolenoidField()
   );
   else {
       cerr<<"AliTPCComparison.C :Can't get AliRun !\n";
       return 1;
   }
   rl->UnloadgAlice();
  
   AliTPCLoader * tpcl = (AliTPCLoader *)rl->GetLoader("TPCLoader");
   if (tpcl == 0x0) {
       cerr<<"AliTPCComparison.C : Can not find TPCLoader\n";
       delete rl;
       return 3;
   }

   Int_t nentr=0,i=0; TObjArray tarray(MAX);
   { /*Load tracks*/
   
     tpcl->LoadTracks();
     
     TTree *tracktree=tpcl->TreeT();
     if (!tracktree) {cerr<<"Can't get a tree with TPC tracks !\n"; return 4;}

     TBranch *tbranch=tracktree->GetBranch("tracks");
     nentr=(Int_t)tracktree->GetEntries();
     AliTPCtrack *iotrack=0;
     for (i=0; i<nentr; i++) {
       iotrack=new AliTPCtrack;
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);
       tarray.AddLast(iotrack);
     }   
     tpcl->UnloadTracks();
   }

   /* Generate a list of "good" tracks */

   GoodTrackTPC gt[MAX];
   Int_t ngood=0;
   ifstream in("good_tracks_tpc");
   if (in) {
      cerr<<"Reading good tracks...\n";
      while (in>>gt[ngood].lab>>gt[ngood].code>>
                 gt[ngood].px>>gt[ngood].py>>gt[ngood].pz>>
                 gt[ngood].x >>gt[ngood].y >>gt[ngood].z) {
	Double_t rin = TMath::Sqrt( gt[ngood].x*gt[ngood].x+gt[ngood].y*gt[ngood].y);
	if (rin<1) continue;
	Double_t theta = gt[ngood].z/rin;
	//theta =0;
	if (TMath::Abs(theta)<thcut){ 
	  ngood++;
	  cerr<<ngood<<'\r';
	  if (ngood==MAX) {
	    cerr<<"Too many good tracks !\n";
	    break;
	  }	
	//ngood++;
        // cerr<<ngood<<'\r';
        // if (ngood==MAX) {
        //    cerr<<"Too many good tracks !\n";
        //    break;
	}
      }
      if (!in.eof()) cerr<<"Read error (good_tracks_tpc) !\n";
   } else {
     cerr<<"Marking good tracks (this will take a while)...\n";
      ngood=good_tracks_tpc(gt,MAX,"COMPARISON");
      ofstream out("good_tracks_tpc");
      if (out) {
         for (Int_t ngd=0; ngd<ngood; ngd++)            
	    out<<gt[ngd].lab<<' '<<gt[ngd].code<<' '<<
                 gt[ngd].px<<' '<<gt[ngd].py<<' '<<gt[ngd].pz<<' '<<
                 gt[ngd].x <<' '<<gt[ngd].y <<' '<<gt[ngd].z <<endl;
      } else cerr<<"Can not open file (good_tracks_tpc) !\n";
      out.close();
   }


   TH1F *hp=new TH1F("hp","PHI resolution",50,-20.,20.); hp->SetFillColor(4);
   TH1F *hl=new TH1F("hl","LAMBDA resolution",50,-20,20);hl->SetFillColor(4);
   TH1F *hpt=new TH1F("hpt","Relative Pt resolution",30,-10.,10.); 
   hpt->SetFillColor(2); 
   TH1F *hmpt=new TH1F("hmpt","Relative Pt resolution (pt>4GeV/c)",30,-60,60); 
   hmpt->SetFillColor(6);

   TH1F *hgood=new TH1F("hgood","Good tracks",30,0.1,6.1);    
   TH1F *hfound=new TH1F("hfound","Found tracks",30,0.1,6.1);
   TH1F *hfake=new TH1F("hfake","Fake tracks",30,0.1,6.1);
   TH1F *hg=new TH1F("hg","",30,0.1,6.1); //efficiency for good tracks
   hg->SetLineColor(4); hg->SetLineWidth(2);
   TH1F *hf=new TH1F("hf","Efficiency for fake tracks",30,0.1,6.1);
   hf->SetFillColor(1); hf->SetFillStyle(3013); hf->SetLineWidth(2);

   TH1F *he =new TH1F("he","dE/dX for pions with 0.4<p<0.5 GeV/c",50,0.,100.);
   TH2F *hep=new TH2F("hep","dE/dX vs momentum",50,0.,2.,50,0.,400.);
   hep->SetMarkerStyle(8);
   hep->SetMarkerSize(0.4);

   //MI change
   Int_t mingood=ngood;
   Int_t track_notfound[MAX], itrack_notfound=0;
   Int_t track_fake[MAX], itrack_fake=0;
   Int_t track_multifound[MAX], track_multifound_n[MAX], itrack_multifound=0;

   while (ngood--) {
      Int_t lab=gt[ngood].lab,tlab=-1;
      Float_t ptg=
      TMath::Sqrt(gt[ngood].px*gt[ngood].px + gt[ngood].py*gt[ngood].py);

      if (ptg<1e-33) continue; // for those not crossing 0 pad row

      hgood->Fill(ptg);

      AliTPCtrack *track=0;
      for (i=0; i<nentr; i++) {
          track=(AliTPCtrack*)tarray.UncheckedAt(i);
          tlab=track->GetLabel();
          if (lab==TMath::Abs(tlab)) break;
      }
      if (i==nentr) {
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
	track_multifound[itrack_multifound]=lab;
	track_multifound_n[itrack_multifound]=micount;
	itrack_multifound++;
      }
      
      Double_t xk=gt[ngood].x;
      track->PropagateTo(xk);

      if (lab==tlab) hfound->Fill(ptg);
      else { 
	track_fake[itrack_fake++]=lab;
	hfake->Fill(ptg); 
      }

      Double_t par[5]; track->GetExternalParameters(xk,par);
      Float_t phi=TMath::ASin(par[2]) + track->GetAlpha();
      if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
      if (phi>=TMath::Pi()) phi-=2*TMath::Pi();
      Float_t lam=TMath::ATan(par[3]); 
      Float_t pt_1=TMath::Abs(par[4]);

      if (TMath::Abs(gt[ngood].code)==11 && ptg>4.) {
         hmpt->Fill((pt_1 - 1/ptg)/(1/ptg)*100.);
      } else {
         Float_t phig=TMath::ATan2(gt[ngood].py,gt[ngood].px);
         hp->Fill((phi - phig)*1000.);

         Float_t lamg=TMath::ATan2(gt[ngood].pz,ptg);
         hl->Fill((lam - lamg)*1000.);

         hpt->Fill((pt_1 - 1/ptg)/(1/ptg)*100.);
      }

      Float_t mom=1./(pt_1*TMath::Cos(lam));
      Float_t dedx=track->GetdEdx();
      hep->Fill(mom,dedx,1.);
      if (TMath::Abs(gt[ngood].code)==211)
	 if (mom>0.4 && mom<0.5) {
            he->Fill(dedx,1.);
         }

   }

   cout<<"\nList of Not found tracks :\n";
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
       cout<<"id.   "<<track_multifound[i]
           <<"     found - "<<track_multifound_n[i]<<"times\n";
   }

   Stat_t ng=hgood->GetEntries(), nf=hfound->GetEntries();
   if (ng!=0) cerr<<"\n\nIntegral efficiency is about "<<nf/ng*100.<<" %\n";
   cout<<"Total number of found tracks ="<<nentr<<endl;
   cout<<"Total number of \"good\" tracks ="
       <<mingood<<"   (selected for comparison: "<<ng<<')'<<endl<<endl;


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

   TLine *line1 = new TLine(0.1,1.0,6.1,1.0); line1->SetLineStyle(4);
   line1->Draw("same");
   TLine *line2 = new TLine(0.1,0.9,6.1,0.9); line2->SetLineStyle(4);
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

   delete rl;
   return 0;
}


Int_t 
good_tracks_tpc(GoodTrackTPC *gt, const Int_t max, const char* evfoldname) {
   Int_t nt=0;

   AliRunLoader *rl = AliRunLoader::GetRunLoader(evfoldname);
   if (rl == 0x0) {
      ::Fatal("AliTPCComparison.C::good_tracks_tpc",
              "Can not find Run Loader in Folder Named %s",
              evfoldname);
   }
   AliTPCLoader *tpcl = (AliTPCLoader *)rl->GetLoader("TPCLoader");
   if (tpcl == 0x0) {
       cerr<<"AliTPCHits2Digits.C : Can not find TPCLoader\n";
       delete rl;
       return 0;
   }
   
   rl->LoadgAlice();
   
   AliTPC *TPC=(AliTPC*)rl->GetAliRun()->GetDetector("TPC");
   Int_t ver = TPC->IsVersion(); 
   cerr<<"TPC version "<<ver<<" has been found !\n";

   rl->CdGAFile();
   AliTPCParamSR *digp=(AliTPCParamSR*)gDirectory->Get("75x40_100x60_150x60");
   if (!digp) { cerr<<"TPC parameters have not been found !\n"; exit(6); }
   TPC->SetParam(digp);

   rl->LoadHeader();

   Int_t np = rl->GetHeader()->GetNtrack();
      

   Int_t nrow_up=digp->GetNRowUp();
   Int_t nrows=digp->GetNRowLow()+nrow_up;
   Int_t zero=digp->GetZeroSup();
   Int_t gap=Int_t(0.125*nrows), shift=Int_t(0.5*gap);
   Int_t good_number=Int_t(0.4*nrows);

   Int_t *good=new Int_t[np];
   Int_t i;
   for (i=0; i<np; i++) good[i]=0;


   switch (ver) {
   case 1:
     {
      tpcl->LoadRecPoints();      
      AliTPCClustersArray *pca=new AliTPCClustersArray, &ca=*pca;
      ca.Setup(digp);
      ca.SetClusterType("AliTPCcluster");
      ca.ConnectTree(tpcl->TreeR());
      Int_t nrows=Int_t(ca.GetTree()->GetEntries());
      for (Int_t n=0; n<nrows; n++) {
          AliSegmentID *s=ca.LoadEntry(n);
          Int_t sec,row;
          digp->AdjustSectorRow(s->GetID(),sec,row);
          AliTPCClustersRow &clrow = *ca.GetRow(sec,row);
          Int_t ncl=clrow.GetArray()->GetEntriesFast();
          while (ncl--) {
              AliTPCcluster *c=(AliTPCcluster*)clrow[ncl];
              Int_t lab=c->GetLabel(0);
              if (lab<0) continue; //noise cluster
              lab=TMath::Abs(lab);

              if (sec>=digp->GetNInnerSector())
                 if (row==nrow_up-1) good[lab]|=0x4000;
              if (sec>=digp->GetNInnerSector())
                 if (row==nrow_up-1-gap) good[lab]|=0x1000;

              if (sec>=digp->GetNInnerSector())
                 if (row==nrow_up-1-shift) good[lab]|=0x2000;
              if (sec>=digp->GetNInnerSector())
                 if (row==nrow_up-1-gap-shift) good[lab]|=0x800;

              good[lab]++;
          }
          ca.ClearRow(sec,row);
      }
      delete pca;
      tpcl->UnloadRecPoints();
     }
     break;
   case 2:
     {
      tpcl->LoadDigits();
      TTree *TD=tpcl->TreeD();
      
      AliSimDigits da, *digits=&da;
      TD->GetBranch("Segment")->SetAddress(&digits);

      Int_t *count = new Int_t[np];
      Int_t i;
      for (i=0; i<np; i++) count[i]=0;

      Int_t sectors_by_rows=(Int_t)TD->GetEntries();
      for (i=0; i<sectors_by_rows; i++) {
          if (!TD->GetEvent(i)) continue;
          Int_t sec,row;
          digp->AdjustSectorRow(digits->GetID(),sec,row);
          cerr<<sec<<' '<<row<<"                                     \r";
          digits->First();
          do { //Many thanks to J.Chudoba who noticed this
              Int_t it=digits->CurrentRow(), ip=digits->CurrentColumn();
              Short_t dig = digits->GetDigit(it,ip);
              Int_t idx0=digits->GetTrackID(it,ip,0); 
              Int_t idx1=digits->GetTrackID(it,ip,1);
              Int_t idx2=digits->GetTrackID(it,ip,2);
              if (idx0>=0 && dig>=zero && idx0<np) count[idx0]+=1;
              if (idx1>=0 && dig>=zero && idx1<np) count[idx1]+=1;
              if (idx2>=0 && dig>=zero && idx2<np) count[idx2]+=1;
          } while (digits->Next());
          for (Int_t j=0; j<np; j++) {
              if (count[j]>1) {
                 if (sec>=digp->GetNInnerSector())
		   if (row==nrow_up-1    ) good[j]|=0x4000;
                 if (sec>=digp->GetNInnerSector())
		   if (row==nrow_up-1-gap) good[j]|=0x1000;

                 if (sec>=digp->GetNInnerSector())
		   if (row==nrow_up-1-shift) good[j]|=0x2000;
                 if (sec>=digp->GetNInnerSector())
		   if (row==nrow_up-1-gap-shift) good[j]|=0x800;
                 good[j]++;
              }
              count[j]=0;
          }
      }
      delete[] count;
      tpcl->UnloadDigits();
     }
      break;
   default:
      cerr<<"Invalid TPC version !\n";
      delete rl;
      exit(7);
   }

   rl->LoadKinematics();
   AliStack* stack = rl->Stack();
   /** select tracks which are "good" enough **/
   for (i=0; i<np; i++) {
      if ((good[i]&0x5000) != 0x5000)
      if ((good[i]&0x2800) != 0x2800) continue;
      if ((good[i]&0x7FF ) < good_number) continue;

      TParticle *p = (TParticle*)stack->Particle(i);
      if (p == 0x0)
       {
         cerr<<"Can not get particle "<<i<<endl;
         continue;
       }
      if (p->Pt()<0.100) continue;
      if (TMath::Abs(p->Pz()/p->Pt())>0.999) continue;

      Int_t j=p->GetFirstMother();
      if (j>=0) {
        TParticle *pp = (TParticle*)stack->Particle(j);
        if (pp == 0x0)
         {
           cerr<<"Can not get particle "<<j<<endl;
           continue;
         }
        if (pp->GetFirstMother()>=0) continue;//only one decay is allowed
	/*  for cascade hyperons only
        Int_t jj=pp->GetFirstMother();
        if (jj>=0) {
          TParticle *ppp = (TParticle*)stack->Particle(jj);
          if (ppp->GetFirstMother()>=0) continue;//two decays are allowed
        }
	*/
      }

      gt[nt].lab=i;
      gt[nt].code=p->GetPdgCode();
      gt[nt].px=0.; gt[nt].py=0.; gt[nt].pz=0.;
      gt[nt].x=0.; gt[nt].y=0.; gt[nt].z=0.;
      nt++;
      if (nt==max) {cerr<<"Too many good tracks !\n"; break;}
      cerr<<np-i<<"                \r";
   }
   rl->UnloadKinematics();

   /** check if there is also information at the entrance of the TPC **/
   
   tpcl->LoadHits();
   TTree *TH=tpcl->TreeH();
   TPC->SetTreeAddress();
   np=(Int_t)TH->GetEntries();
   for (i=0; i<np; i++) {
      TPC->ResetHits();
      TH->GetEvent(i);
      AliTPChit *phit = (AliTPChit*)TPC->FirstHit(-1);
      for ( ; phit; phit=(AliTPChit*)TPC->NextHit() ) {
	 if (phit->fQ !=0. ) continue;

         Double_t px=phit->X(), py=phit->Y(), pz=phit->Z();

         if ((phit=(AliTPChit*)TPC->NextHit())==0) break;
         if (phit->fQ != 0.) continue;

         Double_t x=phit->X(), y=phit->Y(), z=phit->Z();
	 if (TMath::Sqrt(x*x+y*y)>90.) continue;

         Int_t j, lab=phit->Track();
         for (j=0; j<nt; j++) {if (gt[j].lab==lab) break;}
         if (j==nt) continue;         

	 // (px,py,pz) - in global coordinate system, (x,y,z) - in local !
         gt[j].px=px; gt[j].py=py; gt[j].pz=pz;
         Float_t cs,sn; digp->AdjustCosSin(phit->fSector,cs,sn);
         gt[j].x = x*cs + y*sn;
         gt[j].y =-x*sn + y*cs;
         gt[j].z = z;
      }
      cerr<<np-i<<"                \r";
   }

   delete[] good;
   
   tpcl->UnloadHits();
   rl->UnloadgAlice();

   gBenchmark->Stop("AliTPCComparison");
   gBenchmark->Show("AliTPCComparison");   
   return nt;
}
