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
  #include <Riostream.h>
  #include <fstream.h>

  #include "AliRun.h"
  #include "AliHeader.h"
  #include "AliStack.h"
  #include "AliRunLoader.h"
  #include "AliLoader.h"
  #include "AliITSLoader.h"
  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITStrackerV2.h"
  #include "AliITStrackV2.h"
  #include "AliITSclusterV2.h"
  #include "AliMagF.h"
  #include "AliESD.h"

  #include "TFile.h"
  #include "TKey.h"
  #include "TTree.h"
  #include "TH1.h"
  #include "TH2.h"
  #include "TObjArray.h"
  #include "TStyle.h"
  #include "TCanvas.h"
  #include "TLine.h"
  #include "TText.h"
  #include "TParticle.h"
#endif

struct GoodTrackITS {
  Int_t lab;
  Int_t code;
  Float_t px,py,pz;
  Float_t x,y,z;
};

Int_t good_tracks_its(GoodTrackITS *gt, const Int_t max, 
const char* evfoldname = AliConfig::fgkDefaultEventFolderName);

extern AliRun *gAlice;

Int_t AliITSComparisonV2() {
   cerr<<"Doing comparison...\n";
   if (gAlice) {
      delete gAlice->GetRunLoader();
      delete gAlice; 
      gAlice=0;
   }
   
   AliRunLoader *rl = AliRunLoader::Open("galice.root");
   if (!rl) {
       cerr<<"AliITSComparisonV2.C :Can't start sesion !\n";
       return 1;
   }
   rl->LoadgAlice();
   if (rl->GetAliRun())
   AliKalmanTrack::SetConvConst(
     1000/0.299792458/rl->GetAliRun()->Field()->SolenoidField()
   );
   else {
      cerr<<"AliITSComparisonV2.C :Can't get AliRun !\n";
      return 1;
   }
   //rl->UnloadgAlice();
   

   /* Generate a list of "good" tracks */
   const Int_t MAX=15000;
   GoodTrackITS gt[MAX];
   Int_t ngood=0;
   ifstream in("good_tracks_its");
   if (in) {
      cerr<<"Reading good tracks...\n";
      while (in>>gt[ngood].lab>>gt[ngood].code>>
                 gt[ngood].px>>gt[ngood].py>>gt[ngood].pz>>
                 gt[ngood].x >>gt[ngood].y >>gt[ngood].z) {
         ngood++;
         cerr<<ngood<<'\r';
         if (ngood==MAX) {
            cerr<<"Too many good tracks !\n";
            break;
         }
      }
      if (!in.eof()) cerr<<"Read error (good_tracks_its) !\n";
   } else {
      cerr<<"Marking good tracks (this will take a while)...\n";
      ngood=good_tracks_its(gt,MAX,AliConfig::fgkDefaultEventFolderName);
      ofstream out("good_tracks_its");
      if (out) {
	for (Int_t ngd=0; ngd<ngood; ngd++)
	  out<<gt[ngd].lab<<' '<<gt[ngd].code<<' '
             <<gt[ngd].px<<' '<<gt[ngd].py<<' '<<gt[ngd].pz<<' '
             <<gt[ngd].x <<' '<<gt[ngd].y <<' '<<gt[ngd].z <<endl;
      } else cerr<<"Can not open file (good_tracks_its) !\n";
      out.close();
   }

   TObjArray tarray(2000);
   { /*Load tracks*/
   TFile *ef=TFile::Open("AliESDits.root");
   if ((!ef)||(!ef->IsOpen())) {
      ::Fatal("AliITSComparisonV2.C","Can't open AliESDits.root !");
   }
   TKey *key=0;
   TIter next(ef->GetListOfKeys());
   if ((key=(TKey*)next())!=0) {
     AliESD *event=(AliESD*)key->ReadObj();
     Int_t ntrk=event->GetNumberOfTracks();
     for (Int_t i=0; i<ntrk; i++) {
        AliESDtrack *t=event->GetTrack(i);
	UInt_t status=t->GetStatus();
	UInt_t flags=AliESDtrack::kTPCin|AliESDtrack::kITSin;

        if ((status&AliESDtrack::kITSrefit)==0)
	  if ((status&flags)!=flags) continue;

        AliITStrackV2 *iotrack=0;
        iotrack=new AliITStrackV2(*t);
        //if (t->GetConstrainedChi2()>=20) continue;   //  constrained 
        //else iotrack=new AliITStrackV2(*t,kTRUE);    //     track
        if ((status&flags)==flags) {
           iotrack->PropagateTo(3.,0.0028,65.19);
           iotrack->PropagateToVertex();
        }
        tarray.AddLast(iotrack);
     }
     delete event;
   }
   ef->Close();
   }
   Int_t nentr=tarray.GetEntriesFast();

   TH1F *hp=new TH1F("hp","PHI resolution",50,-20.,20.); hp->SetFillColor(4);
   TH1F *hl=new TH1F("hl","LAMBDA resolution",50,-20,20);hl->SetFillColor(4);
   TH1F *hpt=new TH1F("hpt","Relative Pt resolution",30,-10.,10.); 
   hpt->SetFillColor(2); 
   TH1F *hmpt=new TH1F("hmpt","Transverse impact parameter",30,-300,300); 
   hmpt->SetFillColor(6);
   TH1F *hz=new TH1F("hz","Longitudinal impact parameter",30,-300,300); 

   AliITStrackV2 *trk=(AliITStrackV2*)tarray.UncheckedAt(0);
   Double_t pmin=0.1*(100/0.299792458/0.2/trk->GetConvConst());
   Double_t pmax=6.0+pmin;

   TH1F *hgood=new TH1F("hgood","Good tracks",30,pmin,pmax);    
   TH1F *hfound=new TH1F("hfound","Found tracks",30,pmin,pmax);
   TH1F *hfake=new TH1F("hfake","Fake tracks",30,pmin,pmax);
   TH1F *hg=new TH1F("hg","",30,pmin,pmax); //efficiency for good tracks
   hg->SetLineColor(4); hg->SetLineWidth(2);
   TH1F *hf=new TH1F("hf","Efficiency for fake tracks",30,pmin,pmax);
   hf->SetFillColor(1); hf->SetFillStyle(3013); hf->SetLineWidth(2);

   //TH1F *hptw=new TH1F("hptw","Weghted pt",30,pmax,pmin);

   TH1F *he =new TH1F("he","dE/dX for pions with 0.4<p<0.5 GeV/c",50,0.,100.);
   TH2F *hep=new TH2F("hep","dE/dX vs momentum",50,0.,2.,50,0.,200.);
   hep->SetMarkerStyle(8);
   hep->SetMarkerSize(0.4);


   Int_t notf[MAX], nnotf=0;
   Int_t fake[MAX], nfake=0;
   Int_t mult[MAX], numb[MAX], nmult=0;

   Int_t ng;
   for (ng=0; ng<ngood; ng++) {
      Int_t lab=gt[ng].lab, tlab=-1;
      Double_t pxg=gt[ng].px, pyg=gt[ng].py, pzg=gt[ng].pz;
      Double_t ptg=TMath::Sqrt(pxg*pxg+pyg*pyg);


      if (ptg>pmin) hgood->Fill(ptg);

      AliITStrackV2 *track=0;
      Int_t cnt=0;
      Int_t j;
      for (j=0; j<nentr; j++) {
          AliITStrackV2 *trk=(AliITStrackV2*)tarray.UncheckedAt(j);
          Int_t lbl=trk->GetLabel();
          if (lab==TMath::Abs(lbl)) {
	    if (cnt==0) {track=trk; tlab=lbl;}
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

      if (ptg>pmin) {
        if (lab==tlab) hfound->Fill(ptg);
        else {
          fake[nfake++]=lab;
          hfake->Fill(ptg); 
        }
      }

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

      //hptw->Fill(ptg,TMath::Abs(d));

      Double_t z=10000*track->GetZ();
      hz->Fill(z);

      hpt->Fill((pt_1 - 1/ptg)/(1/ptg)*100.);

      Float_t mom=1./(pt_1*TMath::Cos(lam));
      Float_t dedx=track->GetdEdx();
      hep->Fill(mom,dedx,1.);
      if (TMath::Abs(gt[ng].code)==211)
         if (mom>0.4 && mom<0.5) he->Fill(dedx,1.);
   }


   cout<<"\nList of Not found tracks :\n";
   for (ng=0; ng<nnotf; ng++){
     cout<<notf[ng]<<"\t";
     if ((ng%9)==8) cout<<"\n";
   }
   cout<<"\n\nList of fake  tracks :\n";
   for (ng=0; ng<nfake; ng++){
     cout<<fake[ng]<<"\t";
     if ((ng%9)==8) cout<<"\n";
   }
   cout<<"\n\nList of multiple found tracks :\n";
   for (ng=0; ng<nmult; ng++) {
       cout<<"id.   "<<mult[ng]
           <<"     found - "<<numb[ng]<<"times\n";
   }
   cout<<endl;

   cerr<<"Number of good tracks : "<<ngood<<endl;
   cerr<<"Number of found tracks : "<<nentr<<endl;

   ng=(Int_t)hgood->GetEntries(); //cerr<<"Good tracks "<<ng<<endl;
   if (ng!=0)  
   cerr<<"Integral efficiency is about "<<hfound->GetEntries()/ng*100.<<" %\n";
   cerr<<endl;

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

   TLine *line1 = new TLine(pmin,1.0,pmax,1.0); line1->SetLineStyle(4);
   line1->Draw("same");
   TLine *line2 = new TLine(pmin,0.9,pmax,0.9); line2->SetLineStyle(4);
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
   
   return 0;
}

Int_t good_tracks_its(GoodTrackITS *gt, const Int_t max, const char* evfoldname) {
   AliRunLoader* rl = AliRunLoader::GetRunLoader(evfoldname);
   if (rl == 0x0) {
      ::Fatal("AliITSComparisonV2.C::good_tracks_its",
              "Can not find Run Loader in Folder Named %s",
              evfoldname);
   }

   AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0) {
       cerr<<"AliITSComparisonV2.C : Can not find ITSLoader\n";
       delete rl;
       return 3;
   }
   
   rl->LoadgAlice();
   rl->LoadHeader();
   Int_t np = rl->GetHeader()->GetNtrack();

   Int_t *good=new Int_t[np];
   Int_t k;
   for (k=0; k<np; k++) good[k]=0;

   AliITS *ITS=(AliITS*)rl->GetAliRun()->GetDetector("ITS");
   if (!ITS) {
      cerr<<"can't get ITS !\n"; exit(8);
   }
   AliITSgeom *geom=ITS->GetITSgeom();
   if (!geom) {
      cerr<<"can't get ITS geometry !\n"; exit(9);
   }

   itsl->LoadRecPoints();
   TTree *cTree=itsl->TreeR();
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
     cTree->GetEvent(k);
     Int_t ncl=clusters->GetEntriesFast(); if (ncl==0) continue;
     Int_t lay,lad,det;  geom->GetModuleId(k,lay,lad,det);
     if (lay<1 || lay>6) {
	cerr<<"wrong layer !\n"; exit(10);
     }
     while (ncl--) {
        AliITSclusterV2 *pnt=(AliITSclusterV2*)clusters->UncheckedAt(ncl);
        Int_t l0=pnt->GetLabel(0);
	  if (l0>=np) {cerr<<"Wrong label: "<<l0<<endl; continue;}
        Int_t l1=pnt->GetLabel(1);
	  if (l1>=np) {cerr<<"Wrong label: "<<l1<<endl; continue;}
        Int_t l2=pnt->GetLabel(2);
	  if (l2>=np) {cerr<<"Wrong label: "<<l2<<endl; continue;}
        Int_t mask=1<<(lay-1);
        if (l0>=0) good[l0]|=mask; 
        if (l1>=0) good[l1]|=mask; 
        if (l2>=0) good[l2]|=mask;
     }
   }
   clusters->Delete(); delete clusters;
   itsl->UnloadRawClusters();
   
   ifstream in("good_tracks_tpc");
   if (!in) {
     cerr<<"can't get good_tracks_tpc !\n"; exit(11);
   }
   
   rl->LoadKinematics();
   AliStack* stack = rl->Stack();
   Int_t nt=0;
   Double_t px,py,pz,x,y,z;
   Int_t code,lab;
   while (in>>lab>>code>>px>>py>>pz>>x>>y>>z) {
      if (good[lab] != 0x3F) continue;
      TParticle *p = (TParticle*)stack->Particle(lab);
      if (p == 0x0) {
         cerr<<"Can not get particle "<<lab<<endl;
         nt++;
         if (nt==max) {cerr<<"Too many good tracks !\n"; break;}
         continue;
      }
      gt[nt].lab=lab;
      gt[nt].code=p->GetPdgCode();
//**** px py pz - in global coordinate system
      gt[nt].px=p->Px(); gt[nt].py=p->Py(); gt[nt].pz=p->Pz();
      gt[nt].x=gt[nt].y=gt[nt].z=0.;
      nt++;
      if (nt==max) {cerr<<"Too many good tracks !\n"; break;}
   }

   delete[] good;

   rl->UnloadKinematics();
   rl->UnloadHeader();
   rl->UnloadgAlice();

   return nt;
}


