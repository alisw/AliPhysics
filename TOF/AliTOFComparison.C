/****************************************************************************
 *      This macro estimates efficiency of matching with the TOF.           *
 *      TOF "Good" tracks are those originating from the primary vertex,    *
 *      being "good" in the ITS and having at least one digit in the TOF.   * 
 *         (To get the list of "good" tracks one should first run           *
 *          AliTPCComparison.C and AliITSComparisonV2.C macros)             *
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

  #include "AliESD.h"
  #include "AliTOFdigit.h"

  #include "TKey.h"
  #include "TFile.h"
  #include "TTree.h"
  #include "TH1.h"
  #include "TClonesArray.h"
  #include "TStyle.h"
  #include "TCanvas.h"
  #include "TLine.h"
  #include "TText.h"
  #include "TParticle.h"
#endif

struct GoodTrackTOF {
  Int_t lab;
  Int_t code;
  Float_t px,py,pz;
  Float_t x,y,z;
};

extern AliRun *gAlice;

Int_t AliTOFComparison() {
   Int_t good_tracks_tof(GoodTrackTOF *gt, const Int_t max);

   cerr<<"Doing comparison...\n";

   if (gAlice) { 
     delete gAlice->GetRunLoader();
     delete gAlice;//if everything was OK here it is already NULL
     gAlice = 0x0;
   }

   TClonesArray dummy("AliTOFdigit",10000), *digits=&dummy;   

   {
   AliRunLoader *rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"Can not open session"<<endl;
      return 1;
   }
   AliLoader* tofl = rl->GetLoader("TOFLoader");
   if (tofl == 0x0) {
      cerr<<"Can not get the TOF loader"<<endl;
      return 2;
   }
   tofl->LoadDigits("read");

   rl->GetEvent(0);

   TTree *tofTree=tofl->TreeD();
   if (!tofTree) {
      cerr<<"Can't get the TOF cluster tree !\n";
      return 3;
   } 
   TBranch *branch=tofTree->GetBranch("TOF");
   if (!branch) { 
      cerr<<"Can't get the branch with the TOF digits !\n";
      return 4;
   }
   branch->SetAddress(&digits);

   tofTree->GetEvent(0);

   delete rl;
   }

   Int_t nd=digits->GetEntriesFast();
   cerr<<"Number of digits: "<<nd<<endl;

   const Int_t MAX=15000;
   GoodTrackTOF gt[MAX];
   Int_t ngood=0;
   ifstream in("good_tracks_tof");
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
      ngood=good_tracks_tof(gt,MAX);
      ofstream out("good_tracks_tof");
      if (out) {
        for (Int_t ngd=0; ngd<ngood; ngd++)
          out<<gt[ngd].lab<<' '<<gt[ngd].code<<' '
             <<gt[ngd].px<<' '<<gt[ngd].py<<' '<<gt[ngd].pz<<' '
             <<gt[ngd].x <<' '<<gt[ngd].y <<' '<<gt[ngd].z <<endl;
      } else cerr<<"Can not open file (good_tracks_tof) !\n";
      out.close();
   }

   Double_t pmin=0.2;
   Double_t pmax=4.0;

   TH1F *hgood=new TH1F("hgood","Good tracks",30,pmin,pmax);    
   TH1F *hfound=new TH1F("hfound","Matched tracks",30,pmin,pmax);
   TH1F *hfake=new TH1F("hfake","Mismatched tracks",30,pmin,pmax);
   TH1F *hgp=new TH1F("hgp","",30,pmin,pmax); //efficiency for good tracks
   hgp->SetLineColor(4); hgp->SetLineWidth(2);
   TH1F *hfp=new TH1F("hfp","Probability of mismatching",30,pmin,pmax);
   hfp->SetFillColor(1); hfp->SetFillStyle(3013); hfp->SetLineWidth(2);

   TH1F *hgoo=new TH1F("hgoo","Good tracks",30,-1,1);    
   TH1F *hfoun=new TH1F("hfoun","Matched tracks",30,-1,1);
   TH1F *hfak=new TH1F("hfak","Mismatched tracks",30,-1,1);
   TH1F *hgl=new TH1F("hgl","",30,-1,1); //efficiency for good tracks
   hgl->SetLineColor(4); hgl->SetLineWidth(2);
   TH1F *hfl=new TH1F("hfl","Probability of mismatching",30,-1,1);
   hfl->SetFillColor(1); hfl->SetFillStyle(3013); hfl->SetLineWidth(2);

   TFile *ef=TFile::Open("AliESDs.root");
   if (!ef->IsOpen()) {cerr<<"Can't AliESDs.root !\n"; return 1;}

   TIter next(ef->GetListOfKeys());
   TKey *key=0;
   Int_t nev=0;
   while ((key=(TKey*)next())!=0) {
     cerr<<"Processing event number : "<<nev++<<endl;

     AliESD *event=(AliESD*)key->ReadObj();
     Int_t ntrk=event->GetNumberOfTracks();
     cerr<<"Number of ESD tracks : "<<ntrk<<endl; 

     Int_t matched=0;
     Int_t mismatched=0;
     for (Int_t i=0; i<ngood; i++) {
         Int_t lab=gt[i].lab;
         Double_t pxg=gt[i].px, pyg=gt[i].py, pzg=gt[i].pz;
         Double_t ptg=TMath::Sqrt(pxg*pxg+pyg*pyg);

	 if (ptg<0.1) continue;

         Double_t tgl=pzg/ptg; //tan(lambda)

         if (ptg>pmin) { hgood->Fill(ptg); hgoo->Fill(tgl); }

         Int_t j;
	 AliESDtrack *t=0;
         for (j=0; j<ntrk; j++) {
             AliESDtrack *tt=event->GetTrack(j);
             if (lab!=TMath::Abs(tt->GetLabel())) continue;
             t=tt;
             //if ((tt->GetStatus()&AliESDtrack::kTOFpid) == 0) continue;
             if (tt->GetTOFsignal() < 0) continue;
             UInt_t idx=tt->GetTOFcluster();
             if ((Int_t)idx>=nd) {
	       cerr<<"Wrong digit index ! "<<idx<<endl;
               return 5;
             }
             AliTOFdigit *dig=(AliTOFdigit*)digits->UncheckedAt(idx);
             Int_t *label=dig->GetTracks();
             if (label[0]!=lab)
             if (label[1]!=lab)
	       if (label[2]!=lab) {
                 mismatched++; 
                 if (ptg>pmin) { hfake->Fill(ptg); hfak->Fill(tgl); } 
                 break;
               }
             if (ptg>pmin) { hfound->Fill(ptg); hfoun->Fill(tgl); }
             matched++;
             break;
         }
         if (j==ntrk) {
	    cerr<<"Not matched: "<<lab<<"   ";
            if (t) {
               cerr<<(t->GetStatus()&AliESDtrack::kITSout)<<' '
	           <<(t->GetStatus()&AliESDtrack::kTPCout)<<' '
	           <<(t->GetStatus()&AliESDtrack::kTRDout)<<' '
	           <<(t->GetStatus()&AliESDtrack::kTIME);
	    } else cerr<<"No ESD track !";
            cerr<<endl;
         }
     }

     cerr<<"Number of good tracks: "<<ngood<<endl;
     cerr<<"Number of matched tracks: "<<matched<<endl;
     cerr<<"Number of mismatched tracks: "<<mismatched<<endl;
     if (ngood!=0) cerr<<"Efficiency: "<<Float_t(matched)/ngood<<endl;

     hfound->Sumw2(); hgood->Sumw2(); hfake->Sumw2();
     hgp->Divide(hfound,hgood,1,1.,"b");
     hfp->Divide(hfake,hgood,1,1.,"b");
     hgp->SetMaximum(1.4);
     hgp->SetYTitle("Matching efficiency");
     hgp->SetXTitle("Pt (GeV/c)");

     hfoun->Sumw2(); hgoo->Sumw2(); hfak->Sumw2();
     hgl->Divide(hfoun,hgoo,1,1.,"b");
     hfl->Divide(hfak,hgoo,1,1.,"b");
     hgl->SetMaximum(1.4);
     hgl->SetYTitle("Matching efficiency");
     hgl->SetXTitle("Tan(lambda)");

     TCanvas *c1=new TCanvas("c1","",0,0,600,900);
     c1->Divide(1,2);

     c1->cd(1);

     hgp->Draw();
     hfp->Draw("histsame");
     TLine *line1 = new TLine(pmin,1.0,pmax,1.0); line1->SetLineStyle(4);
     line1->Draw("same");
     TLine *line2 = new TLine(pmin,0.9,pmax,0.9); line2->SetLineStyle(4);
     line2->Draw("same");

     c1->cd(2);

     hgl->Draw();
     hfl->Draw("histsame");
     TLine *line3 = new TLine(-1,1.0,1,1.0); line3->SetLineStyle(4);
     line3->Draw("same");
     TLine *line4 = new TLine(-1,0.9,1,0.9); line4->SetLineStyle(4);
     line4->Draw("same");

     break;
   }

   TFile fc("AliTOFComparison.root","RECREATE");
   c1->Write();
   fc.Close();

   return 0;
}

Int_t good_tracks_tof(GoodTrackTOF *gt, const Int_t max) {
   ifstream in("good_tracks_its");
   if (!in) {
     cerr<<"Can't get good_tracks_its !\n"; exit(11);
   }   

   AliRunLoader *rl = AliRunLoader::Open("galice.root","COMPARISON");
   if (!rl) {
       cerr<<"Can't start session !\n";
       exit(1);
   }

   rl->GetEvent(0);


   rl->LoadgAlice();
   rl->LoadHeader();
   Int_t np = rl->GetHeader()->GetNtrack();

   Int_t *good=new Int_t[np];
   Int_t k;
   for (k=0; k<np; k++) good[k]=0;

   AliLoader* tofl = rl->GetLoader("TOFLoader");
   if (tofl == 0x0) {
      cerr<<"Can not get the TOF loader"<<endl;
      exit(2);
   }
   tofl->LoadDigits("read");

   TTree *dTree=tofl->TreeD();
   if (!dTree) {
      cerr<<"Can't get the TOF cluster tree !\n";
      exit(3);
   } 

   TBranch *branch=dTree->GetBranch("TOF");
   if (!branch) { 
     cerr<<"Can't get the branch with the TOF digits !\n";
     exit(4);
   }
   TClonesArray dummy("AliTOFdigit",10000), *digits=&dummy;
   branch->SetAddress(&digits);
   
  dTree->GetEvent(0);
  Int_t nd=digits->GetEntriesFast();
  cerr<<"Number of digits: "<<nd<<endl;

  for (Int_t i=0; i<nd; i++) {
    AliTOFdigit *d=(AliTOFdigit*)digits->UncheckedAt(i);
    Int_t l0=d->GetTrack(0);
       if (l0>=np) {cerr<<"Wrong label: "<<l0<<endl; continue;}
    Int_t l1=d->GetTrack(1);
       if (l1>=np) {cerr<<"Wrong label: "<<l1<<endl; continue;}
    Int_t l2=d->GetTrack(2);
       if (l2>=np) {cerr<<"Wrong label: "<<l2<<endl; continue;}
    if (l0>=0) good[l0]++; 
    if (l1>=0) good[l1]++; 
    if (l2>=0) good[l2]++;
  }

   
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();
  Int_t nt=0;
  Double_t px,py,pz,x,y,z;
  Int_t code,lab;
  while (in>>lab>>code>>px>>py>>pz>>x>>y>>z) {
      if (good[lab] == 0) continue;
      TParticle *p = (TParticle*)stack->Particle(lab);
      if (p == 0x0) {
         cerr<<"Can not get particle "<<lab<<endl;
         exit(1);
      }
      if (TMath::Abs(p->Vx())>0.1) continue;
      if (TMath::Abs(p->Vy())>0.1) continue;
      if (TMath::Abs(p->Vz())>0.1) continue;

      gt[nt].lab=lab;
      gt[nt].code=p->GetPdgCode();
//**** px py pz - in global coordinate system
      gt[nt].px=p->Px(); gt[nt].py=p->Py(); gt[nt].pz=p->Pz();
      gt[nt].x=p->Vx(); gt[nt].y=p->Vy(); gt[nt].z=p->Vz();
      nt++;
      if (nt==max) {cerr<<"Too many good tracks !\n"; break;}
  }

  delete[] good;

  rl->UnloadKinematics();
  rl->UnloadHeader();
  rl->UnloadgAlice();
  delete rl;

  return nt;
}
