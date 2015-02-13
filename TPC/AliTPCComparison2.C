/// \file AliTPCComparison.C
/// \brief Very important, delicate and rather obscure macro
///
/// Creates list of "trackable" tracks, sorts tracks for matching with the ITS,
/// calculates efficiency, resolutions etc.
///
/// There is a possibility to run this macro over several events.
///
/// \author I.Belikov, CERN, Jouri.Belikov@cern.ch, M.Ivanov, GSI, m.ivanov@gsi.de

#if !defined(__CINT__) || defined(__MAKECINT__)
#include<fstream.h>
#include<TROOT.h>
#include<TCanvas.h>
#include<TClassTable.h>
#include<TClonesArray.h>
#include<TFile.h>
#include<TFolder.h>
#include<TH1.h>
#include<TObject.h>
#include<TObjArray.h>
#include<TParticle.h>
#include<TTree.h>
#include <AliMagF.h>
#include <AliRun.h>
#include <alles.h>
#endif

/*
  multievent comparison of reconstructed tracks with "exact tracks" 
  from silulation

*/


struct GoodTrackTPC {
  Int_t fEventN; //event number
  Int_t lab;
  Int_t code;
  Float_t px,py,pz;
  Float_t x,y,z;
};

enum tagprimary {kPrimaryCharged = 0x4000};

Int_t good_tracks(GoodTrackTPC *gt, Int_t max, Int_t firstev, Int_t eventn);

Int_t FindPrimaries(Int_t nparticles);

Int_t AliTPCComparison2(Int_t firstev=0, Int_t eventn=1) {

  const Int_t MAX = 20000;
  const Bool_t kOLD = kFALSE; // if tracking was done with a previous version
  /***********************************************************************/

  TFile *inkin = TFile::Open("rfio:galice.root");
// \file AliTPCComparison2.C
//  if(gAlice)delete gAlice;   COMMENTED BECAUSE OF A BUG (IN COMPILED MODE)

  gAlice = (AliRun*)inkin->Get("gAlice");
  cout<<"AliRun object found on file "<<gAlice<<endl;
  AliKalmanTrack::SetFieldMap(gAlice->Field());
  inkin->Close();
  /*
    delete gAlice;  COMMENTED BECAUSE OF A BUG IN COMPILED MODE
    gAlice = 0;
  */
  /***********************************************************************/
  cerr<<"Doing comparison...\n";
  Int_t i;
  gBenchmark->Start("AliTPCComparison2");

  TFile *cf=0; 
  AliTPCParam *digp=0;
  if(kOLD){
    cf=TFile::Open("AliTPCclusters.root");
    if (!cf->IsOpen()) {cerr<<"Can't open AliTPCclusters.root !\n"; return 1;}
    digp= (AliTPCParam*)cf->Get("75x40_100x60_150x60");
    if (!digp) { cerr<<"TPC parameters have not been found !\n"; return 2; }
  }

  TObjArray tarray(MAX);
  AliTPCtrack *iotrack=0;
  Int_t nentr= 0;
  Int_t eventptr[1000];
  TFile *tf=TFile::Open("AliTPCtracks.root");
  TTree *tracktree=0;

  eventptr[firstev]=0;
  eventptr[firstev+1]=0;
  for (Int_t event=firstev;event<eventn; event++){
    cout<<"================================================"<<endl;
    cout<<"Processing event "<<event<<endl;
    cout<<"================================================"<<endl;
    
    char   tname[100];
    if (eventn==-1) {
      sprintf(tname,"TreeT_TPC");
    }
    else {
      sprintf(tname,"TreeT_TPC_%d",event);
    }
  
    // Load tracks
    if (!tf->IsOpen()) {cerr<<"Can't open AliTPCtracks.root !\n"; return 3;}

    tracktree=(TTree*)tf->Get(tname);
    if (!tracktree) {cerr<<"Can't get a tree with TPC tracks !\n"; return 4;}
    TBranch *tbranch=tracktree->GetBranch("tracks");
    Int_t nentr0=(Int_t)tracktree->GetEntries();
    nentr+=nentr0;
    for (i=0; i<nentr0; i++) {
      iotrack=new AliTPCtrack;
      tbranch->SetAddress(&iotrack);
      tracktree->GetEvent(i);
      tarray.AddLast(iotrack);
    }   
    eventptr[event+1] = nentr;  //store end of the event
    delete tracktree;
  }
  tf->Close();
  if(kOLD)cf->Close();

  GoodTrackTPC gt[MAX];

  Int_t ngood=0;
  ifstream in("good_tracks_tpc");
  if (in) {
    cerr<<"Reading good tracks...\n";
    while (in>>gt[ngood].fEventN>>gt[ngood].lab>>gt[ngood].code>>
           gt[ngood].px>>gt[ngood].py>>gt[ngood].pz>>
           gt[ngood].x >>gt[ngood].y >>gt[ngood].z) {
      ngood++;
      cerr<<ngood<<"\r";
      //cout<<ngood<<"\r";
      if (ngood==MAX) {
        cerr<<"Too many good tracks !\n";
        break;
      }
    }
    cout<<endl;
    if (!in.eof()) cerr<<"Read error (good_tracks_tpc) !\n";
  } else {
    cerr<<"Marking good tracks (this will take a while)...\n";
    ngood=good_tracks(gt,45000,firstev,eventn); 
    printf("Goood %d\n", ngood);
    ofstream out("good_tracks_tpc");
    if (out) {
      cout<<"File good_tracks_tpc opened\n";
      for (Int_t ngd=0; ngd<ngood; ngd++) {
        out<<gt[ngd].fEventN<<' '<<gt[ngd].lab<<' '<<gt[ngd].code<<' '<<
          gt[ngd].px<<' '<<gt[ngd].py<<' '<<gt[ngd].pz<<' '<<
          gt[ngd].x <<' '<<gt[ngd].y <<' '<<gt[ngd].z <<endl;
      }
    } else cerr<<"Can not open file (good_tracks_tpc) !\n";
    out<<flush;
    out.close();

    ofstream out2("good_tracks_tpc_par");

    if (out2) {
      //cout<<"File good_tracks_tpc opened\n";
      for (Int_t ngd=0; ngd<ngood; ngd++) {
	Float_t pt =  TMath::Sqrt(gt[ngd].px*gt[ngd].px+gt[ngd].py*gt[ngd].py);
	Float_t angle = 0;
	if (TMath::Abs(pt)>0.01) angle = gt[ngd].pz/pt;
	out2<<gt[ngd].fEventN<<"\t"<<gt[ngd].lab<<"\t"<<gt[ngd].code<<"\t"<<
	 pt<<"\t"<<angle<<"\t"<<endl;
      }
    } else cerr<<"Can not open file (good_tracks_tpc) !\n";
    out2<<flush;
    out2.close();

  }
  cerr<<"Number of good tracks : "<<ngood<<endl;
  cout<<"Number of good tracks : "<<ngood<<endl;
  if(ngood==0)return 5;
  TH1F *hp=new TH1F("hp","PHI resolution",50,-20.,20.); hp->SetFillColor(4);
  TH1F *hl=new TH1F("hl","LAMBDA resolution",50,-20,20);hl->SetFillColor(4);
  TH1F *hpt=new TH1F("hpt","Relative Pt resolution",30,-10.,10.); 
  hpt->SetFillColor(2); 
  TH1F *hmpt=new TH1F("hmpt","Relative Pt resolution (pt>4GeV/c)",30,-60,60); 
  hmpt->SetFillColor(6);

  AliTPCtrack *trk=(AliTPCtrack*)tarray.UncheckedAt(0);

  Double_t pmin=0.1*(100/0.299792458/0.2/trk->GetConvConst());
  Double_t pmax=6.0+pmin;

  TH1F *hgood=new TH1F("hgood","Good tracks",30,pmin,pmax);    
  TH1F *hfound=new TH1F("hfound","Found tracks",30,pmin,pmax);
  TH1F *hfake=new TH1F("hfake","Fake tracks",30,pmin,pmax);
  TH1F *hg=new TH1F("hg","",30,pmin,pmax); //efficiency for good tracks
  hg->SetLineColor(4); hg->SetLineWidth(2);
  TH1F *hf=new TH1F("hf","Efficiency for fake tracks",30,pmin,pmax);
  hf->SetFillColor(1); hf->SetFillStyle(3013); hf->SetLineWidth(2);

  TH1F *he =new TH1F("he","dE/dX for pions with 0.4<p<0.5 GeV/c",50,0.,100.);
  TH2F *hep=new TH2F("hep","dE/dX vs momentum",50,0.,2.,50,0.,250.);
  hep->SetMarkerStyle(8);
  hep->SetMarkerSize(0.4);

  Int_t mingood=ngood;
  Int_t track_notfound[MAX], itrack_notfound=0;
  Int_t track_fake[MAX], itrack_fake=0;
  Int_t track_multifound[MAX], track_multifound_n[MAX], itrack_multifound=0;
  while (ngood--) {
    Int_t lab=gt[ngood].lab,tlab=-1;
    Float_t ptg=
      TMath::Sqrt(gt[ngood].px*gt[ngood].px + gt[ngood].py*gt[ngood].py);

    if (ptg<pmin) continue;  //  Should be 1e-33 instead pmin?

    hgood->Fill(ptg);
    Int_t ievent = gt[ngood].fEventN;

    AliTPCtrack *track=0;
    for (i=eventptr[ievent]; i<eventptr[ievent+1]; i++) {

      track=(AliTPCtrack*)tarray.UncheckedAt(i);
      tlab=track->GetLabel();
      if (lab==TMath::Abs(tlab)) break;
    }
    if (i==eventptr[ievent+1]) {
      track_notfound[itrack_notfound++]=lab;
      continue;
    }
      
    Int_t micount=0;
    Int_t mi;
    AliTPCtrack * mitrack;
    for (mi=eventptr[ievent]; mi<eventptr[ievent+1]; mi++) {

      mitrack=(AliTPCtrack*)tarray.UncheckedAt(mi);          
      if (lab==TMath::Abs(mitrack->GetLabel())) micount++;
    }
    if (micount>1) {
      track_multifound[itrack_multifound]=lab;
      track_multifound_n[itrack_multifound]=micount;
      itrack_multifound++;
    }
      
    //
    Double_t xk=gt[ngood].x;
    if (!track) continue;
    //    printf("Track =%p\n",track);
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
    } 
    //else 
    {
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
  if (ng!=0) cout<<"\n\nIntegral efficiency is about "<<nf/ng*100.<<" %\n";
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

  gBenchmark->Stop("AliTPCComparison2");
  gBenchmark->Show("AliTPCComparison2");


  return 0;
}


Int_t good_tracks(GoodTrackTPC *gt, Int_t max, Int_t firstev, Int_t eventn) {
  /// eventn  - number of events in file

  TFile *file=TFile::Open("galice.root");
  if (!file->IsOpen()) {cerr<<"Can't open galice.root !\n"; exit(4);}
  //  delete gAlice; gAlice = 0;
  if (!(gAlice=(AliRun*)file->Get("gAlice"))) {
    cerr<<"gAlice have not been found on galice.root !\n";
    exit(5);
  }

  TFile *fdigit = TFile::Open("digits.root");
  file->cd();

  AliTPC *TPC=(AliTPC*)gAlice->GetDetector("TPC");
  Int_t ver = TPC->IsVersion(); 
  cerr<<"TPC version "<<ver<<" has been found !\n";

  AliTPCParam *digp=(AliTPCParam*)file->Get("75x40_100x60");
  if(digp){
    cerr<<"2 pad-lenght geom hits with 3 pad-length geom digits...\n";
    delete digp;
    digp = new AliTPCParamSR();
   }
   else
   {
     digp =(AliTPCParam*)gDirectory->Get("75x40_100x60_150x60");
   }

  if (!digp) { cerr<<"TPC parameters have not been found !\n"; exit(6); }
  TPC->SetParam(digp);

  Int_t nrow_up=digp->GetNRowUp();
  Int_t nrows=digp->GetNRowLow()+nrow_up;
  Int_t zero=digp->GetZeroSup();
  Int_t gap=Int_t(0.125*nrows), shift=Int_t(0.5*gap);
  Int_t good_number=Int_t(0.4*nrows);

  Int_t nt=0;  //reset counter
  char treeName[100];  //declare event identifier
       
  for (Int_t event=firstev;event<eventn;event++){
    Int_t np=gAlice->GetEvent(event);
    Int_t nopr = FindPrimaries(np);
    cout<<"function good_tracks -- event "<<event<<", Charged primaries:"<<nopr<<"\n";

    Int_t *good=new Int_t[np];
    for (Int_t ii=0; ii<np; ii++) good[ii]=0;
     
     
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
        ca->ConnectTree("TreeC_TPC_0");
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
        delete ca;
        cf->Close();
      }
      break;
    case 2:
      {
        sprintf(treeName,"TreeD_75x40_100x60_150x60_%d",event);  
        TD=(TTree*)fdigit->Get(treeName); // To be revised according to
                                          // NewIO schema M.Kowalski
        TD->GetBranch("Segment")->SetAddress(&digits);
        count = new Int_t[np];
        Int_t i;
        for (i=0; i<np; i++) count[i]=0;
        sectors_by_rows=(Int_t)TD->GetEntries();
        for (i=0; i<sectors_by_rows; i++) {
          if (!TD->GetEvent(i)) continue;
          Int_t sec,row;
          digp->AdjustSectorRow(digits->GetID(),sec,row);
          digits->First();
	  digits->ExpandTrackBuffer();
          do { //Many thanks to J.Chudoba who noticed this
            Int_t it=digits->CurrentRow(), ip=digits->CurrentColumn();
	    //            Short_t dig = digits->GetDigit(it,ip);
            Short_t dig = digits->CurrentDigit();

            Int_t idx0=digits->GetTrackIDFast(it,ip,0)-2; 
            Int_t idx1=digits->GetTrackIDFast(it,ip,1)-2;
            Int_t idx2=digits->GetTrackIDFast(it,ip,2)-2;
            if (idx0>=0 && dig>=zero && idx0<np ) count[idx0]+=1;   //background events - higher track ID's
            if (idx1>=0 && dig>=zero && idx1<np ) count[idx1]+=1;
            if (idx2>=0 && dig>=zero && idx2<np ) count[idx2]+=1;
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
        delete TD; //Thanks to Mariana Bondila
      }
      break;
    default:
      cerr<<"Invalid TPC version !\n";
      file->Close();
      exit(7);
    }
  
    /** select tracks which are "good" enough **/
    //printf("\t %d \n",np);
    for (Int_t i=0; i<np; i++) {
      if ((good[i]&0x5000) != 0x5000)
        if ((good[i]&0x2800) != 0x2800) continue;
      if ((good[i]&0x7FF ) < good_number) continue;

      TParticle *p = (TParticle*)gAlice->Particle(i);

      if (p->Pt()<0.100) continue;
      if (TMath::Abs(p->Pz()/p->Pt())>0.999) continue;
      //  THIS IS GOOD IF YOU WANT SECONDARIES
        Int_t j=p->GetFirstMother();
        if (j>=0) {
        TParticle *pp = (TParticle*)gAlice->Particle(j);
        //if (!(pp->TestBit(kPrimaryCharged))) continue; //only one decay is allowed
        }
      
	if(!(p->TestBit(kPrimaryCharged)))continue; // only primaries
	//	printf("1");
  
      gt[nt].fEventN=event;
      gt[nt].lab=i;
      gt[nt].code=p->GetPdgCode();
      gt[nt].px=0.; gt[nt].py=0.; gt[nt].pz=0.;
      gt[nt].x=0.; gt[nt].y=0.; gt[nt].z=0.;
      nt++;
      if (nt==max) {cerr<<"Too many good tracks !\n"; break;}
      cerr<<np-i<<"        \r";    
    }

    /** check if there is also information at the entrance of the TPC **/
    TTree *TH=gAlice->TreeH(); np=(Int_t)TH->GetEntries();
    cout<<endl;
    for (Int_t i=0; i<np; i++) {
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
        for (j=0; j<nt; j++) {if (gt[j].fEventN==event && gt[j].lab==lab) break;}
        if (j==nt) continue;         
	//printf("1-");
        // (px,py,pz) - in global coordinate system, (x,y,z) - in local !
        gt[j].px=px; gt[j].py=py; gt[j].pz=pz;
        Float_t cs,sn; digp->AdjustCosSin(phit->fSector,cs,sn);
        gt[j].x = x*cs + y*sn;
        gt[j].y =-x*sn + y*cs;
        gt[j].z = z;
      }
      cerr<<np-i<<"        \r";
    }
    //printf("\n%d\n",nt);
    cout<<endl;
    delete[] good;
  }             // ///         loop on events
  delete gAlice; gAlice=0;

  file->Close();
  gBenchmark->Stop("AliTPCComparison2");
  gBenchmark->Show("AliTPCComparison2");   
  return nt;

}

Int_t FindPrimaries(Int_t nparticles){

  // cuts:
  Double_t vertcut = 0.001;
  Double_t decacut = 3.;
  Double_t timecut = 0.;
  Int_t nprch1=0;
  TParticle * part = gAlice->Particle(0);
  Double_t xori = part->Vx();
  Double_t yori = part->Vy();
  Double_t zori = part->Vz();
  for(Int_t iprim = 0; iprim<nparticles; iprim++){   //loop on  tracks

    part = gAlice->Particle(iprim);
    Double_t ptot=TMath::Sqrt(part->Px()*part->Px()+part->Py()*part->Py()+part->Pz()*part->Pz());
    if (ptot<0.01) continue;
    char *xxx=strstr(part->GetName(),"XXX");
    if(xxx)continue;

    TParticlePDG *ppdg = part->GetPDG();
    if(TMath::Abs(ppdg->Charge())!=3)continue;  // only charged (no quarks)

    Double_t dist=TMath::Sqrt((part->Vx()-xori)*(part->Vx()-xori)+(part->Vy()-yori)*(part->Vy()-yori)+(part->Vz()-zori)*(part->Vz()-zori));
    if(dist>vertcut)continue;  // cut on the vertex

    if(part->T()>timecut)continue;

    //Double_t ptot=TMath::Sqrt(part->Px()*part->Px()+part->Py()*part->Py()+part->Pz()*part->Pz());
    if(ptot==(TMath::Abs(part->Pz())))continue; // no beam particles

    Bool_t prmch = kTRUE;   // candidate primary track
    Int_t fidau=part->GetFirstDaughter();  // cut on daughters
    Int_t lasdau=0;
    Int_t ndau=0;
    if(fidau>=0){
      lasdau=part->GetLastDaughter();
      ndau=lasdau-fidau+1;
    }
    if(ndau>0){
      for(Int_t j=fidau;j<=lasdau;j++){
        TParticle *dau=gAlice->Particle(j);
        Double_t distd=TMath::Sqrt((dau->Vx()-xori)*(dau->Vx()-xori)+(dau->Vy()-yori)*(dau->Vy()-yori)+(dau->Vz()-zori)*(dau->Vz()-zori));
        if(distd<decacut)prmch=kFALSE;  // eliminate if the decay is near the vertex
      }
    }

    if(prmch){
      nprch1++;
      part->SetBit(kPrimaryCharged);
    }
  }

  return nprch1;
}

