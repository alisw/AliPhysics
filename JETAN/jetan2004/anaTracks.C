// $Id$

#ifndef __CINT__
#include <stdio.h>
#include <iostream.h>
#include <time.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TFile.h>
#include <TMath.h>
#include <TChain.h>
#include <TStopwatch.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <AliStack.h>
#include <AliRunLoader.h>
#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliJetParticle.h>
#include <AliJetParticlesReader.h>
#include <AliJetParticlesReaderKine.h>
#include <AliJetParticlesReaderKineGoodTPC.h>
#include <AliJetParticlesReaderESD.h>
#include <AliJetParticlesReaderHLT.h>
#include <AliJetEventParticles.h>
#endif

#define LABEFF
#define MINHITS 0

void Sort(Int_t na,const Int_t *a,Int_t *index);
void anaTracks(Int_t esd1,Int_t esd2,Float_t qcut,Int_t nMaxEvents=-1);
void anaTracks(Char_t *filename1,Char_t *filename2,Float_t qcut,Int_t nMaxEvents=-1);
void anaTracks(TObjArray *fnames1,TObjArray *fnames2,Float_t qcut,Int_t nMaxEvents=-1);
void anaEtTracks(Int_t esd1,Int_t esd2,Float_t qcut,
		 Int_t etfrom=1, Int_t etto=16, Int_t nMaxEvents=-1);

void anaTracks(Int_t esd1,Int_t esd2,Float_t qcut,Int_t nMaxEvents)
{
  const Int_t nes[5]={50,100,150,200,250};
  const Char_t *dir="/data/aliroot/jets";
  TObjArray *f1=new TObjArray(1);
  TObjArray *f2=new TObjArray(1);
  Char_t filename1[1024];
  Char_t filename2[1024];

  for(Int_t i=0;i<5;i++){
    sprintf(filename1,"%s/%d/aliev-type%d.root",dir,nes[i],esd1);
    f1->Add(new TObjString(filename1));
    sprintf(filename2,"%s/%d/aliev-type%d.root",dir,nes[i],esd2);
    f2->Add(new TObjString(filename2));
  }
  anaTracks(f1,f2,qcut,nMaxEvents);
}

void anaEtTracks(Int_t esd1,Int_t esd2,Float_t qcut,Int_t etfrom, Int_t etto, Int_t nMaxEvents)
{
  const Char_t *dir="/data/aliroot/scan";
  TObjArray *f1=new TObjArray(1);
  TObjArray *f2=new TObjArray(1);
  Char_t filename1[1024];
  Char_t filename2[1024];
  for(Int_t i=etfrom;i<=etto;i++){
    sprintf(filename1,"%s/et%d/aliev-type%d.root",dir,i,esd1);
    f1->Add(new TObjString(filename1));
    sprintf(filename2,"%s/et%d/aliev-type%d.root",dir,i,esd2);
    f2->Add(new TObjString(filename2));
  }
  anaTracks(f1,f2,qcut,nMaxEvents);
}

void anaTracks(Char_t *filename1,Char_t *filename2,Float_t qcut,Int_t nMaxEvents)
/*
  filename1 = monte events
  filename2 = tracked events
*/
{
  TObjArray *f1=new TObjArray(1);
  f1->Add(new TObjString(filename1));
  TObjArray *f2=new TObjArray(1);
  f2->Add(new TObjString(filename2));

  anaTracks(f1,f2,qcut,nMaxEvents);
}

void anaTracks(TObjArray *fnames1,TObjArray *fnames2,Float_t qcut,Int_t nMaxEvents)
/*
  filename1 = monte events
  filename2 = tracked events
*/
{
  gRandom->SetSeed(0); 

 //connect to events
  TChain *theTree1 = new TChain("AJEPtree");
  for(Int_t i=0;i<fnames1->GetEntries();i++){
    theTree1->Add(((TObjString*)fnames1->At(i))->String());
  }
  AliJetEventParticles *ev1=new AliJetEventParticles();
  theTree1->SetBranchAddress("particles",&ev1);

  //connect to events
  TChain *theTree2 = new TChain("AJEPtree");
  for(Int_t i=0;i<fnames2->GetEntries();i++){
    theTree2->Add(((TObjString*)fnames2->At(i))->String());
  }
  AliJetEventParticles *ev2=new AliJetEventParticles();
  theTree2->SetBranchAddress("particles",&ev2);

  Int_t treeentries1=(Int_t)theTree1->GetEntries();
  Int_t treeentries2=(Int_t)theTree2->GetEntries();

  if((nMaxEvents<0) || (nMaxEvents>treeentries1))
    nMaxEvents=treeentries1;

  cout << "Found " << treeentries1 << ", but use " << nMaxEvents << " in chain1" << endl;
  cout << "Found " << treeentries2 << ", but use " << nMaxEvents << " in chain2" << endl;
  
  const Int_t nmax=25000;
  Int_t *ind2=new Int_t[nmax];
  const Float_t ptcut=0.5;
  const Int_t nbinsx=20;
  const Float_t binsx[nbinsx]={1,2,4,6,8,10,12.5,15,17.5,20,25,30,40,45,50,60,70,80,90,100};

  //hinvp
  TH1F *hInvQ2=new TH1F("hInvQ2","",1000,0,1.5);
  TH1F *hInvQ2lab=new TH1F("hInvQ2lab","",1000,0,1.5);
  TH1F *hEffPt1=new TH1F("hEffPt1","",nbinsx-1,binsx);
  TH1F *hEffPt2=new TH1F("hEffPt2","",nbinsx-1,binsx);
  TH1F *hEffEta1=new TH1F("hEffEta1","",40,-1,1);
  TH1F *hEffEta2=new TH1F("hEffEta2","",40,-1,1);
  TH1F *hEffPhi1=new TH1F("hEffPhi1","",62,0,6.2);
  TH1F *hEffPhi2=new TH1F("hEffPhi2","",62,0,6.2);
  TH1F *hFakePt=new TH1F("hFakePt","",nbinsx-1,binsx);
  TH1F *hFakeEta=new TH1F("hFakeEta","",40,-1,1);
  TH1F *hFakePhi=new TH1F("hFakePhi","",62,0,6.2);
  TH1F *hResPt=new TH1F("hResPt","",60,-15,15);
  TH1F *hResEta=new TH1F("hResEta","",200,-0.05,0.05);
  TH1F *hResPhi=new TH1F("hResPhi","",200,-0.05,0.05);
  TProfile *hpResPt=new TProfile("hpResPt","",nbinsx-1,binsx);
  TProfile *hpResEta=new TProfile("hpResEta","",nbinsx-1,binsx);
  TProfile *hpResPhi=new TProfile("hpResPhi","",nbinsx-1,binsx);
  TH1F **hResPtAll=new TH1F*[nbinsx];
  TH1F **hResEtaAll=new TH1F*[nbinsx];
  TH1F **hResPhiAll=new TH1F*[nbinsx];
  for(Int_t i=0;i<nbinsx;i++){
    Char_t fn[100];
    sprintf(fn,"hResPt%d",i);
    hResPtAll[i]=new TH1F(fn,"",60,-15,15);
    sprintf(fn,"hResEta%d",i);
    hResEtaAll[i]=new TH1F(fn,"",200,-0.05,0.05);
    sprintf(fn,"hResPhi%d",i);
    hResPhiAll[i]=new TH1F(fn,"",200,-0.05,0.05);
  }
  //=========================================================================
  // start the event loop
  //=========================================================================
  Int_t nEvent1 = -1;
  Int_t nEvent2 = -1;
  Int_t ntracks1 = 0;
  Int_t ntracks2 = 0;
  Int_t nfakes2 = 0;
  Int_t nlokfakes2 = 0;
  Int_t nlabs = 0;
  while(nEvent1<nMaxEvents && nEvent2<nMaxEvents){
    nEvent1++;
    nEvent2++;
    if((nEvent1 % 100) == 99) {
      cout << "Reading event " << nEvent1 << " " << nlabs << endl;
      nlabs=0;
    }

    ev1->Clear();
    ev2->Clear();
    theTree1->GetEvent(nEvent1);
    theTree2->GetEvent(nEvent2);

    if(ev1->GetEventNr()!=ev2->GetEventNr()){
      cerr << "Need to skip event: " << nEvent1 << " != " << nEvent2 
	   << " (" << ev1->GetEventNr() << " " << ev2->GetEventNr() <<") " << endl;
      if(ev1->GetEventNr()>ev2->GetEventNr()) nEvent2++;
      else nEvent1++;
      continue;
    }

    //ev1->Print();
    //ev2->Print();
    const TClonesArray *p1=ev1->GetParticles();
    if(!p1){
      continue;
    }
    const TClonesArray *p2=ev2->GetParticles();
    if(!p2){
      continue;
    }
    const Int_t n1=p1->GetEntriesFast();
    const Int_t n2=p2->GetEntriesFast();
    if(n1==0 || n2==0){
      continue;
    }
    //clean found particle table
    for(Int_t j=0;j<n2;j++) {
      ind2[j]=1;
    }

    //Compute Q^2 with/without labels
    for(Int_t i=0;i<n1;i++) {
      AliJetParticle *p=(AliJetParticle*)p1->At(i);
      if(p->Pt()<ptcut) continue;
      ntracks1++;
      Float_t px=p->Px();
      Float_t py=p->Py();
      Float_t pz=p->Pz();
      Int_t lab=p->GetLabel();
      if(lab==0) lab=p->GetUID(); //workaround for bug
      for(Int_t j=0;j<n2;j++) {
	AliJetParticle *q=(AliJetParticle*)p2->At(j);
	if(q->Pt()<ptcut) continue;
	if(q->GetNhts()<MINHITS) continue;
	Float_t qx=q->Px();
	Float_t qy=q->Py();
	Float_t qz=q->Pz();
	Float_t Q2=(px-qx)*(px-qx)+(py-qy)*(py-qy)+(pz-qz)*(pz-qz);
	Q2=TMath::Sqrt(Q2)/p->P();
	hInvQ2->Fill(Q2);
	if(q->GetLabel()!=lab) continue;
	hInvQ2lab->Fill(Q2);
      } 
    }

    for(Int_t i=0;i<n1;i++) {
      AliJetParticle *p=(AliJetParticle*)p1->At(i);
      if(p->Pt()<ptcut) continue;
      Float_t px=p->Px();
      Float_t py=p->Py();
      Float_t pz=p->Pz();
      Int_t lab=p->GetLabel();
      if(lab==0) lab=p->GetUID(); //workaround for bug

      Float_t minQ2=1e5;
      Int_t j2=-1,lab2=-1;
      for(Int_t j=0;j<n2;j++) {
	if(ind2[j]!=1) continue; //particle was already found 
	AliJetParticle *q=(AliJetParticle*)p2->At(j);
	if(q->Pt()<ptcut) continue;
	if(q->GetNhts()<MINHITS) continue;
	Float_t qx=q->Px();
	Float_t qy=q->Py();
	Float_t qz=q->Pz();
	Float_t Q2=(px-qx)*(px-qx)+(py-qy)*(py-qy)+(pz-qz)*(pz-qz);
	Q2=TMath::Sqrt(Q2)/p->P();
#ifdef LABEFF
	if(TMath::Abs(q->GetLabel())==lab){
	  if(j2<0){
	    lab2=q->GetLabel();
	    j2=j;
	    minQ2=Q2;
	  } else { //take particle with smaller Q2
	    if(Q2<minQ2){
	      ind2[j2]=2;//mark old found as fake
	      lab2=q->GetLabel();
	      j2=j;
	      minQ2=Q2;
	    } else ind2[j]=2;//mark new found as fake
	  }
	}
#else
	if(Q2<minQ2){ //take particle with smaller Q2
	  minQ2=Q2;
	  j2=j;
	  lab2=q->GetLabel();
	}
#endif
      }
      
      hEffPt1->Fill(p->Pt());
      hEffEta1->Fill(p->Eta());
      hEffPhi1->Fill(p->Phi());

      if(j2>=0){ //take found particle
#ifdef LABEFF
      if(lab!=lab2){
	nlabs++;
	continue;
      }
#else
        if(minQ2>qcut) continue;
#endif
	ntracks2++;
	ind2[j2]=0; //mark as found
	AliJetParticle *q=(AliJetParticle*)p2->At(j2);
	hEffPt2->Fill(q->Pt());
	hEffEta2->Fill(q->Eta());
	hEffPhi2->Fill(q->Phi());
	Float_t ptdiff=p->Pt()-q->Pt();
	Float_t phidiff=p->Phi()-q->Phi();
	Float_t etadiff=p->Eta()-q->Eta();
	hResPt->Fill(ptdiff);
	hResEta->Fill(etadiff);
	hResPhi->Fill(phidiff);
	//cout << ptdiff << etadiff << phidiff << endl;
	hpResPt->Fill(p->Pt(),TMath::Abs(ptdiff)/p->Pt());
	hpResEta->Fill(p->Pt(),TMath::Abs(etadiff));
	hpResPhi->Fill(p->Pt(),TMath::Abs(phidiff));
	Int_t index=hEffPt1->FindBin(p->Pt());
	hResPtAll[index]->Fill(ptdiff);
	hResEtaAll[index]->Fill(etadiff);
	hResPhiAll[index]->Fill(phidiff);
      } 

    } //event loop
    nlokfakes2=0;
    for(Int_t j=0;j<n2;j++) {
      if(!ind2[j]) continue;
      AliJetParticle *q=(AliJetParticle*)p2->At(j);
      if (q->Pt()<ptcut) continue;
      //cout << ind2[j] << " " << q->GetLabel() << endl;
      if(q->GetNhts()<MINHITS) continue;
      hFakePt->Fill(q->Pt());
      hFakeEta->Fill(q->Eta());
      hFakePhi->Fill(q->Phi());
      nlokfakes2++;
    }
    nfakes2+=nlokfakes2;
    if(nlokfakes2/n2>0.1) cout << nEvent1 << " has more than 10% fakes: " << nlokfakes2 << endl;
  } //end of nev loop

  Float_t ffakes=nfakes2*100./ntracks1;
  cout << "Good tracks:  " << ntracks1 << endl;
  cout << "Found tracks: " << ntracks2 << endl;
  cout << "Fake tracks:  " << nfakes2  << endl;
  cout << "Eff:  " << ntracks2*100/ntracks1 << endl;
  cout << "Fake: " << nfakes2*100/ntracks1 << endl;
  TF1 *f1 = new TF1("f1","gaus");
  hResPt->Sumw2();
  Float_t m = hResPt->GetMean();
  Float_t s = hResPt->GetRMS();
  f1->SetParameter(1,m);
  f1->SetParameter(2,s);
  cout << "Res: " << m << " " << s << endl;

  hResPt->Fit(f1,"IN");
  m = f1->GetParameter(1);
  s = f1->GetParameter(2);
  cout << "Res: " << m << " " << s << endl;

  TCanvas *c1=new TCanvas("cInvQ2","",600,400);
  hInvQ2->Sumw2();
  hInvQ2->Draw();
  c1->Update();
  Float_t qq=0,ival;
  while(1){
    ival=hInvQ2->Integral(hInvQ2->FindBin(0),hInvQ2->FindBin(qq));
    if(ival>ntracks1) break;
    qq+=0.001;
  }
  cout << "Qcut: " << qq << " " << ival << " " << ntracks1 << endl;
  TCanvas *c1b=new TCanvas("cInvQ2lab","",600,400);
  hInvQ2lab->Sumw2();
  hInvQ2lab->Draw();
  c1b->Update();

  hEffPt2->Sumw2();
  hEffPt1->Sumw2();
  hEffPt2->Divide(hEffPt1);

  hEffEta2->Sumw2();
  hEffEta1->Sumw2();
  hEffEta2->Divide(hEffEta1);

  hEffPhi2->Sumw2();
  hEffPhi1->Sumw2();
  hEffPhi2->Divide(hEffPhi1);

  hFakePt->Sumw2();
  hFakePt->Divide(hEffPt1);
  hFakeEta->Sumw2();
  hFakeEta->Divide(hEffEta1);
  hFakePhi->Sumw2();
  hFakePhi->Divide(hEffPhi1);

  TGraphErrors *gpts=new TGraphErrors();
  TGraphErrors *getas=new TGraphErrors();
  TGraphErrors *gphis=new TGraphErrors();

//#define fit
  for(Int_t i=0;i<nbinsx;i++){
    hResPtAll[i]->Sumw2();
    m = hResPtAll[i]->GetMean();
    s = hResPtAll[i]->GetRMS();
#ifdef fit
    f1->SetParameter(1,m);
    f1->SetParameter(2,s);
    hResPtAll[i]->Fit(f1,"IN");
    m = f1->GetParameter(1);
    s = f1->GetParameter(2);
    //gpts->SetPointError(i,0,f1->GetParError(2));
#endif
    gpts->SetPoint(i,binsx[i],s/binsx[i]*100);

    hResEtaAll[i]->Sumw2();
    m = hResEtaAll[i]->GetMean();
    s = hResEtaAll[i]->GetRMS();
#ifdef fit
    f1->SetParameter(1,m);
    f1->SetParameter(2,s);
    hResEtaAll[i]->Fit(f1,"QIN");
    m = f1->GetParameter(1);
    s = f1->GetParameter(2);
#endif
    getas->SetPoint(i,binsx[i],s);

    hResPhiAll[i]->Sumw2();
#ifdef fit
    m = hResPhiAll[i]->GetMean();
    s = hResPhiAll[i]->GetRMS();
    f1->SetParameter(1,m);
    f1->SetParameter(2,s);
    hResPhiAll[i]->Fit(f1,"QIN");
    m = f1->GetParameter(1);
    s = f1->GetParameter(2);
#endif
    gphis->SetPoint(i,binsx[i],s);
  }

#if 0
  TCanvas *c2=new TCanvas("cEffPt","",600,400);
  hEffPt2->Draw();
  c2->Update();
  TCanvas *c3=new TCanvas("cEffEta","",600,400);
  hEffEta2->Draw();
  c3->Update();
  TCanvas *c4=new TCanvas("cEffPhi","",600,400);
  hEffPhi2->Draw();
  c4->Update();

  TCanvas *c5=new TCanvas("cResPt","",600,400);
  hpResPt->Draw();
  gpts->SetMarkerSize(1);
  gpts->SetMarkerStyle(20);
  //gpts->Draw("p axis");
  c5->Update();
  TCanvas *c6=new TCanvas("cResPhi","",600,400);
  hpResPhi->Draw();
  gphis->SetMarkerSize(1);
  gphis->SetMarkerStyle(20);
  //gphis->Draw("p axis");
  c6->Update();
  TCanvas *c7=new TCanvas("cResEta","",600,400);
  hpResEta->Draw();
  getas->SetMarkerSize(1);
  getas->SetMarkerStyle(20);
  //getas->Draw("p axis");
  c7->Update();
#else
  TCanvas *c8=new TCanvas("call","",1200,800);
  c8->Divide(3,2);
  c8->cd(1);
  hEffPt2->Draw();
  hFakePt->Draw("same");
  c8->cd(2);
  hEffEta2->Draw();
  hFakeEta->Draw("same");
  c8->cd(3);
  hEffPhi2->Draw();
  hFakePhi->Draw("same");
  c8->cd(4);
  hpResPt->Draw();
  c8->cd(5);
  hpResEta->Draw();
  c8->cd(6);
  hpResPhi->Draw();
#endif

  TFile *f=new TFile("tracks.root","recreate");
  hInvQ2->Write();
  hInvQ2lab->Write();
  hEffPt1->Write();
  hEffPt2->Write();
  hEffEta1->Write();
  hEffEta2->Write();
  hEffPhi1->Write();
  hEffPhi2->Write();
  hResPt->Write();
  hResEta->Write();
  hResPhi->Write();
  hFakePt->Write();
  hFakePhi->Write();
  hFakeEta->Write();
  gpts->Write("gtps");
  getas->Write("gtps");
  gphis->Write("gtps");
  f->mkdir("all");
  f->cd("all");
  for(Int_t i=0;i<nbinsx;i++){
    hResPtAll[i]->Write();
    hResEtaAll[i]->Write();
    hResPhiAll[i]->Write();
  }
  f->Close();
  delete f;
  delete[] ind2;
  delete ev1;
  delete ev2;
  delete theTree1;
  delete theTree2;
}

//_____________________________________________________________________________


void Sort(Int_t na, const Int_t *a, Int_t *index)
{
   if (na <= 0) return;
   Int_t *localArr1 = new Int_t[na];
   Int_t *localArr2 = new Int_t[na];
   Int_t iEl;
   Int_t iEl2;

   for(iEl = 0; iEl < na; iEl++) {
      localArr1[iEl] = a[iEl];
      localArr2[iEl] = iEl;
      //cout << localArr1[iEl] << " " << localArr2[iEl] << endl;
   }

   for (iEl = 0; iEl < na; iEl++) {
      for (iEl2 = na-1; iEl2 > iEl; iEl2--) {
	//cout << iEl << " " << iEl2 << " " << localArr1[iEl2-1] << " " << localArr1[iEl2] <<  endl;

         if (localArr1[iEl2-1] < localArr1[iEl2]) {
            Int_t tmp         = localArr1[iEl2-1];
            localArr1[iEl2-1] = localArr1[iEl2];
            localArr1[iEl2]   = tmp;

            Int_t tmp2        = localArr2[iEl2-1];
            localArr2[iEl2-1] = localArr2[iEl2];
            localArr2[iEl2]   = tmp2;
	    //cout << "Swapped: " << localArr1[iEl2-1] << " " << localArr2[iEl2-1] <<  endl;
	    //cout << "Swapped: " << localArr1[iEl2] << " " << localArr2[iEl2] <<  endl;
         }
	 //cout << iEl << " " << iEl2 << " " << localArr1[iEl2-1] << " " << localArr1[iEl2] <<  endl;
      }
   }

   for (iEl = 0; iEl < na; iEl++) {
      index[iEl] = localArr2[iEl];
      //cout << a[iEl] << " " << index[iEl] <<  endl;
   }
   delete [] localArr2;
   delete [] localArr1;
}
