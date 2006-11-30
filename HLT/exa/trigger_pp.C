// $Id$

/**
   trigger_pp -> Calculate efficiency of trigger and fake tracks in pp pileup events
   plot_pp -> Plot vertex resolution
   plot_pp_vert -> Plot vertex resolution
   If you are looking for the old trigger macro search for 
   OLD_TRIGGER_MACRO.

   Author: C. Loizides
*/


/* ---trigger_pp---
   Calculate the efficiency of found tracks in pileup event 
   and found fake tracks after applying different DCA cuts to 
   the assumed vertex of (0,0,0). Tracks are stored in the ntuple
   root file (see fill_pp.C). Good tracks can either be used from
   the single pp event or the good_tracks_tpc event (the macro
   will detect which ntuple file you specified).

   Script: Track single pp events and store in ntuple:
--------------------------------------------------------
#!/bin/bash

dir=/mnt/data/loizides/alidata/head
wdir=$dir/pp/getrackted-0-375/
mdir=~/work/l3/level3code/exa/
sfile=~/pp-tracks.root
rfile=$dir/PythiaPP_1000Events_TPConly_ev0-375_digits.root

cd $dir
ln -s -f $rfile digitfile.root
cd $mdir;

for i in `seq 0 365`; do
 echo $i
 trigev=$i
 cd $mdir
 aliroot -b -q runtracker_pp.C\($trigev,0,\"$rfile\",\"$wdir\"\)
 aliroot -b -q fill_pp.C\($trigev,\"$wdir\",\"$sfile\"\)
done
--------------------------------------------------------

   Script: Track pileup events and store in ntuple:
-----------------------------------------------------
#!/bin/bash

dir=/mnt/data/loizides/alidata/head
wdir=/mnt/data/loizides/alidata/head/pp/pileup-25
mdir=~/l3/level3code/exa
sfile=~/pileup25-tracks.root

cd $mdir;

for i in `seq 0 100`; do
 echo $i
 pdir=$wdir/pileup-$i/
 cd $pdir
 trigev=`ls digits_*_0_-1.raw | cut -d_ -f 2`
 cd $mdir
 aliroot -b -q runtracker_pp.C\($trigev,\"$pdir\",0,\"$pdir\"\)
 aliroot -b -q fill_pp.C\($trigev,\"$pdir\",\"$sfile\",$i\)
done
-----------------------------------------------------

  Parameter to trigger_pp are:
  Pileupfile is the ntuple root file containing the information of the pileup event
  Ppfile gives the ntuple root file containing the reference information of the single pp event (can either be good_particles or single tracked pp events)
  Evi and Eve give start and end event number to process.
*/

void trigger_pp(Char_t *pileupfile, Char_t *ppfile, Int_t evi=0, Int_t eve=100)
{

  //-----------------------------------------------------------
  const Char_t *LTEXT_="for 25 piles, min. 2 trigger tracks";  
  const Int_t MINTRACKS_=2; //min. tracks of the triggered 
  const Int_t N_=100;       //number of data points
  const Float_t WIDTH_=1; //the width of the xbin
  //-----------------------------------------------------------


  //"pt:phi:eta:xvert:yvert:zvert:imp:nhits:px:py:pz:event:mc")
  Float_t pt1,phi1,eta1,xvert1,yvert1,zvert1,imp1,nhits1,px1,py1,pz1,event1,mc1; //pileup
  Float_t pt2,phi2,eta2,xvert2,yvert2,zvert2,imp2,nhits2,px2,py2,pz2,event2,mc2; //pp

  Int_t countgoodevents=0;
  Int_t nent1,nent2;
  Char_t dummy[1024];

  Int_t tnorm;
  Int_t fnorm;
  Int_t tcounts[N_];
  Int_t fcounts[N_];

  Float_t tallcounts[N_];
  Float_t fallcounts[N_];
  Int_t totnorm=0;

  Float_t thcounts[N_];
  Float_t fhcounts[N_];
  Float_t thmeancounts[N_];
  Float_t fhmeancounts[N_];
  Float_t thsigmacounts[N_];
  Float_t fhsigmacounts[N_];
  Float_t xcut[N_];
  for(Int_t k=0;k<N_;k++){
    tallcounts[k]=0;
    fallcounts[k]=0;
    thmeancounts[k]=0;
    fhmeancounts[k]=0;
    thsigmacounts[k]=0;
    fhsigmacounts[k]=0;
    xcut[k]=(k+1)*WIDTH_; 
  }

  TFile file1=TFile(pileupfile,"READ");
  TFile file2=TFile(ppfile,"READ");
  if(!file1.IsOpen() || !file2.IsOpen()) return;

  for(Int_t i=evi;i<eve;i++){ //loop over pileup-events

    sprintf(dummy,"good-pptracks-%d",i);
    TNtuple *ntuppel1 = (TNtuple*)file1.Get(dummy);
    Float_t nhits1=0,pdg1=0,sector1=0;
    if(ntuppel1){
      ntuppel1->SetBranchAddress("mc",&mc1);
      ntuppel1->SetBranchAddress("pdg",&pdg1);
      ntuppel1->SetBranchAddress("px",&px1);
      ntuppel1->SetBranchAddress("py",&py1);
      ntuppel1->SetBranchAddress("pz",&pz1);
      ntuppel1->SetBranchAddress("xvert",&xvert1);
      ntuppel1->SetBranchAddress("yvert",&yvert1);
      ntuppel1->SetBranchAddress("zvert",&zvert1);
      ntuppel1->SetBranchAddress("nhits",&nhits1);
      ntuppel1->SetBranchAddress("sector",&sector1);
      event1=i;
   } else {
      sprintf(dummy,"pptracks-%d",i);
      TNtuple *ntuppel1 = (TNtuple*)file1.Get(dummy);
      if(!ntuppel1) continue;
      ntuppel1->SetBranchAddress("pt",&pt1);
      ntuppel1->SetBranchAddress("phi",&phi1);
      ntuppel1->SetBranchAddress("eta",&eta1);
      ntuppel1->SetBranchAddress("xvert",&xvert1);
      ntuppel1->SetBranchAddress("yvert",&yvert1);
      ntuppel1->SetBranchAddress("zvert",&zvert1);
      ntuppel1->SetBranchAddress("imp",&imp1);
      ntuppel1->SetBranchAddress("px",&px1);
      ntuppel1->SetBranchAddress("py",&py1);
      ntuppel1->SetBranchAddress("pz",&pz1);
      ntuppel1->SetBranchAddress("event",&event1);
      ntuppel1->SetBranchAddress("mc",&mc1);
    }

    tnorm=0;
    fnorm=0;
    for(Int_t j=0;j<N_;j++){
      tcounts[j]=0;
      fcounts[j]=0;
      thcounts[j]=0;
      fhcounts[j]=0;
    }

    nent1=ntuppel1->GetEntries();
    for(Int_t j=0;j<nent1;j++){
      ntuppel1->GetEvent(j);
      if(nhits1!=0){ // file1 is good-particle file
	pt1=TMath::Sqrt(px1*px1+py1*py1);
	Float_t p=TMath::Sqrt(px1*px1+py1*py1+pz1*pz1);
	eta1=0.5*TMath::log((p+pz1)/(p-pz1));
	phi1=TMath::atan(py1/px1);
	imp1=TMath::Sqrt(xvert1*xvert1+yvert1*yvert1+zvert1*zvert1);
      }
      if(j==0){ //get info of triggered event from file2
	sprintf(dummy,"good-pptracks-%d",event1);
	TNtuple *ntuppel2 = (TNtuple*)file2.Get(dummy);
	if(ntuppel2){
	  Float_t nhits2,pdg2,sector2;
	  ntuppel2->SetBranchAddress("mc",&mc2);
	  ntuppel2->SetBranchAddress("pdg",&pdg2);
	  ntuppel2->SetBranchAddress("px",&px2);
	  ntuppel2->SetBranchAddress("py",&py2);
	  ntuppel2->SetBranchAddress("pz",&pz2);
	  ntuppel2->SetBranchAddress("xvert",&xvert2);
	  ntuppel2->SetBranchAddress("yvert",&yvert2);
	  ntuppel2->SetBranchAddress("zvert",&zvert2);
	  ntuppel2->SetBranchAddress("nhits",&nhits2);
	  ntuppel2->SetBranchAddress("sector",&sector2);
	  event2=event1;
	  nent2=ntuppel2->GetEntries();
	  for(Int_t k=0;k<nent2;k++){
	    ntuppel2->GetEvent(k);
	    if(mc2<0)continue; //dont count fake tracks
	    //if(TMath::fabs(zvert2)>30) continue;
	    pt2=TMath::Sqrt(px2*px2+py2*py2);
	    if(nhits2<63) continue;
	    if(pt2<0.1) continue;
	    if((px2==0)&&(py2==0)&&(pz2==0)) continue;
	    Float_t p=TMath::Sqrt(px2*px2+py2*py2+pz2*pz2);
	    eta2=0.5*TMath::log((p+pz2)/(p-pz2));
	    if(TMath::fabs(eta2)>0.9) continue;
	    phi2=TMath::atan(py2/px2);
	    imp2=TMath::Sqrt(xvert2*xvert2+yvert2*yvert2+zvert2*zvert2);
	    tnorm++;
	  }
	} else {
	  sprintf(dummy,"pptracks-%d",event1);
	  TNtuple *ntuppel2 = (TNtuple*)file2.Get(dummy);
	  if(ntuppel2){
	    ntuppel2->SetBranchAddress("pt",&pt2);
	    ntuppel2->SetBranchAddress("phi",&phi2);
	    ntuppel2->SetBranchAddress("eta",&eta2);
	    ntuppel2->SetBranchAddress("xvert",&xvert2);
	    ntuppel2->SetBranchAddress("yvert",&yvert2);
	    ntuppel2->SetBranchAddress("zvert",&zvert2);
	    ntuppel2->SetBranchAddress("imp",&imp2);
	    ntuppel2->SetBranchAddress("px",&px2);
	    ntuppel2->SetBranchAddress("py",&py2);
	    ntuppel2->SetBranchAddress("pz",&pz2);
	    ntuppel2->SetBranchAddress("event",&event2);
	    ntuppel2->SetBranchAddress("mc",&mc2);

	    nent2=ntuppel2->GetEntries();
	    for(Int_t k=0;k<nent2;k++){
	      ntuppel2->GetEvent(k);
	      if(mc2<0)continue; //dont count fake tracks
	      if(pt2<0.1) continue;
	      //if(TMath::fabs(zvert2)>30) continue;
	      if(TMath::fabs(eta2)>0.9) continue;
	      tnorm++;
	      //cout << k << ": " << pt2 << " " << mc2 <<endl;
	    }
	  }
	}
      }
      if(tnorm<MINTRACKS_){
	tnorm=0;
	break; //triggered event has to have one track at least
      }

      if(pt1<0.1) continue;
      if(TMath::fabs(eta1>0.9)) continue;

      if(mc1<-1) continue; //dont count fake tracks
      else if(mc1==-1) fnorm++;

      for(Int_t k=0;k<N_;k++){ //do the counting
	Float_t cut=xcut[k];
	//Float_t impt=TMath::Sqrt(xvert1*xvert1+yvert1*yvert1);
	Float_t impz=TMath::abs(zvert1);
	if(impz>cut) continue;
	if(mc1==-1) fcounts[k]++;
	else tcounts[k]++;
      }

    }//pp track loop

    if(tnorm==0) continue; //break from above
    countgoodevents++; //for normalization

    totnorm+=tnorm;

    for(Int_t k=0;k<N_;k++){ //do the counting
      if(tcounts[k]> tnorm)tcounts[k]=tnorm; 
      thcounts[k]=(Float_t)tcounts[k]/tnorm;
      fhcounts[k]=(Float_t)fcounts[k]/tnorm;
      tallcounts[k]+=tcounts[k];
      fallcounts[k]+=fcounts[k];
      thmeancounts[k]+=thcounts[k];
      fhmeancounts[k]+=fhcounts[k];
      thsigmacounts[k]+=thcounts[k]*thcounts[k];
      fhsigmacounts[k]+=fhcounts[k]*fhcounts[k];
      if((k==N_-1)&&(thcounts[k]< 0.99)){ 
	//cout << "Warning " << i << " " << event2 << " " << tcounts[k] << " " << tnorm << endl;
      }
    }
  }//pileup loop

  file1.Close();
  file2.Close(); 

  cout << "Events used: " << countgoodevents << endl;
  for(Int_t k=0;k<N_;k++){ //do the counting
    thmeancounts[k]/=countgoodevents;
    fhmeancounts[k]/=countgoodevents;

    tallcounts[k]/=totnorm;
    fallcounts[k]/=totnorm;

    if(countgoodevents>1){
      Float_t stemp=thsigmacounts[k]/countgoodevents;
      Float_t s2temp=fhsigmacounts[k]/countgoodevents;
      Int_t N=countgoodevents-1;
      thsigmacounts[k]=TMath::sqrt((stemp -thmeancounts[k]*thmeancounts[k])/N);
      fhsigmacounts[k]=TMath::sqrt((s2temp-fhmeancounts[k]*fhmeancounts[k])/N);
    }else{
      thsigmacounts[k]=0;
      fhsigmacounts[k]=0;
    }
    cout << k << ": " << thmeancounts[k] << "+-"<< thsigmacounts[k] << " " << fhmeancounts[k] << "+-" << fhsigmacounts[k] << endl;
  }

  TCanvas *c1 = new TCanvas("c1","PileUp",1000,500);
  //c1->SetFillColor(42);
  //c1->SetGrid();
  //c1->GetFrame()->SetFillColor(21);
  //c1->GetFrame()->SetBorderSize(12);

  TGraphErrors *g1=new TGraphErrors(N_,xcut,thmeancounts,0,&thsigmacounts[0]);
  //TGraphErrors *g1=new TGraphErrors(N_,xcut,tallcounts);
  sprintf(dummy,"Good tracks %s",LTEXT_);
  g1->SetTitle(dummy);
  g1->SetMarkerColor(4);
  g1->GetHistogram()->SetXTitle("Z_{cut} [cm]");
  //g1->GetHistogram()->SetYTitle("Efficiency");
  g1->SetMarkerStyle(21);
  
  TGraphErrors *g2=new TGraphErrors(N_,xcut,fhmeancounts,0,fhsigmacounts);
  //TGraphErrors *g2=new TGraphErrors(N_,xcut,fallcounts);
  g2->SetTitle("Fake tracks (relative to good triggered tracks)");
  g2->GetHistogram()->SetXTitle("Z_{cut} [cm]");
  g2->SetMarkerColor(4);
  g2->SetMarkerStyle(21);

  c1->Divide(2,1);
  c1->cd(1);
  g1->Draw("AP");
  c1->cd(2);
  g2->Draw("AP");
  c1->Update();
}

/* Plot the vertex reconstruction resolution */

void plot_pp_vert(Int_t evi=0,Int_t evs=100,Char_t *infile="pp-tracks.root")
{
  const Int_t zcut=3;
  const Int_t zcutplot=5;

  Char_t dummy[1000];
  Char_t fname[1000];
  Char_t fname2[1000];

  Float_t pt,phi,eta,xvert,yvert,zvert,imp,nhits,px,py,pz,event,mc;
  TH1F *vhist = new TH1F("zhist","Vertex reconstruction resolution",50,-zcutplot,zcutplot);

  TFile file = TFile(infile,"READ");
  if(!file.IsOpen()) return;

  TNtuple *ntuppel;

  for(int ev=evi; ev<evs; ev++){
    //sprintf(dummy,"good-pptracks-%d",ev); /*dont use this*/
    sprintf(dummy,"pptracks-%d",ev);

    if(ntuppel) delete ntuppel;
    ntuppel = (TNtuple*)file.Get(dummy);
    if(!ntuppel) continue;
    if(ntuppel.GetEntries()==0) continue;

    sprintf(fname,"(eta >-0.9) && (eta<0.9) && (pt>0.1) && (nhits>63) && (zvert<%d) && (zvert>-%d)",zcut,zcut);
    ntuppel->Draw("zvert>>histo",fname,"groff");
    Float_t mean = histo->GetMean();
    //cout << mean << endl;
    vhist->Fill(mean);
  }
  
  gStyle->SetStatColor(10);
  gStyle->SetOptFit(1);

  TF1 *f1 = new TF1("f1","gaus",-(Float_t)zcutplot/2.,(Float_t)zcutplot/2.);
  vhist->Fit("f1","R");

  vhist->SetXTitle("Z* [cm]");
  vhist->Draw("E1P");
}

/* Plot the impact parameter resolution */
 
void plot_pp(Int_t evi=0,Int_t evs=100,Char_t *infile="pp-tracks.root")
{
  const Int_t zcut=3;
  const Int_t zcutplot=5;

  Char_t dummy[1000];
  Char_t fname[1000];
  Char_t fname2[1000];

  TFile file = TFile(infile,"READ");
  if(!file.IsOpen()) return;
  
  TH1F *vhist = new TH1F("zhist","Impact parameter resolution",50,-zcut,zcut);
  
  TNtuple *ntuppel;
  Int_t a=0;
  for(int ev=evi; ev<evs; ev++){
    //sprintf(dummy,"good-pptracks-%d",ev); /*dont use this*/
    sprintf(dummy,"pptracks-%d",ev);

    if(ntuppel) delete ntuppel;
    ntuppel = (TNtuple*)file.Get(dummy);
    if(!ntuppel) continue;
    if(ntuppel.GetEntries()==0) continue;

    sprintf(fname,"(eta>-0.9) && (eta<0.9) && (pt>0.1) && (nhits>63) && (zvert<%d) && (zvert>-%d)",zcut,zcut);
    ntuppel->Draw("zvert>>histo",fname,"groff");
    Float_t mean = histo->GetMean();
    a++;
    //cout << mean << endl;

    sprintf(fname2,"(zvert-%f)>>+zhist",mean);
    ntuppel->Draw(fname2,fname);
  }
  
  gStyle->SetStatColor(10);
  gStyle->SetOptFit(1);
  cout << a <<endl;
  TF1 *f1 = new TF1("f1","gaus",-(Float_t)zcutplot/2.,(Float_t)zcutplot/2.);
  vhist->Fit("f1","R");

  vhist->SetXTitle("Z_{DCA}-Z* [cm]");
  vhist->DrawCopy("E1P");
  file->Close();
}


#ifdef OLD_TRIGGER_MACRO
void trigger_pp(char *outfile="results.root")
{
  
  TNtuple *ntuppel = new TNtuple("ntuppel","","pt:eta:xvert:yvert:zvert:nhits:px:py:pz:event");
  Float_t meas[10];  
  
  for(int event=0; event<1; event++)
    {
      char fname[256];
      sprintf(fname,"/data1/AliRoot/pp/pileup/tracks_0.raw");
      //sprintf(fname,"aliruntfile.root");

      //Get the tracks:
      /*AliHLTTrackArray *tracks = new AliHLTTrackArray();
      AliHLTFileHandler *file = new AliHLTFileHandler();
      file->SetBinaryInput(fname);
      file->Binary2TrackArray(tracks);
      file->CloseBinaryInput();
      delete file;*/

      sprintf(fname,"/data1/AliRoot/pp/pileup/");
	
      AliHLTEvaluate *eval=new AliHLTEvaluate(fname,63,63);
      eval->LoadData(event,-1);
      eval->AssignIDs();
      
      AliHLTTrackArray *tracks=eval->GetTracks();

      sprintf(fname,"/data1/AliRoot/pp/pileup/");
      //sprintf(fname,"/prog/alice/data/Rawdata/1_patch/pp/recon_%d/",event);

      Int_t ntracks=0;
      Double_t xc,yc,zc;
      Double_t impact;
      AliHLTVertex vertex;

      Int_t mcid = 0;

      AliHLTTrackArray *ftracks = new AliHLTTrackArray();
      
      for(int i=0; i<tracks->GetNTracks(); i++)
	{
	  track = (AliHLTTrack*)tracks->GetCheckedTrack(i);
	  if(!track) continue;
	  
	  //Assign MCid
	  mcid = track->GetMCid();
	  //if (mcid != -1) 
	  //cout << "MCid " << mcid << endl; 
	  
	  track->CalculateHelix();
	  //cout<<"Pt "<<track->GetPt()<<" eta "<<track->GetPseudoRapidity()<<" Nhits "<<track->GetNHits()<<endl;
	  //cout<<"Before second fit; pt "<<track->GetPt()<<" firstpoint "<<track->GetFirstPointX()<<" "<<track->GetFirstPointY()<<" "<<track->GetFirstPointZ()<<endl;

	  track->CalculateHelix();
	  track->GetClosestPoint(&vertex,xc,yc,zc);
	  meas[0]=track->GetPt();
	  meas[1]=track->GetPseudoRapidity();
	  meas[2]=xc;
	  meas[3]=yc;
	  meas[4]=zc;
	  meas[5]=track->GetNHits();
	  meas[6]=track->GetPx();
	  meas[7]=track->GetPy();
	  meas[8]=track->GetPz();
	  meas[9]=event;
	  ntuppel->Fill(meas);
	  //continue;
	  if(fabs(track->GetPseudoRapidity())>0.9) continue;
	  if(track->GetNHits() < 100) continue;
	  if(track->GetPt()<0.2) continue;

	  impact = sqrt(xc*xc+yc*yc+zc*zc);
	  if(fabs(zc)>3) continue;
	  ftracks->AddLast(track);
	  cout<<"-------------------------------------"<<endl;
	  cout<<"Number of hits "<<track->GetNHits()<<endl;
	  cout<<"Transversal impact "<<sqrt(xc*xc+yc*yc)<<endl;
	  cout<<"Longitudinal impact "<<zc<<endl;
	  cout<<"Total impact "<<sqrt(xc*xc+yc*yc+zc*zc)<<endl;
	  cout<<"xc "<<xc<<" yc "<<yc<<" zc "<<zc<<endl;
	  
	  //cout<<"After second fit; pt "<<track->GetPt()<<" firstpoint "<<track->GetFirstPointX()<<" "<<track->GetFirstPointY()<<" "<<track->GetFirstPointZ()<<endl;
	  cout<<"pt "<<track->GetPt()<<" eta "<<track->GetPseudoRapidity()<<endl;
	  
	  ntracks++;
	}
      cout<<endl<<"There was "<<ntracks<<" accepted tracks, out of total "<<tracks->GetNTracks()<<" found tracks"<<endl;
      
      //display(ftracks,fname);
      delete tracks;
      delete ftracks;
    }
  
  TFile *rfile = TFile::Open(outfile,"RECREATE");
  ntuppel->Write();
  rfile->Close();
  delete ntuppel;
}

void display(AliHLTTrackArray *tracks,char *path)
{
  int slice[2]={0,35};
  d = new AliHLTDisplay(slice);
  d->Setup("tracks_0.raw",path);
  d->SetTracks(tracks);
  //d->DisplayClusters();
  d->DisplayAll();
  //d->DisplayTracks();
  
}

void ploteff()
{
  //double x[6]={1,2,4,6,8,10};
  //double y[6]={0.873684,1.01379,1.17751,1.28614,1.31638,1.32022};
  
  double x[6]={0.8,1.6,3.2,4.8,6.4,8};
  double y[6]={0.497006,0.88024,1.19162,1.23952,1.39521,1.40719};//->sigmas

  
  TGraph *gr = new TGraph(6,x,y);
  c1 = new TCanvas("c1","",2);
  SetCanvasOptions(c1);
  gr->SetMarkerStyle(20);
  gr->SetTitle("");
  gr->Draw("APL");
  gr->GetHistogram()->SetXTitle("Z_{CUT} [x #sigma_{Z}]");
  gr->GetHistogram()->SetYTitle("Efficiency");
  SetTGraphOptions(gr);
}

double geteff(char *infile,char *singlefile,double cut)
{
  gStyle->SetStatColor(10);
  gStyle->SetOptFit(1);
  Int_t pileups[25],single[25];
  file = TFile::Open(infile,"READ");
  
  char name[256];
  for(int i=0; i<25; i++)
    {
      sprintf(name,"zvert < %f && zvert > %f && pt>0.2 && nhits>100 && eta>-0.9 && eta<0.9 && event==%d",cut,-1.*cut,i);
      ntuppel->Draw(">>eventlist",name);
      int entries = eventlist->GetN();
      pileups[i]=entries;
    }
  //eventlist->Print("all");
  file->Close();
  
  file = TFile::Open(singlefile,"read");
  for(int i=0; i<25; i++)
    {
      sprintf(name,"zvert < %f && zvert > %f && pt>0.2 && nhits>100 && eta>-0.9 && eta<0.9 && event==%d",3,-3,i);
      ntuppel->Draw(">>eventlist",name);
      int entries = eventlist->GetN();
      single[i]=entries;
    }
  file->Close();
  double totsingle=0,totpileup=0;
  for(int i=0; i<25; i++)
    {
      totsingle += single[i];
      totpileup += pileups[i];
    }
  double toteff = totpileup/totsingle;
  cout<<"Total eff "<<toteff<<endl;
  
  return toteff;
}


void plot(char *infile)
{
  gStyle->SetStatColor(10);
  gStyle->SetOptFit(1);
  file = TFile::Open(infile,"READ");
  
  histo = new TH1F("histo","histo",20,-3,3);
  
  vhist = new TH1F("vhist","",20,-3,3);
  SetTH1Options(vhist);
  
  char fname[256];
  char fname2[256];
  for(int ev=0; ev<25; ev++)
    {
      sprintf(fname,"pt>0.2 && nhits>100 && eta>-0.9 && eta<0.9 && event==%d",ev);
      ntuppel->Draw("zvert>>histo",fname,"goff");
      float mean = histo->GetMean();
      vhist->Fill(mean);
      continue;
      sprintf(fname2,"(zvert-(%f))>>+vhist",mean);
      cout<<fname2<<endl;
      ntuppel->Draw(fname2,fname,"goff");
    }
  
  c1 = new TCanvas("c1","",2);
  SetCanvasOptions(c1);
  vhist->SetXTitle("Z* [cm]");
  vhist->Draw();
  return;
  //ntuppel->Draw("zvert>>histo","pt>0.2");
  TF1 *f1 = new TF1("f1","gaus",-3,3);
  histo->Fit("f1","R");

  //histo->Draw();
  //file->Close();
}

enum tagprimary {kPrimaryCharged=0x4000};
void LoadEvent(Int_t event=0)
{
  //Load the generated particles
  
  gROOT->LoadMacro("$(ALICE_ROOT)/macros/loadlibs.C");
  loadlibs();
  //TFile *rootfile = TFile::Open("/prog/alice/data/pro/25event_pp.root","READ");
  TFile *rootfile = TFile::Open("/prog/alice/pp/pileup/100pileup/alirunfile.root","READ");
  gAlice = (AliRun*)rootfile->Get("gAlice");
  
  //  TNtuple *ntup = new TNtuple(
  
  Int_t nparticles = gAlice->GetEvent(event);
  Int_t nprim = FindPrimaries(nparticles,0.,0.,0.);
  cout<<"Number of primaries "<<nprim<<endl;
  int co=0;
  for(Int_t i=0; i<nparticles; i++)
    {
      TParticle *part = gAlice->Particle(i);
      if(!part->TestBit(kPrimaryCharged)) continue;
      if(fabs(part->Eta())>0.9) continue;
      cout<<part->GetName()<<" pt "<<part->Pt()<<" eta "<<part->Eta()<<" xvert "<<part->Vx()<<" yvert "<<part->Vy()<<" zvert "<<part->Vz()<<endl;
      co++;
    }
  cout<<endl<<"Number of primary tracks in the detector: "<<co<<endl;
  gAlice=0;
  rootfile->Close();
  
}

Int_t FindPrimaries(Int_t nparticles, Double_t xori, Double_t yori, Double_t zori)
{
  //Define primary particles in a pp-event. Code taken from offline.


  // cuts:
  //Double_t vertcut = 0.001;  // cut on the vertex position
  Double_t vertcut = 2.0;
  Double_t decacut = 3.;     // cut if the part. decays close to the vert.
  Double_t timecut = 0.;
  
  TList *listprim = new TList();
  listprim->SetOwner(kFALSE);
  
  Int_t nprch1=0;
  Int_t nprch2=0;
  for(Int_t iprim = 0; iprim<nparticles; iprim++){   //loop on  tracks
    
    TParticle * part = gAlice->Particle(iprim);
    char *xxx=strstr(part->GetName(),"XXX");
    if(xxx)continue;
    
    TParticlePDG *ppdg = part->GetPDG();
    //if(TMath::Abs(ppdg->Charge())!=3)continue;  // only charged (no quarks)
    if(TMath::Abs(ppdg->Charge())<1)continue;  // only charged (no quarks)
    
    Double_t dist=TMath::Sqrt((part->Vx()-xori)*(part->Vx()-xori)+(part->Vy()-yori)*(part->Vy()-yori)+(part->Vz()-zori)*(part->Vz()-zori));
    if(dist>vertcut)continue;  // cut on the vertex

    if(part->T()>timecut)continue;

    Double_t ptot=TMath::Sqrt(part->Px()*part->Px()+part->Py()*part->Py()+part->Pz()*part->Pz());
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
        Double_t rad = TMath::Sqrt((dau->Vx()-xori)*(dau->Vx()-xori)+(dau->Vy()-yori)*(dau->Vy()-yori));
	if(distd<decacut)prmch=kFALSE;  // eliminate if the decay is near the vertex
	//if(rad < 20)prmch=kFALSE;
      }
    }

    if(prmch){
      nprch1++;
      part->SetBit(kPrimaryCharged);
      listprim->Add(part);    // list of primary particles (before cleanup)
    }
  }


  nprch2=0;
  for(Int_t iprim = 0; iprim<nparticles; iprim++){ // cleanup loop
    TParticle * part = gAlice->Particle(iprim);
    if(part->TestBit(kPrimaryCharged)){
      Int_t mothind=part->GetFirstMother();
      if(mothind<0)continue;
      TParticle *moth=gAlice->Particle(mothind);
      TParticle *mothb=(TParticle*)listprim->FindObject(moth);
      if(mothb){
        listprim->Remove(moth);
        moth->ResetBit(kPrimaryCharged);
        nprch2++;
      }
    }
  }

  listprim->Clear("nodelete");
  delete listprim;
  return nprch1-nprch2;
  
}
#endif
