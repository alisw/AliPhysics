// $Id$

void trigger_pp(char *outfile="results.root")
{
  
  TNtuple *ntuppel = new TNtuple("ntuppel","","pt:eta:xvert:yvert:zvert:nhits:px:py:pz:event");
  Float_t meas[10];  
  
  for(int event=0; event<1; event++)
    {
      char fname[256];
      sprintf(fname,"/prog/alice/data/Rawdata/1_patch/pp/pileups/recon_%d/tracks.raw",event);
      //sprintf(fname,"/prog/alice/data/Rawdata/1_patch/pp/recon_%d/tracks.raw",event);

      //Get the tracks:
      AliL3TrackArray *tracks = new AliL3TrackArray();
      AliL3FileHandler *file = new AliL3FileHandler();
      file->SetBinaryInput(fname);
      file->Binary2TrackArray(tracks);
      file->CloseBinaryInput();
      delete file;
      
      sprintf(fname,"/prog/alice/data/Rawdata/1_patch/pp/pileups/recon_%d/",event);
      //sprintf(fname,"/prog/alice/data/Rawdata/1_patch/pp/recon_%d/",event);

      Int_t ntracks=0;
      Double_t xc,yc,zc;
      Double_t impact;
      AliL3Vertex vertex;
      AliL3Fitter *fitter = new AliL3Fitter(&vertex);
      fitter->LoadClusters(fname);
      //fitter->NoVertex();
      
      AliL3TrackArray *ftracks = new AliL3TrackArray();
      
      for(int i=0; i<tracks->GetNTracks(); i++)
	{
	  track = (AliL3Track*)tracks->GetCheckedTrack(i);
	  if(!track) continue;
	  track->CalculateHelix();
	  //cout<<"Pt "<<track->GetPt()<<" eta "<<track->GetPseudoRapidity()<<" Nhits "<<track->GetNHits()<<endl;
	  //cout<<"Before second fit; pt "<<track->GetPt()<<" firstpoint "<<track->GetFirstPointX()<<" "<<track->GetFirstPointY()<<" "<<track->GetFirstPointZ()<<endl;
	  fitter->FitHelix(track);//refit the tracks
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
      
      display(ftracks,fname);
      delete tracks;
      delete fitter;
      delete ftracks;
    }
  
  TFile *rfile = TFile::Open(outfile,"RECREATE");
  ntuppel->Write();
  rfile->Close();
  delete ntuppel;

}

void display(AliL3TrackArray *tracks,char *path)
{
  int slice[2]={0,35};
  d = new AliL3Display(slice);
  d->Setup("tracks.raw",path);
  d->SetTracks(tracks);
  //d->DisplayClusters();
    d->DisplayAll();
    //d->DisplayTracks();
  
}

void ploteff()
{
  gROOT->LoadMacro("XFunct.C");
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
  gROOT->LoadMacro("XFunct.C");
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
  gROOT->LoadMacro("XFunct.C");
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
  TFile *rootfile = TFile::Open("/prog/alice/data/pro/25event_pp.root","READ");
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
  Double_t vertcut = 0.001;  // cut on the vertex position
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
    if(TMath::Abs(ppdg->Charge())!=3)continue;  // only charged (no quarks)
    
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
