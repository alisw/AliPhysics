void trigger()
{
  
  //Get the tracks:
  AliL3TrackArray *tracks = new AliL3TrackArray();
  AliL3FileHandler *file = new AliL3FileHandler();
  file->SetBinaryInput("tracks.raw");
  file->Binary2TrackArray(tracks);
  file->CloseBinaryInput();
  delete file;
  
  Int_t ntracks=0;
  Double_t xc,yc,zc;
  Double_t impact;
  AliL3Vertex vertex;
  AliL3Fitter *fitter = new AliL3Fitter(&vertex);
  fitter->LoadClusters("./");
  //fitter->NoVertex();
  for(int i=0; i<tracks->GetNTracks(); i++)
    {
      track = (AliL3Track*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      if(fabs(track->GetPseudoRapidity())>0.9) continue;
      if(track->GetNHits() < 50) continue;
      if(track->GetPt()<0.1) continue;
      track->CalculateHelix();
      //cout<<"Pt "<<track->GetPt()<<" eta "<<track->GetPseudoRapidity()<<" Nhits "<<track->GetNHits()<<endl;
      //cout<<"Before second fit; pt "<<track->GetPt()<<" firstpoint "<<track->GetFirstPointX()<<" "<<track->GetFirstPointY()<<" "<<track->GetFirstPointZ()<<endl;
      fitter->FitHelix(track);
      track->CalculateHelix();
      
      track->GetClosestPoint(&vertex,xc,yc,zc);
      impact = sqrt(xc*xc+yc*yc+zc*zc);
      if(impact>3) continue;
      cout<<"Number of hits "<<track->GetNHits()<<endl;
      cout<<"Transversal impact "<<sqrt(xc*xc+yc*yc)<<endl;
      cout<<"Longitudinal impact "<<zc<<endl;
      cout<<"Total impact "<<sqrt(xc*xc+yc*yc+zc*zc)<<endl;
      cout<<"xc "<<xc<<" yc "<<yc<<" zc "<<zc<<endl;
      
      //cout<<"After second fit; pt "<<track->GetPt()<<" firstpoint "<<track->GetFirstPointX()<<" "<<track->GetFirstPointY()<<" "<<track->GetFirstPointZ()<<endl;
      cout<<"pt "<<track->GetPt()<<" eta "<<track->GetPseudoRapidity()<<endl;
      
      cout<<"-------------------------------------"<<endl;
      ntracks++;
    }
  cout<<endl<<"There was "<<ntracks<<" found tracks"<<endl;
  delete tracks;
  delete fitter;
}

enum tagprimary {kPrimaryCharged=0x4000};
void LoadEvent(Int_t event=0)
{
  //Load the generated particles

  gROOT->LoadMacro("$(ALICE_ROOT)/macros/loadlibs.C");
  loadlibs();
  TFile *rootfile = TFile::Open("/prog/alice/data/pro/25event_pp.root","READ");
  gAlice = (AliRun*)rootfile->Get("gAlice");
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
        if(distd<decacut)prmch=kFALSE;  // eliminate if the decay is near the vertex
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
