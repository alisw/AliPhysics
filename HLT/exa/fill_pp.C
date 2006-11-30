// $Id$


/**
   fill_pp -> store track info in ntuple root file
   fill_good_tracks_pp -> store good track info in ntuple root file
   print_pp -> print selected events from ntuple root file
   plot_pp -> plot parts of selected events from ntuple root file
   eval_pp -> simple evaluation of impact parameter cut effects
   Author: C. Loizides
*/


/* ---fill_pp---
   Fill ntuple with relevant pp event data 
   from the tracks files:

   Ev is the event number to use
   Path gives the path where you find the digit files (and possibly also l3transform.config)
   Rfile gives the root file where the information will be added
   Saveev will be Ev unless you specify something else (useful for storing the trigger event of pilup events)
*/

void fill_pp(Int_t ev=0,Char_t *path="./",Char_t *rfile="pptracks.root",Int_t saveev=-1)
{
  AliHLTLogger l;
  l.Set(AliHLTLogger::kAll);
  l.UseStdout();
  //l.UseStream();

  if(getenv("TRANSFORMCONFIGPATH")){
    AliHLTTransform::Init(getenv("TRANSFORMCONFIGPATH"));
  } else AliHLTTransform::Init(path);

  TNtuple *ntuppel = new TNtuple("ntuppel","pptracks","pt:phi:eta:xvert:yvert:zvert:imp:nhits:px:py:pz:event:mc");
  Float_t meas[13];  

  TFile file(rfile,"UPDATE");
  Char_t name[1024];

  AliHLTEvaluate *eval=new AliHLTEvaluate(path,63,63);
  eval->LoadData(ev,-1);
  eval->AssignIDs();

  AliHLTTrackArray *a=eval->GetTracks();
  a->Compress();

  AliHLTTrack t;
  Float_t phi,pt,eta;
  Int_t id;
  Int_t nhits;

  AliHLTVertex vertex;
  Double_t xc,yc,zc,impact;

  Int_t ntracks=a.GetNTracks();
  for(Int_t i=0;i<ntracks;i++){
    t.Set(a.GetTrack(i));

    t.CalculateHelix();
    t.GetClosestPoint(&vertex,xc,yc,zc);

    nhits=t.GetNHits();
    if(nhits<63) continue;
    pt=t.GetPt();
    eta=t.GetPseudoRapidity();
    phi=t.GetPsi();
    id=t.GetMCid();

    //impact=TMath::fabs(zc);
    impact=TMath::Sqrt(xc*xc+yc*yc+zc*zc);

    meas[0]=pt;
    meas[1]=phi;
    meas[2]=eta;
    meas[3]=xc;
    meas[4]=yc;   
    meas[5]=zc;
    meas[6]=impact;   
    meas[7]=nhits;
    meas[8]=t.GetPx();
    meas[9]=t.GetPy();
    meas[10]=t.GetPz();
    meas[11]=ev;
    meas[12]=id;

    ntuppel->Fill(meas);

    cout << i << ": pt " << pt << " " << phi << " " << eta << " " << id << " nhits " << nhits << " impact " << impact << endl;

    //if(id!=-1) cout << "mc " << id << " " << t.GetCenterX() << " " << t.GetCenterY() << " " << t.GetRadius() << " fst point " << t.GetFirstPointX() << " " << t.GetFirstPointY() << " " << t.GetFirstPointZ() << " tgl " << t.GetTgl() << endl;
  } 

  if(saveev==-1) saveev=ev;
  sprintf(name,"pptracks-%d",saveev);
  ntuppel->Write(name);
  file.Close();

  delete ntuppel;
  delete eval;
}

/*--fill_good_tracks_pp--

  Fill good tracks taken from offline into ntuple root file for later evaluation.

  Evs is the event to start with
  Eve is the last event to process
  Path gives the path where the offline good_particles_tpc files are
  Rfile specifies the recreated root file where the ntuples will be stored
*/

void fill_good_tracks_pp(Int_t evs=0,Int_t eve=100,Char_t *path="./",Char_t *rfile="pp-good-tracks.root")
{
  TFile file(rfile,"RECREATE");
  if(!file.IsOpen()) return;

  Float_t meas[10];  
  Char_t name[1024];

  for(Int_t i=evs;i<eve;i++){
    TNtuple *ntuppel = new TNtuple("ntuppel-goodtracks","good-pptracks","mc:pdg:px:py:pz:xvert:yvert:zvert:nhits:sector");

    sprintf(name,"%s/good_tracks_tpc_%d",path,i);
    FILE *f=fopen(name,"r");
    cout << "Event " << i << " " << name << endl;

    while(!feof(f)){
      Int_t mc,pdg,nhits,sector;
      Float_t px,py,pz,xv,yv,zv;
      fscanf(f,"%d %d %f %f %f %f %f %f %d %d\n",&mc,&pdg,&px,&py,&pz,&xv,&yv,&zv,&nhits,&sector);
      
      if(nhits<63) continue;
      if((px==0)&&(py==0)&&(pz==0)) continue;
      if((xv==0)&&(yv==0)&&(zv==0)) continue;

      //cout << mc << " " << pdg << " " << px << " " << py << " " << pz << " " << xv << " " << yv << " " << zv << " " << nhits << " " << sector << endl;

      meas[0]=mc;      
      meas[1]=pdg;      
      meas[2]=px;      
      meas[3]=py;      
      meas[4]=pz;      
      meas[5]=xv;      
      meas[6]=yv;      
      meas[7]=zv;      
      meas[8]=nhits;      
      meas[9]=sector;      
      
      ntuppel->Fill(meas);
    }
    fclose(f);
    sprintf(name,"good-pptracks-%d",i);
    ntuppel->Write(name);
    delete ntuppel;
  }

  file.Close();
}

/* ---print_pp---
   Looping over event reading root file containing ntuples 
   and printing wanted information:

   Pfile specifies the root file created with fill_pp
   Evi gives the first event to print
   Eve the last event to print
*/

void print_pp(Char_t *pfile, Int_t evi=0, Int_t eve=-1)
{
  TFile file1=TFile(pfile,"READ");

  if(!file1.IsOpen()) return;

  //"pt:phi:eta:xvert:yvert:zvert:imp:nhits:px:py:pz:event:mc")
  Float_t pt1,phi1,eta1,xvert1,yvert1,zvert1,imp1,nhits1,px1,py1,pz1,event1,mc1; //pileup
  Float_t pt2,phi2,eta2,xvert2,yvert2,zvert2,imp2,nhits2,px2,py2,pz2,event2,mc2; //pp

  Int_t nent1;
  Char_t dummy[1024];
  if(eve==-1) eve=evi;

  for(Int_t i=evi;i<=eve;i++){ //loop over pileup-events
    sprintf(dummy,"pptracks-%d",i);
    TNtuple *ntuppel1 = (TNtuple*)file1.Get(dummy);
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

    nent1=ntuppel1->GetEntries();
    cout << "Event: " << i << " " << nent1 << " tracks" << endl;
    for(Int_t j=0;j<nent1;j++){
      ntuppel1->GetEvent(j);

      Float_t bt=TMath::Sqrt(xvert1*xvert1+yvert1*yvert1);
      Float_t bz=TMath::fabs(zvert1);

      if(mc1!=-1)
      cout << j << ": " << event1 << " " << pt1 << " " << phi1 << " " << eta1 << " mc " << mc1 << " bt " << bt << " bz " << bz << endl;

    }//pp track loop
  }//pileup loop

  file1.Close();
}

/* ---plot_pp---
   Plot some variables from the ntuples (currently impact parameter
   but easily changeable.

   Pfile is the root file containing the ntuples.
   Evi and Eve are start and end event number to process
*/

void plot_pp(Char_t *pfile, Int_t evi=0, Int_t eve=100)
{
  TH1F *hist1=new TH1F("h1","test",100,0,1000);

  TFile file1=TFile(pfile,"READ");

  if(!file1.IsOpen()) return;

  //"pt:phi:eta:xvert:yvert:zvert:imp:nhits:px:py:pz:event:mc")
  Float_t pt1,phi1,eta1,xvert1,yvert1,zvert1,imp1,nhits1,px1,py1,pz1,event1,mc1; //pileup
  Float_t pt2,phi2,eta2,xvert2,yvert2,zvert2,imp2,nhits2,px2,py2,pz2,event2,mc2; //pp

  Int_t nent1;
  Char_t dummy[1024];

  for(Int_t i=evi;i<eve;i++){ //loop over pileup-events
    sprintf(dummy,"pptracks-%d",i);
    TNtuple *ntuppel1 = (TNtuple*)file1.Get(dummy);
    if(!ntuppel1){
      cerr << "Error getting ntuple from " << dummy<<endl;
      continue;
    } 
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

    nent1=ntuppel1->GetEntries();
    cout << "Event: " << i << " " << nent1 << " tracks" << endl;
    for(Int_t j=0;j<nent1;j++){
      ntuppel1->GetEvent(j);

      if(TMath::fabs(eta1)>0.9) continue;
      if(pt1<0.1) continue;

      Float_t bt=TMath::Sqrt(xvert1*xvert1+yvert1*yvert1);
      Float_t bz=TMath::fabs(zvert1);
      hist1->Fill(bz,1);
    }//pp track loop
  }//pileup loop

  file1.Close();

  hist1->Sumw2();
  hist1->Draw("e1p");
}

/* ---eval_pp---
   Little macro to play with pp pileup evaluation and do some
   counting. See real trigger macro trigger.C instead.

   Ev is the event number
   Path is the path to the track arrays.
*/

void eval_pp(Int_t ev=0,Char_t *path="./")
{
  AliHLTLogger l;
  l.Set(AliHLTLogger::kAll);
  l.UseStdout();
  //l.UseStream();

  AliHLTEvaluate *eval=new AliHLTEvaluate(path,63,63);
  eval->LoadData(ev,-1);
  eval->AssignIDs();

  AliHLTTrackArray *a=eval->GetTracks();

  AliHLTTrack t;
  Float_t phi,pt,tgl,eta;
  Int_t id;
  Int_t nhits;
  Int_t abspos=0,pos=0,absneg=0,neg=0;
  //assume vertex around zero
  AliHLTVertex vertex;
  Double_t xc,yc,zc,impact;

  Int_t nents=0;
  Int_t ntracks=a.GetNTracks();
  for(Int_t i=0;i<ntracks;i++){
    t.Set(a.GetTrack(i));

    t.CalculateHelix();
    t.GetClosestPoint(&vertex,xc,yc,zc);

    nhits=t.GetNHits();
    if(nhits<63) continue;
    pt=t.GetPt();
    tgl=t.GetTgl();
    phi=100;
    id=t.GetMCid();
    cout << id << endl;
    if(t.GetPy()!=0) phi=TMath::ATan(t.GetPy()/t.GetPx());

    //impact=TMath::fabs(zc);
    impact=TMath::Sqrt(xc*xc+yc*yc+zc*zc);

    cout << i << ": pt " << pt << " " << phi << " " << tgl << " " << id << " nhits " << nhits << " impact " << impact << endl;

    if(id<0) absneg++;
    else abspos++;
    nents++;
    if(impact>10) continue;
    if(id<0) neg++;
    else pos++;
  } 

  cout << "Result: " << nents << " " << abspos << " " << absneg << endl;
  cout << "After cut: " << pos << " " << neg << endl;
}

void simple_eff_pp(Char_t *pfile, Char_t *reffile, Int_t evi=0, Int_t eve=100,Int_t cut=-1)
{
  //"pt:phi:eta:xvert:yvert:zvert:imp:nhits:px:py:pz:event:mc")
  Float_t pt1,phi1,eta1,xvert1,yvert1,zvert1,imp1,nhits1,px1,py1,pz1,event1,mc1; //pileup
  Float_t pt2,phi2,eta2,xvert2,yvert2,zvert2,imp2,nhits2,px2,py2,pz2,event2,mc2; //pp

  Int_t countgoodevents=0;
  Int_t nent1,nent2;
  Char_t dummy[1024];

  Int_t tnorm;
  Int_t fnorm;
  Int_t pnorm;
  Int_t tcounts=0;
  Int_t fcounts=0;
  Int_t pcounts=0;
  Float_t fmean=0;
  Float_t pmean=0;

  TFile file1=TFile(pfile,"READ");
  TFile file2=TFile(reffile,"READ");
  if(!file1.IsOpen() || !file2.IsOpen()) return;

  for(Int_t i=evi;i<eve;i++){ //loop over pileup-events
    sprintf(dummy,"pptracks-%d",i);
    TNtuple *ntuppel1 = (TNtuple*)file1.Get(dummy);
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

    tnorm=0;
    fnorm=0;
    pnorm=0;
    nent1=ntuppel1->GetEntries();
    for(Int_t j=0;j<nent1;j++){
      ntuppel1->GetEvent(j);

      if(j==0){ //get info of triggered event from file2
	
	sprintf(dummy,"good-pptracks-%d",event1);
	TNtuple *ntuppel2 = (TNtuple*)file2.Get(dummy);
	Float_t nhits2,pdg2,sector2;
	if(ntuppel2){
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
	    }
	  }
	}
      }
      if(tnorm==0){
	tnorm=0;
	break; //triggered event has to have one track at least
      }
      
      if(pt1<0.1) continue;
      if(TMath::fabs(eta1>0.9)) continue;

      if(mc1<-1) continue; //dont count fake tracks
      else if(mc1==-1) fnorm++; //count tracks from pileup
      else {
	if((cut>0)&&(zvert1>cut)) continue;
	pnorm++; //count real tracks
      }
    }//pp track loop

    if(tnorm==0) continue; //break from above
    if(pnorm>tnorm) pnorm=tnorm; //prob. multiple found tracks
    countgoodevents++; //for normalization
    fcounts+=fnorm;
    pcounts+=pnorm;
    tcounts+=tnorm;
    pmean+=(Float_t)pnorm/(Float_t)tnorm;
    fmean+=(Float_t)fnorm/(Float_t)tnorm;
  }//pileup loop

  file1.Close();
  file2.Close(); 

  Float_t peff,feff;
  peff=(Float_t)pcounts/tcounts;
  feff=(Float_t)fcounts/tcounts;
  pmean/=countgoodevents;
  fmean/=countgoodevents;
  cout << "Events used: " << countgoodevents << endl;
  cout << "Good Eff " << peff << " " << pcounts << " " << tcounts << endl;
  cout << "Fake Eff " << feff << " " << fcounts << " " << tcounts << endl;
  cout << "Means are " << pmean << " " << fmean << endl;
}

//////////////////////////
// dont look below here //
//////////////////////////


/** 
    Store certain amount of track 
    arrays in a root file for later
    access - okay needs to many changes
    in AliHLTTrackArray and AliHLTTrack 
    classes so stopped working on it 
    (see fill_pp instead). 
*/

void store_pp(Int_t evs=0,Int_t eve=100,Char_t *path="./",Char_t *rfile="tracks.root")
{
  AliHLTLogger l;
  l.Set(AliHLTLogger::kAll);
  l.UseStdout();
  //l.UseStream();

  if(getenv("TRANSFORMCONFIGPATH")){
    AliHLTTransform::Init(getenv("TRANSFORMCONFIGPATH"));
  }

  TFile file(rfile,"UPDATE");
  Char_t name[1024];

  for(Int_t i=evs;i<eve;i++){
    AliHLTEvaluate *eval=new AliHLTEvaluate(path,63,63);
    eval->LoadData(i,-1);
    eval->AssignIDs();

    AliHLTTrackArray *a=eval->GetTracks();
    sprintf(name,"tracks-%d",i);
    a->Compress();
    a->Write(name);

    delete eval;
  }

  file.Close();
}

/** 
    Loop over stored track arrays 
    and print content. Does not properly
    work (see store_pp).
*/

void load_pp(Int_t evs=0,Int_t eve=100,Char_t *path="./")
{
  AliHLTLogger l;
  l.Set(AliHLTLogger::kAll);
  l.UseStdout();
  //l.UseStream();

  if(getenv("TRANSFORMCONFIGPATH")){
    AliHLTTransform::Init(getenv("TRANSFORMCONFIGPATH"));
  } else {
    AliHLTTransform::Init(path);
  }

  for(Int_t j=evs;j<eve;j++){

    AliHLTEvaluate *eval=new AliHLTEvaluate(path,63,63);
    eval->LoadData(j,-1); //load whole patch
    eval->AssignIDs();

    AliHLTTrackArray *a=eval->GetTracks();

    AliHLTTrack t;
    Float_t phi,pt,tgl,eta;
    Int_t id;
    Int_t nhits;

    AliHLTVertex vertex;
    Double_t xc,yc,zc,impact;

    Int_t nents=0;
    Int_t ntracks=a.GetNTracks();
    for(Int_t i=0;i<ntracks;i++){
      t.Set(a.GetTrack(i));

      t.CalculateHelix();
      t.GetClosestPoint(&vertex,xc,yc,zc);

      nhits=t.GetNHits();
      if(nhits<63) continue;
      pt=t.GetPt();
      tgl=t.GetTgl();
      phi=100;
      id=t.GetMCid();
      if(t.GetPy()!=0) phi=TMath::ATan(t.GetPy()/t.GetPx());

      impact=TMath::fabs(zc);
      //impact=TMath::Sqrt(xc*xc+yc*yc+zc*zc);

      cout << i << ": pt " << pt << " " << phi << " " << tgl << " " << id << " nhits " << nhits << " impact " << impact << endl;

      nents++;
    } 

    delete a;
  }
}
