struct GoodTrack {
  Int_t lab;
  Int_t code;
  Bool_t flag;
  Float_t px,py,pz;
  Float_t x,y,z;
  Float_t pxg,pyg,pzg,ptg;
};


Int_t ITStracks(Int_t evNumber1=0,Int_t evNumber2=0,Int_t nclust=5) {

  cerr<<"Select tracks which have nclust rec points in ITS...\n";

  GoodTrack gt[15000];
  Int_t ngood=0;
  ifstream in("good_tracks");

  cerr<<"Reading good tracks...\n";
  while (in>>gt[ngood].lab>>gt[ngood].code
  >>gt[ngood].px >>gt[ngood].py>>gt[ngood].pz
  >>gt[ngood].x  >>gt[ngood].y >>gt[ngood].z
  >>gt[ngood].pxg  >>gt[ngood].pyg >>gt[ngood].pzg
  >>gt[ngood].ptg >>gt[ngood].flag) {
    ngood++;
    cerr<<ngood<<'\r';
    if (ngood==15000) {
      cerr<<"Too many good tracks !\n";
      break;
    }
  }
  if (!in.eof()) cerr<<"Read error (good_tracks) !\n";


  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } else {
    delete gAlice;
    gAlice=0;
  }

// Connect the Root Galice file containing Geometry, Kine and Hits
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
  if (!file) file = new TFile("galice.root","UPDATE");
  //if (!file) file = new TFile(filename);

// Get AliRun object from file or create it if not on file
  if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
  }


  AliITS* ITS =(AliITS *)gAlice->GetDetector("ITS");
  if (!ITS) return;

  //TClonesArray  *recPoints = ITS->RecPoints();
  //cout<<" recPoints = "<<recPoints<<"\n";
  Int_t numbpoints=0;
  AliITSRecPoint *recp=0;
  Int_t track=-3;


  for (int nev=evNumber1; nev<= evNumber2; nev++) {
    Int_t nparticles = gAlice->GetEvent(nev);
    cout << "nev         " << nev <<endl;
    cout << "nparticles  " << nparticles <<endl;
    if (nev < evNumber1) continue;
    if (nparticles <= 0) return; 
    TClonesArray  *recPoints = ITS->RecPoints();
   TObjArray *particles=gAlice->Particles();     
   // TClonesArray *particles=gAlice->Particles(); 
   // Int_t np=particles->GetEntriesFast();
   Int_t np=gAlice->GetNtrack(); //FCA correction      
    TMatrix goodITS(np,6);
    Int_t modmin[]={0,80,240,324,500,1282};
    Int_t modmax[]={79,239,323,499,1281,2269};
    for(Int_t j=0; j<np; j++){
       for(Int_t j1=0; j1<6; j1++) goodITS(j,j1)=0; 
    }
       
       
    Int_t nent=gAlice->TreeR()->GetEntries();
    printf("Found %d entries in the TreeR (must be one per module per event!)\n",nent); //

    for (Int_t layer=1;layer<=6;layer++) {
      for (Int_t mod=modmin[layer-1];mod<=modmax[layer-1];mod++) {
	     ITS->ResetRecPoints();
	     //gAlice->TreeR()->GetEvent(mod+1); //first entry in TreeR is empty
	     gAlice->TreeR()->GetEvent(mod);   //modificato il 2-3-2001
	     numbpoints = recPoints->GetEntries();
	     if(!numbpoints) continue;
	     for (Int_t irec=0;irec<numbpoints;irec++) {
                 recp=(AliITSRecPoint*)recPoints->UncheckedAt(irec);
		 for (Int_t it=0;it<3;it++) {
		     track=recp->GetLabel(it);
		     if(track<0) continue;
		     if(track > np) { 
		       cout<<" Error on track number \n";
		       getchar();
		     }
		     goodITS(track,layer-1)=1;
	         } //loop over tracks
	     } //loop over points
      } //loop over modules   
    } //loop over layers

    for(Int_t i=0; i<ngood; i++){
       if(gt[i].lab>np) {cout<<" Error on TPC good tracks label \n"; getchar();}    	   
       Int_t  ilab=gt[i].lab;
       Int_t totclust=0;
       for(Int_t j=0; j<6; j++) totclust+=goodITS(ilab,j);	   
	   if(totclust>=nclust) {
	     gt[i].flag=1;
	     printf("label nclusters %d %d\n",gt[i].lab,totclust);  //
	   }
    }
	   
    Int_t itsngood=0;
    ofstream out("itsgood_tracks");
    if (out) {
      for (Int_t ngd=0; ngd<ngood; ngd++) {
	     if(gt[ngd].flag) {
	       out<<gt[ngd].lab<<' '<<gt[ngd].code<<' '
		    <<gt[ngd].px <<' '<<gt[ngd].py<<' '<<gt[ngd].pz<<' '
		    <<gt[ngd].x  <<' '<<gt[ngd].y <<' '<<gt[ngd].z <<' '
		    <<gt[ngd].pxg  <<' '<<gt[ngd].pyg <<' '<<gt[ngd].pzg <<' '
		    <<gt[ngd].ptg<<gt[ngd].flag<<endl;
	       itsngood++;
	     }
      }
    } else cerr<<"Can not open file (itsgood_tracks) !\n";
    out.close();

    cerr<<"Number of good tracks in TPC: "<<ngood<<endl;
    cerr<<"Number of good tracks in ITS: "<<itsngood<<endl;

  }   // event loop 
   
  file->Close();   

  printf("after file close\n");
  
  return 0;
  
}

