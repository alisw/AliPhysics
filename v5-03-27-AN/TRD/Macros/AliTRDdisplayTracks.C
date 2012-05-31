void AliTRDdisplayTracks(Int_t track=157) 
{

  c1 = new TCanvas( "RecPoints", "RecPoints Display", 10, 10, 710, 740);
  c1->SetFillColor(1);
  TView *v=new TView(1);
  
      v->SetRange(-350.,-350.,-400.,350.,350.,400.); // full
  //  v->SetRange(0.,0.,0.,350.,350.,400.);  // top right
  //  v->SetRange(-350.,0.,0.,0.,350.,400.); // top left
  //  v->SetRange(0.,-350.,0.,350.,0.,400.); // bottom right 
  //  v->SetRange(-350.,-350.,0.,0.,0.,400.); // bottom left
  //  v->Side();
  v->Top();
  cerr<<"psi = "<<v->GetPsi()<<endl;

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
    cout << "Loaded shared libraries" << endl;
  }       


// Load tracks

   TFile *tf=TFile::Open("AliTRDtracks.root");
   if (!tf->IsOpen()) {cerr<<"Can't open AliTRDtracks.root !\n"; return 3;}
   TObjArray tarray(2000);
   //  tf->ls();
   TTree *tracktree=(TTree*)tf->Get("TreeT0_TRD");
   //   tracktree->ls();
   TBranch *tbranch=tracktree->GetBranch("tracks");
   Int_t nentr=tracktree->GetEntries();
   cerr<<"Found "<<nentr<<" entries in the track tree"<<endl;
   for (Int_t i=0; i<nentr; i++) {
       AliTRDtrack *iotrack=new AliTRDtrack;
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);
       tarray.AddLast(iotrack);
       cerr<<"got track "<<i<<": index ="<<iotrack->GetLabel()<<endl;
   }
   tf->Close();                 


  // Load clusters
  Char_t *alifile = "TRDclusters.root";
  Int_t   nEvent  = 0;
  TObjArray rparray(2000);
  TObjArray *RecPointsArray = &rparray;
  const TFile *geofile =TFile::Open(alifile);
  //  AliTRDtracker *Tracker = new AliTRDtracker("dummy","dummy");
  AliTRDtracker *Tracker = new AliTRDtracker(geofile);
   Tracker->ReadClusters(RecPointsArray,alifile); 
  Int_t nRecPoints = RecPointsArray->GetEntriesFast();
  cerr<<"Found "<<nRecPoints<<" rec. points"<<endl;


  // Connect the AliRoot file containing Geometry, Kine, Hits, and Digits
      alifile = "galice.root";

   TFile *gafl = (TFile*) gROOT->GetListOfFiles()->FindObject(alifile);
  if (!gafl) {
    cout << "Open the ALIROOT-file " << alifile << endl;
    gafl = new TFile(alifile);
  }
  else {
    cout << alifile << " is already open" << endl;
  }
  

  // Get AliRun object from file or create it if not on file
  gAlice = (AliRun*) gafl->Get("gAlice");
  if (gAlice)
    cout << "AliRun object found on file" << endl;
  else
    gAlice = new AliRun("gAlice","Alice test program");
	

       AliTRDparameter *par =  ( AliTRDparameter *par) gafl->Get("TRDparameter");
       AliTRDv1       *TRD           = (AliTRDv1*) gAlice->GetDetector("TRD");
       AliTRDgeometry *fGeom   = TRD->GetGeometry();   
       
  Int_t i,j,index,det,sector, ti[3];
  Double_t x,y,z, cs,sn,tmp;
  Float_t global[3], local[3];

  // Display all TRD RecPoints
  TPolyMarker3D *pm = new TPolyMarker3D(nRecPoints);

  for (Int_t i = 0; i < nRecPoints; i++) {
    printf("\r point %d out of %d",i,nRecPoints);
    AliTRDcluster *rp = (AliTRDcluster *) RecPointsArray->UncheckedAt(i);    
  
   Int_t  idet=rp->GetDetector();
   Int_t iplan = fGeom->GetPlane(idet);    
   Int_t itt=rp->GetLocalTimeBin();
   Float_t  timeSlice = itt+0.5; 
   Float_t  time0     = par->GetTime0(iplan);

  // calculate (x,y,z) position in rotated chamber
   local[0] = time0 - (timeSlice - par->GetTimeBefore()) 
         * par->GetTimeBinSize();
   local[1]=rp->GetY();
   local[2]=rp->GetZ(); 

    for (j = 0; j < 3; j++) { ti[j] = rp->GetTrackIndex(j); }

    if((track < 0) ||
       ((ti[0]==track)||(ti[1]==track)||(ti[2]==track))) {

        if(fGeom->RotateBack(idet,local,global)) {     
       	x=global[0]; y=global[1]; z=global[2];
       	pm->SetPoint(i,x,y,z);
         }

    }
  }

  pm->SetMarkerSize(1); pm->SetMarkerColor(10); pm->SetMarkerStyle(1);
  pm->Draw();
  
   AliTRDparameter *par =  ( AliTRDparameter *par) gafl->Get("TRDparameter");

  Int_t ntracks = tarray.GetEntriesFast();

  for (i = 0; i < ntracks; i++) {
    AliTRDtrack *pt = (AliTRDtrack *) tarray.UncheckedAt(i), &t=*pt;
    Int_t nclusters = t.GetNumberOfClusters();
    cerr<<"in track "<<i<<" found "<<nclusters<<" clusters"<<endl;
   
    TPolyMarker3D *pm = new TPolyMarker3D(nclusters);
    for(j = 0; j < nclusters; j++) {
      index = t.GetClusterIndex(j);
      AliTRDcluster *rp = (AliTRDcluster *) RecPointsArray->UncheckedAt(index); 

   Int_t  idet=rp->GetDetector();
   Int_t iplan = fGeom->GetPlane(idet);    
   Int_t itt=rp->GetLocalTimeBin();
   Float_t  timeSlice = itt+0.5; 
   Float_t  time0     = par->GetTime0(iplan);


  // calculate (x,y,z) position in rotated chamber
       local[0] = time0 - (timeSlice - par->GetTimeBefore()) 
         * par->GetTimeBinSize();
       local[1]=rp->GetY(); 
       local[2]=rp->GetZ();

            if(fGeom->RotateBack(idet,local,global)) {     
      	x=global[0]; y=global[1]; z=global[2];
      	pm->SetPoint(j,x,y,z);
       }
    } 
    pm->SetMarkerSize(1); pm->SetMarkerColor(i%6+3); pm->SetMarkerStyle(1);
    //    pm->Draw();
  }
  
  TGeometry *geom=(TGeometry*)gafl->Get("AliceGeom");
  geom->Draw("same");

  c1->Modified(); 
  
  c1->Update();
  
  gafl->Close();

}
