void SPDclusterTest(Int_t evNumber1=0,Int_t evNumber2=0){
 //
 //  macro to monitor the SPD digitization and clusterization done with
 //  the Bari/Salerno model
 //
 //  R. Caliandro 15/05/2001 
 //
 //
 //--plots displayed:
 //
 //--pag1:  number of hits     per SPD detector (1-->250)
 //         number of hits     per SPD detector (1-->250)
 //         number of clusters per SPD detector (1-->250)
 //
 //--pag2:  r-phi cluster length layer 1 (red)
 //         z     cluster length layer 1 (red)
 //         r-phi cluster length layer 2 (blue)
 //         z     cluster length layer 2 (blue)
 //
 //--pag3:  r-phi resolution layer 1 (red)
 //         z     resolution layer 1 (red)
 //         r-phi resolution layer 2 (blue)
 //         z     resolution layer 2 (blue)
 //
 //--pag4:  Cluster shape analysis for clusters of 1, 2 and 3 digits
 //         zdim versus xdim for clusters of 4 digits
 //
 // input file name, digitized and clusterized
 char *filein="galice.root";
 // output file name, containing histograms
 char *fileout="SPD_his.root";
 // flag for debugging: 0=no debugging, 1=debugging
 Int_t debug=0;


 // Dynamically link some shared libs
 if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
 } else {
    delete gAlice;
    gAlice=0;
 }

   
 // Connect the Root Galice file containing Geometry, Kine and Hits
 TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filein);
 if (!file) file = new TFile(filein);

 // Get AliRun object from file or create it if not on file
 if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
 }

  //to get the segmentation pointer
  AliITS *ITS  = (AliITS*) gAlice->GetModule("ITS");
  AliITSDetType *iDetType=ITS->DetType(0);
  AliITSsegmentationSPD *seg=(AliITSsegmentationSPD*)iDetType->GetSegmentationModel();

//=======================================================
//--booking of ntuples 

//--ntuple for each detector
  TNtuple *ntuple2 = new TNtuple("ntuple2","","ndet:lay:lad:det:nhits:ndig:nclus");
//--ntuple for each cluster 
  TNtuple *ntuple = new TNtuple("ntuple","","ndet:iclus:ndigclus:xdim:zdim:xdiff:zdiff:anglex:anglez:pmom:errx:errz");

//--booking of histograms 
//layer 1
  TH1F *hist1n1 = new TH1F("hist1n1","xdim",15,0.5,15.5);
  TH1F *hist2n1 = new TH1F("hist2n1","zdim",10,0.5,10.5);
  TH1F *hist3n1 = new TH1F("hist3n1","dig/clus",20,0.5,20.5);
  TH1F *hist4n1 = new TH1F("hist4n1","errx",100,0,0.01);
  TH1F *hist5n1 = new TH1F("hist5n1","errz",500,0,0.05);
  TH2F *hist7n1 = new TH2F("hist7n1","xdim:delx",80,0,800.,15,0.5,15.5);
  TH2F *hist8n1 = new TH2F("hist8n1","zdim:delz",180,0,1800.,10,0.5,10.5);
//layer 2
  TH1F *hist1n2 = new TH1F("hist1n2","xdim",15,0.5,15.5);
  TH1F *hist2n2 = new TH1F("hist2n2","zdim",10,0.5,10.5);
  TH1F *hist3n2 = new TH1F("hist3n2","dig/clus",20,0.5,20.5);
  TH1F *hist4n2 = new TH1F("hist4n2","errx",100,0,0.01);
  TH1F *hist5n2 = new TH1F("hist5n2","errz",500,0,0.05);
  TH2F *hist7n2 = new TH2F("hist7n2","xdim:delx",80,0,800.,15,0.5,15.5);
  TH2F *hist8n2 = new TH2F("hist8n2","zdim:delz",180,0,1800.,10,0.5,10.5);
//--resolution 
  TH1F *hist1 = new TH1F("hist1","xdiff",200,-100,100);
  TH1F *hist3 = new TH1F("hist3","xdiff",200,-100,100);
  TH1F *hist2 = new TH1F("hist2","zdiff",170,-850,850);
  TH1F *hist4 = new TH1F("hist4","zdiff",170,-850,850);
//--momentum 
  TH1F *hist5 = new TH1F("hist5","pmom",200,0,2000);
//--rapidity 
  TH1F *hist6 = new TH1F("hist6","rapidity",60,-3,3);
  TH1F *hist6b= new TH1F("hist6b","rapidity - charged tracks",60,-3,3);
  TH1F *hist6b1= new TH1F("hist6b1","rapidity - charged tracks SPD",60,-3,3);
//--pseudo-rapidity
  TH1F *hist6p = new TH1F("hist6p","eta - charged tracks ",60,-3,3);
  TH1F *hist6p1 = new TH1F("hist6p1","eta - charged tracks SPD ",60,-3,3);
  TH1F *hist6p2 = new TH1F("hist6p2","eta - charged tracks SPD 2 ",60,-3,3);
//--resolution vs angle
  TH1F *hist11n1=new TH1F("hist11n1","anglex - layer 1",180,-90,90);
  TH1F *hist11n2=new TH1F("hist11n2","anglex - layer 2",180,-90,90);
  TH1F *hist12n1=new TH1F("hist12n1","anglez - layer 1",360,-180,180);
  TH1F *hist12n2=new TH1F("hist12n2","anglez - layer 2",360,-180,180);
  TH2F *hist13n1=new TH2F("hist13n1","xidff:anglex",20,-15,15,200,-100,100);
  TH2F *hist13n2=new TH2F("hist13n2","xidff:anglex",20,-30,-15,200,-100,100);
  TH2F *hist14n1=new TH2F("hist14n1","zidff:anglez",18,-90,90,170,-850,850);
  TH2F *hist14n2=new TH2F("hist14n2","zidff:anglez",18,-90,90,170,-850,850);
//--histograms for cluster shape analysis
  TH1F *histsp1=new TH1F("histsp1","Cluster shape (1)",10,0.5,10.5);
  TH2F *histsp2=new TH2F("histsp2","Cluster shape (2)",5,0.5,5.5,5,0.5,5.5);
//=======================================================

 //loop over events
 for (int nev=0; nev<= evNumber2; nev++) {
   Int_t nparticles = gAlice->GetEvent(nev);
   cout << "nev         " <<nev<<endl;
   cout << "nparticles  " <<nparticles<<endl;
   if (nev < evNumber1) continue;
   if (nparticles <= 0) return;

   TTree *TH        = gAlice->TreeH();
   Int_t ntracks    = TH->GetEntries();
   cout << "ntracks  " <<ntracks<<endl;

   // Get pointers to Alice detectors and Digit containers
   AliITS *ITS  = (AliITS *)gAlice->GetModule("ITS");
   TClonesArray *Particles = gAlice->Particles();
   if(!ITS) return;

   // fill modules with sorted by module hits
   Int_t nmodules;
   ITS->InitModules(-1,nmodules);
   ITS->FillModules(nev,evNumber2,nmodules," "," ");

   // get pointer to modules array
   TObjArray *mods = ITS->GetModules();
   AliITShit *itsHit;


   //get the Tree for clusters
   ITS->GetTreeC(nev);
   TTree *TC=ITS->TreeC();
   //TC->Print();
   Int_t nent=TC->GetEntries();
   printf("Found %d entries in the tree of clusters)\n",nent);
   TClonesArray *ITSclusters = ITS->ClustersAddress(0);
   printf("ITSclusters %p\n",ITSclusters);

   //get the Tree for digits
   TTree *TD = gAlice->TreeD();
   //TD->Print();
   Int_t nentd=TD->GetEntries();
   printf("Found %d entries in the tree of digits)\n",nentd);
   TObjArray *fBranches=TD->GetListOfBranches();
   TBranch *branch = (TBranch*)fBranches->UncheckedAt(0);
   printf ("branch %p entries %d \n",branch,branch->GetEntries());
   TClonesArray *ITSdigits  = ITS->DigitsAddress(0);
   printf ("ITSdigits %p \n",ITSdigits);

   //get the Tree for rec points
   TTree *TR = gAlice->TreeR();
   //TR->Print();
   Int_t nentr=TR->GetEntries();
   printf("Found %d entries in the tree of rec points)\n",nentr);
   TClonesArray *ITSrec  = ITS->RecPoints();
   printf ("ITSrec %p \n",ITSrec);
   AliITSRecPoint  *recp;

  // calculus of rapidity distribution for the generated tracks
  gAlice-> ResetHits();
  TParticle *particle;
  for (Int_t track=0; track<ntracks; track++)
  {
    particle  = (TParticle*)gAlice->Particle(track);
    Int_t ikparen   = particle -> GetFirstMother();
    Double_t charge = particle -> GetPDG() ->Charge();
    charge = charge/3.;  //charge is multiplied by 3 in PDG
    Double_t mass   = particle -> GetPDG() -> Mass();
    Double_t eta    = particle -> Eta();
    Int_t pdgcode   = particle -> GetPdgCode();
    char* title     = particle -> GetTitle();
    if (ikparen<0)
    {   
      Double_t part_ene = particle->Energy();
      Double_t part_pz  = particle->Pz();
      Double_t rapid;
      if (part_ene != part_pz) 
      {
        rapid=0.5*TMath::Log((part_ene+part_pz)/(part_ene-part_pz));
      }
      else {
        rapid = 1.e30;
      }
      // filling of the rapidity histogram
      hist6->Fill( (Float_t) rapid);
      if( charge != 0 ) {
        hist6b->Fill( (Float_t) rapid);
        hist6p->Fill( (Float_t) eta);
      }
//          printf("charge= %f, mass = %f , pdg= %d, title = %s\n",
//                    charge,mass,pdgcode,title);
    }
  }

  AliITSgeom *g = ((AliITS *)ITS)->GetITSgeom(); 
  Int_t lay, lad, det;
  //printf("Starts loop on SPD detectors\n");


  //loop over the pixel detectors index=0-79     (1-20)*4 layer 1  
  //                              index=80-239   (1-40)*4 layer 2
  for (Int_t index=g->GetStartSPD();index<=g->GetLastSPD();index++) 
//  for (Int_t index=g->GetStartSPD();index<1;index++)  //debug
  {
  
    g->GetModuleId(index,lay,lad,det); 
    //printf("detector %d (lay=%d lad=%d det=%d)\n",index+1,lay,lad,det);

    AliITSmodule *itsModule = (AliITSmodule*) mods->At(index);
    Int_t numofhits = itsModule->GetNhits();
    //printf("number of hits %d\n",numofhits);
    if(!numofhits) continue;

    //---------- starts test on digits
    ITS->ResetDigits();
    TD->GetEvent(index);
    Int_t ndigits = ITSdigits->GetEntriesFast();
    //if (ndigits) printf("Found %d digits for module %d \n",ndigits,index+1);
    if (!ndigits) printf("no digits found \n");



   if(debug==1) {
    //loop on digits
    for (Int_t digit=0;digit<ndigits;digit++) {
        ITSdigit   = (AliITSdigitSPD*)ITSdigits->UncheckedAt(digit);
        printf("digit=%d fCoord1=%d FCoord2=%d fSignal=%d fTracks=%d fHits=%d \n",digit,ITSdigit->fCoord1,ITSdigit->fCoord2,ITSdigit->fSignal,ITSdigit->fTracks[0],ITSdigit->fHits[0]);
     }
     cout<<"END  test for digits "<<endl;
    }


    //---------- starts test on clusters
    ITS->ResetClusters();
    TC->GetEvent(index);
    Int_t nclust = ITSclusters->GetEntries();
    //printf("Found %d clusters \n",nclust);
    if (!nclust) printf("no clusters found \n");


   if(debug==1) {
     //loop on clusters
     for (Int_t clu=0;clu<nclust;clu++)
     {
      //itsclu = (AliITSRawClusterSPD*) ITSclusters->UncheckedAt(clu);
      itsclu = (AliITSRawClusterSPD*) ITSclusters->At(clu);
      printf("cluster %d nZ=%f nX=%f Q=%f Z=%f X=%f\n",clu+1,itsclu->NclZ(),
                      itsclu->NclX(),itsclu->Q(),itsclu->Z(),itsclu->X());
     }
     cout<<"END  test for clusters "<<endl;
    }



    //---------- starts test on rec points
    ITS->ResetRecPoints();
    TR->GetEvent(index);
    Int_t nrecpoints = ITSrec->GetEntries();
    //printf("Found %d recpoints for module %d \n",nrecpoints,index+1);
    if (!nrecpoints) printf("no recpoints found \n");

   if(debug==1) {
    //loop on rec points
    for (Int_t irec=0;irec<nrecpoints;irec++) {
         recp   = (AliITSRecPoint*)ITSrec->UncheckedAt(irec);
        printf("%d %f %f %f %f  %d %d %d\n",irec+1,recp->GetX(),recp->GetZ(),
            recp->fSigmaX2,recp->fSigmaZ2,
            recp->fTracks[0],recp->fTracks[1],recp->fTracks[2]);
    }
   }

    printf("Detector n.%d (%d hits) (%d digits) (%d clusters)\n",
                                         index+1,numofhits,ndigits,nclust);

           // fill ntuple2
           ntuple2->Fill (   (Float_t) index+1,
		                     (Float_t) lay,
		                     (Float_t) lad,
		                     (Float_t) det,
		                     (Float_t) numofhits,
		                     (Float_t) ndigits,
		                     (Float_t) nclust);

    Int_t xlow; 
    Int_t zlow; 
    Int_t xhigh; 
    Int_t zhigh; 
    Int_t colcenter;
    Int_t rowcenter;

    // loop on clusters in each detector
    for (Int_t i=0; i<nclust; i++)
    {

       irawclu = (AliITSRawClusterSPD*) ITSclusters->UncheckedAt(i);
       irecp   = (AliITSRecPoint*)ITSrec->UncheckedAt(i);

       Int_t xdim = irawclu->NclX();
       Int_t zdim = irawclu->NclZ();
       Float_t errx = TMath::Sqrt(irecp->fSigmaX2);
       Float_t errz = TMath::Sqrt(irecp->fSigmaZ2);
       Float_t xcenter = irawclu->X();
       Float_t zcenter = irawclu->Z();
       Float_t ndigclus = irawclu->Q();
       Int_t itrackclus  = irecp->fTracks[0];

     //Find the hits associated to the main track of the cluster
     // loop on hits in the detector
     for (Int_t hit=0; hit<numofhits; hit++)
     {
          AliITShit *itsHit  = (AliITShit*)itsModule->GetHit(hit);
          Int_t itrackhit = itsHit->GetTrack();
          //Take the same track index of the main track of the cluster
          if (itrackhit == itrackclus) {
               if (itsHit->GetTrackStatus()==66) {
                   Float_t x1l = 10000*itsHit->GetXL(); //in microns
                   Float_t y1l = 10000*itsHit->GetYL();
                   Float_t z1l = 10000*itsHit->GetZL();
                   Float_t p1x = 1000*itsHit->GetPXL(); //in MeV/c
                   Float_t p1y = 1000*itsHit->GetPYL();
                   Float_t p1z = 1000*itsHit->GetPZL();
                }
                else {
                   Float_t x2l = 10000*itsHit->GetXL();
                   Float_t y2l = 10000*itsHit->GetYL();
                   Float_t z2l = 10000*itsHit->GetZL();
                   Float_t p2x = 1000*itsHit->GetPXL();
                   Float_t p2y = 1000*itsHit->GetPYL();
                   Float_t p2z = 1000*itsHit->GetPZL();


                }
          }
     }// end loop on hits on detector

       
     Float_t pmom=TMath::Sqrt(p1x*p1x+p1y*p1y+p1z*p1z); 
     hist5->Fill(pmom);

     Float_t dxhit = TMath::Abs(x2l-x1l);
     Float_t dzhit = TMath::Abs(z2l-z1l);

     Float_t xmidhit = (x1l + x2l)/2;
     Float_t zmidhit = (z1l + z2l)/2;

//   printf("cluster n.%d: x=%f z=%f\n",i,xcenter,zcenter);
//   printf("track n.%d: x1=%f x2=%f z1=%f z2=%f\n",itrackclus,
//                 x1l, x2l, z1l, z2l);

     // analysis of resolution vs angle
     if(index<80)
     {
           Float_t px = -p1x;
           Float_t py = -p1y;
     }
     else{
           Float_t px = p1x;
           Float_t py = p1y;
     }
     Float_t pz = p1z;
     // anglex is the angle in xy plane (local frame)
     Float_t anglex = atan2(px,py); 
     // anglez is the angle in zy plane (local frame)
     Float_t anglez = atan2(pz,py); 
     anglex *= 180.0/TMath::Pi(); // degrees
     anglez *= 180.0/TMath::Pi(); // degrees

     if(xmidhit != 0  || zmidhit != 0)
     {
          Float_t xdiff = (xcenter - xmidhit);
          Float_t zdiff = (zcenter - zmidhit);

          if(index<80)
          {
             // resolution plots
             hist1->Fill(xdiff); 
             hist2->Fill(zdiff); 

             // plots of resolution vs angle
             hist11n1->Fill(anglex);
             hist12n1->Fill(anglez);
             hist13n1->Fill(anglex,xdiff);
             hist14n1->Fill(anglez,zdiff);

          } else {

             // resolution plots
             hist3->Fill(xdiff); 
             hist4->Fill(zdiff); 

             // plots of resolution vs angle
             hist11n2->Fill(anglex);
             hist12n2->Fill(anglez);
             hist13n2->Fill(anglex,xdiff);
             hist14n2->Fill(anglez,zdiff);

          }
       } 
	  
       // fill the ntuple
       ntuple->Fill ( (Float_t) index,
          	          (Float_t) i,
          	          (Float_t) ndigclus,
		              (Float_t) xdim,
		              (Float_t) zdim,
		              (Float_t) xdiff,
		              (Float_t) zdiff,
		              (Float_t) anglex,
		              (Float_t) anglez,
		              (Float_t) pmom,
                                errx,
                                errz);

       // other histograms
       if(index<80)
       {
          hist1n1->Fill((Float_t) xdim);
          hist2n1->Fill((Float_t) zdim);
          hist3n1->Fill(ndigclus);
          hist4n1->Fill(errx);
          hist5n1->Fill(errz);
          hist7n1->Fill(dxhit,(Float_t) xdim);
          hist8n1->Fill(dzhit,(Float_t) zdim);

       } else {
          hist1n2->Fill((Float_t) xdim);
          hist2n2->Fill((Float_t) zdim);
          hist3n2->Fill(ndigclus);
          hist4n2->Fill(errx);
          hist5n2->Fill(errz);
          hist7n2->Fill(dxhit,(Float_t) xdim);
          hist8n2->Fill(dzhit,(Float_t) zdim);
       }

       //histograms for cluster shape analysis
       Int_t xx;
       if(ndigclus<=3) {
          if(ndigclus==1) {
             xx=1;
          } elseif(ndigclus==2){
             if(zdim==2 && xdim==1) xx=2;
             if(zdim==1 && xdim==2) xx=3;
             if(zdim==2 && xdim==2) xx=4;
          } elseif(ndigclus==3){
             if(zdim==2 && xdim==2) xx=5;
             if(zdim==3 && xdim==1) xx=6;
             if(zdim==1 && xdim==3) xx=7;
             if(zdim==3 && xdim==3) xx=8;
             if(zdim==3 && xdim==2) xx=9;
             if(zdim==2 && xdim==3) xx=10;
          }
          histsp1->Fill((Float_t) xx);
       } elseif(ndigclus==4){
          histsp2->Fill((Float_t) xdim,(Float_t) zdim);
       }


    }//end loop on clusters

 } //end loop on the SPD detectors

} //end loop on events 


//================== Plot the results ===================================

gROOT->Reset();
gStyle->SetFillColor(0);
gStyle->SetStatW(0.37);
gStyle->SetStatH(0.22);

//----------------------------------------------------- page 1
gStyle->SetOptLogy(0);
gStyle->SetOptStat(1100);

  TCanvas *c1 = new TCanvas("c1","hits, digits, clusters",200,10,700,780);
c1->SetFillColor(0);
c1->Divide(1,3);
    c1->cd(1);
      ntuple2->SetMarkerColor(kRed);
      ntuple2->SetMarkerStyle(20);
      ntuple2->SetMarkerSize(0.6);
      ntuple2->Draw("nhits:ndet","");
    c1->cd(2);
      ntuple2->SetMarkerColor(kBlue);
      ntuple2->Draw("ndig:ndet","");
    c1->cd(3);
      ntuple2->SetMarkerColor(6);
      ntuple2->Draw("nclus:ndet","");

//----------------------------------------------------- page 2

  TCanvas *c2 = new TCanvas("c2","Cluster Lengths",200,10,700,780);
   //
   // Inside this canvas, we create 4 pads
   //
   pad1 = new TPad("pad1","xdim layer1"     ,0.01,0.51,0.49,0.99,10);
   pad2 = new TPad("pad2","zdim layer1"     ,0.51,0.51,0.99,0.99,10);
   pad3 = new TPad("pad3","xdim layer2"    ,0.01,0.01,0.49,0.49,10);
   pad4 = new TPad("pad4","zdim layer2"    ,0.51,0.01,0.99,0.49,10);
   pad1->Draw();
   pad2->Draw();
   pad3->Draw();
   pad4->Draw();
//
   gStyle->SetStatW(0.40);
   gStyle->SetStatH(0.20);
   gStyle->SetStatColor(42);
   gStyle->SetOptStat(111110);
//  gStyle->SetOptFit(1);
  
   pad1->cd();
   pad1->SetLogy();
   hist1n1->Draw();
   hist1n1->SetLineWidth(2);
   hist1n1->SetLineColor(kRed);
   hist1n1->GetXaxis()->SetNdivisions(110);
   hist1n1->GetYaxis()->SetLabelSize(0.06);
   c2->Update();
   //
   pad2->cd();
   pad2->SetLogy();
   hist2n1->Draw();
   hist2n1->SetLineWidth(2);
   hist2n1->SetLineColor(kRed);
   hist2n1->GetXaxis()->SetNdivisions(110);
   hist2n1->GetYaxis()->SetLabelSize(0.06);
   c2->Update();
   //
   pad3->cd();
   pad3->SetLogy();
   hist1n2->Draw();
   hist1n2->SetLineWidth(2);
   hist1n2->SetLineColor(kBlue);
   hist1n2->GetXaxis()->SetNdivisions(110);
   hist1n2->GetYaxis()->SetLabelSize(0.06);
   c2->Update();
   //
   pad4->cd();
   pad4->SetLogy();
   hist2n2->Draw();
   hist2n2->SetLineColor(kBlue);
   hist2n2->SetLineWidth(2);
   hist2n2->GetXaxis()->SetNdivisions(110);
   hist2n2->GetYaxis()->SetLabelSize(0.06);
   c2->Update();
   //
//----------------------------------------------------- page 3

  TCanvas *c3 = new TCanvas("c3","Resolutions",200,10,700,780);
   //
   // Inside this canvas, we create 4 pads
   //
   pad1 = new TPad("pad1","xdiff layer1"     ,0.01,0.51,0.49,0.99,10);
   pad2 = new TPad("pad2","zdiff layer1"     ,0.51,0.51,0.99,0.99,10);
   pad3 = new TPad("pad3","xdiff layer2"    ,0.01,0.01,0.49,0.49,10);
   pad4 = new TPad("pad4","zdiff layer2"    ,0.51,0.01,0.99,0.49,10);
   pad1->Draw();
   pad2->Draw();
   pad3->Draw();
   pad4->Draw();
//
   gStyle->SetStatW(0.20);
   gStyle->SetStatH(0.20);
   gStyle->SetStatColor(42);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1);
  
   pad1->cd();
   hist1->Draw();
   hist1->SetLineColor(kRed);
   hist1->Fit("gaus");
   c3->Update();
   //
   pad2->cd();
   hist2->Draw();
   hist2->SetLineColor(kRed);
   hist2->Fit("gaus");
   c3->Update();
   //
   pad3->cd();
   hist3->Draw();
   hist3->SetLineColor(kBlue);
   hist3->Fit("gaus");
   c3->Update();
   //
   pad4->cd();
   hist4->Draw();
   hist4->SetLineColor(kBlue);
   hist4->Fit("gaus");
   c3->Update();

//----------------------------------------------------- page 4
  TCanvas *c4 = new TCanvas("c4","Cluster Shape Analysis",200,10,700,780);
//
gStyle->SetOptStat(0);
c4->SetFillColor(0);
c4->Divide(1,2);

    c4->cd(1);
    c4_1->SetLogy();
      histsp1->Draw();
    c4->cd(2);
    c4_2->Divide(2,1);
    c4_2->cd(1);
      histsp2->Draw("box");
    c4_2->cd(2);
     clustershape();




//================== Store the histograms ===================================

  //to write the histograms into a .root file
  TFile outfile(fileout,"RECREATE");

  ntuple->Write();
  ntuple2->Write();
  hist1n1->Write();
  hist2n1->Write();
  hist3n1->Write();
  hist4n1->Write();
  hist5n1->Write();
  hist7n1->Write();
  hist8n1->Write();
  hist1n2->Write();
  hist2n2->Write();
  hist3n2->Write();
  hist4n2->Write();
  hist5n2->Write();
  hist7n2->Write();
  hist8n2->Write();
  hist1->Write();
  hist2->Write();
  hist3->Write();
  hist4->Write();
  hist5->Write();
  hist6->Write();
  hist6b->Write();
  hist6p->Write();
  hist6b1->Write();
  hist6p1->Write();
  hist6p2->Write();
  hist11n1->Write();
  hist11n2->Write();
  hist12n1->Write();
  hist12n2->Write();
  hist13n1->Write();
  hist13n2->Write();
  hist14n1->Write();
  hist14n2->Write();
  histsp1->Write();
  histsp2->Write();

  outfile->Close();
}
//-----------------------------------------------------------------



void clustershape(){
 //
 //macro to display the legend of the cluster shape analysis plot
 //
  
   TPad *pad1 = new TPad("pad1", "This is pad2",0,0,1,1);
   pad1->Draw();
   pad1->cd();
   pad1->Range(0,0.25,1,1);
   pad1->SetFillColor(0);
   pad1->SetBorderSize(1);

//------------------------------------------
   Float_t yfirst= 0.95;
   Float_t ysh   = 0.05;
   Float_t ysiz  = 0.02;
   Float_t ynum  = 0.005;
//------------------------------------------

   //bin 1
   TLatex *tex = new TLatex(0.12,yfirst,"1");
   tex->SetTextSize(0.07);
   tex->SetLineWidth(2);
   tex->Draw();
   TPave *pave = new TPave(0.3,yfirst,0.5,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();

   //bin 2
   yfirst=yfirst-ysh;
   TLatex *tex = new TLatex(0.12,yfirst,"2");
   tex->SetTextSize(0.07);
   tex->SetLineWidth(2);
   tex->Draw();
   TPave *pave = new TPave(0.3,yfirst,0.5,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.5,yfirst,0.7,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();

   //bin 3
   yfirst=yfirst-ysh;
   TLatex *tex = new TLatex(0.12,yfirst-ynum,"3");
   tex->SetTextSize(0.07);
   tex->SetLineWidth(2);
   tex->Draw();
   TPave *pave = new TPave(0.3,yfirst,0.5,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.3,yfirst,0.5,yfirst-ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();

   //bin 4
   yfirst=yfirst-1.8*ysh;
   TLatex *tex = new TLatex(0.12,yfirst+3*ynum,"4");
   tex->SetTextSize(0.07);
   tex->SetLineWidth(2);
   tex->Draw();
   TPave *pave = new TPave(0.3,yfirst,0.5,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.5,yfirst+ysiz,0.7,yfirst+2*ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();

   //bin 5
   yfirst=yfirst-1.5*ysh;
   TLatex *tex = new TLatex(0.12,yfirst+3*ynum,"5");
   tex->SetTextSize(0.07);
   tex->SetLineWidth(2);
   tex->Draw();
   TPave *pave = new TPave(0.3,yfirst,0.5,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.5,yfirst,0.7,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.5,yfirst+ysiz,0.7,yfirst+2*ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();

   //bin 6
   yfirst=yfirst-1.5*ysh;
   TLatex *tex = new TLatex(0.12,yfirst+ynum,"6");
   tex->SetTextSize(0.07);
   tex->SetLineWidth(2);
   tex->Draw();
   TPave *pave = new TPave(0.3,yfirst,0.5,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.5,yfirst,0.7,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.7,yfirst,0.9,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();

   //bin 7
   yfirst=yfirst-2*ysh;
   TLatex *tex = new TLatex(0.12,yfirst+ysiz+ynum,"7");
   tex->SetTextSize(0.07);
   tex->SetLineWidth(2);
   tex->Draw();
   TPave *pave = new TPave(0.3,yfirst,0.5,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.3,yfirst+ysiz,0.5,yfirst+2*ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.3,yfirst+2*ysiz,0.5,yfirst+3*ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();

   //bin 8
   yfirst=yfirst-1.5*ysh;
   TLatex *tex = new TLatex(0.12,yfirst+ynum,"8");
   tex->SetTextSize(0.07);
   tex->SetLineWidth(2);
   tex->Draw();
   TPave *pave = new TPave(0.3,yfirst,0.5,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.5,yfirst+ysiz,0.7,yfirst+2*ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.7,yfirst+2*ysiz,0.9,yfirst+3*ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();

   //bin 9
   yfirst=yfirst-1.5*ysh;
   TLatex *tex = new TLatex(0.12,yfirst+ynum,"9");
   tex->SetTextSize(0.07);
   tex->SetLineWidth(2);
   tex->Draw();
   TPave *pave = new TPave(0.3,yfirst,0.5,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.5,yfirst,0.7,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.7,yfirst+ysiz,0.9,yfirst+2*ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();

   //bin 10
   yfirst=yfirst-1.5*ysh;
   TLatex *tex = new TLatex(0.12,yfirst-ynum,"10");
   tex->SetTextSize(0.07);
   tex->SetLineWidth(2);
   tex->Draw();
   TPave *pave = new TPave(0.3,yfirst-ysiz,0.5,yfirst,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.3,yfirst,0.5,yfirst+ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
   TPave *pave = new TPave(0.5,yfirst+ysiz,0.7,yfirst+2*ysiz,1,"br");
   pave->SetFillColor(18);
   pave->SetLineWidth(2);
   pave->Draw();
}

