void debug(){
 //
 //  macro fro trigger analysis (F. Meddi suggestion)
 //
 Int_t evNumber=0;
 char *filename="galice.root";
 char *fileout="analyse.root";

 // Dynamically link some shared libs
 if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
 } 
   
 // Connect the Root Galice file containing Geometry, Kine and Hits
 TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
 if (!file) file = new TFile(filename);

 // Get AliRun object from file or create it if not on file
 if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
 }

 //loop over events
 for (int nev=0; nev<= evNumber; nev++) {
   Int_t nparticles = gAlice->GetEvent(nev);
   cout << "nev         " <<nev<<endl;
   cout << "nparticles  " <<nparticles<<endl;
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
   ITS->FillModules(nev,-1,nmodules," "," ");

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

  AliITSgeom *g = ((AliITS *)ITS)->GetITSgeom();
  Int_t lay, lad, det;
  printf("Starts loop on SPD detectors\n");


  //loop over the pixel detectors index=0-79     (1-20)*4 layer 1  
  //                              index=80-239   (1-40)*4 layer 2
//  for (Int_t index=g->GetStartSPD();index<=g->GetLastSPD();index++) 
    for (Int_t index=0;index<2;index++)  //debug
//    Int_t index=5;
  {
  
    g->GetModuleId(index,lay,lad,det); 
    printf("detector %d (lay=%d lad=%d det=%d)\n",index+1,lay,lad,det);

    AliITSmodule *itsModule = (AliITSmodule*) mods->At(index);
    Int_t numofhits = itsModule->GetNhits();
    printf("number of hits %d\n",numofhits);
    if(!numofhits) continue;

    //---------- starts test on digits
    ITS->ResetDigits();
    TD->GetEvent(index+1);
    Int_t ndigits = ITSdigits->GetEntriesFast();
    if (ndigits) printf("Found %d digits for module %d \n",ndigits,index+1);
    if (!ndigits) printf("no digits found \n");


    //loop on digits
    for (Int_t digit=0;digit<ndigits;digit++) {
        ITSdigit   = (AliITSdigitSPD*)ITSdigits->UncheckedAt(digit);
        printf("%d %d %d %d \n",ITSdigit->fCoord1,ITSdigit->fCoord2,ITSdigit->fSignal,ITSdigit->fTracks[0]);
     }
     cout<<"END  test for digits "<<endl;


    //---------- starts test on clusters
    ITS->ResetClusters();
    TC->GetEvent(index);
    Int_t nclust = ITSclusters->GetEntries();
    printf("number of clusters %d\n",nclust);


    //loop on clusters
    for (Int_t clu=0;clu<nclust;clu++)
    {
      itsclu = (AliITSRawClusterSPD*) ITSclusters->UncheckedAt(clu);
      printf("cluster %d nZ=%f nX=%f Z=%f X=%f\n",clu+1,itsclu->NclZ(),
                      itsclu->NclX(),itsclu->Z(),itsclu->X());
    }
     cout<<"END  test for clusters "<<endl;


    //---------- starts test on rec points
    ITS->ResetRecPoints();
    TR->GetEvent(index+1);
    Int_t nrecpoints = ITSrec->GetEntries();
    printf("Found %d recpoints for module %d \n",nrecpoints,index+1);

    //loop on rec points
    for (Int_t irec=0;irec<nrecpoints;irec++) {
         recp   = (AliITSRecPoint*)ITSrec->UncheckedAt(irec);
        printf("%d %f %f %f %f  %d %d %d\n",irec+1,recp->GetX(),recp->GetZ(),
            recp->fSigmaX2,recp->fSigmaZ2,
            recp->fTracks[0],recp->fTracks[1],recp->fTracks[2]);
    }

    printf("Detector No. %d (%d total hits) (%d digits) (%d clusters)\n",
                                         index+1,numofhits,ndigits,nclust);

 
  } //end loop on the SPD detectors

} // end loop over events

}
