void  TestReconstruction (Int_t vol=1, const Int_t nRings=128, const Int_t nSectors=20) 
{
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) 
    {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
    }
  // Connect the Root Galice file containing Geometry, Kine and Hits
  char filename[]="galice.root";
  TFile *file =  (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
  if (!file) file = new TFile(filename); 
  // Get AliRun object from file or create it if not on file
  if (!gAlice) 
    {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("\nAliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }
  Int_t nbytes = 0;
  Int_t ipart;
  Int_t volume;
  Int_t np[nRings][nSectors];
 
  
  for (int i=0;i<nRings;i++)
    for(int j=0;j<nSectors;j++)
      np[i][j]=0;
 
  TH1F *hNReal = new TH1F("hNReal","Real number of particles",50,0,50);
  TH1F *hNRec  = new TH1F("hNRec ","Reconst. number of particles",50,0,50);

  Int_t nparticles = gAlice->GetEvent(0);
  if (nparticles <= 0) return;
  printf("\nnparticles=%d\n",nparticles);
   
  gAlice->TreeR()->GetEvent(0);   
  AliFMD *FMD  = (AliFMD*)gAlice->GetDetector("FMD");
  TClonesArray *Particles = gAlice->Particles();
  TParticle *particle;
  AliFMDhit *fmdHit;
  AliFMDReconstParticles *fmdRP;
  if (FMD) 
    {
      TClonesArray *FMDhits   = FMD->Hits();
      TClonesArray *FMDrec    = FMD->ReconParticles(); 
      TTree *TH = gAlice->TreeH();
      Int_t ntracks    = TH->GetEntries();
      if (ntracks<=0) return;

      Int_t  nPads=FMDrec->GetEntries();
       

#ifdef DEBUG
      cout<<"\n(AliFMDReconstParticles*)FMDrec->UncheckedAt(0)="<<(AliFMDReconstParticles*)FMDrec->UncheckedAt(0);	  	   
      cout<<"\nFMDrec->UncheckedAt(0)="<<FMDrec->UncheckedAt(0);	  	   
#endif
    
      for (Int_t track=0; track<ntracks;track++)
	{
	  gAlice->ResetHits();
	  nbytes += TH->GetEvent(track);//?
	  particle=(TParticle*)Particles->UncheckedAt(track);
	  //      Int_t numpart=particle->GetKF();
	  //Float_t eta=particle->GetEta();
	  
	  Int_t  nhits=FMDhits->GetEntriesFast();
	  for (Int_t hit=0;hit<nhits;hit++) 
	    {
	      fmdHit  = (AliFMDhit*)FMDhits->UncheckedAt(hit);
	      volume=fmdHit->Volume();
	      if(volume==vol)
		{
		  np[fmdHit->NumberOfRing()-1][fmdHit->NumberOfSector()-1]++;
		}  
	    }
	}
      //Int_t  nRecPart=FMDrec->GetEntriesFast();      
      Int_t nDeterm=0; Int_t nReal=0; 
      for (Int_t pad=0;pad<nPads;pad++) 
	{
	  fmdRP  = (AliFMDReconstParticles*)FMDrec->UncheckedAt(pad);
	  volume=fmdRP->GetVolume();
	  if(volume==vol)
	    {
;
#ifdef DEBUG
	      fmdDigit  = (AliFMDdigit*)FMDdig->UncheckedAt(pad);
	      cout<<"\nfmdDigit->ADCsignal()="<<fmdDigit->ADCsignal();
	      cout<<"\nfmdDigit->NumberOfRing()="<<fmdDigit->NumberOfRing();
	      cout<<"\nfmdDigit->NumberOfSector()="<<fmdDigit->NumberOfSector();    
#endif  
             nDeterm+=fmdRP->GetNumberOfReconstParticles();
	     nReal+=np[fmdRP->GetNumberOfRing()-1][fmdRP->GetNumberOfSector()-1]; //-1=?
	      Int_t RecRing=fmdRP->GetNumberOfRing()-1;
	      Int_t RecSector=fmdRP->GetNumberOfSector()-1;
	      hNReal->Fill(np[RecRing][RecSector]);
	      hNRec->Fill(fmdRP->GetNumberOfReconstParticles());
	    }
        }  
    }  
  cout<<"\nReal="<<nReal<<
    " nDeterm="<<nDeterm<<
    "\nerror="<<float(nDeterm-nReal)/float(nReal)<<endl;


  TCanvas *c1 = new TCanvas("c1","Alice FMD ",400,10,800,800);
  hNReal->SetFillColor(2);
  hNReal->Draw();
  hNRec->SetFillColor(4);
  hNRec->Draw("same");
  
}






