//
// This macro reads the Hits Tree
//
void AliPMDHitsRead(Int_t nevt = 1)
{
  
  TStopwatch timer;
  timer.Start();
  TH2F *h2 = new TH2F("h2"," Y vs. X",200,-100.,100.,200,-100.,100.); 
  //  FILE *fpw = fopen("alipmdhits.dat","w");

  AliRunLoader *fRunLoader = AliRunLoader::Open("galice.root");

  if (!fRunLoader)
   {
     printf("Can not open session for file ");
   }
  
  if (!fRunLoader->GetAliRun()) fRunLoader->LoadgAlice();
  if (!fRunLoader->TreeE()) fRunLoader->LoadHeader();
  if (!fRunLoader->TreeK()) fRunLoader->LoadKinematics();

  gAlice = fRunLoader->GetAliRun();
  
  if (gAlice)
    {
      printf("Alirun object found\n");
    }
  else
    {
      printf("Could not found Alirun object\n");
    }
  
  fPMD  = (AliPMD*)gAlice->GetDetector("PMD");
  fPMDLoader = fRunLoader->GetLoader("PMDLoader");
  if (fPMDLoader == 0x0)
    {
      printf("Can not find PMDLoader\n");
    }


  fPMDLoader->LoadHits("READ");

  // This reads the PMD Hits tree and assigns the right track number
  // to a cell and stores in the summable digits tree
  //

  const Int_t kPi0 = 111;
  const Int_t kGamma = 22;
  Int_t   npmd;
  Int_t   trackno;
  Int_t   smnumber;
  Int_t   trackpid;
  Int_t   mtrackno;
  Int_t   mtrackpid;
  Int_t   xpad = -1, ypad = -1;
  Float_t edep;
  Float_t vx = -999.0, vy = -999.0, vz = -999.0;
  Float_t xPos, yPos, zPos;
  Float_t xx, yy;

  AliPMDUtility cc;

  for (Int_t ievt = 0; ievt < nevt; ievt++)
    {

      printf("Event Number = %d\n",ievt);
      Int_t nparticles = fRunLoader->GetHeader()->GetNtrack();
      printf("Number of Particles = %d\n",nparticles);
      fRunLoader->GetEvent(ievt);
      // ------------------------------------------------------- //
      // Pointer to specific detector hits.
      // Get pointers to Alice detectors and Hits containers
      
      TTree* treeH = fPMDLoader->TreeH();
      
      Int_t ntracks    = (Int_t) treeH->GetEntries();
      printf("Number of Tracks in the TreeH = %d\n", ntracks);
      
      
      TClonesArray* hits = 0;
      if (fPMD) hits = fPMD->Hits();
      
      // Start loop on tracks in the hits containers
      
      for (Int_t track=0; track<ntracks;track++) 
	{
	  gAlice->ResetHits();
	  treeH->GetEvent(track);
	  if (fPMD) 
	    {
	      npmd = hits->GetEntriesFast();
	      for (int ipmd = 0; ipmd < npmd; ipmd++) 
		{
		  fPMDHit = (AliPMDhit*) hits->UncheckedAt(ipmd);
		  trackno = fPMDHit->GetTrack();

		  //fprintf(fpw,"trackno = %d\n",trackno);

		  //  get kinematics of the particles
		  
		  TParticle* mparticle = gAlice->GetMCApp()->Particle(trackno);
		  trackpid  = mparticle->GetPdgCode();
		  
		  Int_t igatr = -999;
		  Int_t ichtr = -999;
		  Int_t igapid = -999;
		  Int_t imo;
		  Int_t igen = 0;
		  Int_t idmo = -999;
		  
		  Int_t tracknoOld=0, trackpidOld=0, statusOld = 0;
		  if (mparticle->GetFirstMother() == -1)
		    {
		      tracknoOld  = trackno;
		      trackpidOld = trackpid;
		      statusOld   = -1;
		      
		      vx = mparticle->Vx();
		      vy = mparticle->Vy();
		      vz = mparticle->Vz();

		      //fprintf(fpw,"==> Mother ID %5d %5d %5d Vertex: %13.3f %13.3f %13.3f\n", igen, -1, trackpid, vx, vy, vz);
		      
		    }
		  Int_t igstatus = 0;
		  while((imo = mparticle->GetFirstMother()) >= 0)
		    {
		      igen++;
		      
		      mparticle =  gAlice->GetMCApp()->Particle(imo);
		      idmo = mparticle->GetPdgCode();
		      
		      vx = mparticle->Vx();
		      vy = mparticle->Vy();
		      vz = mparticle->Vz();
		      
		      //printf("==> Mother ID %5d %5d %5d Vertex: %13.3f %13.3f %13.3f\n", igen, imo, idmo, vx, vy, vz);
		      //fprintf(fpw,"==> Mother ID %5d %5d %5d Vertex: %13.3f %13.3f %13.3f\n", igen, imo, idmo, vx, vy, vz);
		      
		      if ((idmo == kGamma || idmo == -11 || idmo == 11) && vx == 0. && vy == 0. && vz == 0.)
			{
			  igatr = imo;
			  igapid = idmo;
			  igstatus = 1;
			}
		      if(igstatus == 0)
			{
			  if (idmo == kPi0 && vx == 0. && vy == 0. && vz == 0.)
			    {
			      igatr = imo;
			      igapid = idmo;
			    }
			}
		      ichtr = imo;
		    } // end of while loop
		  
		  if (idmo == kPi0 && vx == 0. && vy == 0. && vz == 0.)
		    {
		      mtrackno = igatr;
		      mtrackpid = igapid;
		    }
		  else
		    {
		      mtrackno  = ichtr;
		      mtrackpid = idmo;
		    }
		  if (statusOld == -1)
		    {
		      mtrackno  = tracknoOld;
		      mtrackpid = trackpidOld;
		    }
		  
		  xPos = fPMDHit->X();
		  yPos = fPMDHit->Y();
		  zPos = fPMDHit->Z();
		  
		  edep       = fPMDHit->GetEnergy();
		  Int_t vol1 = fPMDHit->GetVolume(1); // Column
		  Int_t vol2 = fPMDHit->GetVolume(2); // Row
		  Int_t vol3 = fPMDHit->GetVolume(7); // UnitModule
		  Int_t vol6 = fPMDHit->GetVolume(8); // SuperModule
		  
		  // -----------------------------------------//
		  // For Super Module 1 & 2                   //
		  //  nrow = 96, ncol = 48                    //
		  // For Super Module 3 & 4                   //
		  //  nrow = 48, ncol = 96                    //
		  // -----------------------------------------//
		  
		  smnumber = (vol6-1)*6 + (vol3-1);
		  
		  xpad = vol1 - 1;
		  ypad = vol2 - 1;
		  
		  if(zPos > 361.5)
		    {
		      cc.RectGeomCellPos(smnumber,xpad,ypad,xx,yy);
		      h2->Fill(xx,yy);
		    }
		  
		}
	    }
	} // Track Loop ended
      
    }

  h2->Draw();

  fRunLoader->UnloadgAlice();
  fRunLoader->UnloadHeader();
  fRunLoader->UnloadKinematics();
  fPMDLoader->UnloadHits();

  timer.Stop();
  timer.Print();

}
