//
// Reading macro for the track reference tree
// minimum stuff
//
void AliPMDTrackRefRead()
{

  AliRunLoader*    runLoader =  AliRunLoader::Open("galice.root");

  if (!runLoader)
    { 
      printf("Could not open galice.root");
    }

  runLoader->LoadTrackRefs("READ");

  Int_t nEvents = runLoader->GetNumberOfEvents();

  printf("Total number of Events = %d\n",nEvents);


  FILE *fpw = fopen("trackref.dat","w");



  TClonesArray *arrayTR = 0x0;
  arrayTR = new TClonesArray("AliTrackReference");

  TTree *treeTR;


  AliStack *stack = runLoader->Stack();


  // Starts the event loop

  for (Int_t iev = 0; iev < nEvents; iev++)
    {
      runLoader->GetEvent(iev);
      treeTR = runLoader->TreeTR();
      treeTR->SetBranchAddress("TrackReferences",  &arrayTR);
      
      Int_t nEntries = treeTR->GetEntries();
      
      printf("Total number of Entries = %d \n",nEntries);
      
      for (Int_t i = 0; i < nEntries; i++)
	{
	  Int_t trRead  = treeTR->GetEntry(i);
	  if (trRead <= 0) continue;
	  Int_t nTrackRefs = arrayTR->GetEntries();
	  
	  for (Int_t j = 0; j < nTrackRefs; j++)
	    {
	      AliTrackReference* trackRef = static_cast<AliTrackReference*>(arrayTR->At(j));
	      if (!trackRef) continue;
	      
	      
	      if (trackRef->DetectorId() != AliTrackReference::kPMD) continue;
	      
	      //TParticle* mpart = 0;
	      
	      //Int_t trackno = trackRef->GetTrack();
	      Int_t trackno = trackRef->Label();
	      
	      //printf("Track number = %d \n",trackno);
	      if (iev == 0) fprintf(fpw,"Track number = %d \n",trackno);
	      
	      
	      //mpart = stack->Particle(trackno);

	      //Int_t   MPID   = mpart->GetFirstMother();
	      //Int_t   M2PID  = mpart->GetSecondMother();

	    }
	  
	}
    }
  
  
  
  
  
}
