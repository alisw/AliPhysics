PHOS (Int_t event_number=0,Float_t SignalStep=0.001,Int_t SignalMin=15)
{
// From Int_t Sread ()

   // Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gSystem->Load("libGeant3Dummy.so");   // a dummy version of Geant3
      gSystem->Load("PHOS/libPHOSreconstruction.so");        // the standard Alice classes 
      gSystem->Load("libgalice.so");        // the standard Alice classes 
   }
    
   // Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!file) file = new TFile("galice.root");

   // Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
////////////////////////////////

  AliPHOS *PHOS= gAlice->GetDetector("PHOS");
  if( NULL==PHOS )
  {
    cout << "Can't find PHOS detector!\n";
    exit(1);
  }

  PHOS->SetTreeAddress();

  if( 0==PHOS->fTreePHOS->GetEvent(event_number) )
  {
    printf("Cannot read event number %d\n",event_number);
    return;
  }

  printf("This is event number %d from %d\n",
         event_number,PHOS->fTreePHOS->GetEntries());

  for( int i=0; i<PHOS->fCradles->GetEntries(); i++ )
  {
    AliPHOSCradle &cradle = PHOS->GetCradle(i);
    printf("===============================================================\n");
    printf("Cradle %d\n",i+1);
    cradle.Print();
    //cout.flush();
    

    cradle.Reconstruction(SignalStep,SignalMin);

    printf("It were %d particles in that cradle\n",cradle.GetParticles().GetEntries());
    for( int j=0; j<cradle.GetParticles().GetEntries(); j++ )
    {
      TObjArray &pp = cradle.GetParticles();
      AliPHOSgamma *g = (AliPHOSgamma *) pp.At(j);
      printf("%3d  ",j+1);
      g->Print();
    }

    printf("\nReconstruction: %d gammas\n",cradle.GetGammasReconstructed().GetEntries());
    for( int j=0; j<cradle.GetGammasReconstructed().GetEntries(); j++ )
    {
      TObjArray &pp = cradle.GetGammasReconstructed();
      AliPHOSgamma *g = (AliPHOSgamma *) pp.At(j);
      printf("%3d  ",j+1);
      g->Print();
    }

    //cout.flush();
    
  }

////////////////////////////////////////////////////////////////////////////////

  printf("Done\n");
}
