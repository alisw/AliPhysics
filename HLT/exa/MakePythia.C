//$Id$

/* A little macro filling a root 
   tree with PYTHIA6 events. Same 
   as generate_PYTHIA.C */

typedef struct 
{
  Int_t    no;
  Int_t    part;
  Int_t    total;
  Int_t    npart;
  Int_t    nenergy;
  Int_t    nhard;
  Float_t  b; 
  Float_t  phi;
} THeader;

UInt_t makeSeed(int mode=0) 
{
  switch (mode) {
  case 1: 
    {
      TDatime date;
      return date.GetDate();
    }
  case 2:
    {
      TDatime date;
      UInt_t  seed1 = (date.GetDate() / 1000000 * 100 + 
		       date.GetTime() / 1000000);
      UInt_t  seed2 = (date.GetTime() % 10000);
      TRanMar ranmar(seed1, seed2);
      Int_t   eat = (date.GetDate() - 19980101 + seed2 +
		     gSystem->GetPid()); 
      for (Int_t i = 0; i < (eat + 250000) ; i++) 
	ranmar.Rndm();
      return 2 * Int_t(ranmar.Rndm() * (TMath::Power(2,30) - 1)) + 1;
    }
  case 3:
    {
      TDatime date;
      UInt_t  seed1 = (date.GetDate() / 1000000 * 100 + 
		       date.GetTime() / 1000000 + gSystem->GetPid());
      TRandom rand(seed1);
      return Int_t(rand.Rndm() * (TMath::Power(2,30) - 1)) + 1;
    }

    break;
  }
  return 0;
}

Int_t MakePythia(Int_t nEvents=10,Char_t *rfile="pythia6.root")
{
   gSystem->Load("$ROOTSYS/lib/libPythia6.so");
   gSystem->Load("$ROOTSYS/lib/libEG.so");                   // Root Event Generator interface
   gSystem->Load("$ROOTSYS/lib/libEGPythia6.so");            // Root interface to Pythia6
   //gDebug=1;

   TStopwatch*   stopWatch = new TStopwatch();
   TPythia6* pythia6       = new TPythia6();
   TClonesArray* particles = new TClonesArray("TParticle");
   THeader header;

   Char_t filename[1024];
   sprintf(filename,"%s",rfile);

   TFile* file = new TFile(filename,"RECREATE");
   file->SetCompressionLevel(1);
   TTree* tree = new TTree("pythia","pythia events");
   tree->Branch("particles",&particles);
   //tree->Branch("header",&header, "no/I:part:total:npart:nenergy:nhard:b/F:phi");

   stopWatch->Start();

   //      select Pythia min. bias model, taken from AliPythia.C
   pythia6->SetMSEL(0);
   pythia6->SetMSUB(92,1);      // single diffraction AB-->XB
   pythia6->SetMSUB(93,1);      // single diffraction AB-->AX
   pythia6->SetMSUB(94,1);      // double diffraction
   pythia6->SetMSUB(95,1);      // low pt production
   pythia6->SetMSTP(81,1);      // multiple interactions switched on
   pythia6->SetMSTP(82,3);      // model with varying impact param. & a single Gaussian
   pythia6->SetPARP(82,3.47);   // set value pT_0  for turn-off of the cross section of
   // multiple interaction at a reference energy = 14000 GeV
   pythia6->SetPARP(89,14000.); // reference energy for the above parameter
   pythia6->SetPARP(90,0.174);  // set exponent for energy dependence of pT_0

   // set fragmentation on (default)
   pythia6->SetMSTP(111,1);
  
   // don't smear the primary vertex
   pythia6->SetMSTP(151,0);

   //pythia6->SetMSTP(61,0);
   //pythia6->SetMSTP(64,0);

   pythia6->SetMRPY(1,makeSeed(3));
   pythia6->Initialize("cms","p","p",5500);
   //pythia6->Initialize("cms","p","p",1800);

   stopWatch->Stop();
   stopWatch->Print();

   Int_t i;
   for (i = 0; i < nEvents; i++) {
     stopWatch->Start();
     cout << "Event # " << i << " ... " << flush;

     pythia6->GenerateEvent();

     header.no    = i;
     header.total = pythia6->ImportParticles(particles,"All");
     header.part  = 0;

     header.npart    = 0;
     header.nenergy  = 0;
     header.nhard    = 0;

     header.b     = 0;
     header.phi   = 0;

     //particles->Print();
     //pythia6->Pylist(1);

     tree->Fill();
     
     cout << "done (" << header.total << " particles) " << flush;
     stopWatch->Stop();
     stopWatch->Print();
   }

   file->Write();
   file->Close();

  return 0;
}

