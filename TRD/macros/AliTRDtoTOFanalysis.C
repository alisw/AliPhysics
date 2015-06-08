#ifndef __CINT__
  #include <iostream.h>  

  #include "AliTOF.h"
  #include "AliTOFhit.h"
  #include "AliTPCtrack.h"
  #include "AliTRDtrack.h"
  #include "AliTRDgeometry.h"

  #include "alles.h"  
  #include "TFile.h"
  #include "TParticle.h"
  #include "TStopwatch.h"
#endif     

void AliTRDtoTOFanalysis() 
{

  Int_t nEvent = 0;
  const Int_t maxIndex = 80000;  // max number of primaries to be analysed
  Int_t rtIndex[maxIndex];
  for(Int_t i = 0; i < maxIndex; i++) rtIndex[i] = -1;

  //****** tracking:  Load reconstructed tracks

  TFile *tf=TFile::Open("AliTRDtracks.root");

  if (!tf->IsOpen()) {cerr<<"Can't open AliTRDtracks.root !\n"; return;}
  TObjArray tarray(2000);

  char   tname[100];
  sprintf(tname,"seedsTRDtoTOF2_%d",nEvent);
  //  for other radial positions use different trees:     
  //    sprintf(tname,"seedsTRDtoTOF1_%d",nEvent);     
  //    sprintf(tname,"seedsTRDtoPHOS_%d",nEvent);     
  //    sprintf(tname,"seedsTRDtoRHIC_%d",nEvent);     

  TTree *tracktree=(TTree*)tf->Get(tname);

  TBranch *tbranch=tracktree->GetBranch("tracks");

  Int_t nRecTracks = (Int_t) tracktree->GetEntries();
  cout<<"Found "<<nRecTracks<<" entries in the track tree "<<tname<<endl;

  for (Int_t i=0; i<nRecTracks; i++) {
    AliTPCtrack *iotrack=new AliTPCtrack();
    tbranch->SetAddress(&iotrack);
    tracktree->GetEvent(i);
    tarray.AddLast(iotrack);
    Int_t trackLabel = iotrack->GetLabel();

    //    printf("rt with %d clusters and label %d \n",
    //	   iotrack->GetNumberOfClusters(), trackLabel);

    if(trackLabel < 0) continue;
    if(trackLabel >= maxIndex) continue;
    rtIndex[trackLabel] = i;
  }
  tf->Close();       

  //********

  TH1F *hdx = new TH1F("dx","X(hit) - X(track)",40,-10,10);
  TH1F *hdy = new TH1F("dy","Y(hit) - Y(track)",100,-20,20);
  TH1F *hdz = new TH1F("dz","Z(hit) - Z(track)",100,-20,20);
  TH2F *hdt = new TH2F("dt","T(hit) vs S(track)",100,0,500,500,200,700);


  // Connect the AliRoot file containing Geometry, Kine, Hits, and Digits
  Char_t *alifile = "galice.root";

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


  // Import the Trees for the event nEvent in the file
  const Int_t nparticles = gAlice->GetEvent(nEvent);
  if (nparticles <= 0) return;

   Float_t x,y,z,mass,e;
   Double_t gX, gY, gZ;   // rec. track coordinates in the global system 
   Double_t Px, Py, Pz;   // rec. track momentum components in the global system 
   Float_t X, Y, Z;     // local tracking coordinates for reconstructed track
   Int_t nbytes = 0;
   Int_t j,hit,ipart;
   Int_t nhits;
   Float_t tof;
   TParticle *particle;            
   AliTPCtrack *rt;

// Get pointers to Alice detectors and Hits containers
   AliDetector *TOF  = gAlice->GetDetector("TOF");
   Int_t ntracks    = (Int_t) gAlice->TreeH()->GetEntries();     

// Start loop on tracks in the hits containers

   for (Int_t track=0; track < ntracks; track++) {

     if(TOF) {
       for(AliTOFhit* tofHit = (AliTOFhit*)TOF->FirstHit(track); 
	   tofHit; 
	   tofHit=(AliTOFhit*)TOF->NextHit()) {

         ipart    = tofHit->GetTrack();
	 if(ipart >= maxIndex) continue;
	 if(rtIndex[ipart] < 0) continue; 

	 x = tofHit->X();
	 y = tofHit->Y();
	 z = tofHit->Z();

	 //******* tracking:  extract track coordinates, momentum, etc. 

	 rt = (AliTPCtrack*) tarray.UncheckedAt(rtIndex[ipart]);
	 Double_t tr_length = rt->GetLength();  // meaningfull only with ITS
	                                        // part fully implemented 

	 AliTRDtrack *trd_rt = new AliTRDtrack(*rt, rt->GetAlpha());
	 trd_rt->GetGlobalXYZ(gX,gY,gZ);  // returns running global coordinates 
	 trd_rt->GetPxPyPz(Px,Py,Pz);  // returns momentum in global system
  
	 /*                                             
	 printf("\n hit position - rec. track position:\n");
	 printf(" X: %f - %f = %f \n", x, gX, x - gX);
	 printf(" Y: %f - %f = %f \n", y, gY, y - gY);
	 printf(" Z: %f - %f = %f \n", z, gZ, z - gZ);
	 printf("\n");
	 */

	 delete trd_rt;

	 // conversion from hit coordinates from global coordinate system 
	 // to "local" tracking system 

	 Float_t phi =TMath::ATan2(y,x);
	 if (phi > 2.*TMath::Pi()) phi -= 2.*TMath::Pi();
	 if (phi < 0.            ) phi += 2.*TMath::Pi();
	 Int_t sec = Int_t(phi/AliTRDgeometry::GetAlpha()) % 18;
	 Float_t alpha = (sec + 0.5) * AliTRDgeometry::GetAlpha();

	 Float_t tmp=x*TMath::Cos(alpha) + y*TMath::Sin(alpha);
	 y= - x*TMath::Sin(alpha) + y*TMath::Cos(alpha);
	 x=tmp;  

	 //         particle = (TParticle*)gAlice->Particle(ipart);
	 //        if (particle->GetFirstMother() < 0) continue;

	 X = (Float_t) rt->GetX();   // radial position in the tracking system
	 Y = (Float_t) rt->GetY();   // r*phi position within the sector
	 Z = (Float_t) rt->GetZ();   // Z position


	 if(TMath::Abs(X-x-1) > 2) continue;

	 hdx->Fill(X-x);
	 hdy->Fill(Y-y);
	 hdz->Fill(Z-z);

	 //**********

       }
     }
   }                     
     
  TCanvas* c1 = new TCanvas("c1", "c1", 210, 210, 910, 940);
  c1->SetFillColor(10);
  c1->Divide(2,2);
  c1->cd(1); hdx->Draw();
  c1->cd(2); hdy->Draw();
  c1->cd(3); hdz->Draw();
  c1->cd(4); hdt->Draw();
}
