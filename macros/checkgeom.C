#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TFile.h"

TStopwatch *timer = 0;

Int_t propagate_in_geom(Double_t *, Double_t *);
void score(TGeoVolume *, Int_t, Double_t);
Double_t timing_per_volume(TGeoVolume *);
Double_t *val1 = 0;
Double_t *val2 = 0;
Bool_t   *flags = 0;

void checkgeom()
{
   TGeoManager::Import("geometry.root");
   Int_t nvol = gGeoManager->GetListOfVolumes()->GetEntries();
   Int_t nuid = gGeoManager->GetListOfUVolumes()->GetEntries();
   timer = new TStopwatch();
   Int_t i;
   Double_t value;
   flags = new Bool_t[nuid];
   val1 = new Double_t[nuid];
   val2 = new Double_t[nuid];
   memset(flags, 0, nuid*sizeof(Bool_t));
   memset(val1, 0, nuid*sizeof(Double_t));
   memset(val2, 0, nuid*sizeof(Double_t));
   TGeoVolume *vol;

// STAGE 1: Overlap checking by sampling per volume

   printf("====================================================================\n");
   printf("STAGE 1: Overlap checking by sampling per volume\n");
   printf("====================================================================\n");
   for (i=0; i<nvol; i++) {
      vol = (TGeoVolume*)gGeoManager->GetListOfVolumes()->At(i);
      Int_t uid = vol->GetNumber();
      if (flags[uid]) continue;
      flags[uid] = kTRUE;
      if (!vol->GetNdaughters()) continue;
      vol->CheckOverlaps(0.01, "s"); 
   }

// STAGE 2: Global overlap checking
   printf("====================================================================\n");
   printf("STAGE 2: Global overlap checking\n");
   printf("====================================================================\n");
   gGeoManager->CheckOverlaps(0.01);
   
// STAGE 3: How many crossings per volume in a realistic case ?
   // Ignore volumes with no daughters

   // Generate rays from vertex in phi=[0,2*pi] theta=[0,pi]
   Int_t ntracks = 1000000;
   printf("====================================================================\n");
   printf("STAGE 3: Propagating %i tracks starting from vertex\n and conting number of boundary crossings...\n", ntracks);
   printf("====================================================================\n");
   Int_t nbound = 0;
   Double_t theta, phi;
   Double_t point[3], dir[3];
   new TRandom3();

   timer->Start();
   for (i=0; i<ntracks; i++) {
      phi = 2.*TMath::Pi()*gRandom->Rndm();
      theta= TMath::ACos(1.-2.*gRandom->Rndm());
      point[0] = point[1] = point[2] = 0;
      dir[0]=TMath::Sin(theta)*TMath::Cos(phi);
      dir[1]=TMath::Sin(theta)*TMath::Sin(phi);
      dir[2]=TMath::Cos(theta);
      if ((i%1000)==0) printf("... remaining tracks %i\n", ntracks-i);
      nbound += propagate_in_geom(point,dir);
   }
   Double_t time1 = timer->CpuTime() *1.E6;
   Double_t time2 = time1/ntracks;
   Double_t time3 = time1/nbound;
   printf("Time for crossing %i boundaries: %g [ms]\n", nbound, time1);
   printf("Time per track for full geometry traversal: %g [ms], per crossing: %g [ms]\n", time2, time3);

// STAGE 4: How much time per volume:
   
   printf("====================================================================\n");
   printf("STAGE 4: How much navigation time per volume per next+safety call\n");
   printf("====================================================================\n");
   TGeoIterator next(gGeoManager->GetTopVolume());
   TGeoNode*current;
   TString path;
   vol = gGeoManager->GetTopVolume();
   memset(flags, 0, nuid*sizeof(Bool_t));
   score(vol, 1, timing_per_volume(vol)); 
   while ((current=next())) {
      vol = current->GetVolume();
      Int_t uid = vol->GetNumber();
      if (flags[uid]) continue;
      flags[uid] = kTRUE;
      next.GetPath(path);
      gGeoManager->cd(path.Data());
      score(vol,1,timing_per_volume(vol));
   }   

   // Draw some histos
   Double_t time_tot_pertrack = 0.;
   TCanvas *c1 = new TCanvas("c2","ncrossings",10,10,900,500);
   c1->SetGrid();
   c1->SetTopMargin(0.15);
   TFile *f = new TFile("histos.root", "RECREATE");
   TH1F *h = new TH1F("h","number of boundary crossings per volume",3,0,3);
   h->SetStats(0);
   h->SetFillColor(38);
   h->SetBit(TH1::kCanRebin);
   
   memset(flags, 0, nuid*sizeof(Bool_t));
   for (i=0; i<nuid; i++) {
      vol = gGeoManager->GetVolume(i);
      if (!vol->GetNdaughters()) continue;
      time_tot_pertrack += val1[i]*val2[i];
      h->Fill(vol->GetName(), (Int_t)val1[i]);
   }
   time_tot_pertrack /= ntracks;
   h->LabelsDeflate();
   h->LabelsOption(">","X");
   h->Draw();   


   TCanvas *c2 = new TCanvas("c3","time spent per volume in navigation",10,10,900,500);
   c2->SetGrid();
   c2->SetTopMargin(0.15);
   TH2F *h2 = new TH2F("h2", "time per FNB call vs. ndaughters", 100, 0,100,100,0,15);
   h2->SetStats(0);
   h2->SetMarkerStyle(2);
   TH1F *h1 = new TH1F("h1","percent of time spent per volume",3,0,3);
   h1->SetStats(0);
   h1->SetFillColor(38);
   h1->SetBit(TH1::kCanRebin);
   for (i=0; i<nuid; i++) {
      vol = gGeoManager->GetVolume(i);
      if (!vol->GetNdaughters()) continue;
      value = val1[i]*val2[i]/ntracks/time_tot_pertrack;
      h1->Fill(vol->GetName(), value);
      h2->Fill(vol->GetNdaughters(), val2[i]);
   }     
   h1->LabelsDeflate();
   h1->LabelsOption(">","X");
   h1->Draw();
   TCanvas *c3 = new TCanvas("c4","timing vs. ndaughters",10,10,900,500);
   c3->SetGrid();
   c3->SetTopMargin(0.15);
   h2->Draw();   
   f->Write();
   delete [] flags;
   delete [] val1;
   delete [] val2;
}

Int_t propagate_in_geom(Double_t *start, Double_t *dir)
{
// Propagate from START along DIR from boundary to boundary until exiting 
// geometry. Fill array of hits.
   gGeoManager->InitTrack(start, dir);
   TGeoNode *current = 0;
   Int_t nzero = 0;
   Int_t nhits = 0;
   while (!gGeoManager->IsOutside()) {
      current = gGeoManager->FindNextBoundaryAndStep(TGeoShape::Big(), kFALSE);
      if (!current || gGeoManager->IsOutside()) return nhits;
      Double_t step = gGeoManager->GetStep();
      if (step<2.*TGeoShape::Tolerance()) {
         nzero++;
         continue;
      } 
      else nzero = 0;
      if (nzero>3) {
      // Problems in trying to cross a boundary
         printf("Error in trying to cross boundary of %s\n", current->GetName());
         return nhits;
      }
      // Generate the hit
      nhits++;
      TGeoVolume *vol = current->GetVolume();
      score(vol,0,1.);
      Int_t iup = 1;
      TGeoNode *mother = gGeoManager->GetMother(iup++);
      while (mother && mother->GetVolume()->IsAssembly()) {
         score(mother->GetVolume(), 0, 1.);
         mother = gGeoManager->GetMother(iup++);
      }   
   }
   return nhits;
}      

void score(TGeoVolume *vol, Int_t ifield, Double_t value)
{
// Score something for VOL
   Int_t uid = vol->GetNumber();
   switch (ifield) {
      case 0:
         val1[uid] += value;
         break;
      case 1:
         val2[uid] += value;
   }
}   

Double_t timing_per_volume(TGeoVolume *vol)
{
// Compute timing per "FindNextBoundary" + "Safety" call. Volume must be
// in the current path.
   timer->Reset();
   const TGeoShape *shape = vol->GetShape();
   TGeoBBox *box = (TGeoBBox *)shape;
   Double_t dx = box->GetDX();
   Double_t dy = box->GetDY();
   Double_t dz = box->GetDZ();
   Double_t ox = (box->GetOrigin())[0];
   Double_t oy = (box->GetOrigin())[1];
   Double_t oz = (box->GetOrigin())[2];
   Double_t point[3], dir[3], lpt[3], ldir[3];
   Double_t pstep = 0.;
   pstep = TMath::Max(pstep,dz);
   Double_t theta, phi;
   Int_t idaughter;
   timer->Start();
   Double_t dist;
   Bool_t inside;
   for (Int_t i=0; i<1000000; i++) {
      lpt[0] = ox-dx+2*dx*gRandom->Rndm();
      lpt[1] = oy-dy+2*dy*gRandom->Rndm();
      lpt[2] = oz-dz+2*dz*gRandom->Rndm();
      gGeoManager->GetCurrentMatrix()->LocalToMaster(lpt,point);
      gGeoManager->SetCurrentPoint(point[0],point[1],point[2]);
      phi = 2*TMath::Pi()*gRandom->Rndm();
      theta= TMath::ACos(1.-2.*gRandom->Rndm());
      ldir[0]=TMath::Sin(theta)*TMath::Cos(phi);
      ldir[1]=TMath::Sin(theta)*TMath::Sin(phi);
      ldir[2]=TMath::Cos(theta);
      gGeoManager->GetCurrentMatrix()->LocalToMasterVect(ldir,dir);
      gGeoManager->SetCurrentDirection(dir);
      gGeoManager->SetStep(pstep);
      gGeoManager->ResetState();
      inside = kTRUE;
      dist = TGeoShape::Big();
      if (!vol->IsAssembly()) {
         inside = vol->Contains(lpt);
         if (!inside) {
            dist = vol->GetShape()->DistFromOutside(lpt,ldir,3,pstep); 
//            if (dist>=pstep) continue;
         } else {   
            vol->GetShape()->DistFromInside(lpt,ldir,3,pstep);
         }   
            
         if (!vol->GetNdaughters()) vol->GetShape()->Safety(lpt, inside);
      }   
      if (vol->GetNdaughters()) {
         gGeoManager->Safety();
         gGeoManager->FindNextDaughterBoundary(point,dir,idaughter,kFALSE);
      }   
   }
   timer->Stop();
   Double_t time_per_track = timer->CpuTime();
   Int_t uid = vol->GetNumber();
   Int_t ncrossings = (Int_t)val1[uid];
   if (!vol->GetNdaughters())
      printf("Time for volume %s (shape=%s): %g [ms] ndaughters=%d ncross=%d\n", vol->GetName(), vol->GetShape()->GetName(), time_per_track, vol->GetNdaughters(), ncrossings);
   else
      printf("Time for volume %s (assemb=%d): %g [ms] ndaughters=%d ncross=%d\n", vol->GetName(), vol->IsAssembly(), time_per_track, vol->GetNdaughters(), ncrossings);
   return time_per_track;
}
   
