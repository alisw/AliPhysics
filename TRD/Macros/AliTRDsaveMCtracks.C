#ifndef __CINT__
  #include <iostream.h>
  #include "AliTRDtracker.h"

  #include "AliTRDcluster.h" 
  #include "AliTRDhit.h" 
  #include "AliTRDv1.h"
  #include "AliTRDgeometry.h"    
  #include "AliTRDparameter.h"    
  #include "alles.h"  
  #include "AliTRDmcTrack.h"

  #include "AliTRDtrack.h"
  
  #include "TFile.h"
  #include "TParticle.h"
  #include "TStopwatch.h"
#endif

Int_t Find(Int_t prtcl, TObjArray *mctarray) {

// Returns index of the mct which corresponds to particle <prtcl> 

  Int_t b=0, e=mctarray->GetEntriesFast(), m=(b+e)/2;
  AliTRDmcTrack *mct = 0; 
  while ((b<e)&&(e!=0)) {
    m=(b+e)/2;
    mct = (AliTRDmcTrack*) mctarray->UncheckedAt(m);
    if (prtcl > mct->GetTrackIndex()) b=m+1;
    else e=m;
  }

  if(mct->GetTrackIndex() == prtcl) return m;

  if((m+1) < mctarray->GetEntriesFast()) {
    mct = (AliTRDmcTrack*) mctarray->UncheckedAt(m+1);
    if(mct->GetTrackIndex() == prtcl) return m+1;
  }

  return -1;
}                    

void AliTRDsaveMCtracks() {

  TObjArray mctracks(2000);
  TObjArray *TRDmcTracks = &mctracks; 

  Char_t *alifile = "galice.root";
  Int_t   nEvent  = 0; 
  Float_t low_pt_cut = 0.05;
  Float_t low_p_cut = 0.02;
  Int_t min_no_of_clusters = 10;

  // Connect the AliRoot file containing Geometry, Kine, Hits, and Digits
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

  AliTRDv1       *fTRD     = (AliTRDv1*) gAlice->GetDetector("TRD");   
  
  // Import the Trees for the event nEvent in the file
  const Int_t nparticles = gAlice->GetEvent(nEvent);

  printf("found %d particles in event %d \n", nparticles, nEvent);
    
  // Create TRDmcTracks for each particle
  Int_t label = -1, charge = 0;
  Bool_t primary;
  Float_t mass;
  Int_t pdg_code;

  printf("\n");

  for (Int_t ii=0; ii<nparticles; ii++) {

    printf("\r particle %d out of %d",ii,nparticles);

    TParticle *p = gAlice->Particle(ii);
    if(p->P() < low_p_cut) continue;

    primary = kTRUE;
    if (p->GetFirstMother()>=0) primary = kFALSE;
    
    pdg_code = (Int_t) p->GetPdgCode();

    if ((pdg_code == 10010020) ||
 	(pdg_code == 10010030) || 
 	(pdg_code == 50000050) || 
 	(pdg_code == 50000051) || 
 	(pdg_code == 10020040)) {

      mass = 0.;
      charge = 0;
    }
    else {
      TParticlePDG *pdg = p->GetPDG();
      charge = (Int_t) pdg->Charge(); 
      mass=pdg->Mass();
    }
    if(charge == 0) continue;

    AliTRDmcTrack *t = new AliTRDmcTrack(ii, primary, mass, charge, pdg_code);
    TRDmcTracks->AddLast(t);
  }
  printf("\r\n");

  // Loop through TRD clusters and assign indexes to MC tracks
  
  TFile *geofile =TFile::Open("AliTRDclusters.root");   
  AliTRDtracker *Tracker = new AliTRDtracker(geofile);
  Tracker->SetEventNumber(nEvent);

  AliTRDgeometry *fGeo   = (AliTRDgeometry*) geofile->Get("TRDgeometry"); 
  AliTRDparameter *fPar   = (AliTRDparameter*) geofile->Get("TRDparameter"); 


  Char_t *clusterfile = "AliTRDclusters.root";
  TObjArray carray(2000);
  TObjArray *ClustersArray = &carray;  
  Tracker->ReadClusters(ClustersArray,clusterfile);

  Int_t nClusters = carray.GetEntriesFast();

  Int_t track, index;
  AliTRDmcTrack *tpr = NULL;

  for (Int_t i = 0; i < nClusters; i++) {

    printf("\r assigning cluster %d out of %d", i, nClusters);
    AliTRDcluster *cl = (AliTRDcluster *) ClustersArray->UncheckedAt(i);

    for(Int_t j=0; j<3; j++) {
      track = cl->GetLabel(j);
      if(track < 0) continue;
      if(track >= nparticles) 
	printf("Track index %d is larger than total number of particles %d\n"
               ,track,nparticles);
      else {
	index = Find(track, TRDmcTracks);
	tpr = 0;
	if(index > 0) tpr = (AliTRDmcTrack *) TRDmcTracks->UncheckedAt(index);
	
	if(tpr) {
	  label = tpr->GetTrackIndex();
	  if(label != track) printf("Got label %d while expected %d !\n",
				    label, track);
	  TParticle *p = gAlice->Particle(track);

	  if (p->Pt() > low_pt_cut) tpr->Update(i);	  
	}	  
      }
    }
  }
  
  // Loop through the TRD hits and set XYZ and Pin and Pout

  TTree *hitTree = gAlice->TreeH();
  gAlice->ResetHits();
  Int_t   preamp, amp, det, plane = 0, nBytes = 0;
  Int_t   nTrack = (Int_t) hitTree->GetEntries();
  Float_t  pos[3], rot[3];
  Float_t  prepos[3], prerot[3];

  Bool_t  entrance = kFALSE;
  Bool_t  exit = kFALSE;


  Double_t xin[6], xout[6];
  for(Int_t j = 0; j < 6; j++) {
    xin[j] = Tracker->GetX(0,j,14);
    xout[j] = Tracker->GetX(0,j,0);
  }

  AliTRDhit *hit;

  Int_t current_track = -1;

  printf("\n");
  while(nTrack--) {

    gAlice->ResetHits();
    nBytes += hitTree->GetEvent(nTrack);
    entrance = kTRUE;
    preamp = 1111;
    amp = 1111;
    for(Int_t j = 0; j < 3; j++) { pos[j] = 1111.; rot[j] = 1111.;}

    for(hit = (AliTRDhit *) fTRD->FirstHit(-1); 
	hit; 
	hit = (AliTRDhit *) fTRD->NextHit()) {

      preamp = amp;
      for(Int_t j = 0; j < 3; j++) { prepos[j] = pos[j]; prerot[j] = rot[j];}

      amp = hit->GetCharge();
      if(amp != 0) continue; 
      track   = hit->GetTrack(); 

      TParticle *p = gAlice->Particle(track);
      if (p->Pt() <= low_pt_cut) continue; 

      if(track != current_track) {
	current_track = track;
	//	printf("\n\n");
      }

      det = hit->GetDetector();
      plane    = fGeo->GetPlane(det); 
      pos[0] = hit->X();
      pos[1] = hit->Y();
      pos[2] = hit->Z();
      fGeo->Rotate(det,pos,rot);

      if(TMath::Abs(rot[0] - xin[plane]) < 1.2) entrance = kTRUE;
      else  entrance = kFALSE;
      if(TMath::Abs(rot[0] - xout[plane]) < 1.2) exit = kTRUE;
      else  exit = kFALSE;
      
      if(entrance && (preamp == 0) && (prerot[0] < 200)) {
	index = Find(track, TRDmcTracks);
	if(index > 0) {
	  tpr = 0;
	  tpr = (AliTRDmcTrack *) TRDmcTracks->UncheckedAt(index);
	  if(tpr) {
	    label = tpr->GetTrackIndex();
	    if(label != track) printf("Got label %d while expected %d !\n",
				      label, track);	    
	    else {
	      /*
	      printf("\n momentum at plane %d entrance, Pxyz: %f, %f, %f\n",
		     plane, prepos[0], prepos[1], prepos[2]);
	      printf(" XYZ at plane %d entrance: %f, %f, %f\n",
		     plane, rot[0], rot[1], rot[2]);
	      printf("expected x range: %f .. %f\n",xin[plane]-1.2,xin[plane]+1.2);
	      */
	      tpr->SetPin(plane, prepos[0], prepos[1], prepos[2]);
	      tpr->SetXYZin(plane, rot[0], rot[1], rot[2]);
	    }
	  }
	}
      }


      if(exit && (preamp == 0) && (prerot[0] < 200)) {
	index = Find(track, TRDmcTracks);
	if(index > 0) {
	  tpr = 0;
	  tpr = (AliTRDmcTrack *) TRDmcTracks->UncheckedAt(index);
	  if(tpr) {
	    label = tpr->GetTrackIndex();
	    if(label != track) printf("Got label %d while expected %d !\n",
				      label, track);	    
	    else {
	      /*
	      printf("\n momentum at plane %d exit, Pxyz: %f, %f, %f\n",
		     plane, prepos[0], prepos[1], prepos[2]);
	      printf(" XYZ at plane %d exit: %f, %f, %f\n",
		     plane, rot[0], rot[1], rot[2]);
	      printf("expected x range: %f .. %f\n",xout[plane]-1.2,xout[plane]+1.2);
	      */
	      tpr->SetPout(plane, prepos[0], prepos[1], prepos[2]);
	      tpr->SetXYZout(plane, rot[0], rot[1], rot[2]);
	    }
	  }
	}
      }
    }
  }
  
	
  // Write TRDmcTracks in output file  

  TDirectory *savedir=gDirectory; 
  Char_t *filename   = "AliTRDmcTracks.root";
  TFile *out = new TFile(filename,"RECREATE");

  TTree *tracktree = new TTree("MCtracks","TRD MC tracks");

  AliTRDmcTrack *iotrack=0;
  tracktree->Branch("MCtracks","AliTRDmcTrack",&iotrack,32000,0);
  
  Int_t ntracks = TRDmcTracks->GetEntriesFast();
  
  for (Int_t i=0; i<ntracks; i++) {
    AliTRDmcTrack *pt=(AliTRDmcTrack*)TRDmcTracks->UncheckedAt(i);
    
    Int_t n = pt->GetNumberOfClusters();
    if(n > min_no_of_clusters) { 
      iotrack=pt;
      tracktree->Fill();
      printf("Put track with label %d and %d clusters in the output tree \n",
	     pt->GetTrackIndex(),n);
    }
  }

  tracktree->Write();
  out->Close();     
  savedir->cd(); 

  return;
}

