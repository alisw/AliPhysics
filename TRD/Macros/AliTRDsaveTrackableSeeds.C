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
  #include "AliTrackReference.h"
  
  #include "TFile.h"
  #include "TParticle.h"
  #include "TStopwatch.h"
#endif


void AliTRDsaveTrackableSeeds() {

  TObjArray mctracks(2000);
  TObjArray *TRDmcTracks = &mctracks; 

  TFile *geofile =TFile::Open("AliTRDclusters.root");
  AliTRDtracker *Tracker = new AliTRDtracker(geofile);
  Int_t nEvent = 0;
  Tracker->SetEventNumber(nEvent);

  AliTRDgeometry *fGeo   = (AliTRDgeometry*) geofile->Get("TRDgeometry");
  //AliTRDparameter *fPar = (AliTRDparameter*) geofile->Get("TRDparameter");  

  Char_t *alifile = "galice.root";

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
    
  // Create TRDmcTracks for each tpc seed
  Int_t label = -1, charge = 0;
  Bool_t primary;
  Float_t mass;
  Int_t pdg_code;

  TDirectory *savedir=gDirectory;

  TFile *in=TFile::Open("AliTPCBackTracks.root");
  if (!in->IsOpen()) {
    cerr<<"can't open file AliTPCBackTracks.root  !\n"; return;
  }

  char   tname[100];
  sprintf(tname,"seedsTPCtoTRD_%d",nEvent);
  TTree *seedTree=(TTree*)in->Get(tname);
  if (!seedTree) {
     cerr<<"AliTRDtracker::PropagateBack(): ";
     cerr<<"can't get a tree with seeds from TPC !\n";
  }

  AliTPCtrack *seed=new AliTPCtrack;
  seedTree->SetBranchAddress("tracks",&seed);

  Int_t nSeeds = (Int_t) seedTree->GetEntries();
  for (Int_t is=0; is<nSeeds; is++) {
     seedTree->GetEvent(is);
     Int_t lbl = seed->GetLabel();
     if(TMath::Abs(lbl) > nparticles) { 
       printf("extra high seed label %d \n", lbl);
       continue;
     }
     TParticle *p = gAlice->Particle(TMath::Abs(lbl));

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

     AliTRDmcTrack *t = new AliTRDmcTrack(TMath::Abs(lbl), lbl, primary, 
					  mass, charge, pdg_code);
     TRDmcTracks->AddLast(t);
  }
  delete seed;
  delete seedTree;

  savedir->cd();                         

  Int_t i, mctIndex[nparticles];
  for(i=0; i<nparticles; i++) mctIndex[i]=-1;

  Int_t nMCtracks = TRDmcTracks->GetEntriesFast();
  AliTRDmcTrack *mct = 0;
  for(i = 0; i < nMCtracks; i++) {
    mct = (AliTRDmcTrack*) TRDmcTracks->UncheckedAt(i);
    mctIndex[mct->GetTrackIndex()] = i;
  }

  // Loop through the TRD rec points and set XYZ and Pin and Pout

  Double_t xin[6], xout[6];
  Double_t Px, Py, Pz;

  Float_t  pos[3], rot[3];

  for(Int_t j = 0; j < 6; j++) {
    xin[j] = Tracker->GetX(0,j,0)-3;
    xout[j] = Tracker->GetX(0,j,0);
  }

  TTree *trRefTree = gAlice->TreeTR();
  if (!trRefTree) {
    cout << "<AliTRDreadTrackRef> No TR tree found" << endl;
    return;
  }       

  Int_t nBytes   = 0;
  Bool_t fEntrance, fExit;

  Int_t nTrack = (Int_t) trRefTree->GetEntries();
  cout << "<AliTRDreadTrackRef> Found " << nTrack
       << " primary particles with track refs" << endl;

  for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) {

    printf("TrackRefs for track %d out of %d \n",iTrack,nTrack);
    gAlice->ResetTrackReferences();
    nBytes += trRefTree->GetEvent(iTrack);  

    AliTrackReference *tr = 0;
    tr = (AliTrackReference*) fTRD->FirstTrackReference(-1);
    
    while (tr) {

      Int_t   track = tr->GetTrack();
      if(mctIndex[track] >= 0) { 
	mct = (AliTRDmcTrack*) TRDmcTracks->UncheckedAt(mctIndex[track]);
	
	pos[0] = tr->X();
	pos[1] = tr->Y();
	pos[2] = tr->Z();
	
	Double_t phi = TMath::ATan2(pos[1],pos[0]);
	if(phi < 0) phi = 2 * TMath::Pi() + phi;
	
	Int_t sec = ((Int_t) (phi*180/TMath::Pi())) / 20;
	fGeo->Rotate(fGeo->GetDetector(0,0,17-sec), pos, rot);
	
	fEntrance = kFALSE; fExit = kFALSE;
	for(i = 0; i < 6; i++) {
	  if(TMath::Abs(xin[i] - rot[0]) < 1.4) { fEntrance = kTRUE; break; }
	  if(TMath::Abs(xout[i] - rot[0]) < 1.4) { fExit = kTRUE; break; }
	}
      
	Px = tr->Px();
	Py = tr->Py();
	Pz = tr->Pz();
	
	if(fEntrance) {
	  mct->SetXYZin(i, rot[0], rot[1], rot[2]);
	  mct->SetPin(i, Px, Py, Pz);
	  /*
	    printf("entr plane %d: rp_x = %f vs %f = xin \n",i,rot[0],xin[i]); 
	    printf("             : Px, Py, Pz = %f, %f, %f \n",Px,Py,Pz); 
	    printf("             : y, z = %f, %f\n",rot[1],rot[2]); 
	  */
	}
	if(fExit) {
	  mct->SetXYZout(i, rot[0], rot[1], rot[2]);
	  mct->SetPout(i, Px, Py, Pz);
	  /*
	    printf("exit plane %d: rp_x = %f vs %f = xout \n",i,rot[0],xout[i]); 
	    printf("             : Px, Py, Pz = %f, %f, %f \n",Px,Py,Pz); 
	    printf("             : y, z = %f, %f\n",rot[1],rot[2]);
	  */ 
	}	
      }                                   
      tr = (AliTrackReference *) fTRD->NextTrackReference();
    }
  }  

  // Loop through TRD clusters and assign indexes to MC tracks
  
  Char_t *clusterfile = "AliTRDclusters.root";
  TObjArray carray(2000);
  TObjArray *ClustersArray = &carray;  
  Tracker->ReadClusters(ClustersArray,clusterfile);

  printf("done with ReadClusters()\n");

  Int_t nClusters = carray.GetEntriesFast();

  Int_t track, det, det0, det1, ltb, plane, ind0, ind1, ncl;
  Double_t y, y0, y1, q, q0, q1;

  AliTRDcluster *c0 = NULL, *c1 = NULL;
  printf("nClusters = %d \n", nClusters);

  for (Int_t i = 0; i < nClusters; i++) {

    printf("\r assigning cluster %d out of %d", i, nClusters);
    AliTRDcluster *cl = (AliTRDcluster *) ClustersArray->UncheckedAt(i);

    for(Int_t j=0; j<3; j++) {
      track = cl->GetLabel(j);
      if(track < 0) continue;
      if(track >= nparticles) { 
	printf("Track index %d is larger than total number of particles %d\n"
               ,track,nparticles);
	continue;
      }     
      if(mctIndex[track] < 0) continue;
      mct = (AliTRDmcTrack*) TRDmcTracks->UncheckedAt(mctIndex[track]);

      label = mct->GetTrackIndex();
      if(label != track) { 
	printf("Got label %d while expected %d !\n", label, track);
	continue;
      }

      ncl = mct->GetNumberOfClusters();
      mct->SetNumberOfClusters(ncl+1);

      det=cl->GetDetector();
      plane = fGeo->GetPlane(det);
      ltb=cl->GetLocalTimeBin();
      if((ltb < 0) || (ltb >= kMAX_TB)) continue;

      ind0 = mct->GetClusterIndex(ltb, plane, 0);
      ind1 = mct->GetClusterIndex(ltb, plane, 1);
      if(ind0 < 0) {
	mct->Update(ltb,plane,0,i);	       
      } else {
	c0 = (AliTRDcluster *) ClustersArray->UncheckedAt(ind0);
	det0 = c0->GetDetector();
	y =  cl->GetY();
	y0 = c0->GetY();
	q  = TMath::Abs(cl->GetQ());
	q0 = TMath::Abs(c0->GetQ());
	if((det == det0) && (TMath::Abs(y-y0) < 1.5)) {
	  if(q > q0) mct->Update(ltb,plane,0,i);
	} else {
	  if (ind1 < 0) {
	    mct->Update(ltb,plane,1,i);
	  } 
	  else {
	    c1 = (AliTRDcluster *) ClustersArray->UncheckedAt(ind1);
	    det1 = c1->GetDetector();
	    y1 = c1->GetY();
	    q1 = TMath::Abs(c1->GetQ());
	    if((det == det0) && (TMath::Abs(y-y0) < 1.5)) {
	      if(q > q1) mct->Update(ltb,plane,1,i);
	    }	  
	  }	  
	}
      }
    }
  } 	    
  	
  // Write TRDmcTracks in output file  

  savedir=gDirectory; 
  Char_t *filename   = "AliTRDtrackableSeeds.root";
  TFile *out = new TFile(filename,"RECREATE");

  TTree *tracktree = new TTree("MCtracks","TRD MC tracks");

  AliTRDmcTrack *iotrack=0;
  tracktree->Branch("MCtracks","AliTRDmcTrack",&iotrack,32000,0);
  
  Int_t ntracks = TRDmcTracks->GetEntriesFast();
  
  for (Int_t i=0; i<ntracks; i++) {
    AliTRDmcTrack *pt=(AliTRDmcTrack*)TRDmcTracks->UncheckedAt(i);    
    iotrack=pt;
    tracktree->Fill();
    printf("Put track with label %d and %d clusters in the output tree \n",
	   pt->GetTrackIndex(),pt->GetNumberOfClusters());
  }

  tracktree->Write();
  out->Close();     
  savedir->cd(); 

  return;
}

