#ifndef __CINT__
  #include <iostream.h>

  #include "AliTRDtracker.h"
  #include "AliTRDcluster.h"
  #include "AliTRDhit.h"
  #include "AliTRDv1.h"
  #include "AliTRDgeometry.h"
  #include "AliTRDparameter.h"

  #include "alles.h"
  #include "TFile.h"
  #include "TStopwatch.h"

#endif    

void AliTRDclusterErrors() {

  Int_t No_of_tracks_to_analyze = 40;

  const Int_t tbpp = 15;
  const Int_t nPlanes = 6;
  const Int_t ntb = tbpp * nPlanes; 

  TH1F *hy = new TH1F("delta R*phi","Cluster displacement in R*phi",200,-1.,1.); 
  TH1F *hyp = new TH1F("delta R*phi pos","delta R*phi, positive",200,-1.,1.); 
  TH1F *hym = new TH1F("delta R*phi neg","delta R*phi, negative",200,-1.,1.); 

  TH1F *hyn = new TH1F("Norm., d(R*phi)","Norm. cluster displacement in R*phi",400,-8.,8.); 
  TH1F *hz = new TH1F("delta Z","Cluster displacement in Z",300,-10.,50.); 
  TH2F *hy2 = new TH2F("Amp vs delta R*phi","Amplitude versus delta R*phi",200,-5.,5.,200,0.,600.); 
  TH2F *herr = new TH2F("sigmaY vs delta R*phi","sigmaY vs delta R*phi",200,-1,1,200,0.,0.1); 
  TH2F *hy3 = new TH2F("Position within pad vs delta R*phi","Position within pad vs delta R*phi",200,-1.,1.,200,-0.5,1.5); 
  TH2F *hy4 = new TH2F("local tb vs delta R*phi","local tb vs delta R*phi",200,-1.,1.,20,-2.5,17.5); 

  hy->SetXTitle("Displacement, cm"); 
  hyn->SetXTitle("Displacement, SigmaY"); 
  hy2->SetXTitle("Displacement, cm"); 
  hy2->SetYTitle("Amplitude"); 
  hz->SetXTitle("Displacement, cm"); 
  hy3->SetXTitle("Displacement, cm"); 
  hy3->SetYTitle("Position, cm"); 
  hy4->SetXTitle("Displacement, cm"); 
  hy4->SetYTitle("local time bin"); 

  /*
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
    cout << "Loaded shared libraries" << endl;
  } 
  */      

  // Load clusters
  Char_t *clusterfile = "AliTRDclusters.root";
  Int_t   nEvent  = 0;

  TObjArray carray(2000);
  TObjArray *ClustersArray = &carray;
  TFile *geofile =TFile::Open("AliTRDclusters.root");   

  AliTRDparameter *par = (AliTRDparameter*) geofile->Get("TRDparameter");
  AliTRDgeometry *fGeo = (AliTRDgeometry*) geofile->Get("TRDgeometry");   

  AliTRDtracker *Tracker = new AliTRDtracker(geofile);
  Tracker->SetEventNumber(nEvent);
  Tracker->ReadClusters(ClustersArray,clusterfile);
  Int_t nClusters = carray.GetEntriesFast();

  printf("Total number of clusters %d \n", nClusters);

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

  AliTRDv1       *fTRD           = (AliTRDv1*) gAlice->GetDetector("TRD");

  // Import the Trees for the event nEvent in the file
  Int_t nparticles = gAlice->GetEvent(nEvent);

  TObjArray *particles=gAlice->Particles();

  TTree *hitTree = gAlice->TreeH();
  Int_t  nBytes = 0;

  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrack = (Int_t) hitTree->GetEntries();
  cout << " Found " << nTrack << " primary particles with hits" << endl;   
  No_of_tracks_to_analyze = TMath::Min(No_of_tracks_to_analyze, nTrack);


  // Loop through particles and fill histoes

  Float_t hitY[ntb];
  Float_t hitZ[ntb];
  Float_t hitO[ntb];
  Float_t clusterY[ntb];
  Float_t clusterZ[ntb];
  Float_t clusterQ[ntb];
  Float_t clusterSigmaY[ntb];
  Float_t pos[3];
  Float_t rot[3];
  Float_t global[3];
  Int_t det = 0;
  Int_t track_index[3];

  printf("\n");

  for (Int_t ii=0; ii<No_of_tracks_to_analyze; ii++) {
    printf("track %d out of %d \n", ii+1 , No_of_tracks_to_analyze); 

    for(Int_t plane = 0; plane < nPlanes; plane++) {
      for(Int_t ltb = 14; ltb > -1; ltb--) {

	if(ii >= nTrack) continue;

	TParticle *p = gAlice->Particle(ii);
	if (p->GetFirstMother()>=0) continue;
	TParticlePDG *pdg = p->GetPDG();
	Float_t charge=pdg->Charge(); 
	if(TMath::Abs(charge) < 0.5) continue;

	Int_t gtb = Tracker->GetGlobalTimeBin(0,plane,ltb);
	Double_t x = Tracker->GetX(0,plane,ltb);

	// loop through clusters

	Bool_t cluster_found = kFALSE;
	Int_t nw = 0;

	for (Int_t i = 0; i < nClusters; i++) {

	  AliTRDcluster *cl = (AliTRDcluster *) carray.UncheckedAt(i);

	  nw = cl->GetDetector(); 

	  nw = fGeo->GetPlane(nw);

	  if(nw != plane) continue; 
	  for(Int_t j=0; j<3; j++) track_index[j] = cl->GetLabel(j); 
	  if((track_index[0] != ii) &&
	     (track_index[1] != ii) &&
	     (track_index[2] != ii)) continue;

	  nw = cl->GetLocalTimeBin(); 
	  if(nw != ltb) continue;

	  clusterY[gtb] = cl->GetY();
	  clusterZ[gtb] = cl->GetZ();
	  clusterQ[gtb] = cl->GetQ();

	  clusterSigmaY[gtb] = TMath::Sqrt(cl->GetSigmaY2());
	  cluster_found = kTRUE;
	  break;	  
	}

	if(!cluster_found) continue;

	gAlice->ResetHits();

	//	nBytes += hitTree->GetEvent(nPrimaries - ii - 1);
	nBytes += hitTree->GetEvent(nTrack - ii - 1);

	// Loop through the TRD hits
	Bool_t found_hit = kFALSE;


	for(AliTRDhit *hit = (AliTRDhit *) fTRD->FirstHit(-1); 
	    hit; 
	    hit = (AliTRDhit *) fTRD->NextHit()) {
	  nw = hit->Track();
	  if(nw != ii) continue;
	  det   = hit->GetDetector();
	  nw = fGeo->GetPlane(det);
	  if(nw != plane) continue;           


	  pos[0]=hit->X(); 
	  pos[1]=hit->Y();
	  pos[2]=hit->Z();
	  fGeo->Rotate(det,pos,rot);

	  if(TMath::Abs(rot[0]-x) > 0.01) continue;
	  hitY[gtb] = rot[1];
	  hitZ[gtb] = rot[2];

	  Float_t col0        = par->GetCol0(plane);
	  Float_t colPadSize  = par->GetColPadSize(plane);         
	  Float_t colH = (Int_t ((rot[1] -  col0)/colPadSize)) * colPadSize;  
	  hitO[gtb] = (rot[1] -  col0) - colH;
	  found_hit = kTRUE;
	  break;
	}

	if(!found_hit) continue;

	/*	
	printf("gtb: %d, x: %f, rot[0]: %f, Yhit: %f, Ycl: %f\n",
	       gtb, x, rot[0], rot[1], clusterY[gtb]);
	printf("\n                            Zhit - Zcl = %f - %f = %f\n",
	       rot[2], clusterZ[gtb], rot[2] - clusterZ[gtb]);


		
	printf("found hit within dx = %f - %f \n",rot[0],x);
	printf("pos: %f, %f, %f \n",pos[0],pos[1],pos[2]);
	printf("rot: %f, %f, %f \n",rot[0],rot[1],rot[2]);
	printf("cluster: %d, %f, %f \n",gtb,clusterY[gtb],clusterZ[gtb]);
	*/


	hy->Fill(hitY[gtb]-clusterY[gtb]);
	if(charge > 0) hyp->Fill(hitY[gtb]-clusterY[gtb]);
	else if(charge < 0) hym->Fill(hitY[gtb]-clusterY[gtb]);

	if((clusterQ[gtb]>10)&&(clusterSigmaY[gtb]>0)) 
	  hyn->Fill((hitY[gtb]-clusterY[gtb])/clusterSigmaY[gtb]);
	hz->Fill(hitZ[gtb]-clusterZ[gtb]);
	hy2->Fill(hitY[gtb]-clusterY[gtb],clusterQ[gtb]);
	hy3->Fill(hitY[gtb]-clusterY[gtb],hitO[gtb]);  
	hy4->Fill(hitY[gtb]-clusterY[gtb],(Float_t)(tbpp - 1 - gtb%tbpp));
	herr->Fill(hitY[gtb]-clusterY[gtb],clusterSigmaY[gtb]);
      }
    }
  }
  
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1); 

  TCanvas* c = new TCanvas("c", "c", 110, 110, 810, 840);
  c->SetFillColor(10);
  c->Divide(2,2);
  c->cd(1); hy->SetLineWidth(2); hy->SetFillColor(29); hy->Draw();
  c->cd(2); hz->SetLineWidth(2); hz->SetFillColor(29); hz->Draw();
  c->cd(3); hyn->Draw();
  c->cd(4); hy4->Draw();

  TCanvas* c1 = new TCanvas("c1", "c1", 210, 210, 910, 940);
  c1->SetFillColor(10);
  c1->Divide(2,2);
  c1->cd(1); hyp->SetLineWidth(2); hyp->SetFillColor(29); hyp->Fit("gaus");
  c1->cd(2); hym->SetLineWidth(2); hym->SetFillColor(29); hym->Fit("gaus");
  c1->cd(3); hy3->Draw();
  c1->cd(4); hy2->Draw();
   
}






