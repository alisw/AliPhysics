#include "iostream.h"

void SPDclusterTestDubna (Int_t evNumber1=0,Int_t evNumber2=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

// Dynamically link some shared libs

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   } else {
      delete gAlice;
      gAlice=0;
   }

// Connect the Root Galice file containing Geometry, Kine and Hits

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!file) file = new TFile("galice.root");
   file->ls();

// Get AliRun object from file or create it if not on file

   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }

	  // ------------ Cluster and point analysis histogramms ------------

	  TH1F *Nxpix1 = new TH1F("Nxpix1","Cluster size in x(r*phi) direction for layer 1",20,0.,20.);
	  TH1F *Nxpix2 = new TH1F("Nxpix2","Cluster size in x(r*phi) direction for layer 2",20,0.,20.);
	  TH1F *Nzpix1 = new TH1F("Nzpix1","Cluster size in z direction for layer 1",15,0.,15.);
	  TH1F *Nzpix2 = new TH1F("Nzpix2","Cluster size in z direction for layer 2",15,0.,15.);
	  TH1F *Xpix1 = new TH1F("Xpix1","Local x coordinate (mm) for layer 1",20,-2.,18.);
	  TH1F *Xpix2 = new TH1F("Xpix2","Local x coordinate (mm) for layer 2",20,-2.,18.);
	  TH1F *Zpix1 = new TH1F("Zpix1","Local z coordinate (mm) for layer 1",90,-2.,88.);
	  TH1F *Zpix2 = new TH1F("Zpix2","Lolac z coordinate (mm) for layer 2",90,-2.,88.);

	  TH1F *Xres1 = new TH1F("Xres1","Xrec and Xgen difference (micr) for layers 1",100,-200.,200.);
	  TH1F *Xres2 = new TH1F("Xres2","Xrec and Xgen difference (micr) for layers 2",100,-200.,200.);
	  TH1F *Zres1 = new TH1F("Zres1","Zrec and Zgen difference (micr) for layers 1",100,-800.,800.);
	  TH1F *Zres2 = new TH1F("Zres2","Zrec and Zgen difference (micr) for layers 2",100,-800.,800.);

	  TH1F *Ptot1 = new TH1F("Ptot1","Total momentum (GeV/C) for layers 1",100,0.,5.);
	  TH1F *Pz1 = new TH1F("Pz1","Pz (GeV/C) for layers 1",100,-5.,5.);
	  TH1F *Theta1 = new TH1F("Theta1","Theta angle (rad) for layers 1",100,0.,4.);
	  TH1F *Y1 = new TH1F("Y1","Rapidity for layers 1",100,-4.,4.);
	  TH1F *Eta1 = new TH1F("Eta1","PseudoRapidity for layers 1",100,-4.,4.);
	  TH1F *Y1Den = new TH1F("Y1Den","Rapidity for layers 1",100,-0.5,0.5);
	  TH1F *Eta1Den = new TH1F("Eta1Den","PseudoRapidity for layers 1",100,-0.5,0.5);
	  TH1F *Y1DenA = new TH1F("Y1DenA","Rapidity for layers 1",100,-0.5,0.5);
	  TH1F *Eta1DenA = new TH1F("Eta1DenA","PseudoRapidity for layers 1",100,-0.5,0.5);
	  TH1F *Phi1 = new TH1F("Phi1","Phi angle (rad) for layers 1",100,0.,7.);
	  TH1F *Ptot2 = new TH1F("Ptot2","Total momentum (GeV/C) for layers 2",100,0.,5.);
	  TH1F *Pz2 = new TH1F("Pz2","Pz (GeV/C) for layers 2",100,-5.,5.);
	  TH1F *Theta2 = new TH1F("Theta2","Theta angle (rad) for layers 2",100,0.,4.);
	  TH1F *Y2 = new TH1F("Y2","Rapidity for layers 2",100,-4.,4.);
	  TH1F *Eta2 = new TH1F("Eta2","PseudoRapidity for layers 2",100,-4.,4.);
	  TH1F *Y2Den = new TH1F("Y2Den","Rapidity for layers 2",100,-0.5,0.5);
	  TH1F *Eta2Den = new TH1F("Eta2Den","PseudoRapidity for layers 2",100,-0.5,0.5);
	  TH1F *Y2DenA = new TH1F("Y2DenA","Rapidity for layers 2",100,-0.5,0.5);
	  TH1F *Eta2DenA = new TH1F("Eta2DenA","PseudoRapidity for layers 2",100,-0.5,0.5);
	  TH1F *Phi2 = new TH1F("Phi2","Phi angle (rad) for layers 2",100,0.,7.);

	  // -------------- Create ntuples --------------------
	  //  ntuple structures:

	  struct {
	    Int_t lay;
	    Int_t nx;
	    Int_t nz;
	    Int_t hitprim;
	    Int_t partcode;
	    Int_t ntrover;
	    Float_t dx;
	    Float_t dz;
	    Float_t pmod;
	  } ntuple_st;

	  struct {
	    Int_t lay;
	    Int_t lad;
	    Int_t det;
	    Int_t nx;
	    Int_t nz;
	    Int_t ntrover;
	    Int_t noverlaps;
	    Int_t noverprim;
	    Float_t qcl;
	    Float_t dx;
	    Float_t dz;
	    Float_t x;
	    Float_t z;
	  } ntuple1_st;

	  struct {
	    //	    Int_t lay;
	    Int_t lay;
	    Int_t nx;
	    Int_t nz;
	    Float_t x;
	    Float_t z;
	    Float_t qcl;
	  } ntuple2_st;

	  ntuple = new TTree("ntuple","Demo ntuple");
	  ntuple->Branch("lay",&ntuple_st.lay,"lay/I");
	  ntuple->Branch("nx",&ntuple_st.nx,"nx/I");
	  ntuple->Branch("nz",&ntuple_st.nz,"nz/I");
	  ntuple->Branch("hitprim",&ntuple_st.hitprim,"hitprim/I");
	  ntuple->Branch("partcode",&ntuple_st.partcode,"partcode/I");
	  ntuple->Branch("ntrover",&ntuple_st.ntrover,"ntrover/I");
	  ntuple->Branch("dx",&ntuple_st.dx,"dx/F");
	  ntuple->Branch("dz",&ntuple_st.dz,"dz/F");
	  ntuple->Branch("pmod",&ntuple_st.pmod,"pmod/F");

	  ntuple1 = new TTree("ntuple1","Demo ntuple1");
	  ntuple1->Branch("lay",&ntuple1_st.lay,"lay/I");
	  ntuple1->Branch("lad",&ntuple1_st.lad,"lad/I");
	  ntuple1->Branch("det",&ntuple1_st.det,"det/I");
	  ntuple1->Branch("nx",&ntuple1_st.nx,"nx/I");
	  ntuple1->Branch("nz",&ntuple1_st.nz,"nz/I");
	  ntuple1->Branch("qcl",&ntuple1_st.qcl,"qcl/F");
	  ntuple1->Branch("ntrover",&ntuple1_st.ntrover,"ntrover/I");
	  ntuple1->Branch("noverlaps",&ntuple1_st.noverlaps,"noverlaps/I");
	  ntuple1->Branch("noverprim",&ntuple1_st.noverprim,"noverprim/I");
	  ntuple1->Branch("x",&ntuple1_st.x,"x/F");
	  ntuple1->Branch("z",&ntuple1_st.z,"z/F");
	  ntuple1->Branch("dx",&ntuple1_st.dx,"dx/F");
	  ntuple1->Branch("dz",&ntuple1_st.dz,"dz/F");


	  ntuple2 = new TTree("ntuple2","Demo ntuple2");
	  //	  ntuple2->Branch("lay",&ntuple2_st.lay,"lay/I");
	  ntuple2->Branch("lay",&ntuple2_st.lay,"lay/I");
	  ntuple2->Branch("x",&ntuple2_st.x,"x/F");
	  ntuple2->Branch("z",&ntuple2_st.z,"z/F");
	  ntuple2->Branch("nx",&ntuple2_st.nx,"nx/I");
	  ntuple2->Branch("nz",&ntuple2_st.nz,"nz/I");
	  ntuple2->Branch("qcl",&ntuple2_st.qcl,"qcl/F");

// ------------------------------------------------------------------------
//
//   Loop over events 
//
   for (int nev=0; nev<= evNumber2; nev++) {
     Int_t nparticles = gAlice->GetEvent(nev);
     cout << "nev         " << nev <<endl;
     cout << "nparticles  " << nparticles <<endl;
     if (nev < evNumber1) continue;
     if (nparticles <= 0) return;

     TTree *TH = gAlice->TreeH();
     Int_t ntracks = TH->GetEntries();
     cout<<"ntracks "<<ntracks<<endl;

// Get pointers to Alice detectors and Digits containers
   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   TClonesArray *Particles = gAlice->Particles();

   if (ITS) {
     // fill modules with sorted by module hits
     Int_t nmodules;
     ITS->InitModules(-1,nmodules); 
     //    ITS->FillModules(nev,-1,evNumber2,nmodules," "," ");
     ITS->FillModules(nev,evNumber2,nmodules," "," ");
     //get pointer to modules array
     TObjArray *ITSmodules = ITS->GetModules();
     AliITShit *itsHit;

     // get the Tree for clusters
     ITS->GetTreeC(nev);
     TTree *TC=ITS->TreeC();
     Int_t nent=TC->GetEntries();
     printf("Found %d entries in the tree (must be one per module per event!)\n",nent);
     Int_t lay, lad, det;
     AliITSgeom *geom = ITS->GetITSgeom();
   
     for (Int_t idettype=0;idettype<3;idettype++) {

       TClonesArray *ITSclusters  = ITS->ClustersAddress(idettype);
       //printf ("ITSclusters %p \n",ITSclusters);

          if (idettype != 0) continue;

	  Float_t occup1 = 0;
	  Float_t occup2 = 0;

	  // Module loop
	  for (Int_t mod=0; mod<nent; mod++) {
	      AliITSmodule *itsModule = (AliITSmodule*)ITSmodules->At(mod);
	      geom->GetModuleId(mod,lay,lad,det);

	      Int_t nhits = itsModule->GetNhits();
              //if(nhits) printf("module nhits %d %d\n",mod,nhits);
	      if(!nhits) continue;
     
              ITS->ResetClusters();
              TC->GetEvent(mod);
	      Int_t nclust = ITSclusters->GetEntries();
	      if (!nclust) continue;

	      // cluster/hit loops
	      //cout<<"mod,lay,nclust,nhits ="<<mod<<","<<lay<<","<<nclust<<","<<nhits<<endl;
	for (Int_t clu=0;clu<nclust;clu++) {
		itsclu   = (AliITSRawClusterSPD*)ITSclusters->UncheckedAt(clu);

		Int_t noverlaps = 0;
		Int_t noverprim = 0;

		Int_t clustersizex = itsclu->NclX();
		Int_t clustersizez = itsclu->NclZ();
		//       Int_t xstart = itsclu->XStart();
		//       Int_t xstop = itsclu->XStop();
		Int_t xstart = itsclu->XStartf();
		Int_t xstop = itsclu->XStopf();
		Float_t fxstart = xstart*50;
		Float_t fxstop = (xstop+1)*50;
		Float_t zstart = itsclu->ZStart();
		Float_t zstop = itsclu->ZStop();
		Int_t zend = itsclu->Zend();
		Int_t ntrover = itsclu->NTracks();
		Float_t clusterx = itsclu->X();
		Float_t clusterz = itsclu->Z();
		Float_t clusterQ = itsclu->Q();

		if(lay == 1) occup1 += clusterQ;                
		if(lay == 2) occup2 += clusterQ;                

		ntuple2_st.lay = lay;
		ntuple2_st.x = clusterx/1000.;
		ntuple2_st.z = clusterz/1000.;
		ntuple2_st.nx = clustersizex;
		ntuple2_st.nz = clustersizez;
		ntuple2_st.qcl = clusterQ;

		ntuple2->Fill();

		Int_t icl = 0;
		Float_t dxprimlast = 10.e+6;
		Float_t dzprimlast = 10.e+6;

        	Float_t SPDlength = 83600;	
        	Float_t SPDwidth = 12800;	
                Float_t xhit0 = 1e+5;
                Float_t zhit0 = 1e+5;
		
       for (Int_t hit=0;hit<nhits;hit++) {

		  // Find coordinate differences between the hit and cluster positions
		  // for the resolution determination.

		  itsHit   = (AliITShit*)itsModule->GetHit(hit);

		  Int_t hitlayer = itsHit->GetLayer();
		  Int_t hitladder= itsHit->GetLadder();
		  Int_t hitdet= itsHit->GetDetector();
	          Float_t dEn = 1.0e+6*itsHit->GetIonization(); // hit energy, KeV 
		  Int_t track = itsHit->GetTrack();
		  Int_t dray = 0;
		  Int_t hitstat = itsHit->GetTrackStatus();

		  Float_t zhit = 10000*itsHit->GetZL();
		  Float_t xhit = 10000*itsHit->GetXL();

        Float_t pxsimL = itsHit->GetPXL();  // the momenta at GEANT points
        Float_t pysimL = itsHit->GetPYL();
        Float_t pzsimL = itsHit->GetPZL();
	Float_t psimL = TMath::Sqrt(pxsimL*pxsimL+pysimL*pysimL+pzsimL*pzsimL);

	// Check boundaries
	if(zhit  > SPDlength/2) {
	  //cout<<"!!! z outside ="<<zhit<<endl;
         zhit = SPDlength/2 - 10;
	}
	if(zhit < 0 && zhit < -SPDlength/2) {
	  //cout<<"!!! z outside ="<<zhit<<endl;
         zhit = -SPDlength/2 + 10;
	}
	if(xhit  > SPDwidth/2) {
	  //cout<<"!!! x outside ="<<xhit<<endl;
         xhit = SPDwidth/2 - 10;
	}
	if(xhit  < 0 && xhit < -SPDwidth/2) {
	  //cout<<"!!! x outside ="<<xhit<<endl;
         xhit = -SPDwidth/2 + 10;
	}

		  zhit += SPDlength/2;
		  xhit += SPDwidth/2;
		  Float_t yhit = 10000*itsHit->GetYL();

		  if(hitlayer == 1 && hitstat == 66 && yhit > 71) {
		    xhit0 = xhit;
		    zhit0 = zhit;
		  }
		  if(hitlayer == 2 && hitstat == 66 && yhit < -71) {
		    xhit0 = xhit;
		    zhit0 = zhit;
		  }

		  if(hitstat != 68) continue; // Take only the hit if the last
		  // track point went out from
		  // the detector.

		  if(xhit0 > 9e+4 || zhit0 > 9e+4) continue;

		  Float_t xmed = (xhit + xhit0)/2;
		  Float_t zmed = (zhit + zhit0)/2;

		  Float_t xdif = xmed - clusterx;
		  Float_t zdif = zmed - clusterz;


        // Consider the hits inside of cluster region only

  if((xmed >= fxstart && xmed <= fxstop) && (zmed >= zstart && zmed <= zstop)) {
        icl = 1;

                //        part = (TParticle *)particles.UncheckedAt(track);
                //        Int_t partcode = part->GetPdgCode();
                //              Int_t primery = gAlice->GetPrimary(track);

        Int_t parent = itsHit->GetParticle()->GetFirstMother();
        Int_t partcode = itsHit->GetParticle()->GetPdgCode();

//  partcode (pdgCode): 11 - e-, 13 - mu-, 22 - gamma, 111 - pi0, 211 - pi+
//                      310 - K0s, 321 - K+, 2112 - n, 2212 - p, 3122 - lambda

	Float_t pmod = itsHit->GetParticle()->P(); // total momentum at the
	                                           // vertex
	Float_t energy = itsHit->GetParticle()->Energy(); // energy at the
	                                           // vertex
	Float_t mass = itsHit->GetParticle()->GetMass(); // particle mass 

	Float_t pz = itsHit->GetParticle()->Pz(); // z momentum componetnt  
	                                           // at the vertex
	Float_t px = itsHit->GetParticle()->Px(); // z momentum componetnt  
	                                           // at the vertex
	Float_t py = itsHit->GetParticle()->Py(); // z momentum componetnt  
	                                           // at the vertex
	Float_t phi = itsHit->GetParticle()->Phi(); // Phi angle at the
	                                           // vertex
	Float_t theta = itsHit->GetParticle()->Theta(); // Theta angle at the
	                                           // vertex
	//Float_t eta = itsHit->GetParticle()->Eta(); // Pseudo rapidity at the
	                                           // vertex
	if((energy-pz) > 0) {
	  Float_t y = 0.5*TMath::Log((energy+pz)/(energy-pz));
	}else{
	  cout<<" Warning: energy < pz ="<<energy<<","<<pz<<endl;
	  y = 10;
	}   
        Float_t eta = -TMath::Log(TMath::Tan(theta/2));
        pmod *= 1.0e+3;


        Float_t pxsim = itsHit->GetPXG();  // the momenta at this GEANT point
        Float_t pysim = itsHit->GetPYG();
        Float_t pzsim = itsHit->GetPZG();
	Float_t psim = TMath::Sqrt(pxsim*pxsim+pysim*pysim+pzsim*pzsim);

        Int_t hitprim = 0;

        if(partcode == 11 && pmod < 6) dray = 1; // delta ray is e-
                                                 // at p < 6 MeV/c

        if(dray == 0) noverlaps = noverlaps + 1; // overlapps for all hits but
                                                 // not for delta ray which
                                                 // also went out from the
                                                 // detector and returned
                                                 // again

        if(parent < 0) hitprim = hitprim + 1; // hitprim=1 for the primery
                                              // particles

        if(hitprim > 0) noverprim = noverprim + 1;

        if(hitprim > 0) {
         dxprimlast = xdif;
         dzprimlast = zdif;
        }

        // fill ntuple

         ntuple_st.lay = hitlayer;
         ntuple_st.nx = clustersizex;
         ntuple_st.nz = clustersizez;
         ntuple_st.hitprim = hitprim;
         ntuple_st.partcode = partcode;
         ntuple_st.ntrover = ntrover;
         ntuple_st.dx = xdif;
         ntuple_st.dz = zdif;
         ntuple_st.pmod = pmod;

         ntuple->Fill();

        if(hitlayer == 1) {
           Y1DenA->Fill(y);
           Eta1DenA->Fill(eta);
        }
        if(hitlayer == 2) {
           Y2DenA->Fill(y);
           Eta2DenA->Fill(eta);
        }

      

      if(hitprim > 0) {   // for primary particles

        if(hitlayer == 1) {
           Xres1->Fill(xdif);
           Zres1->Fill(zdif);
           Ptot1->Fill(pmod/1000.);
           Pz1->Fill(pz);
           Theta1->Fill(theta);
           Y1->Fill(y);
           Eta1->Fill(eta);
           Y1Den->Fill(y);
           Eta1Den->Fill(eta);
           Phi1->Fill(phi);
        }
        if(hitlayer == 2) {
           Xres2->Fill(xdif);
           Zres2->Fill(zdif);
           Ptot2->Fill(pmod/1000.);
           Pz2->Fill(pz);
           Theta2->Fill(theta);
           Y2->Fill(y);
           Eta2->Fill(eta);
           Y2Den->Fill(y);
           Eta2Den->Fill(eta);
           Phi2->Fill(phi);
        }
      } // primery particles

     } // end of cluster region
   } // end of hit loop

      if(icl == 1) {

        // fill ntuple1

        //      ntuple1->Fill(clusterlayer,clustersizex,clustersizez,noverlaps,\
noverprim,dx,dz);

        if(noverlaps == 0) noverlaps = 1; // cluster contains one or more
                                          // delta rays only

         ntuple1_st.lay = lay;
         ntuple1_st.lad = lad;
         ntuple1_st.det = det;
         ntuple1_st.x = clusterx*1000.;
         ntuple1_st.z = clusterz*1000.;
         ntuple1_st.nx = clustersizex;
         ntuple1_st.nz = clustersizez;
         ntuple1_st.qcl = clusterQ;
         ntuple1_st.ntrover = ntrover;
         ntuple1_st.noverlaps = noverlaps;
         ntuple1_st.noverprim = noverprim;
         ntuple1_st.dx = dxprimlast;
         ntuple1_st.dz = dzprimlast;

         ntuple1->Fill();

     } // icl = 1
    } // cluster loop
   } // module loop       

     cout<<" Occupancy for layer-1 ="<<occup1<<endl;
     cout<<" Occupancy for layer-2 ="<<occup2<<endl;
     // The real occupancy values are:
     // (for full ALICE event at the full SPD acceptence)
     //   occup1 /= 3932160;    
     //   occup2 /= 7864320;

  } // idettype loop
 } // end if ITS
} // event loop 


   //  Write and Draw Histogramms and ntuples



   TFile fhistos("SPD_his_dubna.root","RECREATE");

   ntuple->Write();
   ntuple1->Write();
   ntuple2->Write();

   Nxpix1->Write();
   Nzpix1->Write();
   Nxpix2->Write();
   Nzpix2->Write();

   Xpix1->Write();
   Zpix1->Write();
   Xpix2->Write();
   Zpix2->Write();

   Xres1->Write();
   Zres1->Write();
   Xres2->Write();
   Zres2->Write();

   Ptot1->Write();
   Pz1->Write();
   Theta1->Write();
   Y1->Write();
   Eta1->Write();
   Y1Den->Write();
   Eta1Den->Write();
   Y1DenA->Write();
   Eta1DenA->Write();
   Phi1->Write();

   Ptot2->Write();
   Pz2->Write();
   Theta2->Write();
   Y2->Write();
   Eta2->Write();
   Y2Den->Write();
   Eta2Den->Write();
   Y2DenA->Write();
   Eta2DenA->Write();
   Phi2->Write();

   fhistos.Close();
   cout<<"!!! Histogramms and ntuples were written"<<endl;

   TCanvas *c1 = new TCanvas("c1","ITS clusters",400,10,600,700);
   c1->Divide(2,2);
   c1->cd(1);
   gPad->SetFillColor(33);
         Xres1->SetFillColor(42);
         Xres1->Draw();
   c1->cd(2);
   gPad->SetFillColor(33);
         Zres1->SetFillColor(46);
         Zres1->Draw();
   c1->cd(3);
   gPad->SetFillColor(33);
         Xres2->SetFillColor(42);
         Xres2->Draw();
   c1->cd(4);
   gPad->SetFillColor(33);
         Zres2->SetFillColor(46);
         Zres2->Draw();

     cout<<"END  test for clusters and hits "<<endl;

//     file->Close();   
}



