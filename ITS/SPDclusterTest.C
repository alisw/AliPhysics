#include "iostream.h"

void SPDclusterTest (Int_t evNumber1=0,Int_t evNumber2=0) 
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
 
//
//   Loop over events 
//
   Int_t Nh=0;
   Int_t Nh1=0;
   for (int nev=0; nev<= evNumber2; nev++) {
     Int_t nparticles = gAlice->GetEvent(nev);
     cout << "nev         " << nev <<endl;
     cout << "nparticles  " << nparticles <<endl;
     if (nev < evNumber1) continue;
     if (nparticles <= 0) return;

     TTree *TH = gAlice->TreeH();
     Int_t ntracks = TH->GetEntries();
     cout<<"ntracks "<<ntracks<<endl;

   Int_t nbytes = 0;

   AliITSRawClusterSPD  *ITSclust;

// Get pointers to Alice detectors and Digits containers
   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   TClonesArray *Particles = gAlice->Particles();

   if (ITS) {
     // fill modules with sorted by module hits
     Int_t nmodules;
     ITS->InitModules(-1,nmodules); 
     ITS->FillModules(nev,-1,evNumber2,nmodules," "," ");
     //get pointer to modules array
     TObjArray *ITSmodules = ITS->GetModules();
     AliITShit *itsHit;

     // get the Tree for clusters
     ITS->GetTreeC(nev);
     TTree *TC=ITS->TreeC();
     Int_t nent=TC->GetEntries();
     printf("Found %d entries in the tree (must be one per module per event!)\n",nent);
   
     for (Int_t idettype=0;idettype<3;idettype++) {

       TClonesArray *ITSclusters  = ITS->ClustersAddress(idettype);
       //printf ("ITSclusters %p \n",ITSclusters);

          if (idettype != 0) continue;


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


	  // -------------- Create ntuples --------------------

	  //  ntuple structures:

	  struct {
	    Int_t lay;
	    Int_t nx;
	    Int_t nz;
	    Int_t hitprim;
	    Int_t partcode;
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
	    Int_t noverlaps;
	    Int_t noverprim;
	    Float_t qcl;
	    Float_t dx;
	    Float_t dz;
	  } ntuple1_st;


	  struct {
	    //	    Int_t lay;
	    Int_t nx;
	    Int_t nz;
	  } ntuple2_st;


	  ntuple = new TTree("ntuple","Demo ntuple");
	  ntuple->Branch("lay",&ntuple_st.lay,"lay/I");
	  ntuple->Branch("nx",&ntuple_st.nx,"nx/I");
	  ntuple->Branch("nz",&ntuple_st.nz,"nz/I");
	  ntuple->Branch("hitprim",&ntuple_st.hitprim,"hitprim/I");
	  ntuple->Branch("partcode",&ntuple_st.partcode,"partcode/I");
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
	  ntuple1->Branch("noverlaps",&ntuple1_st.noverlaps,"noverlaps/I");
	  ntuple1->Branch("noverprim",&ntuple1_st.noverprim,"noverprim/I");
	  ntuple1->Branch("dx",&ntuple1_st.dx,"dx/F");
	  ntuple1->Branch("dz",&ntuple1_st.dz,"dz/F");


	  ntuple2 = new TTree("ntuple2","Demo ntuple2");
	  //	  ntuple2->Branch("lay",&ntuple2_st.lay,"lay/I");
	  ntuple2->Branch("nx",&ntuple2_st.nx,"nx/I");
	  ntuple2->Branch("nz",&ntuple2_st.nz,"nz/I");

// ------------------------------------------------------------------------

	  // Module loop

	  for (Int_t mod=0; mod<nent; mod++) {
	      AliITSmodule *itsModule = (AliITSmodule*)ITSmodules->At(mod);

	      Int_t nhits = itsModule->GetNhits();
              if(nhits) printf("module nhits %d %d\n",mod,nhits);
	      if(!nhits) continue;
     
              ITS->ResetClusters();
              TC->GetEvent(mod);
	      Int_t nclust = ITSclusters->GetEntries();
	      if (nclust) printf("Found %d clust for module %d in det type %d \n",nclust,mod,idettype);
	      if (!nclust) continue;

	      // cluster/hit loops

	for (Int_t clu=0;clu<nclust;clu++) {
		itsclu   = (AliITSRawClusterSPD*)ITSclusters->UncheckedAt(clu);
		printf("%d %d %f %f %f\n",itsclu->NclZ(),itsclu->NclX(),itsclu->Q(),itsclu->X(),itsclu->Z());

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
		Float_t clusterx = itsclu->X();
		Float_t clusterz = itsclu->Z();
		Float_t clusterQ = itsclu->Q();

		ntuple2_st.nx = clustersizex;
		ntuple2_st.nz = clustersizez;

		ntuple2->Fill();

		Int_t icl = 0;
		Float_t dxprimlast = 10.e+6;
		Float_t dzprimlast = 10.e+6;


		//        if(module > 217 && module <  226) {
                cout<<"mod,nclust,clu,Nxpix,Nzpix ="<<mod<<","<<nclust<<","<<clu<<","<<clustersizex<<","<<clustersizez<<endl;
                cout<<"clusx,clusz ="<<clusterx<<","<<clusterz<<endl;
                cout<<"XStartf,XStopf,ZStart,ZStop ="<<fxstart<<","<<fxstop<<","<<zstart<<","<<zstop<<endl;
		//         }
		

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

 		  Int_t clusterlayer = hitlayer;
		  Int_t clusterladder= hitladder;
		  Int_t clusterdetector = hitdet;

		  Int_t track = itsHit->fTrack;
		  Int_t dray = 0;
		  Int_t hitstat = itsHit->GetTrackStatus();


		  Float_t zhit = 10000*itsHit->GetZL();
		  Float_t xhit = 10000*itsHit->GetXL();

		  if(abs(zhit) > SPDlength/2) {
		    if(hitstat == 66) zhit0 = 1e+5;
		    continue;
		  }

		  if(abs(xhit) > SPDwidth/2) {
		    if(hitstat == 66) xhit0 = 1e+5;
		    continue;
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

      cout<<"clu,hit,xmed,fxstart,fxstop,zmed,zstart,zstop ="<<clu<<","<<hit<<","<<xmed<<","<<fxstart<<","<<fxstop<<","<<zmed<<","<<zstart<<","<<zstop<<endl;

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


        Float_t px = itsHit->GetPXL();
        Float_t py = itsHit->GetPYL();
        Float_t pz = itsHit->GetPZL();
        Float_t pmod = 1000*sqrt(px*px+py*py+pz*pz);

              cout<<"track,partcode,pmod,parent ="<<track<<","<<partcode<<","<<pmod<<","<<parent<<endl;

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
         ntuple_st.dx = xdif;
         ntuple_st.dz = zdif;
         ntuple_st.pmod = pmod;

         ntuple->Fill();

      

      if(hitprim > 0) {   // for primary particles
        if(hitlayer == 1) {
     cout<<"!!!!!! lay,hitprim,xdif,zdif ="<<hitlayer<<","<<hitprim<<","<<xdif<<","<<zdif<<endl;
           Xres1->Fill(xdif);
           Zres1->Fill(zdif);
        }
        if(hitlayer == 2) {
     cout<<"!!!!!! lay,hitprim,xdif,zdif ="<<hitlayer<<","<<hitprim<<","<<xdif<<","<<zdif<<endl;
           Xres2->Fill(xdif);
           Zres2->Fill(zdif);
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

         ntuple1_st.lay = clusterlayer;
         ntuple1_st.lad = clusterladder;
         ntuple1_st.det = clusterdetector;
         ntuple1_st.nx = clustersizex;
         ntuple1_st.nz = clustersizez;
         ntuple1_st.qcl = clusterQ;
         ntuple1_st.noverlaps = noverlaps;
         ntuple1_st.noverprim = noverprim;
         ntuple1_st.dx = dxprimlast;
         ntuple1_st.dz = dzprimlast;

         ntuple1->Fill();

     } // icl = 1
    } // cluster loop
   } // module loop       
  } // idettype loop
 } // end if ITS
} // event loop 


   //  Write and Draw Histogramms and ntuples



   TFile fhistos("SPD_his.root","RECREATE");

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

//     cout<<"END  test for clusters and hits "<<endl;

//     file->Close();   
}



