void SSDrecpointTest (Int_t evNumber1=0,Int_t evNumber2=0)
  //void SSDrecpointTest (Int_t evNumber1=0,Int_t evNumber2=999)
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and fill some histograms.
//   
//     Root > .L anal.C   //this loads the macro in memory
//     Root > anal();     //by default process first event   
//     Root > anal(2);    //process third event
//Begin_Html
/*
<img src="gif/anal.gif">
*/
//End_Html
/////////////////////////////////////////////////////////////////////////
    
// Dynamically link some shared libs

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }

// Connect the Root Galice file containing Geometry, Kine and Hits
   TString *str = new TString("galice.root");
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(str->Data());
   if (!file) file = new TFile(str->Data(),"UPDATE");

// Get AliRun object from file or create it if not on file
   //   if (!gAlice) {
     gAlice = (AliRun*)file->Get("gAlice");
     if (gAlice) printf("AliRun object found on file\n");
     if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
     //}


     // -------------- Create ntuples --------------------

     //  ntuple structures:


          struct {
            Int_t lay;
            Int_t nxP;
            Int_t nxN;
            Int_t hitprim;
            Int_t partcode;
            Float_t x;
            Float_t z;
            Float_t dx;
            Float_t dz;
            Float_t pmod;
          } ntuple_st;

          struct {
            Int_t lay;
            Int_t lad;
            Int_t det;
            Int_t nxP;
            Int_t nxN;
            Int_t noverlaps;
            Int_t noverprim;
            Float_t qclP;
            Float_t qclN;
            Float_t qrec;
            Float_t dx;
            Float_t dz;
          } ntuple1_st;

          struct {
            Int_t nxP;
            Int_t nxN;
            Float_t x;
            Float_t z;
          } ntuple2_st;

          ntuple = new TTree("ntuple","Demo ntuple");
          ntuple->Branch("lay",&ntuple_st.lay,"lay/I");
          ntuple->Branch("nxP",&ntuple_st.nxP,"nxP/I");
          ntuple->Branch("nxN",&ntuple_st.nxN,"nxN/I");
          ntuple->Branch("hitprim",&ntuple_st.hitprim,"hitprim/I");
          ntuple->Branch("partcode",&ntuple_st.partcode,"partcode/I");
          ntuple->Branch("x",&ntuple_st.x,"x/F");
          ntuple->Branch("z",&ntuple_st.z,"z/F");
          ntuple->Branch("dx",&ntuple_st.dx,"dx/F");
          ntuple->Branch("dz",&ntuple_st.dz,"dz/F");
          ntuple->Branch("pmod",&ntuple_st.pmod,"pmod/F");

	  ntuple1 = new TTree("ntuple1","Demo ntuple1");
	  ntuple1->Branch("lay",&ntuple1_st.lay,"lay/I");
	  ntuple1->Branch("lad",&ntuple1_st.lad,"lad/I");
	  ntuple1->Branch("det",&ntuple1_st.det,"det/I");
	  ntuple1->Branch("nxP",&ntuple1_st.nxP,"nxP/I");
	  ntuple1->Branch("nxN",&ntuple1_st.nxN,"nxN/I");
	  ntuple1->Branch("qclP",&ntuple1_st.qclP,"qclP/F");
          ntuple1->Branch("qclN",&ntuple1_st.qclN,"qclN/F");
          ntuple1->Branch("qrec",&ntuple1_st.qrec,"qrec/F");
          ntuple1->Branch("dx",&ntuple1_st.dx,"dx/F");
          ntuple1->Branch("dz",&ntuple1_st.dz,"dz/F");
	  ntuple1->Branch("noverlaps",&ntuple1_st.noverlaps,"noverlaps/I");
          ntuple1->Branch("noverprim",&ntuple1_st.noverprim,"noverprim/I");

          ntuple2 = new TTree("ntuple2","Demo ntuple2");
          ntuple2->Branch("nxP",&ntuple2_st.nxP,"nxP/I");
          ntuple2->Branch("nxN",&ntuple2_st.nxN,"nxN/I");
          ntuple2->Branch("x",&ntuple2_st.x,"x/F");
          ntuple2->Branch("z",&ntuple2_st.z,"z/F");


	  // Create Histogramms

	  TH1F *NxP5 = new TH1F("NxP5","P cluster size for layer 5",20,0.,20.);
	  TH1F *NxN5 = new TH1F("NxN5","N cluster size for layer 5",20,0.,20.);
	  TH1F *NxP6 = new TH1F("NxP6","P cluster size for layer 6",20,0.,20.);
	  TH1F *NxN6 = new TH1F("NxN6","N cluster size for layer 6",20,0.,20.);

	  TH1F *Xres5 = new TH1F("Xres5","Xrec and Xgen difference (micr) for layers 5",100,-200.,200.);
	  TH1F *Xres6 = new TH1F("Xres6","Xrec and Xgen difference (micr) for layers 6",100,-200.,200.);
	  TH1F *Zres5 = new TH1F("Zres5","Zrec and Zgen difference (micr) for layers 5",100,-8000.,8000.);
	  TH1F *Zres6 = new TH1F("Zres6","Zrec and Zgen difference (micr) for layers 6",100,-8000.,8000.);
          TH1F *Path5 = new TH1F("Path5","Path length in Si",100,0.,600.);
          TH1F *Path6 = new TH1F("Path6","Path length in Si",100,0.,600.);
          TH1F *dEdX = new TH1F("dEdX","dEdX  (KeV)",100,0.,500.);
          TH2F *adcPadcN5all = new TH2F("adcPadcN5all","adcP/N correlation for lay5",100,0.,200.,100,0.,200.);
          TH2F *adcPadcN6all = new TH2F("adcPadcN6all","adcP/N correlation for lay6",100,0.,200.,100,0.,200.);
          TH2F *adcPadcN5cut = new TH2F("adcPadcN5cut","adcP/N correlation for lay5 and cut of P-N signas",100,0.,200.,100,0.,200.);
          TH2F *adcPadcN6cut = new TH2F("adcPadcN6cut","adcP/N correlation for lay6 and cut of P-N signals",100,0.,200.,100,0.,200.);


   AliITS *ITS  = (AliITS*) gAlice->GetModule("ITS");
   if (!ITS) { cout << "no ITS" << endl; return; }
   
   //AliITSgeom *aliitsgeo = ITS->GetITSgeom();
   AliITSgeom *geom = ITS->GetITSgeom();

   //Int_t cp[8]={0,0,0,0,0,0,0,0};

   cout << "SSD" << endl;

   AliITSDetType *iDetType=ITS->DetType(2);
   AliITSsegmentationSSD *seg2=(AliITSsegmentationSSD*)iDetType->GetSegmentationModel();
   AliITSresponseSSD *res2 = (AliITSresponseSSD*)iDetType->GetResponseModel();
   //res2->SetSigmaSpread(3.,2.);
   AliITSsimulationSSD *sim2=new AliITSsimulationSSD(seg2,res2);
   ITS->SetSimulationModel(2,sim2);

   TClonesArray *dig2  = ITS->DigitsAddress(2);
   TClonesArray *recp2  = ITS->ClustersAddress(2);
   //   AliITSClusterFinderSSD *rec2=new AliITSClusterFinderSSD(seg2,dig2,recp2);
   AliITSClusterFinderSSD *rec2=new AliITSClusterFinderSSD(seg2,dig2);
   ITS->SetReconstructionModel(2,rec2);
   // test
   printf("SSD dimensions %f %f \n",seg2->Dx(),seg2->Dz());
   printf("SSD nstrips %d %d \n",seg2->Npz(),seg2->Npx());

   
//
//   Loop over events
//


   Int_t Nh=0;
   Int_t Nh1=0;
   for (int nev=0; nev<= evNumber2; nev++) {
     Int_t nparticles = 0;
     nparticles = gAlice->GetEvent(nev);
     cout << "nev         " << nev <<endl;
     cout << "nparticles  " << nparticles <<endl;
     if (nev < evNumber1) continue;
     if (nparticles <= 0) return;
     
     AliITShit *itsHit;
     AliITSRecPoint *itsPnt = 0;
     AliITSRawClusterSSD *itsClu = 0;
     
     // Get Hit, Cluster & Recpoints Tree Pointers

     TTree *TH = gAlice->TreeH();
     Int_t nenthit=TH->GetEntries();
     printf("Found %d entries in the Hit tree (must be one per track per event!)\n",nenthit);

     ITS->GetTreeC(nev);
     TTree *TC=ITS->TreeC();
     Int_t nentclu=TC->GetEntries();
     printf("Found %d entries in the Cluster tree (must be one per module per event!)\n",nentclu);

     TTree *TR = gAlice->TreeR();
     Int_t nentrec=TR->GetEntries();
     printf("Found %d entries in the RecPoints tree\n",nentrec);

     // Get Pointers to Clusters & Recpoints TClonesArrays

     TClonesArray *ITSclu  = ITS->ClustersAddress(2); 
     printf ("ITSclu %p \n",ITSclu);
     TClonesArray *ITSrec  = ITS->RecPoints(); 
     printf ("ITSrec %p \n",ITSrec);

     // check recpoints

     //Int_t nbytes = 0;
     Int_t totpoints = 0;
     Int_t totclust = 0;

     // check hits
     
     Int_t nmodules=0;
     Int_t mod;
     
     ITS->InitModules(-1,nmodules); 
     ITS->FillModules(nev,0,nmodules,"","");
     
     TObjArray *fITSmodules = ITS->GetModules();
     
     Int_t first0 = geom->GetStartDet(0);  // SPD
     Int_t last0 = geom->GetLastDet(0);    // SPD
     Int_t first1 = geom->GetStartDet(1);  // SDD
     Int_t last1 = geom->GetLastDet(1);    // SDD
     Int_t first2 = geom->GetStartDet(2);  // SSD
     Int_t last2 = geom->GetLastDet(2);    // SSD

     //  For the SPD: first0 = 0, last0 = 239     (240 modules);  
     //  for the SDD: first1 = 240, last1 = 499   (260 modules);  
     //  for the SSD: first2 = 500, last2 = 2269  (1770 modules).  

     printf("det type %d first0, last0 %d %d \n",0,first0,last0);
     printf("det type %d first1, last1 %d %d \n",1,first1,last1);
     printf("det type %d first2, last2 %d %d \n",2,first2,last2);

     // module loop for the SSD
     for (mod=first2; mod<last2+1; mod++) {  // for the "ALL" option
     //for (mod=0; mod<last2-first2+1; mod++) { //for the "SSD" option

       TTree *TR = gAlice->TreeR();
       Int_t nentrec=TR->GetEntries();
       //printf("Found %d entries in the RecPoints tree\n",nentrec);
      
              //cout << "CLUSTERS: reset" << endl;
       ITS->ResetClusters();
       //cout << "CLUSTERS: get" << endl;
       TC->GetEvent(mod);
       //cout << "RECPOINTS: reset" << endl;
       ITS->ResetRecPoints();
       //cout << "RECPOINTS: get" << endl;
       //TR->GetEvent(mod+1);   // for the V3.04 AliRoot
       TR->GetEvent(mod);       // for the V3.05 AliRoot

       Int_t nrecp = ITSrec->GetEntries();
       totpoints += nrecp;
       if (nrecp) printf("Found %d rec points for module %d\n",nrecp,mod);
       //if (!nrecp) continue;
       Int_t nclusters = ITSclu->GetEntries();
       totclust += nclusters;
       //if (nclusters) printf("Found %d clusters for module %d\n",nrecc,mod);
       
       //AliITSmodule *Mod = (AliITSmodule *)fITSmodules->At(mod+first2);
       // for the "SSD" option

       AliITSmodule *Mod = (AliITSmodule *)fITSmodules->At(mod);
       // for the "ALL" option

       //       printf("Mod: %X\n",Mod);
       Int_t nhits = Mod->GetNhits();
       Float_t epart = 0;
       //cout <<" module,nrecp,nclusters,nhits ="<<mod<<","<<nrecp<<","<<nclusters<<","<<nhits<< endl;

       // ---------------- cluster/hit analysis ---------------------


     Float_t pathInSSD = 300.;

       // ---- Recpoint loop
       for (Int_t pnt=0;pnt<nrecp;pnt++) {
	 itsPnt  = (AliITSRecPoint*)ITSrec->At(pnt);
	 if(!itsPnt) continue;
	 itsClu  = (AliITSRawClusterSSD*)ITSclu->At(pnt);
	 if(!itsClu) continue;

	 Int_t nxP = itsClu->fMultiplicity;
	 Int_t nxN = itsClu->fMultiplicityN;
	 Float_t qclP = itsClu->fSignalP;     // in ADC
	 Float_t qclN = itsClu->fSignalN;     // in ADC
	 //Float_t dq = qclP - qclN;
	 Float_t qcut = itsClu->fQErr;        // abs(dq)/signal,
	                                      // where signal is
	                                      // max of qclP,qclN        
 	 Float_t xrec = 10000*itsPnt->GetX();
	 Float_t zrec = 10000*itsPnt->GetZ();
	 Float_t qrec = itsPnt->GetQ();      // in ADC, maximum from fSignalP/N
	 //Float_t dedx = itsPnt->GetdEdX();   // in KeV (ADC * 2.16)
	 Float_t dedx = itsPnt->fdEdX;   // in KeV (ADC * 2.16)
         Int_t ii = 0;
	 Int_t tr1 = itsPnt->GetLabel(ii);
         Int_t ii = 1;
	 Int_t tr2 = itsPnt->GetLabel(ii);
         Int_t ii = 2;
	 Int_t tr3 = itsPnt->GetLabel(ii);

	 // fill ntuple2
	     ntuple2_st.nxP = nxP;
             ntuple2_st.nxN = nxN;
	     ntuple2_st.x = xrec/1000;
             ntuple2_st.z = zrec/1000;

             if(qcut < 0.18) ntuple2->Fill();


	  Int_t noverlaps = 0;
	  Int_t noverprim = 0;
 	  Int_t flaghit = 0;
          Float_t xhit0 = 1e+7;
          Float_t yhit0 = 1e+7;
          Float_t zhit0 = 1e+7;

       // Hit loop
        for (Int_t hit=0;hit<nhits;hit++) {

	 itsHit   = (AliITShit*)Mod->GetHit(hit);

	 Int_t flagtrack = 0;
	 Int_t hitlayer = itsHit->GetLayer();
	 Int_t hitladder= itsHit->GetLadder();
	 Int_t hitdet= itsHit->GetDetector();

	 Int_t track = itsHit->GetTrack();
	 Int_t dray = 0;
	 Int_t hitstat = itsHit->GetTrackStatus();

 	  Float_t zhit = 10000*itsHit->GetZL();
	  Float_t xhit = 10000*itsHit->GetXL();
	  Float_t yhit = 10000*itsHit->GetYL();
	  Float_t ehit = 1.0e+6*itsHit->GetIonization(); // hit energy, KeV 

	   Int_t parent = itsHit->GetParticle()->GetFirstMother();
	   Int_t partcode = itsHit->GetParticle()->GetPdgCode();

   //  partcode (pdgCode): 11 - e-, 13 - mu-, 22 - gamma, 111 - pi0, 211 - i+
   //  310 - K0s, 321 - K+, 2112 - n, 2212 - p, 3122 - lambda

           Float_t pmod = itsHit->GetParticle()->P(); // the momentum at the
	                                              // vertex
	   pmod *= 1.0e+3;

	  if(hitstat == 66 && yhit < -146.) {  // entering hit
	    xhit0 = xhit;
	    yhit0 = yhit;
	    zhit0 = zhit;
	  }

	  if(hitstat == 66) continue; // Take the not entering hits only 

	  if(xhit0 > 9e+6 || zhit0 > 9e+6 || yhit0 > 9e+6) {
	    //cout<<"default xhit0,zhit0,yhit0 ="<<xhit0<<","<<zhit0<<","<<yhit0<<endl;
	    continue;
	  }



	  // Consider the hits only with the track number equaled to one
	  // of the recpoint
	  if(track == tr1) flagtrack = 1;

         if(flagtrack == 1) {     // the hit corresponds to the recpoint

	   flaghit = 1;

	   //Float_t px = itsHit->GetPXL(); // the momenta at this GEANT point
	   //Float_t py = itsHit->GetPYL();
	   //Float_t pz = itsHit->GetPZL();

         Int_t hitprim = 0;

	 if(partcode == 11 && pmod < 6) dray = 1; // delta ray is e-
	                                          // at p < 6 MeV/c

         if((hitstat == 68 || hitstat == 33) && dray == 0)  noverlaps=noverlaps + 1;
                                                  // overlapps for all hits but
	                                          // not for delta ray which
	                                          // also went out from the
	                                          // detector and returned
	                                          // again


	  // x,z resolution colculation
          if(hitstat == 68 || hitsat == 33) {
  	     Float_t xmed = (xhit + xhit0)/2;
	     Float_t zmed = (zhit + zhit0)/2;
	     Float_t xdif = xmed - xrec;
	     Float_t zdif = zmed - zrec;

            if(parent < 0)  {
	      hitprim = 1; // hitprim=1 for the primery particles
	      noverprim += 1;
	    }
	     pathInSSD = TMath::Sqrt((xhit0-xhit)*(xhit0-xhit)+(yhit0-yhit)*(yhit0-yhit)+(zhit0-zhit)*(zhit0-zhit));

	     //cout<<"lay,pnt,hit,xmed,xrec,xdif,zmed,zrec,zdif ="<<hitlayer<<","<<pnt<<","<<hit<<","<<xmed<<","<<xrec<<","<<xdif<<","<<zmed<<","<<zrec<<","<<zdif<<endl;

	 // fill ntuple
             ntuple_st.lay = hitlayer;
	     ntuple_st.nxP = nxP;
             ntuple_st.nxN = nxN;
	     ntuple_st.hitprim = hitprim;
             ntuple_st.partcode = partcode;
	     ntuple_st.x = xrec/1000;
             ntuple_st.z = zrec/1000;
	     ntuple_st.dx = xdif;
             ntuple_st.dz = zdif;
             ntuple_st.pmod = pmod;

             //if(qcut < 0.18) ntuple->Fill();
             ntuple->Fill();

	     //if(hitlayer == 5 && qcut < 0.18) {
	    if(hitlayer == 5 ) {
             Xres5->Fill(xdif);
             Zres5->Fill(zdif);
             Path5->Fill(pathInSSD);
	    }
            //if(hitlayer == 6 && qcut < 0.18) {
            if(hitlayer == 6) {
             Xres6->Fill(xdif);
             Zres6->Fill(zdif);
             Path6->Fill(pathInSSD);
            }
	  } // hitstat 68/33
	 } else {       // non correspondent hit
	  xhit0 = 1e+7;
	  zhit0 = 1e+7;
	 } // end of hit-recpoint correspondence
	} // hit loop       

	if(flaghit == 1) {

	  if(noverlaps == 0) noverlaps = 1; // cluster contains one or more
	  // delta rays only

	  // fill ntuple1
	  ntuple1_st.lay = hitlayer;
	  ntuple1_st.lad = hitladder;
	  ntuple1_st.det = hitdet;
	  ntuple1_st.nxP = nxP;
	  ntuple1_st.nxN = nxN;
	  ntuple1_st.qclP = qclP*300/pathInSSD; 
	  ntuple1_st.qclN = qclN*300/pathInSSD; 
	  ntuple1_st.qrec = qrec*300/pathInSSD; 
	  ntuple1_st.dx = xdif;
	  ntuple1_st.dz = zdif;
	  noverlaps -= 1;
	  noverprim -= 1;
	  ntuple1_st.noverlaps = noverlaps;
	  ntuple1_st.noverprim = noverprim;

	  //if(qcut < 0.18) ntuple1->Fill();
	  ntuple1->Fill();

          Float_t de = dedx*300./pathInSSD;
          dEdX->Fill(de);
	    if(hitlayer == 5 ) {
             adcPadcN5all->Fill(qclP,qclN);
            }
	    if(hitlayer == 6 ) {
             adcPadcN6all->Fill(qclP,qclN);
            }
	    if(hitlayer == 5 && qcut < 0.18) {
             adcPadcN5cut->Fill(qclP,qclN);
             NxP5->Fill(nxP);
             NxN5->Fill(nxN);
            }
	    if(hitlayer == 6 && qcut < 0.18) {
             adcPadcN6cut->Fill(qclP,qclN);
             NxP6->Fill(nxP);
             NxN6->Fill(nxN);
            }
	} // flaghit = 1
       } //b.b. recpoint loop
     } //b.b. module loop
   } //b.b. evnt loop

   TFile fhistos("SSD_his.root","RECREATE");

   ntuple->Write();
   ntuple1->Write();
   ntuple2->Write();
   NxP5->Write();
   NxN5->Write();
   NxP6->Write();
   NxN6->Write();
   Xres5->Write();
   Zres5->Write();
   Xres6->Write();
   Zres6->Write();
   Path5->Write();
   Path6->Write();
   adcPadcN5all->Write();
   adcPadcN6all->Write();
   adcPadcN5cut->Write();
   adcPadcN6cut->Write();
   dEdX->Write();

   fhistos.Close();

   cout<<"!!! Histogramms and ntuples were written"<<endl;

   TCanvas *c1 = new TCanvas("c1","ITS clusters",400,10,600,700);
   c1->Divide(2,2);
   c1->cd(1);
   gPad->SetFillColor(33);
         Xres5->SetFillColor(42);
         Xres5->Draw();
   c1->cd(2);
   gPad->SetFillColor(33);
         Zres5->SetFillColor(46);
         Zres5->Draw();
   c1->cd(3);
   gPad->SetFillColor(33);
         Xres6->SetFillColor(42);
         Xres6->Draw();
   c1->cd(4);
   gPad->SetFillColor(33);
         Zres6->SetFillColor(46);
         Zres6->Draw();

   cout<<"END  test for clusters and hits "<<endl;

}






