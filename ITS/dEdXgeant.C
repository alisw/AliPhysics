
Int_t max_modules=2269;

void dEdXyy (Int_t evNumber1=0,Int_t evNumber2=0,AliITSPid* pid=0,int mtest=0)
{
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   } else {
      delete gAlice;
      gAlice=0;
   }

// Connect the Root Galice file containing Geometry, Kine and Hits

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   printf("file %p\n",file);
   if (file) file->Close(); 
   file = new TFile("galice.root","UPDATE");
   file->ls();

   printf ("I'm after Map \n");

// Get AliRun object from file or create it if not on file

   if (!gAlice) {
       gAlice = (AliRun*)file->Get("gAlice");
       if (gAlice) printf("AliRun object found on file\n");
       if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
   printf ("I'm after gAlice \n");

 // Create Histogramms -------------------------------------------------

     TH1F *dedxSDD = new TH1F("dedxSDD","Particle energy (KeV) for layer 3 and 4",100,0.,400.);
     //TH1F *dedxSDDnorm = new TH1F("dedxSDDnorm","Particle energy (KeV) for layer 3 and 4",100,0.,400.);
     TH1F *dedxSDDnorm = new TH1F("dedxSDDnorm","Particle energy (KeV) for layer 3 and 4",100,0.,3000.);
     //TH1F *dedxSDDmip = new TH1F("dedxSDDmip","Particle energy (mips) for layer 3 and 4",100,0.,6.);
     TH1F *dedxSDDmip = new TH1F("dedxSDDmip","Particle energy (mips) for layer 3 and 4",100,0.,30.);
     //TH1F *signalSDD = new TH1F("signalSDD","SDD signal (ADC) for layer 3 and 4",100,0.,2000.);
     TH1F *signalSDD = new TH1F("signalSDD","SDD signal (ADC) for layer 3 and 4",100,0.,15000.);
     //TH1F *signalSDDmip = new TH1F("signalSDDmip","SDD signal (mips) for layer 3 and 4",100,0.,6.);
     TH1F *signalSDDmip = new TH1F("signalSDDmip","SDD signal (mips) for layer 3 and 4",100,0.,30.);
     TH1F *dedxSSD = new TH1F("dedxSSD","Particle energy (KeV) for layer 5 and 6",100,0.,400.);
     //TH1F *dedxSSDnorm = new TH1F("dedxSSDnorm","Particle energy (KeV) for layer 5 and 6",100,0.,400.);
     TH1F *dedxSSDnorm = new TH1F("dedxSSDnorm","Particle energy (KeV) for layer 5 and 6",100,0.,3000.);
     //TH1F *dedxSSDmip = new TH1F("dedxSSDmip","Particle energy (mips) for layer 5 and 6",100,0.,6.);
     TH1F *dedxSSDmip = new TH1F("dedxSSDmip","Particle energy (mips) for layer 5 and 6",100,0.,30.);
     //TH1F *signalSSD = new TH1F("signalSSD","SSD signal (ADC) for layer 5 and 6",100,0.,200.);
     TH1F *signalSSD = new TH1F("signalSSD","SSD signal (ADC) for layer 5 and 6",100,0.,2000.);
     //TH1F *signalSSDmip = new TH1F("signalSSDmip","SSD signal (mips) for layer 5 and 6",100,0.,6.);
     TH1F *signalSSDmip = new TH1F("signalSSDmip","SSD signal (mips) for layer 5 and 6",100,0.,30.);



     /*
     TH1F *dedxSDD = new TH1F("dedxSDD","Particle energy (KeV) for layer 3 and 4",100,0.,4000.);
     TH1F *dedxSDDnorm = new TH1F("dedxSDDnorm","Particle energy (KeV) for layer 3 and 4",100,0.,4000.);
     TH1F *dedxSDDmip = new TH1F("dedxSDDmip","Particle energy (mips) for layer 3 and 4",100,0.,40.);
     TH1F *signalSDD = new TH1F("signalSDD","SDD signal (ADC) for layer 3 and 4",100,0.,20000.);
     TH1F *signalSDDmip = new TH1F("signalSDDmip","SDD signal (mips) for layer 3 and 4",100,0.,40.);
     TH1F *dedxSSD = new TH1F("dedxSSD","Particle energy (KeV) for layer 5 and 6",100,0.,4000.);
     TH1F *dedxSSDnorm = new TH1F("dedxSSDnorm","Particle energy (KeV) for layer 5 and 6",100,0.,4000.);
     TH1F *dedxSSDmip = new TH1F("dedxSSDmip","Particle energy (mips) for layer 5 and 6",100,0.,40.);
     TH1F *signalSSD = new TH1F("signalSSD","SSD signal (ADC) for layer 5 and 6",100,0.,2000.);
     TH1F *signalSSDmip = new TH1F("signalSSDmip","SSD signal (mips) for layer 5 and 6",100,0.,40.);
     */

     TH1F *pathSDD = new TH1F("pathSDD","Path length in the SDD",100,0.,1000.);
     TH1F *pathSSD = new TH1F("pathSSD","Path length in the SSD",100,0.,1000.);

 // -------------------------------------------------------------
   AliITS *ITS  = (AliITS*) gAlice->GetModule("ITS");
   if (!ITS) return;
   AliITSgeom *geom = ITS->GetITSgeom();

   // Set the models for cluster finding

   // SPD

   ITS->MakeTreeC();
   Int_t nparticles=gAlice->GetEvent(0);

   AliITSDetType *iDetType=ITS->DetType(0);
   AliITSsegmentationSPD *seg0=(AliITSsegmentationSPD*)iDetType->GetSegmentationModel();
   TClonesArray *dig0  = ITS->DigitsAddress(0);
   TClonesArray *recp0  = ITS->ClustersAddress(0);
   AliITSClusterFinderSPD *rec0=new AliITSClusterFinderSPD(seg0,dig0,recp0);
   ITS->SetReconstructionModel(0,rec0);

   // SDD

   AliITSDetType *iDetType=ITS->DetType(1);
   AliITSgeom *geom = ITS->GetITSgeom();

   AliITSsegmentationSDD *seg1=(AliITSsegmentationSDD*)iDetType->GetSegmentationModel();
   if (!seg1) seg1 = new AliITSsegmentationSDD(geom);
   AliITSresponseSDD *res1 = (AliITSresponseSDD*)iDetType->GetResponseModel();
   if (!res1) res1=new AliITSresponseSDD();

   Float_t baseline,noise;
   res1->GetNoiseParam(noise,baseline);
   Float_t noise_after_el = res1->GetNoiseAfterElectronics();
   Float_t thres = baseline;
   thres += (4.*noise_after_el);  // TB // (4.*noise_after_el);
   printf("thres %f\n",thres);
   //res1->Print();

   TClonesArray *dig1  = ITS->DigitsAddress(1);
   TClonesArray *recp1  = ITS->ClustersAddress(1);
   AliITSClusterFinderSDD *rec1=new AliITSClusterFinderSDD(seg1,res1,dig1,recp1);
   rec1->SetCutAmplitude((int)thres);
   ITS->SetReconstructionModel(1,rec1);
   //rec1->Print();

   // SSD

   AliITSDetType *iDetType=ITS->DetType(2);
   AliITSsegmentationSSD *seg2=(AliITSsegmentationSSD*)iDetType->GetSegmentationModel();
   TClonesArray *dig2  = ITS->DigitsAddress(2);
   AliITSClusterFinderSSD *rec2=new AliITSClusterFinderSSD(seg2,dig2);
   ITS->SetReconstructionModel(2,rec2);
   
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
     AliITShit *itsHitPrev;
     AliITSRecPoint *itsPnt = 0;
     AliITSRawClusterSDD *itsClu = 0;
     
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

     TClonesArray *ITSclu  = ITS->ClustersAddress(1); 
     printf ("ITSclu %p \n",ITSclu);
     //            Int_t ncl = ITSclu->GetEntries();
     //     cout<<"ncluster ="<<ncl<<endl;
     TClonesArray *ITSrec  = ITS->RecPoints(); 
     printf ("ITSrec %p \n",ITSrec);

     // check recpoints

     //Int_t nbytes = 0;
     Int_t totpoints = 0;
     Int_t totclust = 0;

     // check hits
     
     Int_t nmodules=0;
     
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

     Int_t mod;
     printf("det type %d first0, last0 %d %d \n",0,first0,last0);
     printf("det type %d first1, last1 %d %d \n",1,first1,last1);
     printf("det type %d first2, last2 %d %d \n",2,first2,last2);
       Int_t negtrSDD = 0;
       Int_t negtrSSD = 0;

    for (mod=0; mod<last2+1; mod++) {
      if(mod>max_modules)continue;
      if(mod < first1) continue;  // for the SDD/SSD only

       AliITSmodule *Mod = (AliITSmodule *)fITSmodules->At(mod);
       Int_t nhits = Mod->GetNhits();
       Int_t layer = Mod->GetLayer();
       Int_t ladder = Mod->GetLadder();
       Int_t detector = Mod->GetDet();
       Int_t hit;
       Int_t parent;
       Int_t dray;
       Int_t hitstat;
       Int_t partcode;
       Int_t hitprim;
       Float_t epart;
       Float_t pathInSi=300;
       Float_t xhit;
       Float_t yhit;
       Float_t zhit;
       Float_t px;
       Float_t py;
       Float_t pz;
       Float_t pmod;
       Float_t xhit0 = 1e+7;
       Float_t yhit0 = 1e+7;
       Float_t zhit0 = 1e+7;

       TTree *TR = gAlice->TreeR();
       ITS->ResetClusters();
       //cout << "CLUSTERS: get" << endl;
       TC->GetEvent(mod);
       //cout << "RECPOINTS: reset" << endl;
       ITS->ResetRecPoints();
       //cout << "RECPOINTS: get" << endl;
       //TR->GetEvent(mod+1);
       TR->GetEvent(mod);
   
       Int_t nrecp = ITSrec->GetEntries();
       totpoints += nrecp;
       if (!nrecp) continue;

       Int_t trackRecp[3];
       Int_t itr;
       Int_t TrackRecPoint;
       Float_t PmodRecP;
       Float_t signalRP;

//         cout <<" module,layer,ladder,detector,nrecp,nhits ="<<
//  	 mod<<","<<layer<<","<<ladder<<","<<detector<<","<<nrecp<<","<<nhits<< endl;
       cout<<".";

       // ---------- RecPoint signal analysis for the SDD/SSD  ---------

       for (Int_t pnt=0;pnt<nrecp;pnt++) {          // recpoint loop

	 itsPnt  = (AliITSRecPoint*)ITSrec->At(pnt);

	 Int_t RecPointPrim = 0;
	 if(!itsPnt) continue;

	   signalRP = itsPnt->GetQ();
	   trackRecp[0] = itsPnt->GetLabel(0);
	   trackRecp[1] = itsPnt->GetLabel(1);
	   trackRecp[2] = itsPnt->GetLabel(2);

//  	   cout<<"New Recp: pnt,signal,tr0,tr1,tr2 ="<<pnt
//  	       <<","<<signalRP<<","<<trackRecp[0]<<","<<trackRecp[1]<<","<<trackRecp[2]<<endl;

	   /*
            if(trackRecp[0]<0&&trackRecp[1]<0&&trackRecp[2]<0) {
cout<<"pnt,tr0,1,2 ="<<pnt<<","<<trackRecp[0]<<","<<trackRecp[1]<<","<<trackRecp[2]<<endl;
 if(mod<first2) negtrSDD += 1;
 if(mod>=first2) negtrSSD += 1;
	    }
	   */

          xhit0 = 1e+7;
          yhit0 = 1e+7;
          zhit0 = 1e+7;

       // Hit loop
        for (hit=0;hit<nhits;hit++) {

	 itsHit   = (AliITShit*)Mod->GetHit(hit);

	 hitprim = 0;
	 dray = 0;
	 Int_t hitRecp = 0;
	 Int_t trackHit = itsHit->GetTrack();
	 hitstat = itsHit->GetTrackStatus();
 	 zhit = 10000*itsHit->GetZL();
	 xhit = 10000*itsHit->GetXL();
	 yhit = 10000*itsHit->GetYL();
	 Float_t dEn = 1.0e+6*itsHit->GetIonization(); // hit energy, KeV 

	 //	 cout<<" New hit: hit,track,hitstat,yhit,ehit ="<<hit<<","<<trackHit<<","<<hitstat<<","<<yhit<<","<<dEn<<endl;

	 // check the primary particle
	 partcode = itsHit->GetParticle()->GetPdgCode();
	 parent = itsHit->GetParticle()->GetFirstMother();
	 if(parent < 0) hitprim = 1; // primery particle



	 // select the primery hits with the track number equaled to 
	 // the one of the RecPoint
	 for(itr=0;itr<3;itr++) {
	   if(trackRecp[itr] == trackHit) hitRecp = 1; 
	   if(trackRecp[itr] == trackHit && hitprim == 1) {
             RecPointPrim = 1; 
	     TrackRecPoint = trackRecp[itr]; 
	   }
	   if(trackRecp[itr] == trackHit && hitprim == 0) trackRecp[itr]=-100; 
	 }
         
	 if(hitRecp != 1) continue; // hit doesn't correspond to the recpoint

	 //cout<<"Hit Ok: trackRecp0,1,2,hitRecp,RecPointPrim,TrackRecPoint ="<<trackRecp[0]<<","<<trackRecp[1]<<","<<trackRecp[2]<<","<<hitRecp<<","<<RecPointPrim<<","<<TrackRecPoint<<endl;

	 //	 if(hitstat == 66 && yhit < -136.) {  // entering hit
	 if(hitstat == 66 && hitprim == 1) {  // entering hit
	    xhit0 = xhit;
	    yhit0 = yhit;
	    zhit0 = zhit;
	 }
	 //	 	 cout<<" xhit0,zhit0 ="<<xhit0<<","<<zhit0<<","<<endl;
         
	  if(hitstat == 66) continue; // Take only the hits with the not
                                      // zero energy

	  if(xhit0 > 9e+6 || zhit0 > 9e+6) {
	    // cout<<"default xhit0,zhit0 ="<<xhit0<<","<<zhit0<<","<<endl;
	    continue;
	  }

	 // check the particle code (type)
	 //	   Int_t parent = itsHit->GetParticle()->GetFirstMother();
	  //partcode = itsHit->GetParticle()->GetPdgCode();

   //  partcode (pdgCode): 11 - e-, 13 - mu-, 22 - gamma, 111 - pi0, 211 - i+
   //  310 - K0s, 321 - K+, 2112 - n, 2212 - p, 3122 - lambda

	   /*
	   px = itsHit->GetPXL();  // momentum for the hit
	   py = itsHit->GetPYL();
	   pz = itsHit->GetPZL();
	   pmod = 1000*sqrt(px*px+py*py+pz*pz);
	   */
	   
           pmod = itsHit->GetParticle()->P(); // momentum at the vertex
	   pmod *= 1.0e+3;

	   // cout<<"track,partcode,pmod,hitprim ="<<trackHit<<","<<partcode<<","<<pmod<<","<<hitprim<<endl;

	 if(partcode == 11 && pmod < 6) dray = 1; // delta ray is e-
	                                          // at p < 6 MeV/c

  // find the primery particle path length in Si and pmod (MeV/c) 
          if((hitstat == 68 || hitstat == 33) && hitprim == 1)  {
	     pathInSi = TMath::Sqrt((xhit0-xhit)*(xhit0-xhit)+(yhit0-yhit)*(yhit0-yhit)+(zhit0-zhit)*(zhit0-zhit));
             PmodRecP = itsHit->GetParticle()->P();
	     PmodRecP *= 1.0e+3;
	     // cout<<"path in Si="<<pathInSi<<endl;

	 Float_t Signal = signalRP*(300/pathInSi)/38; 
	    if(Signal<0.5) {
	      //cout<<"track,partcode,pmod,hitprim,hitstat,Signal,dEn ="<<trackHit<<","<<partcode<<","<<pmod<<","<<hitprim<<","<<hitstat<<","<<Signal<<","<<dEn<<endl;
	    }
	  }
	} // hit loop

	if(RecPointPrim == 1) {
        
	  //cout<<" SDD/SSD RecPoints: lay,lad,det,track,pmod,signal,path ="<<layer<<","<<ladder<<","<<detector<<","<<TrackRecPoint<<","<<PmodRecP<<","<<signalRP<<","<<pathInSi<<endl;

	 Int_t xpartcode=gAlice->Particle(TrackRecPoint)->GetPdgCode();
         
	 if(mod < first2) {            // SDD
	  signalRP *= (280./pathInSi); // 280 microns of Si thickness    
	  signalSDD->Fill(signalRP); 
	  signalSDDmip->Fill(signalRP/280.); // ADC/280 -> mips
	 InPID(pid,(Int_t)nev,(Int_t)TrackRecPoint, signalRP/280. ,
	       (Float_t)PmodRecP,(Int_t)xpartcode);
	  //cout<<" SDD RecPoints: lay,lad,det,track,pmod,signal,path ="<<layer<<","<<ladder<<","<<detector<<","<<TrackRecPoint<<","<<PmodRecP<<","<<signalRP<<","<<pathInSi<<endl;
	  if(pathInSi < 100 || pathInSi > 2000) cout<<" No normal pathInSi in SDD ="<<pathInSi<<endl;

	  //if(signalRP/280 < 5) cout<<" small signalSDDmip, path ="<<signalRP/280<<","<<pathInSDD<<endl;
	 }else{                        // SSD
	  signalRP *= (300/pathInSi); // 300 microns of Si thickness
	  //cout<<" SSD RecPoints: lay,lad,det,track,pmod,signal,path ="<<layer<<","<<ladder<<","<<detector<<","<<TrackRecPoint<<","<<PmodRecP<<","<<signalRP<<","<<pathInSi<<endl;
	  signalSSD->Fill(signalRP); 
	  signalSSDmip->Fill(signalRP/38.); // ADC/38 -> mips
	  InPID(pid,(Int_t)nev,(Int_t)TrackRecPoint, signalRP/38. ,
		(Float_t)PmodRecP,(Int_t)xpartcode);
	  if(pathInSi < 100 || pathInSi > 2000) cout<<" No normal pathInSi in SSD ="<<pathInSi<<endl;

	  if(signalRP/38 < 0.5) cout<<" small signalSSD (mip), path, module ="<<signalRP/38<<","<<pathInSi<<mod<<endl;
	 }
	} // primery particle
       } // pnt, recpoint loop       
       
       // dEdX hit analysis for the SDD/SSD

         Int_t track = -100;
         Int_t trackPrev = -100;
         parent = -100;
         Int_t parentPrev = -100;
         Int_t flagprim = 0;
	 Int_t TrackPart;
         epart = 0;
         pathInSi = 1.0e+7;
	 Float_t PmodPart;
	 zhit0 = 1.0e+7;
	 xhit0 = 1.0e+7;
	 yhit0 = 1.0e+7;

	 /*
	 if(mod <= 259) {
	 cout<<" SDD hits: nhits ="<<nhits<<endl;
	 }else{
	 cout<<" SSD hits: nhits ="<<nhits<<endl;
	 }
	 */
        
       for (hit=0;hit<nhits;hit++) {
        
	 itsHit   = (AliITShit*)Mod->GetHit(hit);
	 if(hit>0) itsHitPrev   = (AliITShit*)Mod->GetHit(hit-1);
         hitprim = 0;
         Int_t hitprimPrev = 0;
	 track = itsHit->GetTrack();
	 if(hit > 0) Int_t trackPrev = itsHitPrev->GetTrack();
	 dray = 0;
	 hitstat = itsHit->GetTrackStatus();

 	  zhit = 10000*itsHit->GetZL();
	  xhit = 10000*itsHit->GetXL();
	  yhit = 10000*itsHit->GetYL();
	  Float_t ehit = 1.0e+6*itsHit->GetIonization(); // hit energy, KeV 

	  //	  	  cout<<"New hit,lay,hitstat,yhit,ehit ="<<hit<<","<<layer<<","<<hitstat<<","<<yhit<<","<<ehit<<endl;
	  //  cout<<"befor: track,trackPrev ="<<track<<","<<trackPrev<<endl;

	   parent = itsHit->GetParticle()->GetFirstMother();
	   if(hit>0) Int_t parentPrev = itsHitPrev->GetParticle()->GetFirstMother();
	   partcode = itsHit->GetParticle()->GetPdgCode();

   //  partcode (pdgCode): 11 - e-, 13 - mu-, 22 - gamma, 111 - pi0, 211 - i+
   //  310 - K0s, 321 - K+, 2112 - n, 2212 - p, 3122 - lambda

	   px = itsHit->GetPXL();
	   py = itsHit->GetPYL();
	   pz = itsHit->GetPZL();
	   pmod = 1000*sqrt(px*px+py*py+pz*pz);

	   //  cout<<"partcode,pmod,parent,parentPrev,hitprim,flagprim ="<<partcode<<","<<pmod<<","<<parent<<","<<parentPrev<<","<<hitprim<<","<<flagprim<<endl;


	 if(partcode == 11 && pmod < 6) dray = 1; // delta ray is e-
	                                          // at p < 6 MeV/c


          if(parent < 0 && parent > -90) hitprim += 1;
          if(parentPrev < 0 && parentPrev >-90) hitprimPrev += 1; 
         // hitprim=1 for the primery particles
	  //	  if(hit==0&&hitprim==1&&hitstat==66) {
	  if(hitprim==1&&hitstat==66) {
 	    zhit0 = zhit;
 	    xhit0 = xhit;
 	    yhit0 = yhit;
	  }

	  if((hitstat == 68 || hitstat == 33) && hitprim == 1) {
	     pathInSi = TMath::Sqrt((xhit0-xhit)*(xhit0-xhit)+(yhit0-yhit)*(yhit0-yhit)+(zhit0-zhit)*(zhit0-zhit));
	     TrackPart = track;
             PmodPart = itsHit->GetParticle()->P();
	     PmodPart *= 1.0e+3;
	     //  cout<<"xhit0,yhit0,zhit0,pathInSi ="<<xhit0<<","<<yhit0<<","<<zhit0<<","<<pathInSi<<endl;
	  }

	  if(hitprim == 0) track = parent;          
	  if(hitprimPrev == 0) trackPrev = parentPrev;          
	  //  cout<<"after: track,trackPrev ="<<track<<","<<trackPrev<<endl;
          
	  if(hit > 0 && track == trackPrev) epart += ehit;   
	  if((hit > 0 && track == trackPrev) && (hitprim==1)) flagprim +=1;

	  //  cout<<"new hitprim, flagprim ="<<hitprim<<","<<flagprim<<endl;

	  if((hit > 0 && track != trackPrev) || (hit == nhits-1)) {
	    if(flagprim > 0) {
	      if(layer > 2 && layer < 5) {
		if(epart < 40) cout<<" small SDD dedx ="<<epart<<endl;
	        if(pathInSi < 100 || pathInSi > 2000) cout<<" No normal pathSDD ="<<pathInSi<<endl;
		if(epart > 40 && pathInSi > 100 && pathInSi < 1.0e+5) {
                 dedxSDD->Fill(epart); 
	         epart *= (280/pathInSi);
                 dedxSDDnorm->Fill(epart); 
 	         dedxSDDmip->Fill(epart/75); 
		 //InPID(pid,(Int_t)nev,(Int_t)track, epart/75. ,(Float_t)pmod,(Int_t)partcode);
                 pathSDD->Fill(pathInSi); 
		 //  cout<<"SDD hit: lay,lad,det,track,pmod,epart,path ="<<layer<<","<<ladder<<","<<detector<<","<<TrackPart<<","<<PmodPart<<","<<epart<<","<<pathInSi<<endl;
		}
	      }
	      if(layer > 4) {
		if(epart < 40) cout<<" small SSD dedx ="<<epart<<endl;
	        if(pathInSi < 100 || pathInSi > 2000) cout<<" No normal pathSSD ="<<pathInSi<<endl;
		if(epart > 40 && pathInSi > 100 && pathInSi < 1.0e+5) {
                 dedxSSD->Fill(epart); 
		 epart *= (280/pathInSi);
		 dedxSSDnorm->Fill(epart); 
		 dedxSSDmip->Fill(epart/80); 
		 //InPID(pid,(Int_t)nev,(Int_t)track, epart/80. ,(Float_t)pmod,(Int_t)partcode);
                 pathSSD->Fill(pathInSi); 
		 //cout<<"SSD hit: lay,lad,det,track,pmod,epart,path ="<<layer<<","<<ladder<<","<<detector<<","<<TrackPart<<","<<PmodPart<<","<<epart<<","<<pathInSi<<endl;
		}
	      }
              flagprim=0;
              epart = 0;
	    }else{
	      //cout<<"No!, epart for the secondary"<<endl;
              epart = 0;
	    }                 
	  }
       } // SDD/SSD hit loop
    } // module loop
    //cout<<"negtrSDD, negtrSSD ="<<negtrSDD<<","<<negtrSSD<<endl;
   } // event loop
   
   cout<<endl;
     
   TFile fhistos("dEdX_his.root","RECREATE");
  
   dedxSDD->Write();
   dedxSDDnorm->Write();
   dedxSDDmip->Write();
   dedxSSD->Write();
   dedxSSDnorm->Write();
   dedxSSDmip->Write();
   signalSDD->Write();
   signalSDDmip->Write();
   signalSSD->Write();
   signalSSDmip->Write();
   pathSDD->Write();
   pathSSD->Write();

   fhistos.Close();
   cout<<"!!! Histogramms and ntuples were written"<<endl;
   // ------------------------------------------------------------

   TCanvas *c1 = new TCanvas("c1","ITS clusters",400,10,600,700);
   c1->Divide(2,2);


//     c1->cd(1);
//     gPad->SetFillColor(33);
//     signalSDD->SetFillColor(42);
//     signalSDD->Draw();

//     c1->cd(2);
//     gPad->SetFillColor(33);
//     signalSDDmip->SetFillColor(42);
//     signalSDDmip->Draw();

//     c1->cd(3);
//     gPad->SetFillColor(33);
//     dedxSDDnorm->SetFillColor(46);
//     dedxSDDnorm->Draw();

//     c1->cd(4);
//     gPad->SetFillColor(33);
//     dedxSDDmip->SetFillColor(46);
//     dedxSDDmip->Draw();

   c1->cd(1);
   gPad->SetFillColor(33);
   signalSSD->SetFillColor(42);
   signalSSD->Draw();

   c1->cd(2);
   gPad->SetFillColor(33);
   signalSSDmip->SetFillColor(42);
   signalSSDmip->Draw();

   c1->cd(3);
   gPad->SetFillColor(33);
   dedxSSDnorm->SetFillColor(46);
   dedxSSDnorm->Draw();

   c1->cd(4);
   gPad->SetFillColor(33);
   dedxSSDmip->SetFillColor(46);
   dedxSSDmip->Draw();
 

   /*
    c1->cd(1);
   gPad->SetFillColor(33);
   dedxSDD->SetFillColor(42);
   dedxSDD->Draw();

   c1->cd(2);
   gPad->SetFillColor(33);
   pathSDD->SetFillColor(42);
   pathSDD->Draw();

   c1->cd(3);
   gPad->SetFillColor(33);
   dedxSSD->SetFillColor(46);
   dedxSSD->Draw();

   c1->cd(4);
   gPad->SetFillColor(33);
   pathSSD->SetFillColor(46);
   pathSSD->Draw();
   */

//   if(!pid)pid->Tab();
   cout<<"END  test for clusters and hits "<<endl;

}
//============================ InPID() =========================================
int totpid;
void InPID (AliITSPid *pid=0,int     nev,
                             Int_t   track,
			     Float_t signal,
			     Float_t pmom,
			     Int_t   pcode)
{  
// Includes recp in PID  table.
    if(pid){
      //cout<<" In pid: track,sig,pmod,pcode="<<track<<","<<signal<<","<<pmom<<","<<pcode<<endl;                                                                                      
    if( abs(pcode)==321 || abs(pcode)==211 ||abs(pcode)==2212 ||abs(pcode)==11  )
	{
        pid->SetEdep(10000*nev+track,signal);   
	pid->SetPmom(10000*nev+track,pmom);                                     
	pid->SetPcod(10000*nev+track,abs(pcode));
	}                                     
    } 
    totpid++;                            
}
//=======================






