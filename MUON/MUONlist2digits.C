void MUONlist2digits (Int_t evNumber1=0,Int_t evNumber2=0,Int_t nCathode=1) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

// Dynamically link some shared libs


    if (gClassTable->GetID("AliRun") < 0) {
	gSystem->Load("$ALITOP/cern.so/lib/libpdfDUMMY.so");
	gSystem->Load("$ALITOP/cern.so/lib/libPythia.so");
	gSystem->Load("$ROOTSYS/lib/libEG.so");       
	gSystem->Load("$ROOTSYS/lib/libEGPythia.so");       
	gSystem->Load("libGeant3Dummy.so");    //a dummy version of Geant3
	gSystem->Load("PHOS/libPHOSdummy.so"); //the standard Alice classes 
	gSystem->Load("libgalice.so");         // the standard Alice classes 
    }

// Connect the Root Galice file containing Geometry, Kine and Hits

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (file) file->Close(); 
   file = new TFile("galice.root","UPDATE");
   file->ls();
   file->Map();

   printf ("I'm after Map \n");

// Get AliRun object from file or create it if not on file

   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
   printf ("I'm after gAlice \n");

       AliMUON *MUON  = gAlice->GetDetector("MUON");

   //   AliMUONchamber&  iChamber;
   AliMUONchamber*  iChamber;
   AliMUONsegmentation*  segmentation;

   const Float_t rmin[10] = {17.5,17.5,23.5,23.5,33.5,33.5,43,43,50,50};
   const Float_t rmax[10] = {91.5,91.5,122.5,122.5,173,173,221,221,256.5,256.5};


   Int_t Npx[10];
   Int_t Npy[10];


   for ( int i=0;i<10;i++) {
        Npx[i]=(Int_t)rmax[i]/0.75;
        Npy[i]=(Int_t)rmax[i]/0.5;
   }

       Int_t nxmax=1026;
       Int_t nymax=1026;
       AliMUONlist *elem[10][1026][1026];
       TObjArray *obj=new TObjArray;

       Int_t digits[3]; 
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
       Int_t Nc=0;
       //
       Int_t counter=0;

    for (int icat=0;icat<nCathode;icat++) { 

       //initialize the arrray of pointers !!!!!
       for (int i=0; i<10; i++) {
         for (int j=0; j<nymax; j++) {
           for (int k=0; k<nxmax; k++) {
               elem[i][j][k]=0;
           }
         }
       }

     
//
//   Loop over events
//
       for (Int_t track=0; track<ntracks;track++) {
	   gAlice->ResetHits();
	   Int_t nbytes += TH->GetEvent(track);
	   if (MUON)  {
	       for(AliMUONhit* mHit=(AliMUONhit*)MUON->FirstHit(-1); 
		   mHit;
		   mHit=(AliMUONhit*)MUON->NextHit()) 
	       {
		   Int_t   nch   = mHit->fChamber;  // chamber number
		   Float_t x     = mHit->fX;        // x-pos of hit
		   Float_t y     = mHit->fY;        // y-pos

		   iChamber = &(MUON->Chamber(nch-1));
		   response=iChamber->GetResponseModel();

		   Int_t nsec=iChamber->Nsec();
		   //                   Int_t rmin = (Int_t)iChamber->frMin;
		   //                   Int_t rmax = (Int_t)iChamber->frMax;
//
//
		   for (AliMUONcluster* mPad=(AliMUONcluster*)MUON->FirstPad(mHit);
		        mPad;
		        mPad=(AliMUONcluster*)MUON->NextPad())
		   {
	        	Int_t cathode     = mPad->fCathode;   // chamber number
		        Int_t nhit     = mPad->fHitNumber; // hit number
		        Int_t qtot     = mPad->fQ;         // charge
		        Int_t ipx      = mPad->fPadX;      // pad number on X
		        Int_t ipy      = mPad->fPadY;      // pad number on Y
		        Int_t iqpad    = mPad->fQpad;      // charge per pad
		        Int_t rpad     = mPad->fRpad;      // r-pos of pad
//
//
        if (cathode != (icat+1)) continue;
	segmentation=iChamber->GetSegmentationModel(cathode);
	//        Int_t Npx[nch-1]  = segmentation->Npx();
	//        Int_t Npy[nch-1]  = segmentation->Npy();

	if (cathode==1) {
             Int_t npx=Npx[nch-1];
             Int_t npy=Npy[nch-1];
        } else {
             Int_t npx=Npy[nch-1];
             Int_t npy=Npx[nch-1];
        }

	//         printf("nch, cathode, rmin, rmax, npx, npy %d %d %d %d %d %d\n",
	//	                  nch,cathode,rmin,rmax,npx,npy);

	//	       printf("icat, iqpad, ipx, ipy %d %d %d %d \n",
	//			         icat, iqpad, ipx, ipy);


	// check boundaries
        if (rpad < rmin[nch-1] || iqpad ==0 || rpad > rmax[nch-1]) continue;
	//         if (rpad < rmin[nch-1] || iqpad ==0 || rpad > 81.65) continue;
	 //         if (rpad < rmin || iqpad ==0 || rpad > rmax) continue;

	  // fill the info array
          TVector *trinfo;
           trinfo=new TVector(2);

          trinfo(0)=(Float_t)track;
          trinfo(1)=(Float_t)iqpad;

	  //          Int_t digits[3]; 
          digits[0]=ipx;
          digits[1]=ipy;
          digits[2]=iqpad;


          Int_t trk[50];
          Int_t chtrk[50];

    // build the list of fired pads and update the info
    AliMUONlist *pdigit=elem[nch-1][ipy+npy][ipx+npx];
    if (pdigit==0) {
       obj->AddAtAndExpand(new AliMUONlist(rpad,digits),counter);
       counter++;
       //       Int_t nentrobj=obj->GetEntriesFast();
       Int_t last=obj->GetLast();
       pdigit=(AliMUONlist*)obj->At(last);
       elem[nch-1][ipy+npy][ipx+npx]=pdigit;
       //       Int_t q=pdigit->fSignal;
       // list of tracks
       TObjArray *trlist=(TObjArray*)pdigit->TrackList();
       trlist->Add(trinfo);
     } else {
       // update charge
       (*pdigit).fSignal+=iqpad;
       // update list of tracks
       TObjArray* trlist=(TObjArray*)pdigit->TrackList();
       Int_t nptracks=trlist->GetEntriesFast();
        for (Int_t tr=0;tr<nptracks;tr++) {
           TVector *ptrk=(TVector*)trlist->At(tr);
           Int_t entries=ptrk->GetNrows();
           TVector &vtrk = *ptrk;
           trk[tr]=(Int_t)vtrk(0);
           chtrk[tr]=(Int_t)vtrk(1);
           if (trk[tr]==track) {
                chtrk[tr]+=iqpad;
                trlist->RemoveAt(tr);
                trinfo(0)=trk[tr];
                trinfo(1)=chtrk[tr];
                trlist->AddAt(trinfo,tr);
           } else {
                trlist->Add(trinfo);
           }
        } //end loop over list of tracks for one pad
     }  //  end if pdigit
    

		   } //end loop over clust

	       } // hit loop
	   } // if MUON
       } // track loop

    Int_t tracks[10];
    Int_t charges[10];
    cout<<"start filling digits "<<endl;
    //    Int_t digits[3];
    const Int_t zero_supm = 6;
    const Int_t adc_satm = 1024;
    Int_t nd=0;


    // start filling the digits
 for (Int_t id=0;id<10;id++) {
             Int_t nd=0;
	     //             Int_t npx=Npx[id];
	     //             Int_t npy=Npy[id];

	     //	     Int_t nx=2*npx;
	     //	     Int_t ny=2*npx;
   for (Int_t iy=0;iy<nymax;iy++) {
     for (Int_t ix=0;ix<nxmax;ix++) {
        AliMUONlist *address=elem[id][iy][ix];
        if (address==0) continue; 
        // do zero-suppression and signal truncation
        Int_t rpad=address->fRpad;
        Int_t q=address->fSignal;        
	if ( q <= zero_supm  || rpad > 55) continue;
	//        if ( q <= zero_supm ) continue;
        if ( q > adc_satm)  q=adc_satm;
        digits[0]=(*address).fPadX;
        digits[1]=(*address).fPadY;
        digits[2]=address->fSignal;

        TObjArray* trlist=(TObjArray*)address->TrackList();
        Int_t nptracks=trlist->GetEntriesFast();
        // this should be changed to accomodate the real number of tracks
        if (nptracks > 10) {
	    cout<<"Attention - nptracks > 3 "<<nptracks<<endl;
            nptracks=10;
        }
        for (Int_t tr=0;tr<nptracks;tr++) {
           TVector *pp=(TVector*)trlist->At(tr);
           Int_t entries=pp->GetNrows();
           TVector &v2 = *pp;
           tracks[tr]=(Int_t)v2(0);
           charges[tr]=(Int_t)v2(1);
        } //end loop over list of tracks for one pad

	// fill digits
         MUON->AddDigits(id+1,tracks,charges,digits);
         nd++;
     
     }  // end loop over ix
   }  // end loop over iy
    cout<<"I'm out of the loops over pads"<<endl;
    cout<<" id  nd  "<<id<<" "<<nd<<endl;

    //    CathodeIndex[id+icat*10] = nd;     // update cathodes index
   

 } //loop over chambers
    cout<<"I'm out of the loops for digitisation"<<endl;


        gAlice->TreeD()->Fill();
        TTree *TD=gAlice->TreeD();
        int ndig=TD->GetEntries();
        cout<<"number of digits  "<<ndig<<endl;

        int ndig1=(MUON->Dch1())->GetEntriesFast();
        cout<<"number of digits  1   "<<ndig1<<endl;
        int ndig2=(MUON->Dch2())->GetEntriesFast();
        cout<<"number of digits  2   "<<ndig2<<endl;
        int ndig3=(MUON->Dch3())->GetEntriesFast();
        cout<<"number of digits  3   "<<ndig3<<endl;
        int ndig4=(MUON->Dch4())->GetEntriesFast();
        cout<<"number of digits  4   "<<ndig4<<endl;
        int ndig5=(MUON->Dch5())->GetEntriesFast();
        cout<<"number of digits  5   "<<ndig5<<endl;
        int ndig6=(MUON->Dch6())->GetEntriesFast();
        cout<<"number of digits  6   "<<ndig6<<endl;
        int ndig7=(MUON->Dch7())->GetEntriesFast();
        cout<<"number of digits  7   "<<ndig7<<endl;
        int ndig8=(MUON->Dch8())->GetEntriesFast();
        cout<<"number of digits  8   "<<ndig8<<endl;
        int ndig9=(MUON->Dch9())->GetEntriesFast();
        cout<<"number of digits  9   "<<ndig9<<endl;
        int ndig10=(MUON->Dch10())->GetEntriesFast();
        cout<<"number of digits  10   "<<ndig10<<endl;

        MUON->ResetDigits();


	//        char hname[30];
	//        sprintf(hname,"TreeD%d",nev);
	//        gAlice->TreeD()->Write(hname);

    } //end loop over cathodes

        char hname[30];
        sprintf(hname,"TreeD%d",nev);
        gAlice->TreeD()->Write(hname);
    //    MUON->ResetDigits();
    //    if (CathodeIndex) CathodeIndex->Clear();

    //    file->Write();

     file->ls();
     file->Map();
   } // event loop 

   //     delete [] CathodeIndex;

     file->Close();
     cout<<"END  digitisation "<<endl;
   
}
