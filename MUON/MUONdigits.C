void MUONdigits (Int_t evNumber1=0,Int_t evNumber2=0,Int_t nCathode=1) 
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
       gSystem->Load("libGeant3Dummy.so");        //a dummy version of Geant3
       gSystem->Load("PHOS/libPHOSdummy.so");     //the standard Alice classes 
       gSystem->Load("libgalice.so");             //the standard Alice classes 
   }


// Connect the Root Galice file containing Geometry, Kine and Hits

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (file) file->Close(); 
   file = new TFile("galice.root","UPDATE");
   file->ls();
   //   file->Map();

   printf ("I'm after Map \n");

// Get AliRun object from file or create it if not on file

   if (!gAlice) {
       gAlice = (AliRun*)file->Get("gAlice");
       if (gAlice) printf("AliRun object found on file\n");
       if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
   printf ("I'm after gAlice \n");
   
   AliMUON *MUON  = gAlice->GetDetector("MUON");
   
   AliMUONchamber*  iChamber;
   AliMUONsegmentation*  segmentation;
   
   Int_t Npx[10];
   Int_t Npy[10];
   Int_t trk[50];
   Int_t chtrk[50];


   Int_t nxmax=1026;
   Int_t nymax=1026;
   AliMUONlist *elem[10*1026*1026];
   AliMUONlist **ppe=elem;
   TObjArray *obj=new TObjArray;

   Int_t digits[3]; 
   // fill the info array
   TVector *trinfo;
   trinfo=new TVector(2);
   
//
//   Loop over events 
//

   Int_t Nh=0;
   Int_t Nh1=0;
   for (int nev=0; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       cout << "nev         " <<nev<<endl;
       cout << "nparticles  " <<nparticles<<endl;
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;

       TTree *TH = gAlice->TreeH();
       Int_t ntracks = TH->GetEntries();
       Int_t Nc=0;
       //
       Int_t counter=0;

       // loop over cathodes
       for (int icat=0;icat<nCathode;icat++) { 

	   //initialize the arrray of pointers !!!!!
	   ppe=memset(elem,0,sizeof(void*)*10*1026*1026);

	   printf("Start loop over tracks \n");     
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
//
//
		       if (nch >10) continue;
		   
		       iChamber = &(MUON->Chamber(nch-1));
		       response=iChamber->GetResponseModel();

		       Int_t nsec=iChamber->Nsec();
		       Int_t rmin = (Int_t)iChamber->frMin;
		       Int_t rmax = (Int_t)iChamber->frMax;
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
			   Int_t izone    = mPad->fRSec;      // r-pos of pad
//
//
			   if (cathode != (icat+1)) continue;
			   segmentation=iChamber->GetSegmentationModel(cathode);
			   Int_t Npx[nch-1]  = segmentation->Npx();
			   Int_t Npy[nch-1]  = segmentation->Npy();

			   Int_t npx=Npx[nch-1];
			   Int_t npy=Npy[nch-1];
			   Float_t thex, they;
			
			   segmentation->GetPadCxy(ipx,ipy,thex,they);
			   Float_t rpad=TMath::Sqrt(thex*thex+they*they);
			   
			   if (rpad < rmin || iqpad ==0 || rpad > rmax) continue;

			
			   trinfo(0)=(Float_t)track;
			   trinfo(1)=(Float_t)iqpad;

			   digits[0]=ipx;
			   digits[1]=ipy;
			   digits[2]=iqpad;

			   // build the list of fired pads and update the info
			   AliMUONlist *pdigit=
			       elem[(nch-1)*(1026*1026)+(ipy+npy)*1026+(ipx+npx)];
			   if (pdigit==0) {
			       obj->AddAtAndExpand(new AliMUONlist(nch-1,digits),counter);
			       counter++;
			       Int_t last=obj->GetLast();
			       pdigit=(AliMUONlist*)obj->At(last);
			       elem[(nch-1)*(1026*1026)+(ipy+npy)*1026+(ipx+npx)]=pdigit;
			       // list of tracks
			       TObjArray *trlist=(TObjArray*)pdigit->TrackList();
			       trlist->Add(trinfo);
			   } else {
			       // update charge
			       (*pdigit).fSignal+=iqpad;
			       // update list of tracks
			       TObjArray* trlist=(TObjArray*)pdigit->TrackList();
			       Int_t last_entry=trlist->GetLast();
			       TVector *ptrk=(TVector*)trlist->At(last_entry);
			       Int_t last_track=Int_t(ptrk(0));
			       Int_t last_charge=Int_t(ptrk(1));
			       if (last_track==track) {
				   last_charge+=iqpad;
				   trlist->RemoveAt(last_entry);
				   trinfo(0)=last_track;
				   trinfo(1)=last_charge;
				   trlist->AddAt(trinfo,last_entry);
			       } else {
				   trlist->Add(trinfo);
			       }
			       // check the track list
			       Int_t nptracks=trlist->GetEntriesFast();
			       if (nptracks > 2) {
				   printf("Attention - nptracks > 2  %d \n",nptracks);
				   printf("cat,nch,ix,iy %d %d %d %d  \n",icat+1,nch,ipx,ipy);
				   for (Int_t tr=0;tr<nptracks;tr++) {
				       TVector *pptrk=(TVector*)trlist->At(tr);
				       trk[tr]=Int_t(pptrk(0));
				       chtrk[tr]=Int_t(pptrk(1));
				   }
			       } // end if nptracks
			   } //  end if pdigit
		       } //end loop over clusters
		   } // hit loop
	       } // if MUON
	   } // track loop

	   Int_t tracks[10];
	   Int_t charges[10];
	   cout<<"start filling digits "<<endl;
	   const Int_t zero_supm = 6;
	   const Int_t adc_satm = 1024;
	   Int_t nd=0;


	   Int_t nentries=obj->GetEntriesFast();
	   printf(" nentries %d \n",nentries);

	   // start filling the digits
	   
	   for (Int_t nent=0;nent<nentries;nent++) {
	       AliMUONlist *address=(AliMUONlist*)obj->At(nent);
	       if (address==0) continue; 
	       // do zero-suppression and signal truncation
	       Int_t ich=address->fRpad;
	       Int_t q=address->fSignal;        
	       // if ( q <= zero_supm  || rpad > 55) continue;
	       if ( q <= zero_supm ) continue;
	       if ( q > adc_satm)  q=adc_satm;
	       digits[0]=address->fPadX;
	       digits[1]=address->fPadY;
	       digits[2]=q;

	       TObjArray* trlist=(TObjArray*)address->TrackList();
	       Int_t nptracks=trlist->GetEntriesFast();
	       // this was changed to accomodate the real number of tracks
	       if (nptracks > 10) {
		   cout<<"Attention - nptracks > 3 "<<nptracks<<endl;
		   nptracks=10;
	       }
               if (nptracks > 2) {
		   printf("Attention - nptracks > 2  %d \n",nptracks);
		   printf("cat,ich,ix,iy,q %d %d %d %d %d \n",icat,ich,digits[0],digits[1],q);
               }
	       for (Int_t tr=0;tr<nptracks;tr++) {
		   TVector *pp=(TVector*)trlist->At(tr);
		   Int_t entries=pp->GetNrows();
		   tracks[tr]=Int_t(pp(0));
		   charges[tr]=Int_t(pp(1));
	       }      //end loop over list of tracks for one pad
	       
	       // fill digits
	       MUON->AddDigits(ich,tracks,charges,digits);
	   }
	   cout<<"I'm out of the loops for digitisation"<<endl;
	   gAlice->TreeD()->Fill();
	   TTree *TD=gAlice->TreeD();
	   int ndig=TD->GetEntries();
	   cout<<"number of digits  "<<ndig<<endl;
	   TClonesArray *fDch;
	   for (int i=0;i<10;i++) {
	       fDch= MUON->DigitsAddress(i);
	       int ndig=fDch->GetEntriesFast();
	       printf (" i, ndig %d %d \n",i,ndig);
	   }
	   MUON->ResetDigits();
	   obj->Clear();
       } //end loop over cathodes
       char hname[30];
       sprintf(hname,"TreeD%d",nev);
       gAlice->TreeD()->Write(hname);
       file->ls();
   } // event loop 
   file->Close();
//   cout<<"END  digitisation "<<endl;
}
