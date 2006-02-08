/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
// Macro for checking aliroot output and associated files contents
// Gines Martinez, Subatech June 2003
//
#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include "TBranch.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TParticle.h"
#include "TTree.h"
#include "TNtuple.h"

// STEER includes
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliLoader.h"
#include "AliStack.h"

// MUON includes
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONHit.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"

#include "AliMpVSegmentation.h"
#include "AliMpIntPair.h"
#include "AliMpDEManager.h"
#include "AliMpSegFactory.h"
#endif


void MUONkine(Int_t event2Check=0, char * filename="galice.root")
{
  // Stack of particle for each event
  AliStack* stack;
  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }

  RunLoader->LoadKinematics("READ");
  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();

  for(ievent=0; ievent<nevents; ievent++) {  // Event loop
    if (event2Check!=0) ievent=event2Check;
    Int_t iparticle, nparticles;
    // Getting event ievent
    RunLoader->GetEvent(ievent); 
    stack = RunLoader->Stack();
    nparticles = (Int_t) stack->GetNtrack();
    printf(">>> Event %d, Number of particles is %d \n",ievent, nparticles);
    for(iparticle=0; iparticle<nparticles; iparticle++) {
      stack->Particle(iparticle)->Print("");  
    }
    if (event2Check!=0) ievent=nevents;
  }
  RunLoader->UnloadKinematics();
}


void MUONhits(Int_t event2Check=0, char * filename="galice.root")
{
  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  // Loading MUON subsystem
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadHits("READ");  // Loading Tree of hits for MUON
  AliMUONData muondata(MUONLoader,"MUON","MUON");  // Creating MUON data container
  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();

  for(ievent=0; ievent<nevents; ievent++) {  // Event loop
    if (event2Check!=0) ievent=event2Check;
    printf(">>> Event %d \n",ievent);
    // Getting event ievent
    RunLoader->GetEvent(ievent); 
    muondata.SetTreeAddress("H");
    Int_t itrack, ntracks;
    ntracks = (Int_t) muondata.GetNtracks();
    for (itrack=0; itrack<ntracks; itrack++) { // Track loop
      //Getting List of Hits of Track itrack
      muondata.GetTrack(itrack);
      Int_t ihit, nhits;
      nhits = (Int_t) muondata.Hits()->GetEntriesFast();
      printf(">>> Track %d, Number of hits %d \n",itrack,nhits);
      AliMUONHit* mHit;
      for(ihit=0; ihit<nhits; ihit++) {
	mHit = static_cast<AliMUONHit*>(muondata.Hits()->At(ihit));
	Int_t detele   = mHit-> DetElemId(); // Detection element if defined
	Int_t hittrack = mHit->Track();
	Float_t x      = mHit->X();
  	Float_t y      = mHit->Y();
  	Float_t z      = mHit->Z();
  	Float_t elos   = mHit->Eloss();
  	Float_t momentum = mHit->Momentum();
  	printf(">>> >>>  Hit%4d DetEle %4d Track%4d (X,Y,Z)=(%7.2f,%7.2f,%8.2f)cm Elost=%7.2gGeV  P=%6.1fGeV/c\n",
	       ihit, detele, hittrack,x,y,z,elos,momentum);
      }
      muondata.ResetHits();
    } // end track loop
    if (event2Check!=0) ievent=nevents;
  }  // end event loop
  MUONLoader->UnloadHits();
}


void MUONdigits(Int_t event2Check=0, char * filename="galice.root")
{
  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  // Loading MUON subsystem
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadDigits("READ");
  // Creating MUON data container
  AliMUONData muondata(MUONLoader,"MUON","MUON");

  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();
  AliMUONDigit * mDigit;
  
  for(ievent=0; ievent<nevents; ievent++) {
    if (event2Check!=0) ievent=event2Check;
    printf(">>> Event %d \n",ievent);
    RunLoader->GetEvent(ievent);
  
    // Addressing
    Int_t ichamber, nchambers;
    nchambers = AliMUONConstants::NCh(); ;
    muondata.SetTreeAddress("D,GLT");
    //    char branchname[30];    
 
    muondata.GetDigits();
    // Loop on chambers
    for( ichamber=0; ichamber<nchambers; ichamber++) {
      Int_t idigit, ndigits;
      ndigits = (Int_t) muondata.Digits(ichamber)->GetEntriesFast();
      for(idigit=0; idigit<ndigits; idigit++) {
	mDigit = static_cast<AliMUONDigit*>(muondata.Digits(ichamber)->At(idigit));
	Int_t PadX   = mDigit->PadX();     // Pad X number
	Int_t PadY   = mDigit->PadY();     // Pad Y number
	Int_t Signal = mDigit->Signal();   // Physics Signal
	Int_t Physics= mDigit->Physics();  // Physics contribution to signal
	//	  Int_t Hit    = mDigit->Hit();      // iHit
	Int_t Cathode= mDigit->Cathode();  // Cathode
	Int_t Track0 = mDigit->Track(0);
	Int_t Track1 = mDigit->Track(1); 
	//Int_t Track2 = mDigit->Track(2);
	Int_t TCharges0 = mDigit->TrackCharge(0);  //charge per track making this digit (up to 10)
	Int_t TCharges1 = mDigit->TrackCharge(1);
	//Int_t TCharges2 = mDigit->TrackCharge(2);
	Int_t idDE = mDigit->DetElemId();
	//	  printf(">>> Cathode %d\n",Cathode);
	
	printf(">>> DetEle %4d Digit%4d Cath %1d (Ix,Iy)=(%3d,%3d) Signal=%4d Physics=%4d Track0=%4d Charge0=%4d Track1=%4d Charge1=%4d \n",
	       idDE, idigit, Cathode, PadX, PadY, Signal, Physics, Track0, 
	       TCharges0, Track1, TCharges1);
      } // end digit loop
    } // end chamber loop
    muondata.ResetDigits();
    //    } // end cathode loop
    if (event2Check!=0) ievent=nevents;
  }  // end event loop
  MUONLoader->UnloadDigits();
}

void MUONoccupancy(Int_t event2Check=0,  Bool_t perDetEle =kFALSE, char * filename="galice.root") {
  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  // Loading MUON subsystem
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadDigits("READ");
  // Creating MUON data container
  AliMUONData muondata(MUONLoader,"MUON","MUON");
  
  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();
  AliMUONDigit * mDigit =0x0;
  AliMpVSegmentation * segbend = 0x0;
  AliMpVSegmentation * segnonbend = 0x0;
  AliMpIntPair pad(0,0);

  Int_t dEoccupancy_bending[14][26];
  Int_t dEoccupancy_nonbending[14][26];
  Int_t cHoccupancy_bending[14];
  Int_t cHoccupancy_nonbending[14];
  Int_t totaloccupancy_bending =0;
  Int_t totaloccupancy_nonbending =0;

  Int_t dEchannels_bending[14][26];
  Int_t dEchannels_nonbending[14][26];
  Int_t cHchannels_bending[14];
  Int_t cHchannels_nonbending[14];
  Int_t totalchannels_bending =0;
  Int_t totalchannels_nonbending =0;

  Int_t ichamber, nchambers,idetele, detele, ix, iy;
  nchambers = AliMUONConstants::NCh(); ;

  AliMpSegFactory factory;

  for (ichamber=0; ichamber<nchambers; ichamber++) {
    cHchannels_bending[ichamber]=0;
    cHchannels_nonbending[ichamber]=0;
    for (idetele=0; idetele<26; idetele++) {
      detele= 100*(ichamber +1)+idetele;
      dEchannels_bending[ichamber][idetele]=0;
      dEchannels_nonbending[ichamber][idetele]=0;
      dEoccupancy_bending[ichamber][idetele]=0;
      dEoccupancy_nonbending[ichamber][idetele]=0;
      if ( AliMpDEManager::IsValidDetElemId(detele) ) {
	
	segbend    =  factory.CreateMpSegmentation(detele, 0);
	segnonbend =  factory.CreateMpSegmentation(detele, 1);
        if (AliMpDEManager::GetPlaneType(detele, 0) != kBendingPlane ) {
	  AliMpVSegmentation* tmp = segbend;
	  segbend    =  segnonbend;
	  segnonbend =  tmp;
	}  
	  
	for(ix=0; ix<=segbend->MaxPadIndexX(); ix++) {
	  for(iy=0; iy<=segbend->MaxPadIndexY(); iy++) {
	    pad.SetFirst(ix);
	    pad.SetSecond(iy);
	    if( segbend->HasPad(pad) )   {  
	      dEchannels_bending[ichamber][idetele]++;
	      cHchannels_bending[ichamber]++;
	      totalchannels_bending++;
	    }
	  }
	}
	for(ix=0; ix<=segnonbend->MaxPadIndexX(); ix++) {
	  for(iy=0; iy<=segnonbend->MaxPadIndexY(); iy++) {
	    pad.SetFirst(ix);
	    pad.SetSecond(iy);
	    if(segnonbend->HasPad(pad))  {
	      dEchannels_nonbending[ichamber][idetele]++;  
	      cHchannels_nonbending[ichamber]++;
	      totalchannels_nonbending++;
	    }
	  }
	}
	if (perDetEle) printf(">>> Detection element %4d has %5d channels in bending and %5d channels in nonbending \n",
	     detele, dEchannels_bending[ichamber][idetele], dEchannels_nonbending[ichamber][idetele] ); 
      }
    }
    printf(">>> Chamber %2d has %6d channels in bending and %6d channels in nonbending \n",
	   ichamber+1,  cHchannels_bending[ichamber], cHchannels_nonbending[ichamber]);
  }
  printf(">>Spectrometer has  %7d channels in bending and %7d channels in nonbending \n",
	 totalchannels_bending, totalchannels_nonbending);

  factory.DeleteSegmentations();

  ievent=event2Check;
  printf(">>> Event %d \n",ievent);
  RunLoader->GetEvent(ievent);
    
  // Addressing
  muondata.SetTreeAddress("D"); 
  muondata.GetDigits();
  // Loop on chambers
  for( ichamber=0; ichamber<nchambers; ichamber++) {
    cHoccupancy_bending[ichamber]   = 0;
    cHoccupancy_nonbending[ichamber]= 0;
    Int_t idigit, ndigits;
    ndigits = (Int_t) muondata.Digits(ichamber)->GetEntriesFast();
    for(idigit=0; idigit<ndigits; idigit++) {
      mDigit = static_cast<AliMUONDigit*>(muondata.Digits(ichamber)->At(idigit));
      Int_t detele = mDigit->DetElemId();
      Int_t idetele = detele-(ichamber+1)*100;
      if ( mDigit->Cathode() == 0 ) {

	cHoccupancy_bending[ichamber]++;
	dEoccupancy_bending[ichamber][idetele]++;
	totaloccupancy_bending++;
      }
      else {
	cHoccupancy_nonbending[ichamber]++;
	dEoccupancy_nonbending[ichamber][idetele]++;
	totaloccupancy_nonbending++;
      }
    } // end digit loop    

    printf(">>> Chamber %2d  nChannels Bending %5d  nChannels NonBending %5d \n", 
	   ichamber+1, 
	   cHoccupancy_bending[ichamber],
	   cHoccupancy_nonbending[ichamber]);           
    printf(">>> Chamber %2d  Occupancy Bending %5.2f %%  Occupancy NonBending %5.2f %% \n", 
	   ichamber+1, 
	   100.*((Float_t) cHoccupancy_bending[ichamber])/((Float_t) cHchannels_bending[ichamber]),
	   100.*((Float_t) cHoccupancy_nonbending[ichamber])/((Float_t) cHchannels_bending[ichamber])            );


    for(Int_t idetele=0; idetele<26; idetele++) {
      Int_t detele = idetele + 100*(ichamber+1);
      if ( AliMpDEManager::IsValidDetElemId(detele) ) {
	if (perDetEle) {
	  printf(">>> DetEle %4d nChannels Bending %5d  nChannels NonBending %5d \n", 
		 idetele+100*(ichamber+1), 
		 dEoccupancy_bending[ichamber][idetele],
		 dEoccupancy_nonbending[ichamber][idetele]);  
	  printf(">>> DetEle %4d Occupancy Bending %5.2f %%  Occupancy NonBending %5.2f %% \n", 
		 idetele+100*(ichamber+1), 
		 100.*((Float_t) dEoccupancy_bending[ichamber][idetele])/((Float_t) dEchannels_bending[ichamber][idetele]),
		 100.*((Float_t) dEoccupancy_nonbending[ichamber][idetele])/((Float_t) dEchannels_bending[ichamber][idetele]));  
	}
      }
    }
  } // end chamber loop
  printf(">>> Muon Spectrometer  Occupancy Bending %5.2f %%  Occupancy NonBending %5.2f %% \n",  
	   100.*((Float_t) totaloccupancy_bending)/((Float_t) totalchannels_bending),
	 100.*((Float_t) totaloccupancy_nonbending)/((Float_t) totalchannels_nonbending)            );
  muondata.ResetDigits();
  //    } // end cathode loop
  MUONLoader->UnloadDigits();
}

void MUONrecpoints(Int_t event2Check=0, char * filename="galice.root") {

  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  // Getting MUONloader
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadRecPoints("READ");
  // Creating MUON data container
  AliMUONData muondata(MUONLoader,"MUON","MUON");

  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();
  AliMUONRawCluster * mRecPoint = 0;
  
  for(ievent=0; ievent<nevents; ievent++) {
    if (event2Check!=0) ievent=event2Check;
    printf(">>> Event %d \n",ievent);
    RunLoader->GetEvent(ievent);
    // Addressing
    Int_t ichamber, nchambers;
    nchambers = AliMUONConstants::NTrackingCh();
    muondata.SetTreeAddress("RC,TC"); 
    char branchname[30];    
    muondata.GetRawClusters();
    // Loop on chambers
    for( ichamber=0; ichamber<nchambers; ichamber++) {
      sprintf(branchname,"MUONRawClusters%d",ichamber+1);
      //printf(">>>  branchname %s\n",branchname);
      Int_t irecpoint, nrecpoints;
      nrecpoints = (Int_t) muondata.RawClusters(ichamber)->GetEntriesFast();
      // printf(">>> Chamber %2d, Number of recpoints = %6d \n",ichamber+1, nrecpoints);
      for(irecpoint=0; irecpoint<nrecpoints; irecpoint++) {
	mRecPoint = static_cast<AliMUONRawCluster*>(muondata.RawClusters(ichamber)->At(irecpoint));
	Int_t Track0 = mRecPoint->GetTrack(0);
	Int_t Track1 = mRecPoint->GetTrack(1); 
	Int_t Track2 = mRecPoint->GetTrack(2);
	Int_t Q0 = mRecPoint->GetCharge(0);
	Int_t Q1 = mRecPoint->GetCharge(1);
	Float_t x0 = mRecPoint->GetX(0);
	Float_t x1 = mRecPoint->GetX(1);
	Float_t y0 = mRecPoint->GetY(0);
	Float_t y1 = mRecPoint->GetY(1);
	Float_t z0 = mRecPoint->GetZ(0);
	Float_t z1 = mRecPoint->GetZ(1);
	Float_t chi2_0 =  mRecPoint->GetChi2(0);
	//Float_t chi2_1 =  mRecPoint->GetChi2(1);
	Int_t de = mRecPoint->GetDetElemId();
	printf(">>> >>> RecPoint %4d  DetEle %4d (X,Y,Z)=(%7.2f,%7.2f,%8.2f)cm  Q0=%4d  Q1=%4d Hit=%4d Track1=%4d Track2=%4d Chi2=%6.3f \n",
	       irecpoint,de,x0,y0,z0,Q0,Q1,Track0, Track1, Track2, chi2_0);
	if( (x0!=x1) || (y0!=y1) || (z0!=z1) )
	  printf(">>> >>> Warning (X0,Y0,Z0)=(%7.2f, %7.2f, %8.2f)cm != (X1,Y1,Z1)=(%7.2f,%7.2f,%8.2f)cm \n",x0,y0,z0,x1,y1,z1); 
      } // end recpoint loop
    } // end chamber loop
    muondata.ResetRawClusters();
    if (event2Check!=0) ievent=nevents;
  }  // end event loop
  MUONLoader->UnloadRecPoints();
}



void MUONrectrigger (Int_t event2Check=0, char * filename="galice.root"){
  // reads and dumps trigger objects from MUON.RecPoints.root
  TClonesArray * globalTrigger;
  TClonesArray * localTrigger;
  
  // Do NOT print out all the info if the loop runs over all events 
  Int_t PRINTOUT = (event2Check == 0 ) ? 0 : 1 ;  
  
  // Book a ntuple for more detailled studies
  TNtuple *Tgtuple = new TNtuple("Tgtuple","Trigger Ntuple","ev:global:spapt:smapt:undefapt:uplpt:uphpt:upapt:lpapt");
  Int_t WRITE = 0;

  // counters
  Int_t SPLowpt=0,SPHighpt=0,SPAllpt=0;
  Int_t SMLowpt=0,SMHighpt=0,SMAllpt=0;
  Int_t SULowpt=0,SUHighpt=0,SUAllpt=0;
  Int_t USLowpt=0,USHighpt=0,USAllpt=0;
  Int_t LSLowpt=0,LSHighpt=0,LSAllpt=0;

  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadDigits("READ");
  // Creating MUON data container
  AliMUONData muondata(MUONLoader,"MUON","MUON");
  
  
  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();
  
  AliMUONGlobalTrigger *gloTrg(0x0);
  AliMUONLocalTrigger *locTrg(0x0);
  
  for (ievent=0; ievent<nevents; ievent++) {
    if (event2Check!=0) ievent=event2Check;
    if (ievent%100==0 || event2Check) cout << "Processing event " << ievent << endl;
    RunLoader->GetEvent(ievent);
    
    muondata.SetTreeAddress("GLT"); 
    muondata.GetTriggerD();
    
    globalTrigger = muondata.GlobalTrigger();
    localTrigger = muondata.LocalTrigger();
    
    Int_t nglobals = (Int_t) globalTrigger->GetEntriesFast(); // should be 1
    Int_t nlocals  = (Int_t) localTrigger->GetEntriesFast(); // up to 234
    if (PRINTOUT) printf("###################################################\n");
    if (PRINTOUT) {cout << " event " << ievent 
                        << " nglobals nlocals: " << nglobals << " " << nlocals << "\n"; }
    
    for (Int_t iglobal=0; iglobal<nglobals; iglobal++) { // Global Trigger
      gloTrg = static_cast<AliMUONGlobalTrigger*>(globalTrigger->At(iglobal));
      
      SPLowpt+=gloTrg->SinglePlusLpt() ;
      SPHighpt+=gloTrg->SinglePlusHpt() ;
      SPAllpt+=gloTrg->SinglePlusApt() ;
      SMLowpt+=gloTrg->SingleMinusLpt();
      SMHighpt+=gloTrg->SingleMinusHpt();
      SMAllpt+=gloTrg->SingleMinusApt();
      SULowpt+=gloTrg->SingleUndefLpt();
      SUHighpt+=gloTrg->SingleUndefHpt();
      SUAllpt+=gloTrg->SingleUndefApt();
      USLowpt+=gloTrg->PairUnlikeLpt(); 
      USHighpt+=gloTrg->PairUnlikeHpt();
      USAllpt+=gloTrg->PairUnlikeApt();
      LSLowpt+=gloTrg->PairLikeLpt(); 
      LSHighpt+=gloTrg->PairLikeHpt();
      LSAllpt+=gloTrg->PairLikeApt();

      if (PRINTOUT) {
        printf("===================================================\n");
        printf(" Global Trigger output       Low pt  High pt   All\n");
        printf(" number of Single Plus      :\t");
        printf("%i\t%i\t%i\t",gloTrg->SinglePlusLpt(),
               gloTrg->SinglePlusHpt(),gloTrg->SinglePlusApt());
        printf("\n");

        printf(" number of Single Minus     :\t");
        printf("%i\t%i\t%i\t",gloTrg->SingleMinusLpt(),
               gloTrg->SingleMinusHpt(),gloTrg->SingleMinusApt());
        printf("\n");

        printf(" number of Single Undefined :\t"); 
        printf("%i\t%i\t%i\t",gloTrg->SingleUndefLpt(),
               gloTrg->SingleUndefHpt(),gloTrg->SingleUndefApt());
        printf("\n");
        
        printf(" number of UnlikeSign pair  :\t"); 
        printf("%i\t%i\t%i\t",gloTrg->PairUnlikeLpt(),
               gloTrg->PairUnlikeHpt(),gloTrg->PairUnlikeApt());
        printf("\n");
        
        printf(" number of LikeSign pair    :\t");  
        printf("%i\t%i\t%i\t",gloTrg->PairLikeLpt(),
               gloTrg->PairLikeHpt(),gloTrg->PairLikeApt());
        printf("\n");
        
        printf("===================================================\n");
      }
      
    } // end of loop on Global Trigger
    
    for (Int_t ilocal=0; ilocal<nlocals; ilocal++) { // Local Trigger
      if (PRINTOUT) cout << " >>> Output for Local Trigger " << ilocal << "\n";
      
      locTrg = static_cast<AliMUONLocalTrigger*>(localTrigger->At(ilocal));
      
      if (PRINTOUT){
        cout << "Circuit StripX Dev StripY: " 
             << locTrg->LoCircuit() << " " 
             << locTrg->LoStripX() << " " 
             << locTrg->LoDev() << " " 
             << locTrg->LoStripY() 
             << "\n";
        cout << "Lpt Hpt Apt: "     
             << locTrg->LoLpt() << " "   
             << locTrg->LoHpt() << " "  
             << locTrg->LoApt() << "\n";
      }
    } // end of loop on Local Trigger


    // fill ntuple
    //TNtuple *Tgtuple = new TNtuple("Tgtuple","Trigger Ntuple","ev:global:spapt:smapt:undefapt:uplpt:uphpt:upapt:lpapt");
       Tgtuple->Fill(ievent,nglobals,gloTrg->SinglePlusApt(),gloTrg->SingleMinusApt(),gloTrg->SingleUndefApt(),gloTrg->PairUnlikeLpt(),gloTrg->PairUnlikeHpt(),gloTrg->PairUnlikeApt(),gloTrg->PairLikeApt());


    muondata.ResetTrigger();
    if (event2Check!=0) ievent=nevents;
  } // end loop on event  
  
  // Print out summary if loop ran over all event
  if (!event2Check){
    printf("\n");
    printf("===================================================\n");
    printf("===================  SUMMARY  =====================\n");
    printf("\n");
    printf("Total number of events processed %d \n", (event2Check==0) ? nevents : 1);
    printf("\n");
    printf(" Global Trigger output       Low pt  High pt   All\n");
    printf(" number of Single Plus      :\t");
    printf("%i\t%i\t%i\t",SPLowpt,SPHighpt,SPAllpt);
    printf("\n");
    printf(" number of Single Minus     :\t");
    printf("%i\t%i\t%i\t",SMLowpt,SMHighpt,SMAllpt);
    printf("\n");
    printf(" number of Single Undefined :\t"); 
    printf("%i\t%i\t%i\t",SULowpt,SUHighpt,SUAllpt);
    printf("\n");
    printf(" number of UnlikeSign pair  :\t"); 
    printf("%i\t%i\t%i\t",USLowpt,USHighpt,USAllpt);
    printf("\n");
    printf(" number of LikeSign pair    :\t");  
    printf("%i\t%i\t%i\t",LSLowpt,LSHighpt, LSAllpt);
    printf("\n");
    printf("===================================================\n");
  }
  
  if (WRITE){
    TFile *myFile = new TFile("TriggerCheck.root", "RECREATE");
    Tgtuple->Write();
    myFile->Close();
  }


  MUONLoader->UnloadRecPoints();
}

void MUONrectracks (Int_t event2Check=0, char * filename="galice.root"){
// reads and dumps trigger objects from MUON.RecPoints.root
  TClonesArray * RecTracks;
  
  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadTracks("READ");
  // Creating MUON data container
  AliMUONData muondata(MUONLoader,"MUON","MUON");
  
    Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();
  
  //  AliMUONTrack * rectrack;
  
  for (ievent=0; ievent<nevents; ievent++) {
    if (event2Check!=0) ievent=event2Check;
    RunLoader->GetEvent(ievent);
    
    muondata.SetTreeAddress("RT");
    muondata.GetRecTracks();
    RecTracks = muondata.RecTracks();
    
    
    Int_t nrectracks = (Int_t) RecTracks->GetEntriesFast(); //

    printf(">>> Event %d, Number of Recconstructed tracks %d \n",ievent, nrectracks);
    // loop over tracks
 
 
    Int_t nTrackHits;// nPrimary;
    Double_t fitFmin;
    Double_t bendingSlope, nonBendingSlope, inverseBendingMomentum;
    Double_t xRec, yRec, zRec, chi2MatchTrigger;
    Bool_t matchTrigger;
    Double_t Pz,Px,Py,Pt,Ptot,Eta ;

  // setting pointer for tracks, triggertracks & trackparam at vertex
    AliMUONTrack* recTrack = 0;
    AliMUONTrackParam* trackParam = 0;

    for (Int_t iRecTracks = 0; iRecTracks <  nrectracks;  iRecTracks++) {
    // reading info from tracks
      recTrack = (AliMUONTrack*) RecTracks->At(iRecTracks);
      trackParam = (AliMUONTrackParam*) (recTrack->GetTrackParamAtHit())->First();
      trackParam->ExtrapToZ(0.0);
      bendingSlope            = trackParam->GetBendingSlope();
      nonBendingSlope         = trackParam->GetNonBendingSlope();
      inverseBendingMomentum = trackParam->GetInverseBendingMomentum();
      xRec  = trackParam->GetNonBendingCoor();
      yRec  = trackParam->GetBendingCoor();
      zRec  = trackParam->GetZ();

      nTrackHits       = recTrack->GetNTrackHits();
      fitFmin          = recTrack->GetFitFMin();
      matchTrigger     = recTrack->GetMatchTrigger();
      chi2MatchTrigger = recTrack->GetChi2MatchTrigger();
      
      Px = trackParam->Px();
      Py = trackParam->Py(); 
      Pz = trackParam->Pz(); 
      Pt = TMath::Sqrt(Px*Px + Py*Py );
      Ptot = TMath::Sqrt(Px*Px + Py*Py + Pz*Pz);
      Eta =  (Pt!=0) ? 0.5*log( (Ptot+Pz)/(Ptot-Pz) ) : 999999999.999 ;
       
      printf(">>> RecTrack %4d  NofClusters=%2d BendMomentum=%7.2f NonBendSlope=%5.2f  BendSlope=%5.2f Match2Trig=%1d (vertex@z=0)=(%5.2f,%5.2f,%5.1f)cm \n", iRecTracks, nTrackHits, 1/inverseBendingMomentum , nonBendingSlope*180./TMath::Pi(), bendingSlope*180./TMath::Pi(),  matchTrigger, xRec,yRec,zRec);
      printf("    Px=%f  Py =%f  Pz =%f   Pt=%f  Ptot=%f   PseudoRap=%f  \n",Px,Py,Pz,Pt,Ptot,Eta);
    } // end loop tracks

    muondata.ResetRecTracks();
    if (event2Check!=0) ievent=nevents;
  } // end loop on event  
  MUONLoader->UnloadTracks();
}










