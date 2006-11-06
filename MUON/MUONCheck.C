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
#include "AliTracker.h"
#include "AliMagFMaps.h"

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
#include "AliMUONTrackExtrap.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONTriggerCrateStore.h"

#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpIntPair.h"
#include "AliMpDEManager.h"
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
	mHit->Print("full");
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
    
    muondata.GetDigits();
    // Loop on chambers
    for( ichamber=0; ichamber<nchambers; ichamber++) {
      Int_t idigit, ndigits;
      TClonesArray* digits = muondata.Digits(ichamber);
      digits->Sort();
      ndigits = (Int_t)digits->GetEntriesFast();
      for(idigit=0; idigit<ndigits; idigit++) {
        mDigit = static_cast<AliMUONDigit*>(digits->At(idigit));
        mDigit->Print("tracks");
      } // end digit loop
    } // end chamber loop
    muondata.ResetDigits();
    if (event2Check!=0) ievent=nevents;
  }  // end event loop
  MUONLoader->UnloadDigits();
}

void MUONsdigits(Int_t event2Check=0, char * filename="galice.root")
{
  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  // Loading MUON subsystem
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadSDigits("READ");
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
    muondata.SetTreeAddress("S");
    
    muondata.GetSDigits();
    // Loop on chambers
    for( ichamber=0; ichamber<nchambers; ichamber++) {
      Int_t idigit, ndigits;
      TClonesArray* digits = muondata.SDigits(ichamber);
      ndigits = (Int_t)digits->GetEntriesFast();
      for(idigit=0; idigit<ndigits; idigit++) {
        mDigit = static_cast<AliMUONDigit*>(digits->At(idigit));
        mDigit->Print("tracks");
      } // end digit loop
    } // end chamber loop
    muondata.ResetSDigits();
    if (event2Check!=0) ievent=nevents;
  }  // end event loop
  MUONLoader->UnloadSDigits();

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
  const AliMpVSegmentation * segbend = 0x0;
  const AliMpVSegmentation * segnonbend = 0x0;
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
	
	segbend    = AliMpSegmentation::Instance()->GetMpSegmentation(detele, 0);
	segnonbend = AliMpSegmentation::Instance()->GetMpSegmentation(detele, 1);
        if (AliMpDEManager::GetPlaneType(detele, 0) != kBendingPlane ) {
	  const AliMpVSegmentation* tmp = segbend;
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
	mRecPoint->Print("full");
      } // end recpoint loop
    } // end chamber loop
    muondata.ResetRawClusters();
    if (event2Check!=0) ievent=nevents;
  }  // end event loop
  MUONLoader->UnloadRecPoints();
}

void MUONrectrigger (Int_t event2Check=0, char * filename="galice.root", Int_t WRITE = 0, Bool_t readFromRP = 0)
{

  // reads and dumps trigger objects from MUON.RecPoints.root
  TClonesArray * globalTrigger;
  TClonesArray * localTrigger;
  
  // Do NOT print out all the info if the loop runs over all events 
  Int_t PRINTOUT = (event2Check == 0 ) ? 0 : 1 ;  

  // Book a ntuple for more detailled studies
  TNtuple *TgtupleGlo = new TNtuple("TgtupleGlo","Global Trigger Ntuple","ev:global:slpt:shpt:uplpt:uphpt:lplpt:lplpt");
  TNtuple *TgtupleLoc = new TNtuple("TgtupleLoc","Local Trigger Ntuple","ev:LoCircuit:LoStripX:LoDev:StripY:LoLpt:LoHpt:y11:y21:x11");

  // counters
  Int_t SLowpt=0,SHighpt=0;
  Int_t USLowpt=0,USHighpt=0;
  Int_t LSLowpt=0,LSHighpt=0;

  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }

  AliMUONTriggerCrateStore* crateManager = new AliMUONTriggerCrateStore();   
  crateManager->ReadFromFile();

  AliMUONGeometryTransformer* transformer = new AliMUONGeometryTransformer(kFALSE);
  transformer->ReadGeometryData("volpath.dat", "geometry.root");

  TClonesArray*  triggerCircuit = new TClonesArray("AliMUONTriggerCircuit", 234);

  for (Int_t i = 0; i < AliMUONConstants::NTriggerCircuit(); i++)  {
      AliMUONTriggerCircuit* c = new AliMUONTriggerCircuit();
      c->SetTransformer(transformer);
      c->Init(i,*crateManager);
      TClonesArray& circuit = *triggerCircuit;
      new(circuit[circuit.GetEntriesFast()])AliMUONTriggerCircuit(*c);
      delete c;
  }
  
  Char_t fileName[30];
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  if (!readFromRP) {
      cout << " reading from digits \n";
      MUONLoader->LoadDigits("READ");
      sprintf(fileName,"TriggerCheckFromDigits.root");
  } else {
      cout << " reading from RecPoints \n";
      MUONLoader->LoadRecPoints("READ");
      sprintf(fileName,"TriggerCheckFromRP.root");
  }

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
    
    if (!readFromRP) {
	muondata.SetTreeAddress("D,GLT"); 
	muondata.GetTriggerD();
    } else {    
	muondata.SetTreeAddress("RC,TC"); 
	muondata.GetTrigger();
    }

    globalTrigger = muondata.GlobalTrigger();
    localTrigger = muondata.LocalTrigger();
    
    Int_t nglobals = (Int_t) globalTrigger->GetEntriesFast(); // should be 1
    Int_t nlocals  = (Int_t) localTrigger->GetEntriesFast(); // up to 234
    if (PRINTOUT) printf("###################################################\n");
    if (PRINTOUT) printf("event %d nglobal %d nlocal %d \n",ievent,nglobals,nlocals);

    for (Int_t iglobal=0; iglobal<nglobals; iglobal++) { // Global Trigger
      gloTrg = static_cast<AliMUONGlobalTrigger*>(globalTrigger->At(iglobal));
      
      SLowpt+=gloTrg->SingleLpt() ;
      SHighpt+=gloTrg->SingleHpt() ;
      USLowpt+=gloTrg->PairUnlikeLpt(); 
      USHighpt+=gloTrg->PairUnlikeHpt();
      LSLowpt+=gloTrg->PairLikeLpt(); 
      LSHighpt+=gloTrg->PairLikeHpt();
      
      if (PRINTOUT) gloTrg->Print("full");

    } // end of loop on Global Trigger

    for (Int_t ilocal=0; ilocal<nlocals; ilocal++) { // Local Trigger
      locTrg = static_cast<AliMUONLocalTrigger*>(localTrigger->At(ilocal));
      if (PRINTOUT) locTrg->Print("full");
      
      AliMUONTriggerCircuit* circuit = (AliMUONTriggerCircuit*)triggerCircuit->At(locTrg->LoCircuit()-1); 
      
      TgtupleLoc->Fill(ievent,locTrg->LoCircuit(),locTrg->LoStripX(),locTrg->LoDev(),locTrg->LoStripY(),locTrg->LoLpt(),locTrg->LoHpt(),circuit->GetY11Pos(locTrg->LoStripX()),circuit->GetY21Pos(locTrg->LoStripX()+locTrg->LoDev()+1),circuit->GetX11Pos(locTrg->LoStripY()));
    } // end of loop on Local Trigger

    // fill ntuple
    TgtupleGlo->Fill(ievent,nglobals,gloTrg->SingleLpt(),gloTrg->SingleHpt(),gloTrg->PairUnlikeLpt(),gloTrg->PairUnlikeHpt(),gloTrg->PairLikeLpt(),gloTrg->PairLikeHpt());

    muondata.ResetTrigger();
    if (event2Check!=0) ievent=nevents;
  } // end loop on event  
  
  // Print out summary if loop ran over all event
  if (!event2Check){

    printf("\n");
    printf("=============================================\n");
    printf("================  SUMMARY  ==================\n");
    printf("\n");
    printf("Total number of events processed %d \n", (event2Check==0) ? nevents : 1);
    printf("\n");
    printf(" Global Trigger output       Low pt  High pt\n");
    printf(" number of Single           :\t");
    printf("%i\t%i\t",SLowpt,SHighpt);
    printf("\n");
    printf(" number of UnlikeSign pair  :\t"); 
    printf("%i\t%i\t",USLowpt,USHighpt);
    printf("\n");
    printf(" number of LikeSign pair    :\t");  
    printf("%i\t%i\t",LSLowpt,LSHighpt);
    printf("\n");
    printf("=============================================\n");
    fflush(stdout);
  }
  
  if (WRITE){
      TFile *myFile = new TFile(fileName, "RECREATE");
      TgtupleGlo->Write();
      TgtupleLoc->Write();
      myFile->Close();
  }

  MUONLoader->UnloadRecPoints();

  delete crateManager;
  delete transformer;
  delete triggerCircuit;
  
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
    // waiting for mag field in CDB 
  printf("Loading field map...\n");
  if (!AliTracker::GetFieldMap()) {
    AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
    AliTracker::SetFieldMap(field, kFALSE);
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

    // setting pointer for tracks, triggertracks & trackparam at vertex
    AliMUONTrack* recTrack = 0;
    AliMUONTrackParam* trackParam = 0;

    // set the magnetic field for track extrapolations
    AliMUONTrackExtrap::SetField(AliTracker::GetFieldMap());
    for (Int_t iRecTracks = 0; iRecTracks <  nrectracks;  iRecTracks++) {
   //  // reading info from tracks
       recTrack = (AliMUONTrack*) RecTracks->At(iRecTracks);
       trackParam = (AliMUONTrackParam*) (recTrack->GetTrackParamAtHit())->First();
       AliMUONTrackExtrap::ExtrapToZ(trackParam,0.);
      recTrack->Print("full");
    } // end loop tracks

    muondata.ResetRecTracks();
    if (event2Check!=0) ievent=nevents;
  } // end loop on event  
  MUONLoader->UnloadTracks();
}

