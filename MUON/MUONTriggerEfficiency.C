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

/*  */

//
// Macro for checking trigger integrated efficiency for dimuons.
// The efficiency is calculated with respect to the 3/4 coincidence.
// Author: Fabien Guerin, LPC Clermont-Ferrand, Jan. 2006
//

#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include "TBranch.h"
#include "TClonesArray.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TParticle.h"

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
#endif

// Upsilon(1S)


void MUONTriggerEfficiency (char filename[10]="galice.root",  Bool_t readFromRP = 0){
 
// output file
  char digitdat[100];
  char currentfile[100]; 
  sprintf(digitdat,"MUONTriggerEfficiency.out"); 
  
  FILE *fdat=fopen(digitdat,"w");          
  
  Int_t nevents,coincmuon,muonlpt,muonhpt;
  Int_t CoincMuPlus,CoincMuMinus;
  coincmuon=0;
  muonlpt=0;  
  muonhpt=0;
 
// Initialise AliRoot
   // Creating Run Loader and openning file containing Hits
   AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
    
   if (RunLoader ==0x0) {
       printf(">>> Error : Error Opening %s file \n",currentfile);
       return;
   }
             
   nevents = RunLoader->GetNumberOfEvents();          
     
 
   AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
   if (!readFromRP) {
       cout << " reading from digits \n";
       MUONLoader->LoadDigits("READ");
   } else {
       cout << " reading from RecPoints \n";
       MUONLoader->LoadRecPoints("READ");
   }
   MUONLoader->LoadHits("READ");
   
   // Creating MUON data container
   AliMUONData muondata(MUONLoader,"MUON","MUON");

   TClonesArray * globalTrigger;
   AliMUONGlobalTrigger * gloTrg; 
 
 
   for (Int_t ievent=0; ievent<nevents; ievent++) {    // event loop
       CoincMuPlus=0;
       CoincMuMinus=0;
       RunLoader->GetEvent(ievent);
         
     if (ievent%1000==0) printf("\t Event = %d\n",ievent);    

// Hits
     muondata.SetTreeAddress("H");    
        
    Int_t itrack, ntracks, NbHits[2][4];
    Int_t SumNbHits;
    for (Int_t j=0; j<2; j++) {
     for (Int_t jj=0; jj<4; jj++) {
      NbHits[j][jj]=0;
     }
    } 
    ntracks = (Int_t) muondata.GetNtracks();      
    for (itrack=0; itrack<ntracks; itrack++) { // Track loop
      muondata.GetTrack(itrack); 

      Int_t ihit, nhits;
      nhits = (Int_t) muondata.Hits()->GetEntriesFast();   
      AliMUONHit* mHit;
      for(ihit=0; ihit<nhits; ihit++) {
        mHit = static_cast<AliMUONHit*>(muondata.Hits()->At(ihit));
        Int_t Nch        = mHit->Chamber(); 
        Int_t hittrack   = mHit->Track();
        Float_t IdPart     = mHit->Particle();

        for (Int_t j=0;j<4;j++) {
         Int_t kch=11+j;
         if (hittrack==1 || hittrack==2) {
          if (Nch==kch && IdPart==-13 && NbHits[0][j]==0) NbHits[0][j]++; 
          if (Nch==kch && IdPart==13 && NbHits[1][j]==0) NbHits[1][j]++; 
          }       
        }
       }
       
      muondata.ResetHits();
      
    } // end track loop     
    
// 3/4 coincidence
    SumNbHits=NbHits[0][0]+NbHits[0][1]+NbHits[0][2]+NbHits[0][3];
    if (SumNbHits==3 || SumNbHits==4) CoincMuPlus=1;
       
    SumNbHits=NbHits[1][0]+NbHits[1][1]+NbHits[1][2]+NbHits[1][3];
    if (SumNbHits==3 || SumNbHits==4) CoincMuMinus=1;
    
    if (CoincMuPlus==1 && CoincMuMinus==1) coincmuon++;         
           
// Trigger
    if (!readFromRP) {
	muondata.SetTreeAddress("D,GLT"); 
	muondata.GetTriggerD();
    } else {    
	muondata.SetTreeAddress("RC,TC"); 
	muondata.GetTrigger();
    }
   
    globalTrigger = muondata.GlobalTrigger();

    Int_t nglobals = (Int_t) globalTrigger->GetEntriesFast(); // should be 1

    for (Int_t iglobal=0; iglobal<nglobals; iglobal++) { // Global Trigger
      gloTrg = static_cast<AliMUONGlobalTrigger*>(globalTrigger->At(iglobal));
      
       if (gloTrg->PairUnlikeLpt()>=1) muonlpt++;
       if (gloTrg->PairUnlikeHpt()>=1) muonhpt++;
                               
    } // end of loop on Global Trigger         
    
    muondata.ResetTrigger();    

  } // end loop on event  

   MUONLoader->UnloadHits();

  if (!readFromRP) {
      MUONLoader->UnloadDigits();  
  } else {    
      MUONLoader->UnloadRecPoints();
  }  
 
  // calculate efficiency with as a ref. at least 3/4 planes fired
  Float_t efficiencylpt,efficiencyhpt;
  Double_t coincmu,lptmu,hptmu;  
  Float_t error;  
  coincmu=Double_t(coincmuon);
  lptmu=Double_t(muonlpt); 
  hptmu=Double_t(muonhpt); 


  // output
  fprintf(fdat,"\n"); 
  fprintf(fdat," Number of events = %d \n",nevents);  
  fprintf(fdat," Number of events with 3/4 coinc = %d \n",coincmuon);  
  fprintf(fdat," Nomber of dimuons with 3/4 coinc Lpt cut = %d \n",muonlpt);  
  fprintf(fdat," Number of dimuons with 3/4 coinc Hpt cut = %d \n",muonhpt);  
  fprintf(fdat,"\n");
  
  efficiencylpt=lptmu/coincmu;  
  error=efficiencylpt*TMath::Sqrt((lptmu+coincmu)/(lptmu*coincmu));  
  fprintf(fdat," Efficiency Lpt cut = %4.4f +/- %4.4f\n",efficiencylpt,error); 
 
  efficiencyhpt=hptmu/coincmu; 
  error=efficiencyhpt*TMath::Sqrt((hptmu+coincmu)/(hptmu*coincmu));          
  fprintf(fdat," Efficiency Hpt cut = %4.4f +/- %4.4f\n",efficiencyhpt,error);
      

  fclose(fdat);  
}





