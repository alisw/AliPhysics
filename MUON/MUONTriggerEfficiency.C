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

// $Id$

/// \ingroup macros
/// \file MUONTriggerEfficiency.C
/// \brief Macro for checking trigger integrated efficiency for dimuons.
///
/// The efficiency is calculated with respect to the 3/4 coincidence.
///
/// \author Fabien Guerin, LPC Clermont-Ferrand, Jan. 2006

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
#include "AliCDBManager.h"

// MUON includes
#include "AliMUON.h"
#include "AliMUONHit.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONTrack.h"

#include "AliMUONDataInterface.h"
#include "AliMUONMCDataInterface.h"
#include "AliMUONVHitStore.h"
#include "AliMUONVTriggerStore.h"

#endif

// Upsilon(1S)


void MUONTriggerEfficiency(const char* filenameSim="galice_sim.root", 
                           const char* filenameRec="galice.root",  
                           Bool_t readFromRP = 0)
{
  
  // Set default CDB storage
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  // output file
  
  AliMUONMCDataInterface diSim(filenameSim);
  AliMUONDataInterface diRec(filenameRec);
  
  if (!diSim.IsValid())
  {
    cerr << "Cannot access " << filenameSim << endl;
    return;
  }
  
  if (!diRec.IsValid())
  {
    cerr << "Cannot access " << filenameRec << endl;
    return;
  }
  
  FILE* fdat = fopen("MUONTriggerEfficiency.out","w");          
  if (!fdat)
  {
    cerr << "Cannot create output file" << endl;
    return;
  }
  
  Int_t coincmuon,muonlpt,muonhpt;
  Int_t CoincMuPlus,CoincMuMinus;
  coincmuon=0;
  muonlpt=0;  
  muonhpt=0;
  
  Int_t nevents = diSim.NumberOfEvents();          
  
  for (Int_t ievent=0; ievent<nevents; ++ievent) 
  {   
    CoincMuPlus=0;
    CoincMuMinus=0;
    
    if (ievent%1000==0) printf("\t Event = %d\n",ievent);    
    
    // Hits
    
    Int_t NbHits[2][4];
    for (Int_t j=0; j<2; j++) 
    {
      for (Int_t jj=0; jj<4; jj++) 
      {
        NbHits[j][jj]=0;
      }
    } 
    
    Int_t ntracks = (Int_t) diSim.NumberOfTracks(ievent);
    
    for ( Int_t itrack=0; itrack<ntracks; ++itrack ) 
    {
      AliMUONVHitStore* hitStore = diSim.HitStore(ievent,itrack);      
      AliMUONHit* mHit;
      TIter next(hitStore->CreateIterator());
      
      while ( ( mHit = static_cast<AliMUONHit*>(next()) ) )
      {
        Int_t Nch = mHit->Chamber(); 
        Int_t hittrack = mHit->Track();
        Float_t IdPart = mHit->Particle();
        
        for (Int_t j=0;j<4;j++) 
        {
          Int_t kch=11+j;
          if (hittrack==1 || hittrack==2) 
          {
            if (Nch==kch && IdPart==-13 && NbHits[0][j]==0) NbHits[0][j]++; 
            if (Nch==kch && IdPart==13 && NbHits[1][j]==0) NbHits[1][j]++; 
          }       
        }
      }
    } // end track loop     
    
    // 3/4 coincidence
    Int_t SumNbHits=NbHits[0][0]+NbHits[0][1]+NbHits[0][2]+NbHits[0][3];
    if (SumNbHits==3 || SumNbHits==4) CoincMuPlus=1;
    
    SumNbHits=NbHits[1][0]+NbHits[1][1]+NbHits[1][2]+NbHits[1][3];
    if (SumNbHits==3 || SumNbHits==4) CoincMuMinus=1;
    
    if (CoincMuPlus==1 && CoincMuMinus==1) coincmuon++;         
          
    TString tree("D");
    if ( readFromRP ) tree = "R";
    
    AliMUONVTriggerStore* triggerStore = diRec.TriggerStore(ievent,tree.Data());

    AliMUONGlobalTrigger* gloTrg = triggerStore->Global();
    
    if (gloTrg->PairUnlikeLpt()>=1) muonlpt++;
    if (gloTrg->PairUnlikeHpt()>=1) muonhpt++;

  } // end loop on event  
 
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





