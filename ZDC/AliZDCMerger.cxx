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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  			ZDC event merging class                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// --- ROOT system
#include <iostream.h>
#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TSystem.h>

// --- AliRoot header files
#include "AliZDCMerger.h"
#include "AliZDC.h"
#include "AliZDCHit.h"
#include "AliZDCDigit.h"
#include "AliZDCFragment.h"
#include "AliRun.h"
#include "AliDetector.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"

ClassImp(AliZDCMerger)

//____________________________________________________________________________
AliZDCMerger::AliZDCMerger()
{
// Default constructor    
    fMerge       = kDigitize;
    fNEvBgr	 = 0;
    fFnBgr       = 0;
    fBgrFile     = 0;
    fTrHBgr      = 0;
    fTrSDBgr     = 0;
    fImpPar      = 0;
    fSpecn       = 0;
    fSpecp       = 0;
    fFreeSpn     = 0;
    fFreeSpp     = 0;
    fFnSpecn     = 0;
    fSpecnFile   = 0;
    fFnSpecp     = 0;
    fSpecpFile   = 0;
    
}

//____________________________________________________________________________
AliZDCMerger::~AliZDCMerger()
{
// Destructor
    if (fBgrFile)    delete fBgrFile;
    if (fTrHBgr)     delete fTrHBgr;
    if (fTrSDBgr)    delete fTrSDBgr;
    if (fHitsBgr)    delete fHitsBgr;
    if (fSpecnFile)  delete fSpecnFile;
    if (fSpecpFile)  delete fSpecpFile;
}

//____________________________________________________________________________
void AliZDCMerger::InitMerging()
{
    // Hits tree, impact parameter, num. of spectators n & p 
    // 		in background (full Hijing) event
    Float_t b;
    Int_t nspecn, nspecp;
    Background(b, nspecn, nspecp);

    // Production of nuclear fragments -> num. of FREE spectators n & p
    Fragmentation(b, nspecn, nspecp, fFreeSpn, fFreeSpp);
    
    // Extract from spectators distribution the signal events:
    // NFreeSpectatorN spectator n & NFreeSpectatorP spectator p 
    Mixing(fFreeSpn, fFreeSpp);
}

//____________________________________________________________________________
void AliZDCMerger::Background(Float_t &fImpPar, Int_t &fSpecn, Int_t &fSpecp)
{
    
    // --- Open the background file
    if (fMerge && !fBgrFile) fBgrFile = OpenBgrFile();
    
    // --- Read from the TreeE impact parameter (b),
    //     # of spectators n and p (fSpecn, fSpecp)
    fBgrFile->cd();

    // Get AliRun object from file or create it if not on file
    gAlice = (AliRun*)fBgrFile->Get("gAlice");
    if (!gAlice) {
      gAlice = (AliRun*)fBgrFile->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) {
	  printf("\n create new gAlice object");
	  gAlice = new AliRun("gAlice","Alice test program");
      }
    }

//    gAlice = (AliRun*)fBgrFile->Get("gAlice");
    AliHeader *header = gAlice->GetHeader();
    AliGenEventHeader* mcHeader = header->GenEventHeader();
    fImpPar = ((AliGenHijingEventHeader*) mcHeader)->ImpactParameter();
    fSpecn  = ((AliGenHijingEventHeader*) mcHeader)->Spectatorsn();
    fSpecp  = ((AliGenHijingEventHeader*) mcHeader)->Spectatorsp();
//      printf("\n	HIJING simulation - b = %f fm, Nspecn = %d, Nspecp = %d\n",fImpPar,fSpecn,fSpecp);
}

//____________________________________________________________________________
TFile* AliZDCMerger::OpenBgrFile()
{
    // Initialise background event
    TFile *file = new TFile(fFnBgr,"UPDATE");
//      printf("\n AliZDCMerger --- Background event -> %s file opened \n", fFnBgr);
    fHitsBgr = new TClonesArray("AliZDCHit",1000);
    return file;
}

//____________________________________________________________________________
void AliZDCMerger::Fragmentation(Float_t fImpPar, Int_t fSpecn, Int_t fSpecp,
                                 Int_t &fFreeSpn, Int_t &fFreeSpp)
{
    Int_t j, zz[100], nn[100], nAlpha, Ztot, Ntot;
    AliZDCFragment *frag = new AliZDCFragment(fImpPar);
    for(j=0; j<=99; j++){
       zz[j] =0;
       nn[j] =0;
    }

    // Fragments generation
    frag->GenerateIMF(zz, nAlpha);

    // Attach neutrons
    Ztot=0;
    Ntot=0;
    frag->AttachNeutrons(zz, nn, Ztot, Ntot);
    fFreeSpn = fSpecn-Ztot-2*nAlpha;
    fFreeSpp = fSpecp-Ntot-2*nAlpha;
//      printf("\n	Fragmentation -> FreeSpn = %d, FreeSpp = %d\n",fFreeSpn,fFreeSpp);
}

//____________________________________________________________________________
void AliZDCMerger::Mixing(Int_t fFreeSpn, Int_t fFreeSpp)
{

//      printf("\n	AliZDCMerger->Mixing\n");
        
    // ### Background event Hits ###########################################
    fBgrFile->cd();
//    fBgrFile->ls();
    
    AliZDC *ZDC = (AliZDC *)gAlice->GetModule("ZDC");
//    if(ZDC) printf("\n	Ho trovato lo ZDC!\n");

    // Hits tree
    if(fTrHBgr) delete fTrHBgr;
    fTrHBgr = 0;    
    // SDigits tree
    if(fTrSDBgr) delete fTrSDBgr;
    fTrSDBgr = 0;    

//    // Read from Event Tree the # of events of the Hijing file 
//    TTree *thead = (TTree *) fBgrFile->Get("TE");
//    if(!te){
//      printf("\n	ERROR! -> File %s does not contain TreeE\n");
//    }
//    fNEvBgr = (Int_t) te->GetEntries();
    fNEvBgr = 0; // Let's suppose to have 1 full Hijing event per file

    // Hits Tree
    char treeBgrName[20];
    sprintf(treeBgrName,"TreeH%d",fNEvBgr);
    fTrHBgr = (TTree*)gDirectory->Get(treeBgrName); // TreeH
    if(!fTrHBgr){
      printf("\n ERROR -> Can't find TreeH%d in background file\n",fNEvBgr);
    }    
    Int_t ntracks =  (Int_t) fTrHBgr->GetEntries();
//    printf("\n	--- ntracks = %d\n",ntracks);

    // SDigits Tree
    char treeSBgrName[20];
    sprintf(treeSBgrName,"TreeS%d",fNEvBgr);
    fTrSDBgr = (TTree*)gDirectory->Get(treeSBgrName); // TreeH
    if(!fTrSDBgr){
      printf("\n ERROR -> Can't find TreeS%d in background file\n",fNEvBgr);
    }    

    Int_t itrack, i, volume[2], detector, quadrant;
    TClonesArray &sdigits = *fHitsBgr;	// SDigits TCArray
    Float_t hits[10];
    // --- Tracks loop
    for(itrack=0; itrack<ntracks; itrack++){
//       printf("		itrack = %d \n", itrack);
       fTrHBgr->GetEvent(itrack);
       
       fTrack = itrack;
       Int_t NMhits  = 0;
       for(AliZDCHit* zdcHit=(AliZDCHit*)ZDC->FirstHit(-1); 
           zdcHit;
    	   zdcHit = (AliZDCHit*)ZDC->NextHit()){ 

	  for(i=0; i<2; i++) volume[i] = zdcHit->GetVolume(i);
//	  printf("\n		volume[0] = %d volume[1]= %d \n",volume[0], volume[1]);
	  if(volume[0] != 3){
	    // Background hits	
	    hits[7] = zdcHit->GetLightPMQ();
	    hits[8] = zdcHit->GetLightPMC();
//	    printf("\n	Prima ### Background -> PMQ = %f, PMC = %f\n",hits[7],hits[8]);

            //Signal hits
	    detector = volume[0];
	    quadrant = volume[1];
	    ExtractSignal(detector, quadrant, fQuadLight, fComLight);	   
	    hits[7] += fQuadLight;
	    hits[8] += fComLight;
//	    printf("\n	Bckg + Signal -> PMQ = %f, PMC = %f\n",hits[7],hits[8]);

            zdcHit->SetLightPMQ(hits[7]);
            zdcHit->SetLightPMC(hits[8]);
	    hits[7] = zdcHit->GetLightPMQ();
	    hits[8] = zdcHit->GetLightPMC();
//	    printf("\n	Dopo ### Background -> PMQ = %f, PMC = %f\n",hits[7],hits[8]);
	    	    
	    new (sdigits[NMhits++]) AliZDCHit(zdcHit);
//            printf("\n	NMhits = %d\n",NMhits);
	  }
          fBgrFile->cd();
          fTrSDBgr->Fill();
          fTrSDBgr->Write(0,TObject::kOverwrite);
       } //Hits loop
       Digitize();
    } // Tracks loop
        
//    fBgrFile->Close();
        
}

//____________________________________________________________________________
void AliZDCMerger::ExtractSignal(Int_t detector, Int_t quadrant,
                                    Float_t &fQuadLight, Float_t &fComLight)
{

//    printf("\n	Entering ExtractSignal method -> detector = %d quadrant = %d\n",
//           detector, quadrant);

    // Connect spectator n histo's file
    fFnSpecn = gSystem->ExpandPathName("$ALICE/$ALICE_LEVEL/ZDC/ZNsignal.root");
    fSpecnFile = TFile::Open(fFnSpecn,"R");
    fSpecnFile->cd();
//    fSpecnFile->ls();
    TH1F *hPMQ1zn = (TH1F*) gDirectory->Get("hPMQ1zn;1");
    TH1F *hPMQ2zn = (TH1F*) gDirectory->Get("hPMQ2zn;1");
    TH1F *hPMQ3zn = (TH1F*) gDirectory->Get("hPMQ3zn;1");
    TH1F *hPMQ4zn = (TH1F*) gDirectory->Get("hPMQ4zn;1");
    TH1F *hPMC1zn = (TH1F*) gDirectory->Get("hPMC1zn;1");
    TH1F *hPMC2zn = (TH1F*) gDirectory->Get("hPMC2zn;1");
    TH1F *hPMC3zn = (TH1F*) gDirectory->Get("hPMC3zn;1");
    TH1F *hPMC4zn = (TH1F*) gDirectory->Get("hPMC4zn;1");
//    Axis_t x = hPMQ1zn -> GetRandom();
//    printf("	hPMQ1zn -> GetRandom() = %f\n",x);
    
    // Connect spectator p histo's file
    fFnSpecp = gSystem->ExpandPathName("$ALICE/$ALICE_LEVEL/ZDC/ZPsignal.root");
    fSpecpFile = TFile::Open(fFnSpecp,"R");
    fSpecpFile->cd(); 
//    fSpecpFile->ls();
    TH1F *hPMQ1zp = (TH1F*) gDirectory->Get("hPMQ1zp;1");
    TH1F *hPMQ2zp = (TH1F*) gDirectory->Get("hPMQ2zp;1");
    TH1F *hPMQ3zp = (TH1F*) gDirectory->Get("hPMQ3zp;1");
    TH1F *hPMQ4zp = (TH1F*) gDirectory->Get("hPMQ4zp;1");
    TH1F *hPMC1zp = (TH1F*) gDirectory->Get("hPMC1zp;1");
    TH1F *hPMC2zp = (TH1F*) gDirectory->Get("hPMC2zp;1");
    TH1F *hPMC3zp = (TH1F*) gDirectory->Get("hPMC3zp;1");
    TH1F *hPMC4zp = (TH1F*) gDirectory->Get("hPMC4zp;1");

    if(detector == 1){ // --- ZN
      if(quadrant == 1){
        fQuadLight = (Float_t) hPMQ1zn -> GetRandom();
	fComLight = (Float_t) hPMC1zn -> GetRandom();
      }
      else if(quadrant == 2){
        fQuadLight = (Float_t) hPMQ2zn -> GetRandom();
	fComLight = (Float_t) hPMC2zn -> GetRandom();
      }
      else if(quadrant == 3){
        fQuadLight = (Float_t) hPMQ3zn -> GetRandom();
	fComLight = (Float_t) hPMC3zn -> GetRandom();
      }
      else if(quadrant == 4){
        fQuadLight = (Float_t) hPMQ4zn -> GetRandom();
	fComLight = (Float_t) hPMC4zn -> GetRandom();
      }
    }

    else if(detector == 2){ // --- ZP
      fSpecpFile->cd();	// Connect spectator p histo's file
//      fSpecpFile->ls();
      if(quadrant == 1){
        fQuadLight = (Float_t) hPMQ1zp -> GetRandom();
	fComLight = (Float_t) hPMC1zp -> GetRandom();
      }
      else if(quadrant == 2){
        fQuadLight = (Float_t) hPMQ2zp -> GetRandom();
	fComLight = (Float_t) hPMC2zp -> GetRandom();
      }
      else if(quadrant == 3){
        fQuadLight = (Float_t) hPMQ3zp -> GetRandom();
	fComLight = (Float_t) hPMC3zp -> GetRandom();
      }
      else if(quadrant == 4){
        fQuadLight = (Float_t) hPMQ4zp -> GetRandom();
	fComLight = (Float_t) hPMC4zp -> GetRandom();
      }
    }
//    printf("	---	Exiting ExtractSignal -> fQuadLight = %f, fComLight = %f",
//           fQuadLight,fComLight);

}

//____________________________________________________________________________
void AliZDCMerger::Digitize()
{

//  printf("\n	AliZDCMerger->Digitize()");

  AliZDC *ZDC = (AliZDC *)gAlice->GetModule("ZDC");
//  if(ZDC) printf("\n 	Ho trovato lo ZDC!\n");

  Int_t itrack, lightQ, lightC, sector[2], digit;
  Int_t PMCZN = 0, PMCZP = 0, PMQZN[4], PMQZP[4], PMZEM = 0;
  Int_t i;
  for(i=0; i<4; i++){
     PMQZN[i] = 0;
     PMQZP[i] = 0;
  }

  // ### Digitization "on the flight" (no merging) #######################
  if(fMerge == kDigitize){
//    printf("\n		fMerge == kDigitize\n");
    fTrHBgr = gAlice->TreeH(); // TTree
    TTree *treeH = gAlice->TreeH();
    if(!fTrHBgr){
      printf("\n ERROR -> Can't find TreeH%d in background file\n",fNEvBgr);
    }    
    Int_t ntracks =  (Int_t) treeH->GetEntries();
//    printf("\n	--- ntracks = %d\n",ntracks);
        
    // Loop over tracks
    for(itrack=0; itrack<ntracks; itrack++){
       gAlice->ResetHits();
       gAlice->TreeH()->GetEvent(itrack);
       
       // Loop over hits
       for(AliZDCHit* zdcHit=(AliZDCHit*)ZDC->FirstHit(-1); 
           zdcHit;
    	   zdcHit = (AliZDCHit*)ZDC->NextHit()){ 
	  sector[0] = zdcHit->GetVolume(0);
    	  sector[1] = zdcHit->GetVolume(1);
    	  lightQ    = Int_t(zdcHit->GetLightPMQ());
    	  lightC    = Int_t(zdcHit->GetLightPMC());
//	  printf("\n	Digitise -> DET. = %d, quad = %d", sector[0], sector[1]);
//	  printf("		    PMQ = %d, PMC = %d",  lightQ, lightC);
          
    	  if(sector[0] == 1){	//ZN 
    	    PMCZN = PMCZN + lightC;
    	    PMQZN[sector[1]-1] = PMQZN[sector[1]-1] + lightQ;
    	  }
          else if(sector[0] == 2){	//ZP 
    	    PMCZP = PMCZP + lightC;
    	    PMQZP[sector[1]-1] = PMQZP[sector[1]-1] + lightQ;
    	  }
          else if(sector[0] == 3){	//ZEM 
    	    PMZEM = PMZEM + lightC;
    	  }
       } // hits loop
    } // tracks loop
  } // if(fMerge)

  // ### Merging and digitization #######################################
  else if(fMerge == kMerge){
//    printf("\n		fMerge == kMerge\n");

//    Int_t ntracks =  (Int_t) fTrHBgr->GetEntries();
//    printf("\n        --- ntracks = %d\n",ntracks);
//        
//    // Loop over tracks
//    for(itrack=0; itrack<ntracks; itrack++){
       fTrHBgr->GetEvent(fTrack);
       
//       printf("\n\n	Track # %d --- Digitise -> ", fTrack);
       // Loop over hits
       for(AliZDCHit* zdcHit=(AliZDCHit*)ZDC->FirstHit(-1); 
           zdcHit;
    	   zdcHit = (AliZDCHit*)ZDC->NextHit()){ 
	  sector[0] = zdcHit->GetVolume(0);
    	  sector[1] = zdcHit->GetVolume(1);
    	  lightQ    = Int_t(zdcHit->GetLightPMQ());
    	  lightC    = Int_t(zdcHit->GetLightPMC());
//	  printf("\n	DET. = %d, quad = %d,PMQ = %d, PMC = %d", 
//	         sector[0], sector[1],lightQ, lightC);
          
    	  if(sector[0] == 1){	//ZN 
    	    PMCZN = PMCZN + lightC;
    	    PMQZN[sector[1]-1] = PMQZN[sector[1]-1] + lightQ;
    	  }
          else if(sector[0] == 2){	//ZP 
    	    PMCZP = PMCZP + lightC;
    	    PMQZP[sector[1]-1] = PMQZP[sector[1]-1] + lightQ;
    	  }
          else if(sector[0] == 3){	//ZEM 
    	    PMZEM = PMZEM + lightC;
    	  }
       } // hits loop
//    } // tracks loop
  } // if(fMerge)

      
  // Create digits for ZN
  Int_t PedValue;
  sector[0] = 1; // Detector = ZN
  sector[1] = 0; // Common PM ADC
  digit = Phe2ADCch(1, 0, PMCZN);
//  printf("\n\n	ZN ###	PMCZN = %d	ADCZN = %d",PMCZN, digit);
  ZDC->AddDigit(sector, digit);
  Int_t j;
  for(j=0; j<4; j++){
    sector[1] = j+1; // Towers PM ADCs
    digit = Phe2ADCch(1, j+1, PMQZN[j]);
//    printf("\n		PMQZN[%d] = %d	ADCZN[%d] = %d",j,PMQZN[j],j,digit);
    PedValue = AddPedestal();
    digit += PedValue;
//    printf("	PedValue = %d",PedValue);
    ZDC->AddDigit(sector, digit);
  }
//    printf("\n");
  
  // Create digits for ZP
  sector[0] = 2; // Detector = ZP
  sector[1] = 0; // Common PM ADC
  digit = Phe2ADCch(2, 0, PMCZP);
//  printf("\n	ZP --- PMCZP = %d	ADCZP = %d",PMCZP,digit);
  ZDC->AddDigit(sector, digit);
  for(j=0; j<4; j++){
    sector[1] = j+1; // Towers PM ADCs
    digit = Phe2ADCch(2, j+1, PMQZP[j]);
//    printf("\n	       PMQZP[%d] = %d	ADCZP[%d] = %d",j,PMQZP[j],j,digit);
    PedValue = AddPedestal();
    digit += PedValue;
//    printf("	PedValue = %d",PedValue);
    ZDC->AddDigit(sector, digit);
  }
//    printf("\n");
  
  // Create digits for ZEM
  sector[0] = 3; // Detector = ZZEM
  sector[1] = 0; // Single PM ADC
  digit  = Phe2ADCch(3, 0, PMZEM);
//  printf("\n	ZEM *** PMZEM = %d	ADCZEM = %d",PMZEM,digit);
    PedValue = AddPedestal();
    digit += PedValue;
//    printf("	PedValue = %d\n",PedValue);
  ZDC->AddDigit(sector, digit); 

}

//_____________________________________________________________________________
Int_t AliZDCMerger::Phe2ADCch(Int_t Det, Int_t Quad, Int_t Light)
{
  // Evaluation of the ADC channel corresponding to the light yield Light

  if(gAlice->GetDebug() > 0){
//    printf("\n  Phe2ADCch -> Detector = %d, Quadrant = %d, Light = %d\n", Det, Quad, Light);
  }
  
  Int_t ADCch = 0;
  // Parameters for conversion of light yield in ADC channels
  Float_t fPMGain[3][5];      // PM gain
  Float_t fADCRes;            // ADC conversion factor
  
  Int_t j,i;
  for(i=0; i<3; i++){
     for(j=0; j<5; j++){
        fPMGain[i][j]   = 100000.;
     }
  }
  fADCRes   = 0.00000064; // ADC Resolution: 250 fC/ADCch
  
  ADCch = (Int_t) (Light*fPMGain[Det-1][Quad]*fADCRes);
     
  return ADCch;
}

//_____________________________________________________________________________
Int_t AliZDCMerger::AddPedestal()
{
  // --- Pedestal value -> extracted from a gaussian distribution
  // obtained from the beam test on the ZEM prototype (Aug. 2000)
  
  Int_t PedValue;
  Float_t PedMean  = 50.;
  Float_t PedWidth = 10.;
  
  PedValue    = (Int_t) gRandom->Gaus(PedMean,PedWidth);
  
  return PedValue;
}
