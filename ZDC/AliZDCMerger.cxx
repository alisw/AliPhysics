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
#include <TNtuple.h>
#include <TSystem.h>

// --- AliRoot header files
#include "AliZDCMerger.h"
#include "AliZDC.h"
#include "AliZDCHit.h"
#include "AliZDCMergedHit.h"
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
//   fMerge       = kMerge;
    fFnBgr       = 0;
    fBgrFile     = 0;
    fNEvBgr	 = 0;
    fTrHBgr      = 0;
    fHitsBgr     = 0;
    fTrSDBgr     = 0;
    fTrDBgr      = 0;
    fImpPar      = 0;
    fSpecn       = 0;
    fSpecp       = 0;
    fFreeSpn     = 0;
    fFreeSpp     = 0;
    fFnSpecn     = 0;
    fSpecnFile   = 0;
    fFnSpecp     = 0;
    fSpecpFile   = 0;
    fNMhits	 = 0;
     
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
    Mixing();
}

//____________________________________________________________________________
void AliZDCMerger::Background(Float_t &fImpPar, Int_t &fSpecn, Int_t &fSpecp)
{
    
    // --- Open the background file
    if (fMerge && !fBgrFile) fBgrFile = OpenBgrFile();
    
    // --- Read from the TreeE impact parameter (b),
    //     # of spectators n and p (fSpecn, fSpecp)
    fBgrFile->cd();

//    // Get AliRun object from file or create it if not on file
//    gAlice = (AliRun*)fBgrFile->Get("gAlice");
//    if (!gAlice) {
//      gAlice = (AliRun*)fBgrFile->Get("gAlice");
//      if (gAlice) printf("AliRun object found on file\n");
//      if (!gAlice) {
//	  printf("\n create new gAlice object");
//	  gAlice = new AliRun("gAlice","Alice test program");
//      }
//    }

    AliHeader *header = gAlice->GetHeader();
    AliGenEventHeader* mcHeader = header->GenEventHeader();
    fImpPar = ((AliGenHijingEventHeader*) mcHeader)->ImpactParameter();
    fSpecn  = ((AliGenHijingEventHeader*) mcHeader)->Spectatorsn();
    fSpecp  = ((AliGenHijingEventHeader*) mcHeader)->Spectatorsp();
    printf("\n	HIJING simulation - b = %f fm, Nspecn = %d, Nspecp = %d\n",fImpPar,fSpecn,fSpecp);
}

//____________________________________________________________________________
TFile* AliZDCMerger::OpenBgrFile()
{
    // Initialise background event
    TFile *file = new TFile(fFnBgr,"UPDATE");
    printf("\n AliZDCMerger --- Background event -> %s file opened \n", fFnBgr);
    fHitsBgr = new TClonesArray("AliZDCHit",1000);
    fMHits   = new TClonesArray("AliZDCMergedHit",1000);
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
    printf("\n	Fragmentation -> FreeSpn = %d, FreeSpp = %d\n",fFreeSpn,fFreeSpp);
}

//____________________________________________________________________________
void AliZDCMerger::Mixing()
{
    
//    printf("\n	AliZDCMerger->Mixing\n");
       
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

//    fNEvBgr = 0; // Let's suppose to have 1 full Hijing event per file

    // Hits tree
    char treeBgrName[20];
    sprintf(treeBgrName,"TreeH%d",fNEvBgr);
    fTrHBgr = (TTree*)gDirectory->Get(treeBgrName); 
    if(!fTrHBgr){
      printf("\n ERROR -> Can't find TreeH%d in background file\n",fNEvBgr);
    }    
//    fTrHBgr->Print();
    
    // Branch address
    TBranch *branch;
    char branchname[20];
    sprintf(branchname,"%s",ZDC->GetName());
    if(fTrHBgr && fHitsBgr){
//      printf("\n	fTrHBgr!=0 && fHitsBgr!=0\n");
      branch = fTrHBgr->GetBranch(branchname);
      if(branch) branch->SetAddress(&fHitsBgr);
    }
    
    Int_t ntracks =  (Int_t) fTrHBgr->GetEntries();
//    printf("\n	--- ntracks = %d\n\n", ntracks);
    
    Int_t itrack, nhits, ihit, j, sector[2];
    AliZDCHit* zdcHit;
    AliZDCMergedHit *MHit;
    Float_t MHits[7];
    TClonesArray &sdigits = *fMHits;	// SDigits TCArray
    fNMhits = 0;

    // --- Tracks loop
    for(itrack=0; itrack<ntracks; itrack++){
//       printf("               itrack = %d", itrack);
//       gAlice->ResetHits();
       fTrHBgr->GetEvent(itrack);
       
//       for(AliZDCHit* zdcHit=(AliZDCHit*)ZDC->FirstHit(-1); 
//	   zdcHit;
//	   zdcHit = (AliZDCHit*)ZDC->NextHit()){ 

       nhits = fHitsBgr->GetEntries();
//       nhits = ZDCdp->GetNhits();
//       printf(" 	      nhits = %d \n", nhits);
       for(ihit=0; ihit<nhits; ihit++){    
	  zdcHit = (AliZDCHit*) fHitsBgr->UncheckedAt(ihit);

	    for(j=0; j<2; j++) sector[j] = zdcHit->GetVolume(j);
	    MHits[0] = zdcHit->GetPrimKinEn();
	    MHits[1] = zdcHit->GetXImpact();
	    MHits[2] = zdcHit->GetYImpact();
	    MHits[3] = zdcHit->GetSFlag();
	    MHits[4] = zdcHit->GetLightPMQ();
	    MHits[5] = zdcHit->GetLightPMC();
	    MHits[6] = zdcHit->GetEnergy();
	    MHit = new AliZDCMergedHit(sector, MHits);
//	    MHit->Print("");
  	    new((*fMHits)[fNMhits]) AliZDCMergedHit(*MHit);
	    new (sdigits[fNMhits])  AliZDCMergedHit(*MHit);
	    delete MHit;
	    fNMhits++;
	  }//Hits loop
	  
    } // Tracks loop
//    printf("\n	fNMhits (after bckg) = %d, \n",fNMhits); 
//    fMHits->Dump();
//    AliZDCMergedHit *pippo7 = (AliZDCMergedHit*) fMHits->At(7);
//    pippo7->Dump();
    
        
    // ### Signal event Hits ###########################################
    // --- Neutrons
    ExtractSignal(1);
        
    // --- Protons
    ExtractSignal(2);
//    printf("\n	fNMhits (after signal) = %d \n",fNMhits); 
        
}

//____________________________________________________________________________
void AliZDCMerger::ExtractSignal(Int_t SpecType)
{

// printf("\n	Entering in Extract Signal\n");
 
 Int_t NumEvents = 0;
 if(SpecType == 1){		// --- Signal for spectator neutrons
   fFnSpecn = gSystem->ExpandPathName("$ALICE/$ALICE_LEVEL/ZDC/ZNsignalntu.root");
   fSpecnFile = TFile::Open(fFnSpecn,"R");
   fSpecnFile->cd();
//   printf("\n	--- ExtractSignal x n: file %s opened\n", fFnSpecn);
   NumEvents = fFreeSpn;
 }
 else if(SpecType == 2){	// --- Signal for spectator protons
   fFnSpecp = gSystem->ExpandPathName("$ALICE/$ALICE_LEVEL/ZDC/ZPsignalntu.root");
   fSpecpFile = TFile::Open(fFnSpecp,"R");
   fSpecpFile->cd();
//   printf("\n	--- ExtractSignal x p: file %s opened\n", fFnSpecp);
   NumEvents = fFreeSpp;
 }
// printf("\n		# of free spectator = %d\n", NumEvents);
// printf("\n         fNMhits (before adding signal) = %d\n",fNMhits);

  TNtuple *ZDCSignal = (TNtuple*) gDirectory->Get("ZDCSignal");
  Int_t nentries = (Int_t) ZDCSignal->GetEntries();
//  printf("\n   # entries = %d\n", nentries);
  Int_t i, j, k;
  
  AliZDCMergedHit *MHit; 
  Float_t *entry, HitsSpec[7];
  Int_t iev, rnd, Volume[2];
  for(iev=0; iev<NumEvents; iev++){
     rnd = (Int_t) (1000*gRandom->Rndm());
//     printf("\n      rnd = %d\n", rnd);
     for(i=0; i<nentries; i++){  
  	ZDCSignal->GetEvent(i);
  	entry = ZDCSignal->GetArgs();
  	if(entry[0] == rnd){
//        printf("\n   entry[0] = %f\n", entry[0]);
          for(k=0; k<2; k++) Volume[k] = (Int_t) entry[k+1];
          for(j=0; j<7; j++){
             HitsSpec[j] = entry[j+3];
          }
	  MHit = new AliZDCMergedHit(Volume, HitsSpec);
	  new((*fMHits)[fNMhits++]) AliZDCMergedHit(*MHit);
	  delete MHit;
  	}
  	else if(entry[0] > rnd) break;
     }
  }
  
  if(SpecType ==1){
//    printf("\n         fNMhits (after n signal) = %d\n",fNMhits);
    fSpecnFile->Close();
  }
  else if(SpecType == 2){
//    printf("\n         fNMhits (after p signal) = %d\n",fNMhits);
    fSpecpFile->Close();
  }
      
}

//____________________________________________________________________________
void AliZDCMerger::Digitize(Int_t fNMhits, TClonesArray *fMHits)
{

//  printf("\n	AliZDCMerger->Digitize()");

  AliZDC *ZDC = (AliZDC *)gAlice->GetModule("ZDC");
//  if(ZDC) printf("\n 	Ho trovato lo ZDC!\n");

  Int_t lightQ, lightC, sector[2], digit;
  Int_t PMCZN = 0, PMCZP = 0, PMQZN[4], PMQZP[4], PMZEM1 = 0, PMZEM2 = 0;
  Int_t i;
  for(i=0; i<4; i++){
     PMQZN[i] = 0;
     PMQZP[i] = 0;
  }
  
  AliZDCMergedHit *MHit;
  Int_t imhit;
//  printf("\n	   fNMHits = %d\n", fNMhits);
  // Loop over SDigits
  for(imhit=0; imhit<fNMhits; imhit++){
     
     MHit = (AliZDCMergedHit*) fMHits->UncheckedAt(imhit);
     sector[0] = MHit->GetSector(0);
     sector[1] = MHit->GetSector(1);
  if((sector[1]!=1) && (sector[1]!=2) && (sector[1]!=3) && (sector[1]!=4)){
//     printf("\n  *** ERROR!!! sector[0] = %d, sector[1] = %d\n",
// 	      sector[0], sector[1]);
       sector[1] = 0;
   }
     lightQ    = Int_t(MHit->GetLightPMQ());
     lightC    = Int_t(MHit->GetLightPMC());
//     printf("\n    imhit = %d -> DET. = %d, quad = %d,PMQ = %d, PMC = %d", 
//      imhit,sector[0], sector[1],lightQ, lightC);
     
     if(sector[0] == 1){   	   //ZN 
       PMCZN = PMCZN + lightC;
       PMQZN[sector[1]-1] = PMQZN[sector[1]-1] + lightQ;
     }
     else if(sector[0] == 2){	   //ZP 
       PMCZP = PMCZP + lightC;
       PMQZP[sector[1]-1] = PMQZP[sector[1]-1] + lightQ;
     }
     else if(sector[0] == 3){	   //ZEM 
       if(sector[1] ==1) PMZEM1 = PMZEM1 + lightC;
       else              PMZEM2 = PMZEM2 + lightQ;
     }
  } // SDigits loop
      
  // ### Digits creation ###############################################
  // Create digits for ZN
  Int_t PedValue;
  sector[0] = 1; // Detector = ZN
  sector[1] = 0; // Common PM ADC
  digit = Phe2ADCch(1, 0, PMCZN);
//  printf("\n\n	ZN ###	PMCZN = %d	ADCZN = %d",PMCZN, digit);
  PedValue = AddPedestal();
  digit += PedValue;
//  printf("	PedValue = %d",PedValue);
  ZDC->AddDigit(sector, digit);
  Int_t j;
  for(j=0; j<4; j++){
    sector[1] = j+1; // Towers PM ADCs
    digit = Phe2ADCch(1, j+1, PMQZN[j]);
//    printf("\n		PMQZN[%d] = %d	phe	ADCZN[%d] = %d ADCch",j,PMQZN[j],j,digit);
    PedValue = AddPedestal();
    digit += PedValue;
//    printf("	PedValue = %d",PedValue);
    ZDC->AddDigit(sector, digit);
  }
//  printf("\n");
  
  // Create digits for ZP
  sector[0] = 2; // Detector = ZP
  sector[1] = 0; // Common PM ADC
  digit = Phe2ADCch(2, 0, PMCZP);
//  printf("\n	ZP --- PMCZP = %d	phe	ADCZP = %d ADCch",PMCZP,digit);
  PedValue = AddPedestal();
  digit += PedValue;
  printf("	PedValue = %d",PedValue);
  ZDC->AddDigit(sector, digit);
  for(j=0; j<4; j++){
    sector[1] = j+1; // Towers PM ADCs
    digit = Phe2ADCch(2, j+1, PMQZP[j]);
//    printf("\n	       PMQZP[%d] = %d	phe	ADCZP[%d] = %d ADCch",j,PMQZP[j],j,digit);
    PedValue = AddPedestal();
    digit += PedValue;
//    printf("	PedValue = %d",PedValue);
    ZDC->AddDigit(sector, digit);
  }
//  printf("\n");
  
  // Create digits for ZEM
  sector[0] = 3; 
  sector[1] = 1; // Detector = ZEM1
  digit  = Phe2ADCch(3, 1, PMZEM1);
//  printf("\n  ZEM *** PMZEM1 = %d      phe     ADCZEM1 = %d ADCch",PMZEM1,digit);
  PedValue = AddPedestal();
  digit += PedValue;
//  printf("  PedValue = %d\n",PedValue);
  ZDC->AddDigit(sector, digit); 
  sector[1] = 2; // Detector = ZEM2
  digit  = Phe2ADCch(3, 2, PMZEM2);
//  printf("\n  ZEM *** PMZEM2 = %d      phe     ADCZEM2 = %d ADCch",PMZEM2,digit);
  PedValue = AddPedestal();
  digit += PedValue;
//  printf("  PedValue = %d\n",PedValue);
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
