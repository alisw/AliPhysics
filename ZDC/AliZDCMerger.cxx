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
#include <Riostream.h>
#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TNtuple.h>
#include <TSystem.h>
#include <TROOT.h>

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

//int comp(const void *i,const void *j) {return *(int *)i - *(int *)j;}

//____________________________________________________________________________
AliZDCMerger::AliZDCMerger()
{
// Default constructor    
    //fMerge       = kDigitize	-> Only digitization
    //fMerge       = kMerge	-> Digitization + Merging
    fMerge       = kMerge;	
    fFnBgr       = 0;
    fBgrFile     = 0;
    fNEvBgr	 = 0;
    fHitsBgr     = 0;
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
    if (fSpecnFile)  delete fSpecnFile;
    if (fSpecpFile)  delete fSpecpFile;
}

//____________________________________________________________________________
void AliZDCMerger::InitMerging()
{
    // Hits tree, impact parameter, num. of spectators n & p 
    // 		in background (full Hijing) event
    Float_t b;
    Background(b, fSpecn, fSpecp);

    // Production of nuclear fragments -> num. of FREE spectators n & p
    Fragmentation(b, fSpecn, fSpecp, fFreeSpn, fFreeSpp);
    
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
    
    //    gAlice->GetEvent(fNEvBgr);  this is done in the steering macro
    AliHeader *header = gAlice->GetHeader();
    AliGenEventHeader* mcHeader = header->GenEventHeader();
    fImpPar = ((AliGenHijingEventHeader*) mcHeader)->ImpactParameter();
    Int_t dSpecn  = ((AliGenHijingEventHeader*) mcHeader)->Spectatorsn();
    Int_t dSpecp  = ((AliGenHijingEventHeader*) mcHeader)->Spectatorsp();
    // Until there is only 1 ZDC set the # of spectators must be divided by 2!!!
    fSpecn = dSpecn/2;
    fSpecp = dSpecp/2;
    printf("\n	HIJING ev. #%d - b = %f fm, Nspecn = %d, Nspecp = %d\n",
            fNEvBgr,fImpPar,fSpecn,fSpecp);
}

//____________________________________________________________________________
TFile* AliZDCMerger::OpenBgrFile()
{
    // Initialise background event
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(fFnBgr);
  if(!file)cerr<<"AliZDCMerger: background file "<<fFnBgr<<" not found\n";
  //    TFile *file = new TFile(fFnBgr,"UPDATE");
    printf("\n AliZDCMerger --- Background event -> %s file opened \n", fFnBgr);
    fHitsBgr = new TClonesArray("AliZDCHit",1000);
    fMHits   = new TClonesArray("AliZDCMergedHit",1000);
    return file;
}

//____________________________________________________________________________
void AliZDCMerger::Fragmentation(Float_t fImpPar, Int_t fSpecn, Int_t fSpecp,
                                 Int_t &fFreeSpn, Int_t &fFreeSpp)
{
    //printf("\n	Fragmentation -> fSpecn = %d, fSpecp = %d\n",fSpecn,fSpecp);
    Int_t j, zz[100], nn[100], nAlpha, ztot, ntot;
    AliZDCFragment *frag = new AliZDCFragment(fImpPar);
    for(j=0; j<=99; j++){
       zz[j] =0;
       nn[j] =0;
    }

    // Fragments generation
    frag->GenerateIMF(zz, nAlpha);

    // Attach neutrons
    ztot=0;
    ntot=0;
    frag->AttachNeutrons(zz, nn, ztot, ntot);
    fFreeSpn = fSpecn-ntot-2*nAlpha;
    fFreeSpp = fSpecp-ztot-2*nAlpha;
    if(fFreeSpn<0) fFreeSpn=0;
    if(fFreeSpp<0) fFreeSpp=0;
    //printf("\n			2*nAlpha = %d, ztot = %d, ntot = %d\n",2*nAlpha, ztot,    ntot);
    printf("\n	Fragmentation -> FreeSpn = %d, FreeSpp = %d\n",fFreeSpn,fFreeSpp);
}

//____________________________________________________________________________
void AliZDCMerger::Mixing()
{
    
    //printf("\n	AliZDCMerger->Mixing\n");
       
    // ### Background event Hits ###########################################
    fBgrFile->cd();
//    fBgrFile->ls();
    
   AliZDC *zdc = (AliZDC *)gAlice->GetModule("ZDC");
//    if(zdc) printf("\n	Ho trovato lo ZDC!\n");

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
    sprintf(branchname,"%s",zdc->GetName());
    if(fTrHBgr && fHitsBgr){
//      printf("\n	fTrHBgr!=0 && fHitsBgr!=0\n");
      branch = fTrHBgr->GetBranch(branchname);
      if(branch) branch->SetAddress(&fHitsBgr);
    }
    
    Int_t ntracks =  (Int_t) fTrHBgr->GetEntries();
    //printf("\n	--- ntracks = %d\n\n", ntracks);
    
    Int_t itrack, nhits, ihit, j, sector[2];
    AliZDCHit* zdcHit;
    AliZDCMergedHit *mHit;
    Float_t mHits[7];
    TClonesArray &sdigits = *fMHits;	// SDigits TCArray
    fNMhits = 0;

    // --- Tracks loop
    for(itrack=0; itrack<ntracks; itrack++){
       fTrHBgr->GetEvent(itrack);
       
       nhits = fHitsBgr->GetEntries();
//       printf(" 	      nhits = %d \n", nhits);
       for(ihit=0; ihit<nhits; ihit++){    
	  zdcHit = (AliZDCHit*) fHitsBgr->UncheckedAt(ihit);

	    for(j=0; j<2; j++) sector[j] = zdcHit->GetVolume(j);
	    mHits[0] = zdcHit->GetPrimKinEn();
	    mHits[1] = zdcHit->GetXImpact();
	    mHits[2] = zdcHit->GetYImpact();
	    mHits[3] = zdcHit->GetSFlag();
       	    mHits[4] = zdcHit->GetLightPMQ();
	    mHits[5] = zdcHit->GetLightPMC();
	    mHits[6] = zdcHit->GetEnergy();
	    mHit = new AliZDCMergedHit(sector, mHits);
//	    mHit->Print("");
  	    new((*fMHits)[fNMhits]) AliZDCMergedHit(*mHit);
	    new (sdigits[fNMhits])  AliZDCMergedHit(*mHit);
	    delete mHit;
	    fNMhits++;
	  }//Hits loop
	  
    } // Tracks loop
    //printf("	fNMhits (after bckg) = %d, \n",fNMhits);     
        
    // ### Signal event Hits ###########################################
    // --- Neutrons
    ExtractSignal(1);
        
    // --- Protons
    ExtractSignal(2);
    //printf("	fNMhits (after signal) = %d \n",fNMhits); 
        
}

//____________________________________________________________________________
void AliZDCMerger::ExtractSignal(Int_t SpecType)
{

// printf("\n	Entering in Extract Signal\n");
 
 Int_t numEvents = 0;
 if(SpecType == 1){		// --- Signal for spectator neutrons
   fFnSpecn = gSystem->ExpandPathName("$ALICE/$ALICE_LEVEL/ZDC/ZNsignalntu.root");
   fSpecnFile = TFile::Open(fFnSpecn,"R");
   fSpecnFile->cd();
   printf("\n	--- ExtractSignal x n: file %s opened\n", fFnSpecn);
   numEvents = fFreeSpn;
 }
 else if(SpecType == 2){	// --- Signal for spectator protons
   fFnSpecp = gSystem->ExpandPathName("$ALICE/$ALICE_LEVEL/ZDC/ZPsignalntu.root");
   fSpecpFile = TFile::Open(fFnSpecp,"R");
   fSpecpFile->cd();
   printf("\n	--- ExtractSignal x p: file %s opened\n", fFnSpecp);
   numEvents = fFreeSpp;
 }
 //printf("\n		# of free spectators = %d\n", numEvents);
 //printf("\n         fNMhits (before adding signal) = %d\n",fNMhits);

  TNtuple *zdcSignal = (TNtuple*) gDirectory->Get("ZDCSignal");
  Int_t nentries = (Int_t) zdcSignal->GetEntries();
  //printf("\n   # entries = %d\n", nentries);
  
  AliZDCMergedHit *mHit; 
  Float_t *entry, hitsSpec[7];
  Int_t pl, i, j, k, iev=0, rnd[125], volume[2];
  for(pl=0;pl<125;pl++){
     rnd[pl] = 0;
  }
  for(pl=0;pl<numEvents;pl++){
     rnd[pl] = (Int_t) (9999*gRandom->Rndm());
     if(rnd[pl] >= 9998) rnd[pl] = 9997;
     //printf("	rnd[%d] = %d\n",pl,rnd[pl]);     
  }
  // Sorting vector in ascending order with C function QSORT 
  qsort((void*)rnd,numEvents,sizeof(Int_t),comp);
  //for(pl=0;pl<numEvents;pl++){
  ////printf("	rnd[%d] = %d\n",pl,rnd[pl]);     
  //}
  do{
     for(i=0; i<nentries; i++){  
  	zdcSignal->GetEvent(i);
  	entry = zdcSignal->GetArgs();
  	if(entry[0] == rnd[iev]){
          for(k=0; k<2; k++) volume[k] = (Int_t) entry[k+1];
          for(j=0; j<7; j++){
             hitsSpec[j] = entry[j+3];
          }
          //printf("\n   i = %d, iev = %d, entry[0] = %f, rnd[%d] = %d ",i,iev,entry[0],iev,rnd[iev]);
	  mHit = new AliZDCMergedHit(volume, hitsSpec);
	  new((*fMHits)[fNMhits++]) AliZDCMergedHit(*mHit);
	  delete mHit;
  	}
  	else if(entry[0] > rnd[iev]){
	  iev++;
	  continue;
	}
     }
  }while(iev<numEvents);
  
  if(SpecType ==1){
    //printf("\n         fNMhits (after n signal) = %d\n",fNMhits);
    fSpecnFile->Close();
  }
  else if(SpecType == 2){
    //printf("\n         fNMhits (after p signal) = %d\n",fNMhits);
    fSpecpFile->Close();
  }
      
}

//____________________________________________________________________________
void AliZDCMerger::Digitize(Int_t fNMhits, TClonesArray *fMHits)
{
// Digitization

  printf("\n	AliZDCMerger->Digitize()");

  AliZDC *zdc = (AliZDC *)gAlice->GetModule("ZDC");
//  if(zdc) printf("\n 	Ho trovato lo ZDC!\n");
  Int_t lightQ, lightC, sector[2], digit;
  Int_t pmCZN = 0, pmCZP = 0, pmQZN[4], pmQZP[4], pmZEM1 = 0, pmZEM2 = 0;
  Int_t i;
  for(i=0; i<4; i++){
     pmQZN[i] = 0;
     pmQZP[i] = 0;
  }
  
  AliZDCMergedHit *mHit;
  Int_t imhit;
  printf("	   fNMHits = %d\n", fNMhits);
  // Loop over SDigits
  for(imhit=0; imhit<fNMhits; imhit++){
     
     mHit = (AliZDCMergedHit*) fMHits->UncheckedAt(imhit);
     sector[0] = mHit->GetSector(0);
     sector[1] = mHit->GetSector(1);
   // tmp -> c'erano i quadranti cannati!
  if((sector[1]!=1) && (sector[1]!=2) && (sector[1]!=3) && (sector[1]!=4)){
     printf("\n  *** ERROR!!! sector[0] = %d, sector[1] = %d\n",
 	      sector[0], sector[1]);
       sector[1] = 0;
   }
     lightQ    = Int_t(mHit->GetLightPMQ());
     lightC    = Int_t(mHit->GetLightPMC());
//     printf("    imhit = %d -> DET. = %d, quad = %d,PMQ = %d, PMC = %d\n", 
//      imhit,sector[0], sector[1],lightQ, lightC);
     
     if(sector[0] == 1){   	   //ZN 
       pmCZN = pmCZN + lightC;
       pmQZN[sector[1]-1] = pmQZN[sector[1]-1] + lightQ;
     }
     else if(sector[0] == 2){	   //ZP 
       pmCZP = pmCZP + lightC;
       pmQZP[sector[1]-1] = pmQZP[sector[1]-1] + lightQ;
     }
     else if(sector[0] == 3){	   //ZEM 
       if(sector[1] ==1) pmZEM1 = pmZEM1 + lightC;
       else              pmZEM2 = pmZEM2 + lightQ;
     }
  } // SDigits loop
      
  // ### Digits creation ###############################################
  // Create digits for ZN
  Int_t pedValue;
  sector[0] = 1; // Detector = ZN
  sector[1] = 0; // Common PM ADC
  digit = Phe2ADCch(1, 0, pmCZN);
  printf("\n\n	ZN ###	pmCZN = %d	ADCZN = %d",pmCZN, digit);
  pedValue = AddPedestal();
  digit += pedValue;
//  printf("	pedValue = %d",pedValue);
  zdc->AddDigit(sector, digit);
  Int_t j;
  for(j=0; j<4; j++){
    sector[1] = j+1; // Towers PM ADCs
    digit = Phe2ADCch(1, j+1, pmQZN[j]);
    printf("\n		pmQZN[%d] = %d	phe	ADCZN[%d] = %d adcCh",j,pmQZN[j],j,digit);
    pedValue = AddPedestal();
    digit += pedValue;
//    printf("	pedValue = %d",pedValue);
    zdc->AddDigit(sector, digit);
  }
  printf("\n");
  
  // Create digits for ZP
  sector[0] = 2; // Detector = ZP
  sector[1] = 0; // Common PM ADC
  digit = Phe2ADCch(2, 0, pmCZP);
  printf("\n	ZP --- pmCZP = %d	phe	ADCZP = %d adcCh",pmCZP,digit);
  pedValue = AddPedestal();
  digit += pedValue;
//  printf("	pedValue = %d",pedValue);
  zdc->AddDigit(sector, digit);
  for(j=0; j<4; j++){
    sector[1] = j+1; // Towers PM ADCs
    digit = Phe2ADCch(2, j+1, pmQZP[j]);
    printf("\n	       pmQZP[%d] = %d	phe	ADCZP[%d] = %d adcCh",j,pmQZP[j],j,digit);
    pedValue = AddPedestal();
    digit += pedValue;
//    printf("	pedValue = %d",pedValue);
    zdc->AddDigit(sector, digit);
  }
  printf("\n");
  
  // Create digits for ZEM
  sector[0] = 3; 
  sector[1] = 1; // Detector = ZEM1
  digit  = Phe2ADCch(3, 1, pmZEM1);
  printf("\n	ZEM *** pmZEM1 = %d      phe     ADCZEM1 = %d adcCh",pmZEM1,digit);
  pedValue = AddPedestal();
  digit += pedValue;
//  printf("  pedValue = %d\n",pedValue);
  zdc->AddDigit(sector, digit); 
  sector[1] = 2; // Detector = ZEM2
  digit  = Phe2ADCch(3, 2, pmZEM2);
  printf("\n	ZEM *** pmZEM2 = %d      phe     ADCZEM2 = %d adcCh\n",pmZEM2,digit);
  pedValue = AddPedestal();
  digit += pedValue;
//  printf("  pedValue = %d\n",pedValue);
  zdc->AddDigit(sector, digit); 
}

//_____________________________________________________________________________
Int_t AliZDCMerger::Phe2ADCch(Int_t Det, Int_t Quad, Int_t Light)
{
  // Evaluation of the ADC channel corresponding to the light yield Light

  if(gAlice->GetDebug() > 0){
//    printf("\n  Phe2ADCch -> Detector = %d, Quadrant = %d, Light = %d\n", Det, Quad, Light);
  }
  
  Int_t adcCh = 0;
  
  Int_t j,i;
  for(i=0; i<3; i++){
     for(j=0; j<5; j++){
        fPMGain[i][j]   = 100000.;
     }
  }
  fADCRes   = 0.00000064; // ADC Resolution: 250 fC/adcCh
  
  adcCh = (Int_t) (Light*fPMGain[Det-1][Quad]*fADCRes);
     
  return adcCh;
}

//_____________________________________________________________________________
Int_t AliZDCMerger::AddPedestal()
{
  // --- Pedestal value -> extracted from a gaussian distribution
  // obtained from the beam test on the ZEM prototype (Aug. 2000)
  
  Int_t pedValue;
  Float_t pedMean  = 50.;
  Float_t pedWidth = 5.;
  
  pedValue    = (Int_t) gRandom->Gaus(pedMean,pedWidth);
  
  return pedValue;
}
