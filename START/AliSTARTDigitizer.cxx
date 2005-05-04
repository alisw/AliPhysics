
/**************************************************************************
 * Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
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


#include <TTree.h> 
#include <TFile.h>
#include <TDirectory.h>
#include <TRandom.h>
#include <TArrayI.h>
#include <TError.h>
#include <TH1F.h>

#include "AliLog.h"
#include "AliSTARTDigitizer.h"
#include "AliSTART.h"
#include "AliSTARThit.h"
#include "AliSTARTdigit.h"
#include "AliRunDigitizer.h"
#include <AliDetector.h>
#include "AliRun.h"
#include <AliLoader.h>
#include <AliRunLoader.h>
#include <stdlib.h>
#include <Riostream.h>
#include <Riostream.h>

ClassImp(AliSTARTDigitizer)

//___________________________________________
  AliSTARTDigitizer::AliSTARTDigitizer()  :AliDigitizer()
{
// Default ctor - don't use it
  ;
}

//___________________________________________
AliSTARTDigitizer::AliSTARTDigitizer(AliRunDigitizer* manager) 
  :AliDigitizer(manager),
   fSTART(0),
   fHits(0),
   fdigits(0),
   ftimeTDC(0),
   fADC(0),
   ftimeTDCAmp(0),
   fADCAmp(0),
   fSumMult(0),
   fEff(0)
{
// ctor which should be used

  AliDebug(1,"processed");

  fSTART = 0;
  fHits = 0;
  fdigits = 0;

  ftimeTDC = new TArrayI(24); 
  fADC = new TArrayI(24); 
  ftimeTDCAmp = new TArrayI(24); 
  fADCAmp = new TArrayI(24); 
  fSumMult = new TArrayI(6); 
  

  TFile* file = TFile::Open("$ALICE_ROOT/START/PMTefficiency.root");
  fEff = (TH1F*) file->Get("hEff")->Clone();
  fEff->SetDirectory(NULL);
  file->Close();
  delete file;
}

//------------------------------------------------------------------------
AliSTARTDigitizer::~AliSTARTDigitizer()
{
// Destructor

  AliDebug(1,"START");

  delete ftimeTDC;
  delete fADC;
  delete fEff;
  delete ftimeTDCAmp;
  delete  fADCAmp;
  delete fSumMult;
}

 //------------------------------------------------------------------------
Bool_t AliSTARTDigitizer::Init()
{
// Initialization
  AliDebug(1," Init");
 return kTRUE;
}
 

//---------------------------------------------------------------------

void AliSTARTDigitizer::Exec(Option_t* /*option*/)
{

  /*
    Produde digits from hits
        digits is TObject and includes
	We are writing array if left & right  TDC
	left & right  ADC (will need for slow simulation)
	TOF first particle left & right
	mean time and time difference (vertex position)
	
  */

  //output loader 
  AliRunLoader *outRL = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  AliLoader * pOutStartLoader = outRL->GetLoader("STARTLoader");

  AliDebug(1,"start...");
  //input loader
  //
  // From hits to digits
  //
  Int_t hit, nhits;
  Int_t countE[24], sumMult[3];
  Int_t volume,pmt,tr,qt,qtAmp;
  Int_t  bestRightTDC,bestLeftTDC;
  Float_t time[24],besttime[24], timeGaus[24] ;
  Float_t channelWidth=25.; //ps
    //Q->T-> coefficients !!!! should be asked!!!
  Float_t ph2mV = 150./500.;
  Int_t threshold =50; //photoelectrons
  Int_t mV2channel=200000/(25*25);  //5V -> 200ns

  
  AliSTARThit  *startHit;
  TBranch *brHits=0;

  Int_t nFiles=fManager->GetNinputs();
  for (Int_t inputFile=0; inputFile<nFiles;  inputFile++) {
    if (inputFile < nFiles-1) {
      AliWarning(Form("ignoring input stream %d", inputFile));
      continue;
      
    }
    
    Float_t besttimeright=99999.;
    Float_t besttimeleft=99999.;
    Int_t timeDiff, meanTime;
    
    for (Int_t i0=0; i0<24; i0++)
      {
	time[i0]=besttime[i0]=timeGaus[i0]=99999; countE[i0]=0;
      }
    for ( Int_t imu=0; imu<3; imu++) sumMult[imu]=0;
    AliRunLoader * inRL = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inputFile));
    AliLoader * pInStartLoader = inRL->GetLoader("STARTLoader");
    if (!inRL->GetAliRun()) inRL->LoadgAlice();
    AliSTART *fSTART  = (AliSTART*)inRL ->GetAliRun()->GetDetector("START");

       //read Hits 
    pInStartLoader->LoadHits("READ");//probably it is necessary to load them before
    TClonesArray *fHits = fSTART->Hits ();
    TTree *th = pInStartLoader->TreeH();
    brHits = th->GetBranch("START");
    if (brHits) {
      fSTART->SetHitsAddressBranch(brHits);
    }else{
      AliError("Branch START hit not found");
      exit(111);
    } 
    Int_t ntracks    = (Int_t) th->GetEntries();
    
    if (ntracks<=0) return;
    // Start loop on tracks in the hits containers
    for (Int_t track=0; track<ntracks;track++) {
      brHits->GetEntry(track);
      nhits = fHits->GetEntriesFast();
      for (hit=0;hit<nhits;hit++) 
	{
	  startHit   = (AliSTARThit*) fHits->UncheckedAt(hit);
	  if (!startHit) {
 	    AliError("The unchecked hit doesn't exist");
	    break;
	  }
	  pmt=startHit->Pmt();
	  Int_t numpmt=pmt-1;
	  	  Double_t e=startHit->Etot();
	  volume = startHit->Volume();
	  
	  if(e>0 && RegisterPhotoE(e)) {
	    countE[numpmt]++;
	    besttime[numpmt] = startHit->Time();
	    if(besttime[numpmt]<time[numpmt])
	      {
		time[numpmt]=besttime[numpmt];
	      }
	  } 
	} //hits loop
    } //track loop
    
    //spread time right&left by 25ps   && besttime
    for (Int_t ipmt=0; ipmt<12; ipmt++){
      if(countE[ipmt] > threshold) {
	timeGaus[ipmt]=gRandom->Gaus(time[ipmt],25);
	if(timeGaus[ipmt]<besttimeleft) besttimeleft=timeGaus[ipmt]; //timeleft
	sumMult[0] +=  countE[ipmt];
	sumMult[1] +=  countE[ipmt];
      }
    }
    for ( Int_t ipmt=12; ipmt<24; ipmt++){
      if(countE[ipmt] > threshold) {
	timeGaus[ipmt]=gRandom->Gaus(time[ipmt],25); 
	if(timeGaus[ipmt]<besttimeright)  besttimeright=timeGaus[ipmt]; //timeright
	sumMult[0] +=  countE[ipmt];
	sumMult[2] +=  countE[ipmt];
      }	
    }
    //ADC features fill digits
   //folding with experimental time distribution
    Float_t c = 0.0299792; // cm/ps
    Float_t koef=(350.-69.7)/c;
    Float_t  besttimeleftR= besttimeleft;
    besttimeleft=koef+besttimeleft;
    bestLeftTDC=Int_t (besttimeleftR/channelWidth);
    bestRightTDC=Int_t (besttimeright/channelWidth);
    timeDiff=Int_t ((besttimeright-besttimeleftR)/channelWidth);
    meanTime=Int_t (((besttimeright+besttimeleft)/2.)/channelWidth);
    AliDebug(2,Form(" time in ns %f ", besttimeleft));
    for (Int_t i=0; i<24; i++)
      {
	//  fill TDC
	tr= Int_t (timeGaus[i]/channelWidth); 
	if(timeGaus[i]>50000) tr=0;
	
	//fill ADC
       	Int_t al= countE[i]; 
	// QTC procedure:
	// phe -> mV 0.3; 1MIP ->500phe -> ln (amp (mV)) = 5;
	// max 200ns, HIJING  mean 50000phe -> 15000mv -> ln = 15 (s zapasom)
	// channel 25ps
	if (al>threshold) {
	  qt=Int_t (TMath::Log(al*ph2mV) * mV2channel); 
	  qtAmp=Int_t (TMath::Log(al*10*ph2mV) * mV2channel);
	  fADC->AddAt(qt,i);
	  ftimeTDC->AddAt(tr,i);
	  fADCAmp->AddAt(qtAmp,i);
	  ftimeTDCAmp->AddAt(tr,i); //now is the same as non-amplified but will be change
	}
      } //pmt loop
    for (Int_t im=0; im<3; im++)
      { 
	if (sumMult[im]>threshold){
	  Int_t sum=Int_t (TMath::Log(sumMult[im]*ph2mV) * mV2channel);
	  Int_t sumAmp=Int_t (TMath::Log(10*sumMult[im]*ph2mV) * mV2channel);
	  fSumMult->AddAt(sum,im);
	  fSumMult->AddAt(sumAmp,im+3);
	}
      }	

    fSTART->AddDigit(bestRightTDC,bestLeftTDC,meanTime,timeDiff,fSumMult,
		     ftimeTDC,fADC,ftimeTDCAmp,fADCAmp);
    pOutStartLoader->UnloadHits();
  } //input streams loop
  
    //load digits    
    pOutStartLoader->LoadDigits("UPDATE");
    TTree *treeD  = pOutStartLoader->TreeD();
    if (treeD == 0x0) {
      pOutStartLoader->MakeTree("D");
      treeD = pOutStartLoader->TreeD();
    }
    AliSTART *fSTART  = (AliSTART*)outRL ->GetAliRun()->GetDetector("START");
    fSTART->MakeBranch("D");
     treeD->Reset();
     treeD->Fill();
  
  pOutStartLoader->WriteDigits("OVERWRITE");
  
  fSTART->ResetDigits();
  pOutStartLoader->UnloadDigits();
    
}


//------------------------------------------------------------------------
Bool_t AliSTARTDigitizer::RegisterPhotoE(Double_t energy)
{

  
  //  Float_t hc=197.326960*1.e6; //mev*nm
  Double_t hc=1.973*1.e-6; //gev*nm
  Float_t lambda=hc/energy;
  Int_t bin=  fEff->GetXaxis()->FindBin(lambda);
  Float_t eff=fEff->GetBinContent(bin);
  Double_t  p = gRandom->Rndm();
  if (p > eff)
    return kFALSE;
  
  return kTRUE;
}
//----------------------------------------------------------------------------
