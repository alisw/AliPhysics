
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
    :AliDigitizer(manager) 
{
  //	cout<<"AliSTARTDigitizer::AliSTARTDigitizer"<<endl;
// ctor which should be used

  AliDebug(1,"processed");

  ftimeTDC = new TArrayI(24); 
  fADC = new TArrayI(24); 
}

//------------------------------------------------------------------------
AliSTARTDigitizer::~AliSTARTDigitizer()
{
// Destructor

  AliDebug(1,"START"); 
  delete ftimeTDC;
  delete fADC;
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

  AliRunLoader *inRL, *outRL;//in and out Run Loaders
  AliLoader *pInStartLoader, *pOutStartLoader;// in and out STARTLoaders

  outRL = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  pOutStartLoader = outRL->GetLoader("STARTLoader");

  AliDebug(1,"start...");

  //
  // From hits to digits
  //
  Int_t hit, nhits;
  Float_t meanTime;
  Int_t countE[24];
  Int_t volume,pmt,tr,sumRight;
  Int_t  bestRightTDC,bestLeftTDC;
  Float_t time[24]={24*0};
  Float_t besttime[24]={24*0};
  Float_t timeGaus[37]={24*0};
  Float_t channelWidth=25.; //ps
  
  AliSTARThit  *startHit;
  TBranch *brHits=0;

  pOutStartLoader->LoadDigits("UPDATE");//probably it is necessary to load them before
  fdigits= new AliSTARTdigit();
  pOutStartLoader->GetDigitsDataLoader()->GetBaseLoader(0)->Post(fdigits);

  Int_t nFiles=fManager->GetNinputs();
  for (Int_t inputFile=0; inputFile<nFiles;  inputFile++) {
    if (inputFile < nFiles-1) {
      AliWarning(Form("ignoring input stream %d", inputFile));
      continue;

    }

    Float_t besttimeright=9999.;
    Float_t besttimeleft=9999.;
    Float_t timeDiff;
    sumRight=0;
    for (Int_t i0=0; i0<24; i0++)
      {
	time[i0]=9999;  besttime[i0]=9999;	countE[i0]=0;
      }

    inRL = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inputFile));
    if (!inRL->GetAliRun()) inRL->LoadgAlice();
    AliSTART *fSTART  = (AliSTART*) inRL->GetAliRun()->GetDetector("START");
    pInStartLoader = inRL->GetLoader("STARTLoader");
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
	  Float_t e=startHit->Etot();
	  volume = startHit->Volume();
	  if(RegisterPhotoE(e)) countE[numpmt]++;
	  besttime[numpmt] = startHit->Time();
	  if(besttime[numpmt]<time[numpmt])
	    {
	      time[numpmt]=besttime[numpmt];
	    }
	  
	} //hits loop
    } //track loop
    
    //best time right&left   
    for (Int_t ipmt=0; ipmt<12; ipmt++)
      {
	timeGaus[ipmt]=gRandom->Gaus(time[ipmt],0.025);
	if(timeGaus[ipmt]<besttimeleft) besttimeleft=timeGaus[ipmt]; //timeleft
      }
    for ( Int_t ipmt=12; ipmt<24; ipmt++)
      {
   	timeGaus[ipmt]=gRandom->Gaus(time[ipmt],0.025);
	if(timeGaus[ipmt]<besttimeright)  besttimeright=timeGaus[ipmt]; //timeright
      }// besttime
	
    //folding with experimental time distribution
      Float_t c = 29.9792; // mm/ns
    Float_t koef=(350.-69.7)/c;
    Float_t  besttimeleftR= besttimeleft;
    besttimeleft=koef+besttimeleft;
    bestLeftTDC=Int_t (besttimeleftR*1000/channelWidth);
    bestRightTDC=Int_t (besttimeright*1000/channelWidth);
    timeDiff=(c*(besttimeright-besttimeleftR)-(350.-69.7))/(2*c);
    meanTime=(besttimeright+besttimeleftR)/2.;
    Float_t ds=(c*(besttimeright-besttimeleftR)-(350.-69.7))/2;
    AliDebug(2,Form(" timediff in ns %f  z= %f real point%f",timeDiff,timeDiff*c,ds));

 
    // Time to TDC signal
    Int_t iTimeAv=Int_t (meanTime*1000/channelWidth); 
    // time  channel numbres 
    //       fill digits
    fdigits->SetTimeBestLeft(bestLeftTDC);
    fdigits->SetTimeBestRight(bestRightTDC);
    fdigits->SetMeanTime(iTimeAv);
    
    //ADC features
  
    for (Int_t i=0; i<24; i++)
      {

	//  fill TDC
	tr= Int_t (timeGaus[i]*1000/channelWidth); 
	if(timeGaus[i]>100) tr=0;
	ftimeTDC->AddAt(tr,i);
	
	//fill ADC
       	Int_t al= countE[i]; //1024 - channel width shoul be change
 	fADC->AddAt(al,i);
      } //pmt loop

    fdigits->SetTime(*ftimeTDC);
    fdigits->SetADC(*fADC);
     
    pInStartLoader->UnloadHits();
    
  } //input streams loop
     
  pOutStartLoader->WriteDigits("OVERWRITE");
  pOutStartLoader->UnloadDigits();
}


//------------------------------------------------------------------------
Bool_t AliSTARTDigitizer::RegisterPhotoE(Float_t e)
{

  
  //  Float_t hc=197.326960*1.e6; //mev*nm
  Float_t hc=1.973*1.e-6; //gev*nm
  Float_t lambda=hc/e;
  Char_t filename[80];
  TH1F *hEff;
  Char_t *dirname=getenv("ALICE_ROOT");
  sprintf(filename,"%s/START/PMTefficiency.root",dirname);
  TFile *file= new TFile(filename);
  file->cd();
  hEff =  (TH1F*)file->Get("hEff");
  Int_t bin=  hEff->GetXaxis()->FindBin(lambda);
  Float_t eff=hEff->GetBinContent(bin);
  file->Close();
  Double_t  p = gRandom->Rndm();
  if (p > eff)
    return kFALSE;
  
  return kTRUE;
}
//----------------------------------------------------------------------------
