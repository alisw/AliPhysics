
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


#include "AliSTARTDigitizer.h"
#include "AliSTART.h"
#include "AliSTARThit.h"
#include "AliSTARThitPhoton.h"
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
//  fDebug =0;
  if (GetDebug())
    Info("(AliRunDigitizer* manager)" ,"processed");

  ftimeRightTDC = new TArrayI(12); 
  ftimeLeftTDC = new TArrayI(12); 
  fRightADC = new TArrayI(12); 
  fLeftADC = new TArrayI(12); 
}

//------------------------------------------------------------------------
AliSTARTDigitizer::~AliSTARTDigitizer()
{
// Destructor
  if(GetDebug()) Info("dtor","START"); 
  delete ftimeRightTDC;
  delete ftimeLeftTDC;
  delete fRightADC;
  delete fLeftADC;
}

 //------------------------------------------------------------------------
Bool_t AliSTARTDigitizer::Init()
{
// Initialization
// cout<<"AliSTARTDigitizer::Init"<<endl;
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

#ifdef DEBUG
  cout<<"AliSTARTDigitizer::>Hits2Digits start...\n";
#endif
  //
  // From hits to digits
  //
  Int_t hit, nhits;
  Float_t meanTime;
  Int_t countEr[13],countEl[13];
  Int_t volume,pmt,tr,tl,sumRight;
  Int_t  bestRightADC,bestLeftADC;
  Float_t besttimeleftGaus, besttimerightGaus;
  Float_t timeright[13]={13*0};
  Float_t timeleft[13]={13*0};
  Float_t channelWidth=2.5; //ps
  Int_t channelWidthADC=1; //ps
  //  Int_t thresholdAmpl=10;

  
  AliSTARThit  *startHit;
  TBranch *brHits=0;
  TBranch *brHitPhoton=0;
  pOutStartLoader->LoadDigits("UPDATE");//probably it is necessary to load them before
  fdigits= new AliSTARTdigit();
  pOutStartLoader->GetDigitsDataLoader()->GetBaseLoader(0)->Post(fdigits);

  Int_t nFiles=fManager->GetNinputs();
  for (Int_t inputFile=0; inputFile<nFiles;  inputFile++) {
    if (inputFile < nFiles-1) {
      Warning("Exec", "ignoring input stream %d", inputFile);
      continue;
    }

    Float_t besttimeright=9999.;
    Float_t besttimeleft=9999.;
    Int_t iTimeDiff=0;
    Int_t iTimeAv=0;
    Float_t timeDiff,timeAv; 
    sumRight=0;
    for (Int_t i0=0; i0<13; i0++)
      {
	timeright[i0]=0; timeleft[i0]=0;
	countEr[i0]=0;   countEl[i0]=0;
      }

    inRL = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inputFile));
    if (!inRL->GetAliRun()) inRL->LoadgAlice();
    AliSTART *fSTART  = (AliSTART*) inRL->GetAliRun()->GetDetector("START");
    pInStartLoader = inRL->GetLoader("STARTLoader");
    pInStartLoader->LoadHits("READ");//probably it is necessary to load them before
    TClonesArray *fHits = fSTART->Hits ();

    TTree *th = pInStartLoader->TreeH();
    brHits = th->GetBranch("START");
    brHitPhoton = th->GetBranch("STARThitPhoton");
    if (brHits) {
      fSTART->SetHitsAddressBranch(brHits,brHitPhoton);
    }else{
      cerr<<"EXEC Branch START hit not found"<<endl;
      exit(111);
    } 
    Int_t ntracks    = (Int_t) th->GetEntries();
#ifdef DEBUG
   Info("Digitizer",ntracks);
#endif
     if (ntracks<=0) return;
    // Start loop on tracks in the hits containers
    for (Int_t track=0; track<ntracks;track++) {
      brHits->GetEntry(track);
      nhits = fHits->GetEntriesFast();
      for (hit=0;hit<nhits;hit++) {
	startHit   = (AliSTARThit*) fHits->UncheckedAt(hit);
	if (!startHit) {
	  ::Error("Exec","The unchecked hit doesn't exist");
	  break;
	}
	pmt=startHit->Pmt();
	volume = startHit->Volume();
	if(volume==1){
	  timeright[pmt] = startHit->Time();
	  if(timeright[pmt]<besttimeright)
	    {
	      besttimeright=timeright[pmt];
	  } //timeright
	}//time for right shoulder
	if(volume==2){            
	  timeleft[pmt] = startHit->Time();
	  if(timeleft[pmt]<besttimeleft)
	    {
	      besttimeleft=timeleft[pmt];
	    
	  } //timeleftbest
	}//time for left shoulder
      } //hit loop
    } //track loop
  
    // z position

    //folding with experimental time distribution
    
    Float_t koef=69.7/350.;
    besttimeright=koef*besttimeright;
    besttimeleftGaus=gRandom->Gaus(besttimeleft,0.05);
    bestLeftADC=Int_t (besttimeleftGaus*1000/channelWidth);
    besttimerightGaus=gRandom->Gaus(besttimeright,0.05);
    bestRightADC=Int_t (besttimerightGaus*1000/channelWidth);
    timeDiff=besttimerightGaus-besttimeleftGaus;
#ifdef DEBUG
    cout<<" timediff in ns "<<timeDiff<<" z= "<<timeDiff*30<<endl;
#endif
    meanTime=(besttimerightGaus+besttimeleftGaus)/2.;
    if ( TMath::Abs(timeDiff)<TMath::Abs(0.3) ) 
      {
	Float_t t1=1000.*besttimeleftGaus;
	Float_t t2=1000.*besttimerightGaus;
	t1=t1/channelWidth;   //time in ps to channelWidth
	t2=t2/channelWidth;   //time in ps to channelWidth
	timeAv=(t1+t2)/2.;// time  channel numbres
	
	// Time to TDC signal
	// 256 channels for timediff, range 1ns
	iTimeAv=(Int_t)timeAv; 
	timeDiff= 512+1000*timeDiff/channelWidth; // time  channel numbres 
	iTimeDiff=(Int_t)timeDiff;
	//       fill digits
	fdigits->SetTimeBestLeft(bestLeftADC);
	fdigits->SetTimeBestRight(bestRightADC);
	fdigits->SetMeanTime(iTimeAv);
	fdigits->SetTimeDiff(iTimeDiff);
	for (Int_t i=0; i<12; i++)
	  {
	    //  fill TDC
	    timeright[i+1]=gRandom->Gaus(timeright[i+1],0.05);
	    timeleft[i+1]=gRandom->Gaus(timeleft[i+1],0.05);
	    tr= Int_t (timeright[i+1]*1000/channelWidth); 
	    if(tr<200) tr=0;
	    tl= Int_t (timeleft[i+1]*1000/channelWidth); 
	    if(tl<1000) tl=0;
	    
	    ftimeRightTDC->AddAt(tr,i);
	    ftimeLeftTDC->AddAt(tl,i);
	    //fill ADC
	    Int_t al=( Int_t ) countEl[i+1]/ channelWidthADC;
	    Int_t ar=( Int_t ) countEr[i+1]/ channelWidthADC;
	    fRightADC->AddAt(ar,i);
	    fLeftADC ->AddAt(al,i);
	    sumRight+=countEr[i+1];
	  }
	fdigits->SetTimeRight(*ftimeRightTDC);
	fdigits->SetTimeLeft(*ftimeLeftTDC);
	fdigits->SetADCRight(*fRightADC);
	fdigits->SetADCLeft(*fLeftADC);
	fdigits->SetSumADCRight(sumRight);
      }
    else
      {timeAv=999999; timeDiff=99999;}

    pInStartLoader->UnloadHits();
  } //input streams loop

  pOutStartLoader->WriteDigits("OVERWRITE");
  pOutStartLoader->UnloadDigits();
}


//------------------------------------------------------------------------
Bool_t AliSTARTDigitizer::RegisterPhotoE(/*AliSTARThitPhoton *hit*/)
{
    Double_t    pP = 0.2;    
    Double_t    p;
    
    p = gRandom->Rndm();
    if (p > pP)
      return kFALSE;
    
    return kTRUE;
}
//----------------------------------------------------------------------------
