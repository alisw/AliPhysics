
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
#include <TVector.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TRandom.h>
#include <TArrayI.h>
#include <TH1.h>


#include "AliSTARTDigitizer.h"
#include "AliSTART.h"
#include "AliSTARThit.h"
#include "AliSTARThitPhoton.h"
#include "AliSTARTdigit.h"
#include "AliRunDigitizer.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliRun.h"
#include "AliPDG.h"
#include "AliLoader.h"
#include "AliRunLoader.h"
#include "AliSTARTLoader.h"

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
	cout<<"AliSTARTDigitizer::AliSTARTDigitizer"<<endl;
// ctor which should be used
//  fDebug =0;
 // if (GetDebug()>2)
  //  cerr<<"AliSTARTDigitizer::AliSTARTDigitizer"
   //     <<"(AliRunDigitizer* manager) was processed"<<endl;

}

//------------------------------------------------------------------------
AliSTARTDigitizer::~AliSTARTDigitizer()
{
// Destructor


}

 //------------------------------------------------------------------------
Bool_t AliSTARTDigitizer::Init()
{
// Initialization
 cout<<"AliSTARTDigitizer::Init"<<endl;
 return kTRUE;
}
 

//---------------------------------------------------------------------

void AliSTARTDigitizer::Exec(Option_t* /*option*/)
{


  AliRunLoader *inRL, *outRL;//in and out Run Loaders
  AliLoader *ingime, *outgime;// in and out ITSLoaders

  outRL = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  outgime = outRL->GetLoader("STARTLoader");

#ifdef DEBUG
  cout<<"AliSTARTDigitizer::>SDigits2Digits start...\n";
#endif
  //
  // From hits to digits
  //
  Int_t hit, nhits;
  Int_t CountEr[13],CountEl[13];							//!!!
  Int_t volume,pmt,tr,tl,sumRight;
  Float_t timediff,timeav;
  Float_t besttimeright,besttimeleft,meanTime;
  Int_t  bestRightADC,bestLeftADC;
  Float_t besttimeleftGaus, besttimerightGaus;
  Float_t timeright[13]={13*0};
  Float_t timeleft[13]={13*0};
  Float_t channelWidth=2.5; //ps
  Int_t channelWidthADC=1; //ps
  //  Int_t thresholdAmpl=10;

  ftimeRightTDC = new TArrayI(12); 
  ftimeLeftTDC = new TArrayI(12); 
  fRightADC = new TArrayI(12); 
  fLeftADC = new TArrayI(12); 
  
  inRL = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(0));
  inRL->LoadgAlice();
  
  //  fHits = new TClonesArray ("AliSTARThit", 1000);
  fPhotons = new TClonesArray ("AliSTARThitPhoton", 10000);			//!!!
  AliSTART *START  = (AliSTART*) gAlice->GetDetector("START");
  AliSTARThit  *startHit;
  //use if Cherenkov photons
  //  AliSTARThitPhoton  *startHitPhoton;						//!!!
  TBranch *brHits=0;
  TBranch *brHitPhoton=0;
  fdigits= new AliSTARTdigit();

  Int_t nFiles=fManager->GetNinputs();
  for (Int_t inputFile=0; inputFile<nFiles;  inputFile++) {

    besttimeright=9999.;
    besttimeleft=9999.;
    Int_t timeDiff=0;
    Int_t timeAv=0;
    sumRight=0;
    for (Int_t i0=0; i0<13; i0++)
      {
	timeright[i0]=0; timeleft[i0]=0;
	CountEr[i0]=0;   CountEl[i0]=0;
      }

    inRL = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inputFile));
    ingime = inRL->GetLoader("STARTLoader");
    ingime->LoadHits("READ");//probably it is necessary to load them before
    ingime->LoadDigits("UPDATE");//probably it is necessary to load them before
    //use if Cherenkov photons
    //  TClonesArray *STARThitsPhotons = START->Photons ();
    TClonesArray *fHits = START->Hits ();
    //    cout<<" Load  "<<AliSTARTLoader::LoadDigits()<<endl;

    TTree *th = ingime->TreeH();
    brHits = th->GetBranch("START");
    brHitPhoton = th->GetBranch("STARThitPhoton");
    if (brHits) {
      START->SetHitsAddressBranch(brHits,brHitPhoton);
    }else{
      cerr<<"EXEC Branch START hit not found"<<endl;
      exit(111);
    } 
    Int_t ntracks    = (Int_t) th->GetEntries();
    cout<<" ntracks "<<ntracks<<endl;
    if (ntracks<=0) return;
    // Start loop on tracks in the photon hits containers 
    // for amplitude
    /*
    if(brHitPhoton) {
      cout<<"brHitPhoton "<<endl; 
      for (Int_t track=0; track<ntracks;track++) {
	brHitPhoton -> GetEntry(track);;
	nhits = STARThitsPhotons->GetEntriesFast();
	for (hit=0;hit<nhits;hit++) {
	  startHitPhoton   = (AliSTARThitPhoton*) 
	   STARThitsPhotons ->UncheckedAt(hit);
	  pmt=startHitPhoton->fPmt;
	  volume = startHitPhoton->fArray;
	  if(RegisterPhotoE(startHitPhoton)) 
	    {
	      if (volume == 1) CountEr[pmt]++;
	      if (volume == 2) CountEl[pmt]++;
	    }
	} //hit photons
      } //track photons
    } // was photons
    */
    // Start loop on tracks in the hits containers
    for (Int_t track=0; track<ntracks;track++) {
      brHits->GetEntry(track);
      nhits = fHits->GetEntriesFast();
      //  cout<<" brHits hits "<<nhits<<endl;
      for (hit=0;hit<nhits;hit++) {
	startHit   = (AliSTARThit*) fHits->UncheckedAt(hit);
	pmt=startHit->Pmt();
	volume = startHit->Volume();
	if(volume==1){
	  timeright[pmt] = startHit->Time();
	  if(timeright[pmt]<besttimeright)
	    //&&CountEr[pmt-1]>thresholdAmpl)
	    {
	    besttimeright=timeright[pmt];
	  } //timeright
	}//time for right shoulder
	if(volume==2){            
	  timeleft[pmt] = startHit->Time();
	  if(timeleft[pmt]<besttimeleft)
	    //&&CountEl[pmt-1]>thresholdAmpl) 
	    {
	    besttimeleft=timeleft[pmt];
	    
	  } //timeleftbest
	}//time for left shoulder
      } //hit loop
    } //track loop
  
    // z position
    cout<<" right time  "<<besttimeright<<
      " right distance "<<besttimeright*30<<endl;;
    cout<<" left time  "<<besttimeleft<<
      " left distance "<<besttimeleft*30<<endl;;
  

    //folding with experimental time distribution
    
    besttimeleftGaus=gRandom->Gaus(besttimeright,0.05);
    cout<<" besttimeleftGaus "<<besttimeleftGaus<<endl;
    bestLeftADC=Int_t (besttimeleftGaus*1000/channelWidth);
    Float_t koef=69.7/350.;
    besttimeright=koef*besttimeleft;
    besttimerightGaus=gRandom->Gaus(besttimeleft,0.05);
    
    bestRightADC=Int_t (besttimerightGaus*1000/channelWidth);
    timediff=besttimerightGaus-besttimeleftGaus;
    cout<<" timediff in ns "<<timediff<<" z= "<<timediff*30<<endl;
    meanTime=(besttimerightGaus+besttimeleftGaus)/2.;
    if ( TMath::Abs(timediff)<TMath::Abs(0.3) ) 
      {
	Float_t t1=1000.*besttimeleftGaus;
	Float_t t2=1000.*besttimerightGaus;
	t1=t1/channelWidth;   //time in ps to channelWidth
	t2=t2/channelWidth;   //time in ps to channelWidth
	timeav=(t1+t2)/2.;
	
	// Time to TDC signal
	// 256 channels for timediff, range 1ns
	
	timediff=512+1000*timediff/channelWidth; // time in ps 
	
	timeAv = (Int_t)(timeav);   // time  channel numbres
	timeDiff = (Int_t)(timediff); // time  channel numbres
	//       fill digits
	fdigits->SetTimeBestLeft(bestLeftADC);
	fdigits->SetTimeBestRight(bestRightADC);
	fdigits->SetMeanTime(timeAv);
	fdigits->SetTimeDiff(timeDiff);
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
	    Int_t al=( Int_t ) CountEl[i+1]/ channelWidthADC;
	    Int_t ar=( Int_t ) CountEr[i+1]/ channelWidthADC;
	    fRightADC->AddAt(ar,i);
	    fLeftADC ->AddAt(al,i);
	    sumRight+=CountEr[i+1];
	  }
	fdigits->SetTimeRight(*ftimeRightTDC);
	fdigits->SetTimeLeft(*ftimeLeftTDC);
	fdigits->SetADCRight(*fRightADC);
	fdigits->SetADCLeft(*fLeftADC);
	// cout<<" before sum"<<endl;
	fdigits->SetSumADCRight(sumRight);
      }
    else
      {timeAv=999999; timeDiff=99999;}

// trick to find out output dir:


/*
   // trick to find out output dir:
    TTree *outTree = fManager->GetTreeD();
    if (!outTree) {
      cerr<<"something wrong with output...."<<endl;
      exit(111);
    }

    Char_t nameDigits[20];
    TDirectory *wd = gDirectory;
    outTree->GetDirectory()->cd();
    fdigits->Write(nameDigits);
    cout<<nameDigits<<endl;
    wd->cd();
*/  

     Char_t nameDigits[20];
    sprintf(nameDigits,"START_D_%d",fManager->GetOutputEventNr());
    fdigits->Write(nameDigits);

    //    outgime->WriteDigits("OVERWRITE");
  }
}


//------------------------------------------------------------------------
Bool_t AliSTARTDigitizer::RegisterPhotoE(/*AliSTARThitPhoton *hit*/)
{
    Double_t    P = 0.2;    
    Double_t    p;
    
    p = gRandom->Rndm();
    if (p > P)
      return kFALSE;
    
    return kTRUE;
}
//----------------------------------------------------------------------------
