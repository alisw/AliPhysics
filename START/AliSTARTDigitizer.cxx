
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


#include "AliSTARTDigitizer.h"
#include "AliSTART.h"
#include "AliSTARThit.h"
#include "AliSTARTdigit.h"
#include "AliRunDigitizer.h"

#include "AliRun.h"
#include "AliPDG.h"

#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>

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

void AliSTARTDigitizer::Exec(Option_t* option)
{



#ifdef DEBUG
  cout<<"AliSTARTDigitizer::>SDigits2Digits start...\n";
#endif
  //
  // From hits to digits
  //
  Int_t hit;
  Int_t nhits;
  Int_t volume,pmt;
  char nameDigits[20];
  Float_t timediff,timeright,timeleft,timeav;
  Float_t besttimeright,besttimeleft,meanTime;
  Int_t channelWidth=10;
  fHits = new TClonesArray ("AliSTARThit", 1000);
  AliSTART *START  = (AliSTART*) gAlice->GetDetector("START");

  AliSTARThit  *startHit;
  TBranch *brHits=0;
  fdigits= new AliSTARTdigit();

  Int_t nFiles=fManager->GetNinputs();
  for (Int_t inputFile=0; inputFile<nFiles;  inputFile++) {
    sprintf(nameDigits,"START_D_%d",fManager->GetOutputEventNr());
 
    besttimeright=9999.;
    besttimeleft=9999.;
    Int_t timeDiff=0;
    Int_t timeAv=0;
    TClonesArray *STARThits = START->Hits ();

   TTree *th = fManager->GetInputTreeH(inputFile);
    brHits = th->GetBranch("START");
    if (brHits) {
      START->SetHitsAddressBranch(brHits);
    }else{
      cerr<<"EXEC Branch START hit not found"<<endl;
      exit(111);
    } 
    Int_t ntracks    = (Int_t) th->GetEntries();
    if (ntracks<=0) return;
    // Start loop on tracks in the hits containers
    for (Int_t track=0; track<ntracks;track++) {
      brHits->GetEntry(track);
      nhits = STARThits->GetEntriesFast();
      for (hit=0;hit<nhits;hit++) {
	startHit   = (AliSTARThit*) STARThits->UncheckedAt(hit);
	pmt=startHit->fPmt;
	volume = startHit->fVolume;
	if(volume==1){
	  timeright = startHit->fTime;
	  if(timeright<besttimeright) {
	    besttimeright=timeright;
	  } //timeright
	}//time for right shoulder
	if(volume==2){            
	  timeleft = startHit->fTime;
	  if(timeleft<besttimeleft) {
	    besttimeleft=timeleft;
	  } //timeleftbest
	}//time for left shoulder
      } //hit loop
    } //track loop
 
    //folding with experimental time distribution
    Float_t besttimerightGaus=gRandom->Gaus(besttimeright,0.05);
    Float_t koef=69.7/350.;
    besttimeleft=koef*besttimeleft;
    Float_t besttimeleftGaus=gRandom->Gaus(besttimeleft,0.05);
    timediff=besttimerightGaus-besttimeleftGaus;
    meanTime=(besttimerightGaus+besttimeleftGaus)/2.;
    if ( TMath::Abs(timediff)<TMath::Abs(3.) && meanTime<TMath::Abs(5.)) 
     {
       //we assume centre of bunch is 5ns after TTS signal
       //TOF values are relative of the end of bunch
       Float_t ppBunch=25;
    
       ppBunch=ppBunch-10/2;
       Float_t t1=1000.*besttimeleftGaus;
       Float_t t2=1000.*besttimerightGaus;
       t1=t1/channelWidth+ppBunch; //time in ps to channelWidth
       t2=t2/channelWidth+ppBunch; //time in ps to channelWidth
       timeav=(t1+t2)/2.;
     
       // Time to TDC signal
       // 256 channels for timediff, range 1ns
       
       timediff=128+1000*timediff/channelWidth; // time in ps 

       timeAv = (Int_t)(timeav);   // time (ps) channel numbres
       timeDiff = (Int_t)(timediff); // time ( ps) channel numbres
       fdigits->Set(timeAv,timeDiff);
       fdigits->Print();
     }
    else
      {timeAv=999999; timeDiff=99999;}

// trick to find out output dir:
    TTree *outTree = fManager->GetTreeD();
    if (!outTree) {
      cerr<<"something wrong with output...."<<endl;
      exit(111);
    }
    TDirectory *wd = gDirectory;
    outTree->GetDirectory()->cd();
    fdigits->Write(nameDigits);
    wd->cd();
  }

}


