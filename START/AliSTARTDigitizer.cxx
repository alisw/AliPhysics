
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
#include <TGraph.h>

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
#include "AliSTARTParameters.h"
#include "AliCDBLocal.h"
#include "AliCDBStorage.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

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
   ftimeCFD(0),
   ftimeLED(0),
   fADC(0),
   fADC0(0)
{
// ctor which should be used

  AliDebug(1,"processed");

  fSTART = 0;
  fHits = 0;
  fdigits = 0;

  ftimeCFD = new TArrayI(24); 
  fADC = new TArrayI(24); 
  ftimeLED = new TArrayI(24); 
  fADC0 = new TArrayI(24); 
  

}

//------------------------------------------------------------------------
AliSTARTDigitizer::~AliSTARTDigitizer()
{
// Destructor

  AliDebug(1,"START");

  delete ftimeCFD;
  delete fADC;
  delete ftimeLED;
  delete  fADC0;
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
  Int_t countE[24];
  Int_t volume, pmt, trCFD, trLED; 
  Float_t sl, qt;
  Int_t  bestRightTDC, bestLeftTDC;
  Float_t time[24], besttime[24], timeGaus[24] ;
  //  Float_t channelWidth=25.; //ps
    //Q->T-> coefficients !!!! should be asked!!!
  Float_t gain[24],timeDelayCFD[24], timeDelayLED[24];
  Int_t threshold =50; //photoelectrons
  Int_t qtThreshold=13, qtCh;
  Float_t zdetA, zdetC;
  TObjArray slewingLED;
  AliSTARTParameters* param = AliSTARTParameters::Instance();
  Int_t ph2Mip = param->GetPh2Mip();     
  Int_t channelWidth = param->GetChannelWidth() ;  
  
  param->Init();
  for (Int_t i=0; i<24; i++){
    timeDelayCFD[i] = param->GetTimeDelayCFD(i);
    timeDelayLED[i] = param->GetTimeDelayLED(i);
    gain[i] = param->GetGain(i);
    TGraph* gr = param ->GetSlew(i);
    slewingLED.AddAtAndExpand(gr,i);
    TGraph* grEff = param ->GetPMTeff(i);
    fEffPMT.AddAtAndExpand(grEff,i);
  }
  zdetC = param->GetZposition(0);
  zdetA = param->GetZposition(1);
  
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
    Int_t pmtBestRight=9999;
    Int_t pmtBestLeft=9999;
    Int_t timeDiff=99999, meanTime=0;
    Int_t sumMult =0;
    ftimeCFD -> Reset();
    fADC -> Reset();
    fADC0 -> Reset();
    ftimeLED ->Reset();
    cout<<" before reading event "<<endl;
    for (Int_t i0=0; i0<24; i0++)
      {
	time[i0]=besttime[i0]=timeGaus[i0]=99999; countE[i0]=0;
      }
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
	  
	  if(e>0 && RegisterPhotoE(numpmt,e)) {
	    countE[numpmt]++;
	    besttime[numpmt] = startHit->Time();
	    if(besttime[numpmt]<time[numpmt])
	      {
		time[numpmt]=besttime[numpmt];
	      }
	  } //photoelectron accept 
	} //hits loop
    } //track loop
    cout<<" hits was read "<<endl;
    
    //spread time right&left by 25ps   && besttime
     Float_t c = 0.0299792; // cm/ps

   Float_t koef=(zdetA-zdetC)/c; //correction position difference by cable
    for (Int_t ipmt=0; ipmt<12; ipmt++){
      if(countE[ipmt] > threshold) {
	timeGaus[ipmt]=gRandom->Gaus(time[ipmt],25)+koef;
	if(timeGaus[ipmt]<besttimeleft){
	  besttimeleft=timeGaus[ipmt]; //timeleft
	  pmtBestLeft=ipmt;}
      }
    }
    cout<<"left "<<besttimeleft<<" "<<pmtBestLeft<<" koef "<<koef<<endl;
    for ( Int_t ipmt=12; ipmt<24; ipmt++){
      if(countE[ipmt] > threshold) {
	timeGaus[ipmt]=gRandom->Gaus(time[ipmt],25); 
	if(timeGaus[ipmt]<besttimeright) {
	  besttimeright=timeGaus[ipmt]; //timeright
	  pmtBestRight=ipmt;}
      }	
    }
    cout<<"right"<<besttimeright<<" "<<pmtBestRight<<endl;
  
   //folding with alignmentz position distribution  

    //    Float_t koef=(zdetA-zdetC)/c;
    //  Float_t  besttimeleftR= besttimeleft;
    // besttimeleft=koef+besttimeleft;
    cout<<"  besttimeleft "<< besttimeleft<<" delay "<< timeDelayCFD[pmtBestLeft]<<" channelWidth "<<channelWidth<<endl;
    bestLeftTDC=Int_t ((besttimeleft+1000*timeDelayCFD[pmtBestLeft])
		       /channelWidth);
    cout<<" bestLeftTDC "<<bestLeftTDC<<" delay "<<timeDelayCFD[pmtBestLeft]<<endl;
    bestRightTDC=Int_t ((besttimeright+1000*timeDelayCFD[pmtBestLeft])
			/channelWidth);
    timeDiff=Int_t (((besttimeleft+1000*timeDelayCFD[pmtBestLeft])-
		    (besttimeright+1000*timeDelayCFD[pmtBestLeft])
		    )/channelWidth);
    meanTime=Int_t (((besttimeright+1000*timeDelayCFD[pmtBestLeft]+
		      besttimeleft+1000*timeDelayCFD[pmtBestLeft])/2.)
		    /channelWidth);
    cout<<"bestLeftTDC"<<bestLeftTDC<<" bestRightTDC "<< bestRightTDC<<
      " timeDiff "<<timeDiff<<" meanTime= "<<meanTime<<endl;   
    AliDebug(2,Form(" time in ns %f ", besttimeleft));
    for (Int_t i=0; i<24; i++)
      {
       	Float_t  al = countE[i]; 
	if (al>threshold && timeGaus[i]<50000 ) {
	//fill ADC
	// QTC procedure:
	// phe -> mV 0.3; 1MIP ->500phe -> ln (amp (mV)) = 5;
	// max 200ns, HIJING  mean 50000phe -> 15000mv -> ln = 15 (s zapasom)
	// channel 25ps
	  qt= 50.*al*gain[i]/ph2Mip;  // 50mv/Mip amp in mV 
	//  fill TDC
	  trCFD = Int_t (timeGaus[i] + 1000.*timeDelayCFD[i])/channelWidth; 
	  trLED= Int_t (timeGaus[i] + 1000.*timeDelayLED[i]); 
	  sl = ((TGraph*)slewingLED.At(i))->Eval(qt);
	  trLED = Int_t(( trLED + 1000*sl )/channelWidth);
	  qtCh=Int_t (1000.*TMath::Log(qt)) / channelWidth;
	  fADC0->AddAt(0,i);
	  fADC->AddAt(qtCh,i);
	  ftimeCFD->AddAt(Int_t (trCFD),i);
	  ftimeLED->AddAt(trLED,i); 
	  sumMult += Int_t ((al*gain[i]/ph2Mip)*50) ;
	  cout<<al<<" "<<qt<<" trCFD "<<trCFD<<" sl "<<sl<<" trLED "<<trLED<<
	    " qtThreshold "<<qtThreshold<<"  qtCh "<< qtCh<<endl;
	}
      } //pmt loop

	if (sumMult > threshold){
	  fSumMult =  Int_t (TMath::Log(sumMult));
	}
    fSTART->AddDigit(bestRightTDC,bestLeftTDC,meanTime,timeDiff,fSumMult,
		     ftimeCFD,fADC,ftimeLED,fADC0);
    pOutStartLoader->UnloadHits();
  } //input streams loop
  
    //load digits    
    pOutStartLoader->LoadDigits("UPDATE");
    TTree *treeD  = pOutStartLoader->TreeD();
    if (treeD == 0x0) {
      pOutStartLoader->MakeTree("D");
      treeD = pOutStartLoader->TreeD();
    }
    treeD->Reset();
    fSTART  = (AliSTART*)outRL ->GetAliRun()->GetDetector("START");
    // Make a branch in the tree 
    fSTART->MakeBranch("D");
    //     treeD->Print();
     treeD->Fill();
     cout<<"   treeD->Fill(); "<<endl;
  
     pOutStartLoader->WriteDigits("OVERWRITE");
     
     cout<<"  pOutStartLoader->WriteDigits "<<endl;
     fSTART->ResetDigits();
     pOutStartLoader->UnloadDigits();
     
}


//------------------------------------------------------------------------
Bool_t AliSTARTDigitizer::RegisterPhotoE(Int_t ipmt,Double_t energy)
{

  
  //  Float_t hc=197.326960*1.e6; //mev*nm
  Double_t hc=1.973*1.e-6; //gev*nm
  Float_t lambda=hc/energy;
  Float_t eff = ((TGraph*) fEffPMT.At(ipmt))->Eval(lambda);
  Double_t  p = gRandom->Rndm();

  if (p > eff)
    return kFALSE;
  
  return kTRUE;
}

//----------------------------------------------------------------------------
