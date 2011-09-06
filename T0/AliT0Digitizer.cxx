
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

/* $Id$ */

/******************************************************************
 *    Produde digits from hits
 *       digits is TObject and includes
 *	We are writing array if C & A  TDC
 *	C & A  ADC (will need for slow simulation)
 *	TOF first particle C & A
 *	mean time and time difference (vertex position)
 *
 *      Alla.Maevskaya@cern.ch 
 ****************************************************************/


#include <TArrayI.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>
#include <TRandom.h>
#include <TTree.h> 

#include "AliLog.h"
#include "AliT0Digitizer.h"
#include "AliT0.h"
#include "AliT0hit.h"
#include "AliT0digit.h"
#include "AliRunDigitizer.h"
#include "AliRun.h"
#include <AliLoader.h>
#include <AliRunLoader.h>
#include <stdlib.h>
#include "AliT0Parameters.h"

ClassImp(AliT0Digitizer)

//___________________________________________
  AliT0Digitizer::AliT0Digitizer()  :AliDigitizer(),
				     fT0(0),
				     fHits(0),
				     fdigits(0),
				     ftimeCFD(new TArrayI(24)), 
				     ftimeLED (new TArrayI(24)), 
				     fADC(new TArrayI(24)), 
				     fADC0 (new TArrayI(24)),
				     fSumMult(0),
				     fAmpLED(0),
                                     fAmpQTC(0),
				     fParam(0)
{
// Default ctor - don't use it
  ;
}

//___________________________________________
AliT0Digitizer::AliT0Digitizer(AliRunDigitizer* manager) 
  :AliDigitizer(manager),
   fT0(0),
   fHits(0),
   fdigits(0),
   ftimeCFD(new TArrayI(24)), 
   ftimeLED (new TArrayI(24)), 
   fADC(new TArrayI(24)), 
   fADC0 (new TArrayI(24)),
   fSumMult(0),
   fAmpLED(0),
   fAmpQTC(0),
   fParam(0)
{
// ctor which should be used

  AliDebug(1,"processed");
  fParam = AliT0Parameters::Instance();
  fParam->Init();

  for (Int_t i=0; i<24; i++){
    TGraph* gr = fParam ->GetAmpLED(i);
     if(gr) {
      Int_t np = gr->GetN();
      Double_t *x = gr->GetX();
      Double_t *y = gr->GetY();
      
      Double_t *x1 = new Double_t[np];
      Double_t *y1 = new Double_t[np];
      for (Int_t ii=0; ii<np; ii++) {
	y1[ii]=y[np-ii-1]; 
	x1[ii]=x[np-ii-1];
      }
    
    TGraph *grInverse = new TGraph(np,y1,x1);
    fAmpLED.AddAtAndExpand(grInverse,i);
    }
  }
  for (Int_t i=0; i<24; i++){
    TGraph* grq = fParam ->GetQTC(i);
    if(grq){
      Int_t npq = grq->GetN();
      Double_t *xq = grq->GetX();
      Double_t *yq = grq->GetY();
      Double_t *x1q = new Double_t[npq];
      Double_t *y1q = new Double_t[npq];
      for (Int_t ii=1; ii<npq; ii++) {
	y1q[ii]=yq[ii-1]; 
	x1q[ii]=xq[ii-1];
      }
      TGraph *grInverseQTC = new TGraph(npq,y1q,x1q);
      fAmpQTC.AddAtAndExpand(grInverseQTC,i);
    }
  }
}


//------------------------------------------------------------------------
AliT0Digitizer::~AliT0Digitizer()
{
// Destructor

  AliDebug(1,"T0");

  delete ftimeCFD;
  delete fADC;
  delete ftimeLED;
  delete  fADC0;
  //  delete fAmpLED;
}

//------------------------------------------------------------------------
Bool_t AliT0Digitizer::Init()
{
// Initialization
  AliDebug(1," Init");
 return kTRUE;
}
 
//---------------------------------------------------------------------
void AliT0Digitizer::Exec(Option_t* /*option*/)
{

  /*
    Produde digits from hits
        digits is TObject and includes
	We are writing array if C & A  TDC
	C & A  ADC (will need for slow simulation)
	TOF first particle C & A
	mean time and time difference (vertex position)
	
  */

  //output loader 
  AliRunLoader *outRL = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  AliLoader * pOutStartLoader = outRL->GetLoader("T0Loader");

  AliDebug(1,"start...");
  //input loader
  //
  // From hits to digits
  //
  Int_t hit, nhits;
  Int_t countE[24];
  Int_t volume, pmt, trCFD, trLED; 
  Float_t sl=0, qt=0;
  Int_t  bestATDC=0;
  Int_t  bestCTDC=0;
  Float_t qtCh=0;
  Float_t time[24], besttime[24], timeGaus[24] ;
  //Q->T-> coefficients !!!! should be asked!!!
  Float_t timeDelayCFD[24];
  Int_t threshold =50; //photoelectrons
  Float_t zdetA, zdetC;
  Int_t sumMultCoeff = 100;
  Int_t refpoint=0;
   
  Int_t ph2Mip = fParam->GetPh2Mip();     
  Float_t channelWidth = fParam->GetChannelWidth() ;  
  Float_t delayVertex = fParam->GetTimeDelayTVD();
   

  zdetC = TMath::Abs(fParam->GetZPosition("T0/C/PMT1"));
  zdetA  = TMath::Abs(fParam->GetZPosition("T0/A/PMT15"));
  
  //  printf(" !!!!!Z det A = %f C = % f",zdetA,zdetC);
  AliT0hit  *startHit;
  TBranch *brHits=0;
  
  Int_t nFiles=fManager->GetNinputs();
  for (Int_t inputFile=0; inputFile<nFiles;  inputFile++) {
    if (inputFile < nFiles-1) {
      AliWarning(Form("ignoring input stream %d", inputFile));
      continue;
      
    }
    Float_t slew = 0;
    Float_t besttimeC=99999.;
    Float_t besttimeA=99999.;
    Int_t pmtBestC=99999;
    Int_t pmtBestA=99999;
    Int_t timeDiff=99999, meanTime=99999;
    Int_t sumMult =0;   fSumMult=0;
    bestATDC = 99999;  bestCTDC = 99999;
 

    ftimeCFD -> Reset();
    fADC -> Reset();
    fADC0 -> Reset();
    ftimeLED ->Reset();
    for (Int_t i0=0; i0<24; i0++)
      {
	time[i0]=besttime[i0]=timeGaus[i0]=999999; countE[i0]=0;
      }
    AliRunLoader * inRL = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(inputFile));
    AliLoader * pInStartLoader = inRL->GetLoader("T0Loader");
    if (!inRL->GetAliRun()) inRL->LoadgAlice();
    fT0  = (AliT0*)inRL ->GetAliRun()->GetDetector("T0");

       //read Hits 
    pInStartLoader->LoadHits("READ");//probably it is necessary to load them before
    fHits = fT0->Hits ();
    TTree *th = pInStartLoader->TreeH();
    brHits = th->GetBranch("T0");
    if (brHits) {
      fT0->SetHitsAddressBranch(brHits);
    }else{
      AliWarning("Branch T0 hit not found for this event");
      // fT0->AddDigit(bestATDC,bestCTDC,meanTime,timeDiff,fSumMult, refpoint,
      //	       ftimeCFD,fADC0,ftimeLED,fADC);
      continue;      
    } 
    Int_t ntracks    = (Int_t) th->GetEntries();
    
    if (ntracks<=0) return;
    // Start loop on tracks in the hits containers
    for (Int_t track=0; track<ntracks;track++) {
      brHits->GetEntry(track);
      nhits = fHits->GetEntriesFast();
      for (hit=0;hit<nhits;hit++) 
	{
	  startHit   = (AliT0hit*) fHits->UncheckedAt(hit);
	  if (!startHit) {
 	    AliError("The unchecked hit doesn't exist");
	    break;
	  }
	  pmt=startHit->Pmt();
	  Int_t numpmt=pmt-1;
	  Double_t e=startHit->Etot();
	  volume = startHit->Volume();
	  
	  if(e>0 ) {
	    countE[numpmt]++;
	    besttime[numpmt] = startHit->Time();
	    if(besttime[numpmt]<time[numpmt])
	      {
		time[numpmt]=besttime[numpmt];
	      }
	  } //photoelectron accept 
	} //hits loop
    } //track loop
    
    //spread time A&C by 25ps   && besttime
    Float_t c = 0.0299792; // cm/ps
    
    Float_t koef=(zdetA-zdetC)/c; //correction position difference by cable
    for (Int_t ipmt=0; ipmt<12; ipmt++){
      if(countE[ipmt] > threshold) {
	timeGaus[ipmt]=gRandom->Gaus(time[ipmt],25)+koef;
	if(timeGaus[ipmt]<besttimeC){
	  besttimeC=timeGaus[ipmt]; //timeC
	  pmtBestC=ipmt;}
      }
    }
    for ( Int_t ipmt=12; ipmt<24; ipmt++){
      if(countE[ipmt] > threshold) {
	timeGaus[ipmt]=gRandom->Gaus(time[ipmt],25); 
	if(timeGaus[ipmt]<besttimeA) {
	  besttimeA=timeGaus[ipmt]; //timeA
	  pmtBestA=ipmt;}
      }	
    }

    timeDelayCFD[0] = fParam->GetTimeDelayCFD(0);
 
    for (Int_t i=0; i<24; i++)
      {
       	Float_t  al = countE[i]; 
	if (al>threshold && timeGaus[i]<50000 ) {
	  //fill ADC
	  // QTC procedure:
	  // phe -> mV 0.3; 1MIP ->500phe -> ln (amp (mV)) = 5;
	  // max 200ns, HIJING  mean 50000phe -> 15000mv -> ln = 15 (s zapasom)
	  // channel 25ps
	  qt= al/ph2Mip;  // 50mv/Mip amp in mV 
	  // before will we have calibration for high multiplicity 
	  //	  if (qt > 115.) qt =115.; //must be fix!!!
	  //  fill TDC
	  timeDelayCFD[i] = fParam->GetTimeDelayCFD(i);
 	  trCFD = Int_t (timeGaus[i]/channelWidth + timeDelayCFD[i]); 

	  TGraph* gr = ((TGraph*)fAmpLED.At(i));
	  if(gr) sl = gr->Eval(qt);

	  TGraph* gr1 = ((TGraph*)fAmpQTC.At(i));
	  if(gr1)  qtCh = gr1->Eval(qt);
	  fADC0->AddAt(0,i);
	  if(qtCh)
	    fADC->AddAt(Int_t(qtCh),i);
	  //	  sumMult += Int_t ((al*gain[i]/ph2Mip)*50) ;
	  sumMult += Int_t (qtCh/sumMultCoeff)  ;
	 
	  // put slewing 
	  TGraph *fu=(TGraph*) fParam ->GetWalk(i);
	  if(fu)  slew=fu->Eval(Float_t(qtCh));
	  
	  //	  trCFD=trCFD-Int_t(fMaxValue[i]-slew);
	  trCFD = trCFD + Int_t(slew); //for the same channel as cosmic
	  ftimeCFD->AddAt(Int_t (trCFD),i);
	  trLED = Int_t(trCFD  + sl );
	  ftimeLED->AddAt(trLED,i); 
	  AliDebug(1,Form("  pmt %i : delay %f time in ns %f time in channels %i  LEd %i  ",  i, timeDelayCFD[i], timeGaus[i],trCFD, trLED ));
	  AliDebug(1,Form(" qt in MIP %f led-cfd in  %f qt in channels %f   ",qt, sl, qtCh));

	}
      } //pmt loop

    //folding with alignmentz position distribution  
    if( besttimeC > 10000. && besttimeC <15000)
      bestCTDC=Int_t ((besttimeC+timeDelayCFD[pmtBestC])
			 /channelWidth);
 
    if( besttimeA > 10000. && besttimeA <15000)
      bestATDC=Int_t ((besttimeA+timeDelayCFD[pmtBestA])
			/channelWidth);

    if (bestATDC < 99999 && bestCTDC < 99999)
      {
	timeDiff=Int_t (((besttimeC-besttimeA)+1000*delayVertex)
			/channelWidth);
	meanTime=Int_t (((besttimeC+besttimeA)/2. )/channelWidth);
      }

    if (sumMult > threshold){
      fSumMult =  Int_t (1000.* TMath::Log(Double_t(sumMult) / Double_t(sumMultCoeff))
			 /channelWidth);
      AliDebug(10,Form("summult mv %i   mult  in chammens %i in ps %f ", 
		      sumMult, fSumMult, fSumMult*channelWidth));
    }

    fT0->AddDigit(bestATDC,bestCTDC,meanTime,timeDiff,fSumMult, refpoint,
		       ftimeCFD,fADC0,ftimeLED,fADC);
     
    
      AliDebug(10,Form(" Digits wrote refpoint %i bestATDC %i bestCTDC %i  meanTime %i  timeDiff %i fSumMult %i ",refpoint ,bestATDC,bestCTDC,meanTime,timeDiff,fSumMult ));
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
    fT0  = (AliT0*)outRL ->GetAliRun()->GetDetector("T0");
    // Make a branch in the tree 
    fT0->MakeBranch("D");
     treeD->Fill();
  
     pOutStartLoader->WriteDigits("OVERWRITE");
     
     fT0->ResetDigits();
     pOutStartLoader->UnloadDigits();
     
}
