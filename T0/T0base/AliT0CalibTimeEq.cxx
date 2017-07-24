
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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for T0 calibration                       TM-AC-AM_6-02-2006  
// equalize time shift for each time CFD channel
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliT0CalibTimeEq.h"
#include "AliLog.h"
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TH2.h>
#include <TSpectrum.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TFitResult.h>
#include <iostream>

ClassImp(AliT0CalibTimeEq)

//________________________________________________________________
AliT0CalibTimeEq::AliT0CalibTimeEq():TNamed(),
  fMeanVertex(0),        
  fRmsVertex(0),
  fCFDvsTime(24)
{
  //
  for(Int_t i=0; i<24; i++) {
   fTimeEq[i] = 0;	      // Time Equalized for OCDB	 
   fTimeEqRms[i] = -1;	      // RMS of Time Equalized for OCDB	 
   //fCFDvsTime[i] = NULL;
   for (Int_t ih=0; ih<5; ih++)   fCFDvalue[i][ih] = 0;
  }
  fCFDvsTime.SetOwner(kTRUE);
}

//________________________________________________________________
AliT0CalibTimeEq::AliT0CalibTimeEq(const char* name):TNamed(),
						     fMeanVertex(0),        
						     fRmsVertex(0),
						     fCFDvsTime(24)
{
  //constructor

  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  for(Int_t i=0; i<24; i++) {
   fTimeEq[i] = 0;	      // Time Equalized for OCDB	 
   fTimeEqRms[i] = -1;	      // RMS of Time Equalized for OCDB	 
   //fCFDvsTime[i]=NULL;
   for (Int_t ih=0; ih<5; ih++)   fCFDvalue[i][ih] = 0;
  }
  fCFDvsTime.SetOwner(kTRUE);
}

//________________________________________________________________
AliT0CalibTimeEq::AliT0CalibTimeEq(const AliT0CalibTimeEq& calibda):TNamed(calibda),		
								    fMeanVertex(0),        
								    fRmsVertex(0),
								    fCFDvsTime(24)
{
// copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  ((AliT0CalibTimeEq &) calibda).Copy(*this);

}

//________________________________________________________________
AliT0CalibTimeEq &AliT0CalibTimeEq::operator =(const AliT0CalibTimeEq& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
 
  if (this != &calibda) (( AliT0CalibTimeEq &) calibda).Copy(*this);
  return *this;
}

//________________________________________________________________
AliT0CalibTimeEq::~AliT0CalibTimeEq()
{
  //
  // destrictor
}
//________________________________________________________________
void AliT0CalibTimeEq::Reset()
{
  //reset values

  memset(fCFDvalue,0,120*sizeof(Float_t));
  memset(fTimeEq,1,24*sizeof(Float_t));
}


//________________________________________________________________
void  AliT0CalibTimeEq::Print(Option_t*) const
{
  // print time values

  printf("\n	----	PM Arrays	----\n\n");
  printf(" Time delay CFD \n");
  for (Int_t i=0; i<24; i++) {
    printf(" CFD  %f diff %f qt1 %f pedestal %f  \n",fCFDvalue[i][0],fTimeEq[i],fCFDvalue[i][1], fCFDvalue[i][3]);
  printf(" TVDC  %f OrA %f OrC %f  \n",fMeanVertex,fCFDvalue[0][2], fCFDvalue[1][2]);
  TGraph *gr = (TGraph*)fCFDvsTime.At(i);
  gr->Print();
  }
} 


//________________________________________________________________
Bool_t AliT0CalibTimeEq::ComputeOnlineParams(const char* filePhys)
{
  // compute online equalized time
  Float_t meandiff, sigmadiff, meanver, meancfdtime, sigmacfdtime;
  Float_t ora,orc, meanqt, sigmaor,sigmaver, meanEstimate;
  Float_t meanpedold=0, sigmaped=0;

  meandiff = sigmadiff =  meanver = meancfdtime = sigmacfdtime = ora = orc = meanqt =0;
  meanEstimate=sigmaver=sigmaor = 0;
  Int_t maxBin=0;
  Int_t nent=0;
  Int_t okdiff=0;
  Int_t oktime=0;
  Bool_t ok=false;
  //
  gFile = TFile::Open(filePhys);
  if(!gFile) {
    AliError("No input PHYS data found ");
  }
  else
    {
      ok=true;
      printf(" OK %i\n",ok);
      gFile->ls();
      for (Int_t i=0; i<24; i++)  // PMT cycle
	{
	  meandiff = sigmadiff =  meanver = meancfdtime = sigmacfdtime =0;
	  TH1F *hcfd     = (TH1F*) gFile->Get(Form("CFD1minCFD%d",i+1));
	  TH1F *hcfdtime = (TH1F*) gFile->Get(Form("CFD%d",i+1));
	  TH1F *hqt1 = (TH1F*) gFile->Get(Form("QT1%d",i+1));
	  TH1F *hPedOld = (TH1F*) gFile->Get(Form("hPed%i",i+1));
	  TH2F *hCFDvsTime = (TH2F*)  gFile->Get(Form("hCFDvsTimestamp%i_0",i+1) );
	  Float_t meanCFDvsTime[60], timestamp[60];
	  TGraph *cCFDvsTime;
	  for (int ibin=0; ibin<60; ibin++) meanCFDvsTime[ibin]=0;
	  if(!hcfd) {
	    AliWarning(Form("no Diff histograms collected by PHYS DA for channel %i", i));
	    okdiff++;
	    if(okdiff<4) {
	      meandiff = 0;
	      sigmadiff = 0;
	    }
	    else
	      ok = false; 
	  }
	  if(!hcfdtime) {
	    AliWarning(Form("no CFD histograms collected by PHYS DA for channel %i", i));
	    oktime++;
	    if(oktime<4) {
	      meancfdtime = 0;
	    }
	    else
	      ok = false; 
	  }
	  if(hcfd) {
	    nent = Int_t(hcfd->GetEntries());
	    if( nent<=50) {
	      okdiff++;
	      //  printf(" pmt %i nent %i cfd->GetRMS() %f cfd->GetMean() %f \n",
	      //     i, nent, hcfd->GetRMS(), hcfd->GetMean() );
	      if(okdiff<4) {
		meandiff = 0;
		sigmadiff = 0;
	      }
	      else
		{
		  // printf(" OK fsle:: pmt %i nent %i cfd->GetRMS() %f cfd->GetMean() %f \n",
		  // i, nent, cfd->GetRMS(), cfd->GetMean() );
		  AliWarning(Form(" Not  enouph data in PMT %i- PMT1:  %i ", i, nent));
		  ok = false; 
		}
	    }
	    if(nent>50)  { 
	      if(hcfd->GetRMS()>1.5 )
		GetMeanAndSigma(hcfd, meandiff, sigmadiff);
	      if(hcfd->GetRMS()<=1.5) 
		{
		  meandiff = hcfd->GetMean();
		  sigmadiff=hcfd->GetRMS();
		}
	      maxBin = hcfd->GetMaximumBin(); 
	      meanEstimate = hcfd->GetBinCenter( maxBin); 
	      if(TMath::Abs(meanEstimate - meandiff) > 20 ) meandiff = meanEstimate; 
	    }
	  }
	  
	  if(hcfdtime) {
	    nent = Int_t(hcfdtime->GetEntries());
	    if( nent<=50  ) {
	      oktime++;
	      if(oktime<4) {
		meancfdtime = 0;
	      }
	      else
		{
		  AliWarning(Form(" Not  enouph data in PMT %i:  %i ", i, nent));
		  ok = false; 
		}
	    }
	    if(nent > 500  )  { //!!!!!!!!!!
	      if(hcfdtime->GetRMS()>1.5 ) {
		//	TFitResultPtr r = hcfdtime->Fit("gaus","SQ");
		//  if ((Int_t)r==0) 
		//	meancfdtime = r->Parameters()[1];
		GetMeanAndSigma(hcfdtime,meancfdtime, sigmacfdtime);
	      }
	      if(hcfdtime->GetRMS()<=1.5) 
		{
		  if(hcfdtime->GetRMS()==0 ||hcfdtime->GetMean()==0 ) {
		    ok = false;
		  }
		  meancfdtime = hcfdtime->GetMean();
		  sigmacfdtime = hcfdtime->GetRMS();
		}
	      
	      printf ("PMT %i meandiff %f meancfdtime %f meanhist %f entries %f\n",i,meandiff,meancfdtime, hcfdtime->GetMean(), hcfdtime->GetEntries());
	    }
	    if(hqt1) 	GetMeanAndSigma(hqt1,meanqt, sigmaor);
	    //Pedestals
	    if(hPedOld )  meanpedold =hPedOld->GetMean();
	  } //cycle 24 PMT
	  if (hCFDvsTime) {
	    // CFD vs timestamp
	    for (int ibin=1; ibin<59; ibin++) {
	      timestamp[ibin-1] =  hCFDvsTime->GetXaxis()->GetBinCenter(ibin);
	      TH1D* projd = hCFDvsTime->ProjectionY(Form("prY%i_%i",ibin,i),ibin, ibin+1);
	      TH1F* proj = (TH1F*)projd;
	      if(proj->GetEntries()>200 ) {
		Float_t sigmat;
		GetMeanAndSigma(proj, meanCFDvsTime[ibin-1], sigmat);
		//	TFitResultPtr r = proj->Fit("gaus","SQ");
		//	if ((Int_t)r==0) meanCFDvsTime[ibin-1] = r->Parameters()[1];
	      }
	    
	      else 
		meanCFDvsTime[ibin-1] = meancfdtime; 
	    }
	    cCFDvsTime = new TGraph(58, timestamp, meanCFDvsTime);
	    // cCFDvsTime->Print();
	    fCFDvsTime.AddAtAndExpand(cCFDvsTime,i);
	    //  cCFDvsTime->Delete();
	  }
	  SetTimeEq(i,meandiff);
	  SetTimeEqRms(i,sigmadiff);
	  SetCFDvalue(i,0,meancfdtime);
	  SetPedestalOld(i,meanpedold);       
	  SetCFDvalue(i,1,meanqt);
	    
	  if (hcfd) delete hcfd;
	  if (hcfdtime) delete hcfdtime;
	  if(hqt1) delete hqt1;
	  if(hPedOld) delete hPedOld;
	  //	  if(cCFDvsTime) delete cCFDvsTime;
	}
      TH1F *hver = (TH1F*) gFile->Get("hVertex") ;
      TH1F *hora = (TH1F*) gFile->Get("hOrA") ;
      TH1F *horc = (TH1F*) gFile->Get("hOrC") ;
      if(hver) 	GetMeanAndSigma(hver,meanver, sigmaver);
      if(hora) 	GetMeanAndSigma(hora,ora, sigmaor);
      if(horc)  GetMeanAndSigma(horc,orc, sigmaor);
      SetMeanVertex(meanver);
      SetRmsVertex(sigmaver);
      SetOrA(ora);
      SetOrC(orc);

  
      gFile->Close();
      delete gFile;
    }
    return ok; 
}

//________________________________________________________________
Int_t AliT0CalibTimeEq::ComputeOfflineParams(const char* filePhys, Float_t *timecdb, Float_t *cfdvalue, Int_t badpmt)
{
  // fStatus implementation:
  //ok means
  // 1000  - can not open file 
  // for timediff
  // 20 >3 histos are empty
  // -11  low statistics oe empty histos in few channels; for these channels OCDB value will be written back WARNING
  // for cfd
  // 20 >2 histos are empty or with low statistic 
  // -11 if less 3 channels are empty or low statistics in few channels ;for these channels OCDB value will be written back WARNING
   //
  // compute offline equalized time
  Float_t meandiff, sigmadiff, meanver, meancfdtime, sigmacfdtime;
  meandiff = sigmadiff =  meanver = meancfdtime = sigmacfdtime =0;
  Int_t nent=0;
  Int_t ok = 0;
  Int_t okcfd=0;
  TH1F *cfddiff = NULL; 
  TH1F *cfdtime = NULL;
  TH2F *hCFDvsTime = NULL;
  TObjArray * tzeroObj = NULL;
  Float_t qt1[24], ped[24], orA, orC, tvdc;
  Float_t meanCFDvsTime[121], timestamp[121];
  TGraph *cCFDvsTime;
  for (int ibin=0; ibin<121; ibin++) meanCFDvsTime[ibin]=0;
  
  gFile = TFile::Open(filePhys);
  if(!gFile) {
    AliError("No input PHYS data found ");
    ok = 1000;
    return ok;
  }
  else {
      meandiff = sigmadiff =  meanver = meancfdtime = sigmacfdtime =0;
      //      TDirectory *dr = (TDirectory*) gFile->Get("T0Calib");
      tzeroObj = dynamic_cast<TObjArray*>(gFile->Get("T0Calib"));
      for (Int_t i=0; i<24; i++)
	{
	  if (i != badpmt) {	    
	    if(tzeroObj) {
	      cfddiff = (TH1F*) tzeroObj->FindObject(Form("CFD1minCFD%d",i+1));
	      cfdtime = (TH1F*) tzeroObj->FindObject(Form("CFD%d",i+1));
	      hCFDvsTime = (TH2F*) tzeroObj->FindObject(Form("hCFDvsTimestamp%i",i+1) );
	    }
	    else
	      {
		cfddiff = (TH1F*)gFile->Get(Form("CFD1minCFD%d",i+1));
		cfdtime = (TH1F*)gFile->Get(Form("CFD%d",i+1));
		
	      }
	    if(!cfddiff ) {
	      AliWarning(Form("no  histograms collected by pass0 for diff channel %i\n", i));    
		meandiff = timecdb[i];
		sigmadiff = 0;
	      }
	    if(cfddiff) {
	      nent = Int_t(cfddiff->GetEntries());
	      if ( nent == 0  ) {
		AliWarning(Form("no  entries in histogram for diff channel %i\n", i));
		meandiff = timecdb[i];
		sigmadiff = 0;
	      }
	      if(nent<=100 && nent>0) { //!!!!!
		AliWarning(Form(" Not  enouph data in PMT %i- PMT1:  %i \n", i, nent));
		meandiff = timecdb[i];
		sigmadiff = 0; 
	      }
	      if(nent>=100 )  { //!!!!!
		{
		  if(cfddiff->GetRMS()>1.5 )
		    GetMeanAndSigma(cfddiff, meandiff, sigmadiff);
		  if(cfddiff->GetRMS()<=1.5) 
		    {
		      meandiff = cfddiff->GetMean();
		      sigmadiff = cfddiff->GetRMS();
		    }
		  Int_t   maxBin = cfddiff->GetMaximumBin(); 
		  Double_t  meanEstimate = cfddiff->GetBinCenter( maxBin); 
		  if(TMath::Abs(meanEstimate - meandiff) > 20 ) meandiff = meanEstimate; 
		}
	      }	    
	    }
	    if(!cfdtime ) {
	      AliWarning(Form("no  histograms collected by pass0 for time channel %i", i));
	      meancfdtime = cfdvalue[i];
	      okcfd++;
	      if(okcfd<2) {
		meancfdtime = cfdvalue[i];
		//	ok = -11;
	      }
	      else {
		AliError(Form("no  histograms collected by pass0 for time %i channels ", okcfd));
		if (tzeroObj) delete tzeroObj;
		return 20; 
	      }
	    }
	    if(cfdtime) {
	      nent = Int_t(cfdtime->GetEntries());
	      if (nent == 0 || 
		  (cfdtime->GetRMS() == 0 || cfdtime->GetMean() ==0 ) ) 
		{
		  okcfd++;
		  if(okcfd<2) {
		    meancfdtime = cfdvalue[i];
		    //   ok = -11;
		    printf("!!!!bad data:: pmt %i nent%i RMS %f mean %f cdbtime %f \n",
			   i, nent, cfdtime->GetRMS(), cfdtime->GetMean(), cfdvalue[i] );
		  }
		  else
		    {
		    printf("!!!!fatal data:: pmt %i nent%i RMS %f mean %f cdbtime %f \n",
			   i, nent, cfdtime->GetRMS(), cfdtime->GetMean(), cfdvalue[i]);
		      AliError(Form(" histograms collected by pass0 for time %i channels are empty", okcfd));
		      if (tzeroObj) delete tzeroObj;
		      return 20; 
		    }
		}
	      
		  
	      if(nent<=100 && nent>0 ) 
		{
		  okcfd++;
		  AliWarning(Form(" Not  enouph data in PMT in CFD peak %i - %i ", i, nent));
		  meancfdtime = cfdvalue[i];
		  //		  ok = -11;
		  printf("!!!!low statstics:: pmt %i nent%i RMS %f mean %f cdbtime %f \n",
			 i, nent, cfdtime->GetRMS(), cfdtime->GetMean(),  cfdvalue[i]);
		  if (okcfd>2) {
		    ok = -11;
		    if (tzeroObj) delete tzeroObj;
		    return ok;
		  }
		}
	      
	      if( nent>100 )    { //!!!!!
		if(cfdtime->GetRMS()>1.5 )
		  GetMeanAndSigma(cfdtime,meancfdtime, sigmacfdtime);
		if(cfdtime->GetRMS()<=1.5) 
		  {
		    meancfdtime = cfdtime->GetMean();
		    sigmacfdtime=cfdtime->GetRMS();
		  }
		Int_t   maxBin = cfdtime->GetMaximumBin(); 
		Double_t  meanEstimate = cfdtime->GetBinCenter( maxBin); 
		if(TMath::Abs(meanEstimate - meancfdtime) > 20 ) meancfdtime = meanEstimate; 
	      }
	    }
	  
	    SetTimeEq(i,meandiff);
	    SetTimeEqRms(i,sigmadiff);
	    SetCFDvalue(i,0, meancfdtime );
	    qt1[i]=cfdvalue[24+i];
	    SetCFDvalue(i,1,qt1[i]);
	    ped[i]=cfdvalue[52+i];
	    SetCFDvalue(i,3,ped[i]);
	    //    AliInfo(Form(" !!! AliT0CalibTimeEq pmt %i pedestal %f \n ", i, ped[i]) );
	    if (cfddiff) cfddiff->Reset();
	    if (cfdtime) cfdtime->Reset();
	  if (hCFDvsTime) {
	    // CFD vs timestamp
	    for (int ibin=1; ibin<120; ibin++) {
	      timestamp[ibin-1] =  hCFDvsTime->GetXaxis()->GetBinCenter(ibin);
	      TH1D* projd = hCFDvsTime->ProjectionY(Form("prY%i_%i",ibin,i),ibin, ibin+1);
	      TH1F* proj = (TH1F*)projd;
	      if(proj->GetEntries()>200 ) {
		Float_t sigmat;
		GetMeanAndSigma(proj, meanCFDvsTime[ibin-1], sigmat);
	      }
	      else 
		meanCFDvsTime[ibin-1] = meancfdtime; 
	    }
	    cCFDvsTime = new TGraph(58, timestamp, meanCFDvsTime);
	    // cCFDvsTime->Print();
	    fCFDvsTime.AddAtAndExpand(cCFDvsTime,i);
	    //  cCFDvsTime->Delete();
	  }

	    
	  } //bad pmt
	}      
      SetMeanVertex(cfdvalue[48]);
      SetOrA(cfdvalue[49]);
      SetOrC(cfdvalue[50]);
      
      gFile->Close();
      delete gFile;
    }
  if (tzeroObj) delete tzeroObj;
  return ok; 
}

//________________________________________________________________________
void AliT0CalibTimeEq::GetMeanAndSigma(TH1F* hist,  Float_t &mean, Float_t &sigma) {
  
  const double window = 3.;  //fit window 
  double meanEstimate, sigmaEstimate; 
  int maxBin;
  maxBin        =  hist->GetMaximumBin(); //position of maximum
  meanEstimate  =  hist->GetBinCenter( maxBin); // mean of gaussian sitting in maximum
  sigmaEstimate = hist->GetRMS();
  // sigmaEstimate = 10;
  TF1* fit= new TF1("fit","gaus", meanEstimate - window*sigmaEstimate, meanEstimate + window*sigmaEstimate);
  fit->SetParameters(hist->GetBinContent(maxBin), meanEstimate, sigmaEstimate);
  hist->Fit("fit","RQ","");

  mean  = (Float_t) fit->GetParameter(1);
  sigma = (Float_t) fit->GetParameter(2);

 if(TMath::Abs(meanEstimate - mean) > 20 ) mean = meanEstimate; 


  delete fit;
}

