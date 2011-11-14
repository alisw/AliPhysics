
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
#include <TSpectrum.h>
#include <TProfile.h>
#include <iostream>

ClassImp(AliT0CalibTimeEq)

//________________________________________________________________
  AliT0CalibTimeEq::AliT0CalibTimeEq():TNamed(),
				       fMeanVertex(0),        
				       fRmsVertex(0)      
{
  //
  for(Int_t i=0; i<24; i++) {
   fTimeEq[i] = 0;	      // Time Equalized for OCDB	 
   fTimeEqRms[24] = -1;	      // RMS of Time Equalized for OCDB	 
   for (Int_t ih=0; ih<5; ih++)   fCFDvalue[i][ih] = 0;
  }
}

//________________________________________________________________
AliT0CalibTimeEq::AliT0CalibTimeEq(const char* name):TNamed(),
				       fMeanVertex(0),        
				       fRmsVertex(0)      
{
  //constructor

  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  for(Int_t i=0; i<24; i++) {
   fTimeEq[i] = 0;	      // Time Equalized for OCDB	 
   fTimeEqRms[24] = -1;	      // RMS of Time Equalized for OCDB	 
   for (Int_t ih=0; ih<5; ih++)   fCFDvalue[i][ih] = 0;
  }
}

//________________________________________________________________
AliT0CalibTimeEq::AliT0CalibTimeEq(const AliT0CalibTimeEq& calibda):TNamed(calibda),		
				       fMeanVertex(0),        
				       fRmsVertex(0)      
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
  for (Int_t i=0; i<24; i++) 
    printf(" CFD  %f diff %f  \n",fCFDvalue[i][0],fTimeEq[i]);
} 


//________________________________________________________________
Bool_t AliT0CalibTimeEq::ComputeOnlineParams(const char* filePhys)
{
  // compute online equalized time
  Float_t meandiff, sigmadiff, meanver, meancfdtime, sigmacfdtime;
  meandiff = sigmadiff =  meanver = meancfdtime = sigmacfdtime =0;
  Double_t rmsver=0;
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
      for (Int_t i=0; i<24; i++)
	{
	  meandiff = sigmadiff =  meanver = meancfdtime = sigmacfdtime =0;
	  TH1F *cfd     = (TH1F*) gFile->Get(Form("CFD1minCFD%d",i+1));
	  TH1F *cfdtime = (TH1F*) gFile->Get(Form("CFD%d",i+1));
	  if(!cfd) {
	    AliWarning(Form("no Diff histograms collected by PHYS DA for channel %i", i));
	    okdiff++;
	    if(okdiff<4) {
	      meandiff = 0;
	      sigmadiff = 0;
	    }
	    else
	      ok = false; 
	  }
	  if(!cfdtime) {
	    AliWarning(Form("no CFD histograms collected by PHYS DA for channel %i", i));
	    oktime++;
	    if(oktime<4) {
	      meancfdtime = 0;
	    }
	    else
	      ok = false; 
	  }
	  if(cfd) {
	    nent = Int_t(cfd->GetEntries());
	    if( nent<=50) {
	      okdiff++;
	      //	      printf(" pmt %i nent %i cfd->GetRMS() %f cfd->GetMean() %f \n",
	      //     i, nent, cfd->GetRMS(), cfd->GetMean() );
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
	      if(cfd->GetRMS()>1.5 )
		GetMeanAndSigma(cfd, meandiff, sigmadiff);
	      if(cfd->GetRMS()<=1.5) 
		{
		  meandiff = cfd->GetMean();
		  sigmadiff=cfd->GetRMS();
		}
	      Int_t   maxBin = cfd->GetMaximumBin(); 
	      Double_t  meanEstimate = cfd->GetBinCenter( maxBin); 
	      if(TMath::Abs(meanEstimate - meandiff) > 20 ) meandiff = meanEstimate; 
	    }
	  }
	  
	  if(cfdtime) {
	    nent = Int_t(cfdtime->GetEntries());
	    if( nent<=50 ) {
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
	    if(nent > 50  )  { //!!!!!!!!!!
	      if(cfdtime->GetRMS()>1.5 )
		GetMeanAndSigma(cfdtime,meancfdtime, sigmacfdtime);
	      if(cfdtime->GetRMS()<=1.5) 
		{
		  if(cfdtime->GetRMS()==0 ||cfdtime->GetMean()==0 ) {
		    ok = false;
		  }
		  meancfdtime = cfdtime->GetMean();
		  sigmacfdtime = cfdtime->GetRMS();
		}
	      Int_t   maxBin = cfdtime->GetMaximumBin(); 
	      Double_t  meanEstimate = cfdtime->GetBinCenter( maxBin); 
	      if(TMath::Abs(meanEstimate - meancfdtime) > 20 ) meancfdtime = meanEstimate; 
	    }
	  }
	  SetTimeEq(i,meandiff);
	  SetTimeEqRms(i,sigmadiff);
	  SetCFDvalue(i,0,meancfdtime);
	  if (cfd) delete cfd;
	  if (cfdtime) delete cfdtime;
	}
      TH1F *ver = (TH1F*) gFile->Get("hVertex") ;
      if(ver) {
	meanver = ver->GetMean();
	rmsver = ver->GetRMS();
      }
      SetMeanVertex(meanver);
      SetRmsVertex(rmsver);
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
  // 1-24  no histograms with time diff for more than okdiff channenls
  // 100 >3 histos are empty
  // - 1 if less 3 channels are empty ;for these channels OCDB value will be written back WARNING
  // -2  low statistics in few channels; for these channels OCDB value will be written back WARNING
  // for cfd
  // 200 >3 histos are empty 
  // -11 if less 3 channels are empty ;for these channels OCDB value will be written back WARNING
  // -22  low statistics in few channels; for these channels OCDB value will be written back WARNING
   //
  // compute offline equalized time
  Float_t meandiff, sigmadiff, meanver, meancfdtime, sigmacfdtime;
  meandiff = sigmadiff =  meanver = meancfdtime = sigmacfdtime =0;
  Int_t nent=0;
  Int_t ok = 0;
  Int_t okdiff=0;
  Int_t okcfd=0;
  TH1F *cfddiff = NULL; 
  TH1F *cfdtime = NULL;
  TObjArray * tzeroObj = NULL;

  gFile = TFile::Open(filePhys);
  if(!gFile) {
    AliError("No input PHYS data found ");
    ok = 1000;
    return ok;
  }
  else
    {
      meandiff = sigmadiff =  meanver = meancfdtime = sigmacfdtime =0;
      //      TDirectory *dr = (TDirectory*) gFile->Get("T0Calib");
      tzeroObj = dynamic_cast<TObjArray*>(gFile->Get("T0Calib"));
      for (Int_t i=0; i<24; i++)
	{
	  if (i != badpmt) {	    
	    if(tzeroObj) {
	      cfddiff = (TH1F*) tzeroObj->FindObject(Form("CFD1minCFD%d",i+1));
	      cfdtime = (TH1F*)tzeroObj->FindObject(Form("CFD%d",i+1));
	    }
	    else
	      {
		cfddiff = (TH1F*)gFile->Get(Form("CFD1minCFD%d",i+1));
		cfdtime = (TH1F*)gFile->Get(Form("CFD%d",i+1));
		
	      }
	    if(!cfddiff ) {
	      AliWarning(Form("no  histograms collected by pass0 for diff channel %i", i));    
	      okdiff++;
	      if(okdiff<4) {
		meandiff = timecdb[i];
		sigmadiff = 0;
		ok = -1;
	      }
	      else
		return 100+okdiff; 
	    }
	    if(cfddiff) {
	      nent = Int_t(cfddiff->GetEntries());
	      if ( nent == 0) {
		okdiff++;
		if(okdiff<4) {
		  meandiff = timecdb[i];
		  sigmadiff = 0;
		  ok = - 1;
		}
		else
		  return 100; 
	      }	      
	      if(nent<100 && nent>0) { //!!!!!
		AliWarning(Form(" Not  enouph data in PMT %i- PMT1:  %i ", i, nent));
		ok=-2;
		meandiff = timecdb[i];
		sigmadiff = 0; 
	      }
	      if(nent>100 )  { //!!!!!
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
	      if(okcfd<4) {
		meancfdtime = cfdvalue[i];
		ok = -11;
	      }
	      else {
		AliError(Form("no  histograms collected by pass0 for time %i channels ", okcfd));
		return 200; 
	      }
	    }
	    if(cfdtime) {
	      nent = Int_t(cfdtime->GetEntries());
	      if (nent == 0 || 
		  (cfdtime->GetRMS() == 0 || cfdtime->GetMean() ==0 ) ) 
		{
		  okcfd++;
		  if(okcfd<4) {
		    meancfdtime = cfdvalue[i];
		    ok = -11;
		    printf("!!!!bad data:: pmt %i nent%i RMS %f mean %f cdbtime %f \n",
			   i, nent, cfdtime->GetRMS(), cfdtime->GetMean(),  fCFDvalue[i][0]);
		  }
		  else
		    {
		    printf("!!!!fatal data:: pmt %i nent%i RMS %f mean %f cdbtime %f \n",
			   i, nent, cfdtime->GetRMS(), cfdtime->GetMean(),  fCFDvalue[i][0]);
		      AliError(Form(" histograms collected by pass0 for time %i channels are empty", okcfd));
		      return 200; 
		    }
		}
	      
		  
	      if(nent<100 && nent>0 ) 
		{
		AliWarning(Form(" Not  enouph data in PMT in CFD peak %i - %i ", i, nent));
		if(okcfd<4) {
		  meancfdtime = cfdvalue[i];
		  ok = -22;
		  printf("!!!!low statstics:: pmt %i nent%i RMS %f mean %f cdbtime %f \n",
			 i, nent, cfdtime->GetRMS(), cfdtime->GetMean(),  fCFDvalue[i][0]);
		}
		else {
		  printf("!!!!low fatal statstics:: pmt %i nent%i RMS %f mean %f cdbtime %f \n",
			 i, nent, cfdtime->GetRMS(), cfdtime->GetMean(),  fCFDvalue[i][0]);

		  AliError(Form(" Not  enouph data in PMT in CFD peak %i in channels ", okcfd));
		  return 200; 
		}

	      }

	      if( nent>100 && cfdtime->GetRMS() != 0 )    { //!!!!!
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
	  printf(" pmt %i diff %f sigma %f meancfdtime %f cdbtime %f \n",i, meandiff, sigmadiff, meancfdtime, fCFDvalue[i][0]);
	  if (cfddiff) cfddiff->Reset();
	  if (cfdtime) cfdtime->Reset();
	  } //bad pmt
	}      
      gFile->Close();
      delete gFile;
    }
    return ok; 
   }

//________________________________________________________________________
void AliT0CalibTimeEq::GetMeanAndSigma(TH1F* hist,  Float_t &mean, Float_t &sigma) {
  
  const double window = 2.;  //fit window 
  
  double meanEstimate, sigmaEstimate; 
  int maxBin;
  maxBin        =  hist->GetMaximumBin(); //position of maximum
  meanEstimate  =  hist->GetBinCenter( maxBin); // mean of gaussian sitting in maximum
  // sigmaEstimate = hist->GetRMS();
  sigmaEstimate = 10;
  TF1* fit= new TF1("fit","gaus", meanEstimate - window*sigmaEstimate, meanEstimate + window*sigmaEstimate);
  fit->SetParameters(hist->GetBinContent(maxBin), meanEstimate, sigmaEstimate);
  hist->Fit("fit","RQ","Q");

  mean  = (Float_t) fit->GetParameter(1);
  sigma = (Float_t) fit->GetParameter(2);

  delete fit;
}


