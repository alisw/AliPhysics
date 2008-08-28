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




//...
//  Checks the quality assurance. 
//  By comparing with reference data
//  Skeleton for T0
//---------------------------------------------
//checkig without reference data:
//for RAW QA all histograms should have approximatly the same 
//number of entries as RefPoint
//for Rec Points checks 
//  - amplitude measured by 2 methos
//  - online and offline T0 measurements
// for ESD quality of reconstruction ( and measurements):
// RMS of vertex and T0 less than 75ps
//
// Alla.Maevskaya@cern.ch   
//...

// --- ROOT system ---
#include <Riostream.h>
#include <TClass.h>
#include <TH1F.h> 
#include <TH1I.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 
#include <TMath.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliT0QAChecker.h"

ClassImp(AliT0QAChecker)


//__________________________________________________________________
const Double_t AliT0QAChecker::Check(AliQA::ALITASK_t index,TObjArray * list)
{

  // Super-basic check on the QA histograms on the input list:
  // look whether they are empty!


  Double_t test = 10.0  ;
   
  Int_t count = 0 ;
  Double_t nent[200], nentraw[200];
  TString hname[200];
  const char *cname;
  memset(nent,0,200*sizeof(Double_t));
  Double_t w[200];
  memset(w,1,200*sizeof(Double_t));
  TH1 *fhRecLEDAmp[24];  TH1 * fhRecQTC[24];
  TH1 *fhOnlineMean = 0x0;  
  TH1 * fhRecMean = 0x0;
  TH1 *fhESDMean = 0x0;
  TH1 *fhESDVertex = 0x0;
  TString dataType = AliQA::GetAliTaskName(index);

  if (list->GetEntries() == 0){
    test = 1. ; // nothing to check
  }
  else {
    
    TIter next(list) ;
    TH1 * hdata ;
    
    count = 0 ;
    while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
      if (hdata) {
	nent[count] = hdata->GetEntries();
	cname = hdata->GetName();
	hname[count] = cname;
	AliDebug(10,Form("count %i %s -> %f",count, hname[count].Data(),nent[count])) ;
	if(index==2){
	  if(count>23 && count<48)fhRecLEDAmp[count-24] = hdata;
	  if(count>47 && count<72)fhRecQTC[count-48] = hdata;
	  if(count == 72)  fhOnlineMean = hdata; 
	  if(count == 73)  fhRecMean = hdata; 
	}
	

	if(index==3){
	  if(count==0) fhESDMean = hdata;
	  if(count==1) fhESDVertex = hdata;
	  if(count>1){
	    AliWarning("Unknowm ESD QA histograms");
	    return 0;
	  }
	}
	count++ ;
	
        Double_t rv = 0.;
        if(hdata->GetEntries()>0) rv = 1;
	//   AliInfo(Form("%s -> %f", hdata->GetName(), rv)) ;
        test += rv ;
	
     }
      else{
        AliError("Data type cannot be processed") ;
      }

    }
    
    if (count != 0) {
      if (test==0) {
        AliWarning("Histograms are there, but they are all empty: setting flag to kWARNING");
        test = 0.5;  //upper limit value to set kWARNING flag for a task
      }
      else {
	if(index == 2){
	  //rec points
	  if ( TMath::Abs(fhRecMean->GetMean() - fhOnlineMean->GetMean()) > 5) 
	    AliDebug(10,Form("rec mean %f -> online mean %f",fhRecMean->GetMean(), fhOnlineMean->GetMean())) ;
	  Double_t meanLED, meanQTC;
	  for (Int_t idet=0; idet<24; idet++) {
	    meanLED = fhRecLEDAmp[idet]->GetMean();
	    meanQTC = fhRecQTC[idet]->GetMean();
	    if (TMath::Abs(meanLED-meanQTC)> 1.) 
	      AliDebug(1,Form("Amplitude measurements are different in channel %i : Amp LED %f -> Amp QTC %f",idet,meanLED, meanQTC)) ;
	  }
	}	
	 
	
	if (index == 0) {
	  //raw data
	  Float_t realNumber = Float_t(nent[0]);
	  for (Int_t i=77; i<count; i++) 
	    {
	      Double_t diff = TMath::Abs(nent[i]-realNumber);
	      if (diff > 0.1*realNumber )
		AliDebug(1,Form("Problem in Number of entried in hist %s  is %f number of RefPoints %f\n",hname[i].Data() , nent[i],realNumber )) ; 
	    }
	}
 	if (index == 3) {
	  //ESD
	  Double_t rmsMeanTime = fhESDMean->GetRMS();
	  if (rmsMeanTime>3) 		
	    AliDebug(1,Form("Mean time with bad resolution, RMS= %f",rmsMeanTime)) ; 
	  Double_t rmsVertex = fhESDVertex->GetRMS();
	  if (rmsVertex>3) 		
	    AliDebug(1,Form("Vertex with bad resolution, RMS= %f",rmsVertex)) ; 
	}
 
      }
      
    }


  } //  if (list->GetEntries() != 0
  
  AliInfo(Form("Test Result = %f", test)) ;

  return test ;
}

