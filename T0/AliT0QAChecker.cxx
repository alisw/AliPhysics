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
#include <TH2.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 
#include <TMath.h>
#include <TString.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliT0QAChecker.h"

ClassImp(AliT0QAChecker)


//____________________________________________________________________________
Double_t * AliT0QAChecker::Check(AliQAv1::ALITASK_t /*index*/)
{
  Double_t * rv = new Double_t[AliRecoParam::kNSpecies] ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    rv[specie] = 0.0 ; 
  return rv ;  
}

//__________________________________________________________________
Double_t * AliT0QAChecker::Check(AliQAv1::ALITASK_t index,TObjArray ** list)
{

  // Super-basic check on the QA histograms on the input list:
  // look whether they are empty!

  Double_t * test = new Double_t[AliRecoParam::kNSpecies] ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    test[specie]    = 10.0 ; 

  Double_t nent[250];
  TString hname[250];
  const char *cname;
  memset(nent,0,250*sizeof(Double_t));
  Double_t w[250];
  memset(w,1,250*sizeof(Double_t));
  TH2 *fhRecDiff[3];  
  TH2 *fhRawEff[250];
  TH1 *fhESD[2];
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    //  TString dataType = AliQAv1::GetAliTaskName(index);
    if (list[specie]->GetEntries() == 0){
      test[specie] = 1. ; // nothing to check
    }
    else {
    TIter next(list[specie]) ;
    TH1 * hdata ;
    TH2 * h ;
    //  AliInfo(Form(" data type %i %s nentries %i\n",
    //	   index,dataType.Data(),list->GetEntries()));
    for (Int_t ir=0; ir<list[specie]->GetEntries(); ir++) {
      //raw
      if(index==0 ){
        /*
         if(ir < 205) {
         hdata = (TH1*) list[specie]->UncheckedAt(ir);
         if(hdata) {
         cname = hdata->GetName();
         hname[ir] = cname;
         AliDebug(10,Form("count %i %s \n",ir, hname[ir].Data())) ;
         fhRaw[ir] = hdata;
         }
         }*/
        if(ir > 204) {
          //	  else{
          h = (TH2*) list[specie]->UncheckedAt(ir);
          AliInfo(Form(" index %i ir %i \n", index,ir));
          if(h) {
            cname = h->GetName();
            hname[ir] = cname;
            AliDebug(1,Form("count %i %s \n",ir, hname[ir].Data())) ;
            fhRawEff[ir] = h;
          }
        }
      }
     
      //rec
      if(index==2){
        h = (TH2*) list[specie]->UncheckedAt(ir);
        if(h) {
          cname = h->GetName();
          hname[ir] = cname;
          AliDebug(1,Form("count %i %s \n",ir, hname[ir].Data())) ;
          fhRecDiff[ir] = h;
        }
      }
      //esd
      if(index==3){
        cout<<" ir "<<ir<<endl;
        hdata = (TH1*) list[specie]->UncheckedAt(ir);
        if(hdata){
          fhESD[ir] = hdata;
          AliDebug(1,Form("count %i %s ",ir, hname[ir].Data()) );
        }
      }
    }
      if (index == 0) {
        //raw data
	
        for (Int_t icase=205; icase<207; icase++) {
          for (Int_t idet=0; idet<24; idet++) {
            Double_t mean = fhRawEff[icase]->ProjectionY(Form("%s_py_%i_%i",
                                                              fhRawEff[icase]->GetName(), idet,icase),
                                                              idet,idet+1)->GetMean();
            Double_t rms= fhRawEff[icase]->ProjectionY(Form("%s_py%i_%i", 
                                                            fhRawEff[icase]->GetName(), idet,icase),
                                                            idet,idet+1)->GetRMS();
            AliInfo(Form("name %s icase %i idet %i mean %f, rms %f\n",
                    fhRawEff[icase]->GetName(), icase, idet, mean,rms));
            if (mean<1.2 && mean> 0.8 ) {
              test[specie] = 1;
              AliDebug(1,Form("All channels works meane efficieny %f with rms %f test %f",  mean, rms, test[specie])) ; 
            }
            if (mean<=0.8 && mean>= 0.5 ){
              test[specie] = 0.5;
              AliDebug(1,Form("%s problem in channel %i  efficieny %f test %f",
                              fhRawEff[icase]->GetName(), idet,  mean, test[specie])) ; 
            }
            if (mean<0.5 ) { 
              test[specie] = 0.25;
              AliDebug(1,Form("%s big problem in channel %i  efficieny %f test %f",
                              fhRawEff[icase]->GetName(), idet,  mean, test[specie])) ; 
            }
          }
        }
      }
      if(index == 2){
        //rec points
        for (Int_t icase=0; icase<2; icase++) {
          for (Int_t idet=0; idet<24; idet++) {
            Double_t mean = fhRecDiff[icase]->
            ProjectionY(Form("%s_py", fhRecDiff[icase]->GetName()),
                        idet,idet+1)->GetMean();
            Double_t rms= fhRecDiff[icase]->
            ProjectionY(Form("%s_py", fhRecDiff[icase]->GetName()),
                        idet,idet+1)->GetRMS();
            AliInfo(Form("name %s icase %i idet %i mean %f, rms %f\n",
                         fhRecDiff[icase]->GetName(), icase, idet, mean,rms)); 
	  	  
            if(TMath::Abs(mean) >1.5 || rms >1){
              AliDebug(1,Form(" calibration is nor perfect; test=%f", test)) ;
              test[specie]=0.25;
            }
            if(mean>3 || rms >5) {
              test[specie] = 0.1;
              AliDebug(1,Form(" wrong calibration test=%f", test[specie])) ;
            } 
          }
        }	 
      }
      if (index == 3) {
        //ESD
        for (Int_t icase=0; icase<2; icase++) {
          Double_t rmsVertex = fhESD[icase]->GetRMS();
          Double_t meanVertex = fhESD[icase]->GetMean();
          test[specie]=1;
          AliInfo(Form("numentries %d meanVertex %f rmsVertex %f", fhESD[icase]->GetEntries(), meanVertex, rmsVertex));
          if (TMath::Abs(rmsVertex)>3) {
            test[specie]=0.25;
            AliDebug(1,Form("Vertex position resolution not good  , rms= %f test=%f",
                            rmsVertex, test[specie])) ; 
          }
          if (TMath::Abs(meanVertex)>3) {
            test[specie]=0.25;
            AliDebug(1,Form("Vertex position bad calibrated  , Mean= %f test=%f",
                            meanVertex, test[specie])) ; 
          }
        }
      }
    } //  if (list->GetEntries() != 0
    AliInfo(Form("Test Result = %f", test[specie])) ;
  } 
  return test ;
}

